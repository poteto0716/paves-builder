"""Stages of the polypaves adhesion workflow (driven by adhesion.py).

Each stage reads the previous stage's directory under <project>/runs/ and writes its
own: state.xml (OpenMM positions, velocities, box), final.data + in.styles (LAMMPS),
final.pdb, log.csv, summary.json. See README.md for the protocol.
"""
import argparse
import json
import math
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import openmm as mm
import openmm.unit as u

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mdtools as T  # noqa: E402

ATM_TO_BAR = 1.01325
AVOGADRO = 6.02214076e23


# --- plumbing ------------------------------------------------------------------

class Stage:
    def __init__(self, cfg, runs, name):
        self.cfg = cfg
        self.dir = runs / name
        self.dir.mkdir(parents=True, exist_ok=True)
        self.scale = float(cfg.get('step_scale', 1.0))
        self.t0 = time.time()

    def steps(self, n):
        return max(1, int(round(n * self.scale)))

    def say(self, *a):
        print(f'[{self.dir.name} {time.time() - self.t0:8.1f}s]', *a, flush=True)


def integrator(T_K, dt_fs, friction):
    return mm.LangevinMiddleIntegrator(T_K * u.kelvin, friction / u.picosecond, dt_fs * u.femtosecond)


def context_for(cfg, system, integ):
    return T.make_context(system, integ, platform=cfg['platform'], precision=cfg['precision'])


def kcal(e_kj):
    return e_kj / T.KCAL


def potential_kcal(context):
    return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)


def finish(stage, context, system, template, styles, *, top=None, extra=None):
    """state.xml, final.pdb, final.data + in.styles, summary.json."""
    st = T.save_state(context, stage.dir / 'state.xml')
    pos, box = T.positions_A(st), T.box_A(st)
    vel = T.velocities_A_per_fs(st)
    T.write_lammps_data(template, stage.dir / 'final.data', pos, box, velocities=vel,
                        title=f'{stage.dir.name} (OpenMM -> LAMMPS, paves-interface)')
    T.write_styles(styles, stage.dir / 'in.styles', 'final.data')
    if top is None:
        top, _ = T.topology(system)
    T.write_pdb(stage.dir / 'final.pdb', top, pos, box)
    summary = {'stage': stage.dir.name, 'atoms': len(pos), 'box_A': box.tolist(),
               'potential_kcal': potential_kcal(context), 'wall_seconds': time.time() - stage.t0}
    summary.update(extra or {})
    T.dump_json(stage.dir / 'summary.json', summary)
    stage.say('done', json.dumps({k: v for k, v in summary.items() if k != 'box_A'}, default=float))
    return st


def density_g_cm3(system, box_A):
    return T.masses(system).sum() / AVOGADRO / (np.prod(box_A) * 1e-24)


# --- 01 build: PAVES makes the polymer melt ------------------------------------

def silica_cell(cfg):
    """(Lx, Ly, Lz) of the substrate model, from the MOL2's CRYSIN record."""
    lines = Path(cfg['silica_mol2']).read_text().splitlines()
    k = next(i for i, l in enumerate(lines) if l.startswith('@<TRIPOS>CRYSIN'))
    return [float(v) for v in lines[k + 1].split()[:3]]


# Records the workflow sets itself; a user's .paves may not choose them.
OWNED_RECORDS = {'forcefield', 'cell', 'density', 'periodic', 'output', 'openmm_system', 'temperature',
                 'backend', 'platform', 'precision', 'region', 'soft_wall', 'relax', 'minimize'}
REFUSED_RECORDS = {'environment', 'structure', 'structure_format', 'structure_topology'}
CARRIED_RECORDS = {'name', 'seed', 'sequence_seed'}


def split_polypaves(text):
    """A user's .paves -> (chemistry lines kept, {carried record: value}, [owned records dropped]).

    The first word of a line is its record. Lines inside a `composition` ...
    `end` block are kept as they are.
    """
    kept, carried, dropped = [], {}, []
    in_block = False
    for raw in text.splitlines():
        line = raw.rstrip()
        words = line.split('#')[0].split()
        if not words:
            continue
        key = words[0]
        if in_block:
            kept.append(line)
            in_block = key != 'end'
            continue
        if key == 'polypaves-build':
            continue
        if key in REFUSED_RECORDS:
            raise ValueError(f"'{key}' cannot be used here: the workflow builds the polymer from its chemistry "
                             'and adds the silica itself')
        if key in CARRIED_RECORDS:
            carried[key] = ' '.join(words[1:])
        elif key in OWNED_RECORDS:
            dropped.append(key)
        else:
            kept.append(line)
            in_block = key == 'composition' and len(words) == 1
    if not any(l.split()[0] in ('monomer', 'component') for l in kept):
        raise ValueError('no monomer or component record: the file states no chemistry')
    return kept, carried, dropped


def polymer_input(cfg):
    """The .paves text of the build stage, from project.json."""
    y = cfg['system']
    cell = silica_cell(cfg)
    nx, ny, _ = cfg['assemble']['supercell']
    Lx, Ly = cell[0] * nx, cell[1] * ny
    if 'polypaves_lines' in y:
        chemistry = '\n'.join(y['polypaves_lines'])
    else:
        chemistry = (f"monomer       A '{y['monomer']}'\n"
                     f"terminator    '{y.get('terminator', '*C')}'\n"
                     f"degree        {y['dp']}\n"
                     f"chains        {y['chains']}")
    text = f"""polypaves-build 1
# Written by adhesion.py from project.json. Edit project.json, not this file.
# The cell's x and y are the {nx} x {ny} substrate supercell, so the melt fits it.
name          {y.get('name', cfg['name'])}
seed          {y.get('seed', 2026)}
sequence_seed {y.get('sequence_seed', 7)}
forcefield    {cfg['forcefield']}
{chemistry}
cell          {Lx:.6f} {Ly:.6f} density
density       {y.get('build_density', 1.0)}
temperature   {cfg['bulk']['stages'][0][1]}
openmm_system yes
output        . polymer
"""
    return text, (Lx, Ly)


def stage_build(cfg, runs):
    s = Stage(cfg, runs, '01_build')
    y = cfg['system']
    text, (Lx, Ly) = polymer_input(cfg)
    (s.dir / 'polymer.paves').write_text(text)
    s.say(f'building {cfg["name"]}: cell {Lx:.4f} x {Ly:.4f} A, z from {y.get("build_density", 1.0)} g/cm3')
    polypaves_build(cfg, s.dir, 'polymer.paves')
    st = T.load_xml(s.dir / 'polymer.openmm_state.xml')
    info = {'stage': '01_build', 'atoms': int(len(T.positions_A(st))), 'box_A': T.box_A(st).tolist(),
            'wall_seconds': time.time() - s.t0}
    T.dump_json(s.dir / 'summary.json', info)
    s.say('done', json.dumps(info, default=float))


# --- 02 bulk: z-only NPT --------------------------------------------------------

def stage_bulk(cfg, runs):
    s = Stage(cfg, runs, '02_bulk')
    c = cfg['bulk']
    build = runs / '01_build' / 'polymer'
    system = T.load_xml(str(build) + '.openmm_system.xml')
    start = T.load_xml(str(build) + '.openmm_state.xml')
    added, _ = T.configure(system, shake=c['shake'])
    s.say(f'{system.getNumParticles()} atoms, {added} SHAKE constraints')
    P = c['pressure_atm'] * ATM_TO_BAR * u.bar
    T0 = c['stages'][0][1]
    baro = mm.MonteCarloAnisotropicBarostat(mm.Vec3(P, P, P), T0 * u.kelvin, False, False, True,
                                            c['default_barostat_every'])
    system.addForce(baro)
    integ = integrator(T0, c['dt_fs'], cfg['default_friction_per_ps'])
    ctx = context_for(cfg, system, integ)
    # Whole molecules: OpenMM constrains raw coordinates, so an X-H bond split
    # across a face would be dragged across the cell.
    pos0 = T.unwrap(T.positions_A(start), T.box_A(start), T.bonds(system))
    # Not a minimiser: L-BFGS drives a fresh build's worst torsions into a
    # near-linear, singular geometry. A displacement-limited Langevin run
    # (as LAMMPS fix nve/limit) relaxes the build strain instead.
    pr = c['default_prerelax']
    pos1, e = T.prerelax(str(build) + '.openmm_system.xml', pos0, T.box_A(start), steps=s.steps(pr['steps']),
                         T_K=pr['T_K'], friction_per_ps=pr['friction_per_ps'], dt_fs=pr['dt_fs'],
                         max_step_A=pr['max_step_A'], platform=cfg['platform'], precision=cfg['precision'])
    s.say(f'pre-relaxed ({s.steps(pr["steps"])} steps, <= {pr["max_step_A"]} A/step): PE -> '
          + ' -> '.join(f'{x:.0f}' for x in e) + ' kcal/mol')
    ctx.setState(start)
    ctx.setPositions(pos1 * 0.1)
    ctx.applyConstraints(1e-6)
    ctx.setVelocitiesToTemperature(T0 * u.kelvin, 2026)

    def set_T(Tk):
        integ.setTemperature(Tk * u.kelvin)
        ctx.setParameter(mm.MonteCarloAnisotropicBarostat.Temperature(), Tk)
    log = T.Log(s.dir / 'log.csv', ['step', 'time_ps', 'T_target_K', 'T_K', 'PE_kcal', 'Lz_A', 'density_g_cm3'])
    dof = T.degrees_of_freedom(system)
    done = 0
    for kind, ta, tb, n in c['stages']:
        n = s.steps(n)
        base = done

        def report(k, ta=ta, tb=tb, n=n, base=base):
            if k % 1000 and k != n:
                return
            st = ctx.getState(getEnergy=True)
            box = T.box_A(st)
            Tt = ta + (tb - ta) * k / n
            log.row(base + k, (base + k) * c['dt_fs'] / 1000, round(Tt, 2), round(T.temperature(ctx, system, dof), 2),
                    round(kcal(st.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)), 3),
                    round(box[2], 4), round(density_g_cm3(system, box), 5))
        T.run_ramp(integ, n, ta, tb, set_temperature=set_T, chunk=1000, on_chunk=report)
        done += n
        st = ctx.getState()
        s.say(f'{kind} {ta}->{tb} K, {n} steps: Lz {T.box_A(st)[2]:.3f} A, '
              f'density {density_g_cm3(system, T.box_A(st)):.4f} g/cm3')
    log.close()
    st = ctx.getState()
    finish(s, ctx, system, str(build) + '.data', str(build) + '.in.styles',
           extra={'density_g_cm3': density_g_cm3(system, T.box_A(st)), 'steps': done})


# --- 03 surface: vacuum in z, NVT ------------------------------------------------

def stage_surface(cfg, runs):
    s = Stage(cfg, runs, '03_surface')
    c = cfg['surface']
    build = runs / '01_build' / 'polymer'
    system = T.load_xml(str(build) + '.openmm_system.xml')
    prev = T.load_xml(runs / '02_bulk' / 'state.xml')
    box = T.box_A(prev)
    pos = T.unwrap(T.positions_A(prev), box, T.bonds(system))   # whole molecules, never split
    zmin, zmax = pos[:, 2].min(), pos[:, 2].max()
    vac = c['default_vacuum_A']
    pos[:, 2] += vac / 2 - zmin
    new_box = np.array([box[0], box[1], (zmax - zmin) + vac])
    s.say(f'bulk Lz {box[2]:.3f} A; whole molecules span z {zmax - zmin:.3f} A; new Lz {new_box[2]:.3f} A')
    system.setDefaultPeriodicBoxVectors(*[mm.Vec3(*v) for v in np.diag(new_box * 0.1)])
    added, _ = T.configure(system, shake=c['shake'])
    kind, ta, tb, n = c['stages'][0]
    integ = integrator(ta, c['dt_fs'], cfg['default_friction_per_ps'])
    ctx = context_for(cfg, system, integ)
    ctx.setPeriodicBoxVectors(*[mm.Vec3(*v) for v in np.diag(new_box * 0.1)])
    ctx.setPositions(pos * 0.1)
    ctx.applyConstraints(1e-6)
    ctx.setVelocitiesToTemperature(ta * u.kelvin, 2027)
    n = s.steps(n)
    log = T.Log(s.dir / 'log.csv', ['step', 'time_ps', 'T_target_K', 'T_K', 'PE_kcal', 'z_min_A', 'z_max_A'])
    dof = T.degrees_of_freedom(system)

    def report(k):
        if k % 1000 and k != n:
            return
        st = ctx.getState(getEnergy=True, getPositions=True)
        z = T.positions_A(st)[:, 2]
        log.row(k, k * c['dt_fs'] / 1000, round(ta + (tb - ta) * k / n, 2), round(T.temperature(ctx, system, dof), 2),
                round(st.getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole), 3), round(z.min(), 3), round(z.max(), 3))
    T.run_ramp(integ, n, ta, tb, chunk=1000, on_chunk=report)
    log.close()
    finish(s, ctx, system, str(build) + '.data', str(build) + '.in.styles', extra={'steps': n})


# --- 04 assemble: PAVES combines silica supercell + polymer slab ----------------

def stage_assemble(cfg, runs):
    s = Stage(cfg, runs, '04_assemble')
    c = cfg['assemble']
    build = runs / '01_build' / 'polymer'
    polymer_sys = T.load_xml(str(build) + '.openmm_system.xml')
    prev = T.load_xml(runs / '03_surface' / 'state.xml')
    box = T.box_A(prev)
    pos = T.unwrap(T.positions_A(prev), box, T.bonds(polymer_sys))

    # The silica model's own extent, from the file PAVES will replicate.
    silica_z = []
    section = ''
    for line in Path(cfg['silica_mol2']).read_text().splitlines():
        if line.startswith('@<TRIPOS>'):
            section = line[9:].strip()
            continue
        f = line.split()
        if section == 'ATOM' and len(f) >= 6:
            silica_z.append(float(f[4]))
    s_top = max(silica_z)
    nx, ny, nz = c['supercell']
    pos[:, 2] += (s_top + c['default_gap_A']) - pos[:, 2].min()
    Lz = pos[:, 2].max() + c['default_vacuum_above_A']
    s.say(f'silica top {s_top:.3f} A; polymer {pos[:, 2].min():.3f}..{pos[:, 2].max():.3f} A; Lz {Lz:.3f} A')
    cell = [box[0], box[1], Lz]
    T.write_lammps_data(str(build) + '.data', s.dir / 'polymer_slab.data', pos, cell, title='polymer slab for assembly')

    inp = s.dir / 'interface.paves'
    inp.write_text(f"""# Written by pipeline.py assemble. Silica keeps its IFF types and charges;
# the polymer keeps the PCFF types and bond-increment charges it was built with.
name          {cfg['name']}_interface
seed          2026
forcefield    {cfg['forcefield']}
environment   silica file {cfg['silica_mol2']} format mol2 parameterization preserve mobility fixed replicate {nx} {ny} {nz}
environment   polymer file polymer_slab.data format data parameterization preserve mobility fixed
cell          {cell[0]:.6f} {cell[1]:.6f} {cell[2]:.6f}
temperature   550
openmm_system yes
output        out interface
""")
    (s.dir / 'out').mkdir(exist_ok=True)
    polypaves_build(cfg, s.dir, inp.name)

    system = T.load_xml(s.dir / 'out/interface.openmm_system.xml')
    st = T.load_xml(s.dir / 'out/interface.openmm_state.xml')
    p = T.positions_A(st)
    n_si = system.getNumParticles() - len(pos)
    # The polymer must come out exactly where it went in, and after the silica.
    d = p[n_si:] - pos
    d[:, :2] -= box[:2] * np.round(d[:, :2] / box[:2])     # the data file carries x, y wrapped with images
    drift = float(np.abs(d).max())
    if drift > 1e-6:
        sys.exit(f'assembled polymer moved by {drift:.3e} A; atom order or coordinates are not preserved')
    zs = p[:n_si, 2]
    fixed = np.where(zs < zs.min() + cfg['compress']['fix_bottom_A'])[0]
    info = {'silica_atoms': int(n_si), 'polymer_atoms': int(len(pos)), 'atoms': int(len(p)),
            'box_A': T.box_A(st).tolist(), 'silica_z_A': [float(zs.min()), float(zs.max())],
            'polymer_z_A': [float(pos[:, 2].min()), float(pos[:, 2].max())],
            'gap_A': float(pos[:, 2].min() - zs.max()), 'fixed_atoms': int(len(fixed)),
            'polymer_position_check_A': drift}
    np.savetxt(s.dir / 'fixed_atoms.txt', fixed, fmt='%d')
    T.dump_json(s.dir / 'summary.json', info)
    s.say(json.dumps(info))


def interface_files(runs):
    d = runs / '04_assemble'
    info = json.loads((d / 'summary.json').read_text())
    return (d / 'out/interface.openmm_system.xml', d / 'out/interface.data', d / 'out/interface.in.styles',
            np.loadtxt(d / 'fixed_atoms.txt', dtype=int, ndmin=1), info)


def interface_context(cfg, runs, c, T_K, *, state):
    xml, data, styles, fixed, info = interface_files(runs)
    system = T.load_xml(xml)
    added, skipped = T.configure(system, shake=c['shake'], fixed=fixed)
    return system, data, styles, fixed, info, added, skipped


def add_dcd(stage, system, ctx, integ, every):
    """A DCD reporter driven by hand (no Simulation object), plus topology.pdb."""
    if not every:
        return None
    from openmm.app import DCDFile
    top, _ = T.topology(system)
    st = ctx.getState(getPositions=True)
    T.write_pdb(stage.dir / 'topology.pdb', top, T.positions_A(st), T.box_A(st))
    fh = open(stage.dir / 'trajectory.dcd', 'wb')
    dcd = DCDFile(fh, top, integ.getStepSize(), 0, every)

    def write():
        st = ctx.getState(getPositions=True)
        dcd.writeModel(st.getPositions(), periodicBoxVectors=st.getPeriodicBoxVectors())
    return fh, write


# --- 05 compress: wall at 200 MPa --------------------------------------------

def stage_compress(cfg, runs):
    s = Stage(cfg, runs, '05_compress')
    c = cfg['compress']
    system, data, styles, fixed, info, added, skipped = interface_context(cfg, runs, c, c['temperature'], state=None)
    n_si = info['silica_atoms']
    start = T.load_xml(runs / '04_assemble' / 'out/interface.openmm_state.xml')
    box = T.box_A(start)
    area_m2 = box[0] * box[1] * 1e-20
    PA = c['pressure_MPa'] * 1e6 * area_m2 * AVOGADRO / 1000 * 1e-9   # kJ/mol/nm
    kw = c['default_wall_k_kcal_A2'] * T.KCAL * 100                      # kJ/mol/nm^2
    wall = mm.CustomExternalForce('kw*step(z-zw)*(z-zw)^2')
    wall.addGlobalParameter('kw', kw)
    z_top = T.positions_A(start)[n_si:, 2].max()
    wall.addGlobalParameter('zw', (z_top + 1.0) * 0.1)
    for i in range(n_si, system.getNumParticles()):
        wall.addParticle(i, [])
    wall.setForceGroup(31)
    system.addForce(wall)
    s.say(f'{system.getNumParticles()} atoms, {len(fixed)} fixed, SHAKE {added} (+{skipped} skipped at fixed atoms); '
          f'wall target {PA / T.KCAL / 10:.2f} kcal/mol/A = {c["pressure_MPa"]} MPa on {box[0]:.4f} x {box[1]:.4f} A')
    integ = integrator(c['temperature'], c['dt_fs'], cfg['default_friction_per_ps'])
    ctx = context_for(cfg, system, integ)
    xml, *_ = interface_files(runs)
    pr = c['default_prerelax']
    pos1, e = T.prerelax(xml, T.unwrap(T.positions_A(start), box, T.bonds(system)), box, fixed=fixed,
                         steps=s.steps(pr['steps']), T_K=pr['T_K'], friction_per_ps=pr['friction_per_ps'],
                         dt_fs=pr['dt_fs'], max_step_A=pr['max_step_A'],
                         platform=cfg['platform'], precision=cfg['precision'])
    s.say(f'pre-relaxed ({s.steps(pr["steps"])} steps): PE -> ' + ' -> '.join(f'{x:.0f}' for x in e) + ' kcal/mol')
    ctx.setState(start)
    ctx.setPositions(pos1 * 0.1)
    ctx.applyConstraints(1e-6)
    ctx.setVelocitiesToTemperature(c['temperature'] * u.kelvin, 2028)
    zw = (z_top + 1.0) * 0.1
    ctx.setParameter('zw', zw)
    every = c['default_wall_update_every']
    max_move = c['default_wall_max_move_A'] * 0.1
    mobility = max_move / PA
    n = s.steps(c['steps'])
    dcd = add_dcd(s, system, ctx, integ, c['dcd_every'])
    log = T.Log(s.dir / 'log.csv', ['step', 'time_ps', 'T_K', 'PE_kcal', 'wall_z_A', 'wall_force_kcal_A',
                                    'wall_pressure_MPa', 'polymer_z_min_A', 'polymer_z_max_A', 'gap_A'])
    dof = T.degrees_of_freedom(system)
    to_MPa = 1.0 / (area_m2 * AVOGADRO / 1000 * 1e-9) / 1e6
    f_up = 0.0
    for k in range(0, n, every):
        integ.step(min(every, n - k))
        done = min(k + every, n)
        # The polymer's push on the wall: sum of 2 k (z - zw) over atoms beyond it (kJ/mol/nm).
        zp = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(u.nanometer)[n_si:, 2]
        over = zp[zp > zw] - zw
        f_up = float(2 * kw * over.sum())
        zw += float(np.clip(mobility * (f_up - PA), -max_move, max_move))
        ctx.setParameter('zw', zw)
        if dcd and done % c['dcd_every'] == 0:
            dcd[1]()
        if done % 1000 == 0 or done == n:
            st = ctx.getState(getEnergy=True, getPositions=True)
            z = T.positions_A(st)[:, 2]
            log.row(done, done * c['dt_fs'] / 1000, round(T.temperature(ctx, system, dof), 2),
                    round(st.getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole), 3), round(zw * 10, 4),
                    round(f_up / T.KCAL / 10, 3), round(f_up * to_MPa, 3),
                    round(z[n_si:].min(), 3), round(z[n_si:].max(), 3), round(z[n_si:].min() - z[:n_si].max(), 3))
    log.close()
    if dcd:
        dcd[0].close()
    ctx.setParameter('kw', 0.0)       # the wall is not part of what is written
    finish(s, ctx, system, data, styles, extra={'steps': n, 'final_wall_z_A': zw * 10,
                                               'final_wall_pressure_MPa': f_up * to_MPa})


# --- 06 cool / 07 relax: no wall, 0.25 fs -------------------------------------------

def stage_nvt(cfg, runs, name, prev, c, t_start, t_end, steps):
    s = Stage(cfg, runs, name)
    system, data, styles, fixed, info, added, skipped = interface_context(cfg, runs, c, t_start, state=None)
    start = T.load_xml(runs / prev / 'state.xml')
    integ = integrator(t_start, c['dt_fs'], cfg['default_friction_per_ps'])
    ctx = context_for(cfg, system, integ)
    ctx.setState(start)
    if added:
        ctx.applyConstraints(1e-6)
    s.say(f'{system.getNumParticles()} atoms, {len(fixed)} fixed, SHAKE {added}, dt {c["dt_fs"]} fs, '
          f'{t_start}->{t_end} K')
    n = s.steps(steps)
    n_si = info['silica_atoms']
    dcd = add_dcd(s, system, ctx, integ, c.get('dcd_every', 0))
    log = T.Log(s.dir / 'log.csv', ['step', 'time_ps', 'T_target_K', 'T_K', 'PE_kcal', 'polymer_z_min_A',
                                    'polymer_z_max_A', 'gap_A'])
    dof = T.degrees_of_freedom(system)

    def report(k):
        if dcd and k % c['dcd_every'] == 0:
            dcd[1]()
        if k % 4000 and k != n:
            return
        st = ctx.getState(getEnergy=True, getPositions=True)
        z = T.positions_A(st)[:, 2]
        log.row(k, k * c['dt_fs'] / 1000, round(t_start + (t_end - t_start) * k / n, 2),
                round(T.temperature(ctx, system, dof), 2),
                round(st.getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole), 3),
                round(z[n_si:].min(), 3), round(z[n_si:].max(), 3), round(z[n_si:].min() - z[:n_si].max(), 3))
    T.run_ramp(integ, n, t_start, t_end, chunk=1000, on_chunk=report)
    log.close()
    if dcd:
        dcd[0].close()
    finish(s, ctx, system, data, styles, extra={'steps': n})
    # Which atoms LAMMPS should hold still, as LAMMPS ids.
    np.savetxt(s.dir / 'fixed_atoms.lammps_ids.txt', fixed + 1, fmt='%d')


def stage_cool(cfg, runs):
    c = cfg['cool']
    stage_nvt(cfg, runs, '06_cool', '05_compress', c, c['t_start'], c['t_end'], c['steps'])


def stage_relax(cfg, runs):
    c = cfg['relax']
    stage_nvt(cfg, runs, '07_relax', '06_cool', c, c['temperature'], c['temperature'], c['steps'])


# --- 08 interface energy: E_int = E(all) - E(polymer) - E(slab) -----------------

def polypaves_build(cfg, workdir, inp_name):
    """Build a PAVES input in its directory with the installed polypaves package."""
    import contextlib
    import traceback
    from polypaves.config import from_file
    import os
    log = workdir / (inp_name + '.log')
    # The engine writes progress to the process's own stdout/stderr (C level), so
    # both file descriptors go to the log for the duration of the build.
    sys.stdout.flush()
    sys.stderr.flush()
    saved = os.dup(1), os.dup(2)
    with open(log, 'w') as fh:
        os.dup2(fh.fileno(), 1)
        os.dup2(fh.fileno(), 2)
        try:
            with contextlib.redirect_stderr(fh), contextlib.redirect_stdout(fh):
                try:
                    system = from_file(str(workdir / inp_name), write_files=True)
                    print(json.dumps(system.report, default=str))
                    ok = True
                except Exception:
                    traceback.print_exc()
                    ok = False
        finally:
            fh.flush()
            os.dup2(saved[0], 1)
            os.dup2(saved[1], 2)
            os.close(saved[0])
            os.close(saved[1])
    if not ok:
        sys.exit(f'polypaves could not build {workdir / inp_name}; see {log}')


def silica_only_system(cfg, runs):
    """The substrate alone, built by PAVES in the interface cell (same atoms, same order)."""
    d = runs / '04_assemble'
    xml = d / 'out_silica/silica.openmm_system.xml'
    if xml.exists():
        return xml
    info = json.loads((d / 'summary.json').read_text())
    box = info['box_A']
    nx, ny, nz = cfg['assemble']['supercell']
    (d / 'out_silica').mkdir(exist_ok=True)
    (d / 'silica_only.paves').write_text(f"""# Written by pipeline.py: the substrate alone, for interface energies.
name          silica_only
seed          2026
forcefield    {cfg['forcefield']}
environment   silica file {cfg['silica_mol2']} format mol2 parameterization preserve mobility fixed replicate {nx} {ny} {nz}
cell          {box[0]:.6f} {box[1]:.6f} {box[2]:.6f}
temperature   300
openmm_system yes
output        out_silica silica
""")
    polypaves_build(cfg, d, 'silica_only.paves')
    return xml


def nonbonded_groups(system):
    return {f.getForceGroup() for f in system.getForces()
            if isinstance(f, (mm.NonbondedForce, mm.CustomNonbondedForce))}


def energy_parts(ctx, system, positions_A, box_A):
    """(vdW, Coulomb) in kcal/mol from the nonbonded forces only."""
    ctx.setPeriodicBoxVectors(*[mm.Vec3(*v) for v in np.diag(np.asarray(box_A) * 0.1)])
    ctx.setPositions(np.asarray(positions_A) * 0.1)
    out = {}
    for f in system.getForces():
        if isinstance(f, (mm.NonbondedForce, mm.CustomNonbondedForce)):
            e = ctx.getState(getEnergy=True, groups={f.getForceGroup()}).getPotentialEnergy()
            out['coul' if isinstance(f, mm.NonbondedForce) else 'vdw'] = e.value_in_unit(u.kilocalorie_per_mole)
    return out['vdw'], out['coul']


def stage_interface(cfg, runs):
    s = Stage(cfg, runs, '08_interface')
    c = cfg['interface']
    system, data, styles, fixed, info, added, skipped = interface_context(cfg, runs, c, c['temperature'], state=None)
    n_si = info['silica_atoms']
    start = T.load_xml(runs / '07_relax' / 'state.xml')
    box = T.box_A(start)
    area = box[0] * box[1]
    # Evaluation copies, in double precision: the interaction energy is a small
    # difference of large sums. Only the nonbonded terms are evaluated: no bond,
    # angle or torsion crosses the interface, so every bonded term cancels exactly.
    def evaluator(xml):
        sysx = T.load_xml(xml)
        sysx.setDefaultPeriodicBoxVectors(*[mm.Vec3(*v) for v in np.diag(box * 0.1)])
        for f in sysx.getForces():
            if isinstance(f, mm.NonbondedForce):
                f.setEwaldErrorTolerance(c['default_ewald_tolerance'])
        return sysx, T.make_context(sysx, mm.VerletIntegrator(0.001), platform=cfg['platform'], precision='double')
    build = runs / '01_build' / 'polymer'
    ev_all = evaluator(interface_files(runs)[0])
    ev_poly = evaluator(str(build) + '.openmm_system.xml')
    ev_si = evaluator(silica_only_system(cfg, runs))
    if ev_poly[0].getNumParticles() + ev_si[0].getNumParticles() != system.getNumParticles():
        sys.exit('polymer + silica atom counts do not add up to the interface system')
    m_all, m_si = T.masses(ev_all[0]), T.masses(ev_si[0])
    if not np.allclose(m_all[:n_si], m_si):
        sys.exit('silica-only system is not the interface silica in the same order')

    integ = integrator(c['temperature'], c['dt_fs'], cfg['default_friction_per_ps'])
    ctx = context_for(cfg, system, integ)
    ctx.setState(start)
    n = s.steps(c['steps'])
    every = max(1, s.steps(c['snapshot_every']))
    kcal_A2_to_mJ_m2 = T.KCAL * 1000 / AVOGADRO / 1e-20 * 1000
    rows = []
    log = T.Log(s.dir / 'interface_energy.csv', ['step', 'time_ps', 'E_int_kcal', 'E_int_vdw_kcal', 'E_int_coul_kcal',
                                                  'gamma_mJ_m2', 'E_all_nb_kcal', 'E_polymer_nb_kcal', 'E_silica_nb_kcal'])
    for k in range(every, n + 1, every):
        integ.step(every)
        p = T.positions_A(ctx.getState(getPositions=True))
        va, ca = energy_parts(ev_all[1], ev_all[0], p, box)
        vp, cp = energy_parts(ev_poly[1], ev_poly[0], p[n_si:], box)
        vs, cs = energy_parts(ev_si[1], ev_si[0], p[:n_si], box)
        ev, ec = va - vp - vs, ca - cp - cs
        e = ev + ec
        g = e / area * kcal_A2_to_mJ_m2
        rows.append((e, ev, ec, g))
        log.row(k, k * c['dt_fs'] / 1000, round(e, 3), round(ev, 3), round(ec, 3), round(g, 3),
                round(va + ca, 3), round(vp + cp, 3), round(vs + cs, 3))
        s.say(f'step {k}: E_int {e:.1f} kcal/mol (vdW {ev:.1f}, Coulomb {ec:.1f}), {g:.1f} mJ/m^2')
    log.close()
    r = np.array(rows)
    sem = r.std(axis=0, ddof=1) / np.sqrt(len(r)) if len(r) > 1 else np.zeros(4)
    extra = {'snapshots': len(r), 'area_A2': area,
             'E_int_kcal': [float(r[:, 0].mean()), float(sem[0])],
             'E_int_vdw_kcal': [float(r[:, 1].mean()), float(sem[1])],
             'E_int_coul_kcal': [float(r[:, 2].mean()), float(sem[2])],
             'interface_energy_mJ_m2': [float(r[:, 3].mean()), float(sem[3])],
             'work_of_adhesion_mJ_m2': [float(-r[:, 3].mean()), float(sem[3])],
             'note': 'E_int = E_nb(all) - E_nb(polymer) - E_nb(silica), same coordinates and cell, '
                     'double precision; [mean, standard error over snapshots]; gamma = E_int / (Lx Ly)'}
    finish(s, ctx, system, data, styles, extra=extra)


# --- 09 pull: SMD on the top 25 % ------------------------------------------------

def stage_pull(cfg, runs):
    s = Stage(cfg, runs, '09_pull')
    c = cfg['pull']
    system, data, styles, fixed, info, added, skipped = interface_context(cfg, runs, c, c['temperature'], state=None)
    start = T.load_xml(runs / '08_interface' / 'state.xml')
    n_si = info['silica_atoms']
    pos = T.positions_A(start)
    z = pos[n_si:, 2]
    zn = (z - z.min()) / (z.max() - z.min())
    group = (np.where(zn >= 1.0 - c['top_fraction'])[0] + n_si).tolist()
    np.savetxt(s.dir / 'pulled_atoms.txt', group, fmt='%d')
    k_kj = c['k_kcal_A2'] * T.KCAL * 100                     # kJ/mol/nm^2
    com = mm.CustomCentroidBondForce(1, 'z1')                # mass-weighted: the centre of mass
    com.addGroup(group)
    com.addBond([0], [])
    spring = mm.CustomCVForce('0.5*kpull*(zc-zref)^2')
    spring.addCollectiveVariable('zc', com)
    spring.addGlobalParameter('kpull', k_kj)
    m = T.masses(system)
    z0 = float((pos[group, 2] * m[group]).sum() / m[group].sum())
    spring.addGlobalParameter('zref', z0 * 0.1)
    spring.setForceGroup(30)
    system.addForce(spring)
    # 'xy': no friction on the polymer's z velocity (LAMMPS temp/partial 1 1 0).
    # With 'all', Langevin drags the whole film towards v = 0 and, once it has
    # come off, the spring still reads M_polymer * gamma * v (~7 nN for PMMA
    # DP 200 x 20 at 10 m/s and 1/ps).
    mask = np.ones((system.getNumParticles(), 3))
    if c.get('default_thermostat', 'xy') == 'xy':
        mask[n_si:, 2] = 0.0
    integ = T.masked_langevin(c['temperature'], cfg['default_friction_per_ps'], c['dt_fs'], mask)
    ctx = context_for(cfg, system, integ)
    ctx.setState(start)
    n = s.steps(c['steps'])
    dt = c['dt_fs']
    v_nm_fs = c['speed_m_s'] * 1e-6                          # 1 m/s = 1e-5 A/fs = 1e-6 nm/fs
    s.say(f'{len(group)} pulled atoms (normalised z >= {1 - c["top_fraction"]:.2f}), z0 {z0:.4f} A, '
          f'k {c["k_kcal_A2"]} kcal/mol/A^2, v {c["speed_m_s"]} m/s over {n} steps of {dt} fs '
          f'= {v_nm_fs * 10 * n * dt:.4f} A')
    dcd = add_dcd(s, system, ctx, integ, c['dcd_every'])
    fh = open(s.dir / 'pull_force.csv', 'w')
    fh.write('step,time_ps,z_ref_A,z_com_A,extension_A,force_kcal_mol_A,force_nN\n')
    kcal_A_to_nN = T.KCAL * 1000 / AVOGADRO / 1e-10 * 1e9
    rows = []
    every = max(1, int(c.get('force_every', 100)))           # a row every `every` steps (and the last)
    for k in range(n + 1):
        zref = z0 * 0.1 + v_nm_fs * dt * k                   # the reference still moves every step
        ctx.setParameter('zref', zref)
        if k % every == 0 or k == n:
            zc = spring.getCollectiveVariableValues(ctx)[0]
            f = c['k_kcal_A2'] * (zref - zc) * 10              # kcal/mol/A, + = upward on the polymer
            rows.append(f'{k},{k * dt / 1000:.6f},{zref * 10:.6f},{zc * 10:.6f},{(zc * 10 - z0):.6f},{f:.6f},'
                        f'{f * kcal_A_to_nN:.6f}\n')
        if len(rows) >= 10000:
            fh.writelines(rows)
            rows.clear()
        if dcd and k % c['dcd_every'] == 0:
            dcd[1]()
        if k < n:
            integ.step(1)
    fh.writelines(rows)
    fh.close()
    if dcd:
        dcd[0].close()
    ctx.setParameter('kpull', 0.0)
    finish(s, ctx, system, data, styles, extra={'steps': n, 'pulled_atoms': len(group), 'z0_A': z0,
                                               'pull_distance_A': v_nm_fs * 10 * n * dt})


# name -> (run directory, function), in workflow order
STAGES = {
    'build': ('01_build', stage_build),
    'bulk': ('02_bulk', stage_bulk),
    'surface': ('03_surface', stage_surface),
    'assemble': ('04_assemble', stage_assemble),
    'compress': ('05_compress', stage_compress),
    'cool': ('06_cool', stage_cool),
    'relax': ('07_relax', stage_relax),
    'interface': ('08_interface', stage_interface),
    'pull': ('09_pull', stage_pull),
}
