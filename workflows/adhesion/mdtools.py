"""OpenMM helpers for the polypaves adhesion workflow.

The force field always comes from PAVES's export (`openmm_system yes`): a
serialised OpenMM System whose Class II terms were written by C++ (PAVES D069,
D075). Nothing here changes a coefficient. What a stage may change is the
integration setup around it:

* SHAKE: bonds to hydrogen constrained at their Class II r0 (dt = 1 fs), or
  none (dt = 0.25 fs);
* fixed atoms: mass 0, which OpenMM integrators leave in place;
* the periodic box.

Units at the interface of this module: Angstrom, fs, K, kcal/mol, unless a name
says otherwise. OpenMM itself is in nm, ps, kJ/mol.
"""
from __future__ import annotations

import csv
import json
import math
import time
from collections import deque
from pathlib import Path

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u

KCAL = 4.184                      # kJ per kcal
H_MASS_MAX = 1.2                  # amu; a lighter particle is a hydrogen
KB_KJ = 0.0083144626              # kJ/mol/K


# --- System and state -------------------------------------------------------

def load_xml(path):
    return mm.XmlSerializer.deserialize(Path(path).read_text())


def save_state(context, path, *, velocities=True):
    st = context.getState(getPositions=True, getVelocities=velocities, enforcePeriodicBox=False)
    Path(path).write_text(mm.XmlSerializer.serialize(st))
    return st


def bond_force(system):
    """The Class II bond term: CustomBondForce with d = r - r0."""
    for f in system.getForces():
        if isinstance(f, mm.CustomBondForce) and 'r0' in [f.getPerBondParameterName(k)
                                                           for k in range(f.getNumPerBondParameters())]:
            return f
    raise ValueError('no Class II bond force (CustomBondForce with r0) in this System')


def bonds(system):
    """(i, j, r0 in nm) for every bond."""
    f = bond_force(system)
    names = [f.getPerBondParameterName(k) for k in range(f.getNumPerBondParameters())]
    k0 = names.index('r0')
    out = []
    for b in range(f.getNumBonds()):
        i, j, p = f.getBondParameters(b)
        out.append((i, j, p[k0]))
    return out


def configure(system, *, shake, fixed=()):
    """Apply SHAKE (bonds to H at r0) and fixed atoms (mass 0) to a System in place.

    A constraint touching a fixed atom is skipped: OpenMM forbids constraints on
    massless particles, and a fixed atom does not move anyway.
    Returns (constraints added, constraints skipped because of fixed atoms).
    """
    fixed = set(int(i) for i in fixed)
    mass = masses(system)                     # before any is zeroed
    for i in fixed:
        system.setParticleMass(i, 0.0)
    added = skipped = 0
    if shake:
        for i, j, r0 in bonds(system):
            if mass[i] >= H_MASS_MAX and mass[j] >= H_MASS_MAX:
                continue                      # not a bond to hydrogen
            if i in fixed or j in fixed:
                skipped += 1
                continue
            system.addConstraint(i, j, r0)
            added += 1
    return added, skipped


def degrees_of_freedom(system):
    n = sum(1 for i in range(system.getNumParticles()) if system.getParticleMass(i).value_in_unit(u.dalton) > 0)
    return 3 * n - system.getNumConstraints() - 3


def temperature(context, system, dof=None):
    ke = context.getState(getEnergy=True).getKineticEnergy().value_in_unit(u.kilojoule_per_mole)
    return 2 * ke / ((dof or degrees_of_freedom(system)) * KB_KJ)


def minimize(context, *, tolerance=10.0, iterations=5000, rounds=5, rel=1e-4):
    """L-BFGS until the energy stops falling; returns [energy per round] in kcal/mol.

    One call can return at once from a strained start (a failed first line
    search), so it is repeated until a round changes the energy by less than
    `rel` or `rounds` is reached.
    """
    e = [context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)]
    for _ in range(rounds):
        mm.LocalEnergyMinimizer.minimize(context, tolerance, iterations)
        e.append(context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole))
        if abs(e[-2] - e[-1]) <= rel * max(1.0, abs(e[-1])) and len(e) > 2:
            break
    return e


def limited_langevin(T_K, friction_per_ps, dt_fs, max_step_A):
    """Langevin (BAOAB) with every velocity component capped, as LAMMPS fix nve/limit.

    No atom moves more than `max_step_A` per axis per step, however large its
    force. That lets a strained start (a fresh build, whose torsion terms can be
    near-singular) relax without one atom being thrown across the cell. Use it
    only for pre-relaxation: capping velocities is not dynamics. No constraints.
    """
    dt = dt_fs * 1e-3                                   # ps
    vmax = max_step_A * 0.1 / dt                        # nm/ps
    a = math.exp(-friction_per_ps * dt)
    ci = mm.CustomIntegrator(dt)
    ci.addGlobalVariable('a', a)
    ci.addGlobalVariable('b', math.sqrt(1 - a * a))
    ci.addGlobalVariable('kT', KB_KJ * T_K)
    ci.addGlobalVariable('vmax', vmax)
    ci.addUpdateContextState()
    # select(m, ..., 0): a fixed atom (mass 0) keeps v = 0 and so never moves.
    kick = 'select(m, max(-vmax, min(vmax, v + 0.5*dt*f/m)), 0)'
    ci.addComputePerDof('v', kick)
    ci.addComputePerDof('x', 'x + 0.5*dt*v')
    ci.addComputePerDof('v', 'select(m, a*v + b*sqrt(kT/m)*gaussian, 0)')
    ci.addComputePerDof('x', 'x + 0.5*dt*v')
    ci.addComputePerDof('v', kick)
    return ci


def masked_langevin(T_K, friction_per_ps, dt_fs, thermostatted):
    """LangevinMiddleIntegrator, with the thermostat applied only to chosen degrees of freedom.

    `thermostatted` is an (N, 3) array of 0/1: a 1 gets friction and noise, a 0
    is integrated as NVE (it exchanges energy with the others through the
    forces only). Pulling uses it to leave z free on the polymer, as LAMMPS
    `compute temp/partial 1 1 0` does: plain Langevin damps every velocity
    towards 0 in the lab frame, so a film dragged at v feels an extra
    -M*gamma*v that the spring has to supply. Same step sequence as OpenMM's
    LangevinMiddleIntegrator (constraints included); fixed atoms (mass 0) stay put.
    """
    dt = dt_fs * 1e-3                                   # ps
    a = math.exp(-friction_per_ps * dt)
    ci = mm.CustomIntegrator(dt)
    ci.addGlobalVariable('a', a)
    ci.addGlobalVariable('b', math.sqrt(1 - a * a))
    ci.addGlobalVariable('kT', KB_KJ * T_K)
    ci.addPerDofVariable('w', 0)
    ci.addPerDofVariable('x1', 0)
    ci.setPerDofVariableByName('w', [mm.Vec3(*map(float, r)) for r in np.asarray(thermostatted)])
    ci.addUpdateContextState()
    ci.addComputePerDof('v', 'select(m, v + dt*f/m, 0)')
    ci.addConstrainVelocities()
    ci.addComputePerDof('x', 'x + 0.5*dt*v')
    ci.addComputePerDof('v', 'select(m, w*(a*v + b*sqrt(kT/m)*gaussian) + (1-w)*v, 0)')
    ci.addComputePerDof('x', 'x + 0.5*dt*v')
    ci.addComputePerDof('x1', 'x')
    ci.addConstrainPositions()
    ci.addComputePerDof('v', 'v + (x-x1)/dt')
    return ci


def prerelax(system_xml, positions_A, box_A, *, fixed=(), steps=20000, T_K=300.0, friction_per_ps=10.0,
             dt_fs=1.0, max_step_A=0.1, platform='auto', precision='mixed', report=None):
    """Relax a strained structure with limited_langevin; returns (positions_A, [PE per 4000 steps]).

    Runs on a fresh copy of the System: no constraints, no barostat, no extra
    forces; fixed atoms keep mass 0 and stay where they are.
    """
    system = load_xml(system_xml)
    configure(system, shake=False, fixed=fixed)
    integ = limited_langevin(T_K, friction_per_ps, dt_fs, max_step_A)
    ctx = make_context(system, integ, platform=platform, precision=precision)
    ctx.setPeriodicBoxVectors(*[mm.Vec3(*v) for v in np.diag(np.asarray(box_A) * 0.1)])
    ctx.setPositions(np.asarray(positions_A) * 0.1)
    ctx.setVelocitiesToTemperature(T_K * u.kelvin, 1)
    energies = []
    for k in range(0, steps, 4000):
        integ.step(min(4000, steps - k))
        e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
        if not math.isfinite(e):
            raise RuntimeError(f'pre-relaxation diverged after {k} steps')
        energies.append(e)
        if report:
            report(min(k + 4000, steps), e)
    pos = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(u.angstrom)
    return np.asarray(pos), energies


_PLATFORM = {}


def pick_platform(platform='auto'):
    """'auto' takes the fastest platform that actually runs here: CUDA, then OpenCL, then CPU.

    A platform can be listed and still fail (a CUDA build newer than the driver
    gives CUDA_ERROR_UNSUPPORTED_PTX_VERSION), so each is tried on a two-particle
    System before it is chosen. The answer is kept for the rest of the process.
    """
    if platform != 'auto':
        return platform
    if 'auto' in _PLATFORM:
        return _PLATFORM['auto']
    names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    tried = []
    for name in ('CUDA', 'OpenCL', 'CPU', 'Reference'):
        if name not in names:
            continue
        try:
            probe = mm.System()
            probe.addParticle(1.0)
            probe.addParticle(1.0)
            f = mm.CustomBondForce('r^2')
            f.addBond(0, 1, [])
            probe.addForce(f)
            ctx = mm.Context(probe, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName(name))
            ctx.setPositions([mm.Vec3(0, 0, 0), mm.Vec3(0.1, 0, 0)])
            ctx.getState(getEnergy=True)
            del ctx
            _PLATFORM['auto'] = name
            if tried:
                print(f'[platform] using {name}; not usable here: ' + '; '.join(tried), flush=True)
            return name
        except Exception as exc:
            tried.append(f'{name} ({str(exc).splitlines()[0][:80]})')
    raise RuntimeError('no usable OpenMM platform: ' + '; '.join(tried))


def make_context(system, integrator, *, platform='auto', precision='mixed'):
    name = pick_platform(platform)
    plat = mm.Platform.getPlatformByName(name)
    props = {'Precision': precision} if name in ('CUDA', 'OpenCL') else {}
    return mm.Context(system, integrator, plat, props)


# --- Structure --------------------------------------------------------------

def masses(system):
    return np.array([system.getParticleMass(i).value_in_unit(u.dalton) for i in range(system.getNumParticles())])


_ELEMENT_BY_MASS = [(1.2, 'H'), (13.0, 'C'), (14.5, 'N'), (17.0, 'O'), (29.0, 'Si')]


def element_of(mass):
    for limit, symbol in _ELEMENT_BY_MASS:
        if mass < limit:
            return symbol
    return 'X'


def molecules(n, bond_list):
    """Connected components of the bond graph, as a molecule id per atom."""
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for i, j, *_ in bond_list:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
    roots, mol = {}, np.empty(n, int)
    for a in range(n):
        mol[a] = roots.setdefault(find(a), len(roots))
    return mol


def unwrap(positions, box, bond_list, axes=(True, True, True)):
    """Make every molecule whole by walking its bonds (minimum image per bond).

    positions in Angstrom, box = (Lx, Ly, Lz). Only the axes marked True are
    unwrapped; the others are left exactly as given.
    """
    pos = np.array(positions, float)
    n = len(pos)
    adj = [[] for _ in range(n)]
    for i, j, *_ in bond_list:
        adj[i].append(j)
        adj[j].append(i)
    box = np.asarray(box, float)
    mask = np.asarray(axes, bool)
    seen = np.zeros(n, bool)
    for start in range(n):
        if seen[start]:
            continue
        seen[start] = True
        queue = deque([start])
        while queue:
            a = queue.popleft()
            for b in adj[a]:
                if seen[b]:
                    continue
                d = pos[b] - pos[a]
                d[mask] -= box[mask] * np.round(d[mask] / box[mask])
                pos[b] = pos[a] + d
                seen[b] = True
                queue.append(b)
    return pos


def topology(system, bond_list=None):
    """An OpenMM Topology for reporters and viewers: one chain per molecule."""
    m = masses(system)
    bond_list = bonds(system) if bond_list is None else bond_list
    mol = molecules(len(m), bond_list)
    top = app.Topology()
    atoms = []
    chains = {}
    for i, mass in enumerate(m):
        mid = int(mol[i])
        if mid not in chains:
            ch = top.addChain(id=str(mid % 62))
            chains[mid] = top.addResidue('MOL' if mid else 'SUB', ch)
        el = element_of(mass)
        atoms.append(top.addAtom(el, app.Element.getBySymbol(el) if el != 'X' else None, chains[mid]))
    for i, j, *_ in bond_list:
        top.addBond(atoms[i], atoms[j])
    return top, mol


def write_pdb(path, top, positions_A, box_A):
    top.setPeriodicBoxVectors(np.diag(np.asarray(box_A) * 0.1) * u.nanometer)
    with open(path, 'w') as fh:
        app.PDBFile.writeFile(top, np.asarray(positions_A) * 0.1 * u.nanometer, fh, keepIds=False)


# --- Running ---------------------------------------------------------------

class Log:
    """CSV log of a stage: one row every `every` steps."""
    def __init__(self, path, fields):
        self.fh = open(path, 'w', newline='')
        self.w = csv.writer(self.fh)
        self.w.writerow(fields)
        self.t0 = time.time()

    def row(self, *values):
        self.w.writerow(values)
        self.fh.flush()

    def close(self):
        self.fh.close()


def run_ramp(integrator, steps, t_start, t_end, *, set_temperature=None, chunk=1000, on_chunk=None):
    """Integrate `steps`, moving the target temperature linearly from t_start to t_end.

    The target is updated every `chunk` steps, through `set_temperature(T)` if
    given (so a barostat can follow it), else on the integrator. on_chunk(done)
    is called after each chunk.
    """
    set_t = set_temperature or (lambda T: integrator.setTemperature(T * u.kelvin))
    done = 0
    while done < steps:
        n = min(chunk, steps - done)
        set_t(t_start + (t_end - t_start) * (done + 0.5 * n) / steps)
        integrator.step(n)
        done += n
        if on_chunk:
            on_chunk(done)
    return done


def box_A(state):
    return np.diag(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.angstrom)).copy()


def positions_A(state):
    return state.getPositions(asNumpy=True).value_in_unit(u.angstrom)


def velocities_A_per_fs(state):
    return state.getVelocities(asNumpy=True).value_in_unit(u.angstrom / u.femtosecond)


# --- LAMMPS data --------------------------------------------------------------

def write_lammps_data(template, out, positions_A, box, *, velocities=None, title=None, origin=(0.0, 0.0, 0.0)):
    """A LAMMPS data file: the template's types, charges and topology, new coordinates and box.

    `template` is a data file PAVES wrote for this exact system (same atom
    order). Coordinates are wrapped into [origin, origin + L) and the images
    are written as flags, so every molecule is whole on reading. Velocities, if
    given, are in A/fs (LAMMPS real units).
    """
    lines = Path(template).read_text().splitlines()
    pos = np.asarray(positions_A, float)
    box = np.asarray(box, float)
    lo = np.asarray(origin, float)
    img = np.floor((pos - lo) / box).astype(int)
    wrapped = pos - img * box
    out_lines = [title if title is not None else lines[0]]
    k = 1
    while k < len(lines):
        line = lines[k]
        w = line.split('#')[0].split()
        if len(w) == 4 and w[2] in ('xlo', 'ylo', 'zlo'):
            a = 'xyz'.index(w[2][0])
            out_lines.append(f'{lo[a]:.10f} {lo[a] + box[a]:.10f} {w[2]} {w[3]}')
        elif line.startswith('Velocities'):
            k += 2                                    # drop the template's; ours go at the end
            while k < len(lines) and lines[k].strip():
                k += 1
            continue
        elif line.startswith('Atoms'):
            out_lines += [line, '']
            k += 2
            while k < len(lines) and lines[k].strip():
                body, _, comment = lines[k].partition('#')
                f = body.split()
                i = int(f[0]) - 1
                x, y, z = wrapped[i]
                ix, iy, iz = img[i]
                tail = f' #{comment}' if comment else ''
                out_lines.append(f'{f[0]} {f[1]} {f[2]} {f[3]} {x:.10f} {y:.10f} {z:.10f} {ix} {iy} {iz}{tail}')
                k += 1
            continue
        else:
            out_lines.append(line)
        k += 1
    if velocities is not None:
        v = np.asarray(velocities, float)
        out_lines += ['', 'Velocities', '']
        out_lines += [f'{i + 1} {v[i, 0]:.10g} {v[i, 1]:.10g} {v[i, 2]:.10g}' for i in range(len(v))]
    Path(out).write_text('\n'.join(out_lines).rstrip() + '\n')
    return out


def write_styles(template_styles, out_styles, data_name):
    """The PAVES styles file, pointed at another data file."""
    text = Path(template_styles).read_text()
    lines = [(f'read_data {data_name}' if l.startswith('read_data') else l) for l in text.splitlines()]
    Path(out_styles).write_text('\n'.join(lines) + '\n')


def dump_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=1, default=float) + '\n')
