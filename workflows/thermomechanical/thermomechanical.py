#!/usr/bin/env python3
"""PCFF/OpenMM glass-transition, thermal-expansion, and tensile workflow.

Examples:
  python thermomechanical.py new pmma --monomer '*CC(C)(C(=O)OC)*' \
      --atoms-per-chain 1000 --total-atoms 20000
  python thermomechanical.py run projects/pmma
  python thermomechanical.py status projects/pmma
"""
from __future__ import annotations

import argparse
import contextlib
import copy
import csv
import json
import math
import os
import sys
import tempfile
import time
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BUILDER = HERE.parents[1]
sys.path.insert(0, str(HERE))


STAGES = {
    'build': '01_build',
    'cool': '02_cool',
    'elastic': '03_elastic',
    'analyze': '04_analyze',
}
AXIS = {'x': 0, 'y': 1, 'z': 2}
FORCEFIELDS = {'pcff': '@builder/examples/forcefields/pcff.ff'}
ARCHITECTURE_ALIASES = {'random': 'random-exact', 'block': 'blocks'}
OWNED_PAVES_RECORDS = {
    'forcefield', 'cell', 'density', 'periodic', 'output', 'openmm_system', 'temperature',
    'backend', 'platform', 'precision', 'region', 'soft_wall', 'relax', 'minimize',
    'degree', 'chains', 'total_atoms', 'atom_count_mode',
}
REFUSED_PAVES_RECORDS = {'environment', 'structure', 'structure_format', 'structure_topology',
                         'component', 'composition'}


def dump_json(path, value):
    Path(path).write_text(json.dumps(value, indent=1, default=float, allow_nan=False) + '\n')


def resolve(value, project):
    if not isinstance(value, str):
        return value
    for token, root in (('@builder/', BUILDER), ('@workflow/', HERE), ('@project/', project)):
        if value.startswith(token):
            return str(root / value[len(token):])
    return value


def load_project(path):
    source = Path(path).resolve()
    source = source / 'project.json' if source.is_dir() else source
    if not source.exists():
        sys.exit(f'project not found: {source}')
    cfg = json.loads(source.read_text())
    cfg['forcefield'] = resolve(cfg['forcefield'], source.parent)
    if not Path(cfg['forcefield']).exists():
        sys.exit(f'force field not found: {cfg["forcefield"]}')
    return source.parent, cfg


def rep_name(index):
    return f'replica_{index + 1:02d}'


def normalize_block_counts(lines):
    """Turn fixed `block A 70` records into scalable `fraction` records."""
    parsed = []
    fractional = False
    for i, line in enumerate(lines):
        words = line.split('#', 1)[0].split()
        if words and words[0] == 'block':
            if len(words) == 4 and words[2] == 'fraction':
                float(words[3])
                fractional = True
                continue
            if len(words) != 3:
                raise ValueError(f'invalid block record: {line}')
            parsed.append((i, words[1], float(words[2])))
    fixed = [(i, label, value) for i, label, value in parsed]
    if not fixed:
        return lines
    if fractional:
        raise ValueError('do not mix fixed block lengths and block fractions')
    if any(value <= 0 for _, _, value in fixed):
        raise ValueError('fixed block lengths must be positive')
    total = sum(value for _, _, value in fixed)
    out = list(lines)
    for i, label, value in fixed:
        out[i] = f'block         {label} fraction {value / total:.16g}'
    return out


def chemistry_from_file(path):
    """Keep scalable single-polymer chemistry and discard workflow-owned records."""
    kept = []
    sequence_seed = None
    for raw in Path(path).read_text().splitlines():
        line = raw.split('#', 1)[0].strip()
        words = line.split()
        if not words or words[0] == 'polypaves-build':
            continue
        key = words[0]
        if key in REFUSED_PAVES_RECORDS:
            raise ValueError(f"'{key}' is not supported here; provide one polymer/copolymer chemistry")
        if key in OWNED_PAVES_RECORDS or key in ('name', 'seed'):
            continue
        if key == 'sequence_seed':
            if len(words) != 2:
                raise ValueError('sequence_seed must have one integer value')
            sequence_seed = int(words[1])
            continue
        if key == 'sequence':
            raise ValueError('fixed explicit sequences cannot be resized by atom target; use monomer, '
                             'architecture and block/fraction records')
        kept.append(line)
    if not any(line.split()[0] == 'monomer' for line in kept):
        raise ValueError('the PAVES file contains no monomer record')
    kept = normalize_block_counts(kept)
    return kept, sequence_seed


def parse_assignments(values, option):
    result = {}
    for value in values or []:
        if '=' not in value:
            raise ValueError(f'{option} requires NAME=VALUE, got {value!r}')
        name, item = value.split('=', 1)
        if not name or not item or name in result:
            raise ValueError(f'invalid or duplicate {option}: {value!r}')
        result[name] = item
    return result


def chemistry_from_cli(args):
    monomers = {'A': args.monomer}
    monomers.update(parse_assignments(args.comonomer, '--comonomer'))
    architecture = ARCHITECTURE_ALIASES.get(args.architecture, args.architecture)
    if len(monomers) > 1 and not architecture:
        raise ValueError('--architecture is required when --comonomer is used')
    if len(monomers) == 1 and architecture:
        raise ValueError('--architecture requires at least one --comonomer')
    fractions = {key: float(value) for key, value in parse_assignments(args.fraction, '--fraction').items()}
    if architecture in ('random-exact', 'random-probabilistic', 'blocks'):
        if set(fractions) != set(monomers):
            raise ValueError(f'{architecture} requires one --fraction for every monomer: '
                             + ', '.join(monomers))
        if any(value <= 0 for value in fractions.values()) or not math.isclose(
                sum(fractions.values()), 1.0, rel_tol=0, abs_tol=1e-8):
            raise ValueError('--fraction values must be positive and sum to 1')
    elif fractions:
        raise ValueError('--fraction is used only with random or block architectures')
    lines = [f"monomer       {name} '{smiles}'" for name, smiles in monomers.items()]
    lines.append(f"terminator    '{args.terminator}'")
    if architecture:
        lines.append(f'architecture  {architecture}')
    for name, fraction in fractions.items():
        lines.append(f'block         {name} fraction {fraction:.16g}')
    return lines


def chemistry_lines(cfg):
    y = cfg['system']
    if 'chemistry_lines' in y:
        return list(y['chemistry_lines'])
    return [f"monomer       A '{y['monomer']}'", f"terminator    '{y['terminator']}'"]


def chain_metrics(monomer, terminator, dp, forcefield, density, temperature):
    import polypaves
    system = polypaves.build(
        monomer=monomer,
        terminator=terminator,
        dp=int(dp),
        chains=1,
        forcefield=forcefield,
        density=density,
        temperature=temperature,
        seed=17,
    )
    return int(system.n_atoms), float(sum(system.masses))


def chain_atom_count(monomer, terminator, dp, forcefield, density, temperature):
    return chain_metrics(monomer, terminator, dp, forcefield, density, temperature)[0]


def probe_copolymer_chain(cfg, dp, sequence_seed, directory):
    """Build one chain without writing outputs and return its atom count/mass."""
    from polypaves.config import from_file
    y = cfg['system']
    path = Path(directory) / f'probe_dp{dp}_seq{sequence_seed}.paves'
    chemistry = '\n'.join(chemistry_lines(cfg))
    path.write_text(f"""polypaves-build 1
name          sizing_probe
seed          17
sequence_seed {sequence_seed}
forcefield    {cfg['forcefield']}
{chemistry}
degree        {dp}
chains        1
density       {y['initial_density_g_cm3']}
temperature   {y['build_temperature_K']}
output        . probe
""")
    # Native progress is useful for final builds but too noisy for repeated
    # one-chain sizing probes.
    saved = os.dup(1), os.dup(2)
    with open(os.devnull, 'w') as sink:
        os.dup2(sink.fileno(), 1)
        os.dup2(sink.fileno(), 2)
        try:
            with contextlib.redirect_stdout(sink), contextlib.redirect_stderr(sink):
                system = from_file(str(path), write_files=False)
        finally:
            os.dup2(saved[0], 1)
            os.dup2(saved[1], 2)
            os.close(saved[0])
            os.close(saved[1])
    return int(system.n_atoms), float(sum(system.masses))


def resolve_copolymer_dp(cfg, sequence_seed, target, directory):
    """Find a DP whose realized one-chain atom count reaches the target."""
    cache = {}

    def measure(dp):
        if dp not in cache:
            cache[dp] = probe_copolymer_chain(cfg, dp, sequence_seed, directory)
        return cache[dp]

    first = 1
    while True:
        try:
            n1, _ = measure(first)
            break
        except Exception:
            first += 1
            if first > 10000:
                raise RuntimeError('no valid DP found for the supplied copolymer fractions')
    if n1 >= target:
        return first, cache
    second = first + 1
    while True:
        try:
            n2, _ = measure(second)
            break
        except Exception:
            second += 1
    slope = max(1.0, (n2 - n1) / (second - first))
    high = max(second, int(math.ceil(first + (target - n1) / slope)))
    nhigh, _ = measure(high)
    low = first
    while nhigh < target:
        low = high
        high = max(high + 1, high * 2)
        nhigh, _ = measure(high)
        if high > 10_000_000:
            raise RuntimeError('could not bracket a DP for the requested chain atom count')
    while high - low > 1:
        middle = (low + high) // 2
        nmiddle, _ = measure(middle)
        if nmiddle >= target:
            high = middle
        else:
            low = middle
    # A probabilistic architecture may not be perfectly monotonic when DP
    # changes because the full sequence is redrawn. Always verify and advance.
    while measure(high)[0] < target:
        high += 1
    return high, cache


def conservative_random_metrics(cfg, target, directory):
    """Worst-case DP/count/mass for a random-probabilistic copolymer.

    A probabilistic chain can, in principle, contain only the smallest repeat
    unit. Sizing against every corresponding homopolymer guarantees the
    per-chain atom target instead of guaranteeing it only on average.
    """
    lines = chemistry_lines(cfg)
    monomers = [line for line in lines if line.split()[0] == 'monomer']
    shared = [line for line in lines if line.split()[0] == 'terminator']
    candidates = []
    for monomer in monomers:
        variant = copy.deepcopy(cfg)
        variant['system']['chemistry_lines'] = [monomer, *shared]
        dp, _ = resolve_copolymer_dp(variant, 1, target, directory)
        candidates.append((variant, dp))
    safe_dp = max(dp for _, dp in candidates)
    metrics = [probe_copolymer_chain(variant, safe_dp, 1, directory)
               for variant, _ in candidates]
    return safe_dp, min(item[0] for item in metrics), min(item[1] for item in metrics)


def forcefield_cutoff_A(path):
    """Largest numeric pair_style cutoff, when stated in a PAVES .ff file."""
    for raw in Path(path).read_text().splitlines():
        words = raw.split('#', 1)[0].split()
        if words and words[0] == 'pair_style':
            numbers = []
            for word in words[2:]:
                try:
                    numbers.append(float(word))
                except ValueError:
                    pass
            return max(numbers) if numbers else None
    return None


def resolve_size(cfg):
    """Resolve a DP and chain count that meet both all-atom targets."""
    system_cfg = cfg['system']
    target_chain = int(system_cfg['target_atoms_per_chain'])
    if 'chemistry_lines' in system_cfg:
        with tempfile.TemporaryDirectory(prefix='paves-sizing-') as temporary:
            architecture = next((line.split()[1] for line in chemistry_lines(cfg)
                                 if line.split()[0] == 'architecture'), None)
            conservative = (conservative_random_metrics(cfg, target_chain, temporary)
                            if architecture == 'random-probabilistic' else None)
            seeds = [int(system_cfg.get('sequence_seed', 7)) + i
                     for i in range(int(cfg['replicates']))]
            resolved = [resolve_copolymer_dp(cfg, seed, target_chain, temporary)[0]
                        for seed in seeds]
            dp = max([*resolved, conservative[0] if conservative else 1])
            metrics = [probe_copolymer_chain(cfg, dp, seed, temporary) for seed in seeds]
        atoms_by_replica = [item[0] for item in metrics]
        masses_by_replica = [item[1] for item in metrics]
        actual_chain = min(atoms_by_replica)
        chain_mass_Da = min(masses_by_replica)
        if conservative:
            actual_chain = min(actual_chain, conservative[1])
            chain_mass_Da = min(chain_mass_Da, conservative[2])
        if actual_chain < target_chain:
            raise RuntimeError('a copolymer sequence fell below the per-chain atom target at the resolved DP')
        n1 = None
        increment = None
    else:
        args = (system_cfg['monomer'], system_cfg['terminator'])
        kwargs = {
            'forcefield': cfg['forcefield'],
            'density': system_cfg['initial_density_g_cm3'],
            'temperature': system_cfg['build_temperature_K'],
        }
        n1 = chain_atom_count(*args, 1, **kwargs)
        n2 = chain_atom_count(*args, 2, **kwargs)
        increment = n2 - n1
        if increment <= 0:
            raise RuntimeError(f'non-positive repeat-unit atom increment: DP1={n1}, DP2={n2}')
        dp = max(1, 1 + math.ceil((target_chain - n1) / increment))
        actual_chain, chain_mass_Da = chain_metrics(*args, dp, **kwargs)
        while actual_chain < target_chain:
            dp += 1
            actual_chain, chain_mass_Da = chain_metrics(*args, dp, **kwargs)
        atoms_by_replica = [actual_chain] * int(cfg['replicates'])
        masses_by_replica = [chain_mass_Da] * int(cfg['replicates'])
    total_target = int(system_cfg['target_total_atoms'])
    chains_for_target = max(1, math.ceil(total_target / actual_chain))
    cutoff_A = forcefield_cutoff_A(cfg['forcefield'])
    chains_for_cutoff = 1
    if cutoff_A is not None:
        # A cubic density build must exceed twice the nonbonded cutoff for
        # OpenMM's minimum-image periodic nonbonded forces. Add a 0.1% margin.
        min_volume_A3 = (2.002 * cutoff_A) ** 3
        guard_density = float(system_cfg.get('minimum_image_density_g_cm3',
                                              system_cfg['initial_density_g_cm3']))
        mass_per_volume = chain_mass_Da * 1.66053906660 / guard_density
        chains_for_cutoff = max(1, math.ceil(min_volume_A3 / mass_per_volume))
    chains = max(chains_for_target, chains_for_cutoff)
    return {
        'dp': dp,
        'chains': chains,
        'atoms_per_chain': actual_chain,
        'total_atoms': actual_chain * chains,
        'target_atoms_per_chain': target_chain,
        'target_total_atoms': total_target,
        'dp1_atoms': n1,
        'repeat_increment_atoms': increment,
        'chain_mass_Da': chain_mass_Da,
        'pair_cutoff_A': cutoff_A,
        'chains_required_by_atom_target': chains_for_target,
        'chains_required_by_minimum_image': chains_for_cutoff,
        'minimum_image_density_g_cm3': system_cfg.get('minimum_image_density_g_cm3'),
        'atoms_per_chain_by_replica': atoms_by_replica,
        'chain_mass_Da_by_replica': masses_by_replica,
    }


def polypaves_build(input_path):
    """Run PAVES while capturing native stdout/stderr next to the input."""
    from polypaves.config import from_file
    workdir = input_path.parent
    log = workdir / 'build.log'
    saved = os.dup(1), os.dup(2)
    ok = False
    system = None
    sys.stdout.flush()
    sys.stderr.flush()
    with open(log, 'w') as fh:
        os.dup2(fh.fileno(), 1)
        os.dup2(fh.fileno(), 2)
        try:
            with contextlib.redirect_stdout(fh), contextlib.redirect_stderr(fh):
                try:
                    system = from_file(str(input_path), write_files=True)
                    print(json.dumps(system.report, default=str))
                    ok = True
                except Exception:
                    traceback.print_exc()
        finally:
            fh.flush()
            os.dup2(saved[0], 1)
            os.dup2(saved[1], 2)
            os.close(saved[0])
            os.close(saved[1])
    if not ok:
        raise RuntimeError(f'PAVES build failed; see {log}')
    return system


def stage_build(project, cfg, force=False):
    out = project / 'runs' / STAGES['build']
    out.mkdir(parents=True, exist_ok=True)
    sizing_path = out / 'sizing.json'
    if force or not sizing_path.exists():
        sizing = resolve_size(cfg)
        dump_json(sizing_path, sizing)
    else:
        sizing = json.loads(sizing_path.read_text())
    print('[build] resolved sizing:', json.dumps(sizing), flush=True)
    summaries = []
    for index in range(int(cfg['replicates'])):
        directory = out / rep_name(index)
        summary_path = directory / 'summary.json'
        if summary_path.exists() and not force:
            summaries.append(json.loads(summary_path.read_text()))
            continue
        directory.mkdir(parents=True, exist_ok=True)
        seed = int(cfg['base_seed']) + index
        y = cfg['system']
        sequence_seed = int(y.get('sequence_seed', 7)) + index
        chemistry = '\n'.join(chemistry_lines(cfg))
        input_path = directory / 'system.paves'
        input_path.write_text(f"""polypaves-build 1
# Generated by thermomechanical.py. Edit project.json, then rerun build.
name          {cfg['name']}_{rep_name(index)}
seed          {seed}
sequence_seed {sequence_seed}
forcefield    {cfg['forcefield']}
{chemistry}
degree        {sizing['dp']}
chains        {sizing['chains']}
density       {y['initial_density_g_cm3']}
temperature   {y['build_temperature_K']}
openmm_system yes
output        . system
""")
        print(f'[build] {rep_name(index)} seed={seed}', flush=True)
        built = polypaves_build(input_path)
        actual = int(built.n_atoms)
        if actual < sizing['target_total_atoms']:
            raise RuntimeError(f'built {actual} atoms, below target {sizing["target_total_atoms"]}')
        actual_per_chain = actual // int(sizing['chains'])
        if actual_per_chain < sizing['target_atoms_per_chain']:
            raise RuntimeError(f'built {actual_per_chain} atoms/chain, below target '
                               f'{sizing["target_atoms_per_chain"]}')
        summary = {
            'replica': index + 1,
            'seed': seed,
            'sequence_seed': sequence_seed,
            'dp': sizing['dp'],
            'chains': sizing['chains'],
            'atoms_per_chain': actual_per_chain,
            'atoms': actual,
            'initial_density_g_cm3': y['initial_density_g_cm3'],
        }
        dump_json(summary_path, summary)
        summaries.append(summary)
    dump_json(out / 'summary.json', {'stage': 'build', 'sizing': sizing, 'replicas': summaries})


def test_steps(cfg, production, kind):
    if not cfg.get('test_mode'):
        return max(1, int(round(production)))
    return {'prerelax': 2, 'equilibrate': 2, 'cool': 120, 'elastic': 60}.get(kind, 2)


def stage_cool(project, cfg, force=False):
    import openmm as mm
    import openmm.unit as u
    import mdtools as M

    build_dir = project / 'runs' / STAGES['build']
    out = project / 'runs' / STAGES['cool']
    out.mkdir(parents=True, exist_ok=True)
    md = cfg['md']
    cooling = cfg['cooling']
    dt_fs = float(md['dt_fs'])
    duration_ns = ((cooling['temperature_start_K'] - cooling['temperature_end_K']) /
                   cooling['rate_K_per_ns'])
    production_cool_steps = duration_ns * 1e6 / dt_fs
    cool_steps = test_steps(cfg, production_cool_steps, 'cool')
    equil_steps = test_steps(cfg, cooling['equilibrate_ns'] * 1e6 / dt_fs, 'equilibrate')
    prere_steps = test_steps(cfg, cooling['prerelax_steps'], 'prerelax')
    update = 1 if cfg.get('test_mode') else int(cooling['temperature_update_steps'])
    sample = 2 if cfg.get('test_mode') else int(cooling['sample_every_steps'])
    barostat_every = 1 if cfg.get('test_mode') else int(md['barostat_every_steps'])
    summaries = []
    for index in range(int(cfg['replicates'])):
        directory = out / rep_name(index)
        summary_path = directory / 'summary.json'
        if summary_path.exists() and not force:
            summaries.append(json.loads(summary_path.read_text()))
            continue
        directory.mkdir(parents=True, exist_ok=True)
        source = build_dir / rep_name(index)
        xml = source / 'system.openmm_system.xml'
        initial = M.load_xml(source / 'system.openmm_state.xml')
        system = M.load_xml(xml)
        constraints = M.configure(system, md['shake_hydrogen'])
        pressure = md['pressure_atm'] * M.ATM_TO_BAR * u.bar
        barostat = mm.MonteCarloBarostat(pressure, cooling['temperature_start_K'] * u.kelvin,
                                        barostat_every)
        system.addForce(barostat)
        integrator = mm.LangevinMiddleIntegrator(cooling['temperature_start_K'] * u.kelvin,
                                                 md['friction_per_ps'] / u.picosecond,
                                                 dt_fs * u.femtosecond)
        context = M.make_context(system, integrator, cfg['platform'], cfg['precision'])
        prere_cfg = {
            'temperature_K': cooling['temperature_start_K'],
            'friction_per_ps': cooling['prerelax_friction_per_ps'],
            'dt_fs': cooling['prerelax_dt_fs'],
            'max_step_A': cooling['prerelax_max_step_A'],
            'steps': prere_steps,
            'seed': int(cfg['base_seed']) + index,
        }
        print(f'[cool] {rep_name(index)}: prerelax {prere_steps}, equilibrate {equil_steps}, '
              f'cool {cool_steps} steps', flush=True)
        relaxed = M.prerelax(xml, initial, prere_cfg, cfg['platform'], cfg['precision'])
        context.setState(initial)
        context.setPositions(relaxed * 0.1)
        context.applyConstraints(1e-6)
        context.setVelocitiesToTemperature(cooling['temperature_start_K'] * u.kelvin,
                                           int(cfg['base_seed']) + 1000 + index)
        integrator.step(equil_steps)
        t_start = float(cooling['temperature_start_K'])
        t_end = float(cooling['temperature_end_K'])
        step_300 = int(round(cool_steps * (t_start - 300.0) / (t_start - t_end)))
        log_path = directory / 'cooling.csv'
        with open(log_path, 'w', newline='') as fh:
            writer = csv.writer(fh)
            writer.writerow(['step', 'time_ns', 'target_temperature_K', 'measured_temperature_K',
                             'volume_A3', 'density_g_cm3', 'specific_volume_cm3_g'])

            def record(done):
                state = context.getState(getEnergy=True, getPositions=True)
                density = M.density_g_cm3(system, state)
                target = t_start + (t_end - t_start) * done / cool_steps
                writer.writerow([done, done * dt_fs / 1e6, target,
                                 M.measured_temperature(context, system), M.volume_A3(state), density,
                                 1.0 / density])
                fh.flush()

            record(0)
            done = 0
            saved_300 = False
            while done < cool_steps:
                following = min(done + update, cool_steps)
                if done < step_300 < following:
                    following = step_300
                middle = done + 0.5 * (following - done)
                target = t_start + (t_end - t_start) * middle / cool_steps
                integrator.setTemperature(target * u.kelvin)
                context.setParameter(mm.MonteCarloBarostat.Temperature(), target)
                integrator.step(following - done)
                done = following
                if done % sample == 0 or done in (step_300, cool_steps):
                    record(done)
                if done == step_300 and not saved_300:
                    M.save_state(context, directory / 'state_300K.xml')
                    saved_300 = True
            if not saved_300:
                raise RuntimeError('cooling schedule did not produce the required 300 K state')
        final = M.save_state(context, directory / 'state_200K.xml')
        actual_duration_ns = cool_steps * dt_fs / 1e6
        summary = {
            'replica': index + 1,
            'atoms': system.getNumParticles(),
            'constraints': constraints,
            'cooling_steps': cool_steps,
            'nominal_duration_ns': duration_ns,
            'executed_duration_ns': actual_duration_ns,
            'nominal_rate_K_per_ns': cooling['rate_K_per_ns'],
            'executed_rate_K_per_ns': ((t_start - t_end) / actual_duration_ns
                                       if actual_duration_ns else None),
            'test_mode': bool(cfg.get('test_mode')),
            'final_density_g_cm3': M.density_g_cm3(system, final),
        }
        dump_json(summary_path, summary)
        summaries.append(summary)
        del context, integrator
    dump_json(out / 'summary.json', {'stage': 'cool', 'replicas': summaries,
                                     'test_mode': bool(cfg.get('test_mode'))})


def scale_axis(context, axis, factor):
    import openmm as mm
    import openmm.unit as u
    state = context.getState(getPositions=True, enforcePeriodicBox=False)
    pos = np.asarray(state.getPositions(asNumpy=True).value_in_unit(u.nanometer))
    box = np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.nanometer))
    pos[:, axis] *= factor
    box[axis] *= factor
    context.setPeriodicBoxVectors(*[mm.Vec3(*row) for row in box])
    context.setPositions(pos)
    context.applyConstraints(1e-6)


def stage_elastic(project, cfg, force=False):
    import openmm as mm
    import openmm.unit as u
    import mdtools as M

    build_dir = project / 'runs' / STAGES['build']
    cool_dir = project / 'runs' / STAGES['cool']
    out = project / 'runs' / STAGES['elastic']
    out.mkdir(parents=True, exist_ok=True)
    md = cfg['md']
    elastic = cfg['elastic']
    dt_fs = float(md['dt_fs'])
    production_steps = elastic['duration_ns'] * 1e6 / dt_fs
    tensile_steps = test_steps(cfg, production_steps, 'elastic')
    equil_steps = test_steps(cfg, elastic['equilibrate_ns'] * 1e6 / dt_fs, 'equilibrate')
    deform_every = 1 if cfg.get('test_mode') else int(elastic['deform_every_steps'])
    sample_every = 2 if cfg.get('test_mode') else int(elastic['sample_every_steps'])
    barostat_every = 1 if cfg.get('test_mode') else int(md['barostat_every_steps'])
    summaries = []
    for index in range(int(cfg['replicates'])):
        directory = out / rep_name(index)
        summary_path = directory / 'summary.json'
        if summary_path.exists() and not force:
            summaries.append(json.loads(summary_path.read_text()))
            continue
        directory.mkdir(parents=True, exist_ok=True)
        # First equilibrate once at 300 K and 1 atm in all three directions.
        # Every requested tensile axis then branches from exactly this state.
        equil_system = M.load_xml(build_dir / rep_name(index) / 'system.openmm_system.xml')
        M.configure(equil_system, md['shake_hydrogen'])
        pressure = md['pressure_atm'] * M.ATM_TO_BAR * u.bar
        equil_system.addForce(mm.MonteCarloBarostat(
            pressure, elastic['temperature_K'] * u.kelvin, barostat_every))
        equil_integrator = mm.LangevinMiddleIntegrator(
            elastic['temperature_K'] * u.kelvin, md['friction_per_ps'] / u.picosecond,
            dt_fs * u.femtosecond)
        equil_context = M.make_context(equil_system, equil_integrator, cfg['platform'], cfg['precision'])
        equil_context.setState(M.load_xml(cool_dir / rep_name(index) / 'state_300K.xml'))
        equil_context.setVelocitiesToTemperature(elastic['temperature_K'] * u.kelvin,
                                                 int(cfg['base_seed']) + 2000 + index)
        print(f'[elastic] {rep_name(index)}: isotropic 300 K / 1 atm re-equilibration '
              f'{equil_steps} steps', flush=True)
        equil_integrator.step(equil_steps)
        equilibrated_path = directory / 'state_300K_equilibrated.xml'
        M.save_state(equil_context, equilibrated_path)
        del equil_context, equil_integrator

        per_axis = []
        for axis_name in elastic['axes']:
            axis = AXIS[axis_name]
            system = M.load_xml(build_dir / rep_name(index) / 'system.openmm_system.xml')
            M.configure(system, md['shake_hydrogen'])
            pressure = md['pressure_atm'] * M.ATM_TO_BAR * u.bar
            flags = [True, True, True]
            flags[axis] = False
            barostat = mm.MonteCarloAnisotropicBarostat(
                mm.Vec3(pressure, pressure, pressure), elastic['temperature_K'] * u.kelvin,
                flags[0], flags[1], flags[2], barostat_every)
            system.addForce(barostat)
            integrator = mm.LangevinMiddleIntegrator(elastic['temperature_K'] * u.kelvin,
                                                     md['friction_per_ps'] / u.picosecond,
                                                     dt_fs * u.femtosecond)
            context = M.make_context(system, integrator, cfg['platform'], cfg['precision'])
            context.setState(M.load_xml(equilibrated_path))
            print(f'[elastic] {rep_name(index)} {axis_name}: strain '
                  f'{elastic["max_strain"]:.4f} over {tensile_steps} steps', flush=True)
            path = directory / f'stress_{axis_name}.csv'
            with open(path, 'w', newline='') as fh:
                writer = csv.writer(fh)
                writer.writerow(['step', 'time_ns', 'engineering_strain', 'tensile_stress_MPa',
                                 'volume_A3', 'density_g_cm3'])

                def record(done):
                    state = context.getState(getPositions=True)
                    strain = elastic['max_strain'] * done / tensile_steps
                    stress = M.axial_stress_MPa(context, system, axis, elastic['stress_delta'],
                                                elastic['temperature_K'])
                    writer.writerow([done, done * dt_fs / 1e6, strain, stress,
                                     M.volume_A3(state), M.density_g_cm3(system, state)])
                    fh.flush()

                record(0)
                done = 0
                old_strain = 0.0
                while done < tensile_steps:
                    following = min(done + deform_every, tensile_steps)
                    new_strain = elastic['max_strain'] * following / tensile_steps
                    scale_axis(context, axis, (1.0 + new_strain) / (1.0 + old_strain))
                    integrator.step(following - done)
                    done = following
                    old_strain = new_strain
                    if done % sample_every == 0 or done == tensile_steps:
                        record(done)
            M.save_state(context, directory / f'state_{axis_name}_2pct.xml')
            per_axis.append({'axis': axis_name, 'steps': tensile_steps,
                             'nominal_duration_ns': elastic['duration_ns'],
                             'executed_duration_ns': tensile_steps * dt_fs / 1e6,
                             'max_strain': elastic['max_strain']})
            del context, integrator
        summary = {'replica': index + 1, 'equilibrate_300K_steps': equil_steps,
                   'equilibrate_300K_ns': equil_steps * dt_fs / 1e6,
                   'axes': per_axis, 'test_mode': bool(cfg.get('test_mode'))}
        dump_json(summary_path, summary)
        summaries.append(summary)
    dump_json(out / 'summary.json', {'stage': 'elastic', 'replicas': summaries,
                                     'test_mode': bool(cfg.get('test_mode'))})


def stage_analyze(project, cfg, force=False):
    import analysis as A

    out = project / 'runs' / STAGES['analyze']
    out.mkdir(parents=True, exist_ok=True)
    cooling_dir = project / 'runs' / STAGES['cool']
    elastic_dir = project / 'runs' / STAGES['elastic']
    thermal, mechanical, errors = [], [], []
    for index in range(int(cfg['replicates'])):
        name = rep_name(index)
        rows = A.read_csv(cooling_dir / name / 'cooling.csv')
        try:
            fit = A.fit_tg([float(r['target_temperature_K']) for r in rows],
                           [float(r['density_g_cm3']) for r in rows],
                           cfg['cooling']['fit_bin_width_K'])
            fit['replica'] = index + 1
            thermal.append(fit)
        except Exception as exc:
            errors.append(f'{name} thermal: {exc}')
        axis_fits = []
        for axis in cfg['elastic']['axes']:
            rows = A.read_csv(elastic_dir / name / f'stress_{axis}.csv')
            try:
                fit = A.fit_modulus([float(r['engineering_strain']) for r in rows],
                                    [float(r['tensile_stress_MPa']) for r in rows],
                                    cfg['elastic']['fit_min_strain'], cfg['elastic']['fit_max_strain'])
                fit.update({'replica': index + 1, 'axis': axis})
                mechanical.append(fit)
                axis_fits.append(fit['tensile_modulus_GPa'])
            except Exception as exc:
                errors.append(f'{name} elastic {axis}: {exc}')
        if axis_fits:
            mechanical.append({'replica': index + 1, 'axis': 'axis_mean',
                               'tensile_modulus_GPa': float(np.mean(axis_fits)),
                               'tensile_modulus_MPa': float(np.mean(axis_fits) * 1000)})
    if errors and not cfg.get('test_mode'):
        raise RuntimeError('analysis failed:\n  ' + '\n  '.join(errors))
    thermal_keys = {
        'tg_K': [r['tg_K'] for r in thermal],
        'alpha_linear_glass_1_K': [r['alpha_linear_glass_1_K'] for r in thermal],
        'alpha_linear_rubber_1_K': [r['alpha_linear_rubber_1_K'] for r in thermal],
    }
    axis_mean = [r['tensile_modulus_GPa'] for r in mechanical if r['axis'] == 'axis_mean']
    aggregate = {key: A.statistics(value) for key, value in thermal_keys.items()}
    aggregate['tensile_modulus_GPa'] = A.statistics(axis_mean)
    result = {
        'name': cfg['name'],
        'test_mode': bool(cfg.get('test_mode')),
        'replicates_requested': int(cfg['replicates']),
        'aggregate': aggregate,
        'thermal_replicates': thermal,
        'elastic_fits': mechanical,
        'errors': errors,
        'interpretation': {
            'tg': 'centre/intersection of a smooth hyperbola fitted to density versus target temperature',
            'linear_cte': '-(1/(3 rho_Tg)) d(rho)/dT from the low- and high-temperature density asymptotes',
            'tensile_modulus': 'OLS slope of axial tensile stress versus engineering strain; transverse axes at 1 atm',
        },
    }
    dump_json(out / 'results.json', result)
    thermal_rows = [{k: v for k, v in row.items() if k != 'bins'} for row in thermal]
    if thermal_rows:
        A.write_rows(out / 'thermal_replicates.csv', list(thermal_rows[0]), thermal_rows)
    simple_mechanical = [{k: v for k, v in row.items() if k not in ('bins',)} for row in mechanical]
    if simple_mechanical:
        fields = sorted(set().union(*(row.keys() for row in simple_mechanical)))
        A.write_rows(out / 'elastic_fits.csv', fields, simple_mechanical)
    dump_json(out / 'summary.json', {'stage': 'analyze', 'aggregate': aggregate,
                                     'errors': errors, 'results': 'results.json'})
    print('[analyze]', json.dumps(aggregate), flush=True)
    if cfg.get('test_mode'):
        print('[analyze] TEST MODE: fitted values are workflow diagnostics, not physical results.', flush=True)


RUNNERS = {'build': stage_build, 'cool': stage_cool, 'elastic': stage_elastic, 'analyze': stage_analyze}


def cmd_new(args):
    project = Path(args.dir or Path.cwd() / 'projects') / args.name
    target = project / 'project.json'
    if target.exists() and not args.force:
        sys.exit(f'{target} exists; use --force to overwrite it')
    if args.atoms_per_chain < 1 or args.total_atoms < 1 or args.replicates < 1 or args.density <= 0:
        sys.exit('atom targets and replicas must be positive integers, and density must be positive')
    cfg = json.loads((HERE / 'defaults.json').read_text())
    cfg.pop('_comment', None)
    cfg['name'] = args.name
    try:
        if args.polypaves:
            if args.comonomer or args.architecture or args.fraction:
                raise ValueError('--polypaves cannot be combined with copolymer CLI options')
            lines, file_sequence_seed = chemistry_from_file(args.polypaves)
            cfg['system']['chemistry_lines'] = lines
            cfg['system']['polypaves_source'] = str(Path(args.polypaves).resolve())
            sequence_seed = (args.sequence_seed if args.sequence_seed is not None else
                             (file_sequence_seed if file_sequence_seed is not None else 7))
            chemistry_description = Path(args.polypaves).name
        else:
            lines = chemistry_from_cli(args)
            sequence_seed = args.sequence_seed if args.sequence_seed is not None else 7
            if args.comonomer:
                cfg['system']['chemistry_lines'] = lines
                chemistry_description = ARCHITECTURE_ALIASES.get(args.architecture, args.architecture)
            else:
                cfg['system'].update({'monomer': args.monomer, 'terminator': args.terminator})
                chemistry_description = args.monomer
    except (ValueError, OSError) as exc:
        sys.exit(f'invalid chemistry input: {exc}')
    cfg['system'].update({
        'sequence_seed': sequence_seed,
        'target_atoms_per_chain': args.atoms_per_chain,
        'target_total_atoms': args.total_atoms,
        'initial_density_g_cm3': args.density,
    })
    cfg['replicates'] = args.replicates
    cfg['base_seed'] = args.seed
    cfg['platform'] = args.platform
    if args.forcefield:
        cfg['forcefield'] = FORCEFIELDS.get(args.forcefield, str(Path(args.forcefield).resolve()))
    cfg['test_mode'] = bool(args.test)
    project.mkdir(parents=True, exist_ok=True)
    dump_json(target, cfg)
    print(f'created {target}')
    print(f'  chemistry: {chemistry_description}')
    print(f'  {args.replicates} replicas; >= {args.atoms_per_chain} atoms/chain; '
          f'>= {args.total_atoms} atoms/system; density {args.density} g/cm3')
    if args.test:
        print('  TEST MODE: minimal integration steps; results will not be physical')
    print(f'next: python {Path(__file__).name} run {project}')


def cmd_status(args):
    project, cfg = load_project(args.project)
    print(f'{cfg["name"]} ({project})')
    for name, directory in STAGES.items():
        summary = project / 'runs' / directory / 'summary.json'
        if summary.exists():
            print(f'  done  {name:8s} {directory}')
        elif summary.parent.exists():
            count = sum((summary.parent / rep_name(i) / 'summary.json').exists()
                        for i in range(int(cfg['replicates'])))
            print(f'  part  {name:8s} {directory} ({count}/{cfg["replicates"]} replicas)')
        else:
            print(f'  -     {name:8s} {directory}')


def cmd_run(args):
    project, cfg = load_project(args.project)
    names = list(STAGES)
    for value in (args.from_, args.to, args.only):
        if value and value not in names:
            sys.exit(f'unknown stage {value!r}; choose from {", ".join(names)}')
    if args.only:
        todo = [args.only]
        explicit = True
    else:
        first = (names.index(args.from_) if args.from_ else next(
            (i for i, name in enumerate(names)
             if not (project / 'runs' / STAGES[name] / 'summary.json').exists()), len(names)))
        last = names.index(args.to) if args.to else len(names) - 1
        todo = names[first:last + 1]
        explicit = bool(args.from_)
    if not todo:
        print('all stages are complete; use --from STAGE to rerun')
        return
    start = names.index(todo[0])
    if start and not (project / 'runs' / STAGES[names[start - 1]] / 'summary.json').exists():
        sys.exit(f'{todo[0]} requires completed stage {names[start - 1]}')
    print(f'{cfg["name"]}: running {", ".join(todo)}')
    for name in todo:
        before = time.time()
        print(f'=== {name} {time.strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
        RUNNERS[name](project, copy.deepcopy(cfg), force=explicit or bool(args.only))
        print(f'=== {name} complete in {(time.time() - before) / 60:.2f} min', flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command', required=True)
    new = sub.add_parser('new', help='create a project')
    new.add_argument('name')
    source = new.add_mutually_exclusive_group(required=True)
    source.add_argument('--monomer', help="primary repeat-unit SMILES (A) with two '*' connection atoms")
    source.add_argument('--polypaves', metavar='FILE', help='PAVES chemistry file; its size/run records are replaced')
    new.add_argument('--comonomer', action='append', metavar='NAME=SMILES',
                     help='additional repeat unit; repeat for B, C, ...')
    new.add_argument('--architecture',
                     choices=['alternating', 'random', 'random-exact', 'random-probabilistic',
                              'block', 'blocks'],
                     help='copolymer architecture (random means random-exact)')
    new.add_argument('--fraction', action='append', metavar='NAME=VALUE',
                     help='monomer fraction for random/block; one per monomer and sum=1')
    new.add_argument('--sequence-seed', type=int,
                     help='base sequence seed; replica i uses this + i (default file value or 7)')
    new.add_argument('--atoms-per-chain', required=True, type=int, help='minimum all-atom count in each chain')
    new.add_argument('--total-atoms', required=True, type=int, help='minimum all-atom count in the system')
    new.add_argument('--terminator', default='*C', help="end group (default '*C', methyl)")
    new.add_argument('--density', type=float, default=0.5, help='initial density in g/cm3 (default 0.5)')
    new.add_argument('--replicates', type=int, default=4, help='independent seeds (default 4)')
    new.add_argument('--seed', type=int, default=20260929)
    new.add_argument('--forcefield', help='PAVES force-field name or .ff path (default PCFF)')
    new.add_argument('--platform', choices=['auto', 'CPU', 'CUDA', 'OpenCL', 'Reference'], default='auto')
    new.add_argument('--dir', help='project parent directory (default ./projects)')
    new.add_argument('--test', action='store_true', help='minimal CPU-friendly step counts; not physical')
    new.add_argument('--force', action='store_true')
    run = sub.add_parser('run', help='run or resume stages')
    run.add_argument('project')
    run.add_argument('--from', dest='from_', metavar='STAGE')
    run.add_argument('--to', metavar='STAGE')
    run.add_argument('--only', metavar='STAGE')
    status = sub.add_parser('status', help='show completed stages')
    status.add_argument('project')
    args = parser.parse_args(argv)
    {'new': cmd_new, 'run': cmd_run, 'status': cmd_status}[args.command](args)


if __name__ == '__main__':
    main()
