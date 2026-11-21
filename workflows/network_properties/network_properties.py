"""Monomer -> crosslinked network property workflow (epoxy, oxetane, ...).

monomer   pack and relax the monomer at T, thermal history (Tg), 300 K properties,
          then T properties and the relaxed T structure.
crosslink one project per start-site density under <monomer>/crosslinks/; replica N
          crosslinks at T from monomer replica N, then thermal history (Tg) and 300 K.
network   a saved dry network evaluated like a crosslink project.
mechanical  optional elastic constants of a finished crosslink/network replica.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from workflows.network_properties.analysis import analyze_npt, analyze_nvt, aggregate, topology, DA
from workflows.thermomechanical.analysis import fit_tg
from workflows.network_properties.execution import (check_run_signature, configuration_fingerprint,
                                                    file_lock, replica_progress, replica_seed, replica_seeds,
                                                    set_stage)


# Thermal history, 300 K NPT/NVT and solvents use H-bond constraints and 1 fs.
# The elastic branch uses its own unconstrained 0.25 fs (mechanical.py).
DEFAULTS = dict(replicates=4, seed=31415, platform='CUDA', timestep_fs=1.,
    pressure_atm=1., high_K=800., low_K=100., temperature_step_K=20.,
    heat_step_ps=100., high_equilibration_ps=1000., cooling_rate_K_ns=40.,
    reheat_step_ps=100., equilibration_300_ps=1000., production_ps=10000.,
    sample_interval_ps=1., min_effective_samples=100, minimization_iterations=3000,
    volume_interval_steps=50, volume_fraction=.002, solvents=[],
    solvent_molecules=10, solvent_nvt_ps=100., solvent_equilibration_ps=1000.,
    solvent_production_ps=10000., diffusion_min_duration_ps=1000.,
    alchemical_windows=22, alchemical_equilibration_ps=100.,
    alchemical_sampling_ps=500., standard_density_nm3=.602214076,
    nvt_production_ps=10000.,
    elastic_strain=.01, elastic_equilibration_ps=200., elastic_production_ps=1000.,
    cte_offset_K=50., smoke=False)
# Version of the cooling fit stored in dry_result.json; report/collect refit older ones from cooling.csv.
COOLING_FIT_VERSION = 2


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    def convert(x):
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, np.generic):
            return x.item()
        raise TypeError(type(x).__name__)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, default=convert, allow_nan=False) + '\n')
    temporary.replace(path)


def write_csv(path, rows):
    if not rows:
        return
    with Path(path).open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


MODES = ('monomer', 'crosslink', 'saved_network')


def mode(cfg):
    return cfg.get('input_mode', 'saved_network')


def input_digest(project, cfg):
    if mode(cfg) in ('monomer', 'crosslink'):
        from workflows.network_properties.construction import definition_digest
        return definition_digest(project)
    return digest(Path(project) / cfg['network'])


def children(project):
    """Crosslink projects of a monomer project, sorted by name."""
    return sorted(path.parent for path in (Path(project) / 'crosslinks').glob('*/project.json'))


def collection_root(project, cfg):
    """Directory whose project comparison includes this project: the parent of the monomer project."""
    if mode(cfg) == 'crosslink':
        from workflows.network_properties.construction import parent_project
        return parent_project(project, cfg).parent
    return Path(project).resolve().parent


def validate(cfg):
    for key in ('replicates', 'seed', 'solvent_molecules', 'alchemical_windows',
                'min_effective_samples', 'minimization_iterations', 'volume_interval_steps'):
        if type(cfg[key]) is not int or cfg[key] < 1:
            raise ValueError(f'{key} must be a positive integer')
    if cfg['alchemical_windows'] < 4:
        raise ValueError('alchemical_windows must be >= 4')
    for key in ('timestep_fs', 'pressure_atm', 'high_K', 'low_K', 'temperature_step_K',
                'heat_step_ps', 'high_equilibration_ps', 'cooling_rate_K_ns',
                'reheat_step_ps', 'equilibration_300_ps', 'production_ps', 'sample_interval_ps',
                'solvent_nvt_ps', 'solvent_equilibration_ps', 'solvent_production_ps',
                'diffusion_min_duration_ps', 'alchemical_equilibration_ps',
                'alchemical_sampling_ps', 'standard_density_nm3', 'volume_fraction',
                'nvt_production_ps', 'elastic_equilibration_ps', 'elastic_production_ps'):
        if isinstance(cfg[key], bool) or not math.isfinite(cfg[key]) or cfg[key] <= 0:
            raise ValueError(f'{key} must be finite and positive')
    if not cfg['low_K'] < 300 < cfg['high_K']:
        raise ValueError('temperature range must bracket 300 K')
    if mode(cfg) not in MODES:
        raise ValueError(f'input_mode must be one of {MODES}')
    # Projects created before cte_offset_K existed use the default.
    offset = cfg.get('cte_offset_K', DEFAULTS['cte_offset_K'])
    if isinstance(offset, bool) or not isinstance(offset, (int, float)) or not math.isfinite(offset) or offset <= 0:
        raise ValueError('cte_offset_K must be finite and positive')
    if cfg['volume_fraction'] >= .1:
        raise ValueError('volume_fraction must be < 0.1')
    if cfg['production_ps'] / cfg['sample_interval_ps'] < 20:
        raise ValueError('production must have at least 20 samples')
    if cfg['sample_interval_ps'] < cfg['timestep_fs'] / 1000:
        raise ValueError('sample interval is shorter than a time step')
    if not 0 < cfg['elastic_strain'] < .03:
        raise ValueError('elastic_strain must be > 0 and < 0.03')
    names = set()
    for solvent in cfg['solvents']:
        name, smiles = solvent['name'], solvent['smiles']
        if not name or not all(c.isalnum() or c in '_-' for c in name) or name in names:
            raise ValueError('solvent names must be unique, using letters, digits, _ or -')
        if not smiles:
            raise ValueError('empty solvent SMILES')
        names.add(name)
    if mode(cfg) in ('monomer', 'crosslink'):
        from workflows.network_properties.construction import CHEMISTRIES
        if cfg.get('chemistry') not in CHEMISTRIES:
            raise ValueError(f'chemistry must be one of {sorted(CHEMISTRIES)}')
        temperature = cfg.get('temperature_K')
        if isinstance(temperature, bool) or not isinstance(temperature, (int, float)) or not temperature > 0:
            raise ValueError('temperature_K must be positive')
    if mode(cfg) == 'monomer':
        from workflows.network_properties.construction import validate_monomer
        if cfg['solvents']:
            raise ValueError('solvents belong to the crosslink projects, not the monomer project')
        if not cfg['low_K'] < cfg['temperature_K'] < cfg['high_K']:
            raise ValueError('the thermal history range must bracket temperature_K')
        for key in ('t_equilibration_ps', 't_production_ps', 't_nvt_production_ps'):
            if isinstance(cfg[key], bool) or not math.isfinite(cfg[key]) or cfg[key] <= 0:
                raise ValueError(f'{key} must be finite and positive')
        if cfg['t_production_ps'] / cfg['sample_interval_ps'] < 20:
            raise ValueError('the T production must have at least 20 samples')
        cutoffs = cfg['ring_contact_cutoffs_A']
        if not cutoffs or not all(isinstance(c, (int, float)) and not isinstance(c, bool) and c > 0 for c in cutoffs):
            raise ValueError('ring_contact_cutoffs_A must be a list of positive distances')
        validate_monomer(cfg)
    elif mode(cfg) == 'crosslink':
        from workflows.network_properties.construction import validate_crosslink
        validate_crosslink(cfg)


def relaxation_stages_argument(value):
    # Imported here so that status/submit do not load polypaves for argument parsing.
    from workflows.crosslinking.crosslinking import load_relaxation_stages
    return load_relaxation_stages(value)


def base_config(args, smoke):
    cfg = json.loads(json.dumps(DEFAULTS))
    cfg.update(platform=args.platform, seed=args.seed)
    if getattr(args, 'high_temperature', None) is not None:
        cfg['high_K'] = args.high_temperature
    for entry in getattr(args, 'solvent', None) or []:
        name, sep, smiles = entry.partition('=')
        if not sep:
            raise ValueError('--solvent expects NAME=SMILES')
        cfg['solvents'].append(dict(name=name, smiles=smiles))
    if getattr(args, 'solvent_molecules', None) is not None:
        cfg['solvent_molecules'] = args.solvent_molecules
    if smoke:
        cfg.update(smoke=True, timestep_fs=.25,
            heat_step_ps=.002, high_equilibration_ps=.002,
            cooling_rate_K_ns=10_000_000., reheat_step_ps=.002,
            equilibration_300_ps=.002, production_ps=.04, sample_interval_ps=.001,
            minimization_iterations=100, volume_interval_steps=1,
            solvent_nvt_ps=.002, solvent_equilibration_ps=.002, solvent_production_ps=.02,
            alchemical_windows=4, alchemical_equilibration_ps=.001, alchemical_sampling_ps=.004,
            nvt_production_ps=.02)
    return cfg


def _create(project, cfg, definitions=None, network=None):
    validate(cfg)
    project.mkdir(parents=True)
    if definitions is not None:
        from workflows.network_properties.construction import write_definitions
        write_definitions(project, definitions)
    if network is not None:
        (project / cfg['network']).write_bytes(network.read_bytes())
    write_json(project / 'project.json', cfg)
    print(f'Created {project / "project.json"}', flush=True)


def new_monomer(args):
    """Monomer project: build at T, thermal history, 300 K and T properties."""
    from workflows.network_properties import monomer
    from workflows.network_properties.construction import (definitions_for, definition_digest, monomer_options,
                                                           reactive_sites)
    project = args.project
    if project.exists():
        raise FileExistsError(f'Project exists: {project}')
    cfg = base_config(args, args.smoke)
    cfg.update(replicates=args.replicates if args.replicates is not None else (1 if args.smoke else DEFAULTS['replicates']))
    chemistry = args.chemistry or 'epoxy'
    cfg.update(input_mode='monomer', schema_version=3, chemistry=chemistry, **json.loads(json.dumps(monomer.DEFAULTS)))
    cfg['temperature_K'] = float(args.temperature)
    if args.smoke:
        cfg.update(monomer.SMOKE)
    cfg['build'] = monomer_options(args, cfg['temperature_K'])
    definitions = definitions_for(chemistry, args.activation_yaml, args.reactions_yaml)
    cfg['reactive_sites_per_molecule'] = reactive_sites(cfg['build']['smiles'], definitions['activation.yaml'])
    _create(project, cfg, definitions)
    cfg['input_sha256'] = definition_digest(project)
    write_json(project / 'project.json', cfg)
    from workflows.network_properties.directory_summary import enable_collection
    enable_collection(project.resolve().parent)
    print(f'Next: run {project} (or submit), then crosslink {project} --activation-fraction ...', flush=True)


def crosslink_name(count=None, fraction=None):
    return f'n{count}' if count is not None else f'af{fraction:g}'


def new_crosslink(args):
    """One crosslink project per start-site density under <monomer>/crosslinks/."""
    from workflows.network_properties.construction import (crosslink_options, definition_digest, with_selection,
                                                           DEFINITIONS)
    from workflows.network_properties.execution import configuration_fingerprint
    parent = args.monomer.resolve()
    parent_cfg = json.loads((parent / 'project.json').read_text())
    if mode(parent_cfg) != 'monomer':
        raise ValueError(f'{parent} is not a monomer project')
    validate(parent_cfg)
    values = ([dict(count=v) for v in args.activation_count] if args.activation_count is not None
              else [dict(fraction=v) for v in args.activation_fraction])
    if args.name and len(values) != 1:
        raise ValueError('--name names one density; omit it for several')
    if not parent_cfg['reactive_sites_per_molecule'] and any(v.get('count', v.get('fraction')) for v in values):
        raise ValueError(f"the monomer has no site matching the activation target_smarts; only 0 starts are possible")
    replicates = args.replicates if args.replicates is not None else parent_cfg['replicates']
    if not 1 <= replicates <= parent_cfg['replicates']:
        raise ValueError(f'replicates must lie between 1 and the monomer replicates ({parent_cfg["replicates"]}): '
                         'crosslink replica N starts from monomer replica N')
    plans = []
    for value in values:
        project = parent / 'crosslinks' / (args.name or crosslink_name(**value))
        if project.exists():
            raise FileExistsError(f'Project exists: {project}')
        cfg = base_config(args, parent_cfg['smoke'])
        cfg.update(platform=args.platform or parent_cfg['platform'],
                   seed=args.seed if args.seed is not None else parent_cfg['seed'], replicates=replicates)
        cfg.update(input_mode='crosslink', schema_version=3, chemistry=parent_cfg['chemistry'],
                   temperature_K=parent_cfg['temperature_K'], parent='../..',
                   parent_configuration_sha256=configuration_fingerprint(parent_cfg))
        cfg['crosslinking'] = crosslink_options(args, temperature=parent_cfg['temperature_K'],
                                                smoke=parent_cfg['smoke'], **value)
        definitions = with_selection({name: (parent / 'definitions' / name).read_bytes() for name in DEFINITIONS},
                                     **value)
        validate(cfg)
        plans.append((project, cfg, definitions))
    for project, cfg, definitions in plans:
        _create(project, cfg, definitions)
        cfg['input_sha256'] = definition_digest(project)
        write_json(project / 'project.json', cfg)
    print(f'Next: run or submit {parent} --all (or each crosslink directory)', flush=True)


def new_network(args):
    """A saved dry network: thermal history (Tg), 300 K properties and solvents."""
    project = args.project
    if project.exists():
        raise FileExistsError(f'Project exists: {project}')
    cfg = base_config(args, args.smoke)
    cfg.update(replicates=args.replicates if args.replicates is not None else (1 if args.smoke else DEFAULTS['replicates']))
    source = args.network.resolve()
    data = json.loads(source.read_text())
    if not all(k in data for k in ('graphs', 'positions', 'box_vectors')):
        raise ValueError('network must be a PAVES saved-charge structure JSON')
    if data.get('solvent_molecules', 0):
        raise ValueError('input must be a dry network')
    if not all('partial_charge' in a for g in data['graphs'] for a in g['atoms']):
        raise ValueError('all network atoms must have saved partial charges')
    cfg.update(network='input.structure.json', source=str(source), input_sha256=digest(source),
               input_mode='saved_network', schema_version=3)
    _create(project, cfg, network=source)
    from workflows.network_properties.directory_summary import enable_collection
    enable_collection(project.resolve().parent)


def steps(ps, cfg):
    return max(1, int(round(ps * 1000 / cfg['timestep_fs'])))


def set_temperature(runner, temperature):
    from openmm import unit
    runner.temperature = float(temperature)
    runner.integrator.setTemperature(temperature * unit.kelvin)


def phase(runner, directory, cfg, temperature, duration_ps, *, sample=False, minimize=False, base=None, observe=None):
    """Completed phases resume from portable state; interrupted phases rerun.

    base: replica directory for progress records (default: the phase's parent).
    observe: with sample, an object whose add(runner) sees every sample and
    finish(directory, duration_ps) runs before the phase is marked done.
    """
    import openmm as mm
    from openmm import unit
    directory.mkdir(parents=True, exist_ok=True)
    marker = directory / 'done.json'
    set_temperature(runner, temperature)
    if marker.exists():
        runner.context.setState(mm.XmlSerializer.deserialize((directory / 'state.xml').read_text()))
        return json.loads((directory / 'samples.json').read_text())
    print(f'{directory.name}: {temperature:g} K, {duration_ps:g} ps', flush=True)
    stage_base = base or directory.parent
    where = '' if stage_base == directory.parent else f'{directory.parent.name}/'
    set_stage(stage_base, f'thermal history: {where}{directory.name}',
              temperature_K=float(temperature), duration_ps=float(duration_ps))
    if minimize:
        runner.minimize(cfg['minimization_iterations'])
    rows = []
    remaining, elapsed = steps(duration_ps, cfg), 0
    interval = steps(cfg['sample_interval_ps'], cfg)
    mass = sum(runner.system.getParticleMass(i).value_in_unit(unit.dalton)
               for i in range(runner.system.getNumParticles()))
    while remaining:
        chunk = min(remaining, interval)
        runner.md(chunk, ensemble='NPT')
        remaining -= chunk
        elapsed += chunk
        if sample:
            state = runner.context.getState(getEnergy=True)
            volume = state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3)
            rows.append(dict(time_ps=elapsed * cfg['timestep_fs'] / 1000,
                target_temperature_K=float(temperature), volume_nm3=float(volume),
                density_g_cm3=float(mass * DA / (volume * 1e-27) / 1000),
                potential_kJ_mol=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
                kinetic_kJ_mol=state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)))
            if observe is not None:
                observe.add(runner)
    if observe is not None:
        observe.finish(directory, elapsed * cfg['timestep_fs'] / 1000)
    state = runner.context.getState(getPositions=True, getVelocities=True)
    volume = state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3)
    (directory / 'state.xml').write_text(mm.XmlSerializer.serialize(state))
    write_json(directory / 'final.structure.json', runner.structure_data())
    write_json(directory / 'samples.json', rows)
    write_csv(directory / 'samples.csv', rows)
    write_json(marker, dict(temperature_K=temperature, duration_ps=duration_ps,
        actual_duration_ps=elapsed * cfg['timestep_fs'] / 1000,
        final_density_g_cm3=float(mass * DA / (volume * 1e-27) / 1000),
        volume_attempts=runner.volume_attempts, volume_accepts=runner.volume_accepts))
    return rows


def temperature_grid(high, low, increment):
    # Include exactly 300 K even when the chosen step does not divide the range.
    values = list(np.arange(high, low, -increment)) + [low, 300.]
    return sorted(set(float(x) for x in values if low <= x <= high), reverse=True)


def history_densities(directory):
    """Density at the end of every thermal-history phase, and the high-temperature check.

    An uncrosslinked monomer or oligomer liquid can boil or cavitate at the
    high temperature (fixed topology, 1 atm). The voids then survive the cooling
    and corrupt the Tg/CTE fit and the 300 K fluctuations; a large density drop
    at the high temperature is the sign of it.
    """
    phases = []
    for marker in sorted(Path(directory).glob('0*/done.json')):
        done = json.loads(marker.read_text())
        phases.append(dict(phase=marker.parent.name, temperature_K=done['temperature_K'],
                           final_density_g_cm3=done.get('final_density_g_cm3'),
                           volume_acceptance=(done['volume_accepts'] / done['volume_attempts']
                                              if done.get('volume_attempts') else None)))
    density = {p['phase']: p['final_density_g_cm3'] for p in phases}
    start, high = density.get('00_minimize'), density.get('02_high_equilibration')
    ratio = high / start if start and high else None
    return dict(phases=phases, high_to_initial_density_ratio=ratio,
                # Thermal expansion alone gives about 0.8-0.85 for these resins from 300 to 800 K.
                possible_vaporization=bool(ratio is not None and ratio < .7))


def nvt_properties(runner, directory, cfg, mass_Da, *, temperature=300., equilibration_ps=None,
                   production_ps=None, base=None):
    from openmm import unit
    directory.mkdir(exist_ok=True)
    if (directory / 'result.json').exists():
        return json.loads((directory / 'result.json').read_text())
    print(f'{temperature:g} K NVT: direct Cv and dielectric response', flush=True)
    set_stage(base or directory.parent, f'{temperature:g} K NVT: direct Cv and dielectric response')
    set_temperature(runner, temperature)
    # Restart this branch from the completed NPT phase after an interruption.
    runner.md(steps(cfg['equilibration_300_ps'] if equilibration_ps is None else equilibration_ps, cfg),
              ensemble='NVT')
    charges = np.array([a['partial_charge'] for g in runner.host.graphs for a in g['atoms']])
    if abs(charges.sum()) > 1e-4:
        raise ValueError('dielectric estimator requires a neutral network')
    initial = runner.snapshot()
    previous, box = initial['positions'], initial['box']
    inverse = np.linalg.inv(box)
    moment, rows = np.zeros(3), []
    total = steps(cfg.get('nvt_production_ps', cfg['production_ps']) if production_ps is None else production_ps, cfg)
    elapsed = 0
    while elapsed < total:
        chunk = min(steps(cfg['sample_interval_ps'], cfg), total-elapsed)
        runner.md(chunk, ensemble='NVT')
        elapsed += chunk
        state = runner.snapshot()
        displacement = (state['positions']-previous) @ inverse
        displacement -= np.rint(displacement)
        moment += charges @ (displacement @ box)
        previous = state['positions']
        rows.append(dict(time_ps=elapsed*cfg['timestep_fs']/1000,
            volume_nm3=float(np.linalg.det(box)), potential_kJ_mol=state['energy'],
            kinetic_kJ_mol=state['kinetic_energy'],
            **{f'dipole_{axis}_e_nm':float(value) for axis,value in zip('xyz',moment)}))
    result = analyze_nvt(rows, mass_Da, temperature_K=temperature, min_effective_samples=cfg['min_effective_samples'])
    result['temperature_K'] = float(temperature)
    write_json(directory / 'samples.json', rows)
    write_csv(directory / 'samples.csv', rows)
    write_json(directory / 'result.json', result)
    return result


def vaporization_reason(cfg, history):
    """Why the high-temperature density drop invalidates the run, or None."""
    if not history.get('possible_vaporization'):
        return None
    return (f'density at {cfg["high_K"]:g} K fell to {history["high_to_initial_density_ratio"]:.2f} of the start '
            '(boiling/cavitation); lower high_K (--high-temperature)')


def cooling_fit(temperatures, densities, cfg, history):
    """Tg/CTE fit of the cooling trace with its acceptance checks.

    CTEs are the fitted curve's own slope at Tg -/+ cte_offset_K. The checks
    also require both temperatures inside the cooled range (no extrapolation).
    """
    try:
        fit = fit_tg(temperatures, densities, bin_width_K=10.,
                     cte_offset_K=cfg.get('cte_offset_K', DEFAULTS['cte_offset_K']))
        fit['alpha_volume_glass_1_K'] = 3 * fit['alpha_linear_glass_1_K']
        fit['alpha_volume_rubber_1_K'] = 3 * fit['alpha_linear_rubber_1_K']
        t = np.array([b['temperature_K'] for b in fit['bins']])
        rho = np.array([b['density_g_cm3'] for b in fit['bins']])
        linear_sse = np.sum((rho - np.polyval(np.polyfit(t,rho,1),t))**2)
        curved_sse = len(t) * fit['rmse_g_cm3']**2
        fit['AIC_improvement_over_linear'] = float(len(t)*np.log(max(linear_sse,1e-30)/max(curved_sse,1e-30))-6)
        span = cfg['high_K'] - cfg['low_K']
        fit['fit_checks_passed'] = bool(not cfg['smoke'] and fit['r_squared'] is not None
            and fit['r_squared'] > .95 and fit['transition_width_K'] < span / 4
            and fit['AIC_improvement_over_linear'] > 6 and fit['alpha_linear_glass_1_K'] > 0
            and fit['cte_temperatures_within_data']
            and cfg['low_K'] + .1 * span < fit['tg_K'] < cfg['high_K'] - .1 * span)
    except (ValueError, RuntimeError) as exc:
        fit = dict(fit_checks_passed=False, reason=str(exc))
    reason = vaporization_reason(cfg, history)
    if reason:
        # Voids from boiling at high_K remain in the glass: the cooling curve is not the bulk material.
        fit.update(fit_checks_passed=False, vaporization_warning=reason)
    fit.update(fit_version=COOLING_FIT_VERSION, cte_offset_K=cfg.get('cte_offset_K', DEFAULTS['cte_offset_K']))
    return fit


def refresh_cooling(directory, cfg):
    """Refit an older cooling result from cooling.csv, rewrite dry_result.json and the Tg figure.

    report/collect call this, so finished runs get the current CTE definition
    without rerunning MD. Returns True when the result was updated.
    """
    directory = Path(directory)
    source, trace = directory / 'dry_result.json', directory / 'cooling.csv'
    if not source.is_file() or not trace.is_file():
        return False
    dry = json.loads(source.read_text())
    offset = cfg.get('cte_offset_K', DEFAULTS['cte_offset_K'])
    current = dry.get('cooling', {})
    if current.get('fit_version') == COOLING_FIT_VERSION and current.get('cte_offset_K') == offset:
        return False
    with trace.open() as handle:
        rows = list(csv.DictReader(handle))
    dry['cooling'] = cooling_fit([float(r['target_temperature_K']) for r in rows],
                                 [float(r['density_g_cm3']) for r in rows], cfg,
                                 dry.get('thermal_history_density') or {})
    write_json(source, dry)
    try:
        from workflows.network_properties.plots import plot_tg
        plot_tg(directory, dry)
    except Exception as exc:  # a figure must never fail the report
        print(f'WARNING: Tg figure not rewritten for {directory.name}: {exc}', flush=True)
    print(f'{directory.name}: cooling fit updated to version {COOLING_FIT_VERSION} '
          f'(CTE at Tg -/+ {offset:g} K)', flush=True)
    return True


def refresh_cooling_fits(project, cfg):
    return [refresh_cooling(Path(project) / f'replica_{i:02d}', cfg) for i in range(1, cfg['replicates'] + 1)]


def dry_properties(project, cfg, replica):
    from openmm import unit
    from polypaves.sorption import load_network
    from polypaves.sorption_workflow import FixedLoadingSimulation, structural_diagnostics
    directory = project / f'replica_{replica:02d}'
    directory.mkdir(exist_ok=True)
    if (directory / 'dry_result.json').exists():
        return json.loads((directory / 'dry_result.json').read_text())
    from workflows.network_properties.construction import network_path
    host = load_network(network_path(project, cfg, replica))
    # Count zero means no solvent is appended; reusing host avoids typing a dummy solvent.
    runner = FixedLoadingSimulation(host, host, host.positions_nm, 0, box_nm=host.box_nm,
        temperature=300., pressure_bar=cfg['pressure_atm'] * 1.01325,
        timestep_fs=cfg['timestep_fs'], platform=cfg['platform'], seed=replica_seed(project, cfg, replica),
        volume_interval=cfg['volume_interval_steps'], volume_fraction=cfg['volume_fraction'])
    try:
        phase(runner, directory / '00_minimize', cfg, 300., cfg['heat_step_ps'], minimize=True)
        for index, t in enumerate(np.arange(400., cfg['high_K'], 100.)):
            phase(runner, directory / f'01_heat_{index:03d}', cfg, t, cfg['heat_step_ps'])
        phase(runner, directory / '02_high_equilibration', cfg, cfg['high_K'], cfg['high_equilibration_ps'])
        cooling = []
        grid = temperature_grid(cfg['high_K'], cfg['low_K'], cfg['temperature_step_K'])
        previous = grid[0] + cfg['temperature_step_K']
        for index, temperature in enumerate(grid):
            duration = (previous - temperature) / cfg['cooling_rate_K_ns'] * 1000
            # Discard the first half of each plateau for temperature/volume settling.
            phase(runner, directory / f'03_cool_{index:03d}_settle', cfg, temperature, duration / 2)
            rows = phase(runner, directory / f'03_cool_{index:03d}_sample', cfg,
                         temperature, duration / 2, sample=True)
            cooling.extend(rows)
            if temperature == 300.:
                write_json(directory / 'cooling_300K.structure.json', runner.structure_data())
            previous = temperature
        write_json(directory / 'relaxed_100K.structure.json', runner.structure_data())
        write_csv(directory / 'cooling.csv', cooling)
        for index, t in enumerate(sorted(set(np.arange(cfg['low_K'] + cfg['temperature_step_K'],
                                                      300., cfg['temperature_step_K']).tolist() + [300.]))):
            phase(runner, directory / f'04_reheat_{index:03d}', cfg, t, cfg['reheat_step_ps'])
        phase(runner, directory / '05_equilibrate_300K', cfg, 300., cfg['equilibration_300_ps'])
        rows = phase(runner, directory / '06_production_300K', cfg, 300., cfg['production_ps'], sample=True)
        relaxed = runner.structure_data()
        write_json(directory / 'relaxed_300K.structure.json', relaxed)
        mass = sum(host.system.getParticleMass(i).value_in_unit(unit.dalton)
                   for i in range(host.system.getNumParticles()))
        result = dict(smoke=cfg['smoke'], thermal_history=f'{cfg["high_K"]:g} -> {cfg["low_K"]:g} -> 300 K; fixed topology',
            thermal_history_density=history_densities(directory),
            npt=analyze_npt(rows, mass, 300., cfg['pressure_atm'], cfg['min_effective_samples']))
        result['npt']['volume_nm3_mean'] = float(np.mean([r['volume_nm3'] for r in rows]))
        result['cooling'] = cooling_fit([r['target_temperature_K'] for r in cooling],
                                        [r['density_g_cm3'] for r in cooling], cfg,
                                        result['thermal_history_density'])
        snapshot = runner.snapshot()
        result['structure'] = structural_diagnostics(host, host, snapshot['positions'], snapshot['box'],
            0, seed=replica_seed(project, cfg, replica), free_volume_samples=100 if cfg['smoke'] else 10000)
        result['topology'] = topology(host.graphs, snapshot['positions'], snapshot['box'])
        from workflows.network_properties.network_graph import network_structure
        result['network'] = network_structure(host.graphs, snapshot['positions'], snapshot['box'])
        result['npt']['sampling_checks_passed'] &= not cfg['smoke']
        reason = vaporization_reason(cfg, result['thermal_history_density'])
        if reason:
            # Voids from boiling at high_K remain in the glass: neither the cooling
            # curve nor the 300 K fluctuations describe the bulk material.
            result['npt']['sampling_checks_passed'] = False
            result['npt']['vaporization_warning'] = reason
            print(f'WARNING: {reason}', flush=True)
        result['nvt'] = nvt_properties(runner, directory / '07_nvt_300K', cfg, mass)
        write_json(directory / 'dry_result.json', result)
        try:
            from workflows.network_properties.plots import plot_tg
            print(f'Tg figure: {plot_tg(directory, result)}', flush=True)
        except Exception as exc:  # a figure must never fail a finished calculation
            print(f'WARNING: Tg figure not written: {exc}', flush=True)
    finally:
        runner.close()
    return result


def solvent_properties(project, cfg, replica, specification):
    from polypaves import Solvent, pack_network
    from polypaves.sorption import load_network, load_template
    from polypaves.sorption_workflow import FixedLoadingConfig, run_fixed_loading, run_alchemical
    base = project / f'replica_{replica:02d}'
    directory = base / 'solvents' / specification['name']
    directory.mkdir(parents=True, exist_ok=True)
    if (directory / 'result.json').exists():
        return json.loads((directory / 'result.json').read_text())
    seed = replica_seed(project, cfg, replica)
    name = specification['name']
    host = load_network(base / 'relaxed_300K.structure.json')
    packed_path = directory / 'packed'
    if not packed_path.exists():
        set_stage(base, f'solvent {name}: packing', molecules=cfg['solvent_molecules'])
        packed = pack_network(base / 'relaxed_300K.structure.json',
            Solvent(specification['smiles'], name=specification['name']), cfg['solvent_molecules'], seed=seed)
        with tempfile.TemporaryDirectory(prefix='packing_', dir=directory) as temporary:
            complete = Path(temporary) / 'packed'
            packed.save(complete)
            complete.rename(packed_path)
    solvent = load_template(packed_path, 'solvent')
    fe_dir = directory / 'infinite_dilution_free_energy'
    set_stage(base, f'solvent {name}: insertion free energy', windows=cfg['alchemical_windows'])
    # N=0: insert a noninteracting test molecule in the dry relaxed host, then couple.
    # The endpoint difference contains the matched isolated flexible-molecule reference.
    chemical = run_alchemical(host, solvent, host.positions_nm, 0, host.box_nm, fe_dir,
        direction='couple', temperature=300., pressure_bar=cfg['pressure_atm'] * 1.01325,
        platform=cfg['platform'], seed=seed, timestep_fs=cfg['timestep_fs'],
        windows=cfg['alchemical_windows'], equilibration_ps=cfg['alchemical_equilibration_ps'],
        sampling_ps=cfg['alchemical_sampling_ps'], sample_interval_ps=cfg['sample_interval_ps'],
        min_effective_samples=cfg['min_effective_samples'],
        standard_density_nm3=cfg['standard_density_nm3'],
        volume_interval=cfg['volume_interval_steps'], volume_fraction=cfg['volume_fraction'],
        resume=fe_dir.exists())
    md_cfg = FixedLoadingConfig(platform=cfg['platform'], temperature_K=300.,
        pressure_atm=cfg['pressure_atm'], timestep_fs=cfg['timestep_fs'],
        nvt_ps=cfg['solvent_nvt_ps'], npt_equilibration_ps=cfg['solvent_equilibration_ps'],
        production_ps=cfg['solvent_production_ps'], sample_interval_ps=cfg['sample_interval_ps'],
        frame_interval_ps=max(cfg['sample_interval_ps'], cfg['solvent_production_ps'] / 100),
        minimization_iterations=cfg['minimization_iterations'],
        min_effective_samples=cfg['min_effective_samples'],
        diffusion_min_duration_ps=cfg['diffusion_min_duration_ps'],
        volume_interval_steps=cfg['volume_interval_steps'], volume_fraction=cfg['volume_fraction'])
    md_dir = directory / 'diffusion'
    set_stage(base, f'solvent {name}: diffusion MD', production_ps=cfg['solvent_production_ps'])
    initial = json.loads((packed_path / 'packed.structure.json').read_text())
    wet = run_fixed_loading(host, solvent, 0., md_dir, config=md_cfg, seed=seed,
        count_override=cfg['solvent_molecules'], initial_state=initial,
        resume=md_dir.exists())
    result = dict(smoke=cfg['smoke'], solvent=specification,
        infinite_dilution_free_energy=chemical, fixed_loading=wet,
        definitions=dict(free_energy='Gas/isolated molecule -> dry polymer, infinite dilution.',
            standard_density_nm3=cfg['standard_density_nm3'],
            diffusion='Solvent self diffusion at the specified finite loading; host-relative COM MSD.',
            equilibrium_uptake_verified=False))
    write_json(directory / 'result.json', result)
    return result


def report(project):
    project = Path(project).resolve()
    if not (project / 'project.json').is_file():
        raise FileNotFoundError(project / 'project.json')
    with file_lock(project / '.report.lock'):
        result = _report(project)
    cfg = json.loads((project / 'project.json').read_text())
    if mode(cfg) == 'crosslink':
        # The monomer report holds the crosslink series of all densities.
        from workflows.network_properties.construction import parent_project
        parent = parent_project(project, cfg)
        with file_lock(parent / '.report.lock'):
            _report(parent)
    from workflows.network_properties.directory_summary import MARKER, collect_projects
    root = collection_root(project, cfg)
    if (root / MARKER).is_file():
        collect_projects(root)
    return result


NOTES = dict(
    monomer='Each replica packs and relaxes its own monomer liquid at T with an independent seed; '
            'crosslink replica N of every density starts from monomer replica N.',
    crosslink='Each replica creates its start sites in the relaxed T structure of the same-numbered monomer '
              'replica and crosslinks with its own seed.',
    saved_network='Replicas share the supplied topology; only thermal trajectories/seeds differ.')


def _report(project):
    from workflows.network_properties.replica_summary import collect, write_tables, markdown
    from workflows.network_properties.analysis import aggregate as aggregate_values
    cfg = json.loads((project / 'project.json').read_text())
    check_run_signature(project, cfg)
    refresh_cooling_fits(project, cfg)
    replica_ids, properties = collect(project, cfg)
    write_tables(project, replica_ids, properties)
    dry = [json.loads(path.read_text()) for path in sorted(project.glob('replica_*/dry_result.json'))]
    kind = mode(cfg)
    result = dict(smoke=cfg['smoke'], input_mode=kind, completed_replicas=len(dry), requested_replicas=cfg['replicates'],
        replica_statistics=aggregate_values([dict(npt=r['npt']['estimates'], nvt=r.get('nvt', {}),
            cooling={k:v for k,v in r['cooling'].items() if k != 'bins'}) for r in dry]),
        dry_sampling_checks_passed=bool(len(dry) == cfg['replicates'] and not cfg['smoke']
            and all(r['npt']['sampling_checks_passed'] for r in dry)),
        tg_fit_checks_passed=bool(len(dry) == cfg['replicates'] and not cfg['smoke']
            and all(r['cooling']['fit_checks_passed'] for r in dry)),
        note=NOTES[kind],
        scope=dict(thermodynamics='Density, Cp/Cv, CTE, compressibility, Kt/Ks, bulk sound speed.',
            structure='Geometric diagnostics from the existing wheel.',
            optical='Refractive index/Abbe number require electronic polarizabilities; not inferred from fixed charges.',
            dielectric='Fixed-charge NVT dielectric tensor from continuous atomic charge displacements; no electronic polarization.',
            cohesive_energy='A percolating network has no unique isolated-molecule reference; not reported.',
            mechanical='Optional: network_properties.py mechanical <project> after the 300 K properties.'))
    if kind == 'monomer':
        from workflows.network_properties import monomer
        values = [json.loads(p.read_text()) for p in sorted(project.glob(f'replica_*/{monomer.label(cfg)}/result.json'))]
        result['at_temperature'] = dict(temperature_K=cfg['temperature_K'], completed_replicas=len(values),
            replica_statistics=aggregate_values([dict(npt=v['npt']['estimates'], nvt=v['nvt'],
                reactive_ring_density=v['reactive_ring_density'],
                reactive_rings={k: x for k, x in (v['reactive_rings'] or {}).items() if k != 'definitions'})
                for v in values]),
            sampling_checks_passed=bool(len(values) == cfg['replicates'] and not cfg['smoke']
                and all(v['npt']['sampling_checks_passed'] for v in values)))
    from workflows.network_properties.mechanical import result_current
    paths = [p for p in sorted(project.glob('replica_*/mechanical/result.json')) if result_current(p.parent.parent)]
    if paths:
        values = [json.loads(p.read_text()) for p in paths]
        result['mechanical'] = dict(completed_replicas=len(values),
            replica_statistics=aggregate_values([v.get('isotropic', {}) for v in values]),
            sampling_checks_passed=bool(len(values) == cfg['replicates'] and not cfg['smoke']
                and all(v.get('sampling_checks_passed', False) for v in values)))
    result['solvents'] = {}
    for spec in cfg['solvents']:
        paths = sorted(project.glob(f'replica_*/solvents/{spec["name"]}/result.json'))
        values = [json.loads(p.read_text()) for p in paths]
        result['solvents'][spec['name']] = dict(completed_replicas=len(paths),
            results=[str(p.relative_to(project)) for p in paths],
            replica_statistics=aggregate_values([dict(
                coupling_free_energy_kJ_mol=v['infinite_dilution_free_energy']['coupling_free_energy_kJ_mol'],
                chemical_potential_relative_kJ_mol=v['infinite_dilution_free_energy']['chemical_potential_relative_kJ_mol'],
                diffusion_coefficient_m2_s=v['fixed_loading']['diffusion']['diffusion_coefficient_m2_s']) for v in values]),
            free_energy_convergence_verified=bool(len(values) == cfg['replicates'] and not cfg['smoke'] and
                all(v['infinite_dilution_free_energy']['convergence_verified'] for v in values)),
            diffusion_converged=bool(len(values) == cfg['replicates'] and not cfg['smoke'] and
                all(v['fixed_loading']['diffusion']['converged'] for v in values)))
    result.update(replica_ids=replica_ids, property_summary=properties,
        replica_seeds={f'replica_{i:02d}': seed for i, seed in replica_seeds(project, cfg).items()} if 'seed' in cfg else {},
        property_tables=dict(per_replica='replica_properties.csv', aggregate='property_summary.csv'))
    if kind == 'crosslink':
        result['network_builds'] = {replica: json.loads((project / replica / 'network_origin.json').read_text())
            for replica in replica_ids if (project / replica / 'network_origin.json').is_file()}
    elif kind == 'monomer':
        result['monomer_builds'] = {replica: json.loads((project / replica / 'build/done.json').read_text())
            for replica in replica_ids if (project / replica / 'build/done.json').is_file()}
    titles = dict(monomer=f'Monomer properties ({cfg.get("temperature_K", 0):g} K build)',
                  crosslink='Crosslinked network properties', saved_network='Network property calculation')
    text = [f'# {titles[kind]}', '',
        '**SMOKE: numbers are not physical property predictions.**' if cfg['smoke'] else
        'Inspect sampling and fit diagnostics before interpreting estimates.', '',
        f'Completed thermal replicas: {len(dry)}/{cfg["replicates"]}', '',
        f'NPT sampling checks passed: {result["dry_sampling_checks_passed"]}; Tg fit checks passed: {result["tg_fit_checks_passed"]}.', '']
    if kind == 'crosslink':
        from workflows.network_properties.construction import parent_project
        options = cfg['crosslinking']
        starts = (f'{options["activation_count"]} start sites' if options['activation_count'] is not None
                  else f'start-site fraction {options["activation_fraction"]:g}')
        text += [f'Monomer: [{parent_project(project, cfg).name}](../../report.md); {starts}; '
                 f'crosslinked at {cfg["temperature_K"]:g} K from the relaxed monomer replica of the same number.', '']
    text += markdown(replica_ids, properties)
    if result.get('network_builds'):
        text += ['', '**Network construction per replica**', '',
            '| Replica | Seed | Atoms | Start sites | Conversion | Target | Target reached |',
            '|---|---:|---:|---:|---:|---:|---|']
        for replica, build in result['network_builds'].items():
            text.append(f'| {replica} | {build["seed"]} | {build["atoms"]} | {build.get("activation_count")} | '
                        f'{build["conversion"]:.6g} | {build["target_conversion"]} | {build["target_reached"]} |')
    if result.get('monomer_builds'):
        text += ['', '**Monomer build per replica**', '', '| Replica | Seed | Molecules | Atoms | Reactive rings |',
                 '|---|---:|---:|---:|---:|']
        for replica, build in result['monomer_builds'].items():
            text.append(f'| {replica} | {build["seed"]} | {build["molecules"]} | {build["atoms"]} | {build["reactive_sites"]} |')
    if kind == 'monomer' and children(project):
        from workflows.network_properties.series import write_series
        series = write_series(project)
        result['crosslink_series'] = series['csv']
        text += ['', '**Crosslink series**', ''] + series['markdown']
    if result['replica_seeds']:
        text += ['', 'Seeds: ' + ', '.join(f'{r}={seed}' for r,seed in result['replica_seeds'].items()) + '.']
    from workflows.network_properties.plots import render
    figures = render(project)
    for title, path in figures:
        text += ['', f'![{title}]({path})']
    text += ['', result['note'], '', 'Full method scope and limitations: summary.json.']
    write_json(project / 'summary.json', result)
    (project / 'report.md').write_text('\n'.join(text) + '\n')
    print(f'Report: {project / "report.md"}', flush=True)
    return result


def initialize_project(project):
    project = Path(project).resolve()
    if not (project / 'project.json').is_file():
        raise FileNotFoundError(project / 'project.json')
    with file_lock(project / '.project.lock'):
        return _initialize_project(project)


def _initialize_project(project):
    cfg = json.loads((project / 'project.json').read_text())
    validate(cfg)
    if input_digest(project, cfg) != cfg['input_sha256']:
        raise ValueError('input structure was changed')
    check_run_signature(project, cfg)
    write_json(project / 'run_signature.json', dict(sha256=configuration_fingerprint(cfg)))
    if not (project / 'runtime.json').exists():
        import openmm
        import polypaves
        write_json(project / 'runtime.json', dict(python=sys.version, openmm=openmm.__version__,
            polypaves=polypaves.__version__, numpy=np.__version__, platform=cfg['platform']))
    return cfg


def _replicas(cfg, selected):
    if selected is not None and not 1 <= selected <= cfg['replicates']:
        raise ValueError(f'--replica must lie between 1 and {cfg["replicates"]}')
    return [selected] if selected is not None else list(range(1, cfg['replicates'] + 1))


def run(args):
    project = args.project.resolve()
    cfg = initialize_project(project)
    kind = mode(cfg)
    selected = getattr(args, 'replica', None)
    for replica in _replicas(cfg, selected):
        base = project / f'replica_{replica:02d}'
        seed = replica_seed(project, cfg, replica)
        # The lock comes first: a refused duplicate must not overwrite progress/error records.
        with file_lock(base / '.run.lock', blocking=False) as lock, replica_progress(base):
            write_json(base / 'execution.json', dict(replica=replica, seed=seed,
                platform=cfg['platform'], cuda_visible_devices=os.environ.get('CUDA_VISIBLE_DEVICES')))
            if kind == 'monomer':
                from workflows.network_properties.construction import ensure_monomer_build
                set_stage(base, f'monomer build: pack, compress and relax at {cfg["temperature_K"]:g} K', seed=seed)
                ensure_monomer_build(project, cfg, replica, lock)
            elif kind == 'crosslink':
                from workflows.network_properties.construction import ensure_crosslink
                set_stage(base, 'crosslinking: start sites in the relaxed monomer, react', seed=seed)
                ensure_crosslink(project, cfg, replica, lock)
            set_stage(base, 'dry properties: thermal history')
            dry_properties(project, cfg, replica)
            if kind == 'monomer':
                from workflows.network_properties.monomer import temperature_properties
                set_stage(base, f'monomer properties at {cfg["temperature_K"]:g} K')
                temperature_properties(project, cfg, replica)
            for spec in cfg['solvents']:
                solvent_properties(project, cfg, replica, spec)
        if not getattr(args, 'defer_report', False):
            report(project)
    if getattr(args, 'all', False):
        if kind != 'monomer':
            raise ValueError('--all runs a monomer project and its crosslink projects')
        for child in children(project):
            run(argparse.Namespace(project=child, replica=selected, defer_report=getattr(args, 'defer_report', False)))


def mechanical(args):
    """Optional elastic constants of finished crosslink/network replicas (relaxed 300 K structure)."""
    from workflows.network_properties.mechanical import run_mechanical
    project = args.project.resolve()
    cfg = json.loads((project / 'project.json').read_text())
    validate(cfg)
    check_run_signature(project, cfg)
    if mode(cfg) == 'monomer':
        raise ValueError('mechanical properties are computed for crosslink or network projects, not the monomer liquid')
    replicas = _replicas(cfg, getattr(args, 'replica', None))
    missing = [r for r in replicas if not (project / f'replica_{r:02d}/relaxed_300K.structure.json').is_file()]
    if missing:
        raise RuntimeError(f'replicas {missing} have no relaxed 300 K structure yet; run the project first')
    for replica in replicas:
        base = project / f'replica_{replica:02d}'
        with file_lock(base / '.run.lock', blocking=False), replica_progress(base):
            set_stage(base, 'mechanical: elastic constants')
            run_mechanical(base, cfg, replica_seed(project, cfg, replica))
    if not getattr(args, 'defer_report', False):
        report(project)


def _construction_arguments(parser):
    parser.add_argument('--chemistry', choices=('epoxy', 'oxetane'),
                        help='ring-opening preset: activation/reaction YAMLs and default monomer (default epoxy)')
    parser.add_argument('--activation-yaml', type=Path, help='own activation rules (default: the chemistry preset)')
    parser.add_argument('--reactions-yaml', type=Path, help='own reaction rules (default: the chemistry preset)')


def _thermal_arguments(parser, solvents=True):
    parser.add_argument('--high-temperature', type=float,
                        help='highest temperature of the thermal history in K (default 800)')
    if solvents:
        parser.add_argument('--solvent', action='append', default=[], metavar='NAME=SMILES')
        parser.add_argument('--solvent-molecules', type=int, help='molecules per solvent (default 10)')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command', required=True)

    mono = sub.add_parser('monomer', help='create a monomer project (build at T, Tg, 300 K and T properties)')
    mono.add_argument('project', type=Path)
    mono.add_argument('--smiles', help='monomer SMILES; default TMP triglycidyl ether (epoxy) or OXT-221 (oxetane)')
    _construction_arguments(mono)
    mono.add_argument('--temperature', type=float, default=monomer_defaults()['temperature_K'],
                      help='build, T properties and later crosslinking temperature in K (default 363.15 = 90 C)')
    mono.add_argument('--replicates', type=int, help='independent monomer packings (default 4; smoke 1)')
    mono.add_argument('--seed', type=int, default=31415)
    mono.add_argument('--platform', choices=('CUDA', 'CPU'), default='CUDA')
    mono.add_argument('--total-atoms', type=int, help='target atom count (default 20000, rounded above)')
    mono.add_argument('--atom-count-mode', choices=('above', 'nearest', 'floor', 'exact'))
    mono.add_argument('--packing-density', type=float, help='initial packing density in g/cm3 (default 0.5)')
    mono.add_argument('--compression-pressure-mpa', type=float, help='compression NPT pressure in MPa (default 200)')
    mono.add_argument('--compression-steps', type=int, help='compression NPT steps at 1 fs (default 50000; 0 disables)')
    mono.add_argument('--equilibration-steps', type=int, help='1 atm NPT steps at T and 1 fs after compression (default 1000000)')
    _thermal_arguments(mono, solvents=False)
    mono.add_argument('--smoke', action='store_true')
    mono.set_defaults(action=new_monomer)

    cross = sub.add_parser('crosslink', help='create crosslink projects (one per start-site density) of a monomer project')
    cross.add_argument('monomer', type=Path, help='monomer project directory')
    starts = cross.add_mutually_exclusive_group(required=True)
    starts.add_argument('--activation-fraction', type=float, nargs='+', help='fractions of the rings activated (0 <= f < 1)')
    starts.add_argument('--activation-count', type=int, nargs='+', help='numbers of start sites')
    cross.add_argument('--name', help='directory name under crosslinks/ for a single density (default af<f> or n<count>)')
    cross.add_argument('--replicates', type=int, help='replicas (default and maximum: the monomer replicates)')
    cross.add_argument('--seed', type=int, help='base seed (default: the monomer seed, so start sites are nested)')
    cross.add_argument('--platform', choices=('CUDA', 'CPU'), help='default: the monomer platform')
    cross.add_argument('--cycles', type=int, help='maximum crosslinking cycles (default 5000)')
    cross.add_argument('--target-conversion', type=float, help='crosslinking target (default 0.995)')
    cross.add_argument('--max-reactions-per-cycle', type=int)
    cross.add_argument('--equilibration-steps', type=int,
                       help='NPT steps at T and 1 fs after the start sites are added (default 100000 = 100 ps)')
    cross.add_argument('--final-relaxation-steps', type=int,
                       help='after crosslinking: NVT steps at 0.1 fs before the final MD (default 1000)')
    cross.add_argument('--final-md-steps', type=int,
                       help='after crosslinking: NPT steps at 1 fs and 1 atm before saving (default 1000000 = 1 ns)')
    cross.add_argument('--reaction-cutoff', type=float,
                       help='initial reaction search distance in angstrom (default 4, capped by the maximum)')
    cross.add_argument('--max-reaction-cutoff', type=float, help='largest reaction search distance in angstrom (default 6)')
    cross.add_argument('--fixed-reaction-cutoff', action='store_true', help='keep the reaction cutoff fixed')
    cross.add_argument('--reaction-minimization-iterations', type=int,
                       help='minimization iterations after each reaction (default 2000)')
    cross.add_argument('--reaction-minimization-radius', type=float,
                       help='radius in angstrom minimized after reactions (default 12; 0 = whole system)')
    cross.add_argument('--post-reaction-md-steps', type=int,
                       help='NVT steps after the post-reaction minimization, before commit (default 0)')
    cross.add_argument('--relaxation-stages', type=relaxation_stages_argument,
                       help='post-reaction MD stages: JSON text or a JSON/YAML file (see crosslinking/README.md)')
    cross.add_argument('--stop-after-idle-cycles', type=int,
                       help='stop after this many reactionless cycles at the largest cutoff (default 200; 0 disables)')
    cross.add_argument('--template-workers', type=int, help='parallel workers for new reaction templates (default 4)')
    cross.add_argument('--save-trajectory', action='store_true', help='save the OVITO trajectory and cycle structures')
    cross.add_argument('--remove-sol-fraction', action='store_true',
                       help='keep only the network percolating in three dimensions, dropping unreacted monomers '
                            'and unattached oligomers (default OFF)')
    _thermal_arguments(cross)
    cross.set_defaults(action=new_crosslink)

    network = sub.add_parser('network', help='create a project for a saved dry network (Tg, 300 K, solvents)')
    network.add_argument('project', type=Path)
    network.add_argument('--network', type=Path, required=True, help='saved-charge structure JSON')
    network.add_argument('--replicates', type=int)
    network.add_argument('--seed', type=int, default=31415)
    network.add_argument('--platform', choices=('CUDA', 'CPU'), default='CUDA')
    _thermal_arguments(network)
    network.add_argument('--smoke', action='store_true')
    network.set_defaults(action=new_network)

    execute = sub.add_parser('run', help='run replicas locally (one after another)')
    execute.add_argument('project', type=Path)
    execute.add_argument('--replica', type=int, help='run only this replica (1-based)')
    execute.add_argument('--all', action='store_true', help='monomer project: also run every crosslink project')
    execute.add_argument('--defer-report', action='store_true', help=argparse.SUPPRESS)
    execute.set_defaults(action=run)
    elastic = sub.add_parser('mechanical', help='optional elastic constants of finished crosslink/network replicas')
    elastic.add_argument('project', type=Path)
    elastic.add_argument('--replica', type=int, help='only this replica (1-based)')
    elastic.add_argument('--defer-report', action='store_true', help=argparse.SUPPRESS)
    elastic.set_defaults(action=mechanical)
    summary = sub.add_parser('report', help='rewrite report.md and the property tables (no MD)')
    summary.add_argument('project', type=Path)
    summary.set_defaults(action=lambda a: report(a.project))

    from workflows.network_properties.gpu_jobs import reseed, submit, status
    launch = sub.add_parser('submit', help='launch replica jobs on selected local GPUs')
    launch.add_argument('project', type=Path)
    launch.add_argument('--gpus', help='physical GPU indices, e.g. 0,2; default: all visible GPUs')
    launch.add_argument('--dry-run', action='store_true', help='show the GPU/seed plan without launching jobs')
    launch.add_argument('--mechanical', action='store_true',
                        help='queue the optional elastic constants of finished replicas instead')
    launch.add_argument('--monomer-only', action='store_true',
                        help='monomer project: queue only its own replicas, not the crosslink projects')
    launch.set_defaults(action=submit)
    retry = sub.add_parser('reseed', help='give failed replicas a new seed and archive the failed attempt')
    retry.add_argument('project', type=Path)
    retry.add_argument('--replica', type=int, action='append', required=True,
                       help='failed replica to rerun (1-based); repeat for several')
    retry.add_argument('--seed', type=int, help='new seed for one replica (default: automatic, unused)')
    retry.add_argument('--force', action='store_true', help='also reseed a replica that has not failed')
    retry.add_argument('--submit', action='store_true', help='submit the reseeded replicas to the GPUs right away')
    retry.add_argument('--gpus', help='with --submit: physical GPU indices, e.g. 0,2')
    retry.add_argument('--dry-run', action='store_true', help='show the new seeds without changing anything')
    retry.set_defaults(action=reseed, mechanical=False, monomer_only=False)
    check = sub.add_parser('status', help='show replica stages and GPU job states (monomer: with its crosslinks)')
    check.add_argument('project', type=Path)
    check.set_defaults(action=status)
    comparison = sub.add_parser('collect', help='summarize all project directories under a parent directory')
    comparison.add_argument('root', type=Path)
    from workflows.network_properties.directory_summary import collect_projects
    comparison.set_defaults(action=lambda a: collect_projects(a.root))
    args = parser.parse_args(argv)
    try:
        args.action(args)
    except (ValueError, OSError, RuntimeError, ImportError, KeyError) as exc:
        if args.command in ('run', 'mechanical'):
            # Keep the traceback: job logs are the only record of where a detached run stopped.
            import traceback
            traceback.print_exc()
        print(f'network_properties: {exc}', file=sys.stderr)
        return 1
    return 0


def monomer_defaults():
    from workflows.network_properties.monomer import DEFAULTS as values
    return values


if __name__ == '__main__':
    sys.exit(main())
