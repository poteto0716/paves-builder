"""Per-replica epoxy construction -> thermal history -> dry and solvent properties."""
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
from workflows.network_properties.execution import file_lock, replica_progress, set_stage


DEFAULTS = dict(replicates=4, seed=31415, platform='CUDA', timestep_fs=.5,
    pressure_atm=1., high_K=800., low_K=100., temperature_step_K=20.,
    heat_step_ps=100., high_equilibration_ps=1000., cooling_rate_K_ns=40.,
    reheat_step_ps=100., equilibration_300_ps=1000., production_ps=10000.,
    sample_interval_ps=1., min_effective_samples=100, minimization_iterations=3000,
    volume_interval_steps=100, volume_fraction=.002, solvents=[],
    solvent_molecules=10, solvent_nvt_ps=100., solvent_equilibration_ps=1000.,
    solvent_production_ps=10000., diffusion_min_duration_ps=1000.,
    alchemical_windows=22, alchemical_equilibration_ps=100.,
    alchemical_sampling_ps=500., standard_density_nm3=.602214076,
    mechanical=True, nvt_production_ps=10000.,
    elastic_strain=.01, elastic_equilibration_ps=200., elastic_production_ps=1000.,
    smoke=False)


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


def input_digest(project, cfg):
    if cfg.get('input_mode') == 'epoxy':
        from workflows.network_properties.epoxy_pipeline import definition_digest
        return definition_digest(project)
    return digest(Path(project) / cfg['network'])


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
    if cfg.get('input_mode') == 'epoxy':
        from workflows.network_properties.epoxy_pipeline import validate as validate_epoxy
        validate_epoxy(cfg)


def relaxation_stages_argument(value):
    # Imported here so that status/submit do not load polypaves for argument parsing.
    from workflows.crosslinking.crosslinking import load_relaxation_stages
    return load_relaxation_stages(value)


def new(args):
    if args.project.exists():
        raise FileExistsError(f'Project exists: {args.project}')
    cfg = json.loads(json.dumps(DEFAULTS))
    cfg.update(platform=args.platform,
               replicates=args.replicates if args.replicates is not None else (1 if args.smoke else DEFAULTS['replicates']),
               seed=args.seed,
               solvent_molecules=args.solvent_molecules, mechanical=not args.skip_mechanical)
    for entry in args.solvent:
        name, sep, smiles = entry.partition('=')
        if not sep:
            raise ValueError('--solvent expects NAME=SMILES')
        cfg['solvents'].append(dict(name=name, smiles=smiles))
    if args.smoke:
        cfg.update(smoke=True, timestep_fs=.25,
            heat_step_ps=.002, high_equilibration_ps=.002,
            cooling_rate_K_ns=10_000_000., reheat_step_ps=.002,
            equilibration_300_ps=.002, production_ps=.04, sample_interval_ps=.001,
            minimization_iterations=100, volume_interval_steps=1,
            solvent_nvt_ps=.002, solvent_equilibration_ps=.002, solvent_production_ps=.02,
            alchemical_windows=4, alchemical_equilibration_ps=.001, alchemical_sampling_ps=.004,
            nvt_production_ps=.02)
    definitions = None
    if args.network is not None:
        construction_flags = ('total_atoms', 'atom_count_mode', 'cycles', 'target_conversion',
                              'activation_count', 'activation_fraction', 'max_reactions_per_cycle',
                              'activation_yaml', 'reactions_yaml', 'compression_pressure_mpa',
                              'compression_steps', 'reaction_cutoff', 'max_reaction_cutoff',
                              'reaction_minimization_iterations', 'post_reaction_md_steps',
                              'relaxation_stages', 'stop_after_idle_cycles', 'template_workers',
                              'reaction_minimization_radius', 'packing_density',
                              'final_relaxation_steps', 'final_md_steps')
        if (any(getattr(args, name, None) is not None for name in construction_flags)
                or args.save_trajectory or args.fixed_reaction_cutoff):
            raise ValueError('--network cannot be combined with epoxy construction options')
        source = args.network.resolve()
        data = json.loads(source.read_text())
        if not all(k in data for k in ('graphs', 'positions', 'box_vectors')):
            raise ValueError('network must be a PAVES saved-charge structure JSON')
        if data.get('solvent_molecules', 0):
            raise ValueError('input must be a dry network')
        if not all('partial_charge' in a for g in data['graphs'] for a in g['atoms']):
            raise ValueError('all network atoms must have saved partial charges')
        cfg.update(network='input.structure.json', source=str(source), input_sha256=digest(source),
                   input_mode='saved_network', schema_version=1)
    else:
        from workflows.network_properties.epoxy_pipeline import configure
        construction, definitions = configure(args)
        cfg.update(construction)
    validate(cfg)
    args.project.mkdir(parents=True)
    if definitions is None:
        (args.project / cfg['network']).write_bytes(source.read_bytes())
    else:
        (args.project / 'crosslinking').mkdir()
        for name, content in definitions.items():
            (args.project / 'crosslinking' / name).write_bytes(content)
    write_json(args.project / 'project.json', cfg)
    from workflows.network_properties.directory_summary import enable_collection
    enable_collection(args.project.parent)
    print(f'Created {args.project / "project.json"}; edit durations there before run.', flush=True)


def steps(ps, cfg):
    return max(1, int(round(ps * 1000 / cfg['timestep_fs'])))


def set_temperature(runner, temperature):
    from openmm import unit
    runner.temperature = float(temperature)
    runner.integrator.setTemperature(temperature * unit.kelvin)


def phase(runner, directory, cfg, temperature, duration_ps, *, sample=False, minimize=False):
    """Completed phases resume from portable state; interrupted phases rerun."""
    import openmm as mm
    from openmm import unit
    directory.mkdir(parents=True, exist_ok=True)
    marker = directory / 'done.json'
    set_temperature(runner, temperature)
    if marker.exists():
        runner.context.setState(mm.XmlSerializer.deserialize((directory / 'state.xml').read_text()))
        return json.loads((directory / 'samples.json').read_text())
    print(f'{directory.name}: {temperature:g} K, {duration_ps:g} ps', flush=True)
    set_stage(directory.parent, f'thermal history: {directory.name}',
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
    state = runner.context.getState(getPositions=True, getVelocities=True)
    (directory / 'state.xml').write_text(mm.XmlSerializer.serialize(state))
    write_json(directory / 'final.structure.json', runner.structure_data())
    write_json(directory / 'samples.json', rows)
    write_csv(directory / 'samples.csv', rows)
    write_json(marker, dict(temperature_K=temperature, duration_ps=duration_ps,
        actual_duration_ps=elapsed * cfg['timestep_fs'] / 1000,
        volume_attempts=runner.volume_attempts, volume_accepts=runner.volume_accepts))
    return rows


def temperature_grid(high, low, increment):
    # Include exactly 300 K even when the chosen step does not divide the range.
    values = list(np.arange(high, low, -increment)) + [low, 300.]
    return sorted(set(float(x) for x in values if low <= x <= high), reverse=True)


def nvt_properties(runner, directory, cfg, mass_Da):
    from openmm import unit
    directory.mkdir(exist_ok=True)
    if (directory / 'result.json').exists():
        return json.loads((directory / 'result.json').read_text())
    print('300 K NVT: direct Cv and dielectric response', flush=True)
    set_stage(directory.parent, '300 K NVT: direct Cv and dielectric response')
    # Restart this branch from the completed NPT phase after an interruption.
    runner.md(steps(cfg['equilibration_300_ps'], cfg), ensemble='NVT')
    charges = np.array([a['partial_charge'] for g in runner.host.graphs for a in g['atoms']])
    if abs(charges.sum()) > 1e-4:
        raise ValueError('dielectric estimator requires a neutral network')
    initial = runner.snapshot()
    previous, box = initial['positions'], initial['box']
    inverse = np.linalg.inv(box)
    moment, rows = np.zeros(3), []
    total, elapsed = steps(cfg.get('nvt_production_ps', cfg['production_ps']), cfg), 0
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
    result = analyze_nvt(rows, mass_Da, min_effective_samples=cfg['min_effective_samples'])
    write_json(directory / 'samples.json', rows)
    write_csv(directory / 'samples.csv', rows)
    write_json(directory / 'result.json', result)
    return result


def dry_properties(project, cfg, replica):
    from openmm import unit
    from polypaves.sorption import load_network
    from polypaves.sorption_workflow import FixedLoadingSimulation, structural_diagnostics
    directory = project / f'replica_{replica:02d}'
    directory.mkdir(exist_ok=True)
    if (directory / 'dry_result.json').exists():
        return json.loads((directory / 'dry_result.json').read_text())
    from workflows.network_properties.epoxy_pipeline import network_path
    host = load_network(network_path(project, cfg, replica))
    # Count zero means no solvent is appended; reusing host avoids typing a dummy solvent.
    runner = FixedLoadingSimulation(host, host, host.positions_nm, 0, box_nm=host.box_nm,
        temperature=300., pressure_bar=cfg['pressure_atm'] * 1.01325,
        timestep_fs=cfg['timestep_fs'], platform=cfg['platform'], seed=cfg['seed'] + replica,
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
            npt=analyze_npt(rows, mass, 300., cfg['pressure_atm'], cfg['min_effective_samples']))
        try:
            result['cooling'] = fit_tg([r['target_temperature_K'] for r in cooling],
                                      [r['density_g_cm3'] for r in cooling], bin_width_K=10.)
            result['cooling']['alpha_volume_glass_1_K'] = 3 * result['cooling']['alpha_linear_glass_1_K']
            result['cooling']['alpha_volume_rubber_1_K'] = 3 * result['cooling']['alpha_linear_rubber_1_K']
            fit = result['cooling']
            t = np.array([b['temperature_K'] for b in fit['bins']])
            rho = np.array([b['density_g_cm3'] for b in fit['bins']])
            linear_sse = np.sum((rho - np.polyval(np.polyfit(t,rho,1),t))**2)
            curved_sse = len(t) * fit['rmse_g_cm3']**2
            fit['AIC_improvement_over_linear'] = float(len(t)*np.log(max(linear_sse,1e-30)/max(curved_sse,1e-30))-6)
            span = cfg['high_K'] - cfg['low_K']
            fit['fit_checks_passed'] = bool(not cfg['smoke'] and fit['r_squared'] is not None
                and fit['r_squared'] > .95 and fit['transition_width_K'] < span / 4
                and fit['AIC_improvement_over_linear'] > 6 and fit['alpha_linear_glass_1_K'] > 0
                and cfg['low_K'] + .1 * span < fit['tg_K'] < cfg['high_K'] - .1 * span)
        except (ValueError, RuntimeError) as exc:
            result['cooling'] = dict(fit_checks_passed=False, reason=str(exc))
        snapshot = runner.snapshot()
        result['structure'] = structural_diagnostics(host, host, snapshot['positions'], snapshot['box'],
            0, seed=cfg['seed'] + replica, free_volume_samples=100 if cfg['smoke'] else 10000)
        result['topology'] = topology(host.graphs, snapshot['positions'], snapshot['box'])
        from workflows.network_properties.network_graph import network_structure
        result['network'] = network_structure(host.graphs, snapshot['positions'], snapshot['box'])
        result['npt']['sampling_checks_passed'] &= not cfg['smoke']
        result['nvt'] = nvt_properties(runner, directory / '07_nvt_300K', cfg, mass)
        write_json(directory / 'dry_result.json', result)
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
    seed = cfg['seed'] + replica
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
    project = Path(project)
    if not (project / 'project.json').is_file():
        raise FileNotFoundError(project / 'project.json')
    with file_lock(project / '.report.lock'):
        result = _report(project)
    from workflows.network_properties.directory_summary import MARKER, collect_projects
    if (project.parent / MARKER).is_file():
        collect_projects(project.parent)
    return result


def _report(project):
    from workflows.network_properties.replica_summary import collect, write_tables, markdown
    cfg = json.loads((project / 'project.json').read_text())
    replica_ids, properties = collect(project, cfg)
    write_tables(project, replica_ids, properties)
    dry = [json.loads(path.read_text()) for path in sorted(project.glob('replica_*/dry_result.json'))]
    result = dict(smoke=cfg['smoke'], completed_replicas=len(dry), requested_replicas=cfg['replicates'],
        replica_statistics=aggregate([dict(npt=r['npt']['estimates'], nvt=r.get('nvt', {}), cooling={k:v for k,v in r['cooling'].items()
            if k != 'bins'}) for r in dry]),
        dry_sampling_checks_passed=bool(len(dry) == cfg['replicates'] and not cfg['smoke']
            and all(r['npt']['sampling_checks_passed'] for r in dry)),
        tg_fit_checks_passed=bool(len(dry) == cfg['replicates'] and not cfg['smoke']
            and all(r['cooling']['fit_checks_passed'] for r in dry)),
        independent_network_topologies=cfg.get('input_mode') == 'epoxy',
        note=('Each replica packs and crosslinks its own epoxy network using an independent seed.'
              if cfg.get('input_mode') == 'epoxy' else
              'Replicas share the supplied topology; only thermal trajectories/seeds differ.'),
        scope=dict(thermodynamics='Density, Cp/Cv, CTE, compressibility, Kt/Ks, bulk sound speed.',
            structure='Geometric diagnostics from the existing wheel.',
            optical='Refractive index/Abbe number require electronic polarizabilities; not inferred from fixed charges.',
            dielectric='Fixed-charge NVT dielectric tensor from continuous atomic charge displacements; no electronic polarization.',
            cohesive_energy='A percolating network has no unique isolated-molecule reference; not reported.'))
    result['scope']['mechanical'] = 'Requested; see replica outputs.' if cfg['mechanical'] else 'Disabled in project.json.'
    paths = sorted(project.glob('replica_*/mechanical/result.json'))
    if paths:
        values = [json.loads(p.read_text()) for p in paths]
        result['mechanical'] = dict(completed_replicas=len(values),
            replica_statistics=aggregate([v.get('isotropic', {}) for v in values]),
            sampling_checks_passed=bool(len(values) == cfg['replicates'] and not cfg['smoke']
                and all(v.get('sampling_checks_passed', False) for v in values)))
    result['solvents'] = {}
    for spec in cfg['solvents']:
        paths = sorted(project.glob(f'replica_*/solvents/{spec["name"]}/result.json'))
        values = [json.loads(p.read_text()) for p in paths]
        result['solvents'][spec['name']] = dict(completed_replicas=len(paths),
            results=[str(p.relative_to(project)) for p in paths],
            replica_statistics=aggregate([dict(
                coupling_free_energy_kJ_mol=v['infinite_dilution_free_energy']['coupling_free_energy_kJ_mol'],
                chemical_potential_relative_kJ_mol=v['infinite_dilution_free_energy']['chemical_potential_relative_kJ_mol'],
                diffusion_coefficient_m2_s=v['fixed_loading']['diffusion']['diffusion_coefficient_m2_s']) for v in values]),
            free_energy_convergence_verified=bool(len(values) == cfg['replicates'] and not cfg['smoke'] and
                all(v['infinite_dilution_free_energy']['convergence_verified'] for v in values)),
            diffusion_converged=bool(len(values) == cfg['replicates'] and not cfg['smoke'] and
                all(v['fixed_loading']['diffusion']['converged'] for v in values)))
    result.update(replica_ids=replica_ids, property_summary=properties,
        replica_seeds={replica: cfg['seed'] + i for i, replica in enumerate(replica_ids,1)} if 'seed' in cfg else {},
        property_tables=dict(per_replica='replica_properties.csv', aggregate='property_summary.csv'))
    result['network_builds'] = {replica: json.loads((project / replica / 'network_origin.json').read_text())
        for replica in replica_ids if (project / replica / 'network_origin.json').is_file()}
    write_json(project / 'summary.json', result)
    text = ['# Network property calculation', '',
        '**SMOKE: numbers are not physical property predictions.**' if cfg['smoke'] else
        'Inspect sampling and fit diagnostics before interpreting estimates.', '',
        f'Completed thermal replicas: {len(dry)}/{cfg["replicates"]}', '',
        f'NPT sampling checks passed: {result["dry_sampling_checks_passed"]}; Tg fit checks passed: {result["tg_fit_checks_passed"]}.', '']
    text += markdown(replica_ids, properties)
    if result['network_builds']:
        text += ['', '**Epoxy construction per replica**', '',
            '| Replica | Seed | Atoms | Conversion | Target | Target reached |',
            '|---|---:|---:|---:|---:|---|']
        for replica, build in result['network_builds'].items():
            text.append(f'| {replica} | {build["seed"]} | {build["atoms"]} | {build["conversion"]:.6g} | '
                        f'{build["target_conversion"]} | {build["target_reached"]} |')
    if result['replica_seeds']:
        text += ['', 'Seeds: ' + ', '.join(f'{r}={seed}' for r,seed in result['replica_seeds'].items()) + '.']
    from workflows.network_properties.plots import render
    figures = render(project)
    for title, path in figures:
        text += ['', f'![{title}]({path})']
    text += ['', result['note'], '', 'Full method scope and limitations: summary.json.']
    (project / 'report.md').write_text('\n'.join(text) + '\n')
    print(f'Report: {project / "report.md"}', flush=True)


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
    fingerprint = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()
    signature = project / 'run_signature.json'
    if signature.exists() and json.loads(signature.read_text())['sha256'] != fingerprint:
        raise ValueError('configuration changed after starting; create a new project')
    write_json(signature, dict(sha256=fingerprint))
    if not (project / 'runtime.json').exists():
        import openmm
        import polypaves
        write_json(project / 'runtime.json', dict(python=sys.version, openmm=openmm.__version__,
            polypaves=polypaves.__version__, numpy=np.__version__, platform=cfg['platform']))
    return cfg


def run(args):
    project = args.project.resolve()
    cfg = initialize_project(project)
    selected = getattr(args, 'replica', None)
    if selected is not None and not 1 <= selected <= cfg['replicates']:
        raise ValueError(f'--replica must lie between 1 and {cfg["replicates"]}')
    replicas = [selected] if selected is not None else range(1, cfg['replicates'] + 1)
    for replica in replicas:
        base = project / f'replica_{replica:02d}'
        # The lock comes first: a refused duplicate must not overwrite progress/error records.
        with file_lock(base / '.run.lock', blocking=False) as lock, replica_progress(base):
            write_json(base / 'execution.json', dict(replica=replica, seed=cfg['seed'] + replica,
                platform=cfg['platform'], cuda_visible_devices=os.environ.get('CUDA_VISIBLE_DEVICES')))
            if cfg.get('input_mode') == 'epoxy':
                from workflows.network_properties.epoxy_pipeline import ensure_network
                set_stage(base, 'crosslinking: pack, compress, equilibrate and react', seed=cfg['seed'] + replica)
                ensure_network(project, cfg, replica, lock)
            set_stage(base, 'dry properties: thermal history')
            dry_properties(project, cfg, replica)
            if cfg['mechanical']:
                from workflows.network_properties.mechanical import run_mechanical
                set_stage(base, 'mechanical: elastic constants')
                run_mechanical(base, cfg, cfg['seed'] + replica)
            for spec in cfg['solvents']:
                solvent_properties(project, cfg, replica, spec)
        if not getattr(args, 'defer_report', False):
            report(project)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    create = sub.add_parser('new', help='create a full epoxy pipeline, or supply a saved dry network')
    create.add_argument('project', type=Path)
    inputs = create.add_mutually_exclusive_group()
    inputs.add_argument('--network', type=Path, help='reuse one saved network; omit to build epoxy per replica')
    inputs.add_argument('--epoxy-smiles', help='epoxy SMILES; default TMP triglycidyl ether')
    create.add_argument('--total-atoms', type=int, help='epoxy target atom count (default 20000, rounded above)')
    create.add_argument('--atom-count-mode', choices=('above','nearest','floor','exact'))
    create.add_argument('--cycles', type=int, help='maximum crosslinking cycles (default 5000)')
    create.add_argument('--target-conversion', type=float, help='crosslinking target (default 0.995)')
    starts = create.add_mutually_exclusive_group()
    starts.add_argument('--activation-count', type=int)
    starts.add_argument('--activation-fraction', type=float)
    create.add_argument('--max-reactions-per-cycle', type=int)
    create.add_argument('--final-relaxation-steps', type=int,
                        help='after crosslinking: NVT steps at 0.1 fs before the final MD (default 1000)')
    create.add_argument('--final-md-steps', type=int,
                        help='after crosslinking: NPT steps at 1 fs and 1 atm before saving (default 1000000 = 1 ns)')
    create.add_argument('--packing-density', type=float,
                        help='initial epoxy packing density in g/cm3 before compression (default 0.5)')
    create.add_argument('--compression-pressure-mpa', type=float,
                        help='pre-crosslinking compression NPT pressure in MPa (default 200)')
    create.add_argument('--compression-steps', type=int,
                        help='compression NPT steps at 1 fs (default 50000 = 50 ps; 0 disables)')
    create.add_argument('--reaction-cutoff', type=float,
                        help='initial reaction search distance in angstrom (default 5, capped by the maximum)')
    create.add_argument('--max-reaction-cutoff', type=float,
                        help='largest search distance after automatic growth, angstrom (default 7)')
    create.add_argument('--fixed-reaction-cutoff', action='store_true',
                        help='keep the search distance at --reaction-cutoff')
    create.add_argument('--reaction-minimization-iterations', type=int,
                        help='crosslinking minimization iterations, initial and after each reaction (default 2000)')
    create.add_argument('--reaction-minimization-radius', type=float,
                        help='minimize only atoms within this many angstrom of each reaction (default 12; 0 whole system)')
    create.add_argument('--post-reaction-md-steps', type=int,
                        help='NVT steps at 0.1 fs after minimization, before reactions are committed (default 0)')
    create.add_argument('--relaxation-stages', type=relaxation_stages_argument,
                        help='post-reaction MD per cycle: JSON text or a JSON/YAML file (default NVT 0.1 fs x 10000, NPT 1 fs x 1000)')
    create.add_argument('--stop-after-idle-cycles', type=int,
                        help='end crosslinking after this many zero-reaction cycles at the largest search distance (default 200; 0 never)')
    create.add_argument('--template-workers', type=int,
                        help='processes computing new reaction templates in parallel (default 4; 1 serial)')
    create.add_argument('--save-trajectory', action='store_true', help='save crosslinking trajectory (default OFF)')
    create.add_argument('--remove-sol-fraction', action='store_true',
                        help='like development: keep only the network that percolates through the cell, '
                             'dropping unreacted monomers and unattached oligomers (default OFF)')
    create.add_argument('--activation-yaml', type=Path)
    create.add_argument('--reactions-yaml', type=Path)
    create.add_argument('--platform', choices=('CUDA','CPU','Reference'), default='CUDA')
    create.add_argument('--replicates', type=int, help='number of replicas (default 4; smoke default 1)')
    create.add_argument('--seed', type=int, default=31415)
    create.add_argument('--solvent', action='append', default=[], metavar='NAME=SMILES')
    create.add_argument('--solvent-molecules', type=int, default=10)
    create.add_argument('--skip-mechanical', action='store_true')
    create.add_argument('--smoke', action='store_true')
    create.set_defaults(action=new)
    for name, action in [('run',run), ('report',lambda a: report(a.project))]:
        cmd = sub.add_parser(name)
        cmd.add_argument('project', type=Path)
        if name == 'run':
            cmd.add_argument('--replica', type=int, help='run only this replica (1-based)')
            cmd.add_argument('--defer-report', action='store_true', help=argparse.SUPPRESS)
        cmd.set_defaults(action=action)
    from workflows.network_properties.gpu_jobs import submit, status
    launch = sub.add_parser('submit', help='launch replica jobs on selected local GPUs')
    launch.add_argument('project', type=Path)
    launch.add_argument('--gpus', help='physical GPU indices, e.g. 0,2; default: all visible GPUs')
    launch.add_argument('--dry-run', action='store_true', help='show the GPU/seed plan without launching jobs')
    launch.set_defaults(action=submit)
    check = sub.add_parser('status', help='show local GPU job states and logs')
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
        if args.command == 'run':
            # Keep the traceback: job logs are the only record of where a detached run stopped.
            import traceback
            traceback.print_exc()
        print(f'network_properties: {exc}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
