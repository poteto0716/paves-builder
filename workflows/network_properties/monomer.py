"""Monomer properties at the reaction temperature T, after the thermal history and the 300 K stage.

The 300 K NPT production of the thermal history continues (same trajectory)
to T in temperature_step_K steps of reheat_step_ps each, equilibrates for
t_equilibration_ps and samples t_production_ps of NPT. That production gives
the T thermodynamics and the reactive-ring statistics (number density,
coordination, largest contact cluster, ring MSD and diffusion; see
workflows/crosslinking/pre_reaction.py); its last frame is saved as
relaxed_<T>K.structure.json, the start of every crosslink replica with the
same number. An NVT branch at T then gives Cv and the dielectric response.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

RING_LABELS = {'epoxy': 'Epoxide', 'oxetane': 'Oxetane'}
DEFAULTS = dict(temperature_K=363.15, t_equilibration_ps=1000., t_production_ps=2000.,
                t_nvt_production_ps=2000., ring_contact_cutoffs_A=[4., 6.])
SMOKE = dict(t_equilibration_ps=.002, t_production_ps=.04, t_nvt_production_ps=.02)


def label(cfg):
    return f"{float(cfg['temperature_K']):g}K"


def ring_label(cfg):
    return RING_LABELS.get(cfg.get('chemistry'), 'Reactive ring')


def temperature_dir(project, cfg, replica):
    return Path(project) / f'replica_{replica:02d}' / label(cfg)


def relaxed_path(project, cfg, replica):
    return temperature_dir(project, cfg, replica) / f'relaxed_{label(cfg)}.structure.json'


def replica_complete(project, cfg, replica):
    base = Path(project) / f'replica_{replica:02d}'
    return all(path.is_file() for path in (base / 'build/done.json', base / 'dry_result.json',
                                           temperature_dir(project, cfg, replica) / 'result.json',
                                           relaxed_path(project, cfg, replica)))


def heating_grid(start, end, step):
    """Temperatures from start (exclusive) to end (inclusive) in steps of at most step."""
    count = max(1, int(np.ceil(abs(end - start) / step - 1e-9)))
    return [float(start + (end - start) * (i + 1) / count) for i in range(count)]


class RingObserver:
    """Feeds NPT production frames to the reactive-ring sampler; phase() calls add and finish."""

    def __init__(self, groups, masses, cutoffs, temperature_K, pressure_atm):
        from workflows.crosslinking.pre_reaction import RingSampler
        self.sampler = RingSampler(groups, masses, cutoffs)
        self.temperature_K, self.pressure_atm = temperature_K, pressure_atm

    def add(self, runner):
        state = runner.snapshot()
        self.sampler.add(state['positions'] * 10., state['box'] * 10., state['time'])

    def finish(self, directory, duration_ps):
        self.sampler.summary(Path(directory) / 'reactive_rings', temperature_K=self.temperature_K,
                             pressure_atm=self.pressure_atm, sampling_ps=duration_ps)


def temperature_properties(project, cfg, replica):
    import openmm as mm
    from openmm import unit
    from polypaves.sorption import load_network
    from polypaves.sorption_workflow import FixedLoadingSimulation
    from workflows.crosslinking.pre_reaction import PER_NM3_TO_MMOL_CM3, ring_name, smarts_groups
    from workflows.network_properties.analysis import analyze_npt
    from workflows.network_properties.construction import network_path, target_smarts
    from workflows.network_properties.execution import replica_seed
    from workflows.network_properties.network_properties import nvt_properties, phase, write_json
    project = Path(project)
    base = project / f'replica_{replica:02d}'
    directory = temperature_dir(project, cfg, replica)
    if (directory / 'result.json').exists():
        return json.loads((directory / 'result.json').read_text())
    directory.mkdir(exist_ok=True)
    temperature = float(cfg['temperature_K'])
    dry = json.loads((base / 'dry_result.json').read_text())
    host = load_network(network_path(project, cfg, replica))
    runner = FixedLoadingSimulation(host, host, host.positions_nm, 0, box_nm=host.box_nm,
        temperature=300., pressure_bar=cfg['pressure_atm'] * 1.01325,
        timestep_fs=cfg['timestep_fs'], platform=cfg['platform'], seed=replica_seed(project, cfg, replica),
        volume_interval=cfg['volume_interval_steps'], volume_fraction=cfg['volume_fraction'])
    try:
        # Continue the 300 K NPT production of the thermal history (positions, velocities and box).
        runner.context.setState(mm.XmlSerializer.deserialize((base / '06_production_300K/state.xml').read_text()))
        masses = [host.system.getParticleMass(i).value_in_unit(unit.dalton)
                  for i in range(host.system.getNumParticles())]
        mass = float(sum(masses))
        grid = heating_grid(300., temperature, cfg['temperature_step_K'])
        for index, t in enumerate(grid[:-1]):
            phase(runner, directory / f'01_heat_{index:03d}', cfg, t, cfg['reheat_step_ps'], base=base)
        phase(runner, directory / '02_equilibrate', cfg, temperature, cfg['t_equilibration_ps'], base=base)
        groups = smarts_groups(host.graphs, target_smarts(project))
        observer = (RingObserver(groups, masses, cfg['ring_contact_cutoffs_A'], temperature, cfg['pressure_atm'])
                    if groups else None)
        rows = phase(runner, directory / '03_production', cfg, temperature, cfg['t_production_ps'],
                     sample=True, base=base, observe=observer)
        # The crosslink replicas start here: the last NPT frame at T.
        write_json(relaxed_path(project, cfg, replica), runner.structure_data())
        npt = analyze_npt(rows, mass, temperature, cfg['pressure_atm'], cfg['min_effective_samples'])
        npt['volume_nm3_mean'] = float(np.mean([r['volume_nm3'] for r in rows]))
        npt['sampling_checks_passed'] &= not cfg['smoke']
        rings = (json.loads((directory / '03_production/reactive_rings/reactive_rings.json').read_text())
                 if groups else None)
        count = len(groups)
        volumes = {'300K': dry['npt']['volume_nm3_mean'], label(cfg): npt['volume_nm3_mean']}
        density = dict(ring=ring_name(len(groups[0]['heavy'])) if groups else None, rings=count,
            target_smarts=target_smarts(project),
            **{f'number_density_{key}_nm3': count / volume for key, volume in volumes.items()},
            **{f'number_density_{key}_mmol_cm3': count / volume * PER_NM3_TO_MMOL_CM3
               for key, volume in volumes.items()},
            definition='unreacted reactive rings (activation target_smarts) / mean NPT cell volume '
                       'at that temperature; 1 nm^-3 = 1/0.6022 mmol/cm^3')
        nvt = nvt_properties(runner, directory / '04_nvt', cfg, mass, temperature=temperature,
                             equilibration_ps=cfg['t_equilibration_ps'],
                             production_ps=cfg['t_nvt_production_ps'], base=base)
        result = dict(smoke=cfg['smoke'], temperature_K=temperature,
                      history=f'thermal history -> 300 K production -> {temperature:g} K; fixed topology',
                      npt=npt, nvt=nvt, reactive_rings=rings, reactive_ring_density=density,
                      relaxed_structure=str(relaxed_path(project, cfg, replica).relative_to(project)))
        write_json(directory / 'result.json', result)
    finally:
        runner.close()
    return result
