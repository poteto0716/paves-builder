"""Müller-Plathe reverse NEMD with equal-mass atom velocity exchanges.

Natural masses, unconstrained bonded motion, fixed volume, two periodic heat
paths. No thermostat acts during gradient generation or production.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

from workflows.network_properties.analysis import KB, NA
from workflows.network_properties.mechanical import context_for


def conductivity(profile_K, length_nm, area_nm2, transferred_kJ_mol, duration_ps):
    n = len(profile_K)
    if n < 12 or n % 2:
        raise ValueError('use an even slab count >= 12')
    z = (np.arange(n) + .5) * length_nm / n
    fits = []
    for start, end in ((2, n//2-1), (n//2+2, n-1)):
        x, y = z[start:end] * 1e-9, np.asarray(profile_K)[start:end]
        if len(x) < 3 or not np.isfinite(y).all():
            return dict(conductivity_W_m_K=None, reason='Insufficient populated interior slabs.')
        slope, intercept = np.polyfit(x, y, 1)
        residual = np.sum((y - (slope*x+intercept))**2)
        total = np.sum((y-y.mean())**2)
        fits.append(dict(gradient_K_m=float(slope), r_squared=float(1-residual/total) if total > 0 else 0.))
    flux = transferred_kJ_mol * 1000 / NA / (2 * area_nm2 * 1e-18 * duration_ps * 1e-12)
    gradient = np.mean([abs(f['gradient_K_m']) for f in fits])
    opposite = fits[0]['gradient_K_m'] * fits[1]['gradient_K_m'] < 0
    return dict(conductivity_W_m_K=float(flux/gradient) if gradient > 0 and flux > 0 and opposite else None,
                heat_flux_W_m2=float(flux), fits=fits, opposite_gradients=bool(opposite))


def run_transport(base, cfg, seed):
    from openmm import unit
    from polypaves.sorption import load_network
    from workflows.network_properties.network_properties import write_json, write_csv
    directory = Path(base) / 'thermal_conductivity'
    directory.mkdir(exist_ok=True)
    host = load_network(Path(base) / 'relaxed_300K.structure.json')
    if not np.allclose(host.box_nm, np.diag(np.diag(host.box_nm))):
        raise ValueError('RNEMD currently requires an orthorhombic cell')
    slabs = cfg.get('transport_slabs', 20)
    if type(slabs) is not int or slabs < 12 or slabs % 2:
        raise ValueError('transport_slabs must be even and >= 12')
    dt = .25
    eq_ps = .002 if cfg['smoke'] else cfg.get('transport_equilibration_ps', 500.)
    warmup_ps = .004 if cfg['smoke'] else cfg.get('transport_gradient_ps', 1000.)
    production_ps = .02 if cfg['smoke'] else cfg.get('transport_production_ps', 5000.)
    swap_steps = 4 if cfg['smoke'] else cfg.get('transport_swap_steps', 1000)
    if type(swap_steps) is not int or swap_steps < 1:
        raise ValueError('transport_swap_steps must be a positive integer')
    axis_results = []
    for axis in cfg.get('transport_axes', [0,1,2]):
        if axis not in (0,1,2):
            raise ValueError('transport_axes entries must be 0, 1 or 2')
        path = directory / f'axis_{axis}.json'
        if path.exists():
            axis_results.append(json.loads(path.read_text()))
            continue
        print(f'RNEMD heat transport axis {axis}', flush=True)
        system, integrator, context = context_for(host, cfg, seed+axis, dt)
        integrator.step(max(1, round(eq_ps*1000/dt)))
        state = context.getState(getPositions=True, getVelocities=True)
        del context, integrator
        system, integrator, context = context_for(host, cfg, seed+axis, dt, nve=True)
        context.setState(state)
        mass = np.array([system.getParticleMass(i).value_in_unit(unit.dalton)
                         for i in range(system.getNumParticles())])
        unique, counts = np.unique(mass[mass > 0], return_counts=True)
        selected_mass = unique[counts.argmax()]
        selected = mass == selected_mass
        length = float(host.box_nm[axis,axis])
        area = float(np.linalg.det(host.box_nm) / length)
        profiles, log, energy = [], [], []
        transferred, accepted, attempted = 0., 0, 0
        warmup_steps = max(1, round(warmup_ps*1000/dt))
        total_steps = warmup_steps + max(1, round(production_ps*1000/dt))
        elapsed, production_steps = 0, 0
        try:
            while elapsed < total_steps:
                chunk = min(swap_steps, total_steps-elapsed,
                            warmup_steps-elapsed if elapsed < warmup_steps else total_steps-elapsed)
                integrator.step(chunk)
                elapsed += chunk
                in_production = elapsed > warmup_steps
                snapshot = context.getState(getPositions=True, getVelocities=True, getEnergy=True)
                xyz = snapshot.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                velocities = snapshot.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond)
                slab = np.floor((xyz[:,axis] / length % 1) * slabs).astype(int)
                ke = .5 * mass * np.sum(velocities**2, axis=1)
                hot = np.flatnonzero((slab == 0) & selected)
                cold = np.flatnonzero((slab == slabs//2) & selected)
                attempted += 1
                if len(hot) and len(cold):
                    i, j = hot[ke[hot].argmin()], cold[ke[cold].argmax()]
                    exchange = float(ke[j]-ke[i])
                    if exchange > 0:
                        velocities[[i,j]] = velocities[[j,i]]
                        context.setVelocities(velocities * unit.nanometer/unit.picosecond)
                        accepted += 1
                        if in_production:
                            transferred += exchange
                if in_production:
                    production_steps += chunk
                    # Sample before the next swap, with local streaming velocity removed.
                    profile = []
                    for index in range(slabs):
                        take = slab == index
                        m = mass[take]
                        if len(m) < 2:
                            profile.append(float('nan'))
                            continue
                        v = velocities[take]
                        v = v - np.average(v, axis=0, weights=m)
                        kinetic_J = .5 * np.sum(m[:,None]*v*v) * 1000 / NA
                        profile.append(float(2*kinetic_J/((3*len(m)-3)*KB)))
                    profiles.append(profile)
                    total_energy = (snapshot.getPotentialEnergy()+snapshot.getKineticEnergy()).value_in_unit(unit.kilojoule_per_mole)
                    energy.append(total_energy)
                    log.append(dict(time_ps=production_steps*dt/1000, transferred_kJ_mol=transferred,
                                    total_energy_kJ_mol=total_energy))
            with np.errstate(invalid='ignore'):
                mean_profile = np.nanmean(np.asarray(profiles), axis=0)
            estimate = conductivity(mean_profile, length, area, transferred, production_steps*dt/1000)
            block_values = []
            for block in np.array_split(np.arange(len(log)), 5):
                if len(block) < 3:
                    continue
                first, last = int(block[0]), int(block[-1])
                start_q = log[first-1]['transferred_kJ_mol'] if first else 0.
                start_t = log[first-1]['time_ps'] if first else 0.
                profile = np.nanmean(np.asarray(profiles)[block], axis=0)
                value = conductivity(profile, length, area,
                    log[last]['transferred_kJ_mol']-start_q, log[last]['time_ps']-start_t)
                if value['conductivity_W_m_K'] is not None:
                    block_values.append(value['conductivity_W_m_K'])
            estimate.update(axis=axis, smoke=cfg['smoke'], duration_ps=production_steps*dt/1000,
                profile_K=[float(t) if np.isfinite(t) else None for t in mean_profile],
                swaps_accepted=accepted, swaps_attempted=attempted,
                energy_drift_kJ_mol=float(energy[-1]-energy[0]),
                mean_temperature_K=float(np.nanmean(mean_profile)), sampling_checks_passed=False,
                block_conductivities_W_m_K=block_values,
                block_sem_W_m_K=float(np.std(block_values,ddof=1)/np.sqrt(len(block_values))) if len(block_values)>1 else None,
                limitations=['Unconstrained Sage, 0.25 fs NVE, equal-mass velocity swaps.',
                    'Check steady profiles, energy drift, flux dependence and cell-length dependence.',
                    'A short or nonlinear profile cannot establish thermal conductivity.'])
            write_csv(directory / f'axis_{axis}.csv', log)
            write_json(path, estimate)
            axis_results.append(estimate)
        finally:
            del context, integrator
    valid = [r['conductivity_W_m_K'] for r in axis_results if r['conductivity_W_m_K'] is not None]
    result = dict(smoke=cfg['smoke'], axes=axis_results,
        conductivity_direction_mean_W_m_K=float(np.mean(valid)) if len(valid)==len(axis_results) and valid else None,
        sampling_checks_passed=False)
    dry = json.loads((Path(base) / 'dry_result.json').read_text())['npt']['estimates']
    k = result['conductivity_direction_mean_W_m_K']
    cp, rho = dry['Cp_J_kg_K'], dry['density_g_cm3'] * 1000
    result['thermal_diffusivity_m2_s'] = k / (cp * rho) if k is not None and cp > 0 else None
    result['thermal_diffusivity_note'] = 'Uses constrained-model Cp; assess sensitivity to bond constraints.'
    write_json(directory / 'result.json', result)
