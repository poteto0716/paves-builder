"""SI-unit estimators for a fixed-topology periodic network.

Fluctuation relations follow Hayashi et al., npj Comput. Mater. 8, 222 (2022).
Errors between contiguous blocks are diagnostics, not independent replicas.
"""
from __future__ import annotations

import numpy as np

KB = 1.380649e-23
NA = 6.02214076e23
DA = 1e-3 / NA
EPS0 = 8.8541878128e-12
ELEMENTARY_CHARGE = 1.602176634e-19


def effective_samples(values):
    """FFT autocorrelation, truncated at its first nonpositive lag."""
    x = np.asarray(values, float)
    if len(x) < 4 or np.ptp(x) == 0:
        return 0.0
    x = x - x.mean()
    if len(x) < 4 or np.dot(x, x) == 0:
        return 0.0
    f = np.fft.rfft(x, n=2 * len(x))
    ac = np.fft.irfft(f * f.conjugate())[:len(x)]
    ac /= np.arange(len(x), 0, -1)
    ac /= ac[0]
    stop = np.flatnonzero(ac[1:] <= 0)
    end = int(stop[0] + 1) if len(stop) else len(x)
    return float(len(x) / max(1., 1 + 2 * ac[1:end].sum()))


def thermodynamics(rows, mass_Da, temperature_K, pressure_atm):
    """NPT Cp, Cv, expansion, compressibilities, moduli and sound speeds."""
    if len(rows) < 4 or mass_Da <= 0 or temperature_K <= 0:
        raise ValueError('need at least four samples and positive mass/temperature')
    v = np.array([r['volume_nm3'] for r in rows]) * 1e-27
    energy = np.array([r['potential_kJ_mol'] + r['kinetic_kJ_mol'] for r in rows]) * 1000 / NA
    if not np.isfinite(v).all() or not np.isfinite(energy).all() or np.any(v <= 0):
        raise ValueError('nonfinite energies or nonpositive volumes')
    h = energy + pressure_atm * 101325 * v
    mass = mass_Da * DA
    rho = mass / v.mean()
    # Subtract a reference before centering to preserve exactly constant traces.
    dv, dh = v - v[0], h - h[0]
    dv, dh = dv - dv.mean(), dh - dh.mean()
    kappa = float(np.mean(dv * dv) / (KB * temperature_K * v.mean()))
    cp = float(np.mean(dh * dh) / (KB * temperature_K**2 * mass))
    alpha = float(np.mean(dv * dh) / (KB * temperature_K**2 * v.mean()))
    cv = cp - temperature_K * alpha**2 / (rho * kappa) if kappa > 0 else None
    ks = kappa * cv / cp if cv is not None and cp > 0 else None
    return dict(density_g_cm3=float(rho / 1000),
        Cp_J_kg_K=cp, Cv_J_kg_K=cv,
        compressibility_isothermal_1_Pa=kappa,
        bulk_modulus_isothermal_GPa=1e-9 / kappa if kappa > 0 else None,
        alpha_volume_1_K=alpha, alpha_linear_1_K=alpha / 3,
        compressibility_isentropic_1_Pa=ks,
        bulk_modulus_isentropic_GPa=1e-9 / ks if ks and ks > 0 else None,
        bulk_sound_speed_m_s=float(np.sqrt(1 / (ks * rho))) if ks and ks > 0 else None)


def analyze_npt(rows, mass_Da, temperature_K=300., pressure_atm=1., min_effective_samples=100):
    estimate = thermodynamics(rows, mass_Da, temperature_K, pressure_atm)
    blocks = [thermodynamics(list(b), mass_Da, temperature_K, pressure_atm)
              for b in np.array_split(np.asarray(rows, object), 5) if len(b) >= 4]
    errors = {}
    for key in estimate:
        values = [b[key] for b in blocks if b[key] is not None and np.isfinite(b[key])]
        errors[key] = float(np.std(values, ddof=1) / np.sqrt(len(values))) if len(values) > 1 else None
    diagnostics = {}
    for key in ('volume_nm3', 'potential_kJ_mol'):
        x = np.array([r[key] for r in rows])
        drift = abs(x[len(x)//2:].mean() - x[:len(x)//2].mean()) / max(x.std(), 1e-20)
        diagnostics[key] = dict(effective_samples=effective_samples(x), half_drift_in_std=float(drift))
    reliable = (len(blocks) == 5 and
                all(d['effective_samples'] >= min_effective_samples and d['half_drift_in_std'] < .5
                    for d in diagnostics.values()) and
                estimate['Cv_J_kg_K'] is not None and estimate['Cv_J_kg_K'] > 0)
    return dict(estimates=estimate, block_sem=errors, blocks=len(blocks), samples=len(rows),
                diagnostics=diagnostics, sampling_checks_passed=bool(reliable),
                limitations=['Classical force-field heat capacities; no quantum correction.',
                    'Values apply to the constrained Sage model and sampled metastable glass.',
                    'Bulk sound speed is sqrt(Ks/rho), not the solid longitudinal wave speed.',
                    'Block errors omit force-field, finite-size and cooling-rate errors.'])


def analyze_nvt(rows, mass_Da, temperature_K=300., min_effective_samples=100):
    """Fixed-charge dielectric tensor from continuous itinerant polarization.

    Dipoles are cumulative atomic charge displacements at fixed volume, in e nm.
    PME conducting-boundary fluctuation relation; electronic response is absent.
    """
    energy = np.array([r['potential_kJ_mol']+r['kinetic_kJ_mol'] for r in rows])*1000/NA
    moment = np.array([[r[f'dipole_{a}_e_nm'] for a in 'xyz'] for r in rows])*ELEMENTARY_CHARGE*1e-9
    volumes = np.array([r['volume_nm3'] for r in rows])*1e-27
    if len(rows) < 4 or not np.allclose(volumes, volumes[0], rtol=1e-10, atol=0):
        raise ValueError('NVT estimator needs >= 4 samples at fixed volume')
    centered = moment-moment.mean(axis=0)
    tensor = np.eye(3) + centered.T @ centered / (len(rows)*EPS0*KB*temperature_K*volumes[0])
    cv = np.var(energy)/(KB*temperature_K**2*mass_Da*DA)
    return dict(Cv_direct_J_kg_K=float(cv), dielectric_tensor_fixed_charge=tensor.tolist(),
        dielectric_isotropic_fixed_charge=float(np.trace(tensor)/3),
        effective_samples_energy=effective_samples(energy),
        effective_samples_dipole=[effective_samples(moment[:,i]) for i in range(3)],
        sampling_checks_passed=False,
        limitations=['Fixed-charge, conducting PME boundaries; electronic polarization is omitted.',
                    'Continuous charge displacements avoid branch jumps at periodic faces.',
                    'Frozen orientational modes in glass can prevent dielectric convergence.',
                    'Sampling checks alone cannot establish full dielectric relaxation.'])


def topology(graphs, positions_nm, box_nm):
    """Connected components and periodic winding of minimum-image bond cycles."""
    xyz = np.asarray(positions_nm)
    inverse = np.linalg.inv(box_nm)
    adjacent, masses, offset, edges = [[] for _ in xyz], [], 0, 0
    for graph in graphs:
        masses.extend(float(a['mass']) for a in graph['atoms'])
        for bond in graph['bonds']:
            i, j = int(bond[0])+offset, int(bond[1])+offset
            shift = -np.rint((xyz[j]-xyz[i]) @ inverse).astype(int)
            adjacent[i].append((j,shift))
            adjacent[j].append((i,-shift))
            edges += 1
        offset += len(graph['atoms'])
    if offset != len(xyz):
        raise ValueError('graph atom count does not match coordinates')
    lifts, components = {}, []
    for root in range(len(xyz)):
        if root in lifts:
            continue
        lifts[root] = np.zeros(3,dtype=int)
        pending, nodes, winding = [root], [], np.zeros(3,dtype=bool)
        while pending:
            i = pending.pop()
            nodes.append(i)
            for j,shift in adjacent[i]:
                expected = lifts[i]+shift
                if j not in lifts:
                    lifts[j] = expected
                    pending.append(j)
                else:
                    winding |= expected != lifts[j]
        components.append(dict(atoms=len(nodes), mass_Da=float(np.asarray(masses)[nodes].sum()),
                               periodic_winding_axes=winding.tolist()))
    components.sort(key=lambda c:c['mass_Da'],reverse=True)
    return dict(atoms=len(xyz), bonds=edges, connected_components=len(components),
        largest_component_mass_fraction=components[0]['mass_Da']/sum(masses),
        graph_cycle_rank=edges-len(xyz)+len(components), components=components,
        percolates_periodically=any(any(c['periodic_winding_axes']) for c in components),
        limitation='Periodic winding assumes every bond uses its minimum-image displacement.')


def diffusion(time_ps, coordinates_nm, min_duration_ps=1000.):
    """Delegate unwrapped, host-relative COM trajectories to the wheel estimator."""
    from polypaves.sorption_workflow import diffusion_statistics
    return diffusion_statistics(time_ps, coordinates_nm, min_duration_ps=min_duration_ps)


def aggregate(results):
    """Flatten scalar estimates and report sample SD / SEM across thermal seeds."""
    paths = {}
    def visit(value, prefix=''):
        for key, val in value.items():
            path = f'{prefix}.{key}' if prefix else key
            if isinstance(val, dict):
                visit(val, path)
            elif isinstance(val, (int, float)) and not isinstance(val, bool) and np.isfinite(val):
                paths.setdefault(path, []).append(float(val))
    for result in results:
        visit(result)
    return {key: dict(n=len(x), mean=float(np.mean(x)),
                      std=float(np.std(x, ddof=1)) if len(x) > 1 else None,
                      sem=float(np.std(x, ddof=1) / np.sqrt(len(x))) if len(x) > 1 else None)
            for key, x in paths.items()}
