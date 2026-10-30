"""Finite-temperature NVT stress/strain response of the relaxed network.

The network is loaded without X-H constraints, with those bonds as harmonic
terms, so virtual affine strains include every bonded virial consistently.
Full potential energy differences include PME.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np

CONVERSION_GPA = .0016605390671738467  # kJ mol^-1 nm^-3 -> GPa
VOIGT = ('xx', 'yy', 'zz', 'yz', 'xz', 'xy')


def bases():
    matrices = []
    for i, j in ((0,0),(1,1),(2,2),(2,1),(2,0),(1,0)):
        b = np.zeros((3,3))
        b[i,j] = 1.
        matrices.append(b)
    return matrices


def context_for(host, cfg, seed, timestep_fs=.25, nve=False, allow_shear=False):
    import openmm as mm
    from openmm import unit
    system = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(host.system))
    if system.getNumConstraints():
        # OpenFF writes no bond term for a constrained bond; dropping the
        # constraint would leave that hydrogen unbonded (about -0.7 GPa of
        # spurious stress in a 20k-atom epoxy). Load with constrained=False.
        raise ValueError('elastic MD needs a network loaded with load_network(..., constrained=False)')
    if allow_shear:
        # CUDA chooses rectangular/triclinic kernels when the Context is created.
        # Initialize general-cell support, then restore the actual box below.
        default_box = np.array(host.box_nm, copy=True)
        if np.allclose(default_box, np.diag(np.diag(default_box))):
            default_box[1,0] = 1e-4 * default_box[1,1]
        system.setDefaultPeriodicBoxVectors(*(default_box * unit.nanometer))
    if nve:
        # Full-step velocities make an equal-mass swap conserve instantaneous KE.
        integrator = mm.CustomIntegrator(timestep_fs * unit.femtosecond)
        integrator.addUpdateContextState()
        integrator.addComputePerDof('v', 'v+0.5*dt*f/m')
        integrator.addComputePerDof('x', 'x+dt*v')
        integrator.addComputePerDof('v', 'v+0.5*dt*f/m')
    else:
        integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                                timestep_fs * unit.femtosecond)
    if not nve:
        integrator.setRandomNumberSeed(seed)
    properties = {'Precision':'mixed'} if cfg['platform'] == 'CUDA' else {}
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName(cfg['platform']), properties)
    context.setPeriodicBoxVectors(*(host.box_nm * unit.nanometer))
    context.setPositions(host.positions_nm * unit.nanometer)
    context.setVelocitiesToTemperature(300 * unit.kelvin, seed)
    return system, integrator, context


def stress(context, system, delta=1e-4):
    """Tension-positive Cauchy stress in Voigt order (GPa)."""
    from openmm import unit
    state = context.getState(getPositions=True, getVelocities=True)
    xyz = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
    velocities = state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)
    mass = np.array([system.getParticleMass(i).value_in_unit(unit.dalton)
                     for i in range(system.getNumParticles())])
    kinetic = np.einsum('n,ni,nj->ij', mass, velocities, velocities)
    values = []
    try:
        for b in bases():
            energies = []
            for sign in (-1,1):
                transform = np.eye(3) + sign * delta * b
                context.setPeriodicBoxVectors(*((box @ transform) * unit.nanometer))
                context.setPositions((xyz @ transform) * unit.nanometer)
                energies.append(context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
            values.append((energies[1] - energies[0]) / (2 * delta) - np.sum(kinetic * b))
    finally:
        context.setPeriodicBoxVectors(*(box * unit.nanometer))
        context.setPositions(xyz * unit.nanometer)
    return np.asarray(values) / np.linalg.det(box) * CONVERSION_GPA


def isotropic_moduli(stiffness):
    """Voigt/Reuss/Hill bounds for an anisotropic stiffness matrix."""
    c = np.asarray(stiffness, float)
    if np.linalg.eigvalsh(c).min() <= 0:
        return dict(stable=False, reason='Nonpositive stiffness eigenvalue; more sampling/strain checks needed.')
    s = np.linalg.inv(c)
    kv = (np.trace(c[:3,:3]) + 2 * (c[0,1]+c[0,2]+c[1,2])) / 9
    gv = (np.trace(c[:3,:3]) - c[0,1]-c[0,2]-c[1,2] + 3*np.trace(c[3:,3:])) / 15
    kr = 1 / (np.trace(s[:3,:3]) + 2 * (s[0,1]+s[0,2]+s[1,2]))
    gr = 15 / (4*np.trace(s[:3,:3]) - 4*(s[0,1]+s[0,2]+s[1,2]) + 3*np.trace(s[3:,3:]))
    k, g = (kv+kr)/2, (gv+gr)/2
    return dict(stable=True, bulk_Voigt_GPa=float(kv), bulk_Reuss_GPa=float(kr),
        shear_Voigt_GPa=float(gv), shear_Reuss_GPa=float(gr), bulk_Hill_GPa=float(k),
        shear_Hill_GPa=float(g), young_Hill_GPa=float(9*k*g/(3*k+g)),
        poisson_Hill=float((3*k-2*g)/(2*(3*k+g))))


def run_mechanical(base, cfg, seed):
    import openmm as mm
    from openmm import unit
    from polypaves.sorption import load_network
    from workflows.network_properties.network_properties import write_json, write_csv
    directory = Path(base) / 'mechanical'
    directory.mkdir(exist_ok=True)
    if (directory / 'result.json').exists():
        return
    host = load_network(Path(base) / 'relaxed_300K.structure.json', constrained=False)
    strain = cfg.get('elastic_strain', .01)
    if not 0 < strain < .03:
        raise ValueError('elastic_strain must lie between 0 and 0.03')
    timestep = .25
    equil_ps = .001 if cfg['smoke'] else cfg.get('elastic_equilibration_ps', 200.)
    production_ps = .004 if cfg['smoke'] else cfg.get('elastic_production_ps', 1000.)
    interval_ps = cfg['sample_interval_ps']
    means, sems = {}, {}
    for column, basis in enumerate(bases()):
        for sign in (-1,1):
            name = f'{VOIGT[column]}_{"plus" if sign > 0 else "minus"}'
            path = directory / f'{name}.json'
            if path.exists():
                import json
                rows = json.loads(path.read_text())
            else:
                print(f'Elastic response {name}, strain={sign*strain:g}', flush=True)
                system, integrator, context = context_for(host, cfg, seed+column, timestep, allow_shear=True)
                try:
                    transform = np.eye(3) + sign * strain * basis
                    context.setPeriodicBoxVectors(*((host.box_nm @ transform) * unit.nanometer))
                    context.setPositions((host.positions_nm @ transform) * unit.nanometer)
                    integrator.step(max(1, round(equil_ps*1000/timestep)))
                    rows = []
                    for index in range(max(4, round(production_ps/interval_ps))):
                        integrator.step(max(1, round(interval_ps*1000/timestep)))
                        rows.append(dict(time_ps=(index+1)*interval_ps,
                                         **dict(zip(VOIGT, stress(context, system).tolist()))))
                    write_json(path, rows)
                    write_csv(directory / f'{name}.csv', rows)
                finally:
                    del context, integrator
            samples = np.array([[row[k] for k in VOIGT] for row in rows])
            means[column,sign] = samples.mean(axis=0)
            blocks = [b.mean(axis=0) for b in np.array_split(samples, min(5,len(samples)))]
            sems[column,sign] = np.std(blocks, axis=0, ddof=1) / np.sqrt(len(blocks))
    c = np.column_stack([(means[i,1]-means[i,-1])/(2*strain) for i in range(6)])
    errors = np.column_stack([np.sqrt(sems[i,1]**2+sems[i,-1]**2)/(2*strain) for i in range(6)])
    symmetric = (c+c.T)/2
    write_json(directory / 'result.json', dict(smoke=cfg['smoke'], temperature_K=300.,
        voigt_order=list(VOIGT), engineering_shear=True, strain=strain,
        stiffness_GPa=c, symmetric_stiffness_GPa=symmetric, block_sem_GPa=errors,
        asymmetry_relative=float(np.linalg.norm(c-c.T)/max(np.linalg.norm(c),1e-20)),
        isotropic=isotropic_moduli(symmetric), sampling_checks_passed=False,
        limitations=['Unconstrained Sage Hamiltonian, 0.25 fs; differs from constrained NPT heat capacity.',
            'Finite-time isothermal tangent response of the sampled glass.',
            'Repeat with smaller/larger strains and longer trajectories to establish convergence.',
            'Hill moduli assume an isotropic aggregate. No convergence claim is made automatically.']))
