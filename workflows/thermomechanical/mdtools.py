"""Small OpenMM helpers used by the thermomechanical workflow.

All public lengths are Angstrom and times are femtoseconds unless their names
state otherwise.  OpenMM internally uses nm, ps, and kJ/mol.
"""
from __future__ import annotations

import math
from collections import deque
from pathlib import Path

import numpy as np
import openmm as mm
import openmm.unit as u

KB_KJ = 0.00831446261815324
AVOGADRO = 6.02214076e23
ATM_TO_BAR = 1.01325
KJ_MOL_NM3_TO_MPA = 1.6605390671738467
H_MASS_MAX = 1.2
_PLATFORM = {}


def load_xml(path):
    return mm.XmlSerializer.deserialize(Path(path).read_text())


def save_state(context, path):
    state = context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=False)
    Path(path).write_text(mm.XmlSerializer.serialize(state))
    return state


def positions_A(state):
    return np.asarray(state.getPositions(asNumpy=True).value_in_unit(u.angstrom))


def box_vectors_nm(state):
    return np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.nanometer))


def box_A(state):
    vectors = np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.angstrom))
    return np.asarray([np.linalg.norm(v) for v in vectors])


def volume_A3(state):
    vectors = np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.angstrom))
    return float(abs(np.linalg.det(vectors)))


def masses(system):
    return np.asarray([system.getParticleMass(i).value_in_unit(u.dalton)
                       for i in range(system.getNumParticles())])


def density_g_cm3(system, state):
    return float(masses(system).sum() / AVOGADRO / (volume_A3(state) * 1e-24))


def bond_force(system):
    for force in system.getForces():
        if isinstance(force, mm.CustomBondForce):
            names = [force.getPerBondParameterName(k) for k in range(force.getNumPerBondParameters())]
            if 'r0' in names:
                return force
    raise ValueError('no PCFF Class-II CustomBondForce with an r0 parameter')


def bonds(system):
    force = bond_force(system)
    names = [force.getPerBondParameterName(k) for k in range(force.getNumPerBondParameters())]
    r0 = names.index('r0')
    return [(*force.getBondParameters(k)[:2], force.getBondParameters(k)[2][r0])
            for k in range(force.getNumBonds())]


def configure(system, shake_hydrogen=True):
    """Add X-H constraints at the PCFF equilibrium bond lengths."""
    if not shake_hydrogen:
        return 0
    mass = masses(system)
    added = 0
    for i, j, r0 in bonds(system):
        if mass[i] < H_MASS_MAX or mass[j] < H_MASS_MAX:
            system.addConstraint(i, j, r0)
            added += 1
    return added


def degrees_of_freedom(system):
    mobile = sum(system.getParticleMass(i).value_in_unit(u.dalton) > 0
                 for i in range(system.getNumParticles()))
    return 3 * mobile - system.getNumConstraints() - 3


def measured_temperature(context, system):
    energy = context.getState(getEnergy=True).getKineticEnergy().value_in_unit(u.kilojoule_per_mole)
    return 2.0 * energy / (degrees_of_freedom(system) * KB_KJ)


def unwrap(positions, box, bond_list):
    """Make each chain whole in an orthorhombic periodic box."""
    pos = np.array(positions, float)
    box = np.asarray(box, float)
    adjacent = [[] for _ in range(len(pos))]
    for i, j, *_ in bond_list:
        adjacent[i].append(j)
        adjacent[j].append(i)
    seen = np.zeros(len(pos), bool)
    for first in range(len(pos)):
        if seen[first]:
            continue
        seen[first] = True
        queue = deque([first])
        while queue:
            i = queue.popleft()
            for j in adjacent[i]:
                if seen[j]:
                    continue
                delta = pos[j] - pos[i]
                delta -= box * np.round(delta / box)
                pos[j] = pos[i] + delta
                seen[j] = True
                queue.append(j)
    return pos


def limited_langevin(temperature_K, friction_per_ps, dt_fs, max_step_A):
    """BAOAB integrator with capped displacement for fresh packed structures."""
    dt_ps = dt_fs * 1e-3
    vmax = max_step_A * 0.1 / dt_ps
    a = math.exp(-friction_per_ps * dt_ps)
    integrator = mm.CustomIntegrator(dt_ps)
    integrator.addGlobalVariable('a', a)
    integrator.addGlobalVariable('b', math.sqrt(1.0 - a * a))
    integrator.addGlobalVariable('kT', KB_KJ * temperature_K)
    integrator.addGlobalVariable('vmax', vmax)
    integrator.addUpdateContextState()
    kick = 'select(m, max(-vmax, min(vmax, v + 0.5*dt*f/m)), 0)'
    integrator.addComputePerDof('v', kick)
    integrator.addComputePerDof('x', 'x + 0.5*dt*v')
    integrator.addComputePerDof('v', 'select(m, a*v + b*sqrt(kT/m)*gaussian, 0)')
    integrator.addComputePerDof('x', 'x + 0.5*dt*v')
    integrator.addComputePerDof('v', kick)
    return integrator


def pick_platform(requested='auto'):
    """Select a working platform, probing CUDA/OpenCL before falling back."""
    if requested != 'auto':
        return requested
    if requested in _PLATFORM:
        return _PLATFORM[requested]
    available = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    failures = []
    for name in ('CUDA', 'OpenCL', 'CPU', 'Reference'):
        if name not in available:
            continue
        try:
            system = mm.System()
            system.addParticle(1.0)
            system.addParticle(1.0)
            force = mm.HarmonicBondForce()
            force.addBond(0, 1, 0.1, 1.0)
            system.addForce(force)
            context = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName(name))
            context.setPositions([mm.Vec3(), mm.Vec3(0.1, 0, 0)])
            context.getState(getEnergy=True)
            del context
            _PLATFORM[requested] = name
            if failures:
                print(f'[platform] using {name}; unavailable: ' + '; '.join(failures), flush=True)
            return name
        except Exception as exc:  # pragma: no cover - hardware dependent
            failures.append(f'{name} ({str(exc).splitlines()[0]})')
    raise RuntimeError('no working OpenMM platform: ' + '; '.join(failures))


def make_context(system, integrator, platform='auto', precision='mixed'):
    name = pick_platform(platform)
    properties = {'Precision': precision} if name in ('CUDA', 'OpenCL') else {}
    return mm.Context(system, integrator, mm.Platform.getPlatformByName(name), properties)


def prerelax(system_xml, state, cfg, platform='auto', precision='mixed', report=None):
    """Displacement-limited relaxation, returning positions in Angstrom."""
    system = load_xml(system_xml)
    integrator = limited_langevin(cfg['temperature_K'], cfg['friction_per_ps'],
                                  cfg['dt_fs'], cfg['max_step_A'])
    context = make_context(system, integrator, platform, precision)
    vectors = state.getPeriodicBoxVectors()
    context.setPeriodicBoxVectors(*vectors)
    whole = unwrap(positions_A(state), box_A(state), bonds(system))
    context.setPositions(whole * 0.1)
    context.setVelocitiesToTemperature(cfg['temperature_K'] * u.kelvin, cfg.get('seed', 1))
    steps = int(cfg['steps'])
    done = 0
    while done < steps:
        n = min(1000, steps - done)
        integrator.step(n)
        done += n
        if report:
            report(done)
    return positions_A(context.getState(getPositions=True))


def axial_stress_MPa(context, system, axis, delta, temperature_K):
    """Instantaneous tensile stress from a central box-strain energy derivative.

    Scaling coordinates and the periodic vector together includes real-space and
    reciprocal-space PCFF contributions.  The ideal kinetic term is included;
    the returned sign is positive in tension.
    """
    state = context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=False)
    pos = np.asarray(state.getPositions(asNumpy=True).value_in_unit(u.nanometer))
    vectors = box_vectors_nm(state)
    volume = abs(float(np.linalg.det(vectors)))

    def energy_at(strain):
        p = pos.copy()
        b = vectors.copy()
        p[:, axis] *= 1.0 + strain
        b[axis] *= 1.0 + strain
        context.setPeriodicBoxVectors(*[mm.Vec3(*row) for row in b])
        context.setPositions(p)
        return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)

    plus = energy_at(delta)
    minus = energy_at(-delta)
    context.setPeriodicBoxVectors(*[mm.Vec3(*row) for row in vectors])
    context.setPositions(pos)
    derivative = (plus - minus) / (2.0 * delta)
    mobile = sum(system.getParticleMass(i).value_in_unit(u.dalton) > 0
                 for i in range(system.getNumParticles()))
    ideal = mobile * KB_KJ * temperature_K
    return (derivative - ideal) / volume * KJ_MOL_NM3_TO_MPA
