"""Unit/physics checks; installed-wheel integration is opt-in."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from workflows.network_properties.analysis import thermodynamics, analyze_npt, analyze_nvt, aggregate, topology, KB, NA, DA
from workflows.network_properties.mechanical import isotropic_moduli, stress, context_for, CONVERSION_GPA
from workflows.network_properties.network_properties import DEFAULTS, temperature_grid, validate


class ThermodynamicsTests(unittest.TestCase):
    def test_fluctuation_units_and_thermodynamic_identity(self):
        # Orthogonal +/- fluctuations have an exactly specified covariance.
        a = np.tile([-1., -1., 1., 1.], 100)
        b = np.tile([-1., 1., -1., 1.], 100)
        volume = 10e-27 + a*1e-28
        enthalpy = 1e-18 + a*1e-20 + b*2e-20
        energy = enthalpy - 101325*volume
        rows = [dict(volume_nm3=v/1e-27, potential_kJ_mol=e*NA/1000, kinetic_kJ_mol=0.)
                for v,e in zip(volume, energy)]
        result = thermodynamics(rows, 10000, 300, 1.)
        kappa = 1e-56/(KB*300*1e-26)
        cp = 5e-40/(KB*300**2*(10000*DA))
        alpha = 1e-48/(KB*300**2*1e-26)
        self.assertAlmostEqual(result['bulk_modulus_isothermal_GPa'], 1e-9/kappa)
        self.assertAlmostEqual(result['Cp_J_kg_K'], cp)
        self.assertAlmostEqual(result['alpha_volume_1_K'], alpha)
        rho = 10000*DA/1e-26
        self.assertAlmostEqual(result['Cv_J_kg_K'], cp-300*alpha**2/(rho*kappa))
        self.assertAlmostEqual(result['bulk_modulus_isentropic_GPa'],
                               result['bulk_modulus_isothermal_GPa']*cp/result['Cv_J_kg_K'])

    def test_no_volume_fluctuation_does_not_report_infinite_modulus(self):
        rows = [dict(volume_nm3=10., potential_kJ_mol=float(i), kinetic_kJ_mol=1.) for i in range(20)]
        result = analyze_npt(rows, 10000)
        self.assertIsNone(result['estimates']['bulk_modulus_isothermal_GPa'])
        self.assertFalse(result['sampling_checks_passed'])

    def test_nvt_dipole_offsets_do_not_change_dielectric_response(self):
        rows = [dict(volume_nm3=10., potential_kJ_mol=float(i%3), kinetic_kJ_mol=1.,
                     dipole_x_e_nm=float(i%2), dipole_y_e_nm=0., dipole_z_e_nm=0.) for i in range(20)]
        first = analyze_nvt(rows, 10000)
        shifted = copy.deepcopy(rows)
        for r in shifted:
            r['dipole_x_e_nm'] += 123.
        self.assertAlmostEqual(first['dielectric_isotropic_fixed_charge'],
                               analyze_nvt(shifted, 10000)['dielectric_isotropic_fixed_charge'])
        self.assertGreater(first['dielectric_tensor_fixed_charge'][0][0], 1)
        self.assertEqual(first['dielectric_tensor_fixed_charge'][1][1], 1)

    def test_nvt_rejects_changing_volume(self):
        rows = [dict(volume_nm3=float(i+1), potential_kJ_mol=0., kinetic_kJ_mol=1.,
                     dipole_x_e_nm=0., dipole_y_e_nm=0., dipole_z_e_nm=0.) for i in range(5)]
        with self.assertRaises(ValueError):
            analyze_nvt(rows, 10000)

    def test_replica_singleton_has_no_spurious_zero_uncertainty(self):
        self.assertIsNone(aggregate([{'K':2.}])['K']['sem'])
        self.assertAlmostEqual(aggregate([{'K':1.},{'K':3.}])['K']['std'], np.sqrt(2))


class MechanicalTests(unittest.TestCase):
    @unittest.skipUnless(os.environ.get('PAVES_CUDA_TESTS') == '1', 'set PAVES_CUDA_TESTS=1 for CUDA checks')
    def test_cuda_context_allows_virtual_shear_of_rectangular_cell(self):
        import openmm as mm
        from types import SimpleNamespace
        system=mm.System()
        system.addParticle(12.);system.addParticle(12.)
        force=mm.HarmonicBondForce()
        force.addBond(0,1,.2,100.)
        force.setUsesPeriodicBoundaryConditions(True)
        system.addForce(force)
        host=SimpleNamespace(system=system,box_nm=np.eye(3)*3.,
                             positions_nm=np.array([[.1,.1,.1],[.3,.3,.3]]))
        system,integrator,context=context_for(host,{'platform':'CUDA'},31415,allow_shear=True)
        try:
            self.assertTrue(np.isfinite(stress(context,system)).all())
        finally:
            del context,integrator

    @unittest.skipUnless(importlib.util.find_spec('openmm'), 'OpenMM is required')
    def test_virtual_strain_matches_analytic_periodic_bond_and_kinetic_stress(self):
        import openmm as mm
        from openmm import unit
        system = mm.System()
        for _ in range(2):
            system.addParticle(12.)
        force = mm.HarmonicBondForce()
        force.setUsesPeriodicBoundaryConditions(True)
        force.addBond(0,1,.2,100.)
        system.addForce(force)
        box = np.eye(3)*3.
        system.setDefaultPeriodicBoxVectors(*(box*unit.nanometer))
        integrator = mm.VerletIntegrator(.00025)
        context = mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference'))
        positions = np.array([[.1,.1,.1],[2.9,2.9,2.9]])
        velocities = np.array([[1.,2.,3.],[-1.,-2.,-3.]])
        context.setPositions(positions*unit.nanometer)
        context.setVelocities(velocities*unit.nanometer/unit.picosecond)
        try:
            measured = stress(context,system,1e-5)
            r = np.ones(3)*.2
            distance = np.linalg.norm(r)
            tensor = (100*(distance-.2)/distance*np.outer(r,r) -
                      24*np.outer(velocities[0],velocities[0]))/27*CONVERSION_GPA
            expected = tensor[[0,1,2,1,0,0],[0,1,2,2,2,1]]
            np.testing.assert_allclose(measured,expected,rtol=1e-8,atol=1e-10)
            state = context.getState(getPositions=True)
            np.testing.assert_array_equal(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer),positions)
            np.testing.assert_array_equal(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer),box)
        finally:
            del context,integrator

    def test_isotropic_elastic_tensor_and_engineering_shear(self):
        k, g = 4., 1.5
        c = np.eye(6)*g
        c[:3,:3] = k-2*g/3
        np.fill_diagonal(c[:3,:3], k+4*g/3)
        r = isotropic_moduli(c)
        self.assertAlmostEqual(r['bulk_Hill_GPa'], k)
        self.assertAlmostEqual(r['shear_Hill_GPa'], g)
        self.assertAlmostEqual(r['young_Hill_GPa'], 9*k*g/(3*k+g))
        self.assertFalse(isotropic_moduli(-np.eye(6))['stable'])


class ConfigurationTests(unittest.TestCase):
    def test_periodic_winding_distinguishes_network_from_finite_ring(self):
        graphs = [dict(atoms=[{'mass':12.}]*3,bonds=[[0,1,1],[1,2,1],[2,0,1]])]
        spanning = topology(graphs,np.array([[.1,0,0],[.43,0,0],[.76,0,0]]),np.eye(3))
        finite = topology(graphs,np.array([[.1,0,0],[.2,0,0],[.3,0,0]]),np.eye(3))
        self.assertTrue(spanning['percolates_periodically'])
        self.assertEqual(spanning['components'][0]['periodic_winding_axes'],[True,False,False])
        self.assertFalse(finite['percolates_periodically'])
        self.assertEqual(spanning['graph_cycle_rank'],1)

    def test_cooling_grid_brackets_300_and_includes_100(self):
        grid = temperature_grid(800,100,37)
        self.assertEqual(grid[0],800)
        self.assertEqual(grid[-1],100)
        self.assertIn(300,grid)
        self.assertEqual(grid,sorted(set(grid),reverse=True))

    def test_rejects_bad_configuration_before_running(self):
        for patch in ({'replicates':0}, {'timestep_fs':float('nan')},
                      {'elastic_strain':-.01}, {'production_ps':1},
                      {'solvents':[{'name':'../escape','smiles':'O'}]}):
            with self.subTest(patch=patch), self.assertRaises(ValueError):
                validate(dict(DEFAULTS, **patch))


@unittest.skipUnless(os.environ.get('PAVES_NETWORK_INTEGRATION') == '1',
                     'set PAVES_NETWORK_INTEGRATION=1 with the OpenFF / public-wheel environment')
class InstalledWheelIntegration(unittest.TestCase):
    def test_epoxy_all_properties_and_resume(self):
        import shutil
        env = dict(os.environ, OPENMM_CPU_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        def run(args, cwd):
            result = subprocess.run([sys.executable, *map(str,args)], cwd=cwd, env=env,
                                    capture_output=True, text=True, timeout=300)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            for name in ('activation.yaml','reactions.yaml'):
                shutil.copy2(ROOT/'workflows/crosslinking'/name,p/name)
            run([ROOT/'workflows/crosslinking/crosslinking.py','--project','epoxy',
                 '--total-atoms','500','--activation-count','1','--platform','CPU',
                 '--compression-steps','2','--equilibration-steps','2','--cycles','1','--nvt-steps','2','--npt-steps','2',
                 '--max-reactions-per-cycle','1','--reaction-cutoff','10','--max-reaction-cutoff','10',
                 '--minimization-iterations','100'],p)
            source = p/'results/epoxy/final.structure.json'
            original = source.read_bytes()
            script = ROOT/'workflows/network_properties/network_properties.py'
            run([script,'new',p/'evaluation','--network',source,'--platform','CPU',
                 '--solvent','water=O','--solvent-molecules','1','--smoke'],p)
            run([script,'run',p/'evaluation'],p)
            project = p/'evaluation'
            base = project/'replica_01'
            relaxed = json.loads((base/'relaxed_300K.structure.json').read_text())
            self.assertEqual(json.loads(original)['graphs'],relaxed['graphs'])
            self.assertEqual(source.read_bytes(), original)
            for name in ('cooling_300K.structure.json','relaxed_100K.structure.json','dry_result.json',
                         'mechanical/result.json'):
                self.assertTrue((base/name).is_file(),name)
            wet = json.loads((base/'solvents/water/result.json').read_text())
            self.assertFalse(wet['infinite_dilution_free_energy']['convergence_verified'])
            self.assertIsNone(wet['fixed_loading']['diffusion']['diffusion_coefficient_m2_s'])
            before = (project/'summary.json').read_bytes()
            run([script,'run',project],p)
            self.assertEqual((project/'summary.json').read_bytes(),before)
            self.assertTrue((project/'report.md').is_file())


if __name__ == '__main__':
    unittest.main()
