"""Unit/physics checks; installed-wheel integration is opt-in."""
import contextlib
import copy
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from workflows.network_properties.analysis import thermodynamics, analyze_npt, analyze_nvt, aggregate, topology, KB, NA, DA
from workflows.network_properties.mechanical import isotropic_moduli, stress, context_for, CONVERSION_GPA
from workflows.network_properties.network_properties import DEFAULTS, main, report, temperature_grid, validate


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

    def test_kinetic_energy_does_not_enter_the_expansion_covariance(self):
        # K fluctuates but is independent of V: alpha must not move, Cp must include it.
        a = np.tile([-1., -1., 1., 1.], 100)
        c = np.tile([1., -1., -1., 1.], 100)  # orthogonal to a, so Cov(V, K) is exactly zero in total
        rows = [dict(volume_nm3=10.+x*.01, potential_kJ_mol=1000.+x*5., kinetic_kJ_mol=0.) for x in a]
        noisy = [dict(r, kinetic_kJ_mol=500.+k*50.) for r, k in zip(rows, c)]
        rng = np.random.default_rng(1)
        sample = rng.choice(len(rows), 40, replace=False)
        base, with_k = thermodynamics(rows, 10000, 300, 1.), thermodynamics(noisy, 10000, 300, 1.)
        self.assertAlmostEqual(with_k['alpha_volume_1_K'], base['alpha_volume_1_K'])
        self.assertGreater(with_k['Cp_J_kg_K'], base['Cp_J_kg_K'])
        # On a subsample where V and K happen to correlate, alpha still ignores K.
        sub = [rows[i] for i in sample]
        sub_k = [noisy[i] for i in sample]
        self.assertAlmostEqual(thermodynamics(sub_k, 10000, 300, 1.)['alpha_volume_1_K'],
                               thermodynamics(sub, 10000, 300, 1.)['alpha_volume_1_K'])

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


class TgFigureTests(unittest.TestCase):
    def test_tg_figure_with_heating_history_and_without_a_fit(self):
        import csv
        from workflows.thermomechanical.analysis import fit_tg, smooth_hyperbola
        from workflows.network_properties.plots import plot_tg
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            rows = [dict(target_temperature_K=t, volume_nm3=40/d, density_g_cm3=d)
                    for t in np.arange(800, 99, -20.) for d in [smooth_hyperbola(t, 1.15, 420., -2.5e-4, -4e-4, 20.)]]
            with (base/'cooling.csv').open('w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            fit = fit_tg([r['target_temperature_K'] for r in rows], [r['density_g_cm3'] for r in rows])
            history = dict(phases=[dict(phase='01_heat_000', temperature_K=400., final_density_g_cm3=1.1)])
            path = plot_tg(base, dict(cooling=fit, thermal_history_density=history, smoke=False))
            self.assertEqual(path, base/'tg_density_volume.png')
            self.assertGreater(path.stat().st_size, 10000)
            self.assertTrue(plot_tg(base, dict(cooling=dict(reason='too few bins'), smoke=True),
                                    path=base/'none.png').is_file())


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


class NetworkGraphTests(unittest.TestCase):
    def lattice(self):
        """8 six-fold junctions on a periodic 2x2x2 grid (box 6) joined by 24 two-atom strands,
        with a 9-atom loop, a dangling chain, a dangling ring and one sol molecule.
        The grid closes 6-membered cycles only through the boundary: they are not rings."""
        atoms, bonds, xyz = [], [], []
        def atom(position, z=6, aromatic=False):
            atoms.append(dict(atomic_number=z, mass=12., aromatic=aromatic))
            xyz.append(np.asarray(position, float) % 6)
            return len(atoms) - 1
        sites = {(i, j, k): atom((3*i, 3*j, 3*k)) for i in range(2) for j in range(2) for k in range(2)}
        for (i, j, k), a in sites.items():
            for d in np.eye(3, dtype=int):
                first, second = (atom(np.array((3*i, 3*j, 3*k)) + n * d) for n in (1, 2))
                other = sites[tuple((np.array((i, j, k)) + d) % 2)]
                bonds += [[a, first, 1], [first, second, 1], [second, other, 1]]
        hub = sites[0, 0, 0]
        loop = [atom((.1*m, .2, .3)) for m in range(1, 10)]
        bonds += [[hub, loop[0], 1], [loop[-1], hub, 1]] + [[loop[m], loop[m+1], 1] for m in range(8)]
        chain = [atom((3.2, 3.2, 3.4)), atom((3.2, 3.2, 3.6))]
        bonds += [[sites[1, 1, 1], chain[0], 1], [chain[0], chain[1], 1]]
        linker = atom((0.3, 3.3, 3.3))
        ring = [atom((0.5 + .2*np.cos(m), 3.5 + .2*np.sin(m), 3.3), aromatic=True) for m in range(6)]
        bonds += [[sites[0, 1, 1], linker, 1], [linker, ring[0], 1]] + [[ring[m], ring[(m+1) % 6], 4] for m in range(6)]
        network = dict(atoms=atoms, bonds=bonds)
        sol = dict(atoms=[dict(atomic_number=6, mass=12., aromatic=False)] * 2, bonds=[[0, 1, 1]])
        return [network, sol], np.vstack([xyz, [[4.5, 4.5, 4.5], [4.6, 4.5, 4.5]]]), np.eye(3) * 6

    def test_strands_junctions_loops_and_dangling_material(self):
        from workflows.network_properties.network_graph import network_structure
        graphs, xyz, box = self.lattice()
        r = network_structure(graphs, xyz, box, temperature_K=300.)
        total = 76
        self.assertEqual(r['gel_molecules'], 1)
        self.assertEqual(r['junction_functionality'], {6: 8})
        self.assertEqual((r['effective_strands'], r['elastically_ineffective_loops']), (24, 1))
        self.assertAlmostEqual(r['sol_mass_fraction'], 2 / total)
        self.assertAlmostEqual(r['dangling_mass_fraction'], 9 / total)
        self.assertAlmostEqual(r['loop_mass_fraction'], 9 / total)
        self.assertAlmostEqual(r['active_mass_fraction'], 56 / total)
        volume_m3 = 216e-27
        self.assertAlmostEqual(r['strand_density_nu_e_mol_m3'], 24 / NA / volume_m3)
        self.assertAlmostEqual(r['cycle_rank_density_mmol_cm3'], (24 - 8) / NA / volume_m3 / 1000)
        self.assertAlmostEqual(r['Mc_g_mol'], total * 12. / 24)
        self.assertAlmostEqual(r['Mc_active_network_g_mol'], 56 * 12. / 24)
        self.assertAlmostEqual(r['rubbery_shear_modulus_estimates_MPa']['affine'],
                               24 / volume_m3 * KB * 300 / 1e6)
        self.assertEqual(r['polar_groups_count']['aromatic_rings'], 1)
        # 24 strands of two interior atoms; the loop is not a strand.
        self.assertEqual(r['strand_heavy_atoms_histogram'], {2: 24})
        self.assertEqual(r['strand_backbone_bonds_histogram'], {3: 24})
        self.assertEqual(r['strand_molecular_weight_histogram'], {0: 24})  # 2 x 12 g/mol
        self.assertEqual(r['strand_table']['rows'], [[2, 24.0, 3]] * 24)

    def with_loop(self, graphs, xyz, at):
        """Copy of the lattice with a 9-atom loop closing on atom at of the network."""
        network = dict(atoms=list(graphs[0]['atoms']), bonds=list(graphs[0]['bonds']))
        first = len(network['atoms'])
        network['atoms'] += [dict(atomic_number=6, mass=12., aromatic=False)] * 9
        network['bonds'] += [[at, first, 1], [first + 8, at, 1]] + [[first + m, first + m + 1, 1] for m in range(8)]
        added = [(xyz[at] + (.05 * (m + 1), .05, .05)) % 6 for m in range(9)]
        return [network] + graphs[1:], np.vstack([xyz[:first], added, xyz[first:]])

    def test_loops_do_not_change_strands_or_junctions(self):
        from workflows.network_properties.network_graph import network_structure
        graphs, xyz, box = self.lattice()
        base = network_structure(graphs, xyz, box)
        # Atom 8 is inside a strand; atom 60 is inside the loop already on the lattice,
        # which closes a second loop only after the first one is removed.
        for at, loops in ((8, 2), (60, 2)):
            with self.subTest(at=at):
                r = network_structure(*self.with_loop(graphs, xyz, at), box)
                self.assertEqual((r['effective_strands'], r['junctions'], r['elastically_ineffective_loops']),
                                 (24, 8, loops))
                self.assertEqual(r['junction_functionality'], {6: 8})
                self.assertEqual(r['strand_table']['rows'], base['strand_table']['rows'])
                self.assertEqual(r['rubbery_shear_modulus_estimates_MPa'], base['rubbery_shear_modulus_estimates_MPa'])

    def test_mc_definitions_differ_by_which_mass_they_count(self):
        from workflows.network_properties.network_graph import network_structure
        graphs, xyz, box = self.lattice()
        sol = [dict(atoms=[dict(atomic_number=6, mass=12., aromatic=False)], bonds=[]) for _ in range(1000)]
        r = network_structure(graphs + sol, np.vstack([xyz, np.tile([[4.5, 1.5, 1.5]], (1000, 1))]), box)
        self.assertAlmostEqual(r['Mc_g_mol'], 538.)
        self.assertAlmostEqual(r['Mc_active_network_g_mol'], 28.)
        self.assertAlmostEqual(r['strand_molecular_weight_g_mol']['mean'], 24.)

    def test_finite_molecules_are_sol(self):
        from workflows.network_properties.network_graph import network_structure
        graphs, xyz, box = self.lattice()
        r = network_structure(graphs, xyz, box * 100, temperature_K=300.)  # nothing wraps in a large cell
        self.assertEqual((r['gel_molecules'], r['effective_strands'], r['junctions']), (0, 0, 0))
        self.assertAlmostEqual(r['sol_mass_fraction'], 1.)
        self.assertIsNone(r['Mc_g_mol'])

    def test_polar_groups(self):
        from workflows.network_properties.network_graph import polar_groups
        def molecule(heavy, hydrogens, bonds):
            atoms = [dict(atomic_number=z, aromatic=False) for z in heavy] + [dict(atomic_number=1, aromatic=False)] * len(hydrogens)
            return dict(atoms=atoms, bonds=bonds + [[parent, len(heavy) + n, 1] for n, parent in enumerate(hydrogens)])
        methanol = molecule([6, 8], [0, 0, 0, 1], [[0, 1, 1]])
        ether = molecule([6, 8, 6], [0, 0, 0, 2, 2, 2], [[0, 1, 1], [1, 2, 1]])
        oxirane = molecule([6, 6, 8], [0, 0, 1, 1], [[0, 1, 1], [1, 2, 1], [2, 0, 1]])
        acetone = molecule([6, 6, 8, 6], [0, 0, 0, 3, 3, 3], [[0, 1, 1], [1, 2, 2], [1, 3, 1]])
        counts = polar_groups([methanol, ether, oxirane, acetone])
        self.assertEqual(counts, dict(aliphatic_OH=1, ether_O=1, epoxide_rings=1, carbonyl_O=1,
                                      hbond_acceptor_atoms=4, hbond_donor_H=1))

    def test_hbond_acceptors_depend_on_the_functional_group(self):
        from workflows.network_properties.network_graph import polar_groups
        def acceptors(elements, hydrogens, bonds, aromatic=(), charges=None):
            atoms = [dict(atomic_number=z, aromatic=i in aromatic, formal_charge=(charges or {}).get(i, 0))
                     for i, z in enumerate(elements)]
            atoms += [dict(atomic_number=1, aromatic=False, formal_charge=0)] * len(hydrogens)
            graph = dict(atoms=atoms, bonds=bonds + [[p, len(elements) + n, 1] for n, p in enumerate(hydrogens)])
            return polar_groups([graph]).get('hbond_acceptor_atoms', 0)
        ring = [[m, (m + 1) % 5, 4] for m in range(5)]
        self.assertEqual(acceptors([6, 6, 7, 8], [0, 0, 0, 2, 2], [[0, 1, 1], [1, 2, 1], [1, 3, 2]]), 1)  # acetamide
        self.assertEqual(acceptors([6, 6, 8, 8], [0, 0, 0, 3], [[0, 1, 1], [1, 2, 2], [1, 3, 1]]), 1)  # acetic acid
        self.assertEqual(acceptors([6, 7, 8, 8], [0, 0, 0], [[0, 1, 1], [1, 2, 2], [1, 3, 1]], charges={1: 1, 3: -1}), 2)  # nitro
        self.assertEqual(acceptors([6, 7, 6, 6], [], [[0, 1, 1], [1, 2, 1], [1, 3, 1]]), 1)  # tertiary amine
        self.assertEqual(acceptors([7, 6, 6, 6, 6], [0], ring, aromatic=range(5)), 0)  # pyrrole


class ElasticReuseTests(unittest.TestCase):
    CFG = dict(platform='CPU', smoke=False, sample_interval_ps=1.)

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name)
        (self.base/'relaxed_300K.structure.json').write_text('{}')
        self.directory = self.base/'mechanical'
        self.directory.mkdir()

    def branches(self):
        from workflows.network_properties.mechanical import VOIGT
        for name in VOIGT:
            for sign in ('plus', 'minus'):
                (self.directory/f'{name}_{sign}.json').write_text(json.dumps(
                    [dict(time_ps=i, **dict.fromkeys(VOIGT, .1 * i)) for i in range(5)]))

    def run_with_md(self):
        """run_mechanical with MD replaced; True when it had to start MD."""
        from workflows.network_properties import mechanical as m
        class Started(Exception):
            pass
        with patch('polypaves.sorption.load_network', return_value=object()), \
             patch.object(m, 'context_for', side_effect=Started), contextlib.redirect_stdout(io.StringIO()):
            try:
                m.run_mechanical(self.base, self.CFG, 1)
            except Started:
                return True
        return False

    def test_results_without_provenance_are_retired_and_recomputed(self):
        (self.directory/'result.json').write_text(json.dumps(dict(old=True)))
        self.branches()
        self.assertTrue(self.run_with_md())
        retired = self.base/'mechanical.stale-1'
        self.assertEqual(json.loads((retired/'result.json').read_text()), dict(old=True))
        self.assertTrue((retired/'xx_plus.json').is_file())
        self.assertFalse((self.directory/'xx_plus.json').exists())

    def test_branches_of_the_same_method_and_input_are_reused(self):
        self.assertTrue(self.run_with_md())
        self.branches()
        self.assertFalse(self.run_with_md())
        from workflows.network_properties.mechanical import result_current
        self.assertTrue(result_current(self.base))
        (self.base/'relaxed_300K.structure.json').write_text('{"changed": 1}')
        self.assertFalse(result_current(self.base))
        self.assertTrue(self.run_with_md())
        self.assertTrue((self.base/'mechanical.stale-1/result.json').is_file())


class ConfigurationTests(unittest.TestCase):
    def test_long_relaxation_stages_json_is_not_taken_for_a_file_name(self):
        from workflows.crosslinking.crosslinking import load_relaxation_stages
        stages = [dict(ensemble='NVT', steps=10000 + i, timestep=.1) for i in range(10)]
        self.assertEqual(load_relaxation_stages(json.dumps(stages)), stages)

    def test_network_input_refuses_sol_removal(self):
        with tempfile.TemporaryDirectory() as temporary, contextlib.redirect_stderr(io.StringIO()) as error:
            source = Path(temporary)/'source.json'
            source.write_text(json.dumps(dict(graphs=[dict(atoms=[{'partial_charge':0.}])], positions=[[0,0,0]],
                                              box_vectors=[[3,0,0],[0,3,0],[0,0,3]])))
            self.assertEqual(main(['new', str(Path(temporary)/'p'), '--network', str(source),
                                   '--remove-sol-fraction']), 1)
        self.assertIn('--network cannot be combined', error.getvalue())

    def test_report_and_collect_refuse_a_configuration_edited_after_starting(self):
        from workflows.network_properties.execution import configuration_fingerprint
        from workflows.network_properties.directory_summary import collect_projects
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary).resolve()/'p'
            project.mkdir()
            cfg = dict(DEFAULTS, platform='CPU', replicates=1, seed=100, network='input.structure.json',
                       mechanical=False)
            (project/'project.json').write_text(json.dumps(dict(cfg, seed=999)))
            (project/'run_signature.json').write_text(json.dumps(dict(sha256=configuration_fingerprint(cfg))))
            with self.assertRaisesRegex(ValueError, 'changed after starting'):
                report(project)
            with contextlib.redirect_stdout(io.StringIO()):
                collected = collect_projects(project.parent)
            self.assertEqual(collected['projects'], [])
            self.assertIn('changed after starting', collected['errors'][0]['error'])

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
