import importlib.util
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('pre_reaction', ROOT / 'workflows/crosslinking/pre_reaction.py')
pre = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pre)


class PreReactionEpoxides(unittest.TestCase):
    def test_groups_pick_ring_oxygen_and_the_ch2_carbon(self):
        # Epoxide C1(H2)-O2-C3(H): only C1 carries two H and is the attacked carbon.
        atoms = {'c1': 6, 'o2': 8, 'c3': 6, 'h1': 1, 'h2': 1, 'h3': 1, 'x': 6}
        bonds = [('c1', 'o2'), ('o2', 'c3'), ('c1', 'c3'), ('c1', 'h1'), ('c1', 'h2'), ('c3', 'h3'), ('c3', 'x')]
        neighbors = {sid: set() for sid in atoms}
        for a, b in bonds:
            neighbors[a].add(b); neighbors[b].add(a)
        graph = SimpleNamespace(atoms={sid: dict(atomic_number=z) for sid, z in atoms.items()}, neighbors=neighbors)
        order = ('x', 'c1', 'o2', 'c3', 'h1', 'h2', 'h3')
        self.assertEqual(pre.epoxide_groups(graph, {('c1', 'c3', 'o2')}, order),
                         [dict(oxygen=2, carbons=[1], heavy=[1, 3, 2])])

    def test_contacts_use_the_minimum_image_and_both_directions(self):
        box = np.diag([20.] * 3)
        # Group 0 and 1 touch across the x boundary (O0 to C1 is 3 A); group 2 is far away.
        oxygen = np.array([[1., 5, 5], [18., 9, 5], [10., 15, 15]])
        carbons = np.array([[[1., 9, 5]], [[18., 5, 5]], [[10., 18, 15]]])
        near = pre.contacts(oxygen, carbons, box, 4.)
        self.assertEqual(near.tolist(), [[False, True, False], [True, False, False], [False, False, False]])
        stats = pre.cluster_statistics(near)
        self.assertAlmostEqual(stats['coordination_number'], 2 / 3)
        self.assertAlmostEqual(stats['largest_cluster_fraction'], 2 / 3)
        self.assertEqual(stats['clusters'], 2)

    def test_msd_unwraps_across_the_boundary_and_removes_drift(self):
        box = np.diag([10.] * 3)
        groups = [dict(oxygen=0, carbons=[1], heavy=[0, 1])]
        # The ring moves +1 A in x per ps and crosses x = 10; a very heavy third atom pins the
        # centre of mass, so the drift correction is negligible.
        sampler = pre.EpoxideSampler(groups, [1., 1., 1e9], [4.])
        for frame in range(6):
            x = (8.5 + frame) % 10.
            sampler.add(np.array([[x, 1, 1], [(x + 1) % 10, 1, 1], [5., 5, 5]]), box, float(frame))
        lag, msd = pre.mean_squared_displacement(np.array(sampler.centroids), sampler.times)
        np.testing.assert_allclose(lag, [1, 2, 3, 4, 5])
        np.testing.assert_allclose(msd, np.arange(1, 6) ** 2, rtol=1e-6)
        # The whole cell translating together is centre-of-mass drift, not diffusion.
        drifting = pre.EpoxideSampler(groups, [1., 1., 1.], [4.])
        for frame in range(4):
            drifting.add(np.array([[1. + frame, 1, 1], [2. + frame, 1, 1], [5. + frame, 5, 5]]), box, float(frame))
        self.assertAlmostEqual(pre.mean_squared_displacement(np.array(drifting.centroids), drifting.times)[1].max(), 0.)
        fit = pre.diffusion_fit(np.linspace(1, 100, 100), 6 * .5 * np.linspace(1, 100, 100))
        self.assertAlmostEqual(fit['diffusion_coefficient_m2_s'], .5e-8)
        self.assertAlmostEqual(fit['msd_exponent_alpha'], 1.)
        with tempfile.TemporaryDirectory() as temp:
            result = sampler.summary(temp, temperature_K=363.15, pressure_atm=1., sampling_ps=5.)
            self.assertEqual((result['frames'], result['epoxides']), (6, 1))
            self.assertAlmostEqual(result['epoxide_number_density_nm3']['mean'], 1.)
            self.assertAlmostEqual(result['msd_final_A2'], 25., places=4)
            self.assertTrue((Path(temp) / 'epoxide_frames.csv').is_file())


chains_spec = importlib.util.spec_from_file_location('chains', ROOT / 'workflows/crosslinking/chains.py')
chains = importlib.util.module_from_spec(chains_spec)
chains_spec.loader.exec_module(chains)


class KineticChains(unittest.TestCase):
    def test_each_reaction_follows_the_activity_to_its_start(self):
        def move(old, new):
            return {'cycle': 0, 'atom_changes': {old: {'old': {'reaction_state': 'cationic_activity'},
                                                       'new': {'reaction_state': None}},
                                                 new: {'old': {'reaction_state': None},
                                                       'new': {'reaction_state': 'cationic_activity'}}}}
        activation = {'reaction': {'atom_changes': {s: {'old': {'reaction_state': None},
                                                        'new': {'reaction_state': 'cationic_activity'}}
                                                    for s in ('a0', 'b0', 'c0')}}}
        # a grows 3 times, b once, c never.
        history = [move('a0', 'a1'), move('b0', 'b1'), move('a1', 'a2'), move('a2', 'a3')]
        found = chains.kinetic_chains([activation], history)
        self.assertEqual([(c['start_atom'], c['epoxide_units'], c['active_atom']) for c in found],
                         [('a0', 4, 'a3'), ('b0', 2, 'b1'), ('c0', 1, 'c0')])
        stats = chains.chain_statistics(found)
        self.assertAlmostEqual(stats['number_average_epoxide_units'], 7 / 3)
        self.assertAlmostEqual(stats['weight_average_epoxide_units'], 21 / 7)
        self.assertEqual(stats['distribution'], [dict(epoxide_units=1, active_ends=1),
                                                 dict(epoxide_units=2, active_ends=1),
                                                 dict(epoxide_units=4, active_ends=1)])
        with self.assertRaises(ValueError):
            chains.kinetic_chains([activation], [move('zz', 'a1')])


if __name__ == '__main__':
    unittest.main()
