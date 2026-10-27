"""Installed-wheel integration: crosslink -> saved charges -> solvent -> MD/resume."""
import importlib.util
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


class DistributionPolicyTests(unittest.TestCase):
    def test_nested_implementation_source_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            wheel = Path(temp) / 'source-leak.whl'
            with zipfile.ZipFile(wheel, 'w') as archive:
                archive.writestr('polypaves/reaction/engine.py', 'raise RuntimeError()\n')
            result = subprocess.run([str(ROOT / 'scripts/audit_binary_wheel.sh'), str(wheel)],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('source or development artifact', result.stderr)


class CrosslinkingDefaultsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('crosslinking_defaults',
            ROOT / 'workflows/crosslinking/crosslinking.py')
        cls.workflow = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.workflow)

    def test_cli_and_python_defaults_and_overrides(self):
        import inspect
        defaults = dict(total_atoms=20000, atom_count_mode='above',
                        cycles=5000, target_conversion=.995, save_trajectory=False)
        signature = inspect.signature(self.workflow.run)
        for name, value in defaults.items():
            self.assertEqual(signature.parameters[name].default, value)
        for flags, expected in [([], defaults),
                (['--total-atoms', '500', '--atom-count-mode', 'nearest', '--cycles', '1',
                  '--target-conversion', '.9', '--save-trajectory'],
                 dict(total_atoms=500, atom_count_mode='nearest', cycles=1, target_conversion=.9,
                      save_trajectory=True))]:
            with patch.object(sys, 'argv', ['crosslinking.py', '--project', 'check', *flags]), \
                 patch.object(self.workflow, 'run') as run, \
                 patch.object(self.workflow.p, 'configure_logging'):
                self.workflow.main()
                for name, value in expected.items():
                    self.assertEqual(run.call_args.kwargs[name], value)

    def test_default_epoxy_exceeds_twenty_thousand_including_activation(self):
        plan = self.workflow.plan_size(20000, 47, 3, added_per_site=2, fraction=.1)
        self.assertEqual(plan['monomers'], 421)
        self.assertEqual(plan['activation_count'], 126)
        self.assertEqual(plan['actual_total_atoms'], 20039)

    def test_above_is_strict_and_minimal_for_count_and_fraction(self):
        import math
        for target in (500, 4704, 19992, 20000, 20039):
            for selection in ({'count': 2}, {'fraction': .1}):
                with self.subTest(target=target, selection=selection):
                    plan = self.workflow.plan_size(target, 47, 3, added_per_site=2, **selection)
                    self.assertGreater(plan['actual_total_atoms'], target)
                    previous = plan['monomers'] - 1
                    starts = selection.get('count', math.floor(.1 * 3 * previous))
                    self.assertLessEqual(47 * previous + 2 * starts, target)

    def test_explicit_legacy_rounding_modes(self):
        for mode in ('nearest', 'floor', 'exact'):
            plan = self.workflow.plan_size(19992, 47, 3, added_per_site=2, fraction=.1, mode=mode)
            self.assertEqual(plan['actual_total_atoms'], 19992)
        with self.assertRaises(ValueError):
            self.workflow.plan_size(20000, 47, 3, added_per_site=2, fraction=.1, mode='exact')
        with self.assertRaises(ValueError):
            self.workflow.plan_size(20000, 47, 3, added_per_site=2, fraction=.1, mode='invalid')


class PublicReactiveWorkflows(unittest.TestCase):
    def test_installed_implementation_is_compiled(self):
        for module in ('polypaves.api', 'polypaves.reaction.rule', 'polypaves.reaction.local_sage',
                       'polypaves.sorption', 'polypaves.sorption_workflow'):
            spec = importlib.util.find_spec(module)
            self.assertIsNotNone(spec, module)
            self.assertTrue(spec.origin.endswith('.so'), spec.origin)

    def test_identical_molecules_retain_instance_charges(self):
        import numpy as np
        import polypaves as p
        from openmm import unit
        from polypaves.sorption import load_network
        template = p.make_solvent(p.Solvent('CO'))
        graphs = [copy.deepcopy(template.graphs[0]) for _ in range(2)]
        for mi, graph in enumerate(graphs):
            for ai, atom in enumerate(graph['atoms']):
                atom['stable_atom_id'] = f'host-{mi}-{ai}'
        graphs[1]['atoms'][0]['partial_charge'] += .01
        graphs[1]['atoms'][1]['partial_charge'] -= .01
        positions = np.concatenate((template.positions_nm + [2., 2., 2.],
                                    template.positions_nm + [3., 2., 2.])) * 10.
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'final.structure.json'
            path.write_text(json.dumps(dict(graphs=graphs, positions=positions.tolist(),
                                           box_vectors=(np.eye(3) * 60.).tolist())))
            with patch('openff.toolkit.Molecule.assign_partial_charges',
                       side_effect=AssertionError('saved charges must be retained')):
                host = load_network(path)
            self.assertEqual(host.graphs, graphs)
            nb = next(f for f in host.system.getForces() if f.__class__.__name__ == 'NonbondedForce')
            np.testing.assert_array_equal([a['partial_charge'] for g in graphs for a in g['atoms']],
                [nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                 for i in range(host.atom_count)])

    def test_kekule_aromatic_epoxy_is_packed_as_aromatic(self):
        # A Kekulé SMILES used to stop in reaction-template precompute with
        # "non-ring atom N marked aromatic": rings stored as 1/2 bonds were cut.
        kekule = ('C12C=CC=CC=1C1C=CC=CC=1C2(C1C=CC2C=C(OCC3OC3)C=CC=2C=1)'
                  'C1C=C2C=CC(OCC3CO3)=CC2=CC=1')
        env = dict(os.environ, OPENMM_CPU_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        env.pop('PYTHONPATH', None)
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp)
            for name in ('activation.yaml', 'reactions.yaml'):
                shutil.copy2(ROOT / 'workflows/crosslinking' / name, project / name)
            completed = subprocess.run([sys.executable, str(ROOT / 'workflows/crosslinking/crosslinking.py'),
                '--project', 'kekule', '--smiles', kekule, '--total-atoms', '400', '--activation-count', '1',
                '--platform', 'CPU', '--compression-steps', '2', '--equilibration-steps', '2', '--cycles', '1',
                '--nvt-steps', '2', '--npt-steps', '2', '--max-reactions-per-cycle', '1',
                '--reaction-cutoff', '10', '--max-reaction-cutoff', '10', '--minimization-iterations', '100'],
                cwd=project, env=env, capture_output=True, text=True, timeout=600)
            self.assertEqual(completed.returncode, 0, completed.stdout[-3000:] + completed.stderr[-3000:])
            output = project / 'results/kekule'
            settings = json.loads((output / 'settings.json').read_text())
            self.assertEqual(settings['input_smiles'], kekule)
            self.assertNotIn('=', settings['smiles'])
            graph = json.loads((output / 'initial.structure.json').read_text())['graphs'][0]
            self.assertEqual(sum(a['aromatic'] for a in graph['atoms']), 32)
            self.assertGreaterEqual(json.loads((output / 'metrics.json').read_text())['number_of_reactions'], 1)

    def test_configured_post_reaction_relaxation_runs_each_stage(self):
        env = dict(os.environ, OPENMM_CPU_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        env.pop('PYTHONPATH', None)
        stages = [dict(ensemble='NVT', steps=2, timestep=.1), dict(ensemble='NVT', steps=2, timestep=.5),
                  dict(ensemble='NPT', steps=2, timestep=1., pressure_atm=1.)]
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp)
            for name in ('activation.yaml', 'reactions.yaml'):
                shutil.copy2(ROOT / 'workflows/crosslinking' / name, project / name)
            completed = subprocess.run([sys.executable, str(ROOT / 'workflows/crosslinking/crosslinking.py'),
                '--project', 'relax', '--total-atoms', '500', '--activation-count', '1',
                '--platform', 'CPU', '--compression-steps', '2', '--equilibration-steps', '2', '--cycles', '1',
                '--max-reactions-per-cycle', '1', '--reaction-cutoff', '10', '--max-reaction-cutoff', '10',
                '--minimization-iterations', '100', '--post-reaction-md-steps', '2',
                '--relaxation-stages', json.dumps(stages)],
                cwd=project, env=env, capture_output=True, text=True, timeout=600)
            self.assertEqual(completed.returncode, 0, completed.stdout[-3000:] + completed.stderr[-3000:])
            output = project / 'results/relax'
            settings = json.loads((output / 'settings.json').read_text())
            self.assertEqual(settings['relaxation_stages'], stages)
            self.assertEqual(settings['post_reaction_md_steps'], 2)
            self.assertEqual([s['timestep'] for s in settings['cycle_md_stages']], [.1, .5, 1.])
            cycle, = json.loads((output / 'cycles.json').read_text())
            self.assertGreaterEqual(cycle['reactions'], 1)
            self.assertEqual(len(cycle['md_stages']), 3)

    def test_crosslink_water_pack_and_md_resume(self):
        import numpy as np
        env = dict(os.environ, OPENMM_CPU_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        env.pop('PYTHONPATH', None)

        def run(arguments, cwd):
            completed = subprocess.run([sys.executable, *map(str, arguments)], cwd=cwd, env=env,
                capture_output=True, text=True, timeout=300)
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            return completed

        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp)
            for name in ('activation.yaml', 'reactions.yaml'):
                shutil.copy2(ROOT / 'workflows/crosslinking' / name, project / name)
            run([ROOT / 'workflows/crosslinking/crosslinking.py', '--project', 'epoxy',
                '--total-atoms', '500', '--activation-count', '1', '--platform', 'CPU',
                '--compression-steps', '2', '--equilibration-steps', '2', '--cycles', '1', '--nvt-steps', '2', '--npt-steps', '2',
                '--max-reactions-per-cycle', '1', '--reaction-cutoff', '10', '--max-reaction-cutoff', '10',
                '--minimization-iterations', '100'], project)
            network = project / 'results/epoxy/final.structure.json'
            dry = json.loads(network.read_text())
            metrics = json.loads((network.parent / 'metrics.json').read_text())
            self.assertGreaterEqual(metrics['number_of_reactions'], 1)
            self.assertEqual(metrics['frames'], 0)
            self.assertFalse((network.parent / 'cycles').exists())
            self.assertFalse((network.parent / 'frames').exists())
            self.assertFalse((network.parent / 'atoms.lammpstrj.gz').exists())
            self.assertFalse((network.parent / 'bonds.dump.gz').exists())
            self.assertFalse(json.loads((network.parent / 'settings.json').read_text())['save_trajectory'])
            self.assertIn('USER_CHARGES', (network.parent / 'final.mol2').read_text())
            before = network.read_bytes()
            prepared = project / 'packed'
            run([ROOT / 'workflows/sorption/sorption.py', 'pack', '--network', network,
                '--smiles', 'O', '--molecules', '3', '--cell-scale', '1.1', '--output', prepared], project)
            wet = json.loads((prepared / 'packed.structure.json').read_text())
            self.assertEqual(network.read_bytes(), before)
            self.assertEqual(wet['graphs'][:len(dry['graphs'])], dry['graphs'])
            np.testing.assert_allclose(np.array(wet['positions'][:len(dry['positions'])]),
                                       np.array(dry['positions']) * 1.1, atol=1e-12, rtol=0)
            self.assertEqual(wet['solvent_molecules'], 3)
            self.assertEqual(len(wet['positions']), len(dry['positions']) + 9)
            md = project / 'md'
            command = [ROOT / 'workflows/sorption/sorption.py', 'simulate', '--packed', prepared,
                '--output', md, '--platform', 'CPU', '--nvt-ps', '.01', '--npt-equilibration-ps', '.01',
                '--production-ps', '.02', '--sample-interval-ps', '.01', '--frame-interval-ps', '.01',
                '--minimization-iterations', '100']
            run(command, project)
            self.assertTrue((md / 'checkpoint.json').is_file())
            self.assertTrue(list((md / 'frames').glob('*.data.gz')))
            final = json.loads((md / 'final.structure.json').read_text())
            self.assertEqual(final['graphs'][:len(dry['graphs'])], dry['graphs'])
            self.assertEqual(final['solvent_molecules'], 3)
            self.assertTrue(np.isfinite(final['positions']).all())
            original = (md / 'result.json').read_bytes()
            run(command + ['--resume'], project)
            self.assertEqual((md / 'result.json').read_bytes(), original)

            # Opt-in retains visualization output as well as the final network.
            run([ROOT / 'workflows/crosslinking/crosslinking.py', '--project', 'with_trajectory',
                '--total-atoms', '500', '--activation-count', '1', '--platform', 'CPU',
                '--compression-steps', '2', '--equilibration-steps', '2', '--cycles', '1', '--nvt-steps', '2', '--npt-steps', '2',
                '--max-reactions-per-cycle', '1', '--reaction-cutoff', '10', '--max-reaction-cutoff', '10',
                '--minimization-iterations', '100', '--save-trajectory',
                '--template-cache', network.parent / 'reaction_templates'], project)
            visual = project / 'results/with_trajectory'
            self.assertTrue((visual / 'atoms.lammpstrj.gz').is_file())
            self.assertTrue(list((visual / 'frames').glob('*.data.gz')))
            self.assertTrue(list((visual / 'cycles').glob('cycle_*.structure.json')))
            self.assertEqual(json.loads((visual / 'final.structure.json').read_text())['graphs'], dry['graphs'])
            self.assertEqual(json.loads((visual / 'metrics.json').read_text())['frames'], 3)


if __name__ == '__main__':
    unittest.main()
