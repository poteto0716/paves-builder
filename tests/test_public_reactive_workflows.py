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
                '--equilibration-steps', '2', '--cycles', '1', '--nvt-steps', '2', '--npt-steps', '2',
                '--max-reactions-per-cycle', '1', '--reaction-cutoff', '10', '--max-reaction-cutoff', '10',
                '--minimization-iterations', '100'], project)
            network = project / 'results/epoxy/final.structure.json'
            dry = json.loads(network.read_text())
            metrics = json.loads((network.parent / 'metrics.json').read_text())
            self.assertGreaterEqual(metrics['number_of_reactions'], 1)
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


if __name__ == '__main__':
    unittest.main()
