"""Monomer and crosslink projects: defaults, seed/replica routing, stage order and restart provenance."""
from argparse import Namespace
from contextlib import redirect_stderr, redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from workflows.network_properties.network_properties import main, input_digest, run, validate
from workflows.network_properties.construction import (build_command, ensure_crosslink, ensure_monomer_build,
                                                       network_path)
from workflows.network_properties import monomer


def quiet(argv):
    with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
        return main(argv)


def complete_monomer_replica(project, cfg, replica, content='relaxed'):
    """Files a finished monomer replica leaves behind (without running MD)."""
    base = Path(project) / f'replica_{replica:02d}'
    for path, data in ((base / 'build/done.json', {}), (base / 'dry_result.json', {}),
                       (monomer.temperature_dir(project, cfg, replica) / 'result.json', {})):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data))
    monomer.relaxed_path(project, cfg, replica).write_text(json.dumps(dict(graphs=[content])))


class ProjectCreationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.monomer = self.root / 'TMPTGE'
        self.assertEqual(quiet(['monomer', str(self.monomer), '--replicates', '2', '--seed', '100',
                                '--platform', 'CPU', '--smoke']), 0)
        self.cfg = json.loads((self.monomer / 'project.json').read_text())

    def test_monomer_defaults_build_and_relax_at_T_without_reaction(self):
        cfg = self.cfg
        self.assertEqual((cfg['input_mode'], cfg['temperature_K'], cfg['chemistry']), ('monomer', 363.15, 'epoxy'))
        self.assertEqual(cfg['reactive_sites_per_molecule'], 3)
        build = cfg['build']
        self.assertEqual(build['total_atoms'], 500)  # smoke-sized
        for key, value in (('activation_count', 0), ('cycles', 0), ('target_conversion', None),
                           ('reference_npt_steps', 0), ('pre_reaction_sampling_steps', 0),
                           ('final_relaxation_steps', 0), ('final_md_steps', 0), ('temperature', 363.15)):
            self.assertEqual(build[key], value, key)
        self.assertEqual(cfg['solvents'], [])
        self.assertEqual(input_digest(self.monomer, cfg), cfg['input_sha256'])
        self.assertTrue((self.monomer / 'definitions/activation.yaml').is_file())
        full = self.root / 'full'
        self.assertEqual(quiet(['monomer', str(full), '--temperature', '353.15']), 0)
        cfg = json.loads((full / 'project.json').read_text())
        self.assertEqual((cfg['high_K'], cfg['low_K'], cfg['replicates'], cfg['platform']), (800., 100., 4, 'CUDA'))
        self.assertEqual((cfg['t_equilibration_ps'], cfg['t_production_ps'], cfg['t_nvt_production_ps']),
                         (1000., 2000., 2000.))
        self.assertEqual((cfg['production_ps'], cfg['nvt_production_ps']), (10000., 10000.))
        self.assertEqual((cfg['temperature_K'], cfg['build']['temperature']), (353.15, 353.15))
        self.assertEqual((cfg['build']['total_atoms'], cfg['build']['equilibration_steps']), (20000, 1_000_000))
        self.assertEqual(cfg['ring_contact_cutoffs_A'], [4., 6.])
        self.assertEqual(json.loads((self.root / '.network_properties_collection.json').read_text())['kind'],
                         'network_properties')

    def test_crosslink_projects_inherit_the_monomer_and_name_each_density(self):
        self.assertEqual(quiet(['crosslink', str(self.monomer), '--activation-fraction', '0.05', '0.1']), 0)
        names = sorted(p.name for p in (self.monomer / 'crosslinks').iterdir())
        self.assertEqual(names, ['af0.05', 'af0.1'])
        child = self.monomer / 'crosslinks/af0.05'
        cfg = json.loads((child / 'project.json').read_text())
        self.assertEqual((cfg['input_mode'], cfg['parent'], cfg['temperature_K']), ('crosslink', '../..', 363.15))
        self.assertEqual((cfg['seed'], cfg['replicates'], cfg['platform'], cfg['smoke']), (100, 2, 'CPU', True))
        options = cfg['crosslinking']
        self.assertEqual((options['activation_fraction'], options['activation_count'], options['temperature']),
                         (.05, None, 363.15))
        self.assertFalse({'smiles', 'total_atoms', 'density', 'compression_steps'} & set(options))
        self.assertIn('fraction: 0.05', (child / 'definitions/activation.yaml').read_text())
        self.assertEqual(input_digest(child, cfg), cfg['input_sha256'])
        validate(cfg)
        self.assertEqual(quiet(['crosslink', str(self.monomer), '--activation-count', '0']), 0)
        zero = json.loads((self.monomer / 'crosslinks/n0/project.json').read_text())['crosslinking']
        self.assertEqual((zero['cycles'], zero['target_conversion']), (0, None))
        for bad in (['--activation-fraction', '0.05'],                       # exists already
                    ['--activation-fraction', '0.3', '--replicates', '3'],   # more than the monomer replicas
                    ['--activation-fraction', '0.3', '0.4', '--name', 'x']):
            self.assertEqual(quiet(['crosslink', str(self.monomer), *bad]), 1, bad)
        self.assertEqual(quiet(['crosslink', str(child), '--activation-fraction', '0.3']), 1)
        full = self.root / 'full'
        self.assertEqual(quiet(['monomer', str(full)]), 0)
        self.assertEqual(quiet(['crosslink', str(full), '--activation-fraction', '0.1', '--remove-sol-fraction',
                                '--solvent', 'water=O']), 0)
        cfg = json.loads((full / 'crosslinks/af0.1/project.json').read_text())
        options = cfg['crosslinking']
        self.assertEqual((options['equilibration_steps'], options['cycles'], options['target_conversion'],
                          options['final_md_steps']), (100_000, 5000, .995, 1_000_000))
        self.assertTrue(options['remove_sol_fraction'])
        self.assertEqual(cfg['solvents'], [dict(name='water', smiles='O')])
        command = build_command(self.root, 'network', options, seed=1, platform='CUDA', start_structure='s.json')
        self.assertIn('--remove-sol-fraction', command)
        self.assertEqual(command[command.index('--start-structure') + 1], 's.json')

    def test_oxetane_preset(self):
        oxetane = self.root / 'oxetane'
        self.assertEqual(quiet(['monomer', str(oxetane), '--platform', 'CPU', '--chemistry', 'oxetane']), 0)
        cfg = json.loads((oxetane / 'project.json').read_text())
        self.assertEqual((cfg['chemistry'], cfg['build']['smiles']), ('oxetane', 'CCC1(COCC2(CC)COC2)COC1'))
        self.assertEqual(cfg['reactive_sites_per_molecule'], 2)
        self.assertIn('oxetane_OH_relay', (oxetane / 'definitions/reactions.yaml').read_text())
        self.assertEqual(monomer.ring_label(cfg), 'Oxetane')
        plain = self.root / 'plain'
        self.assertEqual(quiet(['monomer', str(plain), '--platform', 'CPU', '--smiles', 'CCCCCCCCO']), 0)
        self.assertEqual(json.loads((plain / 'project.json').read_text())['reactive_sites_per_molecule'], 0)
        # Without reactive rings only 0 starts make sense.
        self.assertEqual(quiet(['crosslink', str(plain), '--activation-count', '3']), 1)
        self.assertEqual(quiet(['crosslink', str(plain), '--activation-count', '0']), 0)

    def test_chemistry_follows_the_ring_of_the_smiles(self):
        # An oxetane SMILES without --chemistry gets the oxetane preset (sites, YAMLs, labels).
        guessed = self.root / 'guessed'
        self.assertEqual(quiet(['monomer', str(guessed), '--platform', 'CPU', '--smiles', 'CCC1(CO)COC1']), 0)
        cfg = json.loads((guessed / 'project.json').read_text())
        self.assertEqual((cfg['chemistry'], cfg['reactive_sites_per_molecule']), ('oxetane', 1))
        self.assertIn('oxetane_OH_relay', (guessed / 'definitions/reactions.yaml').read_text())
        # A preset whose ring the monomer lacks is refused; a monomer with both rings needs --chemistry.
        self.assertEqual(quiet(['monomer', str(self.root / 'wrong'), '--platform', 'CPU',
                                '--smiles', 'CCC1(CO)COC1', '--chemistry', 'epoxy']), 1)
        self.assertFalse((self.root / 'wrong').exists())
        both = 'C1OC1COCC2(CC)COC2'
        self.assertEqual(quiet(['monomer', str(self.root / 'both'), '--platform', 'CPU', '--smiles', both]), 1)
        self.assertEqual(quiet(['monomer', str(self.root / 'both'), '--platform', 'CPU', '--smiles', both,
                                '--chemistry', 'oxetane']), 0)
        self.assertEqual(monomer.ring_label(dict(chemistry='oxetane', build=dict(smiles=both))), 'Reactive ring')
        # The ring label comes from the molecule, not from the preset.
        self.assertEqual(monomer.ring_label(dict(chemistry='epoxy', build=dict(smiles='CCC1(CO)COC1'))), 'Oxetane')
        self.assertEqual(monomer.ring_label(dict(chemistry='oxetane')), 'Oxetane')

    def test_cutoff_and_relaxation_options_reach_crosslinking(self):
        from workflows.crosslinking.crosslinking import load_relaxation_stages
        stages = self.root / 'relax.yaml'
        stages.write_text('stages:\n'
            '  - {ensemble: NVT, steps: 5000, timestep: 0.1}\n'
            '  - {ensemble: NVT, steps: 5000, timestep: 0.5, temperature: 350}\n'
            '  - {ensemble: NPT, steps: 20000, timestep: 1.0, pressure_atm: 1}\n')
        full = self.root / 'full'
        self.assertEqual(quiet(['monomer', str(full), '--platform', 'CPU']), 0)
        self.assertEqual(quiet(['crosslink', str(full), '--activation-fraction', '0.1', '--max-reaction-cutoff', '4',
            '--post-reaction-md-steps', '500', '--reaction-minimization-iterations', '5000',
            '--relaxation-stages', str(stages)]), 0)
        cfg = json.loads((full / 'crosslinks/af0.1/project.json').read_text())
        options = cfg['crosslinking']
        self.assertEqual((options['reaction_cutoff'], options['max_reaction_cutoff']), (4., 4.))
        self.assertEqual(options['post_reaction_md_steps'], 500)
        self.assertEqual(options['minimization_iterations'], 5000)
        self.assertEqual([s['timestep'] for s in options['relaxation_stages']], [.1, .5, 1.])
        command = build_command(self.root, 'network', options, seed=1, platform='CPU')
        self.assertEqual(load_relaxation_stages(command[command.index('--relaxation-stages') + 1]),
                         options['relaxation_stages'])
        for bad in ('[{"ensemble":"NVT","steps":10}]',
                    '[{"ensemble":"NVT","steps":10,"timestep":0.1,"pressure_atm":1}]'):
            self.assertEqual(quiet(['crosslink', str(full), '--activation-fraction', '0.2',
                                    '--relaxation-stages', bad]), 1)
        self.assertEqual(quiet(['crosslink', str(full), '--activation-fraction', '0.2', '--reaction-cutoff', '6',
                                '--max-reaction-cutoff', '5']), 1)
        self.assertFalse((full / 'crosslinks/af0.2').exists())
        cfg['crosslinking']['unknown_option'] = 1
        with self.assertRaises(ValueError):
            validate(cfg)
        del cfg['crosslinking']['unknown_option']
        cfg['crosslinking']['total_atoms'] = 500
        with self.assertRaisesRegex(ValueError, 'starts from the monomer structure'):
            validate(cfg)

    def test_rejects_inconsistent_monomer_settings(self):
        for patch_ in ({'temperature_K': 900.}, {'solvents': [dict(name='water', smiles='O')]},
                       {'t_production_ps': .001}):
            with self.subTest(patch=patch_), self.assertRaises(ValueError):
                validate(dict(self.cfg, **patch_))
        broken = json.loads(json.dumps(self.cfg))
        broken['build']['cycles'] = 3
        with self.assertRaisesRegex(ValueError, 'no start sites'):
            validate(broken)


class BuildRoutingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.monomer = Path(self.temp.name) / 'TMPTGE'
        self.assertEqual(quiet(['monomer', str(self.monomer), '--replicates', '2', '--seed', '100',
                                '--platform', 'CPU', '--smoke']), 0)
        self.assertEqual(quiet(['crosslink', str(self.monomer), '--activation-fraction', '0.1']), 0)
        self.cfg = json.loads((self.monomer / 'project.json').read_text())
        self.child = self.monomer / 'crosslinks/af0.1'
        self.child_cfg = json.loads((self.child / 'project.json').read_text())
        self.commands = []

    def fake_popen(self, command, cwd, *, stderr_lines=(), returncode=0, **options):
        self.commands.append(command)
        if not returncode:
            self.write_build(command, cwd)
        child = MagicMock()
        child.__enter__.return_value = child
        child.stderr = iter(stderr_lines)
        child.returncode = returncode
        return child

    def write_build(self, command, cwd):
        seed = int(command[command.index('--seed') + 1])
        output = Path(cwd) / 'results' / command[command.index('--project') + 1]
        output.mkdir(parents=True, exist_ok=True)
        (output / 'final.structure.json').write_text(json.dumps(dict(
            graphs=[dict(atoms=[dict(partial_charge=0.)], bonds=[])], positions=[[seed, 0, 0]],
            box_vectors=[[10, 0, 0], [0, 10, 0], [0, 0, 10]])))
        start = command[command.index('--start-structure') + 1] if '--start-structure' in command else None
        (output / 'settings.json').write_text(json.dumps(dict(seed=seed, actual_total_atoms=519, monomers=11,
            initial_reactive_sites=33, smiles='C', temperature_K=363.15, activation_count=3,
            actual_activation_fraction=3 / 33,
            start_structure_sha256=hashlib.sha256(Path(start).read_bytes()).hexdigest() if start else None)))
        (output / 'metrics.json').write_text(json.dumps(dict(conversion=.75, cycles_completed=1,
            epoxide_consumption_including_activation=.8)))

    def popen(self, **kwargs):
        return patch('workflows.network_properties.construction.subprocess.Popen',
                     side_effect=lambda command, cwd, **o: self.fake_popen(command, cwd, **kwargs, **o))

    def test_monomer_builds_have_no_starts_and_are_reused_and_verified(self):
        with self.popen(), redirect_stdout(io.StringIO()):
            first = ensure_monomer_build(self.monomer, self.cfg, 1)
            ensure_monomer_build(self.monomer, self.cfg, 2)
            ensure_monomer_build(self.monomer, self.cfg, 1)
        self.assertEqual(len(self.commands), 2)
        command = self.commands[0]
        for flag, value in (('--activation-count', '0'), ('--cycles', '0'), ('--seed', '101'),
                            ('--final-md-steps', '0'), ('--temperature', '363.15'), ('--project', 'monomer')):
            self.assertEqual(command[command.index(flag) + 1], value, flag)
        self.assertNotIn('--start-structure', command)
        self.assertEqual(first, network_path(self.monomer, self.cfg, 1))
        done = json.loads((self.monomer / 'replica_01/build/done.json').read_text())
        self.assertEqual((done['seed'], done['molecules'], done['reactive_sites']), (101, 11, 33))
        first.write_text('{}')
        with self.assertRaisesRegex(ValueError, 'changed'):
            ensure_monomer_build(self.monomer, self.cfg, 1)

    def test_crosslink_replica_starts_from_the_same_numbered_monomer_replica(self):
        with self.assertRaisesRegex(RuntimeError, 'replica_02 has not finished'):
            ensure_crosslink(self.child, self.child_cfg, 2)
        complete_monomer_replica(self.monomer, self.cfg, 1, 'one')
        complete_monomer_replica(self.monomer, self.cfg, 2, 'two')
        with self.popen(), redirect_stdout(io.StringIO()):
            ensure_crosslink(self.child, self.child_cfg, 2)
            ensure_crosslink(self.child, self.child_cfg, 2)
        self.assertEqual(len(self.commands), 1)
        command = self.commands[0]
        start = Path(command[command.index('--start-structure') + 1])
        self.assertEqual(start, monomer.relaxed_path(self.monomer, self.cfg, 2))
        self.assertEqual(command[command.index('--seed') + 1], '102')
        self.assertEqual(command[command.index('--activation-fraction') + 1], '0.1')
        self.assertNotIn('--total-atoms', command)
        origin = json.loads((self.child / 'replica_02/network_origin.json').read_text())
        self.assertEqual(origin['monomer']['replica'], 'replica_02')
        self.assertEqual(origin['monomer']['structure_sha256'], hashlib.sha256(start.read_bytes()).hexdigest())
        self.assertEqual((origin['requested_activation_fraction'], origin['activation_count']), (.1, 3))
        start.write_text(json.dumps(dict(graphs=['changed'])))
        with self.assertRaisesRegex(ValueError, 'changed after this network was built'):
            ensure_crosslink(self.child, self.child_cfg, 2)

    def test_changed_monomer_configuration_is_refused(self):
        complete_monomer_replica(self.monomer, self.cfg, 1)
        (self.monomer / 'project.json').write_text(json.dumps(dict(self.cfg, seed=999)))
        with self.assertRaisesRegex(ValueError, 'changed after this crosslink project was created'):
            ensure_crosslink(self.child, self.child_cfg, 1)

    def test_failed_build_reports_the_crosslinking_error(self):
        lines = ['Traceback (most recent call last):\n', 'ValueError: compression_pressure_mpa must be positive\n']
        with self.popen(stderr_lines=lines, returncode=1), redirect_stderr(io.StringIO()), \
             redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(RuntimeError, 'replica 1: ValueError: compression_pressure_mpa'):
                ensure_monomer_build(self.monomer, self.cfg, 1)
        self.assertFalse((self.monomer / 'replica_01/build/done.json').exists())

    def test_rule_files_are_copied_and_input_changes_are_detectable(self):
        path = self.monomer / 'definitions/reactions.yaml'
        path.write_text(path.read_text() + '\n# modified\n')
        self.assertNotEqual(input_digest(self.monomer, self.cfg), self.cfg['input_sha256'])

    def test_stage_order_of_monomer_and_crosslink_replicas(self):
        events = []
        np = 'workflows.network_properties'
        with patch(f'{np}.network_properties.initialize_project', return_value=self.cfg), \
             patch(f'{np}.construction.ensure_monomer_build', side_effect=lambda p, c, r, lock=None: events.append(('build', r))), \
             patch(f'{np}.network_properties.dry_properties', side_effect=lambda p, c, r: events.append(('properties', r))), \
             patch(f'{np}.monomer.temperature_properties', side_effect=lambda p, c, r: events.append(('T', r))), \
             redirect_stdout(io.StringIO()):
            run(Namespace(project=self.monomer, replica=None, defer_report=True))
        self.assertEqual(events, [('build', 1), ('properties', 1), ('T', 1), ('build', 2), ('properties', 2), ('T', 2)])
        events.clear()
        with patch(f'{np}.network_properties.initialize_project', return_value=self.child_cfg), \
             patch(f'{np}.construction.ensure_crosslink', side_effect=lambda p, c, r, lock=None: events.append(('crosslink', r))), \
             patch(f'{np}.network_properties.dry_properties', side_effect=lambda p, c, r: events.append(('properties', r))), \
             redirect_stdout(io.StringIO()):
            run(Namespace(project=self.child, replica=2, defer_report=True))
        self.assertEqual(events, [('crosslink', 2), ('properties', 2)])

    def test_heating_grid_reaches_T_in_bounded_steps(self):
        grid = monomer.heating_grid(300., 363.15, 20.)
        self.assertEqual(len(grid), 4)
        for value, expected in zip(grid, [315.7875, 331.575, 347.3625, 363.15]):
            self.assertAlmostEqual(value, expected)
        self.assertEqual(grid[-1], 363.15)
        self.assertEqual(monomer.heating_grid(300., 310., 20.), [310.])
        self.assertEqual(monomer.heating_grid(300., 260., 20.), [280., 260.])


if __name__ == '__main__':
    unittest.main()
