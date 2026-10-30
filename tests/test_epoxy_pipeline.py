"""Independent epoxy inputs, seed routing, stage ordering and restart provenance."""
from argparse import Namespace
from contextlib import redirect_stderr
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from workflows.network_properties.network_properties import main, input_digest, run
from workflows.network_properties.epoxy_pipeline import ensure_network, network_path, build_command


class EpoxyPipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.project=Path(self.temp.name)/'project'
        self.assertEqual(main(['new',str(self.project),'--replicates','2','--seed','100',
            '--platform','CPU','--smoke','--skip-mechanical']),0)
        self.cfg=json.loads((self.project/'project.json').read_text())

    def test_new_defaults_to_independent_epoxy_builds(self):
        self.assertEqual(self.cfg['input_mode'],'epoxy')
        self.assertEqual(self.cfg['crosslinking']['total_atoms'],500)  # smoke-sized
        self.assertEqual(self.cfg['crosslinking']['target_conversion'],.995)
        self.assertFalse(self.cfg['crosslinking']['save_trajectory'])
        self.assertNotIn('network',self.cfg)
        self.assertEqual(input_digest(self.project,self.cfg),self.cfg['input_sha256'])
        full=Path(self.temp.name)/'full'
        self.assertEqual(main(['new',str(full),'--platform','CUDA','--replicates','4']),0)
        options=json.loads((full/'project.json').read_text())['crosslinking']
        self.assertEqual(options['total_atoms'],20000)
        self.assertEqual(options['atom_count_mode'],'above')
        self.assertEqual(options['cycles'],5000)

    def fake_popen(self, command, cwd, *, stderr_lines=(), returncode=0, **options):
        if not returncode:
            self.write_build(command, cwd)
        child=MagicMock()
        child.__enter__.return_value=child
        child.stderr=iter(stderr_lines)
        child.returncode=returncode
        return child

    def write_build(self, command, cwd):
        seed=int(command[command.index('--seed')+1])
        output=Path(cwd)/'results/epoxy'
        output.mkdir(parents=True,exist_ok=True)
        (output/'final.structure.json').write_text(json.dumps(dict(
            graphs=[dict(atoms=[dict(partial_charge=0.)],bonds=[])],positions=[[seed,0,0]],
            box_vectors=[[10,0,0],[0,10,0],[0,0,10]])))
        (output/'settings.json').write_text(json.dumps(dict(seed=seed,actual_total_atoms=519)))
        (output/'metrics.json').write_text(json.dumps(dict(conversion=.75,cycles_completed=1,
            epoxide_consumption_including_activation=.8)))

    def test_each_replica_gets_its_own_build_and_completed_build_is_reused(self):
        with patch('workflows.network_properties.epoxy_pipeline.subprocess.Popen',side_effect=self.fake_popen) as build:
            first=ensure_network(self.project,self.cfg,1)
            second=ensure_network(self.project,self.cfg,2)
            self.assertNotEqual(first,second)
            self.assertNotEqual(first.read_bytes(),second.read_bytes())
            ensure_network(self.project,self.cfg,1)
        self.assertEqual(build.call_count,2)
        for replica in (1,2):
            record=json.loads((self.project/f'replica_{replica:02d}/network_origin.json').read_text())
            self.assertEqual(record['seed'],100+replica)
            self.assertFalse(record['target_reached'])
        self.assertEqual(network_path(self.project,self.cfg,2),second)

    def test_incomplete_build_restarts_but_modified_completed_build_is_rejected(self):
        source=network_path(self.project,self.cfg,1)
        source.parent.mkdir(parents=True)
        self.assertIn('--force',build_command(self.project,self.cfg,1))
        with patch('workflows.network_properties.epoxy_pipeline.subprocess.Popen',side_effect=self.fake_popen):
            ensure_network(self.project,self.cfg,1)
        source.write_text('{}')
        with self.assertRaisesRegex(ValueError,'changed'):
            ensure_network(self.project,self.cfg,1)

    def test_failed_build_reports_the_crosslinking_error(self):
        lines=['Traceback (most recent call last):\n','ValueError: compression_pressure_mpa must be positive\n']
        failing=lambda command,cwd,**o:self.fake_popen(command,cwd,stderr_lines=lines,returncode=1,**o)
        with patch('workflows.network_properties.epoxy_pipeline.subprocess.Popen',side_effect=failing), \
             redirect_stderr(io.StringIO()):
            with self.assertRaisesRegex(RuntimeError,'replica 1: ValueError: compression_pressure_mpa'):
                ensure_network(self.project,self.cfg,1)
        self.assertFalse((self.project/'replica_01/crosslinking/done.json').exists())

    def test_rule_files_are_copied_and_input_changes_are_detectable(self):
        original=self.cfg['input_sha256']
        path=self.project/'crosslinking/reactions.yaml'
        path.write_text(path.read_text()+'\n# modified\n')
        self.assertNotEqual(input_digest(self.project,self.cfg),original)

    def test_worker_constructs_before_its_own_physical_properties(self):
        events=[]
        with patch('workflows.network_properties.network_properties.initialize_project',return_value=self.cfg), \
             patch('workflows.network_properties.epoxy_pipeline.ensure_network',side_effect=lambda p,c,r,lock=None:events.append(("build",r))), \
             patch('workflows.network_properties.network_properties.dry_properties',side_effect=lambda p,c,r:events.append(('properties',r))):
            run(Namespace(project=self.project,replica=None,defer_report=True))
        self.assertEqual(events,[('build',1),('properties',1),('build',2),('properties',2)])

    def test_cutoff_and_relaxation_options_reach_crosslinking(self):
        from workflows.crosslinking.crosslinking import load_relaxation_stages
        stages=Path(self.temp.name)/'relax.yaml'
        stages.write_text('stages:\n'
            '  - {ensemble: NVT, steps: 5000, timestep: 0.1}\n'
            '  - {ensemble: NVT, steps: 5000, timestep: 0.5, temperature: 350}\n'
            '  - {ensemble: NPT, steps: 20000, timestep: 1.0, pressure_atm: 1}\n')
        full=Path(self.temp.name)/'careful'
        self.assertEqual(main(['new',str(full),'--platform','CPU','--max-reaction-cutoff','4',
            '--post-reaction-md-steps','500','--reaction-minimization-iterations','5000',
            '--relaxation-stages',str(stages)]),0)
        cfg=json.loads((full/'project.json').read_text())
        options=cfg['crosslinking']
        self.assertEqual((options['reaction_cutoff'],options['max_reaction_cutoff']),(4.,4.))
        self.assertEqual(options['post_reaction_md_steps'],500)
        self.assertEqual(options['minimization_iterations'],5000)
        self.assertEqual([s['timestep'] for s in options['relaxation_stages']],[.1,.5,1.])
        command=build_command(full,cfg,1)
        self.assertEqual(load_relaxation_stages(command[command.index('--relaxation-stages')+1]),
                         options['relaxation_stages'])
        self.assertEqual(command[command.index('--max-reaction-cutoff')+1],'4.0')
        for bad in ('[{"ensemble":"NVT","steps":10}]',
                    '[{"ensemble":"NVT","steps":10,"timestep":0.1,"pressure_atm":1}]'):
            with redirect_stderr(io.StringIO()):
                self.assertEqual(main(['new',str(Path(self.temp.name)/'bad'),'--relaxation-stages',bad]),1)
        with redirect_stderr(io.StringIO()):
            self.assertEqual(main(['new',str(Path(self.temp.name)/'bad'),'--reaction-cutoff','6',
                                   '--max-reaction-cutoff','5']),1)
        self.assertFalse((Path(self.temp.name)/'bad').exists())

    def test_project_without_newer_crosslinking_options_still_validates(self):
        from workflows.network_properties.network_properties import validate
        cfg=json.loads(json.dumps(self.cfg))
        for name in ('relaxation_stages','post_reaction_md_steps','compression_steps'):
            del cfg['crosslinking'][name]
        validate(cfg)
        cfg['crosslinking']['unknown_option']=1
        with self.assertRaises(ValueError):
            validate(cfg)

    def test_saved_network_cannot_silently_ignore_build_options(self):
        output=Path(self.temp.name)/'invalid'
        with redirect_stderr(io.StringIO()) as error:
            result=main(['new',str(output),'--network','saved.structure.json','--total-atoms','500'])
        self.assertEqual(result,1)
        self.assertIn('cannot be combined',error.getvalue())
        self.assertFalse(output.exists())


if __name__=='__main__':
    unittest.main()
