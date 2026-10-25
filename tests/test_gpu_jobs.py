"""Local GPU scheduling tests use real CPU subprocesses and mocked GPU inventory."""
from argparse import Namespace
from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from workflows.network_properties import gpu_jobs as jobs
from workflows.network_properties.execution import file_lock
from workflows.network_properties.network_properties import DEFAULTS, digest, run, main


class GPUJobTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.project = Path(self.temp.name)
        (self.project/'input.structure.json').write_text('{}')
        self.cfg = dict(DEFAULTS, platform='CUDA', replicates=3, seed=100,
            mechanical=False, thermal_conductivity=False, solvents=[], network='input.structure.json',
            input_sha256=digest(self.project/'input.structure.json'))
        (self.project/'project.json').write_text(json.dumps(self.cfg))
        self.devices = [dict(index='0',uuid='GPU-zero'),dict(index='2',uuid='GPU-two')]

    def test_physical_gpu_selection_and_parent_visibility_mask(self):
        inventory = subprocess.CompletedProcess([],0,'0, GPU-zero\n1, GPU-one\n2, GPU-two\n','')
        with patch.dict(os.environ,{'CUDA_VISIBLE_DEVICES':'GPU-two,GPU-zero'}), \
             patch.object(jobs.subprocess,'run',return_value=inventory):
            self.assertEqual(jobs.visible_gpus('2,0'),list(reversed(self.devices)))
            for invalid in ('1','0,0','-1','0;echo bad'):
                with self.subTest(invalid=invalid),self.assertRaises(ValueError):
                    jobs.visible_gpus(invalid)

    def test_dry_run_has_unique_seeds_and_launches_nothing(self):
        args = Namespace(project=self.project,gpus='0,2',dry_run=True)
        with patch.object(jobs,'visible_gpus',return_value=self.devices), \
             patch.object(jobs.subprocess,'Popen') as start, redirect_stdout(io.StringIO()) as output:
            jobs.submit(args)
        plan = json.loads(output.getvalue())
        self.assertEqual([j['seed'] for j in plan['jobs']],[101,102,103])
        self.assertEqual(plan['max_parallel'],2)
        self.assertEqual([d['index'] for d in plan['gpus']],['0','2'])
        start.assert_not_called()
        self.assertFalse((self.project/'run_signature.json').exists())

    def test_explicit_replica_count_is_respected_in_smoke_mode(self):
        source=self.project/'source.json'
        source.write_text(json.dumps(dict(graphs=[dict(atoms=[{'partial_charge':0.}])],
            positions=[[0,0,0]],box_vectors=[[3,0,0],[0,3,0],[0,0,3]])))
        output=self.project/'smoke'
        self.assertEqual(main(['new',str(output),'--network',str(source),'--platform','CUDA',
                               '--replicates','3','--smoke']),0)
        self.assertEqual(json.loads((output/'project.json').read_text())['replicates'],3)

    def test_duplicate_submission_is_rejected_before_launch(self):
        args = Namespace(project=self.project,gpus='0,2',dry_run=False)
        with file_lock(self.project/'jobs/.launcher.lock'), \
             patch.object(jobs,'visible_gpus',return_value=self.devices), \
             patch.object(jobs.subprocess,'Popen') as start:
            with self.assertRaisesRegex(RuntimeError,'already running'):
                jobs.submit(args)
        start.assert_not_called()

    def command(self, fail_replica=None):
        # An exclusive marker detects accidental overlap of workers on the same GPU.
        script = '''import json,os,sys,time
from pathlib import Path
project=Path(sys.argv[1]); replica=int(sys.argv[2]); fail=int(sys.argv[3])
gpu=os.environ['CUDA_VISIBLE_DEVICES']; marker=project/('busy_'+gpu)
with marker.open('x'):
    try:
        cfg=json.loads((project/'project.json').read_text())
        (project/f'observed_{replica}.json').write_text(json.dumps(dict(
            gpu=gpu,seed=cfg['seed']+replica,replica=replica)))
        time.sleep(.05)
    finally:
        marker.unlink()
raise SystemExit(7 if replica==fail else 0)
'''
        return lambda project,replica:[sys.executable,'-c',script,str(project),str(replica),str(fail_replica or 0)]

    def dispatch(self, fail_replica=None):
        plan = jobs.plan(self.project,self.cfg,self.devices)
        path = self.project/'submission.json'
        path.write_text(json.dumps(plan))
        return jobs.dispatch(path,command_factory=self.command(fail_replica),poll_interval=.01,refresh_report=False)

    def test_queue_runs_every_replica_once_without_gpu_overlap(self):
        result = self.dispatch()
        self.assertEqual(result['state'],'completed')
        self.assertTrue(all(j['state']=='completed' for j in result['jobs']))
        observed=[json.loads((self.project/f'observed_{i}.json').read_text()) for i in (1,2,3)]
        self.assertEqual([r['seed'] for r in observed],[101,102,103])
        self.assertEqual(observed[0]['gpu'],'GPU-zero')
        self.assertEqual(observed[1]['gpu'],'GPU-two')
        self.assertIn(observed[2]['gpu'],['GPU-zero','GPU-two'])
        self.assertEqual(json.loads((self.project/'jobs/status.json').read_text())['state'],'completed')
        self.assertTrue(all(Path(j['log']).is_file() for j in result['jobs']))

    def test_failed_job_is_recorded_and_other_jobs_finish(self):
        result = self.dispatch(fail_replica=2)
        self.assertEqual(result['state'],'failed')
        self.assertEqual(result['jobs'][1]['returncode'],7)
        self.assertEqual(result['jobs'][2]['state'],'completed')

    def test_single_replica_run_uses_correct_seed_and_defers_shared_report(self):
        args = Namespace(project=self.project,replica=2,defer_report=True)
        with patch('workflows.network_properties.network_properties.initialize_project',return_value=self.cfg), \
             patch('workflows.network_properties.network_properties.dry_properties') as dry, \
             patch('workflows.network_properties.network_properties.report') as report:
            run(args)
        dry.assert_called_once_with(self.project,self.cfg,2)
        report.assert_not_called()
        self.assertEqual(json.loads((self.project/'replica_02/execution.json').read_text())['seed'],102)
        self.assertFalse((self.project/'replica_01').exists())

    def test_worker_command_and_existing_replica_lock(self):
        command=jobs.worker_command(self.project,3)
        self.assertIn('--defer-report',command)
        self.assertEqual(command[command.index('--replica')+1],'3')
        with file_lock(self.project/'replica_02/.run.lock'), \
             patch('workflows.network_properties.network_properties.initialize_project',return_value=self.cfg), \
             patch('workflows.network_properties.network_properties.dry_properties') as dry:
            with self.assertRaisesRegex(RuntimeError,'Already running'):
                run(Namespace(project=self.project,replica=2,defer_report=True))
        dry.assert_not_called()

    def test_completed_replicas_are_skipped_in_submission_plan(self):
        (self.project/'replica_02').mkdir()
        (self.project/'replica_02/dry_result.json').write_text('{}')
        plan=jobs.plan(self.project,self.cfg,self.devices)
        self.assertEqual([j['state'] for j in plan['jobs']],['queued','completed','queued'])

    def test_changed_configuration_fails_before_starting_any_worker(self):
        plan=jobs.plan(self.project,self.cfg,self.devices)
        path=self.project/'submission.json'
        path.write_text(json.dumps(plan))
        (self.project/'project.json').write_text(json.dumps(dict(self.cfg,seed=999)))
        with patch.object(jobs.subprocess,'Popen') as start, self.assertRaisesRegex(ValueError,'changed'):
            jobs.dispatch(path,refresh_report=False)
        start.assert_not_called()
        self.assertEqual(json.loads((self.project/'jobs/status.json').read_text())['state'],'failed')


if __name__ == '__main__':
    unittest.main()
