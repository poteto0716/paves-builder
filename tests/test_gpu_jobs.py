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
import time
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
            mechanical=False, solvents=[], network='input.structure.json',
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
        with redirect_stdout(io.StringIO()):
            self.assertEqual(main(['network',str(output),'--network',str(source),'--platform','CUDA',
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

    def test_submit_refuses_while_a_replica_lock_is_held(self):
        args = Namespace(project=self.project,gpus='0,2',dry_run=False)
        with file_lock(self.project/'replica_02/.run.lock'), \
             patch.object(jobs,'visible_gpus',return_value=self.devices), \
             patch.object(jobs.subprocess,'Popen') as start:
            with self.assertRaisesRegex(RuntimeError,r'Replicas \[2\]'):
                jobs.submit(args)
        start.assert_not_called()

    def test_worker_stops_when_the_launcher_is_killed(self):
        script = '''import json,sys,time
from pathlib import Path
Path(sys.argv[2]).write_text(sys.argv[1])
time.sleep(60)
'''
        plan = jobs.plan(self.project,self.cfg,self.devices[:1])
        path = self.project/'submission.json'
        path.write_text(json.dumps(plan))
        started = self.project/'worker_started'
        launcher = '''import sys
from workflows.network_properties import gpu_jobs as jobs
worker = [sys.executable,'-c',sys.argv[2],'x',sys.argv[3]]
jobs.dispatch(sys.argv[1],command_factory=lambda p,r:worker,poll_interval=.01,refresh_report=False)
'''
        process = subprocess.Popen([sys.executable,'-c',launcher,str(path),script,str(started)],
                                   cwd=jobs.ROOT,env=dict(os.environ,PYTHONPATH=str(jobs.ROOT)))
        try:
            deadline = time.time()+20
            while not started.exists() and time.time() < deadline:
                time.sleep(.01)
            self.assertTrue(started.exists())
            pid = json.loads((self.project/'jobs/status.json').read_text())['jobs'][0]['pid']
        finally:
            process.kill()
            process.wait()
        deadline = time.time()+10
        while jobs.process_state(pid) not in (None,'Z') and time.time() < deadline:
            time.sleep(.01)
        self.assertIn(jobs.process_state(pid),(None,'Z'))

    def test_refused_duplicate_does_not_mark_a_live_replica_failed(self):
        from workflows.network_properties.execution import read_progress, set_stage
        base = self.project/'replica_01'
        with redirect_stdout(io.StringIO()):
            set_stage(base,'running elsewhere')
        job = dict(replica=1,error='refused')
        with file_lock(base/'.run.lock'):
            jobs.mark_interrupted(self.project,job)
            self.assertEqual(read_progress(base)['state'],'running')
        jobs.mark_interrupted(self.project,job)
        self.assertEqual(read_progress(base)['state'],'failed')

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
if replica==fail:
    print('ValueError: injected failure', file=__import__('sys').stderr)
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
        failed=result['jobs'][1]
        self.assertEqual(failed['error'],'ValueError: injected failure')
        self.assertIn('ValueError: injected failure',Path(failed['error_log']).read_text())
        self.assertNotIn('error',result['jobs'][2])
        with redirect_stdout(io.StringIO()) as output:
            jobs.status(Namespace(project=self.project))
        self.assertIn('cause: ValueError: injected failure',output.getvalue())

    def test_signal_exit_and_progress_records(self):
        from workflows.network_properties.execution import replica_progress, set_stage, read_progress
        self.assertIn('SIGKILL',jobs.failure_summary('',-9))
        self.assertEqual(jobs.failure_summary('x\nRuntimeError: late\n',1),'RuntimeError: late')
        base=self.project/'replica_01'
        with redirect_stdout(io.StringIO()), self.assertRaises(ZeroDivisionError):
            with replica_progress(base):
                set_stage(base,'thermal history: 03_cool_000_sample',temperature_K=800.)
                1/0
        progress=read_progress(base)
        self.assertEqual(progress['state'],'failed')
        self.assertEqual(progress['stage'],'thermal history: 03_cool_000_sample')
        text=(base/'error.log').read_text()
        self.assertIn('during stage: thermal history: 03_cool_000_sample',text)
        self.assertIn('ZeroDivisionError',text)
        with redirect_stdout(io.StringIO()):
            with replica_progress(base):
                set_stage(base,'mechanical: elastic constants')
        progress=read_progress(base)
        self.assertEqual((progress['state'],progress['history'][-1]['stage']),
                         ('completed','mechanical: elastic constants'))

    def test_status_flags_lost_and_zombie_workers(self):
        child=subprocess.Popen([sys.executable,'-c','pass'])
        while jobs.process_state(child.pid)!='Z':
            pass
        self.assertEqual(jobs.process_state(child.pid),'Z')
        child.wait()
        self.assertIsNone(jobs.process_state(child.pid))

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

    def test_reseed_archives_only_the_failed_replica_and_routes_its_new_seed(self):
        from workflows.network_properties.execution import finish_progress, replica_seed, set_stage
        base=self.project/'replica_02'
        with redirect_stdout(io.StringIO()):
            set_stage(base,'thermal history: 02_high_equilibration')
        finish_progress(base,'failed',error='OpenMMException: Particle coordinate is NaN')
        (base/'dry_marker').write_text('old')
        (self.project/'jobs').mkdir()
        (self.project/'jobs/replica_02.error.log').write_text('trace')
        args=lambda **k: Namespace(project=self.project,replica=[2],seed=None,force=False,submit=False,
                                   gpus=None,dry_run=False,**k)
        with redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError,'has not failed'):
                jobs.reseed(Namespace(**{**vars(args()),'replica':[1]}))
            jobs.reseed(args())
        # base seed 100, 3 replicas: first attempts 101..103, the retry of replica 2 is 100+3+2.
        self.assertEqual(replica_seed(self.project,self.cfg,2),105)
        self.assertEqual(replica_seed(self.project,self.cfg,1),101)
        archive=self.project/'failed_attempts/replica_02_seed102_1'
        self.assertEqual((archive/'dry_marker').read_text(),'old')
        self.assertTrue((archive/'job.error.log').is_file())
        self.assertFalse(base.exists())
        self.assertEqual([j['seed'] for j in jobs.plan(self.project,self.cfg,self.devices)['jobs']],[101,105,103])
        with patch('workflows.network_properties.network_properties.initialize_project',return_value=self.cfg), \
             patch('workflows.network_properties.network_properties.dry_properties'), \
             patch('workflows.network_properties.network_properties.report'):
            run(Namespace(project=self.project,replica=2,defer_report=True))
        self.assertEqual(json.loads((base/'execution.json').read_text())['seed'],105)
        # The archive is not a replica directory and is never collected as one.
        self.assertEqual(list(self.project.glob('replica_*/dry_marker')),[])
        self.assertTrue(all(p.is_dir() for p in self.project.glob('replica_*')))
        # A second failure gets yet another unused seed; an explicit used seed is refused.
        with redirect_stdout(io.StringIO()):
            finish_progress(base,'failed',error='again')
            with self.assertRaisesRegex(ValueError,'already used'):
                jobs.reseed(Namespace(**{**vars(args()),'seed':102}))
            jobs.reseed(args())
        self.assertEqual(replica_seed(self.project,self.cfg,2),108)
        self.assertTrue((self.project/'failed_attempts/replica_02_seed105_1').is_dir())
        with file_lock(self.project/'replica_03/.run.lock'), redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(RuntimeError,'running'):
                jobs.reseed(Namespace(**{**vars(args()),'replica':[3],'force':True}))

    def test_monomer_submission_runs_crosslink_replicas_after_their_monomer_replica(self):
        from workflows.network_properties import monomer
        root = self.project/'TMPTGE'
        with redirect_stdout(io.StringIO()):
            self.assertEqual(main(['monomer',str(root),'--platform','CUDA','--replicates','2','--seed','100','--smoke']),0)
            self.assertEqual(main(['crosslink',str(root),'--activation-fraction','0.1']),0)
        cfg = json.loads((root/'project.json').read_text())
        plan = jobs.plan(root,cfg,self.devices)
        self.assertEqual([(j['name'],j['replica'],j['seed']) for j in plan['jobs']],
                         [('.',1,101),('.',2,102),('crosslinks/af0.1',1,101),('crosslinks/af0.1',2,102)])
        self.assertEqual(plan['jobs'][3]['after'],dict(project=str(root),replica=2))
        self.assertEqual(len(plan['configurations']),2)
        alone = jobs.plan(root,cfg,self.devices,monomer_only=True)
        self.assertEqual([(j['name'],j['replica']) for j in alone['jobs']],[('.',1),('.',2)])
        self.assertEqual(list(alone['configurations']),[str(root)])
        with self.assertRaisesRegex(ValueError,'monomer-only'):
            jobs.plan(root/'crosslinks/af0.1',json.loads((root/'crosslinks/af0.1/project.json').read_text()),
                      self.devices,monomer_only=True)
        # Fake workers: a monomer replica writes its completion files (replica 2 fails); a crosslink
        # replica exits 9 if its monomer replica is not complete.
        script = """import json,sys,time
from pathlib import Path
project=Path(sys.argv[1]); replica=int(sys.argv[2]); T=sys.argv[3]
cfg=json.loads((project/'project.json').read_text()); base=project/f'replica_{replica:02d}'
root=project if cfg['input_mode']=='monomer' else (project/cfg['parent']).resolve()
with open(root/'order.txt','a') as fh: fh.write(f"{cfg['input_mode']} {replica}\\n")
time.sleep(.05)
if cfg['input_mode']=='monomer':
    if replica==2: raise SystemExit(7)
    for p in ('build/done.json','dry_result.json',f'{T}/result.json',f'{T}/relaxed_{T}.structure.json'):
        (base/p).parent.mkdir(parents=True,exist_ok=True); (base/p).write_text('{}')
else:
    if not (root/f'replica_{replica:02d}/{T}/result.json').is_file(): raise SystemExit(9)
"""
        factory = lambda project,replica:[sys.executable,'-c',script,str(project),str(replica),monomer.label(cfg)]
        path = root/'submission.json'
        path.write_text(json.dumps(plan))
        result = jobs.dispatch(path,command_factory=factory,poll_interval=.01,refresh_report=False)
        states = [(j['name'],j['replica'],j['state']) for j in result['jobs']]
        self.assertEqual(states,[('.',1,'completed'),('.',2,'failed'),
                                 ('crosslinks/af0.1',1,'completed'),('crosslinks/af0.1',2,'blocked')])
        self.assertIn('replica_02 failed',result['jobs'][3]['error'])
        self.assertEqual(result['state'],'failed')
        order = (root/'order.txt').read_text().split('\n')
        self.assertLess(order.index('monomer 1'),order.index('crosslink 1'))
        self.assertNotIn('crosslink 2',order)
        with redirect_stdout(io.StringIO()) as output:
            jobs.status(Namespace(project=root))
        self.assertIn('crosslinks/af0.1 [crosslink]: 0/2 replicas complete',output.getvalue())
        self.assertIn('blocked',output.getvalue())

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
