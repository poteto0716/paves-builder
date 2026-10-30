"""Detached local GPU queue: one replica per selected GPU, explicit seeds and logs."""
from __future__ import annotations

from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
SCRIPT = Path(__file__).with_name('network_properties.py')


def now():
    return datetime.now(timezone.utc).isoformat()


def visible_gpus(selection=None):
    try:
        result = subprocess.run(['nvidia-smi', '--query-gpu=index,uuid', '--format=csv,noheader,nounits'],
            capture_output=True, text=True, check=True, timeout=10)
    except (OSError, subprocess.SubprocessError) as exc:
        raise RuntimeError(f'Cannot enumerate NVIDIA GPUs: {exc}') from exc
    devices = []
    for line in result.stdout.splitlines():
        if line.strip():
            index, gpu_uuid = [x.strip() for x in line.split(',', 1)]
            devices.append(dict(index=index, uuid=gpu_uuid))
    mask = os.environ.get('CUDA_VISIBLE_DEVICES')
    if mask is not None:
        allowed = [token.strip() for token in mask.split(',') if token.strip()]
        devices = [d for d in devices if any(token == d['index'] or
                   (token.startswith('GPU-') and d['uuid'].startswith(token)) for token in allowed)]
    if selection is not None:
        indices = [token.strip() for token in selection.split(',')]
        if not all(token.isdigit() for token in indices) or len(set(indices)) != len(indices):
            raise ValueError('--gpus must contain distinct GPU indices, e.g. 0,2')
        available = {d['index']: d for d in devices}
        if any(index not in available for index in indices):
            raise ValueError(f'GPU selection {selection} is unavailable; visible physical indices: {list(available)}')
        devices = [available[index] for index in indices]
    if not devices:
        raise ValueError('No visible GPUs; specify an available GPU with --gpus')
    return devices


def replica_complete(project, cfg, replica):
    base = Path(project) / f'replica_{replica:02d}'
    required = [base / 'dry_result.json']
    if cfg.get('input_mode') == 'epoxy':
        required += [base / 'crosslinking/done.json',
                     base / 'crosslinking/results/epoxy/final.structure.json']
    if cfg['mechanical']:
        required.append(base / 'mechanical/result.json')
    required += [base / 'solvents' / s['name'] / 'result.json' for s in cfg['solvents']]
    return all(path.is_file() for path in required)


ERROR_PATTERN = re.compile(r'^(?:[\w.]+(?:Error|Exception|Interrupt)\b|network_properties: ).*')


def log_section(path, offset=0, limit=2_000_000):
    """Text appended to a job log since offset (the current attempt), at most limit bytes."""
    try:
        with open(path, 'rb') as handle:
            size = handle.seek(0, os.SEEK_END)
            handle.seek(max(offset, size - limit))
            return handle.read().decode(errors='replace')
    except OSError:
        return ''


def failure_summary(text, code):
    if code < 0:
        name = signal.Signals(-code).name if -code in signal.valid_signals() else f'signal {-code}'
        hint = ' (check memory: the kernel OOM killer sends SIGKILL)' if -code == signal.SIGKILL else ''
        return f'killed by {name}{hint}'
    errors = [line.strip() for line in text.splitlines() if ERROR_PATTERN.match(line.strip())]
    return errors[-1] if errors else f'exit code {code}; no exception line in the log'


def record_failure(job, code):
    """Store a one-line cause and write jobs/replica_XX.error.log with this attempt's log tail."""
    text = log_section(job['log'], job.get('log_offset', 0))
    job['error'] = failure_summary(text, code)
    path = Path(job['log']).with_suffix('.error.log')
    tail = '\n'.join(text.splitlines()[-300:])
    path.write_text(f'[{now()}] replica={job["replica"]} seed={job["seed"]} returncode={code}\n'
                    f'cause: {job["error"]}\n--- last 300 log lines of this attempt ---\n{tail}\n')
    job['error_log'] = str(path)


def mark_interrupted(project, job):
    """A worker killed by a signal cannot close its own progress record."""
    from workflows.network_properties.execution import read_progress, finish_progress
    base = Path(project) / f'replica_{job["replica"]:02d}'
    progress = read_progress(base)
    if progress and progress.get('state') == 'running':
        finish_progress(base, 'failed', error=job['error'], error_log=job.get('error_log'))


def worker_command(project, replica):
    return [sys.executable, '-u', str(SCRIPT), 'run', str(project), '--replica', str(replica), '--defer-report']


def plan(project, cfg, devices):
    from workflows.network_properties.network_properties import input_digest
    import hashlib
    submission = dict(project=str(project), created_at=now(),
        configuration_sha256=hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest(),
        input_sha256=input_digest(project, cfg),
        gpus=devices, max_parallel=min(len(devices), cfg['replicates']),
        jobs=[dict(replica=i, seed=cfg['seed'] + i,
            state='completed' if replica_complete(project, cfg, i) else 'queued',
            log=str(Path(project) / 'jobs' / f'replica_{i:02d}.log'))
            for i in range(1, cfg['replicates'] + 1)])
    submission['max_parallel'] = min(len(devices), sum(j['state'] == 'queued' for j in submission['jobs']))
    return submission


def submit(args):
    from workflows.network_properties.network_properties import validate, input_digest, initialize_project, write_json, report
    project = args.project.resolve()
    cfg = json.loads((project / 'project.json').read_text())
    validate(cfg)
    if cfg['platform'] != 'CUDA':
        raise ValueError('GPU submission requires platform CUDA in project.json; CPU projects use run')
    if input_digest(project, cfg) != cfg['input_sha256']:
        raise ValueError('input structure was changed')
    submission = plan(project, cfg, visible_gpus(args.gpus))
    if args.dry_run:
        print(json.dumps(submission, indent=2))
        return submission
    jobs = project / 'jobs'
    jobs.mkdir(exist_ok=True)
    handle = (jobs / '.launcher.lock').open('a+')
    try:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f'GPU jobs are already running for {project}; use status') from exc
        initialize_project(project)
        if all(job['state'] == 'completed' for job in submission['jobs']):
            report(project)
            print('All replicas already completed; report refreshed.')
            return submission
        path = jobs / f'submission_{uuid.uuid4().hex}.json'
        write_json(path, submission)
        write_json(jobs / 'plan.json', submission)
        write_json(jobs / 'status.json', dict(submission, state='starting'))
        with (jobs / 'launcher.log').open('ab') as log:
            process = subprocess.Popen([sys.executable, '-u', str(Path(__file__).resolve()),
                'dispatch', str(path), str(handle.fileno())], stdin=subprocess.DEVNULL,
                stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                pass_fds=(handle.fileno(),), cwd=ROOT)
        write_json(jobs / 'launcher.json', dict(pid=process.pid, submitted_at=now(), submission=str(path)))
        print(f'Submitted {sum(j["state"] == "queued" for j in submission["jobs"])} replicas; '
              f'GPUs={",".join(d["index"] for d in submission["gpus"])}; launcher PID={process.pid}')
        print(f'Logs: {jobs}; status: {sys.executable} {SCRIPT} status {project}')
        return submission
    finally:
        # Do not LOCK_UN: the detached controller inherits this locked descriptor.
        # Closing only the parent's copy keeps duplicate submissions excluded.
        handle.close()


def dispatch(path, *, command_factory=worker_command, poll_interval=.5, refresh_report=True):
    """Run the queue and persist status; command injection is used only by CPU tests."""
    from workflows.network_properties.network_properties import write_json, report
    submission = json.loads(Path(path).read_text())
    project = Path(submission['project'])
    cfg = json.loads((project / 'project.json').read_text())
    import hashlib
    fingerprint = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()
    if fingerprint != submission['configuration_sha256']:
        write_json(project / 'jobs/status.json', dict(submission, state='failed',
            error='configuration changed after submission', finished_at=now()))
        raise ValueError('configuration changed after submission')
    status_path = project / 'jobs/status.json'
    state = dict(submission, pid=os.getpid(), state='running', started_at=now())
    pending = [job for job in state['jobs'] if job['state'] == 'queued']
    active, stopped = {}, False
    previous_handlers = {}

    def stop(signum, frame):
        nonlocal stopped
        stopped = True

    def save():
        write_json(status_path, state)

    for sig in (signal.SIGTERM, signal.SIGINT):
        previous_handlers[sig] = signal.signal(sig, stop)
    try:
        save()
        while pending or active:
            if stopped:
                break
            for device in state['gpus']:
                key = device['uuid']
                if key in active or not pending:
                    continue
                job = pending.pop(0)
                environment = dict(os.environ, CUDA_VISIBLE_DEVICES=key, CUDA_DEVICE_ORDER='PCI_BUS_ID')
                for variable in ('OPENMM_CPU_THREADS', 'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
                    environment.setdefault(variable, '1')
                log = Path(job['log'])
                log.parent.mkdir(parents=True, exist_ok=True)
                try:
                    job['log_offset'] = log.stat().st_size if log.exists() else 0
                    with log.open('ab') as output:
                        output.write(f'\n[{now()}] replica={job["replica"]} seed={job["seed"]} GPU={device["index"]}\n'.encode())
                        output.flush()
                        child = subprocess.Popen(command_factory(project, job['replica']),
                            env=environment, cwd=ROOT, stdin=subprocess.DEVNULL,
                            stdout=output, stderr=subprocess.STDOUT, start_new_session=True)
                    job.update(state='running', pid=child.pid, gpu_index=device['index'],
                               gpu_uuid=key, started_at=now())
                    active[key] = (child, job)
                except OSError as exc:
                    job.update(state='failed', error=str(exc), finished_at=now())
                save()
            for key, (child, job) in list(active.items()):
                code = child.poll()
                if code is not None:
                    job.update(state='completed' if code == 0 else 'failed',
                               returncode=code, finished_at=now())
                    if code:
                        try:
                            record_failure(job, code)
                            mark_interrupted(project, job)
                        except OSError as exc:
                            job['error'] = f'exit code {code}; could not write error log: {exc}'
                    del active[key]
                    save()
            if active:
                time.sleep(poll_interval)
        if stopped:
            for child, job in active.values():
                if child.poll() is None:
                    os.killpg(child.pid, signal.SIGTERM)
            for child, job in active.values():
                try:
                    child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait()
                job.update(state='cancelled', returncode=child.returncode, finished_at=now())
            for job in pending:
                job.update(state='cancelled')
        terminal_state = 'cancelled' if stopped else (
            'failed' if any(j['state'] == 'failed' for j in state['jobs']) else 'completed')
        state['state'] = 'collecting' if refresh_report else terminal_state
        save()
        if refresh_report:
            try:
                report(project)
            except Exception as exc:
                state.update(state='failed', report_error=str(exc))
                save()
                raise
        state.update(state=terminal_state, finished_at=now())
        save()
        return state
    except Exception as exc:
        for child, job in active.values():
            if child.poll() is None:
                job.update(state='cancelled')
        for job in pending:
            job.update(state='cancelled')
        state.update(state='failed', error=str(exc), finished_at=now())
        save()
        raise
    finally:
        # Unexpected launcher errors must not leave its own workers unattended.
        for child, job in active.values():
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGTERM)
        for sig, handler in previous_handlers.items():
            signal.signal(sig, handler)


def status(args):
    project = args.project.resolve()
    path = project / 'jobs/status.json'
    if not path.is_file():
        path = project / 'jobs/plan.json'
    if not path.is_file():
        raise FileNotFoundError(f'No GPU submission found: {project}')
    state = json.loads(path.read_text())
    from workflows.network_properties.execution import read_progress
    print(f'Project: {project}\nState: {state.get("state", "starting")}')
    if state.get('state') in ('starting', 'running', 'collecting'):
        launcher = process_state(state.get('pid'))
        if launcher in (None, 'Z'):
            print(f'WARNING: launcher PID {state.get("pid")} is '
                  f'{"a zombie" if launcher == "Z" else "not running"}; job states below are not '
                  'being updated. Check each replica, then submit again to resume unfinished replicas.')
    for job in state['jobs']:
        name = f'replica_{job["replica"]:02d}'
        print(f'{name}: {job["state"]:10s} seed={job["seed"]} GPU={job.get("gpu_index", "—")} log={job["log"]}')
        if job['state'] == 'running':
            alive = process_state(job.get('pid'))
            if alive is None:
                print(f'  WARNING: worker PID {job.get("pid")} no longer exists (lost)')
            elif alive == 'Z':
                print(f'  WARNING: worker PID {job.get("pid")} exited but was not reaped (zombie)')
        progress = read_progress(project / name)
        if progress:
            stage = progress.get('stage', '—')
            since = progress.get('finished_at') if progress.get('state') != 'running' else None
            print(f'  stage: {stage} [{progress.get("state")}] since {progress.get("started_at")} '
                  f'({elapsed(progress.get("started_at"), since)})')
            if progress.get('error'):
                print(f'  error: {progress["error"]}\n  error log: {progress.get("error_log")}')
        text = log_section(job['log'], job.get('log_offset', 0), limit=200_000)
        cycles = CYCLE_PATTERN.findall(text)
        crosslinking = (progress or {}).get('stage', '').startswith('crosslinking')
        if cycles and crosslinking and job['state'] == 'running':
            cycle, conversion = cycles[-1]
            print(f'  crosslinking: cycle {cycle}, conversion {conversion}')
        last = next((line.strip() for line in reversed(text.splitlines()) if line.strip()), None)
        if last and job['state'] == 'running':
            print(f'  last log: {last[:200]}')
        if job.get('error'):
            print(f'  cause: {job["error"]}\n  job error log: {job.get("error_log")}')
    return state


CYCLE_PATTERN = re.compile(r'reaction cycle completed \| cycle=(\d+) .*?conversion=([\d.]+%)')


def process_state(pid):
    """Linux process state letter ('R', 'S', 'Z', ...), or None when the PID is gone."""
    if not pid:
        return None
    try:
        return Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()[0]
    except (OSError, IndexError):
        return None


def elapsed(start, end=None):
    try:
        begin = datetime.fromisoformat(start)
        finish = datetime.fromisoformat(end) if end else datetime.now(timezone.utc)
    except (TypeError, ValueError):
        return '—'
    seconds = int((finish - begin).total_seconds())
    return f'{seconds // 3600}h{seconds % 3600 // 60:02d}m{seconds % 60:02d}s'


if __name__ == '__main__':
    if len(sys.argv) != 4 or sys.argv[1] != 'dispatch':
        raise SystemExit('This helper is started by network_properties.py submit.')
    # Keep the descriptor inherited from submit open for the entire queue lifetime.
    with os.fdopen(int(sys.argv[3]), 'a+'):
        outcome = dispatch(Path(sys.argv[2]))
    raise SystemExit(0 if outcome['state'] == 'completed' else 1)
