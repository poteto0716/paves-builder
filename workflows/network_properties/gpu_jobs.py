"""Detached local GPU queue: one replica per selected GPU, explicit seeds and logs.

Submitting a monomer project also queues its crosslink projects; crosslink
replica N starts when monomer replica N has completed. --mechanical queues
the optional elastic constants of finished crosslink/network replicas.
"""
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
from workflows.network_properties.execution import configuration_fingerprint, die_with_parent, file_lock
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
    """Every default stage of the replica has finished (the optional mechanical stage is not one)."""
    mode = cfg.get('input_mode', 'saved_network')
    if mode == 'monomer':
        from workflows.network_properties.monomer import replica_complete as monomer_complete
        return monomer_complete(project, cfg, replica)
    base = Path(project) / f'replica_{replica:02d}'
    required = [base / 'dry_result.json']
    if mode == 'crosslink':
        from workflows.network_properties.construction import network_path
        required += [base / 'crosslinking/done.json', network_path(project, cfg, replica)]
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
    base = Path(job.get('project', project)) / f'replica_{job["replica"]:02d}'
    try:
        with file_lock(base / '.run.lock', blocking=False):
            pass
    except RuntimeError:
        return  # Refused by the lock: another live worker owns this progress record.
    progress = read_progress(base)
    if progress and progress.get('state') == 'running':
        finish_progress(base, 'failed', error=job['error'], error_log=job.get('error_log'))


def live_replicas(project, cfg):
    """Replicas whose run lock is held by a process, whoever started it."""
    live = []
    for replica in range(1, cfg['replicates'] + 1):
        try:
            with file_lock(Path(project) / f'replica_{replica:02d}/.run.lock', blocking=False):
                pass
        except RuntimeError:
            live.append(replica)
    return live


def worker_command(project, replica):
    return [sys.executable, '-u', str(SCRIPT), 'run', str(project), '--replica', str(replica), '--defer-report']


def mechanical_command(project, replica):
    return [sys.executable, '-u', str(SCRIPT), 'mechanical', str(project), '--replica', str(replica), '--defer-report']


def submission_projects(project, cfg, monomer_only=False):
    """The project and, for a monomer project, its crosslink projects: [(directory, cfg)]."""
    from workflows.network_properties.network_properties import children
    projects = [(Path(project), cfg)]
    if cfg.get('input_mode') == 'monomer' and not monomer_only:
        projects += [(child, json.loads((child / 'project.json').read_text())) for child in children(project)]
    return projects


def plan(project, cfg, devices, mechanical=False, monomer_only=False):
    from workflows.network_properties.network_properties import input_digest
    from workflows.network_properties.execution import replica_seed
    from workflows.network_properties.mechanical import result_current
    project = Path(project)
    jobs, configurations = [], {}
    if monomer_only and (mechanical or cfg.get('input_mode') != 'monomer'):
        raise ValueError('--monomer-only submits the replicas of a monomer project (without --mechanical)')
    for directory, settings in submission_projects(project, cfg, monomer_only):
        configurations[str(directory)] = configuration_fingerprint(settings)
        mode = settings.get('input_mode', 'saved_network')
        if mechanical and mode == 'monomer':
            continue
        name = '.' if directory == project else directory.relative_to(project).as_posix()
        for i in range(1, settings['replicates'] + 1):
            base = directory / f'replica_{i:02d}'
            job = dict(project=str(directory), name=name, replica=i, seed=replica_seed(directory, settings, i),
                       kind='mechanical' if mechanical else 'run',
                       log=str(directory / 'jobs' / f'replica_{i:02d}{".mechanical" if mechanical else ""}.log'))
            if mechanical:
                ready = (base / 'relaxed_300K.structure.json').is_file()
                job['state'] = 'completed' if ready and result_current(base) else 'queued' if ready else 'skipped'
                if not ready:
                    job['error'] = 'no relaxed 300 K structure yet; run the project first'
            else:
                job['state'] = 'completed' if replica_complete(directory, settings, i) else 'queued'
                if mode == 'crosslink':
                    from workflows.network_properties.construction import parent_project
                    job['after'] = dict(project=str(parent_project(directory, settings)), replica=i)
            jobs.append(job)
    submission = dict(project=str(project), created_at=now(),
        configuration_sha256=configuration_fingerprint(cfg), configurations=configurations,
        input_sha256=input_digest(project, cfg), gpus=devices, jobs=jobs)
    submission['max_parallel'] = min(len(devices), sum(j['state'] == 'queued' for j in jobs))
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
    monomer_only = getattr(args, 'monomer_only', False)
    submission = plan(project, cfg, visible_gpus(args.gpus), mechanical=getattr(args, 'mechanical', False),
                      monomer_only=monomer_only)
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
        # Workers stop when the launcher dies, but not instantly; never share a GPU with one.
        for directory, settings in submission_projects(project, cfg, monomer_only):
            live = live_replicas(directory, settings)
            if live:
                raise RuntimeError(f'Replicas {live} of {directory} are still running outside this launcher; '
                                   'wait for them to stop or stop them, then submit again')
            initialize_project(directory)
        if all(job['state'] in ('completed', 'skipped') for job in submission['jobs']):
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
        print(f'Submitted {sum(j["state"] == "queued" for j in submission["jobs"])} replica jobs; '
              f'GPUs={",".join(d["index"] for d in submission["gpus"])}; launcher PID={process.pid}')
        print(f'Logs: {jobs}; status: {sys.executable} {SCRIPT} status {project}')
        return submission
    finally:
        # Do not LOCK_UN: the detached controller inherits this locked descriptor.
        # Closing only the parent's copy keeps duplicate submissions excluded.
        handle.close()


def _dependency(job, jobs):
    """'ready', 'wait' or the reason the job cannot run (its monomer replica did not complete)."""
    after = job.get('after')
    if not after:
        return 'ready'
    parent = Path(after['project'])
    dependency = next((j for j in jobs if Path(j['project']) == parent and j['replica'] == after['replica']
                       and j.get('kind', 'run') == 'run'), None)
    label = f'monomer {parent.name}/replica_{after["replica"]:02d}'
    if dependency is not None and dependency['state'] in ('queued', 'running'):
        return 'wait'
    if dependency is not None and dependency['state'] != 'completed':
        return f'{label} {dependency["state"]}'
    cfg = json.loads((parent / 'project.json').read_text())
    return 'ready' if replica_complete(parent, cfg, after['replica']) else f'{label} has not finished; submit the monomer project'


def dispatch(path, *, command_factory=worker_command, mechanical_factory=mechanical_command,
             poll_interval=.5, refresh_report=True):
    """Run the queue and persist status; command injection is used only by CPU tests."""
    from workflows.network_properties.network_properties import write_json, report
    submission = json.loads(Path(path).read_text())
    project = Path(submission['project'])
    status_path = project / 'jobs/status.json'
    configurations = submission.get('configurations') or {str(project): submission['configuration_sha256']}
    for directory, fingerprint in configurations.items():
        cfg = json.loads((Path(directory) / 'project.json').read_text())
        if configuration_fingerprint(cfg) != fingerprint:
            write_json(status_path, dict(submission, state='failed',
                error=f'configuration changed after submission: {directory}', finished_at=now()))
            raise ValueError(f'configuration changed after submission: {directory}')
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
            for job in list(pending):
                verdict = _dependency(job, state['jobs'])
                if verdict not in ('ready', 'wait'):
                    job.update(state='blocked', error=verdict, finished_at=now())
                    pending.remove(job)
                    save()
            for device in state['gpus']:
                key = device['uuid']
                ready = [job for job in pending if _dependency(job, state['jobs']) == 'ready']
                if key in active or not ready:
                    continue
                job = ready[0]
                pending.remove(job)
                environment = dict(os.environ, CUDA_VISIBLE_DEVICES=key, CUDA_DEVICE_ORDER='PCI_BUS_ID')
                for variable in ('OPENMM_CPU_THREADS', 'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
                    environment.setdefault(variable, '1')
                log = Path(job['log'])
                log.parent.mkdir(parents=True, exist_ok=True)
                try:
                    job['log_offset'] = log.stat().st_size if log.exists() else 0
                    with log.open('ab') as output:
                        output.write(f'\n[{now()}] {job.get("kind", "run")} replica={job["replica"]} seed={job["seed"]} '
                                     f'GPU={device["index"]}\n'.encode())
                        output.flush()
                        factory = mechanical_factory if job.get('kind') == 'mechanical' else command_factory
                        child = subprocess.Popen(factory(Path(job.get('project', project)), job['replica']),
                            env=environment, cwd=ROOT, stdin=subprocess.DEVNULL,
                            stdout=output, stderr=subprocess.STDOUT, start_new_session=True,
                            # A killed launcher must not leave a worker on a GPU it no longer tracks.
                            preexec_fn=die_with_parent)
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
            if active or pending:
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
            'failed' if any(j['state'] in ('failed', 'blocked') for j in state['jobs']) else 'completed')
        state['state'] = 'collecting' if refresh_report else terminal_state
        save()
        if refresh_report:
            try:
                # Crosslink reports also refresh the monomer report and its series.
                children = [Path(d) for d in configurations if Path(d) != project]
                for directory in children:
                    report(directory)
                if not children:
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


def next_seed(cfg, used, replica):
    """First seed + k*replicates + replica (k >= 1) not yet used in this project.

    The initial seeds are seed + 1..replicates, so successive attempts of
    different replicas never collide with each other or with first attempts.
    """
    k = 1
    while (candidate := cfg['seed'] + k * cfg['replicates'] + replica) in used:
        k += 1
    return candidate


def reseed(args):
    """Give failed replicas a new seed and archive their previous attempt.

    The archived directory keeps every output and error log of the old seed.
    The replica then restarts from scratch (in construction mode, from packing),
    because every stage of a replica derives from its seed.
    """
    import shutil
    from workflows.network_properties.execution import (read_progress, replica_seeds, seed_overrides,
                                                        write_seed_overrides)
    project = args.project.resolve()
    cfg = json.loads((project / 'project.json').read_text())
    replicas = sorted(set(args.replica))
    if any(not 1 <= r <= cfg['replicates'] for r in replicas):
        raise ValueError(f'--replica must lie between 1 and {cfg["replicates"]}')
    if args.seed is not None and len(replicas) != 1:
        raise ValueError('--seed sets one replica; omit it to choose new seeds automatically')
    jobs = project / 'jobs'
    jobs.mkdir(exist_ok=True)
    # Exclude a concurrent submit: its launcher holds this lock while replicas run.
    with file_lock(jobs / '.launcher.lock', blocking=False), file_lock(project / '.project.lock'):
        live = set(live_replicas(project, cfg)) & set(replicas)
        if live:
            raise RuntimeError(f'Replicas {sorted(live)} are running; stop them before reseeding')
        overrides = seed_overrides(project)
        current = replica_seeds(project, cfg)
        used = set(current.values()) | {s for o in overrides.values() for s in o.get('previous', [])}
        if args.seed is not None and args.seed in used:
            raise ValueError(f'seed {args.seed} was already used in this project')
        plan = []
        for replica in replicas:
            base = project / f'replica_{replica:02d}'
            progress = read_progress(base) or {}
            if replica_complete(project, cfg, replica) and not args.force:
                raise ValueError(f'replica {replica} has completed; use --force to recompute it with a new seed')
            if progress.get('state') != 'failed' and not args.force:
                raise ValueError(f'replica {replica} has not failed (state: {progress.get("state", "not started")}); '
                                 'use --force to reseed it anyway')
            seed = args.seed if args.seed is not None else next_seed(cfg, used, replica)
            used.add(seed)
            plan.append((replica, base, current[replica], seed, progress))
        for replica, base, old, seed, progress in plan:
            print(f'replica_{replica:02d}: seed {old} -> {seed}'
                  + (f' (failed during: {progress.get("stage")})' if progress.get('state') == 'failed' else ''))
        if args.dry_run:
            return plan
        for replica, base, old, seed, progress in plan:
            archive = None
            if base.exists():
                index = 1
                while (archive := project / 'failed_attempts' / f'replica_{replica:02d}_seed{old}_{index}').exists():
                    index += 1
                archive.parent.mkdir(exist_ok=True)
                shutil.move(base, archive)
                for suffix in ('.error.log',):
                    log = jobs / f'replica_{replica:02d}{suffix}'
                    if log.exists():
                        shutil.move(log, archive / f'job{suffix}')
            entry = overrides.get(str(replica), dict(previous=[]))
            entry.update(seed=seed, previous=entry['previous'] + [old], reseeded_at=now(),
                         previous_error=progress.get('error'))
            entry.setdefault('archives', []).append(None if archive is None else str(archive.relative_to(project)))
            overrides[str(replica)] = entry
            write_seed_overrides(project, overrides)
            print(f'  previous attempt archived: {archive}' if archive else '  (no previous output)')
    if args.submit:
        return submit(args)
    print(f'Run them with: {sys.executable} {SCRIPT} submit {project} [--gpus ...]')
    return plan


def _print_replica(directory, replica, job=None):
    from workflows.network_properties.execution import read_progress
    name = f'replica_{replica:02d}'
    progress = read_progress(Path(directory) / name)
    if job is not None:
        kind = '' if job.get('kind', 'run') == 'run' else f' [{job["kind"]}]'
        print(f'  {name}{kind}: {job["state"]:10s} seed={job["seed"]} GPU={job.get("gpu_index", "—")} log={job["log"]}')
        if job['state'] == 'running':
            alive = process_state(job.get('pid'))
            if alive is None:
                print(f'    WARNING: worker PID {job.get("pid")} no longer exists (lost)')
            elif alive == 'Z':
                print(f'    WARNING: worker PID {job.get("pid")} exited but was not reaped (zombie)')
    else:
        print(f'  {name}: {(progress or {}).get("state", "not started")}')
    if progress:
        stage = progress.get('stage', '—')
        since = progress.get('finished_at') if progress.get('state') != 'running' else None
        print(f'    stage: {stage} [{progress.get("state")}] since {progress.get("started_at")} '
              f'({elapsed(progress.get("started_at"), since)})')
        if progress.get('error'):
            print(f'    error: {progress["error"]}\n    error log: {progress.get("error_log")}')
    if job is None:
        return
    text = log_section(job['log'], job.get('log_offset', 0), limit=200_000)
    cycles = CYCLE_PATTERN.findall(text)
    crosslinking = (progress or {}).get('stage', '').startswith('crosslinking')
    if cycles and crosslinking and job['state'] == 'running':
        cycle, conversion = cycles[-1]
        print(f'    crosslinking: cycle {cycle}, conversion {conversion}')
    last = next((line.strip() for line in reversed(text.splitlines()) if line.strip()), None)
    if last and job['state'] == 'running':
        print(f'    last log: {last[:200]}')
    if job.get('error'):
        print(f'    cause: {job["error"]}')
        if job.get('error_log'):
            print(f'    job error log: {job["error_log"]}')


def status(args):
    """Replica stages of the project (a monomer project: with its crosslinks) and the latest GPU submission."""
    project = args.project.resolve()
    cfg = json.loads((project / 'project.json').read_text())
    path = project / 'jobs/status.json'
    if not path.is_file():
        path = project / 'jobs/plan.json'
    state = json.loads(path.read_text()) if path.is_file() else None
    print(f'Project: {project}')
    if state is None:
        print('GPU submission: none')
    else:
        print(f'GPU submission: {state.get("state", "starting")}')
        if state.get('state') in ('starting', 'running', 'collecting'):
            launcher = process_state(state.get('pid'))
            if launcher in (None, 'Z'):
                print(f'WARNING: launcher PID {state.get("pid")} is '
                      f'{"a zombie" if launcher == "Z" else "not running"}; job states below are not '
                      'being updated. Its workers stop with it; submit again to resume unfinished replicas '
                      '(submit refuses while one is still running).')
    jobs = (state or {}).get('jobs', [])
    for directory, settings in submission_projects(project, cfg):
        mode = settings.get('input_mode', 'saved_network')
        title = directory.name if directory == project else directory.relative_to(project).as_posix()
        done = sum(replica_complete(directory, settings, i) for i in range(1, settings['replicates'] + 1))
        print(f'{title} [{mode}]: {done}/{settings["replicates"]} replicas complete')
        for replica in range(1, settings['replicates'] + 1):
            own = [j for j in jobs if Path(j.get('project', project)) == directory and j['replica'] == replica]
            if own:
                for job in own:
                    _print_replica(directory, replica, job)
            else:
                _print_replica(directory, replica)
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
