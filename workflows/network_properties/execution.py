"""Process locks, child-process lifetime and per-replica progress/error records."""
from contextlib import contextmanager
import ctypes
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import signal
import traceback

PR_SET_PDEATHSIG = 1
# Resolved before fork: the child must not dlopen between fork and exec.
_prctl = ctypes.CDLL(None, use_errno=True).prctl
PROGRESS = 'progress.json'
ERROR_LOG = 'error.log'


@contextmanager
def file_lock(path, *, blocking=True):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a+') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB))
        except BlockingIOError as exc:
            raise RuntimeError(f'Already running: {path.parent}') from exc
        try:
            yield handle
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def die_with_parent():
    """preexec_fn: SIGTERM this child when the starting process dies.

    Without it a killed worker leaves its crosslinking child running, and a
    later run would start a second build in the same directory.
    """
    parent = os.getppid()
    _prctl(PR_SET_PDEATHSIG, signal.SIGTERM)
    if os.getppid() != parent:
        os._exit(1)


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def _write(path, value):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def read_progress(base):
    try:
        return json.loads((Path(base) / PROGRESS).read_text())
    except (OSError, ValueError):
        return None


def set_stage(base, stage, **detail):
    """Record the replica's current stage in progress.json and print it to the job log.

    The previous stage is closed and kept in a short history, so status can show
    where a replica is and where it stopped.
    """
    base = Path(base)
    base.mkdir(parents=True, exist_ok=True)
    progress = read_progress(base) or {}
    history = progress.get('history', [])
    if progress.get('state') == 'running' and progress.get('stage'):
        history.append(dict(stage=progress['stage'], started_at=progress['started_at'], finished_at=now()))
    progress = dict(state='running', stage=stage, detail=detail, started_at=now(), pid=os.getpid(),
                    history=history[-50:])
    _write(base / PROGRESS, progress)
    extra = ' '.join(f'{k}={v}' for k, v in detail.items())
    print(f'[stage {progress["started_at"]}] {stage}' + (f' | {extra}' if extra else ''), flush=True)


def finish_progress(base, state, **fields):
    base = Path(base)
    progress = read_progress(base) or {}
    finished = now()
    if progress.get('state') == 'running' and progress.get('stage'):
        progress.setdefault('history', []).append(dict(stage=progress['stage'], started_at=progress['started_at'],
                                                       finished_at=finished, state=state))
    progress.update(state=state, finished_at=finished, **fields)
    _write(base / PROGRESS, progress)


@contextmanager
def replica_progress(base):
    """Mark the replica completed, or write error.log with the stage and traceback."""
    base = Path(base)
    try:
        yield
    except BaseException as exc:
        progress = read_progress(base) or {}
        stage = progress.get('stage', 'unknown')
        text = (f'[{now()}] replica failed during stage: {stage}\n'
                f'pid={os.getpid()} detail={json.dumps(progress.get("detail", {}))}\n\n'
                + ''.join(traceback.format_exception(exc)))
        (base / ERROR_LOG).write_text(text)
        finish_progress(base, 'failed', error=f'{type(exc).__name__}: {exc}'.strip(),
                        error_log=str(base / ERROR_LOG))
        print(f'[{now()}] replica failed during stage: {stage}; traceback: {base / ERROR_LOG}', flush=True)
        raise
    else:
        finish_progress(base, 'completed', stage='completed')
