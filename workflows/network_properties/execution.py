"""Process locks shared by sequential runs, GPU workers and report generation."""
from contextlib import contextmanager
import ctypes
import fcntl
import os
from pathlib import Path
import signal

PR_SET_PDEATHSIG = 1


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
    ctypes.CDLL(None, use_errno=True).prctl(PR_SET_PDEATHSIG, signal.SIGTERM)
    if os.getppid() != parent:
        os._exit(1)
