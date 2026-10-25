"""Process locks shared by sequential runs, GPU workers and report generation."""
from contextlib import contextmanager
import fcntl
from pathlib import Path


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
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)
