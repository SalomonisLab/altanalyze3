"""Filesystem leases shared by web processes; no scientific workflow changes."""
from contextlib import contextmanager
import fcntl
import math
from pathlib import Path
import time


class JobBusyError(RuntimeError):
    pass


class JobCancelled(RuntimeError):
    pass


class JobLease:
    def __init__(self, path):
        self.handle = Path(path).open('a+')
        try:
            fcntl.flock(self.handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            self.handle.close()
            raise JobBusyError('This session already has an active analysis. Stop it or wait for it to finish.') from None

    def close(self):
        # Closing, rather than LOCK_UN, preserves the inherited child's lease
        # if its supervisor exits. The last descriptor releases the lock.
        self.handle.close()


@contextmanager
def control_job(store, job_id):
    lease = JobLease(store.logs_dir(job_id).parent / '.control.lock')
    try:
        yield
    finally:
        lease.close()


def execution_lease(store, job_id):
    return JobLease(store.logs_dir(job_id).parent / '.execution.lock')


def execution_busy(store, job_id):
    try:
        lease = execution_lease(store, job_id)
    except JobBusyError:
        return True
    lease.close()
    return False


def worker_activity(processes, root):
    """Linux child-tree CPU ticks and logical I/O, including cached file reads."""
    if not Path('/proc/self/stat').exists():
        return None
    children = {}
    for pid, (parent, _) in processes.items():
        children.setdefault(parent, []).append(pid)
    pending, seen, cpu, io = [root], set(), 0, 0
    while pending:
        pid = pending.pop()
        if pid in seen:
            continue
        seen.add(pid)
        pending.extend(children.get(pid, ()))
        try:
            stat = Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()
            cpu += int(stat[11]) + int(stat[12])
            values = dict(line.split(':', 1) for line in Path(f'/proc/{pid}/io').read_text().splitlines())
            io += int(values['rchar']) + int(values['wchar'])
        except (OSError, ValueError, KeyError, IndexError):
            return None  # Missing telemetry is not evidence of a stalled worker.
    return cpu, io


class WorkerWatchdog:
    def __init__(self, idle_seconds=1800, max_seconds=21600, clock=time.monotonic):
        self.clock = clock
        self.started = self.last_progress = clock()
        self.previous = None
        self.idle_seconds = float(idle_seconds)
        self.max_seconds = float(max_seconds)
        if any(not math.isfinite(v) or v < 0 for v in (self.idle_seconds, self.max_seconds)):
            raise ValueError('Worker watchdog time limits must be finite and nonnegative.')

    def check(self, activity, log_signature):
        now = self.clock()
        signature = (activity, log_signature)
        if self.previous != signature:
            self.last_progress = now
            self.previous = signature
        if self.max_seconds > 0 and now - self.started >= self.max_seconds:
            return f'Analysis exceeded the configured {self.max_seconds / 3600:g}-hour runtime limit.'
        if activity is not None and self.idle_seconds > 0 and now - self.last_progress >= self.idle_seconds:
            return (f'Analysis stopped after {self.idle_seconds / 60:g} minutes without worker CPU, '
                    'I/O or log activity. Uploaded data are retained; the analysis can be rerun.')
        return None
