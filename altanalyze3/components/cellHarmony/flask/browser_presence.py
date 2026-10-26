"""Job-scoped browser leases; a fresh page without a job never owns other jobs."""
import json
import math
import os
import re
import tempfile
import time


def record_presence(store, job_id, client_id, state='active', now=None, sequence=0):
    if not re.fullmatch(r'[A-Za-z0-9_-]{12,80}', client_id):
        raise ValueError('Invalid browser client ID.')
    if state not in {'active', 'leaving'}:
        raise ValueError('Invalid browser presence state.')
    if not isinstance(sequence, int) or isinstance(sequence, bool) or not 0 <= sequence <= 2**53:
        raise ValueError('Invalid browser presence sequence.')
    now = time.time() if now is None else now
    seconds = float(os.environ.get('CELLHARMONY_BROWSER_CLOSE_GRACE_SECONDS', '20') if state == 'leaving'
                    else os.environ.get('CELLHARMONY_BROWSER_LEASE_SECONDS', '120'))
    if not math.isfinite(seconds) or seconds <= 0:
        raise ValueError('Browser lease durations must be positive and finite.')
    directory = store.logs_dir(job_id).parent / '.browser_clients'
    directory.mkdir(exist_ok=True)
    path = directory / client_id
    if path.exists():
        previous = json.loads(path.read_text())
        if previous.get('sequence', -1) > sequence:
            return
    temporary = None
    try:
        with tempfile.NamedTemporaryFile('w', dir=directory, delete=False) as handle:
            temporary = handle.name
            json.dump({'expires': now + seconds, 'state': state, 'sequence': sequence}, handle)
        os.replace(temporary, path)
    finally:
        if temporary and os.path.exists(temporary):
            os.unlink(temporary)


def browser_disconnected(store, job_id, now=None):
    directory = store.logs_dir(job_id).parent / '.browser_clients'
    if not directory.exists():
        return False  # CLI/API-only jobs have no browser lease.
    now = time.time() if now is None else now
    records = [p for p in directory.iterdir() if re.fullmatch(r'[A-Za-z0-9_-]{12,80}', p.name)]
    if not records:
        return False
    for path in records:
        try:
            expires = float(json.loads(path.read_text())['expires'])
            if not math.isfinite(expires) or expires > now:
                return False
        except (OSError, ValueError, KeyError, TypeError):
            return False  # Incomplete/missing telemetry never proves abandonment.
    return True
