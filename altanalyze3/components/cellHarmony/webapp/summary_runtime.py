"""Request-scoped summary progress, timings and bounded completion reuse."""
from contextlib import contextmanager
from contextvars import ContextVar
import hashlib
import json
import logging
from pathlib import Path
import tempfile
import time
import urllib.request

_REPORT = ContextVar('summary_report', default=None)
_TIMINGS = ContextVar('summary_timings', default=None)
_CACHE = ContextVar('summary_cache', default=None)
CACHE_SECONDS = 3600
CACHE_ENTRIES = 64
MAX_RESPONSE = 2 * 1024 * 1024


@contextmanager
def reporting(callback, cache_dir=None):
    report_token = _REPORT.set(callback)
    cache_token = _CACHE.set(Path(cache_dir) if cache_dir is not None else None)
    try:
        yield
    finally:
        _CACHE.reset(cache_token)
        _REPORT.reset(report_token)


@contextmanager
def trace():
    records = []
    token = _TIMINGS.set(records)
    try:
        yield records
    finally:
        _TIMINGS.reset(token)


@contextmanager
def phase(stage, detail):
    callback = _REPORT.get()
    if callback:
        callback(stage=stage, detail=detail)
    record = {'stage': stage}
    start = time.monotonic()
    try:
        yield record
    finally:
        record['seconds'] = round(time.monotonic() - start, 4)
        if _TIMINGS.get() is not None:
            _TIMINGS.get().append(record)
        logging.info('AI summary stage=%s seconds=%.3f cache_hit=%s',
                     stage, record['seconds'], record.get('cache_hit', False))


def completion(request, stage, detail):
    """Reuse exact prompts for one job; never persist authentication headers."""
    directory = _CACHE.get()
    # URL, model, full source-bound messages, schema and generation parameters
    # must all match. No evidence is omitted to improve cache hit rates.
    key = hashlib.sha256(request.full_url.encode() + b'\0' + request.data).hexdigest()
    path = directory / (key + '.json') if directory is not None else None
    with phase(stage, detail) as timing:
        if path is not None and path.is_file() and time.time() - path.stat().st_mtime < CACHE_SECONDS:
            raw = path.read_bytes()
            if len(raw) <= MAX_RESPONSE:
                try:
                    result = json.loads(raw)
                except (ValueError, UnicodeDecodeError):
                    pass
                else:
                    timing['cache_hit'] = True
                    timing['usage'] = result.get('usage', {})
                    timing['runtime'] = result.get('scalable_runtime', {})
                    return result
        with urllib.request.urlopen(request, timeout=240) as response:
            raw = response.read(MAX_RESPONSE + 1)
        if len(raw) > MAX_RESPONSE:
            raise ValueError('Summary model response exceeds size limit')
        result = json.loads(raw)
        timing['cache_hit'] = False
        timing['usage'] = result.get('usage', {})
        timing['runtime'] = result.get('scalable_runtime', {})
        if path is not None:
            try:
                directory.mkdir(parents=True, exist_ok=True)
                for old in sorted(directory.glob('*.json'), key=lambda p: p.stat().st_mtime)[:-(CACHE_ENTRIES - 1)]:
                    old.unlink(missing_ok=True)
                with tempfile.NamedTemporaryFile(dir=directory, prefix='.', suffix='.tmp', delete=False) as f:
                    f.write(raw)
                    temporary = Path(f.name)
                temporary.replace(path)
            except OSError:
                logging.warning('AI summary completion cache unavailable; retaining the generated result.')
        return result
