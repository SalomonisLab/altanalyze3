"""Cache small group menus, retaining no expression or per-cell annotations."""
import copy
from pathlib import Path
import threading

from .memory_cache import BoundedCache, CacheBudget
from ..flask import pipeline

_CACHE = BoundedCache(CacheBudget(max_bytes=1024**2, max_entries=64, ttl_seconds=600))
_LOCK = threading.RLock()


def cached_group_fields(path, preferred=None, max_categories=None):
    path = Path(path)
    try:
        stat = path.stat()
    except FileNotFoundError:
        return [], {}
    key = (str(path.resolve()), stat.st_ino, stat.st_size, stat.st_mtime_ns,
           tuple(preferred or ()), max_categories)
    with _LOCK:
        result = _CACHE.get(key)
        if result is None:
            result = pipeline._candidate_group_fields(path, preferred=preferred, max_categories=max_categories)
            _CACHE[key] = result
        # Callers can change their menus without mutating another visitor's cache.
        return copy.deepcopy(result)
