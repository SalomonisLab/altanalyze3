"""Bounded pipeline log snapshots, reused until the file changes."""
from collections import deque
from pathlib import Path
import re
import threading

from .memory_cache import BoundedCache, CacheBudget

_CACHE = BoundedCache(CacheBudget(max_bytes=1024**2, max_entries=64, ttl_seconds=600))
_LOCK = threading.RLock()
_AMBIENT = re.compile(r"Auto-selected rho for library '(.+?)':")


def read_pipeline_log(path):
    path = Path(path)
    with _LOCK:
        try:
            stat = path.stat()
        except FileNotFoundError:
            _CACHE.pop(str(path), None)
            return [], [], []
        signature = (stat.st_ino, stat.st_size, stat.st_mtime_ns)
        cached = _CACHE.get(str(path))
        if cached is not None and cached['signature'] == signature:
            return cached['head'], cached['tail'], cached['progress']
        head, tail, progress = [], deque(maxlen=200), {}
        with path.open(encoding='utf-8') as stream:
            for line in stream:
                if len(head) < 80:
                    head.append(line)
                tail.append(line)
                key = None
                for marker in ('adata shape:', 'Cells remaining after min_genes',
                               'Cells remaining after min_counts', 'Cells remaining after mito-percent',
                               'Applied min_alignment_score='):
                    if marker in line:
                        key = marker
                        break
                if 'Auto-selected rho for library' in line:
                    # Preserve each library's latest correction, not all repetitions.
                    match = _AMBIENT.search(line)
                    if match:
                        key = 'ambient:' + match.group(1)
                if key is not None:
                    progress.pop(key, None)
                    progress[key] = line
        snapshot = dict(signature=signature, head=head, tail=list(tail), progress=list(progress.values()))
        _CACHE[str(path)] = snapshot
        return snapshot['head'], snapshot['tail'], snapshot['progress']
