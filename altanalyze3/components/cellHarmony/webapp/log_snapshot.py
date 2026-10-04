"""Bounded pipeline log snapshots, reused until the file changes."""
from collections import deque
from datetime import datetime, timezone
import math
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
                if line.rstrip().endswith('] Job accepted by worker.'):
                    key = 'analysis_start'
                elif line.rstrip().endswith('] Job completed.'):
                    key = 'analysis_end'
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


def analysis_duration_seconds(meta, progress_lines):
    """Main worker duration, independent of uploads, queue waits and later analyses."""
    if meta.get("status") != "completed":
        return None
    saved = meta.get("analysis_duration_seconds")
    if isinstance(saved, (int, float)) and math.isfinite(saved) and saved >= 0:
        return saved
    # Older jobs have no dedicated timing fields. The bounded log snapshot keeps
    # the latest main-run markers even after many differential/interaction logs.
    start = end = None
    for line in progress_lines:
        match = re.match(r"^\[([^]]+)\] (Job accepted by worker\.|Job completed\.)\s*$", line)
        if not match:
            continue
        try:
            stamp = datetime.fromisoformat(match[1].replace("Z", "+00:00"))
        except ValueError:
            continue
        if stamp.tzinfo is None:
            stamp = stamp.replace(tzinfo=timezone.utc)
        if match[2] == "Job accepted by worker.":
            start, end = stamp, None
        else:
            end = stamp
    if start is not None and end is not None and end >= start:
        return (end - start).total_seconds()
    return None
