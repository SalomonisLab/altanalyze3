"""Shared byte-bounded LRU retention for web results; eviction never closes live views.

The budget bounds retained Python/NumPy buffers, not request-local allocations or
reclaimable mmap pages. Callers may continue using an evicted value safely.
"""
from collections import OrderedDict
import sys
import threading
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp


def retained_bytes(value, seen=None):
    seen = set() if seen is None else seen
    if id(value) in seen:
        return 0
    seen.add(id(value))
    if isinstance(value, np.ndarray):
        base = value
        while isinstance(base, np.ndarray):
            if isinstance(base, np.memmap):
                return 0  # file-backed pages are reclaimable, not owned anonymous buffers
            owner = base
            base = base.base
        if owner is not value:
            if id(owner) in seen:
                return 0
            seen.add(id(owner))
        size = owner.nbytes
        if owner.dtype.hasobject and owner.size:
            sample = owner.flat[:min(owner.size, 32)]
            size += int(sum(sys.getsizeof(v) for v in sample) * owner.size / len(sample))
        return size
    if sp.issparse(value):
        return sum(retained_bytes(getattr(value, k), seen) for k in ('data', 'indices', 'indptr')
                   if hasattr(value, k))
    if isinstance(value, pd.DataFrame):
        return int(value.memory_usage(index=True, deep=True).sum())
    if isinstance(value, (pd.Series, pd.Index)):
        return int(value.memory_usage(deep=True))
    if isinstance(value, dict):
        return sys.getsizeof(value) + sum(retained_bytes(k, seen) + retained_bytes(v, seen)
                                         for k, v in value.items())
    if isinstance(value, (list, tuple, set)):
        return sys.getsizeof(value) + sum(retained_bytes(v, seen) for v in value)
    # AnnData and bundle stand-ins: account for layers/raw/obsm as well as X.
    if hasattr(value, '__dict__') and type(value).__module__.startswith(('anndata.', 'altanalyze3.')):
        return sys.getsizeof(value) + retained_bytes(vars(value), seen)
    return sys.getsizeof(value)


class CacheBudget:
    def __init__(self, max_bytes=2 * 1024**3, max_entries=64, ttl_seconds=600, clock=time.monotonic):
        self.max_bytes = max(0, int(max_bytes))
        self.max_entries = max(0, int(max_entries))
        self.ttl_seconds = max(0, float(ttl_seconds))
        self.clock = clock
        self.lock = threading.RLock()
        self.records = OrderedDict()
        self.bytes = 0
        self.evictions = 0

    def _drop(self, token):
        cache, key, size, _ = self.records.pop(token)
        dict.pop(cache, key, None)
        self.bytes -= size
        self.evictions += 1

    def _expire(self):
        now = self.clock()
        for token, (_, _, _, touched) in list(self.records.items()):
            if now - touched >= self.ttl_seconds:
                self._drop(token)

    def make_room(self, incoming_bytes=0):
        """Release idle retained entries *before* a new whole-matrix read."""
        with self.lock:
            self._expire()
            target = max(0, self.max_bytes - max(0, incoming_bytes))
            while self.records and self.bytes > target:
                self._drop(next(iter(self.records)))

    def snapshot(self):
        with self.lock:
            self._expire()
            return dict(bytes=self.bytes, max_bytes=self.max_bytes,
                        entries=len(self.records), evictions=self.evictions)


class BoundedCache(dict):
    """Dict-compatible LRU sharing one budget with other serving caches.

    Mutated entries must be assigned again to refresh their byte accounting.
    Oversized entries are returned by their builder but are never retained here.
    """
    def __init__(self, budget=None):
        super().__init__()
        self.budget = budget if budget is not None else CacheBudget()

    def __setitem__(self, key, value):
        size = retained_bytes(value)
        budget = self.budget
        with budget.lock:
            self.pop(key, None)
            budget._expire()
            if size > budget.max_bytes or not budget.max_entries or not budget.ttl_seconds:
                return
            budget.make_room(size)
            while len(budget.records) >= budget.max_entries:
                budget._drop(next(iter(budget.records)))
            dict.__setitem__(self, key, value)
            budget.records[(id(self), key)] = (self, key, size, budget.clock())
            budget.bytes += size

    def __getitem__(self, key):
        with self.budget.lock:
            self.budget._expire()
            value = dict.__getitem__(self, key)
            token = (id(self), key)
            if token in self.budget.records:
                cache, k, size, _ = self.budget.records[token]
                self.budget.records[token] = (cache, k, size, self.budget.clock())
                self.budget.records.move_to_end(token)
            return value

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        with self.budget.lock:
            self.budget._expire()
            return dict.__contains__(self, key)

    def __iter__(self):
        with self.budget.lock:
            self.budget._expire()
            return iter(tuple(dict.keys(self)))

    def pop(self, key, *default):
        with self.budget.lock:
            record = self.budget.records.pop((id(self), key), None)
            if record is not None:
                self.budget.bytes -= record[2]
            return dict.pop(self, key, *default)

    def __delitem__(self, key):
        self.pop(key)

    def clear(self):
        with self.budget.lock:
            for key in list(self):
                self.pop(key)

    def update(self, other=(), **kwargs):
        for key, value in dict(other, **kwargs).items():
            self[key] = value

    def __ior__(self, other):
        self.update(other)
        return self

    def setdefault(self, key, default=None):
        with self.budget.lock:
            self.budget._expire()
            if key in self:
                return self[key]
            self[key] = default
            return default
