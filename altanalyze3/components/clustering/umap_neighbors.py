"""Bounded exact correlation search for opt-in landmark UMAP projection.

Centering and unit-normalizing search vectors retains every feature and makes
Pearson correlation a matrix multiplication. The fitted UMAP and its transform
graph/coordinate optimization remain unchanged. No full query-by-fit matrix is
materialized: each search tile has a fixed memory budget.
"""
from __future__ import annotations

import numpy as np


def _search_vectors(values):
    values = np.array(values, dtype=np.float32, order="C", copy=True)
    if values.ndim != 2 or values.shape[1] < 1 or not np.isfinite(values).all():
        raise ValueError("Correlation search requires finite rows with the complete feature panel.")
    # Match correlation's higher precision accumulation; stored search rows stay
    # float32. Constant rows follow UMAP/PyNNDescent's explicit distance rules.
    values -= values.mean(axis=1, keepdims=True, dtype=np.float64)
    norm = np.sqrt(np.einsum("ij,ij->i", values, values, dtype=np.float64))
    constant = norm == 0
    values /= np.where(constant, 1, norm)[:, None]
    return values, constant


class ExactCorrelationIndex:
    """UMAP query interface, in original fitted-row order with Pearson distances."""

    _angular_trees = False

    def __init__(self, fitted_rows, *, working_memory_bytes=64 * 1024**2):
        self._data, self._constant = _search_vectors(fitted_rows)
        if len(self._data) < 1 or working_memory_bytes < len(self._data) * 12:
            raise ValueError("Correlation search needs at least one fit row and one search tile.")
        # Float32 similarity plus int64 partition positions are the two large
        # temporaries. The query input and k-neighbor outputs are separate.
        self.block_rows = max(1, int(working_memory_bytes) // (len(self._data) * 12))

    def query(self, query_data, k=10, epsilon=0.12):
        # Respect small deployment CPU allocations instead of NumPy's host-wide
        # BLAS default. The context restores the previous pool setting.
        from numba import get_num_threads
        from threadpoolctl import threadpool_limits
        with threadpool_limits(limits=min(4, get_num_threads()), user_api="blas"):
            return self._query(query_data, k)

    def _query(self, query_data, k):
        if not 1 <= k <= len(self._data):
            raise ValueError("The neighbor count must be within the fitted cell roster.")
        data, constant = _search_vectors(query_data)
        if data.shape[1] != self._data.shape[1]:
            raise ValueError("Query and fit rows must contain the same complete feature panel.")
        indices = np.empty((len(data), k), dtype=np.int32)
        distances = np.empty((len(data), k), dtype=np.float32)
        for start in range(0, len(data), self.block_rows):
            stop = min(start + self.block_rows, len(data))
            similarity = data[start:stop] @ self._data.T
            # UMAP defines two constant rows as correlation distance zero, and
            # a constant/nonconstant pair as distance one.
            if constant[start:stop].any() and self._constant.any():
                similarity[np.ix_(constant[start:stop], self._constant)] = 1
            positions = np.argpartition(similarity, -k, axis=1)[:, -k:].copy()
            d = np.clip(1 - np.take_along_axis(similarity, positions, axis=1), 0, 2)
            order = np.lexsort((positions, d), axis=1)
            indices[start:stop] = np.take_along_axis(positions, order, axis=1)
            distances[start:stop] = np.take_along_axis(d, order, axis=1)
            del similarity, positions, d, order
        return indices, distances
