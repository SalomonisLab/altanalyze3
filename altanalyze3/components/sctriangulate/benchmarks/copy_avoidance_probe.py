"""
Test three proposed memory changes against the code that is there now.

For each one: how many bytes does it avoid, how long does it take, and is the
result bit-identical to what the current code produces? A change that is not
bit-identical is not a free win, whatever it saves.

Writes /Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/
sctriangulate/benchmarks/copy_avoidance_probe.json
"""

import os
import json
import time
import tracemalloc

import numpy as np
from scipy.sparse import csr_matrix, random as sprandom

HERE = os.path.dirname(os.path.abspath(__file__))


def build(n_cells=100000, n_genes=20000, density=0.05, n_groups=120, seed=0):
    rng = np.random.default_rng(seed)
    X = sprandom(n_cells, n_genes, density=density, format='csr',
                 dtype=np.float32, random_state=seed)
    codes = rng.integers(0, n_groups, n_cells)
    return X, codes, n_groups


def peak(fn):
    tracemalloc.start()
    t0 = time.perf_counter()
    out = fn()
    secs = time.perf_counter() - t0
    _cur, pk = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return out, pk, secs


def main():
    X, codes, n_groups = build()
    n_cells, n_genes = X.shape
    print('matrix {} x {}, nnz {:,}, indices dtype {}'.format(
        n_cells, n_genes, X.nnz, X.indices.dtype))
    print('  data {:.0f} MB, indices {:.0f} MB, indptr {:.1f} MB'.format(
        X.data.nbytes / 1e6, X.indices.nbytes / 1e6, X.indptr.nbytes / 1e6))

    onehot_dense = np.zeros((n_cells, n_groups), dtype=np.float64)
    onehot_dense[np.arange(n_cells), codes] = 1.0
    onehot_sparse = csr_matrix(
        (np.ones(n_cells), (np.arange(n_cells), codes)), shape=(n_cells, n_groups))

    out = {'n_cells': n_cells, 'n_genes': n_genes, 'nnz': int(X.nnz),
           'indices_dtype': str(X.indices.dtype), 'n_groups': n_groups}

    # --- 1. squaring the matrix for the t-test: full copy vs shared indices ---
    def now_copy():
        SQ = X.copy()                       # what metrics.py:565 does today
        SQ.data = SQ.data * SQ.data
        return SQ

    def shared_indices():
        return csr_matrix((X.data * X.data, X.indices, X.indptr),
                          shape=X.shape, copy=False)

    a, pk_a, t_a = peak(now_copy)
    b, pk_b, t_b = peak(shared_indices)
    same = (np.array_equal(a.data, b.data) and np.array_equal(a.indices, b.indices)
            and np.array_equal(a.indptr, b.indptr))
    out['square'] = {'now_peak_MB': pk_a / 1e6, 'now_seconds': t_a,
                     'shared_peak_MB': pk_b / 1e6, 'shared_seconds': t_b,
                     'bit_identical': bool(same),
                     'shares_index_memory': bool(np.shares_memory(b.indices, X.indices))}
    print('\n1. square for the t-test')
    print('   full copy      peak {:7.1f} MB  {:.2f} s'.format(pk_a / 1e6, t_a))
    print('   shared indices peak {:7.1f} MB  {:.2f} s'.format(pk_b / 1e6, t_b))
    print('   bit identical: {}'.format(same))
    del a, b

    # --- 2. nonzero indicator: full copy vs shared indices ---
    def now_ind():
        ind = X.tocsr(copy=True)            # what metrics.py:201 does today
        ind.data = np.ones_like(ind.data, dtype=np.float32)
        return ind

    def shared_ind():
        return csr_matrix((np.ones(X.nnz, dtype=np.float32), X.indices, X.indptr),
                          shape=X.shape, copy=False)

    a, pk_a, t_a = peak(now_ind)
    b, pk_b, t_b = peak(shared_ind)
    same = (np.array_equal(a.data, b.data) and np.array_equal(a.indices, b.indices))
    out['indicator'] = {'now_peak_MB': pk_a / 1e6, 'now_seconds': t_a,
                        'shared_peak_MB': pk_b / 1e6, 'shared_seconds': t_b,
                        'bit_identical': bool(same)}
    print('\n2. nonzero indicator')
    print('   full copy      peak {:7.1f} MB  {:.2f} s'.format(pk_a / 1e6, t_a))
    print('   shared indices peak {:7.1f} MB  {:.2f} s'.format(pk_b / 1e6, t_b))
    print('   bit identical: {}'.format(same))
    del a, b

    # --- 3. one-hot: dense float64 vs sparse, for the group sums ---
    def dense_sums():
        return np.asarray(onehot_dense.T @ X)

    def sparse_sums():
        return np.asarray((onehot_sparse.T @ X).todense())

    a, pk_a, t_a = peak(dense_sums)
    b, pk_b, t_b = peak(sparse_sums)
    delta = float(np.abs(a - b).max())
    bit = bool(np.array_equal(a, b))
    out['onehot'] = {'dense_onehot_MB': onehot_dense.nbytes / 1e6,
                     'sparse_onehot_MB': (onehot_sparse.data.nbytes
                                          + onehot_sparse.indices.nbytes
                                          + onehot_sparse.indptr.nbytes) / 1e6,
                     'dense_matmul_peak_MB': pk_a / 1e6, 'dense_seconds': t_a,
                     'sparse_matmul_peak_MB': pk_b / 1e6, 'sparse_seconds': t_b,
                     'max_abs_diff': delta, 'bit_identical': bit}
    print('\n3. one-hot for the group sums')
    print('   dense  one-hot {:7.1f} MB, matmul peak {:7.1f} MB, {:.2f} s'.format(
        onehot_dense.nbytes / 1e6, pk_a / 1e6, t_a))
    print('   sparse one-hot {:7.1f} MB, matmul peak {:7.1f} MB, {:.2f} s'.format(
        out['onehot']['sparse_onehot_MB'], pk_b / 1e6, t_b))
    print('   max|diff| {:.3e}   bit identical: {}'.format(delta, bit))

    path = os.path.join(HERE, 'copy_avoidance_probe.json')
    json.dump(out, open(path, 'w'), indent=1)
    print('\nwrote {}'.format(path))


if __name__ == '__main__':
    main()
