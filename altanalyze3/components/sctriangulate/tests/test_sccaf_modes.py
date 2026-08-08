"""Guard the two SCCAF modes: same numbers by default, reproducible in both.

The optimized mode is the default. It must not drift from the legacy mode
without someone noticing, and neither mode may return a different number on a
second call.

Run:
  cd /Users/saljh8/Documents/GitHub/altanalyze3
  python3.11 -m pytest altanalyze3/components/sctriangulate/tests/test_sccaf_modes.py -q -s
"""

import os
import sys

import numpy as np
import pytest
import scanpy as sc

HERE = os.path.dirname(os.path.abspath(__file__))
COMPONENT = os.path.dirname(HERE)
REPO = os.path.abspath(os.path.join(COMPONENT, '..', '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from altanalyze3.components.sctriangulate import ScTriangulate      # noqa: E402
from altanalyze3.components.sctriangulate import metrics as M       # noqa: E402

DEMO = os.path.join(COMPONENT, 'benchmarks', 'data', 'demo_pbmc3k.h5ad')
KEYS = ['sctri_rna_leiden_1', 'sctri_rna_leiden_2', 'sctri_rna_leiden_3']


@pytest.fixture(scope='module')
def dense_adata():
    if not os.path.exists(DEMO):
        pytest.skip('demo dataset not present at {}'.format(DEMO))
    adata = sc.read(DEMO)
    sctri = ScTriangulate(dir=os.path.join(COMPONENT, 'benchmarks', 'out_tests'),
                          adata=adata, query=KEYS, predict_doublet=False, verbose=0)
    sctri._to_dense()
    return sctri.adata


def _score(adata, key, mode):
    acc, conf = M.SCCAF_score(adata, key, 'human', 2, False, mode=mode)
    return {str(k): float(v) for k, v in acc.items()}, conf


def test_default_mode_is_optimized():
    assert M.SCCAF_MODE == 'optimized'
    assert M.SCCAF_PER_CLASS is False, \
        'the per-class split changes accuracies; it must stay opt-in'


@pytest.mark.parametrize('key', KEYS)
def test_optimized_matches_legacy_exactly(dense_adata, key):
    """With default settings the sparse path must reproduce the dense numbers."""
    adata_c = M.check_filter_single_cluster(dense_adata, key)
    legacy, conf_legacy = _score(adata_c, key, 'legacy')
    optimized, conf_opt = _score(adata_c, key, 'optimized')

    assert set(legacy) == set(optimized)
    deltas = {c: abs(legacy[c] - optimized[c]) for c in legacy}
    worst = max(deltas.values())
    print('\n{}: {} clusters, max |delta accuracy| = {:.3e}'.format(key, len(legacy), worst))
    assert worst == 0.0, 'largest difference {} at cluster {}'.format(
        worst, max(deltas, key=deltas.get))
    assert np.array_equal(conf_legacy.values, conf_opt.values)


@pytest.mark.parametrize('mode', ['legacy', 'optimized'])
@pytest.mark.parametrize('key', ['sctri_rna_leiden_2'])
def test_mode_is_reproducible(dense_adata, mode, key):
    """A second call on the same input must return the same numbers.

    Upstream failed this: it left random_state=None. An earlier version of the
    optimized mode also failed it, by fitting the per-class problems on threads
    that share liblinear's process-global C random stream.
    """
    adata_c = M.check_filter_single_cluster(dense_adata, key)
    first, _ = _score(adata_c, key, mode)
    second, _ = _score(adata_c, key, mode)
    assert first == second


def test_per_class_is_backend_independent(dense_adata):
    """Opt-in per-class mode must not depend on how many cores are used."""
    key = 'sctri_rna_leiden_1'
    adata_c = M.check_filter_single_cluster(dense_adata, key)
    saved = M.SCCAF_PER_CLASS
    saved_jobs = M.SCCAF_N_JOBS
    try:
        M.SCCAF_PER_CLASS = True
        M.SCCAF_N_JOBS = 1
        serial, _ = _score(adata_c, key, 'optimized')
        M.SCCAF_N_JOBS = 4
        parallel, _ = _score(adata_c, key, 'optimized')
    finally:
        M.SCCAF_PER_CLASS = saved
        M.SCCAF_N_JOBS = saved_jobs
    assert serial == parallel, 'per-class scores changed with the worker count'


def test_sparse_conversion_matches_dense_subset(dense_adata):
    """_csr_from_columns must hold exactly the values the dense subset holds."""
    key = 'sctri_rna_leiden_1'
    adata_c = M.check_filter_single_cluster(dense_adata, key)
    artifact = M._artifact_gene_set('human', 2)
    keep = np.asarray(~adata_c.var_names.isin(artifact))
    dense = M._as_dense(adata_c[:, keep].X)
    sparse = M._csr_from_columns(adata_c.X, keep)
    assert sparse.shape == dense.shape
    assert np.array_equal(sparse.toarray(), dense)
