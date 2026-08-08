"""
check_filter_single_cluster now returns the caller's own object when no cluster
holds one cell, instead of a read-only-ish view. That removes the accidental
write-protection the view gave, so this file is the replacement guarantee: run a
full each_key_run and fail if any byte of X, obs or var moved.

Static reading is not enough. A write can reach the object through a local
alias, through a nested call, or through in-place numpy arithmetic that no
grep finds. These tests checksum the real buffers.
"""

import os
import sys
import hashlib

import numpy as np
import pandas as pd
import pytest
import anndata as ad
from scipy.sparse import csr_matrix, issparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))))))

from altanalyze3.components.sctriangulate import ScTriangulate            # noqa: E402
from altanalyze3.components.sctriangulate import metrics as M             # noqa: E402
from altanalyze3.components.sctriangulate.main_class import each_key_run  # noqa: E402

# the default tfidf5 add-metric that ScTriangulate installs at __init__
ADDED_KWARGS = [{'species': 'human', 'criterion': 2, 'layer': None}]


def digest(adata):
    """Byte fingerprint of everything each_key_run could plausibly touch."""
    h = hashlib.sha256()
    X = adata.X
    if issparse(X):
        for part in (X.data, X.indices, X.indptr):
            h.update(np.ascontiguousarray(part).tobytes())
    else:
        h.update(np.ascontiguousarray(X).tobytes())
    h.update(pd.util.hash_pandas_object(adata.obs, index=True).values.tobytes())
    h.update(pd.util.hash_pandas_object(adata.var, index=True).values.tobytes())
    for name in sorted(adata.layers.keys()):
        L = adata.layers[name]
        h.update(name.encode())
        h.update(np.ascontiguousarray(L.data if issparse(L) else L).tobytes())
    return h.hexdigest()


def build(tmpdir, n_cells=600, n_genes=600, n_clusters=6, singletons=0, seed=0):
    """
    Synthetic object with real cluster structure. Structure matters here:
    reassign_score pools the top marker genes and asks PCA for 30 components, so
    an unstructured matrix yields too few markers and the call raises before the
    test can check anything.
    """
    rng = np.random.default_rng(seed)
    labels = np.array(['c{}'.format(i) for i in rng.integers(0, n_clusters, n_cells)],
                      dtype=object)
    other = np.array(['d{}'.format(i) for i in rng.integers(0, 9, n_cells)], dtype=object)

    dense = ((rng.random((n_cells, n_genes)) < 0.3)
             * rng.random((n_cells, n_genes)) * 2.0).astype(np.float32)
    block = n_genes // n_clusters
    for c in range(n_clusters):                 # each cluster elevates its own block
        rows = np.flatnonzero(labels == 'c{}'.format(c))
        cols = slice(c * block, (c + 1) * block)
        dense[np.ix_(rows, np.arange(n_genes)[cols])] += rng.random(
            (len(rows), block)).astype(np.float32) * 8.0
    X = csr_matrix(dense)

    for i in range(singletons):
        labels[i] = 'solo{}'.format(i)          # clusters of exactly one cell
    genes = ['GENE{}'.format(i) for i in range(n_genes)]
    adata = ad.AnnData(X=X,
                       obs=pd.DataFrame({'anno1': labels, 'anno2': other},
                                        index=['cell{}'.format(i) for i in range(n_cells)]),
                       var=pd.DataFrame(index=genes))
    M.RUN_ENRICHMENT = False
    s = ScTriangulate(dir=str(tmpdir), adata=adata, query=['anno1', 'anno2'],
                      predict_doublet=False, verbose=0)
    s.run_enrichment = False
    return s


@pytest.mark.parametrize('singletons', [0, 3])
def test_each_key_run_does_not_write_to_adata(tmp_path, singletons):
    s = build(tmp_path, singletons=singletons)
    s._install_sparse_sidecar()
    before = digest(s.adata)
    each_key_run(s, 'anno1', scale_sccaf=False, layer=None,
                 added_metrics_kwargs=ADDED_KWARGS)
    assert digest(s.adata) == before, \
        'each_key_run modified adata (singletons={})'.format(singletons)


def test_filter_returns_the_object_itself_when_nothing_is_dropped(tmp_path):
    s = build(tmp_path, singletons=0)
    out = M.check_filter_single_cluster(s.adata, 'anno1')
    assert out is s.adata
    assert not out.is_view


def test_filter_still_drops_singletons(tmp_path):
    s = build(tmp_path, singletons=3)
    out = M.check_filter_single_cluster(s.adata, 'anno1')
    assert out is not s.adata
    assert out.is_view
    assert out.shape[0] == s.adata.shape[0] - 3
    assert not out.obs['anno1'].astype(str).str.startswith('solo').any()
    # the original is untouched
    assert s.adata.shape[0] == 600


def test_shared_index_csr_is_byte_equal_to_the_copy_it_replaces():
    rng = np.random.default_rng(0)
    src = csr_matrix((rng.random((300, 200)) < 0.2).astype(np.float32)
                     * rng.random((300, 200)).astype(np.float32))

    old = src.copy()
    old.data = old.data * old.data
    new = M._shared_index_csr(src, src.data * src.data)

    assert np.array_equal(old.data, new.data)
    assert np.array_equal(old.indices, new.indices)
    assert np.array_equal(old.indptr, new.indptr)
    assert old.shape == new.shape
    assert np.shares_memory(new.indices, src.indices)
    # src itself is untouched by the sharing
    assert not np.shares_memory(new.data, src.data)


def test_shared_index_csr_copies_when_indices_are_unsorted():
    """Sharing unsorted indices would let scipy sort them in place and permute
    them out of step with src.data."""
    data = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32)
    indices = np.array([2, 0, 3, 1], dtype=np.int32)       # deliberately unsorted
    indptr = np.array([0, 2, 4], dtype=np.int32)
    src = csr_matrix((data, indices, indptr), shape=(2, 4))
    assert not src.has_sorted_indices

    out = M._shared_index_csr(src, src.data * src.data)
    assert not np.shares_memory(out.indices, src.indices)
    assert np.array_equal(out.indices, indices)
    assert np.array_equal(out.toarray(), src.toarray() ** 2)
