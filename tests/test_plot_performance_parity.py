"""Faster serving keeps the original statistics, selections and point ordering."""
from importlib import import_module
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from altanalyze3.components.cellHarmony.webapp.job_bundle import _StoreMatrix
from test_plot_transport import unpack

W = import_module("altanalyze3.components.cellHarmony.webapp.app")


@pytest.mark.parametrize("permuted", [False, True])
@pytest.mark.parametrize("dtype", [np.float16, np.float32])
def test_fused_bundle_reductions_are_bitwise_identical(permuted, dtype):
    rng = np.random.default_rng(19)
    n, genes = 90, 24
    matrix = rng.normal(size=(n, genes)).astype(np.float32)
    matrix[rng.random(matrix.shape) < .6] = 0
    order = rng.permutation(n) if permuted else np.arange(n)
    csc = sparse.csc_matrix(matrix[order])
    csc.data = csc.data.astype(dtype)
    owner = SimpleNamespace(n_obs=n, n_vars=genes, _indptr=csc.indptr,
                            _indices=csc.indices, _data=csc.data,
                            _to_h5ad=order if permuted else None, _sparse=True)
    x = _StoreMatrix(owner)
    codes = rng.integers(0, 5, n)
    # Include negative weights and cells outside the selected groups.
    selected = np.flatnonzero(codes < 3)
    indicator = sparse.csr_matrix((rng.normal(size=len(selected)), (codes[selected], selected)), shape=(3, n))
    expected_group = indicator @ x
    expected_total = x.sum(axis=0)
    actual_group, actual_total = x.group_sums_and_total(indicator)
    assert expected_group.tobytes() == actual_group.tobytes()
    assert expected_total.ravel().tobytes() == actual_total.tobytes()


def cache_for(x):
    n = x.shape[0]
    obs = pd.DataFrame({"state": pd.Categorical(np.resize(["B", "A"], n), categories=["B", "A"]),
                        "Library": np.resize(["s1", "s2", "s1"], n)}, index=[f"c{i}α" for i in range(n)])
    a = ad.AnnData(x, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(x.shape[1])]))
    return dict(adata=a, cluster_key="state", populations=obs.state.astype(str).to_numpy(),
                sample_field="Library", obs_names=a.obs_names.to_numpy(), var_names=a.var_names.to_numpy(),
                umap_x=np.arange(n, dtype=float), umap_y=-np.arange(n, dtype=float),
                sample_labels=obs.Library.to_numpy(), obs_filter_values={"Library": obs.Library.to_numpy()})


def test_marker_choices_cached_by_grouping_and_returned_as_copies(monkeypatch):
    cache = cache_for(sparse.csr_matrix(np.random.default_rng(3).random((80, 15)).astype(np.float32)))
    original = W._compute_default_marker_genes
    called = []

    def compute(*args):
        called.append(args[1:])
        return original(*args)

    monkeypatch.setattr(W, "_compute_default_marker_genes", compute)
    expected = original(cache)
    chosen = W._default_marker_genes(cache)
    assert chosen == expected
    chosen.clear()
    assert W._default_marker_genes(cache) == expected
    assert len(called) == 1
    assert W._default_marker_genes(cache, "Library") == original(cache, "Library")
    assert len(called) == 2


@pytest.mark.parametrize("limit", [0, 5, 10, 20, 50])
def test_sparse_batch_combplot_matches_dense_individual_columns(limit):
    rng = np.random.default_rng(7)
    x = rng.normal(size=(145, 130)).astype(np.float32)
    x[rng.random(x.shape) < .4] = 0
    wanted = [f"g{i}" for i in range(129, -1, -1)] + ["absent", "g10"]
    options = dict(group_by="Library", subset_by="state", subset_values=["A"], cells_per_sample=limit)
    assert W._gene_cell_values(cache_for(sparse.csr_matrix(x)), wanted, **options) == W._gene_cell_values(cache_for(x), wanted, **options)


@pytest.mark.parametrize("filters", [None, [("Library", ["s1"])], [("Library", ["absent"])]])
def test_direct_compact_points_match_legacy_with_ties_and_missing_coordinates(monkeypatch, filters):
    cache = cache_for(np.array([[0], [2], [2], [-0.], [np.nan], [1.]], dtype=np.float32))
    cache["umap_x"][1] = np.nan
    cache["umap_y"][2] = np.inf
    monkeypatch.setattr(W, "_get_expression_cache", lambda *a, **k: cache)
    monkeypatch.setattr(W, "_load_reference_adata", lambda *a: None)
    for color_by in ("", "Library"):
        options = dict(display_filters=filters, color_by=color_by)
        legacy = W._build_umap_payload(None, {}, **options)
        assert unpack(W._build_umap_payload(None, {}, compact=True, **options)) == legacy
    for view in ("all", "umap", "violin"):
        legacy = W._build_expression_payload(None, {}, "g0", display_filters=filters, view=view)
        assert unpack(W._build_expression_payload(None, {}, "g0", display_filters=filters, view=view, compact=True)) == legacy
