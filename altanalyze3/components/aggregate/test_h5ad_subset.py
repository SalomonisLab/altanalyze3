"""Prove h5ad_subset copies the selected rows byte-for-byte, for CSR X and CSR layers."""

import os
import tempfile

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")
sparse = pytest.importorskip("scipy.sparse")

from altanalyze3.components.aggregate.h5ad_subset import subset_h5ad


def _toy(n_obs=200, n_vars=40, seed=0):
    import pandas as pd
    rng = np.random.default_rng(seed)
    dense = rng.poisson(0.7, size=(n_obs, n_vars)).astype(np.float64)
    X = sparse.csr_matrix(dense)
    obs = pd.DataFrame(
        {"Library": np.where(np.arange(n_obs) % 4 == 0, "Thymus", "Marrow")},
        index=[f"BC{i:04d}-1" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"Gene{j}" for j in range(n_vars)])
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    adata.layers["counts"] = sparse.csr_matrix(dense * 3)
    adata.obsm["X_umap"] = rng.normal(size=(n_obs, 2))
    adata.uns["lineage_order"] = ["a", "b"]
    return adata


def test_subset_matches_the_source_rows_exactly():
    with tempfile.TemporaryDirectory() as tmp:
        src = _toy()
        src_path = os.path.join(tmp, "src.h5ad")
        out_path = os.path.join(tmp, "out.h5ad")
        src.write(src_path)

        kept, total = subset_h5ad(src_path, out_path, "Library", ["Thymus"], log=lambda *_: None)
        assert (kept, total) == (50, 200)

        out = anndata.read_h5ad(out_path)
        wanted = src[src.obs["Library"] == "Thymus"]
        assert list(out.obs_names) == list(wanted.obs_names)
        assert list(out.var_names) == list(src.var_names)
        assert np.array_equal(out.X.toarray(), wanted.X.toarray())
        assert np.array_equal(out.layers["counts"].toarray(), wanted.layers["counts"].toarray())
        assert np.allclose(out.obsm["X_umap"], wanted.obsm["X_umap"])
        assert list(out.uns["lineage_order"]) == ["a", "b"]
        assert out.X.nnz == wanted.X.nnz


def test_subset_rejects_a_value_the_column_does_not_hold():
    with tempfile.TemporaryDirectory() as tmp:
        src_path = os.path.join(tmp, "src.h5ad")
        _toy().write(src_path)
        with pytest.raises(ValueError, match="holds no value"):
            subset_h5ad(src_path, os.path.join(tmp, "o.h5ad"), "Library", ["Spleen"],
                        log=lambda *_: None)


def test_subset_rejects_an_absent_column():
    with tempfile.TemporaryDirectory() as tmp:
        src_path = os.path.join(tmp, "src.h5ad")
        _toy().write(src_path)
        with pytest.raises(KeyError):
            subset_h5ad(src_path, os.path.join(tmp, "o.h5ad"), "Tissue", ["Thymus"],
                        log=lambda *_: None)
