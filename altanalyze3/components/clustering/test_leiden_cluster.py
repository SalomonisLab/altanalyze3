"""Guard the leiden_cluster CLI additions: --resolution, --normalization, --skip-hvg, --n-pcs.

The defaults must reproduce the behaviour the module had before these flags existed:
resolution 0.5, cp10k-log1p normalization and dispersion-binned highly-variable selection.
"""

import os
import sys
import tempfile

import numpy as np
import pytest

sys.modules.setdefault("tensorflow", None)  # see leiden_cluster._ensure_umap_importable

anndata = pytest.importorskip("anndata")
sparse = pytest.importorskip("scipy.sparse")
pd = pytest.importorskip("pandas")

from altanalyze3.components.clustering import leiden_cluster as LC


def _toy(path, n_cells=600, n_groups=5, n_vars=800, seed=1):
    rng = np.random.default_rng(seed)
    groups = np.repeat(np.arange(n_groups), n_cells // n_groups)
    block = n_vars // n_groups
    base = np.full((len(groups), n_vars), 0.15)
    for g in range(n_groups):
        base[groups == g, g * block:(g + 1) * block] = 6.0
    X = sparse.csr_matrix(rng.poisson(base).astype(np.float32))
    obs = pd.DataFrame({"Library": ["L1"] * len(groups)},
                       index=[f"BC{i:04d}-1.L1" for i in range(len(groups))])
    var = pd.DataFrame(index=[f"Gene{j}" for j in range(n_vars)])
    anndata.AnnData(X=X, obs=obs, var=var).write(path)
    return len(groups)


def _run(tmp, tag, **kwargs):
    cwd = os.getcwd()
    outdir = os.path.join(tmp, tag)
    os.makedirs(outdir, exist_ok=True)
    try:
        LC.combine_and_cluster(
            h5_files=[], h5ad_file=os.path.join(tmp, "src.h5ad"), output_dir=outdir,
            min_genes=0, min_cells=0, min_counts=0, mit_percent=100,
            generate_umap=False, **kwargs)
    finally:
        os.chdir(cwd)
    return pd.read_csv(os.path.join(outdir, "unsupervised_leiden_clusters.tsv"), sep="\t")


def test_higher_resolution_gives_at_least_as_many_clusters():
    with tempfile.TemporaryDirectory() as tmp:
        _toy(os.path.join(tmp, "src.h5ad"))
        low = _run(tmp, "r05", resolution=0.5)
        high = _run(tmp, "r20", resolution=2.0)
        assert high["Leiden"].nunique() >= low["Leiden"].nunique()
        assert len(low) == len(high) == 600


def test_resolution_default_is_the_pre_flag_value():
    import inspect
    sig = inspect.signature(LC.combine_and_cluster)
    assert sig.parameters["resolution"].default == 0.5
    assert sig.parameters["normalization"].default == "cp10k-log1p"
    assert sig.parameters["skip_hvg"].default is False
    assert sig.parameters["n_pcs"].default == 50


def test_skip_hvg_sends_every_feature_to_pca():
    with tempfile.TemporaryDirectory() as tmp:
        _toy(os.path.join(tmp, "src.h5ad"), n_vars=120)
        out = _run(tmp, "skip", resolution=2.0, normalization="log1p", skip_hvg=True, n_pcs=30)
        assert len(out) == 600
        assert out["Leiden"].nunique() >= 2


def test_normalization_none_leaves_the_matrix_alone():
    adata = anndata.AnnData(
        X=sparse.csr_matrix(np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)))
    LC.normalize_adata(adata, mode="none")
    assert np.array_equal(adata.X.toarray(), np.array([[1.0, 2.0], [3.0, 4.0]]))


def test_normalization_log1p_skips_the_depth_step():
    values = np.array([[1.0, 3.0], [10.0, 30.0]], dtype=np.float32)
    adata = anndata.AnnData(X=sparse.csr_matrix(values))
    LC.normalize_adata(adata, mode="log1p")
    assert np.allclose(adata.X.toarray(), np.log1p(values), atol=1e-6)


def test_unknown_normalization_raises():
    with pytest.raises(ValueError, match="Unknown normalization"):
        LC.normalize_adata(anndata.AnnData(X=np.ones((2, 2), dtype=np.float32)), mode="quantile")


def test_non_positive_resolution_raises():
    with tempfile.TemporaryDirectory() as tmp:
        _toy(os.path.join(tmp, "src.h5ad"))
        with pytest.raises(ValueError, match="must be positive"):
            _run(tmp, "bad", resolution=0.0)
