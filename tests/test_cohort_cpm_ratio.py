"""Cohort pseudobulk CPM and isoform ratio (altanalyze3.components.long_read.isoform_matrix).

CPM   = counts summed over a cluster's cells in one sample / that pseudobulk's total over ALL features * 1e6.
ratio = isoform counts / the gene's counts in the same pseudobulk; 0 when the gene is expressed (gene CPM >=
        threshold) and the isoform has no reads; NaN when the gene is below the threshold.
Every expected value is recomputed by hand from the cell-level counts.
"""
import glob
import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from altanalyze3.components.long_read import comparisons as cmp
from altanalyze3.components.long_read import isoform_matrix as iso

FEATS = (["A:iso1", "A:iso2", "A:iso3", "B:iso1", "B:iso2", "C:iso1", "C:iso2", "D:iso1", "D:iso2", "HIGH:iso1"]
         + [f"F{i:04d}:iso1" for i in range(988)] + ["S:iso1", "S:iso2", "S:iso3", "S:iso4"])
FI = {f: i for i, f in enumerate(FEATS)}


def _sample(name, cells):
    M = np.zeros((len(cells), len(FEATS)), dtype=np.float32)
    for r, (_, d) in enumerate(cells):
        for f, c in d.items():
            M[r, FI[f]] = c
    obs = pd.DataFrame({"cluster": [c for c, _ in cells]}, index=[f"{name}_c{i}" for i in range(len(cells))])
    return ad.AnnData(X=sparse.csr_matrix(M), obs=obs, var=pd.DataFrame(index=FEATS))


S1 = _sample("S1", [
    ("X", {"HIGH:iso1": 3_000_000, "A:iso1": 400, "A:iso2": 100, "B:iso1": 1, "B:iso2": 1, "C:iso1": 30,
           "D:iso1": 4, "D:iso2": 1, "S:iso1": 20, "S:iso4": 20}),
    ("Y", {"A:iso1": 50, "A:iso3": 25, "B:iso1": 40, "B:iso2": 10, "C:iso1": 5, "S:iso3": 7})])
S2 = _sample("S2", [
    ("X", {"A:iso1": 80, "A:iso2": 20, "B:iso1": 6, "S:iso1": 3}),
    ("Y", {"A:iso2": 9, "B:iso2": 4, "C:iso1": 2, "C:iso2": 2, "S:iso4": 11})])
assert FI["S:iso1"] < 1000 <= FI["S:iso4"]          # gene S straddles the 1,000-feature chunk boundary


def _hand(col, f, threshold=1.0):
    a, cl = {"S1": S1, "S2": S2}[col.split(".")[1]], col.split(".")[0]
    v = np.asarray(a.X[(a.obs["cluster"] == cl).to_numpy()].sum(axis=0)).ravel().astype(float)
    total = v.sum()
    gene = f.split(":")[0]
    gsum = sum(v[FI[g]] for g in FEATS if g.split(":")[0] == gene)
    ratio = v[FI[f]] / gsum if gsum / total * 1e6 >= threshold else np.nan
    return v[FI[f]] / total * 1e6, ratio


def _combine(tmp_path, monkeypatch, dense=False, threshold=1):
    files = []
    for s, a in (("S1", S1), ("S2", S2)):
        p = tmp_path / f"{s}-isoform.h5ad"
        a.write_h5ad(p)
        files.append(str(p))
    monkeypatch.setenv("ALTANALYZE3_DENSE_COMBINE", "1" if dense else "0")
    monkeypatch.chdir(tmp_path)
    iso.concatenate_h5ad_and_compute_pseudobulks_optimized(files, collection_name="isoform", compute_cpm=True,
                                                           cpm_threshold=threshold, cluster_order=["X", "Y"],
                                                           min_group_size=1)
    cpm = ad.read_h5ad(tmp_path / "isoform_combined_pseudo_cluster_cpm.h5ad")
    rat = ad.read_h5ad(tmp_path / "isoform_combined_pseudo_cluster_ratio.h5ad")
    cols = {c.rsplit("-isoform", 1)[0]: i for i, c in enumerate(cpm.obs_names)}
    return cpm, rat, cols


CASES = [("X.S1", "A:iso1"), ("X.S1", "A:iso3"), ("Y.S1", "A:iso2"), ("X.S1", "B:iso1"), ("X.S1", "B:iso2"),
         ("Y.S1", "B:iso1"), ("X.S2", "C:iso1"), ("X.S1", "S:iso4"), ("Y.S2", "S:iso4"), ("X.S1", "D:iso2")]


@pytest.mark.parametrize("dense", [False, True])
def test_cpm_and_ratio_match_hand_recomputation(tmp_path, monkeypatch, dense):
    cpm, rat, cols = _combine(tmp_path, monkeypatch, dense=dense)
    R = iso.ratio_dense(rat)
    rs = np.asarray(cpm.X.sum(axis=1)).ravel()
    assert np.allclose(rs, 1e6, atol=1.0)
    for col, f in CASES:
        h_cpm, h_ratio = _hand(col, f)
        assert float(cpm.X[cols[col], FI[f]]) == pytest.approx(h_cpm, rel=1e-5, abs=1e-6)
        got = R[cols[col], FI[f]]
        assert (np.isnan(h_ratio) and np.isnan(got)) or got == pytest.approx(h_ratio, abs=1e-6), (col, f)
    assert not glob.glob(str(tmp_path / "*_cpm.txt*")) and not glob.glob(str(tmp_path / "*_ratio.txt*"))


def test_differential_loader_reads_ratio_zero_and_nan(tmp_path, monkeypatch):
    _combine(tmp_path, monkeypatch)
    df = cmp._load_isoform_filtered_matrix("ratio")
    col = [c for c in df.columns if c.startswith("X.S1")][0]
    assert df.at["A:iso3", col] == 0.0             # gene A expressed in X.S1, isoform without reads
    assert np.isnan(df.at["B:iso1", col])          # gene B CPM < 1 in X.S1
    assert df.at["A:iso1", col] == pytest.approx(400 / 500)
    cmp._load_isoform_filtered_matrix("cpm")       # cohort CPM is accepted


@pytest.mark.parametrize("dense", [False, True])
def test_cpm_threshold_reaches_both_paths(tmp_path, monkeypatch, dense):
    # gene D in X.S1: 5 reads of ~3.0e6 -> gene CPM ~1.66: defined at threshold 1, missing at threshold 2
    (tmp_path / "t1").mkdir()
    (tmp_path / "t2").mkdir()
    _, rat1, cols = _combine(tmp_path / "t1", monkeypatch, dense, 1)
    _, rat2, _ = _combine(tmp_path / "t2", monkeypatch, dense, 2)
    r, c = cols["X.S1"], FI["D:iso2"]
    assert iso.ratio_dense(rat1)[r, c] == pytest.approx(0.2)
    assert np.isnan(iso.ratio_dense(rat2)[r, c])


def test_gene_collection_writes_cpm_only(tmp_path, monkeypatch):
    genes = ["G1", "G2", "G3"]
    X = sparse.csr_matrix(np.array([[10, 0, 5], [0, 0, 0], [1, 1, 1]], dtype=float))
    monkeypatch.chdir(tmp_path)
    iso.write_cpm_and_ratio_h5ads("gene", X=X, features=genes, columns=["A.s1", "A.s2", "B.s1"],
                                  cluster_order=["A", "B"], min_group_size=1)
    assert os.path.exists("gene_combined_pseudo_cluster_cpm.h5ad")
    assert not glob.glob("gene_combined_pseudo_cluster_ratio*")
    cpm = ad.read_h5ad("gene_combined_pseudo_cluster_cpm.h5ad").X.toarray()
    assert cpm[0].tolist() == pytest.approx([1e6 * 10 / 15, 0, 1e6 * 5 / 15])


@pytest.mark.parametrize("bad", [np.nan, np.inf, -1.0])
def test_invalid_counts_are_rejected(bad):
    X = sparse.csr_matrix(np.array([[bad, 2.0]]))
    with pytest.raises(ValueError):
        iso.cpm_and_ratio_from_counts(X, ["G:i1", "G:i2"])


def test_per_sample_cpm_is_refused(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError):
        iso.pseudo_cluster_counts_optimized("S1-isoform", S1, compute_cpm=True)
    out = iso.pseudo_cluster_counts_optimized("S1-isoform", S1)
    counts = pd.read_csv(out, sep="\t", index_col=0)
    hand = S1.to_df().groupby(S1.obs["cluster"].to_numpy()).sum().T
    assert np.array_equal(counts.loc[FEATS, ["X", "Y"]].to_numpy(), hand.loc[FEATS, ["X", "Y"]].to_numpy())


def test_loader_refuses_stale_matrices(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    stale = ad.AnnData(X=sparse.csr_matrix(np.array([[2e6, 1e6], [5e5, 5e5]])),
                       obs=pd.DataFrame(index=["A.s1", "B.s1"]), var=pd.DataFrame(index=["G:i1", "G:i2"]))
    stale.write_h5ad("isoform_combined_pseudo_cluster_cpm-filtered.h5ad")
    stale.write_h5ad("isoform_combined_pseudo_cluster_ratio-filtered.h5ad")
    with pytest.raises(ValueError):
        cmp._load_isoform_filtered_matrix("cpm")
    with pytest.raises(ValueError):
        cmp._load_isoform_filtered_matrix("ratio")
