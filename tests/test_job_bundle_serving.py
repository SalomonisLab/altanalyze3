"""A bundled upload job is served exactly as its h5ad files would serve it.

LARGE_DATASET_DESIGN.md steps A and B. flask/pipeline._build_job_bundle writes
outputs/bundle/ through precompute.py; webapp/job_bundle.view stands in for ad.read_h5ad on
the RNA and per-cell modality h5ads. These tests build a small job, bundle it through the
real pipeline step, and compare the stand-in with ad.read_h5ad bit for bit on everything
app.py reads: obs, var, uns, obsm (first two columns), feature slices, X.sum(axis=0), a
group indicator @ X, and dense row selection. A modality whose cells come in a different
order than the RNA h5ad is included on purpose.
"""
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from altanalyze3.components.cellHarmony.flask import pipeline as P
from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
from altanalyze3.components.cellHarmony.webapp import job_bundle as JB

N_CELLS, N_GENES, N_FEAT = 300, 80, 12


def _job(tmp_path):
    rng = np.random.default_rng(7)
    store = JobStore(tmp_path / "jobs")
    meta = store.create_job("human", "ref", None, files=[])
    job_id = meta["job_id"]
    out = store.outputs_dir(job_id)
    states = np.array([f"S{i % 5}" for i in range(N_CELLS)])
    barcodes = [f"cell{i:04d}" for i in range(N_CELLS)]
    obs = pd.DataFrame({"state": pd.Categorical(states, categories=["S3", "S0", "S1", "S2", "S4"]),
                        "Library": pd.Categorical(rng.choice(["A", "B", "C"], N_CELLS)),
                        "pct_counts_mt": rng.random(N_CELLS).astype(np.float32)},
                       index=pd.Index(barcodes))
    x = sp.random(N_CELLS, N_GENES, density=0.2, format="csr", random_state=rng, dtype=np.float32)
    x.data = (x.data * 5).astype(np.float32)
    rna = ad.AnnData(X=x, obs=obs, var=pd.DataFrame(index=[f"G{j}" for j in range(N_GENES)]))
    rna.obsm["X_umap"] = rng.normal(size=(N_CELLS, 2)).astype(np.float32)
    rna.obsm["X_adt"] = rng.random((N_CELLS, 7)).astype(np.float32)
    rna.uns["lineage_order"] = np.array(["S3", "S0", "S1", "S2", "S4"])
    rna_path = out / "combined_with_umap_and_markers.h5ad"
    rna.write_h5ad(rna_path)
    order = rng.permutation(N_CELLS)                    # a modality in another cell order
    dense = rng.normal(size=(N_CELLS, N_FEAT)).astype(np.float32)
    dense[dense < -1.0] = 0.0
    mod = ad.AnnData(X=dense[order], obs=obs.iloc[order].copy(),
                     var=pd.DataFrame(index=[f"TF{j}|T{j}" for j in range(N_FEAT)]))
    mod_path = out / "combined_with_umap_and_markers_grn_edges.h5ad"
    mod.write_h5ad(mod_path)
    artifacts = {"rna": {"h5ad": str(rna_path)}, "grn": {"h5ad": str(mod_path)}}
    payload = {"available": [{"id": "rna", "label": "RNA"},
                             {"id": "grn", "label": "GRN (edges)", "feature_label": "edge"}]}
    return store, job_id, rna_path, mod_path, artifacts, payload


def _flat(x):
    return np.asarray(x.todense() if sp.issparse(x) else x).ravel()


def test_bundle_skipped_below_threshold(tmp_path, monkeypatch):
    store, job_id, rna_path, _m, artifacts, payload = _job(tmp_path)
    monkeypatch.setenv("CELLHARMONY_BUNDLE_MIN_CELLS", str(N_CELLS + 1))
    rec = P._build_job_bundle(store, job_id, rna_path, "state", artifacts, payload)
    assert rec["status"] == "skipped" and "below" in rec["reason"]
    assert not (store.outputs_dir(job_id) / "bundle").exists()
    monkeypatch.setenv("CELLHARMONY_BUNDLE_MIN_CELLS", "many")
    rec = P._build_job_bundle(store, job_id, rna_path, "state", artifacts, payload)
    assert rec["status"] == "skipped" and "not an integer" in rec["reason"]


def test_bundle_view_matches_read_h5ad(tmp_path, monkeypatch):
    store, job_id, rna_path, mod_path, artifacts, payload = _job(tmp_path)
    monkeypatch.setenv("CELLHARMONY_BUNDLE_MIN_CELLS", "0")
    rec = P._build_job_bundle(store, job_id, rna_path, "state", artifacts, payload)
    assert rec["status"] == "completed", rec
    assert set(rec["sources"]) == {"rna", "grn"}
    meta = {"bundle": rec}
    for path, sparse in ((rna_path, True), (mod_path, False)):
        ref = ad.read_h5ad(path)
        got = JB.view(meta, path)
        assert got is not None, JB.refusals()
        assert ref.obs.equals(got.obs) and list(ref.obs.dtypes) == list(got.obs.dtypes)
        assert list(ref.var_names) == list(got.var_names)
        assert list(ref.obsm.keys()) == list(got.obsm.keys())
        for key in ref.obsm.keys():
            assert np.asarray(ref.obsm[key])[:, :2].tobytes() == np.asarray(got.obsm[key]).tobytes()
        mask = (ref.obs["state"] == "S2").to_numpy()
        for g in ref.var_names:
            for rows in (slice(None), mask):
                a, b = ref[rows, g].X, got[rows, g].X
                assert sp.issparse(a) == sp.issparse(b) == sparse
                assert _flat(a).tobytes() == _flat(b).tobytes(), (path.name, g)
        a = np.asarray(ref.X.sum(axis=0)).ravel()
        b = np.asarray(got.X.sum(axis=0)).ravel()
        assert a.dtype == b.dtype and a.tobytes() == b.tobytes()
        codes = pd.Categorical(ref.obs["state"]).codes.astype(np.int64)
        ind = sp.csr_matrix((np.ones(N_CELLS), (codes, np.arange(N_CELLS))), shape=(5, N_CELLS))
        a = ind @ ref.X
        a = np.asarray(a.todense() if sp.issparse(a) else a, dtype=np.float64)
        assert a.tobytes() == np.asarray(ind @ got.X, dtype=np.float64).tobytes()
        # X[rows] then the reductions its callers apply (grn_data, integration_data, GRN net)
        rows = np.flatnonzero(mask)
        cols = [3, 0, 7, 5]
        for sel in (mask, rows):
            a, b = ref.X[sel], got.X[sel]
            assert sp.issparse(a) == sp.issparse(b) == sparse
            assert _flat(a).tobytes() == _flat(b).tobytes()
            assert np.asarray(a.mean(axis=0)).tobytes() == np.asarray(b.mean(axis=0)).tobytes()
            assert (np.asarray(a[:, cols].mean(axis=0)).tobytes()
                    == np.asarray(b[:, cols].mean(axis=0)).tobytes())
        # X[:, j] and X[:, cols] (dot and comb plots)
        for j in (0, 4, N_FEAT - 1):
            a, b = ref.X[:, j], got.X[:, j]
            assert sp.issparse(a) == sp.issparse(b) and np.shape(a) == np.shape(b)
            assert _flat(a).tobytes() == _flat(b).tobytes()
        a, b = ref.X[:, cols], got.X[:, cols]
        assert np.shape(a) == np.shape(b) and _flat(a).tobytes() == _flat(b).tobytes()
        # adata[:, [cols]].X then rows, then a group product (cross_modal._group_means)
        va, vb = ref[:, cols].X, got[:, cols].X
        assert sp.issparse(va) == sp.issparse(vb) and np.shape(va) == np.shape(vb)
        grouper = ind[:, rows]
        ma = grouper @ va[rows].astype(np.float64)
        mb = grouper @ vb[rows].astype(np.float64)
        assert _flat(ma).tobytes() == _flat(mb).tobytes()
        with pytest.raises(NotImplementedError):
            got.X[:]


def test_misaligned_bundle_is_refused(tmp_path, monkeypatch):
    store, job_id, rna_path, mod_path, artifacts, payload = _job(tmp_path)
    monkeypatch.setenv("CELLHARMONY_BUNDLE_MIN_CELLS", "0")
    rec = P._build_job_bundle(store, job_id, rna_path, "state", artifacts, payload)
    # point the RNA source at the modality file: its features are not the RNA store's
    bad = dict(rec, sources={"rna": str(mod_path)})
    assert JB.view({"bundle": bad}, mod_path) is None
    assert any("features differ" in why for why in JB.refusals().values())


def test_float64_h5ad_is_refused(tmp_path, monkeypatch):
    store, job_id, rna_path, mod_path, artifacts, payload = _job(tmp_path)
    monkeypatch.setenv("CELLHARMONY_BUNDLE_MIN_CELLS", "0")
    rec = P._build_job_bundle(store, job_id, rna_path, "state", artifacts, payload)
    ref = ad.read_h5ad(rna_path)
    ref.X = ref.X.astype(np.float64)
    ref.write_h5ad(rna_path)                          # same cells and genes, now float64
    assert JB.view({"bundle": rec}, rna_path) is None
    assert any("float64" in why for why in JB.refusals().values())


def test_unbundled_job_reads_h5ad(tmp_path):
    store, job_id, rna_path, *_ = _job(tmp_path)
    assert JB.view(store.get_job(job_id), rna_path) is None
    assert JB.view({"bundle": {"status": "failed"}}, rna_path) is None
