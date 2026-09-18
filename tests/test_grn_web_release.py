"""GRN release contract: independent inputs, real differentials, and shared views.

Uses synthetic uploaded-job artifacts only; never reruns a published analysis.
"""
from pathlib import Path
from types import SimpleNamespace
import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from fastapi.testclient import TestClient

from altanalyze3.components.cellHarmony import grn_analysis as grn
from altanalyze3.components.cellHarmony.flask import pipeline
from importlib import import_module
web = import_module("altanalyze3.components.cellHarmony.webapp.app")
from altanalyze3.components.cellHarmony.webapp.grn_data import UploadedGrnData
from altanalyze3.components.rna2grn.api import Rna2GrnBundle


def test_activity_sum_preserves_default_mean():
    bundle = object.__new__(Rna2GrnBundle)
    bundle.metadata = {"edge_tf": ["TF1", "TF1", "TF2"]}
    predictions = pd.DataFrame([[1., 3., 8.]], columns=["TF1|A", "TF1|B", "TF2|C"])
    assert bundle.tf_activity(predictions).iloc[0].to_dict() == {"TF1": 2., "TF2": 8.}
    assert bundle.tf_activity(predictions, aggregation="sum").iloc[0].to_dict() == {"TF1": 4., "TF2": 8.}
    with pytest.raises(ValueError):
        bundle.tf_activity(predictions, aggregation="median")


def test_sibling_matching_from_every_modality():
    comparisons = [dict(id=f"{m}::copd::case_vs_control::per_cell_state", modality=m,
                        comparison="case_vs_control", contrast="copd", kind="per_cell_state")
                   for m in ("rna", "grn", "grn_tf")]
    comparisons.append(dict(comparisons[2], id="wrong", contrast="other_groups"))
    ds = SimpleNamespace(deg_manifest=lambda: {"comparisons": comparisons})
    for source in comparisons[:3]:
        for target in comparisons[:3]:
            assert grn._sibling_comparison(ds, source["id"], target["modality"]) == target["id"]
    assert grn._sibling_comparison(ds, "unknown", "grn_tf") == ""


def test_viewer_chat_does_not_substitute_unrelated_tf_comparison():
    from fastapi import FastAPI
    from altanalyze3.components.visualization.scalable_viewer.scalable_app import _install_chat_routes
    comparisons = [dict(id="selected-rna", modality="rna", comparison="case_vs_control",
                        contrast="copd", kind="per_cell_state"),
                   dict(id="unrelated-tf", modality="grn_tf", comparison="male_vs_female",
                        contrast="sex", kind="per_cell_state")]
    ds = SimpleNamespace(states=["AT1"], symbols=["TF1"],
                         deg_manifest=lambda: {"comparisons": comparisons})
    store = SimpleNamespace(dataset=lambda _: ds,
                            get_job=lambda _: {"differential": {"run_id": "selected-rna"}})
    app = FastAPI()
    _install_chat_routes(app, store, {})
    response = TestClient(app).post("/api/jobs/test/chat", json={"question": "Which TF activity changed in AT1?"})
    assert response.status_code == 200
    assert response.json()["status"] == "not_run"
    assert "plot" not in response.json()


def query_data():
    rng = np.random.default_rng(12)
    obs = pd.DataFrame([
        {"Library": f"donor{donor}", "condition": "case" if donor < 4 else "control",
         "cell_type": state}
        for donor in range(8) for state in ("AT1", "AT2") for _ in range(12)
    ])
    obs.index = [f"cell{i}" for i in range(len(obs))]
    x = rng.uniform(1, 2, (len(obs), 5)).astype(np.float32)
    x[obs.condition == "case", 0] += 4
    a = ad.AnnData(x, obs=obs, var=pd.DataFrame(index=["TF1", "TF2", "TARGET1", "TARGET2", "TARGET3"]))
    a.layers["counts"] = x.copy()
    a.obsm["X_umap"] = rng.normal(size=(len(obs), 2))
    a.uns["lineage_order"] = ["AT1", "AT2"]
    return a


class SmallGrnModel:
    """Deterministic predictor; exercises the real pseudobulk and chunk wiring."""
    pseudobulk_statistic = "mean_over_cells_of_log1p_cp10k"
    metadata = {"edge_tf": ["TF1", "TF1", "TF2"], "reference": "synthetic"}
    tf_activity = Rna2GrnBundle.tf_activity

    def __init__(self):
        self.batch_sizes = []

    def predict_from_adata(self, a, *, groupby, **kwargs):
        frame = pd.DataFrame(np.asarray(a.X), index=a.obs_names)
        if groupby:
            frame = frame.groupby(a.obs[groupby], sort=False).mean()
        else:
            self.batch_sizes.append(a.n_obs)
        values = frame.to_numpy()
        edges = pd.DataFrame(np.column_stack([values[:, 0] / 10, values[:, 0] / 20, values[:, 1] / 10]),
                             index=frame.index, columns=["TF1|TARGET1", "TF1|TARGET2", "TF2|TARGET3"])
        return SimpleNamespace(predictions=edges, summary={})


@pytest.fixture
def uploaded(tmp_path, monkeypatch):
    a = query_data()
    bundle = SmallGrnModel()
    monkeypatch.setattr(pipeline, "load_rna2grn_bundle", lambda *args: bundle)
    tf, edges, tf_pb, summary = pipeline._build_imputed_grn_adata(
        a, {"impute_config": {"grn": {"tf_activity_chunk_size": 32}}}, "cell_type")
    app = web.create_app({"JOB_STORAGE": str(tmp_path / "jobs"), "JOB_WORKERS": 1})
    store = app.state.job_store
    client = TestClient(app)
    query_path = tmp_path / "query.h5ad"
    a.write_h5ad(query_path)
    response = client.post("/api/jobs", data={"species": "human", "reference": "synthetic", "sample_names": "query"},
                           files={"files": ("query.h5ad", query_path.read_bytes(), "application/octet-stream")})
    assert response.status_code == 200, response.text
    job = response.json()["job_id"]
    root = store.outputs_dir(job)
    paths = {}
    for name, matrix in (("rna", a), ("grn", edges), ("grn_tf", tf), ("grn_tf_pb", tf_pb)):
        path = root / f"{name}.h5ad"
        matrix.write_h5ad(path)
        paths[name] = str(path)
    coords = root / "coords.tsv"
    pd.DataFrame({"CellBarcode": a.obs_names, "UMAP1": a.obsm["X_umap"][:, 0],
                  "UMAP2": a.obsm["X_umap"][:, 1]}).to_csv(coords, sep="\t", index=False)
    artifacts = {"rna": {"h5ad": paths["rna"]},
                 "grn": {"h5ad": paths["grn"], "differential_h5ad": paths["grn"], "network_h5ad": paths["grn"]},
                 "grn_tf": {"h5ad": paths["grn_tf"], "differential_h5ad": paths["grn_tf_pb"]}}
    store.update_job(job, status="completed", cluster_key="cell_type", reference_cluster_key="cell_type",
                     artifacts={"combined_h5ad": paths["rna"], "umap_coordinates": str(coords)},
                     modality_artifacts=artifacts, modalities=pipeline._modalities_payload(["grn"]),
                     differential_options={"enabled": True, "default_population_col": "cell_type",
                                           "population_columns": [{"value": "cell_type", "label": "cell_type"}],
                                           "default_sample_field": "condition"})
    yield SimpleNamespace(app=app, store=store, job=job, root=root, a=a, tf=tf, edges=edges,
                          tf_pb=tf_pb, bundle=bundle, paths=paths, client=TestClient(app))
    app.state.job_runner.executor.shutdown(wait=True)


def seed_comparisons(u):
    runs = {}
    features = {"rna": ["TARGET1", "TF1"], "grn": ["TF1|TARGET1"], "grn_tf": ["TF1"]}
    for modality, genes in features.items():
        root = u.root / modality
        root.mkdir(exist_ok=True)
        rows = [dict(gene=g, population=s, log2fc=1.5, pval=0.001, fdr=0.01,
                     case_label="case", control_label="control", n_case=4, n_control=4)
                for g in genes for s in ("AT1", "AT2")]
        table = root / "DEG_detailed_case_vs_control.tsv"
        pd.DataFrame(rows).to_csv(table, sep="\t", index=False)
        fold = root / "fold.tsv"
        pd.DataFrame(1.5, index=genes, columns=["AT1", "AT2"]).to_csv(fold, sep="\t")
        run_id = f"run-{modality}"
        runs[run_id] = dict(status="completed", run_id=run_id, case_label="case", control_label="control",
                            comparison_tag="case_vs_control", feature_label="factor" if modality == "grn_tf" else "gene",
                            config=dict(modality=modality, population_col="cell_type", sample_field="condition",
                                        group1_samples=["case"], group2_samples=["control"], comparison_type="pseudobulk"),
                            artifacts={table.stem: str(table), "fold_matrix_tsv": str(fold)})
    u.store.update_job(u.job, differential=runs["run-grn"], differential_history=runs)
    return runs


def test_uploaded_imputation_separates_cells_and_replicates(uploaded):
    u = uploaded
    assert max(u.bundle.batch_sizes) == 32
    assert u.tf.shape == (192, 2)
    assert u.tf_pb.shape == (16, 2)
    assert u.edges.shape == (16, 3)
    np.testing.assert_allclose(u.tf_pb.X[:, 0], u.edges.X[:, :2].sum(axis=1))
    assert np.unique(u.tf.X[:, 0]).size > 100
    assert u.tf.uns["modality"] == "grn_tf"
    assert u.tf_pb.uns["pseudobulk_method"] == "pseudobulk"
    assert set(u.tf_pb.obs.Library) == {f"donor{i}" for i in range(8)}
    meta = u.store.get_job(u.job)
    assert pipeline._modality_differential_h5ad_path(meta, "tf_activity") == Path(u.paths["grn_tf_pb"])
    assert pipeline._modality_differential_h5ad_path(meta, "grn_edges") == Path(u.paths["grn"])


def test_grn_edges_default_factor_state_and_honest_colours(uploaded):
    u = uploaded
    url = f'/api/jobs/{u.job}/grn/network'
    r = u.client.get(url)
    assert r.status_code == 200, r.text
    data = r.json()
    assert data['cell_state'] == 'AT1' and data['genes'] == ['TF1']
    assert data['n_edges'] == 2 and data['node_encoding'] == 'role'
    assert data['n_matching_edges'] == 2 and data['max_edges'] == 25
    assert data['n_aggregates'] == 8 and data['sample_field'] == 'Library'
    limited = u.client.get(url, params={'max_edges': 1}).json()
    assert limited['n_edges'] == 1 and limited['n_matching_edges'] == 2
    strongest = next(e['data'] for e in limited['elements'] if 'source' in e['data'])
    assert strongest['target'] == 'TARGET1'
    expected = np.asarray(u.edges[u.edges.obs.cell_type == 'AT1'].X.mean(axis=0)).ravel()[0]
    assert strongest['score'] == pytest.approx(expected)
    for invalid in (0, -1, 1001):
        assert u.client.get(url, params={'max_edges': invalid}).status_code == 422
    nodes = [e['data'] for e in data['elements'] if 'source' not in e['data']]
    assert all('log2fc' not in n for n in nodes)
    assert {n['role'] for n in nodes} == {'tf', 'target'}
    assert all(e['data']['direction'] == 'neutral' for e in data['elements'] if 'source' in e['data'])
    selected = u.client.get(url, params={'genes': 'TF2', 'cell_state': 'AT2', 'sample': 'donor7'}).json()
    assert selected['genes'] == ['TF2'] and selected['cell_state'] == 'AT2'
    assert selected['n_edges'] == 1
    assert selected['n_aggregates'] == 1
    assert not u.client.get(url, params={'genes': 'TF1', 'threshold': 100}).json()['elements']
    assert not u.client.get(url, params={'cell_state': 'unknown'}).json()['elements']


def test_uploaded_regulatory_endpoints_and_assistant(uploaded, monkeypatch):
    u = uploaded
    seed_comparisons(u)
    base = f"/api/jobs/{u.job}"
    examples = u.client.get(base + "/chat-examples")
    assert examples.status_code == 200, examples.text
    assert any("TFs" in text for text in examples.json()["examples"])
    profile = u.client.get(base + "/grn/tf-activity", params={"cell_state": "AT1"})
    assert profile.status_code == 200
    row = next(r for r in profile.json()["rows"] if r["factor"] == "TF1")
    assert row["log2fc"] == 1.5
    net = u.client.get(base + "/grn/regulator-network", params={"cell_state": "AT1", "features": "TARGET1"})
    assert net.status_code == 200
    assert net.json()["edges"]
    tf = next(n for n in net.json()["nodes"] if n["id"] == "TF1")
    assert tf["activity_log2fc"] == 1.5 and tf["expression_log2fc"] == 1.5
    assert "uploaded cells" in net.json()["arm_scope"]
    assert u.client.get(base + "/grn/tf-activity", params={"contrast": "bad"}).status_code == 404
    monkeypatch.setattr(web.urllib.request, "urlopen", lambda *a, **k: pytest.fail("Regulatory chat must not need the external router"))
    for question, kind in [("Which TFs are most active in AT1 cells?", "barchart"),
                           ("Show regulatory network for TARGET1 in AT1", "integrated_network"),
                           ("Which TF activity changed in AT1?", "volcano")]:
        response = u.client.post(base + "/chat", json={"question": question})
        assert response.status_code == 200, response.text
        assert response.json()["plot"]["kind"] == kind
        if kind == "volcano":
            assert response.json()["table"]["rows"][0][0] == "TF1"
    assert u.store.get_job(u.job)["differential"]["run_id"] == "run-grn"


def test_uploaded_existing_comparison_visualizations(uploaded):
    u = uploaded
    seed_comparisons(u)
    base = f"/api/jobs/{u.job}"
    for modality in ("grn_tf", "grn", "grn_tf"):
        response = u.client.post(base + "/differential/select", params={"contrast": f"run-{modality}"})
        assert response.status_code == 200, response.text
        assert response.json()["selected_modality"] == modality
        assert len(response.json()["completed_comparisons"]) == 3
        for plot, key in (("volcano", "points"), ("summary", "rows"), ("heatmap", "rows")):
            response = u.client.get(base + f"/differential/interactive/{plot}", params={"population": "AT1"})
            assert response.status_code == 200, response.text
            assert response.json()[key]
        gene = "TF1" if modality == "grn_tf" else "TF1|TARGET1"
        response = u.client.get(base + "/differential/interactive/gene", params={"population": "AT1", "gene": gene})
        assert response.status_code == 200, response.text
        response = u.client.get(base + "/differential/interactive/gene/pdf", params={"population": "AT1", "gene": gene})
        assert response.status_code == 200 and response.content.startswith(b"%PDF")
    response = u.client.get(base + "/expression", params={"gene": "TF1", "modality": "grn_tf"})
    assert response.status_code == 200 and response.json()["scatter"]
    for route in ("umap", "dotplot", "combplot", "expression/pdf", "plot-variables", "genes"):
        response = u.client.get(base + "/" + route, params={"gene": "TF1", "genes": "TF1,TF2", "modality": "grn_tf"})
        assert response.status_code == 200, (route, response.text)
        if route.endswith("/pdf"):
            assert response.content.startswith(b"%PDF")


def test_real_tf_and_edge_differentials(uploaded):
    u = uploaded
    runner = u.app.state.job_runner
    for modality in ("grn", "grn_tf"):
        u.store.update_job(u.job, differential={"status": "queued", "config": {
            "modality": modality, "population_col": "cell_type", "sample_field": "condition",
            "group1_samples": ["case"], "group2_samples": ["control"], "comparison_type": "pseudobulk"}})
        runner._run_differential(u.job)
        run = u.store.get_job(u.job)["differential"]
        assert run["status"] == "completed", run.get("message")
        detail = web._get_differential_detail_table(u.app, u.store.get_job(u.job))
        assert not detail.empty
        assert set(detail.n_case) == {4} and set(detail.n_control) == {4}
        assert all(("|" in name) == (modality == "grn") for name in detail.gene)
        assert Path(run["artifacts"]["archive"]).is_file()
    assert len(u.store.get_job(u.job)["differential_history"]) == 2


def test_legacy_enrichment_is_not_predicted_activity(uploaded):
    u = uploaded
    meta = u.store.get_job(u.job)
    meta["modality_artifacts"] = {"grn": {"h5ad": u.paths["grn_tf"], "differential_h5ad": u.paths["grn"]}}
    meta["modalities"]["available"] = [dict(id="rna"), dict(id="grn", label="GRN (TF activity)")]
    assert web._modality_h5ad_path(meta, "grn") == Path(u.paths["grn"])
    assert web._modality_h5ad_path(meta, "grn_tf") == Path(u.paths["grn_tf"])
    assert "grn_tf" not in [entry["id"] for entry in web._differential_options(meta)["modalities"]]
    with pytest.raises(ValueError, match="legacy"):
        pipeline._modality_differential_h5ad_path(meta, "grn_tf")


@pytest.mark.parametrize("render_heatmap", [True, False])
def test_upload_alignment_registers_both_grn_modalities(tmp_path, monkeypatch, render_heatmap):
    from test_flask_pipeline import _write_reference, _write_query_h5ad
    reference = _write_reference(tmp_path)
    registry = tmp_path / "references.json"
    registry.write_text(json.dumps({"species": [{"id": "demo_species", "label": "Demo",
        "references": [{"id": "demo_reference", "label": "Demo", **reference,
                        "impute_modalities": ["grn"]}]}]}))
    monkeypatch.setattr(pipeline, "load_rna2grn_bundle", lambda *args: SmallGrnModel())
    app = web.create_app({"JOB_STORAGE": str(tmp_path / "jobs"), "REFERENCE_REGISTRY": str(registry)})
    client = TestClient(app)
    query = _write_query_h5ad(tmp_path)
    response = client.post("/api/jobs", data={"species": "demo_species", "reference": "demo_reference", "sample_names": "Sample1"},
                           files={"files": ("sample1.h5ad", query.read_bytes(), "application/octet-stream")})
    assert response.status_code == 200, response.text
    job = response.json()["job_id"]
    store = app.state.job_store
    response = client.post(f"/api/jobs/{job}/qc", json={
        "min_genes": 1, "min_counts": 0, "min_cells": 1, "mit_percent": 50,
        "impute_modalities": ["grn"], "marker_render_heatmap": render_heatmap,
        "marker_heatmap_dpi": 60})
    assert response.status_code == 200, response.text
    try:
        app.state.job_runner._run_pipeline(job)
        meta = store.get_job(job)
        assert meta["status"] == "completed", meta.get("message")
        assert {"grn", "grn_tf"} <= set(meta["modality_artifacts"])
        edges = ad.read_h5ad(meta["modality_artifacts"]["grn"]["differential_h5ad"])
        factors = ad.read_h5ad(meta["modality_artifacts"]["grn_tf"]["differential_h5ad"])
        assert all("|" in feature for feature in edges.var_names)
        assert all("|" not in feature for feature in factors.var_names)
        assert "grn_tf" in meta["marker_analysis_by_modality"]
        for modality in ("rna", "grn_tf"):
            marker = meta["marker_analysis_by_modality"][modality]
            if not marker.get("enabled"):
                continue
            assert marker['output_options']['render_heatmap'] == render_heatmap
            assert bool(marker.get('heatmap_pdf')) == render_heatmap
            if render_heatmap:
                assert Path(marker['heatmap_pdf']).is_file()
            # Interactive data and PDF downloads must work with no static figure.
            response = client.get(f"/api/jobs/{job}/marker/heatmap.tsv?modality={modality}")
            assert response.status_code == 200, response.text
            response = client.get(f"/api/jobs/{job}/marker/heatmap.pdf?modality={modality}")
            assert response.status_code == 200, response.text
            assert response.content.startswith(b'%PDF')
        for key in ("imputed_grn_results_zip", "imputed_grn_tf_results_zip"):
            assert Path(meta["artifacts"][key]).is_file()
    finally:
        app.state.job_runner.executor.shutdown(wait=True)
