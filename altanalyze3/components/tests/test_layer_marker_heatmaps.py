"""Layer plumbing with independent software values; no inference models run."""
from importlib import import_module
from io import StringIO

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from fastapi.testclient import TestClient

from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
from altanalyze3.components.cellHarmony.webapp.analysis_workflow import add_imputed_layer_marker_heatmaps
from altanalyze3.components.cellHarmony.webapp.app import create_app
from altanalyze3.components.cellHarmony.flask import pipeline as P

web = import_module("altanalyze3.components.cellHarmony.webapp.app")
discover = import_module("altanalyze3.components.cellHarmony.scalable_discover.pipeline")


def fixture(tmp_path):
    store = JobStore(tmp_path / "jobs")
    job = store.create_job("human", "fixture", None, [])["job_id"]
    obs = pd.DataFrame({"unsupervised_cluster": ["C2", "C1", "Not clustered"],
                        "unsupervised_state": ["state2", "state1", "Not clustered"]}, index=["a", "b", "c"])
    obj = ad.AnnData(np.array([[1., 2.], [3., 4.], [5., 6.]], dtype=np.float32), obs=obs,
                     var=pd.DataFrame(index=["output1", "output2"]))
    path = tmp_path / "prediction.h5ad"; obj.write_h5ad(path)
    result = {"job_id": job, "species": "human", "qc": {}, "cluster_key": "unsupervised_cluster",
              "artifacts": {"combined_h5ad": str(path)}, "modality_artifacts": {"adt": {"h5ad": str(path)}},
              "modalities": {"available": [{"id": "rna"}, {"id": "adt"}]},
              "marker_analysis": {"rna_tag": True}, "icgs3_analysis": {"clusters": ["C1", "C2"]},
              "cell_state_layers": {"names": {"C1": "state1", "C2": "state2"}, "layers": [
                  {"key": "unsupervised_cluster", "marker_analysis": {"rna_tag": True}},
                  {"key": "unsupervised_state", "marker_analysis": {"named_rna_tag": True}}]}}
    return store, job, result, obj, path


def test_standard_modality_hook_receives_all_predictions_and_native_order(tmp_path, monkeypatch):
    store, job, result, obj, source = fixture(tmp_path)
    original_bytes = source.read_bytes()
    calls = []
    def emit(modality, matrix, outputs, key, meta):
        assert modality == "adt" and key == "unsupervised_cluster"
        assert matrix.obs_names.equals(obj.obs_names) and matrix.var_names.equals(obj.var_names)
        np.testing.assert_array_equal(matrix.X, obj.X)
        assert list(matrix.obs[key].astype(str)) == ["C2", "C1", "Not clustered"]
        assert matrix.uns["lineage_order"] == ["C1", "C2", "Not clustered"]
        calls.append(key)
        return {"enabled": True, "modality": modality, "cluster_key": key}
    monkeypatch.setattr(P, "_emit_modality_marker_heatmap", emit)
    monkeypatch.setattr(discover, "relabel_marker_outputs", lambda analysis, names, directory: dict(analysis, names=names))
    add_imputed_layer_marker_heatmaps(job, store, result)
    assert calls == ["unsupervised_cluster"]  # names reuse the same fit
    layers = result["cell_state_layers"]["layers"]
    assert layers[0]["marker_analysis_by_modality"]["adt"]["cluster_key"] == "unsupervised_cluster"
    assert layers[1]["marker_analysis_by_modality"]["adt"]["cluster_key"] == "unsupervised_state"
    assert result["marker_analysis_by_modality"]["adt"]["enabled"]
    assert source.read_bytes() == original_bytes


def test_unknown_barcode_is_rejected_before_scoring(tmp_path, monkeypatch):
    store, job, result, obj, source = fixture(tmp_path)
    bad = obj.copy(); bad.obs_names = ["a", "b", "unreconciled"]
    badpath = tmp_path / "bad.h5ad"; bad.write_h5ad(badpath)
    result["modality_artifacts"]["adt"]["h5ad"] = str(badpath)
    monkeypatch.setattr(P, "_emit_modality_marker_heatmap", lambda *a: pytest.fail("must gate before scoring"))
    with pytest.raises(ValueError, match="identities"):
        add_imputed_layer_marker_heatmaps(job, store, result)


@pytest.mark.parametrize("cells_per_sample", [None, 10])
def test_heatmap_empty_filters_return_200_without_claiming_data(tmp_path, monkeypatch, cells_per_sample):
    app = create_app({"JOB_STORAGE": str(tmp_path / "jobs"), "ISOLATE_JOBS": False})
    store = app.state.job_store
    job = store.create_job("human", "fixture", None, [])["job_id"]
    expression = {"populations": np.array(["A"]), "obs_names": np.array(["a"])}
    base = {"row_ids": np.array(["A:g1"]), "signature": "fixture",
            "matrix": np.array([[1.]], dtype=np.float32), "col_ids": np.array(["A:a"]),
            "col_barcodes": np.array(["a"])}
    monkeypatch.setattr(web, "_get_expression_cache", lambda *a, **k: expression)
    monkeypatch.setattr(web, "_get_marker_heatmap_cache_entry", lambda *a, **k: base)
    monkeypatch.setattr(web, "_apply_display_filter_mask", lambda *a: np.array([False]))
    monkeypatch.setattr(web, "_sample_plot_cells", lambda cache, indices, limit: (indices, {"sample_field": "Library"}))
    try:
        with TestClient(app) as client:
            params = {"filter1_field": "state", "filter1_values": "absent"}
            if cells_per_sample is not None:
                params["cells_per_sample"] = cells_per_sample
            response = client.get(f"/api/jobs/{job}/marker/heatmap.tsv", params=params)
            assert response.status_code == 200
            assert response.headers["X-Marker-Columns"] == "0"
            if cells_per_sample is not None:
                assert response.headers["X-Marker-Rows"] == "1"
            table = pd.read_csv(StringIO(response.text), sep="\t", index_col=0)
            assert table.shape == (1, 0)
            assert list(table.index) == ["A:g1"]
    finally:
        app.state.job_runner.executor.shutdown(wait=True)
