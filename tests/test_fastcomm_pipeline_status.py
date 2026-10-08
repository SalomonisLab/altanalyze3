"""Check workflow status plumbing with explicit independent inputs, without model runs."""
from types import SimpleNamespace
from importlib import import_module

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.cellHarmony.flask import pipeline as P
web = import_module("altanalyze3.components.cellHarmony.webapp.app")


@pytest.mark.parametrize("symbols", [False, True])
@pytest.mark.parametrize("split_fails", [False, True])
def test_optional_split_failure_preserves_global_scores(tmp_path, monkeypatch, symbols, split_fails):
    obj = ad.AnnData(np.ones((10, 2), dtype=np.float32),
                    obs=pd.DataFrame({"state": ["sender"] * 5 + ["receiver"] * 5,
                                      "Library": ["a", "b"] * 5}, index=[f"c{i}" for i in range(10)]),
                    var=pd.DataFrame(index=["ENSG1", "ENSG2"] if symbols else ["TGFB1", "TGFBR1"]))
    if symbols:
        obj.var["gene_symbols"] = ["TGFB1", "TGFBR1"]
    original = obj.copy()
    logs, calls = [], []
    store = SimpleNamespace(update_job=lambda *a, **k: None, append_log=lambda job, text: logs.append(text))
    expected_symbol_col = "gene_symbols" if symbols else None

    def verify(params):
        assert params.state_key == "state" and params.species == "human"
        assert params.gene_symbol_col == expected_symbol_col
        assert params.min_cells == 5 and params.min_lr_expression_score == .2
        assert params.max_lr_candidates_per_state_pair == 5 and not params.include_self_edges
        assert list(params.lr_sources) == ["CellChatDB"] and params.response_matrix is None

    def global_scores(params):
        verify(params)
        assert params.adata is obj
        calls.append("global")
        pd.DataFrame({"sender_state": ["sender"], "receiver_state": ["receiver"],
                      "ligand": ["TGFB1"], "receptor": ["TGFBR1"], "fastcomm_score": [.5]}).to_csv(params.output, sep="\t", index=False)
        params.state_pair_output.write_text("sender_state\treceiver_state\n")
        params.state_expression_output.write_text("gene\tsender\treceiver\n")
        return SimpleNamespace(state_sizes=pd.Series([5, 5], index=["sender", "receiver"]),
                               summary={"n_states": 2, "n_scored_edges": 1, "n_loaded_genes": 2})

    def splits(params):
        verify(params)
        assert params.h5ad == tmp_path / "input.h5ad" and params.split_key == "Library"
        calls.append("split")
        if split_fails:
            raise RuntimeError("explicit split failure")
        params.output_dir.mkdir(parents=True)
        path = params.output_dir / "split_scores_long.tsv"
        path.write_text("split\tsender_state\treceiver_state\n")
        return {"split_scores_long_tsv": str(path)}

    def report(params):
        params.output_tsv.write_text("ligand\treceptor\nTGFB1\tTGFBR1\n")
        params.output_md.write_text("test report")

    monkeypatch.setattr(P, "run_fastcomm", global_scores)
    monkeypatch.setattr(P, "run_fastcomm_benchmark", splits)
    monkeypatch.setattr(P, "write_exemplar_report", report)
    result, artifacts = P._run_fastcomm_analysis("fixture", store, obj, tmp_path / "input.h5ad", tmp_path, "state", {"species": "human"})
    assert result["enabled"] and result["status"] == "completed"
    assert artifacts["fastcomm_scores"].is_file() and artifacts["fastcomm_archive"].is_file()
    assert result["populations"] == ["sender", "receiver"] and result["sample_key"] == "Library"
    assert calls == ["global", "split"]
    np.testing.assert_array_equal(obj.X, original.X)
    pd.testing.assert_frame_equal(obj.obs, original.obs)
    pd.testing.assert_frame_equal(obj.var, original.var)
    if split_fails:
        assert result["per_sample"] == {"status": "failed", "message": "RuntimeError: explicit split failure"}
        assert any("global communication scores remain available" in log for log in logs)
    else:
        assert result["per_sample"]["split_scores_long_tsv"]


def test_global_failure_retains_specific_diagnostic(tmp_path, monkeypatch):
    logs = []
    store = SimpleNamespace(update_job=lambda *a, **k: None, append_log=lambda job, text: logs.append(text))
    obj = ad.AnnData(np.ones((1, 1)), var=pd.DataFrame(index=["TGFB1"]))
    def fail(params):
        raise RuntimeError("explicit global failure")
    monkeypatch.setattr(P, "run_fastcomm", fail)
    result, artifacts = P._run_fastcomm_analysis("fixture", store, obj, tmp_path / "input.h5ad", tmp_path, "state", {"species": "human"})
    assert result == {"enabled": False, "status": "failed", "message": "RuntimeError: explicit global failure", "state_key": "state"}
    assert not artifacts and any("explicit global failure" in log for log in logs)


@pytest.mark.parametrize("split_state", ["completed", "failed", "absent"])
def test_differential_and_split_views_require_successful_saved_splits(tmp_path, split_state):
    path = tmp_path / "split_scores_long.tsv"
    if split_state != "absent":
        path.write_text("split\tsender_state\treceiver_state\n")
    per_sample = {"status": split_state, "message": "explicit split failure", "split_scores_long_tsv": str(path)}
    meta = {"fastcomm_analysis": {"enabled": True, "status": "completed", "sample_key": "Library", "per_sample": per_sample},
            "differential_options": {"upload_profile": {}, "enabled": True}, "files": []}
    success = split_state == "completed"
    assert (web._fastcomm_split_scores_path(meta) is not None) == success
    assert any(m["id"] == "cell_communication" for m in web._differential_options(meta)["modalities"]) == success
    if split_state == "failed":
        with pytest.raises(ValueError, match="explicit split failure"):
            P._fastcomm_split_scores_long_path(meta)
    # A precomputed viewer uses its own already completed contrast, not split files.
    meta["scalable_viewer"] = {"deg_comparisons": [{"modality": "cell_communication"}]}
    assert any(m["id"] == "cell_communication" for m in web._differential_options(meta)["modalities"])
