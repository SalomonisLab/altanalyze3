"""Exercise native rank/NMF behavior with independent synthetic expression."""
import json
from importlib import import_module

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

icgs = import_module("altanalyze3.components.clustering.ICGS")


def fixture():
    rng = np.random.default_rng(73)
    values = rng.uniform(0.01, 0.2, (24, 6)).astype(np.float32)
    values[:12, :3] += 4
    values[12:, 3:] += 4
    return AnnData(values, obs=pd.DataFrame(index=[f"cell{i}" for i in range(24)]),
                   var=pd.DataFrame(index=[f"feature{i}" for i in range(6)]))


def config(path, rank=None, **kwargs):
    return icgs.ICGS3Config(input_paths=[], output_dir=str(path), modality="lipid", rank=rank, **kwargs)


def records(path):
    return json.loads((path / "icgs3_nmf_rank_selection.json").read_text())


@pytest.mark.parametrize("requested", [0, 1, -2, 2.5, True, 100, 10])
def test_invalid_manual_k_matches_real_unbiased_result(tmp_path, requested, monkeypatch):
    monkeypatch.setenv("UDON_NMF_ENGINE", "sklearn")
    monkeypatch.delenv("UDON_NMF_RUNS", raising=False)
    matrix = fixture()
    before = matrix.copy()
    ordinary = tmp_path / "auto"
    recovered = tmp_path / "fallback"
    auto_expr, auto_clusters = icgs._fit_icgs3_nmf(matrix, config(ordinary), str(ordinary), "rowsum")
    manual = config(recovered, requested)
    expr, clusters = icgs._fit_icgs3_nmf(matrix, manual, str(recovered), "rowsum")
    pd.testing.assert_frame_equal(expr, auto_expr)
    pd.testing.assert_frame_equal(clusters, auto_clusters)
    np.testing.assert_array_equal(matrix.X, before.X)
    assert expr.index.equals(before.var_names) and clusters.index.equals(before.obs_names)
    assert manual.rank == requested  # requested configuration is not overwritten
    audit = records(recovered)
    assert [item["status"] for item in audit] == ["failed", "completed"]
    assert audit[1]["mode"] == "auto"
    assert audit[1]["final_rank"] == records(ordinary)[0]["final_rank"]


def test_successful_manual_k_does_not_invoke_estimator(tmp_path, monkeypatch):
    monkeypatch.setattr(icgs, "estimate_udon_nmf_rank", lambda *args: pytest.fail("Manual K should remain manual"))
    expr, clusters = icgs._fit_icgs3_nmf(fixture(), config(tmp_path, 2), str(tmp_path), "rowsum")
    assert expr.shape == (6, 24) and len(clusters) == 24
    assert len(records(tmp_path)) == 1 and records(tmp_path)[0]["final_rank"] == 2


def test_numerical_failure_reselects_features_with_auto_config(tmp_path, monkeypatch):
    seen = []
    def select(matrix, selected):
        seen.append(selected.rank)
        return matrix.var_names[:5] if selected.rank is not None else matrix.var_names
    real_fit = icgs.run_nmf
    def fit(frame, rank, **kwargs):
        if rank == 5 and frame.shape[0] == 5:
            raise icgs.NMFSelectionError("requested factorization was singular")
        return real_fit(frame, rank, **kwargs)
    monkeypatch.setattr(icgs, "_select_icgs3_features", select)
    monkeypatch.setattr(icgs, "run_nmf", fit)
    expr, clusters = icgs._fit_icgs3_nmf(fixture(), config(tmp_path, 5), str(tmp_path))
    assert seen == [5, None] and expr.shape[0] == 6 and len(clusters) == 24
    assert "singular" in records(tmp_path)[0]["error"]


@pytest.mark.parametrize("failure", [MemoryError("memory"), RuntimeError("unrelated"), ValueError("invalid input")])
def test_non_k_failure_is_not_retried(tmp_path, monkeypatch, failure):
    calls = []
    def fit(*args, **kwargs):
        calls.append(kwargs["rank"])
        raise failure
    monkeypatch.setattr(icgs, "run_nmf", fit)
    with pytest.raises(type(failure), match=str(failure)):
        icgs._fit_icgs3_nmf(fixture(), config(tmp_path, 2), str(tmp_path))
    assert calls == [2]


def test_auto_failure_is_not_retried_again_and_minimal_audit_exists(tmp_path, monkeypatch):
    monkeypatch.setattr(icgs, "run_nmf", lambda *args, **kwargs: (_ for _ in ()).throw(icgs.NMFSelectionError("fit failed")))
    with pytest.raises(icgs.NMFSelectionError, match="fit failed"):
        icgs._fit_icgs3_nmf(fixture(), config(tmp_path, 100, minimal_outputs=True), str(tmp_path))
    assert len(records(tmp_path)) == 2
    assert records(tmp_path)[1]["status"] == "failed"


def test_nimfa_fit_error_never_returns_expression_as_basis(tmp_path, monkeypatch):
    monkeypatch.setenv("UDON_NMF_ENGINE", "nimfa")
    def fail(*args, **kwargs):
        raise ValueError("invalid factorization")
    monkeypatch.setattr(icgs._icgs_nimfa, "Snmf", fail)
    frame = pd.DataFrame(fixture().X.T)
    with pytest.raises(icgs.NMFSelectionError, match="invalid factorization"):
        icgs.run_nmf(frame, rank=2)


def test_retry_releases_failed_solver_matrix_before_auto_fit(tmp_path, monkeypatch):
    import weakref
    seen = []
    real_fit = icgs.run_nmf
    def fit(frame, rank, **kwargs):
        if not seen:
            failed_solver_matrix = frame.to_numpy().copy()
            seen.append(weakref.ref(failed_solver_matrix))
            raise icgs.NMFSelectionError("numerical fit failed")
        assert seen[0]() is None, "Failed solver matrices must be freed before retrying"
        return real_fit(frame, rank, **kwargs)
    monkeypatch.setattr(icgs, "run_nmf", fit)
    _, clusters = icgs._fit_icgs3_nmf(fixture(), config(tmp_path, 2), str(tmp_path))
    assert len(clusters)==24 and records(tmp_path)[1]['status']=='completed'
