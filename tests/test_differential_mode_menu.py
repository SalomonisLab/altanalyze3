"""A Differential Explorer view with no data must not be offered.

`supports_differential_network` / `supports_differential_go` are MODALITY capabilities,
not statements about one run. RNA keeps both, so the scALABLE-viewer menu listed
"Network" and "GO Terms" for COPD contrasts that ship neither, and choosing one printed
"No network data are available for this differential run."
"""
from importlib import import_module
from types import SimpleNamespace

import pytest

# `from ...webapp import app` binds the FastAPI instance, not the module.
W = import_module("altanalyze3.components.cellHarmony.webapp.app")


def meta_for(status, networks=(), go_populations=()):
    return {
        "cluster_key": "cell_type",
        "modalities": {"default": "rna",
                       "available": [dict(W._DEFAULT_MODALITY_DEFINITIONS["rna"])]},
        "differential": {
            "status": status,
            "run_id": "cancer_vs_no_cancer",
            "config": {"modality": "rna", "population_col": "cell_type"},
            "networks": [dict(id=n, population=n) for n in networks],
            "artifacts": {},
        },
        "_test_go_populations": list(go_populations),
    }


def labels(monkeypatch, meta):
    monkeypatch.setattr(W, "_differential_go_populations",
                        lambda app, m: list(m.get("_test_go_populations") or []))
    app = SimpleNamespace(state=SimpleNamespace(root_path=""))
    payload = W._build_differential_payload(app, "job", meta)
    return [entry["label"] for entry in payload["visualization_modes"]]


def test_completed_run_without_data_hides_both_views(monkeypatch):
    out = labels(monkeypatch, meta_for("completed"))
    assert "Network" not in out and "GO Terms" not in out
    # The views that always have data stay.
    assert "Volcano" in out and "Heatmap" in out and "Differential counts" in out


def test_completed_run_with_data_keeps_the_view(monkeypatch):
    out = labels(monkeypatch, meta_for("completed", networks=("AT1",), go_populations=("AT1",)))
    assert "Network" in out and "GO Terms" in out


def test_each_view_is_gated_on_its_own_data(monkeypatch):
    only_net = labels(monkeypatch, meta_for("completed", networks=("AT1",)))
    assert "Network" in only_net and "GO Terms" not in only_net
    only_go = labels(monkeypatch, meta_for("completed", go_populations=("AT1",)))
    assert "GO Terms" in only_go and "Network" not in only_go


@pytest.mark.parametrize("status", ["queued", "processing"])
def test_a_run_in_progress_keeps_its_entries(status, monkeypatch):
    """Artifacts arrive later; the menu must not flicker while the run computes."""
    out = labels(monkeypatch, meta_for(status))
    assert "Network" in out and "GO Terms" in out
