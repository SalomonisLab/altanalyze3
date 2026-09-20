"""The volcano still draws when a comparison carried no statistical test.

Cell communication with fewer than 2 samples a side reports an effect and no p-value.
The volcano keyed its y axis on FDR, fell back to p-value, found neither, dropped every
point and the panel read "No volcano data were found", which looks like a broken view
rather than an untested comparison.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import importlib  # noqa: E402

# `from ...webapp import app` binds the FastAPI instance, not the module.
web = importlib.import_module("altanalyze3.components.cellHarmony.webapp.app")


def _frame(**over):
    base = pd.DataFrame({
        "gene": ["CD69->KLRB1", "IL7->IL7R", "CCL5->CCR5"],
        "population": ["NK", "NK", "NK"],
        "log2fc": [-12.34, 3.2, 0.5],
        "pval": [np.nan, np.nan, np.nan],
        "fdr": [np.nan, np.nan, np.nan],
        "abs_delta_score": [0.51, 0.20, 0.02],
        "no_test_reason": ["fewer than 2 cell-communication samples in a group"] * 3,
    })
    for k, v in over.items():
        base[k] = v
    return base


def _payload(frame, monkeypatch):
    monkeypatch.setattr(web, "_get_differential_detail_table", lambda app, meta: frame)
    return web._build_differential_volcano_payload(None, {}, "NK")


def test_an_untested_comparison_still_returns_points(monkeypatch):
    out = _payload(_frame(), monkeypatch)
    assert out["untested"] is True
    assert out["statistic"] == "effect"
    assert len(out["points"]) == 3, "every interaction with an effect must plot"


def test_the_y_axis_is_labelled_as_an_effect_not_a_significance(monkeypatch):
    out = _payload(_frame(), monkeypatch)
    assert out["statistic_label"] == "effect size"
    assert "log" not in out["statistic_label"].lower()


def test_the_reason_reaches_the_browser(monkeypatch):
    out = _payload(_frame(), monkeypatch)
    assert "fewer than 2" in out["reason"]


def test_points_are_ordered_by_effect(monkeypatch):
    out = _payload(_frame(), monkeypatch)
    scores = [p["score"] for p in out["points"]]
    assert scores == sorted(scores, reverse=True)
    assert out["points"][0]["gene"] == "CD69->KLRB1"


def test_a_tested_comparison_is_unchanged(monkeypatch):
    frame = _frame(pval=[1e-8, 0.01, 0.4], fdr=[4.38e-8, 0.02, 0.5])
    out = _payload(frame, monkeypatch)
    assert out["untested"] is False
    assert out["statistic"] == "fdr"
    assert out["statistic_label"] == "FDR"
    assert len(out["points"]) == 3
    # -log10(4.38e-8) is about 7.36, the value the earlier volcano plotted.
    assert 7.0 < out["points"][0]["score"] < 7.7


def test_no_effect_and_no_statistic_gives_no_points(monkeypatch):
    frame = _frame(abs_delta_score=[np.nan] * 3, log2fc=[np.nan] * 3)
    out = _payload(frame, monkeypatch)
    assert out["points"] == []
