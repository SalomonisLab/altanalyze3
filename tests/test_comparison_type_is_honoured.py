"""Selecting "cells" runs a cell-level test, for every modality.

GRN edges used to register the (sample x cell-state) pseudobulk matrix under both
`h5ad` and `differential_h5ad`, and `_differential_runtime_params` returned fixed
pseudobulk settings for grn, grn_tf, metabolite and lipid whatever the reader chose.
A "cells" comparison therefore tested one value per state per group, every cell state
failed `min_cells_per_group`, and the run reported nothing.

Run:  python -m pytest tests/test_comparison_type_is_honoured.py -q
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from altanalyze3.components.cellHarmony.flask import pipeline  # noqa: E402

MODALITIES = ("rna", "adt", "lipids", "lipid", "metabolite", "grn", "grn_tf")


def test_every_modality_reads_comparison_type():
    """No modality may return the same settings for cells and pseudobulk."""
    for modality in MODALITIES:
        cells = pipeline._differential_runtime_params(modality, "cells")
        pseudobulk = pipeline._differential_runtime_params(modality, "pseudobulk")
        assert cells != pseudobulk, f"{modality} ignores comparison_type"


def test_cells_never_uses_the_pseudobulk_minimum():
    """A cell comparison must not fall back to a 2-replicate minimum."""
    for modality in MODALITIES:
        params = pipeline._differential_runtime_params(modality, "cells")
        assert params["min_cells_per_group"] >= 10, modality
        assert params["use_rawp"] is False, modality


def test_pseudobulk_keeps_samples_as_replicates():
    for modality in ("grn", "grn_tf", "metabolite", "lipid"):
        params = pipeline._differential_runtime_params(modality, "pseudobulk")
        assert params["min_cells_per_group"] == 2, modality
        assert params["use_rawp"] is True, modality


def test_the_matrix_follows_the_selection(tmp_path):
    """`pseudobulk_h5ad` is read only when the reader asked for a pseudobulk test."""
    per_cell = tmp_path / "grn_edges.h5ad"
    aggregate = tmp_path / "grn_edges_pseudobulk.h5ad"
    for path in (per_cell, aggregate):
        path.write_bytes(b"")
    meta = {"modality_artifacts": {"grn": {
        "h5ad": str(per_cell), "differential_h5ad": str(per_cell),
        "pseudobulk_h5ad": str(aggregate), "network_h5ad": str(aggregate)}}}

    assert pipeline._modality_differential_h5ad_path(meta, "grn", "cells") == per_cell
    assert pipeline._modality_differential_h5ad_path(meta, "grn", "pseudobulk") == aggregate
    # No selection given keeps the per-cell matrix, which is the documented default.
    assert pipeline._modality_differential_h5ad_path(meta, "grn") == per_cell


def test_grn_edges_no_longer_registers_the_aggregate_as_its_only_matrix():
    source = (REPO / "altanalyze3" / "components" / "cellHarmony" / "flask"
              / "pipeline.py").read_text()
    assert '"differential_h5ad": str(grn_edges_path)' in source
    assert '"pseudobulk_h5ad": str(grn_edges_pb_path)' in source
