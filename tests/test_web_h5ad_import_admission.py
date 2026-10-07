"""The online web wrapper selects disk-backed import without ambient correction."""
from types import SimpleNamespace

import pytest

from altanalyze3.components.cellHarmony import mapped_h5ad
from altanalyze3.components.cellHarmony.flask import pipeline


@pytest.mark.parametrize("cells, matrix_bytes, bounded", [
    (100_000, 100, False), (100_001, 100, True),
    (280_446, 100, True), (4_000, 1024**3, True),
])
def test_web_selects_disk_import_with_ambient_off(tmp_path, monkeypatch, cells, matrix_bytes, bounded):
    logs = []
    store = SimpleNamespace(
        get_job=lambda _: {"species": "human", "reference": "test",
                           "files": [{"filename": "query.h5ad"}],
                           "qc": {"ambient_correction": "no"}},
        uploads_dir=lambda _: tmp_path, outputs_dir=lambda _: tmp_path / "outputs",
        append_log=lambda _, line: logs.append(line),
    )
    monkeypatch.setattr(pipeline, "_lookup_reference", lambda *args: {"states_tsv": "reference.tsv"})
    monkeypatch.setattr(pipeline, "_ensure_reference_fields", lambda _: None)
    monkeypatch.setattr(pipeline, "_first_reference_gene", lambda _: "gene")
    monkeypatch.setattr(pipeline, "_selected_impute_modalities", lambda *args: [])
    monkeypatch.setattr(mapped_h5ad, "inspect_h5ad", lambda _: {
        "cells": cells, "genes": 33_234, "matrix_bytes": matrix_bytes,
    })

    class ReachedImport(Exception):
        pass

    def capture(**kwargs):
        assert kwargs["bounded_h5ad"] is bounded
        assert kwargs["ambient_correct_cutoff"] is None
        assert kwargs["return_adata"] is True
        raise ReachedImport

    monkeypatch.setattr(pipeline.cellHarmony_lite, "combine_and_align_h5", capture)
    with pytest.raises(ReachedImport):
        pipeline.run_cellharmony_pipeline("test", store, tmp_path / "registry.json")
    backend = "disk-backed" if bounded else "in-memory"
    assert any(f"H5AD import backend={backend} cells={cells}" in line for line in logs)
