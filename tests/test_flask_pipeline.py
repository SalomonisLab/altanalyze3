from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.cellHarmony.flask import create_app


def _write_reference(tmp_path: Path) -> dict:
    genes = ["GeneA", "GeneB", "GeneC"]
    states_df = pd.DataFrame(
        {
            "Pop1": [5, 2, 1],
            "Pop2": [1, 4, 3],
        },
        index=genes,
    )
    states_path = tmp_path / "demo_reference.tsv"
    states_df.to_csv(states_path, sep="\t")

    clusters_df = pd.DataFrame(
        {
            "CellBarcode": ["ref_cell_1", "ref_cell_2", "ref_cell_3", "ref_cell_4"],
            "Population": ["Pop1", "Pop2", "Pop1", "Pop2"],
        }
    )
    clusters_path = tmp_path / "demo_reference_clusters.tsv"
    clusters_df.to_csv(clusters_path, sep="\t", index=False)

    coords_df = pd.DataFrame(
        {
            "CellBarcode": clusters_df["CellBarcode"],
            "UMAP1": [0.1, -0.2, 0.3, -0.4],
            "UMAP2": [0.5, -0.6, 0.2, -0.1],
        }
    )
    coords_path = tmp_path / "demo_reference_coords.tsv"
    coords_df.to_csv(coords_path, sep="\t", index=False)

    return {
        "states_tsv": str(states_path),
        "reference_clusters_tsv": str(clusters_path),
        "reference_coords_tsv": str(coords_path),
        "cluster_key": "Population",
    }


def _write_query_h5ad(tmp_path: Path) -> Path:
    cells = [f"cell{i}" for i in range(1, 7)]
    genes = ["GeneA", "GeneB", "GeneC"]
    rng = np.random.default_rng(0)
    counts = rng.integers(0, 6, size=(len(cells), len(genes)))
    counts = counts.astype(float)
    adata = ad.AnnData(X=counts)
    adata.obs_names = cells
    adata.var_names = genes
    adata.obs["sample"] = ["Sample1"] * len(cells)
    adata.obs["group"] = ["Sample1"] * len(cells)
    adata.obs["Library"] = ["Sample1"] * len(cells)
    h5ad_path = tmp_path / "sample1.h5ad"
    adata.write(h5ad_path)
    return h5ad_path


def test_flask_pipeline_end_to_end(tmp_path):
    reference_meta = _write_reference(tmp_path)
    registry = {
        "species": [
            {
                "id": "demo_species",
                "label": "Demo Species",
                "references": [
                    {
                        "id": "demo_reference",
                        "label": "Demo Reference",
                        **reference_meta,
                    }
                ],
            }
        ]
    }
    registry_path = tmp_path / "references.json"
    registry_path.write_text(json.dumps(registry), encoding="utf-8")

    app = create_app(
        {
            "TESTING": True,
            "JOB_STORAGE": str(tmp_path / "jobs"),
            "REFERENCE_REGISTRY": str(registry_path),
            "SECRET_KEY": "test",
        }
    )
    client = app.test_client()

    h5ad_path = _write_query_h5ad(tmp_path)
    with h5ad_path.open("rb") as handle:
        data = {
            "species": "demo_species",
            "reference": "demo_reference",
            "sample_names[]": "Sample1",
            "files[]": (handle, "sample1.h5ad"),
        }
        resp = client.post("/api/jobs", data=data, content_type="multipart/form-data")
    assert resp.status_code == 200
    job_id = resp.json["job_id"]

    qc_payload = {"min_genes": 1, "min_counts": 0, "min_cells": 1, "mit_percent": 50}
    resp = client.post(f"/api/jobs/{job_id}/qc", json=qc_payload)
    assert resp.status_code == 200

    # Run pipeline synchronously via the JobRunner helper.
    app.job_runner._run_pipeline(job_id)

    resp = client.get(f"/api/jobs/{job_id}/status")
    assert resp.status_code == 200
    assert resp.json["status"] == "completed"

    resp = client.get(f"/api/jobs/{job_id}/umap")
    assert resp.status_code == 200
    umap_payload = resp.json
    assert umap_payload["query"], "Query UMAP points missing."
    assert umap_payload["reference"], "Reference UMAP points missing."

    resp = client.get(f"/api/jobs/{job_id}/expression", query_string={"gene": "GeneA"})
    assert resp.status_code == 200
    expr_payload = resp.json
    assert expr_payload["scatter"], "Expression scatter payload missing."
    assert expr_payload["violin"], "Expression violin payload missing."
