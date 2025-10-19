"""
Utility to create a synthetic AnnData object for metacell unit testing.

The generated dataset includes multiple samples, cell types, and donors with
deterministic expression profiles that support exercising the metacell module's
different execution branches (global, boundary-aware, random aggregation, and
sample-restricted modes).
"""

from __future__ import annotations

import datetime as _dt
from pathlib import Path
from typing import Tuple

import anndata as ad
import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent / "data"
DEFAULT_PATH = DATA_DIR / "metacell_synthetic.h5ad"


def _base_expression_profiles(genes: np.ndarray) -> pd.DataFrame:
    """Return a DataFrame of archetypal expression patterns for cell types."""
    profiles = {
        "T_cell": np.linspace(2.0, 10.0, genes.size),
        "B_cell": np.linspace(8.0, 1.0, genes.size),
        "Myeloid": np.concatenate(
            [
                np.linspace(10.0, 4.0, genes.size // 2),
                np.linspace(4.0, 10.0, genes.size - genes.size // 2),
            ]
        ),
    }
    return pd.DataFrame(profiles, index=genes)


def build_synthetic_metacell_adata() -> ad.AnnData:
    """
    Create a deterministic synthetic AnnData object for metacell testing.

    The dataset includes:
    - 3 donors/samples (D1, D2, D3)
    - 3 cell types (T_cell, B_cell, Myeloid)
    - 72 total cells (8 replicates per sample/cell-type combination)
    - 15 genes with distinct expression programs per cell type
    """
    rng = np.random.default_rng(42)
    genes = np.array([f"Gene{i:02d}" for i in range(1, 16)])
    profiles = _base_expression_profiles(genes)

    samples = ["D1", "D2", "D3"]
    cell_types = ["T_cell", "B_cell", "Myeloid"]
    cells_per_combo = 8

    expression_rows = []
    obs_records = []
    var_records = [{"gene_name": g, "category": "protein_coding"} for g in genes]

    for sample_idx, sample in enumerate(samples):
        for cell_type_idx, cell_type in enumerate(cell_types):
            base_profile = profiles[cell_type]
            for replicate in range(cells_per_combo):
                noise = rng.normal(loc=0.0, scale=0.6, size=genes.size)
                counts = np.clip(base_profile + noise, a_min=0.0, a_max=None)
                cell_id = f"{sample}_{cell_type}_{replicate:02d}"
                expression_rows.append(counts)
                obs_records.append(
                    {
                        "cell_id": cell_id,
                        "sample_id": sample,
                        "donor_id": sample,
                        "cell_type": cell_type,
                        "batch": f"batch_{sample_idx % 2}",
                        "replicate": replicate,
                    }
                )

    X = np.vstack(expression_rows)
    obs = pd.DataFrame(obs_records).set_index("cell_id", drop=True)
    var = pd.DataFrame(var_records, index=genes)

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = (adata.X > 0).sum(axis=1)

    # Add a layer with raw counts scaled to integers for testing sum aggregation
    adata.layers["counts"] = np.round(adata.X * 10).astype(np.int32)

    # Create an alternative layer to test layer selection
    adata.layers["normalized"] = adata.X / np.maximum(adata.obs["n_counts"].to_numpy()[:, None], 1.0)

    adata.raw = adata.copy()
    return adata


def write_synthetic_h5ad(path: Path | None = None) -> Tuple[Path, ad.AnnData]:
    """
    Write the synthetic AnnData object to disk and return its path.

    Parameters
    ----------
    path:
        Optional custom path. Defaults to ``tests/data/metacell_synthetic.h5ad``.
    """
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    output_path = path or DEFAULT_PATH
    adata = build_synthetic_metacell_adata()
    adata.uns["generated_at"] = _dt.datetime.utcnow().isoformat()
    adata.write_h5ad(output_path, compression="gzip")
    return output_path, adata


if __name__ == "__main__":
    out_path, _adata = write_synthetic_h5ad()
    print(f"Synthetic metacell AnnData written to {out_path}")
