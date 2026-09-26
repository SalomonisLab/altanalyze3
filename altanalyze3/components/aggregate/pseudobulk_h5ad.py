#!/usr/bin/env python3
"""
Build a pseudobulk AnnData from a single-cell AnnData.

For each (cluster_col, sample_col) pair, this script:
- sums raw counts across all cells in the group
- stores raw summed counts in ``layers["counts"]``
- stores counts normalized by the total pseudobulk counts in ``X``
- carries forward sample-level obs annotations into the export ``.obs``

The output h5ad is written with fast ``lzf`` compression.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def _log_timing(label: str, started_at: float) -> float:
    elapsed = max(0.0, time.perf_counter() - started_at)
    print(f"[timing] {label}={elapsed:.2f}s")
    return elapsed


def _resolve_matrix(adata: ad.AnnData, count_layer: str | None) -> tuple[Any, str]:
    if count_layer and count_layer in adata.layers:
        return adata.layers[count_layer], count_layer
    if count_layer:
        raise KeyError(f"Requested count layer '{count_layer}' was not found in the AnnData object.")
    if "counts" in adata.layers:
        return adata.layers["counts"], "counts"
    return adata.X, "X"


def _to_csr(matrix: Any) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(np.asarray(matrix))


def _mode_or_unique(series: pd.Series) -> Any:
    non_null = series.dropna()
    if non_null.empty:
        return pd.NA
    unique = pd.unique(non_null)
    if len(unique) == 1:
        return unique[0]
    modes = non_null.mode(dropna=True)
    if not modes.empty:
        return modes.iloc[0]
    return non_null.iloc[0]


def _build_sample_metadata(
    obs: pd.DataFrame,
    *,
    sample_col: str,
    exclude_cols: list[str] | None = None,
) -> pd.DataFrame:
    exclude = set(exclude_cols or [])
    columns = [col for col in obs.columns if col not in exclude]
    grouped = obs[columns].groupby(obs[sample_col], sort=False, observed=False)
    sample_meta = grouped.agg(_mode_or_unique)
    sample_meta.insert(0, "sample_n_cells", grouped.size().astype(int))
    return sample_meta


def build_pseudobulk_h5ad(
    h5ad_path: Path,
    *,
    cluster_col: str,
    sample_col: str,
    output_h5ad: Path,
    min_cells: int = 10,
    count_layer: str | None = None,
    extra_layers: list[str] | None = None,
) -> ad.AnnData:
    total_started = time.perf_counter()

    started = time.perf_counter()
    adata = ad.read_h5ad(h5ad_path)
    _log_timing("pseudobulk.load_h5ad", started)
    print(f"[INFO] Loaded AnnData with {adata.n_obs} cells and {adata.n_vars} genes.")
    result = build_pseudobulk_from_adata(
        adata,
        cluster_col=cluster_col,
        sample_col=sample_col,
        output_h5ad=output_h5ad,
        min_cells=min_cells,
        count_layer=count_layer,
        extra_layers=extra_layers,
    )
    _log_timing("pseudobulk.total", total_started)
    return result


def build_pseudobulk_from_adata(
    adata: ad.AnnData,
    *,
    cluster_col: str,
    sample_col: str,
    output_h5ad: Path,
    min_cells: int = 10,
    count_layer: str | None = None,
    extra_layers: list[str] | None = None,
) -> ad.AnnData:
    """Same aggregation as ``build_pseudobulk_h5ad`` on an AnnData already in memory, so a
    caller such as cellHarmony never writes the single-cell object. Each name in
    ``extra_layers`` is summed with the same groups into a pseudobulk layer of that name."""
    extra_layers = [name for name in (extra_layers or []) if name]
    missing_layers = [name for name in extra_layers if name not in adata.layers]
    if missing_layers:
        raise KeyError(f"Requested extra layers were not found in the AnnData object: {missing_layers}")

    if cluster_col not in adata.obs.columns:
        raise KeyError(f"Cluster column '{cluster_col}' was not found in .obs.")
    if sample_col not in adata.obs.columns:
        raise KeyError(f"Sample column '{sample_col}' was not found in .obs.")

    started = time.perf_counter()
    matrix, matrix_name = _resolve_matrix(adata, count_layer)
    counts = _to_csr(matrix)
    print(f"[INFO] Using '{matrix_name}' as the pseudobulk count source.")
    _log_timing("pseudobulk.resolve_counts", started)

    started = time.perf_counter()
    obs = adata.obs.copy()
    cluster_values = obs[cluster_col].astype(str)
    sample_values = obs[sample_col].astype(str)
    group_ids = cluster_values + "|" + sample_values
    codes, unique_groups = pd.factorize(group_ids, sort=False)
    group_sizes = np.bincount(codes)
    valid_group_mask = group_sizes >= int(min_cells)
    n_valid_groups = int(valid_group_mask.sum())
    if n_valid_groups == 0:
        raise ValueError(f"No pseudobulk groups passed min_cells={min_cells}.")
    print(f"[INFO] Retained {n_valid_groups} pseudobulk groups with min_cells={min_cells}.")

    valid_cell_mask = valid_group_mask[codes]
    remap = np.full(len(valid_group_mask), -1, dtype=np.int64)
    remap[np.flatnonzero(valid_group_mask)] = np.arange(n_valid_groups, dtype=np.int64)
    remapped_codes = remap[codes[valid_cell_mask]]
    _log_timing("pseudobulk.build_groups", started)

    started = time.perf_counter()
    design = sparse.csr_matrix(
        (
            np.ones(remapped_codes.shape[0], dtype=np.float32),
            (remapped_codes, np.arange(remapped_codes.shape[0], dtype=np.int64)),
        ),
        shape=(n_valid_groups, remapped_codes.shape[0]),
    )
    retained_counts = counts[valid_cell_mask]
    pseudobulk_counts = (design @ retained_counts).tocsr()
    retained_sum = float(np.sum(retained_counts.data, dtype=np.float64))
    del retained_counts
    pseudobulk_sum = float(np.sum(pseudobulk_counts.data, dtype=np.float64))
    if not np.isclose(pseudobulk_sum, retained_sum, rtol=1e-6, atol=1e-6):
        raise AssertionError(
            f"Pseudobulk counts sum to {pseudobulk_sum:.6g} but the retained cells sum to "
            f"{retained_sum:.6g}; the aggregation lost or invented counts.")
    extra_sums = {}
    for name in extra_layers:
        layer_counts = _to_csr(adata.layers[name])
        extra_sums[name] = (design @ layer_counts[valid_cell_mask]).tocsr()
        if not np.isclose(float(np.sum(extra_sums[name].data, dtype=np.float64)),
                          float(np.sum(layer_counts[valid_cell_mask].data, dtype=np.float64)),
                          rtol=1e-6, atol=1e-6):
            raise AssertionError(f"Pseudobulk layer '{name}' does not conserve its cell counts.")
    print(f"[INFO] Pseudobulks keep {int(valid_cell_mask.sum())} of {adata.n_obs} cells; "
          f"'{matrix_name}' sum conserved ({pseudobulk_sum:.6g}).")
    total_counts = np.asarray(pseudobulk_counts.sum(axis=1)).ravel().astype(np.float64)
    inv_totals = np.divide(
        1.0,
        total_counts,
        out=np.zeros_like(total_counts, dtype=np.float64),
        where=total_counts > 0,
    )
    pseudobulk_normalized = sparse.diags(inv_totals) @ pseudobulk_counts
    _log_timing("pseudobulk.aggregate_counts", started)

    started = time.perf_counter()
    valid_group_labels = unique_groups[np.flatnonzero(valid_group_mask)].astype(str)
    valid_group_sizes = group_sizes[np.flatnonzero(valid_group_mask)].astype(int)
    group_cluster = pd.Index(valid_group_labels).str.split("|", n=1).str[0]
    group_sample = pd.Index(valid_group_labels).str.split("|", n=1).str[1]

    sample_metadata = _build_sample_metadata(
        obs,
        sample_col=sample_col,
        exclude_cols=[cluster_col, sample_col],
    )
    pseudobulk_obs = pd.DataFrame(
        {
            "pseudobulk_id": valid_group_labels,
            cluster_col: group_cluster.to_numpy(dtype=object),
            sample_col: group_sample.to_numpy(dtype=object),
            "n_cells": valid_group_sizes,
            "total_counts": total_counts,
        },
        index=valid_group_labels,
    )
    pseudobulk_obs = pseudobulk_obs.join(sample_metadata, on=sample_col, rsuffix="_sample")
    pseudobulk_obs.index.name = None
    _log_timing("pseudobulk.build_obs", started)

    started = time.perf_counter()
    pseudobulk_adata = ad.AnnData(
        X=pseudobulk_normalized.tocsr(),
        obs=pseudobulk_obs,
        var=adata.var.copy(),
        uns=dict(adata.uns),
    )
    pseudobulk_adata.layers["counts"] = pseudobulk_counts
    for name, summed in extra_sums.items():
        pseudobulk_adata.layers[name] = summed
    pseudobulk_adata.uns["pseudobulk"] = {
        "cluster_col": cluster_col,
        "sample_col": sample_col,
        "min_cells": int(min_cells),
        "count_source": matrix_name,
        "normalization": "counts divided by total counts per pseudobulk",
    }
    if extra_sums:
        pseudobulk_adata.uns["pseudobulk"]["extra_layers"] = list(extra_sums)
    _log_timing("pseudobulk.build_h5ad", started)

    started = time.perf_counter()
    output_h5ad = Path(output_h5ad)
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    pseudobulk_adata.write_h5ad(output_h5ad, compression="lzf")
    _log_timing("pseudobulk.write_h5ad", started)
    print(f"[INFO] Wrote pseudobulk h5ad to: {output_h5ad}")
    return pseudobulk_adata


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a pseudobulk h5ad with raw summed counts and total-normalized expression."
    )
    parser.add_argument("--h5ad", required=True, help="Input single-cell h5ad file.")
    parser.add_argument("--cluster-col", required=True, help="obs column defining cell populations.")
    parser.add_argument("--sample-col", required=True, help="obs column defining donor/sample groups.")
    parser.add_argument("--output-h5ad", required=True, help="Output pseudobulk h5ad path.")
    parser.add_argument(
        "--min-cells",
        type=int,
        default=10,
        help="Minimum number of cells required for a pseudobulk group (default: 10).",
    )
    parser.add_argument(
        "--count-layer",
        default=None,
        help="Layer containing raw counts. Defaults to .layers['counts'] when present, else .X.",
    )
    parser.add_argument(
        "--extra-layers",
        default="",
        help="Comma-separated layers to sum into pseudobulk layers of the same name.",
    )
    args = parser.parse_args()

    build_pseudobulk_h5ad(
        Path(args.h5ad),
        cluster_col=args.cluster_col,
        sample_col=args.sample_col,
        output_h5ad=Path(args.output_h5ad),
        min_cells=args.min_cells,
        count_layer=args.count_layer,
        extra_layers=[name.strip() for name in args.extra_layers.split(",") if name.strip()],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
