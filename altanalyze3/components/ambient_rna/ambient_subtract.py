#!/usr/bin/env python3
"""
Fast ambient RNA correction for filtered-only .h5ad inputs.

This file is intended as a drop-in replacement for the previous SoupX-based
wrapper in altanalyze3/components/ambient_rna/soupx_correct.py.

Design goals
------------
1. Keep the same public entry points:
   - process_h5ad
   - process_anndata
   - run_soupx_correction
   - main
2. Preserve the same default output layers:
   - soupx_raw
   - soupx_corrected
3. Be much faster on filtered-only AnnData objects by:
   - removing the soupx dependency
   - avoiding per-library clustering as part of correction
   - using direct sparse subtraction in CSC format
   - avoiding AnnData concat for reconstruction

Important
---------
This is NOT a full reimplementation of canonical SoupX.
Canonical SoupX expects raw droplets plus filtered cells.
This replacement is an intentionally fast approximation for the
filtered-only case already being used in cellHarmony.

Correction model
----------------
For each library:
- estimate an ambient profile b_g from the lowest-UMI cells
- for each cell j with total count n_j and contamination rho_j:
    ambient_jg = rho_j * n_j * b_g
    corrected_jg = max(0, count_jg - ambient_jg)

Only existing nonzero entries are touched, so runtime is O(nnz).
"""

from __future__ import annotations

import argparse
import csv
import os
import time
import warnings
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


DEFAULT_RHO = 0.15
RAW_LAYER_NAME = "soupx_raw"
CORRECTED_LAYER_NAME = "soupx_corrected"


def validate_rho(value: float) -> float:
    """Ensure rho stays within 0-1 inclusive."""
    if not np.isfinite(value):
        raise ValueError("Contamination fraction (rho) must be finite.")
    if value < 0.0 or value > 1.0:
        raise ValueError("Contamination fraction (rho) must be between 0 and 1.")
    return float(value)


def load_rho_mapping(path: Optional[Path]) -> Dict[str, float]:
    """Load per-library rho values from a CSV/TSV file with columns: library, rho."""
    if path is None:
        return {}

    if not path.exists():
        raise FileNotFoundError("Rho mapping file not found: {0}".format(path))

    with path.open("r", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        reader = csv.DictReader(handle, dialect=dialect)

        if reader.fieldnames is None:
            raise ValueError("Rho mapping file has no header: {0}".format(path))

        required = set(["library", "rho"])
        if not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                "Rho mapping file must contain 'library' and 'rho' columns."
            )

        mapping = {}
        for row in reader:
            lib = row.get("library")
            rho_value = row.get("rho")
            if lib is None or rho_value is None:
                continue
            mapping[str(lib)] = validate_rho(float(rho_value))

    if len(mapping) == 0:
        raise ValueError("No valid library entries were found in {0}".format(path))

    return mapping


def _to_csr(matrix) -> sp.csr_matrix:
    if sp.issparse(matrix):
        return matrix.tocsr()
    return sp.csr_matrix(matrix)


def _to_csc(matrix) -> sp.csc_matrix:
    if sp.issparse(matrix):
        return matrix.tocsc()
    return sp.csc_matrix(matrix)


def _safe_filename(label: str) -> str:
    safe = []
    for char in str(label):
        if char.isalnum() or char in ("-", "_", "."):
            safe.append(char)
        else:
            safe.append("_")
    out = "".join(safe).strip("._")
    if out == "":
        out = "library"
    return out


def _iter_libraries(values: Iterable) -> List[str]:
    seen = []
    seen_set = set()
    for value in values:
        if value is None:
            continue
        if isinstance(value, float) and np.isnan(value):
            continue
        sval = str(value)
        if sval in seen_set:
            continue
        seen.append(sval)
        seen_set.add(sval)
    return seen


def report_status(message: str) -> None:
    print("[STATUS] {0}".format(message), flush=True)


@contextmanager
def status_scope(message: str, done_message: Optional[str] = None):
    report_status(message)
    try:
        yield
    except Exception:
        report_status("FAILED: {0}".format(message))
        raise
    else:
        if done_message is not None:
            report_status(done_message)


def _subset_rows_csr(csr: sp.csr_matrix, row_idx: np.ndarray) -> sp.csr_matrix:
    """Fast row subset helper with explicit CSR output."""
    return csr[row_idx, :].tocsr()


def _compute_totals(csr: sp.csr_matrix) -> np.ndarray:
    return np.asarray(csr.sum(axis=1)).ravel().astype(np.float64, copy=False)


def _estimate_ambient_profile_from_filtered(
    counts_csr: sp.csr_matrix,
    *,
    empty_fraction: float = 0.10,
    min_empty_cells: int = 64,
    max_empty_cells: int = 4000,
    ambient_mode: str = "lowcount",
) -> Tuple[np.ndarray, Dict[str, float]]:
    """
    Estimate ambient profile from filtered cells.

    ambient_mode:
      - "lowcount": use lowest-UMI cells as pseudo-empty droplets
      - "global": use all cells
    """
    n_cells, n_genes = counts_csr.shape
    totals = _compute_totals(counts_csr)

    if ambient_mode not in ("lowcount", "global"):
        raise ValueError("Unsupported ambient_mode: {0}".format(ambient_mode))

    if n_cells == 0 or n_genes == 0:
        raise ValueError("Empty count matrix provided for ambient estimation.")

    if ambient_mode == "global":
        ambient_cells = np.arange(n_cells, dtype=np.int64)
    else:
        k = int(np.ceil(n_cells * float(empty_fraction)))
        if k < min_empty_cells:
            k = min(min_empty_cells, n_cells)
        if k > max_empty_cells:
            k = min(max_empty_cells, n_cells)
        if k < 1:
            k = 1
        ambient_cells = np.argpartition(totals, k - 1)[:k]

    soup = np.asarray(counts_csr[ambient_cells, :].sum(axis=0)).ravel().astype(np.float64)

    if soup.sum() <= 0:
        soup = np.asarray(counts_csr.sum(axis=0)).ravel().astype(np.float64)

    soup_sum = float(soup.sum())
    if soup_sum <= 0:
        soup = np.ones(n_genes, dtype=np.float64) / float(max(n_genes, 1))
    else:
        soup = soup / soup_sum

    meta = {
        "ambient_cells_used": int(ambient_cells.size),
        "ambient_fraction_requested": float(empty_fraction),
        "ambient_mode_lowcount": 1.0 if ambient_mode == "lowcount" else 0.0,
        "median_total_counts": float(np.median(totals)) if totals.size > 0 else 0.0,
        "min_total_counts": float(np.min(totals)) if totals.size > 0 else 0.0,
        "max_total_counts": float(np.max(totals)) if totals.size > 0 else 0.0,
    }
    return soup, meta


def _correct_counts_csc_subtraction(
    counts_csr: sp.csr_matrix,
    soup_profile: np.ndarray,
    rho: float,
    *,
    round_to_int: bool = False,
) -> sp.csr_matrix:
    """
    Fast sparse subtraction on nonzero entries only.

    Input:
        counts_csr: cells x genes
        soup_profile: length n_genes, sums to 1
        rho: scalar contamination fraction for this library

    Output:
        corrected counts in CSR format
    """
    rho = validate_rho(rho)

    if counts_csr.shape[1] != soup_profile.shape[0]:
        raise ValueError("soup_profile length does not match number of genes.")

    counts_csc = counts_csr.tocsc(copy=True)
    indptr = counts_csc.indptr
    indices = counts_csc.indices
    data = counts_csc.data.astype(np.float32, copy=False)

    cell_totals = _compute_totals(counts_csr).astype(np.float32, copy=False)
    soup_profile = np.asarray(soup_profile, dtype=np.float32)

    for gene_idx in range(counts_csc.shape[1]):
        start = indptr[gene_idx]
        end = indptr[gene_idx + 1]
        if start == end:
            continue

        cell_idx = indices[start:end]
        subtract_value = rho * soup_profile[gene_idx] * cell_totals[cell_idx]
        corrected = data[start:end] - subtract_value
        corrected[corrected < 0.0] = 0.0

        if round_to_int:
            corrected = np.rint(corrected)

        data[start:end] = corrected

    counts_csc.data = data
    counts_csc.eliminate_zeros()
    return counts_csc.tocsr()


def _summarize_matrix_change(
    before_csr: sp.csr_matrix,
    after_csr: sp.csr_matrix,
) -> Dict[str, float]:
    before_sum = float(before_csr.sum())
    after_sum = float(after_csr.sum())
    removed = before_sum - after_sum
    frac = 0.0
    if before_sum > 0:
        frac = removed / before_sum

    return {
        "counts_before": before_sum,
        "counts_after": after_sum,
        "counts_removed": removed,
        "fraction_removed": frac,
        "nnz_before": int(before_csr.nnz),
        "nnz_after": int(after_csr.nnz),
    }


def correct_library(
    subset: ad.AnnData,
    *,
    rho: float,
    method: str,
    round_to_int: bool,
    store_raw_layer: bool = True,
    replace_x: bool = True,
    raw_layer: str = RAW_LAYER_NAME,
    corrected_layer: str = CORRECTED_LAYER_NAME,
    library_label: Optional[str] = None,
    clusters=None,
    cluster_source: Optional[str] = None,
    ambient_mode: str = "lowcount",
    empty_fraction: float = 0.10,
    min_empty_cells: int = 64,
    max_empty_cells: int = 4000,
) -> ad.AnnData:
    """
    Apply fast filtered-only ambient correction to a single library AnnData in-place.

    method:
      - "subtraction"
      - "fastsub"
      - "fast_subtraction"
    """
    if subset.X is None or subset.n_obs == 0:
        raise ValueError("Empty library subset provided.")

    supported = set(["subtraction", "fastsub", "fast_subtraction"])
    if method not in supported:
        raise ValueError(
            "Unsupported method '{0}'. Use one of: {1}".format(
                method, ", ".join(sorted(supported))
            )
        )

    counts_csr = _to_csr(subset.X).astype(np.float32)

    if store_raw_layer:
        subset.layers[raw_layer] = counts_csr.copy()

    soup_profile, ambient_meta = _estimate_ambient_profile_from_filtered(
        counts_csr,
        empty_fraction=empty_fraction,
        min_empty_cells=min_empty_cells,
        max_empty_cells=max_empty_cells,
        ambient_mode=ambient_mode,
    )

    corrected_csr = _correct_counts_csc_subtraction(
        counts_csr,
        soup_profile,
        rho,
        round_to_int=round_to_int,
    )

    subset.layers[corrected_layer] = corrected_csr
    if replace_x:
        subset.X = corrected_csr.copy()

    subset.obs["soupx_rho"] = validate_rho(rho)
    if library_label is not None:
        subset.obs["soupx_library"] = str(library_label)

    if clusters is not None:
        subset.obs["soupx_cluster"] = np.asarray(clusters).astype(str)

    change = _summarize_matrix_change(counts_csr, corrected_csr)

    metadata = {
        "library": "" if library_label is None else str(library_label),
        "rho": float(rho),
        "method": str(method),
        "round_to_int": bool(round_to_int),
        "cluster_column": "" if cluster_source is None else str(cluster_source),
        "ambient_mode": str(ambient_mode),
        "empty_fraction": float(empty_fraction),
        "min_empty_cells": int(min_empty_cells),
        "max_empty_cells": int(max_empty_cells),
        "cells": int(subset.n_obs),
        "genes": int(subset.n_vars),
        "ambient_cells_used": int(ambient_meta["ambient_cells_used"]),
        "counts_before": float(change["counts_before"]),
        "counts_after": float(change["counts_after"]),
        "counts_removed": float(change["counts_removed"]),
        "fraction_removed": float(change["fraction_removed"]),
        "nnz_before": int(change["nnz_before"]),
        "nnz_after": int(change["nnz_after"]),
    }
    subset.uns["soupx_correction"] = metadata
    return subset


def _process_library_matrix_only(
    counts_csr: sp.csr_matrix,
    *,
    rho: float,
    method: str,
    round_to_int: bool,
    ambient_mode: str,
    empty_fraction: float,
    min_empty_cells: int,
    max_empty_cells: int,
) -> Tuple[sp.csr_matrix, Dict[str, float]]:
    supported = set(["subtraction", "fastsub", "fast_subtraction"])
    if method not in supported:
        raise ValueError(
            "Unsupported method '{0}'. Use one of: {1}".format(
                method, ", ".join(sorted(supported))
            )
        )

    soup_profile, ambient_meta = _estimate_ambient_profile_from_filtered(
        counts_csr,
        empty_fraction=empty_fraction,
        min_empty_cells=min_empty_cells,
        max_empty_cells=max_empty_cells,
        ambient_mode=ambient_mode,
    )

    corrected_csr = _correct_counts_csc_subtraction(
        counts_csr,
        soup_profile,
        rho,
        round_to_int=round_to_int,
    )

    change = _summarize_matrix_change(counts_csr, corrected_csr)
    meta = {}
    meta.update(ambient_meta)
    meta.update(change)
    return corrected_csr, meta


def process_anndata(
    adata: ad.AnnData,
    *,
    rho: float = DEFAULT_RHO,
    library_col: str = "Library",
    outdir: Path = Path("./soupx_corrected"),
    rho_mapping: Optional[Dict[str, float]] = None,
    store_raw_layer: bool = True,
    replace_x: bool = True,
    write_individual: bool = True,
    merged_filename: str = "soupx_corrected_merged.h5ad",
    cluster_col: Optional[str] = None,
    auto_cluster: bool = True,
    leiden_resolution: float = 0.5,
    leiden_neighbors: int = 15,
    leiden_npcs: int = 30,
    leiden_hvg: int = 3000,
    leiden_random_state: int = 0,
    method: str = "subtraction",
    round_to_int: bool = False,
    input_label: str = "anndata",
    ambient_mode: str = "lowcount",
    empty_fraction: float = 0.10,
    min_empty_cells: int = 64,
    max_empty_cells: int = 4000,
) -> ad.AnnData:
    """
    Process libraries within an AnnData object and write corrected outputs.

    Notes
    -----
    - auto_cluster / Leiden arguments are retained for API compatibility only.
      They are not used in the fast filtered-only correction path.
    - cluster_col is retained and propagated into metadata if present.
    """
    rho_mapping = {} if rho_mapping is None else dict(rho_mapping)

    if library_col not in adata.obs:
        raise ValueError("obs column '{0}' not found in AnnData.".format(library_col))

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    working = adata.copy()
    x_csr = _to_csr(working.X).astype(np.float32)

    if store_raw_layer:
        working.layers[RAW_LAYER_NAME] = x_csr.copy()

    libraries = _iter_libraries(working.obs[library_col].tolist())
    if len(libraries) == 0:
        raise ValueError("No library values found in column '{0}'.".format(library_col))

    report_status(
        "Loaded {0} with {1} cells, {2} genes, {3} libraries.".format(
            input_label,
            working.n_obs,
            working.n_vars,
            len(libraries),
        )
    )

    if auto_cluster:
        warnings.warn(
            "auto_cluster=True was requested, but the fast filtered-only path "
            "does not run clustering. Provide --cluster-col if you want cluster "
            "labels carried into metadata."
        )

    corrected_blocks = []
    row_order_blocks = []
    summary_rows = []

    all_library_values = working.obs[library_col].astype(str).to_numpy()

    for idx, lib in enumerate(libraries, start=1):
        lib_label = str(lib)
        row_idx = np.where(all_library_values == lib_label)[0]
        if row_idx.size == 0:
            report_status("Skipping empty library: {0}".format(lib_label))
            continue

        lib_rho = validate_rho(rho_mapping.get(lib_label, rho))

        report_status(
            "[{0}/{1}] Starting library '{2}' (cells={3}, rho={4:.4f})".format(
                idx,
                len(libraries),
                lib_label,
                row_idx.size,
                lib_rho,
            )
        )

        t0 = time.perf_counter()
        subset_csr = _subset_rows_csr(x_csr, row_idx)

        cluster_source = None
        clusters = None
        if cluster_col is not None:
            if cluster_col not in working.obs:
                raise ValueError(
                    "Requested cluster column '{0}' not found.".format(cluster_col)
                )
            clusters = working.obs.iloc[row_idx][cluster_col].astype(str).to_numpy()
            cluster_source = cluster_col

        corrected_subset_csr, lib_meta = _process_library_matrix_only(
            subset_csr,
            rho=lib_rho,
            method=method,
            round_to_int=round_to_int,
            ambient_mode=ambient_mode,
            empty_fraction=empty_fraction,
            min_empty_cells=min_empty_cells,
            max_empty_cells=max_empty_cells,
        )

        elapsed = time.perf_counter() - t0

        corrected_blocks.append(corrected_subset_csr)
        row_order_blocks.append(row_idx)

        summary_row = {
            "library": lib_label,
            "cells": int(row_idx.size),
            "genes": int(subset_csr.shape[1]),
            "rho": float(lib_rho),
            "method": str(method),
            "ambient_mode": str(ambient_mode),
            "ambient_cells_used": int(lib_meta["ambient_cells_used"]),
            "counts_before": float(lib_meta["counts_before"]),
            "counts_after": float(lib_meta["counts_after"]),
            "counts_removed": float(lib_meta["counts_removed"]),
            "fraction_removed": float(lib_meta["fraction_removed"]),
            "nnz_before": int(lib_meta["nnz_before"]),
            "nnz_after": int(lib_meta["nnz_after"]),
            "cluster_column": "" if cluster_source is None else str(cluster_source),
            "total_time_sec": float(elapsed),
        }
        summary_rows.append(summary_row)

        report_status(
            "[{0}/{1}] Completed library '{2}' in {3:.2f}s".format(
                idx,
                len(libraries),
                lib_label,
                elapsed,
            )
        )

        if write_individual:
            safe_name = _safe_filename(lib_label)
            out_file = outdir / "{0}_soupx_corrected.h5ad".format(safe_name)

            report_status(
                "[{0}/{1}] Writing per-library file to {2}".format(
                    idx,
                    len(libraries),
                    out_file,
                )
            )

            subset_adata = working[row_idx, :].copy()
            if store_raw_layer:
                subset_adata.layers[RAW_LAYER_NAME] = subset_csr.copy()
            subset_adata.layers[CORRECTED_LAYER_NAME] = corrected_subset_csr
            if replace_x:
                subset_adata.X = corrected_subset_csr.copy()

            subset_adata.obs["soupx_rho"] = lib_rho
            subset_adata.obs["soupx_library"] = lib_label
            if clusters is not None:
                subset_adata.obs["soupx_cluster"] = clusters

            subset_adata.uns["soupx_correction"] = {
                "library": lib_label,
                "rho": float(lib_rho),
                "method": str(method),
                "ambient_mode": str(ambient_mode),
                "ambient_cells_used": int(lib_meta["ambient_cells_used"]),
                "cluster_column": "" if cluster_source is None else str(cluster_source),
                "counts_before": float(lib_meta["counts_before"]),
                "counts_after": float(lib_meta["counts_after"]),
                "counts_removed": float(lib_meta["counts_removed"]),
                "fraction_removed": float(lib_meta["fraction_removed"]),
                "nnz_before": int(lib_meta["nnz_before"]),
                "nnz_after": int(lib_meta["nnz_after"]),
                "total_time_sec": float(elapsed),
            }

            subset_adata.write_h5ad(out_file, compression="gzip")

    if len(corrected_blocks) == 0:
        raise RuntimeError("No valid libraries processed.")

    stacked = sp.vstack(corrected_blocks, format="csr")
    stacked_row_order = np.concatenate(row_order_blocks)

    inverse = np.empty(stacked_row_order.shape[0], dtype=np.int64)
    inverse[stacked_row_order] = np.arange(stacked_row_order.shape[0], dtype=np.int64)

    corrected_all = stacked[inverse, :].tocsr()

    working.layers[CORRECTED_LAYER_NAME] = corrected_all
    if replace_x:
        working.X = corrected_all.copy()

    summary_table = pd.DataFrame(summary_rows)

    working.uns.setdefault("soupx_correction", {})
    working.uns["soupx_correction"].update(
        {
            "default_rho": float(rho),
            "library_column": str(library_col),
            "cluster_column": "" if cluster_col is None else str(cluster_col),
            "auto_cluster_requested": bool(auto_cluster),
            "auto_cluster_used": False,
            "leiden_parameters_requested": {
                "resolution": float(leiden_resolution),
                "n_neighbors": int(leiden_neighbors),
                "n_pcs": int(leiden_npcs),
                "n_top_genes": int(leiden_hvg),
                "random_state": int(leiden_random_state),
            },
            "method": str(method),
            "round_to_int": bool(round_to_int),
            "ambient_mode": str(ambient_mode),
            "empty_fraction": float(empty_fraction),
            "min_empty_cells": int(min_empty_cells),
            "max_empty_cells": int(max_empty_cells),
            "libraries": summary_table,
        }
    )

    merged_path = outdir / merged_filename
    report_status("Writing merged corrected dataset to {0}".format(merged_path))
    working.write_h5ad(merged_path, compression="gzip")

    if not summary_table.empty:
        summary_path = outdir / "soupx_summary.tsv"
        summary_table.to_csv(summary_path, sep="\t", index=False)
        report_status("Summary table written to {0}".format(summary_path))

    return working


def process_h5ad(
    h5ad_path: Path,
    *,
    rho: float = DEFAULT_RHO,
    library_col: str = "Library",
    outdir: Path = Path("./soupx_corrected"),
    rho_mapping: Optional[Dict[str, float]] = None,
    store_raw_layer: bool = True,
    replace_x: bool = True,
    write_individual: bool = True,
    merged_filename: str = "soupx_corrected_merged.h5ad",
    cluster_col: Optional[str] = None,
    auto_cluster: bool = True,
    leiden_resolution: float = 0.5,
    leiden_neighbors: int = 15,
    leiden_npcs: int = 30,
    leiden_hvg: int = 3000,
    leiden_random_state: int = 0,
    method: str = "subtraction",
    round_to_int: bool = False,
    ambient_mode: str = "lowcount",
    empty_fraction: float = 0.10,
    min_empty_cells: int = 64,
    max_empty_cells: int = 4000,
) -> ad.AnnData:
    if not h5ad_path.exists():
        raise FileNotFoundError("Input h5ad not found: {0}".format(h5ad_path))

    rho_mapping = {} if rho_mapping is None else dict(rho_mapping)
    adata = ad.read_h5ad(h5ad_path)

    return process_anndata(
        adata,
        rho=rho,
        library_col=library_col,
        outdir=outdir,
        rho_mapping=rho_mapping,
        store_raw_layer=store_raw_layer,
        replace_x=replace_x,
        write_individual=write_individual,
        merged_filename=merged_filename,
        cluster_col=cluster_col,
        auto_cluster=auto_cluster,
        leiden_resolution=leiden_resolution,
        leiden_neighbors=leiden_neighbors,
        leiden_npcs=leiden_npcs,
        leiden_hvg=leiden_hvg,
        leiden_random_state=leiden_random_state,
        method=method,
        round_to_int=round_to_int,
        input_label=h5ad_path.name,
        ambient_mode=ambient_mode,
        empty_fraction=empty_fraction,
        min_empty_cells=min_empty_cells,
        max_empty_cells=max_empty_cells,
    )


def run_soupx_correction(
    h5ad_path: str,
    *,
    rho: float = DEFAULT_RHO,
    library_col: str = "library",
    outdir: str = "./soupx_corrected",
    rho_map_path: Optional[str] = None,
    store_raw_layer: bool = True,
    replace_x: bool = True,
    write_individual: bool = True,
    merged_filename: str = "soupx_corrected_merged.h5ad",
    cluster_col: Optional[str] = None,
    auto_cluster: bool = True,
    leiden_resolution: float = 0.5,
    leiden_neighbors: int = 15,
    leiden_npcs: int = 30,
    leiden_hvg: int = 3000,
    leiden_random_state: int = 0,
    method: str = "subtraction",
    round_to_int: bool = False,
    return_adata: bool = False,
    ambient_mode: str = "lowcount",
    empty_fraction: float = 0.10,
    min_empty_cells: int = 64,
    max_empty_cells: int = 4000,
) -> Optional[ad.AnnData]:
    rho_mapping = load_rho_mapping(Path(rho_map_path)) if rho_map_path else {}

    corrected = process_h5ad(
        Path(h5ad_path),
        rho=rho,
        library_col=library_col,
        outdir=Path(outdir),
        rho_mapping=rho_mapping,
        store_raw_layer=store_raw_layer,
        replace_x=replace_x,
        write_individual=write_individual,
        merged_filename=merged_filename,
        cluster_col=cluster_col,
        auto_cluster=auto_cluster,
        leiden_resolution=leiden_resolution,
        leiden_neighbors=leiden_neighbors,
        leiden_npcs=leiden_npcs,
        leiden_hvg=leiden_hvg,
        leiden_random_state=leiden_random_state,
        method=method,
        round_to_int=round_to_int,
        ambient_mode=ambient_mode,
        empty_fraction=empty_fraction,
        min_empty_cells=min_empty_cells,
        max_empty_cells=max_empty_cells,
    )

    if return_adata:
        return corrected
    return None


def build_cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fast filtered-only ambient RNA correction using direct sparse subtraction."
    )

    parser.add_argument("--h5ad", required=True, help="Input filtered h5ad file.")
    parser.add_argument(
        "--rho",
        type=float,
        default=DEFAULT_RHO,
        help="Default contamination fraction applied to libraries absent from --rho-map.",
    )
    parser.add_argument(
        "--rho-map",
        default=None,
        help="Optional CSV/TSV file with columns: library,rho",
    )
    parser.add_argument(
        "--library-col",
        default="Library",
        help="obs column containing library/sample identifiers.",
    )
    parser.add_argument(
        "--cluster-col",
        default=None,
        help="Optional obs column containing cluster labels. Stored as metadata only.",
    )
    parser.add_argument(
        "--no-auto-cluster",
        dest="auto_cluster",
        action="store_false",
        help="Ignored by the fast path. Retained for CLI compatibility.",
    )
    parser.add_argument(
        "--leiden-resolution",
        type=float,
        default=0.5,
        help="Ignored by the fast path. Retained for CLI compatibility.",
    )
    parser.add_argument(
        "--leiden-neighbors",
        type=int,
        default=15,
        help="Ignored by the fast path. Retained for CLI compatibility.",
    )
    parser.add_argument(
        "--leiden-npcs",
        type=int,
        default=30,
        help="Ignored by the fast path. Retained for CLI compatibility.",
    )
    parser.add_argument(
        "--leiden-hvg",
        type=int,
        default=3000,
        help="Ignored by the fast path. Retained for CLI compatibility.",
    )
    parser.add_argument(
        "--leiden-random-state",
        type=int,
        default=0,
        help="Ignored by the fast path. Retained for CLI compatibility.",
    )
    parser.add_argument(
        "--method",
        default="subtraction",
        choices=["subtraction", "fastsub", "fast_subtraction"],
        help="Fast correction method. 'subtraction' is the intended default.",
    )
    parser.add_argument(
        "--round-to-int",
        action="store_true",
        help="Round corrected values to nearest integers.",
    )
    parser.add_argument(
        "--ambient-mode",
        default="lowcount",
        choices=["lowcount", "global"],
        help="How to estimate the ambient profile from filtered counts.",
    )
    parser.add_argument(
        "--empty-fraction",
        type=float,
        default=0.10,
        help="Fraction of lowest-UMI cells used as pseudo-empty droplets in lowcount mode.",
    )
    parser.add_argument(
        "--min-empty-cells",
        type=int,
        default=64,
        help="Minimum number of pseudo-empty cells used for ambient estimation.",
    )
    parser.add_argument(
        "--max-empty-cells",
        type=int,
        default=4000,
        help="Maximum number of pseudo-empty cells used for ambient estimation.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for corrected h5ad files.",
    )
    parser.add_argument(
        "--merged-name",
        default="soupx_corrected_merged.h5ad",
        help="Filename for the merged corrected h5ad written inside outdir.",
    )
    parser.add_argument(
        "--no-store-raw-layer",
        dest="store_raw_layer",
        action="store_false",
        help="Do not persist original counts in the 'soupx_raw' layer.",
    )
    parser.add_argument(
        "--keep-original-x",
        dest="replace_x",
        action="store_false",
        help="Keep original counts in .X and write corrected counts only to the corrected layer.",
    )
    parser.add_argument(
        "--skip-individual",
        dest="write_individual",
        action="store_false",
        help="Do not emit per-library corrected h5ad files.",
    )

    parser.set_defaults(
        store_raw_layer=True,
        replace_x=True,
        write_individual=True,
        auto_cluster=True,
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_cli_parser()
    args = parser.parse_args(argv)

    run_soupx_correction(
        h5ad_path=args.h5ad,
        rho=args.rho,
        library_col=args.library_col,
        outdir=args.outdir,
        rho_map_path=args.rho_map,
        store_raw_layer=args.store_raw_layer,
        replace_x=args.replace_x,
        write_individual=args.write_individual,
        merged_filename=args.merged_name,
        cluster_col=args.cluster_col,
        auto_cluster=args.auto_cluster,
        leiden_resolution=args.leiden_resolution,
        leiden_neighbors=args.leiden_neighbors,
        leiden_npcs=args.leiden_npcs,
        leiden_hvg=args.leiden_hvg,
        leiden_random_state=args.leiden_random_state,
        method=args.method,
        round_to_int=args.round_to_int,
        return_adata=False,
        ambient_mode=args.ambient_mode,
        empty_fraction=args.empty_fraction,
        min_empty_cells=args.min_empty_cells,
        max_empty_cells=args.max_empty_cells,
    )


if __name__ == "__main__":
    main()