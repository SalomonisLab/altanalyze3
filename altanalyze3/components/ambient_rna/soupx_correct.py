#!/usr/bin/env python3
"""Ambient RNA correction utility using soupx-python.

This script consumes a filtered ``.h5ad`` file that already has an estimated
contamination fraction (rho) and applies SoupX adjustments per library stored in
``adata.obs[library_col]``. By design it does **not** estimate rho internally –
you provide either a global value or an optional mapping with pre-computed
library-specific fractions.

Highlights
~~~~~~~~~~
* Accept either a single rho value or per-library mapping.
* Automatically computes per-library Leiden clusters (Scanpy) to speed SoupX;
  override with `--cluster-col` or disable via CLI if needed.
* Accept optional precomputed cluster labels via `--cluster-col` for reproducible runs.
* Processes each library independently to avoid cross-library bleed.
* Stores the raw counts in a layer (``soupx_raw``) and overwrites ``.X`` with the
  corrected counts unless disabled.
* Persists corrected counts in a ``soupx_corrected`` layer for downstream audit.
* Emits progress updates via ``tqdr`` when available (falls back to console
  prints otherwise).
* When top-of-droplet (raw) matrices are unavailable, derives a fallback ambient
  profile directly from each library's filtered counts.

Example (standalone)::

    python3 soupx_correct.py \
        --h5ad input.h5ad \
        --rho 0.12 \
        --library-col library_id \
        --outdir ./soupx_corrected

Programmatic use::

    from soupx_correct import run_soupx_correction
    corrected = run_soupx_correction("input.h5ad", rho=0.10, return_adata=True)

Dependencies
~~~~~~~~~~~~
``soupx``_, ``anndata``, ``numpy``, ``scipy``

.. _soupx: https://github.com/NiRuff/soupx-python
"""

from __future__ import annotations

import argparse
import csv
import os
import warnings
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, Iterable, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse import SparseEfficiencyWarning
import scanpy as sc

try:  # soupx is mandatory for the core functionality
    import soupx
except ImportError as exc:  # pragma: no cover - explicit error path
    raise ImportError(
        "soupx_correct requires the soupx-python package. Install it with 'pip install soupx'."
    ) from exc

try:  # tqdr is optional but requested for status updates
    import tqdr  # type: ignore
except ImportError:  # pragma: no cover - optional dependency
    tqdr = None


DEFAULT_RHO = 0.15
RAW_LAYER_NAME = "soupx_raw"
CORRECTED_LAYER_NAME = "soupx_corrected"


def _select_tqdr_reporter() -> Optional[callable]:  # pragma: no cover - import-time helper
    if tqdr is None:  # type: ignore[name-defined]
        return None
    for candidate in ("status", "update", "log", "info", "write"):
        reporter = getattr(tqdr, candidate, None)
        if callable(reporter):
            return reporter
    return None


_TQDR_REPORTER = _select_tqdr_reporter()


warnings.filterwarnings(
    "ignore",
    message="`flavor='seurat_v3'` expects raw count data, but non-integers were found.",
    category=UserWarning,
)
warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)


def _try_reporter(message: str, reporter: callable, **kwargs) -> bool:
    """Attempt different call signatures for tqdr reporters."""

    try:
        reporter(message=message, **kwargs)
        return True
    except TypeError:
        pass
    try:
        reporter(message)
        return True
    except TypeError:
        pass
    try:
        reporter(message, **kwargs)
        return True
    except TypeError:
        return False


def report_status(message: str, **kwargs) -> None:
    """Send status updates via tqdr when available, otherwise print."""

    if _TQDR_REPORTER is not None:
        if _try_reporter(message, _TQDR_REPORTER, **kwargs):
            return
    print(f"[STATUS] {message}")


@contextmanager
def status_scope(message: str, done_message: Optional[str] = None, *, status_kwargs: Optional[dict] = None):
    """Emit a start message and optionally a completion message."""

    status_kwargs = status_kwargs or {}
    report_status(message, **status_kwargs)
    try:
        yield
    except Exception:
        report_status(f"FAILED: {message}", **status_kwargs)
        raise
    else:
        if done_message:
            report_status(done_message, **status_kwargs)


def validate_rho(value: float) -> float:
    """Ensure rho stays within 0-1 (inclusive)."""

    if not np.isfinite(value):
        raise ValueError("Contamination fraction (rho) must be a finite number.")
    if value < 0 or value > 1:
        raise ValueError("Contamination fraction (rho) must be between 0 and 1 inclusive.")
    return float(value)


def load_rho_mapping(path: Optional[Path]) -> Dict[str, float]:
    """Load per-library rho values from a CSV/TSV file."""

    if path is None:
        return {}
    if not path.exists():
        raise FileNotFoundError(f"Rho mapping file not found: {path}")

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or {"library", "rho"} - set(reader.fieldnames):
            raise ValueError("Rho mapping file must contain 'library' and 'rho' columns.")
        mapping: Dict[str, float] = {}
        for row in reader:
            library = row.get("library")
            rho_value = row.get("rho")
            if library is None or rho_value is None:
                continue
            try:
                mapping[str(library)] = validate_rho(float(rho_value))
            except ValueError as exc:
                raise ValueError(f"Invalid rho value for library '{library}': {rho_value}") from exc
    if not mapping:
        raise ValueError(f"No valid library entries were found in {path}.")
    return mapping


def _to_csr(matrix) -> sp.csr_matrix:
    """Return matrix as CSR without densifying unnecessarily."""

    if sp.issparse(matrix):
        return matrix.tocsr()
    return sp.csr_matrix(matrix)


def _safe_filename(label: str) -> str:
    """Sanitize a label for filesystem usage."""

    safe_chars = []
    for char in label:
        if char.isalnum() or char in ("-", "_", "."):
            safe_chars.append(char)
        else:
            safe_chars.append("_")
    sanitized = "".join(safe_chars).strip("._")
    return sanitized or "library"


def _run_leiden_clustering(
    subset: ad.AnnData,
    *,
    resolution: float,
    n_neighbors: int,
    n_pcs: int,
    n_top_genes: int,
    random_state: int,
) -> tuple[np.ndarray, str]:
    """Compute Leiden clusters on a copy of the subset and return labels."""

    if subset.n_obs < 3 or subset.n_vars < 3:
        clusters = np.array([f"cell_{i}" for i in range(subset.n_obs)], dtype=str)
        subset.obs["soupx_leiden"] = clusters
        return clusters, "auto_identity"

    work = subset.copy()

    if sp.issparse(work.X):
        work.X = work.X.astype(np.float32)
    else:
        work.X = np.asarray(work.X, dtype=np.float32)

    prev_verbosity = sc.settings.verbosity
    sc.settings.verbosity = 0
    try:
        sc.pp.normalize_total(work, target_sum=1e4, inplace=True)
        sc.pp.log1p(work)

        hvg = min(max(200, n_top_genes), work.n_vars)
        if work.n_vars > hvg:
            sc.pp.highly_variable_genes(work, n_top_genes=hvg, flavor="seurat_v3", subset=True)

        max_pcs = min(n_pcs, max(1, work.n_obs - 1), max(1, work.n_vars - 1))
        if max_pcs < 1:
            max_pcs = 1
        sc.pp.pca(work, n_comps=max_pcs, random_state=random_state)

        effective_neighbors = max(2, min(n_neighbors, work.n_obs - 1))
        sc.pp.neighbors(work, n_neighbors=effective_neighbors, n_pcs=max_pcs, random_state=random_state)
        sc.tl.leiden(
            work,
            resolution=resolution,
            random_state=random_state,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
    finally:
        sc.settings.verbosity = prev_verbosity

    clusters = work.obs["leiden"].astype(str).to_numpy()
    subset.obs["soupx_leiden"] = clusters
    return clusters, f"auto_leiden(res={resolution},k={effective_neighbors},pcs={max_pcs})"


def _configure_soup_profile(sc, counts_csr: sp.csr_matrix) -> str:
    """Ensure SoupChannel has an ambient profile even without raw droplets."""

    # Attempt explicit helper methods provided by soupx-python.
    for method_name in ("setSoupProfileFromTOC", "setSoupProfileFromSoup", "setSoupProfileFromTou"):
        method = getattr(sc, method_name, None)
        if callable(method):
            try:
                method()
                return method_name + "()"
            except TypeError:
                for candidate in (counts_csr, counts_csr.T):
                    try:
                        method(candidate)
                        return f"{method_name}(matrix)"
                    except Exception:
                        continue

    estimate = getattr(sc, "estimateSoupProfile", None)
    setter = getattr(sc, "setSoupProfile", None)
    if callable(estimate) and callable(setter):
        try:
            profile = estimate()
            setter(profile)
            return "estimateSoupProfile()"
        except Exception:
            pass

    # Final fallback: normalize summed gene expression across cells.
    gene_sums = np.asarray(counts_csr.sum(axis=0)).ravel()
    if gene_sums.sum() > 0:
        gene_sums = gene_sums / gene_sums.sum()
    else:
        gene_sums[:] = 1.0 / max(gene_sums.size, 1)

    for attr in ("soup_profile", "soupProfile", "ambientProfile", "_soupProfile"):
        if hasattr(sc, attr):
            setattr(sc, attr, gene_sums)
            return f"setattr:{attr}"

    # Store as a generic attribute; SoupX will read it back if it expects soup_profile.
    sc.soup_profile = gene_sums  # type: ignore[attr-defined]
    return "setattr:soup_profile"


def _adjust_counts(sc, clusters, *, method: str, round_to_int: bool, verbose: int = 0):
    """Call soupx.adjustCounts while tolerating differing signatures."""

    kwargs = {
        "method": method,
        "roundToInt": round_to_int,
        "verbose": verbose,
    }

    try:
        if clusters is None:
            return soupx.adjustCounts(sc, **kwargs)
        return soupx.adjustCounts(sc, clusters=clusters, **kwargs)
    except TypeError:
        if clusters is None:
            return soupx.adjustCounts(sc, **kwargs)
        try:
            return soupx.adjustCounts(sc, clusters, **kwargs)
        except TypeError:
            # Fall back to positional without additional kwargs
            if clusters is None:
                return soupx.adjustCounts(sc, verbose=verbose)
            return soupx.adjustCounts(sc, clusters, verbose=verbose)


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
    clusters = None,
    cluster_source: Optional[str] = None,
) -> ad.AnnData:
    """Apply SoupX correction to a single library AnnData in-place."""

    if subset.X is None or subset.n_obs == 0:
        raise ValueError("Empty library subset provided.")

    counts_csr = _to_csr(subset.X)
    if store_raw_layer:
        subset.layers[raw_layer] = counts_csr.copy()

    try:
        sc = soupx.SoupChannel(tod=counts_csr.T, toc=counts_csr.T)
        soup_channel_mode = "tod=toc"
    except TypeError:
        try:
            sc = soupx.SoupChannel(counts_csr.T)
            soup_channel_mode = "single-arg"
        except TypeError:
            sc = soupx.SoupChannel(toc=counts_csr.T)
            soup_channel_mode = "named-toc"
    rho_value = validate_rho(rho)

    set_contam = getattr(sc, "set_contamination_fraction", None)
    if callable(set_contam):
        set_contam(rho_value)
    else:
        sc.metaData["rho"] = rho_value

    soup_profile_source = _configure_soup_profile(sc, counts_csr)
    corrected = _adjust_counts(sc, clusters, method=method, round_to_int=round_to_int)
    corrected_csr = _to_csr(corrected).T

    subset.layers[corrected_layer] = corrected_csr
    if replace_x:
        subset.X = corrected_csr.copy()

    subset.obs["soupx_rho"] = rho_value
    if library_label is not None:
        subset.obs["soupx_library"] = library_label
    if clusters is not None:
        subset.obs["soupx_cluster"] = clusters

    cluster_count = int(np.unique(clusters).size) if clusters is not None else None

    metadata = {
        "library": library_label or "unknown",
        "rho": float(rho_value),
        "soupx_version": getattr(soupx, "__version__", "unknown"),
        "cells": int(subset.n_obs),
        "genes": int(subset.n_vars),
        "soup_profile_source": soup_profile_source,
        "soup_channel_mode": soup_channel_mode,
        "method": method,
        "round_to_int": bool(round_to_int),
        "n_clusters": int(cluster_count) if cluster_count is not None else 0,
        "cluster_column": cluster_source or "",
    }
    subset.uns["soupx_correction"] = metadata
    return subset


def _iter_libraries(values: Iterable) -> Iterable:
    """Yield sorted unique libraries, skipping NA values."""

    seen = []
    for value in values:
        if value is None:
            continue
        if isinstance(value, float) and np.isnan(value):
            continue
        if value in seen:
            continue
        seen.append(value)
    return seen


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
) -> ad.AnnData:
    """Process libraries within an AnnData object and write corrected outputs."""

    rho_mapping = rho_mapping or {}
    if library_col not in adata.obs:
        raise ValueError(f"obs column '{library_col}' not found in AnnData.")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    working = adata.copy()

    libraries = list(_iter_libraries(working.obs[library_col].tolist()))
    if not libraries:
        raise ValueError(f"No library values found in column '{library_col}'.")

    total_libraries = len(libraries)
    report_status(
        "Loaded {file} with {cells} cells, {genes} genes, {libs} libraries. "
        "Method={method} round_to_int={round_to_int} clustering={cluster_info}".format(
            file=input_label,
            cells=working.n_obs,
            genes=working.n_vars,
            libs=total_libraries,
            method=method,
            round_to_int=round_to_int,
            cluster_info=(
                f"obs:{cluster_col}" if cluster_col
                else ("auto_leiden" if auto_cluster else "cell_level")
            ),
        )
    )

    corrected_libs = []
    summary = []

    for idx, lib in enumerate(libraries, start=1):
        lib_mask = working.obs[library_col] == lib
        subset = working[lib_mask].copy()
        if subset.n_obs == 0:
            report_status(f"Skipping empty library: {lib}")
            continue

        lib_label = str(lib)
        lib_rho = validate_rho(rho_mapping.get(lib_label, rho))
        progress_kwargs = {"progress": {"current": idx, "total": total_libraries}}

        report_status(
            f"[{idx}/{total_libraries}] Starting library '{lib_label}' (cells={subset.n_obs}, rho={lib_rho:.4f})",
            **progress_kwargs,
        )
        library_start = time.perf_counter()

        clusters = None
        cluster_source = None
        cluster_time = 0.0

        if cluster_col:
            if cluster_col not in subset.obs:
                raise ValueError(
                    f"Requested cluster column '{cluster_col}' not found for library '{lib_label}'."
                )
            clusters = subset.obs[cluster_col].astype(str).to_numpy()
            cluster_source = cluster_col
            report_status(
                f"[{idx}/{total_libraries}] Using clusters from column '{cluster_col}'",
                **progress_kwargs,
            )
        elif auto_cluster:
            report_status(
                f"[{idx}/{total_libraries}] Leiden clustering '{lib_label}' (resolution={leiden_resolution})",
                **progress_kwargs,
            )
            cluster_start = time.perf_counter()
            clusters, cluster_source = _run_leiden_clustering(
                subset,
                resolution=leiden_resolution,
                n_neighbors=leiden_neighbors,
                n_pcs=leiden_npcs,
                n_top_genes=leiden_hvg,
                random_state=leiden_random_state,
            )
            cluster_time = time.perf_counter() - cluster_start
            report_status(
                f"[{idx}/{total_libraries}] Clustering completed in {cluster_time:.1f}s",
                **progress_kwargs,
            )

        correction_start = time.perf_counter()
        corrected_subset = correct_library(
            subset,
            rho=lib_rho,
            method=method,
            round_to_int=round_to_int,
            store_raw_layer=store_raw_layer,
            replace_x=replace_x,
            library_label=lib_label,
            clusters=clusters,
            cluster_source=cluster_source,
        )
        correction_time = time.perf_counter() - correction_start
        report_status(
            f"[{idx}/{total_libraries}] SoupX correction completed in {correction_time:.1f}s",
            **progress_kwargs,
        )

        library_metadata = corrected_subset.uns.get("soupx_correction", {}).copy()
        library_metadata.update(
            {
                "cluster_time_sec": round(cluster_time, 3),
                "correction_time_sec": round(correction_time, 3),
            }
        )

        write_time = 0.0
        if write_individual:
            safe_name = _safe_filename(lib_label)
            out_file = outdir / f"{safe_name}_soupx_corrected.h5ad"
            report_status(
                f"[{idx}/{total_libraries}] Writing per-library file to {out_file}",
                **progress_kwargs,
            )
            write_start = time.perf_counter()
            corrected_subset.write_h5ad(out_file, compression="gzip")
            write_time = time.perf_counter() - write_start
            report_status(
                f"[{idx}/{total_libraries}] Per-library file written in {write_time:.1f}s",
                **progress_kwargs,
            )

        library_elapsed = time.perf_counter() - library_start
        library_metadata.update(
            {
                "write_time_sec": round(write_time, 3),
                "total_time_sec": round(library_elapsed, 3),
            }
        )
        corrected_subset.uns["soupx_correction"] = library_metadata

        summary.append(
            {
                "library": lib_label,
                "rho": float(library_metadata.get("rho", lib_rho)),
                "cells": int(library_metadata.get("cells", corrected_subset.n_obs)),
                "genes": int(library_metadata.get("genes", corrected_subset.n_vars)),
                "cluster_column": library_metadata.get("cluster_column", cluster_source or ""),
                "n_clusters": int(library_metadata.get("n_clusters", 0)),
                "method": library_metadata.get("method", method),
                "round_to_int": bool(library_metadata.get("round_to_int", round_to_int)),
                "cluster_time_sec": library_metadata.get("cluster_time_sec", 0.0),
                "correction_time_sec": library_metadata.get("correction_time_sec", 0.0),
                "write_time_sec": library_metadata.get("write_time_sec", write_time),
                "total_time_sec": library_metadata.get("total_time_sec", library_elapsed),
            }
        )

        report_status(
            f"[{idx}/{total_libraries}] Completed library '{lib_label}' in {library_elapsed:.1f}s",
            **progress_kwargs,
        )

        corrected_subset.uns.pop("soupx_correction", None)
        corrected_libs.append(corrected_subset)

    if not corrected_libs:
        raise RuntimeError("No valid libraries processed.")

    if len(corrected_libs) == 1:
        merged = corrected_libs[0]
    else:
        merged = ad.concat(corrected_libs, join="outer", merge="same")

    summary_table = pd.DataFrame(summary) if summary else pd.DataFrame()
    merged.uns.setdefault("soupx_correction", {})
    merged.uns["soupx_correction"].update(
        {
            "default_rho": float(rho),
            "library_column": library_col,
            "libraries": summary_table,
            "cluster_column": cluster_col or "",
            "auto_cluster": bool(auto_cluster),
            "leiden_parameters": {
                "resolution": float(leiden_resolution),
                "n_neighbors": int(leiden_neighbors),
                "n_pcs": int(leiden_npcs),
                "n_top_genes": int(leiden_hvg),
                "random_state": int(leiden_random_state),
            },
            "method": method,
            "round_to_int": bool(round_to_int),
        }
    )

    merged_path = outdir / merged_filename
    report_status(f"Writing merged corrected dataset to {merged_path}")
    merged_write_start = time.perf_counter()
    merged.write_h5ad(merged_path, compression="gzip")
    merged_write_time = time.perf_counter() - merged_write_start
    report_status(f"Merged dataset written in {merged_write_time:.1f}s")

    if not summary_table.empty:
        summary_path = outdir / "soupx_summary.tsv"
        summary_table.to_csv(summary_path, sep="\t", index=False)
        report_status(f"Summary table written to {summary_path}")

    return merged


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
) -> ad.AnnData:
    """Process all libraries within an h5ad file and write corrected outputs."""

    if not h5ad_path.exists():
        raise FileNotFoundError(f"Input h5ad not found: {h5ad_path}")

    rho_mapping = rho_mapping or {}
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
) -> Optional[ad.AnnData]:
    """Convenience wrapper for programmatic usage."""

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
    )
    if return_adata:
        return corrected
    return None


def build_cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Apply SoupX correction using a provided contamination fraction (rho)."
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
        help="Optional CSV/TSV file with columns 'library' and 'rho' for per-library contamination values.",
    )
    parser.add_argument(
        "--library-col",
        default="Library",
        help="Column name in .obs specifying library identifiers.",
    )
    parser.add_argument(
        "--cluster-col",
        help="Optional column in .obs supplying precomputed cluster labels for SoupX.",
    )
    parser.add_argument(
        "--disable-auto-cluster",
        dest="auto_cluster",
        action="store_false",
        help="Disable automatic per-library Leiden clustering.",
    )
    parser.add_argument(
        "--leiden-resolution",
        type=float,
        default=0.5,
        help="Resolution for automatic Leiden clustering (default: 0.5).",
    )
    parser.add_argument(
        "--leiden-neighbors",
        type=int,
        default=15,
        help="Number of neighbors for Leiden clustering (default: 15).",
    )
    parser.add_argument(
        "--leiden-npcs",
        type=int,
        default=30,
        help="Number of principal components for Leiden clustering (default: 30).",
    )
    parser.add_argument(
        "--leiden-hvg",
        type=int,
        default=3000,
        help="Number of highly variable genes for Leiden clustering (default: 3000).",
    )
    parser.add_argument(
        "--leiden-random-state",
        type=int,
        default=0,
        help="Random seed for Leiden clustering (default: 0).",
    )
    parser.add_argument(
        "--method",
        choices=["subtraction", "multinomial", "soupOnly"],
        default="subtraction",
        help="SoupX correction method (default: subtraction).",
    )
    parser.add_argument(
        "--round-to-int",
        action="store_true",
        help="Enable SoupX stochastic rounding to integers.",
    )
    parser.add_argument("--outdir", required=True, help="Output directory for corrected h5ad files.")
    parser.add_argument(
        "--merged-name",
        default="soupx_corrected_merged.h5ad",
        help="Filename for the merged corrected h5ad (written inside outdir).",
    )
    parser.add_argument(
        "--no-store-raw-layer",
        dest="store_raw_layer",
        action="store_false",
        help="Do not persist the original counts in a 'soupx_raw' layer.",
    )
    parser.add_argument(
        "--keep-original-x",
        dest="replace_x",
        action="store_false",
        help="Keep the original counts in .X and only write corrected counts to the layer.",
    )
    parser.add_argument(
        "--skip-individual",
        dest="write_individual",
        action="store_false",
        help="Do not emit per-library corrected h5ad files.",
    )
    parser.set_defaults(store_raw_layer=True, replace_x=True, write_individual=True, auto_cluster=True)
    return parser


def main(argv: Optional[Iterable[str]] = None) -> None:
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
    )


if __name__ == "__main__":
    main()
