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
* Processes each library independently to avoid cross-library bleed.
* Stores the raw counts in a layer (``soupx_raw``) and overwrites ``.X`` with the
  corrected counts unless disabled.
* Persists corrected counts in a ``soupx_corrected`` layer for downstream audit.
* Emits progress updates via ``tqdr`` when available (falls back to console
  prints otherwise).

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
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, Iterable, Optional

import anndata as ad
import numpy as np
import scipy.sparse as sp

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
def status_scope(message: str, done_message: Optional[str] = None):
    """Emit a start message and optionally a completion message."""

    report_status(message)
    try:
        yield
    except Exception:
        report_status(f"FAILED: {message}")
        raise
    else:
        if done_message:
            report_status(done_message)


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


def correct_library(
    subset: ad.AnnData,
    *,
    rho: float,
    store_raw_layer: bool = True,
    replace_x: bool = True,
    raw_layer: str = RAW_LAYER_NAME,
    corrected_layer: str = CORRECTED_LAYER_NAME,
    library_label: Optional[str] = None,
) -> ad.AnnData:
    """Apply SoupX correction to a single library AnnData in-place."""

    if subset.X is None or subset.n_obs == 0:
        raise ValueError("Empty library subset provided.")

    counts_csr = _to_csr(subset.X)
    if store_raw_layer:
        subset.layers[raw_layer] = counts_csr.copy()

    sc = soupx.SoupChannel(toc=counts_csr.T)
    sc.rho = validate_rho(rho)
    corrected = soupx.adjustCounts(sc)
    corrected_csr = _to_csr(corrected).T

    subset.layers[corrected_layer] = corrected_csr
    if replace_x:
        subset.X = corrected_csr.copy()

    subset.obs["soupx_rho"] = rho
    if library_label is not None:
        subset.obs["soupx_library"] = library_label

    metadata = {
        "library": library_label or "unknown",
        "rho": rho,
        "soupx_version": getattr(soupx, "__version__", "unknown"),
        "cells": int(subset.n_obs),
        "genes": int(subset.n_vars),
    }
    history = subset.uns.get("soupx_correction", [])
    if not isinstance(history, list):
        history = [history]
    history.append(metadata)
    subset.uns["soupx_correction"] = history
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


def process_h5ad(
    h5ad_path: Path,
    *,
    rho: float = DEFAULT_RHO,
    library_col: str = "library",
    outdir: Path = Path("./soupx_corrected"),
    rho_mapping: Optional[Dict[str, float]] = None,
    store_raw_layer: bool = True,
    replace_x: bool = True,
    write_individual: bool = True,
    merged_filename: str = "soupx_corrected_merged.h5ad",
) -> ad.AnnData:
    """Process all libraries within an h5ad file and write corrected outputs."""

    if not h5ad_path.exists():
        raise FileNotFoundError(f"Input h5ad not found: {h5ad_path}")

    rho_mapping = rho_mapping or {}
    adata = ad.read_h5ad(h5ad_path)

    if library_col not in adata.obs:
        raise ValueError(f"obs column '{library_col}' not found in h5ad.")

    outdir.mkdir(parents=True, exist_ok=True)

    libraries = list(_iter_libraries(adata.obs[library_col].tolist()))
    if not libraries:
        raise ValueError(f"No library values found in column '{library_col}'.")

    report_status(
        f"Loaded {h5ad_path.name} with {adata.n_obs} cells, {adata.n_vars} genes, {len(libraries)} libraries."
    )

    corrected_libs = []
    summary = []

    for lib in libraries:
        lib_mask = adata.obs[library_col] == lib
        subset = adata[lib_mask].copy()
        if subset.n_obs == 0:
            report_status(f"Skipping empty library: {lib}")
            continue

        lib_label = str(lib)
        lib_rho = rho_mapping.get(lib_label, rho)
        lib_rho = validate_rho(lib_rho)

        with status_scope(
            f"Correcting library '{lib_label}' (rho={lib_rho:.4f})",
            done_message=f"Finished library '{lib_label}'",
        ):
            corrected_subset = correct_library(
                subset,
                rho=lib_rho,
                store_raw_layer=store_raw_layer,
                replace_x=replace_x,
                library_label=lib_label,
            )

        summary.append({
            "library": lib_label,
            "rho": lib_rho,
            "cells": int(corrected_subset.n_obs),
            "genes": int(corrected_subset.n_vars),
        })

        if write_individual:
            safe_name = _safe_filename(lib_label)
            out_file = outdir / f"{safe_name}_soupx_corrected.h5ad"
            corrected_subset.write_h5ad(out_file, compression="gzip")
            report_status(f"Wrote corrected library dataset: {out_file}")

        corrected_subset.uns.pop("soupx_correction", None)
        corrected_libs.append(corrected_subset)

    if not corrected_libs:
        raise RuntimeError("No valid libraries processed.")

    if len(corrected_libs) == 1:
        merged = corrected_libs[0]
    else:
        merged = ad.concat(corrected_libs, join="outer", merge="same")

    merged.uns.setdefault("soupx_correction", {})
    merged.uns["soupx_correction"].update(
        {
            "default_rho": rho,
            "library_column": library_col,
            "libraries": summary,
        }
    )

    merged_path = outdir / merged_filename
    merged.write_h5ad(merged_path, compression="gzip")
    report_status(f"Wrote merged corrected dataset: {merged_path}")

    return merged


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
        default="library",
        help="Column name in .obs specifying library identifiers.",
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
    parser.set_defaults(store_raw_layer=True, replace_x=True, write_individual=True)
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
        return_adata=False,
    )


if __name__ == "__main__":
    main()
