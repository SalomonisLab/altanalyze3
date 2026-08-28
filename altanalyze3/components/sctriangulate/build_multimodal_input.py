#!/usr/bin/env python3
"""Assemble one scTriangulate input object from an RNA h5ad, an optional other-modality h5ad
and any number of external annotation tables.

scTriangulate reconciles several competing annotations of the SAME cells, so every annotation
must reach every cell. Cluster labels arrive from different programs in different files: ICGS3
writes ``icgs3_cell_barcode_clusters.tsv``, leiden_cluster writes
``unsupervised_leiden_clusters.tsv``, cellHarmony-lite writes
``cellHarmony_lite_assignments.txt``, and a reference call may already sit in ``obs``. This
module attaches them all under names you choose, joins the second modality through the existing
``preprocessing.concat_rna_and_other`` with its ``AB_`` convention, and reports what each join
kept.

``--require-complete`` keeps only the cells every named annotation reached. Use it: a cell with a
missing label would otherwise enter scTriangulate under an invented "unassigned" cluster that
competes for cells in every stability metric.

Entry point::

    python -m altanalyze3.components.sctriangulate.build_multimodal_input \
        --rna-h5ad rna.h5ad --adt-h5ad adt.h5ad \
        --annotation Leiden_RNA_r2 leiden/unsupervised_leiden_clusters.tsv Leiden \
        --obs-annotation Ferchen_L4 Mm-MarrowAtlas-L4 \
        --require-complete --out combined.h5ad
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

DEFAULT_ADT_PREFIX = "AB_"
ADT_NORMALIZATIONS = ("log1p", "none")


def read_annotation_table(path: str, column: Optional[str] = None) -> pd.Series:
    """Read a barcode -> label table. Column 1 is the barcode; ``column`` names the label column.

    Without ``column`` the second column is used, which is what a 2-column
    ``barcode<tab>cluster`` export holds. A header is required, because the label column is
    addressed by name everywhere else in this module.
    """
    table = pd.read_csv(path, sep="\t", dtype=str)
    if table.shape[1] < 2:
        raise ValueError(f"{path} holds {table.shape[1]} column(s); a barcode and a label are needed")
    barcode_col = table.columns[0]
    if column is None:
        column = table.columns[1]
    if column not in table.columns:
        raise KeyError(
            f"column '{column}' absent from {path}. Available: {', '.join(map(str, table.columns))}"
        )
    table = table[[barcode_col, column]].dropna()
    duplicated = int(table[barcode_col].duplicated().sum())
    if duplicated:
        raise ValueError(f"{path} repeats {duplicated} barcode(s) in column '{barcode_col}'")
    return pd.Series(table[column].astype(str).values, index=table[barcode_col].astype(str).values)


def build(
    rna_h5ad: str,
    out_path: str,
    *,
    adt_h5ad: Optional[str] = None,
    adt_prefix: str = DEFAULT_ADT_PREFIX,
    adt_normalization: str = "log1p",
    file_annotations: Sequence[Sequence[str]] = (),
    obs_annotations: Sequence[Sequence[str]] = (),
    require_complete: bool = False,
    compression: Optional[str] = "gzip",
    log=print,
) -> Dict[str, object]:
    import anndata as ad
    import scanpy as sc
    from scipy.sparse import csr_matrix, issparse

    if adt_normalization not in ADT_NORMALIZATIONS:
        raise ValueError(
            f"Unknown adt normalization '{adt_normalization}'; choose from "
            f"{', '.join(ADT_NORMALIZATIONS)}"
        )

    rna = sc.read_h5ad(rna_h5ad)
    log(f"[build] RNA {rna.n_obs} cells x {rna.n_vars} features from {rna_h5ad}")
    report: Dict[str, object] = {
        "rna_h5ad": os.path.abspath(rna_h5ad),
        "rna_cells": int(rna.n_obs),
        "rna_features": int(rna.n_vars),
        "annotations": {},
    }

    # ---- attach every annotation onto the RNA cells -------------------------------------------
    names: List[str] = []
    for name, path, column in file_annotations:
        series = read_annotation_table(path, column or None)
        mapped = series.reindex(rna.obs_names.astype(str))
        covered = int(mapped.notna().sum())
        rna.obs[name] = mapped.values
        names.append(name)
        log(f"[build] annotation '{name}': {covered} of {rna.n_obs} RNA cells "
            f"({100.0 * covered / rna.n_obs:.2f}%), {mapped.dropna().nunique()} labels, from {path}")
        report["annotations"][name] = {
            "source": os.path.abspath(path), "column": column or "<second column>",
            "cells_covered": covered, "n_labels": int(mapped.dropna().nunique()),
        }

    for name, obs_column in obs_annotations:
        if obs_column not in rna.obs.columns:
            raise KeyError(
                f"obs column '{obs_column}' absent from {rna_h5ad}. "
                f"Available: {', '.join(map(str, rna.obs.columns))}"
            )
        values = rna.obs[obs_column].astype(str)
        rna.obs[name] = values.values
        names.append(name)
        covered = int(values.notna().sum())
        log(f"[build] annotation '{name}': {covered} of {rna.n_obs} RNA cells from obs['{obs_column}'], "
            f"{values.nunique()} labels")
        report["annotations"][name] = {
            "source": f"obs['{obs_column}']", "column": obs_column,
            "cells_covered": covered, "n_labels": int(values.nunique()),
        }

    # ---- the second modality ------------------------------------------------------------------
    other = None
    if adt_h5ad:
        other = sc.read_h5ad(adt_h5ad)
        log(f"[build] ADT {other.n_obs} cells x {other.n_vars} features from {adt_h5ad}")
        report["adt_h5ad"] = os.path.abspath(adt_h5ad)
        report["adt_cells"] = int(other.n_obs)
        report["adt_features"] = int(other.n_vars)
        if adt_normalization == "log1p":
            # The ADT panel holds linear TotalVI values whose maximum is two orders of magnitude
            # above the log RNA scale. Concatenating them unlogged lets a few bright antibodies
            # dominate every distance the downstream metrics compute. normalize_total is NOT
            # applied: CP10K over a small panel forces each cell's whole panel to one constant total.
            sc.pp.log1p(other)
            log(f"[build] ADT normalization: log1p (no depth step)")
        else:
            log(f"[build] ADT normalization: none")
        report["adt_normalization"] = adt_normalization
        shared = rna.obs_names.intersection(other.obs_names)
        log(f"[build] ADT covers {len(shared)} of {rna.n_obs} RNA cells "
            f"({100.0 * len(shared) / rna.n_obs:.2f}%)")

    # ---- the cell set --------------------------------------------------------------------------
    keep = pd.Series(True, index=rna.obs_names)
    if require_complete:
        for name in names:
            present = rna.obs[name].notna() & ~rna.obs[name].astype(str).isin(["nan", "None", ""])
            lost = int((keep & ~present).sum())
            keep &= present
            log(f"[build] require-complete after '{name}': {int(keep.sum())} cells "
                f"({lost} newly dropped)")
        if other is not None:
            present = rna.obs_names.isin(other.obs_names)
            lost = int((keep & ~present).sum())
            keep &= present
            log(f"[build] require-complete after the ADT join: {int(keep.sum())} cells "
                f"({lost} newly dropped)")
    elif other is not None:
        keep &= rna.obs_names.isin(other.obs_names)
        log(f"[build] restricted to the {int(keep.sum())} cells the ADT panel also holds")

    if int(keep.sum()) == 0:
        raise ValueError("no cell carries every annotation; nothing to write")
    rna = rna[keep.values].copy()
    log(f"[build] cell set: {rna.n_obs} of {report['rna_cells']} RNA cells "
        f"({100.0 * rna.n_obs / report['rna_cells']:.2f}%)")
    report["cells_kept"] = int(rna.n_obs)
    report["require_complete"] = bool(require_complete)

    for name in names:
        rna.obs[name] = pd.Categorical(rna.obs[name].astype(str))
        report["annotations"][name]["n_labels_kept"] = int(rna.obs[name].nunique())

    # ---- concatenate the modalities through the existing scTriangulate helper -------------------
    if other is not None:
        from .preprocessing import concat_rna_and_other
        other = other[rna.obs_names].copy()
        if "X_umap" not in rna.obsm:
            rna.obsm["X_umap"] = np.zeros((rna.n_obs, 2), dtype=np.float32)
        combined = concat_rna_and_other(
            rna, other, umap="rna", umap_key="X_umap", name="adt", prefix=adt_prefix)
        combined.obs = rna.obs.copy()
        n_ab = int(sum(str(v).startswith(adt_prefix) for v in combined.var_names))
        if n_ab != other.n_vars:
            raise ValueError(f"{n_ab} '{adt_prefix}' features in the combined object, "
                             f"expected {other.n_vars}")
        log(f"[build] combined {combined.n_obs} cells x {combined.n_vars} features "
            f"({combined.n_vars - n_ab} RNA + {n_ab} {adt_prefix})")
    else:
        combined = rna

    if not issparse(combined.X):
        combined.X = csr_matrix(combined.X)
    if combined.n_obs != int(keep.sum()):
        raise ValueError(f"combined object holds {combined.n_obs} cells, expected {int(keep.sum())}")
    for name in names:
        if combined.obs[name].isna().any():
            raise ValueError(f"annotation '{name}' still holds NaN after the cell set was fixed")

    report["combined_cells"] = int(combined.n_obs)
    report["combined_features"] = int(combined.n_vars)

    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    combined.write(out_path, compression=compression)
    log(f"[build] wrote {combined.n_obs} x {combined.n_vars} -> {os.path.abspath(out_path)}")

    labels_path = os.path.splitext(out_path)[0] + "_annotations.tsv"
    combined.obs[names].to_csv(labels_path, sep="\t")
    log(f"[build] wrote the annotation matrix -> {os.path.abspath(labels_path)}")
    report["out_h5ad"] = os.path.abspath(out_path)
    report["annotation_matrix"] = os.path.abspath(labels_path)
    report["query"] = ",".join(names)

    report_path = os.path.splitext(out_path)[0] + "_build_report.json"
    with open(report_path, "w") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
    log(f"[build] wrote the build report -> {os.path.abspath(report_path)}")
    log(f"[build] scTriangulate --query {report['query']}")
    return report


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Assemble one scTriangulate input from an RNA h5ad, an optional second "
                    "modality and external annotation tables.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--rna-h5ad", required=True, help="RNA h5ad; its cells define the object.")
    parser.add_argument("--out", required=True, help="Output h5ad.")
    parser.add_argument("--adt-h5ad", default=None,
                        help="Second-modality h5ad joined on the barcode.")
    parser.add_argument("--adt-prefix", default=DEFAULT_ADT_PREFIX,
                        help="Prefix put in front of every second-modality feature.")
    parser.add_argument("--adt-normalization", default="log1p", choices=list(ADT_NORMALIZATIONS),
                        help="log1p brings a linear TotalVI panel onto the log RNA scale. No "
                             "depth step is applied, because CP10K over a small panel forces "
                             "every cell's whole panel to one constant total.")
    parser.add_argument("--annotation", nargs=3, action="append", default=[],
                        metavar=("NAME", "TSV", "COLUMN"),
                        help="Attach a barcode->label TSV as obs[NAME]. COLUMN names the label "
                             "column; pass '' for the second column. Repeatable.")
    parser.add_argument("--obs-annotation", nargs=2, action="append", default=[],
                        metavar=("NAME", "OBS_COLUMN"),
                        help="Copy an existing obs column of the RNA h5ad to obs[NAME]. Repeatable.")
    parser.add_argument("--require-complete", action="store_true",
                        help="Keep only the cells every named annotation reached.")
    parser.add_argument("--compression", default="gzip", choices=["gzip", "lzf", "none"])
    args = parser.parse_args(argv)

    build(
        args.rna_h5ad, args.out,
        adt_h5ad=args.adt_h5ad, adt_prefix=args.adt_prefix,
        adt_normalization=args.adt_normalization,
        file_annotations=args.annotation, obs_annotations=args.obs_annotation,
        require_complete=args.require_complete,
        compression=None if args.compression == "none" else args.compression,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
