#!/usr/bin/env python3
"""Collapse one or more obs columns into named groups, and write an obs sidecar h5ad.

A reference annotation with 51 states cannot be read as a heatmap colour bar, and
``marker_heatmap_h5ad`` drops any covariate above 12 displayed categories. This module maps the
states a reader cares about onto a handful of named groups, sends everything else to one
``Other`` category, and writes the result as a small h5ad the heatmap CLI accepts through
``--covariate-h5ad``. The sidecar carries obs only, so a 308 MB source object yields a file of a
few hundred kilobytes.

The mapping arrives as a 4-column TSV, which makes it a reproducibility artifact rather than a
buried command line::

    out_key        source_key  group_name   source_label
    StJude_group   StJude      DN           DN
    Ferchen_group  Ferchen_L4  CLP-1        CLP1-a

Entry point::

    python -m altanalyze3.components.aggregate.h5ad_group_obs \
        --h5ad <source.h5ad> --mapping <groups.tsv> --out <covariate_groups.h5ad>

Every run prints the cell count of each group and names every source label it did not find.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

MAPPING_COLUMNS = ("out_key", "source_key", "group_name", "source_label")


def read_mapping(path: str) -> pd.DataFrame:
    """Read the 4-column group mapping and reject a duplicate or an incomplete row."""
    table = pd.read_csv(path, sep="\t", dtype=str, comment="#").dropna(how="all")
    missing = [c for c in MAPPING_COLUMNS if c not in table.columns]
    if missing:
        raise ValueError(
            f"{path} lacks the column(s) {missing}. Required: {', '.join(MAPPING_COLUMNS)}"
        )
    table = table[list(MAPPING_COLUMNS)].apply(lambda s: s.astype(str).str.strip())
    blank = table.eq("").any(axis=1)
    if blank.any():
        raise ValueError(f"{path} holds {int(blank.sum())} row(s) with an empty field")
    duplicated = table.duplicated(subset=["out_key", "source_label"], keep=False)
    if duplicated.any():
        offending = table.loc[duplicated].to_string(index=False)
        raise ValueError(
            f"{path} sends one source_label to two groups of the same out_key:\n{offending}"
        )
    for out_key, block in table.groupby("out_key"):
        sources = sorted(block["source_key"].unique())
        if len(sources) != 1:
            raise ValueError(
                f"out_key '{out_key}' draws from more than one source_key: {sources}. "
                "One grouped column must read one obs column."
            )
    return table


def group_obs(
    h5ad_path: str,
    mapping_path: str,
    out_path: str,
    *,
    other_label: str = "Other",
    allow_missing_labels: bool = False,
    compression: Optional[str] = "gzip",
    log=print,
) -> pd.DataFrame:
    """Write an obs sidecar h5ad carrying the source obs plus one column per ``out_key``."""
    import anndata as ad
    from scipy.sparse import csr_matrix

    mapping = read_mapping(mapping_path)
    adata = ad.read_h5ad(h5ad_path, backed="r")
    obs = adata.obs.copy()
    n_obs = int(obs.shape[0])
    log(f"[group] {h5ad_path}: {n_obs} cells, {obs.shape[1]} obs columns")

    absent_total: List[str] = []
    for out_key, block in mapping.groupby("out_key", sort=False):
        source_key = block["source_key"].iloc[0]
        if source_key not in obs.columns:
            raise KeyError(
                f"obs column '{source_key}' absent from {h5ad_path}. "
                f"Available: {', '.join(map(str, obs.columns))}"
            )
        values = obs[source_key].astype(str)
        present = set(values.unique())
        absent = [r.source_label for r in block.itertuples() if r.source_label not in present]
        if absent:
            absent_total.extend(f"{out_key}:{a}" for a in absent)
            message = (f"[group] {out_key}: obs['{source_key}'] holds no cell labelled "
                       f"{absent}; those groups stay empty")
            if not allow_missing_labels:
                raise ValueError(
                    message.replace("[group] ", "")
                    + ". Pass --allow-missing-labels to keep going, or fix the mapping."
                )
            log(message)

        lookup: Dict[str, str] = dict(zip(block["source_label"], block["group_name"]))
        grouped = values.map(lookup).fillna(other_label)
        order = list(dict.fromkeys(block["group_name"])) + [other_label]
        order = [g for g in order if (grouped == g).any()]
        obs[out_key] = pd.Categorical(grouped, categories=order, ordered=True)

        counts = obs[out_key].value_counts().reindex(order)
        log(f"[group] {out_key} from obs['{source_key}'] -> {len(order)} categories")
        for name, count in counts.items():
            log(f"[group]     {name:<16s} {int(count):>6d} of {n_obs} "
                f"({100.0 * int(count) / n_obs:.2f}%)")
        named = int((obs[out_key] != other_label).sum())
        log(f"[group]     named {named} of {n_obs} ({100.0 * named / n_obs:.2f}%), "
            f"{n_obs - named} fell to '{other_label}'")
        if int(counts.sum()) != n_obs:
            raise ValueError(f"{out_key} counts sum to {int(counts.sum())}, expected {n_obs}")

    out_keys = list(dict.fromkeys(mapping["out_key"]))
    sidecar = ad.AnnData(
        X=csr_matrix((n_obs, 1), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(["_placeholder"])),
    )
    if not sidecar.obs_names.equals(adata.obs_names):
        raise ValueError("sidecar barcode order differs from the source h5ad")
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    sidecar.write(out_path, compression=compression)
    log(f"[group] wrote {n_obs} cells x {len(obs.columns)} obs columns -> {os.path.abspath(out_path)}")
    log(f"[group] new columns: {', '.join(out_keys)}")
    if absent_total:
        log(f"[group] source labels not found: {', '.join(absent_total)}")

    table_path = os.path.splitext(out_path)[0] + "_group_counts.tsv"
    rows = []
    for out_key in out_keys:
        for name, count in obs[out_key].value_counts().reindex(
                obs[out_key].cat.categories).items():
            rows.append({"grouped_column": out_key, "group": name, "cells": int(count),
                         "fraction_of_cells": int(count) / n_obs})
    counts_table = pd.DataFrame(rows)
    counts_table.to_csv(table_path, sep="\t", index=False, float_format="%.6g")
    log(f"[group] wrote the group counts -> {os.path.abspath(table_path)}")
    return counts_table


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Collapse obs columns into named groups and write an obs sidecar h5ad.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--h5ad", required=True, help="Source h5ad; only its obs is read.")
    parser.add_argument("--mapping", required=True,
                        help="4-column TSV: out_key, source_key, group_name, source_label.")
    parser.add_argument("--out", required=True, help="Output sidecar h5ad.")
    parser.add_argument("--other-label", default="Other",
                        help="Category every unmapped label falls into.")
    parser.add_argument("--allow-missing-labels", action="store_true",
                        help="Warn instead of failing when the mapping names a label the obs "
                             "column does not hold.")
    parser.add_argument("--compression", default="gzip", choices=["gzip", "lzf", "none"])
    args = parser.parse_args(argv)

    group_obs(
        args.h5ad, args.mapping, args.out,
        other_label=args.other_label,
        allow_missing_labels=args.allow_missing_labels,
        compression=None if args.compression == "none" else args.compression,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
