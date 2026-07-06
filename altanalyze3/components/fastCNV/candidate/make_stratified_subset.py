"""Create a row-stratified AnnData subset for candidate fastCNV evaluation."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd

from altanalyze3.components.fastCNV.candidate.per_gene import _matrix_for_adata, _subset_matrix


def stratified_rows(
    obs: pd.DataFrame,
    group_keys: Sequence[str],
    max_per_group: int,
    *,
    random_state: int,
) -> np.ndarray:
    rng = np.random.default_rng(random_state)
    selected: list[np.ndarray] = []
    indexed = obs.reset_index(drop=True)
    for _, rows in indexed.groupby(list(group_keys), observed=False, sort=True).groups.items():
        row_array = np.asarray(rows, dtype=np.int64)
        if row_array.size > max_per_group:
            row_array = np.sort(rng.choice(row_array, size=max_per_group, replace=False))
        selected.append(row_array)
    if not selected:
        return np.array([], dtype=np.int64)
    return np.sort(np.concatenate(selected).astype(np.int64))


def create_subset(
    h5ad: Path,
    output_h5ad: Path,
    *,
    sample_key: str,
    state_key: str,
    max_per_group: int,
    layer: str,
    random_state: int,
) -> dict:
    source = ad.read_h5ad(h5ad, backed="r")
    try:
        missing = [key for key in (sample_key, state_key) if key not in source.obs.columns]
        if missing:
            raise KeyError(f"Missing obs columns: {missing}")
        rows = stratified_rows(
            source.obs[[sample_key, state_key]].astype(str),
            [sample_key, state_key],
            max_per_group,
            random_state=random_state,
        )
        if rows.size == 0:
            raise ValueError("No rows selected.")
        cols = np.arange(source.n_vars, dtype=np.int64)
        matrix = _subset_matrix(_matrix_for_adata(source, layer), rows, cols)
        subset = ad.AnnData(
            X=matrix,
            obs=source.obs.iloc[rows].copy(),
            var=source.var.copy(),
        )
        subset.uns["fastcnv_subset"] = {
            "source_h5ad": str(h5ad),
            "sample_key": sample_key,
            "state_key": state_key,
            "max_per_group": int(max_per_group),
            "random_state": int(random_state),
            "n_source_cells": int(source.n_obs),
        }
        output_h5ad.parent.mkdir(parents=True, exist_ok=True)
        subset.write_h5ad(output_h5ad, compression="gzip")
        counts = (
            subset.obs.groupby([sample_key, state_key], observed=False)
            .size()
            .rename("n_cells")
            .reset_index()
        )
        counts.to_csv(f"{output_h5ad}.group_counts.tsv", sep="\t", index=False)
        return {
            "source_cells": int(source.n_obs),
            "subset_cells": int(subset.n_obs),
            "genes": int(subset.n_vars),
            "groups": int(counts.shape[0]),
            "output_h5ad": str(output_h5ad),
            "group_counts": f"{output_h5ad}.group_counts.tsv",
        }
    finally:
        source.file.close()


def create_selection(
    h5ad: Path,
    output_tsv: Path,
    *,
    sample_key: str,
    state_key: str,
    max_per_group: int,
    random_state: int,
) -> dict:
    source = ad.read_h5ad(h5ad, backed="r")
    try:
        missing = [key for key in (sample_key, state_key) if key not in source.obs.columns]
        if missing:
            raise KeyError(f"Missing obs columns: {missing}")
        rows = stratified_rows(
            source.obs[[sample_key, state_key]].astype(str),
            [sample_key, state_key],
            max_per_group,
            random_state=random_state,
        )
        if rows.size == 0:
            raise ValueError("No rows selected.")
        selection = source.obs.iloc[rows][[sample_key, state_key]].copy()
        selection.insert(0, "CellBarcode", source.obs_names[rows].astype(str))
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        selection.to_csv(output_tsv, sep="\t", index=False)
        counts = (
            selection.groupby([sample_key, state_key], observed=False)
            .size()
            .rename("n_cells")
            .reset_index()
        )
        counts.to_csv(f"{output_tsv}.group_counts.tsv", sep="\t", index=False)
        return {
            "source_cells": int(source.n_obs),
            "selected_cells": int(selection.shape[0]),
            "groups": int(counts.shape[0]),
            "output_tsv": str(output_tsv),
            "group_counts": f"{output_tsv}.group_counts.tsv",
        }
    finally:
        source.file.close()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Create a stratified H5AD subset.")
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--output-h5ad", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path)
    parser.add_argument("--selection-only", action="store_true")
    parser.add_argument("--sample-key", default="Sample")
    parser.add_argument("--state-key", default="Hs-BM-titrated-reference-centroid")
    parser.add_argument("--max-per-group", type=int, default=200)
    parser.add_argument("--layer", default="X")
    parser.add_argument("--random-state", type=int, default=0)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.selection_only:
        output_tsv = args.output_tsv or Path(f"{args.output_h5ad}.selected_cells.tsv")
        summary = create_selection(
            args.h5ad,
            output_tsv,
            sample_key=args.sample_key,
            state_key=args.state_key,
            max_per_group=args.max_per_group,
            random_state=args.random_state,
        )
    else:
        summary = create_subset(
            args.h5ad,
            args.output_h5ad,
            sample_key=args.sample_key,
            state_key=args.state_key,
            max_per_group=args.max_per_group,
            layer=args.layer,
            random_state=args.random_state,
        )
    print(pd.Series(summary).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
