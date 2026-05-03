"""Downsample a reference h5ad to N cells per group, optionally renaming an obs column.

Typical use: scale down a single-cell atlas before passing it as the
``--control-h5ad`` reference to ``fastCNV.py``. A 4-5x cell-count reduction
roughly halves the fastCNV scoring time without changing the per-cell-state
median+MAD baseline.

Example:
    python3 -m altanalyze3.components.fastCNV.subset_reference \\
      --input /path/to/atlas.h5ad \\
      --output /path/to/atlas.subset50.h5ad \\
      --group-cols Donor "Level 3 Multimodal" \\
      --cells-per-group 50 \\
      --rename "Level 3 Multimodal=Hs-MarrowAtlas-L3M"
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd


LOGGER = logging.getLogger("fastCNV.subset_reference")


def _parse_rename(spec: str) -> Tuple[str, str]:
    if "=" not in spec:
        raise argparse.ArgumentTypeError(
            f"--rename expects 'old=new' (got: {spec!r})"
        )
    old, new = spec.split("=", 1)
    old, new = old.strip(), new.strip()
    if not old or not new:
        raise argparse.ArgumentTypeError(
            f"--rename requires non-empty names on both sides of '=' (got: {spec!r})"
        )
    return old, new


def subset_reference(
    input_path: Path,
    output_path: Path,
    group_cols: Sequence[str],
    cells_per_group: int,
    renames: Sequence[Tuple[str, str]] = (),
    seed: int = 42,
) -> Path:
    LOGGER.info("Reading: %s", input_path)
    t0 = time.perf_counter()
    adata = ad.read_h5ad(input_path)
    LOGGER.info("  loaded %s in %.1fs", adata.shape, time.perf_counter() - t0)

    available = list(adata.obs.columns)
    missing = [c for c in group_cols if c not in adata.obs.columns]
    if missing:
        raise KeyError(
            f"Missing obs columns: {missing}. Available: {available}"
        )
    for old, _ in renames:
        if old not in adata.obs.columns:
            raise KeyError(
                f"--rename source column {old!r} not found. Available: {available}"
            )

    obs_subset = adata.obs[list(group_cols)].copy()
    mask = pd.Series(True, index=obs_subset.index)
    for col in group_cols:
        col_values = obs_subset[col]
        mask &= col_values.notna() & (col_values.astype(str).str.len() > 0)
    n_dropped = int((~mask).sum())
    if n_dropped:
        LOGGER.info("Dropped %d cells missing one or more group columns", n_dropped)
    obs_subset = obs_subset.loc[mask]

    rng = np.random.default_rng(seed)
    keep_indices: List[object] = []
    n_groups = 0
    n_full = 0
    n_truncated = 0
    for _, grp in obs_subset.groupby(list(group_cols), sort=False, observed=True):
        n_groups += 1
        cell_ids = grp.index.to_numpy()
        n = len(cell_ids)
        take = min(cells_per_group, n)
        if take < n:
            chosen = rng.choice(cell_ids, size=take, replace=False)
            n_full += 1
        else:
            chosen = cell_ids
            n_truncated += 1
        keep_indices.extend(chosen)
    LOGGER.info(
        "Selected %d cells across %d groups (%d at full %d, %d truncated by availability)",
        len(keep_indices), n_groups, n_full, cells_per_group, n_truncated,
    )

    keep_index = pd.Index(keep_indices)
    subset = adata[keep_index].copy()
    for old, new in renames:
        subset.obs[new] = subset.obs[old].astype(str)
        LOGGER.info("Added obs column %r (copy of %r as string)", new, old)

    LOGGER.info("Writing: %s", output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    subset.write_h5ad(output_path)
    LOGGER.info(
        "  wrote %s in %.1fs (%.2f GB)",
        subset.shape, time.perf_counter() - t0, output_path.stat().st_size / 1e9,
    )
    return output_path


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Downsample a reference h5ad to N cells per group, optionally renaming an obs column.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input", required=True, type=Path, help="Source h5ad.")
    p.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Destination h5ad. Default: <input>.subset<N>.h5ad next to input.",
    )
    p.add_argument(
        "--group-cols",
        nargs="+",
        required=True,
        help="One or more obs columns defining the group; cells are downsampled within each unique combination.",
    )
    p.add_argument("--cells-per-group", type=int, default=50)
    p.add_argument(
        "--rename",
        action="append",
        type=_parse_rename,
        default=[],
        help="Add a string-typed obs column copied from another. Format: 'old=new'. May be repeated.",
    )
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> Path:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s | %(message)s",
    )
    output = args.output
    if output is None:
        output = args.input.with_name(
            f"{args.input.stem}.subset{args.cells_per_group}.h5ad"
        )
    return subset_reference(
        input_path=args.input,
        output_path=output,
        group_cols=args.group_cols,
        cells_per_group=args.cells_per_group,
        renames=args.rename,
        seed=args.seed,
    )


if __name__ == "__main__":
    main(sys.argv[1:])
