"""Prove the lung port reproduces the known-good mouse whitelist.

``build_whitelist.top_correlates`` and ``build_whitelist.select_features`` are
ports of ``components.rna2adt.mouse.build_whitelist``. This script drives them
with the mouse atlas, the mouse curated map, the mouse seed and no
protein-coding filter, then compares every row against the committed
``mouse/configs/empirical_whitelist.tsv``. Any difference is a port defect.

    python -m altanalyze3.components.rna2adt.lung.validate_port \
      --adata /Users/saljh8/Dropbox/Manuscripts/InProgress/ALRG_Paper_2022/Azimuth-Tests/adata_combined_60k_rna_adt.h5ad
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..mouse.adt_mgi_map import load_curated_adt_mgi_map, strip_prefix as mouse_strip_prefix
from ..mouse.build_whitelist import _is_isotype
from .build_whitelist import entries_to_frame, select_features, top_correlates


MOUSE_ATLAS_DEFAULT = Path(
    "/Users/saljh8/Dropbox/Manuscripts/InProgress/ALRG_Paper_2022/Azimuth-Tests/"
    "adata_combined_60k_rna_adt.h5ad"
)
MOUSE_WHITELIST = Path(__file__).parent.parent / "mouse" / "configs" / "empirical_whitelist.tsv"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--adata", type=Path, default=MOUSE_ATLAS_DEFAULT)
    parser.add_argument("--reference", type=Path, default=MOUSE_WHITELIST)
    parser.add_argument("--max-cells", type=int, default=30000)
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--rna-chunk", type=int, default=512)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    print(f"[load] {args.adata}", flush=True)
    atlas = ad.read_h5ad(args.adata)
    print(f"[load] {atlas.n_obs} cells, {atlas.n_vars} vars", flush=True)
    rng = np.random.default_rng(args.seed)
    if atlas.n_obs > args.max_cells:
        rows = np.sort(rng.choice(atlas.n_obs, args.max_cells, replace=False))
        atlas = atlas[rows].to_memory()
        print(f"[load] subsampled to {atlas.n_obs} cells", flush=True)

    var_names = np.array([str(v) for v in atlas.var_names])
    is_adt = np.array([s.startswith("ADT-") for s in var_names])
    is_isotype = np.array([_is_isotype(s) for s in var_names])
    keep = is_adt & ~is_isotype
    adt_names = var_names[keep].tolist()
    rna_names = var_names[~is_adt].tolist()
    print(f"[load] {len(adt_names)} ADTs kept, {len(rna_names)} RNA genes", flush=True)

    Y = atlas.X[:, np.where(keep)[0]]
    Y = np.asarray(Y.todense()) if sp.issparse(Y) else np.asarray(Y)
    X = atlas.X[:, np.where(~is_adt)[0]]

    top_indices, _ = top_correlates(X, Y, top_n=args.top_n, chunk=args.rna_chunk)
    entries = select_features(adt_names, top_indices, rna_names,
                              load_curated_adt_mgi_map(), fallback_allowed=None)
    for entry in entries:
        entry["adt_clean"] = mouse_strip_prefix(entry["adt_raw"])
    produced = entries_to_frame(entries)

    expected = pd.read_csv(args.reference, sep="\t")
    print(f"\n[compare] produced {produced.shape} vs reference {expected.shape}", flush=True)
    if list(produced.columns) != list(expected.columns):
        print(f"FAIL column mismatch: {list(produced.columns)} vs {list(expected.columns)}")
        return 1
    produced = produced.sort_values("adt_raw").reset_index(drop=True)
    expected = expected.sort_values("adt_raw").reset_index(drop=True)
    if produced.shape != expected.shape:
        print("FAIL row-count mismatch")
        return 1

    def _normalise(series: pd.Series) -> np.ndarray:
        """pandas reads the committed TRUE/FALSE back as booleans, so compare
        on a canonical form rather than on the repr."""
        out = []
        for value in series.astype(str):
            token = value.strip()
            if token.lower() in {"true", "false"}:
                token = token.upper()
            out.append(token)
        return np.array(out, dtype=object)

    mismatches = []
    for column in produced.columns:
        left = pd.Series(_normalise(produced[column]))
        right = pd.Series(_normalise(expected[column]))
        bad = np.where(left.to_numpy() != right.to_numpy())[0]
        for row in bad:
            mismatches.append((column, produced.at[row, "adt_raw"], left[row], right[row]))
    n_rows = int(produced.shape[0])
    n_cells_compared = n_rows * int(produced.shape[1])
    if mismatches:
        print(f"FAIL {len(mismatches)} of {n_cells_compared} cells differ across {n_rows} rows")
        for column, adt, left, right in mismatches[:25]:
            print(f"  {adt}\t{column}\tport={left}\tref={right}")
        return 1
    print(f"PASS all {n_cells_compared} cells identical across {n_rows} ADT rows "
          f"({len(produced.columns)} columns)")
    print(produced.head(10).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
