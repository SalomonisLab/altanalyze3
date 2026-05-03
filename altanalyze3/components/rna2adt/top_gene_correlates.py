"""Compute the top-N RNA gene correlates of every ADT in the atlas.

Spearman correlation between each ADT and each RNA gene's expression
vector across all cells. Sparse-friendly via numpy ranking on densified
chunks. Useful as a sanity check on the curated ADT->RNA map: the
top-5 list for ``Hu.CD19`` should include ``CD19``; for ``Hu.CD3``
it should include ``CD3D``/``CD3E``/``CD3G``/``CD247``; etc.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import List

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def _rank_columns(arr: np.ndarray) -> np.ndarray:
    """Average-rank each column independently. arr shape (n, p) -> (n, p).

    Uses scipy.stats.rankdata vectorised over axis=0 — much faster than a
    per-column pandas loop on wide blocks.
    """
    from scipy.stats import rankdata
    arr = np.asarray(arr, dtype=np.float64)
    return rankdata(arr, axis=0).astype(np.float64, copy=False)


def compute_top_correlates(
    adata_path: Path,
    *,
    output_tsv: Path,
    top_n: int = 5,
    max_cells: int = 30000,
    rna_chunk: int = 256,
    seed: int = 0,
) -> None:
    print(f"[load] {adata_path}", flush=True)
    a = ad.read_h5ad(adata_path)
    print(f"[load] {a.n_obs} cells, {a.n_vars} vars", flush=True)

    rng = np.random.default_rng(seed)
    if a.n_obs > max_cells:
        idx = np.sort(rng.choice(a.n_obs, max_cells, replace=False))
        a = a[idx].to_memory()
        print(f"[load] subsampled to {a.n_obs} cells", flush=True)

    var_names = np.array([str(v) for v in a.var_names])
    # ADT panel = Hu.* (human ADTs) + HuMs.* (cross-reactive Hu/Mouse) +
    # Isotype_* (isotype controls). All three are antibody features and
    # should be treated as ADT outputs, not RNA inputs.
    adt_prefixes = ("Hu.", "HuMs.", "Isotype_")
    is_adt = np.array([s.startswith(adt_prefixes) for s in var_names])
    adt_names_raw = var_names[is_adt]
    rna_names = var_names[~is_adt]
    adt_idx = np.where(is_adt)[0]
    rna_idx = np.where(~is_adt)[0]

    def _strip_prefix(name: str) -> str:
        for prefix in adt_prefixes:
            if name.startswith(prefix):
                return name[len(prefix):]
        return name
    adt_names_clean = np.array([_strip_prefix(n) for n in adt_names_raw])
    print(f"[load] {len(adt_names_raw)} ADTs ({sum(s.startswith('Hu.') for s in adt_names_raw)} Hu, "
          f"{sum(s.startswith('HuMs.') for s in adt_names_raw)} HuMs, "
          f"{sum(s.startswith('Isotype_') for s in adt_names_raw)} Isotype), "
          f"{len(rna_names)} RNA genes", flush=True)

    # Densify ADT block (small) and rank.
    Y = a.X[:, adt_idx]
    if sp.issparse(Y):
        Y = np.asarray(Y.todense())
    Y = np.asarray(Y, dtype=np.float64)
    print("[rank] ranking ADTs", flush=True)
    Yr = _rank_columns(Y)
    Yr_centered = Yr - Yr.mean(axis=0, keepdims=True)
    Yr_norm = np.linalg.norm(Yr_centered, axis=0)
    Yr_norm[Yr_norm == 0] = 1.0
    Yz = Yr_centered / Yr_norm

    # Stream over RNA in chunks; rank each chunk; compute Spearman with each ADT.
    n_adts = len(adt_names_raw)
    top_indices = np.full((n_adts, top_n), -1, dtype=np.int64)
    top_corrs = np.full((n_adts, top_n), -np.inf, dtype=np.float64)

    n_rna = len(rna_idx)
    t0 = time.time()
    for start in range(0, n_rna, rna_chunk):
        stop = min(start + rna_chunk, n_rna)
        cols = rna_idx[start:stop]
        block = a.X[:, cols]
        if sp.issparse(block):
            block = np.asarray(block.todense())
        block = np.asarray(block, dtype=np.float64)
        Xr = _rank_columns(block)
        Xr_centered = Xr - Xr.mean(axis=0, keepdims=True)
        Xr_norm = np.linalg.norm(Xr_centered, axis=0)
        Xr_norm[Xr_norm == 0] = 1.0
        Xz = Xr_centered / Xr_norm
        # corrs[i, j] = Yz[:, i] . Xz[:, j]   shape (n_adts, chunk)
        corrs = Yz.T @ Xz
        # For each ADT, find top_n among current top_n ∪ chunk.
        chunk_size = stop - start
        chunk_global_idx = np.arange(start, stop, dtype=np.int64)
        # Combine current top with new chunk
        for i in range(n_adts):
            combined_corr = np.concatenate([top_corrs[i], corrs[i]])
            combined_idx = np.concatenate([top_indices[i], chunk_global_idx])
            order = np.argsort(combined_corr)[::-1][:top_n]
            top_corrs[i] = combined_corr[order]
            top_indices[i] = combined_idx[order]

        elapsed = time.time() - t0
        pct = 100.0 * stop / n_rna
        print(f"[rank] processed {stop}/{n_rna} ({pct:.1f}%) in {elapsed:.1f}s", flush=True)

    # Build output rows.
    rows: List[dict] = []
    for i, adt_raw in enumerate(adt_names_raw):
        for rank in range(top_n):
            ridx = top_indices[i, rank]
            if ridx < 0:
                continue
            rows.append({
                "adt": adt_names_clean[i],
                "adt_raw": adt_raw,
                "rank": rank + 1,
                "gene": rna_names[ridx],
                "spearman": round(float(top_corrs[i, rank]), 4),
            })
    df = pd.DataFrame(rows)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_tsv, sep="\t", index=False)
    print(f"[done] wrote {len(df)} rows -> {output_tsv}", flush=True)

    # Summary table: one row per ADT.
    from .adt_rna_map import load_curated_adt_rna_map
    curated = load_curated_adt_rna_map()
    summary_rows: List[dict] = []
    for i, adt_raw in enumerate(adt_names_raw):
        top_genes = [rna_names[top_indices[i, r]] for r in range(top_n) if top_indices[i, r] >= 0]
        # Curated partners are keyed by the raw "Hu.*"/"HuMs.*"/"Isotype_*" name
        # in the existing adt_rna_map (which only covers Hu.*). For HuMs and
        # Isotype the curated map has no entry — partner_in_top will be FALSE.
        partners = list(curated.get(adt_raw, []))
        # Find the best (lowest) rank at which any curated partner appears.
        partner_ranks = []
        partners_found = []
        for p in partners:
            if p in top_genes:
                partner_ranks.append(top_genes.index(p) + 1)
                partners_found.append(p)
        if partner_ranks:
            best_rank = min(partner_ranks)
            partner_hit = "TRUE"
            partner_rank_str = str(best_rank)
            partners_in_top_str = ",".join(partners_found)
        else:
            partner_hit = "FALSE"
            partner_rank_str = ""
            partners_in_top_str = ""
        summary_rows.append({
            "adt": adt_names_clean[i],
            "curated_partners": ",".join(partners) if partners else "",
            "partner_in_top": partner_hit,
            "partner_best_rank": partner_rank_str,
            "partners_found": partners_in_top_str,
            "top_genes": "|".join(top_genes),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_path = output_tsv.with_name(output_tsv.stem + "_summary.tsv")
    summary_df.to_csv(summary_path, sep="\t", index=False)
    n_true = (summary_df["partner_in_top"] == "TRUE").sum()
    n_with_partners = (summary_df["curated_partners"] != "").sum()
    print(f"[done] wrote summary {len(summary_df)} rows "
          f"({n_true}/{n_with_partners} with curated partner found in top {top_n}) "
          f"-> {summary_path}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--adata", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--top-n", type=int, default=5)
    p.add_argument("--max-cells", type=int, default=30000)
    p.add_argument("--rna-chunk", type=int, default=256)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    compute_top_correlates(
        args.adata,
        output_tsv=args.out,
        top_n=args.top_n,
        max_cells=args.max_cells,
        rna_chunk=args.rna_chunk,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
