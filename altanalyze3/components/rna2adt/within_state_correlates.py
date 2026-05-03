"""Within-cell-state RNA correlates of each ADT.

For each ADT we (a) find the cell state with the highest mean ADT expression
in the atlas, (b) subset to cells in that state, (c) compute Spearman
correlation between the ADT and every RNA gene using only those cells.

The intent is to remove the dominant lineage-identity signal that drives
the global correlation analysis: a B-lineage marker like ``CD79A`` is the
top global correlate of ``Hu.CD19`` only because both are restricted to
B cells. Within a single B-cell state, however, ``CD79A`` is everywhere
and the correlation collapses, leaving room for genes that **co-vary
cell-by-cell** with CD19 protein -- those are far more likely to share
transcriptional regulators with CD19 itself, and to remain predictive in
cancer cells where lineage identity has been corrupted but the regulon
is preserved.

Top state selection uses mean ADT across cells in each ``cluster_col``
cluster; a minimum of ``min_cells`` cells per state is required, and we
fall back to the next-best state if needed.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import List, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata

from .adt_rna_map import load_curated_adt_rna_map


_ADT_PREFIXES = ("Hu.", "HuMs.", "Isotype_")


def _strip_prefix(name: str) -> str:
    for prefix in _ADT_PREFIXES:
        if name.startswith(prefix):
            return name[len(prefix):]
    return name


def _rank_columns(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr, dtype=np.float64)
    return rankdata(arr, axis=0).astype(np.float64, copy=False)


def _pick_top_state(
    adt_vec: np.ndarray,
    clusters: np.ndarray,
    *,
    min_cells: int,
) -> Tuple[str, np.ndarray, float]:
    """Return (cluster_label, cell_mask, mean_adt) for the highest-mean cluster
    that has at least min_cells cells. If no cluster qualifies, return the
    overall top-mean cluster regardless of size.
    """
    states, inverse = np.unique(clusters, return_inverse=True)
    means = np.zeros(len(states), dtype=np.float64)
    counts = np.zeros(len(states), dtype=np.int64)
    for k in range(len(states)):
        m = inverse == k
        counts[k] = int(m.sum())
        if counts[k] > 0:
            means[k] = float(adt_vec[m].mean())

    # Ranking: prefer states meeting min_cells, sort by mean descending.
    order = np.argsort(-means)
    for k in order:
        if counts[k] >= min_cells:
            mask = inverse == k
            return str(states[k]), mask, float(means[k])
    # Fallback: top mean regardless of size.
    k = int(order[0])
    return str(states[k]), inverse == k, float(means[k])


def compute(
    adata_path: Path,
    *,
    output_tsv: Path,
    cluster_col: str = "Level 3 Multimodal",
    top_n: int = 50,
    min_cells: int = 100,
    max_cells_global: int = 50000,
    rna_chunk: int = 512,
    seed: int = 0,
) -> None:
    print(f"[load] {adata_path}", flush=True)
    a = ad.read_h5ad(adata_path)
    print(f"[load] {a.n_obs} cells, {a.n_vars} vars", flush=True)
    if cluster_col not in a.obs.columns:
        raise KeyError(f"obs[{cluster_col!r}] not present")

    rng = np.random.default_rng(seed)
    if a.n_obs > max_cells_global:
        idx = np.sort(rng.choice(a.n_obs, max_cells_global, replace=False))
        a = a[idx].to_memory()
        print(f"[load] subsampled to {a.n_obs} cells", flush=True)

    var_names = np.array([str(v) for v in a.var_names])
    is_adt = np.array([s.startswith(_ADT_PREFIXES) for s in var_names])
    adt_idx = np.where(is_adt)[0]
    rna_idx = np.where(~is_adt)[0]
    adt_names_raw = var_names[is_adt]
    adt_names_clean = np.array([_strip_prefix(n) for n in adt_names_raw])
    rna_names = var_names[~is_adt]
    print(f"[load] {len(adt_names_raw)} ADTs ({sum(s.startswith('Hu.') for s in adt_names_raw)} Hu, "
          f"{sum(s.startswith('HuMs.') for s in adt_names_raw)} HuMs, "
          f"{sum(s.startswith('Isotype_') for s in adt_names_raw)} Isotype), "
          f"{len(rna_names)} RNA genes", flush=True)

    clusters = a.obs[cluster_col].astype(str).to_numpy()

    # Densify ADT block once.
    Y = a.X[:, adt_idx]
    if sp.issparse(Y):
        Y = np.asarray(Y.todense())
    Y = np.asarray(Y, dtype=np.float64)

    # For each ADT pick the top cell state.
    print("[state] picking top cell state per ADT", flush=True)
    state_info: List[Tuple[str, np.ndarray, float, int]] = []
    for i in range(len(adt_names_raw)):
        cluster_label, mask, mean_adt = _pick_top_state(Y[:, i], clusters, min_cells=min_cells)
        state_info.append((cluster_label, mask, mean_adt, int(mask.sum())))

    curated = load_curated_adt_rna_map()

    # Gather the unique cell-state masks so we can batch correlation work.
    # Each unique mask -> list of ADTs that picked it.
    mask_to_adts: dict[bytes, list[int]] = {}
    mask_storage: dict[bytes, np.ndarray] = {}
    for i, (label, mask, _, _) in enumerate(state_info):
        key = mask.tobytes()
        if key not in mask_to_adts:
            mask_to_adts[key] = []
            mask_storage[key] = mask
        mask_to_adts[key].append(i)
    print(f"[state] {len(mask_to_adts)} unique cell-state masks across {len(adt_names_raw)} ADTs",
          flush=True)

    rows: List[dict] = []
    summary_rows: List[dict] = []
    n_groups = len(mask_to_adts)
    t0 = time.time()
    group_n = 0

    # Cache: for each unique mask compute Yz_within once, plus the X chunks.
    for key, adt_indices in mask_to_adts.items():
        group_n += 1
        mask = mask_storage[key]
        n_cells_in_state = int(mask.sum())
        if n_cells_in_state < 2:
            print(f"[skip] state has {n_cells_in_state} cells", flush=True)
            continue
        # Slice the ADT block for these cells, rank within state.
        Y_state = Y[mask][:, adt_indices]
        Yr = _rank_columns(Y_state)
        Yc = Yr - Yr.mean(axis=0, keepdims=True)
        Yn = np.linalg.norm(Yc, axis=0)
        Yn[Yn == 0] = 1.0
        Yz = Yc / Yn  # (n_cells_in_state, len(adt_indices))

        # Track top-N per ADT in this group.
        n_local = len(adt_indices)
        top_corrs = np.full((n_local, top_n), -np.inf, dtype=np.float64)
        top_global_idx = np.full((n_local, top_n), -1, dtype=np.int64)

        # Stream over RNA chunks using cells in mask.
        cell_indices = np.where(mask)[0]
        Xfull = a.X
        for start in range(0, len(rna_idx), rna_chunk):
            stop = min(start + rna_chunk, len(rna_idx))
            cols = rna_idx[start:stop]
            block = Xfull[cell_indices][:, cols]
            if sp.issparse(block):
                block = np.asarray(block.todense())
            block = np.asarray(block, dtype=np.float64)
            Xr = _rank_columns(block)
            Xc = Xr - Xr.mean(axis=0, keepdims=True)
            Xn = np.linalg.norm(Xc, axis=0)
            Xn[Xn == 0] = 1.0
            Xz = Xc / Xn
            corrs = Yz.T @ Xz                # (n_local, chunk_size)
            chunk_global = np.arange(start, stop, dtype=np.int64)
            for li in range(n_local):
                combined_corr = np.concatenate([top_corrs[li], corrs[li]])
                combined_idx = np.concatenate([top_global_idx[li], chunk_global])
                # Fix NaN safely
                combined_corr = np.where(np.isnan(combined_corr), -np.inf, combined_corr)
                order = np.argsort(combined_corr)[::-1][:top_n]
                top_corrs[li] = combined_corr[order]
                top_global_idx[li] = combined_idx[order]

        elapsed = time.time() - t0
        # Output
        for li, i_global in enumerate(adt_indices):
            adt_raw = adt_names_raw[i_global]
            adt_clean = adt_names_clean[i_global]
            cluster_label, _, mean_adt, n_in_state = state_info[i_global]
            top_genes = []
            for r in range(top_n):
                gi = top_global_idx[li, r]
                if gi < 0:
                    continue
                gene = rna_names[gi]
                top_genes.append(gene)
                rows.append({
                    "adt": adt_clean,
                    "cluster": cluster_label,
                    "n_cells_in_state": n_in_state,
                    "mean_adt_in_state": round(mean_adt, 4),
                    "rank": r + 1,
                    "gene": gene,
                    "spearman": round(float(top_corrs[li, r]), 4),
                })
            partners = list(curated.get(adt_raw, []))
            partner_ranks = []
            partners_found = []
            for p in partners:
                if p in top_genes:
                    partner_ranks.append(top_genes.index(p) + 1)
                    partners_found.append(p)
            if partner_ranks:
                summary_rows.append({
                    "adt": adt_clean,
                    "cluster": cluster_label,
                    "n_cells_in_state": n_in_state,
                    "mean_adt_in_state": round(mean_adt, 4),
                    "curated_partners": ",".join(partners) if partners else "",
                    "partner_in_top": "TRUE",
                    "partner_best_rank": str(min(partner_ranks)),
                    "partners_found": ",".join(partners_found),
                    "top_genes": "|".join(top_genes),
                })
            else:
                summary_rows.append({
                    "adt": adt_clean,
                    "cluster": cluster_label,
                    "n_cells_in_state": n_in_state,
                    "mean_adt_in_state": round(mean_adt, 4),
                    "curated_partners": ",".join(partners) if partners else "",
                    "partner_in_top": "FALSE",
                    "partner_best_rank": "",
                    "partners_found": "",
                    "top_genes": "|".join(top_genes),
                })
        print(f"[group {group_n}/{n_groups}] cluster={cluster_label!r:30s} "
              f"n_cells={n_in_state:6d} adts={n_local:3d} elapsed={elapsed:.1f}s", flush=True)

    df = pd.DataFrame(rows)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_tsv, sep="\t", index=False)
    print(f"[done] wrote {len(df)} long-format rows -> {output_tsv}", flush=True)

    summary_df = pd.DataFrame(summary_rows)
    summary_path = output_tsv.with_name(output_tsv.stem + "_summary.tsv")
    summary_df.to_csv(summary_path, sep="\t", index=False)
    n_with_partners = (summary_df["curated_partners"] != "").sum()
    n_true = (summary_df["partner_in_top"] == "TRUE").sum()
    print(f"[done] wrote summary {len(summary_df)} rows "
          f"({n_true}/{n_with_partners} with curated partner found in top {top_n}) "
          f"-> {summary_path}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--adata", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--cluster-col", default="Level 3 Multimodal")
    p.add_argument("--top-n", type=int, default=50)
    p.add_argument("--min-cells", type=int, default=100)
    p.add_argument("--max-cells-global", type=int, default=50000)
    p.add_argument("--rna-chunk", type=int, default=512)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    compute(
        args.adata,
        output_tsv=args.out,
        cluster_col=args.cluster_col,
        top_n=args.top_n,
        min_cells=args.min_cells,
        max_cells_global=args.max_cells_global,
        rna_chunk=args.rna_chunk,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
