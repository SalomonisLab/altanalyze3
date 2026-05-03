"""Build the mouse empirical RNA->ADT feature whitelist.

Mirrors ``components.rna2adt.build_whitelist`` for the murine panel.
The ADT prefix is ``ADT-`` (vs human ``Hu.``); isotype controls
(``ADT-Mouse_IgG*_isotype_Ctrl``, ``ADT-Rat_IgG*_Isotype_Ctrl``,
``ADT-Armenian_Hamster_IgG_Isotype_Ctrl``) and any ``Hu.IgG.Fc``-style
control are dropped from the output set so they don't pollute the
panel-union with non-informative features.
"""

from __future__ import annotations

import argparse
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata

from .adt_mgi_map import load_curated_adt_mgi_map, strip_prefix


_ADT_PREFIX = "ADT-"
_DROP_TOKENS = ("Isotype_Ctrl", "isotype_Ctrl", "isotype_ctrl")
_PROMISCUITY_LIMIT = 10
_FALLBACK_K = 2


def _is_isotype(name: str) -> bool:
    return any(tok in name for tok in _DROP_TOKENS)


def _rank_columns(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr, dtype=np.float64)
    return rankdata(arr, axis=0).astype(np.float64, copy=False)


def compute(
    adata_path: Path,
    *,
    output_tsv: Path,
    top_n: int = 100,
    fallback_k: int = _FALLBACK_K,
    promiscuity_limit: int = _PROMISCUITY_LIMIT,
    max_cells: int = 30000,
    rna_chunk: int = 512,
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
    is_adt = np.array([s.startswith(_ADT_PREFIX) for s in var_names])
    is_isotype = np.array([_is_isotype(s) for s in var_names])
    is_adt_kept = is_adt & ~is_isotype
    n_dropped = int((is_adt & is_isotype).sum())
    print(f"[load] dropping {n_dropped} isotype/IgG control ADTs from output set", flush=True)
    adt_idx = np.where(is_adt_kept)[0]
    rna_idx = np.where(~is_adt)[0]
    adt_names_raw = var_names[is_adt_kept]
    adt_names_clean = np.array([strip_prefix(n) for n in adt_names_raw])
    rna_names = var_names[~is_adt]
    print(f"[load] {len(adt_names_raw)} ADTs (kept), {len(rna_names)} RNA genes", flush=True)

    Y = a.X[:, adt_idx]
    if sp.issparse(Y):
        Y = np.asarray(Y.todense())
    Y = np.asarray(Y, dtype=np.float64)
    print("[rank] ranking ADTs", flush=True)
    Yr = _rank_columns(Y)
    Yc = Yr - Yr.mean(axis=0, keepdims=True)
    Yn = np.linalg.norm(Yc, axis=0)
    Yn[Yn == 0] = 1.0
    Yz = Yc / Yn

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
        Xc = Xr - Xr.mean(axis=0, keepdims=True)
        Xn = np.linalg.norm(Xc, axis=0)
        Xn[Xn == 0] = 1.0
        Xz = Xc / Xn
        corrs = Yz.T @ Xz
        chunk_global_idx = np.arange(start, stop, dtype=np.int64)
        for i in range(n_adts):
            cc = np.concatenate([top_corrs[i], corrs[i]])
            ci = np.concatenate([top_indices[i], chunk_global_idx])
            cc = np.where(np.isnan(cc), -np.inf, cc)
            order = np.argsort(cc)[::-1][:top_n]
            top_corrs[i] = cc[order]
            top_indices[i] = ci[order]
        elapsed = time.time() - t0
        if (stop // rna_chunk) % 10 == 0 or stop == n_rna:
            print(f"[rank] {stop}/{n_rna} ({100*stop/n_rna:.1f}%) {elapsed:.1f}s", flush=True)

    curated = load_curated_adt_mgi_map()

    # Pass 1: assign curated_in_top100 / curated+fallback / fallback_only
    fallback_usage: Dict[str, int] = defaultdict(int)
    pass1: List[dict] = []
    for i in range(n_adts):
        adt_raw = adt_names_raw[i]
        adt_clean = adt_names_clean[i]
        partners = list(curated.get(adt_raw, []))
        top100_genes = [rna_names[top_indices[i, r]] for r in range(top_n) if top_indices[i, r] >= 0]
        partners_in_top = [g for g in partners if g in top100_genes]
        has_curated = bool(partners)

        if partners_in_top:
            features = list(partners_in_top)
            source = "curated_in_top100"
        else:
            features = list(partners) if partners else []
            source = "fallback_only" if not partners else "curated+fallback"
            n_fallback = fallback_k if partners else (fallback_k + 1)
            n_taken = 0
            for g in top100_genes:
                if g in features:
                    continue
                features.append(g)
                fallback_usage[g] += 1
                n_taken += 1
                if n_taken >= n_fallback:
                    break
        pass1.append({
            "adt_clean": adt_clean,
            "adt_raw": adt_raw,
            "has_curated": has_curated,
            "curated_partners": partners,
            "partners_in_top100": partners_in_top,
            "features": features,
            "source": source,
            "top100_genes": top100_genes,
        })

    # Pass 2: enforce promiscuity
    over = {g for g, c in fallback_usage.items() if c > promiscuity_limit}
    print(f"[dedup] {len(over)} genes are used as fallback by >{promiscuity_limit} ADTs: {sorted(over)[:10]}", flush=True)
    revised_fallback_usage: Dict[str, int] = defaultdict(int)
    for entry in pass1:
        partners = set(entry["curated_partners"])
        kept: List[str] = []
        for g in entry["features"]:
            if g in partners:
                kept.append(g)
        replacements_needed = 0
        for g in entry["features"]:
            if g in partners:
                continue
            if g in over:
                replacements_needed += 1
                continue
            kept.append(g)
            revised_fallback_usage[g] += 1
        if replacements_needed > 0:
            for g in entry["top100_genes"]:
                if replacements_needed == 0:
                    break
                if g in kept or g in over:
                    continue
                if revised_fallback_usage[g] >= promiscuity_limit:
                    continue
                kept.append(g)
                revised_fallback_usage[g] += 1
                replacements_needed -= 1
        entry["features_final"] = kept

    rows = []
    for entry in pass1:
        rows.append({
            "adt_clean": entry["adt_clean"],
            "adt_raw": entry["adt_raw"],
            "has_curated_partner": "TRUE" if entry["has_curated"] else "FALSE",
            "curated_in_top100": "TRUE" if entry["partners_in_top100"] else "FALSE",
            "feature_source": entry["source"],
            "n_features": len(entry["features_final"]),
            "feature_genes": ",".join(entry["features_final"]),
        })
    df = pd.DataFrame(rows)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_tsv, sep="\t", index=False)
    n_curated_only = (df["feature_source"] == "curated_in_top100").sum()
    n_curated_plus = (df["feature_source"] == "curated+fallback").sum()
    n_fallback_only = (df["feature_source"] == "fallback_only").sum()
    unique_genes = sorted({g for entry in pass1 for g in entry["features_final"]})
    print(f"[done] wrote whitelist {len(df)} rows -> {output_tsv}", flush=True)
    print(f"[stats] curated_in_top100: {n_curated_only}", flush=True)
    print(f"[stats] curated+fallback:  {n_curated_plus}", flush=True)
    print(f"[stats] fallback_only:     {n_fallback_only}", flush=True)
    print(f"[stats] unique feature genes: {len(unique_genes)}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--adata", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--top-n", type=int, default=100)
    p.add_argument("--fallback-k", type=int, default=2)
    p.add_argument("--promiscuity-limit", type=int, default=10)
    p.add_argument("--max-cells", type=int, default=30000)
    p.add_argument("--rna-chunk", type=int, default=512)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    compute(
        args.adata,
        output_tsv=args.out,
        top_n=args.top_n,
        fallback_k=args.fallback_k,
        promiscuity_limit=args.promiscuity_limit,
        max_cells=args.max_cells,
        rna_chunk=args.rna_chunk,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
