"""Build the human lung empirical RNA->ADT feature whitelist.

Faithful port of ``components.rna2adt.build_whitelist`` (bone marrow) and
``components.rna2adt.mouse.build_whitelist`` (mouse). The ranking step and the
three-case selection rule are unchanged:

1. Rank every RNA gene against every ADT by Spearman correlation on a
   30,000-cell random subsample; keep the top 100 per ADT.
2. Curated partner(s) in the top-100 -> use exactly those partners.
3. Otherwise -> curated partner(s) plus the top-``fallback_k`` correlates.
4. No curated partner -> top ``fallback_k + 1`` correlates.
5. A fallback gene selected for more than ``promiscuity_limit`` ADTs is
   dropped and replaced by the next eligible correlate. Curated partners are
   exempt.

One deliberate difference from the bone marrow and mouse builds, requested for
this model: **fallback candidates are restricted to protein-coding genes**.
The top-100 ranking itself still runs over all genes, so the
``curated_in_top100`` decision is identical to the other two models. Pass
``--no-protein-coding-filter`` to reproduce the bone marrow / mouse behaviour
exactly.

The lung RNA and ADT live in two files with two cell-name conventions; see
``components.rna2adt.lung.data`` for the bridge.
"""

from __future__ import annotations

import argparse
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata

from . import data as lung_data
from .adt_rna_map import load_curated_adt_rna_map, strip_prefix


PROTEIN_CODING_DEFAULT = Path(
    "/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/"
    "Hs_Ensembl-annotations_simple-PC.txt"
)
_PROMISCUITY_LIMIT = 10
_FALLBACK_K = 2


def _rank_columns(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr, dtype=np.float64)
    return rankdata(arr, axis=0).astype(np.float64, copy=False)


def _zscore_ranks(arr: np.ndarray) -> np.ndarray:
    ranked = _rank_columns(arr)
    centred = ranked - ranked.mean(axis=0, keepdims=True)
    norms = np.linalg.norm(centred, axis=0)
    norms[norms == 0] = 1.0
    return centred / norms


def top_correlates(
    X,
    Y: np.ndarray,
    *,
    top_n: int = 100,
    chunk: int = 512,
    verbose: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """Spearman rank correlation of every column of ``X`` against every column
    of ``Y``. Returns ``(top_indices, top_corrs)`` of shape ``(n_adts, top_n)``,
    ranked best first. ``X`` may be dense or a CSC/CSR sparse matrix."""
    Yz = _zscore_ranks(np.asarray(Y, dtype=np.float64))
    n_adts = Yz.shape[1]
    n_genes = X.shape[1]
    top_indices = np.full((n_adts, top_n), -1, dtype=np.int64)
    top_corrs = np.full((n_adts, top_n), -np.inf, dtype=np.float64)
    if sp.issparse(X) and not sp.isspmatrix_csc(X):
        X = X.tocsc()
    started = time.time()
    for start in range(0, n_genes, chunk):
        stop = min(start + chunk, n_genes)
        block = X[:, start:stop]
        if sp.issparse(block):
            block = np.asarray(block.todense())
        Xz = _zscore_ranks(np.asarray(block, dtype=np.float64))
        corrs = Yz.T @ Xz
        block_global = np.arange(start, stop, dtype=np.int64)
        for i in range(n_adts):
            merged_corrs = np.concatenate([top_corrs[i], corrs[i]])
            merged_idx = np.concatenate([top_indices[i], block_global])
            merged_corrs = np.where(np.isnan(merged_corrs), -np.inf, merged_corrs)
            order = np.argsort(merged_corrs)[::-1][:top_n]
            top_corrs[i] = merged_corrs[order]
            top_indices[i] = merged_idx[order]
        if verbose and ((stop // chunk) % 10 == 0 or stop == n_genes):
            print(f"[rank] {stop}/{n_genes} ({100 * stop / n_genes:.1f}%) "
                  f"{time.time() - started:.1f}s", flush=True)
    return top_indices, top_corrs


def select_features(
    adt_names_raw: Sequence[str],
    top_indices: np.ndarray,
    rna_names: Sequence[str],
    curated: Dict[str, Tuple[str, ...]],
    *,
    fallback_k: int = _FALLBACK_K,
    promiscuity_limit: int = _PROMISCUITY_LIMIT,
    fallback_allowed: Optional[set] = None,
    verbose: bool = True,
) -> List[dict]:
    """The three-case rule plus the promiscuity dedup, ported verbatim from the
    mouse and bone marrow builders. ``fallback_allowed``, when given, restricts
    which genes may be used as a fallback; curated partners are exempt."""
    top_n = top_indices.shape[1]
    rna_names = list(rna_names)
    fallback_usage: Dict[str, int] = defaultdict(int)
    entries: List[dict] = []
    for i, adt_raw in enumerate(adt_names_raw):
        partners = list(curated.get(adt_raw, []))
        top_genes = [rna_names[top_indices[i, r]] for r in range(top_n) if top_indices[i, r] >= 0]
        partners_in_top = [g for g in partners if g in top_genes]
        eligible = ([g for g in top_genes if fallback_allowed is None or g in fallback_allowed])

        if partners_in_top:
            features = list(partners_in_top)
            source = "curated_in_top100"
        else:
            features = list(partners) if partners else []
            source = "fallback_only" if not partners else "curated+fallback"
            n_fallback = fallback_k if partners else (fallback_k + 1)
            taken = 0
            for gene in eligible:
                if gene in features:
                    continue
                features.append(gene)
                fallback_usage[gene] += 1
                taken += 1
                if taken >= n_fallback:
                    break
        entries.append({
            "adt_clean": strip_prefix(adt_raw),
            "adt_raw": adt_raw,
            "has_curated": bool(partners),
            "curated_partners": partners,
            "partners_in_top100": partners_in_top,
            "features": features,
            "source": source,
            "eligible_genes": eligible,
        })

    over = {g for g, c in fallback_usage.items() if c > promiscuity_limit}
    if verbose:
        print(f"[dedup] {len(over)} genes used as fallback by >{promiscuity_limit} ADTs: "
              f"{sorted(over)[:10]}", flush=True)
    revised_usage: Dict[str, int] = defaultdict(int)
    for entry in entries:
        partners = set(entry["curated_partners"])
        kept: List[str] = [g for g in entry["features"] if g in partners]
        replacements_needed = 0
        for gene in entry["features"]:
            if gene in partners:
                continue
            if gene in over:
                replacements_needed += 1
                continue
            kept.append(gene)
            revised_usage[gene] += 1
        if replacements_needed:
            for gene in entry["eligible_genes"]:
                if replacements_needed == 0:
                    break
                if gene in kept or gene in over:
                    continue
                if revised_usage[gene] >= promiscuity_limit:
                    continue
                kept.append(gene)
                revised_usage[gene] += 1
                replacements_needed -= 1
        entry["features_final"] = kept
    return entries


def entries_to_frame(entries: List[dict]) -> pd.DataFrame:
    return pd.DataFrame([{
        "adt_clean": e["adt_clean"],
        "adt_raw": e["adt_raw"],
        "has_curated_partner": "TRUE" if e["has_curated"] else "FALSE",
        "curated_in_top100": "TRUE" if e["partners_in_top100"] else "FALSE",
        "feature_source": e["source"],
        "n_features": len(e["features_final"]),
        "feature_genes": ",".join(e["features_final"]),
    } for e in entries])


def compute(
    *,
    rna_h5ad: Path,
    adt_txt: Path,
    output_tsv: Path,
    protein_coding_tsv: Optional[Path],
    top_n: int = 100,
    fallback_k: int = _FALLBACK_K,
    promiscuity_limit: int = _PROMISCUITY_LIMIT,
    max_cells: int = 30000,
    rna_chunk: int = 512,
    seed: int = 0,
) -> pd.DataFrame:
    print(f"[load] ADT {adt_txt}", flush=True)
    adt = lung_data.load_adt(adt_txt)
    print(f"[load] {adt.obs_names.size} ADT rows, {adt.var_names.size} ADTs "
          f"({adt.n_unparsed} unparsable row names)", flush=True)
    index = lung_data.align_cells(rna_h5ad, adt)
    report = lung_data.report_alignment(index)
    print(f"[align] matched {report['n_matched_cells']} cells "
          f"(RNA retained {100 * report['rna_retained_fraction']:.2f}%, "
          f"ADT retained {100 * report['adt_retained_fraction']:.2f}%)", flush=True)

    rng = np.random.default_rng(seed)
    n_matched = index.rna_rows.size
    if n_matched > max_cells:
        pick = np.sort(rng.choice(n_matched, max_cells, replace=False))
    else:
        pick = np.arange(n_matched)
    rna_rows = index.rna_rows[pick]
    adt_rows = index.adt_rows[pick]
    print(f"[load] subsampled to {rna_rows.size} cells", flush=True)

    rna_names = lung_data.rna_var_names(rna_h5ad)
    print(f"[load] streaming RNA for {rna_rows.size} cells x {rna_names.size} genes", flush=True)
    X = lung_data.read_rna_rows(rna_h5ad, rna_rows, verbose=False)
    Y = adt.values[adt_rows]
    if X.shape[0] != Y.shape[0]:
        raise ValueError(f"cell count mismatch: RNA {X.shape[0]} vs ADT {Y.shape[0]}")
    print(f"[load] X {X.shape} nnz={X.nnz}, Y {Y.shape}", flush=True)

    fallback_allowed = None
    if protein_coding_tsv is not None:
        symbols = lung_data.load_protein_coding_symbols(protein_coding_tsv)
        fallback_allowed = symbols
        in_atlas = sum(1 for g in rna_names if g in symbols)
        print(f"[pc] {len(symbols)} protein-coding symbols; {in_atlas}/{rna_names.size} "
              f"atlas genes are protein-coding", flush=True)

    top_indices, _ = top_correlates(X, Y, top_n=top_n, chunk=rna_chunk)
    entries = select_features(
        list(adt.var_names), top_indices, rna_names, load_curated_adt_rna_map(),
        fallback_k=fallback_k, promiscuity_limit=promiscuity_limit,
        fallback_allowed=fallback_allowed,
    )
    frame = entries_to_frame(entries)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_tsv, sep="\t", index=False)

    unique_genes = sorted({g for e in entries for g in e["features_final"]})
    print(f"[done] wrote whitelist {len(frame)} rows -> {output_tsv}", flush=True)
    for source in ("curated_in_top100", "curated+fallback", "fallback_only"):
        print(f"[stats] {source}: {int((frame['feature_source'] == source).sum())}", flush=True)
    print(f"[stats] unique feature genes: {len(unique_genes)}", flush=True)
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna-h5ad", type=Path, default=lung_data.RNA_H5AD_DEFAULT)
    parser.add_argument("--adt-txt", type=Path, default=lung_data.ADT_TXT_DEFAULT)
    parser.add_argument("--out", type=Path,
                        default=Path(__file__).parent / "configs" / "empirical_whitelist.tsv")
    parser.add_argument("--protein-coding-tsv", type=Path, default=PROTEIN_CODING_DEFAULT)
    parser.add_argument("--no-protein-coding-filter", action="store_true",
                        help="Reproduce the bone marrow / mouse behaviour exactly")
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--fallback-k", type=int, default=_FALLBACK_K)
    parser.add_argument("--promiscuity-limit", type=int, default=_PROMISCUITY_LIMIT)
    parser.add_argument("--max-cells", type=int, default=30000)
    parser.add_argument("--rna-chunk", type=int, default=512)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    compute(
        rna_h5ad=args.rna_h5ad,
        adt_txt=args.adt_txt,
        output_tsv=args.out,
        protein_coding_tsv=None if args.no_protein_coding_filter else args.protein_coding_tsv,
        top_n=args.top_n,
        fallback_k=args.fallback_k,
        promiscuity_limit=args.promiscuity_limit,
        max_cells=args.max_cells,
        rna_chunk=args.rna_chunk,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
