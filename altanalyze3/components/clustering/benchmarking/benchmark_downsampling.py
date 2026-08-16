#!/usr/bin/env python3
"""Benchmark ICGS3 downsampling alternatives for rare-population retention.

This script intentionally benchmarks only the downsampling stage. It reuses the
official ICGS3 input/QC/RNA filtering and ICGS2 hgvfinder code, then evaluates
alternative cell-selection methods against obs annotations such as HLCA.
"""

from __future__ import annotations

import argparse
import json
import os
import threading
import time
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

if __package__ in {None, ""}:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from altanalyze3.components.clustering.ICGS import (  # noqa: E402
    ICGS3Config,
    _icgs2_community_sampling,
    _icgs2_hgvfinder_adata,
    apply_qc,
    apply_rna_unsupervised_gene_filter,
    pagerank_downsample_adata,
    prepare_expression,
    read_inputs,
)


def _rss_mb() -> float:
    try:
        import psutil

        return float(psutil.Process(os.getpid()).memory_info().rss) / (1024.0**2)
    except Exception:
        import resource

        rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        # Linux reports ru_maxrss in KiB; macOS reports bytes.
        if rss > 10_000_000:
            return rss / (1024.0**2)
        return rss / 1024.0


class MemoryMonitor:
    def __init__(self, interval: float = 0.25):
        self.interval = float(interval)
        self.peak_mb = _rss_mb()
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def _run(self) -> None:
        while not self._stop.is_set():
            self.peak_mb = max(self.peak_mb, _rss_mb())
            time.sleep(self.interval)

    def __enter__(self):
        self.start_mb = _rss_mb()
        self.peak_mb = self.start_mb
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb):
        self._stop.set()
        self._thread.join(timeout=2)
        self.end_mb = _rss_mb()
        self.peak_mb = max(self.peak_mb, self.end_mb)


def _audit_selection(
    full: ad.AnnData,
    selected: Sequence[str],
    obs_columns: Sequence[str],
    *,
    min_cells: int = 20,
    max_cells: int = 1000,
) -> pd.DataFrame:
    selected_set = set(map(str, selected))
    rows = []
    for column in obs_columns:
        if column not in full.obs:
            continue
        labels = full.obs[column].astype(str)
        counts = labels.value_counts(dropna=False)
        rare = counts[(counts > int(min_cells)) & (counts < int(max_cells))]
        for label, total in rare.items():
            cells = set(labels.index[labels == label].astype(str))
            retained = len(cells & selected_set)
            rows.append(
                {
                    "obs_column": column,
                    "annotation": label,
                    "original_cells": int(total),
                    "sampled_cells": int(retained),
                    "retained_fraction": float(retained) / float(total),
                    "sampled_below_20": bool(retained < 20),
                }
            )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["obs_column", "original_cells", "annotation"], kind="mergesort")


def _summarize_audit(method: str, audit: pd.DataFrame, selected_n: int, seconds: float, peak_mb: float) -> Dict[str, object]:
    if audit.empty:
        return {
            "method": method,
            "sampled_cells": int(selected_n),
            "seconds": seconds,
            "peak_rss_mb": peak_mb,
            "rare_populations": 0,
            "rare_failures_below20": 0,
        }
    hlca = audit[audit["obs_column"].eq("HLCA")] if "obs_column" in audit else audit
    eval_df = hlca if not hlca.empty else audit
    ionocyte = eval_df[eval_df["annotation"].eq("Ionocyte")]
    total_rare = float(eval_df["original_cells"].sum())
    total_sampled = float(eval_df["sampled_cells"].sum())
    return {
        "method": method,
        "sampled_cells": int(selected_n),
        "seconds": seconds,
        "peak_rss_mb": peak_mb,
        "rare_populations": int(eval_df.shape[0]),
        "rare_failures_below20": int(eval_df["sampled_below_20"].sum()),
        "rare_cells_original": int(total_rare),
        "rare_cells_sampled": int(total_sampled),
        "rare_cell_capture_fraction": total_sampled / total_rare if total_rare else np.nan,
        "mean_rare_retained_fraction": float(eval_df["retained_fraction"].mean()),
        "median_rare_retained_fraction": float(eval_df["retained_fraction"].median()),
        "min_rare_retained_fraction": float(eval_df["retained_fraction"].min()),
        "ionocyte_original": int(ionocyte["original_cells"].iloc[0]) if not ionocyte.empty else np.nan,
        "ionocyte_sampled": int(ionocyte["sampled_cells"].iloc[0]) if not ionocyte.empty else np.nan,
    }


def _ordered_subset(names: Iterable[str], adata: ad.AnnData, target: int) -> List[str]:
    selected = set(map(str, names))
    out = [str(name) for name in adata.obs_names.astype(str) if str(name) in selected]
    return out[: int(target)]


def method_icgs2_louvain_pagerank(adata: ad.AnnData, config: ICGS3Config, target: int) -> Tuple[List[str], pd.DataFrame]:
    cfg = ICGS3Config(**{**asdict(config), "pagerank_cells": int(target), "pre_pagerank_cells": 0})
    sampled, scores = pagerank_downsample_adata(adata, cfg)
    return sampled.obs_names.astype(str).tolist(), scores


def method_louvain_medoid_direct(adata: ad.AnnData, config: ICGS3Config, target: int) -> Tuple[List[str], pd.DataFrame]:
    hvg = _icgs2_hgvfinder_adata(adata, int(config.downsample_var_genes))
    selected, scores = _icgs2_community_sampling(hvg, int(target), pre_pagerank_cells=int(target))
    return _ordered_subset(selected, adata, target), scores


def method_random(adata: ad.AnnData, config: ICGS3Config, target: int) -> Tuple[List[str], pd.DataFrame]:
    rng = np.random.default_rng(int(config.random_state))
    idx = rng.choice(adata.n_obs, size=min(int(target), adata.n_obs), replace=False)
    idx = np.sort(idx)
    names = adata.obs_names[idx].astype(str).tolist()
    return names, pd.DataFrame({"barcode": names, "selected": True})


def _rarity_scores(hvg: ad.AnnData) -> np.ndarray:
    X = hvg.X.tocsr() if sp.issparse(hvg.X) else sp.csr_matrix(np.asarray(hvg.X))
    presence = X.copy()
    presence.data = np.ones_like(presence.data, dtype=np.float32)
    df = np.asarray(presence.sum(axis=0)).ravel()
    idf = -np.log((df + 1.0) / (presence.shape[0] + 1.0)).astype(np.float32)
    scores = np.asarray(presence @ idf).ravel()
    nnz = np.asarray(presence.sum(axis=1)).ravel()
    return scores / np.sqrt(np.maximum(nnz, 1.0))


def method_rarity_quantile(adata: ad.AnnData, config: ICGS3Config, target: int) -> Tuple[List[str], pd.DataFrame]:
    hvg = _icgs2_hgvfinder_adata(adata, int(config.downsample_var_genes))
    scores = _rarity_scores(hvg)
    return _select_rarity_quantiles(hvg, scores, config, target)


def _select_rarity_quantiles(
    adata: ad.AnnData,
    scores: np.ndarray,
    config: ICGS3Config,
    target: int,
) -> Tuple[List[str], pd.DataFrame]:
    rng = np.random.default_rng(int(config.random_state))
    bins = pd.qcut(pd.Series(scores), q=min(20, max(2, int(target) // 50)), duplicates="drop")
    selected = []
    quota = max(1, int(np.ceil(int(target) / max(1, bins.nunique()))))
    for _, idx in pd.Series(np.arange(hvg.n_obs)).groupby(bins, sort=False):
        values = np.asarray(idx.values, dtype=int)
        if values.size <= quota:
            chosen = values
        else:
            local_scores = scores[values]
            top_n = min(values.size, max(1, quota // 2))
            top = values[np.argsort(local_scores)[::-1][:top_n]]
            remaining = np.setdiff1d(values, top, assume_unique=False)
            fill = rng.choice(remaining, size=min(quota - len(top), len(remaining)), replace=False) if quota > len(top) else []
            chosen = np.concatenate([top, np.asarray(fill, dtype=int)])
        selected.extend(chosen.tolist())
    if len(selected) > int(target):
        selected = sorted(selected, key=lambda i: scores[i], reverse=True)[: int(target)]
    selected = sorted(set(selected))
    names = adata.obs_names[selected].astype(str).tolist()
    return names, pd.DataFrame({"barcode": names, "rarity_score": scores[selected], "selected": True})


def method_all_gene_rarity_quantile(adata: ad.AnnData, config: ICGS3Config, target: int) -> Tuple[List[str], pd.DataFrame]:
    scores = _rarity_scores(adata)
    rng = np.random.default_rng(int(config.random_state))
    bins = pd.qcut(pd.Series(scores), q=min(20, max(2, int(target) // 50)), duplicates="drop")
    selected = []
    quota = max(1, int(np.ceil(int(target) / max(1, bins.nunique()))))
    for _, idx in pd.Series(np.arange(adata.n_obs)).groupby(bins, sort=False):
        values = np.asarray(idx.values, dtype=int)
        if values.size <= quota:
            chosen = values
        else:
            local_scores = scores[values]
            top_n = min(values.size, max(1, quota // 2))
            top = values[np.argsort(local_scores)[::-1][:top_n]]
            remaining = np.setdiff1d(values, top, assume_unique=False)
            fill = rng.choice(remaining, size=min(quota - len(top), len(remaining)), replace=False) if quota > len(top) else []
            chosen = np.concatenate([top, np.asarray(fill, dtype=int)])
        selected.extend(chosen.tolist())
    if len(selected) > int(target):
        selected = sorted(selected, key=lambda i: scores[i], reverse=True)[: int(target)]
    selected = sorted(set(selected))
    names = adata.obs_names[selected].astype(str).tolist()
    return names, pd.DataFrame({"barcode": names, "rarity_score": scores[selected], "selected": True})


def _sparse_gene_coverage_select(
    adata: ad.AnnData,
    config: ICGS3Config,
    target: int,
    *,
    top_per_gene: int,
    core_fraction: float,
) -> Tuple[List[str], pd.DataFrame]:
    """Select cells nominated by top expression of many sparse genes.

    This is a cheap rare-state sketch: each expressed gene nominates its top few
    cells, weighted by inverse detection frequency. It uses the full filtered RNA
    feature universe and adds a small random reserve for common-state coverage.
    """
    X = adata.X.tocsc() if sp.issparse(adata.X) else sp.csc_matrix(np.asarray(adata.X))
    n_cells, n_genes = X.shape
    df = np.diff(X.indptr)
    idf = -np.log((df + 1.0) / (n_cells + 1.0)).astype(np.float32)
    score = np.zeros(n_cells, dtype=np.float32)
    nominations = np.zeros(n_cells, dtype=np.int32)
    min_detected = 3
    max_detected = max(10, int(0.25 * n_cells))
    for gene_idx in range(n_genes):
        start, stop = X.indptr[gene_idx], X.indptr[gene_idx + 1]
        detected = stop - start
        if detected < min_detected or detected > max_detected:
            continue
        cells = X.indices[start:stop]
        values = X.data[start:stop]
        k = min(top_per_gene, detected)
        if detected > k:
            local = np.argpartition(values, -k)[-k:]
        else:
            local = np.arange(detected)
        chosen = cells[local]
        weights = np.asarray(values[local], dtype=np.float32)
        max_w = float(weights.max()) if weights.size else 1.0
        if max_w > 0:
            weights = weights / max_w
        score[chosen] += idf[gene_idx] * weights
        nominations[chosen] += 1
    rng = np.random.default_rng(int(config.random_state))
    core_target = int(round(int(target) * float(core_fraction)))
    core = np.argsort(score)[::-1][: min(core_target, n_cells)]
    selected = set(map(int, core))
    remaining = np.setdiff1d(np.arange(n_cells), np.fromiter(selected, dtype=int), assume_unique=False)
    fill_n = min(int(target) - len(selected), remaining.size)
    if fill_n > 0:
        fill = rng.choice(remaining, size=fill_n, replace=False)
        selected.update(map(int, fill))
    selected_idx = np.array(sorted(selected), dtype=int)
    names = adata.obs_names[selected_idx].astype(str).tolist()
    return names, pd.DataFrame(
        {
            "barcode": names,
            "gene_coverage_score": score[selected_idx],
            "gene_nominations": nominations[selected_idx],
            "selected": True,
        }
    )


def method_sparse_gene_coverage(adata: ad.AnnData, config: ICGS3Config, target: int) -> Tuple[List[str], pd.DataFrame]:
    return _sparse_gene_coverage_select(adata, config, target, top_per_gene=3, core_fraction=0.8)


def method_sparse_gene_coverage_aggressive(adata: ad.AnnData, config: ICGS3Config, target: int) -> Tuple[List[str], pd.DataFrame]:
    return _sparse_gene_coverage_select(adata, config, target, top_per_gene=10, core_fraction=0.9)


def method_random_projection_kmeans(adata: ad.AnnData, config: ICGS3Config, target: int) -> Tuple[List[str], pd.DataFrame]:
    from sklearn.cluster import MiniBatchKMeans
    from sklearn.random_projection import SparseRandomProjection

    hvg = _icgs2_hgvfinder_adata(adata, int(config.downsample_var_genes))
    X = hvg.X.tocsr() if sp.issparse(hvg.X) else sp.csr_matrix(np.asarray(hvg.X))
    n_components = min(32, max(2, X.shape[1] - 1))
    projector = SparseRandomProjection(n_components=n_components, random_state=int(config.random_state))
    Z = projector.fit_transform(X)
    Z = Z.astype(np.float32).toarray() if sp.issparse(Z) else np.asarray(Z, dtype=np.float32)
    n_clusters = min(int(target), Z.shape[0])
    model = MiniBatchKMeans(
        n_clusters=n_clusters,
        random_state=int(config.random_state),
        batch_size=8192,
        n_init=1,
        max_iter=40,
        reassignment_ratio=0.01,
    )
    labels = model.fit_predict(Z)
    centers = model.cluster_centers_.astype(np.float32, copy=False)
    dist = np.einsum("ij,ij->i", Z - centers[labels], Z - centers[labels])
    selected = []
    for cluster in range(n_clusters):
        idx = np.where(labels == cluster)[0]
        if idx.size:
            selected.append(int(idx[np.argmin(dist[idx])]))
    names = hvg.obs_names[sorted(selected)].astype(str).tolist()
    return names, pd.DataFrame({"barcode": names, "selected": True})


METHODS = {
    "icgs2_louvain_pagerank": method_icgs2_louvain_pagerank,
    "louvain_medoid_direct": method_louvain_medoid_direct,
    "rarity_quantile": method_rarity_quantile,
    "all_gene_rarity_quantile": method_all_gene_rarity_quantile,
    "sparse_gene_coverage": method_sparse_gene_coverage,
    "sparse_gene_coverage_aggressive": method_sparse_gene_coverage_aggressive,
    "random_projection_kmeans": method_random_projection_kmeans,
    "random": method_random,
}


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark ICGS3 downsampling alternatives.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--methods", default=",".join(METHODS))
    parser.add_argument("--target-cells", type=int, default=5000)
    parser.add_argument("--species", default="Hs")
    parser.add_argument("--modality", default="rna")
    parser.add_argument("--input-normalized", action="store_true")
    parser.add_argument("--min-genes", type=int, default=500)
    parser.add_argument("--min-cells", type=int, default=5)
    parser.add_argument("--min-counts", type=int, default=1000)
    parser.add_argument("--mito-percent", type=float, default=30.0)
    parser.add_argument("--downsample-var-genes", type=int, default=500)
    parser.add_argument("--audit-obs", default="HLCA")
    parser.add_argument("--random-state", type=int, default=0)
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "methods").mkdir(exist_ok=True)
    config = ICGS3Config(
        input_paths=[args.input],
        output_dir=str(outdir),
        modality=args.modality,
        species=args.species,
        input_normalized=bool(args.input_normalized),
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        min_counts=args.min_counts,
        mito_percent=args.mito_percent,
        pagerank_cells=args.target_cells,
        downsample_var_genes=args.downsample_var_genes,
        random_state=args.random_state,
        retention_audit_obs=args.audit_obs,
        generate_umap=False,
        write_h5ad=False,
    )
    (outdir / "benchmark_config.json").write_text(json.dumps(asdict(config), indent=2), encoding="utf-8")

    t0 = time.time()
    adata = read_inputs([args.input])
    adata = apply_qc(
        adata,
        min_genes=config.min_genes,
        min_cells=config.min_cells,
        min_counts=config.min_counts,
        mito_percent=config.mito_percent,
        layer=config.layer,
    )
    adata = prepare_expression(adata, config)
    adata = apply_rna_unsupervised_gene_filter(adata, config, outdir=None)
    preprocess_seconds = time.time() - t0

    obs_columns = [c.strip() for c in args.audit_obs.split(",") if c.strip()]
    selected_methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    summaries = []
    for method in selected_methods:
        if method not in METHODS:
            raise KeyError(f"Unknown method: {method}. Available: {', '.join(METHODS)}")
        print(f"[benchmark] starting {method}", flush=True)
        start = time.time()
        with MemoryMonitor() as monitor:
            selected, scores = METHODS[method](adata, config, int(args.target_cells))
        seconds = time.time() - start
        audit = _audit_selection(adata, selected, obs_columns)
        method_dir = outdir / "methods" / method
        method_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"barcode": selected}).to_csv(method_dir / "selected_barcodes.tsv", sep="\t", index=False)
        scores.to_csv(method_dir / "downsampling_scores.tsv", sep="\t", index=False)
        audit.to_csv(method_dir / "rare_population_retention.tsv", sep="\t", index=False)
        summary = _summarize_audit(method, audit, len(selected), seconds, monitor.peak_mb)
        summary["preprocess_seconds"] = preprocess_seconds
        summaries.append(summary)
        print(f"[benchmark] finished {method}: {json.dumps(summary, default=str)}", flush=True)
        pd.DataFrame(summaries).to_csv(outdir / "downsampling_benchmark_summary.tsv", sep="\t", index=False)

    pd.DataFrame(summaries).to_csv(outdir / "downsampling_benchmark_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
