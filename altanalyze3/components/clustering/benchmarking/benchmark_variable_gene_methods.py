#!/usr/bin/env python3
"""Benchmark variable-gene methods against MarkerFinder population markers.

For a labeled reference such as Adams HLCA, this script:

1. Applies ICGS3 QC and RNA unsupervised gene filtering.
2. Samples up to N cells per population.
3. Runs AltAnalyze3 sparse MarkerFinder on the sampled labeled populations.
4. Stores top MarkerFinder genes per population as a lookup database.
5. Tests variable-gene methods capped at <=3,000 genes.
6. Reports whether each method captures at least K MarkerFinder genes per state.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD

if __package__ in {None, ""}:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from altanalyze3.components.cellHarmony.markerFinder import find_markers_from_adata  # noqa: E402
from altanalyze3.components.clustering.ICGS import (  # noqa: E402
    ICGS3Config,
    _icgs2_hgvfinder_adata,
    apply_qc,
    apply_rna_unsupervised_gene_filter,
    prepare_expression,
    read_inputs,
)
from altanalyze3.components.clustering.benchmark_downsampling import _rss_mb  # noqa: E402


def _sample_cells_per_state(adata: ad.AnnData, label_key: str, cells_per_state: int, seed: int) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    labels = adata.obs[label_key].astype(str)
    selected = []
    for label, idx in labels.groupby(labels).groups.items():
        idx = np.asarray(list(idx), dtype=object)
        if idx.size <= int(cells_per_state):
            chosen = idx
        else:
            chosen = rng.choice(idx, size=int(cells_per_state), replace=False)
        selected.extend([str(x) for x in chosen])
    selected = [name for name in adata.obs_names.astype(str) if name in set(selected)]
    return adata[selected].copy()


def _matrix_stats(adata: ad.AnnData) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    X = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(np.asarray(adata.X))
    means = np.asarray(X.mean(axis=0)).ravel()
    variances = np.asarray(X.multiply(X).mean(axis=0)).ravel() - means**2
    detected = np.asarray((X > 0).sum(axis=0)).ravel()
    return means, variances, detected


def genes_top_dispersion(adata: ad.AnnData, n: int) -> pd.Index:
    hvg = _icgs2_hgvfinder_adata(adata, int(n))
    return pd.Index(hvg.var_names.astype(str))


def genes_top_variance(adata: ad.AnnData, n: int) -> pd.Index:
    _, variances, _ = _matrix_stats(adata)
    idx = np.argsort(variances)[::-1][: int(n)]
    return pd.Index(adata.var_names[idx].astype(str))


def genes_detection_balanced_dispersion(adata: ad.AnnData, n: int) -> pd.Index:
    means, variances, detected = _matrix_stats(adata)
    freq = detected / float(adata.n_obs)
    with np.errstate(divide="ignore", invalid="ignore"):
        dispersion = np.divide(variances, means, out=np.zeros_like(variances, dtype=float), where=means > 0)
    rare_weight = -np.log(np.clip(freq, 1.0 / max(adata.n_obs, 1), 0.95))
    score = dispersion * rare_weight
    usable = np.where((detected >= 3) & np.isfinite(score))[0]
    order = usable[np.argsort(score[usable])[::-1]]
    return pd.Index(adata.var_names[order[: int(n)]].astype(str))


def genes_sparse_gene_coverage_scores(adata: ad.AnnData, n: int, top_per_gene: int = 10) -> pd.Index:
    X = adata.X.tocsc() if sp.issparse(adata.X) else sp.csc_matrix(np.asarray(adata.X))
    n_cells, n_genes = X.shape
    detected = np.diff(X.indptr)
    idf = -np.log((detected + 1.0) / (n_cells + 1.0))
    gene_score = np.zeros(n_genes, dtype=np.float64)
    min_detected = 3
    max_detected = max(10, int(0.25 * n_cells))
    for gene_idx in range(n_genes):
        start, stop = X.indptr[gene_idx], X.indptr[gene_idx + 1]
        d = stop - start
        if d < min_detected or d > max_detected:
            continue
        values = X.data[start:stop]
        k = min(int(top_per_gene), d)
        top = np.argpartition(values, -k)[-k:] if d > k else np.arange(d)
        gene_score[gene_idx] = float(idf[gene_idx]) * float(np.mean(values[top])) * np.sqrt(k)
    order = np.argsort(gene_score)[::-1]
    order = order[gene_score[order] > 0]
    return pd.Index(adata.var_names[order[: int(n)]].astype(str))


def genes_pca_loadings(adata: ad.AnnData, n: int, n_components: int = 50, elbow_fraction: float = 1.0) -> pd.Index:
    X = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(np.asarray(adata.X))
    max_components = max(2, min(int(n_components), adata.n_obs - 1, adata.n_vars - 1))
    svd = TruncatedSVD(n_components=max_components, random_state=0)
    svd.fit(X)
    ratios = svd.explained_variance_ratio_
    # Elbow by maximum distance from first-last line, then threshold multiplier.
    x = np.arange(ratios.size, dtype=float)
    start = np.array([x[0], ratios[0]], dtype=float)
    end = np.array([x[-1], ratios[-1]], dtype=float)
    line = end - start
    denom = np.linalg.norm(line)
    if denom > 0:
        points = np.column_stack([x, ratios])
        elbow = int(np.argmax(np.abs(np.cross(line, start - points)) / denom) + 1)
    else:
        elbow = max_components
    pcs = max(2, min(max_components, int(np.ceil(elbow * float(elbow_fraction)))))
    loadings = np.abs(svd.components_[:pcs, :])
    score = np.max(loadings, axis=0)
    order = np.argsort(score)[::-1][: int(n)]
    return pd.Index(adata.var_names[order].astype(str))


def _marker_database(markers: pd.DataFrame, top_n: int, rho: float) -> pd.DataFrame:
    out = markers.copy()
    out = out[out["direction"].eq("up")]
    out = out[out["pearson_r"] > float(rho)]
    out = out.sort_values(["cluster", "pearson_r", "marker"], ascending=[True, False, True])
    out = out.groupby("cluster", sort=False).head(int(top_n)).copy()
    out["marker_rank"] = out.groupby("cluster")["pearson_r"].rank(method="first", ascending=False).astype(int)
    return out.rename(columns={"cluster": "cell_state"})


def _score_method(method: str, genes: Sequence[str], marker_db: pd.DataFrame, min_hits: int) -> pd.DataFrame:
    genes = set(map(str, genes))
    rows = []
    for state, sub in marker_db.groupby("cell_state", sort=False):
        markers = set(sub["marker"].astype(str))
        hits = sorted(markers & genes)
        rows.append(
            {
                "method": method,
                "cell_state": state,
                "marker_genes_available": len(markers),
                "marker_genes_captured": len(hits),
                "capture_fraction": len(hits) / float(len(markers)) if markers else np.nan,
                "passes_min_10": bool(len(hits) >= int(min_hits)),
                "captured_genes": ",".join(hits),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark variable-gene methods against sparse MarkerFinder state markers.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--label-key", default="HLCA")
    parser.add_argument("--cells-per-state", type=int, default=200)
    parser.add_argument("--marker-top-n", type=int, default=50)
    parser.add_argument("--marker-rho", type=float, default=0.4)
    parser.add_argument("--max-genes", type=int, default=3000)
    parser.add_argument("--min-hits", type=int, default=10)
    parser.add_argument("--input-normalized", action="store_true")
    parser.add_argument("--species", default="Hs")
    parser.add_argument("--random-state", type=int, default=0)
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    config = ICGS3Config(
        input_paths=[args.input],
        output_dir=str(outdir),
        modality="rna",
        species=args.species,
        input_normalized=bool(args.input_normalized),
        min_genes=500,
        min_cells=5,
        min_counts=1000,
        mito_percent=30.0,
        generate_umap=False,
        write_h5ad=False,
    )
    (outdir / "config.json").write_text(json.dumps({**asdict(config), **vars(args)}, indent=2), encoding="utf-8")

    t0 = time.time()
    adata = read_inputs([args.input])
    adata = apply_qc(adata, min_genes=config.min_genes, min_cells=config.min_cells, min_counts=config.min_counts, mito_percent=config.mito_percent, layer=config.layer)
    adata = prepare_expression(adata, config)
    adata = apply_rna_unsupervised_gene_filter(adata, config)
    if args.label_key not in adata.obs:
        raise KeyError(f"Missing label key: {args.label_key}")
    adata = adata[adata.obs[args.label_key].astype(str).ne("")].copy()
    preprocess_seconds = time.time() - t0

    sampled = _sample_cells_per_state(adata, args.label_key, args.cells_per_state, args.random_state)
    sampled.obs["MarkerFinder_state"] = sampled.obs[args.label_key].astype(str)
    pd.DataFrame(
        {
            "cell_state": sampled.obs["MarkerFinder_state"].astype(str).value_counts().index,
            "sampled_cells_for_markerfinder": sampled.obs["MarkerFinder_state"].astype(str).value_counts().values,
        }
    ).to_csv(outdir / "markerfinder_cells_per_state.tsv", sep="\t", index=False)

    marker_dir = outdir / "MarkerFinder"
    t_marker = time.time()
    marker_outputs = find_markers_from_adata(
        sampled,
        "MarkerFinder_state",
        output_dir=str(marker_dir),
        n_markers=args.marker_top_n,
        direction="up",
        rho_threshold=args.marker_rho,
        min_markers_per_cluster=1,
        write_outputs=True,
        heatmap_filename="adams_HLCA_markerfinder_heatmap.pdf",
        marker_table_filename="adams_HLCA_markerfinder_markers.tsv",
        heatmap_table_filename="adams_HLCA_markerfinder_heatmap.tsv",
    )
    marker_seconds = time.time() - t_marker
    marker_db = _marker_database(marker_outputs.markers, args.marker_top_n, args.marker_rho)
    marker_db.to_csv(outdir / "adams_HLCA_markerfinder_top50_gt0.4_marker_database.tsv", sep="\t", index=False)

    methods = {
        "dispersion_500": lambda: genes_top_dispersion(adata, 500),
        "dispersion_1000": lambda: genes_top_dispersion(adata, 1000),
        "dispersion_2000": lambda: genes_top_dispersion(adata, 2000),
        "dispersion_3000": lambda: genes_top_dispersion(adata, args.max_genes),
        "variance_3000": lambda: genes_top_variance(adata, args.max_genes),
        "detection_balanced_dispersion_3000": lambda: genes_detection_balanced_dispersion(adata, args.max_genes),
        "sparse_gene_coverage_3000": lambda: genes_sparse_gene_coverage_scores(adata, args.max_genes, top_per_gene=3),
        "sparse_gene_coverage_aggressive_3000": lambda: genes_sparse_gene_coverage_scores(adata, args.max_genes, top_per_gene=10),
        "pca_elbow_0.5_3000": lambda: genes_pca_loadings(adata, args.max_genes, n_components=50, elbow_fraction=0.5),
        "pca_elbow_1.0_3000": lambda: genes_pca_loadings(adata, args.max_genes, n_components=50, elbow_fraction=1.0),
        "pca_elbow_1.5_3000": lambda: genes_pca_loadings(adata, args.max_genes, n_components=50, elbow_fraction=1.5),
        "pca_50pcs_3000": lambda: genes_pca_loadings(adata, args.max_genes, n_components=50, elbow_fraction=10.0),
    }

    all_scores = []
    summaries = []
    for method, fn in methods.items():
        print(f"[feature-benchmark] starting {method}", flush=True)
        start = time.time()
        mem0 = _rss_mb()
        genes = pd.Index(fn()).drop_duplicates()[: args.max_genes]
        seconds = time.time() - start
        mem1 = _rss_mb()
        pd.DataFrame({"gene": genes}).to_csv(outdir / f"{method}.genes.tsv", sep="\t", index=False)
        score = _score_method(method, genes, marker_db, args.min_hits)
        score["genes_selected"] = len(genes)
        score["method_seconds"] = seconds
        score["rss_delta_mb"] = mem1 - mem0
        all_scores.append(score)
        summaries.append(
            {
                "method": method,
                "genes_selected": len(genes),
                "cell_states_evaluated": score.shape[0],
                "cell_states_passing_min_10": int(score["passes_min_10"].sum()),
                "cell_states_failing_min_10": int((~score["passes_min_10"]).sum()),
                "minimum_marker_hits": int(score["marker_genes_captured"].min()) if not score.empty else 0,
                "median_marker_hits": float(score["marker_genes_captured"].median()) if not score.empty else 0,
                "mean_marker_hits": float(score["marker_genes_captured"].mean()) if not score.empty else 0,
                "method_seconds": seconds,
                "rss_delta_mb": mem1 - mem0,
                "preprocess_seconds": preprocess_seconds,
                "markerfinder_seconds": marker_seconds,
                "marker_database_states": int(marker_db["cell_state"].nunique()),
                "marker_database_genes_unique": int(marker_db["marker"].nunique()),
                "marker_database_rows": int(marker_db.shape[0]),
            }
        )
        print(f"[feature-benchmark] finished {method}: {summaries[-1]}", flush=True)

    scores = pd.concat(all_scores, axis=0, ignore_index=True)
    scores.to_csv(outdir / "variable_gene_method_marker_capture_by_state.tsv", sep="\t", index=False)
    summary = pd.DataFrame(summaries).sort_values(
        ["cell_states_passing_min_10", "minimum_marker_hits", "median_marker_hits", "method_seconds"],
        ascending=[False, False, False, True],
        kind="mergesort",
    )
    summary.to_csv(outdir / "variable_gene_method_marker_capture_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
