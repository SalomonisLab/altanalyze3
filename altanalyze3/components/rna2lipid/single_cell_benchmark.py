from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

from .api import Rna2LipidBundle
from .pipeline import dump_json, json_ready, utc_timestamp


def _require_scanpy():
    try:
        import scanpy as sc  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "scanpy is required for the single-cell benchmark workflow."
        ) from exc
    return sc


@dataclass(frozen=True)
class StrategyResult:
    name: str
    predictions: pd.DataFrame
    metrics: Dict[str, Any]


def _load_single_cell_input(path: Path) -> ad.AnnData:
    sc = _require_scanpy()
    suffix = path.suffix.lower()
    if suffix == ".h5ad":
        adata = sc.read_h5ad(path)
    elif suffix == ".h5":
        adata = sc.read_10x_h5(path)
    else:
        raise ValueError(f"Unsupported single-cell input format: {path}")

    if "gene_symbol" not in adata.var.columns:
        adata.var["gene_symbol"] = adata.var_names.astype(str)
    adata.var["gene_symbol"] = adata.var["gene_symbol"].astype(str).str.strip()
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    return adata


def _prepare_rna_embedding(
    adata: ad.AnnData,
    *,
    random_state: int,
    n_top_genes: int = 3000,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    coarse_clusters: int | None = None,
    metacell_target_size: int = 50,
) -> Dict[str, Any]:
    sc = _require_scanpy()
    work = adata.copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    work.layers["normalized_log1p"] = work.X.copy()
    sc.pp.highly_variable_genes(work, n_top_genes=min(n_top_genes, work.n_vars), flavor="seurat")
    sc.pp.pca(work, n_comps=min(n_pcs, max(2, work.n_obs - 1)), use_highly_variable=True)
    sc.pp.neighbors(work, n_neighbors=min(n_neighbors, max(2, work.n_obs - 1)), n_pcs=min(n_pcs, work.obsm["X_pca"].shape[1]))
    sc.tl.umap(work, random_state=random_state)

    if coarse_clusters is None:
        coarse_clusters = max(8, min(20, int(round(math.sqrt(max(1, work.n_obs / 20.0)) * 1.5))))
    coarse_clusters = min(coarse_clusters, work.n_obs)
    coarse_kmeans = MiniBatchKMeans(n_clusters=coarse_clusters, random_state=random_state, batch_size=1024, n_init="auto")
    work.obs["benchmark_cluster"] = coarse_kmeans.fit_predict(work.obsm["X_pca"]).astype(str)

    eval_clusters = max(coarse_clusters + 8, min(48, coarse_clusters * 2))
    eval_clusters = min(eval_clusters, work.n_obs)
    eval_kmeans = MiniBatchKMeans(n_clusters=eval_clusters, random_state=random_state + 17, batch_size=1024, n_init="auto")
    work.obs["benchmark_eval_cluster"] = eval_kmeans.fit_predict(work.obsm["X_pca"]).astype(str)

    n_metacells = max(2, int(math.ceil(work.n_obs / float(max(10, metacell_target_size)))))
    n_metacells = min(n_metacells, work.n_obs)
    metacell_kmeans = MiniBatchKMeans(n_clusters=n_metacells, random_state=random_state, batch_size=1024, n_init="auto")
    work.obs["benchmark_metacell"] = metacell_kmeans.fit_predict(work.obsm["X_pca"]).astype(str)

    return {
        "adata": work,
        "n_coarse_clusters": int(coarse_clusters),
        "n_eval_clusters": int(eval_clusters),
        "n_metacells": int(n_metacells),
    }


def _predict_from_adata(bundle: Rna2LipidBundle, adata: ad.AnnData) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    gene_labels = adata.var["gene_symbol"].astype(str).str.strip().tolist()
    matched_positions = bundle._matched_model_positions(gene_labels)  # noqa: SLF001
    if not matched_positions:
        raise ValueError("No model genes matched the single-cell input.")

    matrix = adata.layers["normalized_log1p"] if "normalized_log1p" in adata.layers else adata.X
    aligned = bundle._align_matrix_chunk(matrix, matched_positions)  # noqa: SLF001
    predictions = bundle._predict_aligned_matrix(aligned)  # noqa: SLF001
    predictions.index = adata.obs_names
    summary = bundle._build_summary(  # noqa: SLF001
        input_rows=int(adata.n_obs),
        matched_genes=len(matched_positions),
        input_kind="single_cell_benchmark",
    )
    summary["output_rows"] = int(predictions.shape[0])
    return predictions, summary


def _row_normalize_connectivities(connectivities) -> sparse.csr_matrix:
    if not sparse.issparse(connectivities):
        connectivities = sparse.csr_matrix(np.asarray(connectivities))
    csr = connectivities.tocsr().astype(np.float64)
    row_sums = np.asarray(csr.sum(axis=1)).ravel()
    inv = np.divide(1.0, row_sums, out=np.zeros_like(row_sums), where=row_sums > 0)
    return sparse.diags(inv) @ csr


def _cluster_assign(predictions: pd.DataFrame, labels: pd.Series) -> pd.DataFrame:
    cluster_means = predictions.groupby(labels, sort=False).mean()
    assigned = cluster_means.loc[labels].copy()
    assigned.index = predictions.index
    return assigned


def _knn_smooth(predictions: pd.DataFrame, connectivities: sparse.csr_matrix, alpha: float) -> pd.DataFrame:
    weights = _row_normalize_connectivities(connectivities)
    pred = predictions.to_numpy(dtype=float)
    neighbor_mean = weights @ pred
    smoothed = alpha * pred + (1.0 - alpha) * neighbor_mean
    return pd.DataFrame(smoothed, index=predictions.index, columns=predictions.columns)


def _neighbor_pearson_similarity(predictions: pd.DataFrame, connectivities: sparse.csr_matrix) -> float:
    pred = predictions.to_numpy(dtype=float)
    csr = connectivities.tocsr()
    scores: List[float] = []
    for idx in range(csr.shape[0]):
        start, stop = csr.indptr[idx], csr.indptr[idx + 1]
        neighbors = csr.indices[start:stop]
        neighbors = neighbors[neighbors != idx]
        if neighbors.size == 0:
            continue
        centroid = pred[neighbors].mean(axis=0)
        if np.std(pred[idx]) == 0 or np.std(centroid) == 0:
            continue
        scores.append(float(np.corrcoef(pred[idx], centroid)[0, 1]))
    return float(np.mean(scores)) if scores else float("nan")


def _cluster_explained_variance(predictions: pd.DataFrame, labels: pd.Series) -> float:
    pred = predictions.to_numpy(dtype=float)
    total_var = pred.var(axis=0, ddof=0)
    groups = labels.astype(str).to_numpy()
    unique_groups = pd.Index(groups).unique().tolist()
    cluster_means = np.vstack([pred[groups == group].mean(axis=0) for group in unique_groups])
    cluster_sizes = np.asarray([(groups == group).sum() for group in unique_groups], dtype=float)
    grand_mean = pred.mean(axis=0)
    between_var = ((cluster_sizes[:, None] * (cluster_means - grand_mean) ** 2).sum(axis=0)) / max(1.0, cluster_sizes.sum())
    explained = np.divide(between_var, total_var, out=np.zeros_like(between_var), where=total_var > 0)
    return float(np.mean(explained))


def _split_half_reproducibility(predictions: pd.DataFrame, labels: pd.Series, *, random_state: int) -> float:
    rng = np.random.default_rng(random_state)
    pred = predictions
    scores: List[float] = []
    for label, idx in labels.groupby(labels).groups.items():
        cell_ids = np.asarray(list(idx))
        if cell_ids.size < 10:
            continue
        shuffled = cell_ids.copy()
        rng.shuffle(shuffled)
        midpoint = shuffled.size // 2
        left = shuffled[:midpoint]
        right = shuffled[midpoint:]
        if left.size < 3 or right.size < 3:
            continue
        mean_left = pred.loc[left].mean(axis=0).to_numpy(dtype=float)
        mean_right = pred.loc[right].mean(axis=0).to_numpy(dtype=float)
        if np.std(mean_left) == 0 or np.std(mean_right) == 0:
            continue
        scores.append(float(np.corrcoef(mean_left, mean_right)[0, 1]))
    return float(np.mean(scores)) if scores else float("nan")


def _effective_rank(predictions: pd.DataFrame) -> float:
    pred = predictions.to_numpy(dtype=float)
    pred = pred - pred.mean(axis=0, keepdims=True)
    singular_values = np.linalg.svd(pred, full_matrices=False, compute_uv=False)
    singular_values = singular_values[singular_values > 0]
    if singular_values.size == 0:
        return 0.0
    probs = singular_values / singular_values.sum()
    entropy = -np.sum(probs * np.log(probs))
    return float(np.exp(entropy))


def _strategy_metrics(
    predictions: pd.DataFrame,
    *,
    baseline_predictions: pd.DataFrame,
    eval_cluster_labels: pd.Series,
    connectivities: sparse.csr_matrix,
    random_state: int,
) -> Dict[str, Any]:
    values = predictions.to_numpy(dtype=float)
    scaled = StandardScaler().fit_transform(values)
    silhouette = None
    unique_clusters = pd.Index(eval_cluster_labels.astype(str)).nunique()
    if unique_clusters > 1 and predictions.shape[0] > unique_clusters:
        silhouette = float(silhouette_score(scaled, eval_cluster_labels.astype(str).to_numpy()))

    total_variance = float(values.var(axis=0, ddof=0).mean())
    baseline_variance = float(baseline_predictions.to_numpy(dtype=float).var(axis=0, ddof=0).mean())
    variance_retention = total_variance / baseline_variance if baseline_variance > 0 else None

    return {
        "silhouette_score": silhouette,
        "neighbor_pearson_similarity": _neighbor_pearson_similarity(predictions, connectivities),
        "eval_cluster_explained_variance": _cluster_explained_variance(predictions, eval_cluster_labels),
        "split_half_reproducibility": _split_half_reproducibility(predictions, eval_cluster_labels, random_state=random_state),
        "mean_total_variance": total_variance,
        "variance_retention_vs_baseline": variance_retention,
        "effective_rank": _effective_rank(predictions),
    }


def _balanced_score(metrics_df: pd.DataFrame) -> pd.Series:
    weights = {
        "silhouette_score": 0.30,
        "neighbor_pearson_similarity": 0.20,
        "eval_cluster_explained_variance": 0.15,
        "split_half_reproducibility": 0.10,
        "variance_retention_vs_baseline": 0.25,
    }
    normalized_parts = []
    for column, weight in weights.items():
        values = metrics_df[column].astype(float)
        min_value = float(values.min())
        max_value = float(values.max())
        if math.isclose(min_value, max_value):
            normalized = pd.Series(1.0, index=values.index)
        else:
            normalized = (values - min_value) / (max_value - min_value)
        normalized_parts.append(normalized * weight)
    return sum(normalized_parts)


def benchmark_single_cell_assignments(
    *,
    input_path: str | Path,
    output_dir: str | Path,
    bundle_path: str | Path,
    random_state: int = 0,
    knn_alpha: float = 0.5,
    coarse_clusters: int | None = None,
    metacell_target_size: int = 50,
) -> Dict[str, Any]:
    input_path = Path(input_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = _load_single_cell_input(input_path)
    bundle = Rna2LipidBundle.load(bundle_path)
    embedding = _prepare_rna_embedding(
        adata,
        random_state=random_state,
        coarse_clusters=coarse_clusters,
        metacell_target_size=metacell_target_size,
    )
    work = embedding["adata"]
    baseline_predictions, prediction_summary = _predict_from_adata(bundle, work)

    cluster_labels = work.obs["benchmark_cluster"].astype(str)
    eval_cluster_labels = work.obs["benchmark_eval_cluster"].astype(str)
    metacell_labels = work.obs["benchmark_metacell"].astype(str)
    connectivities = work.obsp["connectivities"].tocsr()

    strategies = [
        StrategyResult(
            name="per_cell",
            predictions=baseline_predictions,
            metrics={},
        ),
        StrategyResult(
            name=f"knn_smooth_alpha_{knn_alpha:.2f}",
            predictions=_knn_smooth(baseline_predictions, connectivities, alpha=knn_alpha),
            metrics={},
        ),
        StrategyResult(
            name=f"metacell_mean_{metacell_target_size}",
            predictions=_cluster_assign(baseline_predictions, metacell_labels),
            metrics={},
        ),
        StrategyResult(
            name=f"cluster_mean_{embedding['n_coarse_clusters']}",
            predictions=_cluster_assign(baseline_predictions, cluster_labels),
            metrics={},
        ),
    ]

    metric_rows = []
    strategy_by_name: Dict[str, StrategyResult] = {}
    for strategy in strategies:
        metrics = _strategy_metrics(
            strategy.predictions,
            baseline_predictions=baseline_predictions,
            eval_cluster_labels=eval_cluster_labels,
            connectivities=connectivities,
            random_state=random_state,
        )
        metric_rows.append({"strategy": strategy.name, **metrics})
        strategy_by_name[strategy.name] = StrategyResult(
            name=strategy.name,
            predictions=strategy.predictions,
            metrics=metrics,
        )

    metrics_df = pd.DataFrame(metric_rows).set_index("strategy")
    metrics_df["balanced_score"] = _balanced_score(metrics_df)
    metrics_df = metrics_df.sort_values("balanced_score", ascending=False)
    best_strategy_name = str(metrics_df.index[0])
    best_strategy = strategy_by_name[best_strategy_name]

    metrics_path = output_dir / "benchmark_metrics.tsv"
    metrics_df.to_csv(metrics_path, sep="\t")

    rna_out = work.copy()
    rna_out.obs["benchmark_cluster"] = cluster_labels
    rna_out.obs["benchmark_eval_cluster"] = eval_cluster_labels
    rna_out.obs["benchmark_metacell"] = metacell_labels
    rna_out.write_h5ad(output_dir / "bos_rna_benchmark_reference.h5ad", compression="lzf")

    best_adata = ad.AnnData(
        X=best_strategy.predictions.to_numpy(dtype=np.float32),
        obs=work.obs.copy(),
        var=pd.DataFrame(index=best_strategy.predictions.columns),
        obsm={"X_umap": work.obsm["X_umap"].copy(), "X_pca": work.obsm["X_pca"].copy()},
        uns={
            "benchmark_summary": {
                "strategy": best_strategy_name,
                "metrics": json_ready(best_strategy.metrics),
            }
        },
    )
    best_adata.write_h5ad(output_dir / "best_strategy_lipid_predictions.h5ad", compression="lzf")

    cluster_means = best_strategy.predictions.groupby(cluster_labels, sort=False).mean()
    cluster_means.to_csv(output_dir / "best_strategy_cluster_means.tsv", sep="\t")

    summary = {
        "created_at": utc_timestamp(),
        "input_path": str(input_path),
        "bundle_path": str(Path(bundle_path).resolve()),
        "prediction_summary": prediction_summary,
        "cell_count": int(work.n_obs),
        "gene_count": int(work.n_vars),
        "coarse_cluster_count": int(embedding["n_coarse_clusters"]),
        "evaluation_cluster_count": int(embedding["n_eval_clusters"]),
        "metacell_count": int(embedding["n_metacells"]),
        "benchmark_metrics_path": str(metrics_path),
        "best_strategy": best_strategy_name,
        "best_strategy_metrics": json_ready(best_strategy.metrics),
        "all_metrics": json_ready(metrics_df.reset_index().to_dict(orient="records")),
        "interpretation": {
            "goal": "Improve single-cell lipid assignments by increasing RNA-state coherence while retaining per-cell variation.",
            "note": "Higher silhouette, neighbor Pearson similarity, eval-cluster explained variance, and split-half reproducibility indicate stronger state consistency; variance retention penalizes over-smoothing.",
        },
    }
    dump_json(output_dir / "benchmark_summary.json", json_ready(summary))
    return summary
