"""
Metacell construction module for AltAnalyze3.

This module creates metacells from scRNA-seq AnnData objects with an emphasis
on speed, controllable granularity, and optional restrictions to predefined
cell-state and donor boundaries. It supports three principal strategies:

1. Boundary-aware metacells (per sample and/or cell type).
2. Global metacells across all selected cells.
3. Sample-specific metacells without predefined cell states.

The resulting AnnData contains metacell expression profiles and, in the
unstructured section, a mapping from each metacell to its member cell barcodes.
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.cluster import MiniBatchKMeans
from tqdm.auto import tqdm


###############################################################################
# Data classes and helpers
###############################################################################


@dataclass
class MetacellParams:
    """Collected parameters controlling metacell generation."""

    target_size: int
    min_size: int
    max_size: int
    preserve_small: bool
    aggregation: str
    graph_algorithm: str
    resolution: Optional[float]
    resolution_steps: int
    resolution_tolerance: float
    n_neighbors: int
    neighbor_method: str
    neighbor_metric: str
    n_top_genes: int
    n_pcs: int
    pca_svd_solver: str
    hvg_layer: Optional[str]
    expression_layer: Optional[str]
    use_raw: bool
    random_state: int
    random_metacell_count: Optional[int]
    random_cells_per_metacell: int
    random_sampling_with_replacement: bool


def configure_logging(verbose: bool = False) -> None:
    """Configure root logger with optional verbosity."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(
        description="Construct metacells from an h5ad file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="Input h5ad path.")
    parser.add_argument("--output", required=True, help="Output h5ad path for metacells.")

    # Strategy options
    parser.add_argument(
        "--mode",
        choices=["global", "boundary", "sample"],
        default="global",
        help=(
            "Metacell construction strategy. 'boundary' creates metacells within "
            "each boundary column combination, 'sample' restricts to specific samples."
        ),
    )
    parser.add_argument(
        "--boundary-columns",
        nargs="+",
        default=None,
        help="Columns in obs defining boundaries (e.g., sample cell_type).",
    )
    parser.add_argument("--sample-column", default=None, help="obs column storing sample / donor IDs.")
    parser.add_argument("--cell-type-column", default=None, help="obs column storing cell-state annotations.")
    parser.add_argument(
        "--restrict-sample",
        action="append",
        default=None,
        help="Restrict analysis to the provided sample value (can be repeated).",
    )
    parser.add_argument(
        "--restrict-cell-type",
        action="append",
        default=None,
        help="Restrict analysis to the provided cell-type value (can be repeated).",
    )

    # Graph construction
    parser.add_argument("--graph-algorithm", choices=["leiden", "louvain", "kmeans", "random"], default="leiden")
    parser.add_argument("--resolution", type=float, default=None, help="Resolution for Leiden/Louvain (auto if unset).")
    parser.add_argument("--resolution-steps", type=int, default=5, help="Maximum tuning steps for auto resolution.")
    parser.add_argument(
        "--resolution-tolerance",
        type=float,
        default=0.15,
        help="Acceptable fractional deviation between target and obtained number of clusters.",
    )
    parser.add_argument("--target-metacell-size", type=int, default=50, help="Desired number of cells per metacell.")
    parser.add_argument("--min-metacell-size", type=int, default=80, help="Minimum cells allowed per metacell.")
    parser.add_argument("--max-metacell-size", type=int, default=400, help="Maximum cells allowed per metacell.")
    parser.add_argument(
        "--preserve-small-clusters",
        action="store_true",
        help="Keep clusters smaller than --min-metacell-size instead of merging them.",
    )
    parser.add_argument("--n-top-genes", type=int, default=3000, help="Number of HVGs to use for graph construction.")
    parser.add_argument("--n-pcs", type=int, default=50, help="Number of principal components.")
    parser.add_argument("--pca-svd-solver", choices=["auto", "arpack", "randomized"], default="auto")
    parser.add_argument("--n-neighbors", type=int, default=30, help="Number of neighbors for the KNN graph.")
    parser.add_argument("--neighbor-method", choices=["auto", "umap", "gauss"], default="auto")
    parser.add_argument("--neighbor-metric", default="euclidean", help="Metric for neighbor search.")

    # Expression data handling
    parser.add_argument("--expression-layer", default=None, help="Layer containing expression counts (default: X).")
    parser.add_argument("--hvg-layer", default=None, help="Layer to use for HVG selection (defaults to expression data).")
    parser.add_argument(
        "--use-raw",
        action="store_true",
        help="Use AnnData.raw for graph construction (and expression when no layer provided).",
    )
    parser.add_argument(
        "--aggregation",
        choices=["sum"],
        default="sum",
        help="Aggregation function for metacell expression profiles (sums only).",
    )

    parser.add_argument(
        "--random-metacell-count",
        type=int,
        default=50,
        help="Number of random metacells to generate per group (used when --graph-algorithm random).",
    )
    parser.add_argument(
        "--random-cells-per-metacell",
        type=int,
        default=5,
        help="Number of cells combined for each random metacell.",
    )
    parser.add_argument(
        "--random-sampling-with-replacement",
        action="store_true",
        help="Allow sampling the same cell multiple times when building random metacells.",
    )

    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for Scanpy/Numba where applicable.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging.")

    return parser.parse_args(argv)


###############################################################################
# Data utilities
###############################################################################


def get_obs_matrix(
    adata: ad.AnnData,
    *,
    layer: Optional[str],
    use_raw: bool,
) -> Tuple[sp.spmatrix | np.ndarray, pd.DataFrame]:
    """Return the observation matrix (cells x genes) and corresponding var DataFrame."""
    if use_raw:
        if adata.raw is None:
            raise ValueError("Requested use_raw but AnnData.raw is not set.")
        matrix = adata.raw.X
        var = adata.raw.var.copy()
    elif layer:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in AnnData.")
        matrix = adata.layers[layer]
        var = adata.var.copy()
    else:
        matrix = adata.X
        var = adata.var.copy()

    if sp.issparse(matrix):
        matrix = matrix.tocsr()

    return matrix, var


def subset_adata(
    adata: ad.AnnData,
    sample_column: Optional[str],
    cell_type_column: Optional[str],
    restrict_samples: Optional[List[str]],
    restrict_cell_types: Optional[List[str]],
) -> ad.AnnData:
    """Apply optional sample / cell-type restrictions."""
    mask = np.ones(adata.n_obs, dtype=bool)

    if restrict_samples:
        if not sample_column:
            raise ValueError("--restrict-sample provided but --sample-column is missing.")
        if sample_column not in adata.obs.columns:
            raise KeyError(f"Sample column '{sample_column}' not found in obs.")
        mask &= adata.obs[sample_column].isin(restrict_samples).to_numpy()

    if restrict_cell_types:
        if not cell_type_column:
            raise ValueError("--restrict-cell-type provided but --cell-type-column is missing.")
        if cell_type_column not in adata.obs.columns:
            raise KeyError(f"Cell-type column '{cell_type_column}' not found in obs.")
        mask &= adata.obs[cell_type_column].isin(restrict_cell_types).to_numpy()

    if not mask.any():
        raise ValueError("No cells left after applying restrictions.")

    if mask.all():
        return adata

    logging.info("Restricting analysis to %d of %d cells.", mask.sum(), adata.n_obs)
    return adata[mask].copy()


def group_iterable(
    adata: ad.AnnData,
    boundary_columns: Sequence[str],
    grouped: Optional[pd.core.groupby.generic.DataFrameGroupBy] = None,
) -> Iterable[Tuple[Tuple[str, ...], ad.AnnData]]:
    """Yield (group_key, group_adata) pairs for the requested boundary columns."""
    if not boundary_columns:
        yield ("all",), adata
        return

    missing = [col for col in boundary_columns if col not in adata.obs.columns]
    if missing:
        raise KeyError(f"Boundary columns missing from obs: {', '.join(missing)}")

    if grouped is None:
        grouped = adata.obs.groupby(list(boundary_columns), observed=True, dropna=False)
    for key, group_df in grouped:
        key_tuple = key if isinstance(key, tuple) else (key,)
        subgroup = adata[group_df.index].copy()
        yield key_tuple, subgroup


###############################################################################
# Graph construction and clustering
###############################################################################


def prepare_graph_data(
    adata: ad.AnnData,
    params: MetacellParams,
) -> ad.AnnData:
    """Prepare data (HVGs, PCA, neighbors) used for clustering."""
    logging.debug("Selecting %d HVGs (layer=%s use_raw=%s).", params.n_top_genes, params.hvg_layer, params.use_raw)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=params.n_top_genes,
        flavor="seurat_v3",
        subset=True,
        layer=params.hvg_layer,
    )

    if sp.issparse(adata.X):
        adata.X = adata.X.astype(np.float32)
    else:
        adata.X = np.asarray(adata.X, dtype=np.float32)

    logging.debug("Running PCA (n_pcs=%d).", params.n_pcs)
    solver = params.pca_svd_solver
    if solver == "auto":
        solver = "arpack" if sp.issparse(adata.X) else "randomized"
    elif sp.issparse(adata.X) and solver == "randomized":
        logging.debug("Switching PCA solver to 'arpack' for sparse input to avoid densification.")
        solver = "arpack"
    sc.tl.pca(
        adata,
        n_comps=params.n_pcs,
        svd_solver=solver,
        random_state=params.random_state,
    )

    logging.debug(
        "Building neighbors (n_neighbors=%d, method=%s, metric=%s).",
        params.n_neighbors,
        params.neighbor_method,
        params.neighbor_metric,
    )
    neighbor_kwargs = dict(
        n_pcs=params.n_pcs,
        n_neighbors=params.n_neighbors,
        metric=params.neighbor_metric,
        random_state=params.random_state,
    )
    if params.neighbor_method != "auto":
        neighbor_kwargs["method"] = params.neighbor_method
    sc.pp.neighbors(adata, **neighbor_kwargs)
    return adata


def _run_leiden_or_louvain(
    adata: ad.AnnData,
    params: MetacellParams,
    target_clusters: int,
) -> np.ndarray:
    """Run Leiden/Louvain with optional resolution tuning."""
    algorithm = params.graph_algorithm
    random_state = params.random_state

    if params.resolution is not None:
        res = params.resolution
    else:
        # Heuristic starting point
        res = max(0.1, target_clusters / max(adata.n_obs / 50, 1))

    last_labels: Optional[np.ndarray] = None

    for step in range(params.resolution_steps):
        key_added = f"_tmp_{algorithm}"
        if algorithm == "leiden":
            sc.tl.leiden(
                adata,
                resolution=res,
                random_state=random_state,
                key_added=key_added,
            )
        else:
            sc.tl.louvain(
                adata,
                resolution=res,
                random_state=random_state,
                key_added=key_added,
            )
        labels = adata.obs[key_added].to_numpy().astype(str)
        n_clusters = np.unique(labels).size
        logging.debug("%s resolution %.4f -> %d clusters (target %d).", algorithm, res, n_clusters, target_clusters)
        del adata.obs[key_added]

        last_labels = labels
        if params.resolution is not None:
            break

        if target_clusters <= 1:
            break

        deviation = abs(n_clusters - target_clusters) / target_clusters
        if deviation <= params.resolution_tolerance:
            break

        if n_clusters < target_clusters:
            res *= 1.6
        else:
            res /= 1.6

    if last_labels is None:
        raise RuntimeError("Failed to compute clustering labels.")

    return last_labels


def initial_labels(
    analysis_data: ad.AnnData,
    params: MetacellParams,
) -> np.ndarray:
    """Generate initial cluster labels using the chosen algorithm."""
    n_cells = analysis_data.n_obs
    target_clusters = max(1, int(np.ceil(n_cells / params.target_size)))

    if params.graph_algorithm in {"leiden", "louvain"}:
        return _run_leiden_or_louvain(analysis_data, params, target_clusters)

    # KMeans fallback
    n_clusters = max(1, target_clusters)
    pca = analysis_data.obsm["X_pca"]
    kmeans = MiniBatchKMeans(
        n_clusters=n_clusters,
        batch_size=min(4096, max(2000, params.target_size * 4)),
        random_state=params.random_state,
        n_init="auto",
    )
    labels = kmeans.fit_predict(pca)
    return labels.astype(str)


###############################################################################
# Cluster size enforcement
###############################################################################


def enforce_size_bounds(
    labels: np.ndarray,
    embeddings: np.ndarray,
    params: MetacellParams,
) -> List[np.ndarray]:
    """Split oversized clusters and optionally merge undersized ones."""
    # Map original labels to integer groups
    uniques = np.unique(labels)
    label_to_idx = {label: idx for idx, label in enumerate(uniques)}
    group_indices = [np.where(labels == label)[0] for label in uniques]

    max_size = params.max_size
    target = params.target_size
    min_size = params.min_size
    preserve_small = params.preserve_small

    rng = np.random.default_rng(params.random_state)
    final_groups: List[np.ndarray] = []

    for group_idx, cell_indices in enumerate(group_indices):
        current = cell_indices
        if current.size <= max_size:
            final_groups.append(current)
            continue

        # Split oversized clusters using MiniBatchKMeans on PCA embeddings.
        n_needed = int(np.ceil(current.size / target))
        if n_needed <= 1:
            final_groups.append(current)
            continue

        n_clusters = min(n_needed, max(current.size // max(1, min_size), 1))
        n_clusters = max(2, n_clusters)
        subset_emb = embeddings[current]
        splitter = MiniBatchKMeans(
            n_clusters=n_clusters,
            batch_size=min(4096, max(2000, params.target_size * 4)),
            random_state=params.random_state,
            n_init="auto",
        )
        sub_labels = splitter.fit_predict(subset_emb)
        for sub_label in np.unique(sub_labels):
            member_cells = current[sub_labels == sub_label]
            if member_cells.size == 0:
                continue
            final_groups.append(member_cells)

    if preserve_small or min_size <= 1:
        return final_groups

    # Merge undersized clusters based on centroid similarity.
    centroids = [
        embeddings[indices].mean(axis=0) if indices.size > 0 else np.zeros(embeddings.shape[1])
        for indices in final_groups
    ]
    sizes = [indices.size for indices in final_groups]

    # Keep iteratively merging until all clusters satisfy the minimum size.
    changed = True
    while changed:
        changed = False
        small_indices = [i for i, size in enumerate(sizes) if size < min_size]
        if not small_indices:
            break

        for idx in small_indices:
            if sizes[idx] >= min_size:
                continue

            # Find nearest cluster (excluding itself)
            if len(final_groups) == 1:
                break

            centroid = centroids[idx]
            dists = []
            for j, other_centroid in enumerate(centroids):
                if j == idx or sizes[j] == 0:
                    dists.append(np.inf)
                else:
                    dists.append(np.linalg.norm(centroid - other_centroid))
            partner = int(np.argmin(dists))
            if not np.isfinite(dists[partner]):
                continue

            merged_indices = np.concatenate([final_groups[idx], final_groups[partner]])
            rng.shuffle(merged_indices)
            final_groups[partner] = merged_indices
            centroids[partner] = embeddings[merged_indices].mean(axis=0)
            sizes[partner] = merged_indices.size

            # Clear the merged cluster.
            final_groups[idx] = np.array([], dtype=int)
            sizes[idx] = 0
            changed = True

        # Remove emptied clusters.
        final_groups = [grp for grp in final_groups if grp.size > 0]
        centroids = [embeddings[grp].mean(axis=0) for grp in final_groups]
        sizes = [grp.size for grp in final_groups]

    return final_groups


def generate_random_groups(
    n_cells: int,
    count: int,
    cells_per_group: int,
    with_replacement: bool,
    rng: np.random.Generator,
) -> List[np.ndarray]:
    """Generate random cell index groups for metacell aggregation."""
    if count <= 0 or cells_per_group <= 0 or n_cells == 0:
        return []

    indices = np.arange(n_cells, dtype=int)
    groups: List[np.ndarray] = []

    # If sampling with replacement or requesting more cells than available,
    # draw with replacement for each group.
    if with_replacement or cells_per_group > n_cells:
        for _ in range(count):
            chosen = rng.choice(indices, size=cells_per_group, replace=True)
            groups.append(chosen)
        return groups

    # Otherwise cycle through permutations without replacement per group.
    needed = count * cells_per_group
    pooled = rng.permutation(indices)
    while pooled.size < needed:
        pooled = np.concatenate([pooled, rng.permutation(indices)])

    for start in range(0, count * cells_per_group, cells_per_group):
        segment = pooled[start : start + cells_per_group]
        if segment.size < cells_per_group:
            supplement = rng.choice(indices, size=cells_per_group - segment.size, replace=False)
            segment = np.concatenate([segment, supplement])
        groups.append(segment.astype(int))

    return groups


###############################################################################
# Metacell aggregation
###############################################################################


def aggregate_expression(
    matrix: sp.spmatrix | np.ndarray,
    groups: List[np.ndarray],
    aggregation: str,
) -> sp.spmatrix | np.ndarray:
    """Aggregate expression matrix based on the provided cell index groups."""
    if aggregation != "sum":
        raise ValueError("Only sum aggregation is supported for metacell construction.")

    n_groups = len(groups)
    if n_groups == 0:
        raise ValueError("No metacell groups generated.")

    if sp.issparse(matrix):
        aggregated_rows = []
        for cells in groups:
            cells = np.asarray(cells, dtype=int)
            if cells.size == 0:
                aggregated_rows.append(sp.csr_matrix((1, matrix.shape[1]), dtype=matrix.dtype))
                continue
            subset = matrix[cells]
            agg_row = subset.sum(axis=0)
            aggregated_rows.append(sp.csr_matrix(agg_row))
        aggregated = sp.vstack(aggregated_rows, format="csr")
    else:
        aggregated_rows = []
        for cells in groups:
            cells = np.asarray(cells, dtype=int)
            if cells.size == 0:
                aggregated_rows.append(np.zeros((1, matrix.shape[1]), dtype=matrix.dtype))
                continue
            subset = matrix[cells]
            agg_row = np.sum(subset, axis=0, keepdims=True)
            aggregated_rows.append(agg_row)
        aggregated = np.vstack(aggregated_rows)

    return aggregated


def build_metacell_obs(
    groups: List[np.ndarray],
    metacell_ids: List[str],
    cell_obs: pd.DataFrame,
    sample_column: Optional[str],
    cell_type_column: Optional[str],
    boundary_info: Dict[str, str],
) -> pd.DataFrame:
    """Construct metacell-level metadata."""
    sizes = [cells.size for cells in groups]
    obs = pd.DataFrame({"metacell_id": metacell_ids, "metacell_size": sizes}, index=metacell_ids)

    for column, value in boundary_info.items():
        obs[column] = value

    def _mode(values: pd.Series) -> Optional[str]:
        if values.empty:
            return None
        counts = values.value_counts(dropna=False)
        if counts.empty:
            return None
        return counts.idxmax()

    if sample_column and sample_column in cell_obs.columns:
        dominant_samples = []
        for indices in groups:
            values = cell_obs.iloc[np.asarray(indices, dtype=int)][sample_column]
            dominant_samples.append(_mode(values))
        obs[f"dominant_{sample_column}"] = dominant_samples

    if cell_type_column and cell_type_column in cell_obs.columns:
        dominant_cell_types = []
        for indices in groups:
            values = cell_obs.iloc[np.asarray(indices, dtype=int)][cell_type_column]
            dominant_cell_types.append(_mode(values))
        obs[f"dominant_{cell_type_column}"] = dominant_cell_types

    return obs


def build_membership_df(
    groups: List[np.ndarray],
    metacell_ids: List[str],
    cell_names: Sequence[str],
    cell_obs: pd.DataFrame,
    sample_column: Optional[str],
    cell_type_column: Optional[str],
) -> pd.DataFrame:
    """Return a DataFrame linking metacell IDs to member cell barcodes."""
    records: List[Dict[str, str]] = []
    for metacell_id, indices in zip(metacell_ids, groups):
        for idx in indices:
            idx_int = int(idx)
            entry = {
                "metacell_id": metacell_id,
                "cell_barcode": str(cell_names[idx_int]),
            }
            if sample_column and sample_column in cell_obs.columns:
                entry[sample_column] = cell_obs.iloc[idx_int][sample_column]
            if cell_type_column and cell_type_column in cell_obs.columns:
                entry[cell_type_column] = cell_obs.iloc[idx_int][cell_type_column]
            records.append(entry)
    return pd.DataFrame.from_records(records)


###############################################################################
# Core pipeline
###############################################################################


def process_group(
    group: ad.AnnData,
    params: MetacellParams,
    metacell_start_idx: int,
    boundary_columns: Sequence[str],
    boundary_key: Tuple[str, ...],
    sample_column: Optional[str],
    cell_type_column: Optional[str],
) -> Tuple[ad.AnnData, pd.DataFrame, int]:
    """Compute metacells for a single group and return the aggregated AnnData."""
    rng = np.random.default_rng(params.random_state + metacell_start_idx)
    label_parts = [str(part) for part in boundary_key] if boundary_key else ["all"]
    label = " | ".join(label_parts)
    if len(label) > 48:
        label = label[:45] + "..."

    progress_disable = group.n_obs == 0
    with tqdm(
        total=4,
        desc=f"Group {label}",
        unit="step",
        leave=False,
        disable=progress_disable,
    ) as step_bar:
        if params.graph_algorithm == "random" and params.random_metacell_count:
            groups = generate_random_groups(
                n_cells=group.n_obs,
                count=params.random_metacell_count,
                cells_per_group=params.random_cells_per_metacell,
                with_replacement=params.random_sampling_with_replacement,
                rng=rng,
            )
            if not groups:
                logging.warning(
                    "Group %s has no cells; skipping metacell generation.",
                    boundary_key,
                )
                groups = []
        elif group.n_obs < max(params.min_size, 2):
            logging.warning(
                "Group %s has only %d cells; returning single metacell covering all cells.",
                boundary_key,
                group.n_obs,
            )
            groups = [np.arange(group.n_obs, dtype=int)]
        else:
            # Prepare analysis data for clustering
            analysis_matrix, analysis_var = get_obs_matrix(group, layer=params.hvg_layer, use_raw=params.use_raw)
            analysis = ad.AnnData(X=analysis_matrix, obs=group.obs.copy(), var=analysis_var)
            prepare_graph_data(analysis, params)

            labels = initial_labels(analysis, params)
            embeddings = analysis.obsm["X_pca"]
            groups = enforce_size_bounds(labels, embeddings, params)

        if not groups:
            raise ValueError(f"No metacell groups produced for boundary {boundary_key}.")

        step_bar.update(1)
        if not progress_disable:
            step_bar.set_postfix_str(f"groups={len(groups)}")

        # Aggregate expression using the requested matrix
        expr_matrix, expr_var = get_obs_matrix(group, layer=params.expression_layer, use_raw=params.use_raw)
        aggregated_matrix = aggregate_expression(expr_matrix, groups, params.aggregation)
        step_bar.update(1)
        if not progress_disable:
            step_bar.set_postfix_str(f"metacells={len(groups)}")

        # Assign metacell IDs
        metacell_ids = [f"MC_{idx:06d}" for idx in range(metacell_start_idx, metacell_start_idx + len(groups))]

        boundary_info = dict(zip(boundary_columns, boundary_key)) if boundary_columns else {}

        metacell_obs = build_metacell_obs(
            groups,
            metacell_ids,
            group.obs,
            sample_column=sample_column,
            cell_type_column=cell_type_column,
            boundary_info=boundary_info,
        )
        step_bar.update(1)
        if not progress_disable:
            step_bar.set_postfix_str("obs")

        metacell_var = expr_var.copy()
        metacell_adata = ad.AnnData(X=aggregated_matrix, obs=metacell_obs, var=metacell_var)
        metacell_adata.obs_names = metacell_ids

        membership_df = build_membership_df(
            groups,
            metacell_ids,
            group.obs_names.to_numpy(),
            group.obs,
            sample_column=sample_column,
            cell_type_column=cell_type_column,
        )
        step_bar.update(1)
        if not progress_disable:
            step_bar.set_postfix_str("membership")

    next_index = metacell_start_idx + len(groups)
    return metacell_adata, membership_df, next_index


def assemble_metacells(
    adata: ad.AnnData,
    params: MetacellParams,
    boundary_columns: Sequence[str],
    sample_column: Optional[str],
    cell_type_column: Optional[str],
) -> Tuple[ad.AnnData, pd.DataFrame]:
    """Generate metacells across all groups and assemble the final AnnData."""
    metacell_tables: List[ad.AnnData] = []
    membership_tables: List[pd.DataFrame] = []

    next_metacell_id = 0
    if not boundary_columns:
        total_groups = 1
        grouped = None
    else:
        missing = [col for col in boundary_columns if col not in adata.obs.columns]
        if missing:
            raise KeyError(f"Boundary columns missing from obs: {', '.join(missing)}")
        grouped = adata.obs.groupby(list(boundary_columns), observed=True, dropna=False)
        total_groups = grouped.ngroups

    group_iter = group_iterable(adata, boundary_columns, grouped=grouped)
    group_pbar = tqdm(
        group_iter,
        total=total_groups,
        desc="Metacell groups",
        unit="group",
        disable=total_groups <= 1,
        leave=False,
    )

    for boundary_key, group in group_pbar:
        if hasattr(group_pbar, "set_postfix_str"):
            group_pbar.set_postfix_str(f"cells={group.n_obs}")
        logging.info(
            "Processing group %s with %d cells.",
            boundary_key,
            group.n_obs,
        )
        metacell_table, membership_df, next_metacell_id = process_group(
            group=group,
            params=params,
            metacell_start_idx=next_metacell_id,
            boundary_columns=boundary_columns,
            boundary_key=boundary_key,
            sample_column=sample_column,
            cell_type_column=cell_type_column,
        )
        metacell_tables.append(metacell_table)
        membership_tables.append(membership_df)

    if hasattr(group_pbar, "close"):
        group_pbar.close()

    if not metacell_tables:
        raise RuntimeError("No metacells generated. Check filters and boundaries.")

    logging.info("Concatenating %d metacell batches.", len(metacell_tables))
    combined = ad.concat(metacell_tables, join="outer", merge="first")
    membership = pd.concat(membership_tables, ignore_index=True)
    return combined, membership


###############################################################################
# Entry point
###############################################################################


def resolve_boundary_columns(args: argparse.Namespace) -> List[str]:
    """Resolve boundary columns from CLI arguments and mode selection."""
    if args.boundary_columns:
        return list(args.boundary_columns)

    if args.mode == "boundary":
        inferred = [col for col in [args.sample_column, args.cell_type_column] if col]
        if not inferred:
            raise ValueError(
                "Boundary mode requested but neither --boundary-columns nor --sample-column/--cell-type-column provided."
            )
        return inferred
    return []


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    configure_logging(args.verbose)
    sc.settings.n_jobs = args.threads

    logging.info("Loading AnnData from %s", args.input)
    adata = ad.read_h5ad(args.input)
    logging.info("Loaded %d cells and %d genes.", adata.n_obs, adata.n_vars)

    if args.mode == "sample":
        if not args.sample_column:
            raise ValueError("Sample mode requires --sample-column.")
        if not args.restrict_sample or len(args.restrict_sample) != 1:
            raise ValueError("Sample mode requires exactly one --restrict-sample value.")

    restricted = subset_adata(
        adata,
        sample_column=args.sample_column,
        cell_type_column=args.cell_type_column,
        restrict_samples=args.restrict_sample,
        restrict_cell_types=args.restrict_cell_type,
    )

    boundary_columns = resolve_boundary_columns(args)

    if args.random_cells_per_metacell <= 0:
        raise ValueError("--random-cells-per-metacell must be a positive integer.")

    if args.graph_algorithm == "random":
        if args.random_metacell_count is None or args.random_metacell_count <= 0:
            raise ValueError("--random-metacell-count must be positive when using the random algorithm.")
        random_metacell_count = args.random_metacell_count
    else:
        random_metacell_count = None

    params = MetacellParams(
        target_size=max(1, args.target_metacell_size),
        min_size=max(1, args.min_metacell_size),
        max_size=max(args.min_metacell_size, args.max_metacell_size),
        preserve_small=args.preserve_small_clusters,
        aggregation=args.aggregation,
        graph_algorithm=args.graph_algorithm,
        resolution=args.resolution,
        resolution_steps=max(1, args.resolution_steps),
        resolution_tolerance=max(0.01, args.resolution_tolerance),
        n_neighbors=max(5, args.n_neighbors),
        neighbor_method=args.neighbor_method,
        neighbor_metric=args.neighbor_metric,
        n_top_genes=max(200, args.n_top_genes),
        n_pcs=max(10, args.n_pcs),
        pca_svd_solver=args.pca_svd_solver,
        hvg_layer=args.hvg_layer or args.expression_layer,
        expression_layer=args.expression_layer,
        use_raw=args.use_raw,
        random_state=args.random_state,
        random_metacell_count=random_metacell_count,
        random_cells_per_metacell=args.random_cells_per_metacell,
        random_sampling_with_replacement=args.random_sampling_with_replacement,
    )

    logging.info("Generating metacells (mode=%s, boundaries=%s).", args.mode, boundary_columns or "none")
    metacell_adata, membership_df = assemble_metacells(
        restricted,
        params=params,
        boundary_columns=boundary_columns,
        sample_column=args.sample_column,
        cell_type_column=args.cell_type_column,
    )

    metacell_adata.uns["metacell_membership"] = membership_df
    metacell_adata.uns["metacell_parameters"] = json.loads(json.dumps(vars(args), default=str))

    logging.info("Writing metacell AnnData to %s", args.output)
    metacell_adata.write(args.output, compression="gzip")
    logging.info("Metacell construction finished with %d metacells.", metacell_adata.n_obs)


if __name__ == "__main__":
    main()
