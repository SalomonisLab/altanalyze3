"""
Fast approximate UMAP placement using reference AnnData embeddings.

This streamlined module maps query cells onto a reference UMAP by sampling
reference barcodes from matching clusters and applying a small random offset.
It logs key dataset statistics and performs strict cluster compatibility
checks before proceeding.

CLI usage (excerpt):

    python -m altanalyze3.components.visualization.approximate_umap \
        --query QUERY.h5ad \
        --reference REF.h5ad \
        --query-cluster-key cluster_col \
        --outdir outputs/ \
        [--output-prefix custom_prefix] \
        [--save-updated-h5ad] \
        [--output-h5ad projected_query.h5ad] \
        [--output-pdf custom.pdf]

The `--outdir` argument is required and acts as the base directory for all
outputs. Supplying `--output-h5ad` automatically enables the h5ad export and
treats the value as the filename within `--outdir`; otherwise the export is
skipped unless `--save-updated-h5ad` is explicitly provided.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union
import sys

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import itertools
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.colors import to_hex
from matplotlib.lines import Line2D
from pandas.api.types import is_categorical_dtype

__all__ = ["ApproximateUMAPResult", "approximate_umap", "main"]


AnnDataLike = Union[str, Path, ad.AnnData]


def _str_to_bool(value: Optional[str]) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return True
    text = str(value).strip().lower()
    truthy = {"true", "t", "yes", "y", "1", "on"}
    falsy = {"false", "f", "no", "n", "0", "off"}
    if text in truthy:
        return True
    if text in falsy:
        return False
    raise argparse.ArgumentTypeError(f"Expected a boolean value, got '{value}'.")


@dataclass
class ApproximateUMAPResult:
    query_adata: ad.AnnData
    coordinates: pd.DataFrame
    reference_choices: Dict[str, List[str]]
    cluster_key: str
    cluster_order: List[str]

    def write_text_outputs(self, prefix: Union[str, Path]) -> Dict[str, str]:
        prefix = str(prefix)
        coords_path = f"{prefix}-approximate-umap.tsv"
        self.coordinates.to_csv(coords_path, sep="\t")

        obs_path = f"{prefix}-augmented.tsv"
        obs_df = self.query_adata.obs[[self.cluster_key]].copy()
        obs_df.insert(0, "barcode", self.query_adata.obs_names)
        obs_df.rename(columns={self.cluster_key: "cluster"}, inplace=True)
        obs_df.to_csv(obs_path, sep="\t", index=False)

        placeholder_matrix = pd.DataFrame(
            [
                ["UID"] + self.coordinates.index.tolist(),
                ["UMAP1"] + self.coordinates["umap_0"].astype(str).tolist(),
                ["UMAP2"] + self.coordinates["umap_1"].astype(str).tolist(),
            ]
        )
        placeholder_path = f"{prefix}-placeholder-exp.tsv"
        placeholder_matrix.to_csv(placeholder_path, sep="\t", header=False, index=False)

        return {
            "coordinates": coords_path,
            "augmented": obs_path,
            "placeholder_expression": placeholder_path,
        }

    def write_comparison_pdf(
        self,
        reference_adata: ad.AnnData,
        *,
        umap_key: str,
        reference_cluster_key: str,
        query_cluster_key: str,
        output_path: Union[str, Path],
    ) -> tuple[str, str]:
        output_path = Path(output_path)
        annotated, plain = _save_comparison_pdf(
            reference_adata,
            self.query_adata,
            umap_key=umap_key,
            reference_cluster_key=reference_cluster_key,
            query_cluster_key=query_cluster_key,
            lineage_order=self.cluster_order,
            destination=output_path,
        )
        return annotated, plain


def approximate_umap(
    query: AnnDataLike,
    reference: AnnDataLike,
    *,
    query_cluster_key: str,
    reference_cluster_key: Optional[str] = None,
    umap_key: str = "X_umap",
    jitter: float = 0.05,
    num_reference_cells: int = 1,
    random_state: Optional[int] = None,
) -> ApproximateUMAPResult:
    query_adata = _load_adata(query)
    reference_adata = _load_adata(reference)
    ref_cluster_key = reference_cluster_key or query_cluster_key

    _validate_inputs(
        query_adata,
        reference_adata,
        query_cluster_key=query_cluster_key,
        reference_cluster_key=ref_cluster_key,
        umap_key=umap_key,
        jitter=jitter,
        num_reference_cells=num_reference_cells,
    )

    ref_clusters = reference_adata.obs[ref_cluster_key].astype(str)
    qry_clusters = query_adata.obs[query_cluster_key].astype(str)

    ref_default_order = _determine_cluster_order(ref_clusters)
    cluster_order = list(ref_default_order)
    query_default_order = _determine_cluster_order(qry_clusters)

    query_lineage = _extract_lineage_order(query_adata, query_cluster_key)
    if query_lineage:
        print(
            f"[info] Detected lineage_order in query uns['lineage_order'] (n={len(query_lineage)}). "
            "This will determine category ordering."
        )
        print(f"[debug] Raw query lineage_order list: {list(query_lineage)}")
        sanitized = [str(item) for item in query_lineage if str(item) in query_default_order]
        missing_from_query_lineage = [cluster for cluster in query_default_order if cluster not in sanitized]
        if missing_from_query_lineage:
            print(
                "[warn] lineage_order in query is missing the following clusters; they will be appended "
                "using the observed order: " + ", ".join(missing_from_query_lineage)
            )
            sanitized.extend(missing_from_query_lineage)
        ref_cluster_set = set(ref_default_order)
        extras_in_query_lineage = [cluster for cluster in sanitized if cluster not in ref_cluster_set]
        if extras_in_query_lineage:
            print(
                "[warn] Removing clusters from query lineage_order not present in the reference: "
                + ", ".join(extras_in_query_lineage)
            )
        cluster_order = [cluster for cluster in sanitized if cluster in ref_cluster_set]
        if not cluster_order:
            print(
                "[warn] No clusters remained after aligning query lineage_order with the reference; "
                "reverting to reference ordering."
            )
            cluster_order = list(ref_default_order)
        print(f"[info] Using lineage_order from query AnnData (n={len(cluster_order)}).")
    else:
        print(
            "[warn] Query AnnData missing lineage_order; falling back to reference-derived ordering."
        )
        lineage = _extract_lineage_order(reference_adata, ref_cluster_key)
        if lineage:
            preview = ", ".join(lineage[:10])
            if len(lineage) > 10:
                preview += ", …"
            print(
                f"[info] Detected lineage_order in reference uns['lineage_order'] (n={len(lineage)}): {preview}"
            )
            print(f"[debug] Raw lineage_order list (preserving reference order): {list(lineage)}")
            lineage_order = [str(item) for item in lineage if str(item) in cluster_order]
            print(f"[info] Query clusters matching lineage_order: {len(lineage_order)}")
            missing_from_lineage = [cluster for cluster in cluster_order if cluster not in lineage_order]
            if missing_from_lineage:
                print(
                    "[warn] The following clusters were not found in lineage_order and will be appended: "
                    + ", ".join(missing_from_lineage)
                )
            cluster_order = lineage_order + missing_from_lineage
    if cluster_order:
        sample_preview = ", ".join(cluster_order[:15])
        if len(cluster_order) > 15:
            sample_preview += ", …"
        print(f"[debug] Final cluster_order applied to query AnnData (n={len(cluster_order)}): {sample_preview}")
    query_labels = pd.unique(qry_clusters)

    missing_clusters = sorted(set(query_labels) - set(cluster_order))
    if missing_clusters:
        raise ValueError(
            "The following query clusters are absent from the reference: "
            + ", ".join(missing_clusters)
        )

    logging.info(
        "Reference: %d cells across %d clusters.",
        reference_adata.n_obs,
        len(cluster_order),
    )
    print(f"Reference dataset: {reference_adata.n_obs} cells across {len(cluster_order)} clusters.")

    logging.info(
        "Query: %d cells across %d clusters.",
        query_adata.n_obs,
        len(query_labels),
    )
    print(f"Query dataset: {query_adata.n_obs} cells across {len(query_labels)} clusters.")

    unused_clusters = sorted(set(cluster_order) - set(query_labels))
    if unused_clusters:
        logging.info(
            "Skipping %d reference clusters not present in the query: %s",
            len(unused_clusters),
            ", ".join(unused_clusters),
        )
        print(
            "Skipping reference-only clusters:",
            ", ".join(unused_clusters),
        )

    coords_df = _embedding_to_df(reference_adata, umap_key)
    cells_by_cluster: Dict[str, List[str]] = {}
    cluster_bounds: Dict[str, Tuple[float, float, float, float]] = {}
    coords_by_cluster = coords_df.copy()
    coords_by_cluster["cluster"] = ref_clusters.loc[coords_df.index].astype(str)
    for cluster in cluster_order:
        cluster_key = str(cluster)
        cluster_coords = coords_by_cluster[coords_by_cluster["cluster"] == cluster_key][["umap_0", "umap_1"]]
        if cluster_coords.empty:
            continue
        q1 = cluster_coords.quantile(0.25)
        q3 = cluster_coords.quantile(0.75)
        iqr = q3 - q1
        lower = q1 - 1.5 * iqr
        upper = q3 + 1.5 * iqr
        mask = cluster_coords.apply(
            lambda row: (lower["umap_0"] <= row["umap_0"] <= upper["umap_0"])
            and (lower["umap_1"] <= row["umap_1"] <= upper["umap_1"]),
            axis=1,
        )
        filtered = cluster_coords[mask]
        if filtered.empty:
            filtered = cluster_coords
        cells_by_cluster[cluster_key] = filtered.index.tolist()
        x_low, x_high = float(lower["umap_0"]), float(upper["umap_0"])
        y_low, y_high = float(lower["umap_1"]), float(upper["umap_1"])
        if not np.isfinite(x_low) or not np.isfinite(x_high) or x_high <= x_low:
            span = float(cluster_coords["umap_0"].max() - cluster_coords["umap_0"].min()) or 0.1
            center = float(cluster_coords["umap_0"].median())
            x_low, x_high = center - span / 2.0, center + span / 2.0
        if not np.isfinite(y_low) or not np.isfinite(y_high) or y_high <= y_low:
            span = float(cluster_coords["umap_1"].max() - cluster_coords["umap_1"].min()) or 0.1
            center = float(cluster_coords["umap_1"].median())
            y_low, y_high = center - span / 2.0, center + span / 2.0
        cluster_bounds[cluster_key] = (x_low, x_high, y_low, y_high)

    reference_counts = {cluster: len(indices) for cluster, indices in cells_by_cluster.items()}
    query_counts = Counter(qry_clusters.astype(str))

    rng = np.random.default_rng(random_state)
    records: List[Dict[str, float]] = []
    sampled: Dict[str, List[str]] = {}

    logging.info(
        "Assigning coordinates using up to %d reference cells per query cell; jitter=%.3f.",
        num_reference_cells,
        jitter,
    )
    print(
        f"Assigning coordinates using up to {num_reference_cells} reference cells per query cell "
        f"with jitter {jitter:.3f}."
    )

    for barcode, cluster in zip(query_adata.obs_names, qry_clusters):
        cluster = str(cluster)
        ref_cells = cells_by_cluster.get(cluster)
        if not ref_cells:
            raise ValueError(
                f"Reference cluster '{cluster}' has no cells and cannot be sampled."
            )

        reference_count = reference_counts.get(cluster, 0)
        query_count = query_counts.get(cluster, 0)
        sample_size = num_reference_cells
        if reference_count and query_count > 4 * reference_count:
            sample_size = max(2, num_reference_cells)

        replace = len(ref_cells) < sample_size
        sampled_cells = rng.choice(ref_cells, size=sample_size, replace=replace).tolist()

        if reference_count and query_count > 4 * reference_count and cluster in cluster_bounds:
            coords_subset = coords_df.loc[sampled_cells, ["umap_0", "umap_1"]].to_numpy()
            base = coords_subset.mean(axis=0)
            offset = rng.uniform(-jitter, jitter, size=2) if jitter else 0.0
            approx = base + offset
            x_low, x_high, y_low, y_high = cluster_bounds[cluster]
            approx[0] = float(np.clip(approx[0], x_low, x_high))
            approx[1] = float(np.clip(approx[1], y_low, y_high))
        else:
            base = coords_df.loc[sampled_cells[0], ["umap_0", "umap_1"]].to_numpy()
            offset = rng.uniform(-jitter, jitter, size=2) if jitter else 0.0
            approx = base + offset

        records.append({"barcode": barcode, "umap_0": float(approx[0]), "umap_1": float(approx[1])})
        sampled[barcode] = sampled_cells

    print(f"Processed {len(records)} query cells.")

    coordinates = pd.DataFrame.from_records(records).set_index("barcode")

    query_with_umap = query_adata.copy()
    query_with_umap.obsm[umap_key] = coordinates.loc[query_with_umap.obs_names].to_numpy()
    if cluster_order:
        query_with_umap.obs[query_cluster_key] = pd.Categorical(
            query_with_umap.obs[query_cluster_key],
            categories=cluster_order,
            ordered=True,
        )
        query_with_umap.uns["lineage_order"] = list(cluster_order)
        preview = ", ".join(cluster_order[:15])
        if len(cluster_order) > 15:
            preview += ", …"
        print(f"[debug] Stored query lineage_order (uns['lineage_order']) with n={len(cluster_order)}: {preview}")
        print(f"[info] Stored lineage_order on query AnnData (n={len(cluster_order)}).")

    return ApproximateUMAPResult(
        query_adata=query_with_umap,
        coordinates=coordinates,
        reference_choices=sampled,
        cluster_key=query_cluster_key,
        cluster_order=list(cluster_order),
    )


def _load_adata(obj: AnnDataLike) -> ad.AnnData:
    if isinstance(obj, ad.AnnData):
        return obj
    return ad.read_h5ad(str(obj))


def _create_minimal_adata(barcodes: Sequence[str]) -> ad.AnnData:
    if not barcodes:
        raise ValueError("Cannot build AnnData from an empty barcode list.")
    X = np.zeros((len(barcodes), 1), dtype=float)
    adata = ad.AnnData(X=X)
    adata.var_names = pd.Index(["__dummy__"])
    adata.obs_names = pd.Index([str(bc) for bc in barcodes], name="barcode")
    return adata


def _detect_column(columns: List[str], candidates: Sequence[str]) -> Optional[str]:
    lower_map = {col.lower(): col for col in columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def _load_reference_from_tsv(
    coords_path: Union[str, Path],
    clusters_path: Union[str, Path],
    *,
    umap_key: str,
    cluster_key: str,
) -> ad.AnnData:
    coords_df = pd.read_csv(coords_path, sep="\t")
    cluster_df = pd.read_csv(clusters_path, sep="\t")
    if coords_df.empty:
        raise ValueError("Reference coordinate TSV is empty.")
    if cluster_df.empty:
        raise ValueError("Reference cluster TSV is empty.")

    barcode_col = _detect_column(
        coords_df.columns.tolist(),
        ["barcode", "cellbarcode", "cell_barcode", "cell"],
    )
    if barcode_col is None:
        raise ValueError("Reference coordinate TSV must contain a barcode column.")
    umap1_col = _detect_column(coords_df.columns.tolist(), ["umap1", "umap_1", "umap0", "umap_0", "umapx"])
    umap2_col = _detect_column(coords_df.columns.tolist(), ["umap2", "umap_2", "umapy", "umap1", "umap_1"])
    if umap1_col is None or umap2_col is None:
        raise ValueError(
            "Reference coordinate TSV must contain columns for UMAP1/UMAP2 (e.g. 'UMAP1' and 'UMAP2')."
        )

    cluster_barcode_col = _detect_column(
        cluster_df.columns.tolist(),
        ["barcode", "cellbarcode", "cell_barcode", "cell"],
    )
    if cluster_barcode_col is None:
        raise ValueError("Reference cluster TSV must contain a barcode column.")
    cluster_label_col = cluster_key if cluster_key in cluster_df.columns else _detect_column(
        cluster_df.columns.tolist(),
        ["cluster", "population", "label"],
    )
    if cluster_label_col is None:
        raise ValueError(
            f"Could not locate a cluster column in {clusters_path}; "
            "expected either the provided cluster key or a column named 'cluster'/'population'."
        )

    coords_df = coords_df.rename(
        columns={
            barcode_col: "barcode",
            umap1_col: "umap_0",
            umap2_col: "umap_1",
        }
    )
    cluster_df = cluster_df.rename(
        columns={
            cluster_barcode_col: "barcode",
            cluster_label_col: "cluster",
        }
    )
    merged = pd.merge(coords_df, cluster_df[["barcode", "cluster"]], on="barcode", how="inner")
    if merged.empty:
        raise ValueError("No overlapping barcodes found between coordinate and cluster TSV files.")

    adata = _create_minimal_adata(merged["barcode"].tolist())
    adata.obs[cluster_key] = merged["cluster"].astype(str).values
    adata.obsm[umap_key] = merged[["umap_0", "umap_1"]].to_numpy(dtype=float)
    return adata


def _load_query_from_tsv(
    clusters_path: Union[str, Path],
    *,
    cluster_key: str,
) -> ad.AnnData:
    df = pd.read_csv(clusters_path, sep="\t")
    if df.empty:
        raise ValueError("Query cluster TSV is empty.")
    barcode_col = _detect_column(
        df.columns.tolist(),
        ["barcode", "cellbarcode", "cell_barcode", "cell", "CellBarcode"],
    )
    if barcode_col is None:
        raise ValueError("Query cluster TSV must contain a barcode column.")
    cluster_col = cluster_key if cluster_key in df.columns else _detect_column(
        df.columns.tolist(),
        ["cluster", "population", "label", "alignedpopulation", "reference"],
    )
    if cluster_col is None:
        raise ValueError(
            f"Could not locate a cluster column in {clusters_path}; "
            "expected either '{cluster_key}' or a column named 'cluster'/'population'."
        )
    df = df.rename(columns={barcode_col: "barcode", cluster_col: "cluster"})
    df = df.dropna(subset=["barcode", "cluster"])
    adata = _create_minimal_adata(df["barcode"].astype(str).tolist())
    adata.obs[cluster_key] = df["cluster"].astype(str).values
    return adata


def _validate_inputs(
    query: ad.AnnData,
    reference: ad.AnnData,
    *,
    query_cluster_key: str,
    reference_cluster_key: str,
    umap_key: str,
    jitter: float,
    num_reference_cells: int,
) -> None:
    if query_cluster_key not in query.obs:
        raise KeyError(f"Query AnnData is missing '{query_cluster_key}' in .obs")
    if reference_cluster_key not in reference.obs:
        raise KeyError(f"Reference AnnData is missing '{reference_cluster_key}' in .obs")
    if umap_key not in reference.obsm:
        raise KeyError(f"Reference AnnData is missing '{umap_key}' in .obsm")

    if query.n_obs == 0:
        raise ValueError("Query AnnData contains no cells.")
    if reference.n_obs == 0:
        raise ValueError("Reference AnnData contains no cells.")

    if jitter < 0:
        raise ValueError("Parameter 'jitter' must be non-negative.")
    if num_reference_cells <= 0:
        raise ValueError("Parameter 'num_reference_cells' must be a positive integer.")


def _determine_cluster_order(series: pd.Series) -> List[str]:
    if is_categorical_dtype(series):
        return [str(cat) for cat in series.cat.categories]
    return [str(value) for value in pd.unique(series)]


def _embedding_to_df(adata: ad.AnnData, umap_key: str) -> pd.DataFrame:
    embedding = adata.obsm[umap_key]
    if isinstance(embedding, pd.DataFrame):
        arr = embedding.iloc[:, :2].to_numpy()
    else:
        arr = np.asarray(embedding)
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError(f"Embedding '{umap_key}' must be at least 2-dimensional.")
    return pd.DataFrame(arr, index=adata.obs_names, columns=["umap_0", "umap_1"])


def _save_comparison_pdf(
    reference: ad.AnnData,
    query: ad.AnnData,
    *,
    umap_key: str,
    reference_cluster_key: str,
    query_cluster_key: str,
    lineage_order: Optional[Sequence[str]] = None,
    destination: Path,
) -> tuple[str, str]:
    if umap_key not in reference.obsm:
        raise KeyError(f"Reference AnnData missing '{umap_key}' in .obsm for plotting.")
    if umap_key not in query.obsm:
        raise KeyError(f"Query AnnData missing '{umap_key}' in .obsm for plotting.")

    ref_coords = _embedding_to_df(reference, umap_key).to_numpy()
    qry_coords = _embedding_to_df(query, umap_key).to_numpy()

    ref_labels = reference.obs[reference_cluster_key].astype(str)
    qry_labels = query.obs[query_cluster_key].astype(str)

    categories: List[str] = []
    for lab in pd.unique(ref_labels.tolist() + qry_labels.tolist()):
        if lab not in categories:
            categories.append(lab)
    if lineage_order:
        lineage_order = [str(item) for item in lineage_order]
        lineage_set = set(lineage_order)
        ordered = [label for label in lineage_order if label in categories]
        missing = [label for label in categories if label not in lineage_set]
        print(f"[info] Legend lineage match count: {len(ordered)} / {len(categories)}")
        if missing:
            print(
                "[info] Legend clusters absent from lineage_order: "
                + ", ".join(missing)
            )
        categories = ordered + missing
        if ordered:
            print("[info] Legend will follow lineage_order for matching clusters.")
        else:
            print("[warn] No legend categories matched lineage_order; falling back to detected order.")
        compare_count = min(len(lineage_order), len(categories))
        if compare_count:
            preview_pairs = [
                f"{lineage_order[i]} => {categories[i]}" for i in range(min(compare_count, 15))
            ]
            if compare_count > 15:
                preview_pairs.append("…")
            print(f"[debug] Lineage vs legend alignment preview (first {min(compare_count, 15)} entries): {preview_pairs}")
    if categories:
        legend_preview = ", ".join(categories[:15])
        if len(categories) > 15:
            legend_preview += ", …"
        print(f"[debug] Legend categories passed to Matplotlib (n={len(categories)}): {legend_preview}")

    palette = _resolve_colors(categories, reference, query, reference_cluster_key, query_cluster_key)

    plt.rcParams.update(
        {
            "backend": "Agg",
            "axes.linewidth": 0.5,
            "pdf.fonttype": 42,
            "font.family": "sans-serif",
            "font.sans-serif": "Arial",
            "figure.facecolor": "white",
        }
    )

    legend_handles = [
        Line2D([0], [0], marker="o", color="w", label=cat, markerfacecolor=palette[cat], markersize=3)
        for cat in categories
    ]

    destination.parent.mkdir(parents=True, exist_ok=True)

    fig1, axes1 = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    _scatter(axes1[0], ref_coords, ref_labels, palette, "Reference UMAP", annotate=True)
    _scatter(axes1[1], qry_coords, qry_labels, palette, "Query Approximate UMAP", annotate=True)
    fig1.legend(
        handles=legend_handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        borderaxespad=0.0,
        handlelength=0.6,
        labelspacing=0.3,
        columnspacing=0.6,
        fontsize=6,
    )
    annotated_path = destination
    fig1.savefig(annotated_path, dpi=300, bbox_inches="tight")
    plt.close(fig1)

    plain_path = destination.with_name(destination.stem + "-no-labels" + destination.suffix)
    fig2, axes2 = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    _scatter(axes2[0], ref_coords, ref_labels, palette, "Reference UMAP", annotate=False)
    _scatter(axes2[1], qry_coords, qry_labels, palette, "Query Approximate UMAP", annotate=False)
    fig2.legend(
        handles=legend_handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        borderaxespad=0.0,
        handlelength=0.6,
        labelspacing=0.3,
        columnspacing=0.6,
        fontsize=6,
    )
    fig2.savefig(plain_path, dpi=300, bbox_inches="tight")
    plt.close(fig2)

    return str(annotated_path), str(plain_path)

def _scatter(
    axis: matplotlib.axes.Axes,
    coords: np.ndarray,
    labels: pd.Series,
    palette: Mapping[str, str],
    title: str,
    *,
    annotate: bool,
) -> None:
    colors = [palette.get(label, "#999999") for label in labels]
    point_size = _compute_point_size(coords.shape[0])
    axis.scatter(coords[:, 0], coords[:, 1], c=colors, s=point_size, linewidths=0, alpha=0.85)
    axis.set_title(title, fontsize=12)
    axis.set_xlabel("UMAP1", fontsize=10)
    axis.set_ylabel("UMAP2", fontsize=10)
    axis.tick_params(labelsize=8)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    if annotate:
        _annotate_clusters(axis, coords, labels)


def _extract_lineage_order(adata: ad.AnnData, cluster_key: str) -> Optional[List[str]]:
    if "lineage_order" not in getattr(adata, "uns", {}):
        return None
    raw_order = adata.uns.get("lineage_order")
    try:
        order = [str(item) for item in raw_order if isinstance(item, (str, bytes)) or hasattr(item, "__str__")]
    except Exception:  # pragma: no cover - defensive
        return None
    seen = set()
    deduped = []
    for item in order:
        if item not in seen:
            seen.add(item)
            deduped.append(item)
    return deduped or None


def _resolve_colors(
    categories: Sequence[str],
    reference: ad.AnnData,
    query: ad.AnnData,
    reference_cluster_key: str,
    query_cluster_key: str,
) -> Dict[str, str]:
    candidates = [
        (query, f"{query_cluster_key}_colors"),
        (reference, f"{query_cluster_key}_colors"),
        (reference, f"{reference_cluster_key}_colors"),
    ]
    palette: Optional[List[str]] = None
    for adata_obj, key in candidates:
        if key in adata_obj.uns:
            colors = adata_obj.uns[key]
            if isinstance(colors, Iterable) and not isinstance(colors, (str, bytes)):
                palette = list(dict.fromkeys(colors))  # preserve order, drop duplicates
                break

    if not palette or len(palette) < len(categories):
        base_palette = [
            to_hex(matplotlib.colors.hsv_to_rgb([h, 0.65, 0.95]))
            for h in np.linspace(0, 1, len(categories), endpoint=False)
        ]
        palette = base_palette
    else:
        sanitized: List[str] = []
        for color in palette:
            try:
                sanitized.append(to_hex(mcolors.to_rgb(color)))
            except ValueError:
                continue
        if len(sanitized) < len(categories):
            extra_needed = len(categories) - len(sanitized)
            cmap = plt.get_cmap("tab20", max(len(categories), 1))
            generated = [to_hex(cmap(i)) for i in range(len(categories))]
            sanitized.extend(generated[:extra_needed])
        palette = sanitized

    rng = np.random.default_rng(123)
    rng.shuffle(palette)

    return {str(cat): palette[i % len(palette)] for i, cat in enumerate(categories)}


def _compute_point_size(num_points: int) -> float:
    base = 0.6  # one-tenth of the previous default size (6)
    scale = (1000.0 / max(num_points, 1)) ** 0.5
    size = base * scale
    return float(np.clip(size, 0.1, 2.0))


def _annotate_clusters(axis: matplotlib.axes.Axes, coords: np.ndarray, labels: pd.Series) -> None:
    unique_labels = pd.unique(labels)
    if len(unique_labels) == 0:
        return

    x_span = float(np.ptp(coords[:, 0])) or 1.0
    y_span = float(np.ptp(coords[:, 1])) or 1.0

    for idx, cluster in enumerate(unique_labels):
        mask = labels == cluster
        if not np.any(mask):
            continue
        centroid = coords[mask].mean(axis=0)
        angle = 2 * np.pi * idx / max(len(unique_labels), 1)
        offset = np.array([np.cos(angle) * x_span * 0.06, np.sin(angle) * y_span * 0.06])
        text_pos = centroid + offset

        axis.annotate(
            str(cluster),
            xy=centroid,
            xytext=text_pos,
            textcoords="data",
            ha="center",
            va="center",
            fontsize=7,
            arrowprops={"arrowstyle": "-", "color": "#555555", "lw": 0.6},
            zorder=5,
        )


def _parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Approximate UMAP coordinates for a query AnnData using a reference embedding."
    )
    parser.add_argument("--query", help="Path to the query .h5ad file.")
    parser.add_argument("--reference", help="Path to the reference .h5ad file.")
    parser.add_argument(
        "--query-clusters-tsv",
        help="cellHarmony-lite TSV (CellBarcode + cluster column) to build the query AnnData without .h5ad input.",
    )
    parser.add_argument(
        "--reference-coords-tsv",
        help="TSV with reference barcodes and UMAP coordinates (columns: barcode, UMAP1, UMAP2).",
    )
    parser.add_argument(
        "--reference-clusters-tsv",
        help="TSV with reference barcodes and clusters (columns: barcode, cluster).",
    )
    parser.add_argument("--query-cluster-key", required=True, help="Cluster column in query.obs.")
    parser.add_argument(
        "--reference-cluster-key",
        help="Cluster column in reference.obs. Defaults to the query cluster key.",
    )
    parser.add_argument(
        "--umap-key",
        default="X_umap",
        help="UMAP embedding key in reference.obsm (default: %(default)s).",
    )
    parser.add_argument(
        "--jitter",
        type=float,
        default=0.05,
        help="Uniform random offset applied to sampled coordinates (default: %(default)s).",
    )
    parser.add_argument(
        "--num-reference-cells",
        type=int,
        default=1,
        help="Number of reference barcodes sampled per query cell (default: %(default)s).",
    )
    parser.add_argument("--random-state", type=int, help="Seed for reproducible sampling.")
    parser.add_argument(
        "--output-prefix",
        help="Prefix for outputs (used within --outdir). Defaults to <query>-approximate-umap.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory where outputs will be written.",
    )
    parser.add_argument(
        "--output-h5ad",
        help="Filename for the augmented query AnnData (stored in --outdir; defaults to <prefix>.h5ad).",
    )
    parser.add_argument(
        "--output-pdf",
        help="Filename or path for the comparison PDF (default: <prefix>-comparison.pdf within --outdir).",
    )
    parser.add_argument(
        "--save-updated-h5ad",
        nargs="?",
        const="true",
        default=False,
        type=_str_to_bool,
        help="Write the augmented query AnnData (automatically enabled when --output-h5ad is supplied).",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable INFO logging.")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parse_args(argv)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING)

    if not args.query and not args.query_clusters_tsv:
        print("Error: Provide either --query (h5ad) or --query-clusters-tsv.", file=sys.stderr)
        return 2
    if not args.reference and not (args.reference_coords_tsv and args.reference_clusters_tsv):
        print(
            "Error: Provide either --reference (h5ad) or both --reference-coords-tsv and --reference-clusters-tsv.",
            file=sys.stderr,
        )
        return 2

    if args.reference:
        reference_source: AnnDataLike = ad.read_h5ad(args.reference)
    else:
        ref_cluster_key = args.reference_cluster_key or args.query_cluster_key
        reference_source = _load_reference_from_tsv(
            args.reference_coords_tsv,
            args.reference_clusters_tsv,
            umap_key=args.umap_key,
            cluster_key=ref_cluster_key,
        )

    if args.query:
        query_source: AnnDataLike = args.query
    else:
        query_source = _load_query_from_tsv(
            args.query_clusters_tsv,
            cluster_key=args.query_cluster_key,
        )

    result = approximate_umap(
        query=query_source,
        reference=reference_source,
        query_cluster_key=args.query_cluster_key,
        reference_cluster_key=args.reference_cluster_key,
        umap_key=args.umap_key,
        jitter=args.jitter,
        num_reference_cells=args.num_reference_cells,
        random_state=args.random_state,
    )

    query_path = Path(args.query)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    prefix_name = args.output_prefix or f"{query_path.stem}-approximate-umap"
    prefix = out_dir / prefix_name

    if args.output_pdf:
        output_pdf_candidate = Path(args.output_pdf)
        output_pdf = (
            output_pdf_candidate
            if output_pdf_candidate.is_absolute()
            else out_dir / output_pdf_candidate
        )
    else:
        output_pdf = out_dir / f"{prefix_name}-comparison.pdf"
    output_pdf.parent.mkdir(parents=True, exist_ok=True)

    full_command = " ".join(sys.argv)
    save_h5ad = bool(args.save_updated_h5ad)
    output_h5ad: Optional[Path] = None

    default_h5ad_name = f"{prefix_name}.h5ad"
    if args.output_h5ad:
        output_h5ad_candidate = Path(args.output_h5ad)
        output_h5ad = (
            output_h5ad_candidate
            if output_h5ad_candidate.is_absolute()
            else out_dir / output_h5ad_candidate
        )
        if not save_h5ad:
            print("[info] --output-h5ad supplied; enabling --save-updated-h5ad.")
        save_h5ad = True
    elif save_h5ad:
        output_h5ad = out_dir / default_h5ad_name

    if save_h5ad:
        output_h5ad.parent.mkdir(parents=True, exist_ok=True)
        logging.info("Writing augmented AnnData to %s", output_h5ad)
        print(f"[info] Saving augmented AnnData to {output_h5ad} (command: {full_command})")
        result.query_adata.write(output_h5ad, compression="gzip")
    else:
        logging.info("Skipping augmented AnnData export (enable with --save-updated-h5ad).")
        print(f"[info] Skipping augmented AnnData export (enable with --save-updated-h5ad). Command: {full_command}")

    logging.info("Writing legacy-style exports with prefix %s", prefix)
    result.write_text_outputs(prefix)

    logging.info("Writing comparison PDFs to %s (annotated) and companion without labels", output_pdf)
    reference_for_plot = reference_source if isinstance(reference_source, ad.AnnData) else _load_adata(reference_source)
    annotated_pdf, plain_pdf = result.write_comparison_pdf(
        reference_adata=reference_for_plot,
        umap_key=args.umap_key,
        reference_cluster_key=args.reference_cluster_key or args.query_cluster_key,
        query_cluster_key=args.query_cluster_key,
        output_path=output_pdf,
    )
    logging.info("Annotated UMAP PDF: %s", annotated_pdf)
    logging.info("Label-free UMAP PDF: %s", plain_pdf)

    logging.info("Approximate UMAP mapping completed.")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
