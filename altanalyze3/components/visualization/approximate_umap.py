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
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import anndata as ad
import h5py
import matplotlib
import matplotlib.pyplot as plt
import itertools
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.colors import to_hex
from matplotlib.lines import Line2D
from pandas.api.types import is_categorical_dtype
import scipy.sparse as sp

__all__ = ["ApproximateUMAPResult", "approximate_umap", "main"]


AnnDataLike = Union[str, Path, ad.AnnData]


_PAIRED_COLOR_STOPS = [
    (0.0, (0.6509804129600525, 0.8078431487083435, 0.8901960849761963)),
    (0.09090909090909091, (0.12156862765550613, 0.47058823704719543, 0.7058823704719543)),
    (0.18181818181818182, (0.6980392336845398, 0.8745098114013672, 0.5411764979362488)),
    (0.2727272727272727, (0.20000000298023224, 0.6274510025978088, 0.1725490242242813)),
    (0.36363636363636365, (0.9843137264251709, 0.6039215922355652, 0.6000000238418579)),
    (0.45454545454545453, (0.8901960849761963, 0.10196078568696976, 0.10980392247438431)),
    (0.5454545454545454, (0.9921568632125854, 0.7490196228027344, 0.43529412150382996)),
    (0.6363636363636364, (1.0, 0.49803921580314636, 0.0)),
    (0.7272727272727273, (0.7921568751335144, 0.6980392336845398, 0.8392156958580017)),
    (0.8181818181818182, (0.4156862795352936, 0.239215686917305, 0.6039215922355652)),
    (0.9090909090909091, (1.0, 1.0, 0.6000000238418579)),
    (1.0, (0.6941176652908325, 0.3490196168422699, 0.1568627506494522)),
]


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


def _clean_string(value: object) -> str:
    if value is None:
        return ""
    if pd.isna(value):
        return ""
    return str(value).strip()


def _flatten_expr(values) -> np.ndarray:
    if sp.issparse(values):
        return np.asarray(values.todense()).ravel()
    return np.asarray(values).ravel()


def _interpolate_paired_color(value: float) -> tuple[float, float, float]:
    clipped = max(0.0, min(1.0, float(value)))
    for index in range(1, len(_PAIRED_COLOR_STOPS)):
        right_pos, right_rgb = _PAIRED_COLOR_STOPS[index]
        left_pos, left_rgb = _PAIRED_COLOR_STOPS[index - 1]
        if clipped <= right_pos or index == len(_PAIRED_COLOR_STOPS) - 1:
            span = max(right_pos - left_pos, 1e-9)
            ratio = max(0.0, min(1.0, (clipped - left_pos) / span))
            return tuple(
                left_rgb[channel] + (right_rgb[channel] - left_rgb[channel]) * ratio
                for channel in range(3)
            )
    return _PAIRED_COLOR_STOPS[-1][1]


def _custom_shuffle_indices(indices: list[int]) -> list[int]:
    shuffled: list[int] = []
    for index, value in enumerate(indices):
        if value not in shuffled:
            shuffled.append(value)
        from_end = indices[len(indices) - 1 - index] if indices else value
        if from_end not in shuffled:
            shuffled.append(from_end)
        middle_index = int((index + len(indices)) / 2)
        from_middle = indices[middle_index] if middle_index < len(indices) else indices[-1]
        if from_middle not in shuffled:
            shuffled.append(from_middle)
    return shuffled


def _seeded_shuffle_js(items: list[int], seed: int = 0) -> list[int]:
    state = (int(seed) & 0xFFFFFFFF) or 1

    def _next_rand() -> float:
        nonlocal state
        state = (1664525 * state + 1013904223) & 0xFFFFFFFF
        return state / 4294967296.0

    values = list(items)
    for index in range(len(values) - 1, 0, -1):
        swap_index = int(np.floor(_next_rand() * (index + 1)))
        values[index], values[swap_index] = values[swap_index], values[index]
    return values


def _build_preview_palette(populations: list[str]) -> dict[str, str]:
    ordered = list(populations)
    if len(ordered) <= 4:
        base = ["#ff0000", "#0000ff", "#ffff00", "#00aa00", "#ffffff", "#000000", "#ff00ff"]
        return {population: base[index % len(base)] for index, population in enumerate(ordered)}
    indices = _seeded_shuffle_js(_custom_shuffle_indices(list(range(len(ordered)))), 0)
    denominator = max(len(ordered) - 1, 1)
    colors = [to_hex(_interpolate_paired_color(index / denominator)) for index in indices]
    return {population: colors[index] for index, population in enumerate(ordered)}


def _build_stable_umap_population_order(reference_labels: Sequence[str], query_labels: Sequence[str]) -> list[str]:
    ordered: list[str] = []
    seen: set[str] = set()
    for population in itertools.chain(reference_labels, query_labels):
        label = _clean_string(population)
        if not label or label in seen:
            continue
        seen.add(label)
        ordered.append(label)
    return ordered


def _finite_extent(values: Sequence[float], fallback_min: float = 0.0, fallback_max: float = 1.0) -> tuple[float, float]:
    finite = [float(value) for value in values if np.isfinite(float(value))]
    if not finite:
        return fallback_min, fallback_max
    return float(min(finite)), float(max(finite))


def _median(values: Sequence[float]) -> float:
    arr = np.asarray(list(values), dtype=float)
    if arr.size == 0:
        return 0.0
    return float(np.median(arr))


def _build_population_centroids(coords: np.ndarray, labels: Sequence[str]) -> list[dict[str, float | str]]:
    grouped: dict[str, dict[str, list[float]]] = {}
    for (x, y), population in zip(coords, labels):
        label = _clean_string(population)
        if not label or not np.isfinite(x) or not np.isfinite(y):
            continue
        entry = grouped.setdefault(label, {"xs": [], "ys": []})
        entry["xs"].append(float(x))
        entry["ys"].append(float(y))
    return [
        {"population": population, "x": _median(values["xs"]), "y": _median(values["ys"])}
        for population, values in grouped.items()
    ]


def _relax_umap_labels(
    labels: list[dict[str, float | str]],
    points: np.ndarray,
    *,
    plot_width: int = 900,
    plot_height: int = 520,
) -> list[dict[str, float | str]]:
    if len(labels) <= 1:
        return labels
    x_values = [float(x) for x in points[:, 0] if np.isfinite(x)]
    y_values = [float(y) for y in points[:, 1] if np.isfinite(y)]
    for label in labels:
        x = float(label.get("x", np.nan))
        y = float(label.get("y", np.nan))
        if np.isfinite(x):
            x_values.append(x)
        if np.isfinite(y):
            y_values.append(y)
    x_min, x_max = _finite_extent(x_values, 0.0, 1.0)
    y_min, y_max = _finite_extent(y_values, 0.0, 1.0)
    x_range = max(x_max - x_min, 1.0)
    y_range = max(y_max - y_min, 1.0)
    effective_width = max(plot_width - 80, 480)
    pad_x = (8 / effective_width) * x_range
    pad_y = (5 / plot_height) * y_range
    max_dx = x_range * 0.03
    max_dy = y_range * 0.03
    relaxed: list[dict[str, float | str]] = []
    for label in labels:
        text = _clean_string(label.get("population"))
        width_px = max(42, min(170, len(text) * 6.2))
        height_px = 18
        x = float(label.get("x", 0.0))
        y = float(label.get("y", 0.0))
        relaxed.append(
            {
                **label,
                "x": x,
                "y": y,
                "anchor_x": x,
                "anchor_y": y,
                "half_width": (width_px / effective_width) * x_range * 0.5,
                "half_height": (height_px / plot_height) * y_range * 0.5,
            }
        )
    for _ in range(110):
        for index, current in enumerate(relaxed):
            for compare_index in range(index + 1, len(relaxed)):
                other = relaxed[compare_index]
                dx = float(other["x"]) - float(current["x"])
                dy = float(other["y"]) - float(current["y"])
                overlap_x = float(current["half_width"]) + float(other["half_width"]) + pad_x - abs(dx)
                overlap_y = float(current["half_height"]) + float(other["half_height"]) + pad_y - abs(dy)
                if overlap_x <= 0 or overlap_y <= 0:
                    continue
                shift_x = overlap_x * 0.5 * (1 if dx >= 0 else -1)
                shift_y = overlap_y * 0.5 * (1 if dy >= 0 else -1)
                current["x"] = float(current["x"]) - shift_x
                other["x"] = float(other["x"]) + shift_x
                current["y"] = float(current["y"]) - shift_y
                other["y"] = float(other["y"]) + shift_y
        for item in relaxed:
            item["x"] = float(np.clip(float(item["x"]), float(item["anchor_x"]) - max_dx, float(item["anchor_x"]) + max_dx))
            item["y"] = float(np.clip(float(item["y"]), float(item["anchor_y"]) - max_dy, float(item["anchor_y"]) + max_dy))
    return relaxed


def _square_umap_axes(axis: matplotlib.axes.Axes, coords: np.ndarray, padding_fraction: float = 0.06) -> None:
    if coords.size == 0:
        return
    xs = coords[:, 0]
    ys = coords[:, 1]
    finite_mask = np.isfinite(xs) & np.isfinite(ys)
    if not np.any(finite_mask):
        return
    xs = xs[finite_mask]
    ys = ys[finite_mask]
    min_x, max_x = float(xs.min()), float(xs.max())
    min_y, max_y = float(ys.min()), float(ys.max())
    x_span = max(max_x - min_x, 1.0)
    y_span = max(max_y - min_y, 1.0)
    x_pad = x_span * padding_fraction
    y_pad = y_span * max(0.02, padding_fraction * 0.45)
    axis.set_xlim(min_x - x_pad, max_x + x_pad)
    axis.set_ylim(min_y - y_pad, max_y + y_pad)
    axis.set_xticks([])
    axis.set_yticks([])


def _log_timing(label: str, started_at: float) -> float:
    elapsed = max(0.0, time.perf_counter() - started_at)
    message = f"[timing] {label}={elapsed:.2f}s"
    logging.info(message)
    print(message)
    return elapsed


def ensure_h5ad_compat_for_write(adata: ad.AnnData) -> ad.AnnData:
    """
    Normalize known problematic uns payloads so files remain readable across
    mixed anndata versions (notably older readers that fail on null IOSpec).
    """
    if not isinstance(getattr(adata, "uns", None), dict):
        return adata
    log1p_meta = adata.uns.get("log1p")
    if isinstance(log1p_meta, dict) and log1p_meta.get("base", "__missing__") is None:
        patched = dict(log1p_meta)
        patched["base"] = float(np.e)
        adata.uns["log1p"] = patched
    return adata


def _repair_null_log1p_base(path: Union[str, Path]) -> bool:
    """
    Repair legacy/invalid h5ad payloads where /uns/log1p/base is encoded as null.
    Returns True if a patch was applied.
    """
    patched = False
    with h5py.File(str(path), "r+") as handle:
        if "uns" not in handle or "log1p" not in handle["uns"]:
            return False
        log1p_group = handle["uns"]["log1p"]
        if not isinstance(log1p_group, h5py.Group) or "base" not in log1p_group:
            return False
        base_ds = log1p_group["base"]
        encoding = str(base_ds.attrs.get("encoding-type", "")).strip().lower()
        needs_patch = encoding == "null"
        if not needs_patch:
            try:
                value = base_ds[()]
                if isinstance(value, h5py.Empty):
                    needs_patch = True
            except Exception:
                needs_patch = True
        if not needs_patch:
            return False
        del log1p_group["base"]
        new_ds = log1p_group.create_dataset("base", data=float(np.e))
        new_ds.attrs["encoding-type"] = "numeric-scalar"
        new_ds.attrs["encoding-version"] = "0.2.0"
        patched = True
    return patched


@dataclass
class ApproximateUMAPResult:
    query_adata: ad.AnnData
    coordinates: pd.DataFrame
    reference_choices: Dict[str, List[str]]
    cluster_key: str
    cluster_order: List[str]
    plot_color_map: Optional[Dict[str, str]] = None
    plot_obs_field: Optional[str] = None
    plot_obs_value: Optional[str] = None
    plot_obs_mode: str = "include"
    dot_scale: float = 1.0

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
        custom_color_map: Optional[Mapping[str, str]] = None,
        restrict_obs_field: Optional[str] = None,
        restrict_obs_value: Optional[str] = None,
        restrict_obs_mode: Optional[str] = None,
        dot_scale: Optional[float] = None,
    ) -> tuple[str, str]:
        output_path = Path(output_path)
        annotated, plain = _save_comparison_pdf(
            reference_adata,
            self.query_adata,
            umap_key=umap_key,
            reference_cluster_key=reference_cluster_key,
            query_cluster_key=query_cluster_key,
            lineage_order=self.cluster_order,
            custom_color_map=custom_color_map or self.plot_color_map,
            restrict_obs_field=restrict_obs_field or self.plot_obs_field,
            restrict_obs_value=restrict_obs_value or self.plot_obs_value,
            restrict_obs_mode=restrict_obs_mode or self.plot_obs_mode,
            dot_scale=float(self.dot_scale if dot_scale is None else dot_scale),
            destination=output_path,
        )
        return annotated, plain

    def write_expression_pdfs(
        self,
        *,
        genes: Sequence[str],
        output_dir: Union[str, Path],
        query_cluster_key: str,
        umap_key: str = "X_umap",
        restrict_obs_field: Optional[str] = None,
        restrict_obs_value: Optional[str] = None,
        restrict_obs_mode: Optional[str] = None,
        dot_scale: Optional[float] = None,
    ) -> List[str]:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        written: List[str] = []
        resolved_field = restrict_obs_field or self.plot_obs_field
        resolved_values = restrict_obs_value or self.plot_obs_value
        resolved_mode = restrict_obs_mode or self.plot_obs_mode
        split_values = (
            _parse_obs_values(resolved_values)
            if resolved_field and resolved_mode == "include"
            else []
        )
        for gene in genes:
            gene_text = _clean_string(gene)
            if not gene_text:
                continue
            if len(split_values) > 1:
                for value_text in split_values:
                    pdf_path = output_dir / (
                        f"{_sanitize_filename_component(gene_text)}-"
                        f"{_sanitize_filename_component(value_text)}-expression.pdf"
                    )
                    _save_expression_pdf(
                        self.query_adata,
                        gene_text,
                        query_cluster_key=query_cluster_key,
                        umap_key=umap_key,
                        restrict_obs_field=resolved_field,
                        restrict_obs_value=value_text,
                        restrict_obs_mode=resolved_mode,
                        dot_scale=float(self.dot_scale if dot_scale is None else dot_scale),
                        plot_label=value_text,
                        destination=pdf_path,
                    )
                    written.append(str(pdf_path))
            else:
                pdf_path = output_dir / f"{_sanitize_filename_component(gene_text)}-expression.pdf"
                single_label = split_values[0] if len(split_values) == 1 else None
                _save_expression_pdf(
                    self.query_adata,
                    gene_text,
                    query_cluster_key=query_cluster_key,
                    umap_key=umap_key,
                    restrict_obs_field=resolved_field,
                    restrict_obs_value=resolved_values,
                    restrict_obs_mode=resolved_mode,
                    dot_scale=float(self.dot_scale if dot_scale is None else dot_scale),
                    plot_label=single_label,
                    destination=pdf_path,
                )
                written.append(str(pdf_path))
        return written


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
    custom_color_tsv: Optional[Union[str, Path]] = None,
    restrict_obs_field: Optional[str] = None,
    restrict_obs_value: Optional[str] = None,
    restrict_obs_mode: str = "include",
    dot_scale: float = 1.0,
    copy_query: bool = True,
) -> ApproximateUMAPResult:
    fn_started_at = time.perf_counter()
    if bool(restrict_obs_field) != bool(restrict_obs_value):
        raise ValueError(
            "restrict_obs_field and restrict_obs_value must be supplied together."
        )
    if restrict_obs_mode not in {"include", "exclude"}:
        raise ValueError("restrict_obs_mode must be either 'include' or 'exclude'.")
    step_started_at = time.perf_counter()
    query_adata = _load_adata(query)
    reference_adata = _load_adata(reference)
    _log_timing("approximate_umap.load_inputs", step_started_at)

    step_started_at = time.perf_counter()
    ref_cluster_key = reference_cluster_key or query_cluster_key
    plot_color_map = _load_custom_color_map(custom_color_tsv) if custom_color_tsv else None

    _validate_inputs(
        query_adata,
        reference_adata,
        query_cluster_key=query_cluster_key,
        reference_cluster_key=ref_cluster_key,
        umap_key=umap_key,
        jitter=jitter,
        num_reference_cells=num_reference_cells,
    )
    _log_timing("approximate_umap.validate_inputs", step_started_at)

    step_started_at = time.perf_counter()
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
    _log_timing("approximate_umap.resolve_cluster_order", step_started_at)

    step_started_at = time.perf_counter()
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
    _log_timing("approximate_umap.dataset_summary", step_started_at)

    step_started_at = time.perf_counter()
    coords_df = _embedding_to_df(reference_adata, umap_key)
    cells_by_cluster: Dict[str, np.ndarray] = {}
    coord_arrays_by_cluster: Dict[str, np.ndarray] = {}
    cluster_bounds: Dict[str, Tuple[float, float, float, float]] = {}
    coords_by_cluster = coords_df.copy()
    coords_by_cluster["cluster"] = ref_clusters.loc[coords_df.index].astype(str)
    for cluster in cluster_order:
        cluster_key = str(cluster)
        cluster_coords = coords_by_cluster[coords_by_cluster["cluster"] == cluster_key][["umap_0", "umap_1"]]
        if cluster_coords.empty:
            continue
        cluster_cells = cluster_coords.index.to_numpy(dtype=object, copy=False)
        cluster_coords_np = cluster_coords[["umap_0", "umap_1"]].to_numpy(dtype=np.float32, copy=False)
        q1 = cluster_coords.quantile(0.25)
        q3 = cluster_coords.quantile(0.75)
        iqr = q3 - q1
        lower = q1 - 1.5 * iqr
        upper = q3 + 1.5 * iqr
        mask = (
            (cluster_coords_np[:, 0] >= float(lower["umap_0"]))
            & (cluster_coords_np[:, 0] <= float(upper["umap_0"]))
            & (cluster_coords_np[:, 1] >= float(lower["umap_1"]))
            & (cluster_coords_np[:, 1] <= float(upper["umap_1"]))
        )
        filtered_cells = cluster_cells[mask]
        filtered_coords = cluster_coords_np[mask]
        if filtered_cells.size == 0:
            filtered_cells = cluster_cells
            filtered_coords = cluster_coords_np
        cells_by_cluster[cluster_key] = filtered_cells
        coord_arrays_by_cluster[cluster_key] = filtered_coords
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
    _log_timing("approximate_umap.prepare_reference_clusters", step_started_at)

    rng = np.random.default_rng(random_state)
    qry_clusters_array = qry_clusters.to_numpy(dtype=str, copy=False)
    obs_names_array = query_adata.obs_names.to_numpy(dtype=object, copy=False)
    coordinates_array = np.zeros((query_adata.n_obs, 2), dtype=np.float32)
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

    step_started_at = time.perf_counter()
    for cluster in cluster_order:
        cluster = str(cluster)
        query_indices = np.flatnonzero(qry_clusters_array == cluster)
        if query_indices.size == 0:
            continue

        ref_cells = cells_by_cluster.get(cluster)
        ref_coords = coord_arrays_by_cluster.get(cluster)
        if ref_cells is None or ref_coords is None or ref_cells.size == 0:
            raise ValueError(
                f"Reference cluster '{cluster}' has no cells and cannot be sampled."
            )

        reference_count = int(ref_cells.size)
        query_count = int(query_indices.size)
        sample_size = num_reference_cells
        high_density_cluster = bool(reference_count and query_count > 4 * reference_count)
        if high_density_cluster:
            sample_size = max(2, num_reference_cells)

        sampled_cells_matrix: Optional[np.ndarray] = None
        if sample_size == 1:
            sampled_idx = rng.integers(0, reference_count, size=query_count, endpoint=False)
            base_coords = ref_coords[sampled_idx]
            sampled_cells_matrix = ref_cells[sampled_idx][:, None]
        else:
            replace = reference_count < sample_size
            if replace:
                sampled_idx_matrix = rng.integers(
                    0,
                    reference_count,
                    size=(query_count, sample_size),
                    endpoint=False,
                )
            else:
                sampled_idx_matrix = np.empty((query_count, sample_size), dtype=np.int32)
                for idx in range(query_count):
                    sampled_idx_matrix[idx, :] = rng.choice(
                        reference_count,
                        size=sample_size,
                        replace=False,
                    )
            sampled_coords = ref_coords[sampled_idx_matrix]
            if high_density_cluster:
                base_coords = sampled_coords.mean(axis=1)
            else:
                base_coords = sampled_coords[:, 0, :]
            sampled_cells_matrix = ref_cells[sampled_idx_matrix]

        if jitter:
            offsets = rng.uniform(-jitter, jitter, size=(query_count, 2)).astype(np.float32)
            approx_coords = base_coords + offsets
        else:
            approx_coords = base_coords.copy()

        if high_density_cluster and cluster in cluster_bounds:
            x_low, x_high, y_low, y_high = cluster_bounds[cluster]
            approx_coords[:, 0] = np.clip(approx_coords[:, 0], x_low, x_high)
            approx_coords[:, 1] = np.clip(approx_coords[:, 1], y_low, y_high)

        coordinates_array[query_indices, :] = approx_coords
        if sampled_cells_matrix is not None:
            for local_idx, query_obs_index in enumerate(query_indices):
                barcode = str(obs_names_array[int(query_obs_index)])
                sampled[barcode] = sampled_cells_matrix[local_idx].astype(str).tolist()
    _log_timing("approximate_umap.assign_coordinates", step_started_at)

    step_started_at = time.perf_counter()
    print(f"Processed {query_adata.n_obs} query cells.")

    coordinates = pd.DataFrame(
        coordinates_array,
        index=query_adata.obs_names,
        columns=["umap_0", "umap_1"],
    )

    query_with_umap = query_adata.copy() if copy_query else query_adata
    query_with_umap.obsm[umap_key] = coordinates_array
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
    _log_timing("approximate_umap.finalize_outputs", step_started_at)
    _log_timing("approximate_umap.total", fn_started_at)

    return ApproximateUMAPResult(
        query_adata=query_with_umap,
        coordinates=coordinates,
        reference_choices=sampled,
        cluster_key=query_cluster_key,
        cluster_order=list(cluster_order),
        plot_color_map=plot_color_map,
        plot_obs_field=restrict_obs_field,
        plot_obs_value=restrict_obs_value,
        plot_obs_mode=restrict_obs_mode,
        dot_scale=float(dot_scale),
    )


def _load_adata(obj: AnnDataLike) -> ad.AnnData:
    if isinstance(obj, ad.AnnData):
        return obj
    path = str(obj)
    try:
        _repair_null_log1p_base(path)
    except Exception:
        # Best-effort compatibility repair; fall through to normal loading.
        pass
    try:
        return ad.read_h5ad(path)
    except Exception as exc:
        message = str(exc)
        if "encoding_type='null'" not in message or "/uns/log1p" not in message:
            raise
        if _repair_null_log1p_base(path):
            print(f"[compat] Repaired invalid /uns/log1p/base encoding in {path}")
            return ad.read_h5ad(path)
        raise


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
    custom_color_map: Optional[Mapping[str, str]] = None,
    restrict_obs_field: Optional[str] = None,
    restrict_obs_value: Optional[str] = None,
    restrict_obs_mode: str = "include",
    dot_scale: float = 1.0,
    destination: Path,
) -> tuple[str, str]:
    reference = _subset_adata_for_plot(
        reference,
        obs_field=restrict_obs_field,
        obs_value=restrict_obs_value,
        obs_mode=restrict_obs_mode,
        dataset_label="reference",
    )
    query = _subset_adata_for_plot(
        query,
        obs_field=restrict_obs_field,
        obs_value=restrict_obs_value,
        obs_mode=restrict_obs_mode,
        dataset_label="query",
    )

    if umap_key not in reference.obsm:
        raise KeyError(f"Reference AnnData missing '{umap_key}' in .obsm for plotting.")
    if umap_key not in query.obsm:
        raise KeyError(f"Query AnnData missing '{umap_key}' in .obsm for plotting.")

    ref_coords = _embedding_to_df(reference, umap_key).to_numpy()
    qry_coords = _embedding_to_df(query, umap_key).to_numpy()

    ref_labels = reference.obs[reference_cluster_key].astype(str)
    qry_labels = query.obs[query_cluster_key].astype(str)

    categories = _build_stable_umap_population_order(ref_labels.tolist(), qry_labels.tolist())
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

    palette = _resolve_colors(
        categories,
        reference,
        query,
        reference_cluster_key,
        query_cluster_key,
        custom_color_map=custom_color_map,
    )

    plt.rcParams.update(
        {
            "backend": "Agg",
            "axes.linewidth": 0.5,
            "pdf.fonttype": 42,
            "font.family": "sans-serif",
            "font.sans-serif": "DejaVu Sans",
            "figure.facecolor": "white",
        }
    )

    destination.parent.mkdir(parents=True, exist_ok=True)

    annotated_path = destination
    plain_path = destination.with_name(destination.stem + "-no-labels" + destination.suffix)

    def _render_pdf(path: Path, *, annotate: bool) -> None:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
        try:
            _scatter(axes[0], ref_coords, ref_labels, palette, "Reference UMAP", annotate=annotate, use_reference_grey=False, dot_scale=dot_scale)
            _scatter(axes[1], qry_coords, qry_labels, palette, "Query Approximate UMAP", annotate=annotate, use_reference_grey=False, dot_scale=dot_scale)
            fig.savefig(path, dpi=300, bbox_inches="tight")
        finally:
            plt.close(fig)

    annotated_ok = True
    try:
        _render_pdf(annotated_path, annotate=True)
    except Exception as exc:
        annotated_ok = False
        print(
            "[warn] Annotated approximate UMAP PDF rendering failed; "
            f"falling back to unlabeled output. {type(exc).__name__}: {exc}"
        )

    _render_pdf(plain_path, annotate=False)
    if not annotated_ok:
        _render_pdf(annotated_path, annotate=False)

    return str(annotated_path), str(plain_path)

def _scatter(
    axis: matplotlib.axes.Axes,
    coords: np.ndarray,
    labels: pd.Series,
    palette: Mapping[str, str],
    title: str,
    *,
    annotate: bool,
    use_reference_grey: bool = False,
    dot_scale: float = 1.0,
) -> None:
    colors = ["#e5e7eb" if use_reference_grey else palette.get(label, "#999999") for label in labels]
    point_size = _compute_point_size(coords.shape[0], dot_scale=dot_scale)
    axis.scatter(coords[:, 0], coords[:, 1], c=colors, s=point_size, linewidths=0, alpha=0.5 if not use_reference_grey else 0.3)
    axis.set_title(title, fontsize=12)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["bottom"].set_visible(False)
    axis.spines["left"].set_visible(False)
    _square_umap_axes(axis, coords, 0.06)
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
    custom_color_map: Optional[Mapping[str, str]] = None,
) -> Dict[str, str]:
    generated = _build_preview_palette([str(cat) for cat in categories])
    if custom_color_map:
        sanitized_custom: Dict[str, str] = {}
        for key, value in custom_color_map.items():
            try:
                sanitized_custom[str(key)] = to_hex(mcolors.to_rgb(value))
            except ValueError:
                continue
        resolved = {str(cat): sanitized_custom.get(str(cat), generated.get(str(cat), "#999999")) for cat in categories}
        missing = [str(cat) for cat in categories if str(cat) not in sanitized_custom]
        if missing:
            print(
                "[warn] Custom color TSV is missing colors for: "
                + ", ".join(missing)
                + ". Falling back to web-app palette for those categories."
            )
        return resolved
    return generated


def _load_custom_color_map(path: Union[str, Path]) -> Dict[str, str]:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if df.shape[1] < 2:
        raise ValueError("Custom colors TSV must contain at least two columns: annotation and color.")
    label_column = df.columns[0]
    color_column = df.columns[1]
    color_map: Dict[str, str] = {}
    for _, row in df.iterrows():
        label = _clean_string(row.get(label_column))
        color = _clean_string(row.get(color_column))
        if not label or not color:
            continue
        try:
            color_map[label] = to_hex(mcolors.to_rgb(color))
        except ValueError as exc:
            raise ValueError(f"Invalid color '{color}' for label '{label}' in {path}") from exc
    if not color_map:
        raise ValueError(f"No usable annotation/color pairs found in {path}")
    print(f"[info] Loaded {len(color_map)} custom colors from {path}")
    return color_map


def _parse_obs_values(obs_value: Optional[Union[str, Sequence[str]]]) -> List[str]:
    raw_values: list[str]
    if isinstance(obs_value, (list, tuple, set, np.ndarray, pd.Index)):
        raw_values = [_clean_string(item) for item in obs_value]
    else:
        raw_text = _clean_string(obs_value)
        raw_values = [value for value in re.split(r"[|,]", raw_text) if _clean_string(value)]
    return [value for value in raw_values if value]


def _subset_adata_for_plot(
    adata: ad.AnnData,
    *,
    obs_field: Optional[str],
    obs_value: Optional[Union[str, Sequence[str]]],
    obs_mode: str = "include",
    dataset_label: str,
) -> ad.AnnData:
    if not obs_field:
        return adata
    if obs_field not in adata.obs:
        print(f"[warn] {dataset_label} AnnData missing obs field '{obs_field}'; plotting all cells.")
        return adata
    values = _parse_obs_values(obs_value)
    if not values:
        return adata
    series = adata.obs[obs_field].astype(str)
    include_mask = series.isin(values)
    mask = ~include_mask if obs_mode == "exclude" else include_mask
    if not bool(mask.any()):
        action = "excluded" if obs_mode == "exclude" else "matched"
        raise ValueError(
            f"No {dataset_label} cells {action} obs['{obs_field}'] in {values!r} for plotting."
        )
    subset = adata[mask].copy()
    action_text = "Excluding" if obs_mode == "exclude" else "Restricting"
    print(
        f"[info] {action_text} {dataset_label} plot with obs['{obs_field}'] in {values!r} "
        f"({subset.n_obs}/{adata.n_obs} cells retained)."
    )
    return subset


def _compute_point_size(num_points: int, *, dot_scale: float = 1.0) -> float:
    base = 2.0
    scale = (1000.0 / max(num_points, 1)) ** 0.5
    size = base * scale
    dot_scale = max(float(dot_scale), 0.0)
    return float(np.clip(size, 0.5, 8.0) * dot_scale)


def _annotate_clusters(axis: matplotlib.axes.Axes, coords: np.ndarray, labels: pd.Series) -> None:
    label_points = _build_population_centroids(coords, labels.astype(str).tolist())
    if len(label_points) == 0:
        return
    relaxed = _relax_umap_labels(label_points, coords)
    for label in relaxed:
        axis.text(
            float(label["x"]),
            float(label["y"]),
            str(label["population"]),
            ha="center",
            va="center",
            fontsize=7,
            color="#0f172a",
            zorder=5,
        )


def _sanitize_filename_component(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "").strip()).strip("._")
    return cleaned or "output"


def _save_expression_pdf(
    query: ad.AnnData,
    gene: str,
    *,
    query_cluster_key: str,
    umap_key: str,
    restrict_obs_field: Optional[str],
    restrict_obs_value: Optional[Union[str, Sequence[str]]],
    restrict_obs_mode: str,
    dot_scale: float,
    plot_label: Optional[str],
    destination: Path,
) -> str:
    if gene not in set(query.var_names.astype(str).tolist()):
        raise KeyError(f"Gene '{gene}' was not found in the query AnnData.")
    subset = _subset_adata_for_plot(
        query,
        obs_field=restrict_obs_field,
        obs_value=restrict_obs_value,
        obs_mode=restrict_obs_mode,
        dataset_label="query",
    )
    if umap_key not in subset.obsm:
        raise KeyError(f"Query AnnData missing '{umap_key}' for expression plotting.")
    coords = _embedding_to_df(subset, umap_key).to_numpy()
    values = _flatten_expr(subset[:, gene].X).astype(float)
    zero_mask = np.isfinite(values) & (values <= 0)
    expr_mask = np.isfinite(values) & (values > 0)

    plt.rcParams.update(
        {
            "backend": "Agg",
            "axes.linewidth": 0.5,
            "pdf.fonttype": 42,
            "font.family": "sans-serif",
            "font.sans-serif": "DejaVu Sans",
            "figure.facecolor": "white",
        }
    )
    fig, ax = plt.subplots(figsize=(8.5, 8.5), constrained_layout=True)
    try:
        point_size = _compute_point_size(coords.shape[0], dot_scale=dot_scale)
        if np.any(zero_mask):
            ax.scatter(
                coords[zero_mask, 0],
                coords[zero_mask, 1],
                s=point_size,
                c="#e5e7eb",
                alpha=0.9,
                linewidths=0,
            )
        if np.any(expr_mask):
            expression_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                "expression_grey_red",
                [
                    (0.0, "#f3f4f6"),
                    (0.15, "#fecaca"),
                    (0.35, "#fca5a5"),
                    (0.6, "#ef4444"),
                    (1.0, "#b91c1c"),
                ],
            )
            sc = ax.scatter(
                coords[expr_mask, 0],
                coords[expr_mask, 1],
                s=point_size,
                c=values[expr_mask],
                cmap=expression_cmap,
                alpha=0.9,
                linewidths=0,
            )
            colorbar = fig.colorbar(sc, ax=ax)
            colorbar.set_label(gene)
        title = f"{gene} expression"
        if plot_label:
            title = f"{gene} expression: {plot_label}"
        ax.set_title(title, fontsize=12)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
        _square_umap_axes(ax, coords, 0.06)
        destination.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(destination, dpi=300, bbox_inches="tight")
    finally:
        plt.close(fig)
    return str(destination)


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
        "--custom-colors-tsv",
        help=(
            "TSV with at least two columns: the first is the cell-type annotation "
            "and the second is the color code."
        ),
    )
    parser.add_argument(
        "--restrict-obs-field",
        help="obs column used to restrict which cells are drawn in the PDF exports.",
    )
    parser.add_argument(
        "--restrict-obs-value",
        help="obs value to match when --restrict-obs-field is supplied. Comma- or pipe-delimited values are also accepted.",
    )
    parser.add_argument(
        "--restrict-obs-mode",
        choices=("include", "exclude"),
        default="include",
        help="Whether the supplied obs values are included or excluded from the PDF/expression plots.",
    )
    parser.add_argument(
        "--dot_scale",
        type=float,
        default=1.0,
        help="Multiplier applied to the current auto-computed UMAP point size (default: %(default)s).",
    )
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
        "--expression-gene",
        action="append",
        default=[],
        help="Gene to export as a UMAP expression PDF. Repeat to export multiple genes.",
    )
    parser.add_argument(
        "--expression-pdf-dir",
        help="Directory where gene expression PDFs will be written (default: <prefix>-expression within --outdir).",
    )
    parser.add_argument(
        "--skip-pdf",
        action="store_true",
        help="Skip comparison PDF generation.",
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
    main_started_at = time.perf_counter()
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

    step_started_at = time.perf_counter()
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
    _log_timing("main.load_reference_source", step_started_at)

    step_started_at = time.perf_counter()
    if args.query:
        query_source: AnnDataLike = args.query
    else:
        query_source = _load_query_from_tsv(
            args.query_clusters_tsv,
            cluster_key=args.query_cluster_key,
        )
    _log_timing("main.load_query_source", step_started_at)

    step_started_at = time.perf_counter()
    result = approximate_umap(
        query=query_source,
        reference=reference_source,
        query_cluster_key=args.query_cluster_key,
        reference_cluster_key=args.reference_cluster_key,
        umap_key=args.umap_key,
        jitter=args.jitter,
        num_reference_cells=args.num_reference_cells,
        random_state=args.random_state,
        custom_color_tsv=args.custom_colors_tsv,
        restrict_obs_field=args.restrict_obs_field,
        restrict_obs_value=args.restrict_obs_value,
        restrict_obs_mode=args.restrict_obs_mode,
        dot_scale=args.dot_scale,
        copy_query=False,
    )
    _log_timing("main.compute_approximate_umap", step_started_at)

    step_started_at = time.perf_counter()
    query_path = Path(args.query or args.query_clusters_tsv)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    prefix_name = args.output_prefix or f"{query_path.stem}-approximate-umap"
    prefix = out_dir / prefix_name

    output_pdf: Optional[Path] = None
    expression_pdf_dir: Optional[Path] = None
    if not args.skip_pdf:
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
    if args.expression_gene:
        if args.expression_pdf_dir:
            expression_dir_candidate = Path(args.expression_pdf_dir)
            expression_pdf_dir = (
                expression_dir_candidate
                if expression_dir_candidate.is_absolute()
                else out_dir / expression_dir_candidate
            )
        else:
            expression_pdf_dir = out_dir / f"{prefix_name}-expression"
        expression_pdf_dir.mkdir(parents=True, exist_ok=True)
    _log_timing("main.prepare_output_paths", step_started_at)

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

    step_started_at = time.perf_counter()
    if save_h5ad:
        output_h5ad.parent.mkdir(parents=True, exist_ok=True)
        logging.info("Writing augmented AnnData to %s", output_h5ad)
        print(f"[info] Saving augmented AnnData to {output_h5ad} (command: {full_command})")
        ensure_h5ad_compat_for_write(result.query_adata)
        result.query_adata.write(output_h5ad, compression="gzip")
    else:
        logging.info("Skipping augmented AnnData export (enable with --save-updated-h5ad).")
        print(f"[info] Skipping augmented AnnData export (enable with --save-updated-h5ad). Command: {full_command}")
    _log_timing("main.write_h5ad", step_started_at)

    step_started_at = time.perf_counter()
    logging.info("Writing legacy-style exports with prefix %s", prefix)
    result.write_text_outputs(prefix)
    _log_timing("main.write_text_outputs", step_started_at)

    step_started_at = time.perf_counter()
    if output_pdf is not None:
        logging.info("Writing comparison PDFs to %s (annotated) and companion without labels", output_pdf)
        reference_for_plot = reference_source if isinstance(reference_source, ad.AnnData) else _load_adata(reference_source)
        annotated_pdf, plain_pdf = result.write_comparison_pdf(
            reference_adata=reference_for_plot,
            umap_key=args.umap_key,
            reference_cluster_key=args.reference_cluster_key or args.query_cluster_key,
            query_cluster_key=args.query_cluster_key,
            output_path=output_pdf,
            custom_color_map=result.plot_color_map,
            restrict_obs_field=args.restrict_obs_field,
            restrict_obs_value=args.restrict_obs_value,
            restrict_obs_mode=args.restrict_obs_mode,
            dot_scale=args.dot_scale,
        )
        logging.info("Annotated UMAP PDF: %s", annotated_pdf)
        logging.info("Label-free UMAP PDF: %s", plain_pdf)
    else:
        logging.info("Skipping comparison PDF generation (--skip-pdf).")
    _log_timing("main.write_pdfs", step_started_at)

    step_started_at = time.perf_counter()
    if expression_pdf_dir is not None and args.expression_gene:
        written_expression = result.write_expression_pdfs(
            genes=args.expression_gene,
            output_dir=expression_pdf_dir,
            query_cluster_key=args.query_cluster_key,
            umap_key=args.umap_key,
            restrict_obs_field=args.restrict_obs_field,
            restrict_obs_value=args.restrict_obs_value,
            restrict_obs_mode=args.restrict_obs_mode,
            dot_scale=args.dot_scale,
        )
        for path in written_expression:
            logging.info("Expression UMAP PDF: %s", path)
        print(
            "[info] Wrote expression PDFs for genes: "
            + ", ".join([_clean_string(gene) for gene in args.expression_gene if _clean_string(gene)])
        )
    else:
        logging.info("Skipping expression PDF generation (no --expression-gene supplied).")
    _log_timing("main.write_expression_pdfs", step_started_at)

    logging.info("Approximate UMAP mapping completed.")
    _log_timing("main.total", main_started_at)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
