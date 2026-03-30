#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd


def _coerce_lineage_order(raw_order) -> Optional[List[str]]:
    if raw_order is None:
        return None
    if isinstance(raw_order, dict):
        raw_order = list(raw_order.values())
    if hasattr(raw_order, "tolist"):
        raw_order = raw_order.tolist()
    if isinstance(raw_order, str):
        return None
    try:
        return [str(x) for x in raw_order]
    except TypeError:
        return None


def export_reference_umap_bundle(
    h5ad_path: str | Path,
    *,
    cluster_key: str,
    output_dir: str | Path,
    output_prefix: Optional[str] = None,
    umap_key: str = "X_umap",
) -> Dict[str, str]:
    h5ad_path = Path(h5ad_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not h5ad_path.exists():
        raise FileNotFoundError(f"Reference h5ad file not found: {h5ad_path}")

    adata = ad.read_h5ad(h5ad_path)

    if cluster_key not in adata.obs.columns:
        available_obs = ", ".join(sorted(map(str, adata.obs.columns.tolist())))
        raise KeyError(
            f"Cluster key '{cluster_key}' not found in adata.obs. "
            f"Available obs columns: {available_obs}"
        )
    if umap_key not in adata.obsm:
        raise KeyError(f"UMAP key '{umap_key}' not found in adata.obsm")

    embedding = np.asarray(adata.obsm[umap_key])
    if embedding.ndim != 2 or embedding.shape[1] < 2:
        raise ValueError(f"Embedding '{umap_key}' must be at least two-dimensional.")

    prefix = output_prefix or h5ad_path.stem
    coords_path = output_dir / f"{prefix}_reference_umap.tsv"
    clusters_path = output_dir / f"{prefix}_reference_clusters.tsv"
    metadata_path = output_dir / f"{prefix}_reference_metadata.json"
    config_path = output_dir / f"{prefix}_reference_config_snippet.json"

    barcodes = adata.obs_names.astype(str)
    cluster_values = adata.obs[cluster_key].astype(str).values

    coords_df = pd.DataFrame(
        {
            "barcode": barcodes,
            "UMAP1": embedding[:, 0].astype(float),
            "UMAP2": embedding[:, 1].astype(float),
        }
    )
    coords_df.to_csv(coords_path, sep="\t", index=False)

    clusters_df = pd.DataFrame(
        {
            "barcode": barcodes,
            cluster_key: cluster_values,
            "Population": cluster_values,
        }
    )
    clusters_df.to_csv(clusters_path, sep="\t", index=False)

    lineage_order = _coerce_lineage_order(adata.uns.get("lineage_order"))
    cluster_colors = None
    color_key = f"{cluster_key}_colors"
    if color_key in adata.uns:
        colors = adata.uns[color_key]
        if hasattr(colors, "tolist"):
            colors = colors.tolist()
        try:
            categories = adata.obs[cluster_key].cat.categories.tolist()
        except AttributeError:
            categories = pd.unique(adata.obs[cluster_key].astype(str)).tolist()
        if isinstance(colors, list) and len(colors) == len(categories):
            cluster_colors = dict(zip([str(x) for x in categories], [str(x) for x in colors]))

    metadata = {
        "source_h5ad": str(h5ad_path),
        "cluster_key": cluster_key,
        "umap_key": umap_key,
        "n_cells": int(adata.n_obs),
        "n_features": int(adata.n_vars),
        "reference_coords_tsv": str(coords_path),
        "reference_clusters_tsv": str(clusters_path),
        "lineage_order": lineage_order,
        "cluster_colors": cluster_colors,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")

    config_snippet = {
        "id": f"{prefix}_reference",
        "label": f"{prefix} reference",
        "states_tsv": "/path/to/cellHarmony_reference_states.tsv",
        "reference_clusters_tsv": str(clusters_path),
        "reference_coords_tsv": str(coords_path),
        "cluster_key": cluster_key,
        "soupx_options": ["default"],
    }
    config_path.write_text(json.dumps(config_snippet, indent=2, sort_keys=True), encoding="utf-8")

    return {
        "reference_coords_tsv": str(coords_path),
        "reference_clusters_tsv": str(clusters_path),
        "reference_metadata_json": str(metadata_path),
        "reference_config_snippet_json": str(config_path),
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Export the reference UMAP coordinates and cluster annotation TSVs required "
            "by approximate_umap and the cellHarmony web alignment app."
        )
    )
    parser.add_argument("--h5ad", required=True, help="Input reference .h5ad file.")
    parser.add_argument("--cluster-key", required=True, help="obs column containing the reference labels.")
    parser.add_argument("--outdir", required=True, help="Directory to write the exported reference files.")
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Prefix for the exported files. Defaults to the input h5ad stem.",
    )
    parser.add_argument(
        "--umap-key",
        default="X_umap",
        help="obsm key holding the reference UMAP coordinates (default: X_umap).",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    outputs = export_reference_umap_bundle(
        args.h5ad,
        cluster_key=args.cluster_key,
        output_dir=args.outdir,
        output_prefix=args.output_prefix,
        umap_key=args.umap_key,
    )
    print("[INFO] Exported reference bundle:")
    for key, value in outputs.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
