#!/usr/bin/env python3
"""ICGS3 sparse single-cell clustering and marker discovery.

ICGS3 is the AltAnalyze3 replacement path for the legacy dense/HOPACH-based
ICGS2 workflow.  The official path keeps the matrix sparse through import, QC,
RNA gene filtering, ICGS2-style Louvain/PageRank downsampling, UDON feature
selection, MarkerFinder, and chunked SVM scoring wherever the underlying
algorithm allows it.
It deliberately reuses validated AltAnalyze3 components:

  - cellHarmony-lite import/QC semantics for h5ad, 10x h5, and mtx inputs
  - ICGS2-style PageRank downsampling before expensive NMF
  - UDON NMF, MarkerFinder cluster fitness, and linear SVM reclassification
  - AltAnalyze marker_heatmap_h5ad canonical heatmap outputs
  - AltDatabase BioMarkers gene-set enrichment for cell-state labels
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import sys
import tempfile
import time
import warnings
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MTX_SUFFIXES = (
    "_filtered_matrix.mtx.gz",
    "_counts.mtx.gz",
    "_matrix.mtx.gz",
    "filtered_matrix.mtx.gz",
    "matrix.mtx.gz",
    "matrix.mtx",
    "quants_mat.mtx",
)


@dataclass
class ICGS3Config:
    input_paths: List[str]
    output_dir: str
    modality: str = "rna"
    layer: Optional[str] = None
    input_normalized: bool = False
    normalized_decimals: int = -1
    min_genes: Optional[int] = 500
    min_cells: Optional[int] = 5
    min_counts: Optional[int] = 1000
    mito_percent: Optional[float] = 30.0
    target_cells: Optional[int] = None
    pagerank_cells: int = 5000
    louvain_downsample_cutoff: int = 15000
    pre_pagerank_cells: int = 0
    downsample_var_genes: int = 500
    retention_audit_obs: Optional[str] = "HLCA,TGEN-IPF"
    downsample_key: Optional[str] = None
    max_cells_per_group: Optional[int] = None
    random_state: int = 0
    n_top_features: int = 3000
    max_nmf_variable_genes: int = 0
    n_pcs: int = 0
    max_auto_pcs: int = 50
    n_neighbors: int = 30
    leiden_resolution: float = 0.8
    batch_correction: str = "none"
    batch_key: Optional[str] = None
    technology_batch_key: Optional[str] = None
    tissue_batch_key: Optional[str] = None
    batch_correction_use: str = "graph"
    batch_adjust_expression: bool = False
    expression_batch_key: Optional[str] = None
    cluster_key: str = "ICGS3_cluster"
    presvm_cluster_key: str = "ICGS3_NMF_cluster"
    marker_top_n: int = 60
    marker_rho: Optional[float] = None
    marker_min_per_cluster: Optional[int] = None
    marker_direction: str = "up"
    rank: Optional[int] = None
    nmf_assignment_normalization: str = "auto"
    nmf_k_multiplier: float = 2.0
    max_auto_nmf_k: Optional[int] = None
    nmf_runs: int = 1
    markerfinder_all_genes: bool = True
    resume_sampled_from: Optional[str] = None
    intercorr_threshold: float = 0.4
    corr_n_events: int = 5
    rank_rel_threshold: float = 0.1
    min_group_size: int = 3
    svm_min_decision_score: float = 0.0
    svm_reclassify_all_cells: bool = True
    svm_chunk_size: int = 50000
    species: str = "Hs"
    biomarker_file: Optional[str] = None
    heatmap_cells_per_cluster: int = 0
    write_heatmap_expression_tsv: bool = False
    heatmap_covariates: Optional[str] = None
    heatmap_goelite_terms: bool = False
    heatmap_goelite_max_terms: int = 30
    umap_feature_mode: str = "markerfinder"
    umap_covariates: Optional[str] = None
    umap_genes: Optional[str] = None
    write_h5ad: bool = True
    generate_umap: bool = True


@dataclass
class ICGS3Result:
    adata: ad.AnnData
    clusters: pd.DataFrame
    markers: pd.DataFrame
    biomarker_predictions: pd.DataFrame
    output_dir: str


def _log(message: str) -> None:
    print(f"[ICGS3] {message}", flush=True)


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
            stream.flush()
        return len(data)

    def flush(self):
        for stream in self.streams:
            stream.flush()

    @property
    def encoding(self):
        return "utf-8"

    def close(self):
        self.flush()


def _shell_join(parts: Sequence[object]) -> str:
    return " ".join(shlex.quote(str(p)) for p in parts if p is not None)


def cli_equivalent(config: ICGS3Config) -> str:
    """Reconstruct a reproducible CLI command from the resolved config."""
    cmd: List[object] = ["python3.11", "-m", "altanalyze3.components.clustering.ICGS", "--input"]
    cmd.extend(config.input_paths)
    cmd.extend(["--output-dir", config.output_dir, "--modality", config.modality])
    if config.layer:
        cmd.extend(["--layer", config.layer])
    if config.input_normalized:
        cmd.append("--input-normalized")
    for flag, value in (
        ("--min-genes", config.min_genes),
        ("--normalized-decimals", config.normalized_decimals),
        ("--min-cells", config.min_cells),
        ("--min-counts", config.min_counts),
        ("--mito-percent", config.mito_percent),
        ("--target-cells", config.target_cells),
        ("--pagerank-cells", config.pagerank_cells),
        ("--louvain-downsample-cutoff", config.louvain_downsample_cutoff),
        ("--pre-pagerank-cells", config.pre_pagerank_cells),
        ("--downsample-var-genes", config.downsample_var_genes),
        ("--retention-audit-obs", config.retention_audit_obs),
        ("--max-cells-per-group", config.max_cells_per_group),
        ("--random-state", config.random_state),
        ("--n-top-features", config.n_top_features),
        ("--max-nmf-variable-genes", config.max_nmf_variable_genes),
        ("--n-pcs", config.n_pcs),
        ("--max-auto-pcs", config.max_auto_pcs),
        ("--n-neighbors", config.n_neighbors),
        ("--leiden-resolution", config.leiden_resolution),
        ("--batch-correction", config.batch_correction),
        ("--batch-key", config.batch_key),
        ("--technology-batch-key", config.technology_batch_key),
        ("--tissue-batch-key", config.tissue_batch_key),
        ("--batch-correction-use", config.batch_correction_use),
        ("--expression-batch-key", config.expression_batch_key),
        ("--marker-top-n", config.marker_top_n),
        ("--marker-rho", config.marker_rho),
        ("--marker-min-per-cluster", config.marker_min_per_cluster),
        ("--nmf-k", config.rank),
        ("--nmf-assignment-normalization", config.nmf_assignment_normalization),
        ("--max-auto-nmf-k", config.max_auto_nmf_k),
        ("--nmf-runs", config.nmf_runs),
        ("--resume-sampled-from", config.resume_sampled_from),
        ("--intercorr-threshold", config.intercorr_threshold),
        ("--corr-n-events", config.corr_n_events),
        ("--rank-rel-threshold", config.rank_rel_threshold),
        ("--min-group-size", config.min_group_size),
        ("--svm-min-decision-score", config.svm_min_decision_score),
        ("--svm-chunk-size", config.svm_chunk_size),
        ("--heatmap-cells-per-cluster", config.heatmap_cells_per_cluster),
        ("--heatmap-covariates", config.heatmap_covariates),
        ("--heatmap-goelite-max-terms", config.heatmap_goelite_max_terms),
    ):
        if value is not None:
            cmd.extend([flag, value])
    if config.downsample_key:
        cmd.extend(["--downsample-key", config.downsample_key])
    if not config.svm_reclassify_all_cells:
        cmd.append("--no-svm-reclassify-all-cells")
    if config.batch_adjust_expression:
        cmd.append("--batch-adjust-expression")
    cmd.extend(["--cluster-key", config.cluster_key, "--marker-direction", config.marker_direction, "--species", config.species])
    if config.biomarker_file:
        cmd.extend(["--biomarker-file", config.biomarker_file])
    if config.write_heatmap_expression_tsv:
        cmd.append("--write-heatmap-expression-tsv")
    if config.heatmap_goelite_terms:
        cmd.append("--heatmap-goelite-terms")
    cmd.extend(["--umap-feature-mode", config.umap_feature_mode])
    if config.umap_covariates:
        cmd.extend(["--umap-covariates", config.umap_covariates])
    if config.umap_genes:
        cmd.extend(["--umap-genes", config.umap_genes])
    if not config.write_h5ad:
        cmd.append("--no-h5ad")
    if not config.generate_umap:
        cmd.append("--no-umap")
    return _shell_join(cmd)


def _ensure_sparse_csr(adata: ad.AnnData, layer: Optional[str] = None) -> None:
    target = adata.layers[layer] if layer else adata.X
    target = target.tocsr() if sp.issparse(target) else sp.csr_matrix(target)
    if layer:
        adata.layers[layer] = target
    else:
        adata.X = target


def _feature_paths_for_mtx(path: str) -> Tuple[str, str, str]:
    suffix = next((s for s in MTX_SUFFIXES if path.endswith(s)), None)
    if suffix is None:
        raise ValueError(f"Unsupported mtx filename: {path}")
    prefix = path[: -len(suffix)]
    if path.endswith(".gz"):
        barcodes = prefix + (
            "_barcodes.tsv.gz"
            if any(path.endswith(f"_{s}") for s in ("filtered_matrix.mtx.gz", "counts.mtx.gz", "matrix.mtx.gz"))
            else "barcodes.tsv.gz"
        )
        feature_options = [
            prefix + "_features.tsv.gz",
            prefix + "_genes.tsv.gz",
            prefix + "features.tsv.gz",
            prefix + "genes.tsv.gz",
        ]
    else:
        barcodes = prefix + (
            "_barcodes.tsv"
            if path.endswith("_matrix.mtx")
            else ("quants_mat_rows.txt" if path.endswith("quants_mat.mtx") else "barcodes.tsv")
        )
        feature_options = [
            prefix + "_features.tsv",
            prefix + "_genes.tsv",
            prefix + "features.tsv",
            prefix + "genes.tsv",
            prefix + "quants_mat_cols.txt",
        ]
    features = next((f for f in feature_options if os.path.exists(f)), None)
    if features is None:
        raise FileNotFoundError(f"Missing features/genes file for matrix prefix: {prefix}")
    if not os.path.exists(barcodes):
        raise FileNotFoundError(f"Missing barcodes file for matrix prefix: {prefix}")
    return prefix, barcodes, features


def read_single_input(
    path: str,
    sample_name: Optional[str] = None,
    *,
    append_sample_to_barcodes: bool = False,
) -> ad.AnnData:
    """Read one h5ad/10x-h5/10x-mtx input using cellHarmony-lite-compatible rules."""
    path = os.path.abspath(path)
    if os.path.isdir(path):
        candidates = sorted(
            os.path.join(path, name)
            for name in os.listdir(path)
            if name.endswith(MTX_SUFFIXES)
        )
        if not candidates:
            raise FileNotFoundError(f"No compatible mtx file found in directory: {path}")
        path = candidates[0]

    if path.endswith(".h5ad"):
        adata = sc.read_h5ad(path)
        inferred_sample = os.path.basename(path).replace(".h5ad", "")
    elif path.endswith(".h5"):
        adata = sc.read_10x_h5(path)
        inferred_sample = os.path.basename(path).replace(".h5", "")
    elif path.endswith(MTX_SUFFIXES):
        prefix, barcodes, features = _feature_paths_for_mtx(path)
        inferred_sample = os.path.basename(prefix.rstrip("_")) or os.path.basename(os.path.dirname(path))
        with tempfile.TemporaryDirectory(prefix="icgs3_mtx_") as tmpdir:
            matrix_dest = os.path.join(tmpdir, "matrix.mtx.gz" if path.endswith(".gz") else "matrix.mtx")
            shutil.copy(path, matrix_dest)
            shutil.copy(barcodes, os.path.join(tmpdir, os.path.basename(barcodes)))
            shutil.copy(features, os.path.join(tmpdir, "features.tsv.gz" if features.endswith(".gz") else "features.tsv"))
            try:
                adata = sc.read_10x_mtx(tmpdir, var_names="gene_symbols")
            except Exception:
                from anndata import AnnData
                from scipy.io import mmread

                matrix = mmread(matrix_dest).tocsr()
                if "quants_mat" not in os.path.basename(path):
                    matrix = matrix.T
                barcode_df = pd.read_csv(os.path.join(tmpdir, os.path.basename(barcodes)), header=None)
                feature_df = pd.read_csv(
                    os.path.join(tmpdir, "features.tsv.gz" if features.endswith(".gz") else "features.tsv"),
                    sep="\t",
                    header=None,
                )
                gene_symbols = feature_df[1].astype(str).tolist() if feature_df.shape[1] >= 2 else feature_df[0].astype(str).tolist()
                adata = AnnData(X=matrix)
                adata.obs_names = barcode_df[0].astype(str).tolist()
                adata.var_names = gene_symbols
    else:
        raise ValueError(f"Unsupported ICGS3 input: {path}")

    sample = sample_name or inferred_sample
    _ensure_sparse_csr(adata)
    if "gene_symbols" not in adata.var.columns:
        adata.var["gene_symbols"] = adata.var_names.astype(str)
    adata.var_names_make_unique()
    if sample:
        if append_sample_to_barcodes:
            adata.obs_names = [f"{bc}.{sample}" for bc in adata.obs_names.astype(str)]
        if "sample" not in adata.obs.columns:
            adata.obs["sample"] = str(sample)
        if "Library" not in adata.obs.columns:
            adata.obs["Library"] = str(sample)
    return adata


def read_inputs(input_paths: Sequence[str]) -> ad.AnnData:
    append_sample = len(input_paths) > 1
    adatas = [read_single_input(p, append_sample_to_barcodes=append_sample) for p in input_paths]
    if len(adatas) == 1:
        combined = adatas[0].copy()
    else:
        combined = ad.concat(adatas, label="sample", join="outer", fill_value=0)
    _ensure_sparse_csr(combined)
    combined.var_names_make_unique()
    return combined


def _qc_matrix(adata: ad.AnnData, layer: Optional[str] = None):
    if layer:
        return adata.layers[layer], f"layers['{layer}']"
    if "counts" in adata.layers:
        return adata.layers["counts"], "layers['counts']"
    return adata.X, "X"


def apply_qc(
    adata: ad.AnnData,
    *,
    min_genes: Optional[int],
    min_cells: Optional[int],
    min_counts: Optional[int],
    mito_percent: Optional[float],
    layer: Optional[str] = None,
) -> ad.AnnData:
    """Apply the same core sparse QC filters used by cellHarmony-lite."""
    X, source = _qc_matrix(adata, layer=layer)
    if source != "X":
        _log(f"QC using {source}")
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)

    if min_genes is not None and min_genes > 0:
        keep = np.diff(X.indptr) >= int(min_genes)
        adata = adata[keep].copy()
        X = _qc_matrix(adata, layer=layer)[0].tocsr()
        _log(f"cells after min_genes {min_genes}: {adata.n_obs}")

    if min_cells is not None and min_cells > 0:
        n_cells = np.bincount(X.indices, minlength=X.shape[1])
        keep = n_cells >= int(min_cells)
        adata = adata[:, keep].copy()
        X = _qc_matrix(adata, layer=layer)[0].tocsr()
        _log(f"features after min_cells {min_cells}: {adata.n_vars}")

    if min_counts is not None and min_counts > 0:
        counts = np.asarray(X.sum(axis=1)).ravel()
        keep = counts >= float(min_counts)
        adata = adata[keep].copy()
        X = _qc_matrix(adata, layer=layer)[0].tocsr()
        _log(f"cells after min_counts {min_counts}: {adata.n_obs}")

    if mito_percent is not None:
        mito = pd.Index(adata.var_names.astype(str)).str.upper().str.startswith("MT-")
        if mito.any():
            mito_counts = np.asarray(X[:, mito].sum(axis=1)).ravel()
            total_counts = np.maximum(np.asarray(X.sum(axis=1)).ravel(), 1e-12)
            adata.obs["pct_counts_mt"] = (mito_counts / total_counts) * 100.0
            adata = adata[adata.obs["pct_counts_mt"].values < float(mito_percent)].copy()
            _log(f"cells after mito_percent {mito_percent}: {adata.n_obs}")
    _ensure_sparse_csr(adata)
    return adata


def apply_modality_defaults(config: ICGS3Config) -> ICGS3Config:
    """Avoid RNA-specific QC and marker defaults for small non-RNA panels such as ADT.

    marker_rho and marker_min_per_cluster stay None until this function runs, so an
    explicit user value is distinguishable from an unset one. An earlier version instead
    compared them against the literal defaults of the day (3 and 0.2). Those defaults later
    became 2 and 0.3, which turned the ADT branch into dead code: every ADT run then used
    the RNA gate of rho 0.3, and the gate rejected ADT clusters that carry real markers.
    """
    if config.modality.lower() in {"adt", "grn", "metabolite", "lipid", "psi"}:
        if config.min_genes == 500:
            config.min_genes = 0
        if config.min_counts == 1000:
            config.min_counts = 0
        if config.mito_percent == 30.0:
            config.mito_percent = None
        if config.n_top_features == 3000:
            config.n_top_features = 0
    if config.marker_rho is None:
        config.marker_rho = 0.15 if config.modality.lower() == "adt" else 0.3
    if config.marker_min_per_cluster is None:
        config.marker_min_per_cluster = 1 if config.modality.lower() == "adt" else 2
    return config


def stratified_downsample_adata(
    adata: ad.AnnData,
    *,
    target_cells: Optional[int],
    downsample_key: Optional[str],
    max_cells_per_group: Optional[int],
    random_state: int,
) -> ad.AnnData:
    if (target_cells is None or target_cells <= 0 or adata.n_obs <= target_cells) and not max_cells_per_group:
        return adata

    rng = np.random.default_rng(random_state)
    selected: List[str] = []
    if downsample_key and downsample_key in adata.obs:
        groups = adata.obs[downsample_key].astype(str)
        n_groups = max(groups.nunique(), 1)
        per_group = max_cells_per_group or max(1, int(np.ceil(float(target_cells) / n_groups)))
        for _, idx in groups.groupby(groups).groups.items():
            idx = np.asarray(list(idx), dtype=object)
            n = min(len(idx), per_group)
            selected.extend(rng.choice(idx, size=n, replace=False).tolist())
    else:
        n = int(target_cells or adata.n_obs)
        selected = rng.choice(adata.obs_names.values, size=min(n, adata.n_obs), replace=False).tolist()

    selected = sorted(set(selected), key=lambda x: adata.obs_names.get_loc(x))
    _log(f"downsampled cells: {adata.n_obs} -> {len(selected)}")
    return adata[selected].copy()


def _copy_graph_slots(src: ad.AnnData, dest: ad.AnnData) -> None:
    for key in ("X_pca", "X_pca_uncorrected", "X_pca_harmony"):
        if key in src.obsm:
            dest.obsm[key] = src.obsm[key].copy()
    for key in ("distances", "connectivities"):
        if key in src.obsp:
            dest.obsp[key] = src.obsp[key].copy()
    if "neighbors" in src.uns:
        dest.uns["neighbors"] = dict(src.uns["neighbors"])
    if "icgs3_batch_correction" in src.uns:
        dest.uns["icgs3_batch_correction"] = dict(src.uns["icgs3_batch_correction"])


def _split_key_text(value: Optional[str]) -> List[str]:
    return [k.strip() for k in str(value or "").split(",") if k.strip()]


def _unique_keys(*values: Optional[str]) -> List[str]:
    keys: List[str] = []
    for value in values:
        for key in _split_key_text(value):
            if key not in keys:
                keys.append(key)
    return keys


def _batch_keys(config: ICGS3Config) -> List[str]:
    return _unique_keys(config.batch_key, config.technology_batch_key, config.tissue_batch_key)


def _expression_batch_keys(config: ICGS3Config) -> List[str]:
    base = config.expression_batch_key if config.expression_batch_key else config.batch_key
    return _unique_keys(base, config.technology_batch_key, config.tissue_batch_key)


def _batch_correction_applies(config: ICGS3Config, context: str) -> bool:
    method = str(config.batch_correction or "none").lower()
    use = str(config.batch_correction_use or "graph").lower()
    if method in {"none", "off", "false"}:
        return False
    if context == "graph":
        return use in {"graph", "graph-umap", "all"}
    if context == "umap":
        return use in {"umap", "graph-umap", "all"}
    return False


def apply_batch_corrected_pca(graph_adata: ad.AnnData, config: ICGS3Config, *, n_pcs: int) -> Optional[str]:
    """Optionally replace graph PCA coordinates with Harmony-corrected PCs.

    The correction is intentionally restricted to low-dimensional graph/UMAP
    construction. NMF, MarkerFinder, SVM, and heatmaps continue to use the
    original normalized feature values.
    """
    method = str(config.batch_correction or "none").lower()
    if method not in {"harmony"} or not _batch_correction_applies(config, "graph"):
        return None
    keys = _batch_keys(config)
    if not keys:
        _log("batch correction requested but --batch-key was not supplied; using uncorrected PCA")
        return None
    missing = [k for k in keys if k not in graph_adata.obs]
    if missing:
        _log(f"batch correction skipped; missing obs columns: {', '.join(missing)}")
        return None
    combo = graph_adata.obs[keys].astype(str).agg(" | ".join, axis=1)
    if combo.nunique() < 2:
        _log(f"batch correction skipped; batch keys {keys} define only one group")
        return None
    try:
        import harmonypy as hm
    except Exception as exc:
        _log(f"Harmony correction skipped; harmonypy import failed: {exc}")
        return None
    _log(
        f"Harmony batch correction on {graph_adata.n_obs} cells x {n_pcs} PCs "
        f"using batch key(s): {', '.join(keys)} ({combo.nunique()} groups)"
    )
    ho = hm.run_harmony(graph_adata.obsm["X_pca"], graph_adata.obs[keys].astype(str), vars_use=keys, verbose=False)
    corrected = np.asarray(ho.Z_corr)
    if corrected.shape != graph_adata.obsm["X_pca"].shape and corrected.T.shape == graph_adata.obsm["X_pca"].shape:
        corrected = corrected.T
    if corrected.shape != graph_adata.obsm["X_pca"].shape:
        raise ValueError(
            f"Harmony returned unexpected PCA shape {corrected.shape}; expected {graph_adata.obsm['X_pca'].shape}"
        )
    graph_adata.obsm["X_pca_uncorrected"] = graph_adata.obsm["X_pca"].copy()
    graph_adata.obsm["X_pca_harmony"] = corrected
    graph_adata.obsm["X_pca"] = corrected.copy()
    graph_adata.uns["icgs3_batch_correction"] = {
        "method": "harmony",
        "batch_key": keys,
        "technology_batch_key": _split_key_text(config.technology_batch_key),
        "tissue_batch_key": _split_key_text(config.tissue_batch_key),
        "n_groups": int(combo.nunique()),
        "scope": "PCA/neighbors/PageRank/initial-k/optional-UMAP",
    }
    return "X_pca_harmony"


def apply_expression_batch_adjustment(adata: ad.AnnData, config: ICGS3Config, outdir: str) -> ad.AnnData:
    """Create a temporary non-negative batch-centered analysis matrix.

    This is intended for low-dimensional ADT/intersected-panel workflows where
    panel or capture shifts can dominate NMF. It does not mutate the object used
    for final heatmap rendering or h5ad expression output.
    """
    if not config.batch_adjust_expression:
        return adata
    keys = _expression_batch_keys(config)
    if not keys:
        _log(
            "expression batch adjustment requested but no --expression-batch-key, "
            "--batch-key, --technology-batch-key, or --tissue-batch-key was supplied; using unadjusted matrix"
        )
        return adata
    missing = [k for k in keys if k not in adata.obs]
    if missing:
        _log(f"expression batch adjustment skipped; missing obs columns: {', '.join(missing)}")
        return adata
    combo = adata.obs[keys].astype(str).agg(" | ".join, axis=1)
    if combo.nunique() < 2:
        _log(f"expression batch adjustment skipped; keys {keys} define only one group")
        return adata

    X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
    X = X.astype(np.float32, copy=True)
    global_mean = X.mean(axis=0, keepdims=True)
    rows = []
    for batch, idx in combo.groupby(combo).groups.items():
        pos = adata.obs_names.get_indexer(list(idx))
        batch_mean = X[pos].mean(axis=0, keepdims=True)
        X[pos] = X[pos] - batch_mean + global_mean
        rows.append({"batch": batch, "cells": int(len(pos))})
    negative = int(np.sum(X < 0))
    if negative:
        X[X < 0] = 0.0
    adjusted = adata.copy()
    adjusted.X = sp.csr_matrix(X)
    adjusted.uns["icgs3_expression_batch_adjustment"] = {
        "method": "batch_mean_center_global_mean_clip_zero",
        "batch_key": keys,
        "technology_batch_key": _split_key_text(config.technology_batch_key),
        "tissue_batch_key": _split_key_text(config.tissue_batch_key),
        "n_groups": int(combo.nunique()),
        "negative_values_clipped": negative,
        "scope": "analysis matrix for PageRank/NMF/MarkerFinder/SVM/optional UMAP only",
    }
    pd.DataFrame(rows).to_csv(
        os.path.join(outdir, "sNMF", "icgs3_expression_batch_adjustment_groups.tsv"),
        sep="\t",
        index=False,
    )
    _log(
        f"expression batch adjustment for unsupervised analysis: keys={', '.join(keys)}, "
        f"groups={combo.nunique()}, clipped_negative_values={negative}"
    )
    return adjusted


def _select_elbow_n_pcs(variance_ratio: Sequence[float], *, min_pcs: int = 5) -> int:
    values = np.asarray(variance_ratio, dtype=float)
    if values.size <= 2:
        return max(1, int(values.size))
    x = np.arange(values.size, dtype=float)
    y = values
    start = np.array([x[0], y[0]], dtype=float)
    end = np.array([x[-1], y[-1]], dtype=float)
    line = end - start
    denom = np.linalg.norm(line)
    if denom <= 0:
        selected = int(np.argmax(values < (values[0] * 0.1)) + 1) if np.any(values < (values[0] * 0.1)) else values.size
    else:
        points = np.column_stack([x, y])
        distances = np.abs(np.cross(line, start - points)) / denom
        selected = int(np.argmax(distances) + 1)
    selected = max(int(min_pcs), selected)
    return min(selected, values.size)


def _resolve_graph_n_pcs(adata: ad.AnnData, config: ICGS3Config, *, requested: Optional[int] = None) -> Tuple[int, int, bool]:
    max_allowed = max(1, min(adata.n_obs - 1, adata.n_vars - 1))
    explicit = requested if requested is not None else config.n_pcs
    if explicit and int(explicit) > 0:
        n = min(int(explicit), max_allowed)
        return n, n, False
    candidate = min(max(2, int(config.max_auto_pcs)), max_allowed)
    return candidate, candidate, True


def compute_sparse_graph(
    adata: ad.AnnData,
    config: ICGS3Config,
    *,
    n_pcs: Optional[int] = None,
    use_pca: bool = True,
) -> ad.AnnData:
    graph_adata = adata.copy()
    neighbor_kwargs = {
        "n_neighbors": min(config.n_neighbors, max(2, graph_adata.n_obs - 1)),
        "method": "gauss",
    }
    if use_pca:
        pca_comps, selected_pcs, auto_pcs = _resolve_graph_n_pcs(graph_adata, config, requested=n_pcs)
        sc.pp.pca(graph_adata, n_comps=pca_comps, svd_solver="arpack")
        if auto_pcs:
            selected_pcs = _select_elbow_n_pcs(graph_adata.uns["pca"]["variance_ratio"], min_pcs=min(5, pca_comps))
            selected_pcs = max(1, min(int(selected_pcs), pca_comps))
            graph_adata.obsm["X_pca_full"] = graph_adata.obsm["X_pca"].copy()
            graph_adata.obsm["X_pca"] = graph_adata.obsm["X_pca"][:, :selected_pcs].copy()
            _log(f"automatic PC selection: elbow selected {selected_pcs}/{pca_comps} PCs")
        rep = apply_batch_corrected_pca(graph_adata, config, n_pcs=selected_pcs)
        if rep:
            neighbor_kwargs["use_rep"] = rep
        else:
            neighbor_kwargs["n_pcs"] = selected_pcs
        graph_meta = {
            "n_pcs": int(selected_pcs),
            "pca_computed": int(pca_comps),
            "pc_selection": "elbow" if auto_pcs else "explicit",
            "representation": rep or "X_pca",
        }
    else:
        neighbor_kwargs["use_rep"] = "X"
        graph_meta = {
            "n_pcs": 0,
            "pca_computed": 0,
            "pc_selection": "none",
            "representation": "X",
        }
        _log("graph construction for downsampling uses sparse expression directly (no PCA/PC selection)")
    sc.pp.neighbors(graph_adata, **neighbor_kwargs)
    graph_adata.uns["icgs3_graph"] = {
        **graph_meta,
        "n_neighbors": int(min(config.n_neighbors, max(2, graph_adata.n_obs - 1))),
        "batch_correction": str(config.batch_correction or "none").lower(),
        "batch_key": _batch_keys(config),
    }
    return graph_adata


def _sparse_distance_to_centroid(X, positions: np.ndarray) -> np.ndarray:
    sub = X[positions]
    if sp.issparse(sub):
        sub = sub.tocsr()
        center = np.asarray(sub.mean(axis=0)).ravel()
        center_norm = float(np.dot(center, center))
        row_norm = np.asarray(sub.multiply(sub).sum(axis=1)).ravel()
        dot = np.asarray(sub @ center).ravel()
        return row_norm + center_norm - (2.0 * dot)
    sub = np.asarray(sub, dtype=np.float32)
    center = sub.mean(axis=0)
    diff = sub - center
    return np.einsum("ij,ij->i", diff, diff)


def _icgs2_hgvfinder_adata(adata: ad.AnnData, num_var_genes: int) -> ad.AnnData:
    """ICGS2 hgvfinder equivalent: top dispersion genes with >5 unique values."""
    if num_var_genes is None or int(num_var_genes) <= 0 or adata.n_vars <= int(num_var_genes):
        _log(f"ICGS2 hgvfinder: using all {adata.n_vars} features for downsampling")
        return adata.copy()
    X = adata.X.tocsr() if sp.issparse(adata.X) else np.asarray(adata.X)
    means = np.asarray(X.mean(axis=0)).ravel()
    if sp.issparse(X):
        sq_means = np.asarray(X.multiply(X).mean(axis=0)).ravel()
        variances = sq_means - means**2
        unique_counts = []
        csc = X.tocsc()
        n_obs = X.shape[0]
        for j in range(X.shape[1]):
            values = csc[:, j].data
            unique_counts.append(len(np.unique(values)) + (1 if values.size < n_obs else 0))
        unique_counts = np.asarray(unique_counts)
    else:
        variances = np.var(X, axis=0)
        unique_counts = np.asarray([len(set(X[:, j].tolist())) for j in range(X.shape[1])])
    with np.errstate(divide="ignore", invalid="ignore"):
        dispersion = np.divide(variances, means, out=np.zeros_like(variances, dtype=float), where=means != 0)
    keepable = np.where((unique_counts > 5) & np.isfinite(dispersion))[0]
    if keepable.size == 0:
        keepable = np.arange(adata.n_vars)
    order = keepable[np.argsort(dispersion[keepable])[::-1]]
    keep_idx = order[: min(int(num_var_genes), order.size)]
    _log(f"ICGS2 hgvfinder: selected {len(keep_idx)}/{adata.n_vars} dispersion features for downsampling")
    return adata[:, keep_idx].copy()


def _annoy_neighbor_graph(matrix: np.ndarray, *, n_neighbors: int = 10, n_trees: int = 20, metric: str = "euclidean"):
    """Annoy k-NN graph.

    n_trees was 100. Annoy's own guidance is that trees trade index build time and
    memory against recall, and recall saturates well before 100 trees for small k.
    At k=10 on 280,446 cells the extra 80 trees cost build time and index memory for
    a recall gain that does not change the Louvain partition in practice. Raise it if
    you need higher recall.
    """
    from annoy import AnnoyIndex
    import networkx as nx

    matrix = np.asarray(matrix, dtype=np.float32)
    n_cells, n_features = matrix.shape
    _log(f"Annoy graph build start: {n_cells} cells x {n_features} features, k={n_neighbors}, trees={n_trees}, metric={metric}")
    index = AnnoyIndex(n_features, metric=metric)
    for i in range(n_cells):
        index.add_item(i, matrix[i])
        if (i + 1) % 10000 == 0:
            _log(f"Annoy indexed {i + 1}/{n_cells} cells")
    index.build(int(n_trees))
    _log("Annoy index built; collecting nearest-neighbor graph")
    adjacency = {}
    k = min(int(n_neighbors), n_cells)
    for i in range(n_cells):
        adjacency[i] = index.get_nns_by_item(i, k)
        if (i + 1) % 10000 == 0:
            _log(f"Annoy neighbors collected for {i + 1}/{n_cells} cells")
    return nx.from_dict_of_lists(adjacency), adjacency


def _adata_dense_cells_by_features(adata: ad.AnnData) -> np.ndarray:
    X = adata.X
    return X.toarray().astype(np.float32, copy=False) if sp.issparse(X) else np.asarray(X, dtype=np.float32)


def _mean_pairwise_euclidean_sums(matrix: np.ndarray, positions: Sequence[int]) -> np.ndarray:
    from sklearn.metrics import pairwise_distances_chunked

    positions = np.asarray(positions, dtype=int)
    sub = np.asarray(matrix[positions], dtype=np.float32)
    sums = np.zeros(sub.shape[0], dtype=float)
    start = 0
    for chunk in pairwise_distances_chunked(sub, metric="euclidean", working_memory=512):
        stop = start + chunk.shape[0]
        sums[start:stop] = np.asarray(chunk.sum(axis=1)).ravel()
        start = stop
    return sums / max(sub.shape[0], 1)


def _icgs2_community_sampling(
    adata_hvg: ad.AnnData,
    downsample_cutoff: int,
    pre_pagerank_cells: int = 0,
) -> Tuple[List[str], pd.DataFrame]:
    """ICGS2 community_sampling equivalent using Annoy k=10 and Louvain level 0."""
    try:
        import community
    except Exception as exc:
        raise ImportError("ICGS2-compatible community sampling requires python-louvain.") from exc

    matrix = _adata_dense_cells_by_features(adata_hvg)
    G, _ = _annoy_neighbor_graph(matrix, n_neighbors=10, n_trees=20, metric="euclidean")
    _log("ICGS2 community sampling: running Louvain dendrogram and using level 0 partition")
    dendrogram = community.generate_dendrogram(G)
    partition = community.partition_at_level(dendrogram, 0)
    communities: dict = {}
    for node, group in partition.items():
        communities.setdefault(group, []).append(int(node))
    _log(f"ICGS2 community sampling: identified {len(communities)} level-0 Louvain communities")
    downsample_limit = int(pre_pagerank_cells) if int(pre_pagerank_cells or 0) > 0 else int(downsample_cutoff) * 4
    downsample_limit = max(int(downsample_cutoff), int(downsample_limit))
    per_community = max(1, int(downsample_limit / max(len(communities), 1)))
    _log(
        "ICGS2 community sampling: "
        f"pre-PageRank target={downsample_limit}; selecting up to {per_community} "
        "medoid-proximal cells per community"
    )
    selected_pos: List[int] = []
    rows = []
    for idx, (group, positions) in enumerate(communities.items(), start=1):
        positions = list(positions)
        k = min(per_community, len(positions))
        dist = _mean_pairwise_euclidean_sums(matrix, positions)
        order = np.argsort(dist)[:k]
        chosen = [positions[i] for i in order]
        selected_pos.extend(chosen)
        for pos in chosen:
            rows.append(
                {
                    "barcode": str(adata_hvg.obs_names[pos]),
                    "precommunity": str(group),
                    "preselected": True,
                    "icgs2_community_size": len(positions),
                }
            )
        if idx % 25 == 0 or idx == len(communities):
            _log(f"ICGS2 community sampling: processed {idx}/{len(communities)} communities; selected {len(selected_pos)} cells")
    selected_pos = sorted(set(selected_pos))
    selected = adata_hvg.obs_names[selected_pos].astype(str).tolist()
    return selected, pd.DataFrame(rows)


def _icgs2_caldist(matrix: np.ndarray, current: int, keys: Iterable[int], keylist: Sequence[int]) -> int:
    visited = set(int(k) for k in keylist)
    candidates = [int(k) for k in keys if int(k) != int(current) and int(k) not in visited]
    if not candidates:
        remaining = [int(k) for k in keys if int(k) not in set(keylist)]
        if not remaining:
            return int(current)
        return remaining[0]
    diff = matrix[candidates] - matrix[int(current)]
    dist = np.einsum("ij,ij->i", diff, diff)
    return candidates[int(np.argmin(dist))]


def _icgs2_order_roots(matrix: np.ndarray, roots: Sequence[int]) -> List[int]:
    """Greedy nearest-root ordering equivalent to repeated caldist calls.

    ICGS2 uses this order only for presentation/sample order after PageRank has
    selected cells. Precomputing root-root distances preserves the greedy rule
    while avoiding repeated Python scans over thousands of roots.
    """
    roots = [int(r) for r in roots]
    if len(roots) <= 1:
        return roots
    root_matrix = np.asarray(matrix[roots], dtype=np.float32)
    try:
        from sklearn.metrics import pairwise_distances

        D = pairwise_distances(root_matrix, metric="euclidean", n_jobs=1)
        np.fill_diagonal(D, np.inf)
        ordered_idx = [0]
        remaining = np.ones(len(roots), dtype=bool)
        remaining[0] = False
        while len(ordered_idx) < len(roots):
            row = D[ordered_idx[-1]].copy()
            row[~remaining] = np.inf
            nxt = int(np.argmin(row))
            if not np.isfinite(row[nxt]):
                nxt = int(np.flatnonzero(remaining)[0])
            ordered_idx.append(nxt)
            remaining[nxt] = False
        return [roots[i] for i in ordered_idx]
    except Exception:
        keylist = [roots[0]]
        root_set = set(roots)
        while len(keylist) < len(roots):
            keylist.append(_icgs2_caldist(matrix, keylist[-1], root_set, keylist))
        return keylist


def _icgs2_pagerank_sampling_once(adata_hvg: ad.AnnData, downsample_cutoff: int) -> Tuple[List[str], pd.DataFrame]:
    import networkx as nx

    matrix = _adata_dense_cells_by_features(adata_hvg)
    n = matrix.shape[0]
    downsample_limit = int(downsample_cutoff) * 4
    sampled_names: List[str] = []
    rows = []
    for start in range(0, n, downsample_limit):
        stop = min(start + downsample_limit, n)
        _log(f"ICGS2 PageRank sampling chunk: cells {start + 1}-{stop} of {n}")
        chunk = matrix[start:stop]
        G, adjacency = _annoy_neighbor_graph(chunk, n_neighbors=10, n_trees=20, metric="angular")
        _log("ICGS2 PageRank sampling: computing networkx PageRank")
        pagerank = nx.pagerank(G)
        ranked = sorted(pagerank.items(), key=lambda item: item[1], reverse=True)
        keep_local = [int(node) for node, _ in ranked[: min(int(downsample_cutoff), len(ranked))]]
        _log(f"ICGS2 PageRank sampling: retained {len(keep_local)} top PageRank cells before neighbor grouping")
        keep_set = set(keep_local)
        grouped = []
        grouped_set = set()
        groups = {}
        for key in keep_local:
            if len(grouped) >= len(keep_local):
                break
            if key not in grouped_set:
                groups[key] = []
                grouped.append(key)
                grouped_set.add(key)
                for neighbor in adjacency.get(key, []):
                    if neighbor not in grouped_set and neighbor in keep_set:
                        groups[key].append(neighbor)
                        grouped.append(neighbor)
                        grouped_set.add(neighbor)
            else:
                for root, members in groups.items():
                    if key in members:
                        for neighbor in adjacency.get(key, []):
                            if neighbor not in grouped_set and neighbor in keep_set:
                                groups[root].append(neighbor)
                                grouped.append(neighbor)
                                grouped_set.add(neighbor)
        if not groups:
            continue
        _log(f"ICGS2 PageRank sampling: grouped {sum(len(v) for v in groups.values()) + len(groups)} cells under {len(groups)} roots")
        keylist = _icgs2_order_roots(chunk, list(groups.keys()))
        for key in keylist:
            global_pos = start + key
            sampled_names.append(str(adata_hvg.obs_names[global_pos]))
            rows.append(
                {
                    "barcode": str(adata_hvg.obs_names[global_pos]),
                    "pagerank": float(pagerank.get(key, 0.0)),
                    "selected": True,
                    "icgs2_chunk_start": start,
                    "icgs2_root": True,
                }
            )
            for member in groups.get(key, []):
                global_member = start + int(member)
                sampled_names.append(str(adata_hvg.obs_names[global_member]))
                rows.append(
                    {
                        "barcode": str(adata_hvg.obs_names[global_member]),
                        "pagerank": float(pagerank.get(int(member), 0.0)),
                        "selected": True,
                        "icgs2_chunk_start": start,
                        "icgs2_root": False,
                    }
                )
        _log(f"ICGS2 PageRank sampling chunk complete: selected {len(sampled_names)} cumulative cells")
    ordered = []
    selected_set = set(sampled_names)
    for name in adata_hvg.obs_names.astype(str):
        if name in selected_set:
            ordered.append(str(name))
    return ordered, pd.DataFrame(rows)


def _icgs2_pagerank_sampling(adata_hvg: ad.AnnData, downsample_cutoff: int) -> Tuple[List[str], pd.DataFrame]:
    selected, rows = _icgs2_pagerank_sampling_once(adata_hvg, downsample_cutoff)
    all_rows = [rows]
    iteration = 1
    while len(selected) > int(downsample_cutoff):
        iteration += 1
        _log(f"ICGS2 PageRank recursive pass {iteration}: {len(selected)} cells remain above cutoff {downsample_cutoff}")
        filtered = adata_hvg[selected].copy()
        selected, rows = _icgs2_pagerank_sampling_once(filtered, downsample_cutoff)
        all_rows.append(rows)
    return selected, pd.concat(all_rows, axis=0, ignore_index=True) if all_rows else pd.DataFrame()


def pagerank_downsample_adata(adata: ad.AnnData, config: ICGS3Config) -> Tuple[ad.AnnData, pd.DataFrame]:
    """ICGS2 graph/PageRank downsampling protocol.

    RNA inputs must already be restricted to the unsupervised ICGS gene space:
    protein-coding genes only, excluding mitochondrial, ribosomal, sex-linked and
    other RNASeq.py nuisance genes. Enforce that here as well as in the top-level
    workflow so direct API use cannot select the 500 dispersion genes from an
    invalid RNA feature universe.
    """
    if config.modality.lower() == "rna" and "ICGS_unsupervised_gene_filter" not in adata.var:
        _log("RNA downsampling guard: applying protein-coding/nuisance filter before selecting dispersion genes")
        adata = apply_rna_unsupervised_gene_filter(adata, config, outdir=None)
    if config.pagerank_cells <= 0 or adata.n_obs <= config.pagerank_cells:
        return adata.copy(), pd.DataFrame({"barcode": adata.obs_names, "selected": True})

    work = adata
    score_frames = []
    hvg = _icgs2_hgvfinder_adata(adata, int(config.downsample_var_genes))
    if adata.n_obs > int(config.louvain_downsample_cutoff):
        keep, community_rows = _icgs2_community_sampling(
            hvg,
            int(config.pagerank_cells),
            pre_pagerank_cells=int(config.pre_pagerank_cells or 0),
        )
        score_frames.append(community_rows)
        keep = sorted(set(keep), key=lambda x: adata.obs_names.get_loc(x))
        work = adata[keep].copy()
        _log(f"ICGS2 community sampling: {adata.n_obs} -> {work.n_obs} cells using Annoy k=10 and Louvain level 0")
        hvg = _icgs2_hgvfinder_adata(work, int(config.downsample_var_genes))
    selected, pagerank_rows = _icgs2_pagerank_sampling(hvg, int(config.pagerank_cells))
    score_frames.append(pagerank_rows)
    selected = sorted(set(selected), key=lambda x: adata.obs_names.get_loc(x))
    score_df = pd.concat(score_frames, axis=0, ignore_index=True) if score_frames else pd.DataFrame()
    if not score_df.empty:
        score_df["selected_final"] = score_df["barcode"].astype(str).isin(selected)
    _log(f"ICGS2 PageRank sampling: {work.n_obs} -> {len(selected)} cells")
    return adata[selected].copy(), score_df


def write_retention_audit(
    *,
    full: ad.AnnData,
    candidates: ad.AnnData,
    sampled: ad.AnnData,
    final_clusters: Optional[pd.DataFrame],
    config: ICGS3Config,
    outdir: str,
) -> pd.DataFrame:
    columns = [c.strip() for c in str(config.retention_audit_obs or "").split(",") if c.strip()]
    columns = [c for c in columns if c in full.obs]
    if not columns:
        return pd.DataFrame()
    candidate_names = set(candidates.obs_names.astype(str))
    sampled_names = set(sampled.obs_names.astype(str))
    final_names = set(final_clusters.index.astype(str)) if final_clusters is not None else set()
    rows = []
    for col in columns:
        labels = full.obs[col].astype(str)
        counts = labels.value_counts(dropna=False)
        rare_labels = counts[(counts > 20) & (counts < 1000)]
        for label, total in rare_labels.items():
            cells = set(labels.index[labels == label].astype(str))
            candidate_n = len(cells & candidate_names)
            sampled_n = len(cells & sampled_names)
            final_n = len(cells & final_names) if final_clusters is not None else np.nan
            rows.append(
                {
                    "obs_column": col,
                    "annotation": label,
                    "post_qc_cells": int(total),
                    "sampling_candidate_cells": int(candidate_n),
                    "pagerank_louvain_sampled_cells": int(sampled_n),
                    "final_svm_retained_cells": final_n,
                    "candidate_below_20": bool(candidate_n < 20),
                    "sampled_below_20": bool(sampled_n < 20),
                    "final_below_20": bool(final_clusters is not None and final_n < 20),
                }
            )
    audit = pd.DataFrame(rows)
    if audit.empty:
        return audit
    audit = audit.sort_values(["obs_column", "post_qc_cells", "annotation"], kind="mergesort")
    path = os.path.join(outdir, "sNMF", "icgs3_rare_population_retention_audit.tsv")
    audit.to_csv(path, sep="\t", index=False)
    failures = audit[audit["sampled_below_20"]]
    if not failures.empty:
        _log(
            f"rare-population retention audit: {failures.shape[0]} populations with >20 and <1000 "
            f"post-QC cells have <20 cells after PageRank/Louvain sampling; see {path}"
        )
    else:
        _log(f"rare-population retention audit: all audited >20/<1000 populations retained >=20 sampled cells; see {path}")
    return audit


def round_normalized_expression(adata: ad.AnnData, decimals: Optional[int]) -> ad.AnnData:
    if decimals is None or int(decimals) < 0:
        return adata
    decimals = int(decimals)
    if sp.issparse(adata.X):
        adata.X = adata.X.tocsr(copy=True)
        adata.X.data = np.round(adata.X.data, decimals=decimals).astype(np.float32, copy=False)
        adata.X.eliminate_zeros()
    else:
        adata.X = np.round(np.asarray(adata.X, dtype=np.float32), decimals=decimals)
    adata.uns["icgs3_normalized_decimals"] = decimals
    _log(f"normalized expression rounded to {decimals} decimal places")
    return adata


def prepare_expression(adata: ad.AnnData, config: ICGS3Config) -> ad.AnnData:
    if config.layer:
        if config.layer not in adata.layers:
            raise KeyError(f"Layer '{config.layer}' not found in AnnData.layers.")
        adata.X = adata.layers[config.layer].copy()
        _ensure_sparse_csr(adata)

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    if not config.input_normalized:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.uns["icgs3_normalization"] = "log1p(CP10K)"
    else:
        adata.uns["icgs3_normalization"] = "input treated as normalized"
    _ensure_sparse_csr(adata)
    adata = round_normalized_expression(adata, config.normalized_decimals)
    _ensure_sparse_csr(adata)
    return adata


def select_features(adata: ad.AnnData, n_top_features: int) -> ad.AnnData:
    if n_top_features is None or n_top_features <= 0 or adata.n_vars <= n_top_features:
        adata.var["icgs3_feature"] = True
        return adata
    sc.pp.highly_variable_genes(adata, n_top_genes=int(n_top_features), flavor="seurat")
    if "highly_variable" not in adata.var or int(adata.var["highly_variable"].sum()) == 0:
        variances = np.asarray(adata.X.power(2).mean(axis=0)).ravel() - np.asarray(adata.X.mean(axis=0)).ravel() ** 2
        keep_idx = np.argsort(variances)[::-1][: int(n_top_features)]
        mask = np.zeros(adata.n_vars, dtype=bool)
        mask[keep_idx] = True
        adata.var["highly_variable"] = mask
    adata.var["icgs3_feature"] = adata.var["highly_variable"].astype(bool)
    _log(f"selected features for graph: {int(adata.var['icgs3_feature'].sum())}/{adata.n_vars}")
    return adata


def _add_udon_path() -> str:
    udon_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "udon"))
    if udon_dir not in sys.path:
        sys.path.insert(0, udon_dir)
    return udon_dir


def _adata_to_dense_df(adata: ad.AnnData, genes: Sequence[str], cells: Optional[Sequence[str]] = None) -> pd.DataFrame:
    cells = list(cells) if cells is not None else adata.obs_names.tolist()
    cell_idx = adata.obs_names.get_indexer(cells)
    gene_idx = adata.var_names.get_indexer(pd.Index(genes).astype(str))
    valid_genes = gene_idx >= 0
    if np.any(cell_idx < 0) or not np.any(valid_genes):
        raise ValueError("Could not align cells/genes for ICGS3 matrix conversion.")
    # Slice rows and columns in ONE indexing pass. The previous chained form
    # adata.X[cell_idx, :][:, gene_idx] materialised a full-width intermediate
    # (all vars for the selected cells) before narrowing to the wanted genes,
    # roughly doubling peak memory on the sampled matrix.
    keep_gene_idx = gene_idx[valid_genes]
    if sp.issparse(adata.X):
        X = adata.X[cell_idx, :].tocsc()[:, keep_gene_idx]
        values = X.toarray().astype(np.float32, copy=False)
    else:
        values = np.asarray(adata.X[np.ix_(cell_idx, keep_gene_idx)], dtype=np.float32)
    # np.asfortranarray on the transpose keeps the genes x cells DataFrame from
    # triggering a second full copy inside the pandas block manager.
    return pd.DataFrame(
        np.asfortranarray(values.T),
        index=pd.Index(genes)[valid_genes].astype(str),
        columns=cells,
    )


def _select_icgs3_features(adata: ad.AnnData, config: ICGS3Config) -> pd.Index:
    if config.modality.lower() in {"adt", "grn", "metabolite", "lipid", "psi"}:
        return pd.Index(adata.var_names.astype(str))
    _add_udon_path()
    import fast_feature_selection as FFS

    pc = None
    pc_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "udon")), "ProteinCoding-Hs-Mm.txt")
    if os.path.exists(pc_path):
        pc = load_protein_coding_genes(pc_path)
    keep = legacy_rnaseq_gene_filter(adata, protein_coding_genes=pc, species=config.species)
    if int(np.sum(keep)) == 0:
        keep = legacy_rnaseq_gene_filter(adata, protein_coding_genes=None, species=config.species)
        _log(
            "RNA protein-coding lookup matched 0 features for UDON feature selection; "
            f"using nuisance-filter-only feature universe ({int(np.sum(keep))} features)"
        )
    # Densify ONLY the retained genes. The previous form densified every var first
    # and then took a .loc subset, holding two dense copies of the sampled matrix
    # at once (15,000 cells x 18,100 genes float32 = 1.09 GB each) when the full
    # gene universe is never used beyond this filter. The filtered gene set and
    # its values are identical either way.
    all_vars = pd.Index(adata.var_names.astype(str))
    filtered = _adata_to_dense_df(adata, all_vars[keep])
    _log(f"RNASeq.py-compatible gene filter for NMF feature selection: {all_vars.size} -> {filtered.shape[0]} features")
    # UDON fast_feature_selection stage-2 thresholds. Defaults (0.4 / 5) reproduce the
    # prior behaviour exactly. Relaxing them widens the NMF guide-gene pool: a sweep on
    # the Adams 30k sampling measured 0.4/5 -> 2,370 features capturing 5 of 9 canonical
    # ionocyte markers, and 0.3/3 -> 4,167 features capturing 7 of 9 including CFTR.
    fs = FFS.fast_feature_selection(
        filtered,
        None,
        apply_gene_name_filter=False,
        do_pca=False,
        intercorr_threshold=float(config.intercorr_threshold),
        corr_n_events=int(config.corr_n_events),
    )
    genes = pd.Index(fs["correlated_genes"])
    udon_hvg = pd.Index(fs["udon_hvg"])
    min_required = max(10, int(config.rank or 0) * 2)
    if len(genes) < min_required and len(udon_hvg) > len(genes):
        prior = len(genes)
        genes = udon_hvg
        _log(
            "UDON intercorrelation feature selection returned too few genes for stable NMF "
            f"({prior}); using UDON variance-selected genes ({len(genes)} features)"
        )
    if len(genes) == 0:
        raise ValueError("UDON feature selection returned no RNA features for NMF.")
    max_genes = int(config.max_nmf_variable_genes or 0)
    if max_genes > 0 and len(genes) > max_genes:
        ranked = _rank_features_by_sparse_dispersion(adata, genes)
        genes = pd.Index(ranked[:max_genes])
        _log(f"Capped UDON NMF variable features by sparse dispersion rank: {len(ranked)} -> {len(genes)}")
    return genes


def _rank_features_by_sparse_dispersion(adata: ad.AnnData, genes: Sequence[str]) -> List[str]:
    gene_index = pd.Index(adata.var_names.astype(str))
    positions = gene_index.get_indexer(pd.Index(genes).astype(str))
    valid = positions >= 0
    if not np.any(valid):
        return list(pd.Index(genes).astype(str))
    names = pd.Index(genes).astype(str)[valid].to_numpy()
    X = adata.X[:, positions[valid]]
    if sp.issparse(X):
        X = X.tocsr()
        means = np.asarray(X.mean(axis=0)).ravel()
        second = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    else:
        arr = np.asarray(X, dtype=np.float32)
        means = arr.mean(axis=0)
        second = np.square(arr).mean(axis=0)
    variances = np.maximum(second - np.square(means), 0)
    dispersion = variances / np.maximum(means, 1e-8)
    order = np.argsort(np.nan_to_num(dispersion, nan=-np.inf, posinf=-np.inf))[::-1]
    return names[order].tolist()


def load_protein_coding_genes(path: Optional[str] = None) -> set:
    if path is None:
        path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "udon")), "ProteinCoding-Hs-Mm.txt")
    if not path or not os.path.exists(path):
        return set()
    values = pd.read_csv(path, sep="\t", names=["g"])["g"].dropna().astype(str)
    lowered = set()
    for value in values:
        lowered.add(value.lower())
        lowered.add(value.split("_")[0].lower())
    return lowered


def apply_rna_unsupervised_gene_filter(adata: ad.AnnData, config: ICGS3Config, outdir: Optional[str] = None) -> ad.AnnData:
    """Remove non-protein-coding/nuisance RNA genes before graph, PageRank, NMF, and UMAP."""
    if config.modality.lower() != "rna":
        return adata
    pc = load_protein_coding_genes()
    keep = legacy_rnaseq_gene_filter(adata, protein_coding_genes=pc, species=config.species)
    before = adata.n_vars
    after = int(np.sum(keep))
    if after == 0:
        keep = legacy_rnaseq_gene_filter(adata, protein_coding_genes=None, species=config.species)
        after = int(np.sum(keep))
        _log(
            "RNA protein-coding lookup matched 0 features; falling back to RNASeq.py nuisance-gene "
            f"filter only ({before} -> {after} features)"
        )
    if after == 0:
        raise ValueError("RNA unsupervised gene filtering retained 0 features.")
    filtered = adata[:, keep].copy()
    filtered.var["ICGS_unsupervised_gene_filter"] = True
    _ensure_sparse_csr(filtered)
    _log(
        f"RNA unsupervised gene filter before graph/PageRank/NMF: {before} -> {after} features "
        f"(protein-coding lower-case symbol/id match when available; removed RPL/RPS, MT-, dotted, GM, XIS/TSI, RSP, HLA, *Y)"
    )
    if outdir:
        pd.DataFrame({"feature": filtered.var_names.astype(str)}).to_csv(
            os.path.join(outdir, "sNMF", "icgs3_rna_unsupervised_filtered_genes.tsv"),
            sep="\t",
            index=False,
        )
    return filtered


def _feature_selection_label(config: ICGS3Config, n_features: int, n_total: int) -> str:
    if config.modality.lower() in {"adt", "grn", "metabolite", "lipid", "psi"}:
        return f"{config.modality} panel features ({n_features}/{n_total}; no RNA gene-symbol filter)"
    return f"RNASeq.py nuisance/protein-coding filter + UDON fast_feature_selection ({n_features}/{n_total})"


def _resolve_auto_nmf_rank(udon_estimated_rank: int, config: ICGS3Config) -> Tuple[int, int, Optional[int]]:
    candidate_rank = max(2, int(udon_estimated_rank))
    max_auto_rank = int(config.max_auto_nmf_k) if config.max_auto_nmf_k is not None else None
    final_rank = min(candidate_rank, max_auto_rank) if max_auto_rank is not None else candidate_rank
    return candidate_rank, final_rank, max_auto_rank


def estimate_udon_nmf_rank(expr_sampled: pd.DataFrame, config: ICGS3Config) -> int:
    """ICGS2/UDON auto-rank estimator.

    AltAnalyze2 ICGS2 estimateK returns k from a Tracy-Widom eigenvalue boundary,
    then runs NMF at rank k*2. UDON's determine_nmf_ranks implements that same
    behavior and already returns the multiplied rank. Do not substitute graph
    clustering here without explicit protocol approval.
    """
    _add_udon_path()
    from nmf import determine_nmf_ranks

    return max(2, int(determine_nmf_ranks(expr_sampled, small_feature=False)))


def legacy_rnaseq_gene_filter(
    adata: ad.AnnData,
    *,
    protein_coding_genes: Optional[set] = None,
    species: str = "Hs",
) -> np.ndarray:
    """Match the conventional RNASeq.py ICGS filtering rules for RNA mode.

    RNASeq.py marks non-protein-coding genes as ncRNA and additionally excludes
    highly correlated nuisance symbols: RPL/RPS, MT-, dotted symbols, GM/Gm,
    XIS/TSI, RSP, HLA, and symbols ending in Y.
    """
    var_names = pd.Index(adata.var_names.astype(str))
    if "gene_symbols" in adata.var:
        symbols = adata.var["gene_symbols"].astype(str).values
    else:
        symbols = var_names.astype(str).values
    keep = np.ones(len(var_names), dtype=bool)
    if protein_coding_genes:
        protein_coding_genes = {str(g).lower() for g in protein_coding_genes}
        def is_pc(gene: str, symbol: str) -> bool:
            root = str(gene).split("_")[0].lower()
            sym = str(symbol).split("_")[0].lower()
            gene_l = str(gene).lower()
            symbol_l = str(symbol).lower()
            return (
                root in protein_coding_genes
                or sym in protein_coding_genes
                or gene_l in protein_coding_genes
                or symbol_l in protein_coding_genes
            )

        keep &= np.array([is_pc(g, s) for g, s in zip(var_names, symbols)], dtype=bool)
    lowered = pd.Index([str(s).lower() for s in symbols])
    nuisance = (
        lowered.str.startswith("rpl")
        | lowered.str.startswith("rps")
        | lowered.str.startswith("mt-")
        | lowered.str.contains(r"\.", regex=True)
        | lowered.str.startswith("gm")
        | lowered.str.startswith("xis")
        | lowered.str.startswith("tsi")
        | lowered.str.startswith("rsp")
        | lowered.str.startswith("hla")
        | lowered.str.endswith("y")
    )
    keep &= ~np.asarray(nuisance)
    return keep


def run_nmf_marker_svm(sampled: ad.AnnData, full: ad.AnnData, config: ICGS3Config, outdir: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    _add_udon_path()
    from markerFinder import marker_finder_wrapper
    from nmf import run_nmf
    from linearSVM import generate_train_data

    genes = _select_icgs3_features(sampled, config)
    expr_sampled = _adata_to_dense_df(sampled, genes)
    snmf_dir = os.path.join(outdir, "sNMF")
    marker_dir = os.path.join(outdir, "MarkerFinder")
    os.makedirs(snmf_dir, exist_ok=True)
    os.makedirs(marker_dir, exist_ok=True)
    pd.DataFrame({"feature": pd.Index(expr_sampled.index).astype(str)}).to_csv(
        os.path.join(snmf_dir, "icgs3_nmf_variable_features.tsv"),
        sep="\t",
        index=False,
    )
    _log(f"NMF variable features: {_feature_selection_label(config, expr_sampled.shape[0], sampled.n_vars)}")
    rank = config.rank
    rank_rows = []
    if rank is None:
        udon_rank = estimate_udon_nmf_rank(expr_sampled, config)
        candidate_rank, rank, max_auto_rank = _resolve_auto_nmf_rank(udon_rank, config)
        rank_rows.append(
            {
                "mode": "auto",
                "estimator": "UDON/ICGS2 Tracy-Widom eigenvalue rank",
                "udon_estimated_rank": int(udon_rank),
                "initial_populations": np.nan,
                "nmf_k_multiplier": np.nan,
                "candidate_rank": candidate_rank,
                "max_auto_nmf_k": max_auto_rank,
                "final_rank": int(rank),
            }
        )
        cap_text = str(max_auto_rank) if max_auto_rank is not None else "none"
        _log(
            "auto NMF k: UDON/ICGS2 Tracy-Widom eigenvalue estimator "
            f"returned rank={udon_rank}, max_auto_nmf_k={cap_text} -> rank={rank}"
        )
    else:
        max_auto_rank = int(config.max_auto_nmf_k) if config.max_auto_nmf_k is not None else None
        rank_rows.append(
            {
                "mode": "manual",
                "estimator": "manual override",
                "udon_estimated_rank": np.nan,
                "initial_populations": np.nan,
                "nmf_k_multiplier": np.nan,
                "candidate_rank": np.nan,
                "max_auto_nmf_k": max_auto_rank,
                "final_rank": int(rank),
            }
        )
        _log(f"manual NMF k override: rank={rank}")
    rank = min(int(rank), max(2, sampled.n_obs - 1))
    rank_rows[0]["final_rank"] = int(rank)
    pd.DataFrame(rank_rows).to_csv(os.path.join(snmf_dir, "icgs3_nmf_rank_selection.tsv"), sep="\t", index=False)
    assignment_norm = str(config.nmf_assignment_normalization or "auto").lower()
    if assignment_norm == "auto":
        assignment_norm = "rowsum" if config.modality.lower() in {"adt", "grn", "metabolite", "lipid", "psi"} else "raw"
    if assignment_norm in {"none", "raw"}:
        assignment_norm_arg = None
        assignment_label = "raw"
    else:
        assignment_norm_arg = assignment_norm
        assignment_label = assignment_norm
    _log(
        f"running UDON NMF rank={rank} on {expr_sampled.shape[1]} sampled cells x {expr_sampled.shape[0]} features "
        f"(assignment_normalization={assignment_label}, n_run={int(config.nmf_runs)})"
    )
    _, nmf_clusters = run_nmf(
        expr_sampled,
        rank=rank,
        n_run=int(config.nmf_runs),
        assignment_normalization=assignment_norm_arg,
    )
    nmf_clusters.index = expr_sampled.columns
    nmf_clusters["cluster"] = ["P" + str(x) for x in nmf_clusters["cluster"].astype(int)]
    keep_counts = nmf_clusters["cluster"].value_counts()
    keep_clusters = keep_counts[keep_counts > int(config.min_group_size)].index
    nmf_clusters = nmf_clusters[nmf_clusters["cluster"].isin(keep_clusters)]
    _log(f"pre-SVM NMF clusters after min_group_size>{config.min_group_size}: {nmf_clusters['cluster'].nunique()} clusters across {nmf_clusters.shape[0]} sampled cells")
    nmf_clusters.to_csv(os.path.join(snmf_dir, "icgs3_nmf_presvm_sampled_clusters.tsv"), sep="\t")

    # ICGS2 fidelity: MarkerFinder scores the FULL filtered expression matrix, not the
    # NMF guide genes. AltAnalyze2 NMF_Analysis.py:147 copies filteredInputExpFile to
    # exp.NMF-MarkerFinder.txt and ICGS_NMF.py:1027 runs MarkerFinder on that file.
    # Restricting the pool to the ~2.8k NMF features starves the cluster-fitness gate
    # (ICGS_NMF.py:639 requires >=2 markers above rho) and rejects clusters that ICGS2
    # would keep. Set --markerfinder-nmf-genes-only to restore the narrow pool.
    # Applies to EVERY modality: MarkerFinder scores the full measured feature space,
    # never the NMF guide features. adata is already reduced to that space upstream --
    # for RNA that reduction additionally keeps only protein-coding features and drops
    # RPL/RPS, MT-, HLA, XIST/TSIX and *Y (legacy_rnaseq_gene_filter); for ADT, PSI,
    # GRN, metabolite and lipid the panel is used whole with no such filter.
    use_all_genes = bool(config.markerfinder_all_genes)
    if use_all_genes:
        marker_pool = _adata_to_dense_df(sampled, sampled.var_names)
        extra = (
            " -- protein-coding, minus RPL/RPS, MT-, HLA, XIST/TSIX and *Y"
            if config.modality.lower() == "rna"
            else " -- full feature panel (no protein-coding filter for this modality)"
        )
        _log(
            f"MarkerFinder feature pool (ICGS2-faithful): {marker_pool.shape[0]} features{extra} "
            f"(NMF guide-feature pool was {expr_sampled.shape[0]})"
        )
    else:
        marker_pool = expr_sampled
        _log(
            f"MarkerFinder feature pool restricted to {marker_pool.shape[0]} NMF guide features "
            "(--markerfinder-nmf-genes-only)"
        )

    markers_all, markers_top, heat = marker_finder_wrapper(
        input_df=marker_pool[nmf_clusters.index].T,
        groups=nmf_clusters,
        top_n=config.marker_top_n,
        rho_threshold=config.marker_rho,
        marker_finder_rho=config.marker_rho,
        min_markers_per_cluster=config.marker_min_per_cluster,
    )
    marker_counts = markers_top["top_cluster"].value_counts()
    robust = marker_counts[marker_counts >= int(config.marker_min_per_cluster)].index
    nmf_robust = nmf_clusters[nmf_clusters["cluster"].isin(robust)]
    markers_top = markers_top[markers_top["top_cluster"].isin(robust)]
    if nmf_robust["cluster"].nunique() < 2:
        raise ValueError("ICGS3 found fewer than two marker-robust NMF clusters. Lower marker thresholds or inspect input normalization.")

    # SVM still trains on marker genes only (ICGS2 NMF-SVM behaviour); those markers are
    # now drawn from the full pool, so take their values from that pool.
    expr_markers_sampled = marker_pool.loc[markers_top["marker"].drop_duplicates(), nmf_robust.index]
    del marker_pool
    centroids = generate_train_data(expr_markers_sampled, nmf_robust).dropna(axis=0)
    svm_target = full if bool(config.svm_reclassify_all_cells) else sampled
    target_label = "all post-QC cells" if bool(config.svm_reclassify_all_cells) else "PageRank/Louvain sampled cells only"
    _log(f"SVM reclassification target: {target_label} ({svm_target.n_obs} cells)")
    if bool(config.svm_reclassify_all_cells):
        final_clusters = classify_adata_with_scores_chunked(
            train=centroids,
            adata=svm_target,
            genes=centroids.index,
            cluster_key=config.cluster_key,
            min_decision_score=float(config.svm_min_decision_score),
            chunk_size=int(config.svm_chunk_size),
        )
    else:
        expr_sampled_markers_for_svm = _adata_to_dense_df(svm_target, centroids.index)
        final_clusters = classify_with_scores(
            train=centroids,
            expression=expr_sampled_markers_for_svm,
            cluster_key=config.cluster_key,
            min_decision_score=float(config.svm_min_decision_score),
        )
    final_clusters.to_csv(os.path.join(snmf_dir, "icgs3_svm_reclassification_scores.tsv"), sep="\t")

    # Second MarkerFinder pass, on the same ICGS2-faithful pool. Previously this scored
    # only centroids.index (the first-pass markers), an even narrower pool than the NMF
    # guide genes, so the final cluster-fitness gate could never rescue a cluster whose
    # markers lie outside that set.
    post_svm_pool = svm_target.var_names if use_all_genes else centroids.index
    expr_full_markers = _adata_to_dense_df(svm_target, post_svm_pool, cells=final_clusters.index)
    _log(f"post-SVM MarkerFinder gene pool: {expr_full_markers.shape[0]} genes")
    markers_all2, markers_top2, heat2 = marker_finder_wrapper(
        input_df=expr_full_markers[final_clusters.index].T,
        groups=pd.DataFrame({"cluster": final_clusters[config.cluster_key]}, index=final_clusters.index),
        top_n=config.marker_top_n,
        rho_threshold=config.marker_rho,
        marker_finder_rho=config.marker_rho,
        min_markers_per_cluster=config.marker_min_per_cluster,
    )
    # ICGS2 fidelity FIX (2026-08-18). This pass selects MARKER GENES; it must not delete cells.
    # ICGS_NMF.py:1134-1141 runs markerFinder.analyzeData after Classify, then
    #   markergrps, markerlst = sortFile(allgenesfile, rho_cutoff, name)
    #   ### To plot the heatmap, use the MarkerFinder genes (function pulls those genes out)
    #   ExpandSampleClusters.filterRows(EventAnnot, ...-markers.txt, filterDB=markerlst)
    # filterRows filters ROWS of the expression matrix, that is genes. `name`, the cluster set,
    # is never revised after the SVM, and no cell is ever dropped there.
    # This code previously ran
    #     final_clusters = final_clusters[final_clusters[cluster_key].isin(robust2)]
    # on a table indexed by BARCODE, so a cluster failing the second marker pass took all of its
    # cells with it. On a mouse thymus run that deleted 4,575 of 18,388 cells (24.9%) in two
    # clusters, P48 (4,375) and P16 (200). The cluster-fitness gate that ICGS2 does apply already
    # ran BEFORE the SVM (see `robust` above), where it costs no cell because the SVM reassigns.
    # marker_finder_wrapper already enforces min_markers_per_cluster on markers_top2
    # (markerFinder.py:115), so the old markers_top2 filter below was a no-op.
    counts2 = markers_top2["top_cluster"].value_counts()
    starved = [
        str(c) for c in pd.unique(final_clusters[config.cluster_key].astype(str))
        if int(counts2.get(c, 0)) < int(config.marker_min_per_cluster)
    ]
    if starved:
        n_cells = int(final_clusters[config.cluster_key].astype(str).isin(starved).sum())
        _log(
            f"post-SVM marker pass: {len(starved)} cluster(s) own fewer than "
            f"{config.marker_min_per_cluster} markers at rho >= {config.marker_rho} "
            f"({', '.join(starved[:8])}). ICGS2 keeps them and their {n_cells} cells; they "
            f"contribute no marker row to the heatmap."
        )
    final_clusters, markers_top2, markers_all2 = order_and_renumber_final_clusters(
        final_clusters,
        markers_top2,
        markers_all2,
        expr_full_markers,
        config.cluster_key,
    )
    final_clusters.to_csv(os.path.join(snmf_dir, "icgs3_svm_reclassification_scores.final.tsv"), sep="\t")
    markers_all2.to_csv(os.path.join(marker_dir, "icgs3_markers_all_correlations.tsv"), sep="\t")
    markers_top2.to_csv(os.path.join(marker_dir, "icgs3_markers.tsv"), sep="\t", index=False)
    heat2.to_csv(os.path.join(marker_dir, "icgs3_marker_heatmap_altanalyze_format.tsv"), sep="\t")
    centroids.to_csv(os.path.join(snmf_dir, "icgs3_svm_centroids.tsv"), sep="\t")
    return final_clusters, markers_top2, heat2


def order_and_renumber_final_clusters(
    final_clusters: pd.DataFrame,
    markers_top: pd.DataFrame,
    markers_all: pd.DataFrame,
    expr_full_markers: pd.DataFrame,
    cluster_key: str,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Order surviving clusters by centroid similarity and rename C1..Cn."""
    old_clusters = pd.Index(pd.unique(final_clusters[cluster_key].astype(str)))
    if len(old_clusters) <= 1:
        mapping = {str(c): "C1" for c in old_clusters}
    else:
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import pdist

        centroids = []
        for cluster in old_clusters:
            cells = final_clusters.index[final_clusters[cluster_key].astype(str) == str(cluster)]
            centroids.append(expr_full_markers.loc[:, cells].mean(axis=1).to_numpy(dtype=float))
        C = np.vstack(centroids)
        try:
            order = leaves_list(linkage(pdist(C, metric="correlation"), method="average"))
        except Exception:
            order = np.arange(len(old_clusters))
        ordered = [str(old_clusters[i]) for i in order]
        mapping = {old: f"C{i + 1}" for i, old in enumerate(ordered)}
    final_clusters = final_clusters.copy()
    final_clusters["ICGS3_original_NMF_cluster"] = final_clusters[cluster_key].astype(str)
    cluster_order = [f"C{i + 1}" for i in range(len(mapping))]
    final_clusters[cluster_key] = final_clusters[cluster_key].astype(str).map(mapping)
    final_clusters[cluster_key] = pd.Categorical(final_clusters[cluster_key], categories=cluster_order, ordered=True)
    final_clusters = final_clusters.sort_values([cluster_key, "svm_score"], ascending=[True, False], kind="mergesort")
    final_clusters[cluster_key] = final_clusters[cluster_key].astype(str)
    markers_top = markers_top.copy()
    markers_top.index = pd.RangeIndex(len(markers_top))
    markers_top["original_top_cluster"] = markers_top["top_cluster"].astype(str)
    markers_top["top_cluster"] = markers_top["top_cluster"].astype(str).map(mapping)
    markers_top = markers_top.dropna(subset=["top_cluster"])
    markers_top["top_cluster"] = pd.Categorical(markers_top["top_cluster"], categories=cluster_order, ordered=True)
    score_col = "rho" if "rho" in markers_top.columns else ("correlation" if "correlation" in markers_top.columns else None)
    sort_cols = ["top_cluster"] + ([score_col] if score_col else []) + (["marker"] if "marker" in markers_top.columns else [])
    ascending = [True] + ([False] if score_col else []) + ([True] if "marker" in markers_top.columns else [])
    markers_top = markers_top.sort_values(sort_cols, ascending=ascending, kind="mergesort")
    markers_top["top_cluster"] = markers_top["top_cluster"].astype(str)
    markers_all = markers_all.copy()
    markers_all.index = pd.RangeIndex(len(markers_all))
    markers_all["original_top_cluster"] = markers_all["top_cluster"].astype(str)
    markers_all["top_cluster"] = markers_all["top_cluster"].astype(str).map(mapping)
    markers_all = markers_all.dropna(subset=["top_cluster"])
    markers_all["top_cluster"] = pd.Categorical(markers_all["top_cluster"], categories=cluster_order, ordered=True)
    score_col = "rho" if "rho" in markers_all.columns else ("correlation" if "correlation" in markers_all.columns else None)
    sort_cols = ["top_cluster"] + ([score_col] if score_col else []) + (["marker"] if "marker" in markers_all.columns else [])
    ascending = [True] + ([False] if score_col else []) + ([True] if "marker" in markers_all.columns else [])
    markers_all = markers_all.sort_values(sort_cols, ascending=ascending, kind="mergesort")
    markers_all["top_cluster"] = markers_all["top_cluster"].astype(str)
    return final_clusters, markers_top, markers_all


def classify_with_scores(
    *,
    train: pd.DataFrame,
    expression: pd.DataFrame,
    cluster_key: str,
    min_decision_score: float = 0.0,
) -> pd.DataFrame:
    """Linear SVM reclassification with ICGS2-style positive-margin filtering.

    Returns all cells whose winning class decision score is greater than the
    configured threshold. The row order is cluster then descending score, which
    the canonical MarkerFinder heatmap preserves within each cluster.
    """
    from sklearn.svm import LinearSVC

    clf = LinearSVC(random_state=42, dual=False, max_iter=50000)
    clf.fit(train.T, y=train.columns)
    scores = clf.decision_function(expression.T)
    if np.ndim(scores) == 1:
        pred = np.where(scores > 0, clf.classes_[1], clf.classes_[0])
        win_score = np.abs(scores)
        raw_margin = scores
    else:
        score_array = np.asarray(scores)
        best = np.argmax(score_array, axis=1)
        pred = np.asarray(clf.classes_)[best]
        win_score = score_array[np.arange(score_array.shape[0]), best]
        raw_margin = win_score
    pred = pd.Series(pred, index=expression.columns, dtype=str).str.replace("^U", "", regex=True)
    win_score = pd.Series(np.asarray(win_score, dtype=float), index=expression.columns)
    raw_margin = pd.Series(np.asarray(raw_margin, dtype=float), index=expression.columns)
    keep = win_score > float(min_decision_score)
    n_drop = int((~keep).sum())
    if n_drop:
        _log(
            f"linear SVM: excluding {n_drop} cells with winning decision score <= "
            f"{float(min_decision_score):.3g}"
        )
    out = pd.DataFrame(
        {
            cluster_key: pred,
            "svm_score": win_score,
            "svm_margin": raw_margin,
        },
        index=expression.columns,
    )
    out = out.loc[keep].copy()
    out = out.sort_values([cluster_key, "svm_score"], ascending=[True, False], kind="mergesort")
    return out


def _classification_frame_from_scores(
    *,
    classes: Sequence[object],
    scores,
    index: Sequence[str],
    cluster_key: str,
    min_decision_score: float,
) -> pd.DataFrame:
    if np.ndim(scores) == 1:
        score_values = np.asarray(scores, dtype=float)
        pred = np.where(score_values > 0, classes[1], classes[0])
        win_score = np.abs(score_values)
        raw_margin = score_values
    else:
        score_array = np.asarray(scores, dtype=float)
        best = np.argmax(score_array, axis=1)
        pred = np.asarray(classes)[best]
        win_score = score_array[np.arange(score_array.shape[0]), best]
        raw_margin = win_score
    pred = pd.Series(pred, index=index, dtype=str).str.replace("^U", "", regex=True)
    win_score = pd.Series(np.asarray(win_score, dtype=float), index=index)
    raw_margin = pd.Series(np.asarray(raw_margin, dtype=float), index=index)
    keep = win_score > float(min_decision_score)
    return pd.DataFrame(
        {
            cluster_key: pred.loc[keep],
            "svm_score": win_score.loc[keep],
            "svm_margin": raw_margin.loc[keep],
        },
        index=win_score.index[keep],
    )


def classify_adata_with_scores_chunked(
    *,
    train: pd.DataFrame,
    adata: ad.AnnData,
    genes: Sequence[str],
    cluster_key: str,
    min_decision_score: float = 0.0,
    chunk_size: int = 50000,
) -> pd.DataFrame:
    """Linear SVM reclassification without materializing all cells as one dense matrix."""
    from sklearn.svm import LinearSVC

    gene_index = adata.var_names.get_indexer(pd.Index(genes).astype(str))
    valid = gene_index >= 0
    if not np.any(valid):
        raise ValueError("Could not align SVM centroid genes to AnnData features.")
    train_aligned = train.loc[pd.Index(genes)[valid].astype(str)]
    clf = LinearSVC(random_state=42, dual=False, max_iter=50000)
    clf.fit(train_aligned.T, y=train_aligned.columns)
    idx = gene_index[valid]
    chunk_size = max(1, int(chunk_size))
    frames = []
    dropped = 0
    for start in range(0, adata.n_obs, chunk_size):
        end = min(start + chunk_size, adata.n_obs)
        X = adata.X[start:end, :][:, idx]
        values = X.toarray().astype(np.float32, copy=False) if sp.issparse(X) else np.asarray(X, dtype=np.float32)
        names = adata.obs_names[start:end].astype(str).tolist()
        chunk = pd.DataFrame(values, index=names, columns=train_aligned.index)
        scores = clf.decision_function(chunk)
        frame = _classification_frame_from_scores(
            classes=clf.classes_,
            scores=scores,
            index=names,
            cluster_key=cluster_key,
            min_decision_score=min_decision_score,
        )
        dropped += (end - start) - frame.shape[0]
        frames.append(frame)
    if dropped:
        _log(
            f"linear SVM: excluding {dropped} cells with winning decision score <= "
            f"{float(min_decision_score):.3g}"
        )
    if not frames:
        return pd.DataFrame(columns=[cluster_key, "svm_score", "svm_margin"])
    out = pd.concat(frames, axis=0)
    out = out.sort_values([cluster_key, "svm_score"], ascending=[True, False], kind="mergesort")
    return out


def _import_umap_with_local_retry():
    try:
        import umap

        return umap
    except Exception:
        # Some local installs import ParametricUMAP/TensorFlow from umap.__init__.
        # Blocking tensorflow reproduces the legacy AltAnalyze fallback intent:
        # use standard UMAP without optional parametric extras.
        for mod in list(sys.modules):
            if mod == "umap" or mod.startswith("umap."):
                del sys.modules[mod]
        prior_tf = sys.modules.get("tensorflow", "__missing__")
        sys.modules["tensorflow"] = None
        try:
            import umap

            return umap
        finally:
            if prior_tf == "__missing__":
                sys.modules.pop("tensorflow", None)
            else:
                sys.modules["tensorflow"] = prior_tf


def compute_umap_outputs(
    adata: ad.AnnData,
    config: ICGS3Config,
    outdir: str,
    *,
    marker_features: Optional[Sequence[str]] = None,
) -> None:
    umap_dir = os.path.join(outdir, "UMAPs")
    os.makedirs(umap_dir, exist_ok=True)
    if config.generate_umap:
        feature_mode = str(config.umap_feature_mode or "markerfinder").lower()
        if feature_mode == "variable":
            if "icgs3_feature" in adata.var:
                source_features = adata.var_names[adata.var["icgs3_feature"].values]
            else:
                source_features = adata.var_names
            source_label = "variable_features"
        elif feature_mode == "pca":
            source_features = []
            source_label = "pca"
        else:
            source_features = marker_features or []
            source_label = "final_markerfinder_features"
        features = pd.Index([str(g) for g in source_features]).drop_duplicates()
        features = features[features.isin(adata.var_names.astype(str))]
        if _batch_correction_applies(config, "umap") and str(config.batch_correction or "none").lower() == "harmony":
            graph_source = adata[:, features].copy() if len(features) > 0 else adata.copy()
            _log(
                f"running Harmony-corrected UMAP on {graph_source.n_vars} {source_label} "
                f"via PCA/neighbors (batch key={config.batch_key})"
            )
            graph = compute_sparse_graph(graph_source, config)
            _copy_graph_slots(graph, adata)
            try:
                umap = _import_umap_with_local_retry()
                coords = graph.obsm["X_pca_harmony"] if "X_pca_harmony" in graph.obsm else graph.obsm["X_pca"]
                n_neighbors = min(50, max(2, graph.n_obs - 1))
                model = umap.UMAP(
                    n_neighbors=n_neighbors,
                    min_dist=0.75,
                    metric="correlation",
                    random_state=config.random_state,
                )
                adata.obsm["X_umap"] = model.fit_transform(coords)
                adata.uns["icgs3_umap"] = {
                    "source": f"harmony_corrected_{source_label}_pcs",
                    "n_features": int(graph_source.n_vars),
                    "n_pcs": int(coords.shape[1]),
                    "batch_key": _batch_keys(config),
                    "n_neighbors": int(n_neighbors),
                    "min_dist": 0.75,
                    "metric": "correlation",
                }
                if len(features) > 0:
                    pd.DataFrame({"feature": features}).to_csv(
                        os.path.join(umap_dir, "icgs3_umap_features.tsv"),
                        sep="\t",
                        index=False,
                    )
            except Exception as exc:
                with open(os.path.join(umap_dir, "icgs3_umap_error.txt"), "w", encoding="utf-8") as handle:
                    handle.write(str(exc) + "\n")
                _log(f"UMAP skipped: {exc}")
        elif len(features) == 0:
            _log(f"UMAP using graph/PCA fallback (requested mode={feature_mode})")
            graph = compute_sparse_graph(adata, config)
            _copy_graph_slots(graph, adata)
            try:
                sc.tl.umap(graph, random_state=config.random_state)
                adata.obsm["X_umap"] = graph.obsm["X_umap"].copy()
            except Exception as exc:
                _log(f"Scanpy UMAP skipped: {exc}")
        else:
            feature_path = os.path.join(umap_dir, "icgs3_umap_features.tsv")
            pd.DataFrame({"feature": features}).to_csv(feature_path, sep="\t", index=False)
            X = adata[:, features].X
            X = X.toarray() if sp.issparse(X) else np.asarray(X)
            n_neighbors = min(50, max(2, adata.n_obs - 1))
            _log(f"running final UMAP on {len(features)} {source_label} (metric=correlation, n_neighbors={n_neighbors}, min_dist=0.75)")
            try:
                umap = _import_umap_with_local_retry()
                model = umap.UMAP(
                    n_neighbors=n_neighbors,
                    min_dist=0.75,
                    metric="correlation",
                    random_state=config.random_state,
                )
                adata.obsm["X_umap"] = model.fit_transform(X)
                adata.uns["icgs3_umap"] = {
                    "source": source_label,
                    "n_features": int(len(features)),
                    "feature_file": feature_path,
                    "n_neighbors": int(n_neighbors),
                    "min_dist": 0.75,
                    "metric": "correlation",
                }
            except Exception as exc:
                with open(os.path.join(umap_dir, "icgs3_umap_error.txt"), "w", encoding="utf-8") as handle:
                    handle.write(str(exc) + "\n")
                _log(f"UMAP skipped: {exc}")
    if "X_umap" in adata.obsm:
        umap = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names, columns=["UMAP_1", "UMAP_2"])
        covariate_cols = [c.strip() for c in str(config.umap_covariates or "").split(",") if c.strip()]
        for col in list(dict.fromkeys(["sample", "Library", config.cluster_key, "ICGS3_cell_state_prediction"] + covariate_cols)):
            if col in adata.obs:
                umap[col] = adata.obs[col].astype(str).values
        umap.rename_axis("barcode").reset_index().to_csv(os.path.join(umap_dir, "icgs3_umap.tsv"), sep="\t", index=False)


def write_umap_plots(adata: ad.AnnData, config: ICGS3Config, outdir: str) -> None:
    if "X_umap" not in adata.obsm:
        return
    plt.rcParams.update(
        {
            "backend": "Agg",
            "axes.linewidth": 0.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )
    umap_dir = os.path.join(outdir, "UMAPs")
    restricted_dir = os.path.join(umap_dir, "sample-specific")
    os.makedirs(umap_dir, exist_ok=True)
    os.makedirs(restricted_dir, exist_ok=True)
    coords_all = np.asarray(adata.obsm["X_umap"], dtype=float)
    finite = np.isfinite(coords_all).all(axis=1)
    if np.any(finite):
        x_min, x_max = np.percentile(coords_all[finite, 0], [0.2, 99.8])
        y_min, y_max = np.percentile(coords_all[finite, 1], [0.2, 99.8])
    else:
        x_min, x_max, y_min, y_max = 0.0, 1.0, 0.0, 1.0
    x_pad = max((x_max - x_min) * 0.06, 1e-6)
    y_pad = max((y_max - y_min) * 0.06, 1e-6)
    global_xlim = (float(x_min - x_pad), float(x_max + x_pad))
    global_ylim = (float(y_min - y_pad), float(y_max + y_pad))

    def safe_name(value: str) -> str:
        return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in str(value))

    def ordered_categories(values: pd.Series) -> List[str]:
        if pd.api.types.is_categorical_dtype(values):
            present = set(values.astype(str))
            return [str(c) for c in values.cat.categories if str(c) in present]
        return [str(c) for c in pd.unique(values.astype(str))]

    def point_size(n: int) -> float:
        if n > 60000:
            return 6.0
        if n > 30000:
            return 8.0
        if n > 10000:
            return 10.0
        return 14.0

    def cluster_palette(categories: Sequence[str]) -> dict:
        cmap = plt.get_cmap("tab20", max(1, len(categories)))
        return {str(cat): cmap(i) for i, cat in enumerate(categories)}

    def draw_categorical(
        column: str,
        path: str,
        *,
        mask: Optional[np.ndarray] = None,
        label_centers: bool = False,
        title: Optional[str] = None,
        color_by: Optional[str] = None,
    ) -> None:
        if column not in adata.obs:
            return
        color_col = color_by or column
        if color_col not in adata.obs:
            return
        row_mask = np.ones(adata.n_obs, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
        coords = coords_all[row_mask]
        values = adata.obs[color_col].astype(str).iloc[np.where(row_mask)[0]]
        categories = ordered_categories(adata.obs[color_col])
        palette = cluster_palette(categories)
        fig = plt.figure(figsize=(8.0, 8.4))
        ax = fig.add_axes([0.10, 0.22, 0.86, 0.70])
        for cat in categories:
            cat_mask = values.values == str(cat)
            if not np.any(cat_mask):
                continue
            ax.scatter(
                coords[cat_mask, 0],
                coords[cat_mask, 1],
                s=point_size(coords.shape[0]),
                marker="o",
                linewidths=0,
                edgecolors="none",
                color=palette[str(cat)],
                alpha=1.0,
                label=str(cat),
            )
        if label_centers:
            for cat in categories:
                cat_mask = values.values == str(cat)
                if not np.any(cat_mask):
                    continue
                center = np.median(coords[cat_mask], axis=0)
                ax.text(center[0], center[1], str(cat), ha="center", va="center", fontsize=8, weight="bold")
        ax.set_xlabel("UMAP_1")
        ax.set_ylabel("UMAP_2")
        ax.set_title(title or column)
        ax.set_xlim(global_xlim)
        ax.set_ylim(global_ylim)
        ax.set_aspect("equal", adjustable="box")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if len(categories) <= 30:
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                fig.legend(
                    handles,
                    labels,
                    loc="lower center",
                    ncol=max(1, min(5, len(labels))),
                    frameon=False,
                    fontsize=7,
                    markerscale=2.25,
                )
        fig.savefig(path, dpi=300)
        plt.close(fig)

    #: Continuous gene-expression UMAP. The layering is the point: cells are sorted by
    #: expression and drawn lowest first, so the highest-expressing cells land on top and
    #: are never hidden by the non-expressing background. Drawing in the object's own row
    #: order buries the signal under whichever cells happen to come last.
    #:
    #: Rendering defaults, fixed here so every ICGS3 run matches:
    #:   ramp        white #FFFFFF to red #FF0000
    #:   dot size    half the categorical diameter, 2.1107 points
    #:   outline     none
    #:   raster      none; circles are inline vector ellipses and the colour ramp is
    #:               256 filled rectangles, because fig.colorbar emits a raster image
    GENE_UMAP_RAMP = ("#FFFFFF", "#FF0000")
    GENE_UMAP_DOT_SCALE = 0.5
    GENE_UMAP_RAMP_STEPS = 256

    def draw_continuous_gene(gene: str, path: str, *, title: Optional[str] = None) -> bool:
        from matplotlib.collections import EllipseCollection
        from matplotlib.colors import LinearSegmentedColormap, Normalize
        from matplotlib.patches import Rectangle

        if gene not in adata.var_names:
            _log(f"UMAP gene overlay skipped; {gene} is absent from var_names")
            return False
        column = adata[:, gene].X
        values = np.asarray(column.todense()).ravel() if hasattr(column, "todense") else np.asarray(column).ravel()
        values = np.asarray(values, dtype=float)

        cmap = LinearSegmentedColormap.from_list("icgs3_white_red", list(GENE_UMAP_RAMP))
        vmax = float(np.nanmax(values)) if np.isfinite(values).any() and np.nanmax(values) > 0 else 1.0
        norm = Normalize(vmin=0.0, vmax=vmax)

        order = np.argsort(values, kind="mergesort")   # lowest first, highest drawn last
        diameter = float(2.0 * np.sqrt(point_size(coords_all.shape[0]) / np.pi)) * GENE_UMAP_DOT_SCALE

        fig = plt.figure(figsize=(8.0, 8.4))
        ax = fig.add_axes([0.10, 0.22, 0.86, 0.70])
        ax.add_collection(EllipseCollection(
            widths=diameter, heights=diameter, angles=0.0, units="points",
            offsets=coords_all[order], offset_transform=ax.transData,
            facecolors=cmap(norm(values[order])), edgecolors="none",
            linewidths=0.0, rasterized=False,
        ))
        ax.set_xlabel("UMAP_1")
        ax.set_ylabel("UMAP_2")
        ax.set_title(title or f"{gene} expression")
        ax.set_xlim(global_xlim)
        ax.set_ylim(global_ylim)
        ax.set_aspect("equal", adjustable="box")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        cax = fig.add_axes([0.35, 0.11, 0.30, 0.022])
        edges = np.linspace(0.0, vmax, GENE_UMAP_RAMP_STEPS + 1)
        for i in range(GENE_UMAP_RAMP_STEPS):
            cax.add_patch(Rectangle(
                (edges[i], 0.0), edges[i + 1] - edges[i], 1.0,
                facecolor=cmap(norm(0.5 * (edges[i] + edges[i + 1]))),
                edgecolor="none", linewidth=0.0, rasterized=False,
            ))
        cax.set_xlim(0.0, vmax)
        cax.set_ylim(0.0, 1.0)
        cax.set_yticks([])
        cax.tick_params(labelsize=7, length=2, width=0.5)
        for side in ("top", "right", "left", "bottom"):
            cax.spines[side].set_linewidth(0.5)
            cax.spines[side].set_edgecolor("#000000")
        cax.set_xlabel(gene, fontsize=8)

        fig.savefig(path)
        plt.close(fig)
        _log(f"UMAP gene overlay {gene}: {int((values > 0).sum())} of {values.size} cells above zero -> {path}")
        return True

    draw_categorical(config.cluster_key, os.path.join(umap_dir, "icgs3_umap_clusters.pdf"), label_centers=True, title="ICGS3 clusters")
    draw_categorical(
        "ICGS3_cell_state_prediction",
        os.path.join(umap_dir, "icgs3_umap_cell_state_predictions.pdf"),
        label_centers=True,
        title="ICGS3 cell-state predictions",
    )
    genes = [g.strip() for g in str(config.umap_genes or "").split(",") if g.strip()]
    for gene in dict.fromkeys(genes):
        draw_continuous_gene(gene, os.path.join(umap_dir, f"icgs3_umap_gene_{safe_name(gene)}.pdf"))

    covariates = [c.strip() for c in str(config.umap_covariates or "").split(",") if c.strip()]
    for covariate in covariates:
        draw_categorical(covariate, os.path.join(umap_dir, f"icgs3_umap_{safe_name(covariate)}.pdf"), title=covariate)
        if covariate in adata.obs:
            for level in ordered_categories(adata.obs[covariate]):
                subset = adata.obs[covariate].astype(str).values == str(level)
                if not np.any(subset):
                    continue
                draw_categorical(
                    covariate,
                    os.path.join(restricted_dir, f"icgs3_umap_{safe_name(covariate)}_{safe_name(level)}_clusters.pdf"),
                    mask=subset,
                    color_by=config.cluster_key,
                    label_centers=True,
                    title=f"{covariate}: {level}",
                )


def _heatmap_covariates(config: ICGS3Config) -> List[str]:
    requested = [c.strip() for c in str(config.heatmap_covariates or "").split(",") if c.strip()]
    return list(dict.fromkeys(requested))


def _load_goelite_heatmap_terms(outdir: str, config: ICGS3Config) -> dict:
    if not config.heatmap_goelite_terms:
        return {}
    path = os.path.join(outdir, "GO-Elite", "icgs3_biomarker_enrichment.tsv")
    if not os.path.exists(path):
        _log(f"heatmap GO-Elite terms requested but enrichment file is missing: {path}")
        return {}
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as exc:
        _log(f"heatmap GO-Elite terms skipped; failed reading {path}: {exc}")
        return {}
    needed = {"cluster", "term_name", "p_value"}
    if df.empty or not needed.issubset(df.columns):
        _log("heatmap GO-Elite terms skipped; enrichment table lacks cluster/term_name/p_value columns")
        return {}
    df = df.copy()
    df["cluster"] = df["cluster"].astype(str)
    df["term_name"] = df["term_name"].astype(str)
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    if "fdr" in df.columns:
        df["sort_fdr"] = pd.to_numeric(df["fdr"], errors="coerce")
    else:
        df["sort_fdr"] = np.nan
    df = df.dropna(subset=["p_value"]).sort_values(["cluster", "sort_fdr", "p_value", "term_name"], na_position="last")
    terms = {}
    max_terms = max(1, int(config.heatmap_goelite_max_terms or 30))
    for cluster, sub in df.groupby("cluster", sort=False):
        rows = []
        seen = set()
        for _, row in sub.iterrows():
            term = str(row["term_name"]).strip()
            if not term or term in seen:
                continue
            seen.add(term)
            rows.append((term, float(row["p_value"])))
            if len(rows) >= max_terms:
                break
        if rows:
            terms[str(cluster)] = rows
    _log(f"heatmap GO-Elite term labels: {sum(len(v) for v in terms.values())} terms across {len(terms)} clusters")
    return terms


def run_canonical_heatmap(adata: ad.AnnData, config: ICGS3Config, outdir: str) -> dict:
    from altanalyze3.components.visualization.marker_heatmap_h5ad import generate_marker_heatmap_from_adata

    marker_dir = os.path.join(outdir, "MarkerFinder")
    os.makedirs(marker_dir, exist_ok=True)
    covariates = [c for c in _heatmap_covariates(config) if c in adata.obs]
    missing_covariates = [c for c in _heatmap_covariates(config) if c not in adata.obs]
    if missing_covariates:
        _log(f"heatmap covariates not found in adata.obs: {', '.join(missing_covariates)}")
    go_terms = _load_goelite_heatmap_terms(outdir, config)
    return generate_marker_heatmap_from_adata(
        adata,
        cluster_key=config.cluster_key,
        out=os.path.join(marker_dir, "icgs3_marker_heatmap.pdf"),
        top_n=max(1, int(config.marker_top_n)),
        markers_tsv=os.path.join(marker_dir, "icgs3_marker_heatmap_markers.tsv"),
        heatmap_tsv=os.path.join(marker_dir, "icgs3_marker_heatmap_fold_matrix.tsv"),
        heatmap_column_tsv=os.path.join(marker_dir, "icgs3_marker_heatmap_exp_matrix.tsv") if config.write_heatmap_expression_tsv else None,
        heatmap_cache=os.path.join(marker_dir, "icgs3_marker_heatmap_fold_matrix.npz"),
        marker_method="markerfinder",
        cells_per_cluster=config.heatmap_cells_per_cluster,
        seed=config.random_state,
        species={"Hs": "human", "Mm": "mouse"}.get(config.species, config.species),
        write_expression_tsv=config.write_heatmap_expression_tsv,
        covariate_columns=covariates,
        go_terms=go_terms,
        go_terms_max=config.heatmap_goelite_max_terms,
    )


def _default_biomarker_file(species: str) -> Optional[str]:
    base = "/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart72/goelite"
    candidates = {
        "Hs": os.path.join(base, "Hs", "gene-mapp", "Ensembl-BioMarkers.txt"),
        "Mm": os.path.join(base, "Mm", "gene-mapp", "Ensembl-BioMarkers.txt"),
    }
    path = candidates.get(species)
    return path if path and os.path.exists(path) else None


def clean_biomarker_prediction_labels(labels: pd.DataFrame) -> pd.DataFrame:
    """Apply the legacy RNASeq.py BioMarkers label cleanup rules."""
    if labels.empty or "term_name" not in labels:
        return labels
    labels = labels.copy()
    tissue_by_cluster = {}
    for _, row in labels.iterrows():
        term = str(row["term_name"])
        if any(x in term for x in ("Adult", "Fetal", "Embryo", "Embryonic", "Term")):
            parts = term.split()
            if len(parts) > 1:
                tissue_by_cluster[str(row["cluster"])] = parts[1]
    tissue_preference = "NONE"
    if tissue_by_cluster:
        counts = pd.Series(list(tissue_by_cluster.values())).value_counts()
        if not counts.empty and float(counts.iloc[0]) / max(len(tissue_by_cluster), 1) > 0.33:
            tissue_preference = str(counts.index[0])
            _log(f"likely tissue based on GO-Elite BioMarkers = {tissue_preference}")

    cleaned = []
    for _, row in labels.iterrows():
        cluster = str(row["cluster"])
        original = str(row["term_name"])
        label = original.replace("/", "-")
        for token in (" (", "("):
            if token in label:
                label = label.split(token)[0]
        for token in ("Adult ", "Fetal ", "Embryonic", "Embryo", "Term"):
            label = label.replace(token, "")
        if tissue_preference != "NONE" and tissue_preference != original:
            label = label.replace(tissue_preference, "")
        label = label.strip(" -_")
        if not label:
            label = original.strip() or "UNK"
        cleaned.append(f"{label}_c{cluster}")
    labels["cell_type_prediction"] = cleaned
    return labels


def biomarker_enrichment(markers: pd.DataFrame, background: Sequence[str], config: ICGS3Config, outdir: str) -> pd.DataFrame:
    goelite_dir = os.path.join(outdir, "GO-Elite")
    os.makedirs(goelite_dir, exist_ok=True)
    path = config.biomarker_file or _default_biomarker_file(config.species)
    if not path or not os.path.exists(path) or markers.empty:
        return pd.DataFrame()
    bm = pd.read_csv(path, sep="\t")
    gene_col = "Gene" if "Gene" in bm.columns else ("System" if "System" in bm.columns else bm.columns[1])
    term_col = "Term" if "Term" in bm.columns else ("GeneSet" if "GeneSet" in bm.columns else bm.columns[-1])
    bm = bm[[gene_col, term_col]].dropna()
    bm.columns = ["gene", "term"]
    bm["gene"] = bm["gene"].astype(str).str.upper()
    bg = {str(g).upper() for g in background}
    term_sets = {t: set(g["gene"]) & bg for t, g in bm.groupby("term")}
    term_sets = {t: s for t, s in term_sets.items() if len(s) >= 3}
    if not term_sets or not bg:
        return pd.DataFrame()
    rows = []
    from scipy.stats import hypergeom

    for cluster, grp in markers.groupby("top_cluster"):
        query = {str(g).upper() for g in grp["marker"]} & bg
        if not query:
            continue
        for term, genes in term_sets.items():
            overlap = query & genes
            if not overlap:
                continue
            p = hypergeom(M=len(bg), n=len(genes), N=len(query)).sf(len(overlap) - 1)
            rows.append(
                {
                    "cluster": cluster,
                    "term_name": term,
                    "p_value": p,
                    "overlap": len(overlap),
                    "query_size": len(query),
                    "term_size": len(genes),
                    "overlap_genes": ",".join(sorted(overlap)),
                }
            )
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    out["fdr"] = multipletests(out["p_value"].values, method="fdr_bh")[1]
    out = out.sort_values(["cluster", "fdr", "p_value", "term_name"])
    out.to_csv(os.path.join(goelite_dir, "icgs3_biomarker_enrichment.tsv"), sep="\t", index=False)
    labels = out.groupby("cluster").head(1)[["cluster", "term_name", "fdr", "overlap"]]
    labels = clean_biomarker_prediction_labels(labels)
    labels.to_csv(os.path.join(goelite_dir, "icgs3_cell_state_predictions.tsv"), sep="\t", index=False)
    return labels


def run_icgs3(config: ICGS3Config) -> ICGS3Result:
    config = apply_modality_defaults(config)
    outdir = os.path.abspath(config.output_dir)
    os.makedirs(outdir, exist_ok=True)
    for dirname in ("logs", "sNMF", "MarkerFinder", "UMAPs", "GO-Elite"):
        os.makedirs(os.path.join(outdir, dirname), exist_ok=True)
    log_path = os.path.join(outdir, "logs", f"icgs3_{datetime.now().strftime('%Y%m%d-%H%M%S')}.log")
    original_stdout, original_stderr = sys.stdout, sys.stderr
    log_handle = open(log_path, "a", encoding="utf-8")
    sys.stdout = Tee(original_stdout, log_handle)
    sys.stderr = Tee(original_stderr, log_handle)
    start_time = time.time()
    try:
        with warnings.catch_warnings():
            # Collapse repeated numerical warnings to ONE line per unique source
            # location. The prior run wrote 33,562 duplicate RuntimeWarning lines
            # (6.1 MB of a 6.1 MB log). "once"-per-location keeps every DISTINCT
            # warning visible -- they are results (RULE 4) -- without the flood.
            warnings.simplefilter("default")
            warnings.filterwarnings("once", category=RuntimeWarning)
            warnings.filterwarnings("ignore", category=FutureWarning)
            warnings.filterwarnings("ignore", category=PendingDeprecationWarning)
            warnings.filterwarnings("ignore", message=".*n_jobs value.*", category=UserWarning)
            warnings.filterwarnings("ignore", message=".*default backend for leiden.*", category=FutureWarning)
            result = _run_icgs3_logged(config, outdir, log_path, start_time)
    finally:
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        log_handle.close()
    return result


def _run_icgs3_logged(config: ICGS3Config, outdir: str, log_path: str, start_time: float) -> ICGS3Result:
    def step_time(label: str, t0: float) -> float:
        elapsed = time.time() - t0
        _log(f"{label} completed in {elapsed:.1f}s")
        return time.time()

    _log(f"log file: {log_path}")
    _log(f"cli equivalent: {cli_equivalent(config)}")
    _log(f"parameters: {json.dumps(asdict(config), sort_keys=True)}")
    if float(config.nmf_k_multiplier) != 2.0:
        _log(
            "deprecated parameter ignored: --nmf-k-multiplier is retained for CLI compatibility, "
            "but automatic rank uses the UDON/ICGS2 Tracy-Widom estimator directly"
        )
    Path(os.path.join(outdir, "icgs3_config.json")).write_text(json.dumps(asdict(config), indent=2), encoding="utf-8")

    t = time.time()
    _log(f"loading {len(config.input_paths)} input(s)")
    adata = read_inputs(config.input_paths)
    _log(f"loaded matrix: {adata.n_obs} cells x {adata.n_vars} features")
    t = step_time("input loading", t)
    adata = apply_qc(
        adata,
        min_genes=config.min_genes,
        min_cells=config.min_cells,
        min_counts=config.min_counts,
        mito_percent=config.mito_percent,
        layer=config.layer,
    )
    t = step_time("QC", t)
    adata = prepare_expression(adata, config)
    t = step_time("normalization", t)
    adata = apply_rna_unsupervised_gene_filter(adata, config, outdir=outdir)
    t = step_time("RNA unsupervised gene filtering", t)
    analysis_adata = apply_expression_batch_adjustment(adata, config, outdir)
    if analysis_adata is not adata:
        t = step_time("expression batch adjustment", t)
    if config.downsample_key or config.target_cells:
        # Optional user-facing stratified cap before graph/PageRank. This is not
        # the ICGS2 intelligent sampler; it is a memory guard for extreme inputs.
        adata_for_sampling = stratified_downsample_adata(
            analysis_adata,
            target_cells=config.target_cells,
            downsample_key=config.downsample_key,
            max_cells_per_group=config.max_cells_per_group,
            random_state=config.random_state,
        )
    else:
        adata_for_sampling = analysis_adata
    # Resume entry point: reuse a completed run's PageRank/Louvain selection instead of
    # recomputing it. The sampling is deterministic given the same input and seed, so
    # reusing it isolates downstream changes (rank, marker thresholds) from any change in
    # which cells were sampled. Costs 0 s in place of the 3,068 s that step took on the
    # 280,446-cell input.
    if config.resume_sampled_from:
        src = str(config.resume_sampled_from)
        if not os.path.exists(src):
            raise FileNotFoundError(f"--resume-sampled-from file not found: {src}")
        prior = pd.read_csv(src, sep="\t")
        if "barcode" not in prior.columns:
            raise ValueError(f"{src} lacks a 'barcode' column; expected icgs3_pagerank_downsampling.tsv")
        if "selected_final" in prior.columns:
            wanted = prior.loc[prior["selected_final"].astype(str).str.lower() == "true", "barcode"]
        else:
            wanted = prior["barcode"]
        wanted = pd.Index(wanted.astype(str).unique())
        present = wanted.intersection(pd.Index(adata_for_sampling.obs_names.astype(str)))
        missing = len(wanted) - len(present)
        if len(present) < 2:
            raise ValueError(
                f"--resume-sampled-from matched {len(present)} of {len(wanted)} barcodes in the "
                "current input; the QC/gene-filter settings must match the source run."
            )
        _log(
            f"resuming from prior sampling: {len(present)} of {len(wanted)} barcodes matched "
            f"({missing} absent from this input); PageRank/Louvain SKIPPED -- source {src}"
        )
        if missing:
            _log(f"WARNING resume barcode mismatch: {missing} of {len(wanted)} prior barcodes are not in this input")
        sampled = adata_for_sampling[present].copy()
        pagerank_scores = pd.DataFrame({"barcode": present.astype(str), "selected_final": True})
    else:
        sampled, pagerank_scores = pagerank_downsample_adata(adata_for_sampling, config)
    pagerank_scores.to_csv(os.path.join(outdir, "sNMF", "icgs3_pagerank_downsampling.tsv"), sep="\t", index=False)
    _log(f"downsampling summary: retained {sampled.n_obs} sampled cells for NMF from {adata_for_sampling.n_obs} candidates")
    write_retention_audit(
        full=analysis_adata,
        candidates=adata_for_sampling,
        sampled=sampled,
        final_clusters=None,
        config=config,
        outdir=outdir,
    )
    t = step_time("PageRank/Louvain downsampling", t)

    final_clusters, markers, _ = run_nmf_marker_svm(sampled, analysis_adata, config, outdir)
    write_retention_audit(
        full=analysis_adata,
        candidates=adata_for_sampling,
        sampled=sampled,
        final_clusters=final_clusters,
        config=config,
        outdir=outdir,
    )
    adata.obs[config.cluster_key] = pd.Series(final_clusters[config.cluster_key], index=final_clusters.index).reindex(adata.obs_names)
    adata.obs["ICGS3_SVM_score"] = pd.Series(final_clusters["svm_score"], index=final_clusters.index).reindex(adata.obs_names)
    adata.obs["ICGS3_SVM_margin"] = pd.Series(final_clusters["svm_margin"], index=final_clusters.index).reindex(adata.obs_names)
    adata.obs["ICGS3_original_NMF_cluster"] = pd.Series(final_clusters["ICGS3_original_NMF_cluster"], index=final_clusters.index).reindex(adata.obs_names)
    adata = adata[adata.obs[config.cluster_key].notna()].copy()
    ordered_barcodes = [bc for bc in final_clusters.index.tolist() if bc in set(adata.obs_names)]
    adata = adata[ordered_barcodes].copy()
    cluster_order = [f"C{i + 1}" for i in range(adata.obs[config.cluster_key].astype(str).nunique())]
    observed_order = [c for c in cluster_order if c in set(adata.obs[config.cluster_key].astype(str))]
    adata.obs[config.cluster_key] = pd.Categorical(adata.obs[config.cluster_key].astype(str), categories=observed_order, ordered=True)
    adata.uns["lineage_order"] = observed_order
    adata.obs["icgs3_marker_robust"] = True
    _log(f"marker-robust SVM clusters: {adata.obs[config.cluster_key].nunique()} clusters across {adata.n_obs} cells")
    t = step_time("NMF, MarkerFinder, and SVM reclassification", t)

    biomarker_predictions = biomarker_enrichment(markers, adata.var_names, config, outdir)
    if not biomarker_predictions.empty:
        label_col = "cell_type_prediction" if "cell_type_prediction" in biomarker_predictions.columns else "term_name"
        label_map = biomarker_predictions.set_index("cluster")[label_col].astype(str).to_dict()
        adata.obs["ICGS3_cell_state_prediction"] = [
            label_map.get(str(cluster), f"UNK-c{cluster}")
            for cluster in adata.obs[config.cluster_key].astype(str)
        ]
    else:
        adata.obs["ICGS3_cell_state_prediction"] = [
            f"UNK-c{cluster}" for cluster in adata.obs[config.cluster_key].astype(str)
        ]
    t = step_time("GO-Elite BioMarkers enrichment", t)

    if analysis_adata is not adata:
        umap_adata = analysis_adata[adata.obs_names].copy()
        umap_adata.obs = adata.obs.copy()
        compute_umap_outputs(umap_adata, config, outdir, marker_features=markers["marker"].tolist())
        for key in ("X_umap", "X_pca", "X_pca_harmony", "X_pca_uncorrected"):
            if key in umap_adata.obsm:
                adata.obsm[key] = umap_adata.obsm[key].copy()
        for key in ("icgs3_umap", "icgs3_batch_correction", "icgs3_expression_batch_adjustment"):
            if key in umap_adata.uns:
                adata.uns[key] = dict(umap_adata.uns[key])
    else:
        compute_umap_outputs(adata, config, outdir, marker_features=markers["marker"].tolist())
    write_umap_plots(adata, config, outdir)
    t = step_time("UMAP output", t)

    cluster_cols = [config.cluster_key, "ICGS3_original_NMF_cluster", "ICGS3_SVM_score", "ICGS3_SVM_margin", "ICGS3_cell_state_prediction"]
    clusters = adata.obs[cluster_cols].copy()
    covariate_cols = [c.strip() for c in str(config.umap_covariates or "").split(",") if c.strip()]
    for col in list(dict.fromkeys(["sample", "Library"] + covariate_cols)):
        if col in adata.obs and col not in clusters:
            clusters.insert(0, col, adata.obs[col].astype(str).values)
    clusters.to_csv(os.path.join(outdir, "icgs3_clusters.tsv"), sep="\t")
    clusters.rename_axis("barcode").reset_index().to_csv(
        os.path.join(outdir, "icgs3_cell_barcode_clusters.tsv"),
        sep="\t",
        index=False,
    )

    heatmap_outputs = run_canonical_heatmap(adata, config, outdir)
    adata.uns["icgs3_heatmap_outputs"] = heatmap_outputs
    t = step_time("MarkerFinder heatmap", t)

    adata.uns["icgs3"] = asdict(config)
    if config.write_h5ad:
        adata.write_h5ad(os.path.join(outdir, "icgs3_result.h5ad"), compression="gzip")
        t = step_time("h5ad output", t)
    _log(f"ICGS3 complete in {time.time() - start_time:.1f}s")
    return ICGS3Result(
        adata=adata,
        clusters=clusters,
        markers=markers,
        biomarker_predictions=biomarker_predictions,
        output_dir=outdir,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="ICGS3 sparse clustering, intelligent downsampling, UDON NMF, MarkerFinder, SVM reclassification, UMAP, and GO-Elite summaries.")
    parser.add_argument("--input", nargs="+", required=True, help="Input h5ad, 10x h5, 10x mtx file, or matrix directory.")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--modality", default="rna", choices=["rna", "adt", "grn", "metabolite", "lipid", "psi"])
    parser.add_argument("--layer", default=None)
    parser.add_argument("--input-normalized", action="store_true", help="Treat X/layer as already normalized; skip CP10K log1p.")
    parser.add_argument(
        "--normalized-decimals",
        type=int,
        default=-1,
        help="Optionally round normalized expression values using sparse-safe operations; negative disables. Disabled by default because rounding is not a compute-speed optimization.",
    )
    parser.add_argument("--min-genes", type=int, default=500)
    parser.add_argument("--min-cells", type=int, default=5)
    parser.add_argument("--min-counts", type=int, default=1000)
    parser.add_argument("--mito-percent", type=float, default=30.0)
    parser.add_argument("--target-cells", type=int, default=0, help="Optional memory-guard cap before PageRank; 0 disables.")
    parser.add_argument("--pagerank-cells", type=int, default=5000)
    parser.add_argument("--louvain-downsample-cutoff", type=int, default=15000)
    parser.add_argument(
        "--pre-pagerank-cells",
        type=int,
        default=0,
        help="Optional Louvain pre-reduction target before PageRank; 0 uses ICGS2 default of --pagerank-cells * 4.",
    )
    parser.add_argument("--downsample-var-genes", type=int, default=500, help="ICGS2 hgvfinder variable genes used for Louvain/PageRank downsampling.")
    parser.add_argument(
        "--retention-audit-obs",
        default="HLCA,TGEN-IPF",
        help="Comma-delimited obs columns audited for rare-population retention through sampling and SVM; empty disables.",
    )
    parser.add_argument("--downsample-key", default=None, help="obs column for stratified downsampling.")
    parser.add_argument("--max-cells-per-group", type=int, default=None)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--n-top-features", type=int, default=3000)
    parser.add_argument(
        "--max-nmf-variable-genes",
        type=int,
        default=0,
        help="Optional cap on UDON-selected RNA features used for NMF; 0 keeps all UDON-selected features.",
    )
    parser.add_argument("--n-pcs", type=int, default=0, help="PCs for graph construction; 0 uses automatic elbow selection.")
    parser.add_argument("--max-auto-pcs", type=int, default=50, help="Maximum PCs computed before automatic elbow selection.")
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--leiden-resolution", type=float, default=0.8)
    parser.add_argument(
        "--batch-correction",
        choices=["none", "harmony"],
        default="none",
        help="Optional batch-effect correction applied only to low-dimensional graph/UMAP steps.",
    )
    parser.add_argument(
        "--batch-key",
        default=None,
        help="Comma-delimited obs columns used for batch correction, e.g. Library or ADT_panel,Library. Can be combined with technology/tissue keys.",
    )
    parser.add_argument(
        "--technology-batch-key",
        default=None,
        help="Comma-delimited obs column(s) identifying technology/capture/panel batches, e.g. ADT_panel, Chemistry, or Library.",
    )
    parser.add_argument(
        "--tissue-batch-key",
        default=None,
        help="Comma-delimited obs column(s) identifying tissue/source batches, e.g. TissueGroup for BoneMarrow vs Thymus.",
    )
    parser.add_argument(
        "--batch-correction-use",
        choices=["graph", "umap", "graph-umap", "all"],
        default="graph",
        help="Where to use corrected PCs: graph/PageRank/initial-k, UMAP, or both.",
    )
    parser.add_argument(
        "--batch-adjust-expression",
        action="store_true",
        help=(
            "Use a temporary non-negative batch-centered expression matrix for "
            "PageRank/NMF/MarkerFinder/SVM. Intended for ADT/intersected-panel "
            "batch tests; final heatmap/h5ad expression remains unadjusted."
        ),
    )
    parser.add_argument(
        "--expression-batch-key",
        default=None,
        help=(
            "Comma-delimited obs columns for --batch-adjust-expression. Defaults to --batch-key and also "
            "combines with --technology-batch-key/--tissue-batch-key when supplied."
        ),
    )
    parser.add_argument("--cluster-key", default="ICGS3_cluster")
    parser.add_argument("--marker-top-n", type=int, default=60)
    parser.add_argument(
        "--marker-rho",
        type=float,
        default=None,
        help="MarkerFinder Pearson rho a cluster must reach to survive the marker gate. "
             "Unset resolves by modality: 0.15 for adt, 0.3 otherwise.",
    )
    parser.add_argument(
        "--marker-min-per-cluster",
        type=int,
        default=None,
        help="Markers at or above --marker-rho a cluster must own to survive the marker gate. "
             "Unset resolves by modality: 1 for adt, 2 otherwise.",
    )
    parser.add_argument("--marker-direction", choices=["up", "down", "both"], default="up")
    parser.add_argument("--rank", "--nmf-k", dest="rank", type=int, default=None, help="Override auto-estimated sparse-NMF rank.")
    parser.add_argument(
        "--nmf-assignment-normalization",
        choices=["auto", "raw", "rowsum", "rowmax", "zscore", "quantile95"],
        default="auto",
        help="Normalize SNMF component loadings before assigning cells to their winning component. Auto uses rowsum for ADT/small-feature modalities.",
    )
    parser.add_argument(
        "--nmf-k-multiplier",
        type=float,
        default=2.0,
        help="Deprecated compatibility option; automatic rank now uses the UDON/ICGS2 Tracy-Widom estimator directly.",
    )
    parser.add_argument(
        "--max-auto-nmf-k",
        type=int,
        default=None,
        help="Optional maximum automatically selected NMF rank. By default, auto rank is not artificially capped; explicit --nmf-k can exceed this.",
    )
    parser.add_argument(
        "--nmf-runs",
        type=int,
        default=1,
        help=(
            "Number of SNMF restarts. UDON run_nmf uses seed='nndsvd', a DETERMINISTIC "
            "initialization, so every restart converges to the same factorization and "
            "n_run>1 only multiplies runtime (UDON records run1-vs-run2 ARI=1.000). "
            "Default 1. Raise only if you switch to a stochastic seed."
        ),
    )
    parser.add_argument(
        "--markerfinder-nmf-genes-only",
        action="store_true",
        help=(
            "Restrict MarkerFinder to the NMF guide genes. Default (off) matches ICGS2, "
            "which runs MarkerFinder on the full filtered expression matrix "
            "(protein-coding; RPL/RPS, MT-, HLA, XIST/TSIX and *Y removed). The narrow "
            "pool starves the cluster-fitness gate and rejects clusters ICGS2 would keep."
        ),
    )
    parser.add_argument(
        "--resume-sampled-from",
        type=str,
        default=None,
        help=(
            "Resume at the NMF stage. Path to a completed run's "
            "sNMF/icgs3_pagerank_downsampling.tsv; its selected_final barcodes are reused "
            "as the sampled set and PageRank/Louvain is skipped. QC and gene-filter "
            "settings must match the source run or barcodes will not align."
        ),
    )
    parser.add_argument(
        "--intercorr-threshold", type=float, default=0.4,
        help=("UDON fast_feature_selection gene-gene Pearson threshold for the NMF guide-gene "
              "pool. Default 0.4 reproduces prior behaviour."),
    )
    parser.add_argument(
        "--corr-n-events", type=int, default=5,
        help=("Minimum number of genes a guide gene must correlate with above "
              "--intercorr-threshold (self-correlation counts, as in ICGS2). Default 5. "
              "Lowering this is the cheaper lever for rare populations."),
    )
    parser.add_argument("--rank-rel-threshold", type=float, default=0.1, help="UDON-compatible relative threshold retained for rank reporting/compatibility.")
    parser.add_argument("--min-group-size", type=int, default=3)
    parser.add_argument(
        "--svm-min-decision-score",
        type=float,
        default=0.0,
        help="Minimum winning linear SVM decision score for final cell retention; ICGS2 default is >0.",
    )
    parser.add_argument(
        "--svm-chunk-size",
        type=int,
        default=50000,
        help="Number of cells scored per chunk during all-cell SVM reclassification.",
    )
    svm_group = parser.add_mutually_exclusive_group()
    svm_group.add_argument(
        "--svm-reclassify-all-cells",
        "--svm_reclassify_all_cells",
        dest="svm_reclassify_all_cells",
        action="store_true",
        default=True,
        help="Classify all post-QC cells against surviving NMF/SVM centroids after downsampled NMF training.",
    )
    svm_group.add_argument(
        "--no-svm-reclassify-all-cells",
        "--svm_reclassify_all_cells_FALSE",
        dest="svm_reclassify_all_cells",
        action="store_false",
        help="Testing mode: classify only the sampled cells used for NMF.",
    )
    parser.add_argument("--species", default="Hs", choices=["Hs", "Mm", "human", "mouse"])
    parser.add_argument("--biomarker-file", default=None)
    parser.add_argument(
        "--heatmap-cells-per-cluster",
        type=int,
        default=0,
        help="Cells to display per cluster in the final heatmap; 0 uses all cells (ICGS2-style final output).",
    )
    parser.add_argument("--write-heatmap-expression-tsv", action="store_true", help="Write the large all-cell heatmap expression matrix TSV.")
    parser.add_argument(
        "--heatmap-covariates",
        default=None,
        help="Comma-delimited obs columns to render as compact bottom covariate bars in the MarkerFinder heatmap.",
    )
    parser.add_argument(
        "--heatmap-goelite-terms",
        action="store_true",
        help="Render top GO-Elite/BioMarkers enrichment terms on the left side of the MarkerFinder heatmap.",
    )
    parser.add_argument(
        "--heatmap-goelite-max-terms",
        type=int,
        default=30,
        help="Maximum GO-Elite terms considered per cluster block for MarkerFinder heatmap labels.",
    )
    parser.add_argument(
        "--umap-feature-mode",
        choices=["markerfinder", "variable", "pca"],
        default="markerfinder",
        help="Use final MarkerFinder features, UDON/variable features, or PCA for UMAP.",
    )
    parser.add_argument("--umap-covariates", default=None, help="Comma-delimited obs columns for additional UMAP PDFs.")
    parser.add_argument(
        "--umap-genes",
        default=None,
        help="Comma-delimited genes drawn as continuous expression UMAPs, one PDF per gene. "
             "Cells are drawn lowest expression first and highest last, so the highest sit on top.",
    )
    parser.add_argument("--no-h5ad", action="store_true")
    parser.add_argument("--no-umap", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> ICGS3Result:
    args = build_arg_parser().parse_args(argv)
    config = ICGS3Config(
        input_paths=list(args.input),
        output_dir=args.output_dir,
        modality=args.modality,
        layer=args.layer,
        input_normalized=args.input_normalized,
        normalized_decimals=args.normalized_decimals,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        min_counts=args.min_counts,
        mito_percent=args.mito_percent,
        target_cells=args.target_cells if args.target_cells and args.target_cells > 0 else None,
        pagerank_cells=args.pagerank_cells,
        louvain_downsample_cutoff=args.louvain_downsample_cutoff,
        pre_pagerank_cells=args.pre_pagerank_cells,
        downsample_var_genes=args.downsample_var_genes,
        retention_audit_obs=args.retention_audit_obs,
        downsample_key=args.downsample_key,
        max_cells_per_group=args.max_cells_per_group,
        random_state=args.random_state,
        n_top_features=args.n_top_features,
        max_nmf_variable_genes=args.max_nmf_variable_genes,
        n_pcs=args.n_pcs,
        max_auto_pcs=args.max_auto_pcs,
        n_neighbors=args.n_neighbors,
        leiden_resolution=args.leiden_resolution,
        batch_correction=args.batch_correction,
        batch_key=args.batch_key,
        technology_batch_key=args.technology_batch_key,
        tissue_batch_key=args.tissue_batch_key,
        batch_correction_use=args.batch_correction_use,
        batch_adjust_expression=args.batch_adjust_expression,
        expression_batch_key=args.expression_batch_key,
        cluster_key=args.cluster_key,
        marker_top_n=args.marker_top_n,
        marker_rho=args.marker_rho,
        marker_min_per_cluster=args.marker_min_per_cluster,
        marker_direction=args.marker_direction,
        rank=args.rank,
        nmf_assignment_normalization=args.nmf_assignment_normalization,
        nmf_k_multiplier=args.nmf_k_multiplier,
        max_auto_nmf_k=args.max_auto_nmf_k,
        nmf_runs=args.nmf_runs,
        markerfinder_all_genes=not args.markerfinder_nmf_genes_only,
        resume_sampled_from=args.resume_sampled_from,
        intercorr_threshold=args.intercorr_threshold,
        corr_n_events=args.corr_n_events,
        rank_rel_threshold=args.rank_rel_threshold,
        min_group_size=args.min_group_size,
        svm_min_decision_score=args.svm_min_decision_score,
        svm_reclassify_all_cells=args.svm_reclassify_all_cells,
        svm_chunk_size=args.svm_chunk_size,
        species={"human": "Hs", "mouse": "Mm"}.get(args.species, args.species),
        biomarker_file=args.biomarker_file,
        heatmap_cells_per_cluster=args.heatmap_cells_per_cluster,
        write_heatmap_expression_tsv=args.write_heatmap_expression_tsv,
        heatmap_covariates=args.heatmap_covariates,
        heatmap_goelite_terms=args.heatmap_goelite_terms,
        heatmap_goelite_max_terms=args.heatmap_goelite_max_terms,
        umap_feature_mode=args.umap_feature_mode,
        umap_covariates=args.umap_covariates,
        umap_genes=args.umap_genes,
        write_h5ad=not args.no_h5ad,
        generate_umap=not args.no_umap,
    )
    return run_icgs3(config)


if __name__ == "__main__":
    main()
