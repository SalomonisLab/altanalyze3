"""scALABLE-discover pipeline: scALABLE's QC and ambient RNA correction, then ICGS3.

The stages, in order:

1. `cellHarmony_lite.combine_and_align_h5` with `cellharmony_ref=None` loads the uploads,
   corrects ambient RNA when asked, applies the user's QC thresholds and normalizes. This
   is the same function and the same arguments scALABLE-web uses, stopped before alignment.
2. `ICGS.run_icgs3` clusters the QC-retained counts. QC is not repeated inside ICGS3: its four
   QC thresholds are disabled because step 1 already applied the user's. ICGS3 also writes its
   MarkerFinder marker set, GO-Elite BioMarkers predictions and, with
   `export_marker_networks`, NetPerspective marker networks. The Explore views read those.
3. The combined h5ad holds every QC-retained gene, normalized by the same
   `cellHarmony_lite.normalize_adata` scALABLE-web uses, for the cells ICGS3 assigned to a
   cluster, with ICGS3's cluster columns (without the ICGS3_ prefix) and UMAP attached.
4. Two cell-state layers. `cell_state_predicted`, the default in every view, is ICGS3's GO-Elite
   BioMarkers label of each cluster, one label per cluster; it reads copies of ICGS3's
   MarkerFinder tables in which only the cluster name changes. `cluster` (C1..Cn) is the
   alternative layer and reads ICGS3's own tables.
5. fastComm scores receptor-ligand communication between cell states, once per layer, with
   scALABLE-web's parameters.

No differential expression and no modality imputation run here.
"""
from __future__ import annotations

import gc
import json
import os
import platform
import subprocess
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd

from altanalyze3.components.cellHarmony import cellHarmony_lite
from altanalyze3.components.cellHarmony.flask import pipeline as web_pipeline
from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
from altanalyze3.components.clustering import ICGS
from altanalyze3.components.fastComm.api import FastCommParams, run_fastcomm
from altanalyze3.components.fastComm.benchmark import FastCommBenchmarkParams, run_benchmark as run_fastcomm_benchmark
from altanalyze3.components.fastComm.reporting import ExemplarReportParams, write_exemplar_report
from altanalyze3.components.visualization import marker_heatmap_h5ad

DISCOVER_REFERENCE_ID = "icgs3"
SPECIES_TO_ICGS = {"human": "Hs", "mouse": "Mm"}
# Cell annotations carry no ICGS3_ prefix in the combined h5ad (Nathan, 2026-10-05).
CLUSTER_KEY = "cluster"
PREDICTION_KEY = "cell_state_predicted"
# Every view opens on the predicted names (Nathan, 2026-10-05); the C1..Cn clusters are the
# alternative layer.
DEFAULT_LAYER = PREDICTION_KEY
ICGS3_COLUMN_NAMES = {
    "ICGS3_cluster": CLUSTER_KEY,
    "ICGS3_cell_state_prediction": PREDICTION_KEY,
    "ICGS3_original_NMF_cluster": "original_NMF_cluster",
    "ICGS3_SVM_score": "SVM_score",
    "ICGS3_SVM_margin": "SVM_margin",
}
LAYER_LABELS = {CLUSTER_KEY: "Clusters", PREDICTION_KEY: "Predicted cell states"}
# ICGS3's two downsampling steps, set by scALABLE-discover (Nathan, 2026-10-05). Step 1, Louvain
# community sampling, runs above LOUVAIN_CUTOFF cells and keeps PRE_PAGERANK_CELLS; step 2,
# PageRank, keeps PAGERANK_CELLS, which train NMF. ICGS3's own defaults are 30,000, 4x30,000
# and 30,000; the SVM still classifies every QC-retained cell.
LOUVAIN_CUTOFF = 10000
PRE_PAGERANK_CELLS = 10000
PAGERANK_CELLS = 5000
# What a run exports (Nathan, 2026-10-05). "minimal", the default: the final combined h5ad and
# scALABLE's MarkerFinder output set, with ICGS3 --minimal-outputs and no ICGS3 h5ad; the ICGS3
# input and ambient raw matrices are staged in a temporary folder removed after use. "full":
# every ICGS3 file, the staged inputs, cluster assignments and the ZIP archives.
EXPORT_MODES = ("minimal", "full")
EXPORT_MODE_ENV = "SCALABLE_DISCOVER_EXPORTS"


def export_mode() -> str:
    mode = str(os.environ.get(EXPORT_MODE_ENV, "minimal") or "minimal").strip().lower()
    if mode not in EXPORT_MODES:
        raise ValueError(f"{EXPORT_MODE_ENV}={mode!r}; choose one of {', '.join(EXPORT_MODES)}.")
    return mode


def discover_registry() -> Dict:
    """The species menu. One pseudo-reference per species: ICGS3 needs no atlas."""
    return {"species": [
        {"id": species, "label": species.capitalize(),
         "references": [{"id": DISCOVER_REFERENCE_ID, "label": "ICGS3 unsupervised clustering"}]}
        for species in SPECIES_TO_ICGS
    ]}


def _git_commit(path: Path) -> str:
    try:
        out = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=10)
        dirty = subprocess.run(["git", "-C", str(path), "status", "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True, timeout=10)
        commit = out.stdout.strip() or "unknown"
        return commit + ("+uncommitted-changes" if dirty.stdout.strip() else "")
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def _counts_matrix(adata: ad.AnnData):
    """The QC-retained counts ICGS3 clusters, and whether they are counts at all.

    cellHarmony-lite stores counts in layers['counts'] only when its scale check says X
    held counts (cellHarmony_lite.py, `_scale["x_verdict"] == "counts"`). An h5ad that
    arrived already log-normalized has no counts layer; ICGS3 then reads X as normalized.
    """
    if "counts" in adata.layers:
        return adata.layers["counts"], True
    return adata.X, False


def _network_records(networks) -> List[Dict[str, str]]:
    """ICGS3 stores its NetPerspective networks as a table in uns; the viewer reads dicts."""
    if networks is None:
        return []
    if isinstance(networks, pd.DataFrame):
        return networks.astype(str).to_dict("records")
    return [dict(entry) for entry in networks]


def _write_matrix_h5ad(path: Path, matrix, obs: pd.DataFrame, var: pd.DataFrame, compression) -> None:
    out = ad.AnnData(X=matrix, obs=obs.copy(), var=var.copy())
    out.write_h5ad(path, compression=compression)


def prediction_names(clusters: pd.DataFrame) -> Dict[str, str]:
    """Cluster -> ICGS3 GO-Elite BioMarkers label. Refuses anything but one label per cluster."""
    pairs = clusters[["ICGS3_cluster", "ICGS3_cell_state_prediction"]].astype(str).drop_duplicates()
    if pairs["ICGS3_cluster"].duplicated().any():
        raise RuntimeError("ICGS3 gave one cluster more than one cell-state prediction: "
                           f"{pairs[pairs['ICGS3_cluster'].duplicated(keep=False)].values.tolist()}")
    if pairs["ICGS3_cell_state_prediction"].duplicated().any():
        raise RuntimeError("ICGS3 gave two clusters the same cell-state prediction: "
                           f"{pairs[pairs['ICGS3_cell_state_prediction'].duplicated(keep=False)].values.tolist()}")
    return dict(zip(pairs["ICGS3_cluster"], pairs["ICGS3_cell_state_prediction"]))


def _relabel_tsv_column(src: Path, dest: Path, column: str, names: Dict[str, str]) -> int:
    """Copy a TSV, renaming one column's values. Every other byte stays as written."""
    lines = Path(src).read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    if column not in header:
        raise ValueError(f"{src} has no '{column}' column.")
    index = header.index(column)
    changed = 0
    out = [lines[0]]
    for line in lines[1:]:
        fields = line.split("\t")
        if index < len(fields) and fields[index] in names:
            fields[index] = names[fields[index]]
            changed += 1
        out.append("\t".join(fields))
    Path(dest).write_text("\n".join(out) + "\n", encoding="utf-8")
    return changed


def _relabel_tsv_header(src: Path, dest: Path, names: Dict[str, str]) -> None:
    """Copy a TSV whose columns are cell states, renaming the header only."""
    lines = Path(src).read_text(encoding="utf-8").splitlines()
    header = [names.get(field, field) for field in lines[0].split("\t")]
    Path(dest).write_text("\n".join(["\t".join(header)] + lines[1:]) + "\n", encoding="utf-8")


def relabel_marker_outputs(marker_analysis: Dict, names: Dict[str, str], out_dir: Path) -> Dict:
    """Copies of ICGS3's MarkerFinder set keyed by predicted name instead of cluster id.

    The markers, folds, cells and order are ICGS3's own; only the cluster name changes, by the
    one-to-one map `names`. ICGS3's files stay untouched.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    relabeled = dict(marker_analysis, cluster_key=PREDICTION_KEY, archive="")
    for key in ("markers_tsv", "redundant_markers_tsv"):
        src = str(marker_analysis.get(key) or "")
        if src and Path(src).is_file():
            dest = out_dir / Path(src).name
            _relabel_tsv_column(Path(src), dest, "cluster", names)
            relabeled[key] = str(dest)
    centroids = str(marker_analysis.get("centroids_tsv") or "")
    if centroids and Path(centroids).is_file():
        dest = out_dir / Path(centroids).name
        _relabel_tsv_header(Path(centroids), dest, names)
        relabeled["centroids_tsv"] = str(dest)
    cache = str(marker_analysis.get("heatmap_cache") or "")
    if cache and Path(cache).is_file():
        heatmap_df, row_clusters, column_clusters, _, _ = marker_heatmap_h5ad._read_heatmap_cache(cache)
        dest = out_dir / Path(cache).name
        marker_heatmap_h5ad._write_heatmap_cache(
            str(dest), heatmap_df, [names.get(c, c) for c in row_clusters],
            [names.get(c, c) for c in column_clusters], list(heatmap_df.columns))
        relabeled["heatmap_cache"] = str(dest)
    relabeled["networks"] = [dict(entry, population=names.get(str(entry.get("population")), entry.get("population")))
                             for entry in marker_analysis.get("networks") or []]
    return relabeled


def _run_fastcomm_layer(combined: ad.AnnData, combined_h5ad_path: Path, state_key: str, fastcomm_dir: Path,
                        species: str) -> Tuple[Dict[str, object], Dict[str, Path]]:
    """fastComm between the cell states of one layer, with scALABLE-web's parameters."""
    fastcomm_dir.mkdir(parents=True, exist_ok=True)
    scores_path = fastcomm_dir / "fastcomm_scores.tsv"
    pairs_path = fastcomm_dir / "state_pair_summary.tsv"
    state_expression_path = fastcomm_dir / "state_expression.tsv"
    result = run_fastcomm(FastCommParams(
        adata=combined, output=scores_path, response_matrix=None, lr_sources=("CellChatDB",),
        state_pair_output=pairs_path, state_expression_output=state_expression_path,
        state_key=state_key, species=species, min_cells=5, min_lr_expression_score=0.2,
        max_lr_candidates_per_state_pair=5, include_self_edges=False,
    ))
    significant_path = fastcomm_dir / "significant_interactions.tsv"
    significant_md_path = fastcomm_dir / "significant_interactions.md"
    write_exemplar_report(ExemplarReportParams(
        scores=scores_path, output_tsv=significant_path, output_md=significant_md_path,
        top_n=100000, min_score=0.25, one_per_state_pair=False,
        title="Cell communication significant interactions",
    ))
    sample_key = next((c for c in ("Library", "group", "sample", "Donor") if c in combined.obs.columns), None)
    split_summary: Dict[str, object] = {}
    if sample_key:
        split_summary = run_fastcomm_benchmark(FastCommBenchmarkParams(
            h5ad=combined_h5ad_path, output_dir=fastcomm_dir / "per_sample", state_key=state_key,
            split_key=sample_key, response_matrix=None, lr_sources=["CellChatDB"], species=species,
            min_cells=5, min_lr_expression_score=0.2, max_lr_candidates_per_state_pair=5,
            include_self_edges=False,
        ))
    analysis = {
        "enabled": True, "status": "completed", "state_key": state_key, "sample_key": sample_key,
        "populations": [str(v) for v in result.state_sizes.index.tolist()],
        "scores_tsv": str(scores_path), "state_pair_tsv": str(pairs_path),
        "state_expression_tsv": str(state_expression_path), "significant_tsv": str(significant_path),
        "significant_md": str(significant_md_path), "significance_threshold": 0.25,
        "archive": None, "per_sample": split_summary, "summary": result.summary,
    }
    paths = {"fastcomm_scores": scores_path, "fastcomm_state_pairs": pairs_path,
             "fastcomm_state_expression": state_expression_path,
             "fastcomm_significant_interactions": significant_path,
             "fastcomm_significant_report": significant_md_path}
    return analysis, paths


def run_discover_pipeline(job_id: str, store: JobStore, *, h5ad_compression: Optional[str] = "lzf") -> Dict[str, Path]:
    meta = store.get_job(job_id)
    species = str(meta.get("species") or "").strip().lower()
    if species not in SPECIES_TO_ICGS:
        raise ValueError(f"scALABLE-discover supports human or mouse, not '{species}'.")
    if str(meta.get("reference") or "") != DISCOVER_REFERENCE_ID:
        raise ValueError(f"scALABLE-discover jobs carry reference '{DISCOVER_REFERENCE_ID}', "
                         f"not '{meta.get('reference')}'.")
    compression = web_pipeline._normalize_h5ad_compression(h5ad_compression)
    uploads_dir = store.uploads_dir(job_id)
    outputs_dir = store.outputs_dir(job_id)
    outputs_dir.mkdir(parents=True, exist_ok=True)
    h5_files, h5ad_file = web_pipeline._split_uploads(meta["files"], uploads_dir)
    qc = meta.get("qc", {})
    min_genes = int(qc.get("min_genes", 500))
    min_cells = int(qc.get("min_cells", 0))
    min_counts = int(qc.get("min_counts", 1000))
    mit_percent = int(qc.get("mit_percent", 15))
    ambient_correction = str(qc.get("ambient_correction", "no") or "no").strip().lower()
    ambient_rho = "auto" if ambient_correction == "yes" else None

    # ---- 1. scALABLE QC: the web pipeline's own call, stopped before alignment.
    store.update_job(job_id, progress=20, message="QC: loading samples, ambient RNA correction and QC filters.")
    store.append_log(job_id, "Running cellHarmony_lite QC (no reference; alignment skipped).")
    store.append_log(
        job_id,
        "[params] qc "
        f"input_mode={'single_h5ad' if h5ad_file else 'multi_h5'} n_h5={len(h5_files)} "
        f"h5ad={'yes' if h5ad_file else 'no'} min_genes={min_genes} min_cells={min_cells} "
        f"min_counts={min_counts} mit_percent={mit_percent} ambient_correction={ambient_correction}",
    )
    _, qc_adata = cellHarmony_lite.combine_and_align_h5(
        h5_files=h5_files,
        h5ad_file=h5ad_file,
        cellharmony_ref=None,
        output_dir=str(outputs_dir),
        export_cptt=False,
        export_h5ad=False,
        min_genes=min_genes,
        min_cells=min_cells,
        min_counts=min_counts,
        mit_percent=mit_percent,
        generate_umap=False,
        save_adata=False,
        unsupervised_cluster=False,
        gene_translation_file=None,
        metacell_align=False,
        ambient_correct_cutoff=ambient_rho,
        ambient_memory_efficient=True,
        concat_on_disk=True,
        concat_batch_size=1,
        stream_10x_inputs=True,
        return_adata=True,
    )
    n_qc_cells, n_qc_genes = int(qc_adata.n_obs), int(qc_adata.n_vars)
    if n_qc_cells < 2:
        raise ValueError(f"QC retained {n_qc_cells} cells; ICGS3 needs more. Lower the QC thresholds.")

    # ICGS3 reads counts only. The ambient raw matrix is stored apart and re-attached to the
    # combined h5ad, as scALABLE-web keeps it, without ICGS3 ever loading it.
    exports = export_mode()
    minimal = exports == "minimal"
    staging = tempfile.TemporaryDirectory(prefix="icgs3_staging_", dir=outputs_dir) if minimal else None
    icgs3_input_dir = Path(staging.name) if minimal else outputs_dir / "ICGS3_input"
    icgs3_input_dir.mkdir(parents=True, exist_ok=True)
    icgs3_input_path = icgs3_input_dir / "icgs3_input_counts.h5ad"
    counts, holds_counts = _counts_matrix(qc_adata)
    _write_matrix_h5ad(icgs3_input_path, counts, qc_adata.obs, qc_adata.var, compression)
    soupx_raw_path = None
    if "soupx_raw" in qc_adata.layers:
        soupx_raw_path = (icgs3_input_dir if minimal else outputs_dir / "ambient") / "soupx_raw_counts.h5ad"
        soupx_raw_path.parent.mkdir(parents=True, exist_ok=True)
        _write_matrix_h5ad(soupx_raw_path, qc_adata.layers["soupx_raw"], qc_adata.obs[[]], qc_adata.var[[]], compression)
    store.append_log(job_id, f"QC retained {n_qc_cells} cells x {n_qc_genes} genes; ICGS3 input "
                             f"holds {'counts' if holds_counts else 'normalized X (no counts layer)'}.")
    qc_adata = counts = None
    gc.collect()

    # ---- 2. ICGS3 through its documented top-level function.
    icgs3_dir = outputs_dir / "ICGS3"
    network_jobs = max(1, min(4, os.cpu_count() or 1))
    config = ICGS.ICGS3Config(
        input_paths=[str(icgs3_input_path)],
        output_dir=str(icgs3_dir),
        modality="rna",
        normalization="auto" if holds_counts else "none",
        # QC ran once, above, with the user's thresholds. Zero/None disables each filter in
        # ICGS.apply_qc, so ICGS3 removes no further cell or gene.
        min_genes=0,
        min_cells=0,
        min_counts=0,
        mito_percent=None,
        louvain_downsample_cutoff=LOUVAIN_CUTOFF,
        pre_pagerank_cells=PRE_PAGERANK_CELLS,
        pagerank_cells=PAGERANK_CELLS,
        species=SPECIES_TO_ICGS[species],
        export_marker_networks=True,
        marker_network_top_n=1000,
        marker_network_jobs=network_jobs,
        minimal_outputs=minimal,
        write_h5ad=not minimal,     # the combined h5ad below is the final h5ad
    )
    store.update_job(job_id, progress=35, message="ICGS3: starting.")
    store.append_log(job_id, "Running ICGS3 unsupervised clustering.")
    store.append_log(job_id, f"[params] icgs3 module=altanalyze3.components.clustering.ICGS.run_icgs3 "
                             f"species={config.species} normalization={config.normalization} "
                             "qc=disabled(applied upstream) export_marker_networks=True "
                             f"louvain_downsample_cutoff={LOUVAIN_CUTOFF} pre_pagerank_cells={PRE_PAGERANK_CELLS} "
                             f"pagerank_cells={PAGERANK_CELLS} "
                             f"marker_network_top_n=1000 marker_network_jobs={network_jobs} "
                             f"exports={exports} minimal_outputs={minimal} write_h5ad={not minimal} "
                             "all other parameters=ICGS3Config defaults")
    result = ICGS.run_icgs3(config)
    clusters = result.clusters.copy()
    clusters.index = clusters.index.astype(str)
    if "X_umap" not in result.adata.obsm:
        error_file = icgs3_dir / "UMAPs" / "icgs3_umap_error.txt"
        detail = error_file.read_text().strip() if error_file.exists() else "no error file"
        raise RuntimeError(f"ICGS3 produced no UMAP, which every Explore view needs: {detail}")
    umap = pd.DataFrame(np.asarray(result.adata.obsm["X_umap"], dtype=float)[:, :2],
                        index=result.adata.obs_names.astype(str), columns=["umap_0", "umap_1"])
    cluster_order = [str(v) for v in (result.adata.uns.get("lineage_order") or [])]
    heatmap = dict(result.adata.uns.get("icgs3_heatmap_outputs") or {})
    # The gene universe ICGS3's BioMarkers enrichment tested against (ICGS.biomarker_enrichment).
    biomarker_background = len({str(g).upper() for g in result.adata.var_names})
    icgs3_logs = sorted((icgs3_dir / "logs").glob("icgs3_*.log"))
    cli_line = ICGS.cli_equivalent(config)
    result = None
    gc.collect()
    n_clustered = int(len(clusters))
    store.append_log(job_id, f"ICGS3 assigned {n_clustered} of {n_qc_cells} QC-retained cells "
                             f"({100.0 * n_clustered / n_qc_cells:.1f}%) to {len(cluster_order)} clusters; "
                             f"{n_qc_cells - n_clustered} cells received no ICGS3 cluster and are not shown.")
    missing = [key for key in ("markers_tsv", "heatmap_cache") if not heatmap.get(key)]
    if missing:
        raise RuntimeError(f"ICGS3 MarkerFinder heatmap outputs missing: {', '.join(missing)} ({heatmap})")
    names = prediction_names(clusters)
    name_order = [names.get(c, c) for c in cluster_order]

    # ---- 3. Combined h5ad: every QC gene, scALABLE normalization, ICGS3 cells and labels.
    store.update_job(job_id, progress=81, message="Writing the combined h5ad.")
    combined = ad.read_h5ad(icgs3_input_path)
    combined.obs_names = combined.obs_names.astype(str)
    missing_cells = clusters.index.difference(combined.obs_names)
    if len(missing_cells):
        raise RuntimeError(f"{len(missing_cells)} ICGS3 barcodes are absent from the QC matrix, "
                           f"e.g. {list(missing_cells[:3])}")
    combined = combined[clusters.index].copy()
    if holds_counts:
        combined.layers["counts"] = combined.X.copy()
        cellHarmony_lite.normalize_adata(combined)
    if soupx_raw_path is not None:
        raw = ad.read_h5ad(soupx_raw_path)
        raw.obs_names = raw.obs_names.astype(str)
        combined.layers["soupx_raw"] = raw[clusters.index].X
        raw = None
    if staging is not None:
        staging.cleanup()           # staged ICGS3 input and ambient raw matrix; not exported
    for source, column in ICGS3_COLUMN_NAMES.items():
        if source in clusters.columns:
            combined.obs[column] = clusters[source].values
    combined.obs[CLUSTER_KEY] = pd.Categorical(clusters["ICGS3_cluster"].astype(str).values,
                                               categories=cluster_order, ordered=True)
    combined.obs[PREDICTION_KEY] = pd.Categorical(clusters["ICGS3_cluster"].astype(str).map(names).values,
                                                  categories=name_order, ordered=True)
    combined.obsm["X_umap"] = umap.loc[combined.obs_names, ["umap_0", "umap_1"]].to_numpy()
    combined.uns["lineage_order"] = name_order       # the default layer's state order
    combined.uns["scalable_discover"] = {"default_layer": DEFAULT_LAYER, "cluster_key": CLUSTER_KEY,
                                         "prediction_key": PREDICTION_KEY,
                                         "icgs3_output_dir": str(icgs3_dir),
                                         "qc_cells": n_qc_cells, "clustered_cells": n_clustered}
    combined_h5ad_path = outputs_dir / "combined_with_umap_and_markers.h5ad"
    web_pipeline.approx_mod.ensure_h5ad_compat_for_write(combined)
    combined.write_h5ad(combined_h5ad_path, compression=compression)

    coordinates_path = outputs_dir / "icgs3_umap_coordinates.tsv"
    umap.loc[combined.obs_names].rename_axis("CellBarcode").to_csv(coordinates_path, sep="\t")
    assignments_path = outputs_dir / "icgs3_cluster_assignments.txt"
    if not minimal:
        assignments = combined.obs[[c for c in ("Library", *ICGS3_COLUMN_NAMES.values()) if c in combined.obs]].copy()
        assignments["UMAP1"] = combined.obsm["X_umap"][:, 0]
        assignments["UMAP2"] = combined.obsm["X_umap"][:, 1]
        assignments.rename_axis("CellBarcode").to_csv(assignments_path, sep="\t")

    marker_dir = Path(heatmap["markers_tsv"]).parent
    marker_archive_path = outputs_dir / "icgs3_marker_genes.zip"
    web_pipeline._write_selected_zip(marker_dir, marker_archive_path, suffixes=(".pdf", ".tsv", ".png"))
    icgs3_archive_path = outputs_dir / "icgs3_results.zip"
    if not minimal:
        with zipfile.ZipFile(icgs3_archive_path, "w", zipfile.ZIP_DEFLATED, compresslevel=1) as handle:
            for file_path in sorted(icgs3_dir.rglob("*")):
                if file_path.is_file() and file_path.suffix.lower() != ".h5ad":
                    handle.write(file_path, arcname=str(file_path.relative_to(icgs3_dir)))
    marker_analysis = {
        "enabled": True,
        "source": "ICGS3 MarkerFinder (ICGS.run_canonical_heatmap)",
        "cluster_key": CLUSTER_KEY,
        "markers_tsv": str(heatmap["markers_tsv"]),
        "redundant_markers_tsv": str(heatmap.get("redundant_markers_tsv") or ""),
        "centroids_tsv": str(heatmap.get("centroids_tsv") or ""),
        "heatmap_cache": str(heatmap["heatmap_cache"]),
        "heatmap_pdf": str(heatmap.get("pdf") or ""),
        "networks": _network_records(heatmap.get("networks")),
        "archive": str(marker_archive_path),
        "top_n": int(config.marker_top_n),
    }
    markers = pd.read_csv(heatmap["markers_tsv"], sep="\t")
    default_gene = str(markers["Gene"].iloc[0]) if "Gene" in markers.columns and len(markers) else None
    layer_dir = outputs_dir / "layers" / PREDICTION_KEY
    prediction_marker_analysis = relabel_marker_outputs(marker_analysis, names, layer_dir / "MarkerFinder")
    goelite_dir = icgs3_dir / "GO-Elite"
    goelite = {
        "enrichment_tsv": str(goelite_dir / "icgs3_biomarker_enrichment.tsv"),
        "predictions_tsv": str(goelite_dir / "icgs3_cell_state_predictions.tsv"),
        "background_size": biomarker_background,
        "z_score": "altanalyze3.components.goelite.structures.compute_z_score",
        "available": (goelite_dir / "icgs3_biomarker_enrichment.tsv").is_file(),
    }

    # ---- 4. fastComm between cell states, once per layer, scALABLE-web's parameters.
    # umap_coordinates is the Explore views' UMAP source; the page never lists it as a download.
    artifacts: Dict[str, Path] = {
        "combined_h5ad": combined_h5ad_path,
        "umap_coordinates": coordinates_path,
        "marker_genes_zip": marker_archive_path,
    }
    if not minimal:
        artifacts.update({"assignments": assignments_path, "icgs3_results_zip": icgs3_archive_path})
    fastcomm_by_layer: Dict[str, Dict[str, object]] = {}
    cluster_layer_dir = outputs_dir / "layers" / CLUSTER_KEY
    for layer_key, fastcomm_dir, progress in ((PREDICTION_KEY, outputs_dir / "fastComm", 88),
                                              (CLUSTER_KEY, cluster_layer_dir / "fastComm", 92)):
        try:
            store.update_job(job_id, progress=progress,
                             message=f"fastComm: cell communication between {LAYER_LABELS[layer_key].lower()}.")
            store.append_log(job_id, f"Running fastComm receptor-ligand communication analysis ({layer_key}).")
            analysis, paths = _run_fastcomm_layer(combined, combined_h5ad_path, layer_key, fastcomm_dir, species)
            if layer_key == DEFAULT_LAYER:
                artifacts.update(paths)
                if not minimal:
                    archive = outputs_dir / "cell_communication_fastcomm.zip"
                    web_pipeline._write_selected_zip(fastcomm_dir, archive, suffixes=(".tsv", ".json", ".md"))
                    artifacts["fastcomm_archive"] = archive
                    analysis["archive"] = str(archive)
            fastcomm_by_layer[layer_key] = analysis
            store.append_log(job_id, f"fastComm analysis complete ({layer_key}): "
                                     f"states={analysis['summary'].get('n_states')} "
                                     f"edges={analysis['summary'].get('n_scored_edges')} "
                                     f"genes={analysis['summary'].get('n_loaded_genes')}")
        except Exception as exc:
            # Recorded, not hidden: the job's metadata and log say fastComm failed and why,
            # and the Cell communication view is withheld, exactly as scALABLE-web does.
            fastcomm_by_layer[layer_key] = {
                "enabled": False, "status": "failed", "state_key": layer_key,
                "message": f"{type(exc).__name__}: {exc}" if str(exc) else type(exc).__name__}
            store.append_log(job_id, f"fastComm analysis skipped ({layer_key}): "
                                     f"{fastcomm_by_layer[layer_key]['message']}")
    fastcomm_analysis = fastcomm_by_layer[DEFAULT_LAYER]

    parameters_path = outputs_dir / "scalable_discover_parameters.json"
    repo_root = Path(__file__).resolve().parents[4]
    parameters = {
        "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "job_id": job_id,
        "species": species,
        "exports": exports,
        "inputs": [str(uploads_dir / r["filename"]) for r in meta["files"]],
        "qc": {"module": "altanalyze3.components.cellHarmony.cellHarmony_lite.combine_and_align_h5",
               "cellharmony_ref": None, "min_genes": min_genes, "min_cells": min_cells,
               "min_counts": min_counts, "mit_percent": mit_percent, "ambient_correct_cutoff": ambient_rho,
               "cells_retained": n_qc_cells, "genes_retained": n_qc_genes},
        "icgs3": {"module": "altanalyze3.components.clustering.ICGS.run_icgs3",
                  "config_json": str(icgs3_dir / "icgs3_config.json"),
                  "log": str(icgs3_logs[-1]) if icgs3_logs else None,
                  "cli_equivalent": cli_line,
                  "cli_equivalent_note": "mito_percent=None has no CLI spelling; the CLI default 30 "
                                         "would apply, so rerun through run_icgs3 or icgs3_config.json. "
                                         "Under minimal exports the input path was a temporary file.",
                  "non_default_parameters": {"min_genes": 0, "min_cells": 0, "min_counts": 0,
                                             "mito_percent": None, "normalization": config.normalization,
                                             "louvain_downsample_cutoff": LOUVAIN_CUTOFF,
                                             "pre_pagerank_cells": PRE_PAGERANK_CELLS,
                                             "pagerank_cells": PAGERANK_CELLS,
                                             "species": config.species, "export_marker_networks": True,
                                             "marker_network_top_n": 1000, "marker_network_jobs": network_jobs,
                                             "minimal_outputs": minimal, "write_h5ad": not minimal},
                  "cells_clustered": n_clustered, "clusters": cluster_order},
        "cell_state_layers": {CLUSTER_KEY: cluster_order, PREDICTION_KEY: name_order,
                              "predicted_layer_marker_tables": "copies of ICGS3's MarkerFinder tables with "
                                                               "cluster ids replaced by predicted names"},
        "goelite_biomarkers": goelite,
        "fastcomm": {"module": "altanalyze3.components.fastComm.api.run_fastcomm",
                     "state_keys": [CLUSTER_KEY, PREDICTION_KEY],
                     "lr_sources": ["CellChatDB"], "min_cells": 5, "min_lr_expression_score": 0.2,
                     "max_lr_candidates_per_state_pair": 5,
                     "status": {k: v.get("status", "disabled") for k, v in fastcomm_by_layer.items()}},
        "versions": {"altanalyze3_commit": _git_commit(repo_root), "python": platform.python_version(),
                     "anndata": ad.__version__, "numpy": np.__version__, "pandas": pd.__version__},
    }
    parameters_path.write_text(json.dumps(parameters, indent=2, default=str), encoding="utf-8")
    if not minimal:
        artifacts["parameters_json"] = parameters_path
    for key, path in artifacts.items():
        store.add_artifact(job_id, key, path)

    upload_profile = web_pipeline._upload_profile(meta)
    modalities_payload = web_pipeline._modalities_payload([])
    sample_key = next((c for c in ("Library", "group", "sample") if c in combined.obs.columns), None)
    combined = None
    gc.collect()
    icgs3_analysis = {
        "status": "completed", "output_dir": str(icgs3_dir), "cluster_key": CLUSTER_KEY,
        "clusters": cluster_order, "n_clusters": len(cluster_order),
        "qc_cells": n_qc_cells, "qc_genes": n_qc_genes, "clustered_cells": n_clustered,
        "unassigned_cells": n_qc_cells - n_clustered,
        "log": str(icgs3_logs[-1]) if icgs3_logs else None, "cli_equivalent": cli_line,
        "parameters_json": str(parameters_path), "goelite": goelite,
    }
    # The default layer's values sit at the top level of the job, where scALABLE's readers
    # look. The alternative layer's values sit in cell_state_layers; the app swaps them in
    # for a viewer who picks that layer.
    cell_state_layers = {
        "default": DEFAULT_LAYER,
        "names": names,
        "layers": [
            {"key": PREDICTION_KEY, "label": LAYER_LABELS[PREDICTION_KEY], "states": name_order,
             "marker_analysis": prediction_marker_analysis,
             "fastcomm_analysis": fastcomm_by_layer[PREDICTION_KEY]},
            {"key": CLUSTER_KEY, "label": LAYER_LABELS[CLUSTER_KEY], "states": cluster_order,
             "marker_analysis": marker_analysis,
             "fastcomm_analysis": fastcomm_by_layer[CLUSTER_KEY]},
        ],
    }
    store.update_job(
        job_id,
        cluster_key=DEFAULT_LAYER,
        reference_cluster_key=None,
        default_gene=default_gene,
        marker_analysis=prediction_marker_analysis,
        marker_analysis_by_modality={"rna": dict(prediction_marker_analysis)},
        fastcomm_analysis=fastcomm_analysis,
        fastcnv_analysis={"enabled": False},
        icgs3_analysis=icgs3_analysis,
        cell_state_layers=cell_state_layers,
        modalities=modalities_payload,
        modality_artifacts={"rna": {"h5ad": str(combined_h5ad_path)}},
        selected_impute_modality=None,
        selected_impute_modalities=[],
        differential_options={
            "enabled": False, "upload_profile": upload_profile,
            "sample_names": web_pipeline._job_sample_names(meta),
            "population_columns": [{"value": DEFAULT_LAYER, "label": DEFAULT_LAYER, "n_categories": len(cluster_order)}],
            "default_population_col": DEFAULT_LAYER, "sample_fields": [], "sample_values": {},
            "default_sample_field": sample_key, "comparison_types": ["cells"],
        },
        differential_history={},
        differential=web_pipeline._default_differential_state(False, DEFAULT_LAYER, sample_key),
        message="ICGS3 clustering completed.",
    )
    bundle_record = web_pipeline._build_job_bundle(store, job_id, combined_h5ad_path, DEFAULT_LAYER,
                                                   {"rna": {"h5ad": str(combined_h5ad_path)}}, modalities_payload)
    store.update_job(job_id, bundle=bundle_record, message="ICGS3 clustering completed.")
    store.append_log(job_id, "scALABLE-discover pipeline finished.")
    return artifacts
