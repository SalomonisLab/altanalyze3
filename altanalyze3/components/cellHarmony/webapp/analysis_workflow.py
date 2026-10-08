"""Compose the existing supervised and ICGS3 workflows without changing their methods.

Both branches start from the same QC-normalized matrix. Their accepted barcode
sets define the final union; missing branch assignments are explicit categories.
Expression is taken from shared QC, never averaged across the two pipelines.
"""
from __future__ import annotations

import gc
import contextlib
import json
import shutil
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem

from ..flask import pipeline as P
from ..flask.job_manager import JobStore
from ..mapped_h5ad import Workspace, inspect_h5ad, merge_h5ads, needs_disk_backed_import

MODES = ("supervised", "unsupervised", "both")
SUPERVISED = "supervised_state"
UNSUPERVISED = "unsupervised_cluster"
UNALIGNED = "Unaligned"
UNCLUSTERED = "Not clustered"


class BranchStore(JobStore):
    """Persist branch artifacts separately while reporting progress to the parent."""
    def __init__(self, parent, job_id, branch, prepared):
        self.parent, self.parent_id, self.branch = parent, job_id, branch
        super().__init__(parent.outputs_dir(job_id) / "branches" / branch)
        folder = self._job_dir(job_id)
        folder.mkdir(parents=True, exist_ok=True)
        for sub in ("outputs", "uploads", "logs"):
            (folder / sub).mkdir(exist_ok=True)
        meta = parent.get_job(job_id)
        meta = dict(meta, artifacts={}, bundle={}, cell_state_layers={},
                    files=[{"filename": str(prepared.resolve()), "sample_name": "shared_QC"}],
                    qc=dict(meta.get("qc") or {}, ambient_correction="no"))
        if branch == "unsupervised":
            meta["reference"] = "icgs3"
            meta["qc"].update(impute_modalities=[], impute_modality="none")
        self._write_metadata(job_id, meta)
        self._stage_log = None
        if branch == "unsupervised":
            from ..scalable_discover.tasks import _StageLogStream
            self._stage_log = _StageLogStream(self, job_id)

    def append_log(self, job_id, message):
        super().append_log(job_id, message)
        self.parent.append_log(self.parent_id, f"[{self.branch}] {message}")
        # Reuse discover's native stage reporting in the unified workflow.
        # The analysis method is unchanged; branch updates map into parent progress.
        if self._stage_log is not None:
            self._stage_log._stage(message)

    def update_job(self, job_id, **changes):
        record = super().update_job(job_id, **changes)
        if "message" in changes:
            progress = changes.get("progress", 50)
            mode = (self.parent.get_job(self.parent_id).get("qc") or {}).get("analysis_mode")
            base, width = ((20, 35) if self.branch == "supervised" else (55, 35)) if mode == "both" else (20, 70)
            self.parent.update_job(self.parent_id, message=f"{self.branch.capitalize()}: {changes['message']}",
                                   progress=base + int(width * min(100, progress) / 100))
        return record


def _shared_qc(job_id, store, compression):
    meta = store.get_job(job_id)
    work = store.outputs_dir(job_id) / ".h5ad_work" / "shared"
    work.mkdir(parents=True, exist_ok=True)
    h5, h5ad = P._split_uploads(meta["files"], store.uploads_dir(job_id))
    bounded = False
    if h5ad:
        paths = [p for p, _ in h5ad] if isinstance(h5ad, list) else [h5ad]
        bounded = needs_disk_backed_import([inspect_h5ad(p) for p in paths])
        if isinstance(h5ad, list):
            h5ad = str(merge_h5ads(h5ad, work / "merge"))
    qc = meta.get("qc") or {}
    store.update_job(job_id, progress=18, message="Shared QC and ambient RNA correction.")
    store.append_log(job_id, "Running the standard cellHarmony_lite QC once for the selected analysis branches.")
    _, matrix = P.cellHarmony_lite.combine_and_align_h5(
        h5_files=h5, h5ad_file=h5ad, cellharmony_ref=None, output_dir=str(work),
        export_cptt=False, export_h5ad=False, generate_umap=False, save_adata=False,
        min_genes=int(qc.get("min_genes", 500)), min_cells=int(qc.get("min_cells", 0)),
        min_counts=int(qc.get("min_counts", 1000)), mit_percent=int(qc.get("mit_percent", 15)),
        ambient_correct_cutoff="auto" if qc.get("ambient_correction") == "yes" else None,
        ambient_memory_efficient=True, concat_on_disk=True, concat_batch_size=1,
        stream_10x_inputs=True, return_adata=True, bounded_h5ad=bounded,
    )
    if not matrix.obs_names.is_unique or not matrix.var_names.is_unique:
        raise ValueError("Shared QC requires unique cell barcodes and feature identifiers; no records have been discarded.")
    # Keep the ordered roster independently of either branch for coverage gates.
    cells, genes = matrix.obs_names.copy(), matrix.var_names.copy()
    if len(cells) < 2:
        raise ValueError("Shared QC retained fewer than two cells.")
    matrix.uns["scalable_shared_qc"] = {"analysis_mode": qc.get("analysis_mode"), "qc_applied_once": True}
    prepared = work / "shared_qc_normalized.h5ad"
    P.approx_mod.ensure_h5ad_compat_for_write(matrix)
    matrix.write_h5ad(prepared, compression=compression)
    ws = getattr(matrix, "_matrix_workspace", None)
    if ws:
        ws.release_pages()
    matrix = ws = None
    gc.collect()
    return prepared, cells, genes


def _read_annotations(path):
    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        embeddings = {k: read_elem(v) for k, v in f.get("obsm", {}).items() if "umap" in k.lower()}
    return obs, var.index, embeddings


def enable_differential(meta, original, path, population_keys):
    """The same differential implementation/options for either cell-state layer."""
    fields, values = P._candidate_group_fields(path, preferred=["scalable_upload", "Library", "sample", "group"], max_categories=None)
    generated = set(population_keys) | {"supervised_accepted", "unsupervised_accepted", "original_NMF_cluster"}
    fields = [v for v in fields if v["value"] not in generated]
    values = {k: v for k, v in values.items() if k not in generated}
    sample_field = next((v["value"] for v in fields), None)
    populations = P._candidate_population_columns(path, preferred=population_keys)
    if not P._upload_profile(original).get("allow_alternate_population_fields"):
        populations = [v for v in populations if v["value"] in population_keys]
    enabled = any(len(v) >= 2 for v in values.values())
    meta["differential_options"] = {
        "enabled": enabled, "sample_names": P._job_sample_names(original),
        "sample_fields": fields, "sample_values": values, "default_sample_field": sample_field,
        "population_columns": populations, "default_population_col": meta["cluster_key"],
        "modalities": meta["modalities"]["available"], "default_modality": "rna",
        "comparison_types": ["cells", "pseudobulk"], "upload_profile": P._upload_profile(original),
    }
    meta["differential"] = P._default_differential_state(enabled, meta["cluster_key"], sample_field)
    meta["differential_history"] = {}


def association_payload(obs, filters=None):
    """Complete barcode contingency table; descriptive counts, no inferential test."""
    if not all(k in obs for k in (SUPERVISED, UNSUPERVISED)):
        raise ValueError("Cluster associations require a completed Both analysis.")
    mask = np.ones(len(obs), dtype=bool)
    for field, values in filters or []:
        if field not in obs:
            raise ValueError(f"Unknown annotation field: {field}")
        mask &= obs[field].astype(str).isin(values).to_numpy()
    subset = obs.loc[mask]
    table = pd.crosstab(subset[UNSUPERVISED].astype(str), subset[SUPERVISED].astype(str))
    rows = []
    for cluster in table.index:
        total = int(table.loc[cluster].sum())
        for state in table.columns:
            n = int(table.loc[cluster, state])
            if n:
                rows.append({"cluster": cluster, "state": state, "cells": n,
                             "cluster_cells": total, "percent": 100 * n / total})
    return {"rows": rows, "clusters": table.index.tolist(), "states": table.columns.tolist(),
            "cells": int(len(subset)), "unfiltered_cells": int(len(obs)),
            "description": "Dots count shared cell barcodes. Area shows cell count; color shows the percentage of each unsupervised cluster. Unaligned cells did not pass the reference alignment threshold. Not clustered cells received no ICGS3 assignment. These are descriptive assignments, not independent biological replicates or a statistical test."}


def _combine(job_id, store, prepared, cells, genes, supervised, unsupervised, compression):
    sup = supervised.get_job(job_id)
    unsup = unsupervised.get_job(job_id)
    sup_obs, sup_genes, sup_umaps = _read_annotations(sup["artifacts"]["combined_h5ad"])
    unsup_obs, unsup_genes, unsup_umaps = _read_annotations(unsup["artifacts"]["combined_h5ad"])
    for name, frame, features in (("supervised", sup_obs, sup_genes), ("unsupervised", unsup_obs, unsup_genes)):
        if not frame.index.is_unique or not frame.index.isin(cells).all():
            raise ValueError(f"{name} branch barcodes do not match the shared QC roster.")
        if not features.equals(genes):
            raise ValueError(f"{name} branch changed the shared feature roster; the union has not been written.")
    union = cells[cells.isin(sup_obs.index) | cells.isin(unsup_obs.index)]
    if len(union) != len(sup_obs.index.union(unsup_obs.index)):
        raise ValueError("Cell union integrity check failed.")
    workspace = Workspace(store.outputs_dir(job_id) / ".h5ad_work" / "union")
    combined = workspace.load(prepared)
    combined = workspace.subset(combined, cells.get_indexer(union))
    from ..scalable_discover.pipeline import attach_serving_gene_symbols
    attach_serving_gene_symbols(combined)
    sup_key, unsup_key = sup["cluster_key"], unsup["cluster_key"]
    combined.obs[SUPERVISED] = pd.Categorical(sup_obs[sup_key].astype(str).reindex(union).fillna(UNALIGNED))
    combined.obs[UNSUPERVISED] = pd.Categorical(unsup_obs["cluster"].astype(str).reindex(union).fillna(UNCLUSTERED))
    combined.obs["unsupervised_state"] = pd.Categorical(unsup_obs[unsup_key].astype(str).reindex(union).fillna(UNCLUSTERED))
    combined.obs["supervised_accepted"] = union.isin(sup_obs.index)
    combined.obs["unsupervised_accepted"] = union.isin(unsup_obs.index)
    # Preserve branch scores and ICGS3 annotations without overwriting source obs.
    for frame, prefix in ((sup_obs, "supervised"), (unsup_obs, "unsupervised")):
        for column in frame:
            if column not in combined.obs and column not in (sup_key, "cluster", unsup_key):
                series = frame[column].reindex(union)
                if pd.api.types.is_numeric_dtype(series):
                    combined.obs[column] = series
                else:
                    combined.obs[column] = pd.Categorical(series.astype(object).where(series.notna(), "Unavailable"))
    for name, frame, embeddings in (("supervised", sup_obs, sup_umaps), ("unsupervised", unsup_obs, unsup_umaps)):
        if "X_umap" in embeddings:
            coords = pd.DataFrame(embeddings["X_umap"], index=frame.index).reindex(union).to_numpy(dtype=np.float32)
            combined.obsm[f"X_umap_{name}"] = coords
            combined.obs[f"umap_{name}_x"] = coords[:, 0]
            combined.obs[f"umap_{name}_y"] = coords[:, 1]
    if "X_umap" in combined.obsm:
        combined.obsm["X_umap_input"] = combined.obsm["X_umap"].copy()
    combined.obsm["X_umap"] = combined.obsm["X_umap_unsupervised"].copy()
    # Preserve the standard supervised branch's embedded predictions in the
    # single download as well as its sidecars. Unaligned rows stay missing.
    with h5py.File(sup["artifacts"]["combined_h5ad"], "r") as f:
        sup_uns = read_elem(f["uns"])
        modalities = sup_uns.get("imputed_modalities") or {}
        for modality, info in modalities.items():
            key = info.get("obsm_key")
            if not key or key not in f.get("obsm", {}):
                continue  # some standard large-job modalities use sidecars only
            source = f["obsm"][key]
            if isinstance(source, h5py.Group):
                source = workspace.dataset(source)
            dtype = source.dtype if source.dtype.kind == "f" else np.float32
            values = workspace.array((len(union), source.shape[1]), dtype)
            values[:] = np.nan
            positions = union.get_indexer(sup_obs.index)
            for start in range(0, len(sup_obs), 4096):
                block = source[start:start+4096]
                values[positions[start:start+4096]] = block.toarray() if hasattr(block, "toarray") else block
            if key in combined.obsm:
                combined.obsm[f"source_{key}"] = combined.obsm[key].copy()
            # AnnData registers ndarray writers, not the memmap subclass. This
            # ndarray view keeps the same on-disk buffer without a RAM copy.
            combined.obsm[key] = np.asarray(values)
            features_key = info.get("feature_names_key")
            if features_key in sup_uns:
                if features_key in combined.uns:
                    combined.uns[f"source_{features_key}"] = combined.uns[features_key]
                combined.uns[features_key] = sup_uns[features_key]
        if modalities:
            if "imputed_modalities" in combined.uns:
                # Keep the original predictions' metadata pointing to their
                # preserved source arrays and feature identifiers.
                provenance = {m: dict(v) for m, v in combined.uns["imputed_modalities"].items()}
                for info in provenance.values():
                    for field in ("obsm_key", "feature_names_key"):
                        original_key = info.get(field)
                        source_key = f"source_{original_key}"
                        if source_key in combined.obsm or source_key in combined.uns:
                            info[field] = source_key
                combined.uns["source_imputed_modalities"] = provenance
            combined.uns["imputed_modalities"] = modalities
    combined.uns["scalable_analysis"] = {"mode": "both", "supervised_key": SUPERVISED,
        "unsupervised_key": UNSUPERVISED, "union_cells": len(union), "shared_qc_cells": len(cells),
        "supervised_cells": len(sup_obs), "unsupervised_cells": len(unsup_obs),
        "embedding_note": "Separate coordinate systems; cells not assigned in a branch have missing coordinates in that embedding."}
    path = store.outputs_dir(job_id) / "combined_with_umap_and_markers.h5ad"
    P.approx_mod.ensure_h5ad_compat_for_write(combined)
    combined.write_h5ad(path, compression=compression)
    obs = combined.obs.copy()
    associations = store.outputs_dir(job_id) / "cluster_associations.tsv"
    pd.DataFrame(association_payload(obs)["rows"]).to_csv(associations, sep="\t", index=False)
    assignments = store.outputs_dir(job_id) / "combined_assignments.tsv"
    obs.to_csv(assignments, sep="\t", index_label="CellBarcode")
    workspace.release_pages()
    combined = workspace = None
    gc.collect()
    # The supervised imputation results retain their original coverage. Add the
    # other branch's annotations by barcode, not by row position or fabricated values.
    for modality, record in (sup.get("modality_artifacts") or {}).items():
        if modality == "rna":
            continue
        for key in ("h5ad", "differential_h5ad"):
            source = record.get(key)
            if not source or not Path(source).is_file():
                continue
            with h5py.File(source, "r+") as f:
                frame = read_elem(f["obs"])
                # Pseudobulk outputs have aggregate IDs rather than cell barcodes.
                if not frame.index.isin(obs.index).all():
                    continue
                for col in (SUPERVISED, UNSUPERVISED, "unsupervised_state"):
                    frame[col] = pd.Categorical(obs[col].astype(str).reindex(frame.index))
                del f["obs"]
                write_elem(f, "obs", frame)
    artifacts = dict(sup.get("artifacts") or {})
    artifacts.pop("umap_coordinates", None)  # each selectable embedding lives in obsm
    artifacts.update(combined_h5ad=str(path), assignments=str(assignments), cluster_associations=str(associations),
                     unsupervised_marker_genes_zip=unsup["artifacts"]["marker_genes_zip"])
    layers = [
        {"key": UNSUPERVISED, "label": "Unsupervised clusters", "marker_analysis": unsup["cell_state_layers"]["layers"][1]["marker_analysis"],
         "fastcomm_analysis": unsup["cell_state_layers"]["layers"][1]["fastcomm_analysis"]},
        {"key": "unsupervised_state", "label": "Unsupervised cell states", "marker_analysis": unsup["marker_analysis"], "fastcomm_analysis": unsup["fastcomm_analysis"]},
        {"key": SUPERVISED, "label": "Supervised cell states", "marker_analysis": sup["marker_analysis"], "fastcomm_analysis": sup["fastcomm_analysis"],
         "marker_analysis_by_modality": sup.get("marker_analysis_by_modality", {})},
    ]
    result = dict(sup, cluster_key=UNSUPERVISED, marker_analysis=layers[0]["marker_analysis"],
                  marker_analysis_by_modality={"rna": layers[0]["marker_analysis"]}, fastcomm_analysis=layers[0]["fastcomm_analysis"],
                  cell_state_layers={"default": UNSUPERVISED, "layers": layers,
                                     "names": unsup["cell_state_layers"].get("names", {})}, icgs3_analysis=unsup["icgs3_analysis"],
                  artifacts=artifacts, bundle={}, analysis_mode="both", analysis_summary={
                      "qc_cells": len(cells), "union_cells": len(union), "supervised_cells": len(sup_obs),
                      "unsupervised_cells": len(unsup_obs), "unaligned_clustered_cells": len(unsup_obs.index.difference(sup_obs.index))})
    result["modality_artifacts"]["rna"] = {"h5ad": str(path)}
    return result, path


def run_analysis_workflow(job_id, store, registry_path, *, export_approx_pdfs=False, h5ad_compression="lzf"):
    original = store.get_job(job_id)
    mode = (original.get("qc") or {}).get("analysis_mode", "supervised")
    if mode not in MODES:
        raise ValueError(f"Unknown analysis mode: {mode}")
    if mode == "supervised":
        return P.run_cellharmony_pipeline(job_id, store, registry_path, export_approx_pdfs=export_approx_pdfs,
                                          h5ad_compression=h5ad_compression)
    from ..scalable_discover.pipeline import run_discover_pipeline
    compression = P._normalize_h5ad_compression(h5ad_compression)
    prepared, cells, genes = _shared_qc(job_id, store, compression)
    supervised = None
    if mode == "both":
        supervised = BranchStore(store, job_id, "supervised", prepared)
        P.run_cellharmony_pipeline(job_id, supervised, registry_path, export_approx_pdfs=export_approx_pdfs,
                                  h5ad_compression=h5ad_compression, build_bundle=False, allow_empty_alignment=True)
        gc.collect()
        shutil.rmtree(supervised.outputs_dir(job_id) / ".h5ad_work", ignore_errors=True)
    unsupervised = BranchStore(store, job_id, "unsupervised", prepared)
    from ..flask.tasks import _JobLogStream
    branch_log = _JobLogStream(unsupervised, job_id)
    try:
        with contextlib.redirect_stdout(branch_log), contextlib.redirect_stderr(branch_log):
            run_discover_pipeline(job_id, unsupervised, h5ad_compression=h5ad_compression, build_bundle=False)
    finally:
        branch_log.flush()
    shutil.rmtree(unsupervised.outputs_dir(job_id) / ".h5ad_work", ignore_errors=True)
    if mode == "both":
        result, path = _combine(job_id, store, prepared, cells, genes, supervised, unsupervised, compression)
        keys = [UNSUPERVISED, "unsupervised_state", SUPERVISED]
    else:
        result = unsupervised.get_job(job_id)
        path = Path(result["artifacts"]["combined_h5ad"])
        result["analysis_mode"] = mode
        keys = [e["key"] for e in result["cell_state_layers"]["layers"]]
    enable_differential(result, original, path, keys)
    # Restore upload configuration: branch plumbing must never become user metadata.
    protected = {"job_id", "species", "reference", "files", "qc", "created_at", "updated_at", "status", "progress",
                 "analysis_started_at", "analysis_completed_at", "analysis_duration_seconds", "worker_pid", "message"}
    store.update_job(job_id, **{k: v for k, v in result.items() if k not in protected},
                     message="Analysis complete; preparing Explore results.")
    bundle = P._build_job_bundle(store, job_id, path, result["cluster_key"], result["modality_artifacts"], result["modalities"])
    store.update_job(job_id, bundle=bundle)
    store.append_log(job_id, f"{mode.capitalize()} analysis completed; differential analysis is available for both grouping layers.")
    return result["artifacts"]
