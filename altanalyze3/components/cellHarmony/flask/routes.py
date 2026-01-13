from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from flask import (
    Blueprint,
    current_app,
    jsonify,
    render_template,
    request,
    send_file,
)
from werkzeug.utils import secure_filename

from altanalyze3.components.visualization import approximate_umap as approx_mod

from .job_manager import JobStore
from .tasks import JobRunner


main_bp = Blueprint("main", __name__)
api_bp = Blueprint("api", __name__, url_prefix="/api")


def register_routes(app):
    app.register_blueprint(main_bp)
    app.register_blueprint(api_bp)


def _load_reference_registry() -> Dict:
    registry_path = Path(current_app.config["REFERENCE_REGISTRY"])
    if not registry_path.exists():
        return {"species": []}
    with registry_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _allowed_file(filename: str) -> bool:
    return "." in filename and filename.rsplit(".", 1)[1].lower() in current_app.config["ALLOWED_EXTENSIONS"]


def _job_resources() -> (JobStore, JobRunner):
    return current_app.job_store, current_app.job_runner


def _flatten_expr(values) -> np.ndarray:
    if sp.issparse(values):
        return np.asarray(values.todense()).ravel()
    return np.asarray(values).ravel()


def _load_reference_adata(meta: Dict) -> Optional[ad.AnnData]:
    coords = meta.get("reference_coords_tsv")
    clusters = meta.get("reference_clusters_tsv")
    if not coords or not clusters:
        return None
    cluster_key = meta.get("reference_cluster_key") or meta.get("cluster_key")
    return approx_mod._load_reference_from_tsv(
        coords,
        clusters,
        umap_key="X_umap",
        cluster_key=cluster_key,
    )


@main_bp.route("/")
def index():
    registry = _load_reference_registry()
    return render_template("index.html", registry=json.dumps(registry))


@api_bp.get("/meta/species")
def meta_species():
    return jsonify(_load_reference_registry())


@api_bp.post("/jobs")
def create_job():
    store, _ = _job_resources()
    species = request.form.get("species")
    reference = request.form.get("reference")
    soupx_option = request.form.get("soupx_option")

    if not species or not reference:
        return jsonify({"error": "Species and reference are required."}), 400

    sample_names = request.form.getlist("sample_names[]")
    uploaded_files = request.files.getlist("files[]")

    if not sample_names or not uploaded_files or len(sample_names) != len(uploaded_files):
        return jsonify({"error": "Each uploaded file must include a matching sample name."}), 400
    if len(uploaded_files) > current_app.config["MAX_FILES_PER_JOB"]:
        return jsonify({"error": f"Maximum {current_app.config['MAX_FILES_PER_JOB']} files are allowed per job."}), 400

    normalized_samples: List[str] = []
    for name in sample_names:
        if not name:
            return jsonify({"error": "Sample names cannot be empty."}), 400
        if name in normalized_samples:
            return jsonify({"error": f"Duplicate sample name '{name}' detected."}), 400
        normalized_samples.append(name)

    metadata = store.create_job(species, reference, soupx_option, files=[])
    job_id = metadata["job_id"]
    uploads_dir = store.uploads_dir(job_id)

    records = []
    for sample, file in zip(normalized_samples, uploaded_files):
        if file.filename == "":
            return jsonify({"error": "One of the selected files has no filename."}), 400
        if not _allowed_file(file.filename):
            return jsonify({"error": f"Unsupported file extension for {file.filename}."}), 400
        dest_name = f"{sample}_{secure_filename(file.filename)}"
        dest_path = uploads_dir / dest_name
        file.save(dest_path)
        records.append({"sample_name": sample, "filename": dest_name, "size": os.path.getsize(dest_path)})

    store.update_job(job_id, files=records, message="Upload complete. Configure QC to proceed.")
    return jsonify({"job_id": job_id, "status": "uploaded"})


@api_bp.post("/jobs/<job_id>/qc")
def update_qc(job_id: str):
    store, _ = _job_resources()
    if not store.job_exists(job_id):
        return jsonify({"error": "Job not found."}), 404
    payload = request.get_json(force=True)
    qc_payload = {
        "min_genes": int(payload.get("min_genes", 200)),
        "min_counts": int(payload.get("min_counts", 500)),
        "min_cells": int(payload.get("min_cells", 3)),
        "mit_percent": int(payload.get("mit_percent", 15)),
    }
    meta = store.update_job(job_id, qc=qc_payload, message="QC parameters saved.")
    return jsonify({"job_id": job_id, "qc": meta["qc"]})


@api_bp.post("/jobs/<job_id>/run")
def run_job(job_id: str):
    store, runner = _job_resources()
    if not store.job_exists(job_id):
        return jsonify({"error": "Job not found."}), 404
    runner.submit(job_id)
    store.append_log(job_id, "Job queued by user request.")
    store.update_job(job_id, message="Job submitted to worker.", status="queued", progress=15)
    return jsonify({"job_id": job_id, "status": "queued"})


@api_bp.get("/jobs/<job_id>/status")
def status(job_id: str):
    store, _ = _job_resources()
    if not store.job_exists(job_id):
        return jsonify({"error": "Job not found."}), 404
    meta = store.get_job(job_id)

    log_path = store.logs_dir(job_id) / "pipeline.log"
    log_tail: List[str] = []
    if log_path.exists():
        with log_path.open("r", encoding="utf-8") as handle:
            lines = handle.readlines()
            log_tail = lines[-15:]
    meta["log_tail"] = log_tail
    return jsonify(meta)


@api_bp.get("/jobs/<job_id>/umap")
def umap(job_id: str):
    store, _ = _job_resources()
    if not store.job_exists(job_id):
        return jsonify({"error": "Job not found."}), 404
    meta = store.get_job(job_id)
    artifacts = meta.get("artifacts", {})
    umap_path = artifacts.get("umap_coordinates")
    assignments_path = artifacts.get("assignments")
    if not umap_path or not Path(umap_path).exists():
        return jsonify({"error": "UMAP output not available yet."}), 404
    if not assignments_path or not Path(assignments_path).exists():
        return jsonify({"error": "cellHarmony assignments missing; rerun the job."}), 404

    query_df = pd.read_csv(umap_path, sep="\t")
    barcode_col = query_df.columns[0]
    query_df = query_df.rename(
        columns={
            barcode_col: "CellBarcode",
            "umap_0": "UMAP1",
            "umap_1": "UMAP2",
        }
    )

    assignments_df = pd.read_csv(assignments_path, sep="\t")
    cluster_key = meta.get("cluster_key")
    if not cluster_key or cluster_key not in assignments_df.columns:
        cluster_candidates = [col for col in assignments_df.columns if col not in {"CellBarcode", "Similarity"}]
        if not cluster_candidates:
            return jsonify({"error": "Cluster column missing in assignments file."}), 500
        cluster_key = cluster_candidates[0]
    merged = pd.merge(assignments_df[["CellBarcode", cluster_key]], query_df, on="CellBarcode", how="inner")
    merged["Population"] = merged[cluster_key].astype(str)

    query_points = [
        {
            "barcode": row.CellBarcode,
            "population": row.Population,
            "x": float(row.UMAP1),
            "y": float(row.UMAP2),
        }
        for row in merged.itertuples()
    ]

    reference_points = []
    ref_adata = _load_reference_adata(meta)
    if ref_adata is not None:
        ref_cluster_key = meta.get("reference_cluster_key") or cluster_key
        coords = ref_adata.obsm["X_umap"]
        labels = ref_adata.obs[ref_cluster_key].astype(str).tolist()
        for barcode, (x, y), population in zip(ref_adata.obs_names, coords, labels):
            reference_points.append(
                {
                    "barcode": str(barcode),
                    "population": population,
                    "x": float(x),
                    "y": float(y),
                }
            )

    return jsonify({"reference": reference_points, "query": query_points})


@api_bp.get("/jobs/<job_id>/expression")
def expression(job_id: str):
    gene = request.args.get("gene")
    if not gene:
        return jsonify({"error": "Gene parameter is required."}), 400

    store, _ = _job_resources()
    if not store.job_exists(job_id):
        return jsonify({"error": "Job not found."}), 404
    meta = store.get_job(job_id)
    artifacts = meta.get("artifacts", {})
    h5ad_path = artifacts.get("combined_h5ad")
    if not h5ad_path or not Path(h5ad_path).exists():
        return jsonify({"error": "Combined AnnData output unavailable for this job."}), 404

    adata = ad.read_h5ad(h5ad_path)
    cluster_key = meta.get("cluster_key")
    if not cluster_key or cluster_key not in adata.obs.columns:
        return jsonify({"error": "Cluster assignments missing from AnnData output."}), 500
    if gene not in adata.var_names:
        return jsonify({"error": f"Gene '{gene}' not found in AnnData output."}), 404

    values = _flatten_expr(adata[:, gene].X)
    populations = adata.obs[cluster_key].astype(str).values

    scatter_data = [
        {"population": pop, "value": float(val)}
        for pop, val in zip(populations, values)
    ]
    violin_data = []
    for pop in sorted(pd.unique(populations)):
        mask = populations == pop
        violin_data.append(
            {
                "population": pop,
                "values": [float(v) for v in values[mask]],
            }
        )

    return jsonify({"gene": gene, "scatter": scatter_data, "violin": violin_data})
    if not umap_path.exists():
        return jsonify({"error": "UMAP output not available yet."}), 404

    reference_points = []
    query_points = []
    with umap_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            point = {
                "barcode": row["CellBarcode"],
                "population": row["Population"],
                "x": float(row["UMAP1"]),
                "y": float(row["UMAP2"]),
            }
            if row["Set"] == "reference":
                reference_points.append(point)
            else:
                query_points.append(point)
    return jsonify({"reference": reference_points, "query": query_points})


@api_bp.get("/jobs/<job_id>/download/<artifact>")
def download(job_id: str, artifact: str):
    store, _ = _job_resources()
    if not store.job_exists(job_id):
        return jsonify({"error": "Job not found."}), 404
    meta = store.get_job(job_id)
    artifact_map = meta.get("artifacts", {})
    if artifact not in artifact_map:
        return jsonify({"error": f"Artifact '{artifact}' unavailable."}), 404
    path = Path(artifact_map[artifact])
    if not path.exists():
        return jsonify({"error": "Artifact missing on disk."}), 404
    return send_file(path, as_attachment=True)
