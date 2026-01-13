from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from altanalyze3.components.cellHarmony import cellHarmony_lite
from altanalyze3.components.visualization import approximate_umap as approx_mod

from .job_manager import JobStore


def _flatten_matrix(x) -> np.ndarray:
    if sp.issparse(x):
        return np.asarray(x.todense()).ravel()
    arr = np.asarray(x)
    return arr.ravel()


def _split_uploads(files: List[Dict[str, str]], uploads_dir: Path) -> (List[str], Optional[str]):
    h5_files: List[str] = []
    h5ad_entries: List[str] = []
    for record in files:
        path = uploads_dir / record["filename"]
        suffix = path.suffix.lower()
        if suffix == ".h5ad":
            h5ad_entries.append(str(path))
        else:
            h5_files.append(str(path))
    h5ad_file = None
    if h5ad_entries:
        if len(h5ad_entries) > 1:
            raise ValueError("Multiple .h5ad uploads detected; please upload a single combined .h5ad or 10x files.")
        if h5_files:
            raise ValueError("Mixing .h5ad with other file types is not supported.")
        h5ad_file = h5ad_entries[0]
    if not h5_files and not h5ad_file:
        raise ValueError("No compatible input files detected.")
    return h5_files, h5ad_file


def _lookup_reference(species: str, reference_id: str, registry_path: Path) -> Dict:
    with registry_path.open("r", encoding="utf-8") as handle:
        registry = json.load(handle)
    for entry in registry.get("species", []):
        if entry.get("id") != species:
            continue
        for ref in entry.get("references", []):
            if ref.get("id") == reference_id:
                return ref
    raise ValueError(f"Reference '{reference_id}' for species '{species}' not found in registry.")


def _ensure_reference_fields(reference_entry: Dict) -> None:
    required = ["states_tsv", "reference_clusters_tsv", "reference_coords_tsv"]
    missing = [key for key in required if key not in reference_entry]
    if missing:
        raise ValueError(f"Reference metadata missing required keys: {', '.join(missing)}")


def run_cellharmony_pipeline(job_id: str, store: JobStore, registry_path: Path) -> Dict[str, Path]:
    """
    Execute the production cellHarmony-lite + approximate UMAP workflow for a job.

    Returns
    -------
    Dict[str, Path]
        Mapping of artifact identifiers to their on-disk locations.
    """

    meta = store.get_job(job_id)
    reference_entry = _lookup_reference(meta["species"], meta["reference"], registry_path)
    _ensure_reference_fields(reference_entry)
    store.append_log(job_id, "Reference metadata loaded.")

    uploads_dir = store.uploads_dir(job_id)
    outputs_dir = store.outputs_dir(job_id)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    h5_files, h5ad_file = _split_uploads(meta["files"], uploads_dir)
    qc = meta.get("qc", {})

    store.append_log(job_id, "Running cellHarmony_lite pipeline.")
    ordered_df = cellHarmony_lite.combine_and_align_h5(
        h5_files=h5_files,
        h5ad_file=h5ad_file,
        cellharmony_ref=reference_entry["states_tsv"],
        output_dir=str(outputs_dir),
        export_cptt=False,
        export_h5ad=False,
        min_genes=int(qc.get("min_genes", 200)),
        min_cells=int(qc.get("min_cells", 3)),
        min_counts=int(qc.get("min_counts", 500)),
        mit_percent=int(qc.get("mit_percent", 10)),
        generate_umap=False,
        save_adata=True,
        unsupervised_cluster=False,
        alignment_mode="classic",
        min_alignment_score=None,
        gene_translation_file=None,
        metacell_align=False,
        ambient_correct_cutoff=None,
    )

    assignments_path = outputs_dir / "cellHarmony_lite_assignments.txt"
    combined_h5ad_path = outputs_dir / "combined_with_umap_and_markers.h5ad"
    if not assignments_path.exists():
        raise FileNotFoundError("cellHarmony-lite assignments file missing.")
    if not combined_h5ad_path.exists():
        raise FileNotFoundError("Combined AnnData output missing; ensure save_adata=True.")

    query_cluster_key = Path(reference_entry["states_tsv"]).stem
    reference_cluster_key = reference_entry.get("cluster_key", query_cluster_key)

    store.append_log(job_id, "Running approximate UMAP placement.")
    reference_adata = approx_mod._load_reference_from_tsv(
        reference_entry["reference_coords_tsv"],
        reference_entry["reference_clusters_tsv"],
        umap_key="X_umap",
        cluster_key=reference_cluster_key,
    )
    query_adata = approx_mod._load_query_from_tsv(
        assignments_path,
        cluster_key=query_cluster_key,
    )
    approx_result = approx_mod.approximate_umap(
        query=query_adata,
        reference=reference_adata,
        query_cluster_key=query_cluster_key,
        reference_cluster_key=reference_cluster_key,
        umap_key="X_umap",
        jitter=0.05,
        num_reference_cells=1,
    )

    prefix = outputs_dir / "approximate"
    prefix.mkdir(parents=True, exist_ok=True)
    text_outputs = approx_result.write_text_outputs(prefix)
    annotated_pdf, plain_pdf = approx_result.write_comparison_pdf(
        reference_adata=reference_adata,
        umap_key="X_umap",
        reference_cluster_key=reference_cluster_key,
        query_cluster_key=query_cluster_key,
        output_path=prefix.with_name(prefix.name + "-comparison.pdf"),
    )

    artifacts = {
        "assignments": assignments_path,
        "combined_h5ad": combined_h5ad_path,
        "umap_coordinates": Path(text_outputs["coordinates"]),
        "umap_augmented": Path(text_outputs["augmented"]),
        "umap_placeholder_expression": Path(text_outputs["placeholder_expression"]),
        "umap_pdf": Path(annotated_pdf),
        "umap_pdf_plain": Path(plain_pdf),
    }
    for key, path in artifacts.items():
        store.add_artifact(job_id, key, path)

    store.update_job(
        job_id,
        cluster_key=query_cluster_key,
        reference_cluster_key=reference_cluster_key,
        reference_coords_tsv=reference_entry["reference_coords_tsv"],
        reference_clusters_tsv=reference_entry["reference_clusters_tsv"],
        message="Approximate UMAP completed.",
    )
    store.append_log(job_id, "cellHarmony-lite pipeline finished.")
    return artifacts
