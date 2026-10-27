"""Deterministic methods and notebook export; never load expression or fit a model.

Execution snapshots are captured by workers, not reconstructed from the download
server's environment. Historical runs retain explicit gaps in their provenance.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import re
import subprocess
import shutil
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path


SCHEMA_VERSION = 1
QC_DEFAULTS = {
    "analysis_mode": "supervised", "existing_state_column": None, "existing_counts_source": "auto", "existing_embedding_key": None, "umap_fit_mode": "landmark", "max_k": None,
    "min_genes": 500, "min_counts": 1000, "min_cells": 0, "mit_percent": 15,
    "align_cutoff": 0.4, "ambient_correction": "no", "impute_modality": "none",
    "impute_modalities": None, "marker_render_heatmap": False,
    "marker_write_svg": True, "marker_heatmap_dpi": None, "marker_cells_per_cluster": 100,
}
PACKAGES = ("altanalyze3", "numpy", "scipy", "pandas", "anndata", "h5py",
            "scanpy", "scikit-learn", "umap-learn", "pynndescent", "numba",
            "matplotlib", "plotly", "fastapi", "statsmodels")
CITATIONS = {
    "cellHarmony": {"authors": "DePasquale EAK et al.", "year": 2019,
        "title": "cellHarmony: cell-level matching and holistic comparison of single-cell transcriptomes",
        "journal": "Nucleic Acids Research 47:e138", "doi": "10.1093/nar/gkz789"},
    "ICGS2": {"authors": "Venkatasubramanian M, Chetal K, Schnell DJ, Atluri G, Salomonis N",
        "year": 2020, "title": "Resolving single-cell heterogeneity from hundreds of thousands of cells through sequential hybrid clustering and NMF",
        "journal": "Bioinformatics 36:3773–3780", "doi": "10.1093/bioinformatics/btaa201"},
    "GO-Elite": {"authors": "Zambon AC et al.", "year": 2012,
        "title": "GO-Elite: a flexible solution for pathway and ontology over-representation",
        "journal": "Bioinformatics 28:2209–2210", "doi": "10.1093/bioinformatics/bts366"},
    "ChromLinker": {"authors": "Ferchen K, Zhang X et al.", "year": 2025,
        "title": "A unified multimodal single-cell framework reveals a discrete state model of hematopoiesis in mice",
        "journal": "Nature Immunology 26:2086–2099", "doi": "10.1038/s41590-025-02307-3"},
    "Scanpy": {"authors": "Wolf FA, Angerer P, Theis FJ", "year": 2018,
        "title": "SCANPY: large-scale single-cell gene expression data analysis",
        "journal": "Genome Biology 19:15", "doi": "10.1186/s13059-017-1382-0"},
    "UMAP": {"authors": "McInnes L, Healy J, Melville J", "year": 2018,
        "title": "UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction",
        "journal": "arXiv:1802.03426", "doi": "10.48550/arXiv.1802.03426"},
    "CellChatDB": {"authors": "Jin S et al.", "year": 2021,
        "title": "Inference and analysis of cell-cell communication using CellChat",
        "journal": "Nature Communications 12:1088", "doi": "10.1038/s41467-021-21246-9"},
}


def _json(value):
    return json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False, default=str) + "\n"


def file_identity(path):
    path = Path(path)
    if not path.is_file():
        return {"filename": path.name, "status": "not available at execution"}
    sha = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            sha.update(chunk)
    return {"filename": path.name, "bytes": path.stat().st_size, "sha256": sha.hexdigest()}


def git_execution_identity():
    """Capture the checkout at execution only; never backfill historical jobs."""
    root = Path(__file__).resolve().parents[3]
    def git(*args):
        return subprocess.check_output(["git", "-C", str(root), *args],
            stderr=subprocess.DEVNULL, text=True, timeout=5).strip()
    try:
        return {"commit": git("rev-parse", "HEAD"),
                "commit_date": git("show", "-s", "--format=%cI", "HEAD"),
                "working_tree_modified": bool(git("status", "--porcelain", "--untracked-files=normal"))}
    except (OSError, subprocess.SubprocessError):
        return {"status": "Git identity unavailable at execution"}


def capture_execution(store, job_id, *, application="scALABLE", reference=None,
                      key="pipeline", effective=None):
    """Call at execution, including before a potentially failed stage starts."""
    from ..model_registry.registry import describe_application
    meta = store.get_job(job_id)
    versions = {"python": platform.python_version()}
    for package in PACKAGES:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not installed as a distribution"
    record = {"schema_version": SCHEMA_VERSION,
              "captured_at": datetime.now(timezone.utc).isoformat(),
              "application": application, "runtime_versions": versions,
              "git_checkout": git_execution_identity(),
              "application_method": describe_application(application),
              "species": meta.get("species"), "reference_id": meta.get("reference"),
              "qc": meta.get("qc") or {}, "interface_defaults": dict(QC_DEFAULTS),
              "effective": effective or {}}
    package_root = Path(__file__).resolve().parents[1]
    record["method_files"] = {name: file_identity(package_root / name) for name in (
        "cellHarmony/cellHarmony_lite.py", "cellHarmony/cellHarmony_differential.py",
        "cellHarmony/differential_sampling.py", "cellHarmony/flask/pipeline.py",
        "clustering/ICGS.py", "ambient_rna/ambient_subtract.py", "visualization/approximate_umap.py",
        "visualization/marker_heatmap_h5ad.py", "cellHarmony/analysis_record.py")}
    if reference:
        record["reference"] = dict(reference)
        record["reference_files"] = {name: file_identity(reference[name]) for name in
            ("states_tsv", "reference_clusters_tsv", "reference_coords_tsv") if reference.get(name)}
    if key == "pipeline":
        record["inputs"] = [{"filename": Path(item.get("filename", "")).name,
                              "sample_name": item.get("sample_name"),
                              "bytes": (store.uploads_dir(job_id) / item.get("filename", "")).stat().st_size
                              if (store.uploads_dir(job_id) / item.get("filename", "")).is_file() else None}
                             for item in meta.get("files") or []]
    records = dict(meta.get("analysis_records") or {})
    records[key] = record
    store.update_job(job_id, analysis_records=records)
    return record


def record_effective(store, job_id, key, **effective):
    """Append dispatch decisions to the original execution snapshot."""
    records = dict(store.get_job(job_id).get("analysis_records") or {})
    record = dict(records.get(key) or {})
    record["effective"] = dict(record.get("effective") or {}, **effective)
    records[key] = record
    store.update_job(job_id, analysis_records=records)


def _load_json(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def build_manifest(meta, job_dir):
    """Use persisted evidence only. Never inspect H5AD arrays or model pickles."""
    job_dir = Path(job_dir)
    sources = {"main": meta}
    for path in sorted((job_dir / "outputs" / "branches").glob("*/**/job.json")):
        if not path.resolve().is_relative_to(job_dir.resolve()):
            continue
        branch = _load_json(path)
        if branch:
            sources[path.relative_to(job_dir).as_posix()] = branch
    records = {name: source.get("analysis_records") or {} for name, source in sources.items()}
    pipeline_record = records["main"].get("pipeline") or {}
    models = meta.get("model_versions") or {}
    provenance = {}
    for name, source in sources.items():
        raw = (source.get("artifacts") or {}).get("model_provenance")
        if raw:
            # Export only a job's own artifact; never follow arbitrary metadata paths.
            path = Path(raw).resolve()
            if path.is_relative_to(job_dir.resolve()):
                provenance[name] = _load_json(path)
    logs = sorted(set(job_dir.glob("logs/*.log")) |
                  set(job_dir.glob("outputs/ICGS3/logs/*.log")) |
                  set(job_dir.glob("outputs/branches/*/*/logs/*.log")) |
                  set(job_dir.glob("outputs/branches/*/*/outputs/ICGS3/logs/*.log")))
    logs = [path for path in logs if path.is_file() and path.resolve().is_relative_to(job_dir.resolve())]
    evidence, scale_evidence = [], []
    for path in logs:
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if any(token in line for token in ("[params]", "[ICGS3] parameters:", "[ICGS3] cli equivalent:")):
                    evidence.append({"log": path.relative_to(job_dir).as_posix(), "text": line.rstrip()})
                if any(token in line for token in ("[scale]", "[ambient]", "skipping QC", "[qc]", "[qc-only]",
                        "[mem]", "QC filter", "Cells remaining after", "Auto-selected rho",
                        "Ambient RNA correction", "ambient RNA correction", "Starting library",
                        "Normalization steps:", "adata shape:", "[INFO] Applied min_alignment_score")):
                    scale_evidence.append({"log": path.relative_to(job_dir).as_posix(), "text": line.rstrip()})
    comparisons = dict(meta.get("differential_history") or {})
    current = meta.get("differential") or {}
    if current and (current.get("run_id") or current.get("status") in {"queued", "processing", "completed", "failed"}):
        comparisons[current.get("run_id") or "current"] = current
    version_sources = {"execution_snapshots": records, "model_provenance": provenance}
    # Model runtime versions are execution evidence for those stages, not for all stages.
    if not pipeline_record:
        version_sources["legacy_model_runtime_versions"] = {
            k: v.get("runtime_versions", {}) for k, v in models.items() if isinstance(v, dict)}
    return {
        "schema_version": SCHEMA_VERSION, "generator": "deterministic templates; no LLM",
        "job_id": meta.get("job_id"), "status": meta.get("status"), "message": meta.get("message"),
        "analysis_mode": meta.get("analysis_mode") or (meta.get("qc") or {}).get("analysis_mode") or "supervised",
        "species": meta.get("species"), "reference_id": meta.get("reference"),
        "reference": pipeline_record.get("reference") or {},
        "inputs": pipeline_record.get("inputs") or [{"filename": Path(i.get("filename", "")).name,
                     "sample_name": i.get("sample_name")} for i in meta.get("files") or []],
        "input_memory": meta.get("input_memory") or {},
        "existing_annotations": meta.get("existing_annotations") or {},
        "qc": pipeline_record.get("qc") or meta.get("qc") or {},
        "interface_defaults_at_execution": pipeline_record.get("interface_defaults") or {},
        "model_versions": models, "software_provenance": version_sources,
        "analysis_records": records,
        "analysis_started_at": meta.get("analysis_started_at"),
        "analysis_completed_at": meta.get("analysis_completed_at"),
        "analysis_duration_seconds": meta.get("analysis_duration_seconds"),
        "layers": meta.get("cell_state_layers") or {},
        "marker_analysis": meta.get("marker_analysis") or {},
        "marker_analysis_by_modality": meta.get("marker_analysis_by_modality") or {},
        "fastcomm_analysis": meta.get("fastcomm_analysis") or {},
        "differential_runs": comparisons, "parameter_evidence": evidence,
        "preprocessing_evidence": scale_evidence,
        "provenance_gaps": [] if pipeline_record else [
            "No execution snapshot exists for this historical job. Download-time software versions and current reference metadata are not substituted for execution-time versions.",
            "Historical interface defaults are not recorded; user/default differences cannot be established retrospectively."],
    }, logs


def _display(value):
    if value is None:
        return "Not recorded"
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if isinstance(value, (list, tuple)):
        return "; ".join(_display(item) for item in value) or "None"
    if isinstance(value, dict):
        return "; ".join(f"{key.replace('_', ' ')}: {_display(item)}" for key, item in value.items()) or "None"
    return str(value)


def _table(headers, rows):
    def cell(value):
        return _display(value).replace("\n", " ").replace("|", "&#124;")
    return ["", "| " + " | ".join(headers) + " |",
            "|" + "|".join("---" for _ in headers) + "|"] + [
            "| " + " | ".join(cell(value) for value in row) + " |" for row in rows] + [""]


def _software_methods(m):
    result = []
    for source, records in m["analysis_records"].items():
        for stage, record in records.items():
            versions = record.get("runtime_versions") or {}
            checkout = record.get("git_checkout") or {}
            if versions.get("altanalyze3"):
                result.append(f"{stage} ({source}) used altanalyze3 version {versions['altanalyze3']}.")
            if checkout.get("commit"):
                result.append(f"Git commit `{checkout['commit']}` (commit date {checkout.get('commit_date', 'not recorded')}); "
                              f"working tree modified at execution: {_display(checkout.get('working_tree_modified'))}. "
                              "Execution-time code hashes identify the files used, including local modifications.")
    return result


def model_methods(modality, record):
    component = record.get("component") or {"adt":"rna2adt", "lipids":"rna2lipid", "grn":"rna2grn", "grn_tf":"rna2grn"}.get(modality, modality)
    version = record.get("model_display_id") or record.get("model_version") or record.get("model_version_id") or "not recorded"
    lines = [f"{component} (altanalyze3) used model version **{version}**. These outputs are inferred molecular features."]
    fields = (("model_class", "Model class"), ("architecture", "Architecture"),
        ("model_variant", "Tissue / variant"), ("release_id", "Release identifier"),
        ("model_gene_count", "RNA input features"), ("output_adt_count", "Protein outputs"),
        ("output_lipid_count", "Lipid outputs"), ("n_output_edges", "Regulatory edge outputs"),
        ("n_tfs", "Transcription factors"), ("input_rows", "Input cells"), ("output_rows", "Output cells"),
        ("expression_scale", "Expression scale"), ("prediction_scale", "Prediction scale"),
        ("RNA_input_log_base", "RNA input logarithm base"), ("RNA_model_log_base", "Model RNA logarithm base"),
        ("RNA_log_base_conversion_factor", "RNA log-base conversion factor"),
        ("abundance_units", "Abundance units"), ("pseudobulk_statistic", "Group aggregation statistic"),
        ("n_groups", "Aggregated groups"), ("model_version_id", "Exact model identity"))
    lines += _table(["Model property", "Recorded value"], [(label, record[key]) for key, label in fields if key in record])
    training = [(key.replace("_", " ").capitalize(), record[key]) for key in (
        "training_datasets", "training_provenance", "training_summary", "validation", "validation_summary", "evaluation", "confidence", "declared_method") if key in record]
    if training:
        lines += _table(["Training / validation record", "Recorded value"], training)
    else:
        lines.append("Release-specific training datasets and independent validation metrics were not retained in this job's record. Predictive accuracy cannot be established from the version identifier alone.")
    return "\n".join(lines)


def _preprocessing_methods(m):
    qc = m["qc"]; defaults = m.get("interface_defaults_at_execution") or {}
    labels = {"min_genes":"Minimum detected genes per cell", "min_counts":"Minimum counts per cell",
        "min_cells":"Minimum cells per gene", "mit_percent":"Maximum mitochondrial counts (%)",
        "align_cutoff":"Minimum alignment similarity", "ambient_correction":"Ambient RNA correction"}
    rows = [(label, qc.get(key), defaults.get(key),
             ("yes" if qc[key] != defaults[key] else "no") if key in defaults else "unknown")
            for key, label in labels.items() if key in qc and not (m.get("analysis_mode") == "existing" and key == "align_cutoff")]
    lines = ["Quality-control thresholds were configured as follows." if m.get("analysis_mode") == "existing" else
             "Quality-control thresholds and the reference-alignment cutoff were configured as follows."]
    lines += _table(["Setting", "Recorded value", "Interface default at execution", "Changed from default"], rows)
    evidence = "\n".join(e["text"] for e in m.get("preprocessing_evidence", []))
    if ("skipping QC" in evidence or "Skipping QC filtering" in evidence) and "[qc] Use Existing:" not in evidence:
        lines.append("Count-based QC filtering was skipped because the input was already normalized and log-transformed.")
    if 'using verified layers' in evidence and 'scale_data=True through the standard MarkerFinder hook' in evidence:
        lines.append("Imported expression failed the strict MarkerFinder constant-depth validation on the supplied gene panel. MarkerFinder therefore used the verified counts layer through its standard scale_data=True route (log2 depth normalization with the native median cell-total scale factor). The supplied RNA expression and saved embedding were retained, and scaling validation remained enabled.")
    if 'verified raw.X counts' in evidence:
        lines.append("Raw counts were recovered from the supplied AnnData .raw matrix after confirming the complete matching RNA gene panel; the source raw matrix and normalized RNA were retained.")
    if 'Mitochondrial filtering unavailable' in evidence:
        lines.append("The requested mitochondrial filter could not be evaluated because mitochondrial genes were absent from the supplied gene panel; no zero mitochondrial percentage was fabricated.")
    retained = list(dict.fromkeys(re.findall(r"Cells remaining after (.+?) filtering: (\d+)", evidence)))
    if retained:
        lines += _table(["Applied QC filter", "Cells retained"], [(label, f"{int(count):,}") for label, count in retained])
    shapes = re.findall(r"adata shape: \((\d+), (\d+)\)", evidence)
    if shapes:
        cells, genes = shapes[0]
        lines.append(f"The imported expression matrix contained {int(cells):,} cells and {int(genes):,} genes.")
    rho = list(dict.fromkeys(re.findall(r"Auto-selected rho for library '([^']+)': ([\d.]+)", evidence)))
    if rho or "Ambient RNA correction complete" in evidence:
        lines.append("Ambient RNA correction was performed with AmbientSubtract (altanalyze3; unpublished, manuscript in preparation). The correction fraction was selected separately for each library when automatic selection was used.")
        if rho:
            lines += _table(["Library", "Selected ambient correction fraction (rho)"], rho)
    elif str(qc.get("ambient_correction", "no")).lower() == "yes":
        lines.append("AmbientSubtract (altanalyze3) was requested; completion and library-specific correction fractions were not established by the retained preprocessing record.")
    else:
        lines.append("Ambient RNA correction was not requested.")
    if "skipping normalize_total and log1p" in evidence:
        lines.append("The input expression scale was retained; total-count normalization and log1p were not repeated.")
    elif re.search(r"Normalization steps:.*100%", evidence):
        lines.append("Expression was normalized to 10,000 counts per cell and transformed using the natural logarithm of one plus the normalized count (log1p).")
    else:
        lines.append("The raw-count workflow uses total-count normalization to 10,000 counts per cell followed by natural-log log1p; completion of these transformations was not recorded for this run.")
    alignment = list(dict.fromkeys(re.findall(r"Applied min_alignment_score=([\d.]+)\. Excluded (\d+) cells, kept (\d+)", evidence)))
    for cutoff, excluded, kept in alignment:
        lines.append(f"At an alignment similarity threshold of {cutoff}, {int(kept):,} cells were accepted and {int(excluded):,} cells fell below the threshold." +
            (" Below-threshold cells are retained in the combined Both-mode result as Unaligned." if m["analysis_mode"] == "both" else ""))
    return lines


def render_methods(manifest):
    m = manifest
    mode, qc = m["analysis_mode"], m["qc"]
    paragraphs = ["# scALABLE analysis methods", "",
        "Methods were compiled from the saved analysis record. Parameters and versions refer to this run; the complete execution record and original logs are retained separately in the download archive.",
        "", "## Analysis identity and input", "",
        f"Single-cell expression data from **{m['species']}** were submitted to the **{mode}** scALABLE (altanalyze3) workflow, using reference `{m['reference_id']}`. Recorded analysis status: **{m['status']}**.",
        "Inputs retain their recorded library/sample names. For an annotated H5AD, the selected observation field defines the sample identities used in each differential comparison.",
        *_table(["Library / sample", "Input file"], [(i.get("sample_name"), i.get("filename")) for i in m["inputs"]]), "", "## Selected reference", ""]
    ref = m.get("reference") or {}
    if mode == "existing":
        paragraphs += ["Supplied annotations define the cell populations. Any reference selected below is used only to select optional imputation models; it was not used for alignment."]
    if ref:
        paragraphs += [f"Selected reference: {ref.get('label', m['reference_id'])}."]
        if ref.get("study_url"):
            paragraphs += [f"Reference citation: [{ref.get('study_citation') or ref['study_url']}]({ref['study_url']})."]
        else:
            paragraphs += ["No published study citation is recorded for this reference; no citation to a different reference version has been substituted."]
    elif mode == "existing":
        paragraphs += ["No reference was required to reuse the supplied annotations and embeddings."]
    else:
        paragraphs += ["Reference metadata/version at execution was not recorded. The saved identifier and resource paths in run_manifest.json are the available record."]
    paragraphs += ["Reference resource identities are retained in Software versions and reproducibility and in `run_manifest.json`.",
        "", "## Quality control, ambient RNA correction and expression scale", ""]
    paragraphs += _preprocessing_methods(m)
    paragraphs += ["", "## Cell-state inference and embedding", ""]
    if mode == "existing":
        existing = m.get("existing_annotations") or {}
        paragraphs += [f"Use Existing retained the supplied cell-state/cluster labels from obs[{qc.get('existing_state_column')!r}]. "
            "Requested QC used counts when available, preserving already normalized RNA without a second normalization. "
            "No supervised reference alignment, ICGS3 clustering, UMAP fitting or reference projection was performed. "
            "Supplied annotations and embeddings were retained for QC-passed cells. Missing labels were represented as Unassigned in a separate serving column. "
            f"Annotation column selection was recorded as {existing.get('method', 'user selection')}; "
            "the local model, when used, selected an existing column and did not generate per-cell annotations. "
            f"The requested counts source was {qc.get('existing_counts_source', 'auto')!r}; "
            f"the saved layout source was {qc.get('existing_embedding_key') or 'selected from supplied metadata'!r}. "
            "Downstream analyses used the selected grouping. Differential comparisons retained all selected cells because external annotations need not have native assignment scores."]
    if mode in {"supervised", "both"}:
        paragraphs += ["The supervised branch uses cellHarmony-lite (altanalyze3) reference alignment. Query cells are compared with reference cell-state profiles over the reference marker panel; the native pipeline uses cosine similarity and the requested alignment cutoff. This implementation builds on cellHarmony [cellHarmony], but its centroid alignment and approximate reference placement should not be described as the original full cell-to-cell matching algorithm.",
            "The approximate supervised UMAP places accepted query cells using reference coordinates rather than fitting a new UMAP to all query cells. The recorded placement parameters specify the reference-cell count and jitter; a seed is only claimed if recorded."]
    if mode in {"unsupervised", "both"}:
        paragraphs += ["The unsupervised branch invokes ICGS3 (altanalyze3), through the native `ICGS.run_icgs3` implementation. ICGS3 is unpublished (manuscript in preparation); ICGS2 [ICGS2] is its published methodological precursor, not the version executed here. The workflow includes variable-program discovery, downsampling, NMF, MarkerFinder refinement and SVM assignment. Resolved ICGS3 configuration, including target rank, downsampling, normalization, UMAP features/backend, landmark limits, neighbors and seeds, is retained in the execution snapshots when available. The accelerated UMAP setting concerns embedding; it is not an instruction to discard cells from the returned dataset."]
    if mode in {"unsupervised", "both"}:
        embedding_settings = [("UMAP fitting", qc.get("umap_fit_mode")),
            ("Requested maximum K", "Automatic (no target K)" if "max_k" in qc and qc["max_k"] is None else qc.get("max_k"))]
        paragraphs += _table(["Embedding / clustering setting", "Recorded value"], embedding_settings)
    if mode == "both":
        paragraphs += ["Both branches share upstream QC. The combined H5AD retains the union of accepted barcodes with separate supervised and unsupervised annotation/coordinate layers. Cells without a supervised match are represented as Unaligned rather than silently removed from the combined result. Cluster-association plots summarize barcode intersections; they are not independent biological measurements."]
    rank_selections = {source: record["icgs3"]["effective"]["nmf_rank_selection"]
                       for source, record in m["analysis_records"].items()
                       if (record.get("icgs3", {}).get("effective", {}).get("nmf_rank_selection"))}
    if rank_selections:
        paragraphs += ["An unusable manual K is retried through ICGS3's unbiased automatic feature/rank/NMF selection, retaining upstream QC and the sampled barcode roster. Requested and effective ranks are reported separately."]
        paragraphs += _table(["Branch", "Selection", "Requested K", "Effective K", "Status"],
            [(source, attempt.get("mode"), attempt.get("requested_rank"), attempt.get("final_rank"), attempt.get("status"))
             for source, attempts in rank_selections.items() for attempt in attempts])
    paragraphs += ["", "## Marker identification and visualization", "",
        "MarkerFinder (altanalyze3) runs through the standard native marker-heatmap hook. Recorded parameters distinguish the marker-ranking feature count, network export settings and display-cell sampling from analysis-cell selection. Interactive UMAP, expression, violin, dot, combined and heatmap views read saved results; a selected view does not itself establish a new statistical analysis. Per-modality and per-layer marker status is retained in the manifest.",
        ]
    paragraphs += _table(["Marker display / export setting", "Recorded value"],
        [(label, qc.get(key)) for key, label in (("marker_cells_per_cluster", "Display cells per cluster"),
            ("marker_render_heatmap", "Render heatmap image"), ("marker_write_svg", "Write SVG"),
            ("marker_heatmap_dpi", "Image resolution (dpi)")) if key in qc])
    paragraphs += ["", "## Molecular modality prediction and regulatory provenance", ""]
    if m["model_versions"]:
        paragraphs += ["Model outputs are predictions, not additional measured assays. For each recorded modality, the manifest preserves the exact model/version identifier, available artifact/code checksums, input-gene coverage, expression scale, output dimensions, group aggregation and package versions. No model is loaded or retrained to generate this report. Missing fields are not inferred from a model name."]
        for modality, data in sorted(m["model_versions"].items()):
            paragraphs += [f"### {modality}", "", model_methods(modality, data)]
        if any(k in m["model_versions"] for k in ("grn", "grn_tf")):
            paragraphs += ["Regulatory predictions use the recorded RNA-to-GRN bundle and scale/aggregation conventions. ChromLinker [ChromLinker] is cited as regulatory-reference methodology where applicable, not as evidence that this RNA upload underwent ATAC processing, chromatin model fitting or ChromLinker execution. Exact reference applicability must be established from the bundle provenance. RNA-to-modality models and scALABLE-specific integration components without a published citation are unpublished (manuscript in preparation)."]
    else:
        paragraphs += ["No model-version provenance was retained for this job. Selected imputation settings alone do not establish that inference completed."]
    paragraphs += ["", "## Cell communication", ""]
    comm = m["fastcomm_analysis"]
    paragraphs += [f"Recorded fastComm (altanalyze3) status: {comm.get('status', 'not recorded')}; enabled: {comm.get('enabled', 'not recorded')}. {comm.get('message', '')}"]
    if comm:
        paragraphs += ["The native fastComm (altanalyze3) workflow uses ligand–receptor resources from CellChatDB [CellChatDB]; it does not execute the CellChat inference algorithm. fastComm (altanalyze3) is unpublished (manuscript in preparation). Resolved resource checksums and execution parameters are retained when available. Optional per-sample scoring has its own success/failure status; global results do not imply valid per-sample results."]
    paragraphs += ["", "## Differential comparisons, enrichment and networks", ""]
    paragraphs += ["The scALABLE interface caps differential selection within sample × cell state at 500 cells, unless a smaller cap is recorded. The selection method, score source and tie-breaking rule are reported from each comparison's execution record when available; the current policy is not retroactively attributed to historical analyses. The interface defaults to pseudobulk when more than five selected samples are included; otherwise it defaults to cell-level testing. Complete sample/group membership, cell-state field and cap are retained for each run."]
    paragraphs += ["Differential analysis uses cellHarmony-differential (altanalyze3), with modality-specific thresholds recorded below."]
    if m["differential_runs"]:
        paragraphs += ["Each recorded comparison is listed independently below, including its status, modality, sample observation field, groups, population layer and cell/pseudobulk selection. Only completed runs constitute completed differential analyses. The recorded cell-selection cap also applies before pseudobulk aggregation where recorded. RNA pseudobulk and imputed-prediction aggregation are distinct operations. Recorded native test parameters specify the statistical test, minimum group sizes, fold threshold, significance threshold and raw-p/FDR selection. GO-Elite (altanalyze3) [GO-Elite] is cited for native pathway/ontology enrichment when its execution is recorded, with recorded background/settings retained in the reproducibility record."]
        for key, run in sorted(m["differential_runs"].items()):
            effective = (m["analysis_records"].get("main", {}).get("differential:" + key) or {}).get("effective") or {}
            config = dict(run.get("config") or {}, **(effective.get("config") or {}))
            params = effective.get("statistical_parameters") or {}
            selection = effective.get("cell_selection") or {}
            paragraphs += [f"### {config.get('modality', 'Unspecified modality')} comparison", "",
                f"Comparison {key}; status: **{run.get('status', 'not recorded')}**. "
                f"The comparison groups were {_display(config.get('group1_samples'))} versus {_display(config.get('group2_samples'))}."]
            settings = [("Cell-state annotation", config.get("population_col")),
                ("Sample observation field", config.get("sample_field")),
                ("Analysis unit", config.get("comparison_type")),
                ("Effective test", effective.get("actual_population_test")),
                ("Maximum cells per sample and cell state", config.get("max_cells_per_state_sample")),
                ("Cell-selection method", selection.get("method")),
                ("Selection score", selection.get("metric")),
                ("Equal-score tie-breaking", selection.get("tie_breaker")),
                ("Cells available / selected", f"{_display(selection.get('available'))} / {_display(selection.get('selected'))}" if selection else None),
                ("Minimum cells per group", params.get("min_cells_per_group")),
                ("Fold-change threshold", params.get("fc_thresh")),
                ("Significance threshold", params.get("alpha")),
                ("P-value selection", "raw p value" if params.get("use_rawp") else "Benjamini–Hochberg FDR" if "use_rawp" in params else None)]
            paragraphs += _table(["Parameter", "Recorded setting"], settings)
            if selection.get("method") in {"highest_alignment_score", "highest_icgs3_svm_score"}:
                score_description = ("cellHarmony (altanalyze3) reference alignment score" if selection["method"] == "highest_alignment_score"
                                     else "ICGS3 (altanalyze3) winning SVM assignment score")
                paragraphs += [f"Within each selected sample and cell state, cells were ranked by their saved {score_description} in descending order. Equal scores were resolved by ascending cell barcode. Groups below the cap retained all cells. The selected cells were ordered by barcode before testing or pseudobulk aggregation; their identifiers and scores are retained in selected_cells.tsv in the differential archive. Cell-frequency estimates use all cells in the selected groups."]
        paragraphs += ["Cell-level native expression tests use the recorded Scanpy-compatible ranking method (normally Wilcoxon). Pseudobulk dispatch uses the native limma-like moderated t-test even when the configured ranking-method string is Wilcoxon; the actual dispatch is recorded separately for new comparisons. Benjamini–Hochberg correction is applied to eligible features when FDR is selected. RNA independent filtering and the full-family correction for supported imputed outputs follow the recorded native code; raw p values remain distinct from FDR. Cell-communication comparisons use per-sample scores: replicate-group testing requires at least two samples in each arm. With fewer replicates, scores are descriptive effects, not inferential p values. The applicable correction family and test/effect status must be read from each comparison's output and parameter evidence."]
    else:
        paragraphs += ["No differential comparison was recorded. Availability of the Differential tab does not imply that a test or GO-Elite enrichment was performed."]
    paragraphs += ["", "## Software versions and reproducibility", "",
        "scALABLE (altanalyze3), cellHarmony-lite (altanalyze3), ICGS3 (altanalyze3) and other scALABLE-specific unpublished components are manuscripts in preparation; published predecessors are cited separately. Execution-time software versions and code hashes are included below. Model-stage versions for historical jobs apply to those stages only. Exact replay also requires the original inputs, reference and model resource bytes; the ZIP contains no expression matrices or model files.",
        *_software_methods(m), "```json", _json(m["software_provenance"]).strip(), "```"]
    paragraphs += m["provenance_gaps"]
    paragraphs += ["", "## References", ""]
    for key, citation in CITATIONS.items():
        paragraphs.append(f"- [{key}] {citation['authors']} ({citation['year']}). {citation['title']}. {citation['journal']}. [doi:{citation['doi']}](https://doi.org/{citation['doi']})")
    paragraphs += ["", "References are a methodological catalog. Inclusion does not mean that every cited tool ran; the stage descriptions and execution evidence define applicability.", ""]
    return "\n".join(paragraphs)


def build_notebook(manifest):
    """An honest pseudo-notebook: no fabricated execution counts or outputs."""
    cells = []
    def markdown(text):
        cells.append({"cell_type": "markdown", "metadata": {}, "source": text.splitlines(True)})
    def code(text):
        cells.append({"cell_type": "code", "metadata": {}, "source": text.splitlines(True),
                      "execution_count": None, "outputs": []})
    markdown("# scALABLE recorded analysis\n\nAnnotated pseudo-notebook generated without an LLM. No cells below have been executed. The original result is documented in methods.md and run_manifest.json. Replay is disabled by default; restoring matching code, references and model bytes is a prerequisite.\n")
    code("import json\nfrom pathlib import Path\nrecord = json.loads(Path('run_manifest.json').read_text())\nprint(record['job_id'], record['status'], record['analysis_mode'])\n")
    markdown("## Recorded inputs and settings\n\nInspect all libraries, selected species, reference, QC thresholds, analysis layers and model identities. Input files and trained resources are not inside this ZIP.\n")
    code("for item in record['inputs']:\n    print(item)\nprint(json.dumps(record['qc'], indent=2))\nprint(json.dumps(record['model_versions'], indent=2))\n")
    markdown("## Execution environment and effective parameters\n\nRestore the recorded environment; missing historical details cannot be verified. Do not substitute currently installed versions for the original environment. The native parameter lines and full ICGS3 config snapshots distinguish the actual method from an approximation.\n")
    code("print(json.dumps(record['software_provenance'], indent=2))\nfor entry in record['parameter_evidence']:\n    print(entry['text'])\n")
    markdown("## Native pipeline replay template\n\nFill in local paths only after checking recorded software, application code hashes, reference SHA-256 identities and every model artifact. This calls the same supervised/unsupervised/Both workflow hooks rather than replacing analysis algorithms. The guard deliberately prevents an unattended replay with unverified resources. For historical incomplete records, full equivalence requires additional original provenance.\n")
    code("RUN_REPLAY = False\nRESOURCES_AND_METHOD_VERIFIED = False\nINPUT_PATHS = []  # original files, in recorded order\nREFERENCE_REGISTRY = Path('reference_config.json')\nREPLAY_ROOT = Path('scalable_replay')\n\nif RUN_REPLAY:\n    assert RESOURCES_AND_METHOD_VERIFIED, 'Verify recorded code, environment, reference and complete model identities first.'\n    assert len(INPUT_PATHS) == len(record['inputs'])\n    assert all(Path(p).is_file() for p in INPUT_PATHS)\n    import shutil\n    from altanalyze3.components.cellHarmony.flask.job_manager import JobStore\n    from altanalyze3.components.cellHarmony.webapp.analysis_workflow import run_analysis_workflow\n    from altanalyze3.components.cellHarmony.scalable_discover.pipeline import run_discover_pipeline\n    store = JobStore(REPLAY_ROOT)\n    files = [dict(filename=f'{i}_{Path(p).name}', sample_name=src['sample_name'])\n             for i, (p, src) in enumerate(zip(INPUT_PATHS, record['inputs']))]\n    job = store.create_job(record['species'], record['reference_id'], None, files)\n    for path, item in zip(INPUT_PATHS, files):\n        shutil.copyfile(path, store.uploads_dir(job['job_id']) / item['filename'])\n    store.update_job(job['job_id'], qc=record['qc'])\n    if record['reference_id'] == 'icgs3':\n        run_discover_pipeline(job['job_id'], store)\n    else:\n        run_analysis_workflow(job['job_id'], store, REFERENCE_REGISTRY)\n")
    markdown("## Differential comparisons\n\nEach run retains its own modality and supervised/unsupervised population field. Failed and in-progress comparisons are not silently promoted to completed results. Below are the saved configs; replay uses the standard cellHarmony differential hook after native pipeline replay.\n")
    code("for run_id, run in record['differential_runs'].items():\n    print(run_id, run.get('status'), json.dumps(run.get('config', {}), indent=2))\n\nif RUN_REPLAY:\n    from altanalyze3.components.cellHarmony.flask.pipeline import run_cellharmony_differential\n    for original in record['differential_runs'].values():\n        if original.get('status') == 'completed':\n            store.update_job(job['job_id'], differential={'config': original['config'], 'status': 'processing'})\n            run_cellharmony_differential(job['job_id'], store)\n")
    markdown("## Interpretation and citations\n\nConsult methods.md and references.json. Imputed modalities are predictions; retained GRN/ChromLinker provenance describes reference construction only where established. No cell-communication result should be interpreted when its recorded stage failed.\n")
    return {"nbformat": 4, "nbformat_minor": 5, "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"}}, "cells": [dict(cell, id=f"record-{i}") for i, cell in enumerate(cells)]}


def build_archive(meta, job_dir):
    """Return a seeked spool, with bounded copying/compression and stable ZIP metadata."""
    manifest, logs = build_manifest(meta, job_dir)
    spool = tempfile.SpooledTemporaryFile(max_size=8 * 1024 * 1024, mode="w+b")
    try:
        with zipfile.ZipFile(spool, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=3) as archive:
            def add(name, data=None, source=None):
                info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o644 << 16
                with archive.open(info, "w") as output:
                    if source is not None:
                        with source.open("rb") as handle:
                            shutil.copyfileobj(handle, output, length=64 * 1024)
                    else:
                        output.write(data.encode("utf-8"))
            add("methods.md", render_methods(manifest))
            add("analysis.ipynb", _json(build_notebook(manifest)))
            add("run_manifest.json", _json(manifest))
            add("model_provenance.json", _json({"models": manifest["model_versions"],
                "software_provenance": manifest["software_provenance"],
                "layer_communication": manifest["layers"]}))
            add("references.json", _json(CITATIONS))
            for path in logs:
                if path.resolve().is_relative_to(Path(job_dir).resolve()):
                    add(path.relative_to(job_dir).as_posix(), source=path)
        spool.seek(0)
        return spool
    except BaseException:
        spool.close()
        raise
