"""Deterministic methods and notebook export; never load expression or fit a model.

Execution snapshots are captured by workers, not reconstructed from the download
server's environment. Historical runs retain explicit gaps in their provenance.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import shutil
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path


SCHEMA_VERSION = 1
QC_DEFAULTS = {
    "analysis_mode": "supervised", "umap_fit_mode": "landmark", "max_k": None,
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
              "application_method": describe_application(application),
              "species": meta.get("species"), "reference_id": meta.get("reference"),
              "qc": meta.get("qc") or {}, "interface_defaults": dict(QC_DEFAULTS),
              "effective": effective or {}}
    package_root = Path(__file__).resolve().parents[1]
    record["method_files"] = {name: file_identity(package_root / name) for name in (
        "cellHarmony/cellHarmony_lite.py", "cellHarmony/cellHarmony_differential.py",
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


def render_methods(manifest):
    m = manifest
    mode, qc = m["analysis_mode"], m["qc"]
    paragraphs = ["# scALABLE analysis methods", "",
        "This report was generated deterministically from saved execution metadata, model provenance and logs. It contains no LLM-generated methods. Requested settings are distinguished from stage execution evidence; failed or unfinished stages must not be reported as completed experiments.",
        "", "## Analysis identity and input", "",
        f"Job `{m['job_id']}`; status: **{m['status']}**. Analysis: **{mode}**. Species: **{m['species']}**. Reference identifier: `{m['reference_id']}`.",
        f"Recorded start: {m.get('analysis_started_at')}; completion: {m.get('analysis_completed_at')}; worker duration (seconds): {m.get('analysis_duration_seconds')}.",
        "Inputs retain their recorded library/sample names. For an annotated H5AD, the selected observation field defines the sample identities used in each differential comparison.",
        "```json", _json(m["inputs"]).strip(), "```", "", "## Selected reference", ""]
    ref = m.get("reference") or {}
    if ref:
        paragraphs += [f"Selected reference: {ref.get('label', m['reference_id'])}."]
        if ref.get("study_url"):
            paragraphs += [f"Reference citation: [{ref.get('study_citation') or ref['study_url']}]({ref['study_url']})."]
        else:
            paragraphs += ["No published study citation is recorded for this reference; no citation to a different reference version has been substituted."]
    else:
        paragraphs += ["Reference metadata/version at execution was not recorded. The identifier and reference paths in the parameter evidence below are the available record."]
    paragraphs += ["Reference resource SHA-256 identities, when captured, are in `run_manifest.json`.",
        "", "## Quality control, ambient RNA correction and expression scale", "",
        "The requested QC settings below include minimum detected genes per cell, minimum counts per cell, minimum cells per gene, mitochondrial percentage and the reference alignment score threshold. The execution log determines whether those filters were actually applied. Already normalized input can follow a different preprocessing branch.",
        "", "| Setting | Recorded value | Interface default at execution | Difference |",
        "|---|---|---|---|"]
    defaults = m.get("interface_defaults_at_execution") or {}
    for key in sorted(set(qc) | set(defaults)):
        value = _json(qc[key]).strip().replace("\n", " ") if key in qc else "not recorded"
        default = _json(defaults[key]).strip().replace("\n", " ") if key in defaults else "not recorded"
        difference = "yes" if key in qc and key in defaults and qc[key] != defaults[key] else "no" if key in qc and key in defaults else "unknown"
        paragraphs.append(f"| {key} | {value.replace('|', '&#124;')} | {default.replace('|', '&#124;')} | {difference} |")
    paragraphs += ["", "Ambient correction was requested: **" + str(qc.get("ambient_correction", "not recorded")) + "**. The native implementation is AmbientSubtract (unpublished; manuscript in preparation), not canonical SoupX. Actual execution, counts source, selected correction fractions and any skipped normalization/QC are recorded below. The native raw-count workflow uses total-count normalization to 10,000 counts per cell followed by natural-log `log1p`; already log-normalized X is not normalized a second time. Correction of normalized H5AD input uses the recorded counts-layer branch when applicable. These are implementation descriptions, not evidence that every branch ran for this job.",
        "", "```text"]
    paragraphs += [e["text"] for e in m["preprocessing_evidence"]] or ["No preprocessing execution evidence was retained."]
    paragraphs += ["```", "", "## Cell-state inference and embedding", ""]
    if mode in {"supervised", "both"}:
        paragraphs += ["The supervised branch uses cellHarmony-lite reference alignment. Query cells are compared with reference cell-state profiles over the reference marker panel; the native pipeline uses cosine similarity and the requested alignment cutoff. This implementation builds on cellHarmony [cellHarmony], but its centroid alignment and approximate reference placement should not be described as the original full cell-to-cell matching algorithm.",
            "The approximate supervised UMAP places accepted query cells using reference coordinates rather than fitting a new UMAP to all query cells. The recorded placement parameters specify the reference-cell count and jitter; a seed is only claimed if recorded."]
    if mode in {"unsupervised", "both"}:
        paragraphs += ["The unsupervised branch invokes the native `ICGS.run_icgs3` implementation. ICGS3 is unpublished (manuscript in preparation); ICGS2 [ICGS2] is its published methodological precursor, not the version executed here. The workflow includes variable-program discovery, downsampling, NMF, MarkerFinder refinement and SVM assignment. Resolved ICGS3 configuration, including target rank, downsampling, normalization, UMAP features/backend, landmark limits, neighbors and seeds, is retained in the execution snapshots and parameter appendix when available. The accelerated UMAP setting concerns embedding; it is not an instruction to discard cells from the returned dataset."]
    if mode == "both":
        paragraphs += ["Both branches share upstream QC. The combined H5AD retains the union of accepted barcodes with separate supervised and unsupervised annotation/coordinate layers. Cells without a supervised match are represented as Unaligned rather than silently removed from the combined result. Cluster-association plots summarize barcode intersections; they are not independent biological measurements."]
    paragraphs += ["", "## Marker identification and visualization", "",
        "MarkerFinder runs through the standard native marker-heatmap hook. Recorded parameters distinguish the marker-ranking feature count, network export settings and display-cell sampling from analysis-cell selection. Interactive UMAP, expression, violin, dot, combined and heatmap views read saved results; a selected view does not itself establish a new statistical analysis. Per-modality and per-layer marker status is retained in the manifest.",
        "", "## Molecular modality prediction and regulatory provenance", ""]
    if m["model_versions"]:
        paragraphs += ["Model outputs are predictions, not additional measured assays. For each recorded modality, the manifest preserves the exact model/version identifier, available artifact/code checksums, input-gene coverage, expression scale, output dimensions, group aggregation and package versions. No model is loaded or retrained to generate this report. Missing fields are not inferred from a model name."]
        for modality, data in sorted(m["model_versions"].items()):
            paragraphs += [f"### {modality}", "", "```json", _json(data).strip(), "```"]
        if any(k in m["model_versions"] for k in ("grn", "grn_tf")):
            paragraphs += ["Regulatory predictions use the recorded RNA-to-GRN bundle and scale/aggregation conventions. ChromLinker [ChromLinker] is cited as regulatory-reference methodology where applicable, not as evidence that this RNA upload underwent ATAC processing, chromatin model fitting or ChromLinker execution. Exact reference applicability must be established from the bundle provenance. RNA-to-modality models and scALABLE-specific integration components without a published citation are unpublished (manuscript in preparation)."]
    else:
        paragraphs += ["No model-version provenance was retained for this job. Selected imputation settings alone do not establish that inference completed."]
    paragraphs += ["", "## Cell communication", ""]
    comm = m["fastcomm_analysis"]
    paragraphs += [f"Recorded fastComm status: {comm.get('status', 'not recorded')}; enabled: {comm.get('enabled', 'not recorded')}. {comm.get('message', '')}"]
    if comm:
        paragraphs += ["The native fastComm workflow uses ligand–receptor resources from CellChatDB [CellChatDB]; it does not execute the CellChat inference algorithm. fastComm is unpublished (manuscript in preparation). Resolved resource checksums and execution parameters are retained when available. Optional per-sample scoring has its own success/failure status; global results do not imply valid per-sample results.", "```json", _json(comm).strip(), "```"]
    paragraphs += ["", "## Differential comparisons, enrichment and networks", ""]
    paragraphs += ["The scALABLE interface caps differential selection at 500 randomly sampled cells per sample × cell state (seed 0), unless the user chooses a smaller cap. The interface defaults to pseudobulk when more than five selected samples are included; otherwise it defaults to cell-level testing. These interface policies are distinct from each recorded comparison's actual selection and effective test. Complete sample/group membership, cell-state field and cap are retained for each run."]
    if m["differential_runs"]:
        paragraphs += ["Each recorded comparison is listed independently below, including its status, modality, sample observation field, groups, population layer and cell/pseudobulk selection. Only completed runs constitute completed differential analyses. Cell selection is capped within sample × cell state using the recorded random seed; the cap also applies before pseudobulk aggregation where recorded. RNA pseudobulk and imputed-prediction aggregation are distinct operations. Recorded native test parameters specify the statistical test, minimum group sizes, fold threshold, significance threshold and raw-p/FDR selection. GO-Elite [GO-Elite] is cited for native pathway/ontology enrichment when its execution is recorded, with background/settings in the parameter appendix."]
        for key, run in sorted(m["differential_runs"].items()):
            details = {k: v for k, v in run.items() if k != "artifacts"}
            paragraphs += [f"### Comparison {key}", "", "```json", _json(details).strip(), "```"]
        paragraphs += ["Cell-level native expression tests use the recorded Scanpy-compatible ranking method (normally Wilcoxon). Pseudobulk dispatch uses the native limma-like moderated t-test even when the configured ranking-method string is Wilcoxon; the actual dispatch is recorded separately for new comparisons. Benjamini–Hochberg correction is applied to eligible features when FDR is selected. RNA independent filtering and the full-family correction for supported imputed outputs follow the recorded native code; raw p values remain distinct from FDR. Cell-communication comparisons use per-sample scores: replicate-group testing requires at least two samples in each arm. With fewer replicates, scores are descriptive effects, not inferential p values. The applicable correction family and test/effect status must be read from each comparison's output and parameter evidence."]
    else:
        paragraphs += ["No differential comparison was recorded. Availability of the Differential tab does not imply that a test or GO-Elite enrichment was performed."]
    paragraphs += ["", "## Software versions and reproducibility", "",
        "scALABLE, cellHarmony-lite, ICGS3 and other scALABLE-specific unpublished components are manuscripts in preparation; published predecessors are cited separately. Execution-time software versions and code hashes are included below. Model-stage versions for historical jobs apply to those stages only. Exact replay also requires the original inputs, reference and model resource bytes; the ZIP contains no expression matrices or model files.",
        "```json", _json(m["software_provenance"]).strip(), "```"]
    paragraphs += m["provenance_gaps"]
    paragraphs += ["", "## Native execution parameter appendix", "",
        "The following lines are the recorded effective parameters, not values guessed from current defaults.", "", "```text"]
    paragraphs += [e["text"] for e in m["parameter_evidence"]] or ["No parameter lines were retained."]
    paragraphs += ["```", "", "## References", ""]
    for key, citation in CITATIONS.items():
        paragraphs.append(f"- [{key}] {citation['authors']} ({citation['year']}). {citation['title']}. {citation['journal']}. https://doi.org/{citation['doi']}")
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
