"""scALABLE-discover: scALABLE-web with ICGS3 clustering in place of reference alignment.

This module builds no viewer of its own. `webapp.app.create_app()` builds the whole
scALABLE application: upload, QC, every Explore plot and PDF renderer, Chat, its template
and its front end. Every route reads a job through `app.state.job_store` and
`app.state.job_runner`, so swapping the runner for `DiscoverJobRunner` makes the same
application run the ICGS3 pipeline (`scalable_discover/pipeline.py`). The scalable_viewer
uses the same pattern (`visualization/scalable_viewer/scalable_app.py`).

What this module changes:
  * the runner: ICGS3 instead of cellHarmony alignment, approximate UMAP and imputation,
    with the Run tab's percentage and stage name following ICGS3's major steps
  * the species menu: human or mouse, with no reference atlas
  * upload, QC and configure validate species only and refuse imputation
  * every differential route and the approximate-UMAP tool are removed
  * cell-state layers: every view opens on the predicted cell states; a viewer may switch
    to the C1..Cn clusters (cookie `discover_layer`), read by LayeredJobStore
  * a GO-Elite BioMarkers plot per cell state, drawn from ICGS3's own enrichment table
  * the served page: renamed, Differential tab removed, reference and alignment controls
    hidden, LungMAP page colour, plus static/discover.js. The shared template is not edited;
    the body is rewritten on the way out and every rewrite must match, or the page answers 500.

Launch locally, from the repository root (the directory that holds altanalyze3/):
    python3.11 -m uvicorn altanalyze3.components.cellHarmony.scalable_discover.app:app \\
        --host 127.0.0.1 --port 8010
Online: docker-compose.discover.yml in this directory.
"""
from __future__ import annotations

import contextvars
import json
import os
from importlib import import_module
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from fastapi import FastAPI, File, Form, HTTPException, Query, Request, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
from altanalyze3.components.goelite.structures import compute_z_score

from .pipeline import DISCOVER_REFERENCE_ID, SPECIES_TO_ICGS, discover_registry
from .tasks import DiscoverJobRunner

# The webapp package re-exports a FastAPI instance named `app`; import the module.
W = import_module("altanalyze3.components.cellHarmony.webapp.app")

HERE = Path(__file__).resolve().parent
APP_NAME = "scALABLE-discover"
DEFAULT_JOB_STORAGE = HERE / "jobs"
LAYER_COOKIE = "discover_layer"
_LAYER = contextvars.ContextVar("scalable_discover_layer", default="")
_GITHUB_BLOB = "https://github.com/SalomonisLab/altanalyze3/blob/master/altanalyze3/components/cellHarmony"
EMBEDDING_LABEL = "UMAP"
# scALABLE's GO Terms view highlights a term at FDR <= 0.05 and z > 2
# (webapp/app.py _build_differential_go_payload, `is_positive_sig`); the BioMarkers plot uses
# the same rule. BioMarkers sets carry no ontology, so GO-Elite's parent/child pruning does
# not apply and every term passing the rule is highlighted.
GOELITE_MAX_FDR = 0.05
GOELITE_MIN_Z = 2.0

# (anchor in webapp/templates/index.html, replacement, required count). A count mismatch
# raises, so a template edit that moves an anchor stops the page instead of silently
# serving scALABLE-web's controls under the discover name.
INDEX_REWRITES: Tuple[Tuple[str, str, int], ...] = (
    ("<h1>scALABLE</h1>", f"<h1>{APP_NAME}</h1>", 1),
    ("Align single-cell data to a reference, explore results and perform comprehensive comparisons between groups",
     "Cluster single-cell data without a reference using ICGS3, then explore the clusters and ask questions", 1),
    ('          <button class="workspace-tab-btn" type="button" data-tab="differential">Differential</button>\n', "", 1),
    (f"{_GITHUB_BLOB}/webapp/HOW_TO_USE.md", f"{_GITHUB_BLOB}/scalable_discover/HOW_TO_USE.md", 1),
    (f"{_GITHUB_BLOB}/webapp/README.md", f"{_GITHUB_BLOB}/scalable_discover/README.md", 1),
    ('<label class="field">\n                  <span>Reference</span>',
     '<label class="field hidden">\n                  <span>Reference</span>', 1),
    ("<h2>2. QC and alignment</h2>", "<h2>2. QC and ICGS3 clustering</h2>", 1),
    ('<label class="field">\n                    <span>Minimum cosine similarity score</span>',
     '<label class="field hidden">\n                    <span>Minimum cosine similarity score</span>', 1),
    ("Live counts parsed from alignment log.", "Live counts parsed from the QC log.", 1),
    # The ICGS3 workflow figure fills this panel and carries its own title.
    ("<h2>Reference Preview</h2>", "", 1),
    ('<option value="relative">UMAP broad</option>', "", 2),
    ('UMAP cell types</option>', "UMAP cell states</option>", 2),
)


def rewrite_index(body: str, root_path: str, script_version: int) -> str:
    """Apply INDEX_REWRITES to the rendered scALABLE page and add discover.css and .js."""
    for old, new, expected in INDEX_REWRITES:
        found = body.count(old)
        if found != expected:
            raise RuntimeError(f"scALABLE-discover page rewrite: expected {expected} of {old[:70]!r}, "
                               f"found {found}. webapp/templates/index.html changed; update INDEX_REWRITES.")
        body = body.replace(old, new)
    for tag in ("</head>", "</body>"):
        if body.count(tag) != 1:
            raise RuntimeError(f"scALABLE-discover page rewrite: expected one {tag}.")
    body = body.replace("</head>", f'<link rel="stylesheet" href="{root_path}/discover-static/discover.css'
                                   f'?v={script_version}"></head>')
    tag = (f'<script>window.__SCALABLE_DISCOVER__ = true;</script>'
           f'<script src="{root_path}/discover-static/discover.js?v={script_version}"></script></body>')
    return body.replace("</body>", tag)


# ------------------------------------------------------------------ cell-state layers

def apply_layer(meta: Dict, layer: str) -> Dict:
    """The job as one cell-state layer sees it.

    The default layer's values sit at the top level of the job, where every scALABLE reader
    looks. Another layer swaps in its own cluster key, MarkerFinder set and fastComm run.
    The stored job is never changed: JobStore.update_job reads the raw file.
    """
    layers = meta.get("cell_state_layers") or {}
    default = str(layers.get("default") or meta.get("cluster_key") or "")
    entry = next((e for e in layers.get("layers") or [] if e.get("key") == layer), None)
    if not layer or layer == default or entry is None:
        return dict(meta, active_cell_state_layer=default) if layers else meta
    meta = dict(meta)
    meta["cluster_key"] = entry["key"]
    if entry.get("marker_analysis"):
        meta["marker_analysis"] = entry["marker_analysis"]
        meta["marker_analysis_by_modality"] = dict(meta.get("marker_analysis_by_modality") or {},
                                                   rna=entry["marker_analysis"])
    if entry.get("fastcomm_analysis") is not None:
        meta["fastcomm_analysis"] = entry["fastcomm_analysis"]
    meta["active_cell_state_layer"] = entry["key"]
    return meta


class LayeredJobStore(JobStore):
    """The upload job store, read through the layer the current viewer chose."""

    def get_job(self, job_id: str) -> Dict:
        return apply_layer(super().get_job(job_id), _LAYER.get())


# ------------------------------------------------------------- GO-Elite BioMarkers plot

def _cluster_for_state(meta: Dict, state: str) -> str:
    """A cell state of either layer -> its ICGS3 cluster id."""
    names = (meta.get("cell_state_layers") or {}).get("names") or {}
    reverse = {str(v): str(k) for k, v in names.items()}
    return reverse.get(str(state), str(state))


def _state_label(meta: Dict, cluster: str) -> str:
    """An ICGS3 cluster id as the active layer names it."""
    layers = meta.get("cell_state_layers") or {}
    if meta.get("active_cell_state_layer") == "cluster":
        return cluster
    return str((layers.get("names") or {}).get(cluster, cluster))


def goelite_states(meta: Dict) -> List[str]:
    """The active layer's cell states that hold at least one BioMarkers term."""
    goelite = (meta.get("icgs3_analysis") or {}).get("goelite") or {}
    path = Path(str(goelite.get("enrichment_tsv") or ""))
    if not path.is_file():
        return []
    clusters = set(pd.read_csv(path, sep="\t", usecols=["cluster"])["cluster"].astype(str))
    order = (meta.get("icgs3_analysis") or {}).get("clusters") or sorted(clusters)
    return [_state_label(meta, c) for c in order if c in clusters]


def build_goelite_payload(meta: Dict, state: str) -> Dict:
    """One cell state's BioMarkers terms, in the shape scALABLE's GO Terms view draws.

    p, FDR and overlap are ICGS3's own (ICGS.biomarker_enrichment: hypergeometric p, BH FDR
    over every cluster and term). ICGS3 stores no z-score, so z is GO-Elite's own
    `compute_z_score` (goelite/structures.py) on ICGS3's overlap, query size, term size and
    the gene universe ICGS3 tested against.
    """
    goelite = (meta.get("icgs3_analysis") or {}).get("goelite") or {}
    path = Path(str(goelite.get("enrichment_tsv") or ""))
    background = int(goelite.get("background_size") or 0)
    cluster = _cluster_for_state(meta, state)
    label = _state_label(meta, cluster)
    payload = {"population": label, "cluster": cluster, "terms": [], "labels": [],
               "statistics": {"p_value": "ICGS3 hypergeometric p", "fdr": "ICGS3 BH FDR, all clusters x terms",
                              "z_score": "GO-Elite compute_z_score on ICGS3 overlap counts",
                              "background_size": background,
                              "highlight": f"FDR <= {GOELITE_MAX_FDR} and z > {GOELITE_MIN_Z}"}}
    if not path.is_file() or background <= 0:
        payload["message"] = "This job has no GO-Elite BioMarkers table."
        return payload
    frame = pd.read_csv(path, sep="\t")
    frame = frame.loc[frame["cluster"].astype(str) == cluster]
    if frame.empty:
        payload["message"] = f"No BioMarkers term overlaps the markers of {label}."
        return payload
    # ICGS3 labels each cluster with its first row after sorting by FDR, p and name.
    frame = frame.sort_values(["fdr", "p_value", "term_name"]).reset_index(drop=True)
    terms = []
    for index, row in frame.iterrows():
        z_score = float(compute_z_score(int(row["overlap"]), int(row["query_size"]),
                                        int(row["term_size"]), background))
        fdr = float(row["fdr"])
        genes = [g for g in str(row.get("overlap_genes") or "").split(",") if g]
        significant = fdr <= GOELITE_MAX_FDR and z_score > GOELITE_MIN_Z
        terms.append({
            "term_id": "", "term_name": str(row["term_name"]), "direction": "up",
            "fdr": fdr, "p_value": float(row["p_value"]), "z_score": z_score,
            "fdr_plot": float(min(max(fdr, 1e-300), 1.0)), "score": float(-np.log10(min(max(fdr, 1e-300), 1.0))),
            "overlap": int(row["overlap"]), "query_size": int(row["query_size"]), "term_size": int(row["term_size"]),
            "selected": significant, "is_positive_sig": significant, "is_selected_positive_sig": significant,
            "is_prediction": index == 0, "overlap_genes": genes, "selected_gene": genes[0] if genes else None,
        })
    payload["terms"] = terms
    labelled = [t for t in terms if t["is_selected_positive_sig"]][:4]
    if terms[0] not in labelled:
        labelled = [terms[0]] + labelled[:3]
    payload["labels"] = [{"term_name": t["term_name"], "z_score": t["z_score"], "fdr_plot": t["fdr_plot"],
                          "selected_gene": t["selected_gene"], "overlap_genes": t["overlap_genes"],
                          "label_color": "#1f19c7" if t["is_prediction"] else "#111827",
                          "label_rank": i, "label_role": "prediction" if t["is_prediction"] else "top"}
                         for i, t in enumerate(labelled)]
    return payload


# ------------------------------------------------------------------- wrapped builders

def _relabel_embedding(value):
    """scALABLE-web names its primary embedding "cellHarmony UMAP"; here ICGS3 computed it."""
    if isinstance(value, str):
        return value.replace("cellHarmony UMAP", EMBEDDING_LABEL)
    if isinstance(value, list):
        return [_relabel_embedding(v) for v in value]
    if isinstance(value, dict):
        return {k: (_relabel_embedding(v) if k in {"label", "coords_label", "x_label", "y_label"}
                    or isinstance(v, (dict, list)) and k not in {"query", "reference"} else v)
                for k, v in value.items()}
    return value


def _install_shared_wrappers() -> None:
    """Wrap webapp builders the routes call by module name; once per process.

    As scalable_app.py wraps the same module. Only labels and the live stage message change;
    no value is touched.
    """
    if getattr(W, "_scalable_discover_wrappers_installed", False):
        return
    for name in ("_axis_field_options", "_umap_coordinate_options", "_build_umap_payload"):
        original = getattr(W, name)

        def wrapper(*args, __original=original, **kwargs):
            return _relabel_embedding(__original(*args, **kwargs))

        wrapper.__name__ = original.__name__
        setattr(W, name, wrapper)

    original_examples = W._chat_examples

    def chat_examples(app, meta):
        # The example header reads "Try (<reference label>):". A discover job has no reference.
        payload = original_examples(app, meta)
        if str(meta.get("reference") or "") == DISCOVER_REFERENCE_ID:
            layers = meta.get("cell_state_layers") or {}
            active = meta.get("active_cell_state_layer") or layers.get("default")
            entry = next((e for e in layers.get("layers") or [] if e.get("key") == active), None)
            payload["reference"] = (entry or {}).get("label") or "cell states"
        return payload

    W._chat_examples = chat_examples

    original_live = W._derive_live_pipeline_message

    def live_message(status, log_lines, fallback):
        # The runner keeps the job message on the current ICGS3 stage; scALABLE-web's own
        # log scan would report its last alignment-era marker instead.
        if str(status or "").strip().lower() == "processing" and fallback:
            return str(fallback)
        return original_live(status, log_lines, fallback)

    W._derive_live_pipeline_message = live_message
    W._scalable_discover_wrappers_installed = True


# -------------------------------------------------------------------------- routes

def _pop_route(app: FastAPI, path: str, method: str):
    for index, route in enumerate(list(app.router.routes)):
        if getattr(route, "path", None) == path and method in (getattr(route, "methods", None) or set()):
            return app.router.routes.pop(index)
    raise RuntimeError(f"scALABLE route {method} {path} not found; webapp/app.py changed.")


def _remove_routes(app: FastAPI, predicate) -> List[str]:
    removed = []
    for route in list(app.router.routes):
        path = str(getattr(route, "path", ""))
        if predicate(path):
            app.router.routes.remove(route)
            removed.append(path)
    return removed


def _validate_species(species: str, reference: str) -> str:
    value = str(species or "").strip().lower()
    if value not in SPECIES_TO_ICGS:
        raise HTTPException(status_code=400, detail=f"Species must be human or mouse, not '{species}'.")
    if str(reference or "").strip() != DISCOVER_REFERENCE_ID:
        raise HTTPException(status_code=400, detail=f"scALABLE-discover uses no reference atlas; "
                                                    f"reference must be '{DISCOVER_REFERENCE_ID}'.")
    return value


def create_discover_app(overrides: Optional[Dict] = None) -> FastAPI:
    config = {
        "APP_TITLE": APP_NAME,
        "JOB_STORAGE": os.getenv("SCALABLE_DISCOVER_JOB_STORAGE", str(DEFAULT_JOB_STORAGE)),
        "ROOT_PATH": os.getenv("SCALABLE_DISCOVER_ROOT_PATH", ""),
    }
    config.update(overrides or {})
    app = W.create_app(config)
    cfg = app.state.config
    previous = app.state.job_runner
    app.state.job_store = LayeredJobStore(Path(cfg["JOB_STORAGE"]))
    app.state.job_runner = DiscoverJobRunner(
        app.state.job_store,
        Path(cfg["REFERENCE_REGISTRY"]),
        max_workers=cfg["JOB_WORKERS"],
        export_approx_pdfs=False,
        h5ad_compression=cfg.get("H5AD_COMPRESSION", "lzf"),
        isolate_jobs=cfg.get("ISOLATE_JOBS", True),
        worker_memory_limit_gib=cfg.get("WORKER_MEMORY_LIMIT_GIB", 15),
        total_memory_limit_gib=cfg.get("TOTAL_MEMORY_LIMIT_GIB", 27),
    )
    previous.executor.shutdown(wait=False)
    app.state.discover = True
    _install_shared_wrappers()

    @app.middleware("http")
    async def cell_state_layer(request: Request, call_next):
        # Each viewer's layer choice rides in a cookie, so every scALABLE route reads the job
        # through it without a new query parameter.
        token = _LAYER.set(str(request.cookies.get(LAYER_COOKIE) or ""))
        try:
            return await call_next(request)
        finally:
            _LAYER.reset(token)

    # No comparison runs here, and the approximate-UMAP tool projects onto a reference.
    app.state.removed_routes = _remove_routes(
        app, lambda path: "/differential" in path or path in {"/api/tools/approximate-umap",
                                                              "/api/meta/reference-preview"})
    original_create_job = _pop_route(app, "/api/jobs", "POST").endpoint
    for path, method in (("/", "GET"), ("/api/meta/species", "GET"),
                         ("/api/jobs/{job_id}/qc", "POST"), ("/api/jobs/{job_id}/configure", "POST")):
        _pop_route(app, path, method)

    templates = Jinja2Templates(directory=str(cfg["TEMPLATE_DIR"]))
    static_dir = Path(str(cfg["STATIC_DIR"]))
    discover_static = HERE / "static"
    app.mount("/discover-static", StaticFiles(directory=str(discover_static)), name="discover-static")

    @app.get("/", response_class=HTMLResponse)
    async def index(request: Request):
        css_version = int((static_dir / "styles.css").stat().st_mtime) if (static_dir / "styles.css").exists() else 0
        js_version = max((p.stat().st_mtime_ns for p in (static_dir / "app.js", static_dir / "integrated.js")
                          if p.exists()), default=0)
        response = templates.TemplateResponse(request, str(cfg["INDEX_TEMPLATE"]), {
            "request": request,
            "registry_json": json.dumps(discover_registry()),
            "app_title": cfg["APP_TITLE"],
            "app_root_path": cfg["ROOT_PATH"],
            "styles_version": css_version,
            "app_js_version": js_version,
        })
        body = response.body.decode("utf-8")
        version = max(int(p.stat().st_mtime_ns) for p in discover_static.glob("discover.*"))
        return HTMLResponse(rewrite_index(body, cfg["ROOT_PATH"], version))

    @app.get("/api/meta/species")
    async def meta_species():
        return JSONResponse(discover_registry())

    @app.post("/api/jobs")
    async def create_job(
        species: str = Form(...),
        reference: str = Form(DISCOVER_REFERENCE_ID),
        ambient_option: str | None = Form(None),
        soupx_option: str | None = Form(None),
        sample_names: List[str] = Form(...),
        files: List[UploadFile] = File(...),
    ):
        species = _validate_species(species, reference)
        return await original_create_job(species=species, reference=DISCOVER_REFERENCE_ID,
                                         ambient_option=ambient_option, soupx_option=soupx_option,
                                         sample_names=sample_names, files=files)

    @app.post("/api/jobs/{job_id}/qc")
    async def update_qc(job_id: str, qc: W.QCSettings):
        store, _ = W._job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        requested = [str(v) for v in (qc.impute_modalities or []) if str(v).strip().lower() not in {"", "none"}]
        if requested or str(qc.impute_modality or "none").strip().lower() not in {"", "none"}:
            raise HTTPException(status_code=400, detail="scALABLE-discover runs no modality imputation.")
        values = qc.model_dump()
        values.pop("align_cutoff", None)       # no alignment, so no alignment cutoff is stored
        values.update(impute_modality="none", impute_modalities=[])
        meta = store.update_job(job_id, qc=values, message="QC parameters saved.")
        return JSONResponse({"job_id": job_id, "qc": meta["qc"]})

    @app.post("/api/jobs/{job_id}/configure")
    async def configure_job(job_id: str, settings: W.JobConfigSettings):
        store, _ = W._job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        species = _validate_species(settings.species, settings.reference)
        meta = store.get_job(job_id)
        requested_ambient = (settings.ambient_option if settings.ambient_option is not None
                             else meta.get("ambient_option", meta.get("soupx_option")))
        changed = (str(meta.get("species") or "") != species
                   or str(meta.get("ambient_option", meta.get("soupx_option")) or "") != str(requested_ambient or ""))
        if changed:
            W._invalidate_expression_cache(app, job_id)
            W._invalidate_marker_heatmap_cache(app, job_id)
            W._invalidate_fastcomm_cache(app, job_id)
            W._clear_directory_contents(store.outputs_dir(job_id))
            W._clear_directory_contents(store.logs_dir(job_id))
            meta = store.update_job(
                job_id, species=species, reference=DISCOVER_REFERENCE_ID, ambient_option=requested_ambient,
                status="uploaded", progress=0, message="Species updated. Configure QC and rerun ICGS3.",
                artifacts={}, marker_analysis={}, marker_analysis_by_modality={}, fastcomm_analysis={},
                icgs3_analysis={}, cell_state_layers={}, modality_artifacts={},
                modalities={"default": "rna", "available": [dict(W._DEFAULT_MODALITY_DEFINITIONS["rna"])]},
                differential={},
            )
        return JSONResponse({"job_id": job_id, "status": meta.get("status"), "species": meta.get("species"),
                             "reference": meta.get("reference"),
                             "ambient_option": meta.get("ambient_option", meta.get("soupx_option")),
                             "changed": changed})

    def _job_meta(job_id: str) -> Dict:
        store, _ = W._job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        return store.get_job(job_id)

    @app.get("/api/jobs/{job_id}/biomarkers/goelite/states")
    def goelite_state_list(job_id: str):
        return JSONResponse({"states": goelite_states(_job_meta(job_id))})

    @app.get("/api/jobs/{job_id}/biomarkers/goelite")
    def goelite_plot(job_id: str, population: str = Query(...)):
        return JSONResponse(build_goelite_payload(_job_meta(job_id), population))

    @app.get("/api/jobs/{job_id}/biomarkers/goelite.pdf")
    def goelite_pdf(job_id: str, population: str = Query(...)):
        payload = build_goelite_payload(_job_meta(job_id), population)
        if not payload["terms"]:
            raise HTTPException(status_code=404, detail=payload.get("message") or "No BioMarkers terms.")
        # scALABLE's own GO Terms PDF renderer draws the same plot; its title reads "GO terms:".
        payload = dict(payload, population=f"{payload['population']} (GO-Elite BioMarkers)")
        pdf = W._render_differential_go_pdf(payload)
        name = "".join(ch if ch.isalnum() else "_" for ch in str(population))[:80]
        return StreamingResponse(pdf, media_type="application/pdf",
                                 headers={"Content-Disposition": f'attachment; filename="{job_id}_{name}_goelite_biomarkers.pdf"'})

    return app


app = create_discover_app()
