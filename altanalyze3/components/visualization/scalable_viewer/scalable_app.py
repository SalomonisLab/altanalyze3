"""The scalable_viewer server: scALABLE's own app, served over precomputed bundles.

Design in one sentence: this module does not build a viewer, it *hands bundles to the
viewer that already exists*.

`altanalyze3.components.cellHarmony.webapp.app.create_app()` builds the whole scALABLE
FastAPI application - every route, every payload builder, every matplotlib PDF renderer,
its Jinja template, its stylesheet and its 5,350-line front end. Each of its endpoints
starts with `store.job_exists(job_id)` / `store.get_job(job_id)` (app.py:4499 onward) and
then calls a builder that takes only `(app, meta, ...)`. Swapping the job store for
`bundle_meta.BundleJobStore` therefore makes the entire application serve a precomputed
bundle, with no plot, colour ramp or payload rewritten here.

What this module adds:
  * the bundle-backed job store and cache seeding (bundle_meta.py)
  * `/api/catalog` and a dataset selector, so one process serves many bundles
  * a POST override that switches between the bundle's precomputed DEG contrasts
    instead of launching a differential run
  * a dot plot (`/api/jobs/{id}/dotplot`), which scALABLE does not have, defaulting to
    the top marker of every cell state
  * `/fast/*`, the original binary/memmap endpoints, kept so nothing regresses
  * a bootstrap script that auto-loads the first dataset - there is no Load button
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
from fastapi import HTTPException, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from importlib import import_module

from altanalyze3.components.cellHarmony.webapp.config import BASE_DIR as WEBAPP_DIR

# The webapp package re-exports a FastAPI instance named `app`; import the module.
W = import_module("altanalyze3.components.cellHarmony.webapp.app")

from . import bundle_meta
from . import data_api as da
from .server import create_app as create_fast_app

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_assets(assets_root: Optional[str], catalog: da.Catalog) -> Dict[str, Dict[str, Any]]:
    """Read one `<prefix>_assets.json` per dataset, written by prepare_assets.py."""
    out: Dict[str, Dict[str, Any]] = {}
    if not assets_root:
        return out
    for entry in catalog.entries:
        manifest = Path(assets_root) / entry["id"] / f"{entry['prefix']}_assets.json"
        if not manifest.is_file():
            continue
        with open(manifest, "r") as fh:
            data = json.load(fh)
        out[entry["id"]] = data
    return out


def create_scalable_app(
    catalog: da.Catalog,
    *,
    state_dir: str,
    assets_root: Optional[str] = None,
):
    assets = _load_assets(assets_root, catalog)

    # scALABLE's own template directory and static directory. Nothing is copied: the
    # viewer serves app.js, styles.css and index.html straight out of the webapp.
    app = W.create_app({
        "JOB_STORAGE": state_dir,
        "TEMPLATE_DIR": str(WEBAPP_DIR / "templates"),
        "STATIC_DIR": str(WEBAPP_DIR / "static"),
        "INDEX_TEMPLATE": "index.html",
        "ROOT_PATH": "",
    })

    store = bundle_meta.BundleJobStore(catalog, Path(state_dir), assets)
    app.state.job_store = store
    app.state.catalog = catalog
    app.state.assets = assets
    app.state.bundle_state_dir = state_dir

    # Build meta + seed the caches for every bundle. This is the step that keeps every
    # h5ad closed: `_get_expression_cache` (app.py:1353) and
    # `_get_marker_heatmap_cache_entry` (app.py:1260) both return on a cache hit.
    for entry in catalog.entries:
        store.ensure(app, entry["id"])

    app.mount("/viewer-static", StaticFiles(directory=os.path.join(_HERE, "static")),
              name="viewer-static")
    # The original memmap/binary API, unchanged, so the fast path does not regress.
    app.mount("/fast", create_fast_app(catalog), name="fast")

    _install_bundle_state_colors(app, catalog)
    _install_index_override(app)
    _install_catalog_routes(app, catalog, store)
    _install_differential_select(app, store)
    _install_dotplot_routes(app, store, assets)
    return app


# -------------------------------------------------------------- cell-state colours

def _install_bundle_state_colors(app, catalog: da.Catalog) -> None:
    """Cell-type colours come from the bundle's `cluster_colors`, not a generated ramp.

    scALABLE assigns categorical colours with `_build_preview_palette` (app.py:654), a
    deterministic paired ramp. A bundle ships explicit per-state colours in its metadata
    JSON. This wraps the palette function - it does not edit app.py - so that a request
    whose labels are all known cell states gets the bundle's colours, and every other
    request (sample colours in `Cell frequency`, for example) falls through to
    scALABLE's own ramp unchanged. The same map is served to the browser so the
    on-screen plot and the downloaded PDF agree.

    Every expression, heatmap and volcano gradient stays scALABLE's.
    """
    colors: Dict[str, str] = {}
    for entry in catalog.entries:
        ds = catalog.get(entry["id"])
        for state, value in (ds.meta.get("cluster_colors") or {}).items():
            colors.setdefault(str(state), str(value))
    app.state.bundle_state_colors = colors
    if not colors:
        return

    original = W._build_preview_palette

    def bundle_palette(populations):
        labels = [str(p) for p in populations]
        if labels and all(label in colors for label in labels):
            import matplotlib
            return {label: matplotlib.colors.to_rgb(colors[label]) for label in labels}
        return original(populations)

    W._build_preview_palette = bundle_palette

    @app.get("/api/jobs/{job_id}/state-colors")
    def state_colors(job_id: str):
        ds = catalog.get(job_id)
        return {"cluster_key": ds.cluster_key,
                "colors": {s: ds.colors.get(s) for s in ds.states if ds.colors.get(s)},
                "states": ds.states,
                "source": os.path.join(ds.paths.bundle_dir, ds.paths.prefix + "_metadata.json")}


# ---------------------------------------------------------------------------- index

def _install_index_override(app) -> None:
    """Serve scALABLE's index.html verbatim, plus one bootstrap script tag.

    The template is not edited or copied. The extra script auto-loads the first
    dataset, hides the upload workflow a precomputed bundle cannot use, and adds the
    dataset and contrast selectors.
    """
    bootstrap = os.path.join(_HERE, "static", "viewer_bootstrap.js")

    @app.middleware("http")
    async def inject_bootstrap(request: Request, call_next):
        response = await call_next(request)
        if request.url.path != "/" or response.status_code != 200:
            return response
        chunks = [chunk async for chunk in response.body_iterator]
        body = b"".join(chunks).decode("utf-8")
        version = int(os.path.getmtime(bootstrap))
        tag = (f'<link rel="stylesheet" href="/viewer-static/viewer.css?v={version}">'
               f'<script src="/viewer-static/viewer_bootstrap.js?v={version}"></script></body>')
        return HTMLResponse(body.replace("</body>", tag), status_code=200)


# -------------------------------------------------------------------------- catalog

def _install_catalog_routes(app, catalog: da.Catalog, store) -> None:

    @app.get("/api/catalog")
    def api_catalog():
        rows = []
        for entry in catalog.entries:
            meta = store.get_job(entry["id"])
            diff = meta.get("scalable_viewer", {}).get("deg_comparisons", [])
            rows.append({
                "id": entry["id"], "label": entry["label"],
                "n_cells": entry["n_cells"], "n_genes": entry["n_genes"],
                "n_states": entry["n_states"], "bundle_dir": entry["bundle_dir"],
                "prefix": entry["prefix"], "built_utc": entry["built_utc"],
                "contrasts": diff,
                "has_markers": entry["has_markers"],
                "fastcomm": bool((meta.get("fastcomm_analysis") or {}).get("enabled")),
                "default_gene": meta.get("default_gene"),
            })
        return {"datasets": rows, "load_errors": catalog.load_errors}


# ------------------------------------------------------------- differential contrast

def _install_differential_select(app, store) -> None:
    """Switch between the bundle's precomputed contrasts.

    scALABLE's POST /api/jobs/{id}/differential queues a differential run
    (app.py:4536). A bundle's contrasts are already computed, so this route is
    inserted ahead of it and only re-points meta['differential'] at the requested
    precomputed table. Every downstream view then runs through the unmodified
    scALABLE builders.
    """

    @app.get("/api/jobs/{job_id}/differential/contrasts")
    def contrasts(job_id: str):
        meta = store.get_job(job_id)
        return {"contrasts": meta.get("scalable_viewer", {}).get("deg_comparisons", []),
                "selected": (meta.get("differential") or {}).get("run_id")}

    @app.post("/api/jobs/{job_id}/differential/select")
    def select(job_id: str, contrast: str = Query(...)):
        ds = store.dataset(job_id)
        per_state = [c for c in ds.deg_manifest().get("comparisons", [])
                     if c.get("kind") == "per_cell_state"]
        chosen = next((c for c in per_state if c["id"] == contrast), None)
        if chosen is None:
            raise HTTPException(404, f"unknown contrast: {contrast}. "
                                     f"known: {[c['id'] for c in per_state]}")
        categorical = bundle_meta._categorical_covariates(ds)
        diff_assets = (app.state.assets.get(job_id, {}) or {}).get("differential", {}) or {}
        block = bundle_meta.build_differential_block(ds, chosen, categorical,
                                                     diff_assets.get(chosen["id"]))
        W._invalidate_differential_cache(app, job_id)
        store.update_job(job_id, differential=block)
        meta = store.get_job(job_id)
        return JSONResponse(W._build_differential_payload(app, job_id, meta, root_path=""))

    # Intercept the run POST: a precomputed bundle has nothing to queue.
    from fastapi.routing import APIRoute

    async def run_differential_blocked(job_id: str):
        meta = store.get_job(job_id)
        return JSONResponse(W._build_differential_payload(app, job_id, meta, root_path=""))

    for index, route in enumerate(list(app.router.routes)):
        if isinstance(route, APIRoute) and route.path == "/api/jobs/{job_id}/differential" \
                and "POST" in route.methods:
            app.router.routes.pop(index)
            app.post("/api/jobs/{job_id}/differential")(run_differential_blocked)
            break


# --------------------------------------------------------------------------- dotplot

def _dotplot_default_genes(assets: Dict[str, Any], ds: da.Dataset) -> List[Dict[str, str]]:
    default = (assets or {}).get("dotplot_default") or {}
    pairs = default.get("pairs") or []
    if pairs:
        return pairs
    # Fall back to the bundle's own marker table, top marker per state, canonical order.
    by_state: Dict[str, str] = {}
    for row in ds.markers():
        by_state.setdefault(str(row["cluster"]), str(row["gene"]))
    return [{"state": s, "gene": by_state[s]} for s in ds.states if s in by_state]


def _install_dotplot_routes(app, store, assets: Dict[str, Dict[str, Any]]) -> None:

    @app.get("/api/jobs/{job_id}/dotplot")
    def dotplot(job_id: str, genes: str = Query("")):
        """Mean and detected fraction per (gene, cell state), from the bundle's
        precomputed stats matrices. Default gene set = top marker of every state."""
        ds = store.dataset(job_id)
        pairs = _dotplot_default_genes(assets.get(job_id, {}), ds)
        requested = [g.strip() for g in genes.split(",") if g.strip()]
        wanted = requested or [p["gene"] for p in pairs]
        rows, missing, labels = [], [], []
        seen = set()
        for gene in wanted:
            if gene in seen:
                continue
            seen.add(gene)
            row = ds.resolve_gene(gene)
            if row is None:
                missing.append(gene)
                continue
            rows.append(row)
            labels.append(gene)
        if not rows:
            raise HTTPException(404, f"none of the {len(wanted)} requested genes are in this dataset")
        mean = np.asarray(ds.stats_mean[rows, :], dtype=np.float32)
        frac = np.asarray(ds.stats_frac[rows, :], dtype=np.float32)
        return {
            "genes": labels, "states": ds.states, "state_n": ds.state_n,
            "state_colors": [ds.colors.get(s, "#BBBBBB") for s in ds.states],
            "mean": mean.tolist(), "frac": frac.tolist(),
            "layer": ds.sv.get("layer"),
            "default_pairs": pairs, "is_default": not requested,
            "n_requested": len(wanted), "n_returned": len(labels),
            "n_missing": len(missing), "missing": missing,
        }
