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

import io
import json
import os
import re
import sqlite3
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

# The served page carries this name. scALABLE is the analysis tool; this deployment
# serves precomputed bundles, so it is named for what it is.
VIEWER_NAME = "scALABLE-viewer"


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
    _install_goelite_significance_tiers()
    _install_index_override(app)
    _install_catalog_routes(app, catalog, store)
    _install_differential_select(app, store)
    _install_dotplot_routes(app, store, assets)
    _install_study_route(app, catalog)
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


# -------------------------------------------------------- GO-Elite colour tiers
#
# The differential GO scatter drew two colours: blue for `is_selected_positive_sig`
# (app.py:1728 = the GO-Elite `selected` flag AND FDR<=0.05 AND Z>2) and grey for
# everything else, with `showlegend: false` (app.js:3687).
#
# `selected` is GO-Elite's DAG pruning, not a significance call. prio.py:34-35 keeps a
# term only when |Z|>=1.96, FDR<=0.1 AND at least 3 query genes overlap it, and
# prio.py:46-54 then drops it when a parent or child term already represents it. So the
# blue set is "representative term", and an FDR=0.0145 term carried by a single gene is
# correctly not blue. The plot was not wrong; it was unlabelled, and grey collapsed
# "significant but not representative" together with "not significant".
#
# This adds a third colour and a legend. It rewrites neither webapp/app.py nor
# webapp/static/app.js: the payload builder and the PDF renderer are wrapped on the
# module, the same way `_build_preview_palette` is wrapped above.

_GO_SIG_THRESHOLD = 0.05

# One hue family, dark to light, so the evidence order reads before the legend does.
_GO_TIERS = (
    ("representative", "#1f19c7", "GO-Elite representative"),
    ("significant", "#60a5fa", "Significant, not representative"),
    ("other", "#d1d5db", "Not significant"),
)
_GO_TIER_COLORS = {key: color for key, color, _ in _GO_TIERS}

_GO_SIG_NOTE = (
    "Significant = Fisher FDR < {threshold}. A GO-Elite representative term also needs "
    ">=3 overlapping query genes, |Z| >= 1.96,\nand no parent or child term that already "
    "represents it (goelite/prio.py:34-54); {n_one_gene} of the {n_sig} significant terms "
    "rest on a single gene."
)


def _go_term_tier(term: Dict[str, Any]) -> str:
    """representative / significant / other for one payload term."""
    if bool(term.get("is_selected_positive_sig")):
        return "representative"
    for key in ("fdr_plot", "fdr", "p_value"):
        value = term.get(key)
        if value is None:
            continue
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if number != number:  # NaN
            continue
        return "significant" if number <= _GO_SIG_THRESHOLD else "other"
    return "other"


def _install_goelite_significance_tiers() -> None:
    # The wrap is on the shared webapp module, so a second create_scalable_app in the
    # same process must not stack a second copy of it.
    if getattr(W, "_scalable_viewer_go_tiers_installed", False):
        return
    original_payload = W._build_differential_go_payload

    def go_payload(app, meta, population):
        payload = original_payload(app, meta, population)
        terms = payload.get("terms") or []
        counts = {key: 0 for key, _, _ in _GO_TIERS}
        n_one_gene = 0
        for term in terms:
            tier = _go_term_tier(term)
            term["color_tier"] = tier
            n_genes = len(term.get("overlap_genes") or [])
            term["n_overlap_genes"] = n_genes
            counts[tier] += 1
            if tier != "other" and n_genes <= 1:
                n_one_gene += 1
        n_sig = counts["representative"] + counts["significant"]
        payload["significance"] = {
            "fdr_threshold": _GO_SIG_THRESHOLD,
            "n_terms": len(terms),
            "n_significant": n_sig,
            "counts": counts,
            "tiers": [
                {"key": key, "color": color, "label": label, "n": counts[key]}
                for key, color, label in _GO_TIERS
            ],
            "note": _GO_SIG_NOTE.format(
                threshold=_GO_SIG_THRESHOLD, n_one_gene=n_one_gene, n_sig=n_sig
            ),
        }
        return payload

    W._build_differential_go_payload = go_payload
    W._render_differential_go_pdf = _render_go_pdf
    W._scalable_viewer_go_tiers_installed = True


def _render_go_pdf(payload: Dict[str, Any]) -> io.BytesIO:
    """The scALABLE GO scatter with three colour tiers and a legend.

    Axes, scale, labels, annotation offsets and the term-label arrows are the ones
    webapp/app.py:4027 `_render_differential_go_pdf` draws. Only the point colouring and
    the legend differ. Output stays vector: scatter markers and 2-point leader lines.
    """
    from matplotlib.lines import Line2D

    plt = W.plt
    W._configure_matplotlib_pdf_style()
    terms = payload.get("terms", []) or []
    if not terms:
        raise HTTPException(status_code=404, detail="No differential GO terms were available.")

    groups: Dict[str, List[tuple]] = {key: [] for key, _, _ in _GO_TIERS}
    x_values: List[float] = []
    y_values: List[float] = []
    for term in terms:
        z_score = float(term.get("z_score", np.nan) or np.nan)
        fdr_plot = float(term.get("fdr_plot", np.nan) or np.nan)
        if not (W._is_finite_number(z_score) and W._is_finite_number(fdr_plot) and fdr_plot > 0):
            continue
        tier = str(term.get("color_tier") or "")
        if tier not in groups:
            tier = _go_term_tier(term)
        x_values.append(z_score)
        y_values.append(fdr_plot)
        groups[tier].append((z_score, fdr_plot))
    if not x_values or not y_values:
        raise HTTPException(status_code=404, detail="No differential GO terms were available.")

    fig, ax = plt.subplots(figsize=(9.2, 7.8))
    # Weakest evidence first, so the representative terms are never hidden underneath.
    for zorder, (key, color, _) in enumerate(reversed(_GO_TIERS), start=2):
        points = groups[key]
        if not points:
            continue
        ax.scatter(
            [entry[0] for entry in points],
            [entry[1] for entry in points],
            s=48 if key == "representative" else 42,
            c=color,
            alpha=0.98 if key == "representative" else 0.95,
            linewidths=0,
            zorder=zorder,
        )

    labels = payload.get("labels", []) or []
    annotation_offsets = [(-6, 38), (-2, 10), (2, -10), (6, -34), (18, -60)]
    for index, label in enumerate(labels):
        z_score = float(label.get("z_score", np.nan) or np.nan)
        fdr_plot = float(label.get("fdr_plot", np.nan) or np.nan)
        if not (W._is_finite_number(z_score) and W._is_finite_number(fdr_plot) and fdr_plot > 0):
            continue
        dx, dy = annotation_offsets[min(index, len(annotation_offsets) - 1)]
        ax.annotate(
            str(label.get("term_name", "")),
            xy=(z_score, fdr_plot),
            xytext=(dx, dy),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=9,
            color=str(label.get("label_color", "#111827")),
            arrowprops={
                "arrowstyle": "-",
                "color": str(label.get("label_color", "#111827")),
                "linewidth": 1.0,
                "alpha": 0.9,
                "shrinkA": 0,
                "shrinkB": 0,
            },
            zorder=6,
        )

    x_min = min(-10.0, float(np.floor(min(x_values) - 0.5)))
    x_max = max(20.0, float(np.ceil(max(x_values) + 2.5)))
    y_min = max(min(y_values) * 0.5, 1e-300)
    ax.set_yscale("log")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, 1.0)
    ax.set_xlabel("Z-Score")
    ax.set_ylabel("Fishers FDR p")
    ax.set_title(f"GO terms: {payload.get('population', '')}")
    ax.axvline(0.0, color="#111827", linewidth=1.2, alpha=0.95, zorder=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    fig.tight_layout()

    significance = payload.get("significance") or {}
    counts = significance.get("counts") or {}
    total = int(significance.get("n_terms") or len(terms))
    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markersize=7,
            markerfacecolor=color,
            markeredgecolor="none",
            label=f"{label} ({int(counts.get(key, 0))} of {total})",
        )
        for key, color, label in _GO_TIERS
    ]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.085),
        ncol=3,
        frameon=False,
        fontsize=8,
        handletextpad=0.4,
        columnspacing=1.4,
    )
    note = str(significance.get("note") or "")
    if note:
        ax.text(
            0.5,
            -0.145,
            note,
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=7,
            color="#374151",
        )

    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


# ---------------------------------------------------------------------------- index

def _install_index_override(app) -> None:
    """Serve scALABLE's index.html verbatim, plus one bootstrap script tag.

    The template is not edited or copied. The extra script auto-loads the first
    dataset, hides the upload workflow a precomputed bundle cannot use, and adds the
    dataset and contrast selectors.

    The page is also renamed. This deployment serves precomputed bundles, so it is the
    *viewer*, not the analysis tool that shares the name. The response body is rewritten
    on the way out - `<title>` (webapp/templates/index.html:6, filled by app_title) and
    `<h1>` (index.html:19) - so neither app.py nor the template changes, and the served
    HTML already carries the new name before any script runs.
    """
    bootstrap = os.path.join(_HERE, "static", "viewer_bootstrap.js")

    @app.middleware("http")
    async def inject_bootstrap(request: Request, call_next):
        response = await call_next(request)
        if request.url.path != "/" or response.status_code != 200:
            return response
        chunks = [chunk async for chunk in response.body_iterator]
        body = b"".join(chunks).decode("utf-8")
        body = body.replace("<title>scALABLE</title>", f"<title>{VIEWER_NAME}</title>")
        body = body.replace("<h1>scALABLE</h1>", f"<h1>{VIEWER_NAME}</h1>")
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


# ---------------------------------------------------------------------- study record
#
# The Study tab reads the LungMAP site database instead of carrying hard-coded prose.
# The connection is opened read-only (`mode=ro` URI plus `PRAGMA query_only`), one
# connection per request, and it is closed in a `finally`. Nothing writes to the file.
#
# Schema used (breath.sqlite):
#   entity(id, namespace, class, graph, source_table, label, comment)
#   entity_value(subject_id, predicate, object_id, ordinal, value)   -- a view
# The dataset row carries rdfs:label, rdfs:comment, lmdb:is_data_type, lmdb:in_taxon,
# lmdb:has_sample_type, lmdb:has_dataset_sample_count, lmdb:has_ingest_stage.
# Linked rows point AT the dataset with lmdb:applies_to_dataset; entity.class names the
# kind (experiment_tool, experiment_file, experiment_sample, ...). A link row holds one
# `lmdb:has_<kind>` predicate whose object_id is the real record, so the real record's
# own properties (a file's URL and size, a sample's age and sex) are read in a second
# pass.

SITE_DB = os.environ.get(
    "LUNGMAP_SITE_DB",
    "/Users/saljh8/Dropbox/LungMAP/refactored_website/build/breath.sqlite")

# EXACTLY ONE id. There is deliberately no fallback to another study: showing a
# different study's title, description, tools, files and samples is presenting wrong
# data, whatever the notice says. When this record is absent from the site database the
# endpoint reads the SOURCE TABLE row instead (see read_study_record_from_tables), which
# is this study's own pending record, and reports source="source_tables".
STUDY_CANDIDATES = [s.strip() for s in os.environ.get(
    "LUNGMAP_STUDY_IDS", "lmdata:LMEX0000004416").split(",") if s.strip()]

# The site database is built from these TSVs, so a row here is this study's real record
# before publication. Columns are declared by the file's own #predicate header row.
SOURCE_TABLES = os.environ.get(
    "LUNGMAP_SOURCE_TABLES",
    "/Users/saljh8/Dropbox/LungMAP/refactored_website/build/lungmap-data-new/data/metadata")


# Every source TSV carries three header rows. `#predicate` names the RDF predicate each
# column writes and is the only header this reader needs; `#object_type` and
# `#search_weight` belong to the site's own loader. A fourth `row_type` line names the
# columns in plain words and is ignored. Every real row starts with the literal `data`.
def _tsv_rows(path, where=None, value=None, limit=None):
    """Rows of one source TSV as {predicate: value} dicts, in file order.

    The file is streamed and a row is kept only when its `where` predicate equals
    `value`, so a 5 MB link table costs one pass and holds only the matching rows.
    Returns [] when the file is absent, so a missing table degrades to "no rows"
    instead of raising inside a request.
    """
    if not os.path.exists(path):
        return []
    out, header = [], None
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "#predicate":
                header = parts
                continue
            if parts[0] != "data" or header is None:
                continue
            row = {header[i]: (parts[i] if i < len(parts) else "")
                   for i in range(len(header))}
            if where is not None and row.get(where, "") != value:
                continue
            out.append(row)
            if limit is not None and len(out) >= limit:
                break
    return out


def _tsv_by_id(path):
    """{identifier: {predicate: value}} for one entity or vocabulary table.

    The identifier column is the one whose `#predicate` header cell is `NA`.
    """
    index = {}
    for row in _tsv_rows(path):
        key = row.get("NA", "")
        if key:
            index[key] = row
    return index


# external_db/*.tsv holds one row per external accession, all of them under the same
# `lmdb:has_resource_id` predicate. external_db/external_api.tsv turns an accession into
# a URL: `lmdb:has_id_type` names the pipe-separated table types the template serves and
# `lmdb:has_link` is the template with `[ID]` where the accession goes. The reader uses
# that mechanism rather than hard-coding any host, so a new accession type needs no code.
_EXTERNAL_DB_FILES = ("geo_dataset.tsv", "geo_sample.tsv", "doi.tsv", "pmid.tsv",
                      "dbgap_dataset.tsv", "ega_dataset.tsv", "massive_dataset.tsv",
                      "proteome_xchange_dataset.tsv")


def _external_links(tables_root):
    """{xref_id: {"db", "accession", "name", "url", "icon"}} for the small tables."""
    templates = []
    for row in _tsv_rows(os.path.join(tables_root, "external_db", "external_api.tsv")):
        link = row.get("lmdb:has_link", "")
        if "[ID]" not in link:
            continue
        for kind in row.get("lmdb:has_id_type", "").split("|"):
            kind = kind.strip()
            if kind:
                templates.append((kind, row.get("rdfs:label", ""), link,
                                  row.get("lmdb:has_icon", "")))
    by_kind = {}
    for kind, label, link, icon in templates:
        by_kind.setdefault(kind, (label, link, icon))

    out = {}
    for name in _EXTERNAL_DB_FILES:
        for row in _tsv_rows(os.path.join(tables_root, "external_db", name)):
            xref_id = row.get("NA", "")
            kind = row.get("rdf:type", "")
            accession = row.get("lmdb:has_resource_id", "")
            if not xref_id or not accession:
                continue
            label, link, icon = by_kind.get(kind, ("", "", ""))
            out[xref_id] = {"db": kind, "accession": accession, "name": label or kind,
                            "url": link.replace("[ID]", accession) if link else "",
                            "icon": icon}
    return out


# One entry per link table this reader follows: the field it fills, the link table, the
# column naming the target record, and the entity or vocabulary table that holds the
# target's own label and description.
_SOURCE_LINKS = (
    ("tools", "lungmap_xref/dataset__tool.tsv", "lmdb:has_tool",
     "lungmap_data/tool.tsv"),
    ("files", "lungmap_xref/dataset__file.tsv", "lmdb:has_file",
     "lungmap_data/supporting_file.tsv"),
    ("samples", "lungmap_xref/dataset__sample.tsv", "lmdb:has_sample",
     "lungmap_data/sample.tsv"),
    ("researchers", "lungmap_xref/dataset__researcher.tsv", "lmdb:has_researcher",
     "lungmap_data/researcher.tsv"),
    ("publications", "lungmap_xref/dataset__publication.tsv", "lmdb:has_publication",
     "lungmap_data/publication.tsv"),
    ("technologies", "lungmap_xref/dataset__technology.tsv", "lmdb:has_technology",
     "lungmap_data/technology.tsv"),
    ("age_ranges", "lungmap_xref/dataset__age_range.tsv", "lmdb:has_age_range",
     "lungmap_vocabulary/age_range.tsv"),
)

_URL_IN_TEXT = re.compile(r"https?://[^\s]+")


def _link_order(row):
    try:
        return float(row.get("lmdb:has_display_order", ""))
    except (TypeError, ValueError):
        return 1e9


def read_study_record_from_tables(study_id, tables_root=SOURCE_TABLES):
    """This study's own record from the source TSVs, for when the site database lacks it.

    The site database is BUILT from these tables, so a row here is this study's real
    record before publication. This reader follows the same shape the database reader
    returns - tools, files, samples, researchers, publications, technologies,
    age_ranges, reference and raw_data - so the Study tab renders one way whichever
    source answered. It never invents a value: a field with no row stays empty.
    """
    bare = study_id.split(":", 1)[-1]
    head = _tsv_rows(os.path.join(tables_root, "lungmap_data", "dataset.tsv"),
                     where="NA", value=bare, limit=1)
    if not head:
        return None
    by_pred = head[0]

    def pretty(term):
        t = str(term or "")
        for prefix in ("data_type_", "sample_type_", "ingest_stage_"):
            if t.startswith(prefix):
                t = t[len(prefix):]
        return t.replace("_", " ").strip()

    out = {}
    for field, link_rel, head_pred, target_rel in _SOURCE_LINKS:
        links = _tsv_rows(os.path.join(tables_root, *link_rel.split("/")),
                          where="lmdb:applies_to_dataset", value=bare)
        if not links:
            out[field] = []
            continue
        targets = _tsv_by_id(os.path.join(tables_root, *target_rel.split("/")))
        entries = []
        for link in sorted(links, key=_link_order):
            target_id = link.get(head_pred, "")
            target = targets.get(target_id, {})
            record = {
                "record_id": link.get("NA", ""),
                "target_id": target_id,
                # A target with no row of its own still shows its identifier, so a
                # dangling reference is visible instead of blank.
                "name": target.get("rdfs:label", "") or target_id,
                "description": target.get("rdfs:comment", ""),
            }
            path = link.get("lmdb:has_path", "")
            if path:
                record["path"] = path
                record["url"] = (SITE_BASE + path) if path.startswith("/") else path
            for extra_pred, key in (("lmdb:has_role", "role"),
                                    ("lmdb:has_site", "site"),
                                    ("lmdb:has_icon", "icon")):
                if link.get(extra_pred):
                    record[key] = link[extra_pred]
            if target.get("lmdb:has_icon"):
                record.setdefault("icon", target["lmdb:has_icon"])
            entries.append(record)
        out[field] = entries

    external = _external_links(tables_root)

    # The paper's URL. A publication row carries no URL column, so the link comes from
    # publication__db_xref -> external_db, exactly as the other 14 publications do. The
    # publisher URL written into the label wins when it is there, because it is the page
    # a reader is sent to; the DOI resolver is the fallback.
    pub_xrefs = _tsv_rows(os.path.join(tables_root, "lungmap_xref",
                                       "publication__db_xref.tsv"))
    by_publication = {}
    for row in pub_xrefs:
        by_publication.setdefault(row.get("lmdb:applies_to_publication", ""), []).append(
            row.get("oboInOwl:hasDbXref", ""))
    for pub in out.get("publications", []):
        found = _URL_IN_TEXT.search(pub.get("name", ""))
        urls = []
        for xref_id in by_publication.get(pub.get("target_id", ""), []):
            entry = external.get(xref_id)
            if entry and entry["url"]:
                urls.append(entry)
                pub.setdefault(entry["db"], entry["accession"])
        pub["url"] = (found.group(0).rstrip(".,;)") if found
                      else (urls[0]["url"] if urls else ""))
        pub["links"] = [{"name": e["name"], "url": e["url"]} for e in urls]
        # The raw label is kept under "label". "name" is the citation a reader sees, so
        # the URL is taken out of it once it has become the link target.
        pub["label"] = pub.get("name", "")
        if found:
            pub["name"] = _URL_IN_TEXT.sub("", pub["label"]).strip().rstrip(".").strip()

    # Raw data: every external accession this dataset declares, resolved to a URL.
    accessions = []
    for row in _tsv_rows(os.path.join(tables_root, "lungmap_xref",
                                      "dataset__db_xref.tsv"),
                         where="lmdb:applies_to_dataset", value=bare):
        entry = external.get(row.get("oboInOwl:hasDbXref", ""))
        if entry and entry["url"]:
            accessions.append({"name": "%s %s" % (entry["name"], entry["accession"]),
                               "url": entry["url"], "db": entry["db"],
                               "accession": entry["accession"]})
    out["accessions"] = accessions

    def pick(entries, hosts):
        for entry in entries:
            url = (entry.get("url") or entry.get("path") or "").lower()
            if any(host in url for host in hosts):
                return {"name": entry.get("name", ""), "url": entry.get("url", "")}
        return None

    raw_data = (accessions[0] if accessions
                else pick(out["files"], _RAW_DATA_HOSTS) or pick(out["tools"], _RAW_DATA_HOSTS))
    if raw_data:
        raw_data = {"name": raw_data.get("name", ""), "url": raw_data.get("url", "")}
    reference = None
    if out["publications"]:
        first = out["publications"][0]
        reference = {"name": first.get("name", ""), "url": first.get("url", "")}
    if not (reference and reference["url"]):
        hit = pick(out["files"], _REFERENCE_HOSTS) or pick(out["tools"], _REFERENCE_HOSTS)
        if hit:
            reference = hit if not reference else {"name": reference["name"],
                                                   "url": hit["url"]}

    # The database reader returns these two as plain label lists. Match it, so the
    # front end joins the same shape whichever source answered.
    out["technologies"] = [t["name"] for t in out["technologies"]]
    out["age_ranges"] = [a["name"] for a in out["age_ranges"]]

    out.update({
        "dataset_id": bare,
        "study_id": bare,
        "title": by_pred.get("rdfs:label", ""),
        "description": by_pred.get("rdfs:comment", ""),
        "assay": pretty(by_pred.get("lmdb:is_data_type")),
        "sample_type": pretty(by_pred.get("lmdb:has_sample_type")),
        "organism": "Homo sapiens" if "9606" in str(by_pred.get("lmdb:in_taxon", "")) else "",
        "cell_count": by_pred.get("lmdb:has_dataset_sample_count", ""),
        "ingest_stage": pretty(by_pred.get("lmdb:has_ingest_stage")),
        "release_date": by_pred.get("lmdb:has_release_date", ""),
        "reference": reference,
        "raw_data": raw_data,
        "site_base": SITE_BASE,
        "tables_root": tables_root,
        "n_tools": len(out["tools"]), "n_files": len(out["files"]),
        "n_samples": len(out["samples"]), "n_researchers": len(out["researchers"]),
        "source": "source_tables",
        "published": False,
        "notice": ("This study is not published to the site database yet. The record below "
                   "is its own pending entry from the LungMAP source tables."),
    })
    return out

# Several tool rows store a site-relative path ("/cell-cards/"). The public host is the
# one the site itself uses: app/lungmap/web/data_routes.py:650 and api/routes.py:5.
SITE_BASE = os.environ.get("LUNGMAP_SITE_BASE", "https://www.lungmap.net").rstrip("/")

_LINK_CLASS_FIELD = {
    "experiment_tool": "tools",
    "experiment_file": "files",
    "experiment_sample": "samples",
    "experiment_researcher": "researchers",
    "experiment_publication": "publications",
    "experiment_technology": "technologies",
    "experiment_age_range": "age_ranges",
}
_LINK_CLASS_HEAD = {
    "experiment_tool": "lmdb:has_tool",
    "experiment_file": "lmdb:has_file",
    "experiment_sample": "lmdb:has_sample",
    "experiment_researcher": "lmdb:has_researcher",
    "experiment_publication": "lmdb:has_publication",
    "experiment_technology": "lmdb:has_technology",
    "experiment_age_range": "lmdb:has_age_range",
}
_RAW_DATA_HOSTS = ("ncbi.nlm.nih.gov/geo", "ncbi.nlm.nih.gov/sra", "ncbi.nlm.nih.gov/gap",
                   "data-browser.lungmap.net", "ega-archive.org", "dbgap")
_REFERENCE_HOSTS = ("doi.org", "pubmed", "nature.com", "sciencedirect", "cell.com")


def _study_connect(path: str) -> sqlite3.Connection:
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=15.0)
    con.execute("PRAGMA query_only = ON")
    return con


def _short_key(predicate: str) -> str:
    key = predicate.split(":", 1)[-1]
    for prefix in ("has_", "in_", "is_"):
        if key.startswith(prefix):
            return key[len(prefix):]
    return key


def _props_for(cur, ids: List[str]) -> Dict[str, List[tuple]]:
    """(predicate, object_id, ordinal, value) for every requested subject id."""
    out: Dict[str, List[tuple]] = {}
    unique = list(dict.fromkeys(ids))
    for start in range(0, len(unique), 400):
        chunk = unique[start:start + 400]
        query = ("SELECT subject_id, predicate, object_id, ordinal, value FROM entity_value "
                 "WHERE subject_id IN (%s)" % ",".join("?" * len(chunk)))
        for sid, pred, obj, ordinal, value in cur.execute(query, chunk):
            out.setdefault(sid, []).append((pred, obj, ordinal, value))
    return out


def _first(props: List[tuple], predicate: str):
    for pred, _obj, _ordinal, value in props:
        if pred == predicate:
            return value
    return None


def _object_of(props: List[tuple], predicate: str):
    for pred, obj, _ordinal, _value in props:
        if pred == predicate:
            return obj
    return None


def _order_of(props: List[tuple]) -> float:
    raw = _first(props, "lmdb:has_display_order")
    try:
        return float(raw)
    except (TypeError, ValueError):
        return 1e9


def read_study_record(db_path: str, candidates: List[str]) -> Dict[str, Any]:
    """One read-only pass over breath.sqlite for one dataset record."""
    if not os.path.isfile(db_path):
        return {"ok": False, "error": f"site database not found: {db_path}",
                "db_path": db_path, "candidates": candidates}
    con = _study_connect(db_path)
    try:
        cur = con.cursor()
        found = None
        for candidate in candidates:
            row = cur.execute("SELECT id, label, comment FROM entity WHERE id = ?",
                              (candidate,)).fetchone()
            if row is not None:
                found = row
                break
        if found is None:
            return {"ok": False, "error": "none of the candidate study ids is in the database",
                    "db_path": db_path, "candidates": candidates,
                    "missing": candidates}
        study_id = found[0]
        head = _props_for(cur, [study_id]).get(study_id, [])

        # Link rows: everything that points AT this dataset.
        link_rows = cur.execute(
            "SELECT e.class, ev.subject_id FROM entity_value ev "
            "JOIN entity e ON e.id = ev.subject_id "
            "WHERE ev.predicate = 'lmdb:applies_to_dataset' AND ev.object_id = ?",
            (study_id,)).fetchall()
        by_class: Dict[str, List[str]] = {}
        for cls, sid in link_rows:
            by_class.setdefault(cls, []).append(sid)
        link_props = _props_for(cur, [sid for _cls, sid in link_rows])

        # Second pass: the records the link rows point to (file URLs, sample attributes).
        targets: List[str] = []
        for cls, ids in by_class.items():
            predicate = _LINK_CLASS_HEAD.get(cls)
            if not predicate:
                continue
            for sid in ids:
                obj = _object_of(link_props.get(sid, []), predicate)
                if obj:
                    targets.append(obj)
        target_props = _props_for(cur, targets)

        out: Dict[str, Any] = {field: [] for field in _LINK_CLASS_FIELD.values()}
        for cls, ids in by_class.items():
            field = _LINK_CLASS_FIELD.get(cls)
            predicate = _LINK_CLASS_HEAD.get(cls)
            if not field or not predicate:
                continue
            entries = []
            for sid in sorted(ids, key=lambda s: _order_of(link_props.get(s, []))):
                props = link_props.get(sid, [])
                name = _first(props, predicate)
                if name is None:
                    continue
                record: Dict[str, Any] = {"name": str(name), "record_id": sid}
                target = _object_of(props, predicate)
                if target:
                    record["target_id"] = target
                # Every other lmdb: property on the link row (role, site, path, ...).
                for pred, _obj, _ordinal, value in props:
                    if pred in (predicate, "lmdb:applies_to_dataset", "lmdb:has_display_order",
                                "rdf:type") or value is None:
                        continue
                    record.setdefault(_short_key(pred), str(value))
                # The target record's own properties (URL, size, age, sex, ...).
                for pred, _obj, _ordinal, value in target_props.get(target, []):
                    if pred in ("rdf:type", "rdfs:label") or value is None:
                        continue
                    if pred == "rdfs:comment":
                        record.setdefault("description", str(value))
                        continue
                    record.setdefault(_short_key(pred), str(value))
                entries.append(record)
            out[field] = entries

        # tools[] must expose {name, path}; files[] a URL the browser can follow.
        for tool in out["tools"]:
            tool["path"] = tool.get("path") or tool.get("display_url") or ""
            tool["url"] = (SITE_BASE + tool["path"]) if tool["path"].startswith("/") \
                else tool["path"]
        for entry in out["files"]:
            entry["url"] = entry.get("display_url") or ""
            entry["size"] = entry.get("file_size") or ""
        for pub in out["publications"]:
            pub["url"] = pub.get("display_url") or pub.get("doi") or ""
        out["technologies"] = [t["name"] for t in out["technologies"]]
        out["age_ranges"] = [a["name"] for a in out["age_ranges"]]

        def pick(urls_from, hosts):
            for entry in urls_from:
                url = entry.get("url") or entry.get("path") or ""
                if any(host in url.lower() for host in hosts):
                    return {"name": entry.get("name", ""), "url": url}
            return None

        reference = None
        if out["publications"]:
            first_pub = out["publications"][0]
            reference = {"name": first_pub["name"], "url": first_pub.get("url", "")}
        if not (reference and reference["url"]):
            hit = pick(out["files"], _REFERENCE_HOSTS) or pick(out["tools"], _REFERENCE_HOSTS)
            if hit:
                reference = hit if not reference else {"name": reference["name"], "url": hit["url"]}
        raw_data = pick(out["files"], _RAW_DATA_HOSTS) or pick(out["tools"], _RAW_DATA_HOSTS)

        out.update({
            "ok": True,
            "db_path": db_path,
            "db_mtime": os.path.getmtime(db_path),
            "site_base": SITE_BASE,
            "candidates": candidates,
            "study_id": study_id,
            "requested_id": candidates[0] if candidates else None,
            "is_fallback": bool(candidates) and study_id != candidates[0],
            "dataset_id": study_id.split(":", 1)[-1],
            "title": _first(head, "rdfs:label") or found[1] or study_id,
            "description": _first(head, "rdfs:comment") or found[2] or "",
            "assay": _first(head, "lmdb:is_data_type") or "",
            "organism": _first(head, "lmdb:in_taxon") or "",
            "sample_type": _first(head, "lmdb:has_sample_type") or "",
            "cell_count": _first(head, "lmdb:has_dataset_sample_count") or "",
            "ingest_stage": _first(head, "lmdb:has_ingest_stage") or "",
            "reference": reference,
            "raw_data": raw_data,
            "n_tools": len(out["tools"]), "n_files": len(out["files"]),
            "n_samples": len(out["samples"]), "n_researchers": len(out["researchers"]),
        })
        return out
    finally:
        con.close()


def viewer_meta_samples(ds: da.Dataset, column: str = "meta_sample") -> Dict[str, Any]:
    """The viewer's own sample list, read from the bundle's obs covariates.

    Used when the site database record carries no sample rows yet. One row per level of
    `meta_sample`, with the number of metacells behind it and the single value of every
    covariate that is constant within that meta-sample (sex, copd_status, Group, ...).
    A covariate that varies inside a meta-sample is reported as "mixed", never as one of
    its values.
    """
    covariates = ds.covariate_names()
    if column not in covariates:
        return {"column": column, "available": False, "rows": [], "n": 0,
                "columns": [], "error": f"'{column}' is not a covariate of this bundle"}
    _kind, codes, categories = ds.covariate_values(column)
    codes = np.asarray(codes, dtype=np.int64)
    n_levels = len(categories)
    counts = np.bincount(codes[codes >= 0], minlength=n_levels)

    extras = [name for name in ("sex", "copd_status", "cancer_status", "Group", "clincluster",
                                "Collection.Method", "meta_sample_n_donors", "n_cells_total")
              if name in covariates]
    values: Dict[str, List[Any]] = {}
    for name in extras:
        kind, arr, cats = ds.covariate_values(name)
        arr = np.asarray(arr)
        column_values: List[Any] = []
        for level in range(n_levels):
            sel = arr[codes == level]
            if sel.size == 0:
                column_values.append("")
                continue
            unique = np.unique(sel)
            if unique.size != 1:
                column_values.append("mixed")
            elif kind == "numeric":
                value = float(unique[0])
                column_values.append(int(value) if value == int(value) else round(value, 3))
            else:
                index = int(unique[0])
                column_values.append(cats[index] if 0 <= index < len(cats) else "")
        values[name] = column_values

    rows = []
    for level, label in enumerate(categories):
        row = {"name": str(label), "metacells": int(counts[level])}
        for name in extras:
            row[name] = values[name][level]
        rows.append(row)
    return {"column": column, "available": True, "n": len(rows), "rows": rows,
            "columns": ["name", "metacells"] + extras,
            "source": os.path.join(ds.paths.bundle_dir, ds.paths.prefix + "_metadata.json")}


def _install_study_route(app, catalog: da.Catalog) -> None:

    @app.get("/api/study")
    def api_study(dataset: str = Query(""), study_id: str = Query("")):
        """The LungMAP study record for the Study tab, read-only from breath.sqlite.

        `study_id` overrides the candidate list. `dataset` names the viewer bundle whose
        own meta-samples are returned as `viewer_samples`, used by the browser when the
        site record has no sample rows yet.
        """
        candidates = [study_id] if study_id else list(STUDY_CANDIDATES)
        try:
            record = read_study_record(SITE_DB, candidates)
        except sqlite3.Error as err:
            # The site database is rebuilt in place by another process. A locked or
            # half-written file is reported to the browser, never swallowed.
            record = {"ok": False, "db_path": SITE_DB, "candidates": candidates,
                      "error": f"sqlite3.{type(err).__name__}: {err}"}
        # Not in the site database yet. Read THIS study's own pending row from the
        # source tables. Never answer with a different study: a record that is not the
        # requested one is wrong data, and a notice does not make it right.
        if not record.get("ok"):
            db_error = record.get("error", "")
            pending = read_study_record_from_tables(candidates[0])
            if pending:
                pending["ok"] = True
                pending["requested_id"] = candidates[0]
                # The provenance line must name the file the record actually came from.
                # Reporting breath.sqlite here would credit a database that does not
                # hold this row. The database path stays beside it, with the reason it
                # did not answer.
                pending["site_db"] = SITE_DB
                pending["site_db_error"] = db_error
                pending["db_path"] = pending.get("tables_root", SOURCE_TABLES)
                record = pending
        else:
            record.setdefault("source", "site_database")
            record.setdefault("published", True)

        ds_id = dataset or (catalog.entries[0]["id"] if catalog.entries else "")
        record["dataset_bundle"] = ds_id
        if ds_id:
            try:
                record["viewer_samples"] = viewer_meta_samples(catalog.get(ds_id))
            except Exception as err:                      # noqa: BLE001 - reported, not hidden
                record["viewer_samples"] = {"available": False, "rows": [], "n": 0,
                                            "columns": [], "error": f"{type(err).__name__}: {err}"}
            entry = next((e for e in catalog.entries if e["id"] == ds_id), None)
            if entry:
                record["viewer_stats"] = {"n_cells": entry["n_cells"], "n_genes": entry["n_genes"],
                                          "n_states": entry["n_states"], "label": entry["label"]}
        if not record.get("ok"):
            return JSONResponse(record, status_code=200)
        return record
