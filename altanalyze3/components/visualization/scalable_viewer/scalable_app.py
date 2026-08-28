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

import contextvars
import io
import json
import os
import re
import sqlite3
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from fastapi import Body, HTTPException, Query, Request
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

# The header's "How to Use" and "ReadMe" links. webapp/templates/index.html:30-31 points
# both at scALABLE's own documents, which describe the upload-and-run workflow this
# deployment removes. The viewer's own pair lives at
# altanalyze3/components/visualization/scalable_viewer/HOW_TO_USE.md and README.md, on the
# master branch of https://github.com/SalomonisLab/altanalyze3. The shared template is not
# edited: _install_index_override rewrites the served body, and viewer_bootstrap.js
# repeats the swap in the DOM so a cached page also lands on the right document.
_GITHUB_BLOB = "https://github.com/SalomonisLab/altanalyze3/blob/master/altanalyze3/components"
DOC_LINKS = {
    f"{_GITHUB_BLOB}/cellHarmony/webapp/HOW_TO_USE.md":
        f"{_GITHUB_BLOB}/visualization/scalable_viewer/HOW_TO_USE.md",
    f"{_GITHUB_BLOB}/cellHarmony/webapp/README.md":
        f"{_GITHUB_BLOB}/visualization/scalable_viewer/README.md",
}


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
    _install_violin_covariate(app)
    _install_index_override(app)
    _install_catalog_routes(app, catalog, store)
    _install_differential_select(app, store)
    # scALABLE builds its examples from job metadata, which for a bundle yields
    # placeholder contrast names such as "Group 1 versus Group 2". The bundle
    # knows its real cell states and comparisons, so it answers this instead.
    # scALABLE gained its own chat route, which computes from a job's h5ad via
    # adata.X. A bundle has no adata, so the official handler raised
    # "'BundleAnnData' object has no attribute 'X'" for every question the
    # bundle-backed executors answer. FastAPI matches the first route
    # registered, so the official one is removed here.
    _drop_official_route(app, "/api/jobs/{job_id}/chat", "POST")
    _drop_official_route(app, "/api/jobs/{job_id}/chat-examples")
    _drop_official_route(app, "/api/jobs/{job_id}/dotplot")
    _drop_official_route(app, "/api/jobs/{job_id}/combplot")
    _install_dotplot_routes(app, store, assets)
    _install_combplot_routes(app, store, assets)
    _install_chat_routes(app, store, assets)
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


# ------------------------------------------------------ violin grouping covariate
#
# scALABLE groups the violin by one column only: the cluster key. `_build_expression_payload`
# (app.py:3552-3564) walks `cache_entry["populations"]`, which `_get_expression_cache`
# fills from `adata.obs[cluster_key]` (app.py:1396). The same cache entry already holds
# every other obs column, keyed by name, in `obs_filter_values` (app.py:1432, and for a
# bundle bundle_meta.py:156-166) - the exact list the "Annotation 1" selector offers, so
# the new "Covariate" control needs no new precomputation and no new bundle field.
#
# The contract is a `covariate` query parameter:
#   * absent            -> scALABLE's own builder runs untouched, so today's payload is
#                          returned byte for byte,
#   * `covariate=<col>` -> the violin, and ONLY the violin, is regrouped by that column.
#                          Scatter, UMAP points, UMAP colouring and the global range stay
#                          the cluster-key payload scALABLE built.
#
# The default the browser sends is the cluster key itself, which takes the regroup path.
# That is deliberate and testable: `?covariate=cell_state` must return exactly the payload
# the parameter-free request returns, which proves the regrouping loop below is scALABLE's
# loop and not a lookalike.
#
# Neither webapp/app.py nor webapp/static/app.js is edited. The two expression routes are
# re-registered ahead of the originals purely to declare the extra query parameter; each
# then delegates to scALABLE's own endpoint object, so the error mapping, the store lookup
# and the response headers are still scALABLE's.

# Group ORDER and CAP are scALABLE's own rule, kept as-is: order by mean expression
# descending, draw at most 10 (app.py:3564). Keeping the cap identical is what makes the
# default grouping identical. A covariate with more groups than the cap is NOT silently
# truncated: the payload carries the totals below and the plot title names them.
_VIOLIN_MAX_GROUPS = 10

# Set by the expression routes, read by the wrapped payload builder. It is set inside the
# route coroutine that awaits the original endpoint, so it is visible for the whole of
# that request and for nothing else.
_VIOLIN_COVARIATE: contextvars.ContextVar[str] = contextvars.ContextVar(
    "scalable_viewer_violin_covariate", default=""
)


def _violin_note(covariate: str, n_total: int, n_shown: int, is_default: bool) -> str:
    """The truncation statement that goes on the plot.

    The default covariate keeps scALABLE's own wording verbatim, so the default plot is
    unchanged. Every other covariate states the denominator as well as the cap.
    """
    if is_default:
        return f"top {n_shown} states by mean"
    if n_total > n_shown:
        return f"top {n_shown} of {n_total} groups by mean"
    return f"all {n_total} groups by mean"


def _regroup_violin(app, meta, payload, covariate, modality, display_filters):
    """Rebuild `payload['violin']` grouped by `covariate`, with scALABLE's own loop."""
    cache_entry = W._get_expression_cache(app, meta, modality=W._normalize_modality_id(modality))
    cluster_key = str(cache_entry.get("cluster_key") or "")
    labels_by_field = cache_entry.get("obs_filter_values") or {}
    labels = labels_by_field.get(covariate)
    if labels is None:
        payload["violin_covariate"] = covariate
        payload["violin_covariate_error"] = (
            f"'{covariate}' is not one of this dataset's annotation columns "
            f"({len(labels_by_field)} available); the violin is grouped by '{cluster_key}'."
        )
        return payload

    gene = str(payload.get("resolved_gene") or payload.get("gene") or "")
    resolved = W._resolve_gene_name(cache_entry["var_names"], gene)
    if not resolved:
        # A reference-centroid payload (app.py:3369) has no per-cell matrix to regroup.
        payload["violin_covariate"] = covariate
        payload["violin_covariate_error"] = (
            f"'{gene}' has no per-cell values in this dataset, so the violin cannot be "
            f"regrouped by '{covariate}'."
        )
        return payload

    values = W._flatten_expr(cache_entry["adata"][:, resolved].X).astype(float)
    display_mask = W._apply_display_filter_mask(cache_entry, display_filters)
    labels = np.asarray(labels, dtype=str)

    # app.py:3552-3563 verbatim, with `labels` in place of `populations`. The one addition
    # is the empty-string group: an obs column can leave a cell unannotated, and "" is not
    # a biological group. The cluster key never carries it, which is why the default
    # payload is unaffected.
    groups: List[Dict[str, Any]] = []
    n_unlabeled = 0
    for label in sorted(pd.unique(labels)):
        mask = (labels == label) & display_mask
        if label == "":
            n_unlabeled = int(np.count_nonzero(mask))
            continue
        group_values = values[mask]
        finite_values = group_values[np.isfinite(group_values)]
        groups.append(
            {
                "population": label,
                "values": [float(v) for v in finite_values],
                "mean": float(np.mean(finite_values)) if len(finite_values) else 0.0,
            }
        )
    groups = sorted(groups, key=lambda entry: entry["mean"], reverse=True)

    is_default = covariate == cluster_key
    payload["violin"] = groups[:_VIOLIN_MAX_GROUPS]
    payload["violin_covariate"] = covariate
    payload["violin_cluster_key"] = cluster_key
    payload["violin_covariate_is_default"] = bool(is_default)
    payload["violin_n_groups_total"] = len(groups)
    payload["violin_n_groups_shown"] = len(payload["violin"])
    payload["violin_n_cells_unlabeled"] = n_unlabeled
    payload["violin_note"] = _violin_note(
        covariate, len(groups), len(payload["violin"]), is_default
    )
    return payload


def _install_violin_covariate(app) -> None:
    """Wrap the payload builder and the violin PDF, then declare the query parameter."""
    if not getattr(W, "_scalable_viewer_violin_covariate_installed", False):
        original_payload = W._build_expression_payload

        # **kwargs, deliberately: this shim stands in for
        # `_build_expression_payload`, and every parameter that function gains
        # upstream must pass straight through. Naming them one by one meant the
        # new `violin_limit` raised "unexpected keyword argument" the moment it
        # was added.
        def expression_payload(app_, meta, gene, modality="rna",
                               display_filters=None, **extra):
            payload = original_payload(
                app_, meta, gene, modality=modality,
                display_filters=display_filters, **extra
            )
            covariate = str(_VIOLIN_COVARIATE.get("") or "").strip()
            if not covariate:
                return payload
            return _regroup_violin(app_, meta, payload, covariate, modality, display_filters)

        original_pdf = W._render_expression_pdf

        def expression_pdf(payload, mode):
            # The default covariate goes through scALABLE's own renderer, so the default
            # PDF stays exactly the PDF this route produced before. Only a non-default
            # covariate needs the covariate title and the x-axis label.
            if (
                mode == "violin"
                and payload.get("violin_covariate")
                and payload.get("violin")
                and not payload.get("violin_covariate_is_default")
            ):
                return _render_violin_pdf(payload)
            return original_pdf(payload, mode)

        W._build_expression_payload = expression_payload
        W._render_expression_pdf = expression_pdf
        W._scalable_viewer_violin_covariate_installed = True

    _install_expression_covariate_routes(app)


def _install_expression_covariate_routes(app) -> None:
    """Re-register the two expression routes with a `covariate` parameter in front of
    scALABLE's own, and delegate the body to scALABLE's own endpoint objects."""
    from fastapi.routing import APIRoute

    targets = {}
    for index, route in enumerate(list(app.router.routes)):
        if not isinstance(route, APIRoute) or "GET" not in route.methods:
            continue
        if route.path in ("/api/jobs/{job_id}/expression", "/api/jobs/{job_id}/expression/pdf"):
            targets.setdefault(route.path, (index, route.endpoint))
    if len(targets) != 2:
        raise RuntimeError(
            "scalable_viewer: expected the two scALABLE expression routes, found "
            f"{sorted(targets)}"
        )

    json_index, json_endpoint = targets["/api/jobs/{job_id}/expression"]

    async def expression_with_covariate(
        job_id: str,
        gene: str = Query(...),
        modality: str = Query("rna"),
        covariate: str = Query(""),
        filter1_field: Optional[str] = Query(None),
        filter1_values: List[str] = Query([]),
        filter2_field: Optional[str] = Query(None),
        filter2_values: List[str] = Query([]),
        violin_limit: int = Query(10),
    ):
        token = _VIOLIN_COVARIATE.set(str(covariate or "").strip())
        try:
            return await json_endpoint(
                job_id,
                gene=gene,
                modality=modality,
                # Every parameter the wrapped route declares must be passed
                # explicitly. Calling it directly bypasses FastAPI, so an
                # unpassed argument arrives as the Query object itself, and the
                # route raised "int() argument must be ... not 'Query'" the
                # moment `violin_limit` was added upstream.
                violin_limit=violin_limit,
                filter1_field=filter1_field,
                filter1_values=filter1_values,
                filter2_field=filter2_field,
                filter2_values=filter2_values,
            )
        finally:
            _VIOLIN_COVARIATE.reset(token)

    pdf_index, pdf_endpoint = targets["/api/jobs/{job_id}/expression/pdf"]

    async def expression_pdf_with_covariate(
        job_id: str,
        gene: str = Query(...),
        mode: str = Query("umap"),
        modality: str = Query("rna"),
        covariate: str = Query(""),
        filter1_field: Optional[str] = Query(None),
        filter1_values: List[str] = Query([]),
        filter2_field: Optional[str] = Query(None),
        filter2_values: List[str] = Query([]),
    ):
        token = _VIOLIN_COVARIATE.set(str(covariate or "").strip())
        try:
            return await pdf_endpoint(
                job_id,
                gene=gene,
                mode=mode,
                modality=modality,
                filter1_field=filter1_field,
                filter1_values=filter1_values,
                filter2_field=filter2_field,
                filter2_values=filter2_values,
            )
        finally:
            _VIOLIN_COVARIATE.reset(token)

    # Registered last, then moved in front of the original: the router matches the first
    # route whose path and method match, so position is what makes the override effective.
    for path, handler, index in (
        ("/api/jobs/{job_id}/expression", expression_with_covariate, json_index),
        ("/api/jobs/{job_id}/expression/pdf", expression_pdf_with_covariate, pdf_index),
    ):
        app.get(path)(handler)
        app.router.routes.insert(index, app.router.routes.pop())


def _render_violin_pdf(payload: Dict[str, Any]) -> io.BytesIO:
    """The scALABLE violin PDF (app.py:3826-3858) with the covariate in the title.

    Reached only for a NON-default covariate; the default covariate is drawn by
    scALABLE's own `_render_expression_pdf`, unchanged.

    Figure size, violin bodies, colours, jitter seed, marker size, tick rotation and the
    y-limit padding are the ones `_render_expression_pdf` draws. Only the title text and
    the x-axis label differ, so a covariate PDF and the on-screen plot say the same thing.
    """
    plt = W.plt
    W._configure_matplotlib_pdf_style()
    fig, ax = plt.subplots(figsize=(8.5, 8.5))
    gene = payload["gene"]
    global_min = float(payload.get("global_min", 0.0) or 0.0)
    global_max = float(payload.get("global_max", 0.0) or 0.0)
    if global_max <= global_min:
        global_max = global_min + 1e-9
    violin_data = payload["violin"]
    positions = np.arange(1, len(violin_data) + 1)
    parts = ax.violinplot(
        [entry["values"] for entry in violin_data],
        positions=positions,
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )
    for body in parts["bodies"]:
        body.set_facecolor("#94a3b8")
        body.set_edgecolor("#475569")
        body.set_alpha(0.55)
    for idx, entry in enumerate(violin_data, start=1):
        vals = np.asarray(entry["values"], dtype=float)
        jitter = np.random.default_rng(0).normal(0, 0.035, size=len(vals))
        ax.scatter(np.full(len(vals), idx) + jitter, vals, s=4, c="#0f172a", alpha=0.35, linewidths=0)
    ax.set_xticks(positions)
    ax.set_xticklabels([entry["population"] for entry in violin_data], rotation=45, ha="right")
    covariate = str(payload.get("violin_covariate") or "")
    note = str(payload.get("violin_note") or "")
    ax.set_title(f"{gene} expression by {covariate} ({note})")
    ax.set_xlabel(covariate)
    ax.set_ylabel("Expression")
    pad = max((global_max - global_min) * 0.04, 0.05)
    ax.set_ylim(global_min - pad, global_max + pad)
    fig.tight_layout()
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

    The two documentation links (index.html:30-31) are repointed the same way. The shared
    template sends them to scALABLE's own HOW_TO_USE.md and README.md, which describe the
    upload-and-run workflow this deployment removes. The viewer ships its own pair, so the
    served page carries those instead. The template is not edited: the swap is a body
    rewrite here, and viewer_bootstrap.js repeats it in the DOM for a cached page.
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
        for old, new in DOC_LINKS.items():
            body = body.replace(old, new)
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


def split_gene_list(text: str) -> List[str]:
    """Gene symbols out of whatever a user pasted.

    A user pastes a column out of Excel, which arrives newline separated, or a
    row, which arrives tab separated. Others type commas, spaces or semicolons.
    Splitting on commas alone accepted only one of those: `SFTPC AGER` resolved
    to a single symbol named "SFTPC AGER" and the request 404ed.

    Duplicates are dropped and order is kept, so the plot reads in the order the
    user listed the genes.
    """
    parts = re.split(r"[\s,;|]+", str(text or ""))
    seen, out = set(), []
    for part in parts:
        gene = part.strip().strip('"').strip("'")
        if gene and gene not in seen:
            seen.add(gene)
            out.append(gene)
    return out


def _donor_covariate(ds) -> Optional[str]:
    """The covariate that names the donor, or None.

    Tried in the order a bundle is likely to spell it. `meta_sample` holds 178
    donors in the COPD atlas; other bundles use `Library`, `donor` or `sample`.
    """
    covariates = ds.covariate_names() or {}
    for name in ("meta_sample", "donor", "Donor", "Library", "sample", "sample_id", "pool"):
        info = covariates.get(name)
        if info and info.get("kind") == "categorical":
            return name
    return None


def _bundle_subset_mask(ds, subset_by: str, subset_values: List[str]):
    """Cells a bundle-backed figure is restricted to, or None for all of them.

    Restricting to one cell state and then grouping by a clinical variable is
    what turns "COPD versus control" from an atlas-wide average into the
    question the user meant: within AT2, what differs between disease and
    control.
    """
    if not subset_by or not subset_values:
        return None
    keep = set(str(v) for v in subset_values)
    cluster_key = ds.sv.get("cluster_key") or "cell_state"
    if subset_by in (cluster_key, "cell_state"):
        codes = np.asarray(ds.state_code, dtype=np.int64)
        wanted = {i for i, s in enumerate(ds.states) if s in keep}
        if not wanted:
            return None
        return np.isin(codes, list(wanted))
    try:
        kind, values, labels = ds.covariate_values(subset_by)
    except KeyError:
        return None
    if kind != "categorical" or labels is None:
        return None
    wanted = {i for i, s in enumerate(labels) if str(s) in keep}
    if not wanted:
        return None
    return np.isin(np.asarray(values, dtype=np.int64), list(wanted))


def _bundle_subset_masks(ds, subset_by: str, subset_values: List[str],
                         subset2_by: str = "", subset2_values: List[str] = None):
    """Cells that satisfy BOTH display filters, or None when neither restricts.

    The viewer's "Filter data to display" block carries two annotation rows, so
    a figure can be held to one cell state AND one clinical level at the same
    time. The two masks combine with AND: a cell must satisfy every filter the
    user set to stay on the plot.
    """
    first = _bundle_subset_mask(ds, subset_by, list(subset_values or []))
    second = _bundle_subset_mask(ds, subset2_by, list(subset2_values or []))
    if first is None:
        return second
    if second is None:
        return first
    return first & second


def _bundle_group_axis(ds, group_by: str = ""):
    """The grouping variable for a bundle-backed DotPlot or CombPlot.

    Cell state is the default and is the fast path: the bundle already carries
    per-state means and fractions. Any other variable is a covariate, so the
    values are grouped per cell instead.
    """
    column = str(group_by or "").strip()
    if not column or column == (ds.sv.get("cluster_key") or "cell_state"):
        return "cell_state", list(ds.states), np.asarray(ds.state_code, dtype=np.int64), True
    kind, values, labels = ds.covariate_values(column)
    if kind != "categorical" or labels is None:
        raise HTTPException(400, f"{column} is not categorical, so it cannot group this plot")
    return column, [str(v) for v in labels], np.asarray(values, dtype=np.int64), False


def _bundle_default_genes(ds, groups) -> List[str]:
    """One marker gene per group, for a blank gene set.

    The bundle ships a marker table keyed by cell state, so that is used when the
    grouping is cell state. Any other variable has no such table, and the top
    markers of the states are still the informative genes to open on.
    """
    by_group: Dict[str, str] = {}
    for row in ds.markers():
        by_group.setdefault(str(row["cluster"]), str(row["gene"]))
    ordered = [by_group[g] for g in groups if g in by_group]
    if ordered:
        return ordered[:12]
    return [by_group[s] for s in ds.states if s in by_group][:12]


#: Donor covariates the CombPlot draws as annotation bands under the gene rows,
#: matching the binary bands of the ICGS CombPlot. Only the names a bundle
#: actually records are used, so a bundle without them simply draws no bands.
COMBPLOT_TRACK_DEFAULTS = ("copd_status", "Group", "sex", "Smoking Status")

#: A covariate with more levels than this needs one band per level and would
#: bury the figure, so it is reported as skipped instead of drawn.
COMBPLOT_TRACK_MAX_LEVELS = 40


def _combplot_tracks(ds, names, group, keep, n_groups):
    """The value of each covariate for every CombPlot column.

    A column is one (cell state, donor) pair, so a donor-level covariate holds
    one value across it. This function does not assume that. It counts the
    cells of every (column, level) pair, reports the most common level, and
    returns the smallest fraction that level covers over the columns, so a
    covariate that varies inside a column shows up as a purity below 1.0
    instead of hiding behind a single label.
    """
    valid = group >= 0
    at = group[valid]
    values: Dict[str, List[str]] = {}
    levels: Dict[str, List[str]] = {}
    purity: Dict[str, float] = {}
    used: List[str] = []
    skipped: List[str] = []
    for name in names:
        try:
            kind, codes, cats = ds.covariate_values(name)
        except KeyError:
            skipped.append(name)
            continue
        if kind != "categorical" or cats is None or not len(cats):
            skipped.append(name)
            continue
        if len(cats) > COMBPLOT_TRACK_MAX_LEVELS:
            skipped.append(name)
            continue
        codes = np.asarray(codes, dtype=np.int64)
        n_levels = len(cats)
        counts = np.bincount(at * n_levels + codes[valid],
                             minlength=n_groups * n_levels).reshape(n_groups, n_levels)
        sub = counts[keep] if len(keep) else np.zeros((0, n_levels), dtype=np.int64)
        totals = sub.sum(axis=1)
        top = sub.argmax(axis=1) if sub.shape[0] else np.zeros(0, dtype=np.int64)
        best = sub.max(axis=1) if sub.shape[0] else np.zeros(0, dtype=np.int64)
        share = np.where(totals > 0, best / np.maximum(totals, 1), 1.0)
        per_column = [str(cats[int(t)]) for t in top]
        seen = set(per_column)
        values[name] = per_column
        levels[name] = [str(c) for c in cats if str(c) in seen]
        purity[name] = round(float(share.min()) if share.size else 1.0, 4)
        used.append(name)
    return used, values, levels, purity, skipped


def _install_combplot_routes(app, store, assets: Dict[str, Dict[str, Any]]) -> None:

    @app.get("/api/jobs/{job_id}/combplot")
    def combplot(job_id: str, genes: str = Query(""), donor_key: str = Query(""),
                 min_cells: int = Query(5), group_by: str = Query(""),
                 groups: List[str] = Query([]), subset_by: str = Query(""),
                 subset_values: List[str] = Query([]), subset2_by: str = Query(""),
                 subset2_values: List[str] = Query([]), tracks: str = Query("")):
        """Per-donor pseudobulk for each gene, grouped by cell state.

        One bar per (cell state, donor). The bar is the mean of the gene over
        that donor's cells in that state, so a state with 178 donors draws 178
        bars, coloured by the state. Cell-level values would draw 123,076 bars
        and hide the donor structure the plot exists to show.

        States run in the bundle's canonical order, which is the centroid
        ordering `lineage_order` records, so the x axis reads the same way as
        every other plot in the viewer.

        The mean is over the cells present, so a donor contributing no cells to
        a state produces no bar rather than a zero. `min_cells` raises that
        floor; the response reports how many groups it removed.
        """
        ds = store.dataset(job_id)
        group_column, group_names, group_code, _ = _bundle_group_axis(ds, group_by)
        if groups:
            wanted_groups = [g for g in group_names if g in set(groups)]
            if wanted_groups:
                group_names = wanted_groups
        # Blank means the marker gene of every group, matching the DotPlot.
        wanted = split_gene_list(genes) or _bundle_default_genes(ds, group_names)
        if not wanted:
            raise HTTPException(400, "give at least one gene")

        key = donor_key.strip() or _donor_covariate(ds)
        if not key:
            raise HTTPException(404, "this bundle records no donor covariate")
        try:
            kind, donor_code, donor_labels = ds.covariate_values(key)
        except KeyError:
            raise HTTPException(404, f"no covariate named {key}")
        if kind != "categorical" or donor_labels is None:
            raise HTTPException(400, f"{key} is not categorical, so it cannot name donors")

        donor_code = np.asarray(donor_code, dtype=np.int64)
        keep_names = {g: i for i, g in enumerate(group_names)}
        all_names, _ = list(group_names), None
        # Re-index the per-cell codes onto the groups that survived the filter.
        original = _bundle_group_axis(ds, group_by)[1]
        remap = np.full(len(original), -1, dtype=np.int64)
        for i, name in enumerate(original):
            if name in keep_names:
                remap[i] = keep_names[name]
        state_code = remap[np.asarray(group_code, dtype=np.int64)]
        restrict = _bundle_subset_masks(ds, subset_by, list(subset_values),
                                        subset2_by, list(subset2_values))
        if restrict is not None:
            state_code = np.where(restrict, state_code, -1)
        n_states, n_donors = len(group_names), len(donor_labels)
        group = np.where(state_code >= 0, state_code * n_donors + donor_code, -1)
        n_groups = n_states * n_donors
        valid = group >= 0
        cells_per_group = np.bincount(group[valid], minlength=n_groups)

        # Columns are the groups that hold cells, ordered state-major so the
        # plot reads left to right in canonical state order.
        keep = np.nonzero(cells_per_group >= max(1, int(min_cells)))[0]
        keep = keep[np.argsort(keep, kind="stable")]
        dropped = int(np.count_nonzero((cells_per_group > 0)
                                       & (cells_per_group < max(1, int(min_cells)))))

        columns = [{"group": group_names[int(g) // n_donors],
                    "state": group_names[int(g) // n_donors],
                    "donor": donor_labels[int(g) % n_donors],
                    "n_cells": int(cells_per_group[int(g)])} for g in keep]
        colors = [ds.colors.get(c["group"], "#BBBBBB") for c in columns]

        # The annotation bands under the gene rows. Blank asks for the default
        # donor covariates this bundle records.
        asked = [t.strip() for t in str(tracks).split(",") if t.strip()]
        if not asked:
            recorded = ds.covariate_names() or {}
            asked = [t for t in COMBPLOT_TRACK_DEFAULTS if t in recorded]
        (track_names, track_values, track_levels,
         track_purity, track_skipped) = _combplot_tracks(ds, asked, group, keep, n_groups)

        series, labels, missing = [], [], []
        for gene in wanted:
            row = ds.resolve_gene(gene)
            if row is None:
                missing.append(gene)
                continue
            idx, val = ds.gene_column(row)
            sums = np.zeros(n_groups, dtype=np.float64)
            if idx.size:
                cells = idx.astype(np.int64)
                inside = group[cells] >= 0
                if inside.any():
                    np.add.at(sums, group[cells][inside], val.astype(np.float64)[inside])
            with np.errstate(invalid="ignore", divide="ignore"):
                means = np.where(cells_per_group > 0, sums / np.maximum(cells_per_group, 1), 0.0)
            labels.append(gene)
            series.append([round(float(v), 5) for v in means[keep]])

        if not series:
            raise HTTPException(404, f"none of the {len(wanted)} requested genes are in this dataset")

        return {
            "genes": labels, "values": series,
            "columns": columns, "colors": colors,
            "states": group_names, "groups": group_names,
            "group_by": group_column, "group_label": group_column,
            "subset_by": subset_by or "", "subset_values": list(subset_values),
            "subset2_by": subset2_by or "", "subset2_values": list(subset2_values),
            "n_cells_kept": int(np.count_nonzero(restrict)) if restrict is not None
                            else int(np.asarray(ds.state_code).shape[0]),
            "track_names": track_names, "tracks": track_values,
            "track_levels": track_levels, "track_purity": track_purity,
            "track_skipped": track_skipped,
            "donor_key": key,
            "n_donors": n_donors, "n_states": n_states,
            "n_columns": len(columns), "n_groups_dropped": dropped,
            "layer": ds.sv.get("layer"),
            "n_requested": len(wanted), "n_returned": len(labels),
            "n_missing": len(missing), "missing": missing,
        }


#: Where the viewer sends a chat question to be read. The model runs in the
#: LungMAP site process; this viewer holds none. Override with
#: SCALABLE_ASSISTANT_URL when the site is not on the same host.
ASSISTANT_URL = os.environ.get(
    "SCALABLE_ASSISTANT_URL", "http://127.0.0.1:8001/api/assistant/viewer-intent")


def _markers_for_states(ds, states: List[str], limit: int = 25) -> List[Dict[str, Any]]:
    """Marker genes for these cell states, from the table or computed.

    The bundle's marker table covers 14 of this atlas's 39 states, so DC1, DC2,
    NK cells, Mast cells and 21 others returned an empty table. The per-state
    mean and detected fraction exist for all 39, so where the table is silent
    the markers are computed as a one-versus-rest contrast on those: the gap
    between a gene's mean inside the state and its mean everywhere else. The
    `source` field on every row says which of the two produced it, because a
    stored fold and a computed gap are not the same statistic.
    """
    wanted = [s for s in states if s]
    if not wanted:
        return []
    rows = [dict(r, source="marker table") for r in ds.markers()
            if str(r.get("cluster")) in set(wanted)]
    covered = {str(r["cluster"]) for r in rows}
    missing = [s for s in wanted if s not in covered]
    if not missing:
        return rows[:limit]

    mean = np.asarray(ds.stats_mean, dtype=np.float32)
    counts = np.asarray(ds.state_n, dtype=np.float64)
    total = counts.sum()
    per_state = max(1, limit // max(1, len(wanted)))
    for state in missing:
        index = ds.states.index(state)
        inside = mean[:, index]
        # Mean everywhere else, weighted by how many cells each state holds.
        others = (mean * counts).sum(axis=1) - inside * counts[index]
        outside = others / max(1.0, total - counts[index])
        gap = inside - outside
        for row in np.argsort(-gap)[:per_state]:
            rows.append({"gene": ds.symbols[int(row)], "cluster": state,
                         "fold": round(float(gap[int(row)]), 4),
                         "p": None, "source": "computed one-vs-rest"})
    return rows[:limit]


def _donor_axis(ds):
    """Per-cell donor codes and the donor labels, or (None, None)."""
    key = _donor_covariate(ds)
    if not key:
        return None, None
    kind, values, labels = ds.covariate_values(key)
    if kind != "categorical" or labels is None:
        return None, None
    return np.asarray(values, dtype=np.int64), [str(v) for v in labels]


def _donor_pseudobulk(ds, rows: List[int], state: str = "",
                      min_cells: int = 5):
    """Mean of each gene per donor, optionally inside one cell state.

    The matrix every per-donor protocol needs: genes down, donors across. A
    donor contributing fewer than `min_cells` cells is dropped rather than
    represented by a mean over one or two cells.

    Returns (matrix, donors, n_cells) with matrix shape (len(rows), n_donors).
    """
    donor_code, donor_labels = _donor_axis(ds)
    if donor_code is None:
        return None, [], []
    keep = np.ones(donor_code.shape[0], dtype=bool)
    if state:
        try:
            keep = np.asarray(ds.state_code, dtype=np.int64) == ds.states.index(state)
        except ValueError:
            return None, [], []
    counts = np.bincount(donor_code[keep], minlength=len(donor_labels))
    usable = np.nonzero(counts >= max(1, int(min_cells)))[0]
    if not usable.size:
        return None, [], []

    index = np.full(len(donor_labels), -1, dtype=np.int64)
    index[usable] = np.arange(usable.size)
    per_cell = np.where(keep, index[donor_code], -1)

    matrix = np.zeros((len(rows), usable.size), dtype=np.float64)
    for position, row in enumerate(rows):
        idx, val = ds.gene_column(row)
        if idx.size:
            cells = idx.astype(np.int64)
            inside = per_cell[cells] >= 0
            if inside.any():
                np.add.at(matrix[position], per_cell[cells][inside],
                          val.astype(np.float64)[inside])
    matrix /= np.maximum(counts[usable], 1)
    return matrix, [donor_labels[i] for i in usable], [int(counts[i]) for i in usable]


def _candidate_rows(ds, limit: int = 3000) -> List[int]:
    """A spread of gene rows to scan, so a whole-transcriptome sweep stays quick."""
    step = max(1, len(ds.symbols) // limit)
    return list(range(0, len(ds.symbols), step))


def _follow_ups(intent: str, ds, state: str = "", covariate: str = "",
                genes: Optional[List[str]] = None) -> List[str]:
    """What to ask next, given what was just answered.

    Each suggestion is a question this viewer can answer, built from the same
    vocabulary the answer used, so a follow-up never lands on a protocol that
    has no recipe. The point is to carry an association forward: a gene that
    tracks age in one cell state raises the question of whether it does so in
    another, and whether the disease contrast moves the same gene.
    """
    genes = [g for g in (genes or []) if g]
    other = next((s for s in ds.states if s != state), "")
    contrasts = [(c.get("id") or "") for c in (ds.deg_manifest() or {}).get("comparisons", [])]
    case = ""
    per_state = next((c for c in contrasts if "per_cell_state" in c), "")
    if per_state:
        head = per_state.split("::")[0].replace("_", " ")
        case = head.split(" vs ")[0] if " vs " in head else head

    out: List[str] = []
    if intent == "severity_gradient" and state and covariate:
        if genes:
            out.append(f"Which genes co-vary with {genes[0]} across donors in {state}?")
        if other:
            out.append(f"Which genes in {other} track {covariate}?")
        if case:
            out.append(f"Which genes are significant in {case} in {state} cells?")
        out.append(f"What processes are up in {case or 'disease'} {state} cells?")
    elif intent in ("state_contrast", "donor_heterogeneity"):
        if state:
            out.append(f"Is the {state} {case} signature present in all donors or a subset?")
            out.append(f"What processes are up in {case} {state} cells?")
            out.append(f"Which genes in {state} track Age?")
        out.append(f"Which cell type is most affected in {case}?")
    elif intent == "cell_identity" and state:
        out += [f"What distinguishes {state} from {other}?",
                f"Which genes are significant in {case} in {state} cells?",
                f"Which genes in {state} track Age?"]
    elif intent == "composition_shift":
        out += [f"Which cell type is most affected in {case}?",
                f"Which genes are significant in {case} in {state or ds.states[0]} cells?"]
    elif intent == "coexpression_module" and state:
        out += [f"Which genes in {state} track Age?",
                f"What processes are up in {case} {state} cells?"]
    elif intent == "most_affected_state":
        out += [f"Which genes are significant in {case} in {state or ds.states[0]} cells?",
                f"Which cell states are depleted in {case}?"]
    seen, unique = set(), []
    for question in out:
        if question and question not in seen:
            seen.add(question)
            unique.append(question)
    return unique[:4]


def _install_chat_routes(app, store, assets: Dict[str, Dict[str, Any]]) -> None:

    @app.get("/api/jobs/{job_id}/chat-examples")
    def chat_examples(job_id: str):
        """Example questions built from THIS dataset's own vocabulary.

        The tab used to ship eight fixed lung sentences, so a bone-marrow job
        offered AT2 and COPD questions its data cannot answer. These are
        assembled from the states and contrasts the bundle actually holds, so
        every chip is answerable by the dataset that is loaded.
        """
        ds = store.dataset(job_id)
        states = list(ds.states)
        # A state with enough cells to answer a per-donor question.
        counts = list(ds.state_n)
        ranked = [s for _, s in sorted(zip(counts, states), reverse=True)]
        first = ranked[0] if ranked else ""
        second = ranked[1] if len(ranked) > 1 else first

        contrasts = [(c.get("id") or "") for c in (ds.deg_manifest() or {}).get("comparisons", [])]
        per_state = next((c for c in contrasts if "per_cell_state" in c), "")
        sides = per_state.split("::")[0].replace("_", " ") if per_state else ""
        case = sides.split(" vs ")[0] if " vs " in sides else sides

        covariates = ds.covariate_names() or {}
        numeric = next((n for n, i in covariates.items()
                        if i.get("kind") == "numeric" and not n.endswith("__n_obs")), "")
        categorical = next((n for n, i in covariates.items()
                            if i.get("kind") == "categorical"
                            and 1 < len(i.get("categories") or []) <= 8), "")

        examples = [f"What are the best marker genes of {first} cells?",
                    f"What distinguishes {first} from {second}?"]
        if case:
            examples += [
                f"Which genes are significant in {case} in {first} cells?",
                f"Is the {first} {case} signature present in all donors or a subset?",
                f"Which cell type is most affected in {case}?",
                f"What processes are up in {case} {first} cells?",
            ]
        if numeric:
            examples.append(f"Which genes in {first} track {numeric}?")
        if categorical:
            examples.append(f"Do the {categorical} groups differ molecularly in {first}?")
        examples.append(f"Show me the transcriptional targets of a regulator in {first}")

        return {"reference": ds.sv.get("label") or job_id,
                "examples": examples[:8],
                "placeholder": f"e.g. What are the best marker genes of {first} cells?"}

    @app.post("/api/jobs/{job_id}/chat")
    def chat(job_id: str, body: Dict[str, Any] = Body(...)):
        """Answer one question about this dataset, with numbers from the bundle.

        Two steps, deliberately separated:

        1. The LungMAP site reads the sentence and returns which supported query
           to run and on which cell state, contrast or genes. It sees no data.
        2. This route runs that query against the bundle and returns the result.

        So the model chooses the question and the bundle answers it. A model of
        this size asked for a fold change would invent one; asked only to pick
        between five readings it is reliable. Every number below is read off the
        precomputed matrices.
        """
        ds = store.dataset(job_id)
        question = str(body.get("question") or "").strip()
        if not question:
            raise HTTPException(400, "question is required")

        contrasts = [c.get("id") or c.get("label") or ""
                     for c in (ds.deg_manifest() or {}).get("comparisons", [])]
        contrasts = [c for c in contrasts if c]

        reading: Dict[str, Any]
        try:
            import urllib.request
            covariates = [name for name, info in (ds.covariate_names() or {}).items()
                          if info.get("kind") in ("numeric", "categorical")]
            modalities = [m.get("id") for m in
                          ((ds.sv.get("modalities") or {}).get("available") or [])] or ["rna"]
            payload = json.dumps({"question": question, "states": ds.states,
                                  "contrasts": contrasts,
                                  # Without these the router cannot fill a
                                  # covariate slot and every severity or
                                  # composition question returned `clarify`.
                                  "covariates": covariates,
                                  "modalities": modalities}).encode()
            req = urllib.request.Request(
                ASSISTANT_URL, data=payload,
                headers={"Content-Type": "application/json"})
            # The model runs on CPU: about 24 s on the first call, which loads
            # the 1.07 GB file, and about 15 s warm. A 30 s ceiling failed the
            # cold call outright, so the window is wide enough for both.
            # One retry. Routing is deterministic now, so a failure here is a
            # dropped connection rather than a slow answer, and it cost one
            # question in a 68-question sweep. The second attempt uses a fresh
            # request object because a urllib Request cannot be replayed.
            last_error = None
            for attempt in (1, 2):
                try:
                    retry = urllib.request.Request(
                        ASSISTANT_URL, data=payload,
                        headers={"Content-Type": "application/json"})
                    with urllib.request.urlopen(retry, timeout=120) as resp:
                        reading = json.loads(resp.read().decode())
                    last_error = None
                    break
                except Exception as exc:            # noqa: BLE001
                    last_error = exc
                    if attempt == 1:
                        time.sleep(0.25)
            if last_error is not None:
                raise last_error
        except Exception as exc:                  # noqa: BLE001 - the site may be down
            raise HTTPException(
                503,
                f"the assistant at {ASSISTANT_URL} did not answer ({exc}). "
                "The scALABLE viewer holds no model of its own.")

        intent = reading.get("intent")
        state = reading.get("cell_state") or ""
        state2 = reading.get("cell_state_2") or ""
        genes = reading.get("genes") or []
        contrast = reading.get("contrast") or ""
        covariate = reading.get("covariate") or ""

        result: Dict[str, Any] = {"question": question, "reading": reading,
                                  "intent": intent}

        if intent == "clarify":
            result["answer"] = (
                f"I need to know {reading.get('missing') or 'a little more'}. "
                f"This dataset has {len(ds.states)} cell states and "
                f"{len(contrasts)} precomputed comparisons.")
            result["choices"] = {"states": ds.states, "contrasts": contrasts}
            return result

        if intent == "unsupported":
            result["answer"] = (
                "I can answer four things about this dataset: the marker genes of a "
                "cell state, genes differing between two groups within a cell state, "
                "where named genes are expressed, and what separates two cell states.")
            return result

        # The router speaks protocol names. Several of them are answered by the
        # executors already written here; the rest are specified in
        # lungmap/assistant/viewer_protocols.py but have no endpoint yet, and
        # this says so rather than returning a neighbouring protocol's answer,
        # which is how "transcriptional targets of RUNX1" came back as the
        # marker genes of alveolar macrophages.
        PROTOCOL_ALIAS = {
            # New protocols that the executors below already answer.
            "cell_identity": "markers",
            "state_comparison": "compare",
            "expression_lookup": "expression",
            "state_contrast": "differential",
            "patient_stratification": "differential",
            "shared_vs_state_specific": "differential",
            "contrast_specificity": "differential",
            "pathway_program": "goelite",
            "regulatory_driver": "network",
            "communication_rewiring": "ccc",
        }
        NOT_YET = {
            "severity_gradient": "a Spearman correlation of per-donor pseudobulk against a clinical variable",
            "dose_response": "a monotonic trend test across ordered stages",
            "composition_shift": "per-donor cell-count composition",
            "coexpression_module": "per-donor co-expression around a seed gene",
            "annotation_concordance": "a cross-tabulation of the two annotations",
        }
        # --- protocols computed here, from the bundle -----------------------
        if intent == "severity_gradient" and state and covariate:
            answer = _run_severity_gradient(ds, state, covariate)
            result.update(answer)
            top = [r[0] for r in ((answer.get("table") or {}).get("rows") or [])[:1]]
            result["follow_ups"] = _follow_ups(intent, ds, state, covariate, top)
            return result
        if intent == "coexpression_module" and state and genes:
            result.update(_run_coexpression(ds, state, genes[0]))
            result["follow_ups"] = _follow_ups(intent, ds, state, covariate, genes)
            return result
        if intent == "composition_shift" and covariate:
            result.update(_run_composition(ds, covariate, question))
            result["follow_ups"] = _follow_ups(intent, ds, state, covariate, genes)
            return result
        if intent == "annotation_concordance":
            result.update(_run_concordance(ds, state))
            return result
        if intent == "dose_response" and state:
            result.update(_run_dose_response(ds, state, covariate or "Group", genes))
            return result

        if intent in NOT_YET:
            result["answer"] = (
                f"That is the {intent} protocol. It needs {NOT_YET[intent]}, which this "
                "viewer does not compute yet, so I am not going to answer it with a "
                "different analysis. The protocol is specified and the endpoint is the "
                "remaining work.")
            result["status"] = "not_implemented"
            return result
        # GO-Elite is shipped for every comparison, so this is answered here
        # rather than pointing at another tab.
        # Where the disease acts, not which genes move. Aliasing this to the
        # generic differential answered a cell-type question with a gene volcano.
        # Who carries the signature, not which genes define it. Aliasing this to
        # the generic differential answered a donor question with a gene volcano.
        if intent == "donor_heterogeneity" and state:
            result.update(_run_donor_heterogeneity(
                ds, state, contrast or (contrasts[0] if contrasts else "")))
            result["follow_ups"] = _follow_ups(intent, ds, state, covariate, genes)
            return result

        if intent == "most_affected_state":
            result.update(_run_most_affected_state(ds, contrast or (contrasts[0] if contrasts else "")))
            result["follow_ups"] = _follow_ups(intent, ds, state, covariate, genes)
            return result

        if intent == "pathway_program" and state:
            answer = _run_pathway_program(assets.get(job_id, {}), ds, state,
                                          contrast or (contrasts[0] if contrasts else ""),
                                          reading.get("direction") or "both")
            result.update(answer)
            result["follow_ups"] = _follow_ups(intent, ds, state, covariate, genes)
            return result

        if intent in ("regulatory_driver", "communication_rewiring", "pathway_program"):
            target = PROTOCOL_ALIAS[intent]
            label = {"network": "differential network",
                     "ccc": "cell-cell communication",
                     "goelite": "GO-Elite enrichment"}[target]
            result["answer"] = (
                f"That is the {intent} protocol, answered by this dataset's {label} "
                f"results for {state or 'the selected cell state'}. Open the "
                f"{'Differential' if target != 'ccc' else 'Explore'} tab for the figure; "
                "the chat does not yet inline it.")
            result["status"] = "use_existing_view"
            result["plot"] = {"kind": "network"}
            return result
        intent = PROTOCOL_ALIAS.get(intent, intent)


        if intent == "markers":
            rows = _markers_for_states(ds, [state], 25)
            result["answer"] = (f"{len(rows)} top marker genes of {state}, ranked as the "
                                f"bundle's marker table ranks them.")
            result["table"] = {
                "columns": ["gene", "cluster", "fold", "p", "source"],
                "rows": [[r.get("gene"), r.get("cluster"),
                          round(float(r.get("fold") or 0), 4),
                          (float(r["p"]) if r.get("p") is not None else None),
                          r.get("source", "marker table")] for r in rows],
            }
            result["plot"] = {"kind": "dotplot",
                              "genes": [r.get("gene") for r in rows[:12]]}
            result["follow_ups"] = _follow_ups(reading.get("intent") or "", ds, state,
                                               covariate, [r.get("gene") for r in rows[:1]])
            return result

        if intent == "expression":
            rows_idx, labels, missing = [], [], []
            for gene in genes:
                row = ds.resolve_gene(gene)
                if row is None:
                    missing.append(gene)
                else:
                    rows_idx.append(row)
                    labels.append(gene)
            if not rows_idx:
                result["answer"] = f"None of {', '.join(genes)} are in this dataset."
                return result
            mean = np.asarray(ds.stats_mean[rows_idx, :], dtype=np.float32)
            table_rows = []
            for i, gene in enumerate(labels):
                order = np.argsort(-mean[i])[:5]
                for j in order:
                    table_rows.append([gene, ds.states[int(j)],
                                       round(float(mean[i][int(j)]), 3),
                                       round(float(ds.stats_frac[rows_idx[i]][int(j)]), 3)])
            result["answer"] = (f"Highest-expressing cell states for "
                                f"{', '.join(labels)}, by mean of the {ds.sv.get('layer')} layer."
                                + (f" Not in this dataset: {', '.join(missing)}." if missing else ""))
            result["table"] = {"columns": ["gene", "cell state", "mean", "fraction"],
                               "rows": table_rows}
            result["plot"] = {"kind": "dotplot", "genes": labels}
            return result

        if intent == "differential":
            manifest = ds.deg_manifest() or {}
            comparisons = manifest.get("comparisons", [])
            chosen = None
            for comp in comparisons:
                cid = comp.get("id") or comp.get("label") or ""
                if contrast and cid == contrast:
                    chosen = comp
                    break
            if chosen is None and comparisons:
                chosen = comparisons[0]
            if chosen is None:
                result["answer"] = "This bundle carries no precomputed comparisons."
                return result
            table = ds.deg_table(chosen.get("id"), 25, 0.05, state or None)
            rows_out = (table.get("rows") if isinstance(table, dict) else table) or []
            if not rows_out and state:
                # The contrast is not computed for every cell state: the COPD
                # differential covers 18 of this atlas's 39. An empty table read
                # as "nothing changes here", which is a different and wrong
                # claim, so the absence is stated instead.
                covered = sorted({str(r.get("population") or r.get("cluster") or "")
                                  for r in (ds.deg_table(chosen.get("id"), 100000, None, None)
                                            or {}).get("rows", [])} - {""})
                result["answer"] = (
                    f"{chosen.get('label') or chosen.get('id')} is not computed for {state}. "
                    f"It covers {len(covered)} of the {len(ds.states)} cell states in this "
                    f"dataset. This is missing data, not an absence of change.")
                result["status"] = "not_covered"
                result["table"] = {"columns": ["cell state covered"],
                                   "rows": [[s] for s in covered]}
                result["plot"] = None
                return result
            result["answer"] = (
                f"Top genes for {chosen.get('label') or chosen.get('id')}"
                + (f" in {state}" if state else "")
                + ", FDR below 0.05, from the bundle's precomputed differential.")
            result["table"] = table if isinstance(table, dict) else {"rows": table}
            result["follow_ups"] = _follow_ups(reading.get("intent") or "", ds,
                                               state, covariate, genes)
            result["plot"] = {"kind": "volcano", "comparison": chosen.get("id")}
            return result

        if intent == "compare":
            want = {state, state2}
            rows = _markers_for_states(ds, [state, state2], 30)
            result["answer"] = f"Marker genes separating {state} and {state2}."
            result["table"] = {
                "columns": ["gene", "cluster", "fold", "p", "source"],
                "rows": [[r.get("gene"), r.get("cluster"),
                          round(float(r.get("fold") or 0), 4),
                          (float(r["p"]) if r.get("p") is not None else None),
                          r.get("source", "marker table")] for r in rows],
            }
            result["plot"] = {"kind": "dotplot",
                              "genes": [r.get("gene") for r in rows[:12]]}
            return result

        result["answer"] = "I did not understand that."
        return result


def _drop_official_route(app, path: str, method: str = "GET") -> None:
    """Remove scALABLE's own handler for `path`, so a bundle-backed one can serve it.

    scALABLE computes the DotPlot and CombPlot from a job's h5ad. A precomputed
    bundle has no h5ad: it ships the same statistics already reduced, which is
    why it can answer for 123,076 cells instantly. FastAPI matches the first
    route registered, so the official handler would win and then fail on a
    bundle. The route is removed and re-registered against the bundle here.
    """
    from fastapi.routing import APIRoute

    for index, route in enumerate(list(app.router.routes)):
        if isinstance(route, APIRoute) and route.path == path and method in route.methods:
            app.router.routes.pop(index)
            return


def _install_dotplot_routes(app, store, assets: Dict[str, Dict[str, Any]]) -> None:

    @app.get("/api/jobs/{job_id}/plot-variables")
    def plot_variables(job_id: str):
        """The variables this bundle's DotPlot and CombPlot may group or filter by.

        Cell state first, then every categorical covariate with a workable number
        of levels. A covariate with one level per cell would draw one column per
        cell, so it is left out.
        """
        ds = store.dataset(job_id)
        cluster_key = ds.sv.get("cluster_key") or "cell_state"
        variables = [{"field": cluster_key, "values": list(ds.states), "n": len(ds.states)}]
        for name, info in (ds.covariate_names() or {}).items():
            if info.get("kind") != "categorical":
                continue
            levels = [str(v) for v in (info.get("categories") or [])]
            if 1 < len(levels) <= 60:
                variables.append({"field": name, "values": levels, "n": len(levels)})
        return {"cluster_key": cluster_key, "variables": variables}

    @app.get("/api/jobs/{job_id}/dotplot")
    def dotplot(job_id: str, genes: str = Query(""), group_by: str = Query(""),
                groups: List[str] = Query([]), subset_by: str = Query(""),
                subset_values: List[str] = Query([]), subset2_by: str = Query(""),
                subset2_values: List[str] = Query([])):
        """Mean and detected fraction per (gene, cell state), from the bundle's
        precomputed stats matrices. Default gene set = top marker of every state."""
        ds = store.dataset(job_id)
        group_column, group_names, group_code, is_states = _bundle_group_axis(ds, group_by)
        if groups:
            chosen = [g for g in group_names if g in set(groups)]
            if chosen:
                group_names = chosen
        pairs = _dotplot_default_genes(assets.get(job_id, {}), ds)
        requested = split_gene_list(genes)
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
        restrict = _bundle_subset_masks(ds, subset_by, list(subset_values),
                                        subset2_by, list(subset2_values))
        if restrict is None and is_states and len(group_names) == len(ds.states):
            # Fast path: the bundle already holds the per-state statistics.
            mean = np.asarray(ds.stats_mean[rows, :], dtype=np.float32).tolist()
            frac = np.asarray(ds.stats_frac[rows, :], dtype=np.float32).tolist()
            counts = list(ds.state_n)
        else:
            codes = np.asarray(group_code, dtype=np.int64)
            index = {g: i for i, g in enumerate(group_names)}
            original = _bundle_group_axis(ds, group_by)[1]
            remap = np.array([index.get(n, -1) for n in original], dtype=np.int64)
            per_cell = remap[codes]
            if restrict is not None:
                per_cell = np.where(restrict, per_cell, -1)
            counts = [int((per_cell == i).sum()) for i in range(len(group_names))]
            mean, frac = [], []
            n_cells = int(per_cell.shape[0])
            for row in rows:
                idx, val = ds.gene_column(row)
                sums = np.zeros(len(group_names), dtype=np.float64)
                hits = np.zeros(len(group_names), dtype=np.int64)
                if idx.size:
                    cells = idx.astype(np.int64)
                    keep_cells = per_cell[cells] >= 0
                    np.add.at(sums, per_cell[cells][keep_cells], val.astype(np.float64)[keep_cells])
                    np.add.at(hits, per_cell[cells][keep_cells], 1)
                mean.append([float(sums[i] / counts[i]) if counts[i] else 0.0
                             for i in range(len(group_names))])
                frac.append([float(hits[i] / counts[i]) if counts[i] else 0.0
                             for i in range(len(group_names))])
        return {
            "genes": labels, "states": group_names, "groups": group_names,
            "group_by": group_column, "group_label": group_column,
            "subset_by": subset_by or "", "subset_values": list(subset_values),
            "subset2_by": subset2_by or "", "subset2_values": list(subset2_values),
            "state_n": counts,
            "state_colors": [ds.colors.get(s, "#BBBBBB") for s in group_names],
            "colors": [ds.colors.get(s, "#BBBBBB") for s in group_names],
            "mean": mean, "frac": frac,
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

# EXACTLY ONE id, and it belongs to the DATASET, not to this module. There is
# deliberately no fallback to another study: showing a different study's title,
# description, tools, files and samples is presenting wrong data, whatever the notice
# says. When this record is absent from the site database the endpoint reads the SOURCE
# TABLE row instead (see read_study_record_from_tables), which is this study's own
# pending record, and reports source="source_tables".
#
# `study_ids_for_dataset` resolves the id per bundle. This constant is only the
# deployment-wide override, and it now defaults to EMPTY.
#
# Until 2026-08-12 it defaulted to `lmdata:LMEX0000004416`, and every bundle inherited
# that id. The COPD bundle then served the Study tab of `LMEX0000004416`, "Lipidomics
# Imaging of Human Postnatal Lung in Health and Bronchopulmonary Dysplasia", a different
# study with a different title, abstract, sample list and file list. The COPD record is
# `lmdata:LMEX0000009416`. A shared default cannot be right for more than one dataset, so
# there is no longer a default at all: a bundle states its own id, or the Study tab says
# that no id is configured and names the two flags that set one.
STUDY_CANDIDATES = [s.strip() for s in os.environ.get(
    "LUNGMAP_STUDY_IDS", "").split(",") if s.strip()]


def study_ids_for_dataset(assets: Dict[str, Dict[str, Any]], catalog: da.Catalog,
                          dataset_id: str) -> List[str]:
    """The LungMAP study id of one bundle, most specific source first.

    1. `study_id` in the dataset's asset manifest (`prepare_assets.py --study-id`)
    2. `scalable_viewer.study_id` in the bundle metadata (`precompute.py --study-id`)
    3. `LUNGMAP_STUDY_IDS`, the deployment-wide override
    4. nothing - the Study tab reports that no id is configured

    The asset manifest wins because rebuilding it costs minutes, while rebuilding a
    bundle costs hours. Step 4 never guesses: a wrong study record is wrong data.
    """
    entry = (assets or {}).get(dataset_id) or {}
    from_assets = str(entry.get("study_id") or "").strip()
    if from_assets:
        return [from_assets]
    try:
        ds = catalog.get(dataset_id)
        from_bundle = str((ds.sv or {}).get("study_id") or "").strip()
    except Exception:                                   # unknown dataset id, reported below
        from_bundle = ""
    if from_bundle:
        return [from_bundle]
    return list(STUDY_CANDIDATES)

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

    # A controlled-vocabulary term carries its own label. Deriving one by stripping the
    # prefix and swapping underscores invents a string the database does not use:
    # `data_type_single_nucleus_rna_seq` became "single nucleus rna seq" instead of the
    # real label "Single-nucleus RNA-seq". Read the label; never manufacture it.
    _VOCAB_FILES = {
        "data_type_": "lungmap_vocabulary/data_type.tsv",
        "sample_type_": "lungmap_vocabulary/sample_type.tsv",
        "ingest_stage_": "lungmap_vocabulary/ingest_stage.tsv",
        "age_range_": "lungmap_vocabulary/age_range.tsv",
    }

    def pretty(term):
        t = str(term or "").strip()
        if not t:
            return ""
        for prefix, rel in _VOCAB_FILES.items():
            if not t.startswith(prefix):
                continue
            path = os.path.join(tables_root, rel)
            if not os.path.exists(path):
                break
            with open(path, "r", encoding="utf-8") as fh:
                header = None
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if parts[0] == "#predicate":
                        header = parts
                    elif header and len(parts) > 1 and parts[1] == t:
                        if "rdfs:label" in header:
                            i = header.index("rdfs:label")
                            if i < len(parts) and parts[i].strip():
                                return parts[i].strip()
                        break
            break
        return t

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
    # A technology is named manufacturer + label, as the site does it:
    # "10x Genomics Chromium Single Cell Gene Expression".
    def _label_map(rel, pred):
        path, out_map, head = os.path.join(tables_root, rel), {}, None
        if not os.path.exists(path):
            return out_map
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if p[0] == "#predicate":
                    head = p
                elif head and p[0] == "data" and len(p) > 1 and pred in head:
                    i = head.index(pred)
                    if i < len(p):
                        out_map[p[1]] = p[i].strip()
        return out_map

    # manufacturer.tsv carries no rdfs:label; its name lives in lmdb:short_label
    # ("10x Genomics"), with lmdb:long_label as the fallback.
    makers = _label_map("lungmap_data/manufacturer.tsv", "lmdb:short_label")
    if not makers:
        makers = _label_map("lungmap_data/manufacturer.tsv", "lmdb:long_label")
    tech_maker = _label_map("lungmap_data/technology.tsv", "lmdb:has_manufacturer")
    tech_names = []
    for t in out["technologies"]:
        label = str(t.get("name") or "")
        maker = makers.get(tech_maker.get(str(t.get("target_id") or ""), ""), "")
        tech_names.append(f"{maker} {label}".strip()
                          if maker and not label.startswith(maker) else label)
    out["technologies"] = tech_names
    out["age_ranges"] = [_ucfirst(a["name"]) for a in out["age_ranges"]]

    # The analysis protocol as its own field, matching the database reader. A
    # protocol is the supporting_file whose file__file_type row names
    # `file_type_protocol`; nothing else in this schema links a study to a protocol
    # (15_add_analysis_protocol.py:5-16 measured that).
    protocol_files = set()
    ft_path = os.path.join(tables_root, "lungmap_xref", "file__file_type.tsv")
    if os.path.exists(ft_path):
        with open(ft_path, "r", encoding="utf-8") as fh:
            head_row = None
            for line in fh:
                cells = line.rstrip("\n").split("\t")
                if cells[0] == "#predicate":
                    head_row = cells
                elif head_row and cells[0] == "data" \
                        and "lmdb:applies_to_file" in head_row \
                        and "lmdb:has_file_type" in head_row:
                    file_col = head_row.index("lmdb:applies_to_file")
                    type_col = head_row.index("lmdb:has_file_type")
                    if max(file_col, type_col) < len(cells) \
                            and cells[type_col].strip() == "file_type_protocol":
                        protocol_files.add(cells[file_col].strip())
    protocol = None
    for entry in out["files"]:
        target = str(entry.get("target_id") or "").split(":")[-1]
        if target in protocol_files:
            protocol = {"name": str(entry.get("name") or "Analysis protocol"),
                        "url": str(entry.get("url") or "")}
            break

    out.update({
        "dataset_id": bare,
        "study_id": bare,
        "title": by_pred.get("rdfs:label", ""),
        "description": by_pred.get("rdfs:comment", ""),
        "assay": pretty(by_pred.get("lmdb:is_data_type")),
        "sample_type": pretty(by_pred.get("lmdb:has_sample_type")),
        # The COMMON name, which is what the site prints (pages.py:1052). The raw
        # ontology label of obo:NCBITaxon_9606 is "homo sapiens" and belongs in the
        # vocabulary, not on a study page.
        "organism": "Human" if "9606" in str(by_pred.get("lmdb:in_taxon", "")) else "",
        "cell_count": by_pred.get("lmdb:has_dataset_sample_count", ""),
        "ingest_stage": pretty(by_pred.get("lmdb:has_ingest_stage")),
        "release_date": by_pred.get("lmdb:has_release_date", ""),
        "reference": reference,
        "raw_data": raw_data,
        "protocol": protocol,
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
# Tool paths in the source tables are site-relative ("/breath-omics-analysis-page/?...")
# and belong to the LungMAP SITE, which is a DIFFERENT SERVER from this viewer. Serving
# them relative makes the browser resolve them against the viewer's own origin, so
# every tool link 404s here. They must carry the site's origin.
#
# That origin is the LOCAL redesign instance, not the public host: this machine holds
# the data, and a public hostname can change. Override with LUNGMAP_SITE_BASE only to
# point at a different deployment on purpose.
SITE_BASE = os.environ.get("LUNGMAP_SITE_BASE", "http://127.0.0.1:8001").rstrip("/")

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


def _ucfirst(text: str) -> str:
    """Capitalise the first letter and leave every other letter alone.

    The site prints a vocabulary label this way (app/lungmap/pages.py:1052 and the
    `(ucfirst)` filter named at pages.py:2459), so `adult` reads `Adult` there. The
    viewer printed the raw label, so the two surfaces disagreed on the word. This
    capitalises the vocabulary's OWN label. It never rebuilds a label from a slug.
    """
    text = str(text or "")
    return text[:1].upper() + text[1:] if text else text


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
            # A file URL gets the same SITE_BASE prefix the tool branch above
            # applies. Without it a site-relative display_url such as
            # /static/protocols/LMEX0000004416_analysis_protocol.html was handed
            # to the browser unchanged, so the viewer resolved it against its own
            # origin (port 8062) and every such link answered 404.
            path = entry.get("display_url") or ""
            entry["url"] = (SITE_BASE + path) if path.startswith("/") else path
            entry["size"] = entry.get("file_size") or ""
        # A publication row may carry its URL only inside its label, which is how the
        # COPD record was written. Pull the first http(s) token out rather than showing
        # a Reference with no link.
        for pub in out["publications"]:
            pub["url"] = pub.get("display_url") or pub.get("doi") or ""
            if not pub["url"]:
                hit = re.search(r"https?://\S+", str(pub.get("name") or ""))
                if hit:
                    pub["url"] = hit.group(0).rstrip(".,;)")
        # A technology is named manufacturer + product, the way the site names it:
        # "10x Genomics Chromium Single Cell Gene Expression". The site builds that
        # in pages.py:147 (DbAPITechnologies) by joining the technology to its
        # manufacturer and reading `lmdb:short_label`. manufacturer.tsv carries NO
        # rdfs:label, so `entity.label` is NULL for every manufacturer and reading it
        # returns nothing. The SOURCE-TABLE reader above was fixed on 2026-08-11 and
        # this one was not, and this study serves from the database, so the viewer
        # printed the bare product name beside the site's full one.
        tech_ids = [str(t.get("target_id") or "") for t in out["technologies"]]
        tech_ids = [t for t in tech_ids if t]
        makers: Dict[str, str] = {}
        if tech_ids:
            for tid, maker in cur.execute(
                "SELECT tech.id, COALESCE(short.object_literal, long.object_literal) "
                "FROM entity tech "
                "LEFT JOIN entity_property mlink ON mlink.subject_id = tech.id "
                "       AND mlink.predicate = 'lmdb:has_manufacturer' "
                "LEFT JOIN entity_property short ON short.subject_id = mlink.object_id "
                "       AND short.predicate = 'lmdb:short_label' "
                "LEFT JOIN entity_property long ON long.subject_id = mlink.object_id "
                "       AND long.predicate = 'lmdb:long_label' "
                "WHERE tech.id IN (%s)" % ",".join("?" * len(tech_ids)), tech_ids):
                makers[str(tid)] = str(maker or "").strip()
        technology_names = []
        for tech in out["technologies"]:
            label = str(tech.get("name") or "")
            maker = makers.get(str(tech.get("target_id") or ""), "")
            technology_names.append("%s %s" % (maker, label)
                                    if maker and not label.startswith(maker) else label)
        out["technologies"] = technology_names
        out["age_ranges"] = [_ucfirst(a["name"]) for a in out["age_ranges"]]

        # External accessions (GEO, dbGaP, PMID, DOI). The site resolves these through
        # a record entity whose `lmdb:has_resource_id` holds the accession and an
        # external_api whose `lmdb:has_id_type` names the record's class, then
        # substitutes [ID] into `lmdb:has_link` (see refactored_website
        # app/lungmap/pages.py:475). Reading only files and tools missed them entirely,
        # so Raw data rendered as a dash.
        out["db_xrefs"] = []
        for row in cur.execute(
            """
            SELECT COALESCE(rid.object_literal, REPLACE(x.object_id,'lmdata:','')) AS acc,
                   api.label AS api_label,
                   link.object_literal AS template
            FROM entity_property x
            JOIN entity target ON target.id = x.object_id
            LEFT JOIN entity_property rid
                   ON rid.subject_id = x.object_id AND rid.predicate = 'lmdb:has_resource_id'
            LEFT JOIN entity_property idt
                   ON idt.predicate = 'lmdb:has_id_type'
                  AND REPLACE(idt.object_id,'lmdata:','') = target.class
            LEFT JOIN entity api ON api.id = idt.subject_id
            LEFT JOIN entity_property link
                   ON link.subject_id = idt.subject_id AND link.predicate = 'lmdb:has_link'
            WHERE x.predicate IN ('oboInOwl:hasDbXref','lmdb:has_db_xref')
              AND x.subject_id IN (
                  SELECT subject_id FROM entity_property
                  WHERE predicate IN ('lmdb:part_of_experiment','lmdb:applies_to_dataset')
                    AND object_id = ?)
            GROUP BY x.object_id
            """, (study_id,)).fetchall():
            # This cursor yields plain tuples, not sqlite3.Row, so index positionally:
            # 0 = accession, 1 = api label, 2 = url template.
            acc = str(row[0] or "")
            api_label = str(row[1] or "")
            template = str(row[2] or "")
            url = template.replace("[ID]", acc) if template and acc else ""
            out["db_xrefs"].append({"name": ("%s %s" % (api_label, acc)).strip(),
                                    "accession": acc, "url": url})

        def pick(urls_from, hosts):
            for entry in urls_from:
                url = entry.get("url") or entry.get("path") or ""
                if any(host in url.lower() for host in hosts):
                    return {"name": entry.get("name", ""), "url": url}
            return None

        reference = None
        if out["publications"]:
            first_pub = out["publications"][0]
            # The citation and its URL are one string in some publication rows. Show
            # the citation as the link TEXT and keep the URL in the href; printing the
            # raw address inside the link text wraps over four lines and reads badly.
            pub_name = re.sub(r"\s*https?://\S+\s*$", "", str(first_pub["name"] or "")).strip()
            reference = {"name": pub_name or str(first_pub["name"] or ""),
                         "url": first_pub.get("url", "")}
        if not (reference and reference["url"]):
            hit = pick(out["files"], _REFERENCE_HOSTS) or pick(out["tools"], _REFERENCE_HOSTS)
            if hit:
                reference = hit if not reference else {"name": reference["name"], "url": hit["url"]}
        # Accessions first: a GEO series is the raw data, and the files/tools scan only
        # ever found it by accident.
        raw_data = (pick(out["db_xrefs"], _RAW_DATA_HOSTS)
                    or pick(out["files"], _RAW_DATA_HOSTS)
                    or pick(out["tools"], _RAW_DATA_HOSTS))
        if not (reference and reference.get("url")):
            hit = pick(out["db_xrefs"], _REFERENCE_HOSTS)
            if hit:
                reference = hit if not reference else {"name": reference["name"],
                                                       "url": hit["url"]}

        # Organism. `lmdb:in_taxon` carries the ontology term's own label, "homo
        # sapiens". The site shows the COMMON name: it reads `species_common_name`
        # out of ontology_property and capitalises it (pages.py:88 and pages.py:1052),
        # which gives "Human". The viewer printed the raw term label.
        taxon = _object_of(head, "lmdb:in_taxon")
        organism = ""
        if taxon:
            row = cur.execute("SELECT value FROM ontology_property WHERE term_id = ? "
                              "AND property = 'species_common_name'", (taxon,)).fetchone()
            if row and row[0]:
                organism = _ucfirst(str(row[0]))
        if not organism:
            organism = _ucfirst(str(_first(head, "lmdb:in_taxon") or ""))

        # The analysis protocol, as its OWN field. It is a supporting_file whose file
        # type is `protocol` (15_add_analysis_protocol.py registers it), so it reached
        # the payload only buried inside files[] and the Study tab had no Protocol row.
        # This is the query the site's Protocols tab runs, verbatim (pages.py:1315).
        protocol = None
        has_view = cur.execute("SELECT COUNT(*) FROM sqlite_master WHERE name = "
                               "'v_supporting_file'").fetchone()[0]
        if has_view:
            row = cur.execute(
                "SELECT f.label, f.display_url "
                "FROM entity_property x "
                "JOIN v_supporting_file f ON f.file_id = x.object_id "
                "WHERE x.predicate = 'lmdb:has_file' "
                "  AND x.subject_id IN (SELECT subject_id FROM entity_property "
                "      WHERE predicate IN ('lmdb:part_of_experiment', "
                "                          'lmdb:applies_to_dataset') "
                "        AND object_id = ?) "
                "  AND f.file_type_id = 'protocol' ORDER BY f.label", (study_id,)).fetchone()
            if row:
                path = str(row[1] or "")
                protocol = {"name": str(row[0] or "Analysis protocol"),
                            "url": (SITE_BASE + path) if path.startswith("/") else path}

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
            "organism": organism,
            "sample_type": _first(head, "lmdb:has_sample_type") or "",
            "cell_count": _first(head, "lmdb:has_dataset_sample_count") or "",
            "ingest_stage": _first(head, "lmdb:has_ingest_stage") or "",
            "reference": reference,
            "raw_data": raw_data,
            "protocol": protocol,
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
        site record has no sample rows yet, and whose own study id is used when the
        caller names no `study_id`.
        """
        assets = getattr(app.state, "assets", {}) or {}
        resolved_dataset = dataset or (catalog.entries[0]["id"] if catalog.entries else "")
        candidates = ([study_id] if study_id
                      else study_ids_for_dataset(assets, catalog, resolved_dataset))
        if not candidates:
            # Never answer with another study's record. Say what is missing instead.
            return JSONResponse({
                "ok": False,
                "dataset_bundle": resolved_dataset,
                "candidates": [],
                "error": (
                    f"no LungMAP study id is configured for the bundle {resolved_dataset!r}. "
                    f"Set one with `prepare_assets.py --study-id lmdata:LMEX...`, or with "
                    f"`precompute.py --study-id`, or with the LUNGMAP_STUDY_IDS environment "
                    f"variable. The Study tab shows no record rather than another study's."),
            }, status_code=200)
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

        ds_id = resolved_dataset
        record["dataset_bundle"] = ds_id
        record["study_id_source"] = ("query" if study_id else "dataset")
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


# --------------------------------------------------------------- protocol executors
#
# Five protocols that previously answered "not implemented". Each computes from
# the bundle: per-donor pseudobulk for the three that need a donor axis, cell
# counts for composition, and a covariate cross-tabulation for concordance.

def _run_severity_gradient(ds, state: str, covariate: str) -> Dict[str, Any]:
    """Genes whose per-donor level tracks a numeric clinical variable.

    A case-control contrast says whether a gene differs between two labelled
    groups. It cannot say whether the gene follows severity. Spearman against a
    measured variable separates a step at the diagnosis boundary from a decline
    with lung function, and only the second is a graded mechanism.
    """
    from scipy import stats as _st

    try:
        kind, values, _ = ds.covariate_values(covariate)
    except KeyError:
        return {"answer": f"{covariate} is not recorded in this dataset.",
                "status": "not_covered"}
    if kind != "numeric":
        return {"answer": f"{covariate} is categorical, so it has no gradient. "
                          "Ask for a group comparison instead.",
                "status": "not_covered"}

    rows = _candidate_rows(ds)
    matrix, donors, counts = _donor_pseudobulk(ds, rows, state)
    if matrix is None or len(donors) < 8:
        return {"answer": f"Too few donors contribute cells to {state} to correlate anything.",
                "status": "not_covered"}

    # One covariate value per donor: the mean over that donor's cells, which is
    # exact here because these covariates are donor-level and constant within a
    # donor.
    donor_code, donor_labels = _donor_axis(ds)
    per_donor = []
    values = np.asarray(values, dtype=np.float64)
    for name in donors:
        mask = donor_code == donor_labels.index(name)
        per_donor.append(float(np.nanmean(values[mask])) if mask.any() else np.nan)
    covariate_values = np.asarray(per_donor, dtype=np.float64)
    good = np.isfinite(covariate_values)
    if int(good.sum()) < 8:
        return {"answer": f"{covariate} is missing for too many donors in {state}.",
                "status": "not_covered"}

    # spearmanr returns a scalar for one gene, a matrix for many. Indexing the
    # scalar case as a matrix raised "IndexError: invalid index to scalar
    # variable" for EC general capillary, where only one gene survived the
    # per-donor filter.
    # spearmanr returns a scalar for a single pair and a correlation matrix for
    # many, so the RESULT decides how to read it, not the input shape. Indexing
    # a scalar as a matrix raised IndexError for EC general capillary, where the
    # per-donor filter left one gene.
    # Spearman is Pearson on ranks, so the whole scan is one matrix product.
    # `spearmanr` on the stacked matrix builds the full gene-by-gene correlation
    # matrix, 3022 x 3022 here, to use one column of it: that cost 5.07 s.
    scanned = matrix[:, good]
    n = int(good.sum())

    def rank_rows(block: np.ndarray) -> np.ndarray:
        """Average ranks along each row, ties shared, as Spearman requires."""
        order = np.argsort(block, axis=1, kind="mergesort")
        ranks = np.empty_like(order, dtype=np.float64)
        rows_index = np.arange(block.shape[0])[:, None]
        ranks[rows_index, order] = np.arange(block.shape[1], dtype=np.float64)
        # Ties: average the ranks of equal values, or a flat gene ranks 0..n-1
        # and correlates perfectly with anything.
        sorted_block = np.take_along_axis(block, order, axis=1)
        for row in range(block.shape[0]):
            values = sorted_block[row]
            start = 0
            for index in range(1, block.shape[1] + 1):
                if index == block.shape[1] or values[index] != values[start]:
                    if index - start > 1:
                        ranks[row, order[row, start:index]] = (start + index - 1) / 2.0
                    start = index
        return ranks

    gene_ranks = rank_rows(scanned)
    covariate_ranks = rank_rows(covariate_values[good][None, :])[0]
    gene_centred = gene_ranks - gene_ranks.mean(axis=1, keepdims=True)
    covariate_centred = covariate_ranks - covariate_ranks.mean()
    denominator = (np.linalg.norm(gene_centred, axis=1)
                   * np.linalg.norm(covariate_centred))
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = np.where(denominator > 0,
                       gene_centred @ covariate_centred / np.maximum(denominator, 1e-12),
                       0.0)
    rho = np.clip(np.nan_to_num(rho, nan=0.0), -1.0, 1.0)
    # The same t approximation scipy uses for Spearman.
    with np.errstate(invalid="ignore", divide="ignore"):
        tstat = rho * np.sqrt((n - 2) / np.maximum(1e-12, 1 - rho ** 2))
    pval = np.nan_to_num(2 * _st.t.sf(np.abs(tstat), max(1, n - 2)), nan=1.0)

    varying = scanned.std(axis=1) > 0
    rho = np.where(varying[:len(rho)], rho, 0.0)
    order = [int(i) for i in np.argsort(-np.abs(rho)) if varying[int(i)]][:25]
    # A linear fit alongside the rank correlation. Spearman says a relationship
    # is monotonic; the slope says how much the gene moves per unit of the
    # variable, which is what "regulated with age" actually claims. Both are
    # reported, because a high rho with a slope near zero is a real ordering of
    # noise.
    x = covariate_values[good]
    x_centred = x - x.mean()
    denominator = float((x_centred ** 2).sum())
    table_rows = []
    points: Dict[str, Any] = {}
    for i in order:
        y = scanned[int(i)]
        slope = float((x_centred * (y - y.mean())).sum() / denominator) if denominator else 0.0
        intercept = float(y.mean() - slope * x.mean())
        predicted = slope * x + intercept
        ss_res = float(((y - predicted) ** 2).sum())
        ss_tot = float(((y - y.mean()) ** 2).sum())
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        gene = ds.symbols[rows[int(i)]]
        table_rows.append([gene, round(float(rho[int(i)]), 4), float(pval[int(i)]),
                           round(slope, 6), round(r2, 4), int(good.sum()),
                           "up with " + covariate if slope > 0 else "down with " + covariate])
        if len(points) < 6:
            points[gene] = {"x": [round(float(v), 4) for v in x],
                            "y": [round(float(v), 5) for v in y],
                            "slope": round(slope, 6), "intercept": round(intercept, 6),
                            "rho": round(float(rho[int(i)]), 4), "r2": round(r2, 4)}
    if not table_rows:
        # Every scanned gene was flat across these donors, so no correlation is
        # defined. An empty table would read as "nothing correlates", which is a
        # different and unsupported claim.
        return {
            "answer": (f"No gene varies across the {int(good.sum())} donors that "
                       f"contribute cells to {state}, so no correlation with "
                       f"{covariate} is defined. Too few cells per donor in this "
                       "state, not an absence of signal."),
            "status": "not_covered",
            "table": {"columns": ["donors with cells in " + state],
                      "rows": [[d] for d in donors[:40]]},
        }
    return {
        "answer": (f"Genes in {state} whose per-donor level tracks {covariate}, "
                   f"Spearman across {int(good.sum())} donors, "
                   f"{len(rows)} genes scanned."),
        "status": "answered",
        "table": {"columns": ["gene", "rho", "p", "slope", "r2", "n_donors", "direction"],
                  "rows": table_rows},
        # One panel per gene: each donor's level against the variable, with the
        # fitted line. A CombPlot here drew all 39 cell states for a correlation
        # computed inside one of them, which showed the wrong thing entirely.
        "plot": {"kind": "gradient", "covariate": covariate, "state": state,
                 "donors": donors, "points": points},
    }


def _run_coexpression(ds, state: str, seed: str) -> Dict[str, Any]:
    """Genes co-varying with a seed gene across donors, inside one cell state."""
    seed_row = ds.resolve_gene(seed)
    if seed_row is None:
        return {"answer": f"{seed} is not in this dataset.", "status": "not_covered"}
    rows = _candidate_rows(ds)
    if seed_row not in rows:
        rows = [seed_row] + rows
    matrix, donors, _ = _donor_pseudobulk(ds, rows, state)
    if matrix is None or len(donors) < 8:
        return {"answer": f"Too few donors contribute cells to {state}.",
                "status": "not_covered"}

    seed_index = rows.index(seed_row)
    seed_values = matrix[seed_index]
    centred = matrix - matrix.mean(axis=1, keepdims=True)
    seed_centred = seed_values - seed_values.mean()
    denominator = (np.linalg.norm(centred, axis=1) * np.linalg.norm(seed_centred))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(denominator > 0, centred @ seed_centred / np.maximum(denominator, 1e-12), 0.0)
    r[seed_index] = 1.0
    order = [int(i) for i in np.argsort(-r) if int(i) != seed_index][:25]
    return {
        "answer": (f"Genes co-varying with {seed} across {len(donors)} donors in {state}. "
                   "Correlation is over per-donor pseudobulk, the same axis disease varies along."),
        "status": "answered",
        "table": {"columns": ["gene", "r with " + seed, "n_donors"],
                  "rows": [[ds.symbols[rows[i]], round(float(r[i]), 4), len(donors)]
                           for i in order]},
        "plot": {"kind": "combplot",
                 "genes": [seed] + [ds.symbols[rows[i]] for i in order[:5]]},
    }


def _categorical_twin(ds, covariate: str) -> str:
    """A categorical covariate standing for a numeric one, or "".

    Asking which states are depleted in GOLD IV routes the covariate slot to
    `gold_ordinal`, which is numeric, and composition needs groups. The same
    grading is also stored categorically in `Group`, whose levels are named
    GOLD I, II / GOLD III / GOLD IV. Matching on the words shared between the
    numeric covariate's name and a categorical covariate's levels finds it.
    """
    wanted = {w for w in bundle_meta._label_tokens(covariate) if len(w) > 2}
    if not wanted:
        return ""
    for name, info in (ds.covariate_names() or {}).items():
        if (info or {}).get("kind") != "categorical":
            continue
        try:
            kind, _v, labels = ds.covariate_values(name)
        except KeyError:
            continue
        if kind != "categorical" or labels is None or len(labels) < 2:
            continue
        words = {w for label in labels for w in bundle_meta._label_tokens(str(label))}
        words |= set(bundle_meta._label_tokens(name))
        if wanted & words:
            return name
    return ""


def _named_level(question: str, labels) -> int:
    """Index of the level the question names, or -1.

    `Group` carries three GOLD levels. A question about GOLD IV must compare
    GOLD IV, not whichever two levels happen to be encoded first.
    """
    asked = set(bundle_meta._label_tokens(question or ""))
    best, best_score = -1, 0
    for index, label in enumerate(labels):
        tokens = set(bundle_meta._label_tokens(str(label)))
        if tokens and tokens <= asked and len(tokens) > best_score:
            best, best_score = index, len(tokens)
    return best


def _run_composition(ds, covariate: str, question: str = "") -> Dict[str, Any]:
    """Which cell states change in abundance between the levels of a variable."""
    try:
        kind, values, labels = ds.covariate_values(covariate)
    except KeyError:
        return {"answer": f"{covariate} is not recorded here.", "status": "not_covered"}
    swapped = ""
    if kind != "categorical" or labels is None:
        twin = _categorical_twin(ds, covariate)
        if twin:
            swapped, covariate = covariate, twin
            kind, values, labels = ds.covariate_values(covariate)
    if kind != "categorical" or labels is None:
        return {"answer": f"{covariate} is numeric; composition needs groups.",
                "status": "not_covered"}

    donor_code, donor_labels = _donor_axis(ds)
    if donor_code is None:
        return {"answer": "This bundle records no donor column.", "status": "not_covered"}
    states = np.asarray(ds.state_code, dtype=np.int64)
    values = np.asarray(values, dtype=np.int64)

    # Fraction of each donor's cells in each state, then the mean per group.
    per_donor = np.zeros((len(donor_labels), len(ds.states)), dtype=np.float64)
    np.add.at(per_donor, (donor_code, states), 1.0)
    totals = per_donor.sum(axis=1, keepdims=True)
    fractions = np.divide(per_donor, np.maximum(totals, 1))
    donor_group = np.full(len(donor_labels), -1, dtype=np.int64)
    for donor in range(len(donor_labels)):
        mask = donor_code == donor
        if mask.any():
            donor_group[donor] = int(np.bincount(values[mask]).argmax())

    present = [g for g in range(len(labels)) if (donor_group == g).any()]
    if len(present) < 2:
        return {"answer": f"{covariate} has fewer than two groups with donors.",
                "status": "not_covered"}
    # Reference group first, so it takes the sky-blue bar and the disease the
    # red one. Stored order put COPD first and inverted the colours.
    named = _named_level(question, labels)
    if named in present and len(present) > 2:
        # Three GOLD levels: compare the one asked about against the mildest
        # remaining level rather than an arbitrary pair.
        rest = [g for g in present if g != named]
        first, second = _control_first(labels, [rest[0], named])
    else:
        first, second = _control_first(labels, present)
    a = fractions[donor_group == first].mean(axis=0)
    b = fractions[donor_group == second].mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.log2((a + 1e-6) / (b + 1e-6))
    # A Mann-Whitney per state: donor fractions are bounded and skewed, so a
    # rank test is the defensible one, and the two groups have unequal sizes.
    from scipy import stats as _st

    a_donors = fractions[donor_group == first]
    b_donors = fractions[donor_group == second]
    pvals = []
    for index in range(len(ds.states)):
        x, y = a_donors[:, index], b_donors[:, index]
        if x.size < 3 or y.size < 3 or (x.std() == 0 and y.std() == 0):
            pvals.append(1.0)
            continue
        try:
            pvals.append(float(_st.mannwhitneyu(x, y, alternative="two-sided").pvalue))
        except ValueError:
            pvals.append(1.0)
    # Rank by the test, not by the log ratio. The ratio carries a pseudocount,
    # so a state present at 0.00012 against 0.00000 topped the list on a
    # difference of one ten-thousandth of a donor's cells, which is noise. The
    # rank test asks whether the donors actually separate.
    order = sorted(range(len(ds.states)),
                   key=lambda i: (pvals[i], -abs(float(a[i] - b[i]))))[:25]
    return {
        "answer": ((f"{swapped} is numeric, so the grouped variable {covariate} "
                    f"was used instead. " if swapped else "")
                   + f"Cell-state abundance, {labels[first]} against {labels[second]}, "
                   f"as each donor's share of their own cells. "
                   f"{int((donor_group == first).sum())} and "
                   f"{int((donor_group == second).sum())} donors. "
                   "Abundance and expression are confounded: a depleted state looks "
                   "changed in any pooled expression comparison."),
        "status": "answered",
        "table": {"columns": ["cell state", f"mean fraction {labels[first]}",
                              f"mean fraction {labels[second]}", "log2 ratio", "p"],
                  "rows": [[ds.states[int(i)], round(float(a[int(i)]), 5),
                            round(float(b[int(i)]), 5), round(float(ratio[int(i)]), 4),
                            float(pvals[int(i)])]
                           for i in order]},
        # A paired bar per cell state, one bar per group, so the two frequencies
        # are read directly. A log ratio alone hides whether a state is common
        # in both groups or rare in both.
        # Every donor's own fraction travels with the summary, so the figure can
        # show the biological replicates behind each bar rather than only the
        # mean. A mean of 141 donors and a mean of 4 look identical otherwise.
        "plot": {"kind": "frequency",
                 "groups": [str(labels[first]), str(labels[second])],
                 "n_donors": [int((donor_group == first).sum()),
                              int((donor_group == second).sum())],
                 "states": [ds.states[int(i)] for i in order],
                 "a": [round(float(a[int(i)]), 6) for i in order],
                 "b": [round(float(b[int(i)]), 6) for i in order],
                 # Standard error of the mean, which is what an error bar on a
                 # group mean should show.
                 "a_sem": [round(float(a_donors[:, int(i)].std(ddof=1)
                                       / max(1.0, np.sqrt(a_donors.shape[0]))), 6)
                           for i in order],
                 "b_sem": [round(float(b_donors[:, int(i)].std(ddof=1)
                                       / max(1.0, np.sqrt(b_donors.shape[0]))), 6)
                           for i in order],
                 "a_points": [[round(float(v), 6) for v in a_donors[:, int(i)]]
                              for i in order],
                 "b_points": [[round(float(v), 6) for v in b_donors[:, int(i)]]
                              for i in order],
                 "p": [float(pvals[int(i)]) for i in order]},
    }


def _run_concordance(ds, state: str = "") -> Dict[str, Any]:
    """Where the two cell annotations agree, and where they do not."""
    second = next((name for name in ("TGEN-IPF", "Population", "celltype_level3")
                   if name in (ds.covariate_names() or {})), "")
    if not second:
        return {"answer": "This bundle carries only one cell annotation.",
                "status": "not_covered"}
    kind, values, labels = ds.covariate_values(second)
    if kind != "categorical" or labels is None:
        return {"answer": f"{second} is not categorical.", "status": "not_covered"}

    states = np.asarray(ds.state_code, dtype=np.int64)
    values = np.asarray(values, dtype=np.int64)
    wanted = range(len(ds.states)) if not state else [ds.states.index(state)]
    table_rows = []
    for index in wanted:
        mask = states == index
        total = int(mask.sum())
        if not total:
            continue
        counts = np.bincount(values[mask], minlength=len(labels))
        for other in np.argsort(-counts)[:3]:
            if counts[int(other)] == 0:
                continue
            table_rows.append([ds.states[int(index)], str(labels[int(other)]),
                               int(counts[int(other)]),
                               round(float(counts[int(other)]) / total, 4)])
    scope = state or f"all {len(ds.states)} cell states"
    return {
        "answer": (f"How {ds.sv.get('cluster_key') or 'cell_state'} maps onto {second}, "
                   f"for {scope}. A state whose top match holds well under 1.0 is a "
                   "boundary case, and a differential there may be an artefact of "
                   "which labelling was used."),
        "status": "answered",
        "table": {"columns": ["cell state", second, "n cells", "fraction"],
                  "rows": table_rows[:60]},
        "plot": {"kind": "heatmap"},
    }


def _run_dose_response(ds, state: str, covariate: str, genes: List[str]) -> Dict[str, Any]:
    """Whether a gene changes stepwise across the levels of an ordered variable."""
    from scipy import stats as _st

    try:
        kind, values, labels = ds.covariate_values(covariate)
    except KeyError:
        return {"answer": f"{covariate} is not recorded here.", "status": "not_covered"}
    if kind != "categorical" or labels is None:
        # "across GOLD stages" resolves to `gold_ordinal`, which is numeric.
        # The same variable exists categorically as `Group`, and an ordered
        # question wants the levels, so swap rather than refuse.
        swapped = ""
        for name, info in (ds.covariate_names() or {}).items():
            if info.get("kind") != "categorical":
                continue
            stem = covariate.lower().split("_")[0][:4]
            if stem and (stem in name.lower() or name.lower()[:4] == stem):
                swapped = name
                break
        if not swapped and "Group" in (ds.covariate_names() or {}):
            swapped = "Group"
        if not swapped:
            return {"answer": f"{covariate} is numeric; ask for a gradient instead.",
                    "status": "not_covered"}
        covariate = swapped
        kind, values, labels = ds.covariate_values(covariate)

    wanted = [g for g in (genes or []) if ds.resolve_gene(g) is not None]
    if not wanted:
        wanted = _bundle_default_genes(ds, [state])[:6]
    rows = [ds.resolve_gene(g) for g in wanted]
    matrix, donors, _ = _donor_pseudobulk(ds, rows, state)
    if matrix is None or len(donors) < 8:
        return {"answer": f"Too few donors contribute cells to {state}.",
                "status": "not_covered"}

    donor_code, donor_labels = _donor_axis(ds)
    values = np.asarray(values, dtype=np.int64)
    donor_level = []
    for name in donors:
        mask = donor_code == donor_labels.index(name)
        donor_level.append(int(np.bincount(values[mask]).argmax()) if mask.any() else -1)
    level = np.asarray(donor_level, dtype=np.float64)
    good = level >= 0

    table_rows = []
    for position, gene in enumerate(wanted):
        per_level = []
        for index in sorted(set(level[good].astype(int))):
            selected = good & (level == index)
            per_level.append(round(float(matrix[position][selected].mean()), 4)
                             if selected.any() else None)
        rho, pval = _st.spearmanr(level[good], matrix[position][good])
        clean = [v for v in per_level if v is not None]
        monotonic = (all(x <= y for x, y in zip(clean, clean[1:]))
                     or all(x >= y for x, y in zip(clean, clean[1:])))
        table_rows.append([gene, str(per_level), round(float(rho), 4),
                           float(pval), "yes" if monotonic else "no"])
    order = [str(l) for l in labels]
    return {
        "answer": (f"{', '.join(wanted)} across {covariate} in {state}, "
                   f"{int(good.sum())} donors. Levels in order: {', '.join(order)}. "
                   "A stepwise change behaves like a progression marker; a change only "
                   "at the extreme behaves like an end-stage marker."),
        "status": "answered",
        "table": {"columns": ["gene", "mean per level", "trend rho", "p", "monotonic"],
                  "rows": table_rows},
        "plot": {"kind": "combplot", "genes": wanted, "group_by": covariate},
    }


def _run_pathway_program(assets: Dict[str, Any], ds, state: str, contrast: str,
                         direction: str = "both", limit: int = 25) -> Dict[str, Any]:
    """Enriched processes for one cell state, read from the shipped GO-Elite table.

    A gene list is not a mechanism. GO-Elite is precomputed for each comparison,
    so the enriched terms are read directly rather than re-derived, and the
    genes driving each term come with them.

    The table keys its rows as `<cell state>__<direction>`, so a question about
    what is up in AT2 reads `AT2__up`.
    """
    import csv as _csv

    # The assets record the file directly; an earlier version guessed at the
    # directory layout with rglob and found nothing.
    head = (contrast or "").split("::")[0]
    differential = (assets or {}).get("differential") or {}
    entry = differential.get(contrast) or {}
    if not entry:
        entry = next((v for k, v in differential.items() if k.startswith(head)), {})
    path = str((entry or {}).get("goelite_tsv") or "")
    if not path or not Path(path).exists():
        return {"answer": f"No GO-Elite results are shipped for {head}.",
                "status": "not_covered"}
    candidates = [Path(path)]

    wanted = set()
    if direction in ("up", "both"):
        wanted.add(f"{state}__up")
    if direction in ("down", "both"):
        wanted.add(f"{state}__down")

    rows = []
    with candidates[0].open() as handle:
        for row in _csv.DictReader(handle, delimiter="\t"):
            if row.get("population") in wanted:
                try:
                    z = float(row.get("z_score") or 0)
                    fdr = float(row.get("fdr") or 1)
                except ValueError:
                    continue
                rows.append((z, fdr, row))
    if not rows:
        return {
            "answer": (f"GO-Elite is shipped for {head} but holds no terms for {state}. "
                       "It covers the cell states with enough differential genes to test."),
            "status": "not_covered",
        }

    # A term overlapping one gene can reach Z=27 and says nothing: the Z score
    # rewards a small expected count. Terms are ranked by Z but must overlap at
    # least two genes, and the count that was dropped is reported rather than
    # hidden, because a state whose only terms are one-gene terms has no
    # enrichment worth reading.
    single = [item for item in rows if int(float(item[2].get("overlap") or 0)) < 2]
    rows = [item for item in rows if int(float(item[2].get("overlap") or 0)) >= 2]
    if not rows:
        return {
            "answer": (f"GO-Elite found {len(single)} terms for {state} in {head}, and every "
                       "one overlaps a single gene. A one-gene term reaches a high Z score "
                       "because the expected count is tiny, so there is no enrichment here "
                       "worth reporting."),
            "status": "not_covered",
        }
    rows.sort(key=lambda item: -item[0])
    table_rows, chart = [], []
    for z, fdr, row in rows[:limit]:
        way = "up" if str(row.get("population", "")).endswith("__up") else "down"
        genes = str(row.get("overlap_genes") or "").replace("|", ", ")[:70]
        table_rows.append([row.get("term_name"), way, round(z, 3), fdr,
                           int(float(row.get("overlap") or 0)), genes])
        chart.append([f"{row.get('term_name')} ({way})", z if way == "up" else -z])

    return {
        "answer": (f"Enriched processes in {state} for {head}: {len(rows)} terms "
                   f"overlapping two or more genes, ranked by Z score. "
                   f"{len(single)} one-gene terms were left out, because a term with one "
                   "overlapping gene reaches a high Z on a tiny expected count. The genes "
                   "driving each term are listed beside it."),
        "status": "answered",
        "table": {"columns": ["term", "direction", "z", "fdr", "n overlap", "genes"],
                  "rows": table_rows},
        # Name the column explicitly: the last column here is the gene list, a
        # string, so a bar chart guessing "rightmost column" drew nothing.
        "plot": {"kind": "barchart", "value_column": "z",
                 "label_column": "term", "sign_column": "direction"},
        "chart": chart,
    }


#: Words that mark the reference side of a comparison. A level carrying one of
#: these is drawn first, in sky blue, so the disease bar sits beside it in red.
_REFERENCE_WORDS = {"non", "no", "not", "control", "controls", "ctrl", "healthy",
                    "normal", "never", "none", "negative", "unaffected", "ref"}


def _control_first(labels, present, reference: str = ""):
    """Order two levels so the reference group comes first.

    The bundle stores levels in whatever order the categories were encoded, and
    for `copd_status` that is COPD before non-COPD. Drawn as stored, the disease
    takes the control colour. `reference`, when given, is the right-hand side of
    a contrast id (`COPD_vs_non-COPD` -> `non-COPD`) and decides it outright;
    otherwise a level whose label carries a negation or health word is the
    reference. Neither test firing leaves the stored order alone.
    """
    first, second = present[0], present[1]
    words = lambda index: {w.lower() for w in
                           str(labels[index]).replace("-", " ").replace("_", " ").split()}
    if reference:
        wanted = {w.lower() for w in
                  reference.replace("-", " ").replace("_", " ").split() if w}
        score_first = len(wanted & words(first))
        score_second = len(wanted & words(second))
        if score_second > score_first:
            return second, first
        if score_first > score_second:
            return first, second
    if (words(second) & _REFERENCE_WORDS) and not (words(first) & _REFERENCE_WORDS):
        return second, first
    return first, second


def _frequency_payload(ds, covariate: str, states_in_order: Optional[List[str]] = None,
                       limit: int = 25, case_label: str = "",
                       control_label: str = "") -> Optional[Dict[str, Any]]:
    """Per-donor cell-state frequency for the two groups of a categorical variable.

    Shared by `composition_shift`, which ranks by the rank test, and by
    `most_affected_state`, which keeps its own transcriptional ranking and uses
    this only for the figure. Returns None when the variable has fewer than two
    groups with donors.
    """
    from scipy import stats as _st

    try:
        kind, values, labels = ds.covariate_values(covariate)
    except KeyError:
        return None
    if kind != "categorical" or labels is None:
        return None
    donor_code, donor_labels = _donor_axis(ds)
    if donor_code is None:
        return None

    states = np.asarray(ds.state_code, dtype=np.int64)
    values = np.asarray(values, dtype=np.int64)
    per_donor = np.zeros((len(donor_labels), len(ds.states)), dtype=np.float64)
    np.add.at(per_donor, (donor_code, states), 1.0)
    totals = per_donor.sum(axis=1, keepdims=True)
    fractions = np.divide(per_donor, np.maximum(totals, 1))

    donor_group = np.full(len(donor_labels), -1, dtype=np.int64)
    for donor in range(len(donor_labels)):
        mask = donor_code == donor
        if mask.any():
            donor_group[donor] = int(np.bincount(values[mask]).argmax())
    present = [g for g in range(len(labels)) if (donor_group == g).any()]
    if len(present) < 2:
        return None
    # When the caller knows which two levels the contrast compared, use exactly
    # those. `gold_ordinal` carries four levels, and taking the first two
    # present drew GOLD III for a GOLD IV contrast.
    first = second = None
    if case_label and control_label:
        by_key = {bundle_meta._label_tokens(str(labels[g])): g for g in present}
        first = by_key.get(bundle_meta._label_tokens(control_label))
        second = by_key.get(bundle_meta._label_tokens(case_label))
    if first is None or second is None or first == second:
        first, second = _control_first(labels, present)

    a_donors = fractions[donor_group == first]
    b_donors = fractions[donor_group == second]
    a, b = a_donors.mean(axis=0), b_donors.mean(axis=0)
    pvals = []
    for index in range(len(ds.states)):
        x, y = a_donors[:, index], b_donors[:, index]
        if x.size < 3 or y.size < 3 or (x.std() == 0 and y.std() == 0):
            pvals.append(1.0)
            continue
        try:
            pvals.append(float(_st.mannwhitneyu(x, y, alternative="two-sided").pvalue))
        except ValueError:
            pvals.append(1.0)

    if states_in_order:
        order = [ds.states.index(s) for s in states_in_order if s in ds.states][:limit]
    else:
        order = sorted(range(len(ds.states)),
                       key=lambda i: (pvals[i], -abs(float(a[i] - b[i]))))[:limit]

    sem = lambda block, i: float(block[:, i].std(ddof=1) / max(1.0, np.sqrt(block.shape[0])))
    return {
        "kind": "frequency",
        "groups": [str(labels[first]), str(labels[second])],
        "n_donors": [int((donor_group == first).sum()), int((donor_group == second).sum())],
        "states": [ds.states[int(i)] for i in order],
        "a": [round(float(a[int(i)]), 6) for i in order],
        "b": [round(float(b[int(i)]), 6) for i in order],
        "a_sem": [round(sem(a_donors, int(i)), 6) for i in order],
        "b_sem": [round(sem(b_donors, int(i)), 6) for i in order],
        "a_points": [[round(float(v), 6) for v in a_donors[:, int(i)]] for i in order],
        "b_points": [[round(float(v), 6) for v in b_donors[:, int(i)]] for i in order],
        "p": [float(pvals[int(i)]) for i in order],
        "_table": {"a": a, "b": b, "p": pvals, "order": order,
                   "labels": [str(labels[first]), str(labels[second])]},
    }


def _run_most_affected_state(ds, contrast: str, limit: int = 39) -> Dict[str, Any]:
    """Rank cell states by how strongly they respond to a contrast.

    The question is where the disease acts, not which genes move. Ranking the
    states by the number of genes passing FDR, and by the median effect among
    them, says which compartment to look at before any gene list is opened.

    This used to be aliased to the generic differential, so it answered with a
    gene volcano: the right numbers for a different question.
    """
    manifest = ds.deg_manifest() or {}
    comparisons = manifest.get("comparisons", [])
    chosen = next((c for c in comparisons if (c.get("id") or "") == contrast), None)
    if chosen is None:
        chosen = next((c for c in comparisons if "per_cell_state" in (c.get("id") or "")), None)
    if chosen is None:
        return {"answer": "This bundle carries no per-cell-state comparison.",
                "status": "not_covered"}

    table = ds.deg_table(chosen.get("id"), 1000000, 0.05, None) or {}
    rows = table.get("rows") if isinstance(table, dict) else table
    if not rows:
        return {"answer": f"{chosen.get('id')} returned no genes below FDR 0.05.",
                "status": "not_covered"}

    per_state: Dict[str, List[float]] = {}
    top_gene: Dict[str, tuple] = {}
    for row in rows:
        state = str(row.get("population") or row.get("cluster") or "")
        if not state:
            continue
        try:
            fold = float(row.get("log2fc") or 0.0)
        except (TypeError, ValueError):
            continue
        per_state.setdefault(state, []).append(abs(fold))
        best = top_gene.get(state)
        if best is None or abs(fold) > best[1]:
            top_gene[state] = (row.get("gene"), abs(fold), fold)

    ranked = sorted(per_state.items(), key=lambda kv: (-len(kv[1]), -float(np.median(kv[1]))))
    table_rows, chart = [], []
    for state, folds in ranked[:limit]:
        gene, _, signed = top_gene.get(state, ("", 0.0, 0.0))
        table_rows.append([state, len(folds), round(float(np.median(folds)), 4),
                           gene, round(float(signed), 4)])
        chart.append([state, len(folds)])

    covered = len(per_state)
    # The abundance figure must split donors on the same variable the
    # differential compared, so it is derived from the contrast id.
    frequency = None
    try:
        categorical = bundle_meta._categorical_covariates(ds)
        grouping, case_value, control_value = bundle_meta._contrast_group_field(
            ds, chosen, categorical)
    except Exception:                                    # noqa: BLE001
        grouping, case_value, control_value = "", "", ""
    if grouping:
        frequency = _frequency_payload(ds, grouping, [r[0] for r in table_rows],
                                       case_label=case_value,
                                       control_label=control_value)
        if frequency:
            frequency.pop("_table", None)
    return {
        "answer": (f"Cell states ranked by how strongly they respond to "
                   f"{chosen.get('id')}, counting genes below FDR 0.05. "
                   f"{covered} of the {len(ds.states)} states have any, and the "
                   "contrast is only computed for the states with enough cells. "
                   + (f"The table ranks by differential genes. The bars show the "
                      f"same states' abundance in each group, "
                      f"{frequency['groups'][0]} against {frequency['groups'][1]}, "
                      f"as each donor's share of their own cells, with one dot per "
                      f"donor and standard error. A state can change in expression, "
                      f"in abundance, or in both, and the two are confounded: a "
                      f"depleted state looks changed in any pooled comparison."
                      if frequency else
                      "Where the disease acts, before any gene list is opened.")),
        "status": "answered",
        "table": {"columns": ["cell state", "n significant", "median |log2fc|",
                              "top gene", "its log2fc"],
                  "rows": table_rows},
        # The question spans two things: which states change transcriptionally,
        # and which change in abundance. The table ranks by differential genes;
        # the figure shows each group's cell frequency for those same states, in
        # the same order, so both are visible at once.
        "plot": frequency or {"kind": "barchart", "value_column": "n significant",
                              "label_column": "cell state"},
    }


def _run_donor_heterogeneity(ds, state: str, contrast: str,
                             n_genes: int = 40) -> Dict[str, Any]:
    """Whether a signature is carried by most donors or by a few.

    A significant fold change is a statement about group means. With 178 donors
    a gene can clear FDR because a handful carry it strongly while the rest show
    nothing, and that is a subtype rather than a disease mechanism.

    So the signature is scored per donor: the genes up in the contrast minus the
    genes down, each standardised across donors, inside the one cell state. A
    heatmap of those genes with donors ordered by their score shows at a glance
    whether the case donors form a block or scatter through the controls.

    This used to answer with a gene volcano, which shows the group means the
    question is asking to look past.
    """
    manifest = ds.deg_manifest() or {}
    comparisons = manifest.get("comparisons", [])
    chosen = next((c for c in comparisons if (c.get("id") or "") == contrast), None)
    if chosen is None:
        chosen = next((c for c in comparisons if "per_cell_state" in (c.get("id") or "")), None)
    if chosen is None:
        return {"answer": "This bundle carries no per-cell-state comparison.",
                "status": "not_covered"}

    table = ds.deg_table(chosen.get("id"), 100000, 0.05, state) or {}
    rows = (table.get("rows") if isinstance(table, dict) else table) or []
    if not rows:
        return {"answer": f"{chosen.get('id')} is not computed for {state}.",
                "status": "not_covered"}

    ups, downs = [], []
    for row in sorted(rows, key=lambda r: -abs(float(r.get("log2fc") or 0))):
        gene = row.get("gene")
        row_index = ds.resolve_gene(str(gene))
        if row_index is None:
            continue
        if float(row.get("log2fc") or 0) > 0 and len(ups) < n_genes // 2:
            ups.append((gene, row_index))
        elif float(row.get("log2fc") or 0) < 0 and len(downs) < n_genes // 2:
            downs.append((gene, row_index))
        if len(ups) >= n_genes // 2 and len(downs) >= n_genes // 2:
            break
    signature = ups + downs
    if len(signature) < 4:
        return {"answer": f"Too few significant genes in {state} to score a signature.",
                "status": "not_covered"}

    indices = [i for _, i in signature]
    matrix, donors, counts = _donor_pseudobulk(ds, indices, state)
    if matrix is None or len(donors) < 8:
        return {"answer": f"Too few donors contribute cells to {state}.",
                "status": "not_covered"}

    # Standardise each gene across donors so one loud gene cannot set the score.
    centre = matrix.mean(axis=1, keepdims=True)
    spread = matrix.std(axis=1, keepdims=True)
    z = np.divide(matrix - centre, np.where(spread > 0, spread, 1.0))
    up_rows = np.arange(len(ups))
    down_rows = np.arange(len(ups), len(signature))
    score = z[up_rows].mean(axis=0) - (z[down_rows].mean(axis=0) if down_rows.size else 0.0)

    donor_code, donor_labels = _donor_axis(ds)
    group_of = {}
    case_name = (chosen.get("id") or "").split("::")[0].split("_vs_")[0].replace("_", " ")
    for name in ("copd_status", "dx_category", "Group"):
        if name in (ds.covariate_names() or {}):
            kind, values, labels = ds.covariate_values(name)
            if kind == "categorical" and labels is not None:
                values = np.asarray(values, dtype=np.int64)
                for donor in donors:
                    mask = donor_code == donor_labels.index(donor)
                    if mask.any():
                        group_of[donor] = str(labels[int(np.bincount(values[mask]).argmax())])
            break

    order = np.argsort(-score)
    ordered = [donors[int(i)] for i in order]
    groups = [group_of.get(d, "") for d in ordered]
    levels = sorted({g for g in groups if g})
    # How far up the ranking the case donors sit: if the signature is universal
    # they are spread through it, if it is a subtype they cluster at the top.
    summary = ""
    if len(levels) == 2:
        case = next((l for l in levels if case_name.lower().split()[0] in l.lower()), levels[0])
        positions = [i for i, g in enumerate(groups) if g == case]
        n_case = len(positions)
        top_half = sum(1 for i in positions if i < len(ordered) / 2)
        summary = (f" {top_half} of the {n_case} {case} donors fall in the top half of the "
                   f"ranking; an even split would put about {n_case // 2} there.")

    return {
        "answer": (f"The {chosen.get('id')} signature in {state}, scored per donor: "
                   f"{len(ups)} up genes minus {len(downs)} down genes, each standardised "
                   f"across the {len(donors)} donors with cells in this state.{summary} "
                   "A fold change is a group mean; this is who carries it."),
        "status": "answered",
        "table": {"columns": ["donor", "signature score", "group", "n cells"],
                  "rows": [[donors[int(i)], round(float(score[int(i)]), 4),
                            group_of.get(donors[int(i)], ""), int(counts[int(i)])]
                           for i in order]},
        "plot": {"kind": "signature",
                 "genes": [g for g, _ in signature],
                 "n_up": len(ups),
                 "donors": ordered,
                 "groups": groups,
                 "score": [round(float(score[int(i)]), 4) for i in order],
                 "z": [[round(float(z[r][int(c)]), 3) for c in order]
                       for r in range(len(signature))]},
    }
