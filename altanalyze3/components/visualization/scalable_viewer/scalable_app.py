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
from . import grn_network as gnet
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
        # prepare_assets writes the manifest FLAT, as <assets_root>/<prefix>_assets.json.
        # This looked only in <assets_root>/<id>/, found nothing, and returned an empty
        # dict, so every asset-driven feature reported itself unavailable: the marker
        # heatmap, marker networks, fastComm, GRN networks, and the GO-Elite and network
        # panels of the Differential tab. The contrast dropdown kept working because
        # precompute embeds the DEG tables in the bundle rather than the manifest, which
        # is why the failure looked partial. Both layouts are accepted now.
        candidates = [
            Path(assets_root) / f"{entry['prefix']}_assets.json",
            Path(assets_root) / entry["id"] / f"{entry['prefix']}_assets.json",
        ]
        manifest = next((c for c in candidates if c.is_file()), None)
        if manifest is None:
            continue
        with open(manifest, "r") as fh:
            data = json.load(fh)
        out[entry["id"]] = _absolutise_asset_paths(data)
    return out


def _absolutise_asset_paths(data: Dict[str, Any]) -> Dict[str, Any]:
    """Resolve every relative path in a manifest against the project root.

    prepare_assets writes paths relative to the project it ran in, but the server
    runs from wherever launchd starts it. `Path(rel).exists()` then answers False
    and the feature reports itself unavailable: fastComm returned "fastComm scores
    are unavailable", the marker heatmap returned "Marker heatmap matrix
    unavailable", and marker networks returned an empty element list. Every one of
    those files existed.

    `bundle_dir` is absolute, and every other path in the manifest is relative to
    the project root three levels above it (<root>/scalable_viewer/bundles_*/<id>).
    A path that already resolves is left alone, so nothing is rewritten twice.
    """
    bundle_dir = str(data.get("bundle_dir") or "").strip()
    if not bundle_dir or not os.path.isabs(bundle_dir):
        return data
    root = os.path.abspath(os.path.join(bundle_dir, "..", "..", ".."))

    def fix(value):
        if isinstance(value, str):
            if not value or os.path.isabs(value) or os.path.exists(value):
                return value
            candidate = os.path.join(root, value)
            return candidate if os.path.exists(candidate) else value
        if isinstance(value, dict):
            return {k: fix(v) for k, v in value.items()}
        if isinstance(value, list):
            return [fix(v) for v in value]
        return value

    return fix(data)


def create_scalable_app(
    catalog: da.Catalog,
    *,
    state_dir: str,
    assets_root: Optional[str] = None,
):
    assets = _load_assets(assets_root, catalog)

    # scALABLE's own template directory and static directory. Nothing is copied: the
    # viewer serves app.js, styles.css and index.html straight out of the webapp.
    # SCALABLE_ROOT_PATH is the path prefix the viewer answers on, empty at an origin
    # of its own and "/scalable-viewer" behind the LungMAP proxy. The web app already
    # carries the prefix through: FastAPI routes under it, the template publishes it as
    # `__APP_ROOT_PATH__`, and app.js `apiPath()` prefixes every call. The proxy must
    # pass the prefix on rather than strip it, as it does for /cellharmony/.
    app = W.create_app({
        "JOB_STORAGE": state_dir,
        "TEMPLATE_DIR": str(WEBAPP_DIR / "templates"),
        "STATIC_DIR": str(WEBAPP_DIR / "static"),
        "INDEX_TEMPLATE": "index.html",
        "ROOT_PATH": os.environ.get("SCALABLE_ROOT_PATH", ""),
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

    # Static JS and CSS went out with an ETag and NO Cache-Control. A browser
    # given that combination applies HEURISTIC caching: it may reuse the file for
    # a fraction of its age without ever asking the server. Nathan reported the
    # same two viewer bugs three times while the server was already serving the
    # fixed files, because his browser kept replaying old JavaScript.
    # `no-cache` does not disable caching. It requires REVALIDATION, so the ETag
    # above decides, and a 304 still costs nothing when nothing changed.
    @app.middleware("http")
    async def revalidate_static(request, call_next):
        response = await call_next(request)
        path = request.url.path
        if path.endswith((".js", ".css", ".map")) or "/static" in path:
            response.headers["Cache-Control"] = "no-cache, must-revalidate"
        return response

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
    _install_feature_name_routes(app, store)
    _install_differential_feature_name_route(app, store)
    _install_grn_routes(app, store)
    _install_grn_edge_adata(app, store)
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
        if not len(finite_values):
            continue
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
        x_field: str = Query(""),
        y_field: str = Query(""),
    ):
        token = _VIOLIN_COVARIATE.set(str(covariate or "").strip())
        try:
            response = await json_endpoint(
                job_id,
                gene=_to_feature_key(app, job_id, modality, gene),
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
                # The expression UMAP takes the same coordinate pair as every
                # other UMAP panel. Omitting them here would silently drop the
                # reader's choice, the way violin_limit once was.
                x_field=x_field,
                y_field=y_field,
            )
            # THE TITLE NAMES THE FEATURE THE WAY THE READER DOES.
            #
            # `_build_expression_payload` sets `gene` to the resolved KEY, so a lipid plot
            # was titled "CE(18:2) expression" even when the reader picked
            # "18:2 Cholesterol ester". `resolved_gene` still carries the key, so nothing
            # that needs the key loses it. A store with no display column returns the key
            # from `display_of`, which makes this a no-op everywhere else.
            return _retitle_with_display_name(app, job_id, modality, response)
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
                gene=_to_feature_key(app, job_id, modality, gene),
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
        # Under a prefix the shell arrives as "/scalable-viewer/", not "/": routing
        # strips the root path, `request.url.path` keeps it. Comparing the raw path
        # left the page un-renamed and the bootstrap script off the only page that
        # loads it.
        root = W._normalize_root_path(request.scope.get("root_path") or "")
        path = request.url.path
        if root and path.startswith(root):
            path = path[len(root):] or "/"
        if path != "/" or response.status_code != 200:
            return response
        chunks = [chunk async for chunk in response.body_iterator]
        body = b"".join(chunks).decode("utf-8")
        body = body.replace("<title>scALABLE</title>", f"<title>{VIEWER_NAME}</title>")
        body = body.replace("<h1>scALABLE</h1>", f"<h1>{VIEWER_NAME}</h1>")
        for old, new in DOC_LINKS.items():
            body = body.replace(old, new)
        version = int(os.path.getmtime(bootstrap))
        tag = (f'<link rel="stylesheet" href="{root}/viewer-static/viewer.css?v={version}">'
               f'<script src="{root}/viewer-static/viewer_bootstrap.js?v={version}"></script></body>')
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

    _drop_official_route(app, "/api/jobs/{job_id}/differential/select", "POST")

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
        # THE NAMES FOLLOW THE CONTRAST. Updating only `differential` left
        # `feature_display` as it was built on first load, which is empty for RNA, so
        # switching to the lipid differential kept showing `PE(16:0/22:4)`.
        store.update_job(job_id, differential=block,
                         feature_display=bundle_meta.feature_display_map(
                             ds, str(chosen.get("modality") or "rna")))
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


def split_gene_list(text: str, keep_pipe: bool = False, known=None) -> List[str]:
    """Gene symbols out of whatever a user pasted.

    A user pastes a column out of Excel, which arrives newline separated, or a
    row, which arrives tab separated. Others type commas, spaces or semicolons.
    Splitting on commas alone accepted only one of those: `SFTPC AGER` resolved
    to a single symbol named "SFTPC AGER" and the request 404ed.

    Duplicates are dropped and order is kept, so the plot reads in the order the
    user listed the genes.

    `keep_pipe` stops the vertical bar being a separator. A GRN edge is named
    `TF|target`, so splitting on the bar turned `TEAD3|MRPL33` into two names that
    the store does not hold, and the request 404ed.
    """
    pattern = r"[\s,;]+" if keep_pipe else r"[\s,;|]+"
    seen, out = set(), []

    def keep(raw: str) -> None:
        gene = raw.strip().strip('"').strip("'")
        if gene and gene not in seen:
            seen.add(gene)
            out.append(gene)

    # A FEATURE NAME MAY CONTAIN A SPACE, SO A CHUNK IS TRIED WHOLE BEFORE SPLITTING.
    #
    # Nathan, 2026-09-08: lipids display as "18:2 Cholesterol ester". Splitting on
    # whitespace first turned that one name into three tokens and the request 404ed with
    # "none of the 3 requested genes are in this dataset". Splitting on commas alone is
    # not the fix either: the docstring above records that `SFTPC AGER` must still resolve
    # to two genes.
    #
    # So the separators a name can never contain -- newline, tab, semicolon, comma, and
    # the bar unless keep_pipe -- cut the text into chunks, and each chunk is offered to
    # `known` intact. A chunk that names a real feature is kept whole; one that does not
    # falls back to whitespace splitting, which is the old behaviour exactly. Without
    # `known` nothing changes, so every existing caller is unaffected.
    hard = r"[\n\r\t;,]+" if keep_pipe else r"[\n\r\t;,|]+"
    for chunk in re.split(hard, str(text or "")):
        chunk = chunk.strip().strip('"').strip("'")
        if not chunk:
            continue
        if known is not None and " " in chunk:
            try:
                if known(chunk) is not None:
                    keep(chunk)
                    continue
            except Exception:                      # a resolver that raises is not a match
                pass
        for part in re.split(pattern, chunk):
            keep(part)
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
        return np.isin(codes, list(wanted))
    try:
        kind, values, labels = ds.covariate_values(subset_by)
    except KeyError:
        return None
    if kind != "categorical" or labels is None:
        return None
    wanted = {i for i, s in enumerate(labels) if str(s) in keep}
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


def _feature_source(ds, modality: str = ""):
    """The store the DotPlot and CombPlot read: the RNA store, or one modality's.

    A `da.ModalityStore` answers `resolve_gene`, `gene_column`, `symbols`, `stats_mean`
    and `stats_frac` exactly as the Dataset does, so the two plots need no other change.
    """
    name = str(modality or "").strip().lower()
    if not name or name == "rna":
        return ds
    try:
        return ds.modality(name)
    except (KeyError, FileNotFoundError) as exc:
        raise HTTPException(404, f"this bundle carries no modality '{name}': {exc}")


def _bundle_default_genes(ds, groups, source=None) -> List[str]:
    """One marker gene per group, for a blank gene set.

    The bundle ships a marker table keyed by cell state, so that is used when the
    grouping is cell state. Any other variable has no such table, and the top
    markers of the states are still the informative genes to open on.
    """
    # An imputed modality ships no marker table, so its opening set is its own first
    # features, in the order the prediction table named them.
    if source is not None and source is not ds:
        return list(source.symbols[:12])
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
                 min_cells: int = Query(5), cells_per_sample: int = Query(10), group_by: str = Query(""),
                 groups: List[str] = Query([]), subset_by: str = Query(""),
                 subset_values: List[str] = Query([]), subset2_by: str = Query(""),
                 subset2_values: List[str] = Query([]), tracks: str = Query(""),
                 modality: str = Query(""), unit: str = Query("cells", pattern="^(cells|donor)$")):
        """Individual stored observations by default, or explicit per-donor means."""
        ds = store.dataset(job_id)
        group_column, group_names, group_code, _ = _bundle_group_axis(ds, group_by)
        if groups:
            wanted_groups = [g for g in group_names if g in set(groups)]
            group_names = wanted_groups
        # Blank means the marker gene of every group, matching the DotPlot.
        features = _feature_source(ds, modality)
        wanted = (split_gene_list(genes, keep_pipe=features.names_have_pipe)
                  or _bundle_default_genes(ds, group_names, features))
        if not wanted:
            raise HTTPException(400, "give at least one gene")

        if unit == "cells":
            codes = np.asarray(group_code, dtype=np.int64)
            valid = (codes >= 0) & (codes < len(_bundle_group_axis(ds, group_by)[1]))
            original = _bundle_group_axis(ds, group_by)[1]
            selected_codes = [i for i, name in enumerate(original) if name in group_names]
            valid &= np.isin(codes, selected_codes)
            restrict = _bundle_subset_masks(ds, subset_by, list(subset_values),
                                            subset2_by, list(subset2_values))
            if restrict is not None:
                valid &= restrict
            keep = np.flatnonzero(valid)
            keep = keep[np.argsort(codes[keep], kind="stable")]
            names = bundle_meta._obs_names(ds)
            from types import SimpleNamespace
            sample_field = _donor_covariate(ds) or ""
            obs = pd.DataFrame(index=names)
            if sample_field:
                _, sample_codes, sample_names = ds.covariate_values(sample_field)
                obs[sample_field] = [str(sample_names[int(c)]) if int(c) >= 0 else "" for c in sample_codes]
            sampling_cache = {"adata": SimpleNamespace(obs=obs), "obs_names": names,
                              "sample_field": sample_field, "populations": np.asarray(ds.states)[np.asarray(ds.state_code)]}
            keep, sampling = W._sample_plot_cells(sampling_cache, keep, cells_per_sample)
            columns = [{"cell": str(names[i]), "group": original[codes[i]],
                        "state": original[codes[i]], "n_cells": 1} for i in keep]
            labels, series, missing = [], [], []
            for gene in wanted:
                row = features.resolve_gene(gene)
                if row is None:
                    missing.append(gene)
                    continue
                idx, val = features.gene_column(row)
                values = np.zeros(ds.n_cells, dtype=float)
                values[np.asarray(idx, dtype=np.int64)] = val
                labels.append(gene)
                series.append(np.round(values[keep], 5).tolist())
            if not labels:
                raise HTTPException(404, "none of the requested genes are in this dataset")
            # Each track column is now an observation, so its annotation is exact.
            observation_codes = np.arange(ds.n_cells, dtype=np.int64)
            asked = [t.strip() for t in tracks.split(",") if t.strip()]
            if not asked:
                asked = [t for t in COMBPLOT_TRACK_DEFAULTS if t in (ds.covariate_names() or {})]
            track_names, track_values, track_levels, purity, skipped = _combplot_tracks(
                ds, asked, observation_codes, keep, ds.n_cells)
            return {"unit": "cells", "sampling": sampling, "observation_unit": ds.sv.get("observation_unit", "cells"),
                    "genes": labels, "values": series, "columns": columns,
                    "colors": [ds.colors.get(c["group"], "#BBBBBB") for c in columns],
                    "groups": group_names, "states": group_names,
                    "group_by": group_column, "group_label": group_column,
                    "n_columns": len(columns), "n_cells_kept": len(columns), "n_groups_dropped": 0,
                    "track_names": track_names, "tracks": track_values, "track_levels": track_levels,
                    "track_purity": purity, "track_skipped": skipped,
                    "n_requested": len(wanted), "n_returned": len(labels),
                    "n_missing": len(missing), "missing": missing}

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
            row = features.resolve_gene(gene)
            if row is None:
                missing.append(gene)
                continue
            idx, val = features.gene_column(row)
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
            "unit": "donor",
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

        from altanalyze3.components.cellHarmony.webapp.cross_modal import examples as cross_examples
        available = (store.get_job(job_id).get("modalities") or {}).get("available", [])
        examples = cross_examples([m["id"] for m in available]) + [f"What pathways have the best cross-modality {first} representation?", f"What is the best modality marker of {first}?"] + examples
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

        from altanalyze3.components.cellHarmony.webapp.cross_pathways import answer_if_requested as pathway_answer
        pathways = pathway_answer(app, store.get_job(job_id), question)
        if pathways is not None:
            return pathways
        from altanalyze3.components.cellHarmony.webapp.modality_markers import answer_if_requested as marker_answer
        markers = marker_answer(app, store.get_job(job_id), question)
        if markers is not None:
            return markers
        from altanalyze3.components.cellHarmony.webapp.cross_modal import answer_if_requested
        cross_answer = answer_if_requested(app, store.get_job(job_id), question)
        if cross_answer is not None:
            return cross_answer
        from altanalyze3.components.cellHarmony.webapp.integration_chat import read_question as read_integrated, answer as integrated_answer
        integrated = read_integrated(question, ds.states, ds.symbols)
        if integrated is not None:
            return integrated_answer(app, store.get_job(job_id), question, integrated)
        from altanalyze3.components.cellHarmony.webapp.chat_service import read_question as read_protocol, execute as execute_protocol
        local=read_protocol(question,ds.states,ds.symbols,getattr(ds,"covariate_names",lambda:{})())
        if local is not None:
            answer=execute_protocol(app,store.get_job(job_id),question,local)
            if answer is not None:return answer
        reading = gnet.read_regulatory_question(question, ds.states, ds.symbols)
        if reading is None and local is not None:
            reading=local
            from altanalyze3.components.cellHarmony.webapp.chat_service import resolve_contrast
            reading['contrast']=resolve_contrast(ds,store.get_job(job_id),question)

        if reading is not None and reading.get("router") != "local_protocol":
            selected = (store.get_job(job_id).get("differential") or {}).get("run_id", "")
            modality = reading["modality"]
            candidates = [c for c in ds.deg_manifest().get("comparisons", [])
                          if c.get("modality", "rna") == modality and c.get("kind") == "per_cell_state"]
            question_words = re.sub(r"[^a-z0-9]+", " ", question.lower())
            named = next((c["id"] for c in candidates if c.get("comparison") and
                          re.sub(r"[^a-z0-9]+", " ", c.get("comparison", "").lower()) in question_words), "")
            reading["contrast"] = named or gnet._sibling_comparison(ds, selected, modality)
            if not reading["contrast"] and reading["intent"] != "differential":
                reading["contrast"] = selected
        elif reading is None:
            try:
                import urllib.request
                covariates = [name for name, info in (ds.covariate_names() or {}).items()
                              if info.get("kind") in ("numeric", "categorical")]
                # THE ROUTER MUST BE TOLD EVERY MODALITY THIS BUNDLE CARRIES.
                #
                # `ds.sv` is the bundle's STORED metadata block, where `modalities` is the raw
                # manifest keyed by id -- {"adt": {...}, "grn": {...}, "grn_tf": {...}} -- and
                # carries no "available" key. So `.get("available")` was always None and this
                # sent ["rna"] for every dataset. The router then had no modality slot to fill,
                # and no GRN, TF-activity, ADT, lipid or cell-communication question could
                # route at all. Measured on the COPD-metacells bundle 2026-09-08: the runtime
                # block holds ['rna','adt','grn','grn_tf','lipid'] and this line sent ['rna'].
                #
                # `bundle_meta._modalities_block(ds)` is the runtime shape, {"default", ...,
                # "available": [...]}, and it is what /api/catalog already serves.
                modalities = [m.get("id") for m in
                              (bundle_meta._modalities_block(ds).get("available") or [])] or ["rna"]
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

        protocol_result=execute_protocol(app,store.get_job(job_id),question,reading)
        if protocol_result is not None:return protocol_result
        intent = reading.get("intent")
        state = reading.get("cell_state") or ""
        state2 = reading.get("cell_state_2") or ""
        genes = reading.get("genes") or []
        contrast = reading.get("contrast") or ""
        covariate = reading.get("covariate") or ""

        result: Dict[str, Any] = {"question": question, "reading": reading,
                                  "intent": intent}

        if reading.get("router") == "local_regulatory" and intent == "differential" and not contrast:
            result.update(status="not_run", answer="No precomputed comparison for this modality matches the selected groups.")
            return result

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
        # --- the two regulatory answers, computed here -----------------------
        #
        # Nathan, 2026-09-07: "Both GRN and TF-activity should be able to have separate
        # types of plots supported in Chat." Before this, `regulatory_driver` fell into the
        # branch below and told the reader to open the differential-network tab, which
        # answers a different question with a different feature space.
        #
        # WHICH OF THE TWO. A question about a factor's own activity wants the ranked
        # profile; a question about what a factor regulates wants the graph. The router
        # names the modality now that it is told the bundle carries one, so `grn_tf` picks
        # the profile and anything else picks the network.
        _grn_modality = str(reading.get("modality") or "").strip().lower()
        if intent in ("tf_activity", "regulator_activity") or (
                intent == "regulatory_driver" and _grn_modality == "grn_tf"):
            answer = gnet.tf_activity_profile(
                ds, cell_state=state,
                contrast=contrast,
                factors=genes or None, limit=int(reading.get("limit") or 25))
            result.update(answer)
            if answer.get("by_cell_state"):
                result.update(gnet.tf_activity_state_chat(answer))
                return result
            rows = answer.get("rows") or []
            # `plot` IS A SPEC OBJECT, NOT A NAME. app.js:6924 reads `result.plot` and
            # dispatches on `spec.kind`; a bare string matches no branch and the panel
            # prints "This answer has no figure". A ranked activity chart is the bar chart
            # the chat already draws, so this reuses it rather than adding a renderer.
            #
            # ACTIVITY IS UNSIGNED AND ITS DIRECTION IS A SEPARATE COLUMN. A summed
            # regulatory activity is always positive, so `sign_column` carries the tested
            # log2 fold change and colours the bar by which way the factor moved.
            # app.js:7309 documents that exact case for a GO Z score.
            result["table"] = {
                "columns": ["factor", "activity", "log2fc", "fdr"],
                "rows": [[r["factor"], r["activity"], r["log2fc"], r["fdr"]]
                         for r in rows],
            }
            result["plot"] = {"kind": "barchart", "label_column": "factor",
                              "value_column": "activity", "sign_column": "log2fc"}
            if not rows:
                result["answer"] = (answer.get("note")
                                    or "This bundle carries no regulatory activity.")
            else:
                moved = [r for r in rows if r.get("tested")]
                lead = rows[0]
                result["answer"] = (
                    f"In {answer.get('cell_state')}, {lead['factor']} leads this ranking "
                    f"of {answer.get('n_factors')} modelled factors (reported change first, then activity); "
                    f"its summed regulatory activity is {lead['activity']}. "
                    + (f"{len(moved)} of the {len(rows)} shown changed in "
                       f"{answer.get('tf_activity_contrast') or 'the comparison'}."
                       if moved else
                       "None of the factors shown passed that comparison's own gates, so "
                       "the ranking is by activity alone."))
            result["follow_ups"] = _follow_ups(intent, ds, state, covariate, genes)
            return result

        if intent == "regulatory_driver" and state:
            answer = gnet.regulator_network(
                ds, state, contrast=contrast,
                features=genes or None,
                limit=int(reading.get("limit") or gnet.DEFAULT_LIMIT))
            result.update(answer)
            # A SPEC OBJECT, for the same reason as above. `kind: "network"` is new, and
            # app.js draws it from `result.nodes` and `result.edges`, which this answer
            # already carries, so no second request is made.
            result["plot"] = {"kind": "network",
                              "legend": answer.get("legend") or "",
                              "edge_score_range": answer.get("edge_score_range")}
            if not (answer.get("nodes") or []):
                result["answer"] = answer.get("note") or "No regulator reaches those features."
            else:
                result["answer"] = (
                    f"{answer['n_regulators']} factors reach "
                    f"{answer['n_targets_reached']} of the {answer['n_targets_requested']} "
                    f"features through {answer['n_edges_drawn']} edges in {state}. "
                    f"The edge floor is the {answer['edge_percentile']:.0f}th percentile of "
                    f"{state}'s own edge scores, {answer['edge_cut']}, and a factor had to "
                    f"clear {answer['expression_cut']} to be drawn, which dropped "
                    f"{answer['n_regulators_dropped_as_silent']} silent factors. "
                    f"{answer['n_factors_with_tested_activity']} factors carry a tested "
                    f"activity change here.")
            result["follow_ups"] = _follow_ups(intent, ds, state, covariate, genes)
            return result

        if intent in ("communication_rewiring", "pathway_program"):
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
        # Numeric obs fields, so the UMAP panel can plot ANY pair of float columns as
        # X and Y. app.js:6733 adds "obs columns (pick X and Y)" once more than one is
        # offered. Nathan, 2026-09-01: any float-valued obs field qualifies, excluding
        # counts and scaled counts. A field with a missing value is skipped, because an
        # axis needs a coordinate for every cell.
        COUNT_LIKE = ("n_counts", "n_cells", "n_cells_total", "n_genes_detected",
                      "metacell", "n_donors", "meta_sample_n_donors")
        numeric_variables = []
        for name, info in (ds.covariate_names() or {}).items():
            if info.get("kind") != "numeric":
                continue
            if name in COUNT_LIKE or name.endswith("__n_obs"):
                continue
            try:
                if int(float(info.get("n_missing") or 0)) != 0:
                    continue
            except (TypeError, ValueError):
                continue
            numeric_variables.append({"field": name, "min": info.get("min"),
                                      "max": info.get("max")})
        numeric_variables.sort(key=lambda v: v["field"])
        coords = [{"key": "", "label": "cellHarmony UMAP"}]
        for key in (ds.sv.get("embeddings") or []):
            if str(key) in ("X_umap", ""):
                continue
            coords.append({"key": str(key),
                           "label": str(key).replace("X_", "").replace("_", " ")})
        return {"cluster_key": cluster_key, "variables": variables,
                "color_variables": variables,
                "numeric_variables": numeric_variables,
                "coords": coords}

    @app.get("/api/jobs/{job_id}/dotplot")
    def dotplot(job_id: str, genes: str = Query(""), group_by: str = Query(""),
                groups: List[str] = Query([]), subset_by: str = Query(""),
                subset_values: List[str] = Query([]), subset2_by: str = Query(""),
                subset2_values: List[str] = Query([]), modality: str = Query("")):
        """Mean and detected fraction per (gene, cell state), from the bundle's
        precomputed stats matrices. Default gene set = top marker of every state."""
        ds = store.dataset(job_id)
        group_column, group_names, group_code, is_states = _bundle_group_axis(ds, group_by)
        if groups:
            chosen = [g for g in group_names if g in set(groups)]
            group_names = chosen
        features = _feature_source(ds, modality)
        pairs = ([] if features is not ds
                 else _dotplot_default_genes(assets.get(job_id, {}), ds))
        requested = split_gene_list(genes, keep_pipe=features.names_have_pipe,
                                    known=features.resolve_gene)
        wanted = (requested or [p["gene"] for p in pairs]
                  or _bundle_default_genes(ds, group_names, features))
        rows, missing, labels = [], [], []
        seen = set()
        for gene in wanted:
            if gene in seen:
                continue
            seen.add(gene)
            row = features.resolve_gene(gene)
            if row is None:
                missing.append(gene)
                continue
            rows.append(row)
            # THE AXIS SHOWS THE READER'S NAME WHATEVER THEY TYPED. A store without a
            # display column returns the key, so this is a no-op for every other modality.
            shown = features.display_of(row) if hasattr(features, "display_of") else ""
            labels.append(shown or gene)
        if not rows:
            raise HTTPException(404, f"none of the {len(wanted)} requested genes are in this dataset")
        restrict = _bundle_subset_masks(ds, subset_by, list(subset_values),
                                        subset2_by, list(subset2_values))
        if restrict is None and is_states and len(group_names) == len(ds.states):
            # Fast path: the bundle already holds the per-state statistics.
            mean = np.asarray(features.stats_mean[rows, :], dtype=np.float32).tolist()
            counts = list(ds.state_n)
            # Older bundles count stored nonzeros (including negative predictions).
            # Derive the displayed >0 fraction consistently with filtered views.
            frac = []
            state_codes = np.asarray(ds.state_code, dtype=np.int64)
            for row in rows:
                idx, val = features.gene_column(row)
                positive_cells = np.asarray(idx, dtype=np.int64)[np.asarray(val) > 0]
                hits = np.bincount(state_codes[positive_cells], minlength=len(ds.states))
                frac.append((hits / np.maximum(np.asarray(counts), 1)).tolist())
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
                idx, val = features.gene_column(row)
                sums = np.zeros(len(group_names), dtype=np.float64)
                hits = np.zeros(len(group_names), dtype=np.int64)
                if idx.size:
                    cells = idx.astype(np.int64)
                    keep_cells = per_cell[cells] >= 0
                    np.add.at(sums, per_cell[cells][keep_cells], val.astype(np.float64)[keep_cells])
                    np.add.at(hits, per_cell[cells][keep_cells], (np.asarray(val)[keep_cells] > 0).astype(np.int64))
                mean.append([float(sums[i] / counts[i]) if counts[i] else 0.0
                             for i in range(len(group_names))])
                frac.append([float(hits[i] / counts[i]) if counts[i] else 0.0
                             for i in range(len(group_names))])
        return {
            "modality": modality or "rna", "genes": labels, "states": group_names, "groups": group_names,
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
    "/Users/saljh8/Dropbox/LungMAP/refactored_website/lungmap-data/site/breath.sqlite")

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


def _to_feature_key(app, job_id: str, modality: str, gene: str) -> str:
    """The store's own key for whatever name the reader sent.

    WHY THIS IS NEEDED AND WHY IT IS NOT OPTIONAL. `_build_expression_payload` resolves
    against `cache_entry["var_names"]`, which holds the KEYS, and an unresolved name does
    not raise: it SILENTLY FALLS BACK TO FEATURE 0. Measured on this bundle 2026-09-08,
    `?gene=NOT_A_LIPID_AT_ALL` returned `resolved_gene=CE(18:2)`. So once lipids began
    showing "18:2 Cholesterol ester", every pick became an unresolved name and every plot
    drew CE(18:2) -- which is exactly what Nathan reported.

    Translating here, before the shipped route sees it, keeps the fallback untouched for
    every caller that was already sending a key. A name this cannot resolve is passed
    through unchanged, so the behaviour for a genuine typo is the old behaviour.
    """
    name = str(gene or "").strip()
    if not name:
        return gene
    try:
        store, _ = W._job_resources(app)
        ds = store.dataset(job_id)
        features = _feature_source(ds, modality)
        if features is ds or not hasattr(features, "display_of"):
            return gene
        row = features.resolve_gene(name)
        if row is None:
            return gene
        key = features.symbols[row]
        return str(key or gene)
    except Exception:                              # never fail a plot over a lookup
        return gene


def _retitle_with_display_name(app, job_id: str, modality: str, response):
    """Rewrite an expression payload's `gene` to the name a reader recognises.

    Returns the response untouched unless the modality carries a display column and the
    resolved feature has a different display name, so every other modality and every
    bundle built without one behaves exactly as before. A failure to rewrite returns the
    original payload rather than an error: a title is not worth a 500.
    """
    try:
        body = getattr(response, "body", None)
        if not body:
            return response
        payload = json.loads(body)
        key = str(payload.get("resolved_gene") or payload.get("gene") or "")
        if not key:
            return response
        # THE INSTALLER FOR THESE ROUTES IS HANDED ONLY `app`.
        # Taking `store` as an argument raised NameError on every expression request,
        # RNA included, because _install_expression_covariate_routes captures no store.
        store, _ = W._job_resources(app)
        ds = store.dataset(job_id)
        features = _feature_source(ds, modality)
        if features is ds or not hasattr(features, "display_of"):
            return response
        row = features.resolve_gene(key)
        if row is None:
            return response
        shown = features.display_of(row)
        if not shown or shown == payload.get("gene"):
            return response
        payload["gene"] = shown
        payload["feature_key"] = key
        return JSONResponse(payload)
    except Exception:                              # a title is never worth a failure
        return response


def _install_differential_feature_name_route(app, store) -> None:
    """The differential detail panel takes the name the volcano shows.

    WHY THIS IS REQUIRED, NOT A POLISH. Once the differential tables render lipids as
    "18:0/20:4 Phosphatidylcholine", clicking a volcano point sends THAT to
    `/differential/interactive/gene`, which looks the feature up in the differentials
    object by KEY and raised
        KeyError: Gene '18:0/20:4 Phosphatidylcholine' not found in the aligned AnnData
    So renaming the tables without this breaks every click-through.

    Translated in, relabelled out, exactly as the expression routes are.
    """
    from fastapi.routing import APIRoute

    path = "/api/jobs/{job_id}/differential/interactive/gene"
    original = None
    for route in list(app.router.routes):
        if isinstance(route, APIRoute) and route.path == path and "GET" in route.methods:
            original = route.endpoint
            break
    if original is None:
        return
    _drop_official_route(app, path)

    @app.get(path)
    def differential_gene_named(job_id: str, population: str = Query(...),
                                gene: str = Query(...),
                                feature: Optional[str] = Query(None)):
        # build_differential_block records the chosen contrast's modality under
        # `config`, not at the top of the block.
        meta = store.get_job(job_id)
        modality = str(((meta.get("differential") or {}).get("config") or {})
                       .get("modality") or "rna")
        key = _to_feature_key(app, job_id, modality, gene)
        response = original(job_id, population=population, gene=key, feature=feature)
        return _retitle_with_display_name(app, job_id, modality, response)


def _install_feature_name_routes(app, store) -> None:
    """The feature picker offers the name a reader recognises.

    Nathan, 2026-09-08: lipids must read "18:2 Cholesterol ester" rather than "CE(18:2)",
    with the abbreviation still searchable.

    WHY THIS REPLACES THE SHIPPED ROUTE RATHER THAN CHANGING THE CACHE. The official
    `/genes` builds its list from `cache_entry["var_names"]`, and those same strings are
    the expression cache's lookup keys. Renaming them there would rename the key as well
    as the label, which is the one thing this change must not do. So the key stays put and
    only the suggestion list is answered from the bundle's display column.

    A modality with no display column returns exactly what it returned before, because
    `ModalityStore.display` falls back to the feature name.
    """
    _drop_official_route(app, "/api/jobs/{job_id}/genes")

    @app.get("/api/jobs/{job_id}/genes")
    def job_genes(job_id: str, modality: str = Query("rna")):
        ds = store.dataset(job_id)
        features = _feature_source(ds, modality)
        shown = list(getattr(features, "display", None) or features.symbols)
        feature_label = "gene"
        if features is not ds:
            feature_label = getattr(features, "feature_label", "feature") or "feature"
        return {"genes": shown,
                "modality": (modality or "rna").strip().lower() or "rna",
                "feature_label": feature_label,
                # Both names resolve, so a caller may send either back.
                "keys": list(features.symbols) if shown != list(features.symbols) else None}


def _install_grn_routes(app, store) -> None:
    """The two regulatory views Chat draws, and the HTTP routes behind them.

    Nathan, 2026-09-07: "Both GRN and TF-activity should be able to have separate types of
    plots supported in Chat." So there are two routes, not one with a flag: a network
    answers who regulates a set of genes, and an activity profile answers which factors are
    most active in a cell state and which of them moved. The rules are the site's own, and
    `grn_network.py` records where each one comes from.

    These are ADDED, never substituted. The shipped
    `/api/jobs/{job_id}/grn/network` keeps its absolute-threshold behaviour, so nothing
    that already calls it changes.
    """

    _drop_official_route(app, "/api/jobs/{job_id}/grn/regulator-network")
    _drop_official_route(app, "/api/jobs/{job_id}/grn/tf-activity")

    @app.get("/api/jobs/{job_id}/grn/regulator-network")
    def grn_regulator_network(
        job_id: str,
        cell_state: str = Query(""),
        contrast: str = Query(""),
        features: str = Query(""),
        limit: int = Query(gnet.DEFAULT_LIMIT),
        edge_percentile: float = Query(gnet.DEFAULT_EDGE_PERCENTILE),
        expression_percentile: float = Query(gnet.DEFAULT_EXPRESSION_PERCENTILE),
    ):
        """The factors that regulate the features one differential put on screen.

        `limit` is the "Show" number: how many of the contrast's features become targets,
        ranked by absolute log2 fold change among those at or below FDR 0.05. It defaults
        to 200 and a reader may raise or lower it.
        """
        ds = store.dataset(job_id)
        want = [f.strip() for f in str(features).replace(",", " ").split() if f.strip()]
        return gnet.regulator_network(
            ds, cell_state, contrast=contrast, features=want or None,
            limit=int(limit), edge_percentile=float(edge_percentile),
            expression_percentile=float(expression_percentile))

    @app.get("/api/jobs/{job_id}/grn/tf-activity")
    def grn_tf_activity(
        job_id: str,
        cell_state: str = Query(""),
        contrast: str = Query(""),
        factors: str = Query(""),
        limit: int = Query(25),
    ):
        """Per-factor regulatory activity in one cell state, with its tested change."""
        ds = store.dataset(job_id)
        want = [f.strip() for f in str(factors).replace(",", " ").split() if f.strip()]
        return gnet.tf_activity_profile(
            ds, cell_state=cell_state, contrast=contrast,
            factors=want or None, limit=int(limit))


# =================================================================================
# GRN edges from the bundle's own per-cell-state store
# =================================================================================
#
# webapp/app.py `_grn_edges_adata` opens
# `meta['modality_artifacts']['grn']['network_h5ad']`, the edge-level object
# cellHarmony-differential writes beside each contrast's DEG tables. A precomputed bundle
# has no such object. Measured on the COPD atlas: 0 files match
# `differentials_only_*.h5ad` anywhere under the study tree, so `network_h5ad` was never
# set and every GRN edges request answered "GRN edge output is unavailable for this job."
#
# The edges were present the whole time, as the bundle's own `grn` modality store:
# `<prefix>_grn_genes.tsv` names 63,647 `TF|target` edges and `<prefix>_grn_stats_mean.npy`
# holds each edge's mean score in each of the 50 cell states. This assembles those two
# arrays into the cell-state x edge AnnData the payload builder already knows how to read.
# No score is recomputed and nothing is written: the values are the ones precompute stored.
#
# The store is per cell state, so the AnnData carries no sample column, `disp_col` stays
# None and the panel's Sample menu stays empty. That is the store's granularity, not a
# missing field.

_GRN_EDGE_ADATA_CACHE: Dict[str, Any] = {}


def _bundle_grn_edge_adata(ds, cluster_key: str):
    """Cell-state x edge AnnData from one bundle's `grn` store."""
    import anndata as ad

    store_grn = ds.modality("grn")
    mean = np.asarray(store_grn.stats_mean)          # (n_edges, n_states)
    states = [str(s) for s in ds.states]
    names = [str(n) for n in store_grn.features]
    if mean.shape != (len(names), len(states)):
        raise ValueError(
            f"GRN store is {mean.shape} but names {len(names)} edges and "
            f"{len(states)} cell states.")
    key = str(cluster_key or "cell_type")
    obs = pd.DataFrame({key: states}, index=pd.Index(states, name="cell_state"))
    var = pd.DataFrame(index=pd.Index(names, name="edge"))
    adata = ad.AnnData(X=np.ascontiguousarray(mean.T, dtype=np.float32), obs=obs, var=var)
    if adata.shape != (len(states), len(names)):
        raise ValueError(f"assembled GRN AnnData is {adata.shape}, expected "
                         f"{(len(states), len(names))}")
    return adata


def _install_grn_edge_adata(app, store) -> None:
    """Serve `_grn_edges_adata` from the bundle when the job is one.

    The wrap is on the shared webapp module, the same way the GO payload builder is
    wrapped above. An uploaded cellHarmony job still reaches the original, so the
    contrast-written edge h5ad remains the source there.
    """
    original = W._grn_edges_adata
    if getattr(original, "_bundle_aware", False):
        return

    def bundle_aware(meta):
        sv = (meta or {}).get("scalable_viewer") or {}
        bundle_dir = str(sv.get("bundle_dir") or "")
        if not bundle_dir:
            return original(meta)
        cached = _GRN_EDGE_ADATA_CACHE.get(bundle_dir)
        if cached is not None:
            return cached
        catalog = getattr(app.state, "catalog", None)
        ids = list(catalog.ids()) if catalog is not None else []
        job_id = next((jid for jid in ids
                       if str(store.dataset(jid).paths.bundle_dir) == bundle_dir), "")
        if not job_id:
            return original(meta)
        cluster_key = str(meta.get("cluster_key") or meta.get("reference_cluster_key") or "")
        adata = _bundle_grn_edge_adata(store.dataset(job_id), cluster_key)
        _GRN_EDGE_ADATA_CACHE[bundle_dir] = adata
        return adata

    bundle_aware._bundle_aware = True
    W._grn_edges_adata = bundle_aware


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
                                          "n_states": entry["n_states"], "label": entry["label"],
                                          "observation_unit": catalog.get(ds_id).sv.get("observation_unit", "cells")}
        if not record.get("ok"):
            return JSONResponse(record, status_code=200)
        return record


# --------------------------------------------------------------- protocol executors
#
# Five protocols that previously answered "not implemented". Each computes from
# the bundle: per-donor pseudobulk for the three that need a donor axis, cell
# counts for composition, and a covariate cross-tabulation for concordance.

def _run_severity_gradient(ds, *args, **kwargs):
    from altanalyze3.components.cellHarmony import chat_protocols
    from altanalyze3.components.cellHarmony.webapp.chat_data import protocol_adapter
    return chat_protocols._run_severity_gradient(protocol_adapter(ds), *args, **kwargs)


def _run_coexpression(ds, *args, **kwargs):
    from altanalyze3.components.cellHarmony import chat_protocols
    from altanalyze3.components.cellHarmony.webapp.chat_data import protocol_adapter
    return chat_protocols._run_coexpression(protocol_adapter(ds), *args, **kwargs)


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


def _run_composition(ds, *args, **kwargs):
    from altanalyze3.components.cellHarmony import chat_protocols
    from altanalyze3.components.cellHarmony.webapp.chat_data import protocol_adapter
    return chat_protocols._run_composition(protocol_adapter(ds), *args, **kwargs)


def _run_concordance(ds, *args, **kwargs):
    from altanalyze3.components.cellHarmony import chat_protocols
    from altanalyze3.components.cellHarmony.webapp.chat_data import protocol_adapter
    return chat_protocols._run_concordance(protocol_adapter(ds), *args, **kwargs)


def _run_dose_response(ds, *args, **kwargs):
    from altanalyze3.components.cellHarmony import chat_protocols
    from altanalyze3.components.cellHarmony.webapp.chat_data import protocol_adapter
    return chat_protocols._run_dose_response(protocol_adapter(ds), *args, **kwargs)


def _run_pathway_program(ds, *args, **kwargs):
    from altanalyze3.components.cellHarmony import chat_protocols
    from altanalyze3.components.cellHarmony.webapp.chat_data import protocol_adapter
    return chat_protocols._run_pathway_program(protocol_adapter(ds), *args, **kwargs)


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


def _run_most_affected_state(ds, *args, **kwargs):
    from altanalyze3.components.cellHarmony import chat_protocols
    from altanalyze3.components.cellHarmony.webapp.chat_data import protocol_adapter
    return chat_protocols._run_most_affected_state(protocol_adapter(ds), *args, **kwargs)


def _run_donor_heterogeneity(ds, *args, **kwargs):
    from altanalyze3.components.cellHarmony import chat_protocols
    from altanalyze3.components.cellHarmony.webapp.chat_data import protocol_adapter
    return chat_protocols._run_donor_heterogeneity(protocol_adapter(ds), *args, **kwargs)
