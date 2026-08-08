"""FastAPI server for the scalable_viewer.

One process serves many precomputed bundles. Datasets load lazily: the catalog reads
only the small metadata JSON of each bundle, and the first request for a dataset
memory-maps its arrays. No endpoint opens an h5ad.

Binary endpoints (*.bin) exist because the embedding of 161,432 cells is 1.3 MiB as
raw float32 and about 6 MiB as JSON. The client parses them with a DataView.
"""
from __future__ import annotations

import io
import json
import os
import time
from typing import List, Optional

import numpy as np
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import HTMLResponse, PlainTextResponse, Response
from fastapi.staticfiles import StaticFiles

from . import data_api as da

_HERE = os.path.dirname(os.path.abspath(__file__))

# Morpheus is loaded from the Broad CDN, the same source the validated cellHarmony
# webapp uses (app.py lines 349 and 377-378). No Morpheus asset ships in this repo,
# so the heatmap tab needs network access to software.broadinstitute.org.
MORPHEUS_CSS = "https://software.broadinstitute.org/morpheus/css/morpheus-latest.min.css"
MORPHEUS_JS_EXT = "https://software.broadinstitute.org/morpheus/js/morpheus-external-latest.min.js"
MORPHEUS_JS = "https://software.broadinstitute.org/morpheus/js/morpheus-latest.min.js"

# Two colour schemes, both given as RGB hex.
#  altanalyze : the scheme the cellHarmony webapp already uses (app.py line 585).
#  cyan_yellow: the linear cyan -> yellow ramp this project prefers for heatmaps.
COLOR_SCHEMES = {
    "altanalyze": {"values": [-2, 0, 2], "colors": ["#00F0FF", "#000000", "#FFFF00"]},
    "cyan_yellow": {"values": [-2, 2], "colors": ["#00FFFF", "#FFFF00"]},
}
FALLBACK_COLOR = "#BBBBBB"


def _ds(app: FastAPI, ds_id: str) -> da.Dataset:
    try:
        return app.state.catalog.get(ds_id)
    except KeyError:
        raise HTTPException(404, f"unknown dataset: {ds_id}. "
                                 f"known: {app.state.catalog.ids()}")
    except (FileNotFoundError, ValueError) as exc:
        raise HTTPException(500, f"dataset {ds_id} failed to load: {exc}")


def create_app(catalog: da.Catalog) -> FastAPI:
    app = FastAPI(title="AltAnalyze3 scALABLE viewer", version="1.0")
    app.state.catalog = catalog
    app.mount("/static", StaticFiles(directory=os.path.join(_HERE, "static")), name="static")
    index_path = os.path.join(_HERE, "templates", "index.html")

    # ------------------------------------------------------------------ shell

    @app.get("/", response_class=HTMLResponse)
    def index():
        with open(index_path, "r") as fh:
            return HTMLResponse(fh.read())

    @app.get("/healthz")
    def healthz():
        return {"ok": True, "n_datasets": len(catalog.entries),
                "loaded": [e["id"] for e in catalog.entries if e["loaded"]],
                "load_errors": catalog.load_errors}

    @app.get("/api/catalog")
    def api_catalog():
        """Every dataset this server can serve. Reads metadata JSON only; nothing is loaded."""
        return {"datasets": catalog.entries, "load_errors": catalog.load_errors,
                "color_schemes": list(COLOR_SCHEMES.keys())}

    # ---------------------------------------------------------------- dataset

    @app.get("/api/dataset/{ds_id}/meta")
    def meta(ds_id: str):
        d = _ds(app, ds_id)
        return {
            "id": d.id, "label": d.label, "n_cells": d.n_cells, "n_genes": d.n_genes,
            "cluster_key": d.cluster_key, "layer": d.sv.get("layer"),
            "states": d.states, "state_n": d.state_n,
            "state_colors": [d.colors.get(s, FALLBACK_COLOR) for s in d.states],
            "states_missing_color": [s for s in d.states if s not in d.colors],
            "canonical_order_source": d.sv.get("canonical_order_source"),
            "n_states_in_canonical_order": len(d.meta.get("lineage_order", [])),
            "states_in_canonical_order_not_present": d.sv.get(
                "states_in_canonical_order_not_present", []),
            "covariates": d.covariate_names(),
            "embedding_method": d.sv.get("embedding_method"),
            "stats_method": d.sv.get("stats_method"),
            "hvg_method": d.sv.get("hvg_method"),
            "source_h5ad": d.meta.get("source_h5ad"),
            "bundle_dir": d.paths.bundle_dir, "prefix": d.paths.prefix,
            "built_utc": d.sv.get("built_utc"), "nnz": d.sv.get("nnz"),
            "markers": {k: v for k, v in d.sv.get("markers", {}).items() if k != "clusters"},
            "n_deg_tables": len(d.deg_manifest().get("comparisons", [])),
            "ccc": d.sv.get("ccc", {"available": False}),
            "warnings": d.sv.get("warnings", []),
        }

    @app.get("/api/dataset/{ds_id}/genes")
    def genes(ds_id: str, q: str = Query(""), limit: int = 25):
        d = _ds(app, ds_id)
        return {"matches": d.search_genes(q, limit=limit), "n_genes": d.n_genes}

    # ------------------------------------------------------------- embeddings

    @app.get("/api/dataset/{ds_id}/embedding.bin")
    def embedding_bin(ds_id: str):
        """Raw float32 [x0,y0,x1,y1,...]. Length is 8 * n_cells bytes."""
        d = _ds(app, ds_id)
        buf = np.ascontiguousarray(d.embedding, dtype="<f4").tobytes()
        return Response(buf, media_type="application/octet-stream",
                        headers={"X-N-Cells": str(d.n_cells), "Cache-Control": "no-store"})

    @app.get("/api/dataset/{ds_id}/colorby")
    def colorby(ds_id: str, key: str = Query("cell_state")):
        """Legend for a colouring. The values themselves come from colorby.bin."""
        d = _ds(app, ds_id)
        if key in ("cell_state", d.cluster_key, "__state__"):
            return {"key": "cell_state", "kind": "categorical", "dtype": "int16",
                    "categories": d.states,                       # CANONICAL order
                    "colors": [d.colors.get(s, FALLBACK_COLOR) for s in d.states],
                    "counts": d.state_n, "order_source": d.sv.get("canonical_order_source")}
        try:
            kind, vals, cats = d.covariate_values(key)
        except KeyError:
            raise HTTPException(404, f"unknown covariate: {key}. "
                                     f"known: {sorted(d.covariate_names().keys())}")
        if kind == "numeric":
            finite = np.isfinite(vals)
            return {"key": key, "kind": "numeric", "dtype": "float32",
                    "min": float(vals[finite].min()) if finite.any() else None,
                    "max": float(vals[finite].max()) if finite.any() else None,
                    "n_missing": int((~finite).sum()), "n_cells": int(vals.size)}
        counts = np.bincount(vals[vals >= 0].astype(np.int64), minlength=len(cats)).tolist()
        return {"key": key, "kind": "categorical", "dtype": "int32", "categories": cats,
                "colors": [_palette(i) for i in range(len(cats))], "counts": counts,
                "n_missing": int((vals < 0).sum())}

    @app.get("/api/dataset/{ds_id}/colorby.bin")
    def colorby_bin(ds_id: str, key: str = Query("cell_state")):
        d = _ds(app, ds_id)
        if key in ("cell_state", d.cluster_key, "__state__"):
            buf = np.ascontiguousarray(d.state_code, dtype="<i2").tobytes()
            return Response(buf, media_type="application/octet-stream",
                            headers={"X-Dtype": "int16", "Cache-Control": "no-store"})
        try:
            kind, vals, _ = d.covariate_values(key)
        except KeyError:
            raise HTTPException(404, f"unknown covariate: {key}")
        dt, tag = ("<f4", "float32") if kind == "numeric" else ("<i4", "int32")
        buf = np.ascontiguousarray(vals, dtype=dt).tobytes()
        return Response(buf, media_type="application/octet-stream",
                        headers={"X-Dtype": tag, "Cache-Control": "no-store"})

    @app.get("/api/dataset/{ds_id}/expression.bin")
    def expression_bin(ds_id: str, gene: str):
        """[n:int32][cell index:uint32 * n][value:float32 * n]. Cells not listed are zero."""
        d = _ds(app, ds_id)
        row = d.resolve_gene(gene)
        if row is None:
            raise HTTPException(404, f"gene not in this dataset: {gene}")
        idx, val = d.gene_column(row)
        return Response(da.pack_sparse(idx, val), media_type="application/octet-stream",
                        headers={"X-Gene": d.symbols[row], "X-N-Detected": str(int(idx.size)),
                                 "X-Max": str(float(val.max()) if val.size else 0.0),
                                 "Cache-Control": "no-store"})

    # ------------------------------------------------------------------ plots

    @app.get("/api/dataset/{ds_id}/violin")
    def violin(ds_id: str, gene: str, bins: int = 40):
        d = _ds(app, ds_id)
        row = d.resolve_gene(gene)
        if row is None:
            raise HTTPException(404, f"gene not in this dataset: {gene}")
        t0 = time.time()
        res = d.violin(row, n_bins=bins)
        res["state_colors"] = [d.colors.get(s, FALLBACK_COLOR) for s in d.states]
        res["elapsed_ms"] = round((time.time() - t0) * 1000, 1)
        return res

    @app.get("/api/dataset/{ds_id}/dotplot")
    def dotplot(ds_id: str, genes: str):
        d = _ds(app, ds_id)
        want = [g.strip() for g in genes.split(",") if g.strip()]
        rows, missing = [], []
        for g in want:
            r = d.resolve_gene(g)
            (rows.append(r) if r is not None else missing.append(g))
        if not rows:
            raise HTTPException(404, f"none of the {len(want)} requested genes are in this dataset")
        res = d.dotplot(rows)
        res.update({"n_requested": len(want), "n_returned": len(rows),
                    "n_missing": len(missing), "missing": missing,
                    "state_colors": [d.colors.get(s, FALLBACK_COLOR) for s in d.states]})
        return res

    # ---------------------------------------------------------------- markers

    @app.get("/api/dataset/{ds_id}/markers")
    def markers(ds_id: str, state: Optional[str] = None, top: int = 25):
        d = _ds(app, ds_id)
        rows = d.markers()
        n_all = len(rows)
        if state:
            rows = [r for r in rows if r["cluster"] == state]
        rows = sorted(rows, key=lambda r: -(r["fold"] if r["fold"] is not None else -1e9))[:top]
        return {"n_rows_total": n_all, "n_rows_returned": len(rows), "state": state,
                "rows": rows, "source": d.sv.get("markers", {}).get("source")}

    @app.get("/api/dataset/{ds_id}/heatmap.tsv", response_class=PlainTextResponse)
    def heatmap_tsv(ds_id: str, states: Optional[str] = None, top: int = 5,
                    row_scale: bool = True):
        """Morpheus-ready TSV. Columns are cell states in CANONICAL order."""
        d = _ds(app, ds_id)
        sel = [s.strip() for s in states.split(",")] if states else None
        text, info = d.heatmap_tsv(sel, top, row_scale)
        return PlainTextResponse(text, media_type="text/tab-separated-values",
                                 headers={"X-Rows": str(info["n_genes_written"]),
                                          "X-Cols": str(info["n_states"]),
                                          "X-Missing": str(info["n_genes_missing"]),
                                          "Cache-Control": "no-store"})

    @app.get("/api/dataset/{ds_id}/morpheus", response_class=HTMLResponse)
    def morpheus(ds_id: str, states: Optional[str] = None, top: int = 5,
                 row_scale: bool = True, scheme: str = "altanalyze"):
        """The Morpheus heatmap page, loaded in an iframe by the main shell.

        Reuses the embed the cellHarmony webapp already ships: CDN scripts, fetch the
        TSV, wrap it in a Blob, hand the blob URL to new morpheus.HeatMap({dataset: ...}).
        """
        d = _ds(app, ds_id)
        cs = COLOR_SCHEMES.get(scheme, COLOR_SCHEMES["altanalyze"])
        qs = f"?top={int(top)}&row_scale={'true' if row_scale else 'false'}"
        if states:
            qs += "&states=" + states.replace(" ", "%20")
        url = f"/api/dataset/{ds_id}/heatmap.tsv{qs}"
        html = _MORPHEUS_PAGE.replace("__CSS__", MORPHEUS_CSS) \
                             .replace("__JS_EXT__", MORPHEUS_JS_EXT) \
                             .replace("__JS__", MORPHEUS_JS) \
                             .replace("__DATASET_URL__", url) \
                             .replace("__VALUES__", json.dumps(cs["values"])) \
                             .replace("__COLORS__", json.dumps(cs["colors"])) \
                             .replace("__TITLE__", f"{d.label} markers")
        return HTMLResponse(html)

    # -------------------------------------------------------------------- DEG

    @app.get("/api/dataset/{ds_id}/deg")
    def deg_list(ds_id: str):
        d = _ds(app, ds_id)
        comps = d.deg_manifest().get("comparisons", [])
        return {"n_comparisons": len(comps),
                "comparisons": [{k: c[k] for k in ("id", "comparison", "kind", "n_rows",
                                                   "columns", "source")} for c in comps]}

    @app.get("/api/dataset/{ds_id}/deg/table")
    def deg_table(ds_id: str, id: str, max_rows: int = 500,
                  fdr_max: Optional[float] = None, state: Optional[str] = None):
        d = _ds(app, ds_id)
        try:
            return d.deg_table(id, max_rows=max_rows, fdr_max=fdr_max, state=state)
        except KeyError:
            raise HTTPException(404, f"unknown DEG table: {id}")

    # -------------------------------------------------------------------- CCC

    @app.get("/api/dataset/{ds_id}/ccc")
    def ccc(ds_id: str, limit: int = 5000):
        d = _ds(app, ds_id)
        return d.ccc(limit=limit)

    return app


def _palette(i: int) -> str:
    """Fixed RGB hex ring for covariate categories. No named colours, no rainbow."""
    ring = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860",
            "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD", "#1B4400", "#B79762",
            "#6A3A4C", "#0000A6", "#997D87", "#FF34FF", "#008941", "#A30059"]
    return ring[i % len(ring)]


_MORPHEUS_PAGE = """<!doctype html>
<html><head><meta charset="utf-8"><title>__TITLE__</title>
<link rel="stylesheet" href="__CSS__" />
<style>html,body{margin:0;height:100%;font-family:Arial,Helvetica,sans-serif;}
#t{position:absolute;inset:0;} #err{padding:14px;color:#C44E52;font-size:13px;}</style>
</head><body>
<div id="t"></div><div id="err" hidden></div>
<script src="__JS_EXT__"></script>
<script src="__JS__"></script>
<script>
(async function(){
  const fail = (m)=>{const e=document.getElementById('err');e.hidden=false;e.textContent=m;};
  try{
    const r = await fetch("__DATASET_URL__", {cache:"no-store"});
    if(!r.ok){ fail("heatmap.tsv returned HTTP "+r.status); return; }
    const text = await r.text();
    if(!window.morpheus || !window.morpheus.HeatMap){
      fail("Morpheus did not load from the Broad CDN. This tab needs network access to "
         + "software.broadinstitute.org. Rows="+text.split("\\n").length);
      return;
    }
    const blobUrl = URL.createObjectURL(new Blob([text],{type:"text/tab-separated-values"}));
    new window.morpheus.HeatMap({
      el: document.getElementById("t"),
      dataset: blobUrl,
      rowSize: 12,
      columnSize: 12,
      drawGrid: false,
      rows: [{field:"id", display:["text"]}],
      columns: [{field:"id", display:["text"]}],
      colorScheme: {scalingMode:"fixed", stepped:false, values:__VALUES__, colors:__COLORS__}
    });
  }catch(e){ fail("Morpheus embed failed: "+(e&&e.message?e.message:e)); }
}());
</script></body></html>
"""
