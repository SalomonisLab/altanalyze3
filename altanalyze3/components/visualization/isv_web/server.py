"""FastAPI server for the interactive ISV web app.

The single RunContext (loaded once at startup) holds the catalog + caches; every endpoint reads from
it. The heavy isoform query is memoized per (gene + selection + thresholds) signature so repeated /
threshold-tweak requests are instant.
"""
from __future__ import annotations

import json
import os
from typing import List, Optional

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel

from . import data_api as da

_HERE = os.path.dirname(os.path.abspath(__file__))
_TEMPLATES = Jinja2Templates(directory=os.path.join(_HERE, "templates"))


class IsoformQuery(BaseModel):
    gene: str
    samples: Optional[List[str]] = None
    groups: Optional[List[str]] = None
    cell_types: Optional[List[str]] = None
    combine_by: str = "group"
    cluster_similarity_threshold: float = 0.85
    min_split_fraction: float = 0.05
    min_count: int = 1
    max_isoforms: int = 1500
    cluster_strategy: str = "subsequence"
    cluster_mode: str = "block"
    filter_junctions: Optional[List[str]] = None
    include_introns: bool = True


class MoleculeQuery(BaseModel):
    gene: str
    scope: str = "combined"               # 'combined' (pool the selection) or 'sample' (one library)
    sample: Optional[str] = None          # required when scope == 'sample'
    samples: Optional[List[str]] = None
    groups: Optional[List[str]] = None
    cell_types: Optional[List[str]] = None
    min_count: int = 1
    max_isoforms: int = 400


class ReadsQuery(BaseModel):
    gene: str
    cell_types: Optional[List[str]] = None   # cell-state selection; None == all
    conditions: Optional[List[str]] = None   # covariates; None == all groups
    max_isoforms: int = 300                  # per-panel molecule cap (engine default)
    panel_by: str = "covariate"              # 'covariate' (panel per covariate) | 'cell_type' (panel per cell type)


def create_app(ctx: da.RunContext) -> FastAPI:
    app = FastAPI(title="AltAnalyze3 ISV Viewer", version="1.0")
    app.state.ctx = ctx
    app.state.query_cache = {}   # signature -> response dict (simple in-process memo)

    static_dir = os.path.join(_HERE, "static")
    app.mount("/static", StaticFiles(directory=static_dir), name="static")

    _index_path = os.path.join(_HERE, "templates", "index.html")

    @app.get("/", response_class=HTMLResponse)
    def index():
        # index.html is a static shell (no template vars); serve it directly to avoid Starlette
        # TemplateResponse signature differences across versions.
        with open(_index_path, "r") as fh:
            return HTMLResponse(fh.read())

    @app.get("/healthz")
    def healthz():
        return {"ok": True, "genes": len(ctx.all_genes), "samples": len(ctx.samples)}

    @app.get("/api/catalog")
    def catalog():
        return {
            "samples": [{"library": s, "group": ctx.samples[s]["group"]} for s in ctx.sample_order],
            "groups": sorted(ctx.groups.keys()),
            "cell_types": ctx.cell_types,
            "gene_count": len(ctx.all_genes),
            "defaults": {
                "cluster_similarity_threshold": 0.85, "min_split_fraction": 0.05,
                "min_count": 1, "max_isoforms": 1500,
                "cluster_strategy": "subsequence", "cluster_mode": "block",
            },
        }

    @app.get("/api/genes")
    def genes(q: str = Query("", min_length=0), limit: int = 25):
        return {"matches": ctx.search_genes(q, limit=limit)}

    @app.get("/api/junctions")
    def junctions(gene: str):
        g = ctx.resolve_gene(gene)
        if not g:
            raise HTTPException(404, f"gene not found: {gene}")
        return {"gene": g, "junctions": da.gene_junctions(ctx, g)}

    @app.post("/api/isoforms")
    def isoforms(q: IsoformQuery):
        g = ctx.resolve_gene(q.gene)
        if not g:
            raise HTTPException(404, f"gene not found in this run: {q.gene}")
        sig = json.dumps({**q.dict(), "gene": g}, sort_keys=True, default=str)
        cache = app.state.query_cache
        if sig in cache:
            return cache[sig]
        res = da.query_isoforms(
            ctx, g, samples=q.samples, groups=q.groups, cell_types=q.cell_types,
            combine_by=q.combine_by,
            cluster_similarity_threshold=q.cluster_similarity_threshold,
            min_split_fraction=q.min_split_fraction, min_count=q.min_count,
            max_isoforms=q.max_isoforms, cluster_strategy=q.cluster_strategy,
            cluster_mode=q.cluster_mode, filter_junctions=q.filter_junctions,
            include_introns=q.include_introns,
        )
        if len(cache) > 256:   # bound memory
            cache.clear()
        cache[sig] = res
        return res

    @app.post("/api/molecules")
    def molecules(q: MoleculeQuery):
        g = ctx.resolve_gene(q.gene)
        if not g:
            raise HTTPException(404, f"gene not found in this run: {q.gene}")
        sig = "mol:" + json.dumps({**q.dict(), "gene": g}, sort_keys=True, default=str)
        cache = app.state.query_cache
        if sig in cache:
            return cache[sig]
        res = da.query_molecules(
            ctx, g, scope=q.scope, sample=q.sample, samples=q.samples, groups=q.groups,
            cell_types=q.cell_types, min_count=q.min_count, max_isoforms=q.max_isoforms,
        )
        if len(cache) > 256:
            cache.clear()
        cache[sig] = res
        return res

    @app.post("/api/reads")
    def reads(q: ReadsQuery):
        """Read-level (molecule) view: per-covariate panels of individual molecules from the engine's
        own *_isoform_ids.tsv output. First call for a (gene, selection) may take ~10-30s (the engine
        builds counts caches); subsequent calls are served from the on-disk + in-process cache."""
        g = ctx.resolve_gene(q.gene)
        if not g:
            raise HTTPException(404, f"gene not found in this run: {q.gene}")
        sig = "reads:" + json.dumps({**q.dict(), "gene": g, "grp": ctx.molecule_grouping},
                                    sort_keys=True, default=str)
        cache = app.state.query_cache
        if sig in cache:
            return cache[sig]
        res = da.query_reads(ctx, g, cell_types=q.cell_types, conditions=q.conditions,
                             max_isoforms=q.max_isoforms, panel_by=q.panel_by)
        if len(cache) > 256:
            cache.clear()
        cache[sig] = res
        return res

    @app.get("/api/grouping")
    def grouping_get():
        """Current molecule-view grouping mode."""
        return {"molecule_grouping": ctx.molecule_grouping, "modes": ["final_isoform", "structure"]}

    @app.post("/api/grouping")
    def grouping_set(mode: str = Query(..., description="final_isoform | structure")):
        """Backend switch for the molecule-view grouping (final-isoform vs legacy live structure-clustering).
        Flushes the query cache so the next read render uses the new mode.
            curl -X POST 'http://HOST:PORT/api/grouping?mode=structure'
        """
        if mode not in ("final_isoform", "structure"):
            raise HTTPException(400, "mode must be 'final_isoform' or 'structure'")
        ctx.molecule_grouping = mode
        app.state.query_cache.clear()
        return {"molecule_grouping": ctx.molecule_grouping}

    @app.get("/api/isoform/{isoform_id:path}/protein", response_class=PlainTextResponse)
    def protein(isoform_id: str):
        fa = ctx.protein_fasta(isoform_id)
        if not fa:
            raise HTTPException(404, f"no protein sequence for isoform: {isoform_id}")
        return PlainTextResponse(fa, headers={
            "Content-Disposition": f'attachment; filename="{_safe(isoform_id)}.fasta"',
        })

    @app.get("/api/isoform/{isoform_id:path}/mrna", response_class=PlainTextResponse)
    def mrna(isoform_id: str):
        fa = ctx.mrna_fasta(isoform_id)
        if not fa:
            raise HTTPException(404, f"no mRNA sequence for isoform: {isoform_id}")
        return PlainTextResponse(fa, headers={
            "Content-Disposition": f'attachment; filename="{_safe(isoform_id)}.mrna.fasta"',
        })

    @app.get("/api/isoform/{isoform_id:path}/orf", response_class=PlainTextResponse)
    def orf(isoform_id: str):
        fa = ctx.orf_fasta(isoform_id)
        if not fa:
            raise HTTPException(404, f"no ORF/CDS sequence for isoform: {isoform_id}")
        return PlainTextResponse(fa, headers={
            "Content-Disposition": f'attachment; filename="{_safe(isoform_id)}.orf.fasta"',
        })

    @app.post("/api/proteins")
    def proteins(ids: List[str]):
        """Export multiple isoform protein sequences as one FASTA (e.g. a whole cluster)."""
        chunks = [fa for i in ids if (fa := ctx.protein_fasta(i))]
        if not chunks:
            raise HTTPException(404, "no protein sequences found for the requested isoforms")
        body = "".join(c if c.endswith("\n") else c + "\n" for c in chunks)
        return PlainTextResponse(body, headers={
            "Content-Disposition": 'attachment; filename="isoform_proteins.fasta"',
        })

    return app


def _safe(name):
    return "".join(c if (c.isalnum() or c in "._-") else "_" for c in str(name))
