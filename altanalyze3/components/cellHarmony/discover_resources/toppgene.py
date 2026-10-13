"""ToppGene enrichment for the Discover tab.

Nathan, 2026-09-11: "Add a ToppGene tab which should use the toppgene API (as LungChat does).
First, no analysis should be run when the user selects it, but the user will have a drop down
to select the major gene set categories ... then have a run option."

NOTHING RUNS UNTIL ASKED. Every other Discover view draws as soon as its tab is opened,
because its data is already in the index. This one posts a gene list to a server at CCHMC and
waits, so opening the tab only offers the choice; `/api/discover/toppgene` is called by the
Run button alone.

THE API, verified against the live service on 2026-09-11 rather than from documentation:

    POST https://toppgene.cchmc.org/API/lookup
         {"Symbols": ["SPP1", ...]}
      -> {"Genes": [{"OfficialSymbol", "Entrez", "Submitted", "Description"}, ...]}

    POST https://toppgene.cchmc.org/API/enrich
         {"Genes": [6696, ...], "Categories": [{"Type", "PValue", "MinGenes",
                                                "MaxGenes", "MaxResults", "Correction"}]}
      -> {"Annotations": [{"Category", "ID", "Name", "PValue", "QValueFDRBH",
                           "QValueFDRBY", "QValueBonferroni", "TotalGenes", "GenesInTerm",
                           "GenesInQuery", "GenesInTermInQuery", "Genes", "Source", "URL"}]}

Symbols must become Entrez ids first: the enrich endpoint takes ids, not names.

THE CERTIFICATE. toppgene.cchmc.org presents a self-signed certificate in its chain, so a
default verified request fails with SSLError on a machine that does not carry the CCHMC root.
Verification stays ON. Set LUNGMAP_TOPPGENE_CA to the CA bundle, or -- knowing what it costs
-- LUNGMAP_TOPPGENE_INSECURE=1 to skip verification. Defaulting to insecure would have made
every install silently accept any certificate for a request carrying the reader's gene list.
"""

from __future__ import annotations

import json
import os
import urllib.request
import ssl
from typing import Any

BASE = os.environ.get("LUNGMAP_TOPPGENE_URL", "https://toppgene.cchmc.org/API")
TIMEOUT = float(os.environ.get("LUNGMAP_TOPPGENE_TIMEOUT", "60"))

# The categories the dropdown offers. Value is the API's own `Type`; the label is what a
# reader recognises. Grouped so Gene Ontology reads as one idea with three branches.
CATEGORIES: list[dict[str, str]] = [
    # THE LIST IS LUNGCHAT'S, NOT ONE I ASSEMBLED. Nathan pointed at
    # src/lungchat/agents/research_agent/tools/gene_functional_enrichment_tool.py:566 in the LungChat repository,
    # which names all nineteen types the API accepts. My first version invented
    # "ComputationalGeneSet" (the API calls it "Computational") and omitted ToppCell,
    # Cytoband and GeneFamily -- ToppCell being the one Nathan asked for by name.
    {"value": "GeneOntologyBiologicalProcess", "label": "GO: Biological Process",
     "group": "Gene Ontology"},
    {"value": "GeneOntologyMolecularFunction", "label": "GO: Molecular Function",
     "group": "Gene Ontology"},
    {"value": "GeneOntologyCellularComponent", "label": "GO: Cellular Component",
     "group": "Gene Ontology"},
    {"value": "Pathway", "label": "Pathway", "group": "Pathways and interactions"},
    {"value": "Interaction", "label": "Protein interaction",
     "group": "Pathways and interactions"},
    {"value": "Domain", "label": "Protein domain", "group": "Pathways and interactions"},
    {"value": "ToppCell", "label": "ToppCell", "group": "Expression"},
    {"value": "Coexpression", "label": "Co-expression", "group": "Expression"},
    {"value": "CoexpressionAtlas", "label": "Co-expression atlas", "group": "Expression"},
    {"value": "Computational", "label": "Computational gene set", "group": "Expression"},
    {"value": "TFBS", "label": "Transcription factor binding site", "group": "Regulation"},
    {"value": "MicroRNA", "label": "microRNA target", "group": "Regulation"},
    {"value": "HumanPheno", "label": "Human phenotype", "group": "Phenotype and disease"},
    {"value": "MousePheno", "label": "Mouse phenotype", "group": "Phenotype and disease"},
    {"value": "Disease", "label": "Disease", "group": "Phenotype and disease"},
    {"value": "Drug", "label": "Drug", "group": "Phenotype and disease"},
    {"value": "GeneFamily", "label": "Gene family", "group": "Gene structure"},
    {"value": "Cytoband", "label": "Cytoband", "group": "Gene structure"},
    {"value": "Pubmed", "label": "PubMed co-citation", "group": "Literature"},
]
VALID = {c["value"] for c in CATEGORIES}


def _context() -> ssl.SSLContext | None:
    ca = os.environ.get("LUNGMAP_TOPPGENE_CA")
    if ca and os.path.exists(ca):
        return ssl.create_default_context(cafile=ca)
    if os.environ.get("LUNGMAP_TOPPGENE_INSECURE") == "1":
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        return ctx
    return None                      # the default verified context


def _post(path: str, body: dict) -> dict:
    req = urllib.request.Request(
        f"{BASE}/{path}", data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=TIMEOUT, context=_context()) as fh:
        return json.loads(fh.read().decode())


def categories() -> list[dict[str, str]]:
    return list(CATEGORIES)


# The measurements whose features ARE gene symbols. Nathan, 2026-09-11: "No ToppGene for
# imputed modalities (except grn_tfs)." A transcription factor is a gene, so grn_tf qualifies;
# an antibody name, a lipid species, a "TF|target" edge and a cell-cell interaction are not
# symbols, and posting them returns an empty enrichment that reads like a biological result
# rather than the category error it is.
GENE_MODALITIES = {"rna", "grn_tf"}


def modality_allowed(modality: str) -> bool:
    return (modality or "rna").strip().lower() in GENE_MODALITIES


def enrich(symbols: list[str], category: str,
           max_results: int = 50, p_value: float = 0.05,
           min_genes: int = 2, max_genes: int = 1500) -> dict[str, Any]:
    """Run one category over one gene list, and say plainly what was asked and answered."""
    symbols = [s.strip() for s in symbols if s and s.strip()]
    if not symbols:
        return {"terms": [], "error": "No genes to test."}
    if category not in VALID:
        return {"terms": [], "error": f"Unknown category {category!r}."}
    try:
        found = _post("lookup", {"Symbols": symbols}).get("Genes") or []
    except Exception as exc:                        # noqa: BLE001
        return {"terms": [], "error": f"ToppGene lookup failed: {exc}"}
    ids = [g["Entrez"] for g in found if g.get("Entrez")]
    resolved = {str(g.get("Submitted") or g.get("OfficialSymbol")) for g in found}
    missing = [s for s in symbols if s not in resolved]
    if not ids:
        return {"terms": [], "error": "ToppGene recognised none of these gene symbols.",
                "n_submitted": len(symbols), "n_resolved": 0}
    try:
        out = _post("enrich", {
            "Genes": ids,
            "Categories": [{"Type": category, "PValue": float(p_value),
                            "MinGenes": int(min_genes), "MaxGenes": int(max_genes),
                            "MaxResults": int(max_results), "Correction": "FDR"}]})
    except Exception as exc:                        # noqa: BLE001
        return {"terms": [], "error": f"ToppGene enrichment failed: {exc}",
                "n_submitted": len(symbols), "n_resolved": len(ids)}
    terms = []
    for a in out.get("Annotations") or []:
        # THE SCORE IS -log10 OF THE CORRECTED p, WHICH IS WHAT THE PLOT DRAWS. Nathan asked
        # for a GO-Elite-like figure using "the score data for the stats"; GO-Elite ranks by
        # the same quantity, so a reader moving between the two is reading one scale.
        q = a.get("QValueFDRBH")
        p = a.get("PValue")
        import math
        score = -math.log10(max(q if q is not None else (p or 1.0), 1e-300))
        terms.append({
            "id": a.get("ID"), "name": a.get("Name"), "category": a.get("Category"),
            "p_value": p, "q_fdr_bh": q,
            "q_bonferroni": a.get("QValueBonferroni"), "q_fdr_by": a.get("QValueFDRBY"),
            "genes_in_term": a.get("GenesInTerm"), "genes_in_query": a.get("GenesInQuery"),
            "hits": a.get("GenesInTermInQuery"), "total_genes": a.get("TotalGenes"),
            "score": score,
            "url": a.get("URL"), "source": a.get("Source"),
            "gene_symbols": [g.get("Symbol") for g in (a.get("Genes") or [])
                             if g.get("Symbol")],
        })
    terms.sort(key=lambda t: -t["score"])
    return {
        "terms": terms, "category": category,
        "n_submitted": len(symbols), "n_resolved": len(ids), "n_unresolved": len(missing),
        "unresolved": missing[:25],
        "note": (f"{len(terms)} terms at FDR (BH) below {p_value} over {len(ids)} of "
                 f"{len(symbols)} symbols ToppGene could resolve. Score is "
                 f"-log10 of the BH-corrected q."),
    }


def tsv(result: dict[str, Any]) -> str:
    """The terms as a table, with the statistic behind each one named in its header."""
    head = ["term_id", "term_name", "category", "score_minus_log10_q_fdr_bh",
            "p_value", "q_fdr_bh", "q_bonferroni", "q_fdr_by",
            "genes_in_term", "genes_in_query", "hits_in_term", "total_genes",
            "url", "hit_genes"]
    lines = [f"# ToppGene, category {result.get('category', '')}",
             f"# {result.get('note', '')}",
             "# score = -log10(q_fdr_bh); rows are ordered by it",
             "\t".join(head)]
    for t in result.get("terms") or []:
        lines.append("\t".join(str(x if x is not None else "") for x in [
            t["id"], t["name"], t["category"], round(t["score"], 4),
            t["p_value"], t["q_fdr_bh"], t["q_bonferroni"], t["q_fdr_by"],
            t["genes_in_term"], t["genes_in_query"], t["hits"], t["total_genes"],
            t["url"], ",".join(t.get("gene_symbols") or []),
        ]))
    return "\n".join(lines) + "\n"
