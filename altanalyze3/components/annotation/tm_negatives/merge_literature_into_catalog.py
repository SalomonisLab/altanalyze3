"""
Merge per-gene literature evidence (from WebSearches) into the filtered
isoform catalog. Produces uniprot_isoform_catalog_with_literature.tsv
where every row carries the gene-level mislocalization literature call
for its parent gene.
"""

from __future__ import annotations

import csv
from pathlib import Path

HERE = Path(__file__).resolve()
DAEDALUS_ROOT = HERE.parents[2] / "Daedalus" / "phase_a"
IN_DIR = DAEDALUS_ROOT / "data" / "interim" / "tm_negatives"
CATALOG = IN_DIR / "uniprot_isoform_catalog.tsv"
LIT = IN_DIR / "gene_literature_evidence.tsv"
OUT = IN_DIR / "uniprot_isoform_catalog_with_literature.tsv"


def main():
    lit_by_gene = {}
    for r in csv.DictReader(open(LIT), delimiter="\t"):
        lit_by_gene[r["gene_name"]] = r

    rows = list(csv.DictReader(open(CATALOG), delimiter="\t"))
    extra = [
        "literature_evidence_found",
        "literature_confidence",
        "literature_pubmed_ids",
        "literature_source_urls",
        "literature_evidence_snippet",
        "literature_search_query",
    ]
    header = list(rows[0].keys()) + extra
    n_yes = n_no = n_amb = n_missing = 0
    with open(OUT, "w") as h:
        h.write("\t".join(header) + "\n")
        for r in rows:
            lit = lit_by_gene.get(r["gene_name"], {})
            if lit:
                ef = lit.get("evidence_found", "")
                if ef == "YES": n_yes += 1
                elif ef == "NO": n_no += 1
                elif ef == "AMBIGUOUS": n_amb += 1
            else:
                n_missing += 1
            r["literature_evidence_found"] = lit.get("evidence_found", "MISSING")
            r["literature_confidence"] = lit.get("confidence", "")
            r["literature_pubmed_ids"] = lit.get("pubmed_ids", "")
            r["literature_source_urls"] = lit.get("source_urls", "")
            r["literature_evidence_snippet"] = lit.get("evidence_snippet", "")
            r["literature_search_query"] = lit.get("search_query_used", "")
            h.write("\t".join(
                str(r.get(k, "")).replace("\t", " ").replace("\n", " ")
                for k in header
            ) + "\n")

    print(f"Wrote: {OUT}")
    print(f"Rows by literature evidence: YES={n_yes}  AMBIGUOUS={n_amb}  NO={n_no}  MISSING={n_missing}")


if __name__ == "__main__":
    main()
