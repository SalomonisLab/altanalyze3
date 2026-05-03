"""
Annotate the filtered isoform catalog with literature citations drawn
from UniProt's own evidence attachments.

For each isoform row in uniprot_isoform_catalog.tsv we look up the same
UniProt entry in the raw JSONL and extract:

    * direct_pubmed_ids      PubMed IDs attached to the isoform's
                             molecule-tagged subcellular location
                             comment or its isoform note
    * similar_iso_pubmed_ids PubMed IDs attached to OTHER non-canonical
                             isoforms of the same gene (reasonably
                             similar isoform evidence per the user's
                             instruction)
    * gene_level_pubmed_ids  PubMed IDs attached to entry-level
                             trafficking / subcellular sentences
                             (not isoform-attributed)
    * evidence_quotes        Quoted sentences from UniProt comments
                             that mention trafficking, retention,
                             mislocalization, or this isoform's name
    * evidence_level         DIRECT_ISOFORM | SIMILAR_ISOFORM |
                             GENE_LEVEL | NONE

This does NOT fetch or read the actual papers. It records the
UniProt-curator-supplied provenance for each annotation. The user is
expected to spot-check the PubMed IDs before treating them as confirmed
mislocalization literature.

Run:
    python verify_catalog_citations.py
"""

from __future__ import annotations

import csv
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
DAEDALUS_ROOT = HERE.parents[2] / "Daedalus" / "phase_a"
RAW_JSONL = DAEDALUS_ROOT / "data" / "raw" / "uniprot_reviewed_human_mouse.jsonl.gz"
OUT_DIR = DAEDALUS_ROOT / "data" / "interim" / "tm_negatives"
IN_TSV = OUT_DIR / "uniprot_isoform_catalog.tsv"
OUT_TSV = OUT_DIR / "uniprot_isoform_catalog_with_citations.tsv"

ISO_INLINE_RE = re.compile(r"\[Isoform\s+([^\]]+)\]", re.IGNORECASE)
SENT_SPLIT_RE = re.compile(r"(?<=[.!?])\s+(?=[A-Z\[\(])")
TRAFFICKING_RE = re.compile(
    r"retain|retention|mislocal|misfold|mistraffic|trafficking|"
    r"fails?\s+to\s+reach|fails?\s+to\s+traffic|fails?\s+to\s+localize|"
    r"does\s+not\s+reach|does\s+not\s+traffic|not\s+expressed\s+at|"
    r"not\s+detected\s+at|absent\s+from|unable\s+to\s+(reach|traffic|localize)|"
    r"sequester|intracellularly\s+retained|endo[- ]?h\s+sensitive|"
    r"aberrant|abnormal|dominant[- ]negative|lacks?\s+(the|a)\s+(signal|transmembrane)|"
    r"truncated|deletion|ER\s+retention|retained\s+in\s+the",
    re.IGNORECASE,
)

NON_PM_KEYWORDS = [
    "endoplasmic reticulum", "sarcoplasmic reticulum", "golgi",
    "cytoplasm", "cytosol", "nucleus", "nuclear", "mitochondri",
    "peroxisome", "lysosome", "endosome", "intracellular",
]


def _pubmeds(evidences):
    """Extract PubMed IDs from a UniProt evidences list."""
    out = []
    for ev in evidences or []:
        if not isinstance(ev, dict):
            continue
        if (ev.get("source") or "").lower() == "pubmed":
            pid = ev.get("id")
            if pid:
                out.append(str(pid))
    return out


def _resolve_molecule(mol, iso_map):
    if not mol:
        return None
    bare = re.sub(r"^[Ii]soform\s+", "", mol).strip()
    for iid, meta in iso_map.items():
        if meta.get("name", "").strip().lower() == bare.lower():
            return iid
        for syn in meta.get("synonyms") or []:
            if syn.strip().lower() == bare.lower():
                return iid
    m = re.match(r"^\s*[Ii]soform\s+(\w+)\s*$", mol)
    if m:
        tail = m.group(1)
        for iid in iso_map:
            if iid.endswith(f"-{tail}"):
                return iid
    return None


def _collect_iso_map(comments):
    out = {}
    for c in comments:
        if c.get("commentType") != "ALTERNATIVE PRODUCTS":
            continue
        for iso in c.get("isoforms") or []:
            ids = iso.get("isoformIds") or []
            name_obj = iso.get("name") or {}
            name = name_obj.get("value", "") if isinstance(name_obj, dict) else ""
            syns = [s.get("value", "") for s in iso.get("synonyms") or [] if isinstance(s, dict)]
            for iid in ids:
                out[iid] = {"name": name, "synonyms": syns}
    return out


def _attribute_sentence(sent, iso_map):
    """Try to resolve a sentence to a specific iso id (or None)."""
    m = ISO_INLINE_RE.search(sent)
    if m:
        iid = _resolve_molecule(f"Isoform {m.group(1).strip()}", iso_map)
        if iid:
            return iid
    low = sent.lower()
    for iid, meta in iso_map.items():
        for nm in [meta.get("name", "")] + (meta.get("synonyms") or []):
            if nm and len(nm) > 2 and nm.lower() in low:
                return iid
    m2 = re.search(r"\bisoform\s+(\w+)\b", low)
    if m2:
        return _resolve_molecule(f"Isoform {m2.group(1)}", iso_map)
    return None


def _extract_evidence(entry):
    """Return per-iso dict of:
        {iso_id: {
            "direct_pmids": [...],
            "quotes": [ (pmids, sentence, source_comment_type) ],
        }}
       plus entry-level:
        "_GENE_": {
            "gene_level_pmids": [...],
            "gene_level_quotes": [ (pmids, sentence, ctype) ],
            "all_non_canonical_pmids_by_iso": {iid: [pmids]}
        }
    """
    comments = entry.get("comments") or []
    iso_map = _collect_iso_map(comments)
    per_iso = defaultdict(lambda: {"direct_pmids": set(), "quotes": []})
    gene_level_pmids = set()
    gene_level_quotes = []

    for c in comments:
        ctype = c.get("commentType", "")
        # Molecule-tagged subcellular location: direct isoform evidence
        if ctype == "SUBCELLULAR LOCATION":
            mol = c.get("molecule")
            if mol:
                iid = _resolve_molecule(mol, iso_map)
                if iid:
                    for sl in c.get("subcellularLocations") or []:
                        loc_obj = sl.get("location") or {}
                        pmids = _pubmeds(loc_obj.get("evidences"))
                        per_iso[iid]["direct_pmids"].update(pmids)
                        loc = loc_obj.get("value") or ""
                        if loc and pmids:
                            per_iso[iid]["quotes"].append(
                                (pmids, f"[{iid}] subcellular location: {loc}", ctype)
                            )

        # Note inside comments: may have inline [Isoform X] markers
        note_obj = c.get("note")
        note_texts = note_obj.get("texts") if isinstance(note_obj, dict) else []
        for nt in note_texts or []:
            text = nt.get("value") if isinstance(nt, dict) else None
            if not text:
                continue
            pmids = _pubmeds(nt.get("evidences")) if isinstance(nt, dict) else []
            for sent in SENT_SPLIT_RE.split(text):
                hit_traffic = TRAFFICKING_RE.search(sent)
                hit_non_pm = any(kw in sent.lower() for kw in NON_PM_KEYWORDS)
                if not (hit_traffic or hit_non_pm):
                    continue
                iid = _attribute_sentence(sent, iso_map)
                if iid:
                    per_iso[iid]["direct_pmids"].update(pmids)
                    if sent.strip():
                        per_iso[iid]["quotes"].append(
                            (pmids, sent.strip(), ctype)
                        )
                else:
                    gene_level_pmids.update(pmids)
                    gene_level_quotes.append((pmids, sent.strip(), ctype))

        # Plain comment texts: FUNCTION, MISCELLANEOUS, TISSUE
        for nt in c.get("texts") or []:
            text = nt.get("value") if isinstance(nt, dict) else None
            if not text:
                continue
            pmids = _pubmeds(nt.get("evidences")) if isinstance(nt, dict) else []
            for sent in SENT_SPLIT_RE.split(text):
                hit_traffic = TRAFFICKING_RE.search(sent)
                hit_non_pm = any(kw in sent.lower() for kw in NON_PM_KEYWORDS)
                if not (hit_traffic or hit_non_pm):
                    continue
                iid = _attribute_sentence(sent, iso_map)
                if iid:
                    per_iso[iid]["direct_pmids"].update(pmids)
                    per_iso[iid]["quotes"].append(
                        (pmids, sent.strip(), ctype)
                    )
                else:
                    gene_level_pmids.update(pmids)
                    gene_level_quotes.append((pmids, sent.strip(), ctype))

        # Alternative sequence feature evidences: isoform-attributed
        # (not via attribution here but via iso_tokens downstream; captured
        # via the sentence mining above if description is in a comment).

    # ALTERNATIVE PRODUCTS isoform notes
    for c in comments:
        if c.get("commentType") != "ALTERNATIVE PRODUCTS":
            continue
        for iso in c.get("isoforms") or []:
            ids = iso.get("isoformIds") or []
            nobj = iso.get("note")
            if not isinstance(nobj, dict):
                continue
            for nt in nobj.get("texts") or []:
                text = nt.get("value") if isinstance(nt, dict) else None
                if not text:
                    continue
                pmids = _pubmeds(nt.get("evidences") if isinstance(nt, dict) else None)
                for sent in SENT_SPLIT_RE.split(text):
                    if not (TRAFFICKING_RE.search(sent) or
                            any(kw in sent.lower() for kw in NON_PM_KEYWORDS)):
                        continue
                    for iid in ids:
                        per_iso[iid]["direct_pmids"].update(pmids)
                        per_iso[iid]["quotes"].append(
                            (pmids, sent.strip(), "ALT_PRODUCTS_NOTE")
                        )

    return per_iso, gene_level_pmids, gene_level_quotes


def main():
    if not IN_TSV.exists():
        sys.exit(f"missing: {IN_TSV}")
    if not RAW_JSONL.exists():
        sys.exit(f"missing: {RAW_JSONL}")

    rows = list(csv.DictReader(open(IN_TSV), delimiter="\t"))
    primary_of = {r["uniprotkb_id"]: r["isoform_id"].split("-")[0] for r in rows}
    wanted = set(primary_of.values())

    # Build lookup: accession -> (per_iso_evidence, gene_level_pmids, quotes)
    evidence_by_acc = {}
    with gzip.open(RAW_JSONL, "rt") as h:
        for line in h:
            try:
                d = json.loads(line)
            except Exception:
                continue
            acc = d.get("primaryAccession")
            if acc not in wanted:
                continue
            evidence_by_acc[acc] = _extract_evidence(d)

    header = list(rows[0].keys()) + [
        "evidence_level", "direct_pubmed_ids",
        "similar_iso_pubmed_ids", "gene_level_pubmed_ids",
        "evidence_quotes",
    ]
    with open(OUT_TSV, "w") as out:
        out.write("\t".join(header) + "\n")
        n_direct = n_similar = n_gene = n_none = 0
        for r in rows:
            acc = r["isoform_id"].split("-")[0]
            iid = r["isoform_id"]
            ev = evidence_by_acc.get(acc)
            if ev is None:
                direct, gl_pmids, gl_quotes = {}, set(), []
            else:
                direct, gl_pmids, gl_quotes = ev

            direct_pmids = sorted(direct.get(iid, {}).get("direct_pmids") or set())
            direct_quotes = direct.get(iid, {}).get("quotes") or []

            similar_pmids = set()
            similar_quotes = []
            for other_iid, other_ev in direct.items():
                if other_iid == iid:
                    continue
                similar_pmids.update(other_ev.get("direct_pmids") or set())
                similar_quotes.extend(other_ev.get("quotes") or [])
            similar_pmids = sorted(similar_pmids)

            gene_pmids = sorted(gl_pmids)

            if direct_pmids:
                level = "DIRECT_ISOFORM"
                n_direct += 1
            elif similar_pmids:
                level = "SIMILAR_ISOFORM"
                n_similar += 1
            elif gene_pmids:
                level = "GENE_LEVEL"
                n_gene += 1
            else:
                level = "NONE"
                n_none += 1

            quotes = []
            for pmids, sent, ctype in direct_quotes[:6]:
                tag = f"PMID:{','.join(pmids)}" if pmids else "no-pmid"
                quotes.append(f"[DIRECT][{ctype}][{tag}] {sent}")
            for pmids, sent, ctype in similar_quotes[:4]:
                tag = f"PMID:{','.join(pmids)}" if pmids else "no-pmid"
                quotes.append(f"[SIMILAR][{ctype}][{tag}] {sent}")
            if level == "GENE_LEVEL":
                for pmids, sent, ctype in gl_quotes[:3]:
                    tag = f"PMID:{','.join(pmids)}" if pmids else "no-pmid"
                    quotes.append(f"[GENE][{ctype}][{tag}] {sent}")

            r["evidence_level"] = level
            r["direct_pubmed_ids"] = ",".join(direct_pmids)
            r["similar_iso_pubmed_ids"] = ",".join(similar_pmids)
            r["gene_level_pubmed_ids"] = ",".join(gene_pmids)
            r["evidence_quotes"] = (" || ".join(quotes))[:6000].replace(
                "\t", " ").replace("\n", " ")

            out.write("\t".join(
                str(r.get(k, "")).replace("\t", " ").replace("\n", " ")
                for k in header
            ) + "\n")

    print(f"Wrote: {OUT_TSV}")
    print(f"  DIRECT_ISOFORM   : {n_direct}")
    print(f"  SIMILAR_ISOFORM  : {n_similar}")
    print(f"  GENE_LEVEL       : {n_gene}")
    print(f"  NONE             : {n_none}")


if __name__ == "__main__":
    main()
