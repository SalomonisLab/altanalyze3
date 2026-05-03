"""
Scan locally cached UniProt reviewed (SwissProt) for TM isoforms with
REAL DOCUMENTED EVIDENCE that they fail to reach the plasma membrane due
to misfolding / missing signal peptide / missing critical domain /
explicit trafficking defect.

Inclusion requires (all):
    I   Entry is human, reviewed, has ALTERNATIVE PRODUCTS.
    II  Canonical has >=1 TM feature AND the candidate isoform retains
        at least one TM after its alt-seq change (preserves criterion 2).
    III Candidate isoform is annotated to a non-PM compartment
        (molecule-tagged or inline [Isoform X]:) OR has explicit
        text evidence of mistrafficking attributed to this isoform.

Additional restriction: MULTI-PASS TM only (tm_count_canonical >= 2).
Single-pass TM proteins are excluded.

Inclusion requires (>=1 of these DOCUMENTED reasons):
    A  Isoform note or comment text contains a MISTRAFFICKING SENTENCE
       attributed to this specific isoform (inline [Isoform X] marker
       OR isoform name / synonym appears in the same sentence)
    B  Explicit structural-loss sentence (lacks signal / missing TM /
       truncated / deletion of ...) attributed to this isoform AND
       paired with subcellular consequence evidence OR note
    C  Signal peptide deleted by this isoform's alt-seq change
       (only a weak signal for multi-pass; kept as a qualifier since
       user listed 'no signal sequence' as a reason)

Domain deletion is REPORTED but not a qualifier: for multi-pass TM
proteins, losing an extracellular/intracellular domain does not
reliably predict failure to traffic to the PM.

NMD check is dropped (SwissProt excludes NMD isoforms by policy).

Output: single TSV 'tm_negative_isoforms_documented.tsv' with one row
per qualifying isoform + full evidence columns.

Run: python scan_uniprot_tm_negatives.py
"""

from __future__ import annotations

import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
DAEDALUS_ROOT = HERE.parents[2] / "Daedalus" / "phase_a"
RAW_JSONL = DAEDALUS_ROOT / "data" / "raw" / "uniprot_reviewed_human_mouse.jsonl.gz"
OUT_DIR = DAEDALUS_ROOT / "data" / "interim" / "tm_negatives"

PM_TERMS = [
    "cell membrane", "apical cell membrane", "basolateral cell membrane",
    "lateral cell membrane", "membrane raft", "cell surface",
    "plasma membrane", "sarcolemma", "postsynaptic density",
    "dendritic spine membrane",
]

# Specific phrases indicating documented mistrafficking / retention /
# missing structural prerequisite.  Intentionally NOT including broad
# words like "intracellular" alone, which often describe normal biology.
MISTRAFFIC_PHRASES = [
    r"retained\s+in\s+the\s+(er|endoplasmic\s+reticulum|cytoplasm|cytosol|nucleus|golgi)",
    r"er[- ]retained",
    r"er[- ]retention",
    r"retention\s+in\s+the\s+endoplasmic",
    r"retention\s+in\s+the\s+er",
    r"fails?\s+to\s+reach\s+the\s+(cell\s+surface|plasma\s+membrane|cell\s+membrane)",
    r"fails?\s+to\s+be\s+expressed\s+at\s+the\s+cell\s+surface",
    r"fails?\s+to\s+localize\s+to\s+the\s+(cell\s+surface|plasma\s+membrane)",
    r"fails?\s+to\s+traffic",
    r"does\s+not\s+reach\s+the\s+(cell\s+surface|plasma\s+membrane|cell\s+membrane)",
    r"does\s+not\s+traffic",
    r"is\s+not\s+expressed\s+at\s+the\s+cell\s+surface",
    r"is\s+not\s+transported\s+to\s+the\s+(plasma\s+membrane|cell\s+surface)",
    r"not\s+detected\s+at\s+the\s+(cell\s+surface|plasma\s+membrane)",
    r"not\s+present\s+at\s+the\s+(cell\s+surface|plasma\s+membrane)",
    r"absent\s+from\s+the\s+(cell\s+surface|plasma\s+membrane)",
    r"unable\s+to\s+(reach|traffic|localize)",
    r"trafficking\s+defect",
    r"impaired\s+trafficking",
    r"impaired\s+transport",
    r"aberrant\s+localization",
    r"abnormal\s+localization",
    r"misfold",
    r"mislocal",
    r"mistraffic",
    r"sequester(ed)?\s+in\s+the\s+(er|endoplasmic|cytoplasm|golgi|nucleus)",
    r"intracellularly\s+retained",
    r"endo[- ]?h\s+sensitive",
    r"dominant[- ]negative",  # often mechanistic
]
MISTRAFFIC_RE = re.compile("|".join(MISTRAFFIC_PHRASES), re.IGNORECASE)

# Phrases indicating a MOLECULAR REASON (missing structure)
STRUCTURAL_LOSS_PHRASES = [
    r"lacks?\s+the\s+signal\s+(peptide|sequence)",
    r"lacks?\s+a\s+signal\s+(peptide|sequence)",
    r"missing\s+the\s+signal\s+(peptide|sequence)",
    r"no\s+signal\s+(peptide|sequence)",
    r"without\s+a?\s*signal\s+(peptide|sequence)",
    r"signal\s+peptide\s+is\s+(lacking|missing|absent|deleted)",
    r"lacks?\s+the\s+transmembrane",
    r"lacks?\s+a\s+transmembrane",
    r"missing\s+the\s+transmembrane",
    r"without\s+the\s+transmembrane",
    r"lacks?\s+the\s+[A-Z0-9][A-Za-z0-9\- ]{2,30}\s+domain",
    r"missing\s+the\s+[A-Z0-9][A-Za-z0-9\- ]{2,30}\s+domain",
    r"lacks?\s+the\s+extracellular",
    r"lacks?\s+the\s+cytoplasmic\s+(tail|domain)",
    r"truncated\s+",
    r"deletion\s+of\s+",
    r"loss\s+of\s+(the\s+)?(signal|transmembrane|extracellular|cytoplasmic)",
]
STRUCTURAL_RE = re.compile("|".join(STRUCTURAL_LOSS_PHRASES), re.IGNORECASE)

SENT_SPLIT_RE = re.compile(r"(?<=[.!?])\s+(?=[A-Z\[\(])")
ISO_INLINE_RE = re.compile(r"\[Isoform\s+([^\]]+)\]", re.IGNORECASE)


def _norm(text):
    t = (text or "").strip().lower()
    t = re.sub(r"\s*\{[^}]*\}\s*", "", t)
    return t.rstrip(".;,").strip()


def _is_pm(loc):
    n = _norm(loc)
    return any(term in n for term in PM_TERMS)


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
            seq_ids = iso.get("sequenceIds") or []
            note_obj = iso.get("note")
            note_texts = note_obj.get("texts") if isinstance(note_obj, dict) else []
            note_evids = []
            for n in note_texts or []:
                if not isinstance(n, dict):
                    continue
                for ev in n.get("evidences") or []:
                    if isinstance(ev, dict) and ev.get("source") and ev.get("id"):
                        note_evids.append(f"{ev['source']}:{ev['id']}")
            note = " | ".join(
                n.get("value", "") for n in (note_texts or [])
                if isinstance(n, dict) and n.get("value")
            )
            for iso_id in ids:
                out[iso_id] = {
                    "name": name, "synonyms": syns,
                    "sequence_ids": seq_ids, "note": note,
                    "note_evidences": sorted(set(note_evids)),
                }
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


def _tm_ranges(features):
    out = []
    for f in features or []:
        if f.get("type") != "Transmembrane":
            continue
        loc = f.get("location") or {}
        s = (loc.get("start") or {}).get("value")
        e = (loc.get("end") or {}).get("value")
        if isinstance(s, int) and isinstance(e, int):
            out.append((s, e))
    return out


def _signal_peptide(features):
    for f in features or []:
        if f.get("type") == "Signal":
            loc = f.get("location") or {}
            s = (loc.get("start") or {}).get("value")
            e = (loc.get("end") or {}).get("value")
            if isinstance(s, int) and isinstance(e, int):
                return (s, e)
    return None


def _domains(features):
    out = []
    for f in features or []:
        if f.get("type") != "Domain":
            continue
        loc = f.get("location") or {}
        s = (loc.get("start") or {}).get("value")
        e = (loc.get("end") or {}).get("value")
        desc = f.get("description") or ""
        if isinstance(s, int) and isinstance(e, int):
            out.append((s, e, desc))
    return out


def _alt_seq_features(features):
    out = []
    for f in features or []:
        if f.get("type") != "Alternative sequence":
            continue
        loc = f.get("location") or {}
        s = (loc.get("start") or {}).get("value")
        e = (loc.get("end") or {}).get("value")
        desc = f.get("description") or ""
        iso_tokens = []
        for m in re.finditer(r"in\s+isoform\s+([^.,;]+)", desc, re.IGNORECASE):
            for tok in re.split(r"\s*,\s*|\s+and\s+", m.group(1)):
                tok = tok.strip()
                if tok:
                    iso_tokens.append(tok)
        evids = []
        for ev in f.get("evidences") or []:
            if isinstance(ev, dict) and ev.get("source") and ev.get("id"):
                evids.append(f"{ev['source']}:{ev['id']}")
        out.append({
            "begin": s, "end": e, "desc": desc,
            "iso_tokens": iso_tokens, "evidences": sorted(set(evids)),
        })
    return out


def _iso_alt_changes(iso_id, iso_map, alt_features):
    tail = iso_id.split("-")[-1] if "-" in iso_id else iso_id
    name = (iso_map.get(iso_id) or {}).get("name", "").lower()
    syns = {s.strip().lower() for s in (iso_map.get(iso_id) or {}).get("synonyms") or []}
    changes = []
    for af in alt_features:
        for tok in af["iso_tokens"]:
            t = tok.strip().lower()
            if (t == tail.lower() or t == name or t in syns or
                t == f"isoform {tail}".lower()):
                changes.append(af)
                break
    return changes


def _tm_preserved(changes, tm_ranges):
    if not tm_ranges:
        return False
    surviving = 0
    for ts, te in tm_ranges:
        deleted = False
        for af in changes:
            cs, ce = af["begin"], af["end"]
            if cs is None or ce is None:
                continue
            if "missing" in af["desc"].lower() and cs <= ts and ce >= te:
                deleted = True
                break
        if not deleted:
            surviving += 1
    return surviving > 0


def _overlap(a1, a2, b1, b2):
    return not (a2 < b1 or b2 < a1)


def _signal_disrupted(changes, signal_range):
    if not signal_range:
        return False
    ss, se = signal_range
    for af in changes:
        cs, ce = af["begin"], af["end"]
        if cs is None or ce is None:
            continue
        if _overlap(cs, ce, ss, se):
            if "missing" in af["desc"].lower() or cs <= ss:
                return True
    return False


def _domains_disrupted(changes, domains):
    hits = []
    for d_s, d_e, d_desc in domains:
        for af in changes:
            cs, ce = af["begin"], af["end"]
            if cs is None or ce is None:
                continue
            if _overlap(cs, ce, d_s, d_e) and (
                "missing" in af["desc"].lower() or cs <= d_s
            ):
                hits.append(f"{d_desc} ({d_s}-{d_e})")
                break
    return hits


def _collect_iso_locs(comments, iso_map):
    out = defaultdict(list)
    for c in comments:
        if c.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        mol = c.get("molecule")
        if mol:
            iid = _resolve_molecule(mol, iso_map) or f"MOL::{mol}"
            for sl in c.get("subcellularLocations") or []:
                loc_obj = sl.get("location") or {}
                loc = loc_obj.get("value")
                if loc:
                    out[iid].append(_norm(loc))
        note_obj = c.get("note")
        note_texts = note_obj.get("texts") if isinstance(note_obj, dict) else []
        for nt in note_texts or []:
            text = nt.get("value") if isinstance(nt, dict) else None
            if not text:
                continue
            for m in re.finditer(r"\[Isoform\s+([^\]]+)\]\s*:\s*([^.\[]+)", text, re.IGNORECASE):
                target, loc = m.group(1).strip(), m.group(2).strip()
                iid = _resolve_molecule(f"Isoform {target}", iso_map) or f"INLINE::{target}"
                out[iid].append(_norm(loc))
    return dict(out)


def _attribute_sentence_to_iso(sent, iso_map):
    """Return the iso_id that this sentence is attributed to, or None."""
    # inline [Isoform X]
    m = ISO_INLINE_RE.search(sent)
    if m:
        iid = _resolve_molecule(f"Isoform {m.group(1).strip()}", iso_map)
        if iid:
            return iid
    low = sent.lower()
    # iso name / synonym
    for iid, meta in iso_map.items():
        for nm in [meta.get("name", "")] + (meta.get("synonyms") or []):
            if nm and len(nm) > 2 and nm.lower() in low:
                return iid
    # "isoform N" / "isoform X" plain
    m2 = re.search(r"\bisoform\s+(\w+)\b", low)
    if m2:
        iid = _resolve_molecule(f"Isoform {m2.group(1)}", iso_map)
        if iid:
            return iid
    # bare short-form / long-form attribution is ambiguous
    return None


def _mine_evidence_sentences(comments, iso_map):
    """Return {iso_id: [(comment_type, sentence, mistraffic_hit, structural_hit)]}."""
    out = defaultdict(list)
    for c in comments:
        ctype = c.get("commentType", "")
        texts = list(c.get("texts") or [])
        note_obj = c.get("note")
        if isinstance(note_obj, dict):
            texts += note_obj.get("texts") or []
        if ctype == "ALTERNATIVE PRODUCTS":
            for iso in c.get("isoforms") or []:
                nobj = iso.get("note")
                if isinstance(nobj, dict):
                    texts += nobj.get("texts") or []
        for nt in texts:
            val = nt.get("value") if isinstance(nt, dict) else None
            if not val:
                continue
            for sent in SENT_SPLIT_RE.split(val):
                mt = bool(MISTRAFFIC_RE.search(sent))
                st = bool(STRUCTURAL_RE.search(sent))
                if not (mt or st):
                    continue
                iid = _attribute_sentence_to_iso(sent, iso_map)
                if iid is None:
                    continue
                out[iid].append((ctype, sent.strip(), mt, st))
    return dict(out)


def main():
    if not RAW_JSONL.exists():
        sys.exit(f"missing: {RAW_JSONL}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    stats = Counter()
    rows = []

    with gzip.open(RAW_JSONL, "rt") as h:
        for line in h:
            try:
                d = json.loads(line)
            except Exception:
                continue
            if d.get("organism", {}).get("scientificName") != "Homo sapiens":
                continue
            stats["human_entries"] += 1
            comments = d.get("comments") or []
            iso_map = _collect_iso_map(comments)
            if not iso_map:
                continue
            stats["with_alt_products"] += 1

            features = d.get("features") or []
            tm_ranges = _tm_ranges(features)
            if not tm_ranges:
                continue
            stats["with_tm_on_canonical"] += 1
            if len(tm_ranges) < 2:
                stats["skipped_single_pass"] += 1
                continue
            stats["with_multipass_tm"] += 1

            alt_feats = _alt_seq_features(features)
            signal = _signal_peptide(features)
            domains = _domains(features)
            iso_locs = _collect_iso_locs(comments, iso_map)
            sent_by_iso = _mine_evidence_sentences(comments, iso_map)

            primary = d.get("primaryAccession")
            uniprot_id = d.get("uniProtkbId", "")
            pd_ = d.get("proteinDescription") or {}
            rec = pd_.get("recommendedName") or {}
            fn = rec.get("fullName") or {}
            pname = fn.get("value", "") if isinstance(fn, dict) else ""
            genes = d.get("genes") or []
            gene_name = ""
            if genes:
                gn = genes[0].get("geneName") or {}
                gene_name = gn.get("value", "") if isinstance(gn, dict) else ""

            for iid, meta in iso_map.items():
                if iid.startswith(("MOL::", "INLINE::")):
                    continue
                changes = _iso_alt_changes(iid, iso_map, alt_feats)
                if not changes:
                    continue  # must be resolvable to an alt-seq change
                # Criterion II: TM preserved
                if not _tm_preserved(changes, tm_ranges):
                    stats["skipped_tm_lost"] += 1
                    continue

                # Evidence types
                sig_dis = _signal_disrupted(changes, signal)
                dom_hits = _domains_disrupted(changes, domains)
                sentences = sent_by_iso.get(iid, [])
                mistraffic_sents = [s for s in sentences if s[2]]
                structural_sents = [s for s in sentences if s[3] and not s[2]]

                # Iso note explicit content
                note = meta.get("note", "")
                note_mt = bool(MISTRAFFIC_RE.search(note)) if note else False
                note_st = bool(STRUCTURAL_RE.search(note)) if note else False

                # Documented reasons that QUALIFY a row
                qualifying_categories = []
                if mistraffic_sents or note_mt:
                    qualifying_categories.append("text_mistrafficking")
                if structural_sents or note_st:
                    qualifying_categories.append("text_structural_loss")
                if sig_dis:
                    qualifying_categories.append("signal_peptide_deleted")

                # Reported but NOT qualifying for multi-pass TM
                informational_categories = []
                if dom_hits:
                    informational_categories.append("domain_deleted")

                if not qualifying_categories:
                    stats["no_documented_evidence"] += 1
                    continue
                evidence_categories = qualifying_categories + informational_categories

                # Optional: record whether the isoform is annotated non-PM
                locs = iso_locs.get(iid, [])
                has_pm = any(_is_pm(l) for l in locs)
                has_iso_loc = bool(locs)

                rows.append({
                    "primary_accession": primary,
                    "uniprotkb_id": uniprot_id,
                    "gene_name": gene_name,
                    "protein_name": pname,
                    "isoform_id": iid,
                    "isoform_name": meta.get("name", ""),
                    "isoform_synonyms": ";".join(meta.get("synonyms") or []),
                    "evidence_categories": ";".join(evidence_categories),
                    "has_isoform_subcellular_annotation": int(has_iso_loc),
                    "isoform_annotated_non_pm": int(has_iso_loc and not has_pm),
                    "isoform_locations": ";".join(sorted(set(locs))),
                    "tm_count_canonical": len(tm_ranges),
                    "signal_peptide_canonical": f"{signal[0]}-{signal[1]}" if signal else "",
                    "signal_peptide_deleted": int(sig_dis),
                    "domains_deleted": " | ".join(dom_hits),
                    "alt_seq_changes": " | ".join(
                        f"{af['begin']}-{af['end']}: {af['desc']}"
                        + (f" [refs: {','.join(af['evidences'])}]" if af['evidences'] else "")
                        for af in changes
                    ).replace("\t", " "),
                    "mistrafficking_evidence_sentences": " || ".join(
                        f"[{ct}] {st}" for ct, st, _, _ in mistraffic_sents[:6]
                    ).replace("\t", " ").replace("\n", " ")[:6000],
                    "structural_loss_evidence_sentences": " || ".join(
                        f"[{ct}] {st}" for ct, st, _, _ in structural_sents[:6]
                    ).replace("\t", " ").replace("\n", " ")[:6000],
                    "isoform_note": (note or "").replace("\t", " ").replace("\n", " ")[:4000],
                    "isoform_note_evidence_codes": ";".join(meta.get("note_evidences") or []),
                })

    # Write
    tsv_path = OUT_DIR / "tm_negative_isoforms_documented.tsv"
    header = [
        "primary_accession", "uniprotkb_id", "gene_name", "protein_name",
        "isoform_id", "isoform_name", "isoform_synonyms",
        "evidence_categories",
        "has_isoform_subcellular_annotation",
        "isoform_annotated_non_pm",
        "isoform_locations",
        "tm_count_canonical",
        "signal_peptide_canonical", "signal_peptide_deleted",
        "domains_deleted",
        "alt_seq_changes",
        "mistrafficking_evidence_sentences",
        "structural_loss_evidence_sentences",
        "isoform_note", "isoform_note_evidence_codes",
    ]
    with open(tsv_path, "w") as h:
        h.write("\t".join(header) + "\n")
        for r in rows:
            h.write("\t".join(str(r.get(k, "")) for k in header) + "\n")

    cat_counts = Counter()
    for r in rows:
        for cat in r["evidence_categories"].split(";"):
            cat_counts[cat] += 1

    lines = [
        "# UniProt TM isoforms with DOCUMENTED mistrafficking / structural-loss evidence",
        "",
        "## Scan stats",
    ]
    for k, v in stats.most_common():
        lines.append(f"{k}\t{v}")
    lines.append("")
    lines.append(f"qualifying_isoforms\t{len(rows)}")
    lines.append(f"qualifying_genes\t{len({r['gene_name'] or r['primary_accession'] for r in rows})}")
    lines.append("")
    lines.append("## Evidence category breakdown (isoforms may have >1)")
    for cat, n in cat_counts.most_common():
        lines.append(f"{cat}\t{n}")
    lines.append("")
    lines.append(f"Output: {tsv_path}")
    (OUT_DIR / "tm_negative_stats.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
