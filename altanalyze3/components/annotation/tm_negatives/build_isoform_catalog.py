"""
Build a gene-level UniProt isoform catalog for human + mouse.

For every UniProt reviewed (SwissProt) entry we emit one row per isoform
with the following fields plus the same fields (without redundant
gene/uniprotkb_id columns) for the canonical reference isoform.

Per-isoform fields computed:
    * isoform_locations         molecule-tagged subcellular locations
                                (falls back to entry-level if no
                                molecule-tagged annotation exists)
    * has_cell_membrane         1 if any location mentions 'cell membrane'
                                (covers 'Cell membrane', 'Apical cell
                                membrane', 'Basolateral cell membrane',
                                etc.)
    * has_transmembrane         1 if the isoform retains >= 1 TM helix
                                after applying its alt-seq changes
    * tm_count                  number of canonical TM helices surviving
                                in this isoform
    * signal_peptide            canonical signal-peptide range in the
                                form '1-20', or 'deleted' if the iso's
                                alt-seq change removes it, or '' if
                                canonical has no signal peptide
    * alt_seq_changes           pipe-separated alt-seq events attributed
                                to this isoform (begin-end: description)
    * protein_length            length after applying alt-seq changes

Run:
    python build_isoform_catalog.py
"""

from __future__ import annotations

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
OUT_TSV = OUT_DIR / "uniprot_isoform_catalog.tsv"
OUT_TSV_FULL = OUT_DIR / "uniprot_isoform_catalog_full.tsv"

ORGANISMS = {"Homo sapiens", "Mus musculus"}

ISO_INLINE_RE = re.compile(r"\[Isoform\s+([^\]]+)\]\s*:\s*([^.\[]+)", re.IGNORECASE)


def _norm(text):
    t = (text or "").strip()
    t = re.sub(r"\s*\{[^}]*\}\s*", "", t)
    return t.rstrip(".;,").strip()


def _has_cell_membrane(locs):
    for l in locs:
        ll = l.lower()
        if "cell membrane" in ll or "plasma membrane" in ll:
            return True
    return False


def _collect_iso_map(comments, primary_accession, canonical_length):
    """Return ({iso_id: meta}, canonical_iso_id).

    If the entry has no ALTERNATIVE PRODUCTS comment, we synthesize a
    single-isoform map keyed by the primary accession.
    """
    iso_map = {}
    canonical = None
    for c in comments:
        if c.get("commentType") != "ALTERNATIVE PRODUCTS":
            continue
        for iso in c.get("isoforms") or []:
            ids = iso.get("isoformIds") or []
            if not ids:
                continue
            name_obj = iso.get("name") or {}
            name = name_obj.get("value", "") if isinstance(name_obj, dict) else ""
            syns = [s.get("value", "") for s in iso.get("synonyms") or [] if isinstance(s, dict)]
            status = iso.get("isoformSequenceStatus", "") or iso.get("sequenceStatus", "")
            note_obj = iso.get("note")
            note_texts = note_obj.get("texts") if isinstance(note_obj, dict) else []
            note = " | ".join(
                n.get("value", "") for n in (note_texts or [])
                if isinstance(n, dict) and n.get("value")
            )
            for iso_id in ids:
                iso_map[iso_id] = {
                    "name": name,
                    "synonyms": syns,
                    "status": status,
                    "note": note,
                }
                if str(status).lower() == "displayed" and canonical is None:
                    canonical = iso_id
    if not iso_map:
        iso_map[primary_accession] = {
            "name": "", "synonyms": [], "status": "Displayed", "note": "",
        }
        canonical = primary_accession
    if canonical is None:
        # Fallback: guess canonical as ACC-1 if present, else first key
        guess = f"{primary_accession}-1"
        canonical = guess if guess in iso_map else next(iter(iso_map))
    return iso_map, canonical


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


def _collect_iso_locs(comments, iso_map):
    out = defaultdict(list)
    for c in comments:
        if c.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        mol = c.get("molecule")
        if mol:
            iid = _resolve_molecule(mol, iso_map)
            if iid:
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
            for m in ISO_INLINE_RE.finditer(text):
                target = m.group(1).strip()
                loc = m.group(2).strip()
                iid = _resolve_molecule(f"Isoform {target}", iso_map)
                if iid:
                    out[iid].append(_norm(loc))
    return dict(out)


def _entry_subcell(comments):
    out = []
    for c in comments:
        if c.get("commentType") != "SUBCELLULAR LOCATION" or c.get("molecule"):
            continue
        for sl in c.get("subcellularLocations") or []:
            loc_obj = sl.get("location") or {}
            loc = loc_obj.get("value")
            if loc:
                out.append(_norm(loc))
    return out


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
        alt = f.get("alternativeSequence") or {}
        orig_seq = alt.get("originalSequence") or ""
        alt_seqs = alt.get("alternativeSequences") or []
        alt_seq = alt_seqs[0] if alt_seqs else ""
        out.append({
            "begin": s, "end": e, "desc": desc, "iso_tokens": iso_tokens,
            "original": orig_seq, "alt": alt_seq,
            "is_missing": "missing" in desc.lower(),
        })
    return out


def _iso_changes(iso_id, iso_map, alt_features):
    tail = iso_id.split("-")[-1] if "-" in iso_id else ""
    meta = iso_map.get(iso_id, {})
    name = meta.get("name", "").strip().lower()
    syns = {s.strip().lower() for s in meta.get("synonyms") or []}
    out = []
    for af in alt_features:
        for tok in af["iso_tokens"]:
            t = tok.strip().lower()
            if (t == tail.lower() or t == name or t in syns or
                t == f"isoform {tail}".lower()):
                out.append(af)
                break
    return out


def _apply_len_delta(canonical_len, changes):
    delta = 0
    for af in changes:
        if af["begin"] is None or af["end"] is None:
            continue
        if af["is_missing"]:
            delta -= (af["end"] - af["begin"] + 1)
        else:
            orig_len = af["end"] - af["begin"] + 1
            if af["original"]:
                orig_len = len(af["original"])
            delta += len(af["alt"]) - orig_len
    return max(0, canonical_len + delta)


def _tm_count_after(tm_ranges, changes):
    surviving = 0
    for ts, te in tm_ranges:
        deleted = False
        for af in changes:
            if not af["is_missing"]:
                continue
            cs, ce = af["begin"], af["end"]
            if cs is None or ce is None:
                continue
            if cs <= ts and ce >= te:
                deleted = True
                break
        if not deleted:
            surviving += 1
    return surviving


def _signal_status(signal, changes):
    if not signal:
        return ""
    ss, se = signal
    for af in changes:
        if af["begin"] is None or af["end"] is None:
            continue
        # Disrupted if the isoform is "Missing" over a range that
        # overlaps or precedes the signal peptide.
        if af["is_missing"] and af["begin"] <= ss and af["end"] >= ss:
            return "deleted"
    return f"{ss}-{se}"


def build_row(entry):
    organism = entry.get("organism", {}).get("scientificName", "")
    if organism not in ORGANISMS:
        return []

    primary = entry.get("primaryAccession")
    uniprot_id = entry.get("uniProtkbId", "")
    seq_len = (entry.get("sequence") or {}).get("length") or 0

    genes = entry.get("genes") or []
    gene_name = ""
    if genes:
        gn = genes[0].get("geneName") or {}
        gene_name = gn.get("value", "") if isinstance(gn, dict) else ""

    comments = entry.get("comments") or []
    features = entry.get("features") or []

    iso_map, canonical_iid = _collect_iso_map(comments, primary, seq_len)
    tm_ranges = _tm_ranges(features)
    signal = _signal_peptide(features)
    alt_feats = _alt_seq_features(features)
    iso_locs = _collect_iso_locs(comments, iso_map)
    entry_locs = _entry_subcell(comments)

    # Pre-compute canonical
    can_changes = _iso_changes(canonical_iid, iso_map, alt_feats)
    can_len = seq_len  # canonical is the displayed sequence by definition
    can_tm = _tm_count_after(tm_ranges, [])  # canonical has full set
    can_sig = f"{signal[0]}-{signal[1]}" if signal else ""
    can_locs = iso_locs.get(canonical_iid, []) or entry_locs
    can_alt = " | ".join(
        f"{af['begin']}-{af['end']}: {af['desc']}" for af in can_changes
    )

    rows = []
    for iid, meta in iso_map.items():
        changes = _iso_changes(iid, iso_map, alt_feats) if iid != canonical_iid else []
        iso_len = seq_len if iid == canonical_iid else _apply_len_delta(seq_len, changes)
        iso_tm = can_tm if iid == canonical_iid else _tm_count_after(tm_ranges, changes)
        iso_sig = can_sig if iid == canonical_iid else _signal_status(signal, changes)
        iso_locs_for_iid = iso_locs.get(iid, []) or entry_locs
        iso_alt = " | ".join(
            f"{af['begin']}-{af['end']}: {af['desc']}" for af in changes
        )
        rows.append({
            "organism": organism,
            "gene_name": gene_name,
            "uniprotkb_id": uniprot_id,
            "isoform_id": iid,
            "isoform_name": meta.get("name", ""),
            "isoform_locations": ";".join(sorted(set(iso_locs_for_iid))),
            "has_cell_membrane": int(_has_cell_membrane(iso_locs_for_iid)),
            "has_transmembrane": int(iso_tm > 0),
            "tm_count": iso_tm,
            "signal_peptide": iso_sig,
            "alt_seq_changes": iso_alt,
            "protein_length": iso_len,
            "canonical_isoform_id": canonical_iid,
            "canonical_locations": ";".join(sorted(set(can_locs))),
            "canonical_has_cell_membrane": int(_has_cell_membrane(can_locs)),
            "canonical_has_transmembrane": int(can_tm > 0),
            "canonical_tm_count": can_tm,
            "canonical_signal_peptide": can_sig,
            "canonical_alt_seq_changes": can_alt,
            "canonical_protein_length": can_len,
            "is_canonical": int(iid == canonical_iid),
        })
    return rows


def main():
    if not RAW_JSONL.exists():
        sys.exit(f"missing: {RAW_JSONL}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    header = [
        "organism", "gene_name", "uniprotkb_id",
        "isoform_id", "isoform_name",
        "isoform_locations", "has_cell_membrane",
        "has_transmembrane", "tm_count",
        "signal_peptide", "alt_seq_changes", "protein_length",
        "canonical_isoform_id", "canonical_locations",
        "canonical_has_cell_membrane", "canonical_has_transmembrane",
        "canonical_tm_count", "canonical_signal_peptide",
        "canonical_alt_seq_changes", "canonical_protein_length",
        "is_canonical",
    ]

    n_entries = 0
    n_rows = 0
    n_by_org = defaultdict(int)

    n_full_rows = 0
    with (
        gzip.open(RAW_JSONL, "rt") as inp,
        open(OUT_TSV, "w") as out,
        open(OUT_TSV_FULL, "w") as out_full,
    ):
        out.write("\t".join(header) + "\n")
        out_full.write("\t".join(header) + "\n")
        for line in inp:
            try:
                d = json.loads(line)
            except Exception:
                continue
            rows = build_row(d)
            if not rows:
                continue
            # Gate on a TM+CellMembrane canonical so the catalog only
            # contains isoforms belonging to bona fide surface-protein genes.
            if not rows[0]["canonical_has_transmembrane"]:
                continue
            if not rows[0]["canonical_has_cell_membrane"]:
                continue
            n_entries += 1
            for r in rows:
                full_line = "\t".join(
                    str(r.get(k, "")).replace("\t", " ").replace("\n", " ")
                    for k in header
                ) + "\n"
                out_full.write(full_line)
                n_full_rows += 1
                # Negative catalog (legacy file): non-canonical isoforms
                # that retain TM but lose Cell membrane.
                if r["is_canonical"]:
                    continue
                if r["has_cell_membrane"]:
                    continue
                out.write(full_line)
                n_rows += 1
                n_by_org[r["organism"]] += 1

    print(f"Wrote: {OUT_TSV}")
    print(f"Wrote: {OUT_TSV_FULL}  (all isoforms in TM+CellMembrane genes: {n_full_rows})")
    print(f"Entries processed: {n_entries}")
    print(f"Total negative isoform rows: {n_rows}")
    for org, n in n_by_org.items():
        print(f"  {org}: {n}")


if __name__ == "__main__":
    main()
