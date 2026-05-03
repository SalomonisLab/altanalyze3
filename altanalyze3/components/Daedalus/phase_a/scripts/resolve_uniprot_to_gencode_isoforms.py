#!/usr/bin/env python3
"""Resolve each UniProt isoform suffix (e.g. O14788-2) to a gencode ENSP/ENST
by sequence alignment.

For each gene/accession in the TM-negative catalog:
  1. Reconstruct full protein sequence for every UniProt isoform by applying
     the entry's `Alternative sequence` (VSP) features to the canonical
     sequence.
  2. Collect every gencode ENSP/ENST whose UniProt-Gencode map links it to
     this accession.
  3. For each (isoform, ENSP) pair, score by exact-match first, then by
     pairwise global alignment identity. Pick the top-scoring ENSP per
     isoform with the constraint that an ENSP is assigned at most once.

Gencode is the authoritative transcript source: every gencode protein-coding
transcript on the gene is a candidate, regardless of whether UniProt's
`Ensembl` cross-reference table happens to list it. This avoids missing
matches when UniProt's xrefs are stale.

Inputs:
    phase_a/data/raw/uniprot_reviewed_human_mouse.jsonl.gz
    phase_a/data/raw/gencode/gencode.v49.pc_translations.fa.gz
    phase_a/data/raw/gencode/gencode.vM38.pc_translations.fa.gz
    phase_a/data/interim/tm_negatives/uniprot_isoform_catalog_full.tsv

Output:
    phase_a/data/interim/tm_negatives/uniprot_isoform_to_ensp_resolved.tsv
        columns: species, gene_name, primary_accession, isoform_id,
                 isoform_protein_length, ensp_id, enst_id,
                 ensp_protein_length, identity_pct, length_match,
                 assignment_method
"""
from __future__ import annotations

import argparse
import gzip
import json
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio.Align import PairwiseAligner

PHASE_A = Path(__file__).resolve().parents[1] / "data"
RAW_JSONL = PHASE_A / "raw" / "uniprot_reviewed_human_mouse.jsonl.gz"
GENCODE_HUMAN = PHASE_A / "raw" / "gencode" / "gencode.v49.pc_translations.fa.gz"
GENCODE_MOUSE = PHASE_A / "raw" / "gencode" / "gencode.vM38.pc_translations.fa.gz"
CATALOG_FULL = PHASE_A / "interim" / "tm_negatives" / "uniprot_isoform_catalog_full.tsv"
OUT_PATH = PHASE_A / "interim" / "tm_negatives" / "uniprot_isoform_to_ensp_resolved.tsv"


def _norm_species(value: str) -> str:
    v = str(value or "").strip().lower()
    if v in {"homo sapiens", "human", "hs"}:
        return "human"
    if v in {"mus musculus", "mouse", "mm"}:
        return "mouse"
    return v


def _apply_vsp_patches(canonical: str, patches: list[tuple[int, int, str]]) -> str:
    """Apply VSP patches in coordinate order. Each patch = (start_1based, end_1based, replacement_str).
    An empty replacement means deletion."""
    if not patches:
        return canonical
    # Apply right-to-left so earlier coordinates remain valid
    seq = canonical
    for start, end, replacement in sorted(patches, key=lambda x: -x[0]):
        # 1-based inclusive coordinates
        s, e = start - 1, end
        seq = seq[:s] + replacement + seq[e:]
    return seq


def _extract_isoform_sequences(entry: dict) -> dict[str, str]:
    """Return {isoform_id: full_protein_sequence}. The canonical (Displayed)
    isoform takes the entry's `sequence.value`. Each non-canonical isoform's
    sequence is obtained by applying its VSP features to the canonical."""
    canonical_seq = entry.get("sequence", {}).get("value", "")
    if not canonical_seq:
        return {}

    # Map VSP feature_id -> (start, end, replacement)
    vsps: dict[str, tuple[int, int, str]] = {}
    for feat in entry.get("features", []) or []:
        if feat.get("type") != "Alternative sequence":
            continue
        fid = feat.get("featureId", "")
        if not fid:
            continue
        loc = feat.get("location", {})
        try:
            start = int(loc.get("start", {}).get("value"))
            end = int(loc.get("end", {}).get("value"))
        except (TypeError, ValueError):
            continue
        alt = feat.get("alternativeSequence") or {}
        replacement = alt.get("originalSequence") if False else (alt.get("alternativeSequences") or [""])[0] if alt.get("alternativeSequences") else alt.get("sequence", "") or ""
        vsps[fid] = (start, end, str(replacement or ""))

    isoforms: dict[str, str] = {}
    canonical_iso_id: str | None = None
    for comment in entry.get("comments", []) or []:
        if comment.get("commentType") != "ALTERNATIVE PRODUCTS":
            continue
        for iso in comment.get("isoforms", []) or []:
            iso_ids = iso.get("isoformIds") or []
            if not iso_ids:
                continue
            iso_id = str(iso_ids[0])
            status = (iso.get("isoformSequenceStatus") or "").lower()
            seq_ids = iso.get("sequenceIds") or []
            if status == "displayed":
                canonical_iso_id = iso_id
                isoforms[iso_id] = canonical_seq
            else:
                patches = [vsps[s] for s in seq_ids if s in vsps]
                isoforms[iso_id] = _apply_vsp_patches(canonical_seq, patches)
    if canonical_iso_id is None:
        # Single-isoform entry — emit accession-1 as canonical
        acc = entry.get("primaryAccession", "")
        isoforms[f"{acc}-1"] = canonical_seq
    return isoforms


def _load_uniprot_isoform_sequences(target_accessions: set[str]) -> dict[str, dict[str, str]]:
    """Return {accession: {isoform_id: sequence}} for the requested accessions."""
    out: dict[str, dict[str, str]] = {}
    with gzip.open(RAW_JSONL, "rt") as f:
        for line in f:
            try:
                d = json.loads(line)
            except Exception:
                continue
            acc = d.get("primaryAccession")
            if acc not in target_accessions:
                continue
            seqs = _extract_isoform_sequences(d)
            if seqs:
                out[acc] = seqs
    return out


def _load_gencode_protein_sequences(target_genes_by_species: dict[str, set[str]]) -> dict[tuple[str, str], dict[str, tuple[str, str]]]:
    """Return {(species, gene): {ensp: (enst, sequence)}} for genes of interest."""
    out: dict[tuple[str, str], dict[str, tuple[str, str]]] = defaultdict(dict)
    for species, fasta_path in (("human", GENCODE_HUMAN), ("mouse", GENCODE_MOUSE)):
        if species not in target_genes_by_species:
            continue
        wanted = target_genes_by_species[species]
        if not wanted:
            continue
        with gzip.open(fasta_path, "rt") as f:
            current_ensp = None
            current_enst = None
            current_gene = None
            buf: list[str] = []
            for line in f:
                if line.startswith(">"):
                    # flush previous record
                    if current_ensp is not None and current_gene in wanted:
                        out[(species, current_gene)][current_ensp] = (current_enst, "".join(buf))
                    parts = line.strip()[1:].split("|")
                    current_ensp = parts[0] if parts else None
                    current_enst = parts[1] if len(parts) > 1 else None
                    current_gene = parts[6] if len(parts) > 6 else None
                    buf = []
                else:
                    buf.append(line.strip())
            if current_ensp is not None and current_gene in wanted:
                out[(species, current_gene)][current_ensp] = (current_enst, "".join(buf))
    return dict(out)


def _identity_pct(a: str, b: str, aligner: PairwiseAligner) -> float:
    """Global alignment identity in percent (matches / max(len(a), len(b)))."""
    if not a or not b:
        return 0.0
    if a == b:
        return 100.0
    try:
        score = aligner.score(a, b)
    except Exception:
        return 0.0
    return 100.0 * score / max(len(a), len(b))


def _local_identity_pct(a: str, b: str, aligner: PairwiseAligner) -> float:
    """Local alignment identity = matches / min(len(a), len(b))."""
    if not a or not b:
        return 0.0
    try:
        score = aligner.score(a, b)
    except Exception:
        return 0.0
    return 100.0 * score / max(min(len(a), len(b)), 1)


def _resolve_gene(
    isoforms: dict[str, str],
    ensp_seqs: dict[str, tuple[str, str]],
    global_aligner: PairwiseAligner,
    local_aligner: PairwiseAligner,
    min_identity: float,
    min_local_identity: float,
) -> list[dict]:
    """Greedy assignment: each isoform_id picks its best-scoring ENSP among
    those still unassigned. Score = max(global_identity, local_identity).
    Local identity rescues truncated isoforms (alt seq is a strict subset
    of the canonical / ENSP); global identity preserves preference for
    exact-length matches."""
    rows: list[dict] = []
    if not isoforms or not ensp_seqs:
        return rows
    scored: list[tuple[float, float, float, int, str, str, str]] = []
    for iso_id, iso_seq in isoforms.items():
        for ensp_id, (enst_id, ensp_seq) in ensp_seqs.items():
            if iso_seq == ensp_seq:
                scored.append((100.0, 100.0, 100.0, 0, iso_id, ensp_id, "exact_sequence"))
                continue
            global_id = _identity_pct(iso_seq, ensp_seq, global_aligner)
            local_id = _local_identity_pct(iso_seq, ensp_seq, local_aligner)
            score = max(global_id, local_id)
            length_diff = abs(len(iso_seq) - len(ensp_seq))
            scored.append((score, global_id, local_id, length_diff, iso_id, ensp_id, "pairwise_align"))
    scored.sort(key=lambda x: (-x[0], x[3]))
    used_iso: set[str] = set()
    used_ensp: set[str] = set()
    for score, global_id, local_id, length_diff, iso_id, ensp_id, method in scored:
        if iso_id in used_iso or ensp_id in used_ensp:
            continue
        if score < min_identity and local_id < min_local_identity:
            continue
        used_iso.add(iso_id)
        used_ensp.add(ensp_id)
        enst_id, ensp_seq = ensp_seqs[ensp_id]
        rows.append({
            "isoform_id": iso_id,
            "isoform_protein_length": len(isoforms[iso_id]),
            "ensp_id": ensp_id,
            "enst_id": enst_id,
            "ensp_protein_length": len(ensp_seq),
            "identity_pct": round(global_id, 3),
            "local_identity_pct": round(local_id, 3),
            "length_match": int(len(isoforms[iso_id]) == len(ensp_seq)),
            "assignment_method": method,
        })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scope", default="negatives_and_canonicals",
                    choices=["negatives_and_canonicals", "all_tm_cm"],
                    help="Which accessions to resolve")
    ap.add_argument("--min-identity", type=float, default=80.0,
                    help="Reject assignments below this %% global identity")
    ap.add_argument("--min-local-identity", type=float, default=95.0,
                    help="Accept truncated isoforms with local identity above this %%")
    ap.add_argument("--out", type=Path, default=OUT_PATH)
    args = ap.parse_args()

    cat = pd.read_csv(CATALOG_FULL, sep="\t")
    cat["species"] = cat["organism"].map(_norm_species)
    cat["primary_accession"] = cat["isoform_id"].astype(str).str.split("-").str[0]
    cat["canonical_primary_accession"] = cat["canonical_isoform_id"].astype(str).str.split("-").str[0]

    if args.scope == "negatives_and_canonicals":
        neg = cat[
            (cat["canonical_has_cell_membrane"] == 1)
            & (cat["canonical_has_transmembrane"] == 1)
            & (cat["is_canonical"] == 0)
            & (cat["has_cell_membrane"] == 0)
            & (cat["has_transmembrane"] == 1)
        ]
        target_acc = set(neg["primary_accession"]) | set(neg["canonical_primary_accession"])
    else:
        pos = cat[
            (cat["canonical_has_cell_membrane"] == 1)
            & (cat["canonical_has_transmembrane"] == 1)
            & (cat["has_cell_membrane"] == 1)
            & (cat["has_transmembrane"] == 1)
        ]
        target_acc = set(pos["primary_accession"]) | set(cat["canonical_primary_accession"])

    target_genes_by_species: dict[str, set[str]] = defaultdict(set)
    for _, r in cat[cat["primary_accession"].isin(target_acc)].iterrows():
        target_genes_by_species[r["species"]].add(str(r["gene_name"]))

    print(f"[info] target accessions: {len(target_acc)}")
    print(f"[info] target genes: human={len(target_genes_by_species.get('human', set()))} "
          f"mouse={len(target_genes_by_species.get('mouse', set()))}")

    print("[info] reading UniProt sequences ...")
    iso_seqs = _load_uniprot_isoform_sequences(target_acc)
    n_iso = sum(len(v) for v in iso_seqs.values())
    print(f"[info] loaded {len(iso_seqs)} accessions / {n_iso} isoforms")

    print("[info] reading gencode ENSP sequences ...")
    ensp_seqs = _load_gencode_protein_sequences(target_genes_by_species)
    n_ensp = sum(len(v) for v in ensp_seqs.values())
    print(f"[info] loaded {len(ensp_seqs)} (species, gene) groups / {n_ensp} ENSPs")

    # Map gene -> accessions in catalog (for pairing genes with the right UniProt entry)
    gene_to_acc: dict[tuple[str, str], set[str]] = defaultdict(set)
    for _, r in cat.iterrows():
        gene_to_acc[(r["species"], str(r["gene_name"]))].add(str(r["primary_accession"]))

    global_aligner = PairwiseAligner()
    global_aligner.mode = "global"
    global_aligner.match_score = 1.0
    global_aligner.mismatch_score = 0.0
    global_aligner.open_gap_score = 0.0
    global_aligner.extend_gap_score = 0.0

    local_aligner = PairwiseAligner()
    local_aligner.mode = "local"
    local_aligner.match_score = 1.0
    local_aligner.mismatch_score = -1.0
    local_aligner.open_gap_score = -2.0
    local_aligner.extend_gap_score = -1.0

    rows: list[dict] = []
    for (species, gene), accessions in sorted(gene_to_acc.items()):
        ensps = ensp_seqs.get((species, gene), {})
        if not ensps:
            continue
        merged_iso: dict[str, str] = {}
        for acc in accessions:
            for iso_id, seq in iso_seqs.get(acc, {}).items():
                merged_iso[iso_id] = seq
        if not merged_iso:
            continue
        for r in _resolve_gene(
            merged_iso, ensps, global_aligner, local_aligner,
            args.min_identity, args.min_local_identity,
        ):
            rows.append({
                "species": species,
                "gene_name": gene,
                "primary_accession": str(r["isoform_id"]).split("-")[0],
                **r,
            })

    out_df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"[ok] wrote {args.out}  rows={len(out_df)}")
    if not out_df.empty:
        print("    exact_sequence:", int((out_df["assignment_method"] == "exact_sequence").sum()))
        print("    pairwise_align:", int((out_df["assignment_method"] == "pairwise_align").sum()))
        print(f"    median identity_pct: {out_df['identity_pct'].median():.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
