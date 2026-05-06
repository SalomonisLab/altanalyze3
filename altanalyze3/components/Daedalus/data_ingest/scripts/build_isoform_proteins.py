#!/usr/bin/env python3
"""Reconstruct full-length protein sequences for every UniProt isoform that
belongs to a TM + Cell-Membrane gene (per the isoform catalog).

For each catalog isoform we apply its `Alternative sequence` (VSP) features to
the canonical sequence to obtain the isoform's full sequence. The output is
the per-isoform sequence table consumed by build_uniprot_isoform_features.py.

Inputs:
    data_ingest/data/raw/uniprot_reviewed_human_mouse.jsonl.gz
    data_ingest/data/interim/tm_negatives/uniprot_isoform_catalog_full.tsv

Output:
    data_ingest/data/interim/uniprot_isoform_proteins.tsv
        species, gene_name, primary_accession, isoform_id,
        canonical_isoform_id, is_canonical, has_cell_membrane,
        has_transmembrane, isoform_locations, protein_length,
        protein_sequence
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
RAW_JSONL = DATA_DIR / "raw" / "uniprot_reviewed_human_mouse.jsonl.gz"
CATALOG_FULL = DATA_DIR / "interim" / "tm_negatives" / "uniprot_isoform_catalog_full.tsv"
OUT_PATH = DATA_DIR / "interim" / "uniprot_isoform_proteins.tsv"


def _norm_species(value: str) -> str:
    v = str(value or "").strip().lower()
    if v in {"homo sapiens", "human", "hs"}:
        return "human"
    if v in {"mus musculus", "mouse", "mm"}:
        return "mouse"
    return v


def _apply_vsp_patches(canonical: str, patches: list[tuple[int, int, str]]) -> str:
    if not patches:
        return canonical
    seq = canonical
    for start, end, replacement in sorted(patches, key=lambda x: -x[0]):
        s, e = start - 1, end
        seq = seq[:s] + replacement + seq[e:]
    return seq


def _extract_isoform_sequences(entry: dict) -> tuple[str | None, dict[str, str]]:
    """Return (canonical_isoform_id, {iso_id: full_sequence})."""
    canonical_seq = entry.get("sequence", {}).get("value", "")
    if not canonical_seq:
        return None, {}

    vsps: dict[str, tuple[int, int, str]] = {}
    for feat in entry.get("features", []) or []:
        if feat.get("type") != "Alternative sequence":
            continue
        fid = feat.get("featureId", "")
        if not fid:
            continue
        loc = feat.get("location", {}) or {}
        try:
            start = int(loc.get("start", {}).get("value"))
            end = int(loc.get("end", {}).get("value"))
        except (TypeError, ValueError):
            continue
        alt = feat.get("alternativeSequence") or {}
        replacement = ""
        if alt.get("alternativeSequences"):
            replacement = str(alt["alternativeSequences"][0])
        elif alt.get("sequence"):
            replacement = str(alt["sequence"])
        vsps[fid] = (start, end, replacement)

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
        canonical_iso_id = f"{acc}-1"
        isoforms[canonical_iso_id] = canonical_seq
    return canonical_iso_id, isoforms


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=OUT_PATH)
    args = ap.parse_args()

    cat = pd.read_csv(CATALOG_FULL, sep="\t")
    cat["species"] = cat["organism"].map(_norm_species)
    cat["primary_accession"] = cat["isoform_id"].astype(str).str.split("-").str[0]
    target_acc = set(cat["primary_accession"].astype(str))
    target_iso_ids = set(cat["isoform_id"].astype(str))
    print(f"[info] target accessions: {len(target_acc)}  isoforms in catalog: {len(target_iso_ids)}")

    # Map for catalog metadata lookups (some isoform_ids span multi-rows; keep first)
    cat_idx = cat.drop_duplicates(subset=["isoform_id"]).set_index("isoform_id")

    rows: list[dict] = []
    with gzip.open(RAW_JSONL, "rt") as f:
        for line in f:
            try:
                d = json.loads(line)
            except Exception:
                continue
            acc = d.get("primaryAccession")
            if acc not in target_acc:
                continue
            canonical_iso, isoforms = _extract_isoform_sequences(d)
            if not isoforms:
                continue
            for iso_id, seq in isoforms.items():
                if iso_id not in cat_idx.index:
                    continue
                meta = cat_idx.loc[iso_id]
                rows.append({
                    "species": _norm_species(meta["organism"]),
                    "gene_name": str(meta["gene_name"]),
                    "primary_accession": acc,
                    "isoform_id": iso_id,
                    "canonical_isoform_id": str(meta["canonical_isoform_id"]),
                    "is_canonical": int(meta["is_canonical"]),
                    "has_cell_membrane": int(meta["has_cell_membrane"]),
                    "has_transmembrane": int(meta["has_transmembrane"]),
                    "tm_count": int(meta["tm_count"]) if pd.notna(meta["tm_count"]) else 0,
                    "isoform_locations": str(meta["isoform_locations"] or ""),
                    "alt_seq_changes": str(meta["alt_seq_changes"] or ""),
                    "protein_length": len(seq),
                    "protein_sequence": seq,
                })

    out_df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"[ok] wrote {args.out}  rows={len(out_df)}")
    print(f"     unique accessions: {out_df['primary_accession'].nunique()}")
    print(f"     unique isoform_ids: {out_df['isoform_id'].nunique()}")
    print(f"     genes: {out_df['gene_name'].nunique()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
