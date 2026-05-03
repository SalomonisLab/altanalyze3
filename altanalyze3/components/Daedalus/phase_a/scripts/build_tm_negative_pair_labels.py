#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


PHASE_A = Path(__file__).resolve().parents[1]
INTERIM = PHASE_A / "data" / "interim"
TM_NEG = INTERIM / "tm_negatives"

NEG_PATH = TM_NEG / "tm_negative_isoforms_documented.tsv"
LIT_PATH = TM_NEG / "gene_literature_evidence.tsv"
MAP_PATH = INTERIM / "uniprot_gencode_map.tsv"
PAIR_PATH = INTERIM / "isoform_pair_candidates.with_splits.tsv"

OUT_TSV = TM_NEG / "tm_negative_pair_labels.tsv"
OUT_JSON = TM_NEG / "tm_negative_pair_label_summary.json"


def _normalize_species(value: str) -> str:
    value = str(value or "").strip().lower()
    if value in {"homo sapiens", "human", "hs"}:
        return "human"
    if value in {"mus musculus", "mouse", "mm"}:
        return "mouse"
    return value


def _parse_shared(value: object) -> set[str]:
    if pd.isna(value):
        return set()
    return {token for token in str(value).split(";") if token}


def main() -> int:
    neg = pd.read_csv(NEG_PATH, sep="\t")
    lit = pd.read_csv(LIT_PATH, sep="\t")
    mapping = pd.read_csv(MAP_PATH, sep="\t")
    pairs = pd.read_csv(PAIR_PATH, sep="\t", low_memory=False)

    lit["species"] = lit["organism"].map(_normalize_species)
    lit = lit.rename(
        columns={
            "confidence": "literature_confidence",
            "evidence_found": "literature_evidence_found",
            "source_urls": "literature_source_urls",
            "evidence_snippet": "literature_evidence_snippet",
        }
    )
    lit = lit[
        [
            "gene_name",
            "species",
            "literature_evidence_found",
            "literature_confidence",
            "pubmed_ids",
            "literature_source_urls",
            "literature_evidence_snippet",
        ]
    ].drop_duplicates()

    mapping["species"] = mapping["gencode_species"].map(_normalize_species)
    acc_to_tx = (
        mapping.dropna(subset=["primary_accession", "gencode_transcript_id"])
        .groupby(["primary_accession", "species"])["gencode_transcript_id"]
        .apply(lambda s: sorted(set(str(x) for x in s if str(x))))
        .to_dict()
    )

    rows: list[dict[str, object]] = []
    pairs["species"] = pairs["species"].map(_normalize_species)

    for neg_row in neg.itertuples(index=False):
        gene_name = str(neg_row.gene_name)
        accession = str(neg_row.primary_accession)
        species_candidates = [
            species
            for (acc, species) in acc_to_tx.keys()
            if acc == accession
        ]
        for species in species_candidates:
            mapped_txs = set(acc_to_tx.get((accession, species), []))
            if not mapped_txs:
                continue
            sub = pairs[
                (pairs["species"] == species)
                & (pairs["gene_name"].astype(str) == gene_name)
                & (pairs["alternative_transcript_id"].astype(str).isin(mapped_txs))
            ].copy()
            if sub.empty:
                continue
            sub = sub[
                sub["shared_uniprot_accessions"].map(lambda v: accession in _parse_shared(v))
            ].copy()
            if sub.empty:
                continue
            lit_sub = lit[(lit["gene_name"].astype(str) == gene_name) & (lit["species"] == species)]
            lit_row = lit_sub.iloc[0].to_dict() if not lit_sub.empty else {}
            for pair_row in sub.itertuples(index=False):
                rows.append(
                    {
                        "species": species,
                        "gene_name": gene_name,
                        "primary_accession": accession,
                        "uniprotkb_id": getattr(neg_row, "uniprotkb_id", ""),
                        "reference_transcript_id": str(pair_row.reference_transcript_id),
                        "alternative_transcript_id": str(pair_row.alternative_transcript_id),
                        "split": str(pair_row.split),
                        "label_task": "cell_surface_localization",
                        "label_binary": 0,
                        "label_state": "negative",
                        "label_source": "curated_tm_negative",
                        "evidence_categories": getattr(neg_row, "evidence_categories", ""),
                        "has_isoform_subcellular_annotation": getattr(neg_row, "has_isoform_subcellular_annotation", ""),
                        "isoform_annotated_non_pm": getattr(neg_row, "isoform_annotated_non_pm", ""),
                        "isoform_locations": getattr(neg_row, "isoform_locations", ""),
                        "tm_count_canonical": getattr(neg_row, "tm_count_canonical", ""),
                        "signal_peptide_canonical": getattr(neg_row, "signal_peptide_canonical", ""),
                        "signal_peptide_deleted": getattr(neg_row, "signal_peptide_deleted", ""),
                        "domains_deleted": getattr(neg_row, "domains_deleted", ""),
                        "alt_seq_changes": getattr(neg_row, "alt_seq_changes", ""),
                        "mistrafficking_evidence_sentences": getattr(neg_row, "mistrafficking_evidence_sentences", ""),
                        "structural_loss_evidence_sentences": getattr(neg_row, "structural_loss_evidence_sentences", ""),
                        "isoform_note": getattr(neg_row, "isoform_note", ""),
                        "literature_evidence_found": lit_row.get("literature_evidence_found", ""),
                        "literature_confidence": lit_row.get("literature_confidence", ""),
                        "pubmed_ids": lit_row.get("pubmed_ids", ""),
                        "literature_source_urls": lit_row.get("literature_source_urls", ""),
                        "literature_evidence_snippet": lit_row.get("literature_evidence_snippet", ""),
                    }
                )

    out = pd.DataFrame(rows).drop_duplicates(
        subset=[
            "species",
            "gene_name",
            "primary_accession",
            "reference_transcript_id",
            "alternative_transcript_id",
            "label_task",
        ]
    )
    out = out.sort_values(
        ["species", "gene_name", "primary_accession", "reference_transcript_id", "alternative_transcript_id"]
    ).reset_index(drop=True)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)

    summary = {
        "rows": int(len(out)),
        "genes": int(out["gene_name"].nunique()) if not out.empty else 0,
        "species_counts": out["species"].value_counts().to_dict() if not out.empty else {},
        "split_counts": out["split"].value_counts().to_dict() if not out.empty else {},
        "source_files": {
            "negatives": str(NEG_PATH),
            "literature": str(LIT_PATH),
            "uniprot_gencode_map": str(MAP_PATH),
            "pair_candidates": str(PAIR_PATH),
        },
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"[ok] Wrote {OUT_TSV} rows={len(out)}")
    print(f"[ok] Wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
