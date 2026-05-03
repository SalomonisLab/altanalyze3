#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Daedalus.phase_b.query_inference import DEFAULT_CHECKPOINT, format_prediction_json, load_query_protein_sequence, predict_query_isoform


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Predict reference-conditioned preservation and channel evidence for a query alternative isoform."
    )
    parser.add_argument("--species", required=True, choices=["human", "mouse"])
    parser.add_argument("--gene-id", default=None)
    parser.add_argument("--gene-name", default=None)
    parser.add_argument("--reference-transcript-id", default=None)
    parser.add_argument("--alternative-transcript-id", default="query_alt")
    parser.add_argument("--alt-protein-seq", default=None)
    parser.add_argument("--alt-protein-fasta", default=None)
    parser.add_argument("--alt-transcript-length", type=int, default=None)
    parser.add_argument("--checkpoint", default=str(DEFAULT_CHECKPOINT))
    parser.add_argument("--output-json", default=None)
    parser.add_argument("--include-biogrid-partners", action="store_true", help="Enumerate specific BioGRID partners from the raw BioGRID archive.")
    args = parser.parse_args()

    protein_seq = load_query_protein_sequence(
        path=Path(args.alt_protein_fasta) if args.alt_protein_fasta else None,
        inline_sequence=args.alt_protein_seq,
    )
    result = predict_query_isoform(
        species=args.species,
        gene_id=args.gene_id,
        gene_name=args.gene_name,
        reference_transcript_id=args.reference_transcript_id,
        alternative_transcript_id=args.alternative_transcript_id,
        alt_protein_seq=protein_seq,
        alt_transcript_length=args.alt_transcript_length,
        checkpoint_path=Path(args.checkpoint),
        include_biogrid_partners=args.include_biogrid_partners,
    )
    payload = format_prediction_json(result)
    if args.output_json:
        Path(args.output_json).write_text(payload + "\n", encoding="utf-8")
    print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
