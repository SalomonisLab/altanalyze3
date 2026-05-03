"""
validate_ppi_impact.py: Batch PPI gain/loss analysis for isoform pairs.

Combines three analyses for each reference-alternative pair:
  1. Domain disruption: which UniProt domains are affected
  2. PPI impact: which specific BioGRID partners are predicted lost/retained
  3. Functional category impacts: DNA-binding (TFs), kinase, signaling

Also supports APPRIS reference auto-selection: if only novel transcript IDs
are provided (no explicit reference), the MANE SELECT / APPRIS Principal 1
isoform for that gene is automatically selected as the reference.

Input TSV columns:
  Option A (explicit pairs):
    pair_id, reference_protein_id, reference_gene_name,
    alternative_protein_seq  [optional: reference_protein_seq, junction_peptide]

  Option B (novel transcripts only — APPRIS auto-selection):
    alternative_transcript_id, [alternative_protein_seq], [junction_peptide]
    Requires --daedalus_interim to find the reference.

Output TSV:
  pair_id, reference_gene, n_ref_domains, n_alt_domains_estimated,
  disrupted_domains, retained_domains, n_ppis_lost, n_ppis_retained,
  lost_ppi_partners, retained_ppi_partners, ppi_interface_disrupted,
  dna_binding_disrupted, kinase_domain_disrupted, summary

Usage
-----
# Explicit pairs (reference already known)
python scripts/validate_ppi_impact.py \
    --input  altanalyze3_pairs.tsv \
    --daedalus_interim data/interim \
    --output ppi_impact.tsv

# Novel transcripts only (auto-select APPRIS reference)
python scripts/validate_ppi_impact.py \
    --input  altanalyze3_novel_transcripts.tsv \
    --daedalus_interim data/interim \
    --auto_reference \
    --output ppi_impact.tsv

# With 3DID domain-domain interaction database (higher accuracy)
python scripts/validate_ppi_impact.py \
    --input altanalyze3_pairs.tsv \
    --daedalus_interim data/interim \
    --ddi_flat 3did_flat.gz \
    --output ppi_impact.tsv

# Only report pairs where ≥1 PPI is disrupted
python scripts/validate_ppi_impact.py \
    --input altanalyze3_pairs.tsv \
    --daedalus_interim data/interim \
    --output ppi_impact.tsv \
    --disrupted_only
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Batch PPI impact analysis for isoform pairs")
    p.add_argument("--input", required=True, help="Input TSV of isoform pairs or novel transcripts")
    p.add_argument("--daedalus_interim", required=True, help="Daedalus Phase A data/interim/ directory")
    p.add_argument("--output", required=True, help="Output TSV path")
    p.add_argument("--ddi_flat", default=None,
                   help="Optional 3DID flat file (3did_flat.gz) for domain-domain interactions. "
                        "If absent, DDIs are inferred from BioGRID + UniProt (lower confidence).")
    p.add_argument("--auto_reference", action="store_true",
                   help="Auto-select MANE SELECT / APPRIS reference for each novel transcript. "
                        "Requires alternative_transcript_id column in input.")
    p.add_argument("--pair_id_col", default="pair_id")
    p.add_argument("--ref_protein_col", default="reference_protein_id")
    p.add_argument("--ref_gene_col", default="reference_gene_name",
                   help="Column with reference gene symbol (HGNC). "
                        "Fallback: gene_name column.")
    p.add_argument("--alt_seq_col", default="alternative_protein_seq")
    p.add_argument("--ref_seq_col", default="reference_protein_seq")
    p.add_argument("--junction_col", default="junction_peptide")
    p.add_argument("--has_protein_col", default="has_alt_protein")
    p.add_argument("--alt_transcript_col", default="alternative_transcript_id",
                   help="Used when --auto_reference is set.")
    p.add_argument("--disrupted_only", action="store_true",
                   help="Only output pairs with ≥1 predicted PPI lost or domain disrupted.")
    p.add_argument("--min_lost_ppis", type=int, default=0,
                   help="Minimum number of lost PPIs to include in output (default: 0 = all).")
    return p.parse_args()


def main():
    args = parse_args()

    script_dir = Path(__file__).parent
    proteus_root = script_dir.parent
    if str(proteus_root) not in sys.path:
        sys.path.insert(0, str(proteus_root))

    try:
        import pandas as pd
    except ImportError:
        print("ERROR: pandas required. pip install pandas", file=sys.stderr)
        sys.exit(1)

    from proteus.validation.ppi_impact import PPIImpactAnalyzer

    # ---- Load input ----
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(input_path, sep="\t", low_memory=False)
    print(f"[validate_ppi_impact] Loaded {len(df):,} rows from {input_path.name}")

    interim_dir = Path(args.daedalus_interim)

    # ---- Auto-reference selection ----
    if args.auto_reference:
        from proteus.data.reference_selector import ReferenceSelector
        selector = ReferenceSelector()
        selector.load(interim_dir)

        # Build pairs from novel transcript IDs
        alt_tids = df[args.alt_transcript_col].astype(str).tolist()
        alt_seqs = df[args.alt_seq_col].astype(str).tolist() if args.alt_seq_col in df.columns else [""] * len(df)
        pairs_meta = selector.build_pairs(alt_tids, alt_seqs)
        pairs_df = pd.DataFrame(pairs_meta)

        # Merge with original for sequence / junction columns
        pairs_df = pairs_df.merge(
            df.rename(columns={args.alt_transcript_col: "alternative_transcript_id"}),
            on="alternative_transcript_id",
            how="left",
        )
        print(f"  Auto-selected references for {len(pairs_df):,} pairs.")
    else:
        pairs_df = df.copy()

    # ---- Build pair dicts ----
    pairs = []
    for _, row in pairs_df.iterrows():
        pair_id = str(row.get(args.pair_id_col, f"pair_{len(pairs)}"))
        ref_pid = str(row.get(args.ref_protein_col, "") or "")
        ref_gene = str(
            row.get(args.ref_gene_col, row.get("gene_name", row.get("reference_gene_name", ""))) or ""
        )
        alt_seq = str(row.get(args.alt_seq_col, "") or "")
        ref_seq = str(row.get(args.ref_seq_col, "") or "")
        junction = str(row.get(args.junction_col, "") or "")

        has_prot_raw = row.get(args.has_protein_col, 1)
        if isinstance(has_prot_raw, str):
            has_protein = has_prot_raw.strip().lower() not in ("0", "false", "no", "nmd")
        else:
            try:
                has_protein = bool(int(has_prot_raw))
            except (TypeError, ValueError):
                has_protein = True

        pairs.append({
            "pair_id": pair_id,
            "reference_protein_id": ref_pid,
            "reference_gene_name": ref_gene,
            "alternative_protein_seq": alt_seq,
            "reference_protein_seq": ref_seq,
            "junction_peptide": junction,
            "has_alt_protein": has_protein,
        })

    # ---- Load analyzer and run ----
    analyzer = PPIImpactAnalyzer()
    analyzer.load(
        interim_dir,
        ddi_flat_path=Path(args.ddi_flat) if args.ddi_flat else None,
    )

    print(f"[validate_ppi_impact] Analyzing {len(pairs):,} pairs...")
    results = analyzer.analyze_batch(pairs)

    # ---- Apply filters ----
    if args.disrupted_only:
        results = [r for r in results if r.ppi_interface_disrupted or r.disrupted_domains]
        print(f"  After --disrupted_only filter: {len(results):,} pairs.")

    if args.min_lost_ppis > 0:
        results = [r for r in results if r.n_lost >= args.min_lost_ppis]
        print(f"  After --min_lost_ppis {args.min_lost_ppis} filter: {len(results):,} pairs.")

    # ---- Write output ----
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    analyzer.to_tsv(results, output_path)

    # Summary stats
    n_dna = sum(1 for r in results if r.dna_binding_disrupted)
    n_kinase = sum(1 for r in results if r.kinase_domain_disrupted)
    n_ppi_lost = sum(r.n_lost for r in results)
    n_ppi_retained = sum(r.n_retained for r in results)

    print(f"\n[Summary]")
    print(f"  Pairs analyzed:            {len(results):,}")
    print(f"  Pairs with PPI disruption: {sum(1 for r in results if r.ppi_interface_disrupted):,}")
    print(f"  Total PPIs predicted lost: {n_ppi_lost:,}")
    print(f"  Total PPIs retained:       {n_ppi_retained:,}")
    print(f"  DNA-binding disrupted:     {n_dna:,}")
    print(f"  Kinase domain disrupted:   {n_kinase:,}")
    print(f"  Output: {output_path}")


if __name__ == "__main__":
    main()
