"""
validate_junction_topology.py

Batch script: given AltAnalyze3 novel isoform output (a TSV of isoform pairs
with junction-spanning peptides), predict the TM topology of each alternative
protein and classify whether the unique junction peptide falls in a
transmembrane, extracellular, cytoplasmic, or signal-peptide region.

This is the correct surfaceome validation for novel isoforms:
  - The input is the *junction-spanning peptide* — the short AA sequence unique
    to the alternative isoform at the novel exon junction
  - We locate that peptide in the full alternative protein and classify its
    position using Kyte-Doolittle sliding-window TM prediction (same algorithm
    as Daedalus build_gencode_protein_sequence_features.py)
  - A junction landing in a TM helix or extracellular loop is a candidate
    cell-surface epitope / CAR-T target

Input TSV columns (required):
    pair_id                 unique identifier for this isoform pair
    alternative_protein_seq full translated AA sequence of the alternative
    junction_peptide        6-30 AA sequence spanning the novel junction
    has_alt_protein         1/0 or True/False (0 = predicted NMD)

Optional columns (passed through to output):
    gene_name, reference_protein_id, alternative_protein_id,
    reference_transcript_id, alternative_transcript_id

Output TSV columns (all input columns + validation results):
    junction_start_aa, junction_end_aa, junction_found,
    n_tm_helices, tm_segments, alt_seq_signal_peptide,
    topology_of_junction, is_surface_accessible, confidence,
    evidence_summary

Usage
-----
# Basic
python scripts/validate_junction_topology.py \
    --input  altanalyze3_novel_junctions.tsv \
    --output junction_topology_validation.tsv

# With protein FASTA (if alternative_protein_seq not in TSV)
python scripts/validate_junction_topology.py \
    --input  altanalyze3_novel_junctions.tsv \
    --fasta  novel_isoform_proteins.fa \
    --output junction_topology_validation.tsv

# Only surface-accessible hits
python scripts/validate_junction_topology.py \
    --input  altanalyze3_novel_junctions.tsv \
    --output junction_topology_validation.tsv \
    --surface_only

AltAnalyze3 integration note
-----------------------------
AltAnalyze3 long-read pipeline produces GFF files and junction BED/TSV files.
The junction_peptide column should contain the 6-30 AA translation of the
reads/bases spanning the novel exon junction that are ABSENT from the reference
protein.  If AltAnalyze3 provides a nucleotide junction sequence, you can
translate it first using --nt_junction_col with --reading_frame 0.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _translate(nt: str, frame: int = 0) -> str:
    """Translate nucleotide sequence in given frame to amino acids."""
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    seq = nt.upper()[frame:]
    aa = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        a = codon_table.get(codon, 'X')
        if a == '*':
            break
        aa.append(a)
    return ''.join(aa)


def _read_fasta(fasta_path: Path) -> dict:
    """Return {sequence_id: sequence} from FASTA file."""
    seqs = {}
    current_id = None
    buf = []
    opener = open
    if str(fasta_path).endswith('.gz'):
        import gzip
        opener = lambda p: gzip.open(p, 'rt')  # noqa: E731
    with opener(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(buf)
                current_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if current_id is not None:
            seqs[current_id] = ''.join(buf)
    return seqs


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Validate AltAnalyze3 junction peptides against predicted TM topology"
    )
    p.add_argument('--input', required=True,
                   help='Input TSV from AltAnalyze3 with pair_id, junction_peptide, '
                        'alternative_protein_seq columns')
    p.add_argument('--output', required=True,
                   help='Output TSV path')
    p.add_argument('--fasta', default=None,
                   help='Optional FASTA of alternative protein sequences (if not in TSV). '
                        'Sequence IDs should match alternative_protein_id column.')
    p.add_argument('--pair_id_col', default='pair_id',
                   help='Column name for pair identifier (default: pair_id)')
    p.add_argument('--seq_col', default='alternative_protein_seq',
                   help='Column containing alternative protein AA sequence (default: alternative_protein_seq). '
                        'If absent, uses --fasta with alternative_protein_id column.')
    p.add_argument('--junction_col', default='junction_peptide',
                   help='Column containing junction-spanning AA peptide (default: junction_peptide)')
    p.add_argument('--nt_junction_col', default=None,
                   help='If set, translate this nucleotide column to get the junction peptide '
                        'instead of --junction_col')
    p.add_argument('--reading_frame', type=int, default=0, choices=[0, 1, 2],
                   help='Reading frame for --nt_junction_col translation (default: 0)')
    p.add_argument('--has_protein_col', default='has_alt_protein',
                   help='Column for NMD flag (0/False = NMD predicted, default: has_alt_protein)')
    p.add_argument('--tm_threshold', type=float, default=1.6,
                   help='Kyte-Doolittle score threshold for TM helix calling (default: 1.6)')
    p.add_argument('--signal_score_min', type=float, default=2.0,
                   help='Min signal peptide score to call as signal (default: 2.0)')
    p.add_argument('--surface_only', action='store_true',
                   help='Only output rows where is_surface_accessible=1')
    p.add_argument('--min_confidence', type=float, default=0.0,
                   help='Minimum topology confidence score to include (default: 0.0)')
    p.add_argument('--min_peptide_len', type=int, default=6,
                   help='Minimum junction peptide length (shorter peptides skipped, default: 6)')
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # ------------------------------------------------------------------ #
    # Import after arg parsing so --help works without torch/pandas        #
    # ------------------------------------------------------------------ #
    try:
        import pandas as pd
    except ImportError:
        print("ERROR: pandas is required. Install with: pip install pandas", file=sys.stderr)
        sys.exit(1)

    # Add Proteus to path if running from scripts/
    script_dir = Path(__file__).parent
    proteus_root = script_dir.parent
    if str(proteus_root) not in sys.path:
        sys.path.insert(0, str(proteus_root))

    from proteus.validation.junction_topology import JunctionTopologyValidator

    # ------------------------------------------------------------------ #
    # Load input                                                           #
    # ------------------------------------------------------------------ #
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(input_path, sep='\t', low_memory=False)
    print(f"[validate_junction_topology] Loaded {len(df):,} rows from {input_path.name}")

    # ------------------------------------------------------------------ #
    # Optionally load FASTA for sequences                                  #
    # ------------------------------------------------------------------ #
    fasta_seqs: dict = {}
    if args.fasta:
        fasta_path = Path(args.fasta)
        if fasta_path.exists():
            fasta_seqs = _read_fasta(fasta_path)
            print(f"  FASTA: {len(fasta_seqs):,} sequences loaded from {fasta_path.name}")
        else:
            print(f"WARNING: FASTA not found: {fasta_path}", file=sys.stderr)

    # ------------------------------------------------------------------ #
    # Build pair list                                                      #
    # ------------------------------------------------------------------ #
    validator = JunctionTopologyValidator(
        tm_threshold=args.tm_threshold,
        signal_score_min=args.signal_score_min,
    )

    pairs = []
    skipped_no_pep = 0
    skipped_short = 0

    for _, row in df.iterrows():
        pair_id = str(row.get(args.pair_id_col, f"row_{len(pairs)}"))

        # Resolve alternative protein sequence
        seq = ""
        if args.seq_col in df.columns:
            seq = str(row.get(args.seq_col, "") or "").strip()
        if not seq and fasta_seqs:
            alt_pid = str(row.get("alternative_protein_id", row.get("alt_protein_id", "")) or "")
            seq = fasta_seqs.get(alt_pid, "")

        # Resolve junction peptide
        if args.nt_junction_col and args.nt_junction_col in df.columns:
            nt = str(row.get(args.nt_junction_col, "") or "").strip()
            junction = _translate(nt, args.reading_frame)
        else:
            junction = str(row.get(args.junction_col, "") or "").strip()

        if not junction:
            skipped_no_pep += 1
            continue
        if len(junction) < args.min_peptide_len:
            skipped_short += 1
            continue

        # NMD flag
        has_protein_raw = row.get(args.has_protein_col, 1)
        if isinstance(has_protein_raw, str):
            has_protein = has_protein_raw.strip().lower() not in ('0', 'false', 'no', 'nmd')
        else:
            try:
                has_protein = bool(int(has_protein_raw))
            except (TypeError, ValueError):
                has_protein = True

        pairs.append({
            "pair_id": pair_id,
            "alternative_protein_seq": seq,
            "junction_peptide": junction,
            "has_alt_protein": has_protein,
        })

    if skipped_no_pep:
        print(f"  Skipped {skipped_no_pep:,} rows with missing junction peptide")
    if skipped_short:
        print(f"  Skipped {skipped_short:,} rows with peptide shorter than {args.min_peptide_len} aa")
    print(f"  Validating {len(pairs):,} pairs...")

    # ------------------------------------------------------------------ #
    # Run validation                                                       #
    # ------------------------------------------------------------------ #
    results = validator.validate(pairs)

    # ------------------------------------------------------------------ #
    # Build output DataFrame                                               #
    # ------------------------------------------------------------------ #
    result_df = validator.to_dataframe(results)

    # Merge pass-through columns from input
    passthrough_cols = [
        c for c in df.columns
        if c not in (
            args.pair_id_col, args.seq_col, args.junction_col,
            args.has_protein_col,
            'junction_start_aa', 'junction_end_aa', 'junction_found',
            'n_tm_helices', 'tm_segments', 'alt_seq_signal_peptide',
            'topology_of_junction', 'is_surface_accessible', 'confidence',
            'evidence_summary',
        )
    ]
    passthrough_cols = [c for c in passthrough_cols if c != 'pair_id']

    if passthrough_cols:
        pass_df = df[[args.pair_id_col] + passthrough_cols].copy()
        pass_df = pass_df.rename(columns={args.pair_id_col: 'pair_id'})
        pass_df['pair_id'] = pass_df['pair_id'].astype(str)
        result_df = result_df.merge(pass_df, on='pair_id', how='left')

    # Apply filters
    if args.surface_only:
        result_df = result_df[result_df['is_surface_accessible'] == 1]
        print(f"  After --surface_only filter: {len(result_df):,} rows")

    if args.min_confidence > 0:
        result_df = result_df[result_df['confidence'] >= args.min_confidence]
        print(f"  After --min_confidence {args.min_confidence} filter: {len(result_df):,} rows")

    # ------------------------------------------------------------------ #
    # Write output                                                         #
    # ------------------------------------------------------------------ #
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(output_path, sep='\t', index=False)

    n_surface = int(result_df['is_surface_accessible'].sum()) if 'is_surface_accessible' in result_df else 0
    n_tm      = int((result_df['topology_of_junction'] == 'transmembrane_helix').sum()) if 'topology_of_junction' in result_df else 0
    n_ec      = int((result_df['topology_of_junction'] == 'extracellular').sum()) if 'topology_of_junction' in result_df else 0
    n_nmd     = int((result_df['topology_of_junction'] == 'nmd_absent').sum()) if 'topology_of_junction' in result_df else 0

    print(f"\n[Results] {len(result_df):,} pairs evaluated")
    print(f"  Surface-accessible (TM + extracellular + signal): {n_surface:,}")
    print(f"  In transmembrane helix:  {n_tm:,}")
    print(f"  In extracellular loop:   {n_ec:,}")
    print(f"  NMD predicted (absent):  {n_nmd:,}")
    print(f"  Written to: {output_path}")


if __name__ == '__main__':
    main()
