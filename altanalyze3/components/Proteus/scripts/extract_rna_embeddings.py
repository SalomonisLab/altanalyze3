"""
extract_rna_embeddings.py: Batch extract Orthrus RNA embeddings for all
transcripts in the Proteus train/val/test splits.

Reads proteus_train/val/test.tsv to collect unique transcript IDs,
parses sequences from a GENCODE FASTA file, and calls OrthrusRNAEncoder
to produce cached .pt files in data/interim/rna_embeddings/.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional


def parse_fasta(fasta_path: Path) -> Dict[str, str]:
    """
    Parse a FASTA file (plain or gzipped) into a dict of {id: sequence}.

    Transcript IDs are taken from the first word of the header line.
    """
    sequences: Dict[str, str] = OrderedDict()
    opener = gzip.open if str(fasta_path).endswith(".gz") else open

    current_id: Optional[str] = None
    current_seq: List[str] = []

    with opener(str(fasta_path), "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                header = line[1:].split()[0]
                # Strip version suffix (ENST00000123456.5 → ENST00000123456)
                current_id = header.split(".")[0]
                current_seq = []
            else:
                current_seq.append(line.upper().replace("T", "U"))

    if current_id is not None:
        sequences[current_id] = "".join(current_seq)

    return sequences


def get_unique_transcript_ids(data_processed_dir: Path) -> List[str]:
    """Collect all unique transcript IDs from Proteus split TSVs."""
    import pandas as pd

    ids = set()
    for split in ("train", "val", "test"):
        tsv = data_processed_dir / f"proteus_{split}.tsv"
        if not tsv.exists():
            continue
        df = pd.read_csv(tsv, sep="\t", low_memory=False)
        for col in ("reference_transcript_id", "alternative_transcript_id"):
            if col in df.columns:
                ids.update(df[col].dropna().astype(str).tolist())

    # Remove empty/nan strings
    ids.discard("")
    ids.discard("nan")
    ids.discard("None")
    return sorted(ids)


def find_fasta(data_raw_dir: Path) -> Optional[Path]:
    """Search for GENCODE transcript FASTA in data/raw/gencode/."""
    for pattern in [
        "gencode.*.transcripts.fa.gz",
        "gencode.*.lncRNA_transcripts.fa.gz",
        "gencode.*.pc_transcripts.fa.gz",
        "*.transcripts.fa.gz",
        "*.fa.gz",
        "*.fasta.gz",
        "*.fa",
        "*.fasta",
    ]:
        candidates = sorted((data_raw_dir / "gencode").glob(pattern))
        if candidates:
            return candidates[0]
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract Orthrus RNA embeddings for all Proteus transcripts."
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=None,
        help=(
            "Path to GENCODE transcript FASTA (.fa or .fa.gz). "
            "If not provided, searches data/raw/gencode/."
        ),
    )
    parser.add_argument(
        "--data_processed_dir",
        type=Path,
        default=Path("data/processed"),
        help="Directory containing proteus_train/val/test.tsv. Default: data/processed",
    )
    parser.add_argument(
        "--cache_dir",
        type=Path,
        default=Path("data/interim/rna_embeddings"),
        help="Output cache directory for .pt files. Default: data/interim/rna_embeddings",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=32,
        help="Batch size for encoding. Default: 32",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="auto",
        help="Device (auto/cuda/mps/cpu). Default: auto",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Print statistics without running inference.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Find FASTA
    fasta_path = args.fasta
    if fasta_path is None:
        data_raw = Path("data/raw")
        fasta_path = find_fasta(data_raw)
        if fasta_path is None:
            print("ERROR: No FASTA file found. Please provide --fasta path.")
            sys.exit(1)
    print(f"Using FASTA: {fasta_path}")

    # Get unique transcript IDs
    transcript_ids = get_unique_transcript_ids(args.data_processed_dir)
    print(f"Unique transcript IDs across all splits: {len(transcript_ids)}")

    if args.dry_run:
        # Count already cached
        cache_dir = Path(args.cache_dir)
        cached = sum(1 for tid in transcript_ids if (cache_dir / f"{tid}.pt").exists())
        print(f"Already cached: {cached}")
        print(f"To encode: {len(transcript_ids) - cached}")
        print("Dry run — exiting.")
        return

    # Parse FASTA
    print(f"Parsing FASTA sequences...")
    all_sequences = parse_fasta(fasta_path)
    print(f"  FASTA contains {len(all_sequences)} sequences.")

    # Filter to only needed transcripts
    sequences_needed = {
        tid: all_sequences[tid]
        for tid in transcript_ids
        if tid in all_sequences
    }
    missing_from_fasta = set(transcript_ids) - set(all_sequences.keys())
    print(f"  Found in FASTA: {len(sequences_needed)}")
    if missing_from_fasta:
        print(
            f"  Not found in FASTA (will use zero embeddings): {len(missing_from_fasta)}"
        )
        if len(missing_from_fasta) <= 10:
            print(f"    {sorted(missing_from_fasta)}")

    # Import and run encoder
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from proteus.encoders.rna import OrthrusRNAEncoder

    encoder = OrthrusRNAEncoder(device=args.device)

    ids_list = list(sequences_needed.keys())
    seqs_list = [sequences_needed[tid] for tid in ids_list]

    cached = encoder.extract_and_cache(
        transcript_ids=ids_list,
        sequences=seqs_list,
        output_dir=args.cache_dir,
        device=args.device,
        batch_size=args.batch_size,
    )

    print(f"\nExtraction complete.")
    print(f"  Total cached: {len(cached)}")
    print(f"  Cache directory: {args.cache_dir}")


if __name__ == "__main__":
    main()
