"""
extract_protein_embeddings.py: Batch extract ESM2 protein embeddings for all
proteins in the Proteus train/val/test splits.

Reads proteus_train/val/test.tsv to collect unique protein IDs,
parses sequences from a GENCODE protein FASTA file, and calls
ESM2ProteinEncoder to produce cached .pt or .absent files in
data/interim/protein_embeddings/.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional


def parse_fasta(fasta_path: Path, strip_version: bool = True) -> Dict[str, str]:
    """
    Parse a protein FASTA file into {protein_id: sequence}.

    For GENCODE protein FASTAs, headers look like:
    >ENSP00000123456.3|ENST00000123456.5|ENSG00000123456.10|...

    The protein ID is taken as the first pipe-delimited field,
    with optional version stripping.
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
                header = line[1:]
                # Try pipe-delimited (GENCODE protein format)
                first_field = header.split("|")[0].split()[0]
                if strip_version and "." in first_field:
                    pid = first_field.split(".")[0]
                else:
                    pid = first_field
                current_id = pid
                current_seq = []
            else:
                current_seq.append(line.upper())

    if current_id is not None:
        sequences[current_id] = "".join(current_seq)

    return sequences


def get_unique_protein_ids(data_processed_dir: Path) -> List[str]:
    """Collect all unique protein IDs from Proteus split TSVs."""
    import pandas as pd

    ids = set()
    for split in ("train", "val", "test"):
        tsv = data_processed_dir / f"proteus_{split}.tsv"
        if not tsv.exists():
            continue
        df = pd.read_csv(tsv, sep="\t", low_memory=False)
        for col in ("reference_protein_id", "alternative_protein_id"):
            if col in df.columns:
                ids.update(df[col].dropna().astype(str).tolist())

    ids.discard("")
    ids.discard("nan")
    ids.discard("None")
    return sorted(ids)


def find_fasta(data_raw_dir: Path) -> Optional[Path]:
    """Search for GENCODE protein FASTA in data/raw/gencode/."""
    for pattern in [
        "gencode.*.pc_translations.fa.gz",
        "gencode.*.proteins.fa.gz",
        "*translations*.fa.gz",
        "*proteins*.fa.gz",
        "*.fa.gz",
        "*.fasta.gz",
    ]:
        subdir = data_raw_dir / "gencode"
        if subdir.exists():
            candidates = sorted(subdir.glob(pattern))
            if candidates:
                return candidates[0]
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract ESM2 protein embeddings for all Proteus proteins."
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=None,
        help=(
            "Path to GENCODE protein FASTA (.fa or .fa.gz). "
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
        default=Path("data/interim/protein_embeddings"),
        help="Output cache directory for .pt/.absent files. Default: data/interim/protein_embeddings",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=16,
        help="Batch size for encoding. Default: 16",
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
            print("ERROR: No protein FASTA file found. Please provide --fasta path.")
            sys.exit(1)
    print(f"Using FASTA: {fasta_path}")

    # Get unique protein IDs
    protein_ids = get_unique_protein_ids(args.data_processed_dir)
    print(f"Unique protein IDs across all splits: {len(protein_ids)}")

    if args.dry_run:
        cache_dir = Path(args.cache_dir)
        cached = sum(
            1 for pid in protein_ids
            if (cache_dir / f"{pid}.pt").exists() or (cache_dir / f"{pid}.absent").exists()
        )
        print(f"Already cached: {cached}")
        print(f"To encode: {len(protein_ids) - cached}")
        print("Dry run — exiting.")
        return

    # Parse FASTA
    print("Parsing protein FASTA sequences...")
    all_sequences = parse_fasta(fasta_path)
    print(f"  FASTA contains {len(all_sequences)} sequences.")

    sequences_map: Dict[str, Optional[str]] = {}
    for pid in protein_ids:
        sequences_map[pid] = all_sequences.get(pid, None)

    found = sum(1 for v in sequences_map.values() if v)
    missing = len(sequences_map) - found
    print(f"  Found in FASTA: {found}")
    print(f"  Not found (will be marked absent): {missing}")

    # Import and run encoder
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from proteus.encoders.protein import ESM2ProteinEncoder

    encoder = ESM2ProteinEncoder(device=args.device)

    ids_list = list(sequences_map.keys())
    seqs_list = [sequences_map[pid] for pid in ids_list]

    cached = encoder.extract_and_cache(
        protein_ids=ids_list,
        sequences=seqs_list,
        output_dir=args.cache_dir,
        device=args.device,
        batch_size=args.batch_size,
    )

    print(f"\nExtraction complete.")
    pt_count = sum(1 for p in cached.values() if str(p).endswith(".pt"))
    absent_count = sum(1 for p in cached.values() if str(p).endswith(".absent"))
    print(f"  Embedded (.pt): {pt_count}")
    print(f"  Absent (.absent): {absent_count}")
    print(f"  Cache directory: {args.cache_dir}")


if __name__ == "__main__":
    main()
