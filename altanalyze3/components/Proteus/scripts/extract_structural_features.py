"""
extract_structural_features.py: Compile structural/functional feature vectors
for all proteins in the Proteus splits from Daedalus Phase A UniProt tables.

CPU-only.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List


DEFAULT_DAEDALUS_INTERIM = (
    "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/"
    "Daedalus/phase_a/data/interim"
)


def get_unique_protein_gene_pairs(data_processed_dir: Path) -> List[tuple]:
    """Return list of (protein_id, gene_id) pairs from all split TSVs."""
    import pandas as pd

    pairs = {}  # protein_id -> gene_id
    for split in ("train", "val", "test"):
        tsv = data_processed_dir / f"proteus_{split}.tsv"
        if not tsv.exists():
            continue
        df = pd.read_csv(tsv, sep="\t", low_memory=False)

        for _, row in df.iterrows():
            gene_id = str(row.get("gene_id", ""))
            for col_prefix in ("reference", "alternative"):
                pid_col = f"{col_prefix}_protein_id"
                if pid_col in df.columns:
                    pid = str(row.get(pid_col, ""))
                    if pid and pid.lower() not in ("nan", "none", ""):
                        pairs[pid] = gene_id

    return [(pid, gid) for pid, gid in pairs.items()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract structural/functional features for all Proteus proteins."
    )
    parser.add_argument(
        "--daedalus_interim",
        type=Path,
        default=Path(DEFAULT_DAEDALUS_INTERIM),
        help=f"Path to Daedalus Phase A interim directory. Default: {DEFAULT_DAEDALUS_INTERIM}",
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
        default=Path("data/interim/structural_features"),
        help="Output cache directory for .pt files. Default: data/interim/structural_features",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Print statistics without running computation.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    cache_dir = args.cache_dir
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Get protein/gene pairs
    pairs = get_unique_protein_gene_pairs(args.data_processed_dir)
    print(f"Unique proteins to process: {len(pairs)}")

    if args.dry_run:
        cached = sum(
            1 for pid, _ in pairs
            if (cache_dir / f"{pid}.pt").exists()
        )
        print(f"Already cached: {cached}")
        print(f"To compute: {len(pairs) - cached}")
        print("Dry run — exiting.")
        return

    sys.path.insert(0, str(Path(__file__).parent.parent))
    from proteus.encoders.structural import StructuralFeatureEncoder

    encoder = StructuralFeatureEncoder()

    protein_ids = [p[0] for p in pairs]
    gene_ids = [p[1] for p in pairs]

    encoder.extract_and_cache(
        protein_ids=protein_ids,
        gene_ids=gene_ids,
        daedalus_interim_dir=args.daedalus_interim,
        output_dir=cache_dir,
    )

    print(f"\nStructural feature extraction complete.")
    cached_count = sum(1 for pid in protein_ids if (cache_dir / f"{pid}.pt").exists())
    print(f"  Total cached: {cached_count}")
    print(f"  Cache directory: {cache_dir}")


if __name__ == "__main__":
    main()
