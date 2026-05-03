"""
generate_pretrain_events.py: Generate synthetic exon-skip events for
Proteus SSL pretraining using SyntheticExonSkipGenerator.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


DEFAULT_DAEDALUS_INTERIM = (
    "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/"
    "Daedalus/phase_a/data/interim"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate synthetic exon-skip events for Proteus SSL pretraining."
    )
    parser.add_argument(
        "--gtf",
        type=Path,
        required=True,
        help="Path to GENCODE GTF file (.gtf or .gtf.gz).",
    )
    parser.add_argument(
        "--daedalus_interim",
        type=Path,
        default=Path(DEFAULT_DAEDALUS_INTERIM),
        help=f"Path to Daedalus Phase A interim directory. Default: {DEFAULT_DAEDALUS_INTERIM}",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=Path("data/interim/pretrain_events"),
        help="Output directory for synthetic events TSV. Default: data/interim/pretrain_events",
    )
    parser.add_argument(
        "--max_genes",
        type=int,
        default=None,
        help="Limit to N genes for debugging. Default: None (all genes)",
    )
    parser.add_argument(
        "--n_workers",
        type=int,
        default=4,
        help="Number of parallel worker processes. Default: 4",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    sys.path.insert(0, str(Path(__file__).parent.parent))
    from proteus.pretraining.synthetic import SyntheticExonSkipGenerator

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "synthetic_exon_skips.tsv"

    print(f"Generating synthetic exon-skip events...")
    print(f"  GTF: {args.gtf}")
    print(f"  Daedalus interim: {args.daedalus_interim}")
    print(f"  Output: {output_path}")
    print(f"  Max genes: {args.max_genes or 'all'}")
    print(f"  Workers: {args.n_workers}")

    gen = SyntheticExonSkipGenerator()
    df = gen.generate(
        gtf_path=args.gtf,
        daedalus_interim=args.daedalus_interim,
        output_path=output_path,
        max_genes=args.max_genes,
        n_workers=args.n_workers,
    )

    print(f"\nGeneration complete: {len(df):,} synthetic pairs saved to {output_path}")


if __name__ == "__main__":
    main()
