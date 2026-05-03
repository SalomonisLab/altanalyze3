"""
init_proteus.py: Initialize Proteus data directory from Daedalus Phase A pairs.

Reads isoform_pair_candidates.with_splits.tsv and priority_pair_subsets.tsv,
creates Proteus data/processed/ directory, and writes train/val/test TSVs
with embedding cache status columns.
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
        description="Initialize Proteus data directory from Daedalus Phase A isoform pairs."
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
        default=Path("data/processed"),
        help="Output directory for Proteus TSV files. Default: data/processed",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Print statistics without writing files.",
    )
    return parser.parse_args()


def main() -> None:
    import pandas as pd

    args = parse_args()

    daedalus_interim = args.daedalus_interim
    output_dir = args.output_dir

    # Load main pair candidates TSV
    pair_candidates_path = daedalus_interim / "isoform_pair_candidates.with_splits.tsv"
    if not pair_candidates_path.exists():
        print(f"ERROR: Pair candidates file not found: {pair_candidates_path}")
        print("Please ensure Daedalus Phase A has been run first.")
        sys.exit(1)

    print(f"Loading pair candidates from: {pair_candidates_path}")
    df = pd.read_csv(pair_candidates_path, sep="\t", low_memory=False)
    print(f"  Loaded {len(df)} rows with columns: {list(df.columns)}")

    # Load priority pair subsets if available
    priority_path = daedalus_interim / "priority_pair_subsets.tsv"
    priority_df = None
    if priority_path.exists():
        print(f"Loading priority pair subsets from: {priority_path}")
        priority_df = pd.read_csv(priority_path, sep="\t", low_memory=False)
        print(f"  Loaded {len(priority_df)} priority pair rows.")

        # Merge priority annotations into main df — only add columns not already present
        merge_keys = ["reference_transcript_id", "alternative_transcript_id"]
        merge_keys = [k for k in merge_keys if k in priority_df.columns and k in df.columns]

        if merge_keys:
            new_cols = [c for c in priority_df.columns if c not in df.columns]
            if new_cols:
                df = df.merge(
                    priority_df[merge_keys + new_cols],
                    on=merge_keys,
                    how="left",
                )
                print(f"  Merged priority columns: {new_cols}")
            else:
                print("  No new priority columns to merge (already present in main df).")
    else:
        print(f"Priority pair subsets not found at {priority_path}. Skipping.")

    # Add cache status columns (initialized to 0 = not cached)
    df["rna_ref_cached"] = 0
    df["rna_alt_cached"] = 0
    df["protein_ref_cached"] = 0
    df["protein_alt_cached"] = 0

    # Determine split column
    split_col = None
    for candidate in ("split", "data_split", "partition", "set"):
        if candidate in df.columns:
            split_col = candidate
            break

    if split_col is None:
        print("WARNING: No split column found. Assigning all rows to 'train'.")
        df["split"] = "train"
        split_col = "split"
    elif split_col != "split":
        df["split"] = df[split_col]
        print(f"  Renamed split column '{split_col}' -> 'split'")

    # Normalize split values
    df["split"] = df["split"].fillna("train").astype(str).str.lower().str.strip()
    split_map = {
        "train": "train", "training": "train",
        "val": "val", "validation": "val", "valid": "val", "dev": "val",
        "test": "test", "holdout": "test", "eval": "test",
    }
    df["split"] = df["split"].map(lambda x: split_map.get(x, "train"))

    # Split
    train_df = df[df["split"] == "train"].reset_index(drop=True)
    val_df = df[df["split"] == "val"].reset_index(drop=True)
    test_df = df[df["split"] == "test"].reset_index(drop=True)

    print(f"\nSplit statistics:")
    print(f"  Train: {len(train_df):>8,} pairs")
    print(f"  Val:   {len(val_df):>8,} pairs")
    print(f"  Test:  {len(test_df):>8,} pairs")
    print(f"  Total: {len(df):>8,} pairs")

    # Print additional statistics
    for col in ("is_membrane_pair", "is_kinase_pair", "is_tf_pair", "is_surface_pair"):
        if col in df.columns:
            n = df[col].astype(bool).sum()
            print(f"  {col}: {n:,} ({100*n/len(df):.1f}%)")

    if "weak_label" in df.columns:
        label_counts = df["weak_label"].value_counts(dropna=False)
        print(f"\n  Weak label distribution:")
        for label, count in label_counts.items():
            print(f"    {label!r}: {count:,}")

    if args.dry_run:
        print("\nDry run — no files written.")
        return

    # Write output files
    output_dir.mkdir(parents=True, exist_ok=True)

    train_path = output_dir / "proteus_train.tsv"
    val_path = output_dir / "proteus_val.tsv"
    test_path = output_dir / "proteus_test.tsv"

    train_df.to_csv(train_path, sep="\t", index=False)
    val_df.to_csv(val_path, sep="\t", index=False)
    test_df.to_csv(test_path, sep="\t", index=False)

    print(f"\nFiles written to {output_dir}:")
    print(f"  {train_path} ({len(train_df):,} rows)")
    print(f"  {val_path} ({len(val_df):,} rows)")
    print(f"  {test_path} ({len(test_df):,} rows)")
    print("\nInit complete. Next step: run extract_rna_embeddings.py")


if __name__ == "__main__":
    main()
