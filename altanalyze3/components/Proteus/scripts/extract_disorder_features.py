"""
extract_disorder_features.py: Batch extract disorder feature vectors
for all proteins in the Proteus train/val/test splits.

CPU-only, parallelizable with --n_workers.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def parse_protein_fasta(fasta_path: Path) -> Dict[str, str]:
    """Parse protein FASTA (plain or gzipped) into {protein_id: sequence}."""
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
                # Preserve full versioned ID (e.g. ENSP00000360644.5)
                current_id = header.split("|")[0].split()[0]
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


def _worker_compute_features(
    chunk: List[Tuple[str, str]],
    output_dir: str,
    backend: str,
) -> Tuple[int, int]:
    """Worker function: compute and save disorder features for a chunk of proteins."""
    import numpy as np
    from pathlib import Path

    # Import disorder encoder in worker
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from proteus.encoders.disorder import DisorderEncoder

    encoder = DisorderEncoder(backend=backend)
    output_path = Path(output_dir)
    done, failed = 0, 0

    for pid, seq in chunk:
        out_file = output_path / f"{pid}.npy"
        if out_file.exists():
            done += 1
            continue
        try:
            feats = encoder.compute_features(seq)
            np.save(str(out_file), feats)
            done += 1
        except Exception as exc:
            import warnings
            warnings.warn(f"Failed for {pid}: {exc}")
            np.save(str(out_file), np.zeros(50, dtype=np.float32))
            failed += 1

    return done, failed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract disorder features for all Proteus proteins."
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        required=True,
        help="Path to protein FASTA file (.fa or .fa.gz).",
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
        default=Path("data/interim/disorder_features"),
        help="Output cache directory for .npy files. Default: data/interim/disorder_features",
    )
    parser.add_argument(
        "--n_workers",
        type=int,
        default=4,
        help="Number of parallel worker processes. Default: 4",
    )
    parser.add_argument(
        "--backend",
        type=str,
        default="auto",
        choices=["auto", "starling", "metapredict"],
        help="Disorder prediction backend. Default: auto (STARLING > metapredict)",
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

    # Get unique protein IDs
    protein_ids = get_unique_protein_ids(args.data_processed_dir)
    print(f"Unique protein IDs across all splits: {len(protein_ids)}")

    # Parse FASTA
    print(f"Parsing protein FASTA: {args.fasta}")
    all_sequences = parse_protein_fasta(args.fasta)
    print(f"  FASTA contains {len(all_sequences)} sequences.")

    to_process: List[Tuple[str, str]] = []
    cached_count = 0

    for pid in protein_ids:
        out_file = cache_dir / f"{pid}.npy"
        if out_file.exists():
            cached_count += 1
        else:
            seq = all_sequences.get(pid, "")
            to_process.append((pid, seq))

    print(f"Already cached: {cached_count}")
    print(f"To compute: {len(to_process)}")

    if args.dry_run:
        print("Dry run — exiting.")
        return

    if not to_process:
        print("All features already cached.")
        return

    # Chunk for parallel processing
    n_workers = min(args.n_workers, len(to_process))
    chunk_size = max(1, len(to_process) // n_workers)
    chunks = [
        to_process[i : i + chunk_size]
        for i in range(0, len(to_process), chunk_size)
    ]

    print(
        f"Processing {len(to_process)} proteins with {n_workers} workers "
        f"({chunk_size} per chunk)..."
    )

    total_done, total_failed = 0, 0

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(
                _worker_compute_features,
                chunk,
                str(cache_dir),
                args.backend,
            ): i
            for i, chunk in enumerate(chunks)
        }
        for future in as_completed(futures):
            try:
                done, failed = future.result()
                total_done += done
                total_failed += failed
                print(f"  Chunk done: {done} succeeded, {failed} failed")
            except Exception as exc:
                print(f"  Chunk failed: {exc}")

    print(f"\nDisorder feature extraction complete.")
    print(f"  Succeeded: {total_done}")
    print(f"  Failed: {total_failed}")
    print(f"  Cache directory: {cache_dir}")


if __name__ == "__main__":
    main()
