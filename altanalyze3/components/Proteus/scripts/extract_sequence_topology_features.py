"""
extract_sequence_topology_features.py: Compute sequence-based TM topology
and signal peptide predictions for novel AltAnalyze3 isoforms.

Reads a protein FASTA (plain or gzipped) and writes per-protein .npy files
of shape [9] to data/interim/sequence_topology_features/.

Unlike the cached gencode_protein_sequence_features.tsv in Daedalus (which
only covers GENCODE proteins), this script processes any input FASTA including
novel AltAnalyze3 output sequences not present in GENCODE.

CPU-only. Can be called at inference time for novel isoforms.

Usage
-----
# Standard GENCODE proteins:
python scripts/extract_sequence_topology_features.py \
    --fasta ../Daedalus/phase_a/data/raw/gencode/gencode.v49.pc_translations.fa.gz

# Novel AltAnalyze3 isoforms:
python scripts/extract_sequence_topology_features.py \
    --fasta /path/to/altanalyze3_novel_isoforms.fa \
    --cache_dir data/interim/sequence_topology_features
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import OrderedDict
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract sequence-based TM topology features for novel isoforms."
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        required=True,
        help="Path to protein FASTA file (.fa, .fasta, .fa.gz).",
    )
    parser.add_argument(
        "--cache_dir",
        type=Path,
        default=Path("data/interim/sequence_topology_features"),
        help="Output cache directory for .npy files. Default: data/interim/sequence_topology_features",
    )
    parser.add_argument(
        "--n_workers",
        type=int,
        default=4,
        help="Number of parallel worker processes. Default: 4",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Print statistics without running computation.",
    )
    return parser.parse_args()


def _worker_compute(chunk: List[Tuple[str, str]], output_dir: str) -> Tuple[int, int]:
    """Worker: compute and save topology features for a protein chunk."""
    import numpy as np
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from proteus.encoders.topology import SequenceTopologyEncoder

    enc = SequenceTopologyEncoder()
    out = Path(output_dir)
    done, failed = 0, 0
    for pid, seq in chunk:
        f = out / f"{pid}.npy"
        if f.exists():
            done += 1
            continue
        try:
            feats = enc.compute_features(seq)
            np.save(str(f), feats)
            done += 1
        except Exception as exc:
            import warnings
            warnings.warn(f"Failed for {pid}: {exc}")
            np.save(str(f), np.zeros(enc.FEATURE_DIM, dtype=np.float32))
            failed += 1
    return done, failed


def main() -> None:
    from concurrent.futures import ProcessPoolExecutor, as_completed

    args = parse_args()
    cache_dir = args.cache_dir
    cache_dir.mkdir(parents=True, exist_ok=True)

    print(f"Parsing FASTA: {args.fasta}")
    sequences = parse_protein_fasta(args.fasta)
    print(f"  {len(sequences):,} sequences loaded.")

    to_process: List[Tuple[str, str]] = [
        (pid, seq)
        for pid, seq in sequences.items()
        if not (cache_dir / f"{pid}.npy").exists()
    ]
    cached = len(sequences) - len(to_process)
    print(f"  Cached: {cached:,}   To compute: {len(to_process):,}")

    if args.dry_run:
        print("Dry run — exiting.")
        return

    if not to_process:
        print("All features cached.")
        return

    n_workers = min(args.n_workers, len(to_process))
    chunk_size = max(1, len(to_process) // n_workers)
    chunks = [to_process[i: i + chunk_size] for i in range(0, len(to_process), chunk_size)]

    print(f"Processing {len(to_process):,} proteins with {n_workers} workers...")
    total_done, total_failed = 0, 0
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(_worker_compute, chunk, str(cache_dir)): i
            for i, chunk in enumerate(chunks)
        }
        for future in as_completed(futures):
            try:
                done, failed = future.result()
                total_done += done
                total_failed += failed
                print(f"  Chunk done: {done} succeeded, {failed} failed")
            except Exception as exc:
                print(f"  Chunk error: {exc}")

    print(f"\nDone. Succeeded: {total_done}, Failed: {total_failed}")
    print(f"Cache: {cache_dir}")


if __name__ == "__main__":
    main()
