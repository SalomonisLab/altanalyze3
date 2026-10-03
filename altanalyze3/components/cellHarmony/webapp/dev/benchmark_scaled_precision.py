"""Measure rounding versus float16 on a saved RNA serving store, in bounded blocks."""
import argparse
import json
import time
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path, help="float32 job_expr_data.npy")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--compact-data", type=Path, help="optional float16 .npy output")
    args = parser.parse_args()
    source = np.load(args.data, mmap_mode="r")
    if source.dtype != np.float32 or source.ndim != 1:
        parser.error("expected a one-dimensional float32 expression value array")
    if args.compact_data and args.compact_data.resolve() == args.data.resolve():
        parser.error("the output must differ from the source")
    if args.compact_data and args.compact_data.exists():
        parser.error("the compact output must be a new file")
    report = {"nonzeros": source.size, "float32_bytes": source.nbytes, "alternatives": {}}
    for name in ("round_3_decimals", "float16"):
        started = time.perf_counter()
        maximum = total = 0.0
        zeros = nonfinite = 0
        dest = None
        if name == "float16" and args.compact_data:
            dest = np.lib.format.open_memmap(args.compact_data, mode="w+", dtype=np.float16,
                                            shape=source.shape)
        for start in range(0, source.size, 4_000_000):
            block = source[start:start + 4_000_000]
            compact = np.round(block, 3) if name == "round_3_decimals" else block.astype(np.float16)
            if dest is not None:
                dest[start:start + block.size] = compact
            restored = compact.astype(np.float32, copy=False)
            delta = np.abs(restored - block)
            maximum = max(maximum, float(delta.max(initial=0)))
            total += float(delta.sum(dtype=np.float64))
            zeros += int(np.count_nonzero((block != 0) & (restored == 0)))
            nonfinite += int(np.count_nonzero(~np.isfinite(restored)))
        if dest is not None:
            dest.flush()
            del dest
        report["alternatives"][name] = {
            "value_bytes": source.size * (2 if name == "float16" else 4),
            "seconds": time.perf_counter() - started,
            "max_absolute_error": maximum,
            "mean_absolute_error": total / max(1, source.size),
            "new_zeros": zeros, "nonfinite": nonfinite,
        }
    args.output.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
