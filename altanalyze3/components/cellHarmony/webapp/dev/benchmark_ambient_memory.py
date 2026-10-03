"""Measure ambient correction in a fresh process; stdout ends with JSON metrics.

Run before/after separately so ru_maxrss cannot carry over between modes.
This measures correction only, not file import, alignment, or viewer memory.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import subprocess
import sys
import tempfile
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["before", "after"])
    parser.add_argument("--baseline-ref", default="5d4dfd1")
    parser.add_argument("--cells", type=int, default=218000)
    parser.add_argument("--genes", type=int, default=32000)
    parser.add_argument("--nnz-per-cell", type=int, default=80)
    parser.add_argument("--libraries", type=int, default=6)
    args = parser.parse_args()
    if min(args.cells, args.genes, args.nnz_per_cell, args.libraries) < 1:
        parser.error("matrix dimensions and library count must be positive")
    package_root = Path(__file__).resolve().parents[4]
    with tempfile.TemporaryDirectory() as tmp:
        module_path = package_root / "components/ambient_rna/ambient_subtract.py"
        if args.mode == "before":
            code = subprocess.check_output([
                "git", "-C", str(package_root), "show",
                f"{args.baseline_ref}:altanalyze3/components/ambient_rna/ambient_subtract.py",
            ])
            module_path = Path(tmp) / "ambient_before.py"
            module_path.write_bytes(code)
        spec = importlib.util.spec_from_file_location("ambient_benchmark", module_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        rng = np.random.default_rng(19)
        indices = rng.integers(0, args.genes, size=args.cells * args.nnz_per_cell, dtype=np.int32)
        data = rng.integers(1, 20, size=indices.size, dtype=np.int32).astype(np.float32)
        indptr = np.arange(0, indices.size + 1, args.nnz_per_cell, dtype=np.int64)
        x = sp.csr_matrix((data, indices, indptr), shape=(args.cells, args.genes))
        x.sum_duplicates()
        nnz = x.nnz
        obs = pd.DataFrame({"Library": [f"L{i * args.libraries // args.cells}" for i in range(args.cells)]},
                           index=[f"c{i}" for i in range(args.cells)])
        query = ad.AnnData(x, obs=obs)
        del x, data, indices, indptr
        kwargs = dict(rho="auto", outdir=Path(tmp), write_individual=False, write_merged=False)
        if args.mode == "after":
            kwargs.update(inplace=True, store_corrected_layer=False)
        start = time.perf_counter()
        result = module.process_anndata(query, **kwargs)
        elapsed = time.perf_counter() - start
        digest = hashlib.sha256()
        for array in (result.X.indptr, result.X.indices, result.X.data):
            digest.update(memoryview(array))
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        print(json.dumps(dict(
            mode=args.mode, cells=args.cells, genes=args.genes, libraries=args.libraries,
            nnz=nnz, elapsed_sec=round(elapsed, 3),
            peak_rss_mb=round(rss / (1024**2 if sys.platform == "darwin" else 1024), 1),
            corrected_sha256=digest.hexdigest(),
        )))


if __name__ == "__main__":
    main()
