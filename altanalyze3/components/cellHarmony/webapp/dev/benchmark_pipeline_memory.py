"""Reproducible large 10x/h5ad benchmark of import, correction, QC and alignment.

Prepare inputs once, then run each mode in a separate process. This intentionally
does not call MarkerFinder, networks, UMAP placement, or viewer endpoints.
"""
import argparse
import os
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import subprocess
import sys
import tempfile
import threading
import time

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


def load_revision_module(name, relative_path, revision, root, tmp):
    path = Path(tmp) / (name.rsplit(".", 1)[-1] + ".py")
    path.write_bytes(subprocess.check_output(["git", "-C", str(root), "show", f"{revision}:{relative_path}"]))
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare(root, cells, genes, entries, samples):
    root.mkdir(parents=True, exist_ok=True)
    if any(root.glob("sample*.h5")) or (root / "merged.h5ad").exists():
        raise ValueError("Use a fresh fixture directory to avoid mixing benchmark inputs.")
    rng = np.random.default_rng(19)
    names = np.array([f"g{i}".encode() for i in range(genes)])
    for sample, rows in enumerate(np.array_split(np.arange(cells), samples)):
        size = len(rows)
        indices = rng.integers(0, genes, size=size * entries, dtype=np.int32)
        data = rng.integers(1, 20, size=indices.size, dtype=np.int32)
        # Variable sequencing depth exercises count retention and normalization;
        # near-constant totals would be classified as pre-scaled model output.
        data *= np.repeat(rng.integers(1, 10, size=size, dtype=np.int32), entries)
        matrix = sp.csr_matrix((data, indices, np.arange(0, indices.size + 1, entries)), shape=(size, genes))
        matrix.sum_duplicates()
        with h5py.File(root / f"sample{sample}.h5", "w") as handle:
            group = handle.create_group("matrix")
            for key, values in dict(data=matrix.data, indices=matrix.indices, indptr=matrix.indptr).items():
                group.create_dataset(key, data=values, compression="lzf")
            group.create_dataset("shape", data=[genes, size])
            group.create_dataset("barcodes", data=[f"c{i}".encode() for i in rows])
            features = group.create_group("features")
            for key, values in dict(id=names, name=names, feature_type=np.repeat(b"Gene Expression", genes),
                                    genome=np.repeat(b"GRCh38", genes)).items():
                features.create_dataset(key, data=values)
        del matrix, indices, data
    reference = pd.DataFrame(rng.random((min(512, genes), 20)), index=[f"g{i}" for i in range(min(512, genes))],
                             columns=[f"state{i}" for i in range(20)])
    reference.to_csv(root / "reference.tsv", sep="\t")
    (root / "fixture.json").write_text(json.dumps(dict(cells=cells, genes=genes, nnz_per_cell=entries, samples=samples)))


def guard_memory(limit_gib):
    """Terminate the benchmark if peak RSS exceeds the requested local guard."""
    if limit_gib <= 0:
        return
    def watch():
        while True:
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            gib = rss / (1024**3 if sys.platform == "darwin" else 1024**2)
            if gib > limit_gib:
                print(json.dumps(dict(error="memory guard exceeded", peak_rss_gib=gib,
                                      limit_gib=limit_gib)), file=sys.stderr, flush=True)
                os._exit(70)
            time.sleep(0.1)
    threading.Thread(target=watch, daemon=True).start()


def prepare_h5ad(root):
    import scanpy as sc
    from altanalyze3.components.cellHarmony.merge_inputs import concat_matching_10x
    def load(path, sample_name_override=None):
        sample = sc.read_10x_h5(path)
        name = sample_name_override or Path(path).stem
        sample.var_names_make_unique()
        sample.obs_names = [f"{cell}.{name}" for cell in sample.obs_names]
        sample.obs["sample"] = name
        sample.obs["Library"] = name
        sample.obs["group"] = ""
        return sample, name
    merged = concat_matching_10x([str(p) for p in sorted(root.glob("sample*.h5"))], load)
    if merged is None:
        raise ValueError("Prepare matching 10x benchmark inputs first.")
    merged.write_h5ad(root / "merged.h5ad", compression="lzf")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["prepare", "prepare-h5ad", "baseline", "memory", "disk", "stream", "h5ad"])
    parser.add_argument("root", type=Path)
    parser.add_argument("--baseline-ref", default="5d4dfd1")
    parser.add_argument("--cells", type=int, default=218000)
    parser.add_argument("--genes", type=int, default=32000)
    parser.add_argument("--nnz-per-cell", type=int, default=1000)
    parser.add_argument("--samples", type=int, default=6)
    parser.add_argument("--rss-limit-gib", type=float, default=28)
    args = parser.parse_args()
    if min(args.cells, args.genes, args.nnz_per_cell, args.samples) < 1 or args.samples > args.cells:
        parser.error("positive dimensions and at least one cell per sample are required")
    guard_memory(args.rss_limit_gib)
    if args.mode == "prepare":
        prepare(args.root, args.cells, args.genes, args.nnz_per_cell, args.samples)
        return
    if args.mode == "prepare-h5ad":
        prepare_h5ad(args.root)
        return
    from altanalyze3.components.cellHarmony import cellHarmony_lite as module
    package_root = Path(__file__).resolve().parents[4]
    with tempfile.TemporaryDirectory(prefix="cellharmony_bench_") as tmp:
        if args.mode == "baseline":
            import altanalyze3.components.ambient_rna as ambient_package
            legacy = load_revision_module("legacy_ambient", "altanalyze3/components/ambient_rna/ambient_subtract.py",
                                          args.baseline_ref, package_root, tmp)
            ambient_package.ambient_subtract = legacy
            module = load_revision_module("legacy_cellharmony", "altanalyze3/components/cellHarmony/cellHarmony_lite.py",
                                          args.baseline_ref, package_root, tmp)
        disk = args.mode in {"disk", "stream"}
        options = {} if args.mode == "baseline" else dict(ambient_memory_efficient=True,
                    concat_on_disk=disk, concat_batch_size=1 if disk else None, stream_10x_inputs=args.mode == "stream")
        start = time.perf_counter()
        assignments, result = module.combine_and_align_h5(
            h5_files=[] if args.mode == "h5ad" else [str(p) for p in sorted(args.root.glob("sample*.h5"))],
            h5ad_file=str(args.root / "merged.h5ad") if args.mode == "h5ad" else None,
            cellharmony_ref=str(args.root / "reference.tsv"), output_dir=tmp,
            min_genes=0, min_counts=0, min_cells=0, mit_percent=100,
            generate_umap=False, save_adata=False, export_h5ad=False, export_cptt=False,
            ambient_correct_cutoff="auto", return_adata=True, **options,
        )
        elapsed = time.perf_counter() - start
        digest = hashlib.sha256()
        for matrix in (result.X, result.layers.get("counts", result.X), result.layers["soupx_raw"]):
            for values in (matrix.indptr, matrix.indices, matrix.data):
                digest.update(memoryview(values))
        labels = hashlib.sha256(assignments[["CellBarcode", "reference"]].to_csv(index=False).encode()).hexdigest()
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        print(json.dumps(dict(mode=args.mode, elapsed_sec=round(elapsed, 3), cells=result.n_obs, genes=result.n_vars,
            input_files=1 if args.mode == "h5ad" else len(list(args.root.glob("sample*.h5"))),
            input_nnz=int(result.layers["soupx_raw"].nnz),
            peak_rss_mb=round(rss / (1024**2 if sys.platform == "darwin" else 1024), 1),
            counts_layer="counts" in result.layers, normalized="log1p" in result.uns,
            expression_sha256=digest.hexdigest(), assignments_sha256=labels)))


if __name__ == "__main__":
    main()
