"""Benchmark output options on saved scALABLE artifacts, without rerunning imputation.

Usage: python tests/benchmark_marker_outputs.py JOB_JSON NEW_OUTPUT_DIRECTORY
Outputs must be outside the job. Large benchmark artifacts may be moved to delete/.
Input loading and network export are excluded from marker-stage wall times.
"""
import argparse
import contextlib
import gc
import hashlib
import json
import platform
import time
from pathlib import Path

import anndata as ad
from altanalyze3.components.visualization.marker_heatmap_h5ad import generate_marker_heatmap_from_adata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("job_json", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--modalities", nargs="+", default=["rna", "adt", "metabolite", "lipid", "grn_tf"])
    parser.add_argument("--variants", nargs="+", default=None,
                        choices=["default", "no_static", "markers_only", "25_cells", "600_dpi", "no_svg"])
    args = parser.parse_args()
    source = args.job_json.resolve()
    destination = args.output.resolve()
    if destination == source.parent or source.parent in destination.parents:
        parser.error("Benchmark output must be outside the saved job.")
    destination.mkdir(parents=True, exist_ok=False)
    meta = json.loads(source.read_text())
    cases = {
        "default": {},
        "no_static": {"render_heatmap": False},
        "markers_only": {"render_heatmap": False, "write_heatmap_cache": False},
        "25_cells": {"cells_per_cluster": 25},
        "600_dpi": {"heatmap_dpi": 600},
        "no_svg": {"write_svg": False},
    }
    report = {"job_id": meta["job_id"], "platform": platform.platform(),
              "scope": "Marker scoring and output generation; no imputation, input loading or networks in stage timings",
              "rows": []}
    for modality in args.modalities:
        path = Path(meta["modality_artifacts"][modality]["h5ad"])
        before = (path.stat().st_size, path.stat().st_mtime_ns)
        started = time.perf_counter()
        data = ad.read_h5ad(path)
        load_seconds = time.perf_counter() - started
        analysis = meta["marker_analysis_by_modality"][modality]
        expected = None
        for name in args.variants or cases:
            if not args.variants and modality != "rna" and name not in ("default", "no_static", "markers_only"):
                continue
            overrides = cases[name]
            out = destination / modality / name
            out.mkdir(parents=True)
            opts = dict(cluster_key=analysis["cluster_key"], out=str(out / "markers.pdf"),
                        top_n=50 if modality == "rna" else 5, marker_method="markerfinder",
                        cells_per_cluster=100, seed=0, export_networks=False,
                        validate_scaling=modality == "rna",
                        centroid_method="log2_cp10k" if modality == "rna" else "mean",
                        write_heatmap_tsv=False, write_expression_tsv=False,
                        write_heatmap_cache=True, heatmap_dpi=2400)
            opts.update(overrides)
            with (out / "benchmark.log").open("w") as log, contextlib.redirect_stdout(log):
                started = time.perf_counter()
                result = generate_marker_heatmap_from_adata(data, **opts)
                elapsed = time.perf_counter() - started
            hashes = {k: hashlib.sha256(Path(result[k]).read_bytes()).hexdigest()
                      for k in ("markers_tsv", "redundant_markers_tsv", "centroids_tsv")}
            if expected is None:
                expected = hashes
            if hashes != expected:
                raise AssertionError(f"Marker/centroid results changed for {modality}/{name}")
            row = dict(modality=modality, variant=name, cells=data.n_obs, features=data.n_vars,
                       load_seconds=load_seconds, seconds=elapsed, timings=result["timings"],
                       bytes=sum(p.stat().st_size for p in out.iterdir() if p.name != "benchmark.log"),
                       numerical_outputs_identical=True, hashes=hashes)
            report["rows"].append(row)
            (destination / "results.json").write_text(json.dumps(report, indent=2))
            print(f"{modality}/{name}: {elapsed:.2f}s; {row['bytes']/1e6:.2f} MB; identical marker tables/centroids", flush=True)
            gc.collect()
        assert before == (path.stat().st_size, path.stat().st_mtime_ns)
        del data
        gc.collect()


if __name__ == "__main__":
    main()
