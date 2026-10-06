"""Compare verified identical inputs; change only the explicitly authorized embedding."""
import argparse
from importlib.metadata import version
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from altanalyze3.components.clustering.benchmarking.benchmark_landmark_umap import (
    digest, load_verified_input,
)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--method", choices=["legacy", "scanpy", "landmark"], required=True)
    args = parser.parse_args()
    x, cells, features, labels, manifest = load_verified_input(args.input)
    if manifest.get("pca_benchmark_authorized") is not True:
        raise ValueError("PCA embedding comparison requires the user's explicit authorization.")
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    if args.method == "legacy":
        from altanalyze3.components.clustering.ICGS import _import_umap_with_local_retry
        from altanalyze3.components.clustering.umap_fit import fit_umap
        model = _import_umap_with_local_retry().UMAP(**manifest["umap_parameters"])
        coords, info, selected = fit_umap(
            model, lambda rows: np.asarray(x[rows]), len(cells), mode="landmark",
            labels=labels, max_fit_cells=30000, batch_cells=50000, seed=0, log=print)
    else:
        from altanalyze3.components.clustering.accelerated_umap import pca_umap
        matrix = sp.csr_matrix(x)
        coords, info, selected = pca_umap(matrix, method=args.method, labels=labels, log=print)
    seconds = time.perf_counter() - started
    np.save(args.output / "coordinates.npy", coords)
    np.save(args.output / "selected_rows.npy", selected)
    report = {"method": args.method, "embedding_seconds": seconds, "parameters": info,
              "cells": len(cells), "features": len(features),
              "input_manifest_sha256": digest(args.input / "manifest.json"),
              "cell_sha256": digest(args.input / "cells.txt"),
              "feature_sha256": digest(args.input / "features.txt"),
              "cluster_sha256": digest(args.input / "states.npy"),
              "finite_coordinates": bool(np.isfinite(coords).all()),
              "peak_rss_gib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (
                  1024**3 if sys.platform == "darwin" else 1024**2),
              "versions": {name: version(name) for name in (
                  "scanpy", "anndata", "umap-learn", "scikit-learn", "numpy", "scipy")}}
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
