"""Reproducible large synthetic serving benchmark; no biological analysis/model changes.

Prepare once, then run baseline/optimized in separate fresh processes. All source
cells and features, full marker-selection statistics and point values are retained.
This is a serving/component benchmark, not a full workflow or the online visitor job.
"""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import resource
import sys
import time
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from fastapi.responses import JSONResponse

from altanalyze3.components.cellHarmony.webapp.job_bundle import _StoreMatrix
from altanalyze3.components.cellHarmony.webapp.plot_payload import compact_plot_payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["prepare", "baseline", "optimized", "serve"])
    parser.add_argument("root", type=Path)
    parser.add_argument("--cells", type=int, default=600000)
    parser.add_argument("--genes", type=int, default=1024)
    parser.add_argument("--density", type=float, default=.04)
    parser.add_argument("--port", type=int, default=8012)
    args = parser.parse_args()
    w = importlib.import_module("altanalyze3.components.cellHarmony.webapp.app")
    root = args.root.resolve()
    if args.mode == "prepare":
        root.mkdir(parents=True, exist_ok=False)
        rng = np.random.default_rng(71)
        n, genes = args.cells, args.genes
        # Build one gene at a time to keep preparation memory bounded, too.
        nnz_per_gene = max(1, int(n * args.density))
        total = nnz_per_gene * genes
        idx = np.lib.format.open_memmap(root / "indices.npy", mode="w+", dtype=np.int32, shape=(total,))
        data = np.lib.format.open_memmap(root / "data.npy", mode="w+", dtype=np.float32, shape=(total,))
        for j in range(genes):
            cells = np.sort(rng.choice(n, nnz_per_gene, replace=False))
            start = j * nnz_per_gene
            idx[start:start + nnz_per_gene] = cells
            # Distinct group-specific signals and ties; exact values in both modes.
            data[start:start + nnz_per_gene] = rng.random(nnz_per_gene).astype(np.float32) + (cells % 40 == j % 40) * 5
        idx.flush(); data.flush()
        np.save(root / "indptr.npy", np.arange(genes + 1, dtype=np.int64) * nnz_per_gene)
        obs = pd.DataFrame({"state": pd.Categorical([f"c{i % 40}" for i in range(n)]),
                            "Library": pd.Categorical([f"s{i % 127}" for i in range(n)])},
                           index=[f"performance-cell-{i:07}" for i in range(n)])
        coords = rng.normal(size=(n, 2)).astype(np.float32)
        coords[:, 0] += np.arange(n) % 40
        np.save(root / "coords.npy", coords)
        obs.to_pickle(root / "obs.pkl")
        (root / "dimensions.json").write_text(json.dumps(dict(cells=n, genes=genes, nnz=total)))
        print(json.dumps(dict(root=str(root), cells=n, genes=genes, nnz=total)), flush=True)
        return

    dims = json.loads((root / "dimensions.json").read_text())
    n, genes = dims["cells"], dims["genes"]
    idx, data, indptr = [np.load(root / f"{name}.npy", mmap_mode="r") for name in ("indices", "data", "indptr")]
    owner = SimpleNamespace(n_obs=n, n_vars=genes, _indices=idx, _data=data, _indptr=indptr, _to_h5ad=None, _sparse=True)
    owner._dense = lambda j: np.bincount(idx[indptr[j]:indptr[j + 1]], weights=data[indptr[j]:indptr[j + 1]], minlength=n).astype(np.float32)
    matrix = _StoreMatrix(owner)
    obs = pd.read_pickle(root / "obs.pkl")
    # AnnData holds metadata and column reads only; never instantiate n×genes dense.
    a = ad.AnnData(sp.csr_matrix((n, genes), dtype=np.float32), obs=obs,
                   var=pd.DataFrame(index=[f"G{j}" for j in range(genes)]))
    a.__dict__["_X"] = matrix
    a.obsm["X_umap"] = np.load(root / "coords.npy", mmap_mode="r")
    cache = dict(adata=a, obs_names=a.obs_names.to_numpy(), var_names=a.var_names.to_numpy(),
                 populations=obs.state.astype(str).to_numpy(), cluster_key="state", sample_field="Library",
                 sample_labels=obs.Library.astype(str).to_numpy(), umap_x=a.obsm["X_umap"][:, 0], umap_y=a.obsm["X_umap"][:, 1],
                 obs_filter_values={"Library": obs.Library.astype(str).to_numpy(), "state": obs.state.astype(str).to_numpy()})
    cache["display_filters_meta"] = {
        "fields": [{"value": "Library", "label": "Library"}, {"value": "state", "label": "state"}],
        "values": {"Library": obs.Library.cat.categories.tolist(), "state": obs.state.cat.categories.tolist()},
        "default_primary_field": "Library", "default_secondary_field": "state", "show_secondary": True}
    w._get_expression_cache = lambda *args, **kwargs: cache
    w._load_reference_adata = lambda *args: None
    if args.mode == "serve":
        import uvicorn
        application = w.create_app({"JOB_STORAGE": str(root / "jobs"), "APP_TITLE": "scALABLE — synthetic performance test"})
        store = application.state.job_store
        meta = store.create_job("human", "hs_lung_cellref2_reference", None, files=[])
        coords = store.outputs_dir(meta["job_id"]) / "coords.tsv"
        pd.DataFrame({"CellBarcode": a.obs_names, "UMAP1": cache["umap_x"], "UMAP2": cache["umap_y"]}).to_csv(coords, sep="\t", index=False)
        # Header-only source for asynchronous admission; all serving values above
        # come from the complete generated store, including explicit zero entries.
        source = store.outputs_dir(meta["job_id"]) / "performance.h5ad"
        ad.AnnData(sp.csr_matrix((n, genes), dtype=np.float32), obs=obs, var=a.var).write_h5ad(source)
        store.update_job(meta["job_id"], status="completed", cluster_key="state", progress=100,
                         message="Synthetic serving performance test; no biological analysis.",
                         artifacts={"combined_h5ad": str(source), "umap_coordinates": str(coords)},
                         modality_artifacts={"rna": {"h5ad": str(source)}},
                         modalities={"default":"rna", "available":[{"id":"rna","label":"RNA"}]})
        print(f"http://127.0.0.1:{args.port}/?job_id={meta['job_id']}", flush=True)
        uvicorn.run(application, host="127.0.0.1", port=args.port)
        return

    report = dict(mode=args.mode, scope="Synthetic serving components only", **dims)
    if args.mode == "baseline":
        # Same production statistic, using the previous two-pass bundle adapter.
        matrix.group_sums_and_total = None
    start = time.perf_counter()
    selected = w._default_marker_genes(cache)
    report["default_genes_cold_seconds"] = time.perf_counter() - start
    start = time.perf_counter()
    repeated = w._compute_default_marker_genes(cache) if args.mode == "baseline" else w._default_marker_genes(cache)
    report["default_genes_repeat_seconds"] = time.perf_counter() - start
    assert repeated == selected
    report["selected_genes"] = selected
    start = time.perf_counter()
    cells = w._gene_cell_values(cache, selected, cells_per_sample=10)
    report["combplot_cells_seconds"] = time.perf_counter() - start
    report["combplot_cell_count"] = len(cells["columns"])
    report["combplot_sha256"] = hashlib.sha256(JSONResponse(cells).body).hexdigest()
    start = time.perf_counter()
    payload = w._build_umap_payload(None, {}, compact=args.mode == "optimized")
    if args.mode == "baseline":
        payload = compact_plot_payload(payload)
    report["compact_umap_build_seconds"] = time.perf_counter() - start
    start = time.perf_counter()
    body = JSONResponse(payload).body
    report["umap_json_seconds"] = time.perf_counter() - start
    report["umap_bytes"] = len(body)
    report["umap_sha256"] = hashlib.sha256(body).hexdigest()
    report["peak_rss_gib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024**3 if sys.platform == "darwin" else 1024**2)
    (root / f"{args.mode}.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
