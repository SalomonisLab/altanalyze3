"""Stress bounded statistics with replicated dense H5AD rows, without copying the input.

This is a computational component benchmark, not a complete biological workflow.
HDF5 virtual datasets preserve every source value while expanding the logical matrix.
"""
import argparse
import json
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from altanalyze3.components.cellHarmony.flask.pipeline import _read_differential_h5ad
from altanalyze3.components.cellHarmony.markerFinder import marker_finder
from altanalyze3.components.cellHarmony.cellHarmony_differential import _moderated_t_test
from altanalyze3.components.cellHarmony.webapp.dev.benchmark_full_workflow import monitor


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--cluster-key', required=True)
    parser.add_argument('--cells', type=int, default=1_000_000)
    parser.add_argument('--rss-limit-gib', type=float, default=28)
    args = parser.parse_args()
    if args.cells < 4 or args.output.exists():
        parser.error('Use at least four cells and a fresh output directory')
    try:
        from anndata.io import read_elem
    except ImportError:
        from anndata.experimental import read_elem
    with h5py.File(args.source) as source:
        if not isinstance(source['X'], h5py.Dataset):
            parser.error('This benchmark requires a dense imputed-modality H5AD')
        source_obs, var = read_elem(source['obs']), read_elem(source['var'])
        if args.cluster_key not in source_obs:
            parser.error('Cluster key is missing from the source observations')
        shape, dtype = source['X'].shape, source['X'].dtype
        source_uns = read_elem(source['uns']) if 'uns' in source else {}
    args.output.mkdir(parents=True)
    positions = np.arange(args.cells) % shape[0]
    obs = source_obs.iloc[positions][[args.cluster_key]].copy()
    obs.index = [f'benchmark_{i}' for i in range(args.cells)]
    obs['condition'] = np.where(np.arange(args.cells) % 2, 'case', 'control')
    path = args.output / 'replicated.h5ad'
    uns = {key: source_uns[key] for key in ('expression_scale', 'log1p') if key in source_uns}
    ad.AnnData(X=sp.csr_matrix((args.cells, len(var)), dtype=dtype), obs=obs, var=var,
               uns=uns).write_h5ad(path)
    with h5py.File(path, 'r+') as handle:
        del handle['X']
        layout = h5py.VirtualLayout(shape=(args.cells, len(var)), dtype=dtype)
        source = h5py.VirtualSource(str(args.source.resolve()), 'X', shape=shape)
        for start in range(0, args.cells, shape[0]):
            count = min(shape[0], args.cells - start)
            layout[start:start + count, :] = source[:count, :]
        handle.create_virtual_dataset('X', layout)
        handle['X'].attrs.update({'encoding-type': 'array', 'encoding-version': '0.2.0'})
    measured, stop = monitor(args.output, args.rss_limit_gib)
    disk = _read_differential_h5ad(path, disk_backed=True)
    report = {}
    try:
        started = time.perf_counter()
        r, p = marker_finder(disk.X, disk.obs[args.cluster_key], disk.var_names)
        report.update(markerfinder_seconds=time.perf_counter() - started, marker_shape=r.shape)
        started = time.perf_counter()
        differential, _ = _moderated_t_test(disk, 'condition', 'case', 'control', 'all', store_full_log=False)
        report.update(pooled_differential_seconds=time.perf_counter() - started,
                      differential_rows=len(differential))
    finally:
        disk._analysis_h5_handle.close()
        stop.set()
    report.update(measured, n_cells=args.cells, n_features=len(var),
                  logical_matrix_gib=args.cells * len(var) * np.dtype(dtype).itemsize / 1024**3,
                  source=str(args.source.resolve()),
                  scope='Replicated dense-matrix statistics only; excludes upload, alignment, imputation and web views')
    (args.output / 'results.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
