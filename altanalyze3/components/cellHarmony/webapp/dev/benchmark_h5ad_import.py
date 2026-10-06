"""Synthetic import/QC/alignment/MarkerFinder stress test; no biological claims.

Run in a disposable process. The exact full feature/cell roster is generated here;
no user data or production models are substituted. Large scratch files are removed.
"""
import argparse
import json
import resource
import sys
import shutil
import tempfile
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from altanalyze3.components.cellHarmony.cellHarmony_lite import combine_and_align_h5


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cells', type=int, default=400_000)
    parser.add_argument('--nnz-per-cell', type=int, default=1000)
    parser.add_argument('--files', type=int, default=1)
    parser.add_argument('--min-counts', type=int, default=500)
    parser.add_argument('--mode', choices=['baseline', 'bounded'], default='bounded')
    parser.add_argument('--markers', action='store_true')
    parser.add_argument('--umap', action='store_true', help='Also run approximate UMAP and save the H5AD.')
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    n, p, k = args.cells, 32_738, args.nnz_per_cell
    assert 32 < k < p and n > 0 and 1 <= args.files <= n
    if args.files > 1 and args.mode != 'bounded':
        parser.error('Multi-file validation uses the new disk merge backend.')
    result = {'synthetic': True, 'cells': n, 'genes': p, 'stored_values': n * k,
              'mode': args.mode, 'scope': 'import, configured QC, normalization, cosine alignment; optional standard MarkerFinder'}
    with tempfile.TemporaryDirectory(prefix='scalable-h5ad-benchmark-') as scratch:
        root = Path(scratch)
        required = n * k * ((28 if args.files > 1 else 20) if args.mode == 'bounded' else 8) + 3 * 1024**3
        if shutil.disk_usage(root).free < required:
            raise RuntimeError(f'Benchmark needs at least {required / 1024**3:.1f} GiB free scratch space.')
        path = root / 'synthetic.h5ad'
        template_rows = 1024
        variable = (np.arange(template_rows)[:, None] * 997 + np.arange(k - 32)) % (p - 32) + 32
        indices = np.concatenate([np.tile(np.arange(32), (template_rows, 1)), np.sort(variable, axis=1)], axis=1).astype(np.int32)
        values = np.random.default_rng(71).integers(1, 11, (template_rows, k)).astype(np.float32)
        values[np.arange(template_rows), np.arange(template_rows) % 2] = 100
        obs = pd.DataFrame({'Library': pd.Categorical(np.where(np.arange(n) % 2, 'B', 'A'))},
                           index=[f'cell_{i}' for i in range(n)])
        var = pd.DataFrame(index=[f'G{i}' for i in range(p)])
        print(f'Preparing {args.files} synthetic H5AD file(s): {n:,} cells, {n * k:,} stored values.', flush=True)
        entries, expected_ids = [], []
        boundaries = np.linspace(0, n, args.files + 1, dtype=int)
        for i, (lo, hi) in enumerate(zip(boundaries[:-1], boundaries[1:])):
            path = root / f'synthetic-{i}.h5ad'
            entries.append((path, f'upload{i}'))
            expected_ids.extend(obs.index[lo:hi] if args.files == 1 else
                                obs.index[lo:hi] + f'::upload{i}')
            ad.AnnData(X=sp.csr_matrix((hi - lo, p), dtype=np.float32),
                       obs=obs.iloc[lo:hi].copy(), var=var).write_h5ad(path)
            with h5py.File(path, 'r+') as f:
                x = f['X']
                for key in ('data', 'indices', 'indptr'):
                    del x[key]
                data = x.create_dataset('data', ((hi - lo) * k,), dtype=np.float32, compression='lzf')
                cols = x.create_dataset('indices', ((hi - lo) * k,), dtype=np.int32, compression='lzf')
                x.create_dataset('indptr', data=np.arange(hi - lo + 1, dtype=np.int64) * k)
                for start in range(0, hi - lo, template_rows):
                    count = min(template_rows, hi - lo - start)
                    template = (np.arange(count) + lo + start) % template_rows
                    data[start * k:(start + count) * k] = values[template].ravel()
                    cols[start * k:(start + count) * k] = indices[template].ravel()
        print('Synthetic inputs prepared; starting measured analysis.', flush=True)
        reference = root / 'ref.tsv'
        ref = pd.DataFrame(np.ones((32, 2)), index=var.index[:32], columns=['A', 'B'])
        ref.iloc[0, 0] = ref.iloc[1, 1] = 100
        ref.to_csv(reference, sep='\t')
        started = time.perf_counter()
        if args.files > 1:
            from altanalyze3.components.cellHarmony.mapped_h5ad import merge_h5ads
            path = merge_h5ads(entries, root / 'merge')
            result['merge_seconds'] = time.perf_counter() - started
        result['files'] = args.files
        _, output = combine_and_align_h5([], reference, h5ad_file=path, output_dir=root / 'output',
                                         min_genes=500, min_cells=0, min_counts=args.min_counts, mit_percent=15,
                                         min_alignment_score=None, return_adata=True,
                                         bounded_h5ad=args.mode == 'bounded')
        result['alignment_seconds'] = time.perf_counter() - started
        retained = np.flatnonzero(values.sum(axis=1)[np.arange(n) % template_rows] >= args.min_counts)
        result['min_counts'] = args.min_counts
        result['retained_cells'] = len(retained)
        assert output.obs_names.equals(pd.Index(expected_ids)[retained]) and output.var_names.equals(var.index)
        assert output.X.nnz == output.layers['counts'].nnz == len(retained) * k
        np.testing.assert_array_equal(output.obs['ref'].astype(str), obs['Library'].iloc[retained].astype(str))
        probe = np.array([0, len(retained) // 2, len(retained) - 1])
        for row in probe:
            source = values[retained[row] % template_rows]
            np.testing.assert_array_equal(output.layers['counts'][row].data, source)
            np.testing.assert_allclose(output.X[row].data, np.log1p(source / source.sum() * 1e4), rtol=1e-6)
        result['roster_counts_and_probe_checks'] = 'passed'
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2) + '\n')
        if args.markers:
            from altanalyze3.components.visualization.marker_heatmap_h5ad import generate_marker_heatmap_from_adata
            started = time.perf_counter()
            marker = generate_marker_heatmap_from_adata(output, cluster_key='ref',
                         out=str(root / 'markers.pdf'), top_n=50, marker_method='markerfinder',
                         seed=0, render_heatmap=False, write_heatmap_tsv=False,
                         write_expression_tsv=False, write_heatmap_cache=True)
            result['markerfinder_seconds'] = time.perf_counter() - started
            result['marker_artifacts_created'] = bool(marker.get('markers_tsv'))
        if args.umap:
            from altanalyze3.components.visualization.approximate_umap import approximate_umap
            reference_adata = ad.AnnData(X=sp.csr_matrix((2, 0)),
                 obs=pd.DataFrame({'state': ['A', 'B']}, index=['rA', 'rB']))
            reference_adata.obsm['X_umap'] = np.array([[-1., 0.], [1., 0.]])
            started = time.perf_counter()
            mapped = approximate_umap(query=output, reference=reference_adata,
                     query_cluster_key='ref', reference_cluster_key='state', umap_key='X_umap',
                     jitter=0.05, num_reference_cells=1, copy_query=False)
            assert mapped.query_adata is output
            assert output.obsm['X_umap'].shape == (len(retained), 2)
            assert np.isfinite(output.obsm['X_umap']).all()
            result['approximate_umap_seconds'] = time.perf_counter() - started
            started = time.perf_counter()
            output.write_h5ad(root / 'result.h5ad', compression='lzf')
            result['save_h5ad_seconds'] = time.perf_counter() - started
            result['scope'] += '; approximate UMAP and H5AD serialization'
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        result['peak_rss_gib'] = peak / (1024**3 if sys.platform == 'darwin' else 1024**2)
        result['compressed_input_gib'] = sum(path.stat().st_size for path, _ in entries) / 1024**3
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2) + '\n')
        print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
