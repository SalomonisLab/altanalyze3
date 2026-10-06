"""Verify completed standard web outputs without loading their full RNA matrices."""
import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
from anndata.io import read_elem, sparse_dataset
import h5py
import numpy as np
import pandas as pd

from altanalyze3.components.rna2lipid.api import load_bundle
from altanalyze3.components.cellHarmony.markerFinder import detect_input_scaling


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--job-dir', type=Path, required=True)
    parser.add_argument('--model-identity', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    meta = json.loads((args.job_dir / 'job.json').read_text())
    assert meta['status'] == 'completed', 'A partial/failed job cannot satisfy validation'
    assert meta.get('bundle', {}).get('status') == 'completed', 'Large-job serving bundle must succeed'
    assert set(meta['bundle']['sources']) >= {'rna', 'lipids'}
    assert meta.get('fastcomm_analysis', {}).get('enabled'), 'Communication analysis must succeed'
    identity = json.loads(args.model_identity.read_text())
    bundle = load_bundle()
    with bundle.bundle_path.open('rb') as model_file:
        assert hashlib.file_digest(model_file, 'sha256').hexdigest() == identity['model_sha256']
    assert list(bundle.input_genes) == identity['inputs']
    assert list(bundle.output_lipids) == identity['outputs']
    outputs = args.job_dir / 'outputs'
    with h5py.File(args.source, 'r') as source, \
            h5py.File(outputs / 'combined_with_umap_and_markers.h5ad', 'r') as rna, \
            h5py.File(outputs / 'qc_passed_unaligned_cells.h5ad', 'r') as excluded, \
            h5py.File(outputs / 'combined_with_umap_and_markers_lipids.h5ad', 'r') as lipid:
        original_obs, original_var = read_elem(source['obs']), read_elem(source['var'])
        rna_obs, rna_var = read_elem(rna['obs']), read_elem(rna['var'])
        excluded_obs, excluded_var = read_elem(excluded['obs']), read_elem(excluded['var'])
        lipid_obs, lipid_var = read_elem(lipid['obs']), read_elem(lipid['var'])
        assert original_obs.index.is_unique and rna_obs.index.is_unique and excluded_obs.index.is_unique
        assert rna_var.index.equals(original_var.index) and excluded_var.index.equals(original_var.index)
        assert not rna_obs.index.intersection(excluded_obs.index).size
        assert rna_obs.index.union(excluded_obs.index).sort_values().equals(original_obs.index.sort_values())
        for field in ('Library', 'sample'):
            for actual in (rna_obs, excluded_obs):
                pd.testing.assert_series_equal(actual[field].astype(str),
                                               original_obs.loc[actual.index, field].astype(str))
        assert lipid_obs.index.equals(rna_obs.index)
        communication = meta['fastcomm_analysis']
        assert communication['summary']['n_cells'] == len(rna_obs)
        assert communication['per_sample']['n_splits'] == rna_obs[communication['sample_key']].astype(str).nunique()
        assert list(lipid_var.index) == identity['outputs']
        assert lipid['X'].shape == (len(rna_obs), len(identity['outputs']))
        for start in range(0, len(rna_obs), 8192):
            values = lipid['X'][start:start + 8192]
            assert np.isfinite(values).all() and (values >= 0).all()
        assert np.isfinite(read_elem(rna['obsm/X_umap'])).all()
        rows = np.sort(np.random.default_rng(61).choice(len(rna_obs), min(128, len(rna_obs)), replace=False))
        expression = sparse_dataset(rna['X'])[rows]
        counts = sparse_dataset(rna['layers/counts'])[rows]
        raw = sparse_dataset(rna['layers/soupx_raw'])[rows]
        source_rows = original_obs.index.get_indexer(rna_obs.index[rows])
        original_counts = sparse_dataset(source['layers/counts'])[source_rows].astype(np.float32)
        assert (raw != original_counts).nnz == 0
        expected = counts.copy()
        totals = np.asarray(expected.sum(axis=1)).ravel()
        divisors = totals / 10000
        divisors[divisors == 0] = 1
        np.divide(expected.data, np.repeat(divisors, np.diff(expected.indptr)), out=expected.data)
        np.log1p(expected.data, out=expected.data)
        np.testing.assert_array_equal(expression.indptr, expected.indptr)
        np.testing.assert_array_equal(expression.indices, expected.indices)
        np.testing.assert_allclose(expression.data, expected.data, rtol=1e-6, atol=1e-6)
        scaling = detect_input_scaling(expression)
        assert scaling['status'] == 'ok', scaling
        probe = ad.AnnData(expression, obs=rna_obs.iloc[rows].copy(), var=rna_var.copy())
        optimized = bundle.predict_from_adata(probe).predictions
        bundle._lipidwise_linear_parameters = lambda: None
        baseline = bundle.predict_from_adata(probe).predictions
        np.testing.assert_allclose(optimized, baseline, rtol=1e-12, atol=1e-12)
        expected_lipids = np.maximum(baseline.to_numpy(dtype=np.float32), 0)
        np.testing.assert_array_equal(lipid['X'][rows], expected_lipids)
        report = dict(job_id=meta['job_id'], status='verified', input_cells=len(original_obs),
                      aligned_cells=len(rna_obs), below_cutoff_cells=len(excluded_obs),
                      full_gene_panel=len(original_var), source_libraries=original_obs['Library'].nunique(),
                      lipid_outputs=len(lipid_var), all_lipid_values_finite_nonnegative=True,
                      source_cell_roster_preserved=True, gene_order_preserved=True,
                      model_identity_unchanged=True, raw_counts_probe_exact=True,
                      normalized_rna_probe_matches_corrected_counts=True,
                      lipid_probe_matches_original_estimators=True, probe_cells=len(rows),
                      communication_and_serving_bundle_completed=True,
                      max_lipid_prediction_difference=float(np.max(np.abs(optimized.to_numpy()-baseline.to_numpy()))),
                      normalization=scaling)
    args.report.write_text(json.dumps(report, indent=2, default=int) + '\n')
    print(json.dumps(report, indent=2, default=int))


if __name__ == '__main__':
    main()
