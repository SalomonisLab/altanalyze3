"""Train the explicitly authorized CPTT counterpart through the integrity gate.

This does not activate defaults, replace CPM weights or publish a model.
Run with the baseline scikit-learn 1.6.1 environment, from the package parent.
"""
from pathlib import Path
import json
import pickle

import numpy as np
import pandas as pd

from .candidate_integrity import (
    STATE, read_json, sha256, validate_run_manifest, validate_tables,
    training_RNA_for_manifest, fit_verified_candidate,
)

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / 'artifacts/lung_CPTT_model_20261008'
MANIFEST = OUTPUT / 'run_manifest.json'
BUNDLE = HERE / 'rna2lipid_hs_lung_CPTT_native_log2_20261008_bundle.pkl'


def main():
    contract, paths = validate_run_manifest(MANIFEST)
    if BUNDLE.exists():
        raise RuntimeError('CPTT release artifact already exists; do not overwrite immutable weights.')
    x = pd.read_csv(paths['RNA'], index_col=0).T
    x = x.T.groupby(level=0).mean().T
    x.index = x.index.str.strip()
    x.columns = x.columns.str.strip()
    x = x.loc[contract['samples'], contract['genes']].apply(pd.to_numeric, errors='raise')
    x = x.fillna(x.median())
    cptt = training_RNA_for_manifest(x, read_json(MANIFEST), read_json(STATE))
    y = pd.read_csv(paths['corrected_targets'], index_col=0)
    validate_tables(cptt, y, contract)
    # Verify the depth conversion by an independent inverse before fitting.
    recovered = np.log2(1 + 100 * np.expm1(cptt.to_numpy() * np.log(2)))
    roundtrip_error = float(np.max(abs(recovered - x.to_numpy())))
    if roundtrip_error > 1e-12:
        raise RuntimeError('CPM/CPTT conversion failed its round-trip check.')
    predictions, bundle, seconds = fit_verified_candidate(MANIFEST)
    if bundle['X_columns'] != contract['genes'] or bundle['Y_columns'] != contract['lipids']:
        raise RuntimeError('Fitted model differs from the complete instructed panel.')
    if not np.isfinite(predictions.to_numpy()).all():
        raise RuntimeError('Nonfinite fitted predictions; do not fill, filter or publish.')
    bundle['training_samples'] = contract['samples']
    bundle['final_training_seconds'] = seconds
    bundle['metadata'] = {
        'RNA_input_scale': 'log2(CP10k + 1)',
        'RNA_normalization_target_sum': 10000,
        'intended_modality': 'single_cell',
        'target_scale': 'native_relative_log2',
        'linear_inverse': '2**prediction; no subtraction or added pseudocount',
        'target_scaling': {'mode': 'standard'},
        'candidate_only': True,
        'training_manifest_sha256': sha256(MANIFEST),
        'RNA_scale_change': read_json(MANIFEST)['approved_RNA_scale_change'],
        'training_algorithm': 'Identical pinned lipid-wise ElasticNetCV; only authorized RNA representation changes',
    }
    # Exact target scalers are an independent check that the lipid space did not change.
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        bulk = pickle.loads((HERE / 'rna2lipid_hs_lung_native_log2_20261005_bundle.pkl').read_bytes())
    for field in ('mean_', 'scale_', 'var_'):
        np.testing.assert_array_equal(getattr(bundle['scaler_y'], field), getattr(bulk['scaler_y'], field))
    cptt.to_csv(OUTPUT / 'training_RNA_log2_CP10k1.csv')
    predictions.to_csv(OUTPUT / 'training_predictions_native_log2.csv')
    bundle['summary'].to_csv(OUTPUT / 'training_fit_summary.csv', index=False)
    bundle['candidate_models'].to_csv(OUTPUT / 'all_808_candidate_fits.csv', index=False)
    bundle['coefficients'].to_csv(OUTPUT / 'selected_coefficients.csv', index=False)
    with BUNDLE.open('xb') as handle:
        pickle.dump(bundle, handle, protocol=4)
    audit = {
        'training_profiles': len(contract['samples']), 'RNA_genes': len(contract['genes']),
        'lipid_targets': len(contract['lipids']), 'estimators': len(bundle['models']),
        'RNA_source_scale': 'log2(CPM + 1)', 'RNA_model_scale': 'log2(CP10k + 1)',
        'conversion_roundtrip_max_abs_error': roundtrip_error,
        'lipid_target_scalers_byte_equivalent_arrays': True,
        'RNA_source_sha256': sha256(paths['RNA']),
        'target_source_sha256': sha256(paths['corrected_targets']),
        'trainer_sha256': sha256(HERE / 'training.py'),
        'training_manifest_sha256': sha256(MANIFEST),
        'CPTT_bundle_sha256': sha256(BUNDLE),
        'CPM_bundle_sha256_unchanged': sha256(HERE / 'rna2lipid_hs_lung_native_log2_20261005_bundle.pkl'),
        'training_seconds': seconds,
        'validation_type': 'Training/inference/scale integrity; not external disease validation',
    }
    (OUTPUT / 'training_integrity.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == '__main__':
    main()
