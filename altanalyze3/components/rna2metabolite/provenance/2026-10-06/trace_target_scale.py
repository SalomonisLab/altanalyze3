"""Read-only trace of packaged AML models to saved targets and source workbook.

No fitting, missing-value filling, feature/sample exclusion or model writes.
Reuses the original target extraction functions to verify their saved outputs.
"""
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import pickle

import numpy as np
import openpyxl
import pandas as pd

ROOT = Path('/Users/saljh8/Dropbox/Collaborations/Grimes/Human-MS-impute')
A3 = Path(__file__).resolve().parents[3]
WORKBOOK = Path('/Users/saljh8/Downloads/43018_2026_1175_MOESM3_ESM.xlsx')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    path = ROOT / 'code/build_unique_ms_tables.py'
    spec = importlib.util.spec_from_file_location('original_ms_extraction', path)
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)  # main() is not executed.
    workbook = openpyxl.load_workbook(WORKBOOK, read_only=True, data_only=True)
    matched = json.loads((ROOT / 'data/matched_case_ids.json').read_text())
    report = {'source_root': str(ROOT), 'model_or_dataset_changed': False,
              'fit_executed': False, 'missing_values_filled': False, 'targets_excluded': False,
              'workbook': {'path': str(WORKBOOK), 'sha256': sha(WORKBOOK)}, 'modalities': {}}
    for name, sheets, column in [
        ('lipid', [('Table26', 'pos'), ('Table27', 'neg')], 'Lipid'),
        ('metabolite', [('Table21', 'HILIC'), ('Table22', 'RP')], 'Metabolite'),
    ]:
        records = []
        def lipid_transform(raw):
            logged = np.log2(raw.clip(lower=0) + 1)
            return logged.sub(logged.median(axis=0), axis=1)
        for sheet, assay in sheets:
            records.extend(source.candidates(source.sheet_df(workbook, sheet), column, assay,
                lipid_transform if name == 'lipid' else lambda raw: raw,
                detect_on_raw=(name == 'lipid')))
        regenerated, annotation = source.assemble(source.pick_best(records))
        basename = 'lipidomics' if name == 'lipid' else 'metabolomics'
        target_path = ROOT / f'data/{basename}_unique_matched.txt.gz'
        saved = pd.read_csv(target_path, sep='\t', index_col=0)
        expected_cases = matched['lipid_rna' if name == 'lipid' else 'metab_rna']
        assert list(saved.columns) == expected_cases, 'Training sample roster/order mismatch'
        assert list(saved.index) == list(regenerated.index), 'Target roster/order mismatch'
        # .loc requires every saved sample; there is no permissive intersection.
        reconstructed = regenerated.loc[saved.index, saved.columns]
        np.testing.assert_allclose(reconstructed, saved, rtol=0, atol=1e-12, equal_nan=True)
        saved_annotation = pd.read_csv(ROOT / f'data/{basename}_unique_annotation.tsv', sep='\t', index_col=0)
        assert list(saved_annotation.index) == list(saved.index)
        assert list(saved_annotation.protocol) == list(annotation.loc[saved.index].protocol)
        source_model_path = ROOT / f'models/rna2{name}_aml.pkl'
        component = A3 / ('rna2lipid/aml' if name == 'lipid' else 'rna2metabolite')
        bundle_path = component / f'artifacts/rna2{name}_aml_bundle.pkl.gz'
        with source_model_path.open('rb') as handle:
            original = pickle.load(handle)
        with gzip.open(bundle_path, 'rb') as handle:
            bundle = pickle.load(handle)
        assert original['universe_genes'] == bundle['X_columns'], 'RNA input roster mismatch'
        assert original['targets'] == bundle['Y_columns'] == list(saved.index), 'Model target roster mismatch'
        assert original['n_cases'] == saved.shape[1], 'Training sample count mismatch'
        for key in ('mu', 'sd', 'intercept'):
            np.testing.assert_array_equal(np.asarray(original[key], np.float32), bundle[key])
        for original_values, packaged_values in zip(original['coef'], bundle['coef']):
            np.testing.assert_array_equal(np.asarray(original_values, np.float32), packaged_values)
        assert len(original['coef']) == len(bundle['coef']) == len(saved.index)
        assert len(original['sel_idx']) == len(bundle['sel_idx']) == len(saved.index)
        for original_values, packaged_values in zip(original['sel_idx'], bundle['sel_idx']):
            np.testing.assert_array_equal(original_values, packaged_values)
        difference = np.abs(reconstructed.to_numpy() - saved.to_numpy())
        valid = np.array([len(c) > 0 and np.isfinite(i) for c, i in zip(bundle['coef'], bundle['intercept'])])
        result = {
            'targets': len(saved), 'training_cases': saved.shape[1], 'RNA_genes': len(bundle['X_columns']),
            'target_identities_and_sample_roster_verified': True, 'source_model_matches_packaged_arrays': True,
            'source_workbook_matches_all_saved_targets': True,
            'maximum_target_reconstruction_error': float(np.nanmax(difference)),
            'target_file': {'path': str(target_path), 'sha256': sha(target_path)},
            'source_model': {'path': str(source_model_path), 'sha256': sha(source_model_path)},
            'packaged_model': {'path': str(bundle_path), 'sha256': sha(bundle_path)},
            'fitted_targets': int(valid.sum()), 'unfitted_targets': int((~valid).sum()),
            'unfitted_target_ids': [t for t, fitted in zip(bundle['Y_columns'], valid) if not fitted],
            'saved_target_transform': ('log2(max(raw_intensity, 0) + 1) - sample_median_within_ion_mode'
                                       if name == 'lipid' else 'identity import of workbook Tables21/22'),
            'upstream_metabolite_scale_evidence': ('Original extraction script declares log2/median-centered; '
                'workbook values imported verbatim. This audit does not establish an upstream absolute-intensity inverse.'
                if name == 'metabolite' else None),
            'training_target_transform_after_loading': 'none; ridge fitted directly to saved Y',
            'normalized_linear_representation': '2**saved_target',
            'subtract_one_after_exponentiating_centered_targets': False,
        }
        report['modalities'][name] = result
        print(name, {k: result[k] for k in ['targets', 'training_cases', 'RNA_genes',
              'maximum_target_reconstruction_error', 'fitted_targets', 'unfitted_targets']}, flush=True)
    workbook.close()
    report['scripts'] = {str(ROOT / f'code/{name}'): sha(ROOT / f'code/{name}') for name in [
        'extract_ms_matrices.py', 'build_unique_ms_tables.py', 'evaluate_imputation.py',
        'build_imputation_model.py', 'make_altanalyze3_bundles.py']}
    out = Path(__file__).with_name('target_scale_trace.json')
    out.write_text(json.dumps(report, indent=2) + '\n')
    print(out, flush=True)


if __name__ == '__main__':
    main()
