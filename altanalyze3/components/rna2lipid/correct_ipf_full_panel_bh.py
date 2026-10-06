"""Apply the explicitly approved BH correction to saved full-panel IPF statistics."""
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

from candidate_integrity import CONTRACT, STATE, checked_file, read_json, require_review, sha256, validate_run_manifest
from full_panel_bh import KEYS, correct_full_panel_bh
import evaluate_full_panel_ipf as evaluation

HERE = Path(__file__).resolve().parent
SOURCE = HERE / 'artifacts/LungMAP_full202_native_log2_20261005'
OUT = HERE / 'artifacts/LungMAP_full202_BH_corrected_20261006'
PHASE = 'full202_lipid_BH_correction'
MODELS = [evaluation.PRIOR, evaluation.LABEL]


def read_csv(path, **kwargs):
    return pd.read_csv(path, float_precision='round_trip', **kwargs)


def run():
    # The original source/method gate remains binding; this is not a bypass or new model.
    validate_run_manifest(SOURCE / 'run_manifest.json')
    manifest_path = OUT / 'BH_correction_manifest.json'
    manifest = read_json(manifest_path)
    state = read_json(STATE)
    if 'BH_correction_approval' not in state:
        raise RuntimeError('No recorded user authorization for the BH correction.')
    reviewed = dict(state, approval=state['BH_correction_approval'])
    require_review(reviewed, sha256(manifest_path), PHASE)
    if manifest['analysis_phase'] != PHASE or manifest['BH_policy'] != 'full_panel_no_abundance_filter':
        raise RuntimeError('Unexpected approved correction policy.')
    for name, record in manifest['sources'].items():
        checked_file(record, name)
    c = read_json(CONTRACT)
    if manifest['required_lipids'] != c['lipids'] or manifest['models'] != MODELS:
        raise RuntimeError('Approved correction identity contract differs from baseline.')
    before = read_csv(SOURCE / 'all_cellHarmony_tests.csv')
    original_joined = read_csv(SOURCE / 'all_measured_vs_imputed_tests.csv')
    protected = {str(p): sha256(p) for p in [SOURCE / 'candidate_bundle.pkl', Path(c['production_bundle']['path'])]}
    corrected = correct_full_panel_bh(before, c['lipids'], MODELS)
    lookup = corrected.set_index(KEYS + ['gene'])
    joined = original_joined.copy()
    joined['fdr_legacy_filtered'] = joined.fdr
    indices = pd.MultiIndex.from_frame(joined[KEYS + ['gene']])
    if not indices.isin(lookup.index).all():
        raise RuntimeError('Experimental comparison has unresolved test identities.')
    changes = ['fdr', 'significant_BH_005', 'significant_BH_010', 'BH_panel_size', 'BH_policy']
    for field in changes:
        joined[field] = lookup[field].reindex(indices).to_numpy()
    untouched = [f for f in original_joined if f not in changes]
    pd.testing.assert_frame_equal(original_joined[untouched], joined[untouched])
    corrected.to_csv(OUT / 'all_cellHarmony_tests.csv', index=False)
    joined.to_csv(OUT / 'all_measured_vs_imputed_tests.csv', index=False)
    # Check the serialized corrected results too, including every unchanged original field.
    serialized = read_csv(OUT / 'all_cellHarmony_tests.csv')
    unchanged = [f for f in before if f not in changes]
    pd.testing.assert_frame_equal(before[unchanged], serialized[unchanged])
    comparison_keys = ['contrast', 'unit', 'model']
    for key, group in corrected.groupby(comparison_keys, sort=False):
        contrast, unit, label = key
        folder = OUT / 'differentials' / contrast / unit / label
        folder.mkdir(parents=True, exist_ok=True)
        old_folder = SOURCE / 'differentials' / contrast / unit / label
        group.to_csv(folder / 'all_tested_lipids.csv', index=False)
        group[group.pval.lt(.05)].to_csv(folder / 'raw_p_005_differentials.csv', index=False)
        for cutoff, suffix in [(.05, '005'), (.1, '010')]:
            group[group.fdr.lt(cutoff)].to_csv(folder / f'BH_{suffix}_differentials.csv', index=False)
        counts = group.groupby('population', sort=False).agg(
            n_case=('n_case', 'first'), n_control=('n_control', 'first'), tested_lipids=('gene', 'size'),
            raw_p_005_calls=('significant_raw_p_005', 'sum'), BH_005_calls=('significant_BH_005', 'sum'),
            BH_010_calls=('significant_BH_010', 'sum')).reset_index()
        counts.to_csv(folder / 'differential_counts.csv', index=False)
        parameters = read_json(old_folder / 'parameters.json')
        parameters.update(BH_scope='All 202 required lipid outputs per population/contrast',
                          BH_policy=manifest['BH_policy'], abundance_filter=False,
                          raw_statistics_source=str(old_folder / 'all_tested_lipids.csv'),
                          approval_manifest_sha256=sha256(manifest_path))
        evaluation.dump(folder / 'parameters.json', parameters)
        shutil.copyfile(old_folder / 'population_replication_census.csv', folder / 'population_replication_census.csv')
    for filename in ['population_replication_census.csv', 'lipid_matching_audit.csv',
                     'unmatched_experimental_identity_audit.csv', 'experimental_lung_IPF_lipids_all_544.csv']:
        shutil.copyfile(SOURCE / filename, OUT / filename)
    maps = read_csv(SOURCE / 'lipid_matching_audit.csv')
    measured = read_csv(SOURCE / 'experimental_lung_IPF_lipids_all_544.csv', index_col=0)
    registry = read_csv(SOURCE / 'provenance/IPF_vs_healthy_registry.tsv', sep='\t')
    audit = {'user_authorized': True, 'user_message': state['BH_correction_approval']['user_message'],
             'BH_policy': manifest['BH_policy'], 'panel_size': len(c['lipids']),
             'rows_before': len(before), 'rows_after': len(corrected),
             'raw_p_values_unchanged': True, 'folds_and_group_sizes_unchanged': True,
             'experimental_statistics_unchanged': True, 'raw_test_functions_not_rerun': True,
             'retraining_or_reimputation_performed': False, 'source_artifacts_modified': False,
             'complete_model_comparisons': corrected.groupby(KEYS).ngroups,
             'BH_values_before': int(before.fdr.notna().sum()), 'BH_values_after': int(corrected.fdr.notna().sum()),
             'production_and_candidate_bundle_sha256': protected}
    evaluation.dump(OUT / 'BH_correction_audit.json', audit)
    evaluation.OUT = OUT
    evaluation.summarize(corrected, joined, maps, measured, registry,
                         bh_policy=manifest['BH_policy'], reference_out=SOURCE)
    old_summary = read_csv(SOURCE / 'cohort_concordant_union_summary.csv')
    new_summary = read_csv(OUT / 'cohort_concordant_union_summary.csv')
    pd.testing.assert_frame_equal(old_summary[old_summary.overlap_gate.eq('both_rawp005')].reset_index(drop=True),
                                  new_summary[new_summary.overlap_gate.eq('both_rawp005')].reset_index(drop=True))
    if protected != {p: sha256(p) for p in protected}:
        raise RuntimeError('Protected model changed during correction.')
    audit.update(raw_p_concordance_counts_unchanged=True, serialized_raw_statistics_unchanged=True,
                 protected_bundles_unchanged=True)
    evaluation.dump(OUT / 'BH_correction_audit.json', audit)
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    run()
