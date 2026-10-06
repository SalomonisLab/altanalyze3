"""Build the user-authorized target-filtered AML bundle, retaining fitted weights."""
import gzip
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd

from altanalyze3.components.rna2metabolite.missingness_release import (
    RELEASE_PATH, filter_bundle, sha256, target_policy, validate_release,
)

ROOT = Path('/Users/saljh8/Dropbox/Collaborations/Grimes/Human-MS-impute')
COMPONENT = RELEASE_PATH.parent
OUT = COMPONENT / 'artifacts/rna2metabolite_aml_NA30_20261006_bundle.pkl.gz'
STATE = COMPONENT.parent / 'rna2lipid/integrity/decision_state.json'


def main():
    state = json.loads(STATE.read_text())
    approval = state.get('AML_target_missingness_approval', {})
    if approval.get('granted_by') != 'user' or approval.get('user_message') != 'If >30% values are NA, remove that feature':
        raise RuntimeError('Actual user missingness-rule approval required')
    questions = {q['id']: q for q in state['questions']}
    for decision in ('AML_target_encoding', 'AML_ten_all_missing_target_sources'):
        if questions[decision]['status'] != 'resolved' or not questions[decision].get('user_answer'):
            raise RuntimeError(f'Unresolved required AML decision: {decision}')
    source_trace = json.loads(Path(__file__).with_name('target_scale_trace.json').read_text())
    trace = source_trace['modalities']['metabolite']
    for script, digest in source_trace['scripts'].items():
        if sha256(script) != digest:
            raise RuntimeError(f'Verified original method script changed: {script}')
    for key in ('target_file', 'source_model', 'packaged_model'):
        if sha256(trace[key]['path']) != trace[key]['sha256']:
            raise RuntimeError(f'Verified baseline source changed: {key}')
    original_path = Path(trace['packaged_model']['path'])
    with gzip.open(original_path, 'rb') as handle:
        original = pickle.load(handle)
    table = pd.read_csv(trace['target_file']['path'], sep='\t', index_col=0)
    case_source = ROOT / 'data/matched_case_ids.json'
    cases = json.loads(case_source.read_text())['metab_rna']
    if len(cases) != 84 or len(original['X_columns']) != 12416 or len(original['Y_columns']) != 2533:
        raise RuntimeError('Original baseline dimensions changed')
    if len(set(original['X_columns'])) != 12416:
        raise RuntimeError('Duplicate original RNA identities')
    if original['metadata']['alpha'] != 100 or original['metadata']['nfeat'] != 1000:
        raise RuntimeError('Original ridge settings changed')
    policy = target_policy(table, original['Y_columns'], cases)
    updated = filter_bundle(original, policy)
    indices = np.flatnonzero(policy.retained.to_numpy())
    # Verify every retained learned value, plus every original RNA input scaler.
    for key in ('mu', 'sd'):
        np.testing.assert_array_equal(updated[key], original[key])
    for key in ('sel_idx', 'coef'):
        for new, i in zip(updated[key], indices):
            np.testing.assert_array_equal(new, original[key][i])
    np.testing.assert_array_equal(updated['intercept'], original['intercept'][indices])
    if updated['X_columns'] != original['X_columns']:
        raise RuntimeError('Original RNA roster changed')
    manifest = {'schema_version': 1, 'user_authorization': approval,
                'baseline_bundle': str(original_path), 'baseline_bundle_sha256': sha256(original_path),
                'target_source': trace['target_file'],
                'training_case_source': {'path': str(case_source), 'sha256': sha256(case_source)},
                'original_method_scripts': source_trace['scripts'],
                'training_cases': cases, 'RNA_genes': original['X_columns'],
                'original_targets': original['Y_columns'], 'retained_targets': updated['Y_columns'],
                'removed_targets': policy.index[~policy.retained].tolist(),
                'original_estimator': original['metadata']['estimator'],
                'retained_weights_identical': True, 'RNA_scalers_identical': True,
                'training_roster_unchanged': True, 'new_fitting_executed': False,
                'new_CV_executed': False, 'target_values_filled': False,
                'original_targets_count': len(policy), 'retained_targets_count': len(indices),
                'removed_targets_count': int((~policy.retained).sum()),
                'bundle_path': str(OUT), 'all_retained_models_fitted': True}
    # Validate all inputs before writing; publish via an immutable versioned path.
    OUT.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUT.with_suffix(OUT.suffix + '.tmp')
    with temporary.open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as handle:
            pickle.dump(updated, handle, protocol=4)
    manifest['bundle_sha256'] = sha256(temporary)
    if OUT.exists() and sha256(OUT) != manifest['bundle_sha256']:
        temporary.unlink()
        raise RuntimeError('Versioned release already exists with different bytes')
    temporary.replace(OUT)
    RELEASE_PATH.write_text(json.dumps(manifest, indent=2) + '\n')
    policy.to_csv(Path(__file__).with_name('target_missingness_policy.tsv'), sep='\t')
    validate_release(OUT, updated)
    if sha256(original_path) != trace['packaged_model']['sha256']:
        raise RuntimeError('Original rollback bundle changed')
    print(json.dumps({key: manifest[key] for key in ['original_targets_count', 'retained_targets_count',
                    'removed_targets_count', 'retained_weights_identical', 'new_fitting_executed', 'bundle_sha256']}))


if __name__ == '__main__':
    main()
