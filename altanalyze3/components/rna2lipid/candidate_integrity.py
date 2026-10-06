"""Fail-closed guards for the user-paused LungMAP corrected-model comparison.

No scientific operation is performed by the status command. Local guards do
not authenticate a conversation or prevent someone from editing Python code.
Approval records must be grounded in an actual user response.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONTRACT = HERE / 'integrity/baseline_contract.json'
STATE = HERE / 'integrity/decision_state.json'
BASELINE_CONTRACT_SHA256 = 'b2620153c969a566329d241003f54b3631e656d9acbf9677743d22bb47a0529c'


class IntegrityError(RuntimeError):
    pass


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    try:
        value = json.loads(Path(path).read_text())
    except (OSError, ValueError) as exc:
        raise IntegrityError(f'Missing or unreadable integrity document: {path}. Ask the user; do not substitute.') from exc
    if not isinstance(value, dict):
        raise IntegrityError(f'Integrity document must be an object: {path}')
    return value


def exact_ids(expected, actual, label):
    expected, actual = list(expected), list(actual)
    duplicates = [key for key, count in Counter(actual).items() if count > 1]
    missing = [key for key in expected if key not in set(actual)]
    extra = [key for key in actual if key not in set(expected)]
    if duplicates or missing or extra:
        raise IntegrityError(f'{label} mismatch: missing={missing}, extra={extra}, duplicates={duplicates}. '
                             'STOP and ask for the source/mapping; do not intersect or discard records.')
    if actual != expected:
        raise IntegrityError(f'{label} order differs from the reviewed baseline; resolve the explicit mapping before execution.')


def require_method(expected, proposed):
    if proposed != expected:
        keys = sorted(set(expected) | set(proposed)) if isinstance(proposed, dict) else ['method']
        changed = [k for k in keys if not isinstance(proposed, dict) or expected.get(k) != proposed.get(k)]
        raise IntegrityError('Baseline method differs: ' + ', '.join(changed) +
                             '. STOP and discuss the concrete change with the user before proceeding.')


def checked_file(record, label):
    if not isinstance(record, dict) or not record.get('path') or not record.get('sha256'):
        raise IntegrityError(f'{label}: source path and SHA-256 are required; ask the user for the source.')
    path = Path(record['path'])
    if not path.is_file():
        raise IntegrityError(f'{label}: source is unavailable at {path}. Ask the user where it is; do not guess.')
    if sha256(path) != record['sha256']:
        raise IntegrityError(f'{label}: source hash changed after review: {path}. Reconcile and review before use.')
    return path


def validate_contract(contract):
    if contract.get('schema_version') != 1:
        raise IntegrityError('Unsupported or missing baseline contract schema.')
    for label in ('samples', 'lipids', 'genes'):
        ids = contract.get(label)
        if not isinstance(ids, list) or not ids or len(ids) != len(set(ids)):
            raise IntegrityError(f'Invalid baseline {label}; do not derive it from a candidate subset.')
    checked_file(contract.get('production_bundle'), 'Production bundle')
    for label, record in contract.get('pinned_code', {}).items():
        checked_file(record, 'Pinned method code: ' + label)
    if set(contract.get('pinned_code', {})) != {'trainer', 'api', 'delivered_trainer'}:
        raise IntegrityError('The baseline must pin trainer, API and delivered trainer code.')


def require_review(state, manifest_digest):
    if state.get('schema_version') != 1:
        raise IntegrityError('Missing or unsupported decision-state schema.')
    pending = [x for x in state.get('questions', []) if x.get('status') != 'resolved']
    if state.get('analysis_status') != 'approved' or pending:
        question = pending[0].get('question', '') if pending else 'The user has not approved resuming this analysis.'
        raise IntegrityError('ANALYSIS BLOCKED. ' + question +
                             ' Wait for the user. Safeguard tests, silence and generic proceed instructions do not resolve it.')
    approval = state.get('approval')
    if not isinstance(approval, dict) or approval.get('manifest_sha256') != manifest_digest:
        raise IntegrityError('No user approval for this exact run manifest; do not reuse approval for a different run.')
    for field in ('user_message', 'conversation_reference'):
        if not isinstance(approval.get(field), str) or not approval[field].strip():
            raise IntegrityError('Approval requires the actual user message and conversation reference; never fabricate them.')
    if approval.get('granted_by') != 'user':
        raise IntegrityError('Only the user can authorize this paused analysis.')
    for item in state.get('questions', []):
        if not item.get('user_answer') or not item.get('resolution_evidence'):
            raise IntegrityError('A resolved question lacks its user answer or supporting evidence.')


def validate_provenance(provenance, contract):
    if provenance.get('schema_version') != 1 or provenance.get('data_role') != 'corrected_training_targets':
        raise IntegrityError('Corrected-target provenance is required; legacy reconstruction is not corrected evidence.')
    if provenance.get('unresolved') != []:
        raise IntegrityError('Corrected-target provenance contains unresolved dependencies. Ask the user and wait.')
    for kind, expected in [('sample_sources', contract['samples']), ('feature_sources', contract['lipids'])]:
        records = provenance.get(kind)
        if not isinstance(records, list):
            raise IntegrityError(f'Missing {kind}; every required identity needs provenance.')
        exact_ids(expected, [r.get('id') for r in records], kind)
        for record in records:
            if record.get('status') != 'verified' or not record.get('evidence') or not record.get('source_identifier'):
                raise IntegrityError(f'{kind}/{record.get("id")}: unresolved source or mapping. Ask the user before reconstruction.')
            checked_file(record.get('source'), f'{kind}/{record["id"]}')
    if not provenance.get('transformation_description'):
        raise IntegrityError('The corrected target transformation must be documented and reviewed.')
    if 'imputed_entries' not in provenance or not isinstance(provenance['imputed_entries'], list):
        raise IntegrityError('Explicitly enumerate target imputation, even when the list is empty.')
    if provenance['imputed_entries'] and not provenance.get('imputation_user_decision'):
        raise IntegrityError('Target imputation requires its specifically discussed user decision.')


def validate_tables(x, y, contract):
    exact_ids(contract['samples'], x.index, 'RNA training samples')
    exact_ids(contract['samples'], y.index, 'Lipid training samples')
    exact_ids(contract['genes'], x.columns, 'RNA input features')
    exact_ids(contract['lipids'], y.columns, 'Lipid output features')
    for name, frame in [('RNA', x), ('lipid targets', y)]:
        for row in frame.itertuples(index=False, name=None):
            try:
                valid = all(math.isfinite(float(v)) for v in row)
            except (TypeError, ValueError):
                valid = False
            if not valid:
                raise IntegrityError(f'Nonfinite or nonnumeric {name}. Ask about the missing source; do not fill or drop it automatically.')


def validate_run_manifest(manifest_path):
    """Validate authorization, exact intended method, files and provenance.

    Expected contract and decision-state paths are fixed. No force/skip flag.
    """
    contract, state = read_json(CONTRACT), read_json(STATE)
    if sha256(CONTRACT) != BASELINE_CONTRACT_SHA256:
        raise IntegrityError('Baseline contract changed; discuss and review it instead of narrowing it to fit a candidate.')
    manifest = read_json(manifest_path)
    require_review(state, sha256(manifest_path))
    validate_contract(contract)
    if manifest.get('schema_version') != 1 or manifest.get('baseline_contract_sha256') != sha256(CONTRACT):
        raise IntegrityError('Run manifest does not reference the reviewed baseline contract.')
    require_method(contract['method'], manifest.get('method'))
    inputs = manifest.get('inputs', {})
    if set(inputs) != {'RNA', 'corrected_targets', 'provenance'}:
        raise IntegrityError('Run inputs must include RNA, corrected targets and complete provenance.')
    if inputs['RNA'] != contract['RNA_source']:
        raise IntegrityError('Training RNA source differs from the pinned production source.')
    paths = {key: checked_file(value, key) for key, value in inputs.items()}
    if inputs['corrected_targets'].get('data_role') != 'corrected_training_targets':
        raise IntegrityError('Legacy targets cannot stand in for corrected targets.')
    provenance = read_json(paths['provenance'])
    if provenance.get('target_sha256') != inputs['corrected_targets']['sha256']:
        raise IntegrityError('Provenance is not bound to this target table.')
    validate_provenance(provenance, contract)
    return contract, paths


def fit_verified_candidate(manifest_path):
    """The candidate fitting gateway: all checks precede estimator fitting.

    Returns the normal trainer result without writing/deploying a bundle.
    Only an explicitly reviewed future full-panel run can reach this call.
    """
    contract, paths = validate_run_manifest(manifest_path)
    import pandas as pd
    import sklearn
    if sklearn.__version__ != contract['method']['sklearn_version']:
        raise IntegrityError('scikit-learn version differs from the deployed baseline.')
    x = pd.read_csv(paths['RNA'], index_col=0).T
    # These are the production source preprocessing operations, not a
    # candidate intersection. Missing required keys raise before fitting.
    x = x.T.groupby(level=0).mean().T
    x.index = x.index.str.strip()
    x.columns = x.columns.str.strip()
    try:
        x = x.loc[contract['samples'], contract['genes']].apply(pd.to_numeric, errors='raise')
    except (KeyError, ValueError) as exc:
        raise IntegrityError('Pinned RNA source no longer resolves every baseline identity. Ask the user.') from exc
    # Explicit baseline feature-median policy only; target filling is not done.
    x = x.fillna(x.median())
    y = pd.read_csv(paths['corrected_targets'], index_col=0)
    validate_tables(x, y, contract)
    try:
        from .training import fit_sparse_lipidwise_elasticnet
    except ImportError:
        import sys
        sys.path.insert(0, str(HERE.parents[2]))
        from altanalyze3.components.rna2lipid.training import fit_sparse_lipidwise_elasticnet
    import inspect
    if sha256(inspect.getsourcefile(fit_sparse_lipidwise_elasticnet)) != contract['pinned_code']['trainer']['sha256']:
        raise IntegrityError('Loaded trainer is not the pinned baseline implementation.')
    # Recheck source hashes immediately before fitting.
    validate_run_manifest(manifest_path)
    return fit_sparse_lipidwise_elasticnet(x, y, x, **contract['method']['fit_settings'])


def reject_retired_workflow(name):
    """Old narrowed/alternative candidate drivers are not a permitted fallback."""
    raise IntegrityError(f'RETIRED WORKFLOW BLOCKED: {name}. This driver does not satisfy the requested '
        '202-lipid / 50-profile / 1303-gene production ElasticNetCV comparison. '
        'Do not bypass this guard or drop records. Resolve pending questions with the user; '
        'a future reviewed run must use fit_verified_candidate and its complete provenance manifest.')


def require_diagnostic_authorization(stage):
    """A pending source question also stops further reconstruction diagnostics."""
    state = read_json(STATE)
    for authorization in state.get('diagnostic_authorizations', []):
        if (authorization.get('stage') == stage and authorization.get('granted_by') == 'user'
                and authorization.get('user_message') and authorization.get('conversation_reference')):
            return
    raise IntegrityError(f'DIAGNOSTIC BLOCKED: {stage}. A source question is pending. '
        'Ask the user and wait before further searching or reconstruction; '
        'safeguard implementation is not authorization to resume scientific diagnostics.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--status', action='store_true', required=True)
    parser.parse_args()
    try:
        if sha256(CONTRACT) != BASELINE_CONTRACT_SHA256:
            raise IntegrityError('Baseline contract changed; explicit review is required.')
        validate_contract(read_json(CONTRACT))
        require_review(read_json(STATE), '')
    except IntegrityError as exc:
        print(json.dumps({'execution_allowed': False, 'reason': str(exc)}, indent=2))
        return 2
    print(json.dumps({'execution_allowed': False, 'reason': 'A specific reviewed run manifest is still required.'}))
    return 2


if __name__ == '__main__':
    raise SystemExit(main())
