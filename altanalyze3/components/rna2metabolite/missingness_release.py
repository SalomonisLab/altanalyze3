"""Apply the explicitly authorized >30% training-target missingness rule.

The estimator and all retained weights are unchanged. Missingness is measured
before filling, across the full original matched training-case roster.
"""
from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd

RELEASE_PATH = Path(__file__).with_name('missingness_release.json')


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def target_policy(table, target_ids, case_ids):
    if not table.index.is_unique or not table.columns.is_unique:
        raise ValueError('Duplicate target or training-case identity')
    if list(table.index) != list(target_ids):
        raise ValueError('Original target roster/order mismatch')
    if list(table.columns) != list(case_ids) or not case_ids:
        raise ValueError('Original training-case roster/order mismatch')
    values = table.to_numpy(dtype=float)
    if np.isinf(values).any():
        raise ValueError('Infinite source values are not NA')
    missing = np.isnan(values).sum(axis=1)
    # Integer comparison keeps exactly 30% and removes strictly greater than 30%.
    retained = missing * 10 <= len(case_ids) * 3
    return pd.DataFrame({'missing_count': missing, 'training_cases': len(case_ids),
                         'missing_fraction': missing / len(case_ids),
                         'retained': retained}, index=table.index)


def filter_bundle(original, policy):
    if list(policy.index) != list(original['Y_columns']):
        raise ValueError('Policy/model target roster mismatch')
    indices = np.flatnonzero(policy.retained.to_numpy())
    if not len(indices):
        raise ValueError('Missingness rule retains no targets')
    result = deepcopy(original)
    result['Y_columns'] = [original['Y_columns'][i] for i in indices]
    for key in ('sel_idx', 'coef'):
        if len(original[key]) != len(policy):
            raise ValueError(f'Incomplete original {key} roster')
        result[key] = [original[key][i].copy() for i in indices]
    if len(original['intercept']) != len(policy):
        raise ValueError('Incomplete original intercept roster')
    result['intercept'] = original['intercept'][indices].copy()
    for idx, coef, intercept in zip(result['sel_idx'], result['coef'], result['intercept']):
        if not len(idx) or len(idx) != len(coef) or not np.isfinite(intercept) or not np.isfinite(coef).all():
            raise ValueError('A retained target lacks a finite fitted model')
    metadata = result.setdefault('metadata', {})
    for key, values in metadata.get('per_target', {}).items():
        metadata['per_target'][key] = {t: values[t] for t in result['Y_columns'] if t in values}
    scores = np.asarray(list(metadata.get('per_target', {}).get('heldout_spearman', {}).values()))
    scores = scores[np.isfinite(scores)]
    metadata.update(heldout_median_spearman=round(float(np.median(scores)), 3) if scores.size else None,
                    n_imputable_sp_gt_0p3=int((scores > .3).sum()),
                    evaluation_status='Historical per-target CV for unchanged models; no new CV',
                    expression_scale='log2', log_base=2, log_pseudocount=0,
                    target_scale='Relative normalized log2 signal; linear representation is 2**Y',
                    missingness_release='AML_metabolite_NA30_20261006')
    return result


def validate_release(path, bundle):
    """Fail closed for the authorized release; historical explicit paths stay usable."""
    if bundle.get('metadata', {}).get('missingness_release') != 'AML_metabolite_NA30_20261006':
        if Path(path).name == 'rna2metabolite_aml_NA30_20261006_bundle.pkl.gz':
            raise ValueError('Release metadata missing')
        return
    manifest = json.loads(RELEASE_PATH.read_text())
    if sha256(path) != manifest['bundle_sha256']:
        raise ValueError('Missingness release bundle hash mismatch')
    if bundle['X_columns'] != manifest['RNA_genes'] or bundle['Y_columns'] != manifest['retained_targets']:
        raise ValueError('Missingness release identity mismatch')
    if bundle['metadata'].get('n_train_cases') != len(manifest['training_cases']):
        raise ValueError('Missingness release training roster mismatch')
    if bundle['metadata'].get('alpha') != 100 or bundle['metadata'].get('nfeat') != 1000:
        raise ValueError('Original ridge settings changed')
    if not np.isfinite(bundle['intercept']).all() or any(not np.isfinite(c).all() or not len(c) for c in bundle['coef']):
        raise ValueError('Missingness release includes unfitted/nonfinite targets')
