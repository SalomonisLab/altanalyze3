"""Verify the approved target rule and unchanged learned predictions."""
import gzip
import json
import pickle

import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.rna2metabolite.missingness_release import (
    RELEASE_PATH, filter_bundle, sha256, target_policy, validate_release,
)


def test_strict_missingness_boundary_and_observed_zero_negative_values():
    values = np.full((4, 100), -2.0)
    values[0, :30] = np.nan
    values[1, :31] = np.nan
    values[2, :] = np.nan
    values[3, :] = 0.0
    table = pd.DataFrame(values, index=['at30', 'over30', 'all_missing', 'zero'])
    policy = target_policy(table, list(table.index), list(table.columns))
    assert policy.retained.tolist() == [True, False, False, True]
    assert policy.missing_count.tolist() == [30, 31, 100, 0]


def test_training_denominator_and_identity_gates():
    table = pd.DataFrame([[np.nan] * 25 + [1.] * 59], index=['A'])
    assert target_policy(table, ['A'], list(table.columns)).retained.tolist() == [True]
    table.iloc[0, 25] = np.nan
    assert target_policy(table, ['A'], list(table.columns)).retained.tolist() == [False]
    with pytest.raises(ValueError, match='target roster'):
        target_policy(table, ['B'], list(table.columns))
    with pytest.raises(ValueError, match='case roster'):
        target_policy(table, ['A'], list(reversed(table.columns)))
    table.columns = [0] * 84
    with pytest.raises(ValueError, match='Duplicate'):
        target_policy(table, ['A'], list(table.columns))


def test_filter_retains_learned_values_and_rejects_unfitted_retained_target():
    model = dict(X_columns=['g1'], Y_columns=['A', 'B'],
                 mu=np.array([.5]), sd=np.array([2.]),
                 sel_idx=[np.array([0]), np.array([], dtype=int)],
                 coef=[np.array([1.25]), np.array([])], intercept=np.array([3., np.nan]),
                 metadata={'per_target': {'heldout_spearman': {'A': .4, 'B': .2}}})
    policy = pd.DataFrame({'retained': [True, False]}, index=['A', 'B'])
    result = filter_bundle(model, policy)
    assert result['Y_columns'] == ['A']
    np.testing.assert_array_equal(result['mu'], model['mu'])
    np.testing.assert_array_equal(result['sd'], model['sd'])
    np.testing.assert_array_equal(result['coef'][0], model['coef'][0])
    assert model['Y_columns'] == ['A', 'B']  # baseline is not edited
    assert result['metadata']['per_target']['heldout_spearman'] == {'A': .4}
    policy.iloc[1, 0] = True
    with pytest.raises(ValueError, match='finite fitted'):
        filter_bundle(model, policy)


def test_complete_real_release_matches_original_predictions_and_authorized_rosters(tmp_path):
    manifest = json.loads(RELEASE_PATH.read_text())
    assert len(manifest['training_cases']) == 84
    assert len(manifest['RNA_genes']) == 12416
    assert len(manifest['original_targets']) == 2533
    assert len(manifest['retained_targets']) == 2023
    assert len(manifest['removed_targets']) == 510
    baseline_path = RELEASE_PATH.parent / 'artifacts/rna2metabolite_aml_bundle.pkl.gz'
    bundle_path = RELEASE_PATH.parent / 'artifacts/rna2metabolite_aml_NA30_20261006_bundle.pkl.gz'
    with gzip.open(baseline_path, 'rb') as handle:
        old = pickle.load(handle)
    with gzip.open(bundle_path, 'rb') as handle:
        new = pickle.load(handle)
    assert sha256(baseline_path) == manifest['baseline_bundle_sha256']
    validate_release(bundle_path, new)
    positions = {t: i for i, t in enumerate(old['Y_columns'])}
    z = np.random.default_rng(17).normal(size=(3, len(old['X_columns'])))
    for j, target in enumerate(new['Y_columns']):
        i = positions[target]
        np.testing.assert_array_equal(new['sel_idx'][j], old['sel_idx'][i])
        np.testing.assert_array_equal(new['coef'][j], old['coef'][i])
        np.testing.assert_array_equal(z[:, new['sel_idx'][j]] @ new['coef'][j] + new['intercept'][j],
                                      z[:, old['sel_idx'][i]] @ old['coef'][i] + old['intercept'][i])
    altered = tmp_path / 'altered.pkl.gz'
    altered.write_bytes(b'changed')
    with pytest.raises(ValueError, match='hash mismatch'):
        validate_release(altered, new)


def test_api_and_scalable_use_finite_release_and_verified_inverse():
    import anndata as ad
    from altanalyze3.components.rna2metabolite.api import load_bundle
    from altanalyze3.components.cellHarmony.imputed_scale import prediction_encoding, uses_imputed_scale
    from altanalyze3.components.cellHarmony.flask.pipeline import _finalize_imputed_adata
    model = load_bundle()
    assert len(model.targets) == 2023
    query = ad.AnnData(np.array([[1.], [2.]], dtype=np.float32),
                      obs=pd.DataFrame(index=['cell1', 'cell2']),
                      var=pd.DataFrame(index=[model.input_genes[0]]))
    result = model.predict_from_adata(query)
    assert result.predictions.shape == (2, 2023)
    assert np.isfinite(result.predictions.to_numpy()).all()
    exported = model.impute_anndata(query)
    assert uses_imputed_scale(exported)
    assert prediction_encoding(exported) == ('log2', 2.0, 0.0)
    np.testing.assert_allclose(exported.X, result.predictions, rtol=1e-6, atol=1e-6)
    config = json.loads((RELEASE_PATH.parent.parent / 'cellHarmony/flask/reference_config.json').read_text())
    selections = [ref['impute_config']['metabolite'] for species in config['species']
                  for ref in species['references'] if 'metabolite' in ref.get('impute_config', {})]
    assert len(selections) == 1
    cfg = selections[0]
    assert cfg['bundle_path'].endswith(model.bundle_path.name)
    assert cfg['log_pseudocount'] == 0
    scalable, _ = _finalize_imputed_adata(query, result.predictions, modality_id='metabolite',
        feature_label='Metabolite', feature_type='Metabolite', expression_scale=cfg['expression_scale'],
        log_base=cfg['log_base'], log_pseudocount=cfg['log_pseudocount'], base_summary=result.summary)
    np.testing.assert_allclose(scalable.layers['counts'], np.exp2(scalable.X), rtol=1e-6)
    assert scalable.n_vars == 2023
