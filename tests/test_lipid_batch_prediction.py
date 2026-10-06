"""Exact deployed model coverage and prediction parity for batched evaluation."""
import numpy as np
import pandas as pd
import pytest
from threadpoolctl import threadpool_limits

from altanalyze3.components.rna2lipid.api import load_bundle


@pytest.mark.parametrize('threads', [1, 4])
@pytest.mark.parametrize('rows', [1, 31, 2048])
def test_all_production_lipids_match_the_original_estimator_calls(monkeypatch, threads, rows):
    bundle = load_bundle()
    assert len(bundle.input_genes) == 1303 and len(bundle.output_lipids) == 202
    assert set(bundle.models) == set(bundle.output_lipids)
    panel = (bundle.input_genes, bundle.output_lipids)
    values = pd.DataFrame(np.random.default_rng(61).uniform(.1, 8, (rows, len(bundle.input_genes))),
                          columns=bundle.input_genes, index=[f'cell{i}' for i in range(rows)])
    pack = bundle._lipidwise_linear_parameters()
    assert pack is not None and pack[0].shape == (1303, 202)
    with threadpool_limits(limits=threads):
        optimized = bundle._predict_aligned_matrix(values)
        # The original per-lipid estimator.predict implementation remains the
        # executable numerical baseline and fallback for unsupported models.
        monkeypatch.setattr(bundle, '_lipidwise_linear_parameters', lambda: None)
        baseline = bundle._predict_aligned_matrix(values)
    pd.testing.assert_frame_equal(optimized, baseline, atol=1e-12, rtol=1e-12)
    assert tuple(optimized.columns) == panel[1]
    assert optimized.index.equals(values.index)
    assert panel == (bundle.input_genes, bundle.output_lipids)


def test_linear_pack_does_not_cache_updated_model_coefficients():
    bundle = load_bundle()
    first = bundle.output_lipids[0]
    weights, _ = bundle._lipidwise_linear_parameters()
    estimator = bundle.models[first]['model']
    estimator.coef_[0] += 1.
    updated, _ = bundle._lipidwise_linear_parameters()
    assert not np.array_equal(weights, updated)


def test_mismatched_feature_names_keep_the_original_validation(monkeypatch):
    bundle = load_bundle()
    estimator = bundle.models[bundle.output_lipids[0]]['model']
    if not hasattr(estimator, 'feature_names_in_'):
        estimator.feature_names_in_ = np.asarray(bundle.models[bundle.output_lipids[0]]['genes'])
    estimator.feature_names_in_ = estimator.feature_names_in_.copy()
    estimator.feature_names_in_[0] = 'unresolved_test_feature'
    assert bundle._lipidwise_linear_parameters() is None
    with pytest.raises(RuntimeError, match='Prediction failed specifically'):
        bundle._predict_aligned_matrix(np.zeros((1, 1303)))


@pytest.mark.parametrize('panel', ['input_genes', 'output_lipids'])
def test_duplicate_model_panel_names_use_original_prediction_validation(panel):
    bundle = load_bundle()
    names = getattr(bundle, panel)
    setattr(bundle, panel, (names[0], names[0], *names[2:]))
    assert bundle._lipidwise_linear_parameters() is None
