"""Predicted modalities retain their scale; RNA input guards still apply."""
import numpy as np
import pandas as pd
import anndata as ad
import pytest
from altanalyze3.components.cellHarmony.flask import pipeline


@pytest.mark.parametrize('modality', ['lipids', 'adt', 'grn_tf'])
def test_signed_predictions_make_markers_without_count_normalization(tmp_path, modality):
    rng = np.random.default_rng(13)
    groups = np.repeat(['A', 'B'], 40)
    values = rng.normal(0, .05, (80, 2))
    values[:40] += [-1., -5.]
    values[40:] += [-5., -1.]
    a = ad.AnnData(values.copy(), obs=pd.DataFrame({'state': groups}, index=[str(i) for i in range(80)]),
                   var=pd.DataFrame(index=['PC(16:0/18:1)', 'PE(18:0/20:4)']))
    # The RNA default must still reject this signed matrix.
    with pytest.raises(ValueError, match='input rejected'):
        pipeline.marker_mod._markerfinder_stats(a, 'state', False, None)
    result = pipeline._emit_modality_marker_heatmap(modality, a, tmp_path, 'state', {'species': 'human'})
    assert result['enabled']
    centroids = pd.read_csv(result['centroids_tsv'], sep='\t', index_col=0)
    expected = pd.DataFrame({'A': values[:40].mean(axis=0), 'B': values[40:].mean(axis=0)}, index=a.var_names)
    np.testing.assert_allclose(centroids.loc[expected.index, expected.columns], expected, rtol=5e-4)
    np.testing.assert_array_equal(a.X, values)
    markers = pd.read_csv(result['markers_tsv'], sep='\t')
    assert set(markers.Gene) == set(a.var_names)
    assert (markers.rho > .99).all()
    aggregates = pipeline.marker_mod._prepare_marker_stats_aggregates(a, 'state', list(a.var_names), False, None)
    # Default RNA/count centroid export keeps its guard as well.
    with pytest.raises(ValueError, match='Centroid matrix cannot be built'):
        pipeline.marker_mod._build_marker_centroids_from_aggregates(list(a.var_names), ['A', 'B'], aggregates)


@pytest.mark.parametrize('modality', ['lipids', 'adt', 'grn_tf', 'metabolite', 'lipid'])
def test_empty_predicted_markers_are_nonfatal(tmp_path, modality):
    a = ad.AnnData(np.ones((40, 2), dtype=np.float32),
                   obs=pd.DataFrame({'state': ['A'] * 20 + ['B'] * 20}, index=list(map(str, range(40)))),
                   var=pd.DataFrame(index=['constant_A', 'constant_B']))
    result = pipeline._emit_modality_marker_heatmap(modality, a, tmp_path, 'state', {})
    assert result['enabled'] is False
    assert 'differential views remain available' in result['message']
    assert 'heatmap_pdf' not in result


def test_unrelated_marker_errors_are_not_suppressed(tmp_path, monkeypatch):
    def broken(*args, **kwargs):
        raise ValueError('invalid cluster data')
    monkeypatch.setattr(pipeline.marker_mod, 'generate_marker_heatmap_from_adata', broken)
    with pytest.raises(ValueError, match='invalid cluster data'):
        pipeline._emit_modality_marker_heatmap('lipids', None, tmp_path, 'state', {})


def test_marker_correlations_match_scipy_independent_of_cell_order():
    from scipy.stats import pearsonr
    from altanalyze3.components.cellHarmony.markerFinder import marker_finder
    rng = np.random.default_rng(9)
    groups = np.repeat(['A', 'B'], 4000)
    signal = np.r_[np.zeros(4000), np.ones(4000)]
    values = np.column_stack([
        20 + .001 * signal + rng.normal(0, .0001, 8000),
        np.full(8000, 20.1234),
        signal,
        rng.normal(size=8000),
    ]).astype(np.float32)
    names = ['low_variance', 'constant', 'perfect', 'noise']
    r, p = marker_finder(values, groups, names)
    reverse_r, reverse_p = marker_finder(values[::-1], groups[::-1], names)
    assert 'constant' not in r.index
    assert np.isfinite(p).all().all()
    assert (r.abs() <= 1).all().all()
    np.testing.assert_allclose(r, reverse_r.loc[r.index, r.columns], atol=1e-8)
    np.testing.assert_allclose(p, reverse_p.loc[p.index, p.columns], atol=1e-8)
    for name in r.index:
        expected_r, expected_p = pearsonr(values[:, names.index(name)].astype(float), signal)
        assert r.loc[name, 'B'] == pytest.approx(expected_r, abs=1e-8)
        assert p.loc[name, 'B'] == pytest.approx(expected_p, abs=1e-8)


def test_missing_pvalue_does_not_erase_other_fdrs():
    adjusted = pipeline.marker_mod._bh_fdr([.01, np.nan, .04, .03])
    np.testing.assert_allclose(adjusted, [.03, np.nan, .04, .04], equal_nan=True)


@pytest.mark.parametrize('architecture', ['lipidwise', 'single_model'])
def test_lipid_prediction_preserves_model_feature_names(architecture):
    import warnings
    from pathlib import Path
    from sklearn.linear_model import LinearRegression
    from sklearn.preprocessing import StandardScaler
    from altanalyze3.components.rna2lipid.api import Rna2LipidBundle
    frame = pd.DataFrame({'A': [1., 4., 2., 3.], 'B': [2., 1., 3., 4.]})
    scaler = StandardScaler().fit(frame)
    transformed = pd.DataFrame(scaler.transform(frame), columns=frame.columns)
    target = np.array([1., 2., 4., 3.])
    estimator = LinearRegression().fit(transformed, target)
    kwargs = ({'models': {'PC': {'genes': ['A', 'B'], 'model': estimator}}}
              if architecture == 'lipidwise' else {'model': estimator})
    bundle = Rna2LipidBundle(bundle_path=Path('test.pkl'), scaler_x=scaler,
                            input_genes=['A', 'B'], output_lipids=['PC'], **kwargs)
    with warnings.catch_warnings(record=True) as emitted:
        result = bundle.predict_from_dataframe(frame)
    assert not [w for w in emitted if 'valid feature names' in str(w.message)]
    np.testing.assert_allclose(result.predictions['PC'], estimator.predict(transformed))


def test_lipid_predictions_are_floored_before_export_and_analysis(monkeypatch):
    from types import SimpleNamespace
    query = ad.AnnData(np.ones((2, 2)), obs=pd.DataFrame(index=['cell1', 'cell2']))
    predictions = pd.DataFrame([[-3.5, 20.], [1.5, -0.01]], index=query.obs_names,
                               columns=['PC(16:0/18:1)', 'PE(18:0/20:4)'])
    fake = SimpleNamespace(predict_from_adata=lambda a: SimpleNamespace(predictions=predictions, summary={}))
    monkeypatch.setattr(pipeline, 'load_rna2lipid_bundle', lambda *args: fake)
    result, summary = pipeline._build_imputed_lipid_adata(query)
    np.testing.assert_array_equal(result.X, [[0., 20.], [1.5, 0.]])
    assert (result.layers['counts'] >= 0).all()
    assert summary['clipped_negative_values'] == 2
    assert result.uns['prediction_summary']['clipped_negative_values'] == 2
    assert predictions.iloc[0, 0] == -3.5  # source/model predictions are not mutated
    assert list(result.var_names) == list(predictions.columns)


def test_single_cluster_has_no_defined_marker_correlations():
    from altanalyze3.components.cellHarmony.markerFinder import marker_finder
    values = np.random.default_rng(17).normal(20, .01, (4000, 3)).astype(np.float32)
    r, p = marker_finder(values, ['A'] * len(values), ['a', 'b', 'c'])
    assert r.empty and p.empty
