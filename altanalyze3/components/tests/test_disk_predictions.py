import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.cellHarmony.flask.disk_predictions import PredictionWriter


@pytest.mark.parametrize('negative', [False, True])
def test_prediction_blocks_are_normal_h5ad(tmp_path, negative):
    values = np.arange(39, dtype=np.float32).reshape(13, 3)
    if negative:
        values[9, 1] = -1
    obs = pd.DataFrame({'state': ['a'] * 13}, index=[f'c{i}' for i in range(13)])
    var = pd.DataFrame(index=['f1', 'f2', 'f3'])
    path = tmp_path / 'predictions.h5ad'
    writer = PredictionWriter(path, obs, var, {'expression_scale': 'linear'}, {})
    for start in range(0, 13, 4):
        writer.append(values[start:start + 4])
    writer.close()
    loaded = ad.read_h5ad(path)
    np.testing.assert_array_equal(loaded.X, values)
    np.testing.assert_array_equal(loaded.layers['counts'], np.maximum(values, 0))
    assert list(loaded.obs_names) == list(obs.index)
    with h5py.File(path) as handle:
        assert (handle['X'].id == handle['layers/counts'].id) is (not negative)


def test_incomplete_prediction_fails_and_closes_file(tmp_path):
    path = tmp_path / 'partial.h5ad'
    writer = PredictionWriter(path, pd.DataFrame(index=['a', 'b']),
                              pd.DataFrame(index=['f']), {}, {})
    writer.append(np.ones((1, 1), dtype=np.float32))
    with pytest.raises(ValueError, match='incomplete'):
        writer.close()
    assert not writer.handle.id.valid
    writer.abort()
    assert not path.exists()

@pytest.mark.parametrize('scale,base', [('log2', '2'), ('log1p', 'e'), ('log1p', '10')])
def test_logged_prediction_counts_match_metadata_conversion(tmp_path, scale, base):
    from altanalyze3.components.cellHarmony.flask.pipeline import _attach_imputed_expression_metadata
    values = np.linspace(-0.2, 3, 77, dtype=np.float32).reshape(11, 7)
    obs, var = pd.DataFrame(index=[f'c{i}' for i in range(11)]), pd.DataFrame(index=[f'f{i}' for i in range(7)])
    expected = ad.AnnData(X=values.copy(), obs=obs, var=var)
    _attach_imputed_expression_metadata(expected, values, expression_scale=scale, log_base=base)
    path = tmp_path / 'logged.h5ad'
    writer = PredictionWriter(path, obs, var, dict(expected.uns), {}, expression_scale=scale, log_base=base)
    for start in range(0, 11, 3):
        writer.append(values[start:start + 3])
    writer.close()
    actual = ad.read_h5ad(path)
    np.testing.assert_array_equal(actual.layers['counts'], expected.layers['counts'])

@pytest.mark.parametrize('scale,base', [('linear', None), ('log2', '2'), ('log1p', 'e')])
def test_streamed_broadcast_profiles_preserve_cell_level_values(tmp_path, scale, base):
    from altanalyze3.components.cellHarmony.flask.pipeline import _stream_pseudobulk_viewer, _pseudobulk_group_key, _finalize_imputed_adata
    from altanalyze3.components.cellHarmony.disk_differential import read_rows, bounded_materialize
    obs = pd.DataFrame({'Library': np.tile(['s1', 's2', 's3'], 39),
                        'state': np.tile(['a', 'b', 'b', 'a', 'c', 'a', 'b', 'c', 'b'], 13)},
                       index=[f'c{i}' for i in range(117)])
    query = ad.AnnData(X=np.ones((117, 1), dtype=np.float32), obs=obs)
    keys = _pseudobulk_group_key(obs, 'state')
    groups = pd.Index(keys.unique())
    values = np.random.default_rng(19).uniform(0, 3, (len(groups), 37)).astype(np.float32)
    predictions = pd.DataFrame(values, index=groups, columns=[f'f{i}' for i in range(37)])
    frame = predictions.reindex(keys)
    frame.index = obs.index
    expected, _ = _finalize_imputed_adata(query, frame, modality_id='metabolite', feature_label='metabolite',
                                         feature_type='metabolite', expression_scale=scale, log_base=base, base_summary={})
    path = tmp_path / 'broadcast.h5ad'
    disk, _ = _stream_pseudobulk_viewer(query, predictions, 'state', {'expression_scale': scale, 'log_base': base},
                                       'metabolite', {}, path, 'lzf')
    try:
        physical = ad.read_h5ad(path)
        np.testing.assert_array_equal(physical.X, expected.X)
        np.testing.assert_array_equal(physical.layers['counts'], expected.layers['counts'])
        selected = np.array([9, 1, 2, 67, 31])
        np.testing.assert_array_equal(read_rows(disk.X, selected), expected.X[selected])
        np.testing.assert_array_equal(read_rows(disk.layers['counts'], selected), expected.layers['counts'][selected])
        memory = bounded_materialize(disk)
        np.testing.assert_array_equal(memory.X, expected.X)
        np.testing.assert_array_equal(memory.layers['counts'], expected.layers['counts'])
    finally:
        disk._analysis_h5_handle.close()
