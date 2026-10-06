"""Regression coverage for selective H5AD reads and bounded split aggregation."""
import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from altanalyze3.components.fastComm.api import FastCommParams, _matrix_inputs_from_h5ad
from altanalyze3.components.fastComm.benchmark import _batched_state, _score_subset, FastCommBenchmarkParams
from altanalyze3.components.fastComm.scoring import make_state_pseudobulk


@pytest.mark.parametrize('dtype', [np.float32, np.float64])
@pytest.mark.parametrize('as_sparse', [True, False])
def test_batched_means_detection_preserve_dense_baseline(dtype, as_sparse):
    rng = np.random.default_rng(18)
    values = rng.normal(size=(113, 7)).astype(dtype)
    values[rng.random(values.shape) < .7] = 0
    values[2, 1] = np.nan
    values[3, 2:4] = [2, -2]  # Duplicate symbols must be averaged before detection.
    columns = ['A', 'B', 'C', 'C', 'D', 'E', 'F']
    index = [f'cell{i}' for i in range(len(values))]
    metadata = pd.DataFrame({'state': [' a ', 'b', '', 'rare'] + ['a', 'b', 'c'] * 36 + ['a']}, index=index)
    from altanalyze3.components.fastComm.api import _deduplicate_columns
    dense = _deduplicate_columns(pd.DataFrame(values, index=index, columns=columns))
    metadata = metadata.iloc[np.random.default_rng(9).permutation(len(metadata))[:100]]
    expected = make_state_pseudobulk(dense, metadata, state_key='state', min_cells=3)
    matrix = sparse.csr_matrix(values) if as_sparse else values
    actual = _batched_state(matrix, index, columns, metadata, state_key='state', min_cells=3, block_genes=2)
    pd.testing.assert_frame_equal(actual.expression, expected.expression, check_exact=True, check_names=False)
    pd.testing.assert_frame_equal(actual.detection, expected.detection, check_names=False)
    pd.testing.assert_series_equal(actual.sizes, expected.sizes)


@pytest.mark.parametrize('encoding', ['dense', 'csr', 'csc'])
@pytest.mark.parametrize('layer', [None, 'counts'])
def test_selective_reader_never_reads_unrequested_layers_or_raw(tmp_path, monkeypatch, encoding, layer):
    values = np.arange(60, dtype=np.float32).reshape(12, 5)
    matrix = values if encoding == 'dense' else getattr(sparse, encoding + '_matrix')(values)
    obj = ad.AnnData(matrix, obs=pd.DataFrame({'state': ['a', 'b'] * 6}, index=[f'c{i}' for i in range(12)]),
                     var=pd.DataFrame({'symbol': ['A', 'B', 'B', 'D', 'E']}, index=[f'g{i}' for i in range(5)]))
    obj.layers['counts'] = matrix * 2
    obj.layers['unrelated'] = matrix * 3
    obj.raw = obj.copy()
    path = tmp_path / 'input.h5ad'
    obj.write_h5ad(path)
    # A backed AnnData reader would load these slots. Poison them to enforce the
    # selective-read boundary rather than testing just a small RAM measurement.
    with h5py.File(path, 'a') as handle:
        del handle['raw']; handle.create_group('raw').attrs['encoding-type'] = 'do-not-read'
        del handle['layers/unrelated']; handle.create_group('layers/unrelated').attrs['encoding-type'] = 'do-not-read'
    monkeypatch.setattr(ad, 'read_h5ad', lambda *a, **k: pytest.fail('Eager AnnData reader invoked'))
    matrix, index, genes, metadata, info = _matrix_inputs_from_h5ad(
        FastCommParams(h5ad=path, gene_symbol_col='symbol', layer=layer), required_genes=['E', 'B'])
    np.testing.assert_array_equal((matrix.toarray() if sparse.issparse(matrix) else matrix), values[:, [1, 2, 4]] * (2 if layer else 1))
    assert genes == ['B', 'B', 'E']
    assert index == list(obj.obs_names) and metadata.equals(obj.obs)
    assert info['n_loaded_columns'] == 3



def test_batched_state_never_densifies_all_features(monkeypatch):
    import altanalyze3.components.fastComm.benchmark as benchmark
    original = benchmark.to_dense_frame
    seen = []
    def bounded(matrix, **kwargs):
        assert matrix.shape[1] <= 8
        seen.extend(kwargs['columns'])
        return original(matrix, **kwargs)
    monkeypatch.setattr(benchmark, 'to_dense_frame', bounded)
    index = [f'c{i}' for i in range(61)]
    columns = [f'g{i}' for i in range(29)]
    state = _batched_state(sparse.csr_matrix(np.ones((61, 29))), index, columns,
                           pd.DataFrame({'state': ['a'] * 61}, index=index),
                           state_key='state', min_cells=1, block_genes=8)
    assert seen == columns
    assert state.expression.eq(1).all().all() and state.detection.eq(1).all().all()
