import anndata as ad
import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.visualization.scalable_viewer import bundle as B, precompute as P
from altanalyze3.components.cellHarmony.flask.pipeline import _attach_imputed_expression_metadata


@pytest.mark.parametrize('layout', ['ordered', 'reordered_missing', 'all_zero', 'per_state'])
def test_dense_h5_block_builder_matches_reference(tmp_path, monkeypatch, layout):
    rng = np.random.default_rng(19)
    matrix = rng.normal(size=(211, 37)).astype(np.float32)
    matrix[rng.random(matrix.shape) < .6] = 0
    if layout == 'all_zero':
        matrix[:] = 0
    rows = [f'cell{i}' for i in range(len(matrix))]
    cells = rows.copy()
    states = ['a', 'b', 'c']
    codes = np.arange(len(cells)) % 3
    counts = np.bincount(codes, minlength=3)
    if layout == 'reordered_missing':
        cells = cells[::-1]
        cells[0] = 'missing'
    if layout == 'per_state':
        matrix = matrix[:3]
        rows = states
    obj = ad.AnnData(X=matrix, obs=pd.DataFrame(index=rows),
                     var=pd.DataFrame(index=[f'f{i}' for i in range(matrix.shape[1])]))
    # These unused fields must not enter the block builder.
    obj.layers['counts'] = np.ones_like(matrix)
    obj.obsm['unused'] = np.zeros((len(rows), 5), dtype=np.float32)
    source = tmp_path / 'input.h5ad'
    obj.write_h5ad(source)
    outputs = []
    for name in ['fast', 'reference']:
        folder = tmp_path / name
        folder.mkdir()
        paths = B.BundlePaths(str(folder), 'test')
        if name == 'reference':
            original = P._read_feature_matrix
            monkeypatch.setattr(P, '_read_feature_matrix', lambda path, **kwargs: original(path))
        info = P.ingest_modality('test', str(source), paths=paths, barcodes=cells,
                                states=states, state_code=codes, state_n=counts)
        outputs.append((info, paths.modality('test')))
    assert outputs[0][0] == outputs[1][0]
    for key in ['stats_mean', 'stats_frac'] + ([] if layout == 'per_state' else ['expr_indptr', 'expr_indices', 'expr_data']):
        arrays = [np.load(getattr(paths, key)) for _, paths in outputs]
        # Empty stores may carry an unused sentinel in the legacy builder.
        if key in ['expr_indices', 'expr_data']:
            arrays = [x[:outputs[0][0]['nnz']] for x in arrays]
        np.testing.assert_array_equal(*arrays)


@pytest.mark.parametrize('scale,base', [('linear', None), ('log2', 2), ('log1p', 'e'), ('log1p', 10)])
def test_count_conversion_matches_full_float64_transform(scale, base):
    rng = np.random.default_rng(19)
    matrix = rng.normal(size=(101, 53)).astype(np.float32)
    obj = ad.AnnData(X=matrix)
    _attach_imputed_expression_metadata(obj, matrix, expression_scale=scale, log_base=base)
    values = matrix.astype(np.float64)
    expected = values if scale == 'linear' else (np.exp2(values) - 1 if scale == 'log2'
                  else np.expm1(values) if base == 'e' else np.power(base, values) - 1)
    np.testing.assert_array_equal(obj.layers['counts'], np.maximum(expected, 0).astype(np.float32))
