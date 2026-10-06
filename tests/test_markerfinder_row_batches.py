"""Changing storage-based batch sizes preserves every cell's marker statistics."""
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from altanalyze3.components.cellHarmony.markerFinder import _marker_finder_rows


@pytest.mark.parametrize('encoding', ['csr', 'csc', 'dense'])
def test_row_batches_match_original_statistics_and_marker_ranking(encoding):
    rng = np.random.default_rng(124)
    cells, genes = 6000, 2101
    values = rng.random((cells, genes), dtype=np.float32)
    values[rng.random(values.shape) < .8] = 0
    values[:, -1] = 0
    groups = np.asarray([f'group{i % 17}' for i in range(cells)])
    groups[-3:] = 'rare'
    names = pd.Index([f'gene{i}' for i in range(genes)])
    matrix = {'csr': sparse.csr_matrix, 'csc': sparse.csc_matrix,
              'dense': np.asarray}[encoding](values)
    original_rows = max(1, min(8192, 8_000_000 // genes))
    if original_rows >= 256:
        original_rows = original_rows // 256 * 256
    expected = _marker_finder_rows(matrix, groups, names, row_block_size=original_rows)
    actual = _marker_finder_rows(matrix, groups, names)
    for before, after in zip(expected, actual):
        pd.testing.assert_frame_equal(before, after, rtol=1e-10, atol=1e-12)
    assert list(actual[0].columns) == list(expected[0].columns)
    assert 'rare' in actual[0].columns
    for group in actual[0]:
        assert actual[0][group].sort_values(ascending=False).index.tolist() == expected[0][group].sort_values(ascending=False).index.tolist()


def test_empty_sparse_rows_keep_their_group_membership():
    matrix = sparse.csr_matrix(([1., 2., 3., 4.], [0, 1, 0, 1], [0, 2, 2, 4]), shape=(3, 2))
    expected = _marker_finder_rows(matrix, ['A', 'A', 'B'], ['x', 'y'], row_block_size=1)
    actual = _marker_finder_rows(matrix, ['A', 'A', 'B'], ['x', 'y'])
    for before, after in zip(expected, actual):
        pd.testing.assert_frame_equal(before, after, rtol=1e-10, atol=1e-12)
