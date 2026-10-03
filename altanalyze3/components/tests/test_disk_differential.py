import anndata as ad
import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.cellHarmony import cellHarmony_differential as D
from altanalyze3.components.cellHarmony.flask.pipeline import _read_differential_h5ad


@pytest.mark.parametrize('sparse', [False, True])
@pytest.mark.parametrize('scale', ['logged', 'linear', 'raw_large'])
def test_disk_rank_and_fold_match_memory(tmp_path, scale, sparse):
    rng = np.random.default_rng(41)
    x = rng.uniform(0, 3, size=(97, 701)).astype(np.float32)
    x[:48, :51] *= 2
    if scale == 'raw_large':
        x *= 20
    obj = ad.AnnData(X=x, obs=pd.DataFrame({'group': ['case'] * 48 + ['control'] * 49},
                                          index=[f'c{i}' for i in range(len(x))]),
                     var=pd.DataFrame(index=[f'g{i}' for i in range(x.shape[1])]))
    obj.layers['counts'] = np.expm1(x) if scale == 'logged' else x.copy()
    if scale == 'logged':
        obj.uns['log1p'] = {'base': None}
    if sparse:
        import scipy.sparse as sp
        obj.X = sp.csr_matrix(obj.X)
        obj.layers['counts'] = sp.csr_matrix(obj.layers['counts'])
    path = tmp_path / 'source.h5ad'
    obj.write_h5ad(path)
    disk = _read_differential_h5ad(path, disk_backed=True)
    try:
        for selected in [np.ones(len(x), dtype=bool), np.arange(len(x)) % 3 != 0]:
            frames = []
            folds = []
            for matrix in [obj[selected].copy(), disk[selected]]:
                names, fdr, fc, p = D._rank_genes_scanpy(matrix, 'group', 'case', 'control', 'wilcoxon')
                frames.append(pd.DataFrame({'fdr': fdr.values, 'fc': fc.values, 'p': p.values}, index=names).sort_index())
                folds.append(D._compute_pseudobulk_log2fc(matrix, 'group', 'case', 'control').sort_index())
            np.testing.assert_allclose(frames[0], frames[1], rtol=1e-12, atol=1e-12)
            np.testing.assert_allclose(folds[0], folds[1], rtol=1e-12, atol=1e-12)
    finally:
        disk._analysis_h5_handle.close()


@pytest.mark.parametrize('sparse', [False, True])
def test_disk_differential_export_preserves_selected_values(tmp_path, sparse):
    x = np.arange(60, dtype=np.float32).reshape(12, 5)
    obj = ad.AnnData(X=x, obs=pd.DataFrame(index=[f'c{i}' for i in range(12)]),
                     var=pd.DataFrame(index=[f'g{i}' for i in range(5)]), layers={'counts': x * 2})
    if sparse:
        import scipy.sparse as sp
        obj.X = sp.csr_matrix(obj.X)
        obj.layers['counts'] = sp.csr_matrix(obj.layers['counts'])
    source, target = tmp_path / 'source.h5ad', tmp_path / 'result.h5ad'
    obj.write_h5ad(source)
    disk = _read_differential_h5ad(source, disk_backed=True)
    try:
        store = {'detailed_deg': pd.DataFrame({'gene': ['g1', 'g3']})}
        D.write_differentials_only_h5ad(disk[np.arange(12) % 2 == 0], store, target)
        loaded = ad.read_h5ad(target)
        np.testing.assert_array_equal(loaded.X.toarray() if sparse else loaded.X, x[::2][:, [1, 3]])
        np.testing.assert_array_equal(loaded.layers['counts'].toarray() if sparse else loaded.layers['counts'], x[::2][:, [1, 3]] * 2)
    finally:
        disk._analysis_h5_handle.close()


def test_moderated_statistics_match_full_matrix_formula(tmp_path):
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    rng = np.random.default_rng(73)
    x = rng.uniform(0, 8, (47, 701)).astype(np.float32)
    obj = ad.AnnData(X=x, obs=pd.DataFrame({'group': ['case'] * 23 + ['control'] * 24},
                                          index=[f'c{i}' for i in range(47)]),
                     var=pd.DataFrame(index=[f'g{i}' for i in range(701)]))
    case, control = x[:23], x[23:]
    diff = case.mean(0) - control.mean(0)
    variance = (case.var(0, ddof=1) + control.var(0, ddof=1)) / 2
    shrunk = .8 * variance + .2 * np.median(variance)
    t = diff / np.maximum(np.sqrt(shrunk / 23 + shrunk / 24), 1e-12)
    p = 2 * stats.t.sf(np.abs(t), df=45)
    keep = D._independent_filter_mask(np.vstack([case, control])) & np.isfinite(p)
    fdr = np.ones_like(p)
    fdr[keep] = multipletests(p[keep], method='fdr_bh')[1]
    path = tmp_path / 'moderated.h5ad'
    obj.write_h5ad(path)
    disk = _read_differential_h5ad(path, disk_backed=True)
    try:
        view = obj[:]
        for matrix in (obj, disk, view):
            result, _ = D._moderated_t_test(matrix, 'group', 'case', 'control', 'all', store_full_log=False)
            result = result.set_index('gene').reindex(obj.var_names)
            np.testing.assert_allclose(result['log2fc'], diff, rtol=1e-6, atol=1e-6)
            np.testing.assert_allclose(result['t'], t, rtol=1e-6, atol=1e-6)
            np.testing.assert_allclose(result['pval'], p, rtol=1e-5, atol=1e-8)
            np.testing.assert_allclose(result['fdr'], fdr, rtol=1e-5, atol=1e-8)
        assert view.is_view
    finally:
        disk._analysis_h5_handle.close()


def test_markerfinder_disk_blocks_preserve_statistics(tmp_path):
    import h5py
    import scipy.sparse as sp
    from altanalyze3.components.cellHarmony.markerFinder import marker_finder
    from altanalyze3.components.visualization.marker_heatmap_h5ad import _prepare_marker_stats_aggregates
    rng = np.random.default_rng(74)
    x = rng.uniform(0, 8, (1203, 71)).astype(np.float32)
    groups = np.array(['b', 'a', 'c'])[np.arange(len(x)) % 3]
    names = [f'g{i}' for i in range(x.shape[1])]
    obj = ad.AnnData(X=x, obs=pd.DataFrame({'group': groups}, index=[f'c{i}' for i in range(len(x))]),
                     var=pd.DataFrame(index=names))
    path = tmp_path / 'markers.h5ad'
    obj.write_h5ad(path)
    expected = marker_finder(sp.csr_matrix(x), groups, names)
    blocked = marker_finder(sp.csr_matrix(x), groups, names, feature_block_size=13)
    with h5py.File(path) as handle:
        actual = marker_finder(handle['X'], groups, names)
        for a, b, c in zip(expected, blocked, actual):
            np.testing.assert_allclose(a, b, rtol=1e-11, atol=1e-12)
            np.testing.assert_allclose(a, c, rtol=1e-11, atol=1e-12)
        disk = _read_differential_h5ad(path, disk_backed=True)
        try:
            selected = ['g12', 'g3', 'g65', 'missing', 'g12']
            a = _prepare_marker_stats_aggregates(obj, 'group', selected, False, None)
            b = _prepare_marker_stats_aggregates(disk, 'group', selected, False, None)
            for key in ('sums', 'counts', 'total_sum'):
                np.testing.assert_allclose(a[key], b[key], rtol=1e-6, atol=1e-6)
        finally:
            disk._analysis_h5_handle.close()


def test_streaming_markerfinder_large_sparse_preserves_small_variances():
    import scipy.sparse as sp
    from altanalyze3.components.cellHarmony.markerFinder import marker_finder
    rng = np.random.default_rng(75)
    x = (1000 + rng.normal(0, .01, (50_001, 11))).astype(np.float32)
    x[:, 0] = 1000
    groups = np.where(np.arange(len(x)) % 2, 'a', 'b')
    matrix = sp.csr_matrix(x)
    expected = marker_finder(matrix, groups, feature_block_size=11)
    actual = marker_finder(matrix, groups)
    for a, b in zip(expected, actual):
        np.testing.assert_allclose(a, b, rtol=1e-6, atol=2e-8)
    assert '0' not in actual[0].index


def test_sparse_disk_export_no_degs_stays_sparse(tmp_path):
    import scipy.sparse as sp
    x = sp.random(73, 701, density=.05, format='csr', random_state=71, dtype=np.float32)
    obj = ad.AnnData(X=x, obs=pd.DataFrame(index=[f'c{i}' for i in range(73)]),
                     var=pd.DataFrame(index=[f'g{i}' for i in range(701)]))
    source, target = tmp_path / 'source.h5ad', tmp_path / 'all.h5ad'
    obj.write_h5ad(source)
    disk = _read_differential_h5ad(source, disk_backed=True)
    try:
        D.write_differentials_only_h5ad(disk, {'detailed_deg': pd.DataFrame(columns=['gene'])}, target)
        actual = ad.read_h5ad(target)
        assert sp.issparse(actual.X)
        np.testing.assert_array_equal(actual.X.toarray(), x.toarray())
    finally:
        disk._analysis_h5_handle.close()


def test_anonymous_row_store_preserves_selection_and_shared_counts(tmp_path):
    import h5py
    from altanalyze3.components.cellHarmony.disk_differential import bind_row_store, read_rows
    x = np.random.default_rng(95).uniform(0, 3, (117, 701)).astype(np.float32)
    path = tmp_path / 'cache.h5ad'
    ad.AnnData(X=x, obs=pd.DataFrame(index=[f'c{i}' for i in range(117)]),
               var=pd.DataFrame(index=[f'g{i}' for i in range(701)])).write_h5ad(path)
    with h5py.File(path, 'r+') as handle:
        handle['layers/counts'] = handle['X']
    disk = _read_differential_h5ad(path, disk_backed=True)
    try:
        bind_row_store(disk)
        assert disk.X._analysis_row_store is disk.layers['counts']._analysis_row_store
        selected = np.array([31, 6, 102, 1, 7])
        np.testing.assert_array_equal(read_rows(disk.X, selected, slice(19, 67)), x[selected, 19:67])
        np.testing.assert_array_equal(read_rows(disk.X, slice(None)), x)
    finally:
        disk._analysis_h5_handle.close()


def test_unfiltered_export_keeps_values_metadata_and_shared_counts(tmp_path):
    import h5py
    x = np.random.default_rng(96).uniform(0, 3, (91, 701)).astype(np.float32)
    obj = ad.AnnData(X=x, obs=pd.DataFrame(index=[f'c{i}' for i in range(91)]),
                     var=pd.DataFrame(index=[f'g{i}' for i in range(701)]), uns={'note': 'preserved'},
                     obsm={'X_umap': np.ones((91, 2), dtype=np.float32)})
    obj.layers['counts'] = obj.X
    path = tmp_path / 'unfiltered.h5ad'
    D.write_differentials_only_h5ad(obj, {'detailed_deg': pd.DataFrame(columns=['gene'])}, path)
    actual = ad.read_h5ad(path)
    np.testing.assert_array_equal(actual.X, x)
    np.testing.assert_array_equal(actual.layers['counts'], x)
    np.testing.assert_array_equal(actual.obsm['X_umap'], obj.obsm['X_umap'])
    assert actual.uns['note'] == 'preserved'
    assert 'cellHarmony_DE' not in obj.uns
    with h5py.File(path) as handle:
        assert handle['X'].id == handle['layers/counts'].id
