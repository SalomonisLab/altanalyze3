import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from altanalyze3.components.clustering.umap_input import dense_umap_input, sparse_umap_input


@pytest.mark.parametrize('kind', ['csr','csc','dense'])
def test_umap_consumes_identical_full_matrix(kind):
    rng = np.random.default_rng(1)
    x = rng.normal(size=(71,31))
    x[x < .1] = 0
    if kind != 'dense': x = getattr(sp, kind + '_matrix')(x)
    a = ad.AnnData(X=x, var=pd.DataFrame(index=[f'G{i}' for i in range(31)]))
    features = a.var_names[::-2]
    old = a[:,features].X
    expected = np.asarray(old.toarray() if sp.issparse(old) else old, dtype=np.float32, order='C')
    actual = dense_umap_input(a, features, block_rows=7)
    assert actual.dtype == np.float32 and actual.flags.c_contiguous
    assert actual.tobytes() == expected.tobytes()
    with pytest.raises(ValueError, match='identities'):
        dense_umap_input(a, ['unresolved'])


def test_reordered_cells_match_original_subset_and_input(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from altanalyze3.components.clustering import ICGS
    x = np.random.default_rng(8).normal(size=(31, 13))
    a = ad.AnnData(X=sp.csr_matrix(x), obs=pd.DataFrame(index=[f'C{i}' for i in range(31)]),
                   var=pd.DataFrame(index=[f'G{i}' for i in range(13)]))
    cells = a.obs_names[[8,2,24,0,12]]
    features = a.var_names[[3,1,7,8]]
    expected = a[cells,features].X.toarray().astype(np.float32)
    assert dense_umap_input(a,features,cells=cells,block_rows=2).tobytes() == expected.tobytes()
    seen = []
    class Model:
        def __init__(self, **kwargs): seen.append(kwargs)
        def fit_transform(self, matrix):
            assert matrix.tobytes() == expected.tobytes()
            return matrix[:,:2].copy()
    monkeypatch.setattr(ICGS, '_import_umap_with_local_retry', lambda:SimpleNamespace(UMAP=Model))
    config = ICGS.ICGS3Config(input_paths=[], output_dir=str(tmp_path), generate_umap=True, minimal_outputs=True)
    target = a[cells].copy()
    ICGS.compute_umap_outputs(target,config,str(tmp_path),marker_features=list(features),matrix_source=a)
    np.testing.assert_array_equal(target.obsm['X_umap'],expected[:,:2])
    assert target.obs_names.equals(cells) and target.var_names.equals(a.var_names)
    assert list(pd.read_csv(tmp_path/'UMAPs/icgs3_umap_features.tsv',sep='\t')['feature']) == list(features)
    assert target.uns['icgs3_umap']['feature_file'] is not None
    assert seen[0]['metric'] == 'correlation' and seen[0]['random_state'] == config.random_state
    with pytest.raises(ValueError,match='cell identities'):
        dense_umap_input(a,features,cells=['unresolved'])


@pytest.mark.parametrize('kind', ['csr', 'csc', 'dense', 'backed'])
def test_pca_input_preserves_all_ordered_features_cells_and_values(kind, tmp_path):
    x = np.random.default_rng(12).normal(size=(71, 13))
    x[x < .1] = 0
    matrix = x if kind == 'dense' else getattr(sp, ('csc' if kind == 'csc' else 'csr') + '_matrix')(x)
    a = ad.AnnData(matrix, obs=pd.DataFrame(index=[f'C{i}' for i in range(71)]),
                   var=pd.DataFrame(index=[f'G{i}' for i in range(13)]))
    if kind == 'backed':
        a.write_h5ad(tmp_path / 'input.h5ad')
        a = ad.read_h5ad(tmp_path / 'input.h5ad', backed='r')
    try:
        features, cells = a.var_names[::-2], a.obs_names[::-1]
        expected = x[::-1, ::-2].astype(np.float32)
        actual = sparse_umap_input(a, features, cells=cells, block_rows=7)
        np.testing.assert_array_equal(actual.toarray(), expected)
        assert sp.isspmatrix_csr(actual) and actual.dtype == np.float32
        np.testing.assert_array_equal(dense_umap_input(a, a.var_names), x.astype(np.float32))
        with pytest.raises(ValueError, match='feature identities'):
            sparse_umap_input(a, ['unresolved'])
        with pytest.raises(ValueError, match='cell identities'):
            sparse_umap_input(a, features, cells=['unresolved'])
    finally:
        if a.isbacked: a.file.close()


def test_accelerated_changes_only_embedding_and_records_exact_panel(tmp_path, monkeypatch):
    from altanalyze3.components.clustering import ICGS, accelerated_umap
    x = np.random.default_rng(17).normal(size=(71, 13)).astype(np.float32)
    a = ad.AnnData(sp.csr_matrix(x), obs=pd.DataFrame(index=[f'C{i}' for i in range(71)]),
                   var=pd.DataFrame(index=[f'G{i}' for i in range(13)]))
    a.obs['ICGS3_cluster'] = ['C1'] * 70 + ['C2']
    a.layers['counts'] = a.X.copy()
    target = a[a.obs_names[::-1]].copy()
    features = ['G8', 'G4', 'G1']
    def pca(matrix, **kwargs):
        np.testing.assert_array_equal(matrix.toarray(), x[::-1][:, [8, 4, 1]])
        assert kwargs['labels'][0] == 'C2'
        assert kwargs['method'] == 'landmark' and kwargs['n_neighbors'] == 15
        return np.zeros((71, 2), dtype=np.float32), {'fit_mode': 'full', 'total_cells': 71}, np.arange(71)
    monkeypatch.setattr(accelerated_umap, 'pca_umap', pca)
    config = ICGS.ICGS3Config(input_paths=[], output_dir=str(tmp_path), minimal_outputs=True,
                            umap_fit_mode='pca_landmark')
    ICGS.compute_umap_outputs(target, config, str(tmp_path), marker_features=features, matrix_source=a)
    np.testing.assert_array_equal(target.X.toarray(), x[::-1])
    np.testing.assert_array_equal(target.layers['counts'].toarray(), x[::-1])
    assert list(target.var_names) == list(a.var_names)
    assert list(target.obs['ICGS3_cluster']) == list(a.obs['ICGS3_cluster'])[::-1]
    assert target.uns['icgs3_umap']['requested_fit_mode'] == 'pca_landmark'
    assert list(pd.read_csv(tmp_path / 'UMAPs/icgs3_umap_features.tsv', sep='\t').feature) == features
    with pytest.raises(ValueError, match='feature identities'):
        ICGS.compute_umap_outputs(target, config, str(tmp_path), marker_features=['unresolved'], matrix_source=a)
