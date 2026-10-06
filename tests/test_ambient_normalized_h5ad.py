"""Normalized RNA uses its actual counts for standard ambient correction."""
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from altanalyze3.components.ambient_rna.ambient_subtract import process_anndata
from altanalyze3.components.cellHarmony.cellHarmony_lite import combine_and_align_h5, normalize_adata
from altanalyze3.components.cellHarmony.markerFinder import detect_input_scaling
from altanalyze3.components.cellHarmony.mapped_h5ad import Workspace


@pytest.mark.parametrize('bounded', [False, True])
@pytest.mark.parametrize('rho', [0.2, 'auto'])
def test_normalized_counts_routing_matches_standard_raw_count_workflow(tmp_path, bounded, rho):
    rng = np.random.default_rng(932)
    counts = sparse.csr_matrix(rng.poisson(3, (140, 70)).astype(np.float32))
    obs = pd.DataFrame({'Library': ['a'] * 70 + ['b'] * 70}, index=[f'cell{i}' for i in range(140)])
    var = pd.DataFrame(index=[f'gene{i}' for i in range(70)])
    expected = process_anndata(ad.AnnData(counts.copy(), obs=obs.copy(), var=var.copy()),
                              rho=rho, library_col='Library', outdir=tmp_path / 'expected',
                              inplace=True, write_individual=False, write_merged=False,
                              store_corrected_layer=False)
    expected_counts = expected.X.copy()
    normalize_adata(expected)
    source = ad.AnnData(counts.copy(), obs=obs.copy(), var=var.copy())
    source.layers['counts'] = counts.copy()
    normalize_adata(source)
    source.raw = source.copy()
    path = tmp_path / 'normalized.h5ad'
    source.write_h5ad(path)
    _, actual = combine_and_align_h5(h5_files=[], h5ad_file=str(path), cellharmony_ref=None,
                                    output_dir=str(tmp_path / 'actual'), return_adata=True,
                                    ambient_correct_cutoff=str(rho), ambient_memory_efficient=True,
                                    bounded_h5ad=bounded, min_genes=99999, min_counts=999999,
                                    min_cells=99999, mit_percent=0)
    assert actual.obs_names.equals(source.obs_names)
    assert actual.var_names.equals(source.var_names)
    pd.testing.assert_frame_equal(actual.obs, source.obs)
    np.testing.assert_allclose(actual.X.toarray(), expected.X.toarray(), rtol=1e-6, atol=1e-6)
    np.testing.assert_array_equal(actual.layers['counts'].toarray(), expected_counts.toarray())
    np.testing.assert_array_equal(actual.layers['soupx_raw'].toarray(), counts.toarray())
    np.testing.assert_array_equal(actual.raw.X.toarray(), source.raw.X.toarray())
    assert detect_input_scaling(actual.X)['status'] == 'ok'
    assert actual.uns['soupx_correction']['input_source'] == 'layers/counts'


def test_logged_rna_without_raw_counts_fails_before_subtraction(tmp_path, monkeypatch):
    source = ad.AnnData(sparse.csr_matrix(np.full((10, 5), 3.5)))
    path = tmp_path / 'no_counts.h5ad'
    source.write_h5ad(path)
    def unexpected(*args, **kwargs):
        raise AssertionError('Logged RNA must never reach count subtraction')
    monkeypatch.setattr('altanalyze3.components.ambient_rna.ambient_subtract.process_anndata', unexpected)
    with pytest.raises(ValueError, match='raw counts'):
        combine_and_align_h5(h5_files=[], h5ad_file=str(path), cellharmony_ref=None,
                             output_dir=str(tmp_path / 'actual'), ambient_correct_cutoff='auto')


def test_scale_probe_preserves_exact_classifier_sample(tmp_path):
    rng = np.random.default_rng(932)
    source = ad.AnnData(sparse.random(23519, 90, density=.08, random_state=rng, format='csr'))
    source.layers['counts'] = (source.X * 19).astype(np.int64)
    path = tmp_path / 'source.h5ad'
    source.write_h5ad(path)
    workspace = Workspace(tmp_path / 'work')
    mapped = workspace.load(path)
    actual = workspace.scale_probe(mapped, random_state=14)
    rows = np.sort(np.random.default_rng(14).choice(source.n_obs, 20000, replace=False))
    for expected, observed in [(source.X[rows], actual.X),
                               (source.layers['counts'][rows], actual.layers['counts'])]:
        for a, b in [(expected.data, observed.data), (expected.indices, observed.indices),
                     (expected.indptr, observed.indptr)]:
            np.testing.assert_array_equal(a, b)
