"""Serving optimizations preserve complete catalogs and exact baseline statistics."""
from importlib import import_module
import h5py
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from test_plot_performance_parity import cache_for
from altanalyze3.components.cellHarmony.webapp.feature_lookup import FeatureLookup

W = import_module('altanalyze3.components.cellHarmony.webapp.app')


@pytest.mark.parametrize("n_features", [4, 32000])
def test_catalog_reads_only_var_and_preserves_categorical_aliases(tmp_path, monkeypatch, n_features):
    app = W.create_app({'JOB_STORAGE': str(tmp_path / 'jobs')})
    source = tmp_path / 'source.h5ad'
    var = pd.DataFrame({'gene_symbols': pd.Categorical(['TP53', 'SHARED', 'SHARED', ''] + [f'SYM{i}' for i in range(4, n_features)])},
                       index=['ENSG000001.5', 'ENSG000002', 'ENSG000003', 'other'] + [f'ENSG{i:011}' for i in range(4,n_features)])
    ad.AnnData(sparse.csr_matrix((10000, n_features)), var=var).write_h5ad(source)
    # This intentionally metadata-only fixture has no X, obs or embeddings.
    # Reading an expression cache cannot work; the complete var catalog can.
    with h5py.File(source, 'a') as handle:
        del handle['X']; del handle['obs']
    meta = {'status': 'completed', 'artifacts': {'combined_h5ad': str(source)}}
    monkeypatch.setattr(W, '_get_expression_cache', lambda *a, **k: pytest.fail('matrix loaded for gene catalog'))
    lookup = FeatureLookup(var.index, var)
    result = W._build_gene_suggestions_payload(app, meta)
    assert result['genes'] == lookup.display_names
    assert result['keys'] == lookup.suggestions
    assert set(var.index) <= set(result['keys'])
    assert 'TP53' in result['keys'] and 'SHARED' not in result['keys']
    assert W._build_gene_suggestions_payload(app, meta) is result
    # A new source timestamp invalidates the cached response.
    with h5py.File(source, 'a') as handle:
        handle.attrs['source-revision'] = 2
    assert W._build_gene_suggestions_payload(app, meta) == result
    assert W._build_gene_suggestions_payload(app, meta) is not result


@pytest.mark.parametrize('limit', [0, 1, 10, 30])
@pytest.mark.parametrize('filter_values', [None, ['s1'], ['missing']])
def test_violin_matches_previous_all_state_mean_order_and_values(monkeypatch, limit, filter_values):
    rng = np.random.default_rng(41)
    x = rng.normal(size=(510, 1)).astype(np.float32)
    x[::9] = 0; x[::13] = np.nan; x[::41] = np.inf
    cache = cache_for(x)
    cache['populations'] = np.array([f'c{i % 33}' for i in range(len(x))])
    # Equal means across two states exercise the previous alphabetical tie order.
    for state in ('c2', 'c11'):
        x[cache['populations'] == state] = 0
    monkeypatch.setattr(W, '_get_expression_cache', lambda *a, **k: cache)
    filters = [('Library', filter_values)] if filter_values else None
    mask = W._apply_display_filter_mask(cache, filters)
    previous = []
    for state in sorted(pd.unique(cache['populations'])):
        values = x[:, 0][(cache['populations'] == state) & mask].astype(float)
        values = values[np.isfinite(values)]
        if len(values):
            previous.append(dict(population=state, values=values.tolist(), mean=float(np.mean(values))))
    previous.sort(key=lambda p: p['mean'], reverse=True)
    actual = W._build_expression_payload(None, {}, 'g0', view='violin', violin_limit=limit, display_filters=filters)
    assert actual['violin'] == previous[:limit]


@pytest.mark.parametrize('dtype', [np.float32, np.float64])
def test_dotplot_batch_preserves_individual_column_reduction(dtype):
    rng = np.random.default_rng(17)
    x = rng.normal(size=(451, 23)).astype(dtype)
    x[rng.random(x.shape) < .75] = 0
    rows = list(range(22, -1, -1))
    cache = cache_for(sparse.csr_matrix(x))
    actual = W._gene_state_stats(cache, [f'g{i}' for i in rows], group_by='Library', subset_by='state', subset_values=['A'])
    _, groups, labels = W._group_axis(cache, 'Library')
    masks = [(labels == group) & (cache['populations'] == 'A') for group in groups]
    assert actual['mean'] == [[float(x[:, row][mask].mean()) for mask in masks] for row in rows]
    assert actual['frac'] == [[float((x[:, row][mask] > 0).mean()) for mask in masks] for row in rows]


def test_morpheus_uses_single_dataset_get_and_visible_failure():
    html = W._build_marker_heatmap_viewer_html('test', '')
    assert 'method: "HEAD"' not in html
    assert html.count('fetch(datasetUrl') == 1
    assert 'marker-heatmap-dataset' in html
    assert 'marker-heatmap-error' in html
    assert 'const viewerUrl = new URL(document.baseURI);' in html
    assert '${viewerUrl.search}' in html
    assert 'event.origin !== viewerUrl.origin' in html
    assert 'render().catch(err => showFallback' in html
