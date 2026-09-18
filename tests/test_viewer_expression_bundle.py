"""Explore must use published bundle arrays after source H5AD cleanup."""
from types import SimpleNamespace
import numpy as np
from fastapi.testclient import TestClient
from altanalyze3.components.visualization.scalable_viewer import bundle_meta as B


def test_bundle_explore_without_source_h5ad(tmp_path, monkeypatch):
    clusters = tmp_path / 'clusters.tsv'
    clusters.write_text('barcode\tcluster\nc1\tAT1\nc2\tAT2\n')
    ds = SimpleNamespace(
        id='bundle', paths=SimpleNamespace(bundle_dir=str(tmp_path), clusters=str(clusters)),
        n_cells=2, symbols=['GENE1'], states=['AT1', 'AT2'], state_code=np.array([0, 1]),
        cluster_key='cell_type', cells={}, covariate_names=lambda: {},
        embedding=np.array([[1., 2.], [3., 4.]]),
        resolve_gene=lambda gene: 0 if gene == 'GENE1' else None,
        gene_column=lambda row: (np.array([0, 1]), np.array([2., 5.])),
    )
    app = B.W.create_app({'JOB_STORAGE': str(tmp_path / 'jobs')})
    store = B.BundleJobStore(SimpleNamespace(get=lambda job: ds, ids=lambda: ['bundle']), tmp_path, {})
    app.state.job_store = store
    meta = dict(job_id='bundle', cluster_key='cell_type', status='completed',
                artifacts={'combined_h5ad': str(tmp_path / 'removed.h5ad')})
    store._meta['bundle'] = meta
    def forbidden(*args, **kwargs):
        raise AssertionError('A bundle must never open or validate the source H5AD')
    monkeypatch.setattr(B.W, '_modality_h5ad_path', forbidden)
    monkeypatch.setattr(B.W.ad, 'read_h5ad', forbidden)
    try:
        client = TestClient(app)
        for endpoint in ('status', 'differential/status'):
            response = client.get('/api/jobs/bundle/' + endpoint)
            assert response.status_code == 200, response.text
        for cold in (True, False):
            if cold:
                app.state.expression_cache.clear()
            variables = client.get('/api/jobs/bundle/plot-variables')
            assert variables.status_code == 200, variables.text
            assert variables.json()['cluster_key'] == 'cell_type'
            assert variables.json()['coords']
            filters = client.get('/api/jobs/bundle/display-filters')
            assert filters.status_code == 200, filters.text
            assert filters.json()['values']['cell_type'] == ['AT1', 'AT2']
            expression = client.get('/api/jobs/bundle/expression', params={'gene': 'GENE1'})
            assert expression.status_code == 200, expression.text
            assert [(r['x'], r['y'], r['value']) for r in expression.json()['umap']] == [(1., 2., 2.), (3., 4., 5.)]
            umap = client.get('/api/jobs/bundle/umap')
            assert umap.status_code == 200, umap.text
        B.W._invalidate_expression_cache(app, 'bundle')
        assert B.W._get_expression_cache(app, meta)['bundle_path'] == str(tmp_path)
    finally:
        app.state.job_runner.executor.shutdown(wait=True)
