"""Retention, live-request safety and real serving parity under cache eviction."""
import gc
from types import SimpleNamespace
import weakref
from concurrent.futures import ThreadPoolExecutor

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from fastapi.testclient import TestClient

from altanalyze3.components.cellHarmony.webapp.memory_cache import BoundedCache, CacheBudget, retained_bytes
from altanalyze3.components.cellHarmony.webapp.app import create_app, _read_cached_h5ad, _get_cache_lock
from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
from altanalyze3.components.cellHarmony.flask.tasks import JobRunner


def test_shared_budget_evicts_lru_and_releases_arrays():
    budget = CacheBudget(max_bytes=20000)
    a, b = BoundedCache(budget), BoundedCache(budget)
    value = np.ones(1000)
    ref = weakref.ref(value)
    a['first'] = value
    del value
    b['second'] = np.ones(1000)
    a.get('first')  # second is now the least recently used
    b['third'] = np.ones(1000)
    assert 'second' not in b and 'first' in a
    b['fourth'] = np.ones(1000)
    gc.collect()
    assert ref() is None and 'first' not in a
    assert budget.snapshot()['bytes'] <= budget.max_bytes


def test_expiration_oversize_and_live_request_safety():
    now = [0.0]
    budget = CacheBudget(max_bytes=10000, ttl_seconds=10, clock=lambda: now[0])
    cache = BoundedCache(budget)
    live = np.arange(1000)
    cache['live'] = live
    now[0] = 11
    assert cache.get('live') is None
    np.testing.assert_array_equal(live, np.arange(1000))
    cache['too_large'] = np.ones(2000)
    assert 'too_large' not in cache and budget.snapshot()['bytes'] == 0


def test_layers_raw_sparse_and_mmap_accounting(tmp_path):
    x = sp.eye(100, format='csr', dtype=np.float32)
    obj = ad.AnnData(x)
    before = retained_bytes(obj)
    obj.layers['counts'] = x.copy()
    obj.obsm['wide'] = np.ones((100, 1000), dtype=np.float32)
    obj.raw = obj.copy()
    assert retained_bytes(obj) >= before + 400000
    mapped = np.lib.format.open_memmap(tmp_path / 'matrix.npy', mode='w+', dtype=np.float32, shape=(100, 1000))
    assert retained_bytes(np.asarray(mapped)[:, :3]) == 0
    anonymous = np.ones((100, 1000), dtype=np.float32)
    assert retained_bytes(anonymous[:, :3]) == anonymous.nbytes


def test_concurrent_retention_stays_bounded():
    budget = CacheBudget(max_bytes=50000, max_entries=4)
    cache = BoundedCache(budget)
    def access(i):
        cache[str(i)] = np.ones(1000)
        cache.get(str(i))
        assert budget.snapshot()['bytes'] <= budget.max_bytes
    with ThreadPoolExecutor(max_workers=8) as ex:
        list(ex.map(access, range(100)))
    assert len(cache) <= 4
    cache.clear()
    assert budget.snapshot()['bytes'] == 0


def test_serving_reader_keeps_values_counts_and_coordinates(tmp_path):
    obj = ad.AnnData(sp.eye(10, format='csr', dtype=np.float32),
                     obs=pd.DataFrame({'state': ['A'] * 5 + ['B'] * 5}, index=[f'c{i}' for i in range(10)]))
    obj.layers['counts'] = obj.X * 3
    obj.layers['soupx_raw'] = obj.X * 4
    obj.raw = obj.copy()
    obj.obsm['X_adt'] = np.arange(1000, dtype=np.float32).reshape(10, 100)
    obj.uns['lineage_order'] = ['B', 'A']
    path = tmp_path / 'data.h5ad'
    obj.write_h5ad(path)
    got = _read_cached_h5ad(path)
    np.testing.assert_array_equal(got.X.toarray(), obj.X.toarray())
    np.testing.assert_array_equal(got.layers['counts'].toarray(), obj.layers['counts'].toarray())
    np.testing.assert_array_equal(got.obsm['X_adt'], obj.obsm['X_adt'][:, :2])
    assert got.obs.equals(obj.obs) and got.uns['lineage_order'].tolist() == ['B', 'A']
    assert got.raw is None and list(got.layers) == ['counts']


def test_serving_estimate_accepts_null_metadata(tmp_path):
    import h5py
    from altanalyze3.components.cellHarmony.webapp.app import _h5ad_memory_estimate
    obj = ad.AnnData(np.ones((4, 3), dtype=np.float32))
    obj.uns['optional'] = None
    path = tmp_path / 'null_metadata.h5ad'
    obj.write_h5ad(path)
    with h5py.File(path, 'r') as fh:
        assert fh['uns/optional'].size is None
    assert _h5ad_memory_estimate(path) >= obj.X.nbytes
    got = _read_cached_h5ad(path, CacheBudget())
    np.testing.assert_array_equal(got.X, obj.X)


def test_multiple_jobs_reload_after_eviction_without_changing_results(tmp_path):
    app = create_app(dict(JOB_STORAGE=str(tmp_path / 'jobs'), CACHE_MAX_GIB=0.00002))
    obj = ad.AnnData(sp.eye(30, format='csr', dtype=np.float32),
                     obs=pd.DataFrame({'state': ['A'] * 15 + ['B'] * 15}, index=[f'c{i}' for i in range(30)]),
                     var=pd.DataFrame(index=[f'G{i}' for i in range(30)]))
    obj.obsm['X_umap'] = np.zeros((30, 2), dtype=np.float32)
    store = app.state.job_store
    ids = []
    for i in range(5):
        job = store.create_job('human', 'test', None, [])['job_id']
        path = store.outputs_dir(job) / 'combined.h5ad'
        obj.write_h5ad(path)
        store.update_job(job, cluster_key='state', status='completed', artifacts={'combined_h5ad': str(path)})
        ids.append(job)
    try:
        with TestClient(app) as client:
            first = client.get(f'/api/jobs/{ids[0]}/expression', params={'gene': 'G0'}).json()
            for job in ids[1:]:
                response = client.get(f'/api/jobs/{job}/expression', params={'gene': 'G0'})
                assert response.status_code == 200
                assert app.state.result_cache_budget.snapshot()['bytes'] <= 0.00002 * 1024**3
            again = client.get(f'/api/jobs/{ids[0]}/expression', params={'gene': 'G0'}).json()
            assert first == again
    finally:
        app.state.job_runner.executor.shutdown(wait=True)


def test_cache_locks_survive_waiters_but_not_finished_requests():
    app = SimpleNamespace(state=SimpleNamespace(expression_cache_locks=weakref.WeakValueDictionary()))
    first = _get_cache_lock(app, 'expression_cache_locks', 'job:rna')
    assert _get_cache_lock(app, 'expression_cache_locks', 'job:rna') is first
    ref = weakref.ref(first)
    del first
    gc.collect()
    assert ref() is None


@pytest.mark.parametrize('task', ['pipeline', 'differential'])
def test_child_worker_kill_is_terminal_failure(tmp_path, monkeypatch, task):
    store = JobStore(tmp_path / 'jobs')
    job = store.create_job('human', 'test', None, [])['job_id']
    store.update_job(job, differential={'status': 'queued'})
    runner = JobRunner(store, tmp_path / 'registry.json', isolate_jobs=True)
    monkeypatch.setattr(runner, '_memory_wait_reason', lambda: None)
    monkeypatch.setattr('altanalyze3.components.cellHarmony.flask.tasks.subprocess.Popen',
                        lambda *a, **k: SimpleNamespace(pid=123456, wait=lambda timeout=None: -9))
    try:
        if task == 'pipeline':
            runner.submit_pipeline(job)
        else:
            runner.submit_differential(job)
        runner._futures[f'{task}:{job}'].result(timeout=5)
        meta = store.get_job(job)
        state = meta if task == 'pipeline' else meta['differential']
        assert state['status'] == 'failed' and state['worker_pid'] is None
        assert 'killed by the operating system' in state['message']
        assert 'uploaded files remain available' in state['message']
    finally:
        runner.executor.shutdown(wait=True)


@pytest.mark.parametrize('modality', ['metabolite', 'lipid', 'grn', 'grn_tf'])
def test_cell_comparison_reads_cells_from_legacy_metadata(tmp_path, modality):
    from altanalyze3.components.cellHarmony.flask.pipeline import _modality_differential_h5ad_path
    cells, pseudobulk = tmp_path / 'cells.h5ad', tmp_path / 'pb.h5ad'
    cells.touch()
    pseudobulk.touch()
    meta = {'modality_artifacts': {modality: {'h5ad': str(cells), 'differential_h5ad': str(pseudobulk),
                                             'pseudobulk_h5ad': str(pseudobulk)}}}
    if modality == 'grn':
        # Without this key the established legacy adapter interprets grn.h5ad
        # as enrichment and differential_h5ad as the edge matrix.
        meta['modality_artifacts']['grn_tf'] = {'h5ad': str(cells)}
    assert _modality_differential_h5ad_path(meta, modality, 'cells') == cells
    assert _modality_differential_h5ad_path(meta, modality, 'pseudobulk') == pseudobulk


def test_differential_reader_preserves_raw_and_counts_without_wide_copies(tmp_path):
    from altanalyze3.components.cellHarmony.flask.pipeline import _read_differential_h5ad
    obj = ad.AnnData(sp.eye(10, format='csr', dtype=np.float32))
    obj.layers['counts'] = obj.X * 3
    obj.obsm['X_adt'] = np.ones((10, 100), dtype=np.float32)
    obj.raw = obj.copy()
    path = tmp_path / 'data.h5ad'
    obj.write_h5ad(path)
    got = _read_differential_h5ad(path)
    np.testing.assert_array_equal(got.X.toarray(), obj.X.toarray())
    np.testing.assert_array_equal(got.raw.X.toarray(), obj.raw.X.toarray())
    np.testing.assert_array_equal(got.layers['counts'].toarray(), obj.layers['counts'].toarray())
    assert list(got.obsm) == []


def test_many_sample_names_have_short_distinct_output_tags():
    from altanalyze3.components.cellHarmony.flask.pipeline import _comparison_tag
    names = [f'replica{i}_a_long_sample_name' for i in range(30)]
    tag = _comparison_tag(names[:15], names[15:])
    assert len(tag.encode()) <= 120
    assert tag == _comparison_tag(names[:15], names[15:])
    assert tag != _comparison_tag(names[:15], names[15:-1] + ['different'])
    assert _comparison_tag(['A'], ['B']) == 'A_vs_B'


@pytest.mark.parametrize('dense', [False, True])
def test_bundle_counts_sum_matches_in_memory_across_blocks(tmp_path, dense):
    from altanalyze3.components.cellHarmony.webapp.integration_data import _sample_counts
    rng = np.random.default_rng(9)
    matrix = rng.random((1200, 11), dtype=np.float32)
    obj = ad.AnnData(matrix if dense else sp.csr_matrix(matrix))
    obj.layers['counts'] = obj.X.copy()
    path = tmp_path / 'counts.h5ad'
    obj.write_h5ad(path)
    proxy = SimpleNamespace(n_vars=11)  # bundle view exposes no in-memory layers
    rows = np.arange(0, 1200, 2)
    np.testing.assert_array_equal(_sample_counts(proxy, rows, path), _sample_counts(obj, rows))


@pytest.mark.parametrize('artifact_key', ['DEG_summary_A_vs_B', 'fold_matrix_tsv'])
def test_empty_completed_differential_has_empty_detail_table(tmp_path, artifact_key):
    from altanalyze3.components.cellHarmony.webapp.app import _get_differential_detail_table
    app = create_app(dict(JOB_STORAGE=str(tmp_path / 'jobs')))
    try:
        meta = {'job_id':'empty', 'differential': {'status':'completed',
                'artifacts': {artifact_key: str(tmp_path / 'summary.tsv')}}}
        frame = _get_differential_detail_table(app, meta)
        assert frame.empty and {'gene','population','fdr','log2fc'} <= set(frame.columns)
        meta['differential']['artifacts'] = {'DEG_detailed_A_vs_B': str(tmp_path / 'missing.tsv')}
        with pytest.raises(FileNotFoundError):
            _get_differential_detail_table(app, meta)
    finally:
        app.state.job_runner.executor.shutdown(wait=True)


def test_bundle_default_includes_medium_jobs(monkeypatch):
    from altanalyze3.components.cellHarmony.flask.pipeline import _bundle_min_cells
    monkeypatch.delenv('CELLHARMONY_BUNDLE_MIN_CELLS', raising=False)
    assert _bundle_min_cells() == (10000, '')


def test_metadata_options_read_annotations_without_loading_matrices(tmp_path, monkeypatch):
    from altanalyze3.components.cellHarmony.flask import pipeline
    obs = pd.DataFrame({'state': ['A'] * 6 + ['B'] * 6,
                        'Library': ['one', 'two'] * 6, 'pct_counts_mt': np.zeros(12)},
                       index=[f'c{i}' for i in range(12)])
    obj = ad.AnnData(sp.eye(12, format='csr', dtype=np.float32), obs=obs)
    obj.layers['counts'] = obj.X.copy()
    obj.obsm['X_adt'] = np.ones((12, 1000), dtype=np.float32)
    path = tmp_path / 'data.h5ad'
    obj.write_h5ad(path)
    monkeypatch.setattr(pipeline.ad, 'read_h5ad', lambda *a, **k: pytest.fail('loaded the matrix'))
    assert pipeline._read_h5ad_obs(path).equals(obj.obs)
    assert pipeline._candidate_population_columns(path, ['state']) == [dict(value='state', label='state', n_categories=2)]
    fields, values = pipeline._candidate_group_fields(path, ['Library'])
    assert fields[0]['value'] == 'Library' and values['Library'] == ['one', 'two']
