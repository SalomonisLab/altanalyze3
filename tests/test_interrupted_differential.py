from concurrent.futures import Future
from importlib import import_module
import os

import pytest
from fastapi.testclient import TestClient
from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
from altanalyze3.components.cellHarmony.flask.tasks import JobRunner


@pytest.fixture
def saved_job(tmp_path):
    store = JobStore(tmp_path / 'jobs')
    job = store.create_job('human', 'reference', None, [])['job_id']
    config = dict(modality='rna', population_col='state', sample_field='Library',
                  comparison_type='cells', group1_samples=['case'], group2_samples=['control'])
    store.update_job(job, status='completed', differential={'status': 'processing', 'progress': 15,
                     'config': config, 'artifacts': {'partial': 'keep.tsv'}},
                     differential_history={'prior': {'status': 'completed'}})
    runner = JobRunner(store, tmp_path / 'registry.json', max_workers=1)
    yield store, runner, job, config
    runner.executor.shutdown(wait=True)


def test_orphaned_processing_is_retryable_without_losing_results(saved_job):
    store, runner, job, config = saved_job
    result = runner.recover_interrupted_differential(job)
    assert result['status'] == 'completed'
    assert result['differential']['status'] == 'failed'
    assert result['differential']['config'] == config
    assert result['differential']['artifacts'] == {'partial': 'keep.tsv'}
    assert result['differential_history']['prior']['status'] == 'completed'
    assert runner.recover_interrupted_differential(job) == result
    assert (store.logs_dir(job) / 'pipeline.log').read_text().count('was interrupted') == 1


def test_active_worker_is_not_interrupted(saved_job):
    store, runner, job, _ = saved_job
    pending = Future()
    runner._futures[f'differential:{job}'] = pending
    assert runner.recover_interrupted_differential(job)['differential']['status'] == 'processing'
    pending.set_result(None)
    assert runner.recover_interrupted_differential(job)['differential']['status'] == 'failed'


def test_other_live_process_is_not_interrupted(saved_job, monkeypatch):
    store, runner, job, _ = saved_job
    meta = store.get_job(job)
    store.update_job(job, differential={**meta['differential'], 'worker_pid': os.getpid() + 10000})
    monkeypatch.setattr(os, 'kill', lambda pid, signal: None)
    assert runner.recover_interrupted_differential(job)['differential']['status'] == 'processing'


def test_dead_process_is_retryable(saved_job, monkeypatch):
    store, runner, job, _ = saved_job
    meta = store.get_job(job)
    store.update_job(job, differential={**meta['differential'], 'worker_pid': os.getpid() + 10000})
    def missing(pid, signal):
        raise ProcessLookupError()
    monkeypatch.setattr(os, 'kill', missing)
    assert runner.recover_interrupted_differential(job)['differential']['status'] == 'failed'


@pytest.mark.parametrize('read_status_first', [True, False])
def test_web_retry_recovers_stale_flag_and_blocks_active_run(saved_job, monkeypatch, read_status_first):
    store, runner, job, config = saved_job
    web = import_module('altanalyze3.components.cellHarmony.webapp.app')
    app = web.create_app({'JOB_STORAGE': str(store.root), 'JOB_WORKERS': 1})
    app.state.job_runner.executor.shutdown(wait=True)
    app.state.job_store, app.state.job_runner = store, runner
    monkeypatch.setattr(web, '_build_differential_payload', lambda app, job_id, meta, **kw: meta['differential'])
    monkeypatch.setattr(web, '_validate_differential_request', lambda meta, payload: payload.model_dump())
    monkeypatch.setattr(web, '_differential_options', lambda meta: {'default_population_col': 'state'})
    submissions = []
    def submit(job_id):
        submissions.append(job_id)
        runner._futures[f'differential:{job_id}'] = Future()
    monkeypatch.setattr(runner, 'submit_differential', submit)
    with TestClient(app) as client:
        if read_status_first:
            response = client.get(f'/api/jobs/{job}/differential/status')
            assert response.json()['status'] == 'failed'
        response = client.post(f'/api/jobs/{job}/differential', json=config)
        assert response.status_code == 200, response.text
        assert response.json()['status'] == 'queued'
        assert submissions == [job]
        response = client.post(f'/api/jobs/{job}/differential', json=config)
        assert response.status_code == 409
        assert submissions == [job]


def test_differential_working_copy_retains_test_inputs():
    import anndata as ad
    import numpy as np
    import pandas as pd
    from altanalyze3.components.cellHarmony.flask.pipeline import _trim_differential_working_data
    a = ad.AnnData(np.array([[1., 2.], [3., 4.]]),
                   obs=pd.DataFrame({'state': ['A', 'B']}, index=['c1', 'c2']),
                   var=pd.DataFrame(index=['g1', 'g2']))
    a.layers['counts'] = a.X * 10
    a.layers['soupx_raw'] = a.X * 11
    a.obsm['X_umap'] = np.array([[0., 0.], [1., 1.]])
    a.obsm['X_lipid'] = np.ones((2, 100))
    a.uns['imputed_modalities'] = {'lipid': {'obsm_key': 'X_lipid'}}
    a.uns['lineage_order'] = ['A', 'B']
    working = _trim_differential_working_data(a.copy())
    np.testing.assert_array_equal(working.X, a.X)
    np.testing.assert_array_equal(working.layers['counts'], a.layers['counts'])
    np.testing.assert_array_equal(working.obsm['X_umap'], a.obsm['X_umap'])
    assert working.obs.equals(a.obs) and working.var.equals(a.var)
    assert working.uns['lineage_order'] == ['A', 'B']
    assert set(working.layers) == {'counts'}
    assert set(working.obsm) == {'X_umap'}
    assert 'imputed_modalities' not in working.uns
    assert 'soupx_raw' in a.layers and 'X_lipid' in a.obsm
