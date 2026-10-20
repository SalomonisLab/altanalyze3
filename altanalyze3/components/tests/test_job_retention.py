from datetime import datetime, timedelta, timezone

import pytest

from altanalyze3.components.cellHarmony.flask.job_manager import JobStore


def old_job(store, hours, **changes):
    job = store.create_job('mouse', 'mm_bm_reference', None, [])['job_id']
    meta = store.get_job(job)
    meta.update(status='completed', updated_at=(datetime.now(timezone.utc)-timedelta(hours=hours)).isoformat(), **changes)
    store._write_metadata(job,meta)
    return job


def test_default_retains_yesterdays_results_but_expires_after_three_days(tmp_path, monkeypatch):
    monkeypatch.delenv('CELLHARMONY_JOB_RETENTION_HOURS', raising=False)
    store = JobStore(tmp_path)
    yesterday = old_job(store,24)
    expired = old_job(store,73)
    store.create_job('human','ref',None,[])
    assert store.job_exists(yesterday)
    assert not store.job_exists(expired)


def test_user_visit_renews_retention_without_changing_analysis_metadata(tmp_path):
    store = JobStore(tmp_path)
    job = old_job(store,74)
    before = store.get_job(job)
    store.touch_job(job)
    assert store.get_job(job) == before
    assert store.purge_old_jobs() == 0
    assert store.job_exists(job)


def test_configured_retention_preserves_active_pipeline_and_differential(tmp_path, monkeypatch):
    monkeypatch.setenv('CELLHARMONY_JOB_RETENTION_HOURS','8')
    store = JobStore(tmp_path)
    expired = old_job(store,9)
    store.update_job(expired,status='processing')
    # Restore old timestamp to exercise status protection rather than freshness.
    meta=store.get_job(expired);meta['updated_at']=(datetime.now(timezone.utc)-timedelta(hours=100)).isoformat();store._write_metadata(expired,meta)
    de = old_job(store,9,differential={'status':'processing'})
    assert store.purge_old_jobs() == 0
    assert store.job_exists(expired) and store.job_exists(de)


@pytest.mark.parametrize('value',['0','-1','nan','inf'])
def test_invalid_retention_rejected(tmp_path, monkeypatch,value):
    monkeypatch.setenv('CELLHARMONY_JOB_RETENTION_HOURS',value)
    with pytest.raises(ValueError):JobStore(tmp_path)


def test_status_visit_renews_real_web_job(tmp_path):
    from fastapi.testclient import TestClient
    from altanalyze3.components.cellHarmony.webapp.app import create_app
    app = create_app({'JOB_STORAGE': str(tmp_path/'jobs'),
                      'REFERENCE_REGISTRY': str(tmp_path/'registry.json')})
    store = app.state.job_store
    job = old_job(store,24)
    before = store.get_job(job)['updated_at']
    try:
        with TestClient(app) as client:
            assert client.get(f'/api/jobs/{job}/status').status_code == 200
        assert (store.root/job/'.last_access').exists()
        assert store.get_job(job)['updated_at'] == before
    finally:
        app.state.job_runner.executor.shutdown(wait=True)
