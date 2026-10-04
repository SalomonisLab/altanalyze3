from altanalyze3.components.cellHarmony.webapp.log_snapshot import read_pipeline_log, _CACHE


def test_unchanged_log_is_cached_and_append_refreshes(tmp_path):
    path=tmp_path/'pipeline.log'
    lines=[f'message {i}\n' for i in range(500)]
    path.write_text(''.join(lines))
    head,tail,progress=read_pipeline_log(path)
    assert head==lines[:80] and tail==lines[-200:] and not progress
    assert read_pipeline_log(path)[0] is head
    with path.open('a') as fh:
        fh.write('Cells remaining after min_genes 500 filtering: 350\n')
    assert read_pipeline_log(path)[2]==['Cells remaining after min_genes 500 filtering: 350\n']
    assert _CACHE.budget.snapshot()['bytes']<=1024**2


def test_qc_survives_later_logs_and_uses_latest_counts(tmp_path):
    path=tmp_path/'pipeline.log'
    relevant=['adata shape: (100, 40)\n','Cells remaining after min_genes 500 filtering: 90\n',
              'Cells remaining after min_genes 500 filtering: 85\n',
              "Auto-selected rho for library 'A:1': 0.1\n","Auto-selected rho for library 'A:2': 0.2\n",
              '[INFO] Applied min_alignment_score=0.4. Excluded 5 cells, kept 80.\n']
    path.write_text(''.join(relevant)+('later work\n'*500))
    _,tail,progress=read_pipeline_log(path)
    assert relevant[-1] in progress and relevant[0] in progress
    assert relevant[1] not in progress and relevant[2] in progress
    assert relevant[3] in progress and relevant[4] in progress
    assert relevant[-1] not in tail


def test_replaced_truncated_deleted_and_partial_logs(tmp_path):
    path=tmp_path/'pipeline.log'
    path.write_text('old\n'*100)
    read_pipeline_log(path)
    path.write_text('new partial')
    assert read_pipeline_log(path)[:2]==(['new partial'],['new partial'])
    other=tmp_path/'replacement'
    other.write_text('replacement\n')
    other.replace(path)
    assert read_pipeline_log(path)[0]==['replacement\n']
    path.unlink()
    assert read_pipeline_log(path)==([],[],[])


def test_main_duration_ignores_later_differentials_and_reruns(tmp_path):
    from altanalyze3.components.cellHarmony.webapp.log_snapshot import analysis_duration_seconds

    path = tmp_path / 'pipeline.log'
    path.write_text('[2026-10-03T01:18:22.054320] Job accepted by worker.\n'
                    '[2026-10-03T01:19:28.925823] Job completed.\n'
                    + 'later work\n' * 600
                    + '[2026-10-03T01:32:19.266739] Differential analysis completed.\n')
    meta = {'status': 'completed', 'created_at': '2026-10-03T01:18:15.332648Z',
            'updated_at': '2026-10-03T01:32:19.262797Z'}
    progress = read_pipeline_log(path)[2]
    assert round(analysis_duration_seconds(meta, progress)) == 67
    assert analysis_duration_seconds(dict(meta, analysis_duration_seconds=12.5), progress) == 12.5
    with path.open('a') as stream:
        stream.write('[2026-10-03T02:00:00] Job accepted by worker.\n')
    assert analysis_duration_seconds(meta, read_pipeline_log(path)[2]) is None
    with path.open('a') as stream:
        stream.write('[2026-10-03T02:00:20] Job completed.\n')
    assert analysis_duration_seconds(meta, read_pipeline_log(path)[2]) == 20
    assert analysis_duration_seconds({'status': 'processing', 'analysis_duration_seconds': 67}, progress) is None
    assert analysis_duration_seconds(meta, []) is None


def test_worker_duration_survives_differential_metadata_updates(tmp_path, monkeypatch):
    from altanalyze3.components.cellHarmony.flask import tasks
    from altanalyze3.components.cellHarmony.flask.job_manager import JobStore

    store = JobStore(tmp_path / 'jobs')
    job = store.create_job('human', 'test', None, files=[])['job_id']
    runner = tasks.JobRunner(store, tmp_path / 'registry.json', max_workers=1, isolate_jobs=False)
    calls = []
    def fake_pipeline(job_id, *args, **kwargs):
        state = store.get_job(job_id)
        assert state['analysis_duration_seconds'] is None
        assert state['analysis_completed_at'] is None
        assert state['analysis_started_at']
        calls.append(job_id)
    monkeypatch.setattr(tasks, 'run_cellharmony_pipeline', fake_pipeline)
    try:
        runner._run_pipeline(job)
        first = store.get_job(job)
        assert first['status'] == 'completed' and first['analysis_duration_seconds'] >= 0
        store.update_job(job, differential={'status': 'completed'})
        assert store.get_job(job)['analysis_duration_seconds'] == first['analysis_duration_seconds']
        runner._run_pipeline(job)
        assert len(calls) == 2
    finally:
        runner.executor.shutdown(wait=True)
