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
