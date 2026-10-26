"""Session lifecycle tests with actual isolated workers and cross-process leases."""
import io
import os
from pathlib import Path
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import anndata as ad
from fastapi.testclient import TestClient
import numpy as np
import pytest

from altanalyze3.components.cellHarmony.flask import tasks
from altanalyze3.components.cellHarmony.flask.job_control import (
    JobBusyError, WorkerWatchdog, control_job, execution_busy, execution_lease)
from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
from altanalyze3.components.cellHarmony.webapp.app import create_app


def wait_for(predicate):
    deadline = time.monotonic()+10
    while time.monotonic()<deadline:
        if predicate(): return
        time.sleep(.02)
    pytest.fail('Worker state did not change')


@pytest.fixture
def server(tmp_path,monkeypatch):
    app=create_app({'JOB_STORAGE':str(tmp_path/'jobs'),'ISOLATE_JOBS':True,'JOB_WORKERS':2})
    runner=app.state.job_runner
    monkeypatch.setattr(runner,'_memory_wait_reason',lambda:None)
    original=tasks.subprocess.Popen
    code='''
import os,sys,time
from pathlib import Path
from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
s=JobStore(Path(sys.argv[1])); j,t=sys.argv[2:4]
if t=='pipeline': s.update_job(j,status='processing',worker_pid=os.getpid(),progress=30)
else: s.update_job(j,differential={'status':'processing','worker_pid':os.getpid()})
s.append_log(j,'TEST worker entered QC')
while True: time.sleep(.1)
'''
    def launch(cmd,**kw):
        if cmd[1:3]==['-m',runner.WORKER_MODULE]:
            return original([sys.executable,'-c',code,cmd[3],cmd[5],cmd[6]],**kw)
        return original(cmd,**kw)
    monkeypatch.setattr(tasks.subprocess,'Popen',launch)
    with TestClient(app) as client:
        yield app,client
        for job in app.state.job_store.list_jobs():
            runner.cancel_job(job['job_id'])
    runner.executor.shutdown(wait=True)


def new_job(app):
    return app.state.job_store.create_job('human','fixture',None,[])['job_id']


def test_duplicate_run_and_mutation_guards_and_real_stop(server):
    app,client=server;store=app.state.job_store;runner=app.state.job_runner
    job=new_job(app)
    assert client.post(f'/api/jobs/{job}/run').status_code==200
    wait_for(lambda:store.get_job(job)['status']=='processing')
    pid=next(iter(runner._active_workers))
    for suffix,payload in [('run',None),('qc',{}),('configure',{'species':'human','reference':'fixture'}),
                           ('differential/select?contrast=test',None)]:
        r=client.post(f'/api/jobs/{job}/{suffix}',json=payload)
        assert r.status_code==409,r.text
    before=store.get_job(job)['qc']
    assert client.post(f'/api/jobs/{job}/cancel').status_code==200
    assert store.get_job(job)['status']=='cancelled'
    assert store.get_job(job)['qc']==before
    assert not execution_busy(store,job)
    with pytest.raises(ProcessLookupError):os.kill(pid,0)
    # A stopped session can be rerun with its uploads; logs do not repeat QC.
    assert client.post(f'/api/jobs/{job}/run').status_code==200
    wait_for(lambda:store.get_job(job)['status']=='processing')
    assert len(list((store.logs_dir(job)/'history').glob('*/pipeline.log')))==1
    assert (store.logs_dir(job)/'pipeline.log').read_text().count('Job queued by user request.')==1


def test_two_simultaneous_requests_start_only_one_worker(server):
    app,client=server;job=new_job(app)
    with ThreadPoolExecutor(2) as pool:
        responses=list(pool.map(lambda _:client.post(f'/api/jobs/{job}/run'),range(2)))
    assert sorted(r.status_code for r in responses)==[200,409]
    wait_for(lambda:app.state.job_store.get_job(job)['status']=='processing')
    assert len(app.state.job_runner._active_workers)==1


def test_cancel_memory_queued_job_does_not_launch(server,monkeypatch):
    app,client=server;runner=app.state.job_runner;job=new_job(app)
    monkeypatch.setattr(runner,'_memory_wait_reason',lambda:'Waiting for memory.')
    assert client.post(f'/api/jobs/{job}/run').status_code==200
    assert client.post(f'/api/jobs/{job}/cancel').status_code==200
    assert app.state.job_store.get_job(job)['status']=='cancelled'
    assert not runner._active_workers


def test_pipeline_blocked_while_differential_runs(server):
    app,client=server;runner=app.state.job_runner;store=app.state.job_store;job=new_job(app)
    store.update_job(job,status='completed')
    runner.submit_differential(job)
    wait_for(lambda:(store.get_job(job).get('differential') or {}).get('status')=='processing')
    assert client.post(f'/api/jobs/{job}/run').status_code==409
    assert client.post(f'/api/jobs/{job}/qc',json={}).status_code==409
    assert client.post(f'/api/jobs/{job}/cancel').status_code==200
    assert store.get_job(job)['status']=='completed'
    assert store.get_job(job)['differential']['status']=='cancelled'


def test_new_upload_stops_previous_run_and_always_has_fresh_id(server,tmp_path):
    app,client=server;store=app.state.job_store;job=new_job(app)
    client.post(f'/api/jobs/{job}/run')
    wait_for(lambda:store.get_job(job)['status']=='processing')
    path=tmp_path/'small.h5ad';ad.AnnData(np.ones((8,4),dtype=np.float32)).write_h5ad(path)
    data={'species':'human','reference':'fixture','sample_names':'new','previous_job_id':job}
    r=client.post('/api/jobs',data=data,files={'files':('small.h5ad',path.read_bytes(),'application/octet-stream')})
    assert r.status_code==200,r.text
    assert r.json()['job_id']!=job
    assert store.get_job(job)['status']=='cancelled'
    assert store.get_job(r.json()['job_id'])['status']=='uploaded'


def test_file_lease_blocks_independent_web_processes(tmp_path):
    store=JobStore(tmp_path);job=store.create_job('human','fixture',None,[])['job_id']
    lease=execution_lease(store,job)
    code='''
import sys
from pathlib import Path
from altanalyze3.components.cellHarmony.flask.job_control import execution_busy
from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
sys.exit(0 if execution_busy(JobStore(Path(sys.argv[1])),sys.argv[2]) else 1)
'''
    try: assert subprocess.run([sys.executable,'-c',code,str(tmp_path),job]).returncode==0
    finally:lease.close()
    assert not execution_busy(store,job)


def test_watchdog_does_not_treat_quiet_active_computation_as_stalled():
    now=[0.];w=WorkerWatchdog(idle_seconds=10,max_seconds=100,clock=lambda:now[0])
    for i in range(8):
        now[0]=i*9
        assert w.check((i,0),(1,100)) is None
    now[0]=73
    assert 'without worker CPU' in w.check((7,0),(1,100))


def test_watchdog_bounds_busy_hangs_and_never_guesses_missing_telemetry():
    now=[0.];w=WorkerWatchdog(10,100,clock=lambda:now[0])
    assert w.check(None,()) is None
    now[0]=50
    assert w.check(None,()) is None
    now[0]=101
    assert 'runtime limit' in w.check((50,0),())


def test_watchdog_releases_real_stalled_worker(server,monkeypatch):
    app,client=server;runner=app.state.job_runner;store=app.state.job_store;job=new_job(app)
    monkeypatch.setenv('CELLHARMONY_WORKER_IDLE_SECONDS','1')
    monkeypatch.setattr(tasks,'worker_activity',lambda *_:(0,0))
    client.post(f'/api/jobs/{job}/run')
    wait_for(lambda:store.get_job(job)['status']=='failed')
    assert 'without worker CPU' in store.get_job(job)['message']
    wait_for(lambda:not execution_busy(store,job))
    assert not runner._active_workers


def test_discover_rejects_mutating_active_sessions(tmp_path):
    from altanalyze3.components.cellHarmony.scalable_discover.app import create_discover_app
    app=create_discover_app({'JOB_STORAGE':str(tmp_path)})
    store=app.state.job_store;job=store.create_job('human','icgs3',None,[])['job_id']
    store.update_job(job,status='processing')
    try:
        with TestClient(app) as client:
            assert client.post(f'/api/jobs/{job}/qc',json={}).status_code==409
            assert client.post(f'/api/jobs/{job}/configure',json={'species':'human','reference':'icgs3'}).status_code==409
    finally:app.state.job_runner.executor.shutdown(wait=True)


def test_termination_includes_child_that_ignores_sigterm(tmp_path):
    ready=tmp_path/'child_pid'
    child="import os,signal,time;from pathlib import Path;signal.signal(signal.SIGTERM,signal.SIG_IGN);Path("+repr(str(ready))+").write_text(str(os.getpid()));time.sleep(100)"
    parent="import subprocess,sys,time;subprocess.Popen([sys.executable,'-c',"+repr(child)+"]);time.sleep(100)"
    process=subprocess.Popen([sys.executable,'-c',parent],start_new_session=True)
    try:
        wait_for(ready.exists)
        pid=int(ready.read_text())
        tasks.JobRunner._terminate_worker(process)
        assert process.poll() is not None
        def stopped():
            result=subprocess.run(['ps','-p',str(pid),'-o','stat='],capture_output=True,text=True)
            return result.returncode!=0 or result.stdout.strip().startswith('Z')
        wait_for(stopped)
    finally:
        try:os.killpg(process.pid,9)
        except ProcessLookupError:pass
        process.wait()


def test_cancellation_never_signals_the_web_process(tmp_path,monkeypatch):
    runner=tasks.JobRunner(JobStore(tmp_path),tmp_path/'registry.json')
    signals=[]
    monkeypatch.setattr(tasks.os,'killpg',lambda *args:signals.append(args))
    try:
        runner._stop_recorded_worker(os.getpid(),'test')
        assert not signals
    finally:runner.executor.shutdown(wait=True)


@pytest.mark.parametrize('value',['nan','inf','-1'])
def test_watchdog_rejects_invalid_limits(value):
    with pytest.raises(ValueError):WorkerWatchdog(value,100)


def test_fresh_visitor_cannot_cancel_another_running_job(server):
    app,client=server;store=app.state.job_store;job=new_job(app)
    client.post(f'/api/jobs/{job}/run')
    wait_for(lambda:store.get_job(job)['status']=='processing')
    assert client.get('/').status_code==200
    app.state.job_runner.cancel_disconnected_jobs()
    assert store.get_job(job)['status']=='processing'
    assert client.post('/api/jobs/not-a-real-job/browser-presence',json={
        'client_id':'visitor_new_page','state':'leaving'}).status_code==404
    assert store.get_job(job)['status']=='processing'


def test_all_tabs_must_leave_before_cancellation(server,monkeypatch):
    from altanalyze3.components.cellHarmony.flask.browser_presence import record_presence
    app,client=server;store=app.state.job_store;runner=app.state.job_runner;job=new_job(app)
    clients=['browser_client_first','browser_client_second']
    for key in clients:
        assert client.post(f'/api/jobs/{job}/browser-presence',json={
            'client_id':key,'state':'active','sequence':1}).status_code==200
    client.post(f'/api/jobs/{job}/run')
    wait_for(lambda:store.get_job(job)['status']=='processing')
    record_presence(store,job,clients[0],'leaving',now=0,sequence=2)
    runner.cancel_disconnected_jobs()
    assert store.get_job(job)['status']=='processing'
    # Same tab returning on a reload renews the lease.
    record_presence(store,job,clients[0],'active',sequence=3)
    runner.cancel_disconnected_jobs()
    assert store.get_job(job)['status']=='processing'
    for key in clients:
        record_presence(store,job,key,'leaving',now=0,sequence=4)
    runner.cancel_disconnected_jobs()
    assert store.get_job(job)['status']=='cancelled'
    assert 'browser session disconnected' in store.get_job(job)['message']


def test_disconnect_is_scoped_to_one_id_and_cancels_queued_work(server,monkeypatch):
    from altanalyze3.components.cellHarmony.flask.browser_presence import record_presence
    app,client=server;store=app.state.job_store;runner=app.state.job_runner
    first,other=new_job(app),new_job(app)
    monkeypatch.setattr(runner,'_memory_wait_reason',lambda:'Waiting for memory.')
    for job in (first,other):client.post(f'/api/jobs/{job}/run')
    record_presence(store,first,'browser_client_first','leaving',now=0)
    runner.cancel_disconnected_jobs()
    assert store.get_job(first)['status']=='cancelled'
    assert store.get_job(other)['status']=='queued'


def test_close_notification_and_timeout_are_ordered_and_reload_safe(tmp_path):
    from altanalyze3.components.cellHarmony.flask.browser_presence import record_presence,browser_disconnected
    store=JobStore(tmp_path);job=store.create_job('human','fixture',None,[])['job_id']
    key='browser_client_first'
    record_presence(store,job,key,'active',now=100,sequence=1)
    assert not browser_disconnected(store,job,now=219)
    assert browser_disconnected(store,job,now=221)  # Browser crashed; no beacon received.
    record_presence(store,job,key,'leaving',now=101,sequence=3)
    record_presence(store,job,key,'active',now=102,sequence=2)  # Delayed, older heartbeat.
    assert browser_disconnected(store,job,now=122)
    record_presence(store,job,'browser_reloaded_page','active',now=110,sequence=1)
    assert not browser_disconnected(store,job,now=122)


def test_invalid_presence_does_not_register_a_lease(server):
    app,client=server;job=new_job(app)
    for payload in ({'client_id':'../../outside'}, {'client_id':'browser_client_first','state':'kill-all'},[] ):
        assert client.post(f'/api/jobs/{job}/browser-presence',json=payload).status_code==400


def test_new_supervisor_recovers_persisted_browser_leases(server):
    from altanalyze3.components.cellHarmony.flask.browser_presence import record_presence
    app,client=server;store=app.state.job_store;job=new_job(app)
    client.post(f'/api/jobs/{job}/run')
    wait_for(lambda:store.get_job(job)['status']=='processing')
    record_presence(store,job,'browser_client_first','leaving',now=0)
    other_supervisor=tasks.JobRunner(store,Path(app.state.config['REFERENCE_REGISTRY']),isolate_jobs=True)
    try:
        assert not other_supervisor._futures
        other_supervisor.cancel_disconnected_jobs()
        assert store.get_job(job)['status']=='cancelled'
    finally:other_supervisor.executor.shutdown(wait=True)


@pytest.mark.parametrize('finishing_worker', [False, True])
def test_browser_close_preserves_completed_session_and_results(server, finishing_worker):
    from concurrent.futures import Future
    from altanalyze3.components.cellHarmony.flask.browser_presence import record_presence
    app,client=server;store=app.state.job_store;runner=app.state.job_runner;job=new_job(app)
    result=store.logs_dir(job).parent/'outputs'/'completed_result.tsv'
    result.write_text('cell\tstate\ncell1\tAM\n')
    store.update_job(job,status='completed',progress=100,message='Job finished successfully.',
                     differential={'status':'completed','run_id':'saved-comparison'})
    before=store.get_job(job)
    record_presence(store,job,'browser_client_finished','leaving',now=0)
    lease=execution_lease(store,job) if finishing_worker else None
    # Cover the interval after completion is saved, before the future and lock finish.
    pending=Future()
    if finishing_worker:runner._futures[f'pipeline:{job}']=pending
    try:
        runner.cancel_disconnected_jobs()
        runner.cancel_job(job,only_if_disconnected=True)
        assert store.get_job(job)==before
        assert result.read_text()=='cell\tstate\ncell1\tAM\n'
        assert not (result.parent.parent/'.cancel_requested').exists()
    finally:
        if lease:lease.close()
        if finishing_worker:
            pending.set_result(None)
            runner._futures.pop(f'pipeline:{job}')


def test_browser_close_stops_active_differential_but_preserves_completed_analysis(server):
    from altanalyze3.components.cellHarmony.flask.browser_presence import record_presence
    app,client=server;store=app.state.job_store;runner=app.state.job_runner;job=new_job(app)
    result=store.logs_dir(job).parent/'outputs'/'completed_result.tsv'
    result.write_text('retained analysis')
    store.update_job(job,status='completed',progress=100)
    runner.submit_differential(job)
    wait_for(lambda:(store.get_job(job).get('differential') or {}).get('status')=='processing')
    record_presence(store,job,'browser_client_finished','leaving',now=0)
    runner.cancel_disconnected_jobs()
    assert store.get_job(job)['status']=='completed'
    assert store.get_job(job)['differential']['status']=='cancelled'
    assert result.read_text()=='retained analysis'
