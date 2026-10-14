"""Summary progress and exact-prompt reuse preserve provider/source boundaries."""
import json
import os
import threading
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from fastapi.responses import JSONResponse
from altanalyze3.components.cellHarmony.webapp import summary_runtime as runtime
from altanalyze3.components.cellHarmony.webapp import biology_literature as literature
from altanalyze3.components.cellHarmony.webapp.plot_build import PlotBuildQueue


def request(*, content='complete source-bound evidence',model='model-a',url='http://model/completions'):
    return urllib.request.Request(url,data=json.dumps({'model':model,'messages':[{'role':'user','content':content}]}).encode(),
                                  headers={'Authorization':'Bearer fixture-secret'})


class Reply:
    def __enter__(self):return self
    def __exit__(self,*args):pass
    def read(self,*args):return json.dumps({'choices':[{'message':{'content':'{}'}}],'usage':{'completion_tokens':3}}).encode()


def test_exact_cache_survives_reporting_contexts_and_never_stores_credentials(tmp_path,monkeypatch):
    calls=[]
    monkeypatch.setattr(runtime.urllib.request,'urlopen',lambda *a,**k:calls.append(1) or Reply())
    with runtime.reporting(lambda **event:None,tmp_path/'job-a'),runtime.trace() as timings:
        first=runtime.completion(request(),'model.first','Interpreting…')
    with runtime.reporting(lambda **event:None,tmp_path/'job-a'),runtime.trace() as cached_timings:
        second=runtime.completion(request(),'model.first','Interpreting…')
    assert first==second and len(calls)==1
    assert timings[0]['cache_hit'] is False and cached_timings[0]['cache_hit'] is True
    assert cached_timings[0]['usage']['completion_tokens']==3
    for p in tmp_path.rglob('*.json'):
        assert 'fixture-secret' not in p.read_text() and 'Authorization' not in p.read_text()
        assert p.stat().st_mode & 0o077==0
    # A source, model, endpoint or job change never reuses another result.
    with runtime.reporting(None,tmp_path/'job-a'):
        for req in [request(content='changed evidence'),request(model='other'),request(url='http://other/completions')]:
            runtime.completion(req,'model.first','Interpreting…')
    with runtime.reporting(None,tmp_path/'job-b'):
        runtime.completion(request(),'model.first','Interpreting…')
    assert len(calls)==5
    for p in (tmp_path/'job-a').glob('*.json'):os.utime(p,(0,0))
    with runtime.reporting(None,tmp_path/'job-a'):
        runtime.completion(request(),'model.first','Interpreting…')
    assert len(calls)==6


def test_completion_cache_is_bounded_and_corrupt_files_are_rebuilt(tmp_path,monkeypatch):
    monkeypatch.setattr(runtime.urllib.request,'urlopen',lambda *a,**k:Reply())
    monkeypatch.setattr(runtime,'CACHE_ENTRIES',2)
    with runtime.reporting(None,tmp_path):
        for i in range(4):runtime.completion(request(content=str(i)),'model','Interpreting…')
        assert len(list(tmp_path.glob('*.json')))==2
        for p in tmp_path.glob('*.json'):p.write_text('invalid JSON')
        assert runtime.completion(request(content='3'),'model','Interpreting…')['choices']


def test_queue_pending_reports_current_stage_and_separate_queued_work():
    queue=PlotBuildQueue(max_entries=2)
    reached=threading.Event();release=threading.Event()
    def build():
        queue.report_progress('a',stage='model.section.P1',detail='Interpreting pathway 1 of 6…')
        reached.set();assert release.wait(5)
        return JSONResponse({'done':True})
    try:
        assert queue.request('a',build).status_code==202
        assert reached.wait(5)
        progress=json.loads(queue.request('a',build).body)
        assert progress['stage']=='model.section.P1' and '1 of 6' in progress['detail']
        assert progress['elapsed_seconds']>=0
        other=json.loads(queue.request('b',lambda:JSONResponse({})).body)
        assert other['stage']=='queued' and '1 of 6' not in other['detail']
    finally:
        release.set();queue.close()


def test_parallel_literature_checks_preserve_exact_results_and_order(monkeypatch):
    d={'comparisons':[{'id':'r','comparison':'Disease vs Control'}],'modalities':['rna'],
       'species':'human','evidence':[{'id':f'E{i}','kind':'differential_feature','modality':'rna',
           'run_id':'r','state':f'state{i}','feature':f'gene{i}','title':f'finding{i}','facts':{'log2fc':1}}
           for i in range(8)]}
    lock=threading.Lock();active=[0];peak=[0]
    def search(*args):
        with lock:active[0]+=1;peak[0]=max(peak[0],active[0])
        time.sleep(.025)
        with lock:active[0]-=1
        return {'status':'searched','papers':[],'query':str(args),'note':'fixture transport'}
    monkeypatch.setattr(literature,'search',search)
    monkeypatch.setattr(literature,'ThreadPoolExecutor',lambda max_workers:ThreadPoolExecutor(max_workers=1))
    expected=literature.review(d)
    peak[0]=0
    monkeypatch.setattr(literature,'ThreadPoolExecutor',ThreadPoolExecutor)
    actual=literature.review(d)
    assert actual==expected and peak[0]>1
    assert [r['feature'] for r in actual['novel_predictions']]==[f'gene{i}' for i in range(8)]
