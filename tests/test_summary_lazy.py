"""Progressive summaries retain evidence and infer only explicitly opened results."""
import importlib
import json
import threading
import time

from fastapi.testclient import TestClient
from altanalyze3.components.cellHarmony.webapp import biology_summary as summary
from altanalyze3.components.cellHarmony.webapp import biology_interpretation as interpretation
from altanalyze3.components.cellHarmony.webapp import biology_literature as literature
from test_biology_summary import fixture
from test_biology_interpretation import digest, response


def test_surface_has_identical_computed_evidence_and_no_remote_calls(tmp_path,monkeypatch):
    store,job=fixture(tmp_path);meta=store.get_job(job)
    calls=[]
    monkeypatch.setattr(summary,'interpret',lambda d:calls.append('model') or (None,''))
    monkeypatch.setattr(interpretation,'interpret_sections',lambda d:calls.append('sections') or {})
    monkeypatch.setattr(literature,'review',lambda d:calls.append('literature') or {})
    full=summary.build_summary(meta)
    assert calls==['literature','model','sections']
    calls.clear()
    surface=summary.build_surface(meta)
    assert calls==[]
    for key in ('evidence','observations','comparisons','synthesis','question_framework',
                'regulatory_review','integrated_pathway_review','missing_reference_states','methods'):
        assert surface['digest'][key]==full['digest'][key],key
    assert surface['overview']==full['overview'] and surface['hypotheses']==full['hypotheses']
    assert surface['digest']['section_interpretation']['status']=='deferred'
    assert surface['digest']['deferred_sections']['pathways'][0]['source_id']=='P1'


def test_single_source_interpretation_preserves_exact_prompt_and_guards(monkeypatch):
    d=digest();d['regulatory_review']['rows'].append(dict(d['regulatory_review']['rows'][0],state='AT2'))
    text='Could the predicted matrix-remodeling program in AM affect tissue repair rather than reflect a change in cell identity?'
    response(monkeypatch,text)
    captured=[];original=interpretation.urllib.request.urlopen
    def capture(req,*args,**kw):
        captured.append(json.loads(req.data));return original(req,*args,**kw)
    monkeypatch.setattr(interpretation.urllib.request,'urlopen',capture)
    selected=interpretation.interpret_sections(d,section='grn',source_id='G1')
    assert len(captured)==1 and set(selected['sources'])=={'G1'}
    first_prompt=captured[0]
    assert selected['sections']['grn'][0]['state']=='AM'
    captured.clear()
    interpretation.interpret_sections(d)
    assert len(captured)==2 and captured[0]==first_prompt
    assert 'AT2' not in json.loads(first_prompt['messages'][1]['content'])['findings'][0]['state']


def test_detail_routes_are_cached_source_bound_and_do_not_block_surface(tmp_path,monkeypatch):
    store,job=fixture(tmp_path)
    web=importlib.import_module('altanalyze3.components.cellHarmony.webapp.app')
    app=web.create_app({'JOB_STORAGE':str(store.root)})
    reached=threading.Event();release=threading.Event();calls=[]
    def detail(surface,section,source_id=None):
        calls.append((section,source_id));reached.set();assert release.wait(5)
        return {'kind':'summary_interpretation','review':{'status':'completed','sections':{},'sources':{}}}
    monkeypatch.setattr(summary,'build_summary_detail',detail)
    def poll(client,url):
        for _ in range(300):
            r=client.post(url)
            if r.status_code!=202:return r
            time.sleep(.01)
        raise AssertionError('Request did not finish')
    try:
        with TestClient(app) as client:
            base=f'/api/jobs/{job}/chat-summary'
            surface=poll(client,base);assert surface.status_code==200,surface.text
            ident=surface.json()['source_identity']
            url=base+'-section?section=pathways&source_id=P1&expected_identity='+ident
            assert client.post(url).status_code==202 and reached.wait(5)
            start=time.monotonic()
            assert client.post(base).content==surface.content
            assert time.monotonic()-start<.5
            assert client.post(base+'-section?section=grn&source_id=G999').status_code==422
            assert client.post(base+'-section?section=pathways&source_id=P1&expected_identity=old').status_code==409
            release.set();result=poll(client,url);assert result.status_code==200,result.text
            assert client.post(url).content==result.content and calls==[('pathways','P1')]
            store.update_job(job,status='processing')
            assert client.post(url).status_code==409
            assert client.post('/api/jobs/missing/chat-summary-section?section=pathways&source_id=P1').status_code==404
    finally:
        release.set()


def test_failed_detail_can_be_retried_without_caching_transport_errors():
    from fastapi.responses import JSONResponse
    from altanalyze3.components.cellHarmony.webapp.plot_build import PlotBuildQueue
    queue=PlotBuildQueue(cache_errors=False);calls=[]
    def build():
        calls.append(1)
        return JSONResponse({'detail':'fixture outage'} if len(calls)==1 else {'done':True},
                            status_code=503 if len(calls)==1 else 200)
    def poll():
        for _ in range(200):
            result=queue.request('source',build)
            if result.status_code!=202:return result
            time.sleep(.005)
        raise AssertionError('Background detail did not finish')
    try:
        assert poll().status_code==503
        assert poll().status_code==200
        assert poll().status_code==200 and len(calls)==2
    finally:queue.close()
