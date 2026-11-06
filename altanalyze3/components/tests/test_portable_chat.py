import io
import json
from types import SimpleNamespace

import pytest
from altanalyze3.components.cellHarmony.webapp import local_router, local
from importlib import import_module
web = import_module('altanalyze3.components.cellHarmony.webapp.app')


@pytest.fixture
def routing(monkeypatch):
    from altanalyze3.components.cellHarmony.webapp import feature_lookup
    monkeypatch.setattr(feature_lookup, 'expression_lookup', lambda cache: SimpleNamespace(suggestions=['SFTPC']))
    monkeypatch.setattr(web, '_chat_states_by_size', lambda cache: [('AT1',3),('AT2',3)])
    monkeypatch.setattr(web, '_groupable_columns', lambda cache: [{'field':'disease'}])
    monkeypatch.setattr(web, '_chat_genes_in_question', lambda *args: [])
    monkeypatch.setattr(web, 'CHAT_ASSISTANT_URL', '')
    monkeypatch.delenv('SCALABLE_SUMMARY_URL', raising=False)
    monkeypatch.delenv('SCALABLE_SUMMARY_MODEL', raising=False)
    return {'adata':SimpleNamespace(obs=SimpleNamespace(columns=['disease']))}


def test_builtin_chat_requires_no_site_or_model(routing,monkeypatch):
    monkeypatch.setattr(web.urllib.request,'urlopen',lambda *args,**kw:pytest.fail('No external dependency'))
    result=web._chat_read_question('What are the marker genes of AT2?',routing,{})
    assert result['intent']=='markers' and result['cell_state']=='AT2'
    result=web._chat_read_question('Help me understand this dataset',routing,{})
    assert result['intent']=='clarify'


@pytest.mark.parametrize('bad', [{'cell_state':'invented'}, {'genes':['invented']},
                                 {'intent':'run_shell'}, {'covariate':'invented'}, {'genes':'SFTPC'}])
def test_model_cannot_invent_dataset_slots(monkeypatch,bad):
    monkeypatch.setenv('SCALABLE_SUMMARY_URL','http://localhost:9/v1/chat/completions')
    monkeypatch.setenv('SCALABLE_SUMMARY_MODEL','fixture')
    answer={'intent':'markers','cell_state':'AT2','genes':[],**bad}
    reply={'choices':[{'message':{'content':json.dumps(answer)}}]}
    monkeypatch.setattr(local_router.urllib.request,'urlopen',lambda *args,**kw:io.BytesIO(json.dumps(reply).encode()))
    assert local_router.read_with_model('question',states=['AT1','AT2'],genes=['SFTPC'],covariates=['disease'],modalities=['rna']) is None


def test_model_uses_same_summary_service_for_unhandled_question(routing,monkeypatch):
    monkeypatch.setenv('SCALABLE_SUMMARY_URL','http://localhost:9/v1/chat/completions')
    monkeypatch.setenv('SCALABLE_SUMMARY_MODEL','fixture')
    calls=[]
    def request(req,timeout):
        calls.append((req.full_url,json.loads(req.data),timeout))
        return io.BytesIO(json.dumps({'choices':[{'message':{'content':json.dumps({'intent':'markers','cell_state':'AT2'})}}]}).encode())
    monkeypatch.setattr(local_router.urllib.request,'urlopen',request)
    result=web._chat_read_question('Tell me how AT2 is identified',routing,{})
    assert result['intent']=='markers' and result['router']=='local_summary_model'
    assert len(calls)==1 and calls[0][0].endswith('/v1/chat/completions') and calls[0][2]==30


def test_split_model_requires_all_shards_and_selects_first(tmp_path):
    first=tmp_path/'model-00001-of-00002.gguf'; second=tmp_path/'model-00002-of-00002.gguf'
    first.write_bytes(b'GGUF')
    with pytest.raises(ValueError,match='incomplete'):local.validate_model(first)
    second.write_bytes(b'GGUF')
    assert local.validate_model(second)==first
    assert local.model_environment(8090)['SCALABLE_SUMMARY_URL']=='http://127.0.0.1:8090/v1/chat/completions'


def test_launcher_supplies_environment_and_stops_its_model(tmp_path,monkeypatch):
    model=tmp_path/'fixture.gguf';model.write_bytes(b'GGUF')
    monkeypatch.setattr(local.sys,'argv',['local','--model',str(model),'--port','18000','--model-port','18086'])
    monkeypatch.setattr(local.importlib.util,'find_spec',lambda key:object())
    monkeypatch.setattr(local.urllib.request,'urlopen',lambda *args,**kw:io.BytesIO(b'{}'))
    class Socket:
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def bind(self,address):pass
    monkeypatch.setattr(local.socket,'socket',Socket)
    calls=[]
    class Process:
        def __init__(self,cmd,**kwargs):self.cmd=cmd;self.kwargs=kwargs;self.stopped=False;calls.append(self)
        def poll(self):return 0 if self.stopped else None
        def wait(self,timeout=None):
            if 'uvicorn' in self.cmd:self.stopped=True
            return 0
        def terminate(self):self.stopped=True
        def kill(self):self.stopped=True
    monkeypatch.setattr(local.subprocess,'Popen',Process)
    monkeypatch.setattr(local.signal,'signal',lambda *args:None)
    local.main()
    assert len(calls)==2
    assert calls[1].kwargs['env']['SCALABLE_SUMMARY_URL']=='http://127.0.0.1:18086/v1/chat/completions'
    assert calls[1].kwargs['env']['SCALABLE_SUMMARY_MODEL']=='Qwen2.5-7B-Instruct-Q4_K_M'
    assert all(p.stopped for p in calls)


def test_existing_smoke_fixture_passes_default_qc_and_inspection(tmp_path):
    import h5py
    from altanalyze3.components.cellHarmony.webapp.dev.create_existing_fixture import create
    from altanalyze3.components.cellHarmony.webapp.existing_annotations import inspect_annotations
    from altanalyze3.components.anndata_compat import read_elem, sparse_dataset
    path=create(tmp_path/'synthetic.h5ad')
    report=inspect_annotations([path])
    assert report['selected']=='cell_type' and report['embedding_key']=='X_umap'
    assert report['inputs'][0]['count_source']=='layers/counts'
    with h5py.File(path,'r') as f:
        counts=sparse_dataset(f['layers/counts']).to_memory()
        assert ((counts>0).sum(axis=1)>=500).all() and (counts.sum(axis=1)>=1000).all()
        assert read_elem(f['uns'])['synthetic_fixture']
    with pytest.raises(FileExistsError):create(path)
