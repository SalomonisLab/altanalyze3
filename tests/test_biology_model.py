"""Native inference runs in a subprocess; an abort cannot kill the web server."""
import io
import json
import subprocess
from altanalyze3.components.cellHarmony.webapp import biology_model as model


class Process:
    def __init__(self, response):
        self.stdin=io.StringIO();self.stdout=io.StringIO(response);self.killed=False
    def poll(self):return -9 if self.killed else None
    def kill(self):self.killed=True
    def wait(self,timeout):return -9


def test_worker_recovers_after_native_process_exit(monkeypatch):
    failed=Process('');healthy=Process(json.dumps({'choices':[{'message':{'content':'{}'}}]}))
    created=[]
    def spawn(*a,**k):
        process=[failed,healthy][len(created)];created.append(process);return process
    monkeypatch.setattr(subprocess,'Popen',spawn)
    monkeypatch.setattr(model.select,'select',lambda streams,*a:(streams,[],[]))
    worker=model.ModelWorker({'model_path':'installed.gguf','idle_seconds':300})
    import pytest
    with pytest.raises(RuntimeError,match='process exited'):worker.complete({'messages':[]})
    assert worker.process is None and failed.killed
    assert worker.complete({'messages':[]})['choices']
    assert len(created)==2 and not healthy.killed


def test_worker_timeout_releases_native_memory(monkeypatch):
    process=Process('')
    monkeypatch.setattr(subprocess,'Popen',lambda *a,**k:process)
    monkeypatch.setattr(model.select,'select',lambda *a:([],[],[]))
    worker=model.ModelWorker({'model_path':'installed.gguf','idle_seconds':300})
    import pytest
    with pytest.raises(TimeoutError):worker.complete({'messages':[]})
    assert process.killed and worker.process is None


def test_invalid_schema_is_reported_before_generation(monkeypatch):
    import sys,types
    generated=[]
    class Engine:
        def __init__(self,**kwargs):pass
        def reset(self):pass
        def create_chat_completion(self,**kwargs):generated.append(kwargs);return {}
    class Grammar:
        @staticmethod
        def from_json_schema(schema,verbose=False):raise ValueError('invalid fixture schema')
    monkeypatch.setitem(sys.modules,'llama_cpp',types.SimpleNamespace(Llama=Engine,LlamaGrammar=Grammar))
    options={'model_path':'installed.gguf','ctx_size':4096,'threads':2,'gpu_layers':0}
    request={'options':options,'payload':{'messages':[{'role':'user','content':'fixture'}],
        'response_format':{'type':'json_object','schema':{'invalid':'fixture'}}}}
    monkeypatch.setattr(sys,'stdin',io.StringIO(json.dumps(request)+'\n'))
    output=io.StringIO();monkeypatch.setattr(sys,'stdout',output)
    model.worker_main()
    result=json.loads(output.getvalue())
    assert 'invalid fixture schema' in result['error'] and generated==[]
