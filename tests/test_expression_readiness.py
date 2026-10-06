"""A partial analytical output must never trigger a whole-H5AD serving read."""
from types import SimpleNamespace
import importlib
import pytest

W=importlib.import_module('altanalyze3.components.cellHarmony.webapp.app')


@pytest.mark.parametrize('status',['queued','processing','running'])
def test_processing_requests_do_not_open_matrices(status,monkeypatch):
    app=SimpleNamespace(state=SimpleNamespace(job_store=object()))
    monkeypatch.setattr(W._job_bundle,'view',lambda *a,**k:pytest.fail('Bundle opened before readiness'))
    monkeypatch.setattr(W,'_read_cached_h5ad',lambda *a,**k:pytest.fail('Whole matrix opened before readiness'))
    with pytest.raises(FileNotFoundError,match='being finalized'):
        W._get_expression_cache(app,{'job_id':'large','status':status,'cluster_key':'state'})


def test_new_bundle_invalidates_previous_h5ad_cache(tmp_path,monkeypatch):
    source=tmp_path/'output.h5ad';source.write_bytes(b'fixture')
    stat=source.stat()
    old={'h5ad_path':str(source),'umap_path':'','cluster_key':'state','modality':'rna',
         'source_stamp':(stat.st_mtime_ns,stat.st_size),'bundle_signature':(None,)*4}
    app=SimpleNamespace(state=SimpleNamespace(job_store=object(),expression_cache={'large:rna':old},expression_cache_locks={}))
    monkeypatch.setattr(W,'_modality_h5ad_path',lambda *a:str(source))
    class BundleReached(Exception):pass
    def view(*a):raise BundleReached()
    monkeypatch.setattr(W._job_bundle,'view',view)
    monkeypatch.setattr(W,'_read_cached_h5ad',lambda *a,**k:pytest.fail('Whole H5AD used despite ready bundle'))
    with pytest.raises(BundleReached):
        W._get_expression_cache(app,{'job_id':'large','status':'completed','cluster_key':'state',
                                     'bundle':{'status':'completed','dir':'bundle','prefix':'job','built_utc':'today'}})


def test_viewer_store_keeps_its_own_readiness_handling():
    expected=object()
    store=SimpleNamespace(get_expression_cache=lambda *a:expected)
    app=SimpleNamespace(state=SimpleNamespace(job_store=store))
    assert W._get_expression_cache(app,{'job_id':'viewer','status':'processing'}) is expected
