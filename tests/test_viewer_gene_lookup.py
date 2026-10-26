"""Unknown viewer genes are normal UI outcomes; aliases read the native gene values."""
import threading
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from fastapi.testclient import TestClient

from altanalyze3.components.visualization.scalable_viewer import bundle_meta as B
from altanalyze3.components.visualization.scalable_viewer import data_api as D
from altanalyze3.components.visualization.scalable_viewer import scalable_app as S
from altanalyze3.components.visualization.scalable_viewer.server import create_app as fast_app


def dataset(tmp_path, names=('MRC1', 'AGER')):
    genes = tmp_path / 'genes.tsv'
    genes.write_text('index\tgene_id\tsymbol\n' + ''.join(
        f'{i}\tENSG{i+1:011d}\t{name}\n' for i,name in enumerate(names)))
    clusters = tmp_path / 'clusters.tsv'
    clusters.write_text('barcode\tcluster\n' + ''.join(f'c{i}\t{("AM" if i<3 else "AT2")}\n' for i in range(6)))
    ds = D.Dataset.__new__(D.Dataset)
    ds.id = 'fixture'; ds.label = 'fixture'; ds.n_cells = 6; ds.n_genes = len(names)
    ds.cluster_key = 'cell_type'; ds.states = ['AM', 'AT2']; ds.state_n = [3,3]
    ds.paths = SimpleNamespace(genes=str(genes), clusters=str(clusters), bundle_dir=str(tmp_path))
    ds._lock = threading.Lock(); ds._symbols = ds._gene_ids = ds._sym_index = None
    ds._cells = {'state_code':np.array([0,0,0,1,1,1]), 'embedding':np.arange(12).reshape(6,2)}
    ds.sv = {'covariates':{}}; ds.meta = {}; ds.colors = {}; ds._modality_stores = {}
    ds._mean = np.tile([2.,5.], (len(names),1)); ds._frac = np.ones((len(names),2))
    ds.gene_column = lambda row: (np.arange(6,dtype=np.uint32), np.arange(1,7,dtype=np.float32)*(row+1))
    return ds


@pytest.fixture
def viewer(tmp_path):
    ds = dataset(tmp_path)
    catalog = SimpleNamespace(get=lambda job: ds, ids=lambda:['fixture'], entries=[])
    app = B.W.create_app({'JOB_STORAGE':str(tmp_path/'jobs')})
    store = B.BundleJobStore(catalog, tmp_path, {})
    app.state.job_store = store
    store._meta['fixture'] = {'job_id':'fixture', 'status':'completed', 'cluster_key':'cell_type',
                             'species':'Hs', 'reference':'fixture', 'scalable_viewer':{}, 'artifacts':{}}
    S._install_violin_covariate(app)
    try:
        with TestClient(app) as client:
            yield app, client, ds
    finally:
        app.state.job_runner.executor.shutdown(wait=True)


@pytest.mark.parametrize('view', ['all','umap','violin'])
@pytest.mark.parametrize('compact', [False,True])
def test_cd206_expression_matches_mrc1_in_viewer_routes(viewer,view,compact):
    _,client,ds = viewer
    options = {'view':view, 'compact':compact, 'covariate':'cell_type'}
    alias = client.get('/api/jobs/fixture/expression',params={**options,'gene':'CD206'})
    canonical = client.get('/api/jobs/fixture/expression',params={**options,'gene':'MRC1'})
    assert alias.status_code == canonical.status_code == 200,alias.text
    a,b = alias.json(),canonical.json()
    assert a['resolved_gene'] == b['resolved_gene'] == 'MRC1'
    for field in ('umap','violin','scatter','global_min','global_max'):
        assert a[field] == b[field]
    assert ds.symbols == ['MRC1','AGER']


@pytest.mark.parametrize('view', ['all','umap','violin'])
def test_missing_gene_returns_a_message_without_reference_or_axis_work(viewer,view,monkeypatch):
    _,client,_ = viewer
    def forbidden(*args,**kwargs):
        pytest.fail('An absent bundle gene must not read a different reference or build axes')
    monkeypatch.setattr(B.W,'_build_reference_expression_payload',forbidden)
    monkeypatch.setattr(B.W,'_parent_plot_context',forbidden)
    for _ in range(3):
        response = client.get('/api/jobs/fixture/expression',params={
            'gene':'THIS_GENE_IS_ABSENT','view':view,'covariate':'cell_type','compact':True})
        assert response.status_code == 200,response.text
        body = response.json()
        assert body['source'] == 'missing' and body['resolved_gene'] is None
        assert 'was not found in this dataset' in body['message']
        assert body['umap'] == body['violin'] == body['scatter'] == []


@pytest.mark.parametrize('gene', ['CD206','THIS_GENE_IS_ABSENT'])
@pytest.mark.parametrize('mode', ['umap','violin'])
def test_expression_pdf_handles_aliases_and_missing_features(viewer,gene,mode):
    _,client,_ = viewer
    response = client.get('/api/jobs/fixture/expression/pdf',params={'gene':gene,'mode':mode,'covariate':'cell_type'})
    assert response.status_code == 200,response.text[:100]
    assert response.content.startswith(b'%PDF-')


def test_ambiguous_alias_is_a_client_message_not_a_500(viewer):
    app,client,_ = viewer
    cache = B.W._get_expression_cache(app,app.state.job_store.get_job('fixture'))
    cache['adata'].var = pd.DataFrame({'gene_symbols':['SHARED','SHARED']},index=cache['var_names'])
    for route in ('expression','expression/pdf'):
        response = client.get('/api/jobs/fixture/'+route,params={'gene':'SHARED'})
        assert response.status_code == 400,response.text
        assert 'multiple primary IDs' in response.json()['detail']


def test_native_bundle_resolver_search_and_binary_api_share_cd206_alias(tmp_path):
    ds = dataset(tmp_path)
    assert ds.resolve_gene('CD206') == ds.resolve_gene('MRC1') == 0
    assert ds.resolve_gene(' CD-206 ') == 0
    assert ds.search_genes('CD206') == ['MRC1']
    assert ds.resolve_gene('ABSENT') is None
    app = fast_app(SimpleNamespace(get=lambda job:ds,ids=lambda:['fixture']))
    with TestClient(app) as client:
        a = client.get('/api/dataset/fixture/expression.bin',params={'gene':'CD206'})
        b = client.get('/api/dataset/fixture/expression.bin',params={'gene':'MRC1'})
        assert a.status_code == b.status_code == 200
        assert a.content == b.content
        missing = client.get('/api/dataset/fixture/expression.bin',params={'gene':'ABSENT'})
        assert missing.status_code == 404 and 'gene not in this dataset' in missing.json()['detail']


def test_native_bundle_alias_does_not_override_features_or_guess_between_duplicates(tmp_path):
    ds = dataset(tmp_path, ('CD206','MRC1'))
    assert ds.resolve_gene('CD206') == 0 and ds.resolve_gene('MRC1') == 1
    ds = dataset(tmp_path, ('MRC1','MRC1'))
    assert ds.resolve_gene('CD206') is None
