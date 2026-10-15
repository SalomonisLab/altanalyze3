"""Scientific/UI integration checks: native method parity and retained state identities."""
import ast
import importlib
import json
import hashlib
from pathlib import Path
import time

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr
from fastapi.testclient import TestClient

from altanalyze3.components.cellHarmony.webapp import biology_sections as sections,drug_chat,biology_literature
from altanalyze3.components.cellHarmony.discover_resources import drug_matching as dm
from test_biology_summary import fixture


@pytest.fixture(autouse=True)
def literature_transport(monkeypatch):
    monkeypatch.setattr(biology_literature,'search',lambda *a:{'status':'fixture','papers':[],'note':'fixture'})


def test_methods_execution_version_not_current_registry(tmp_path):
    store,job=fixture(tmp_path);m=store.get_job(job)
    m['model_versions']={'lipids':{'model_version':'historic-test','training_datasets':['Retained source DOI'], 'validation_summary':{'external_holdout':False}}}
    r=sections.methods(m)
    assert r['models'][0]['execution_record']['model_version']=='historic-test'
    assert r['models'][0]['training_and_validation']['training_datasets']==['Retained source DOI']
    assert r['sections'][0]['title']=='Analysis identity and input'
    assert r['provenance_gaps']


def test_native_drug_method_matches_pinned_source_ast():
    source=Path(dm.__file__)
    manifest=json.loads(source.with_name('source_manifest.json').read_text())
    nodes={n.name:n for n in ast.parse(source.read_text()).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
    for name,expected in manifest['native_method_ast_sha256'].items():
        assert hashlib.sha256(ast.dump(nodes[name],include_attributes=False).encode()).hexdigest()==expected,name


def test_literature_uses_documented_state_names_without_gene_synonyms():
    query=biology_literature.query_terms('SPI1',('AM','alveolar macrophage'),'IPF','human')
    assert 'SPI1' in query and 'SP1' not in query and 'alveolar macrophage' in query


def test_verified_claim_requires_exact_identity_direction_and_audited_source():
    row=dict(species='human',state='AT2',disease='IPF',feature='SPI1',modality='rna',effect=2)
    record=dict(row,status='reviewed',direction='up',reviewed_by='test fixture',reviewed_at='2026-10-11',claim='Fixture only',evidence_excerpt='Fixture only',citation_url='https://example.org/fixture')
    assert biology_literature.verified_sources(row,[record])==[record]
    for key,value in [('feature','SP1'),('species','mouse'),('state','AM'),('direction','down'),('status','pending')]:
        assert not biology_literature.verified_sources(row,[dict(record,**{key:value})])


def test_fact_summary_retains_retrieval_uncertainty():
    row={'compound':'Example','rho':-.4,'cell_line':'A549','dose':'1 uM','time':'24h'}
    text=drug_chat.fact_summary(row,{'fda':{'fda_approved':None},'trials':{'errors':['HTTP 503'],'n_in_condition':0}},'IPF','AT2','opposing')
    assert 'unresolved' in text and 'incomplete' in text and '0 registered trials' not in text


def test_concentrated_missing_state_and_diffuse_negative(tmp_path):
    n=200;obs=pd.DataFrame({'supervised_state':['Unaligned']*80+['AM']*120,
        'unsupervised_cluster':['c18']*80+['c2']*120,
        'unsupervised_state':['Fibroblast_c18']*80+['Macrophage_c2']*120},index=[str(i) for i in range(n)])
    a=ad.AnnData(np.zeros((n,1)),obs=obs);p=tmp_path/'a.h5ad';a.write_h5ad(p)
    with h5py.File(p,'a') as f:del f['X']
    r=sections.missing_states({'analysis_mode':'both'},str(p))
    assert r['rows'][0]['cluster']=='c18' and r['rows'][0]['concentrated']
    assert r['rows'][0]['names']=={'Fibroblast_c18':80}
    assert sum(v['cells'] for v in r['rows'])==200
    # Same total unaligned cells scattered evenly across states must not nominate a state.
    obs['supervised_state']=['Unaligned' if i%5<2 else 'AM' for i in range(n)]
    a=ad.AnnData(np.zeros((n,1)),obs=obs);a.write_h5ad(p)
    r=sections.missing_states({'analysis_mode':'both'},str(p));assert not any(v['candidate'] for v in r['rows'])


def test_native_spearman_ties_directions_neutral_and_one_row_per_signature():
    q=np.array([1,1,-2,3,5,0,7,8,9,10.])
    block=np.column_stack([q,-q,q[::-1],np.zeros(len(q))])
    got,shared=dm._spearman_block(q,block)
    for i in range(3):assert got[i]==pytest.approx(spearmanr(q,block[:,i]).statistic)
    counts=dm._direction_counts(q,block)
    assert np.isfinite(got[:3]).all()
    assert dm.METHOD_VERSION=='cellref2-dart-like-signature-v1'
    assert counts['drug_zero'][-1]==len(q)


def test_drug_routing_does_not_guess_state_or_contrast(tmp_path):
    store,job=fixture(tmp_path);m=store.get_job(job)
    assert drug_chat.read(m,'Which drugs reverse my signature?')['kind']=='drug_match_choices'
    r=drug_chat.read(m,'Which drugs mimic the AT2 signature with iLINCS?')
    assert r['cell_state']=='AT2' and r['run_id']=='run_rna' and r['gene_space']=='measured' and r['tail']=='mimicking'
    prepared,_=drug_chat.prepare(m,'run_rna','AT2')
    assert prepared.genes==['SFTPC']
    m['differential_history']['duplicate']=dict(m['differential_history']['run_rna'],run_id='duplicate')
    assert drug_chat.read(m,'drugs reverse AT2 signature')['kind']=='drug_match_choices'


def test_native_network_hubs_all_states_and_counts(tmp_path):
    store,job=fixture(tmp_path);m=store.get_job(job);path=Path(m['artifacts']['combined_h5ad'])
    obs=pd.DataFrame({'supervised_state':['AT2']*4+['Unaligned']*4,'disease':['IPF','IPF','Healthy','Healthy']*2,'Library':['D1','D1','D2','D2']*2},index=[str(i) for i in range(8)])
    a=ad.AnnData(np.tile(np.array([[2,10,8]],dtype='f4'),(8,1)),obs=obs,var=pd.DataFrame(index=['STAT3','COL1A1','CTGF']));a.layers['counts']=a.X.copy();a.write_h5ad(path)
    for mod,features in [('rna',['COL1A1','CTGF','STAT3']),('grn',['STAT3|COL1A1','TEAD2|CTGF']),('grn_tf',['STAT3','TEAD2'])]:
        run=m['differential_history']['run_'+mod];p=run['artifacts']['DEG_detailed_fixture'];pd.DataFrame([{'population':'AT2','gene':g,'log2fc':-2,'fdr':.01,'pval':.001} for g in features]).to_csv(p,sep='\t',index=False)
    ds=sections.SavedInterpretData(m,'run_rna');case,control=ds.arm_expression('run_rna','AT2',['STAT3'])
    assert case['STAT3']==pytest.approx(np.log2(1+10000*2/20));assert case==control
    from altanalyze3.components.cellHarmony.webapp.biology_summary import collect_evidence
    d=collect_evidence(m);r=sections.regulatory_sections(m,d)
    assert r['states_examined']==2
    hit=next(v for v in r['rows'] if v['state']=='AT2')
    assert hit['status']=='completed' and hit['hubs'][0]['factor']=='STAT3'
    assert hit['hubs'][0]['target_RNA_down']==1 and hit['hubs'][0]['TF_log2fc']==-2
    assert next(v for v in r['rows'] if v['state']=='Unaligned')['n_edges']==0
    # Different arm identities cannot supply a TF comparison.
    m['differential_history']['run_grn_tf']['config']['group2_samples']=['Different']
    assert sections.regulatory_sections(m,d)['states_examined']==2  # digest identity is immutable
    without_tf=sections.regulatory_sections(m,collect_evidence(m))
    assert without_tf['states_examined']==2
    hit=next(v for v in without_tf['rows'] if v['state']=='AT2')
    assert hit['hubs'] and hit['tf_run_id']==''
    assert hit['hubs'][0]['TF_log2fc'] is None
    assert 'No matching TF-activity differential' in hit['tf_activity_status']


def test_pubmed_title_hits_are_not_declared_verified_or_proven_novel(monkeypatch):
    d={'comparisons':[{'id':'r','comparison':'IPF vs Healthy'}],'modalities':['rna'],'species':'human','evidence':[
        {'kind':'differential_feature','modality':'rna','run_id':'r','state':'AT2','feature':'SPI1','id':'E1','title':'SPI1 higher'}]}
    monkeypatch.setattr(biology_literature,'search',lambda *a:{'status':'searched','papers':[{'title':'SP1 in fibrosis','url':'https://pubmed.ncbi.nlm.nih.gov/1/'}]})
    r=biology_literature.review(d);assert not r['verified_known_biology'];assert r['novel_predictions'][0]['feature']=='SPI1'
    assert 'unconfirmed' in r['verification_note']


def test_pathway_review_retains_states_with_unavailable_results(tmp_path,monkeypatch):
    store,job=fixture(tmp_path);m=store.get_job(job)
    from altanalyze3.components.cellHarmony.webapp.biology_summary import collect_evidence
    d=collect_evidence(m)
    def ranked(ds,**params):
        if params['cell_state']=='AT2':raise ValueError('Fixture resource unavailable')
        return {'pathways':[]}
    monkeypatch.setattr(sections.discover_integration,'pathways',ranked)
    r=sections.pathway_sections(m,d)
    assert r['coverage'] and any(v['status']=='unavailable' for v in r['coverage'])
    assert {v['state'] for v in r['coverage']}==set(d['observations']['layers']['supervised_state']['counts'])


def test_async_drug_and_pathway_api(tmp_path,monkeypatch):
    store,job=fixture(tmp_path);web=importlib.import_module('altanalyze3.components.cellHarmony.webapp.app')
    app=web.create_app({'JOB_STORAGE':str(store.root)})
    monkeypatch.setattr(drug_chat,'screen',lambda meta,**params:{'kind':'drug_match','rows':[],'run_id':params['run_id']})
    with TestClient(app) as c:
        response=c.post(f'/api/jobs/{job}/chat',json={'question':'which drugs reverse my signature?'})
        assert response.status_code==200 and response.json()['kind']=='drug_match_choices'
        for path in [f'/api/jobs/{job}/chat-drugs?run_id=run_rna&cell_state=AT2',f'/api/jobs/{job}/chat-interpret?kind=pathway&run_id=run_rna&cell_state=AT2&feature=Surfactant%20metabolism']:
            for _ in range(100):
                response=c.post(path)
                if response.status_code!=202:break
                time.sleep(.01)
            assert response.status_code==200,response.text
        assert c.post(f'/api/jobs/{job}/chat-drugs?run_id=r&cell_state=AT2&tail=garbage').status_code==422


def test_missing_drug_library_returns_clear_503_without_setup_details(tmp_path,monkeypatch,caplog):
    store,job=fixture(tmp_path);web=importlib.import_module('altanalyze3.components.cellHarmony.webapp.app')
    app=web.create_app({'JOB_STORAGE':str(store.root)})
    def missing(*args,**kwargs):raise FileNotFoundError('fixture missing library')
    monkeypatch.setattr(drug_chat,'screen',missing)
    with TestClient(app) as client:
        for _ in range(100):
            response=client.post(f'/api/jobs/{job}/chat-drugs?run_id=run_rna&cell_state=AT2')
            if response.status_code!=202:break
            time.sleep(.01)
        assert response.status_code==503
        detail=response.json()['detail']
        assert 'reference library is missing or incomplete' in detail
        assert 'other results remain available' in detail
        assert 'CELLREF2_DRUGSTORE' not in detail
        assert 'configured directory=' in caplog.text
