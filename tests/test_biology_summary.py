"""Read-only summary contracts, all modalities/layers, model grounding and API polls."""
import io
import json
from pathlib import Path
import time
from importlib import import_module

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
from fastapi.testclient import TestClient
from altanalyze3.components.cellHarmony.webapp import biology_summary as summary
from altanalyze3.components.cellHarmony.flask.job_manager import JobStore


@pytest.fixture(autouse=True)
def no_external_literature(monkeypatch):
    from altanalyze3.components.cellHarmony.webapp import biology_literature
    monkeypatch.setattr(biology_literature,'search',lambda *a: {'status':'test_fixture','papers':[],'note':'Fixture literature transport'})


def fixture(tmp_path):
    store=JobStore(tmp_path/'jobs'); job=store.create_job('human','hs_lung_cellref2_reference',None,[])['job_id']
    p=store.outputs_dir(job);p.mkdir(exist_ok=True)
    obs=pd.DataFrame({'supervised_state':['AT2']*4+['Unaligned']*4,
                      'unsupervised_state':['AT2_c1']*4+['Fibroblast_c2']*4,
                      'disease':pd.Categorical(['IPF','IPF','Healthy','Healthy']*2),
                      'Library':pd.Categorical(['D1','D1','D2','D2']*2), 'age':np.arange(8)+30},index=[f'cell{i}' for i in range(8)])
    a=ad.AnnData(np.ones((8,2),dtype='f4'),obs=obs,var=pd.DataFrame(index=['SFTPC','ACTA2']))
    a.write_h5ad(p/'rna.h5ad')
    # Delete expression to enforce metadata-only reading; RNA values cannot be read.
    with h5py.File(p/'rna.h5ad','a') as f: del f['X']
    mods=['rna','adt','lipids','metabolite','grn','grn_tf']; markers={}; runs={}; manifest={'rna':{'h5ad':str(p/'rna.h5ad')}}
    features=['SFTPC','CD36','PC(16:0/18:1)','Citrate','STAT3|COL1A1','STAT3']
    for mod,feature in zip(mods,features):
        marker=p/f'{mod}_markers.tsv'; pd.DataFrame([{'Gene':feature,'Fold':3,'cluster':'AT2','rho':.8}]).to_csv(marker,sep='\t',index=False)
        markers[mod]={'markers_tsv':str(marker)}
        detail=p/f'{mod}_de.tsv'; pd.DataFrame([{'population':'AT2','gene':feature,'log2fc':2,'fdr':.01,'pval':.001},
                       {'population':'Unaligned','gene':feature,'log2fc':-2,'fdr':.02,'pval':.002},
                       {'population':'AT2','gene':'nonsignificant','log2fc':4,'fdr':.8,'pval':.5}]).to_csv(detail,sep='\t',index=False)
        runs[f'run_{mod}']={'status':'completed','run_id':f'run_{mod}','case_label':'IPF','control_label':'Healthy',
            'config':{'modality':mod,'population_col':'supervised_state','sample_field':'disease','comparison_type':'pseudobulk',
                      'group1_samples':['IPF'],'group2_samples':['Healthy']},'artifacts':{'DEG_detailed_fixture':str(detail)}}
    go=p/'go.tsv';pd.DataFrame([{'population':'AT2','term_name':'Surfactant metabolism','fdr':.001,'overlap_genes':'SFTPC'}]).to_csv(go,sep='\t',index=False)
    runs['run_rna']['artifacts']['goelite_tsv']=str(go)
    comm=p/'comm.tsv'; pd.DataFrame([{'sender_state':'AT2','receiver_state':'Unaligned','ligand':'TGFB1','receptor':'TGFBR2',
                                     'fastcomm_score':.7,'fastcomm_empirical_p':.01}]).to_csv(comm,sep='\t',index=False)
    aggregate=ad.AnnData(np.array([[.4,.9]],dtype='f4'),obs=pd.DataFrame({'supervised_state':['AT2'],'sample':['D1']},index=['D1_AT2']),
                        var=pd.DataFrame(index=['STAT3|COL1A1','TEAD2|CTGF']))
    aggregate.write_h5ad(p/'grn_network.h5ad');manifest['grn']={'network_h5ad':str(p/'grn_network.h5ad')}
    store.update_job(job,status='completed',analysis_mode='both',cluster_key='supervised_state',
        artifacts={'combined_h5ad':str(p/'rna.h5ad')},modality_artifacts=manifest,
        modalities={'available':[{'id':m,'label':m} for m in mods]},marker_analysis_by_modality=markers,
        cell_state_layers={'default':'supervised_state','layers':[
            {'key':'supervised_state','marker_analysis':markers['rna'],'marker_analysis_by_modality':markers,
             'fastcomm_analysis':{'scores_tsv':str(comm)}}, {'key':'unsupervised_state'}]},
        differential_history=runs)
    return store,job


def test_all_modalities_arms_layers_counts_and_directions_without_expression(tmp_path, monkeypatch):
    store,job=fixture(tmp_path)
    monkeypatch.setattr(ad,'read_h5ad',lambda *a,**k:pytest.fail('Full AnnData loaded'))
    r=summary.build_summary(store.get_job(job));d=r['digest']
    assert d['observations']['cells']==8
    assert d['observations']['layers']['supervised_state']['counts']=={'AT2':4,'Unaligned':4}
    assert d['observations']['covariates']['disease']['counts']=={'IPF':4,'Healthy':4}
    assert len(d['comparisons'])==6
    assert {h['modality'] for h in r['hypotheses']}==set(d['modalities'])
    for e in d['evidence']:
        if e['kind']=='differential_patterns':
            assert e['facts']['directions_by_state']=={'AT2':{'up':1},'Unaligned':{'down':1}}
            assert e['facts']['scanned_rows']==3
            assert e['action']['run_id']==f"run_{e['modality']}"
            assert e['action']['layer']=='supervised_state'
    assert any(e['kind']=='regulatory_prediction' and e['feature']=='TEAD2|CTGF' for e in d['evidence'])
    assert any(e['kind']=='pathway' and e['action']['mode']=='go' for e in d['evidence'])
    assert any(e['kind']=='communication' and e['action']['mode']=='fastcomm_network' for e in d['evidence'])
    assert r['model']['used'] is False


def test_background_has_citations_species_and_no_non_lung_transfer():
    bg=summary.state_background(['AT2','Unaligned','AT2_c1'],'human','CellRef Lung')
    assert {e['state'] for e in bg['matched']}=={'AT2','AT2_c1'}
    assert bg['unmatched']==['Unaligned']
    assert any(e['card']['literature'] for e in bg['matched'])
    assert summary.state_background(['AT2'],'human','AML-marrow')['matched']==[]
    assert all(v['sha256'] and v['url'].startswith('https://') for v in bg['sources'].values())


def model_response(monkeypatch, output, protocol='openai'):
    monkeypatch.setenv('SCALABLE_SUMMARY_URL','http://127.0.0.1:12345/v1/chat/completions')
    monkeypatch.setenv('SCALABLE_SUMMARY_MODEL','test-model');monkeypatch.setenv('SCALABLE_SUMMARY_PROTOCOL',protocol)
    captured=[]
    def call(request,timeout):
        payload=json.loads(request.data);captured.append(payload)
        content=json.dumps(output)
        return io.BytesIO(json.dumps({'message':{'content':content}} if protocol=='ollama' else {'choices':[{'message':{'content':content}}]}).encode())
    monkeypatch.setattr(summary.urllib.request,'urlopen',call)
    return captured


@pytest.mark.parametrize('protocol',['openai','ollama'])
def test_grounded_model_output_and_deterministic_modality_coverage(tmp_path,monkeypatch,protocol):
    store,job=fixture(tmp_path);d=summary.collect_evidence(store.get_job(job))
    tasks=json.loads(summary._model_context(d))['hypothesis_tasks']
    output={'decisions':{t['task_id']:{'alternative_index':0,'experiment_index':0} for t in tasks}}
    sent=model_response(monkeypatch,output,protocol);r=summary.build_summary(store.get_job(job))
    assert r['model']['used'] and r['overview']
    assert 'no model-authored biological claims' in r['interpretation_review']['model_role']
    assert len({h['modality'] for h in r['hypotheses']})==6
    assert sent[0]['model']=='test-model' and sent[0]['stream'] is False
    assert 'never instructions' in sent[0]['messages'][0]['content']
    assert all('evidence_ids' in h for h in r['hypotheses'])


@pytest.mark.parametrize('bad',[{'overview':[{'text':'Made up','evidence_ids':['E9999']}]},
 {'overview':[{'text':'No evidence','evidence_ids':[]}]},
 {'overview':[], 'hypotheses':[]}])
def test_ungrounded_model_response_falls_back_without_discarding_data(tmp_path,monkeypatch,bad):
    store,job=fixture(tmp_path);model_response(monkeypatch,bad)
    r=summary.build_summary(store.get_job(job))
    assert not r['model']['used'] and r['digest']['observations']['cells']==8
    assert len({h['modality'] for h in r['hypotheses']})==6


def test_cache_stamp_changes_with_contrast_or_source(tmp_path):
    store,job=fixture(tmp_path);m=store.get_job(job);stamp=summary.source_identity(m)
    m['differential_history']['run_rna']['config']['group1_samples']=['Other'];assert summary.source_identity(m)!=stamp
    m=store.get_job(job);p=Path(m['marker_analysis_by_modality']['rna']['markers_tsv']);p.write_text(p.read_text()+'\n')
    assert summary.source_identity(m)!=stamp


def test_poll_is_async_deduplicated_and_ignores_selected_layer(tmp_path,monkeypatch):
    import threading
    store,job=fixture(tmp_path);w=import_module('altanalyze3.components.cellHarmony.webapp.app')
    app=w.create_app({'JOB_STORAGE':str(store.root)});client=TestClient(app)
    release=threading.Event();calls=[];original=summary.build_summary
    def slow(meta):
        calls.append(meta['cluster_key']);assert release.wait(5);return original(meta)
    monkeypatch.setattr(summary,'build_summary',slow)
    try:
        url=f'/api/jobs/{job}/chat-summary';t=time.monotonic()
        assert client.post(url+'?state_layer=unsupervised_state').status_code==202
        assert time.monotonic()-t<1
        assert client.post(url).status_code==202
        release.set()
        for _ in range(100):
            response=client.post(url)
            if response.status_code!=202:break
            time.sleep(.01)
        assert response.status_code==200,response.text
        assert calls==['supervised_state']
        assert client.post(url).content==response.content
        store.update_job(job,status='processing');assert client.post(url).status_code==409
        assert client.post('/api/jobs/missing/chat-summary').status_code==404
    finally:
        release.set();app.state.summary_build_queue.close();app.state.plot_build_queue.close()


def test_zero_differentials_no_invented_disease_effect(tmp_path):
    store,job=fixture(tmp_path);store.update_job(job,differential_history={})
    r=summary.build_summary(store.get_job(job))
    assert r['digest']['comparisons']==[]
    assert len([e for e in r['digest']['evidence'] if e['kind'].startswith('differential')])==0
    assert any('no completed covariate comparison' in gap for gap in r['digest']['coverage'])


def test_chunked_scan_retains_rare_and_unaligned_cells(tmp_path,monkeypatch):
    store,job=fixture(tmp_path);monkeypatch.setattr(summary,'CHUNK',3)
    d=summary.collect_evidence(store.get_job(job))
    assert sum(d['observations']['layers']['supervised_state']['counts'].values())==8
    assert d['observations']['layers']['supervised_state']['counts']['Unaligned']==4


def test_missing_categorical_and_nullable_developmental_covariate(tmp_path):
    store,job=fixture(tmp_path);m=store.get_job(job)
    path=m['artifacts']['combined_h5ad']
    from anndata.io import write_elem
    obs=pd.DataFrame({'supervised_state':pd.Categorical(['AT2',None]),'developmental_age':pd.array([35,None],dtype='Int64'),
                      'missing_annotation':pd.Categorical([None,None])},index=['a','b'])
    with h5py.File(path,'a') as f:del f['obs'];write_elem(f,'obs',obs)
    d=summary.collect_evidence(m)
    assert d['observations']['cells']==2
    assert d['observations']['layers']['supervised_state']['counts']=={'AT2':1,'(missing)':1}
    assert d['observations']['covariates']['developmental_age']['mean']==35
    assert d['observations']['covariates']['missing_annotation']['counts']=={'(missing)':2}


def test_llm_rejects_cross_covariate_evidence_mix(tmp_path,monkeypatch):
    store,job=fixture(tmp_path);m=store.get_job(job)
    m['differential_history']['run_grn']['config']['group1_samples']=['Male']
    d=summary.collect_evidence(m)
    rna=next(e for e in d['evidence'] if e['kind']=='differential_feature' and e['modality']=='rna')
    grn=next(e for e in d['evidence'] if e['kind']=='differential_feature' and e['modality']=='grn')
    model_response(monkeypatch,{'overview':[{'text':'Saved RNA contrast','evidence_ids':[rna['id']]}],
        'hypotheses':[{'modality':'rna','title':'Invalid cross-covariate link','rationale':'','validation':'',
                      'evidence_ids':[rna['id'],grn['id']]}]})
    # Exercise the comparison gate independently of context prioritization.
    monkeypatch.setattr(summary,'_model_context',lambda digest:json.dumps({'evidence':digest['evidence']}))
    with pytest.raises(ValueError,match='comparison arms'):summary.interpret(d)


def test_model_unavailable_preserves_report_and_never_echoes_credentials(tmp_path,monkeypatch):
    store,job=fixture(tmp_path)
    monkeypatch.setenv('SCALABLE_SUMMARY_URL','http://127.0.0.1:1/v1/chat/completions')
    monkeypatch.setenv('SCALABLE_SUMMARY_MODEL','no-model')
    monkeypatch.setenv('SCALABLE_SUMMARY_API_KEY','PRIVATE_SECRET')
    def unavailable(*a,**k):raise RuntimeError('PRIVATE_SECRET')
    monkeypatch.setattr(summary.urllib.request,'urlopen',unavailable)
    r=summary.build_summary(store.get_job(job))
    assert r['model']['used'] is False
    assert 'PRIVATE_SECRET' not in json.dumps(r)
    assert r['digest']['observations']['cells']==8


def test_empty_optional_tables_preserve_dataset_and_modality_inventory(tmp_path):
    store,job=fixture(tmp_path);m=store.get_job(job)
    Path(m['differential_history']['run_rna']['artifacts']['DEG_detailed_fixture']).write_text('')
    Path(m['differential_history']['run_rna']['artifacts']['goelite_tsv']).write_text('')
    r=summary.build_summary(m)
    assert r['digest']['observations']['cells']==8
    assert len(r['digest']['comparisons'])==6
    assert any(h['modality']=='rna' and ' in AT2' in h['title'] for h in r['hypotheses'])


def test_llm_cannot_cite_records_not_sent_or_mix_state_layers(tmp_path,monkeypatch):
    store,job=fixture(tmp_path);d=summary.collect_evidence(store.get_job(job))
    marker=next(e for e in d['evidence'] if e['kind']=='marker' and e['modality']=='rna')
    extras=[dict(marker,id=f'E{900+i}') for i in range(7)]
    d['evidence'].extend(extras)
    context=json.loads(summary._model_context(d));context.pop('hypothesis_tasks',None)
    monkeypatch.setattr(summary,'_model_context',lambda digest:json.dumps(context))
    model_response(monkeypatch,{'overview':[{'text':'Unsupported citation','evidence_ids':[extras[-1]['id']]}]})
    with pytest.raises(ValueError,match='unsupported evidence'):summary.interpret(d)
    other=dict(marker,id='E1000',layer='unsupervised_state');d['evidence'].append(other)
    model_response(monkeypatch,{'overview':[{'text':'Saved marker evidence','evidence_ids':[marker['id']]}],
        'hypotheses':[{'modality':'rna','title':'Layer mismatch','evidence_ids':[marker['id'],other['id']]}]})
    monkeypatch.setattr(summary,'_model_context',lambda digest:json.dumps({'evidence':digest['evidence']}))
    with pytest.raises(ValueError,match='cell-state layers'):summary.interpret(d)


def test_bundle_store_protocol_without_source_h5ad_returns_explicit_coverage(tmp_path):
    store,job=fixture(tmp_path);m=store.get_job(job)
    genes=tmp_path/'genes.tsv';genes.write_text('SFTPC\nACTA2\n')
    m['artifacts']={};m['modality_artifacts']['rna']['h5ad']=str(genes)
    class BundleStore:
        def job_exists(self,key):return key==job
        def get_job(self,key):return m
    w=import_module('altanalyze3.components.cellHarmony.webapp.app')
    app=w.create_app({'JOB_STORAGE':str(store.root)});app.state.job_store=BundleStore()
    with TestClient(app) as client:
        for _ in range(100):
            r=client.post(f'/api/jobs/{job}/chat-summary')
            if r.status_code!=202:break
            time.sleep(.01)
        assert r.status_code==200,r.text
        assert r.json()['digest']['observations']['cells'] is None
        assert any('observation metadata unavailable' in s for s in r.json()['digest']['coverage'])
        assert len(r.json()['digest']['comparisons'])==6


def test_synthesis_keeps_untested_communication_and_links_exact_regulatory_chain(tmp_path):
    store,job=fixture(tmp_path);meta=store.get_job(job);out=store.outputs_dir(job)
    pd.DataFrame([
        {'population':'AT2','gene':'SFTPC','log2fc':-1,'fdr':.01},
        {'population':'AT2','gene':'COL1A1','log2fc':2,'fdr':.04},
    ]).to_csv(out/'rna_de.tsv',sep='\t',index=False)
    go=out/'go.tsv';pd.DataFrame([
        {'population':'AT2__down','term_name':'cholesterol biosynthetic process','fdr':.001,'overlap_genes':'SFTPC'},
        {'population':'coreg_global','term_name':'Global genes','fdr':.001,'overlap_genes':'COL1A1'},
    ]).to_csv(go,sep='\t',index=False)
    pd.DataFrame([{'population':'AT2','gene':'CE(22:4)','log2fc':2,'fdr':.01}]).to_csv(out/'lipids_de.tsv',sep='\t',index=False)
    comm=out/'contrast.tsv';pd.DataFrame([
        {'population':'AT2','gene':'Fibroblast->AT2:FN1->CD44','delta_score':-.6,'log2fc':-20,'fdr':np.nan,'test_applied':False},
        {'population':'AT2','gene':'Basal->AT2:APP->CD74','delta_score':.4,'log2fc':20,'fdr':np.nan,'test_applied':False},
    ]).to_csv(comm,sep='\t',index=False)
    run=dict(meta['differential_history']['run_rna']);run['config']=dict(run['config'],modality='cell_communication')
    run['artifacts']={'DEG_detailed_fixture':str(comm)}
    meta['differential_history']['run_comm']=run
    r=summary.build_summary(meta);p=next(p for p in r['digest']['synthesis']['programs'] if p['state']=='AT2')
    assert p['regulatory_chains'][0]['regulator']=='STAT3'
    assert p['regulatory_chains'][0]['target']=='COL1A1'
    assert p['regulatory_chains'][0]['TF_prediction_change']=='higher'
    assert 'coreg_global' not in {p['state'] for p in r['digest']['synthesis']['programs']}
    effects=[e for e in r['digest']['evidence'] if e['kind']=='communication_contrast']
    assert len(effects)==2 and all(e['facts']['test_applied'] is False for e in effects)
    assert effects[0]['facts']['delta_score']==-.6
    assert any('discordance' in i['text'] for i in p['integrative_insights'])
    assert any('opposite directions' in i['text'] for i in p['integrative_insights'])
    assert any(h['modality']=='cell_communication' and 'run_comm' in
               {e['run_id'] for e in r['digest']['evidence'] if e['id'] in h['evidence_ids']} for h in r['hypotheses'])
    assert r['digest']['synthesis']['design_limits']


def test_synthesis_does_not_join_regulatory_predictions_from_different_arms(tmp_path):
    store,job=fixture(tmp_path);meta=store.get_job(job);out=store.outputs_dir(job)
    pd.DataFrame([{'population':'AT2','gene':'COL1A1','log2fc':2,'fdr':.01}]).to_csv(out/'rna_de.tsv',sep='\t',index=False)
    meta['differential_history']['run_grn']['config']['group1_samples']=['Other']
    d=summary.collect_evidence(meta)
    assert all(not p['regulatory_chains'] for p in d['synthesis']['programs'])


def test_differential_summary_ranking_preserves_smaller_states(tmp_path):
    store,job=fixture(tmp_path);meta=store.get_job(job);out=store.outputs_dir(job)
    pd.DataFrame([{'population':'AT2','gene':f'G{i}','log2fc':10,'fdr':1e-9} for i in range(40)]+
                 [{'population':'Unaligned','gene':'rare_state_feature','log2fc':1,'fdr':.04}]).to_csv(out/'rna_de.tsv',sep='\t',index=False)
    d=summary.collect_evidence(meta)
    assert any(e['feature']=='rare_state_feature' for e in d['evidence'])
    assert next(e for e in d['evidence'] if e['kind']=='differential_patterns' and e['modality']=='rna')['facts']['directions_by_state']['AT2']['up']==40


def test_compact_model_context_keeps_native_directions_and_all_state_inventory(tmp_path):
    store,job=fixture(tmp_path);d=summary.collect_evidence(store.get_job(job))
    c=json.loads(summary._model_context(d))
    assert c['observations']['layers']['supervised_state']['counts']=={'AT2':4,'Unaligned':4}
    assert {p['modality'] for p in c['global_patterns']}==set(d['modalities'])
    rna=next(e for e in c['evidence'] if e['modality']=='rna' and e['kind']=='differential_feature')
    assert rna['direction_in_case']=='higher'
    assert 'log2fc' not in rna['facts'] and 'action' not in rna
    assert len(c['evidence'])<=len(d['evidence'])
    assert {e['id'] for e in c['evidence']}<={e['id'] for e in d['evidence']}


def test_truncated_live_model_output_is_not_presented_as_completed_interpretation(tmp_path,monkeypatch):
    store,job=fixture(tmp_path)
    monkeypatch.setenv('SCALABLE_SUMMARY_URL','http://127.0.0.1:12345/v1/chat/completions')
    monkeypatch.setenv('SCALABLE_SUMMARY_MODEL','test-model')
    output={'choices':[{'finish_reason':'length','message':{'content':'{}'}}]}
    monkeypatch.setattr(summary.urllib.request,'urlopen',lambda *a,**k:io.BytesIO(json.dumps(output).encode()))
    r=summary.build_summary(store.get_job(job))
    assert not r['model']['used']
    assert 'Model interpretation rejected:' in r['model']['message']
    assert r['digest']['observations']['cells']==8


def test_question_framework_records_design_and_requires_pathway_target_membership(tmp_path):
    store,job=fixture(tmp_path);d=summary.collect_evidence(store.get_job(job));f=d['question_framework']
    assert len(f['analysis_plan'])==1
    assert set(f['analysis_plan'][0]['completed_modalities'])==set(d['modalities'])
    assert f['analysis_plan'][0]['missing_modality_comparisons']==[]
    q=next(q for q in f['questions'] if q['state']=='AT2')
    assert q['domain']=='surfactant' and q['RNA_anchor_ids']
    # STAT3->COL1A1 has no overlap with the saved SFTPC pathway: co-occurrence
    # alone must not be promoted to a regulator of surfactant metabolism.
    assert q['regulatory_chains']==[] and q['replication_established'] is False
    assert any(r['state']=='Unaligned' for r in f['state_coverage'])
    assert not next(r for r in f['state_coverage'] if r['state']=='Unaligned')['RNA_functional_question_available']


def test_model_cannot_invent_a_biological_question_even_with_valid_evidence_ids(tmp_path,monkeypatch):
    store,job=fixture(tmp_path);d=summary.collect_evidence(store.get_job(job))
    e=next(e for e in d['evidence'] if e['modality']=='rna' and e['kind']=='differential_feature')
    model_response(monkeypatch,{'overview':[{'text':'RNA observation','evidence_ids':[e['id']]}],
        'hypotheses':[{'modality':'rna','question_id':'invented','title':'A new mechanism',
            'rationale':'Unsupported question','validation':'Experiment','evidence_ids':[e['id']]}]})
    r=summary.build_summary(store.get_job(job));assert not r['model']['used']
    assert 'Model interpretation rejected:' in r['model']['message']


def test_task_map_keeps_code_owned_scope_and_requires_every_modality(tmp_path,monkeypatch):
    store,job=fixture(tmp_path);d=summary.collect_evidence(store.get_job(job))
    tasks=json.loads(summary._model_context(d))['hypothesis_tasks']
    assert len(tasks)==6 and {t['modality'] for t in tasks}==set(d['modalities'])
    output={'decisions':{t['task_id']:{'alternative_index':1,'experiment_index':2,
        'title':'Invented title','evidence_ids':['E9999']} for t in tasks}}
    model_response(monkeypatch,output);narrative,_=summary.interpret(d)
    assert {h['modality'] for h in narrative['hypotheses']}==set(d['modalities'])
    assert all(h['title']!='Invented title' and 'E9999' not in h['evidence_ids'] for h in narrative['hypotheses'])
    assert all(h['validation']==t['experiment_options'][2] and h['priority_alternative']==t['alternatives'][1]
               for h,t in zip(narrative['hypotheses'],tasks))
    del output['decisions'][tasks[-1]['task_id']]
    with pytest.raises(ValueError,match='did not complete'):summary.interpret(d)


@pytest.mark.parametrize('text,modality',[
 ('BOS alters the AM state.','rna'),
 ('NR4A3 is predicted to repress CEBPB.','grn'),
 ('Predicted CD326 indicates increased antigen presentation.','adt'),
 ('Reduced FN1 signaling changes AM.','cell_communication'),
 ('CD169 is higher on macrophages.','adt')])
def test_semantic_guards_reject_observed_failure_modes(text,modality):
    with pytest.raises(ValueError):summary._check_interpretation_scope(text,modality)


def test_closed_mode_rejects_free_prose_without_discarding_computed_evidence(tmp_path,monkeypatch):
    store,job=fixture(tmp_path);d=summary.collect_evidence(store.get_job(job))
    tasks=json.loads(summary._model_context(d))['hypothesis_tasks']
    output={'overview':[{'text':'An invented causal explanation','task_ids':[tasks[0]['task_id']]}],
        'hypotheses':{t['task_id']:{'rationale':'STAT3 is predicted to repress COL1A1.',
            'validation':'Directly measure the phenotype.'} for t in tasks}}
    model_response(monkeypatch,output);report=summary.build_summary(store.get_job(job))
    assert not report['model']['used']
    assert 'required closed experiment-priority decisions' in report['model']['message']
    assert report['digest']['observations']['cells']==8
    assert len(report['hypotheses'])==6
    assert all('predicted to repress' not in h['rationale'] for h in report['hypotheses'])


@pytest.mark.parametrize('index',[3,-1,'0',True,None])
def test_closed_mode_rejects_out_of_contract_experiment_choices(tmp_path,monkeypatch,index):
    store,job=fixture(tmp_path);d=summary.collect_evidence(store.get_job(job))
    tasks=json.loads(summary._model_context(d))['hypothesis_tasks']
    output={'decisions':{t['task_id']:{'alternative_index':index,'experiment_index':0} for t in tasks}}
    model_response(monkeypatch,output)
    with pytest.raises(ValueError,match='outside the defined alternatives'):summary.interpret(d)


def test_question_overview_leads_with_named_relationships_and_original_comparison(tmp_path):
    from altanalyze3.components.cellHarmony.webapp.biology_questions import overview
    store,job=fixture(tmp_path);d=summary.collect_evidence(store.get_job(job))
    path=next(e for e in d['evidence'] if e['kind']=='pathway')
    f=d['question_framework'];f['cross_state_domain_questions']=[{'contrast':d['comparisons'][0]['contrast'],
        'opposing_directions':True,'state_processes':[{'state':'AT2','direction':'up','evidence_id':path['id']},
        {'state':'Unaligned','direction':'down','evidence_id':path['id']}]}]
    rows=overview(d)
    assert rows[0]['text'].startswith('IPF vs Healthy: AT2: Surfactant metabolism')
    assert 'Unaligned' in rows[0]['text'] and 'cell-to-cell coupling' in rows[0]['text']
    assert rows[0]['evidence_ids']==[path['id']]
    assert f['questions'][0]['background_role'].startswith('Established state function')
