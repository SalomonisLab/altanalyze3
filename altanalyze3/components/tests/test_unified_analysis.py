"""Independent barcode/expression fixtures: no scientific models are run here."""
from types import SimpleNamespace
from pathlib import Path
from importlib import import_module

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from fastapi.testclient import TestClient

from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
from altanalyze3.components.cellHarmony.webapp.analysis_workflow import _combine, association_payload, enable_differential
from altanalyze3.components.cellHarmony.flask import pipeline as pipeline_mod
from altanalyze3.components.cellHarmony.webapp.state_layers import ACTIVE_LAYER, AnalysisJobStore
from altanalyze3.components.cellHarmony.webapp.app import create_app
web_mod = import_module('altanalyze3.components.cellHarmony.webapp.app')


@pytest.fixture()
def combined_fixture(tmp_path):
    store = JobStore(tmp_path / 'jobs')
    job = store.create_job('human', 'fixture', None, [])['job_id']
    # Independently defined baseline: every source feature and row has known values.
    baseline = ad.AnnData(sp.csr_matrix(np.array([[1.,2.,3.],[4.,5.,6.],[7.,8.,9.],[10.,11.,12.]], dtype=np.float32)),
                         obs=pd.DataFrame({'Library':['a','a','b','b']}, index=['s1','s2','s3','neither']),
                         var=pd.DataFrame(index=['gA','gB','gC']))
    baseline.layers['counts'] = baseline.X.copy()
    path = tmp_path / 'qc.h5ad'; baseline.write_h5ad(path)
    return store, job, baseline, path


def branch(tmp_path, baseline, names, key, values, unsup=False):
    matrix = baseline[names].copy(); matrix.obs[key] = values
    matrix.obsm['X_umap'] = np.arange(2*len(names),dtype=np.float32).reshape(-1,2)
    if not unsup:
        matrix.obsm['X_fixture_modality'] = np.arange(2*len(names),dtype=np.float64).reshape(-1,2) + 0.123456789
        matrix.uns['imputed_modalities'] = {'fixture': {'obsm_key':'X_fixture_modality','feature_names_key':'fixture_names'}}
        matrix.uns['fixture_names'] = ['measured_output_A','measured_output_B']
    if unsup: matrix.obs['cluster'] = ['c1']*len(names)
    path = tmp_path / (key+'.h5ad'); matrix.write_h5ad(path)
    markers = {'enabled':False}; comm = {'enabled':False}
    meta = {'cluster_key':key,'artifacts':{'combined_h5ad':str(path),'marker_genes_zip':'fixture.zip'},
            'marker_analysis':markers,'marker_analysis_by_modality':{'rna':markers},
            'fastcomm_analysis':comm,'modalities':{'available':[{'id':'rna'}]},'modality_artifacts':{'rna':{'h5ad':str(path)}},
            'icgs3_analysis':{},'cell_state_layers':{'layers':[{}, {'marker_analysis':markers,'fastcomm_analysis':comm}]}}
    return SimpleNamespace(get_job=lambda job:meta)


def test_union_retains_sup_only_and_unsup_only_with_exact_source_values(combined_fixture,tmp_path):
    store,job,base,path = combined_fixture
    sup = branch(tmp_path,base,['s1','s2'],'ref',['A','B'])
    unsup = branch(tmp_path,base,['s2','s3'],'cell_state_predicted',['c1','c1'],True)
    meta,output = _combine(job,store,path,base.obs_names,base.var_names,sup,unsup,'lzf')
    result=ad.read_h5ad(output)
    assert list(result.obs_names)==['s1','s2','s3'] # union, not intersection or concatenation duplicates
    assert result.var_names.equals(base.var_names)
    np.testing.assert_array_equal(result.X.toarray(),base.X[:3].toarray())
    np.testing.assert_array_equal(result.layers['counts'].toarray(),base.layers['counts'][:3].toarray())
    assert list(result.obs['supervised_state'].astype(str))==['A','B','Unaligned']
    assert list(result.obs['unsupervised_cluster'].astype(str))==['Not clustered','c1','c1']
    assert np.isnan(result.obsm['X_umap_unsupervised'][0]).all()
    assert np.isnan(result.obsm['X_umap_supervised'][2]).all()
    np.testing.assert_array_equal(result.obsm['X_fixture_modality'][:2],np.arange(4,dtype=np.float64).reshape(-1,2)+0.123456789)
    assert np.isnan(result.obsm['X_fixture_modality'][2]).all()
    assert list(result.uns['fixture_names']) == ['measured_output_A','measured_output_B']
    payload=association_payload(result.obs)
    assert sum(row['cells'] for row in payload['rows'])==3
    assert 'Unaligned' in payload['states'] and 'Not clustered' in payload['clusters']
    filtered=association_payload(result.obs,[('Library',['b'])])
    assert filtered['cells']==1 and filtered['rows'][0]['percent']==100


def test_empty_supervised_roster_keeps_every_clustered_cell(combined_fixture,tmp_path):
    store,job,base,path=combined_fixture
    sup=branch(tmp_path,base,[],'ref',[])
    unsup=branch(tmp_path,base,['s1','s2','s3'],'cell_state_predicted',['c1']*3,True)
    _,output=_combine(job,store,path,base.obs_names,base.var_names,sup,unsup,'lzf')
    result=ad.read_h5ad(output)
    assert len(result)==3 and set(result.obs.supervised_state.astype(str))=={'Unaligned'}


def test_changed_feature_roster_is_rejected_before_union_write(combined_fixture,tmp_path):
    store,job,base,path=combined_fixture
    sup=branch(tmp_path,base,['s1'],'ref',['A'])
    bad=base[:,:2].copy()
    unsup=branch(tmp_path,bad,['s2'],'cell_state_predicted',['c1'],True)
    with pytest.raises(ValueError,match='feature roster'):
        _combine(job,store,path,base.obs_names,base.var_names,sup,unsup,'lzf')
    assert not (store.outputs_dir(job)/'combined_with_umap_and_markers.h5ad').exists()


def test_layer_selection_is_request_local_and_does_not_mutate_disk(tmp_path):
    store=AnalysisJobStore(tmp_path)
    job=store.create_job('human','fixture',None,[])['job_id']
    store.update_job(job,cluster_key='unsupervised_cluster',marker_analysis={'tag':'unsup'},cell_state_layers={
        'default':'unsupervised_cluster','layers':[{'key':'supervised_state','marker_analysis':{'tag':'sup'}}]})
    token=ACTIVE_LAYER.set('supervised_state')
    try:
        assert store.get_job(job)['marker_analysis']['tag']=='sup'
        assert store.get_job(job)['cluster_key']=='supervised_state'
        store.update_job(job,message='independent update')
    finally: ACTIVE_LAYER.reset(token)
    assert store.get_job(job)['cluster_key']=='unsupervised_cluster'
    assert store.get_job(job)['marker_analysis']['tag']=='unsup'


def test_invalid_mode_and_running_configuration_changes_are_rejected(tmp_path):
    app=create_app({'JOB_STORAGE':str(tmp_path),'ISOLATE_JOBS':False})
    job=app.state.job_store.create_job('human','fixture',None,[])['job_id']
    with TestClient(app) as client:
        assert client.post(f'/api/jobs/{job}/qc',json={'analysis_mode':'silently_fallback'}).status_code==422
        app.state.job_store.update_job(job,status='processing')
        assert client.post(f'/api/jobs/{job}/qc',json={'analysis_mode':'both'}).status_code==409
        assert client.post(f'/api/jobs/{job}/run').status_code==409
    app.state.job_runner.executor.shutdown(wait=True)


def test_h5ad_differential_keeps_original_population_annotations(combined_fixture):
    store,job,base,path=combined_fixture
    base.obs['original_state']=['A','A','B','B']
    base.obs['cluster']=['c1','c2','c1','c2']
    base.write_h5ad(path)
    meta={'cluster_key':'cluster','modalities':{'available':[{'id':'rna'}]}}
    original={'files':[{'filename':'original.h5ad','sample_name':'original'}]}
    enable_differential(meta,original,path,['cluster'])
    assert {'original_state','cluster'} <= {v['value'] for v in meta['differential_options']['population_columns']}
    assert meta['differential_options']['default_population_col']=='cluster'


def test_communication_worker_selects_explicit_population_branch(tmp_path,monkeypatch):
    store=JobStore(tmp_path)
    job=store.create_job('human','fixture',None,[])['job_id']
    store.update_job(job,cluster_key='unsupervised_cluster',fastcomm_analysis={'tag':'unsupervised'},
        cell_state_layers={'default':'unsupervised_cluster','layers':[
            {'key':'supervised_state','fastcomm_analysis':{'tag':'supervised'}}]},
        differential={'config':{'modality':'cell_communication','population_col':'supervised_state',
            'group1_samples':['a'],'group2_samples':['b'],'sample_field':'Library'}})
    monkeypatch.setattr(pipeline_mod,'_run_cell_communication_differential',lambda **kw: kw['meta']['fastcomm_analysis'])
    assert pipeline_mod.run_cellharmony_differential(job,store)=={'tag':'supervised'}
    assert store.get_job(job)['fastcomm_analysis']=={'tag':'unsupervised'}


def test_association_endpoint_counts_union_and_recomputes_filtered_denominator(tmp_path):
    app=create_app({'JOB_STORAGE':str(tmp_path/'jobs'),'ISOLATE_JOBS':False})
    store=app.state.job_store
    job=store.create_job('human','fixture',None,[])['job_id']
    matrix=ad.AnnData(obs=pd.DataFrame({
        'Library':['a','a','b'],
        'unsupervised_cluster':['c1','c1','c2'],
        'supervised_state':['A','Unaligned','A']},index=['s1','s2','s3']))
    path=tmp_path/'fixture.h5ad';matrix.write_h5ad(path)
    store.update_job(job,status='completed',analysis_mode='both',artifacts={'combined_h5ad':str(path)})
    with TestClient(app) as client:
        url=f'/api/jobs/{job}/cluster-associations'
        response=client.get(url)
        assert response.status_code==200
        assert sum(r['cells'] for r in response.json()['rows'])==3
        assert 'Unaligned' in response.json()['states']
        filtered=client.get(url,params={'subset_by':'Library','subset_values':'b'})
        assert filtered.json()['cells']==1 and filtered.json()['rows'][0]['percent']==100
        assert client.get(url,params={'subset_by':'unknown','subset_values':'b'}).status_code==400
        store.update_job(job,status='processing')
        assert client.get(url).status_code==404
    app.state.job_runner.executor.shutdown(wait=True)


@pytest.mark.parametrize('mode,coords,overlay',[
    ('supervised','',True),('unsupervised','',False),('both','',False),
    ('both','X_umap_supervised',True),('both','X_umap_unsupervised',False)])
def test_reference_overlay_only_uses_the_supervised_coordinate_space(monkeypatch,mode,coords,overlay):
    cache={'obs_names':np.array(['s1']), 'populations':np.array(['c1']),
        'sample_field':'Library','sample_labels':np.array(['a']), 'cluster_key':'cluster',
        'umap_x':np.array([0.]),'umap_y':np.array([1.]),
        'adata':SimpleNamespace(obsm={'X_umap_supervised':np.array([[2.,3.]]),
            'X_umap_unsupervised':np.array([[0.,1.]])})}
    called=[]
    monkeypatch.setattr(web_mod,'_get_expression_cache',lambda *a,**kw:cache)
    monkeypatch.setattr(web_mod,'_load_reference_adata',lambda *a:called.append(True))
    payload=web_mod._build_umap_payload(None,{'analysis_mode':mode},coords_key=coords)
    assert bool(called)==overlay
    assert payload['reference_hidden']==(not overlay)


@pytest.mark.parametrize('mode',['supervised','unsupervised','both'])
def test_upload_saves_analysis_choice_before_qc_and_recovers_missing_storage(tmp_path,mode):
    matrix=ad.AnnData(sp.csr_matrix(np.array([[1,2],[3,4]],dtype=np.int32)),
        obs=pd.DataFrame(index=['s1','s2']),var=pd.DataFrame(index=['gA','gB']))
    path=tmp_path/'counts.h5ad';matrix.write_h5ad(path)
    root=tmp_path/'jobs'
    app=create_app({'JOB_STORAGE':str(root),'ISOLATE_JOBS':False})
    root.rmdir() # reproduce the removed-directory error in the screenshot
    with TestClient(app) as client:
        with path.open('rb') as f:
            response=client.post('/api/jobs',data={'species':'human','reference':'fixture' if mode!='unsupervised' else '',
                'analysis_mode':mode,'sample_names':'library'},files={'files':('counts.h5ad',f,'application/octet-stream')})
        assert response.status_code==200,response.text
        meta=app.state.job_store.get_job(response.json()['job_id'])
        assert meta['analysis_mode']==mode and meta['qc']['analysis_mode']==mode
        assert meta['status']=='uploaded' # registration only, no scientific inference
        with path.open('rb') as f:
            bad=client.post('/api/jobs',data={'species':'human','reference':'fixture','analysis_mode':'invalid','sample_names':'library'},
                files={'files':('counts.h5ad',f,'application/octet-stream')})
        assert bad.status_code==422
        html=client.get('/').text
        assert html.index('id="job-form"') < html.index('id="analysis-mode"') < html.index('id="qc-form"')
        assert html.count('id="analysis-mode"')==1
    app.state.job_runner.executor.shutdown(wait=True)


@pytest.mark.parametrize('mode', ['unsupervised', 'both'])
def test_unified_icgs3_stage_logs_update_parent_progress(tmp_path, mode):
    from altanalyze3.components.cellHarmony.flask.tasks import _JobLogStream
    from altanalyze3.components.cellHarmony.webapp.analysis_workflow import BranchStore
    store = JobStore(tmp_path / 'jobs')
    job = store.create_job('human', 'fixture', None, [])['job_id']
    store.update_job(job, qc={'analysis_mode': mode})
    branch_store = BranchStore(store, job, 'unsupervised', tmp_path / 'prepared.h5ad')
    stream = _JobLogStream(branch_store, job)
    stream.write('[ICGS3] loading 1 input(s)\n')
    loading = store.get_job(job)
    assert 'step 1 of 10' in loading['message']
    stream.write('[ICGS3] running final UMAP on 619 features')
    stream.flush()  # ICGS3 Tee flushes partial lines before their newline.
    umap = store.get_job(job)
    assert umap['progress'] > loading['progress']
    assert 'step 10 of 10: UMAP' in umap['message']
    branch_store.append_log(job, '[ICGS3] UMAP transformed 50,000 remaining cells (50,000/100,000 completed;')
    projected = store.get_job(job)
    assert projected['progress'] >= umap['progress']
    assert '50,000 of 100,000' in projected['message']
    branch_store.append_log(job, '[ICGS3] loading repeated earlier stage')
    assert store.get_job(job)['progress'] == projected['progress']


def test_unified_workflow_routes_icgs3_stdout_into_live_parent_status(tmp_path, monkeypatch):
    workflow = import_module('altanalyze3.components.cellHarmony.webapp.analysis_workflow')
    discover = import_module('altanalyze3.components.cellHarmony.scalable_discover.pipeline')
    store = JobStore(tmp_path / 'jobs')
    job = store.create_job('human', 'fixture', None, [])['job_id']
    store.update_job(job, qc={'analysis_mode': 'unsupervised'})
    monkeypatch.setattr(workflow, '_shared_qc', lambda *a: (tmp_path / 'prepared.h5ad', [], []))
    def stage_fixture(*args, **kwargs):
        print('[ICGS3] running final UMAP on software test fixture', flush=True)
        raise RuntimeError('stop software fixture before scientific analysis')
    monkeypatch.setattr(discover, 'run_discover_pipeline', stage_fixture)
    with pytest.raises(RuntimeError, match='stop software fixture'):
        workflow.run_analysis_workflow(job, store, tmp_path / 'registry.json')
    assert 'step 10 of 10: UMAP' in store.get_job(job)['message']
    assert '[unsupervised] [ICGS3]' in (store.logs_dir(job) / 'pipeline.log').read_text()


@pytest.mark.parametrize('mode', ['unsupervised', 'both'])
def test_status_api_preserves_each_native_icgs3_step_after_shared_normalization(tmp_path, mode):
    import os
    from altanalyze3.components.cellHarmony.scalable_discover.tasks import STAGES
    app = create_app({'JOB_STORAGE': str(tmp_path), 'ISOLATE_JOBS': False})
    store = app.state.job_store
    job = store.create_job('human', 'icgs3', None, [])['job_id']
    # A live external owner keeps recovery from mistaking this status fixture for an orphan.
    store.update_job(job, status='processing', worker_pid=os.getppid(), analysis_mode=mode,
                     qc={'analysis_mode': mode}, cell_state_layers={})
    store.append_log(job, 'Normalization steps: 100%')
    with TestClient(app) as client:
        for fragment, progress, stage in STAGES:
            if not stage.startswith('ICGS3'):
                continue
            message = f'Unsupervised: {stage}'
            store.update_job(job, message=message, progress=progress)
            store.append_log(job, f'[unsupervised] [ICGS3] {fragment}')
            response = client.get(f'/api/jobs/{job}/status')
            assert response.status_code == 200
            assert response.json()['message'] == message
            assert response.json()['progress'] == progress
    app.state.job_runner.executor.shutdown(wait=True)


def test_status_message_keeps_supervised_log_updates_and_terminal_errors():
    derive = web_mod._derive_live_pipeline_message
    assert derive('processing', ['Normalization steps: 100%', 'Running rna2lipid lipid imputation.'],
                  'Preparing inputs…') == 'Running rna2lipid lipid imputation.'
    assert derive('failed', ['Normalization steps: 100%'], 'worker failed') == 'worker failed'
    assert derive('queued', ['Normalization steps: 100%'], 'Analysis queued.') == 'Analysis queued.'


def test_both_default_explore_serves_native_branch_embeddings_and_expression(combined_fixture, tmp_path):
    store,job,base,path = combined_fixture
    sup = branch(tmp_path,base,['s1','s2'],'ref',['A','B'])
    unsup = branch(tmp_path,base,['s2','s3'],'cell_state_predicted',['c1','c1'],True)
    meta,_ = _combine(job,store,path,base.obs_names,base.var_names,sup,unsup,'lzf')
    app = create_app({'JOB_STORAGE': str(store.root), 'ISOLATE_JOBS': False})
    app.state.job_store.update_job(job, **meta, status='completed', progress=100)
    with TestClient(app) as client:
        for key,expected in [('unsupervised_cluster', ['s2','s3']), ('unsupervised_state',['s2','s3']),
                             ('supervised_state',['s1','s2'])]:
            client.cookies.set(f'scalable_layer_{job}',key)
            response=client.get(f'/api/jobs/{job}/umap',params={'x_field':'__umap_1','y_field':'__umap_2'})
            assert response.status_code==200, response.text
            assert [p['barcode'] for p in response.json()['query']]==expected
            if key == 'supervised_state':
                assert response.json()['unplaced_population_counts']=={'Unaligned':1}
                # Display the retained unaligned barcode using its own recorded
                # unsupervised coordinates, without inventing a supervised point.
                alternate=client.get(f'/api/jobs/{job}/umap',params={
                    'color_by':'supervised_state','coords':'X_umap_unsupervised'})
                assert alternate.status_code==200
                assert {p['barcode']:p['population'] for p in alternate.json()['query']}=={'s2':'B','s3':'Unaligned'}
                assert alternate.json()['n_cells_selected']==3
            variables=client.get(f'/api/jobs/{job}/plot-variables').json()
            labels={v['label'] for v in variables['numeric_variables']}
            assert {'Supervised UMAP 1','Supervised UMAP 2','Unsupervised UMAP 1','Unsupervised UMAP 2'} <= labels
            expression=client.get(f'/api/jobs/{job}/expression',params={'gene':'gA','view':'umap',
                                  'x_field':'__umap_1','y_field':'__umap_2'})
            assert expression.status_code==200, expression.text
            points=expression.json()['umap']
            assert {p['barcode'] for p in points}==set(expected)
            assert {p['barcode']:p['value'] for p in points}=={s:float(base.X[base.obs_names.get_loc(s),0]) for s in expected}
            assert expression.json()['message'] is None
    app.state.job_runner.executor.shutdown(wait=True)


def test_empty_coordinate_expression_does_not_claim_present_gene_is_missing(combined_fixture):
    store,job,base,path = combined_fixture
    base.obsm['X_umap']=np.full((len(base),2),np.nan)
    base.obs['state']=['a','b','a','b'];base.write_h5ad(path)
    app=create_app({'JOB_STORAGE':str(store.root),'ISOLATE_JOBS':False})
    app.state.job_store.update_job(job,status='completed',cluster_key='state',artifacts={'combined_h5ad':str(path)})
    with TestClient(app) as client:
        response=client.get(f'/api/jobs/{job}/expression',params={'gene':'gA','view':'umap'})
        assert response.status_code==200
        assert response.json()['resolved_gene']=='gA'
        assert 'feature is present' in response.json()['message']
    app.state.job_runner.executor.shutdown(wait=True)


def test_supervised_progress_reports_loading_qc_normalization_and_alignment(tmp_path):
    from altanalyze3.components.cellHarmony.flask.tasks import _JobLogStream
    store=JobStore(tmp_path/'jobs');job=store.create_job('human','fixture',None,[])['job_id']
    stream=_JobLogStream(store,job)
    for line,step in [('adata shape: (20, 30)',1),('Running ambient RNA correction',2),
                      ('Cells remaining after min_genes filtering: 20',3),('Normalization steps: 100%',4),
                      ('Aligning cells to reference... cosine',5)]:
        stream.write(line);stream.flush()
        meta=store.get_job(job)
        assert f'cellHarmony step {step} of 10' in meta['message']
        assert web_mod._derive_live_pipeline_message('processing',['Normalization steps: 100%'],meta['message'])==meta['message']
    progress=store.get_job(job)['progress']
    stream.write('Normalization steps: 100%\n')
    assert store.get_job(job)['progress']==progress


def test_both_imputed_sidecars_keep_barcode_coordinates_and_differential_layers(combined_fixture,tmp_path):
    store,job,base,path=combined_fixture
    sup=branch(tmp_path,base,['s1','s2'],'ref',['A','B'])
    unsup=branch(tmp_path,base,['s2','s3'],'cell_state_predicted',['c1','c1'],True)
    sup_meta=sup.get_job(job)
    modalities=['lipids','adt','metabolite','lipid','grn','grn_tf']
    for modality in modalities:
        # Independently assigned software values, deliberately reversed barcode order.
        prediction=ad.AnnData(np.array([[20.,2.],[10.,1.]],dtype=np.float32),
                             obs=base.obs.loc[['s2','s1']].copy(),var=pd.DataFrame(index=['F1','F2']))
        output=tmp_path/(modality+'.h5ad');prediction.write_h5ad(output)
        sup_meta['modality_artifacts'][modality]={'h5ad':str(output)}
        sup_meta['modalities']['available'].append({'id':modality,'label':modality,'supports_differential':True})
    meta,output=_combine(job,store,path,base.obs_names,base.var_names,sup,unsup,'lzf')
    enable_differential(meta,{'files':[]},output,['unsupervised_cluster','unsupervised_state','supervised_state'])
    app=create_app({'JOB_STORAGE':str(store.root),'ISOLATE_JOBS':False})
    app.state.job_store.update_job(job,**meta,status='completed',progress=100)
    with TestClient(app) as client:
        for layer,expected in [('unsupervised_cluster',{'s2':20.}),('unsupervised_state',{'s2':20.}),
                               ('supervised_state',{'s1':10.,'s2':20.})]:
            client.cookies.set(f'scalable_layer_{job}',layer)
            for modality in modalities:
                response=client.get(f'/api/jobs/{job}/expression',params={'gene':'F1','modality':modality,'view':'umap'})
                assert response.status_code==200,response.text
                assert {p['barcode']:p['value'] for p in response.json()['umap']}==expected
                assert next(p for p in response.json()['umap'] if p['barcode']=='s2')['x']==(2. if layer=='supervised_state' else 0.)
            status=client.get(f'/api/jobs/{job}/status').json()['differential_ui']
            assert status['default_population_col']==layer
            assert {p['value'] for p in status['population_columns']}=={'unsupervised_cluster','unsupervised_state','supervised_state'}
            assert status['enabled']
            assert not {'supervised_state','unsupervised_cluster','supervised_accepted','original_NMF_cluster'} & {p['value'] for p in status['sample_fields']}
    app.state.job_runner.executor.shutdown(wait=True)


def test_log_download_contains_full_provenance_without_mutating_pipeline_log(tmp_path):
    import json
    app=create_app({'JOB_STORAGE':str(tmp_path),'ISOLATE_JOBS':False})
    store=app.state.job_store;job=store.create_job('human','fixture',None,[])['job_id']
    store.append_log(job,'Job completed.')
    provenance={'fixture-models':{'lipids':{'version':'fixture-L','inputs':['gA','gB']},
                                 'adt':{'version':'fixture-P','outputs':['F1','F2']}}}
    source=store.outputs_dir(job)/'model_provenance.json';source.write_text(json.dumps(provenance))
    store.update_job(job,status='completed',artifacts={'model_provenance':str(source)})
    logpath=store.logs_dir(job)/'pipeline.log';before=logpath.read_bytes()
    with TestClient(app) as client:
        response=client.get(f'/api/jobs/{job}/log')
        assert response.status_code==200
        assert response.content.startswith(before)
        assert 'Model versions and provenance' in response.text
        assert json.loads(response.text.split('=== Model versions and provenance ===\n')[1])==provenance
        assert 'attachment;' in response.headers['content-disposition']
    assert logpath.read_bytes()==before
    app.state.job_runner.executor.shutdown(wait=True)
