"""Use Existing preserves supplied identities and never runs a clustering branch."""
import json
from pathlib import Path
from unittest.mock import Mock

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from fastapi.testclient import TestClient

from altanalyze3.components.cellHarmony.webapp.existing_annotations import inspect_annotations, suggest_annotation
from altanalyze3.components.cellHarmony.webapp.existing_workflow import prepare_annotations, supplied_embedding_result
from altanalyze3.components.cellHarmony.webapp.app import create_app
from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
from altanalyze3.components.cellHarmony.webapp.analysis_workflow import run_analysis_workflow
from altanalyze3.components.cellHarmony.flask import pipeline as P


def fixture(path, normalized=False, name='cell_type', n=40):
    rng = np.random.default_rng(34)
    counts = rng.poisson(3, (n, 30)).astype(np.float32)
    counts[:n//2, :10] += 30
    counts[n//2:, 10:20] += 30
    counts[-1] = 0
    obs = pd.DataFrame({name: ['Macrophage']*(n//2) + ['Epithelial']*(n-n//2),
                        'Library': ['case', 'control']*(n//2)}, index=[f'cell{i}' for i in range(n)])
    matrix = ad.AnnData(sp.csr_matrix(counts), obs=obs,
                        var=pd.DataFrame(index=['MRC1','EPCAM','MT-CO1','CXCL12','CXCR4','TGFB1','TGFBR1','TGFBR2','SPP1','CD44',
                                                  'ITGAV','ITGB1','IL6','IL6R','IL6ST','CSF1','CSF1R','VEGFA','KDR',
                                                  'JAG1','NOTCH1','DLL4','HGF','MET','FGF7','FGFR2','WNT5A','FZD2','IL1B','IL1R1']))
    matrix.obsm['X_umap'] = rng.normal(size=(n,2))
    matrix.obsm['X_pca'] = rng.normal(size=(n,4))
    matrix.layers['counts'] = matrix.X.copy()
    matrix.uns['source_tool'] = 'independent synthetic fixture'
    if normalized:
        matrix.X = sp.csr_matrix(np.log1p(counts / np.maximum(counts.sum(axis=1, keepdims=True), 1) * 10000))
    matrix.write_h5ad(path)
    return matrix


def test_clear_selection_does_not_call_model(tmp_path, monkeypatch):
    path=tmp_path/'a.h5ad'; fixture(path)
    report=inspect_annotations([path])
    assert report['selected']=='cell_type'
    monkeypatch.setattr('urllib.request.urlopen', Mock(side_effect=AssertionError('not ambiguous')))
    assert suggest_annotation(report)['method']=='column_metadata'


def test_ambiguous_local_model_is_constrained_to_existing_columns(tmp_path, monkeypatch):
    path=tmp_path/'a.h5ad'; matrix=fixture(path, name='annotation_A')
    matrix.obs['annotation_B']=matrix.obs.annotation_A; matrix.write_h5ad(path)
    report=inspect_annotations([path]); assert report['selected'] is None
    monkeypatch.setenv('SCALABLE_SUMMARY_URL','http://localhost:9/v1/chat/completions')
    monkeypatch.setenv('SCALABLE_SUMMARY_MODEL','fixture-model')
    class Response:
        def __enter__(self): return self
        def __exit__(self,*_): pass
        def read(self,*_): return json.dumps({'choices':[{'message':{'content':json.dumps(answer)}}]}).encode()
    answer={'column':'annotation_B','reason':'Contains named cell states'}
    monkeypatch.setattr('urllib.request.urlopen',lambda request,timeout: Response())
    selected=suggest_annotation(report)
    assert selected['selected']=='annotation_B' and selected['method']=='local_model'
    answer['column']='invented_cells'
    assert suggest_annotation(report)['selected'] is None
    answer['column']='Library'
    assert suggest_annotation(report)['selected'] is None


def test_unavailable_model_warns_and_numeric_cluster_is_detected(tmp_path, monkeypatch):
    path=tmp_path/'a.h5ad'; matrix=fixture(path, name='unusual_label')
    monkeypatch.delenv('SCALABLE_SUMMARY_URL',raising=False)
    assert suggest_annotation(inspect_annotations([path]))['model_status']=='unavailable'
    matrix.obs['leiden']=pd.Categorical(['0']*20+['1']*20); matrix.write_h5ad(path)
    report=inspect_annotations([path])
    assert report['selected']=='leiden'
    assert next(c for c in report['candidates'] if c['key']=='leiden')['kind']=='cluster'


@pytest.mark.parametrize('normalized',[False,True])
def test_requested_qc_preserves_supplied_obs_embeddings_and_expression(tmp_path, normalized):
    path=tmp_path/'a.h5ad'; source=fixture(path,normalized)
    _, result=P.cellHarmony_lite.combine_and_align_h5(h5_files=[],h5ad_file=str(path),cellharmony_ref=None,
        output_dir=str(tmp_path/'qc'), export_cptt=False,generate_umap=False,save_adata=False,export_h5ad=False,
        min_genes=1,min_cells=0,min_counts=1,mit_percent=100,return_adata=True,
        bounded_h5ad=True,qc_normalized_input=True)
    assert list(result.obs_names)==list(source.obs_names[:-1])
    assert result.var_names.equals(source.var_names)
    assert result.obs['cell_type'].tolist()==source.obs['cell_type'].tolist()[:-1]
    for key in source.obsm:
        np.testing.assert_array_equal(result.obsm[key],source.obsm[key][:-1])
    np.testing.assert_array_equal(result.layers['counts'].toarray(), source.layers['counts'][:-1].toarray())
    if normalized:
        np.testing.assert_array_equal(result.X.toarray(),source.X[:-1].toarray())


def test_missing_counts_does_not_treat_log_values_as_counts(tmp_path):
    path=tmp_path/'a.h5ad'; matrix=fixture(path,True); del matrix.layers['counts']; matrix.write_h5ad(path)
    with pytest.raises(ValueError,match='requires counts'):
        P.cellHarmony_lite.combine_and_align_h5(h5_files=[],h5ad_file=str(path),cellharmony_ref=None,
            output_dir=str(tmp_path/'qc'),export_cptt=False,generate_umap=False,
            min_genes=1,min_counts=1000,qc_normalized_input=True)


def test_upload_annotation_inspection_and_qc_validation(tmp_path,monkeypatch):
    path=tmp_path/'a.h5ad'; fixture(path)
    app=create_app({'JOB_STORAGE':str(tmp_path/'jobs'),'ISOLATE_JOBS':False})
    with TestClient(app) as client:
        assert 'Use Existing' in client.get('/').text
        response=client.post('/api/jobs',data={'species':'human','analysis_mode':'existing','sample_names':'a'},
                             files={'files':('a.h5ad',path.read_bytes(),'application/octet-stream')})
        assert response.status_code==200,response.text
        job=response.json()['job_id']
        report=client.get(f'/api/jobs/{job}/existing-annotations')
        assert report.status_code==200,report.text
        assert report.json()['selected']=='cell_type'
        assert client.post(f'/api/jobs/{job}/qc',json={'analysis_mode':'existing','existing_state_column':'fake'}).status_code==400
        response=client.post(f'/api/jobs/{job}/qc',json={'analysis_mode':'existing','existing_state_column':'cell_type'})
        assert response.status_code==200,response.text
        assert client.get(f'/api/jobs/{job}/status').json()['qc']['existing_state_column']=='cell_type'
    app.state.job_runner.executor.shutdown(wait=True)


def test_multiple_inputs_require_a_common_cell_state_column(tmp_path):
    paths=[tmp_path/'a.h5ad',tmp_path/'b.h5ad']
    fixture(paths[0],name='cell_type'); fixture(paths[1],name='leiden')
    report=inspect_annotations(paths)
    assert report['selected'] is None
    assert not any(c['key']=='cell_type' for c in report['candidates'])
    assert any('No unambiguous' in w for w in report['warnings'])


def test_missing_labels_retained_and_source_column_not_overwritten(tmp_path):
    path=tmp_path/'a.h5ad'; source=fixture(path)
    source.obs.loc['cell0','cell_type']=None; source.write_h5ad(path)
    store=JobStore(tmp_path/'jobs'); job=store.create_job('human','existing',None,[])['job_id']
    meta=store.update_job(job,files=[{'filename':str(path),'sample_name':'a'}],qc={'existing_state_column':'cell_type'})
    key=prepare_annotations(source,meta,store,job)
    assert pd.isna(source.obs.loc['cell0','cell_type'])
    assert source.obs.loc['cell0',key]=='Unassigned'
    assert len(source)==40
    del source.obsm['X_umap']
    with pytest.raises(ValueError, match='UMAP or t-SNE coordinates missing'):
        supplied_embedding_result(source,key,store,job)
    assert 'X_umap' not in source.obsm


def test_pipeline_uses_standard_downstream_hooks_without_alignment_or_icgs(tmp_path,monkeypatch):
    path=tmp_path/'a.h5ad'; source=fixture(path,True)
    store=JobStore(tmp_path/'jobs'); job=store.create_job('human','existing',None,[])['job_id']
    store.update_job(job,analysis_mode='existing', files=[{'filename':str(path),'sample_name':'a'}],qc={
        'analysis_mode':'existing','existing_state_column':'cell_type','min_genes':1,'min_counts':1,
        'min_cells':0,'mit_percent':100,'ambient_correction':'no','impute_modalities':[]})
    # Sentinels prohibit both reference placement and clustering. MarkerFinder and
    # communication are separate standard downstream hooks, exercised below in integration validation.
    monkeypatch.setattr(P.approx_mod,'approximate_umap',Mock(side_effect=AssertionError('must not align')))
    monkeypatch.setattr(P.approx_mod,'_load_reference_from_tsv',Mock(side_effect=AssertionError('must not load reference')))
    from altanalyze3.components.clustering import ICGS
    monkeypatch.setattr(ICGS,'run_icgs3',Mock(side_effect=AssertionError('must not cluster')))
    marker=Mock(return_value={'markers_tsv':path,'centroids_tsv':path,'networks':[]})
    monkeypatch.setattr(P.marker_mod,'generate_marker_heatmap_from_adata',marker)
    comm=Mock(return_value=({'enabled':False},{})); monkeypatch.setattr(P,'_run_fastcomm_analysis',comm)
    monkeypatch.setattr(P,'_build_job_bundle',lambda *a,**k:{})
    artifacts=run_analysis_workflow(job,store,tmp_path/'nonexistent-registry.json')
    output=ad.read_h5ad(artifacts['combined_h5ad'])
    np.testing.assert_array_equal(output.X.toarray(),source.X[:-1].toarray())
    np.testing.assert_array_equal(output.obsm['X_umap'],source.obsm['X_umap'][:-1])
    assert marker.call_count==1 and comm.call_count==1
    assert store.get_job(job)['cluster_key'] in output.obs
    assert store.get_job(job)['differential_options']['enabled']


def test_embedded_prediction_reuse_keeps_values_features_and_provenance(tmp_path,monkeypatch):
    from altanalyze3.components.cellHarmony.webapp.existing_workflow import reuse_embedded_modalities
    source=fixture(tmp_path/'a.h5ad',True)
    values=np.arange(source.n_obs*3,dtype=np.float64).reshape(-1,3)/7
    source.obsm['X_lipids']=values
    source.uns['lipid_feature_names']=['lipid-A','lipid-B','lipid-C']
    source.uns['imputed_modalities']={'lipids':{'obsm_key':'X_lipids','feature_names_key':'lipid_feature_names',
        'summary':{'model_version_id':'original-production-id','expression_scale':'linear'}}}
    marker=Mock(return_value={'enabled':True})
    monkeypatch.setattr(P,'_emit_modality_marker_heatmap',marker)
    artifacts,markers,warnings=reuse_embedded_modalities(source,tmp_path,[],'cell_type',{},'lzf')
    output=ad.read_h5ad(artifacts['lipids']['h5ad'])
    np.testing.assert_array_equal(output.X,values)
    assert list(output.var_names)==['lipid-A','lipid-B','lipid-C']
    assert output.uns['prediction_summary']['model_version_id']=='original-production-id'
    assert not warnings and 'lipids' in markers


def test_prior_scalable_output_remembers_its_existing_layer(tmp_path):
    path=tmp_path/'a.h5ad'; matrix=fixture(path)
    matrix.obs['scalable_existing_state']=matrix.obs.cell_type
    matrix.uns['scalable_existing_analysis']={'effective_column':'scalable_existing_state'}
    matrix.write_h5ad(path)
    assert inspect_annotations([path])['selected']=='scalable_existing_state'


def test_user_override_records_actual_selected_column(tmp_path):
    path=tmp_path/'input.h5ad'; matrix=fixture(path)
    matrix.obs['alternate_state']=matrix.obs.cell_type; matrix.write_h5ad(path)
    app=create_app({'JOB_STORAGE':str(tmp_path/'jobs'),'ISOLATE_JOBS':False})
    with TestClient(app) as client:
        response=client.post('/api/jobs',data={'species':'human','analysis_mode':'existing','sample_names':'a'},
            files={'files':('input.h5ad',path.read_bytes(),'application/octet-stream')})
        job=response.json()['job_id']
        client.get(f'/api/jobs/{job}/existing-annotations')
        response=client.post(f'/api/jobs/{job}/qc',json={'analysis_mode':'existing','existing_state_column':'alternate_state'})
        assert response.status_code==200,response.text
        record=app.state.job_store.get_job(job)['existing_annotations']
        assert record['selected']=='alternate_state' and record['method']=='user_selection'
    app.state.job_runner.executor.shutdown(wait=True)


def test_unidentified_annotation_warns_and_cannot_start_with_empty_selection(tmp_path,monkeypatch):
    path=tmp_path/'a.h5ad'; matrix=fixture(path)
    matrix.obs=matrix.obs[['Library']]; matrix.write_h5ad(path)
    monkeypatch.delenv('SCALABLE_SUMMARY_URL',raising=False)
    app=create_app({'JOB_STORAGE':str(tmp_path/'jobs'),'ISOLATE_JOBS':False})
    with TestClient(app) as client:
        response=client.post('/api/jobs',data={'species':'human','analysis_mode':'existing','sample_names':'a'},
            files={'files':('a.h5ad',path.read_bytes(),'application/octet-stream')})
        job=response.json()['job_id']
        report=client.get(f'/api/jobs/{job}/existing-annotations').json()
        assert report['selected'] is None and report['warnings']
        response=client.post(f'/api/jobs/{job}/qc',json={'analysis_mode':'existing'})
        assert response.status_code==400 and 'Select an existing' in response.text
    app.state.job_runner.executor.shutdown(wait=True)


def raw_tsne_fixture(path):
    matrix=fixture(path,True)
    raw=ad.AnnData(X=matrix.layers['counts'].copy(),obs=matrix.obs.copy(),var=matrix.var.copy())
    matrix.raw=raw
    del matrix.layers['counts']
    matrix.obsm['X_tSNE']=matrix.obsm.pop('X_umap')
    matrix.uns['default_embedding']='X_tSNE'
    matrix.write_h5ad(path)
    return matrix


def test_raw_counts_and_tsne_detected_without_false_missing_warning(tmp_path):
    path=tmp_path/'input.h5ad'; raw_tsne_fixture(path)
    report=inspect_annotations([path])
    assert report['inputs'][0]['count_source']=='raw/X'
    assert report['embedding_key']=='X_tSNE' and report['embedding_label']=='t-SNE'
    assert not report['warnings']


def test_logged_raw_is_not_misidentified_as_counts(tmp_path):
    path=tmp_path/'input.h5ad'; matrix=fixture(path,True)
    matrix.raw=matrix.copy();del matrix.layers['counts'];matrix.write_h5ad(path)
    report=inspect_annotations([path])
    assert report['inputs'][0]['count_source'] is None
    assert any('Raw counts were not identified' in w for w in report['warnings'])


def test_raw_counts_qc_preserves_normalized_rna_and_all_source_features(tmp_path):
    path=tmp_path/'input.h5ad'; source=raw_tsne_fixture(path)
    _,result=P.cellHarmony_lite.combine_and_align_h5(h5_files=[],h5ad_file=str(path),cellharmony_ref=None,
        output_dir=str(tmp_path/'qc'),export_cptt=False,generate_umap=False,save_adata=False,export_h5ad=False,
        min_genes=1,min_cells=0,min_counts=1,mit_percent=100,return_adata=True,bounded_h5ad=True,qc_normalized_input=True)
    assert result.var_names.equals(source.var_names) and result.raw.var_names.equals(source.raw.var_names)
    np.testing.assert_array_equal(result.X.toarray(),source.X[:-1].toarray())
    np.testing.assert_array_equal(result.layers['counts'].toarray(),source.raw.X[:-1].toarray())
    np.testing.assert_array_equal(result.raw.X.toarray(),source.raw.X[:-1].toarray())
    np.testing.assert_array_equal(result.obsm['X_tSNE'],source.obsm['X_tSNE'][:-1])
    assert 'X_umap' not in result.obsm


def test_native_existing_pipeline_uses_tsne_without_inventing_umap(tmp_path):
    path=tmp_path/'input.h5ad';source=raw_tsne_fixture(path)
    app=create_app({'JOB_STORAGE':str(tmp_path/'jobs'),'ISOLATE_JOBS':False})
    store=app.state.job_store
    job=store.create_job('human','existing',None,[])['job_id']
    store.update_job(job,analysis_mode='existing',files=[{'filename':str(path),'sample_name':'fixture'}],qc={
        'analysis_mode':'existing','existing_state_column':'cell_type','min_genes':1,'min_counts':1,
        'min_cells':0,'mit_percent':100,'ambient_correction':'no','impute_modalities':[]})
    artifacts=run_analysis_workflow(job,store,tmp_path/'unused.json')
    result=ad.read_h5ad(artifacts['combined_h5ad'])
    assert 'X_umap' not in result.obsm and 'umap_coordinates' not in artifacts
    np.testing.assert_array_equal(result.obsm['X_tSNE'],source.obsm['X_tSNE'][:-1])
    store.update_job(job,status='completed',progress=100)
    with TestClient(app) as client:
        response=client.get(f'/api/jobs/{job}/umap')
        assert response.status_code==200,response.text
        payload=response.json()
        assert 't-SNE' in payload['x_label'] and 't-SNE' in payload['y_label']
        assignments=pd.read_csv(artifacts['assignments'],sep='\t')
        assert 't-SNE 1' in assignments and 'UMAP1' not in assignments
    app.state.job_runner.executor.shutdown(wait=True)


def test_counts_override_preserves_original_layers_and_drives_qc(tmp_path):
    path = tmp_path/'input.h5ad'; source = fixture(path, True)
    # Different sources give different QC outcomes. The requested source, rather
    # than the conventional counts name, must determine the retained barcodes.
    alternative = source.layers['counts'].copy()
    alternative[0, :] = 0
    source.layers['raw_counts'] = alternative
    source.write_h5ad(path)
    report = inspect_annotations([path])
    assert {'layers/counts', 'layers/raw_counts'} <= {c['key'] for c in report['count_candidates']}
    _, result = P.cellHarmony_lite.combine_and_align_h5(h5_files=[], h5ad_file=str(path), cellharmony_ref=None,
        output_dir=str(tmp_path/'qc'), export_cptt=False, generate_umap=False, save_adata=False, export_h5ad=False,
        min_genes=1, min_cells=0, min_counts=1, mit_percent=100, return_adata=True,
        bounded_h5ad=True, qc_normalized_input=True, existing_counts_source='layers/raw_counts')
    assert result.obs_names.equals(source.obs_names[1:-1])
    assert result.var_names.equals(source.var_names)
    np.testing.assert_array_equal(result.X.toarray(), source.X[1:-1].toarray())
    np.testing.assert_array_equal(result.layers['counts'].toarray(), alternative[1:-1].toarray())
    np.testing.assert_array_equal(result.layers['source_counts'].toarray(), source.layers['counts'][1:-1].toarray())
    assert result.uns['scalable_count_source'] == 'layers/raw_counts'


def test_counts_in_X_retained_before_normalization(tmp_path):
    path = tmp_path/'input.h5ad'; source = fixture(path)
    del source.layers['counts']; source.write_h5ad(path)
    _, result = P.cellHarmony_lite.combine_and_align_h5(h5_files=[], h5ad_file=str(path), cellharmony_ref=None,
        output_dir=str(tmp_path/'qc'), export_cptt=False, generate_umap=False, save_adata=False, export_h5ad=False,
        min_genes=0, min_cells=0, min_counts=0, mit_percent=100, return_adata=True,
        bounded_h5ad=True, qc_normalized_input=True, existing_counts_source='X')
    np.testing.assert_array_equal(result.layers['counts'].toarray(), source.X.toarray())
    assert not np.array_equal(result.X.toarray(), source.X.toarray())


@pytest.mark.parametrize('axes', [('UMAP_1','UMAP_2'), ('tSNE_x','tSNE_y'), ('umap_0','umap_1')])
def test_obs_layout_slots_are_used_without_overwriting_source(tmp_path, axes):
    path = tmp_path/'input.h5ad'; source = fixture(path)
    coordinates = source.obsm.pop('X_umap')
    source.obs[axes[0]] = coordinates[:, 0]; source.obs[axes[1]] = coordinates[:, 1]
    source.write_h5ad(path)
    report = inspect_annotations([path])
    key = report['embedding_key']
    assert key.startswith('obs:')
    store=JobStore(tmp_path/'jobs'); job=store.create_job('human','existing',None,[])['job_id']
    store.update_job(job, qc={'existing_embedding_key': key})
    result = supplied_embedding_result(source, 'cell_type', store, job)
    np.testing.assert_array_equal(result.coordinates.to_numpy(), coordinates)
    np.testing.assert_array_equal(source.obs[list(axes)].to_numpy(), coordinates)


def test_missing_layout_blocked_before_run_or_output_invalidation(tmp_path, monkeypatch):
    path=tmp_path/'input.h5ad'; source=fixture(path)
    del source.obsm['X_umap']; source.write_h5ad(path)
    app=create_app({'JOB_STORAGE':str(tmp_path/'jobs'),'ISOLATE_JOBS':False})
    store=app.state.job_store
    job=store.create_job('human','existing',None,[])['job_id']
    store.update_job(job, analysis_mode='existing', files=[{'filename':str(path),'sample_name':'a'}],
                     qc={'analysis_mode':'existing','existing_state_column':'cell_type'},artifacts={'prior':'retained'})
    submit=Mock(side_effect=AssertionError('No worker should launch')); monkeypatch.setattr(app.state.job_runner,'submit',submit)
    with TestClient(app) as client:
        response=client.post(f'/api/jobs/{job}/qc',json={'analysis_mode':'existing','existing_state_column':'cell_type'})
        assert response.status_code==400 and 'UMAP or t-SNE coordinates missing' in response.text
        response=client.post(f'/api/jobs/{job}/run')
        assert response.status_code==400 and 'UMAP or t-SNE coordinates missing' in response.text
    assert store.get_job(job)['artifacts']=={'prior':'retained'}
    with pytest.raises(ValueError, match='UMAP or t-SNE coordinates missing'):
        run_analysis_workflow(job, store, tmp_path/'unused.json')
    submit.assert_not_called()
    app.state.job_runner.executor.shutdown(wait=True)


def test_layout_and_counts_choices_validated_and_recorded(tmp_path):
    path=tmp_path/'input.h5ad'; source=raw_tsne_fixture(path)
    source.obsm['X_umap']=source.obsm['X_tSNE']+10; source.write_h5ad(path)
    app=create_app({'JOB_STORAGE':str(tmp_path/'jobs'),'ISOLATE_JOBS':False})
    store=app.state.job_store; job=store.create_job('human','existing',None,[])['job_id']
    store.update_job(job, files=[{'filename':str(path),'sample_name':'a'}])
    settings={'analysis_mode':'existing','existing_state_column':'cell_type',
              'existing_embedding_key':'X_umap','existing_counts_source':'raw/X'}
    with TestClient(app) as client:
        assert client.post(f'/api/jobs/{job}/qc',json={**settings,'existing_embedding_key':'X_pca'}).status_code==400
        assert client.post(f'/api/jobs/{job}/qc',json={**settings,'existing_counts_source':'X'}).status_code==400
        response=client.post(f'/api/jobs/{job}/qc',json=settings)
        assert response.status_code==200,response.text
    assert store.get_job(job)['qc']['existing_embedding_key']=='X_umap'
    result=supplied_embedding_result(source,'cell_type',store,job)
    np.testing.assert_array_equal(result.coordinates.to_numpy(),source.obsm['X_umap'])
    app.state.job_runner.executor.shutdown(wait=True)


def test_nonfinite_and_latent_coordinates_are_not_valid_layouts(tmp_path):
    path=tmp_path/'input.h5ad'; matrix=fixture(path)
    matrix.obsm['X_umap'][0,0]=np.nan; matrix.write_h5ad(path)
    report=inspect_annotations([path])
    assert not report['embedding_key'] and not report['embedding_candidates']


def test_imported_normalization_mismatch_uses_standard_counts_route_and_preserves_X(tmp_path, monkeypatch):
    path=tmp_path/'input.h5ad'; source=fixture(path,True)
    # Imported logs can describe normalization over an upstream panel, whose
    # excluded contributions cause >0.1% variation on the supplied panel.
    counts=source.layers['counts'].toarray()
    multiplier=1 + np.linspace(0, .008, source.n_obs)[:, None]
    source.X=sp.csr_matrix(np.log1p(counts / np.maximum(counts.sum(axis=1,keepdims=True),1) * 10000 / multiplier))
    source.write_h5ad(path)
    store=JobStore(tmp_path/'jobs'); job=store.create_job('human','existing',None,[])['job_id']
    store.update_job(job,analysis_mode='existing',files=[{'filename':str(path),'sample_name':'a'}],qc={
        'analysis_mode':'existing','existing_state_column':'cell_type','existing_embedding_key':'X_umap',
        'min_genes':1,'min_counts':1,'min_cells':0,'mit_percent':100,'impute_modalities':[]})
    hook=Mock(wraps=P.marker_mod.generate_marker_heatmap_from_adata)
    monkeypatch.setattr(P.marker_mod,'generate_marker_heatmap_from_adata',hook)
    artifacts=run_analysis_workflow(job,store,tmp_path/'unused.json')
    assert hook.call_count==2
    assert 'scale_data' not in hook.call_args_list[0].kwargs
    assert hook.call_args.kwargs['layer']=='counts' and hook.call_args.kwargs['scale_data'] is True
    assert hook.call_args.kwargs.get('validate_scaling',True) is True
    result=ad.read_h5ad(artifacts['combined_h5ad'])
    assert result.var_names.equals(source.var_names) and result.obs_names.equals(source.obs_names[:-1])
    np.testing.assert_array_equal(result.X.toarray(),source.X[:-1].toarray())
    np.testing.assert_array_equal(result.layers['counts'].toarray(),source.layers['counts'][:-1].toarray())
    np.testing.assert_array_equal(result.obsm['X_umap'],source.obsm['X_umap'][:-1])
    assert store.get_job(job)['marker_analysis']['input_options']=={'layer':'counts','scale_data':True}
    assert Path(store.get_job(job)['marker_analysis']['markers_tsv']).is_file()
    assert Path(store.get_job(job)['marker_analysis']['heatmap_cache']).is_file()


def test_imported_normalization_mismatch_without_counts_still_requires_real_counts(tmp_path):
    path=tmp_path/'input.h5ad'; source=fixture(path,True)
    values=source.X.toarray()
    source.X=sp.csr_matrix(np.log1p(np.expm1(values)/(1+np.linspace(0,.008,source.n_obs)[:,None])))
    del source.layers['counts']; source.write_h5ad(path)
    store=JobStore(tmp_path/'jobs'); job=store.create_job('human','existing',None,[])['job_id']
    store.update_job(job,analysis_mode='existing',files=[{'filename':str(path),'sample_name':'a'}],qc={
        'analysis_mode':'existing','existing_state_column':'cell_type','existing_embedding_key':'X_umap',
        'min_genes':1,'min_counts':0,'min_cells':0,'mit_percent':100,'impute_modalities':[]})
    with pytest.raises(P.marker_mod._load_marker_finder().MarkerFinderInputError):
        run_analysis_workflow(job,store,tmp_path/'unused.json')
    assert source.n_obs==40 and source.n_vars==30
