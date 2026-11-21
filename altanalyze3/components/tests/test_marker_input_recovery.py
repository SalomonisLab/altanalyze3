"""Imported-RNA recovery preserves data and keeps the native statistical gate."""
import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path
import pytest
from scipy import sparse

from altanalyze3.components.cellHarmony import markerFinder as native
from altanalyze3.components.cellHarmony.marker_input_recovery import prepare_marker_input, RecoveryUnavailable
from altanalyze3.components.cellHarmony.mapped_h5ad import Workspace
from altanalyze3.components.visualization.marker_heatmap_h5ad import generate_marker_heatmap_from_adata


def data(base=2., dense=False):
    rng = np.random.default_rng(7)
    shape = np.linspace(.1, 5, 120)[:, None]
    counts = rng.poisson(rng.gamma(shape, 10 / shape, (120, 60))).astype(float)
    counts[:60, :20] += 35
    counts[60:, 20:40] += 35
    abundance = counts / counts.sum(axis=1)[:, None] * 10000 / np.linspace(.99, 1.01, 120)[:, None]
    x = np.log1p(abundance) / np.log(base)
    obj = ad.AnnData(x if dense else sparse.csr_matrix(x),
                     obs=pd.DataFrame({'state': ['a']*60 + ['b']*60}, index=[f'c{i}' for i in range(120)]),
                     var=pd.DataFrame(index=[f'g{i}' for i in range(60)]))
    obj.obsm['X_umap'] = rng.normal(size=(120, 2))
    return obj, counts


@pytest.mark.parametrize('base', [2., np.e])
@pytest.mark.parametrize('dense', [False, True])
def test_inferred_recovery_native_hook_and_source_identity(tmp_path, base, dense):
    obj, counts = data(base, dense)
    original = obj.X.copy()
    with pytest.raises(native.MarkerFinderInputError):
        native.enforce_input_scaling(obj.X)
    work, options, record = prepare_marker_input(obj)
    assert record['log_base'] == pytest.approx(base)
    assert record['raw_counts_reconstructed'] is False
    assert work.obs_names.equals(obj.obs_names) and work.var_names.equals(obj.var_names)
    np.testing.assert_array_equal(obj.X if dense else obj.X.data, original if dense else original.data)
    assert not obj.layers
    assert native.detect_input_scaling(work.X)['status'] == 'ok'
    expected, _ = native.scale_expression_matrix(np.expm1((original if dense else original.toarray()) * np.log(base)))
    np.testing.assert_allclose(work.X if dense else work.X.toarray(), expected, atol=1e-12)
    outputs = generate_marker_heatmap_from_adata(work, cluster_key='state', out=str(tmp_path/'markers.pdf'),
        marker_method='markerfinder', render_heatmap=False, write_heatmap_cache=True,
        write_heatmap_tsv=False, write_expression_tsv=False, **options)
    assert Path(outputs['markers_tsv']).is_file() and Path(outputs['heatmap_cache']).is_file()


@pytest.mark.parametrize('source', ['counts', 'raw_counts', 'raw'])
def test_real_counts_precede_log_inference(source):
    obj, counts = data()
    if source == 'raw':
        raw = ad.AnnData(sparse.csr_matrix(counts[:, ::-1]), obs=obj.obs.copy(), var=obj.var.iloc[::-1].copy())
        obj.raw = raw
    else:
        obj.layers[source] = sparse.csr_matrix(counts)
    work, options, record = prepare_marker_input(obj)
    assert record['method'] == 'verified_linear_source' and options['scale_data'] is True
    matrix = work.layers[options['layer']] if 'layer' in options else work.X
    np.testing.assert_array_equal(matrix.toarray(), counts)


def test_mislabeled_counts_not_used_and_incomplete_raw_not_intersected():
    obj, _ = data()
    obj.layers['counts'] = obj.X.copy()
    obj.raw = obj[:, :-1].copy()
    work, options, record = prepare_marker_input(obj)
    assert work.shape == obj.shape and record['method'] == 'processed_abundance_log_inverse'
    assert record['inspected_sources'][-1]['missing_genes'] == ['g59']


def test_invalid_counts_source_does_not_hide_valid_imported_expression():
    obj, _ = data()
    obj.layers['counts'] = sparse.csr_matrix(obj.shape)
    assert prepare_marker_input(obj)[2]['method'] == 'processed_abundance_log_inverse'


def test_fortran_dense_input_and_raw_superset_preserve_all_rna_axes():
    obj, counts = data(dense=True)
    obj.X = np.asfortranarray(obj.X)
    assert native.detect_input_scaling(prepare_marker_input(obj)[0].X)['status'] == 'ok'
    extra = np.column_stack([counts, np.ones(obj.n_obs)])
    obj.raw = ad.AnnData(extra, obs=obj.obs.copy(), var=pd.DataFrame(index=obj.var_names.tolist()+['extra']))
    work, _, record = prepare_marker_input(obj)
    np.testing.assert_array_equal(work.X, counts)
    assert work.var_names.equals(obj.var_names) and obj.raw.n_vars == obj.n_vars+1
    assert record['source'] == 'raw.X'


@pytest.mark.parametrize('bad', ['negative', 'nan', 'ambiguous'])
def test_unverified_input_is_not_guessed(bad):
    obj, _ = data(dense=True)
    if bad == 'negative': obj.X[0, 0] = -1
    elif bad == 'nan': obj.X[0, 0] = np.nan
    else: obj.X = np.random.default_rng(10).uniform(0, 6, obj.shape)
    with pytest.raises(RecoveryUnavailable): prepare_marker_input(obj)


def test_recorded_base_and_overflow():
    obj, _ = data(base=10.)
    obj.uns['log1p'] = {'base': 10.}
    assert prepare_marker_input(obj)[2]['log_base'] == 10.
    obj.uns['log1p']['base'] = 1e300
    with pytest.raises(RecoveryUnavailable): prepare_marker_input(obj)


@pytest.mark.parametrize('dense', [False, True])
def test_disk_backed_native_normalization_matches_in_memory(tmp_path, dense):
    obj, _ = data(dense=dense)
    obj.X[-1] = 0
    expected, info = native.scale_expression_matrix(obj.X)
    workspace = Workspace(tmp_path/'mapped')
    mapped, actual = native.scale_expression_matrix(obj.X, workspace=workspace)
    np.testing.assert_allclose(mapped if dense else mapped.toarray(), expected if dense else expected.toarray(), atol=0, rtol=0)
    assert actual == info
    obj._matrix_workspace = workspace
    work, _, record = prepare_marker_input(obj)
    assert native.detect_input_scaling(work.X)['status'] == 'ok'
    assert record['method'] == 'processed_abundance_log_inverse'


def test_supervised_recovery_completes_alignment_and_explore(tmp_path):
    import json
    from fastapi.testclient import TestClient
    from altanalyze3.components.cellHarmony.webapp.app import create_app
    from altanalyze3.components.cellHarmony.flask.pipeline import run_cellharmony_pipeline
    obj, counts = data()
    obj.var_names = ['EPOR', 'CSF2RB'] + obj.var_names.tolist()[2:]
    original = obj.X.copy()
    source = tmp_path/'atlas.h5ad'; obj.write_h5ad(source)
    normalized, _ = native.scale_expression_matrix(counts)
    reference = pd.DataFrame({'a': normalized[:60].mean(axis=0), 'b': normalized[60:].mean(axis=0)}, index=obj.var_names)
    reference.to_csv(tmp_path/'states.tsv', sep='\t')
    pd.DataFrame({'barcode': ['r1','r2','r3','r4'], 'states': ['a','a','b','b']}).to_csv(tmp_path/'clusters.tsv',sep='\t',index=False)
    pd.DataFrame({'barcode': ['r1','r2','r3','r4'], 'UMAP1':[0,1,5,6], 'UMAP2':[0,1,5,6]}).to_csv(tmp_path/'coords.tsv',sep='\t',index=False)
    registry = tmp_path/'registry.json'
    registry.write_text(json.dumps({'species':[{'id':'human','references':[{'id':'fixture',
        'states_tsv':'states.tsv','reference_clusters_tsv':'clusters.tsv','reference_coords_tsv':'coords.tsv'}]}]}))
    app=create_app({'JOB_STORAGE':str(tmp_path/'jobs'),'ISOLATE_JOBS':False})
    store=app.state.job_store
    job=store.create_job('human','fixture',None,[])['job_id']
    store.update_job(job,analysis_mode='supervised',files=[{'filename':str(source),'sample_name':'atlas'}],qc={
        'analysis_mode':'supervised','min_genes':1,'min_counts':0,'min_cells':0,'mit_percent':100,
        'align_cutoff':0,'ambient_correction':'no','impute_modalities':[]})
    try:
        artifacts=run_cellharmony_pipeline(job,store,registry,export_approx_pdfs=False,build_bundle=False)
        store.update_job(job,status='completed',progress=100)
        result=ad.read_h5ad(artifacts['combined_h5ad'])
        np.testing.assert_array_equal(result.X.toarray(), original.toarray())
        assert result.obs_names.equals(obj.obs_names) and result.var_names.equals(obj.var_names)
        assert np.isfinite(result.obsm['X_umap']).all()
        marker=store.get_job(job)['marker_analysis']
        assert marker['enabled'] and marker['input_recovery']['method']=='processed_abundance_log_inverse'
        assert json.loads(result.uns['scalable_marker_input_recovery']) == marker['input_recovery']
        cache=ad.read_h5ad(marker['heatmap_cache'])
        assert min(cache.shape)>0 and np.isfinite(cache.X).all()
        with TestClient(app) as client:
            for route in ('umap','expression?gene=EPOR','expression?gene=CSF2RB','genes','marker/heatmap/viewer'):
                prefix='/jobs' if route=='marker/heatmap/viewer' else '/api/jobs'
                response=client.get(f'{prefix}/{job}/{route}')
                assert response.status_code==200, (route,response.text[:300])
                assert response.content
                if route.startswith('expression?'):
                    payload=response.json()
                    gene=route.split('=')[1]
                    assert payload['source']=='query' and payload['resolved_gene']==gene
                    assert len(payload['umap'])==obj.n_obs and payload['violin']
                    expected=dict(zip(obj.obs_names, original[:,obj.var_names.get_loc(gene)].toarray().ravel()))
                    for point in payload['umap']:
                        assert point['value']==expected[point['barcode']]
                elif route=='genes':
                    assert {'EPOR','CSF2RB'} <= set(response.json()['genes'])
        from altanalyze3.components.cellHarmony.analysis_record import build_manifest, _preprocessing_methods
        manifest,_=build_manifest(store.get_job(job),store._job_dir(job))
        assert 'not reconstructed raw counts' in '\n'.join(_preprocessing_methods(manifest))
    finally:
        app.state.job_runner.executor.shutdown(wait=True)
