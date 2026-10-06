import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import pytest

from altanalyze3.components.cellHarmony.mapped_h5ad import Workspace, inspect_h5ad, merge_h5ads
from altanalyze3.components.cellHarmony.cellHarmony_lite import combine_and_align_h5


@pytest.mark.parametrize('whole', [True, False])
def test_cache_advice_preserves_flushed_scratch_values(tmp_path, monkeypatch, whole):
    import os
    import mmap
    calls = []
    monkeypatch.setattr(os, 'POSIX_FADV_DONTNEED', 4, raising=False)
    monkeypatch.setattr(os, 'posix_fadvise',
                        lambda fd, offset, length, advice: calls.append((offset, length, advice)),
                        raising=False)
    workspace = Workspace(tmp_path / 'work')
    values = workspace.array((mmap.PAGESIZE * 4,), np.float64)
    values[:] = np.arange(values.size)
    if whole:
        workspace.release_pages()
        assert calls == [(0, 0, 4)]
    else:
        workspace.release_range(np.asarray(values), 0, mmap.PAGESIZE, write=True)
        assert calls == [(0, mmap.PAGESIZE * 8, 4)]
    np.testing.assert_array_equal(np.fromfile(values.filename, dtype=np.float64), np.arange(values.size))
    np.testing.assert_array_equal(values, np.arange(values.size))


def test_unsupported_cache_advice_does_not_corrupt_scratch(tmp_path, monkeypatch):
    import os
    def unsupported(*args):
        raise OSError('unsupported filesystem advice')
    monkeypatch.setattr(os, 'POSIX_FADV_DONTNEED', 4, raising=False)
    monkeypatch.setattr(os, 'posix_fadvise', unsupported, raising=False)
    workspace = Workspace(tmp_path / 'work')
    values = workspace.array((16,), np.float32)
    values[:] = np.arange(16)
    workspace.release_pages()
    np.testing.assert_array_equal(np.fromfile(values.filename, dtype=np.float32), np.arange(16))


def test_range_release_flushes_only_the_written_mapping_region(tmp_path):
    import mmap
    workspace = Workspace(tmp_path / 'ranges')
    mapped = workspace.array((mmap.PAGESIZE * 8,), np.float64)
    mapped[:] = np.arange(mapped.size)
    workspace.release_pages()
    original = mapped._mmap
    calls = []

    class ObserveMapping:
        def flush(self, offset, length):
            calls.append(('flush', offset, length))
            return original.flush(offset, length)

        def madvise(self, option, offset, length):
            calls.append(('release', offset, length))
            return original.madvise(option, offset, length)

    mapped._mmap = ObserveMapping()
    view = np.asarray(mapped[mmap.PAGESIZE:2 * mmap.PAGESIZE])
    view[123:456] += 1
    workspace.release_range(view, 123, 456, write=True)
    assert calls[0][0] == 'flush' and calls[1][0] == 'release'
    assert calls[0][1:] == calls[1][1:]
    assert calls[0][1] % mmap.PAGESIZE == 0
    assert calls[0][2] < mapped.nbytes
    calls.clear()
    workspace.release_range(view, 123, 456)
    assert [entry[0] for entry in calls] == ['release']
    expected = np.arange(mapped.size, dtype=np.float64)
    expected[mmap.PAGESIZE + 123:mmap.PAGESIZE + 456] += 1
    np.testing.assert_array_equal(mapped, expected)
    mapped._mmap = original


@pytest.mark.parametrize('encoding', ['csr', 'csc', 'dense'])
@pytest.mark.parametrize('rows', [[2, 0, 2, 1], [], [-1, 0]])
@pytest.mark.parametrize('columns', [None, [0, 1, 2, 3], [3, 1, 3]])
def test_bounded_selection_preserves_values_order_and_explicit_zeros(tmp_path, encoding, rows, columns):
    source = sp.csr_matrix((np.array([5., 0., 2., 7., 8.]),
                            np.array([3, 0, 1, 2, 0]),
                            np.array([0, 3, 3, 5])), shape=(3, 4))
    matrix = {'csr': lambda: source, 'csc': source.tocsc,
              'dense': source.toarray}[encoding]()
    rows = np.asarray(rows, dtype=np.int64)
    cols = None if columns is None else np.asarray(columns, dtype=np.int64)
    expected = matrix[rows, :]
    if cols is not None:
        expected = expected[:, cols]
    actual = Workspace(tmp_path / 'selection').select(matrix, rows, cols)
    if sp.issparse(expected):
        expected = expected.tocsr()
        np.testing.assert_array_equal(actual.indptr, expected.indptr)
        np.testing.assert_array_equal(actual.indices, expected.indices)
        np.testing.assert_array_equal(actual.data, expected.data)
    else:
        np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize('cells,size,number,expected', [
    (100_000, 1024**3 - 1, 1, False),
    (100_001, 100, 1, True),
    (65_662, 1024**3, 1, True),
    (3, 100, 2, True),
])
def test_disk_import_decision_uses_uncompressed_size_cells_and_file_count(cells, size, number, expected):
    from altanalyze3.components.cellHarmony.mapped_h5ad import needs_disk_backed_import
    assert needs_disk_backed_import([{'cells': cells, 'matrix_bytes': size}] * number) is expected


def test_skipped_qc_notice_survives_long_analysis_logs(tmp_path):
    from altanalyze3.components.cellHarmony.webapp.log_snapshot import read_pipeline_log
    path = tmp_path / 'pipeline.log'
    skipped = '...skipping QC: input is already scaled and log-transformed\n'
    path.write_text('work\n' * 90 + skipped + 'UMAP work\n' * 300)
    head, tail, progress = read_pipeline_log(path)
    assert skipped not in head and skipped not in tail
    assert skipped in progress


@pytest.mark.parametrize('encoding', ['csr', 'csc', 'dense'])
@pytest.mark.parametrize('logged', [False, True])
def test_import_qc_normalization_alignment_match_baseline(tmp_path, encoding, logged):
    rng = np.random.default_rng(31)
    values = rng.poisson(2, (81, 17)).astype(np.float32)
    values[0] = 0
    values[1, :15] = 0
    values[2, 0] = 300
    values[:, -1] = 0
    values[3, -1] = 1
    counts = values.copy()
    if logged:
        values = np.log1p(values / np.maximum(values.sum(axis=1, keepdims=True), 1) * 1e4)
    convert = {'csr': sp.csr_matrix, 'csc': sp.csc_matrix, 'dense': lambda v: v}[encoding]
    source = ad.AnnData(convert(values), obs=pd.DataFrame({'Library': ['lib'] * 81},
                       index=[f'c{i}' for i in range(81)]),
                       var=pd.DataFrame(index=['MT-a'] + [f'g{i}' for i in range(16)]))
    source.layers['counts'] = convert(counts)
    source.obsm['user_embedding'] = rng.random((81, 2))
    source.raw = source.copy()
    source.uns['user_metadata'] = {'value': 'preserve'}
    path = tmp_path / 'query.h5ad'
    source.write_h5ad(path)
    reference = tmp_path / 'ref.tsv'
    pd.DataFrame(rng.random((17, 3)), index=source.var_names, columns=['A', 'B', 'C']).to_csv(reference, sep='\t')
    kwargs = dict(h5_files=[], h5ad_file=path, cellharmony_ref=reference, min_genes=3,
                  min_cells=0, min_counts=10, mit_percent=30, min_alignment_score=.45,
                  return_adata=True, concat_batch_size=0)
    _, expected = combine_and_align_h5(**kwargs, output_dir=str(tmp_path / 'baseline'))
    _, actual = combine_and_align_h5(**kwargs, output_dir=str(tmp_path / 'bounded'), bounded_h5ad=True)
    def dense(x):
        return x.toarray() if sp.issparse(x) else np.asarray(x)
    assert list(actual.obs_names) == list(expected.obs_names)
    assert list(actual.var_names) == list(expected.var_names)
    pd.testing.assert_frame_equal(actual.obs, expected.obs)
    np.testing.assert_allclose(dense(actual.X), dense(expected.X), rtol=1e-6)
    np.testing.assert_array_equal(dense(actual.layers['counts']), dense(expected.layers['counts']))
    np.testing.assert_array_equal(dense(actual.raw.X), dense(expected.raw.X))
    np.testing.assert_array_equal(actual.obsm['user_embedding'], expected.obsm['user_embedding'])
    assert actual.uns['user_metadata'] == expected.uns['user_metadata']
    assert getattr(actual, '_matrix_workspace', None) is not None


@pytest.mark.parametrize('encoding', ['csr', 'csc', 'dense', 'mixed'])
def test_disk_merge_preserves_union_and_sample_annotations(tmp_path, encoding):
    entries = []
    for i, genes in enumerate((['a', 'b'], ['b', 'c'])):
        data = ad.AnnData(sp.csr_matrix([[1, 2], [3, 4]], dtype=np.float32),
                         obs=pd.DataFrame({'donor': ['D', 'E']}, index=['same1', 'same2']),
                         var=pd.DataFrame({'gene_symbols': [gene.upper() for gene in genes]}, index=genes))
        if encoding == 'csc':
            data.X = data.X.tocsc()
        elif encoding == 'dense' or (encoding == 'mixed' and i == 1):
            data.X = data.X.toarray()
        data.layers['counts'] = data.X.copy()
        data.raw = data.copy()
        data.obsp['neighbors'] = sp.csr_matrix([[0, i + 1], [i + 1, 0]], dtype=np.float32)
        path = tmp_path / f'{i}.h5ad'
        data.write_h5ad(path)
        entries.append((path, f'lib{i}'))
    result = Workspace(tmp_path / 'work').load(merge_h5ads(entries, tmp_path / 'merge'))
    assert result.shape == (4, 3)
    assert list(result.var_names) == ['a', 'b', 'c']
    assert list(result.var['gene_symbols']) == ['A', 'B', 'C']
    assert list(result.raw.var['gene_symbols']) == ['A', 'B', 'C']
    assert list(result.obs_names) == ['same1::lib0', 'same2::lib0', 'same1::lib1', 'same2::lib1']
    assert list(result.obs['donor']) == ['D', 'E', 'D', 'E']
    assert list(result.obs['Library']) == ['lib0', 'lib0', 'lib1', 'lib1']
    values = result.X.toarray() if sp.issparse(result.X) else result.X
    raw_values = result.raw.X.toarray() if sp.issparse(result.raw.X) else result.raw.X
    np.testing.assert_array_equal(values, [[1, 2, 0], [3, 4, 0], [0, 1, 2], [0, 3, 4]])
    np.testing.assert_array_equal(raw_values, values)
    layer = result.layers['counts']
    np.testing.assert_array_equal(layer.toarray() if sp.issparse(layer) else layer, values)
    np.testing.assert_array_equal(result.obsp['neighbors'].toarray(),
                                  [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 2], [0, 0, 2, 0]])


def test_min_cells_uses_global_counts_after_min_genes(tmp_path):
    source = ad.AnnData(sp.csr_matrix([[1, 0, 1], [1, 0, 1], [0, 1, 0]], dtype=np.float32),
                       var=pd.DataFrame(index=['a', 'b', 'c']))
    result = Workspace(tmp_path).qc(source, 2, 2, 0, None)
    assert list(result.obs_names) == ['0', '1']
    assert list(result.var_names) == ['a', 'c']
    np.testing.assert_array_equal(result.X.toarray(), [[1, 1], [1, 1]])


def test_header_reports_uncompressed_matrices(tmp_path):
    source = ad.AnnData(sp.eye(20, format='csr', dtype=np.float64))
    source.layers['counts'] = source.X.copy()
    path = tmp_path / 'input.h5ad'
    source.write_h5ad(path, compression='gzip')
    info = inspect_h5ad(path)
    assert info['cells'] == info['genes'] == 20
    assert info['matrix_bytes'] == 2 * info['x_bytes']


def test_mapped_marker_outputs_match_standard_hook(tmp_path):
    from altanalyze3.components.visualization.marker_heatmap_h5ad import generate_marker_heatmap_from_adata
    from altanalyze3.components.cellHarmony.cellHarmony_lite import normalize_adata
    rng = np.random.default_rng(712)
    counts = rng.poisson(3, (240, 300)).astype(np.float32)
    counts[:120, :5] += 25
    counts[120:, 5:10] += 25
    source = ad.AnnData(sp.csr_matrix(counts),
              obs=pd.DataFrame({'state': ['A'] * 120 + ['B'] * 120}, index=[f'c{i}' for i in range(240)]),
              var=pd.DataFrame(index=[f'G{i}' for i in range(300)]))
    normalize_adata(source)
    path = tmp_path / 'markers.h5ad'
    source.write_h5ad(path)
    # Exercise the mixed-index-width input produced by some H5AD writers.
    source.X.indptr = source.X.indptr.astype(np.int64)
    outputs = []
    for name, obj in [('base', source), ('mapped', Workspace(tmp_path / 'maps').load(path))]:
        outputs.append(generate_marker_heatmap_from_adata(obj, cluster_key='state', marker_method='markerfinder',
                       out=str(tmp_path / f'{name}.pdf'), render_heatmap=False,
                       write_heatmap_tsv=False, write_expression_tsv=False, write_heatmap_cache=True))
    for key in ('markers_tsv', 'centroids_tsv'):
        left, right = [pd.read_csv(result[key], sep='\t') for result in outputs]
        pd.testing.assert_frame_equal(left, right, check_exact=False, rtol=1e-5, atol=1e-6)


@pytest.mark.parametrize('rho', [0.2, 'auto'])
def test_ambient_preserves_raw_with_mapped_work_arrays(tmp_path, rho):
    from altanalyze3.components.ambient_rna.ambient_subtract import process_anndata
    rng = np.random.default_rng(17)
    obj = ad.AnnData(sp.csr_matrix(rng.poisson(3, (140, 70)).astype(np.float32)),
                    obs=pd.DataFrame({'Library': ['a'] * 70 + ['b'] * 70}, index=[f'c{i}' for i in range(140)]))
    path = tmp_path / 'source.h5ad'
    obj.write_h5ad(path)
    options = dict(rho=rho, inplace=True, write_individual=False, write_merged=False, store_corrected_layer=False)
    expected = process_anndata(obj.copy(), outdir=tmp_path / 'base', **options)
    actual = process_anndata(Workspace(tmp_path / 'work').load(path), outdir=tmp_path / 'mapped', **options)
    np.testing.assert_array_equal(actual.X.toarray(), expected.X.toarray())
    np.testing.assert_array_equal(actual.layers['soupx_raw'].toarray(), obj.X.toarray())


def test_upload_multiple_h5ad_and_memory_summary(tmp_path):
    from fastapi.testclient import TestClient
    from altanalyze3.components.cellHarmony.webapp.app import create_app
    from altanalyze3.components.cellHarmony.flask.pipeline import _split_uploads, _upload_profile
    source = tmp_path / 'query.h5ad'
    ad.AnnData(sp.eye(3, format='csr')).write_h5ad(source)
    app = create_app({'JOB_STORAGE': str(tmp_path / 'jobs')})
    with TestClient(app) as client:
        resp = client.post('/api/jobs', data={'species': 'human', 'reference': 'demo',
                           'sample_names': ['a', 'b']}, files=[('files', ('one.h5ad', source.read_bytes())),
                                                              ('files', ('two.h5ad', source.read_bytes()))])
    assert resp.status_code == 200, resp.text
    assert resp.json()['input_memory']['bounded_h5ad']
    assert resp.json()['input_memory']['cells'] == 6
    meta = app.state.job_store.get_job(resp.json()['job_id'])
    h5, h5ad = _split_uploads(meta['files'], app.state.job_store.uploads_dir(meta['job_id']))
    assert not h5 and len(h5ad) == 2
    assert _upload_profile(meta)['differential_eligible']


def test_memory_guard_terminates_real_child(tmp_path, monkeypatch):
    import os
    import subprocess
    import sys
    from altanalyze3.components.cellHarmony.flask import tasks
    from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
    runner = tasks.JobRunner(JobStore(tmp_path / 'jobs'), tmp_path / 'registry.json')
    process = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(30)'], start_new_session=True)
    runner._active_workers[process.pid] = process
    monkeypatch.setattr(tasks, 'container_memory', lambda: (28 * 1024**3, 30 * 1024**3))
    try:
        with pytest.raises(MemoryError, match='safety limit'):
            runner._wait_for_worker(process)
        assert process.poll() is not None
    finally:
        if process.poll() is None:
            process.kill()
        process.wait()
        runner.executor.shutdown()


def test_four_h5ads_complete_pipeline_and_expose_views(tmp_path):
    import json
    import shutil
    from test_flask_pipeline import _write_reference, _write_query_h5ad
    from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
    from altanalyze3.components.cellHarmony.flask.tasks import JobRunner
    from altanalyze3.components.cellHarmony.webapp.app import create_app
    from fastapi.testclient import TestClient
    reference = _write_reference(tmp_path)
    registry = tmp_path / 'registry.json'
    registry.write_text(json.dumps({'species': [{'id': 'demo_species', 'label': 'Demo',
        'references': [{'id': 'demo_reference', 'label': 'Demo', **reference}]}]}))
    store = JobStore(tmp_path / 'jobs')
    job = store.create_job('demo_species', 'demo_reference', None, [])['job_id']
    source = _write_query_h5ad(tmp_path)
    records = []
    for i in range(4):
        name = f'input{i}.h5ad'
        shutil.copyfile(source, store.uploads_dir(job) / name)
        records.append({'filename': name, 'sample_name': f'upload{i}'})
    store.update_job(job, files=records, qc={'min_genes': 1, 'min_cells': 0,
                     'min_counts': 0, 'mit_percent': 50, 'align_cutoff': 0})
    runner = JobRunner(store, registry, export_approx_pdfs=False)
    try:
        runner._run_pipeline(job)
    finally:
        runner.executor.shutdown()
    meta = store.get_job(job)
    assert meta['status'] == 'completed', meta.get('message')
    result = ad.read_h5ad(meta['artifacts']['combined_h5ad'])
    assert result.n_obs == 24 and result.n_vars == 3
    assert result.obs_names.is_unique
    assert set(result.obs['scalable_upload']) == {f'upload{i}' for i in range(4)}
    assert result.obsm['X_umap'].shape == (24, 2)
    options = meta['differential_options']
    assert options['default_sample_field'] == 'scalable_upload'
    assert set(options['sample_values']['scalable_upload']) == {f'upload{i}' for i in range(4)}
    assert not (store.outputs_dir(job) / '.h5ad_work').exists()
    app = create_app({'JOB_STORAGE': str(store.root), 'REFERENCE_REGISTRY': str(registry)})
    with TestClient(app) as client:
        for suffix in ('status', 'umap', 'expression?gene=GeneA'):
            response = client.get(f'/api/jobs/{job}/{suffix}')
            assert response.status_code == 200, response.text


@pytest.mark.parametrize('encoding', ['csr', 'dense'])
def test_scale_probe_matches_original_random_sample(tmp_path, encoding):
    from altanalyze3.components.cellHarmony.cellHarmony_lite import assess_expression_scale
    rng = np.random.default_rng(456)
    values = rng.poisson(3, (25001, 13)).astype(np.float32)
    obj = ad.AnnData(sp.csr_matrix(values) if encoding == 'csr' else values)
    obj.layers['counts'] = obj.X.copy()
    path = tmp_path / 'scale.h5ad'
    obj.write_h5ad(path)
    loaded = Workspace(tmp_path / 'work').load(path)
    assert assess_expression_scale(loaded) == assess_expression_scale(obj)
    # Releasing dead scratch arrays must not invalidate views retained by AnnData.
    import gc
    gc.collect()
    loaded._matrix_workspace.release_pages()
    saved = tmp_path / 'saved.h5ad'
    loaded.write_h5ad(saved)
    result = ad.read_h5ad(saved)
    matrix = result.X.toarray() if sp.issparse(result.X) else result.X
    np.testing.assert_array_equal(matrix, values)


def test_merge_rejects_mixed_scales_before_creating_output(tmp_path):
    entries = []
    for name, matrix in [('counts', [[1., 2.], [2., 3.]]),
                         ('logged', [[1.2, 1.9], [1.8, 2.1]])]:
        path = tmp_path / (name + '.h5ad')
        ad.AnnData(sp.csr_matrix(matrix)).write_h5ad(path)
        entries.append((path, name))
    with pytest.raises(ValueError, match='different expression scales'):
        merge_h5ads(entries, tmp_path / 'merge')
    assert not (tmp_path / 'merge' / 'merged.h5ad').exists()


def test_supervised_memory_failure_records_status_and_cleans_scratch(tmp_path, monkeypatch):
    import subprocess
    import sys
    from altanalyze3.components.cellHarmony.flask import tasks
    from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
    store = JobStore(tmp_path / 'jobs')
    job = store.create_job('human', 'reference', None, [])['job_id']
    upload = store.uploads_dir(job) / 'keep.h5ad'
    upload.write_bytes(b'unchanged input fixture')
    scratch = store.outputs_dir(job) / '.h5ad_work'
    scratch.mkdir(parents=True)
    (scratch / 'array.bin').write_bytes(b'scratch')
    runner = tasks.JobRunner(store, tmp_path / 'registry.json')
    original_popen = subprocess.Popen
    def child(*args, **kwargs):
        return original_popen([sys.executable, '-c', 'import time; time.sleep(30)'], **kwargs)
    monkeypatch.setattr(tasks.subprocess, 'Popen', child)
    monkeypatch.setattr(runner, '_memory_wait_reason', lambda: None)
    monkeypatch.setattr(tasks, 'container_memory', lambda: (28 * 1024**3, 30 * 1024**3))
    # Avoid the process-table subprocess being replaced by the worker stub.
    monkeypatch.setattr(tasks, 'process_memory', lambda: {})
    try:
        runner._run_isolated(job, 'pipeline')
    finally:
        runner.executor.shutdown()
    meta = store.get_job(job)
    assert meta['status'] == 'failed'
    assert 'safety limit' in meta['message']
    assert not scratch.exists()
    assert upload.read_bytes() == b'unchanged input fixture'


def test_reuploaded_job_preserves_existing_upload_annotation(tmp_path):
    entries = []
    for i in range(2):
        obj = ad.AnnData(sp.csr_matrix([[1., 2.]]),
                        obs=pd.DataFrame({'scalable_upload': ['original'],
                                          'source_scalable_upload': ['older']}, index=['cell']))
        path = tmp_path / f'{i}.h5ad'
        obj.write_h5ad(path)
        entries.append((path, f'new{i}'))
    result = ad.read_h5ad(merge_h5ads(entries, tmp_path / 'merge'))
    assert list(result.obs['scalable_upload']) == ['new0', 'new1']
    assert list(result.obs['source_scalable_upload']) == ['older', 'older']
    assert list(result.obs['source_source_scalable_upload']) == ['original', 'original']


def test_conflicting_feature_identity_stops_merge():
    from altanalyze3.components.cellHarmony.mapped_h5ad import _merge_feature_annotations
    tables = [pd.DataFrame({'gene_symbols':['A']}, index=['ENSG1']),
              pd.DataFrame({'gene_symbols':['B']}, index=['ENSG1'])]
    with pytest.raises(ValueError, match='conflicting gene_symbols'):
        _merge_feature_annotations(tables, pd.Index(['ENSG1']))
