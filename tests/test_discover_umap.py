from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from altanalyze3.components.clustering import ICGS
from altanalyze3.components.clustering.umap_fit import fit_umap, select_landmarks


def test_core_defaults_remain_full_and_cli_can_opt_in():
    config = ICGS.ICGS3Config(input_paths=[], output_dir="unused")
    assert config.umap_fit_mode == "full"
    parser = ICGS.build_arg_parser()
    args = parser.parse_args(["--input", "demo.h5ad", "--output-dir", "unused"])
    assert args.umap_fit_mode == "full"
    config.umap_fit_mode = "landmark"
    assert "--umap-fit-mode landmark" in ICGS.cli_equivalent(config)


def test_landmarks_preserve_small_states_and_every_cell_in_order():
    labels = np.array(["large"] * 1700 + ["small"] * 280 + ["rare"] * 2)
    x = np.arange(len(labels) * 3, dtype=np.float32).reshape(len(labels), 3)
    batches = []
    class Model:
        n_neighbors = 10
        def fit_transform(self, block):
            batches.append(("fit", block.copy()))
            return block[:, :2]
        def transform(self, block):
            batches.append(("transform", block.copy()))
            return block[:, :2]
    coords, info, selected = fit_umap(Model(), lambda rows: x[rows], len(x),
        mode="landmark", labels=labels, max_fit_cells=1000, batch_cells=31, seed=42)
    np.testing.assert_array_equal(coords, x[:, :2])
    assert {1980, 1981}.issubset(selected)
    assert np.sum(labels[selected] == "small") >= 200
    assert np.sum(labels[selected] == "large") >= 200
    assert len(selected) == 1000 and info["fit_mode"] == "landmark"
    assert info["min_landmarks_per_state"] == 200
    assert all(len(block) <= 31 for kind, block in batches if kind == "transform")
    np.testing.assert_array_equal(selected, select_landmarks(labels, 1000, 42))
    assert np.concatenate([block[:, 0] for _, block in batches]).size == len(x)


def test_landmark_budget_expands_for_complete_state_coverage():
    labels = np.repeat(np.arange(6).astype(str), 200)
    selected = select_landmarks(labels, 100, 0)
    assert len(selected) == len(labels)


def test_stratified_budget_fills_exactly_and_keeps_every_rare_cell():
    labels = np.repeat(["large", "medium", "rare"], [40000, 12000, 199])
    selected = select_landmarks(labels, 30000, 0)
    assert len(selected) == len(np.unique(selected)) == 30000
    assert set(np.flatnonzero(labels == "rare")).issubset(selected)
    assert all(np.sum(labels[selected] == s) >= min(200, np.sum(labels == s))
               for s in np.unique(labels))
    np.testing.assert_array_equal(selected, select_landmarks(labels, 30000, 0))


def test_small_jobs_fit_all_cells_and_nonfinite_coordinates_fail():
    x = np.ones((20, 3), dtype=np.float32)
    model = SimpleNamespace(fit_transform=lambda block: block[:, :2])
    coords, info, selected = fit_umap(model, lambda rows: x[rows], 20, mode="landmark")
    assert info["fit_mode"] == "full" and len(selected) == 20
    model.fit_transform = lambda block: np.full((len(block), 2), np.nan)
    with pytest.raises(ValueError, match="finite coordinates"):
        fit_umap(model, lambda rows: x[rows], 20)


def test_landmark_reads_only_bounded_blocks_and_preserves_features(tmp_path, monkeypatch):
    from altanalyze3.components.clustering import umap_input
    x = np.random.default_rng(19).normal(size=(600, 8))
    source = ad.AnnData(sp.csr_matrix(x), obs=pd.DataFrame(index=[f"C{i}" for i in range(600)]),
                       var=pd.DataFrame(index=[f"G{i}" for i in range(8)]))
    target = ad.AnnData(sp.csr_matrix((600, 8)), obs=source.obs.iloc[::-1].copy(), var=source.var.copy())
    target.obs["ICGS3_cluster"] = ["C1"] * 600
    features = ["G6", "G1", "G3"]
    original = umap_input.dense_umap_input
    calls = []
    def bounded(*args, **kwargs):
        calls.append(len(kwargs["cells"]))
        assert calls[-1] <= 220
        return original(*args, **kwargs)
    monkeypatch.setattr(umap_input, "dense_umap_input", bounded)
    class Model:
        def __init__(self, **kwargs): self.n_neighbors = kwargs["n_neighbors"]
        def fit_transform(self, block): return block[:, :2]
        def transform(self, block): return block[:, :2]
    monkeypatch.setattr(ICGS, "_import_umap_with_local_retry", lambda: SimpleNamespace(UMAP=Model))
    config = ICGS.ICGS3Config(input_paths=[], output_dir=str(tmp_path), minimal_outputs=True,
                            umap_fit_mode="landmark", umap_fit_cells=220, umap_transform_batch_cells=30)
    ICGS.compute_umap_outputs(target, config, str(tmp_path), marker_features=features, matrix_source=source)
    np.testing.assert_array_equal(target.obsm["X_umap"], x[::-1][:, [6, 1]].astype(np.float32))
    assert list(pd.read_csv(tmp_path / "UMAPs/icgs3_umap_features.tsv", sep="\t").feature) == features
    landmarks = pd.read_csv(tmp_path / "UMAPs/icgs3_umap_landmarks.tsv", sep="\t").barcode
    assert len(landmarks) == 220 and set(landmarks).issubset(target.obs_names)


@pytest.mark.parametrize("holds_counts", [False, True])
@pytest.mark.parametrize("fit_mode", [None, "full"])
@pytest.mark.parametrize("max_k", [None, 12])
def test_discover_handoff_preserves_scale_and_passes_option(tmp_path, monkeypatch, holds_counts, fit_mode, max_k):
    from altanalyze3.components.cellHarmony.scalable_discover import pipeline
    x = np.random.default_rng(1).integers(0, 100, size=(40, 10)).astype(np.float32)
    normalized = np.log1p(x)
    a = ad.AnnData(sp.csr_matrix(normalized))
    if holds_counts: a.layers["counts"] = sp.csr_matrix(x)
    if not holds_counts:
        # The interface fix must not weaken ICGS3's scientific scale guard.
        with pytest.raises(ICGS.ExpressionScaleError, match="input-normalized"):
            ICGS.report_expression_scale(a, ICGS.ICGS3Config(input_paths=[], output_dir=str(tmp_path)))
    qc = {} if fit_mode is None else {"umap_fit_mode": fit_mode}
    if max_k is not None:
        qc["max_k"] = max_k
    a.write_h5ad(tmp_path / "source.h5ad")
    store = SimpleNamespace(get_job=lambda _: {"species": "human", "reference": "icgs3",
        "files": [{"filename": "source.h5ad"}], "qc": qc},
        uploads_dir=lambda _: tmp_path, outputs_dir=lambda _: tmp_path / "outputs",
        update_job=lambda *a, **kw: None, append_log=lambda *a: None)
    monkeypatch.setattr(pipeline.cellHarmony_lite, "combine_and_align_h5", lambda **kw: (None, a))
    class ReachedICGS(Exception): pass
    def capture(config):
        staged = ad.read_h5ad(config.input_paths[0])
        np.testing.assert_array_equal(staged.X.toarray(), x if holds_counts else normalized)
        assert config.input_normalized is not holds_counts
        assert config.rank == max_k
        assert config.max_auto_nmf_k is None
        if max_k is not None:
            assert "--nmf-k 12" in ICGS.cli_equivalent(config)
        else:
            assert "--nmf-k" not in ICGS.cli_equivalent(config)
        assert config.umap_fit_mode == ("full" if fit_mode == "full" else "landmark")
        assert config.umap_n_neighbors == (0 if fit_mode == "full" else 15)
        assert config.umap_transform_backend == ("umap" if fit_mode == "full" else "exact_correlation")
        assert config.umap_min_dist == 0.75
        report = ICGS.report_expression_scale(staged, config)
        assert report["X"]["verdict"] == ("counts" if holds_counts else "log")
        assert ICGS.resolve_normalization(config) == ("cp10k-log1p" if holds_counts else "none")
        if not holds_counts:
            prepared = ICGS.prepare_expression(staged, config)
            np.testing.assert_array_equal(prepared.X.toarray(), normalized)
            assert "counts" not in prepared.layers
        raise ReachedICGS
    monkeypatch.setattr(ICGS, "run_icgs3", capture)
    with pytest.raises(ReachedICGS): pipeline.run_discover_pipeline("test", store)


@pytest.mark.parametrize("matrix_bytes, expected", [(1024**3 - 1, False), (1024**3, True)])
def test_discover_uses_advertised_disk_import_below_cell_threshold(tmp_path, monkeypatch, matrix_bytes, expected):
    from altanalyze3.components.cellHarmony.scalable_discover import pipeline
    from altanalyze3.components.cellHarmony import mapped_h5ad
    store = SimpleNamespace(get_job=lambda _: {"species": "human", "reference": "icgs3",
        "files": [{"filename": "source.h5ad"}], "qc": {}},
        uploads_dir=lambda _: tmp_path, outputs_dir=lambda _: tmp_path / "outputs",
        update_job=lambda *a, **kw: None, append_log=lambda *a: None)
    monkeypatch.setattr(mapped_h5ad, "inspect_h5ad", lambda _: {"cells": 65_662, "matrix_bytes": matrix_bytes})
    class ReachedImport(Exception): pass
    def capture(**kwargs):
        assert kwargs["bounded_h5ad"] is expected
        assert kwargs["min_genes"] == 500 and kwargs["min_counts"] == 1000
        raise ReachedImport
    monkeypatch.setattr(pipeline.cellHarmony_lite, "combine_and_align_h5", capture)
    with pytest.raises(ReachedImport): pipeline.run_discover_pipeline("test", store)


def test_discover_multiple_h5ads_preserve_union_samples_and_count_handoff(tmp_path, monkeypatch):
    from altanalyze3.components.cellHarmony.scalable_discover import pipeline
    records = []
    for i, genes in enumerate((["a", "b"], ["b", "c"])):
        values = np.array([[1, 2], [3, 4]], dtype=np.float32) + i * 4
        source = ad.AnnData(sp.csr_matrix(values), obs=pd.DataFrame(index=["cell1", "cell2"]),
                           var=pd.DataFrame(index=genes))
        source.write_h5ad(tmp_path / f"{i}.h5ad")
        records.append({"filename": f"{i}.h5ad", "sample_name": f"lib{i}"})
    store = SimpleNamespace(get_job=lambda _: {"species": "human", "reference": "icgs3",
        "files": records, "qc": {"min_genes": 0, "min_cells": 0, "min_counts": 0}},
        uploads_dir=lambda _: tmp_path, outputs_dir=lambda _: tmp_path / "outputs",
        update_job=lambda *a, **kw: None, append_log=lambda *a: None)
    class ReachedICGS(Exception): pass
    def capture(config):
        staged = ad.read_h5ad(config.input_paths[0])
        assert list(staged.var_names) == ["a", "b", "c"]
        assert list(staged.obs_names) == ["cell1::lib0", "cell2::lib0", "cell1::lib1", "cell2::lib1"]
        assert list(staged.obs["Library"]) == ["lib0", "lib0", "lib1", "lib1"]
        np.testing.assert_array_equal(staged.X.toarray(), [[1, 2, 0], [3, 4, 0], [0, 5, 6], [0, 7, 8]])
        assert config.input_normalized is False
        raise ReachedICGS
    monkeypatch.setattr(ICGS, "run_icgs3", capture)
    with pytest.raises(ReachedICGS): pipeline.run_discover_pipeline("test", store)


def test_discover_umap_progress_advances_without_changing_analysis():
    from altanalyze3.components.cellHarmony.scalable_discover.tasks import _StageLogStream
    updates = []
    store = SimpleNamespace(append_log=lambda *a: None,
                            update_job=lambda job, **fields: updates.append(fields))
    stream = _StageLogStream(store, "job")
    for line in ["running final UMAP", "UMAP fit mode=landmark", "UMAP mapping 32,498 remaining cells",
                 "UMAP fit completed in 120.0 s; transform in 60.0 s"]:
        stream.write(line + "\n")
    assert [item["progress"] for item in updates] == [68, 69, 72, 75]
    assert "mapping remaining cells" in updates[2]["message"]


def test_discover_ui_and_qc_api_round_trip(tmp_path):
    from fastapi.testclient import TestClient
    from altanalyze3.components.cellHarmony.scalable_discover.app import create_discover_app
    app = create_discover_app({"JOB_STORAGE": str(tmp_path), "TESTING": True})
    job_dir = tmp_path / "test"
    job_dir.mkdir()
    app.state.job_store._write_metadata("test", {"job_id": "test", "status": "uploaded"})
    try:
        with TestClient(app) as client:
            page = client.get("/")
            assert page.status_code == 200
            assert '<option value="landmark" selected>' in page.text
            assert '<span>Max K</span>' in page.text
            assert 'type="number" min="2" step="1"' in page.text
            assert 'placeholder="—"' in page.text
            for choice in ["full", "landmark"]:
                response = client.post("/api/jobs/test/qc", json={"umap_fit_mode": choice})
                assert response.status_code == 200 and response.json()["qc"]["umap_fit_mode"] == choice
                assert app.state.job_store.get_job("test")["qc"]["umap_fit_mode"] == choice
            assert client.post("/api/jobs/test/qc", json={"umap_fit_mode": "unknown"}).status_code == 422
            assert client.post("/api/jobs/test/qc", json={}).json()["qc"]["umap_fit_mode"] == "landmark"
            assert app.state.job_store.get_job("test")["qc"]["max_k"] is None
            for max_k in [2, 12, 100, None]:
                response = client.post("/api/jobs/test/qc", json={"max_k": max_k})
                assert response.status_code == 200
                assert response.json()["qc"]["max_k"] == max_k
                assert app.state.job_store.get_job("test")["qc"]["max_k"] == max_k
            for invalid in [0, 1, -2, 12.5, True, "twelve"]:
                assert client.post("/api/jobs/test/qc", json={"max_k": invalid}).status_code == 422
                assert app.state.job_store.get_job("test")["qc"]["max_k"] is None
    finally:
        app.state.job_runner.executor.shutdown(wait=True)


def test_projection_logs_all_cells_and_timings_without_changing_outputs():
    x = np.arange(2000 * 3, dtype=np.float32).reshape(2000, 3)
    labels = np.array(['c1'] * 1000 + ['c2'] * 1000)
    class Model:
        n_neighbors = 15
        def fit_transform(self, values): return values[:, :2]
        def transform(self, values): return values[:, :2]
    logs = []
    coords, _, selected = fit_umap(Model(), lambda rows: x[rows], len(x), mode='landmark',
                                   labels=labels, max_fit_cells=500, batch_cells=500, log=logs.append)
    np.testing.assert_array_equal(coords, x[:, :2])
    assert len(selected) == 500
    progress = [line for line in logs if line.startswith('UMAP transformed')]
    assert len(progress) == 3
    assert '(500/1,500 completed;' in progress[0]
    assert '(1,500/1,500 completed;' in progress[-1]
    assert all('input ' in line and 'projection ' in line for line in progress)
    assert len([line for line in logs if line.startswith('UMAP projection input loaded:')]) == 3


def test_discover_projection_progress_changes_message_for_each_completed_block():
    from altanalyze3.components.cellHarmony.scalable_discover.tasks import _StageLogStream
    updates = []
    store = SimpleNamespace(append_log=lambda *a: None,
                            update_job=lambda job, **fields: updates.append(fields))
    stream = _StageLogStream(store, 'job')
    for line in ['UMAP mapping 42,083 remaining cells',
                 'UMAP transformed 42,083 remaining cells (42,083/210,414 completed; input 1.0 s; projection 2.0 s)',
                 'UMAP transformed 42,083 remaining cells (84,166/210,414 completed; input 1.0 s; projection 2.0 s)',
                 'UMAP transformed 42,082 remaining cells (210,414/210,414 completed; input 1.0 s; projection 2.0 s)',
                 'UMAP fit completed in 45.6 s; transform in 1427.1 s']:
        stream.write(line+'\n')
    assert [update['progress'] for update in updates] == [72, 72, 73, 75, 75]
    assert '42,083 of 210,414' in updates[1]['message']
    assert 'writing UMAP' in updates[-1]['message']
    count = len(updates)
    stream.write('UMAP transformed 1 remaining cells (211,000/210,414 completed; invalid)\n')
    assert len(updates) == count


def test_exact_projection_search_matches_full_panel_correlation_and_bounds_tiles():
    from altanalyze3.components.clustering.umap_neighbors import ExactCorrelationIndex
    from scipy.spatial.distance import cdist
    rng = np.random.default_rng(53)
    train = rng.normal(size=(73, 93)).astype(np.float32)
    query = rng.normal(size=(41, 93)).astype(np.float32)
    before = train.copy(), query.copy()
    search = ExactCorrelationIndex(train, working_memory_bytes=73 * 12 * 3)
    assert search.block_rows == 3
    indices, distances = search.query(query, k=11)
    reference = cdist(query.astype(np.float64), train.astype(np.float64), metric='correlation')
    np.testing.assert_array_equal(indices, np.argsort(reference, axis=1)[:, :11])
    np.testing.assert_allclose(distances, np.take_along_axis(reference, indices, axis=1), atol=3e-7)
    np.testing.assert_array_equal(train, before[0])
    np.testing.assert_array_equal(query, before[1])


def test_exact_correlation_preserves_negative_correlations_and_constant_rules():
    from altanalyze3.components.clustering.umap_neighbors import ExactCorrelationIndex
    train = np.array([[0, 1, 2], [2, 1, 0], [4, 4, 4], [0, 0, 0]], dtype=np.float32)
    index = ExactCorrelationIndex(train)
    indices, distances = index.query(train, k=4)
    matrix = np.empty((4, 4))
    np.put_along_axis(matrix, indices, distances, axis=1)
    np.testing.assert_allclose(matrix, [[0, 2, 1, 1], [2, 0, 1, 1],
                                       [1, 1, 0, 0], [1, 1, 0, 0]], atol=3e-7)
    with pytest.raises(ValueError, match='complete feature panel'):
        index.query(np.ones((3, 2)), k=2)
    with pytest.raises(ValueError, match='finite'):
        index.query(np.array([[np.nan, 0, 0]]), k=2)


def test_projection_replacement_restores_fitted_index_on_success_and_error():
    from altanalyze3.components.clustering.umap_fit import projection_search
    from altanalyze3.components.clustering.umap_neighbors import ExactCorrelationIndex
    original = object()
    model = SimpleNamespace(metric='correlation', _small_data=False,
                            _knn_search_index=original, _raw_data=np.eye(7, dtype=np.float32))
    with projection_search(model, 'exact_correlation') as backend:
        assert backend == 'exact_correlation'
        assert isinstance(model._knn_search_index, ExactCorrelationIndex)
        assert model._knn_search_index._data.shape == model._raw_data.shape
    assert model._knn_search_index is original
    with pytest.raises(RuntimeError, match='deliberate'):
        with projection_search(model, 'exact_correlation'):
            raise RuntimeError('deliberate')
    assert model._knn_search_index is original
    model.metric = 'euclidean'
    with projection_search(model, 'exact_correlation') as backend:
        assert backend == 'umap'
        assert model._knn_search_index is original
    model.metric = 'correlation'
    model._sparse_data = True
    with projection_search(model, 'exact_correlation') as backend:
        assert backend == 'umap'
        assert model._knn_search_index is original


def test_full_fit_keeps_established_search_even_when_projection_backend_requested():
    x = np.arange(30, dtype=np.float32).reshape(10, 3)
    class Model:
        def fit_transform(self, rows):
            return rows[:, :2]
        def transform(self, rows):
            raise AssertionError('Full fit must not project cells')
    coords, info, selected = fit_umap(Model(), lambda rows: x[rows], len(x),
                                     transform_backend='exact_correlation')
    np.testing.assert_array_equal(coords, x[:, :2])
    np.testing.assert_array_equal(selected, np.arange(len(x)))
    assert info['transform_search_backend'] == 'umap'
    config = ICGS.ICGS3Config(input_paths=[], output_dir='unused')
    assert config.umap_transform_backend == 'umap'
    config.umap_transform_backend = 'exact_correlation'
    assert '--umap-transform-backend exact_correlation' in ICGS.cli_equivalent(config)


def test_landmark_projection_installs_exact_search_and_preserves_complete_roster():
    from altanalyze3.components.clustering.umap_neighbors import ExactCorrelationIndex
    rng = np.random.default_rng(8)
    x = rng.normal(size=(2002, 17)).astype(np.float32)
    labels = np.array(['large'] * 1700 + ['small'] * 300 + ['rare'] * 2)
    original = object()
    class Model:
        n_neighbors = 15
        metric = 'correlation'
        _small_data = False
        _knn_search_index = original
        calls = 0
        def fit_transform(self, rows):
            self._raw_data = rows.copy()
            return rows[:, :2]
        def transform(self, rows):
            assert isinstance(self._knn_search_index, ExactCorrelationIndex)
            assert self._knn_search_index._data.shape == self._raw_data.shape
            positions, distances = self._knn_search_index.query(rows, self.n_neighbors)
            assert positions.shape == distances.shape == (len(rows), self.n_neighbors)
            self.calls += 1
            return rows[:, :2]
    model = Model()
    coords, info, selected = fit_umap(model, lambda rows: x[rows], len(x),
        mode='landmark', labels=labels, max_fit_cells=600, batch_cells=501,
        transform_backend='exact_correlation')
    assert info['transform_search_backend'] == 'exact_correlation'
    assert model._knn_search_index is original and model.calls == 3
    np.testing.assert_array_equal(coords, x[:, :2])
    np.testing.assert_array_equal(selected, select_landmarks(labels, 600, 0, 200))
