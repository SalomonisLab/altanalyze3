"""Regression checks for correction memory options and orphaned pipeline jobs."""
import os
import threading
from concurrent.futures import Future

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from altanalyze3.components.ambient_rna import ambient_subtract as ambient


def _legacy_subtract(counts, profile, rho, *, round_to_int=False):
    """Independent reference for the former CSC algorithm."""
    csc = counts.tocsc(copy=True)
    csc.data = csc.data.astype(np.float32, copy=False)
    totals = np.asarray(counts.sum(axis=1)).ravel().astype(np.float32)
    profile = np.asarray(profile, dtype=np.float32)
    for gene in range(csc.shape[1]):
        lo, hi = csc.indptr[gene:gene + 2]
        values = csc.data[lo:hi] - rho * profile[gene] * totals[csc.indices[lo:hi]]
        values[values < 0] = 0
        csc.data[lo:hi] = np.rint(values) if round_to_int else values
    csc.eliminate_zeros()
    return csc.tocsr()


@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.int32])
@pytest.mark.parametrize("rho", [0.0, 0.35, 1.0])
def test_csr_kernel_matches_legacy_for_input_dtypes(dtype, rho):
    rng = np.random.default_rng(12)
    values = rng.random((9, 30)) * 30
    values[0] = 0
    counts = sp.csr_matrix(values.astype(dtype))
    profile = rng.random(30)
    profile /= profile.sum()
    expected = _legacy_subtract(counts, profile, rho)
    actual = ambient._correct_counts_csc_subtraction(counts, profile, rho)
    np.testing.assert_array_equal(actual.toarray(), expected.toarray())


@pytest.mark.parametrize("rounded", [False, True])
@pytest.mark.parametrize("interleaved", [False, True])
def test_six_libraries_match_legacy_correction_and_auto_rho(tmp_path, monkeypatch, rounded, interleaved):
    rng = np.random.default_rng(21)
    counts = sp.random(600, 120, density=0.2, random_state=rng, format="csr",
                       data_rvs=lambda n: rng.integers(1, 30, n).astype(np.float32), dtype=np.float32)
    libs = np.tile(np.arange(6), 100) if interleaved else np.repeat(np.arange(6), 100)
    source = ad.AnnData(counts.copy(), obs=pd.DataFrame({"Library": libs}, index=[f"c{i}" for i in range(600)]))
    expected = np.zeros(counts.shape, dtype=np.float32)
    rhos = []
    for lib in range(6):
        rows = np.where(libs == lib)[0]
        subset = counts[rows].copy()
        profile, _ = ambient._estimate_ambient_profile_from_filtered(subset)
        with monkeypatch.context() as patch:
            patch.setattr(ambient, "_correct_counts_csc_subtraction", _legacy_subtract)
            rho, _ = ambient._auto_select_rho(subset, profile)
        rhos.append(rho)
        expected[rows] = _legacy_subtract(subset, profile, rho, round_to_int=rounded).toarray()
    result = ambient.process_anndata(source, rho="auto", round_to_int=rounded, inplace=True,
                                    store_corrected_layer=False, write_individual=False,
                                    write_merged=False, outdir=tmp_path)
    np.testing.assert_array_equal(result.X.toarray(), expected)
    np.testing.assert_array_equal(result.layers["soupx_raw"].toarray(), counts.toarray())
    np.testing.assert_array_equal(result.uns["soupx_correction"]["libraries"]["rho"], rhos)


def test_contiguous_library_buffers_share_memory():
    source = sp.csr_matrix(np.ones((600, 120), dtype=np.float32))
    subset = ambient._subset_rows_csr(source, np.arange(100))
    assert np.shares_memory(source.data, subset.data)
    assert np.shares_memory(source.indices, subset.indices)


def test_duplicate_sparse_coordinates_match_legacy(tmp_path):
    counts = sp.csr_matrix((np.array([2, 3, 4, 1] * 4, dtype=np.float32),
                            np.array([0, 0, 1, 2] * 4), np.arange(0, 17, 4)), shape=(4, 3))
    source = ad.AnnData(counts.copy(), obs=pd.DataFrame({"Library": ["a", "a", "b", "b"]},
                                                      index=[f"c{i}" for i in range(4)]))
    expected = np.zeros(counts.shape, dtype=np.float32)
    for rows in (np.array([0, 1]), np.array([2, 3])):
        subset = counts[rows].copy()
        profile, _ = ambient._estimate_ambient_profile_from_filtered(subset)
        expected[rows] = _legacy_subtract(subset, profile, 0.2).toarray()
    result = ambient.process_anndata(source, rho=0.2, inplace=True, store_corrected_layer=False,
                                    write_individual=False, write_merged=False, outdir=tmp_path)
    np.testing.assert_array_equal(result.X.toarray(), expected)


@pytest.mark.parametrize("sparse", [False, True])
def test_cosine_blocks_preserve_scores_and_ties(sparse):
    from altanalyze3.components.cellHarmony.cellHarmony_lite import _cosine_best_matches
    rng = np.random.default_rng(9)
    query = rng.random((23, 17)).astype(np.float32)
    query[0] = 0
    ref = rng.random((7, 17))
    ref[1] = ref[0]  # Preserve first-match tie behavior.
    if sparse:
        query = sp.csr_matrix(query)
        norms = np.sqrt(query.multiply(query).sum(axis=1)).A1
        norms[norms == 0] = 1
        normalized = sp.diags(1 / norms).dot(query)
    else:
        norms = np.linalg.norm(query, axis=1)
        norms[norms == 0] = 1
        normalized = query / norms[:, None]
    expected = np.asarray(normalized.dot((ref / np.linalg.norm(ref, axis=1)[:, None]).T))
    matches, scores = _cosine_best_matches(query, ref, chunk_size=4)
    np.testing.assert_array_equal(matches, expected.argmax(axis=1))
    np.testing.assert_allclose(scores, expected.max(axis=1), rtol=1e-7, atol=1e-7)


def test_alignment_subset_omits_unused_layers():
    from altanalyze3.components.cellHarmony.cellHarmony_lite import subset_to_reference_genes
    source = ad.AnnData(sp.csr_matrix(np.eye(5)), var=pd.DataFrame(index=[f"g{i}" for i in range(5)]))
    source.layers["counts"] = source.X.copy()
    selected, present, missing = subset_to_reference_genes(source, ["g3", "absent", "g1"], copy_layers=False)
    assert present == ["g3", "g1"] and missing == ["absent"]
    assert not selected.layers
    np.testing.assert_array_equal(selected.X.toarray(), source.X[:, [3, 1]].toarray())


@pytest.mark.parametrize("interleaved", [False, True])
@pytest.mark.parametrize("rho", [0.2, "auto"])
def test_owned_correction_matches_default(tmp_path, interleaved, rho):
    rng = np.random.default_rng(14)
    counts = sp.csr_matrix(rng.poisson(2, (120, 70)).astype(np.float32))
    libs = np.tile(["a", "b", "c"], 40) if interleaved else np.repeat(["a", "b", "c"], 40)
    source = ad.AnnData(counts, obs=pd.DataFrame({"Library": libs}, index=[f"c{i}" for i in range(120)]))
    expected = ambient.process_anndata(
        source, rho=rho, outdir=tmp_path / "default",
        write_individual=False, write_merged=False,
    )
    owned = source.copy()
    result = ambient.process_anndata(
        owned, rho=rho, outdir=tmp_path / "owned", inplace=True,
        store_corrected_layer=False, write_individual=False, write_merged=False,
    )
    assert result is owned
    assert "soupx_corrected" not in result.layers
    np.testing.assert_array_equal(source.X.toarray(), counts.toarray())
    np.testing.assert_array_equal(result.layers["soupx_raw"].toarray(), counts.toarray())
    np.testing.assert_allclose(result.X.toarray(), expected.X.toarray())
    pd.testing.assert_index_equal(result.obs_names, source.obs_names)
    cols = [c for c in expected.uns["soupx_correction"]["libraries"] if c != "total_time_sec"]
    pd.testing.assert_frame_equal(result.uns["soupx_correction"]["libraries"][cols],
                                  expected.uns["soupx_correction"]["libraries"][cols])
    # Later normalization must not modify the retained raw counts.
    result.X.data *= 2
    np.testing.assert_array_equal(result.layers["soupx_raw"].toarray(), counts.toarray())


def test_invalid_correction_options_do_not_mutate_input(tmp_path):
    source = ad.AnnData(sp.eye(3, format="csr"), obs=pd.DataFrame({"Library": ["a"] * 3}))
    with pytest.raises(ValueError, match="requires replace_x"):
        ambient.process_anndata(source, outdir=tmp_path, store_corrected_layer=False, replace_x=False)
    assert not source.layers


def test_orphaned_pipeline_fails_without_deleting_outputs(tmp_path):
    from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
    from altanalyze3.components.cellHarmony.flask.tasks import JobRunner

    store = JobStore(tmp_path)
    job = store.create_job("human", "ref", None, files=[])
    job_id = job["job_id"]
    output = store.outputs_dir(job_id) / "saved.txt"
    output.write_text("safe")
    runner = JobRunner(store, tmp_path / "registry.json", max_workers=1)
    try:
        # Same PID is deliberately used to simulate PID reuse after a Docker restart.
        store.update_job(job_id, status="queued", progress=15, worker_pid=os.getpid())
        recovered = runner.recover_interrupted_pipeline(job_id)
        assert recovered["status"] == "failed"
        assert "interrupted" in recovered["message"]
        assert output.read_text() == "safe"
        store.update_job(job_id, status="processing")
        pending = Future()
        runner._futures[f"pipeline:{job_id}"] = pending
        assert runner.recover_interrupted_pipeline(job_id)["status"] == "processing"
        pending.set_result(None)
        store.update_job(job_id, status="completed")
        assert runner.recover_interrupted_pipeline(job_id)["status"] == "completed"
    finally:
        runner.executor.shutdown(wait=True)


def test_submission_cannot_overwrite_fast_completion(tmp_path):
    from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
    from altanalyze3.components.cellHarmony.flask.tasks import JobRunner

    store = JobStore(tmp_path)
    job_id = store.create_job("human", "ref", None, files=[])["job_id"]
    runner = JobRunner(store, tmp_path / "registry.json", max_workers=1)
    runner._run_pipeline = lambda jid: store.update_job(jid, status="completed", progress=100)
    try:
        runner.submit(job_id)
        runner._futures[f"pipeline:{job_id}"].result(timeout=5)
        assert runner.recover_interrupted_pipeline(job_id)["status"] == "completed"
    finally:
        runner.executor.shutdown(wait=True)


@pytest.mark.parametrize("second_task", ["pipeline", "differential"])
def test_second_visitor_queues_without_false_failure(tmp_path, second_task):
    from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
    from altanalyze3.components.cellHarmony.flask.tasks import JobRunner
    store = JobStore(tmp_path)
    first = store.create_job("human", "ref", None, files=[])["job_id"]
    second = store.create_job("human", "ref", None, files=[])["job_id"]
    runner = JobRunner(store, tmp_path / "registry.json")
    started, release = threading.Event(), threading.Event()
    def run(job_id):
        store.update_job(job_id, status="processing")
        if job_id == first:
            started.set()
            assert release.wait(5)
        store.update_job(job_id, status="completed", message="Done.")
    runner._run_pipeline = run
    def run_differential(job_id):
        store.update_job(job_id, differential={"status": "completed", "message": "Done."})
    runner._run_differential = run_differential
    try:
        runner.submit(first)
        assert started.wait(5)
        if second_task == "pipeline":
            runner.submit(second)
        else:
            store.update_job(second, status="completed", differential={"config": {"modality": "rna"}})
            runner.submit_differential(second)
        assert runner.recover_interrupted_pipeline(first)["status"] == "processing"
        queued = store.get_job(second)
        if second_task == "differential":
            queued = queued["differential"]
            assert queued["config"] == {"modality": "rna"}
        else:
            assert runner.recover_interrupted_pipeline(second)["status"] == "queued"
        assert queued["status"] == "queued"
        assert "Other analyses" in queued["message"]
        assert "start automatically" in queued["message"]
        release.set()
        runner._futures[f"{second_task}:{second}"].result(timeout=5)
        completed = store.get_job(second)
        if second_task == "differential":
            completed = completed["differential"]
        assert completed["status"] == "completed"
        assert "queued" not in completed.get("message", "")
    finally:
        release.set()
        runner.executor.shutdown(wait=True)


def test_other_live_worker_is_preserved(tmp_path, monkeypatch):
    from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
    from altanalyze3.components.cellHarmony.flask.tasks import JobRunner
    store = JobStore(tmp_path)
    job_id = store.create_job("human", "ref", None, files=[])["job_id"]
    store.update_job(job_id, status="processing", worker_pid=os.getpid() + 1)
    runner = JobRunner(store, tmp_path / "registry.json")
    monkeypatch.setattr(os, "kill", lambda pid, signal: None)
    try:
        assert runner.recover_interrupted_pipeline(job_id)["status"] == "processing"
    finally:
        runner.executor.shutdown(wait=True)


def test_status_endpoint_reports_interrupted_main_job(tmp_path):
    from fastapi.testclient import TestClient
    from altanalyze3.components.cellHarmony.webapp.app import create_app
    app = create_app({"JOB_STORAGE": str(tmp_path / "jobs"), "REFERENCE_REGISTRY": str(tmp_path / "registry.json")})
    store = app.state.job_store
    job_id = store.create_job("human", "ref", None, files=[])["job_id"]
    store.update_job(job_id, status="queued", progress=15, worker_pid=os.getpid())
    try:
        with TestClient(app) as client:
            response = client.get(f"/api/jobs/{job_id}/status")
        assert response.status_code == 200
        payload = response.json()
        assert payload["status"] == "failed"
        assert "interrupted" in payload["message"]
        assert payload["progress"] == 100
    finally:
        app.state.job_runner.executor.shutdown(wait=True)


def test_interrupted_metadata_write_preserves_previous_record(tmp_path, monkeypatch):
    from altanalyze3.components.cellHarmony.flask import job_manager
    store = job_manager.JobStore(tmp_path)
    job_id = store.create_job("human", "ref", None, files=[])["job_id"]
    def partial_write(data, handle, **kwargs):
        handle.write('{"status":')
        raise OSError("simulated interrupted write")
    monkeypatch.setattr(job_manager.json, "dump", partial_write)
    with pytest.raises(OSError, match="interrupted write"):
        store.update_job(job_id, status="processing")
    assert store.get_job(job_id)["status"] == "uploaded"
    assert not list((tmp_path / job_id).glob(".job-*.tmp"))
