"""Startup thread policy and the actual isolated-worker launch environment."""
from types import SimpleNamespace

import pytest

from altanalyze3.components.cellHarmony.flask.worker_runtime import worker_environment


def test_mac_policy_preserves_parallel_graph_and_analysis_settings():
    original = {"OPENBLAS_NUM_THREADS": "4", "OMP_NUM_THREADS": "4", "NUMBA_NUM_THREADS": "4",
                "UDON_NMF_RUNS": "1", "UDON_NMF_ENGINE": "sklearn", "PYTHONPATH": "/existing"}
    env = worker_environment(original, platform="darwin")
    assert env["OPENBLAS_NUM_THREADS"] == env["VECLIB_MAXIMUM_THREADS"] == "1"
    assert env["MKL_NUM_THREADS"] == env["BLIS_NUM_THREADS"] == "1"
    assert {k: env[k] for k in original if k != "OPENBLAS_NUM_THREADS"} == {
        k: v for k, v in original.items() if k != "OPENBLAS_NUM_THREADS"}
    assert original["OPENBLAS_NUM_THREADS"] == "4"


@pytest.mark.parametrize("platform", ["linux", "win32"])
def test_other_platforms_keep_configured_thread_budgets(platform):
    original = {"OPENBLAS_NUM_THREADS": "4", "OMP_NUM_THREADS": "4"}
    assert worker_environment(original, platform=platform) == original
    assert worker_environment({}, platform=platform) == {}


@pytest.mark.parametrize("task", ["pipeline", "differential"])
@pytest.mark.parametrize("discover", [False, True])
def test_isolated_launch_sets_blas_before_interpreter_imports(tmp_path, monkeypatch, task, discover):
    from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
    from altanalyze3.components.cellHarmony.flask.tasks import JobRunner
    from altanalyze3.components.cellHarmony.flask import tasks
    store = JobStore(tmp_path / "jobs")
    job = store.create_job("human", "test", None, [])["job_id"]
    store.update_job(job, status="completed", differential={"status": "completed"})
    monkeypatch.setattr(tasks.sys, "platform", "darwin")
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "4")
    monkeypatch.setenv("NUMBA_NUM_THREADS", "4")
    launched = []

    def launch(cmd, **kwargs):
        launched.append((cmd, kwargs))
        return SimpleNamespace(pid=12345)

    monkeypatch.setattr(tasks.subprocess, "Popen", launch)
    runner_type = JobRunner
    if discover:
        from altanalyze3.components.cellHarmony.scalable_discover.tasks import DiscoverJobRunner
        runner_type = DiscoverJobRunner
    runner = runner_type(store, tmp_path / "registry.json", isolate_jobs=True)
    monkeypatch.setattr(runner, "_memory_wait_reason", lambda: None)
    monkeypatch.setattr(runner, "_wait_for_worker", lambda process: 0)
    try:
        runner._run_isolated(job, task)
    finally:
        runner.executor.shutdown(wait=True)
    assert len(launched) == 1
    cmd, kwargs = launched[0]
    assert "-m" in cmd and kwargs["start_new_session"] is True
    assert runner.WORKER_MODULE in cmd
    assert kwargs["env"]["OPENBLAS_NUM_THREADS"] == "1"
    assert kwargs["env"]["NUMBA_NUM_THREADS"] == "4"
    assert str(tmp_path / "jobs") in cmd
