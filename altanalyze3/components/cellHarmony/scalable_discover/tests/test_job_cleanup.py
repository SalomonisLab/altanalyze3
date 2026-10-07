"""Automatic retention uses the shared web policy and runs without new uploads."""
import json
import threading
from datetime import datetime, timedelta, timezone
from unittest.mock import Mock

from fastapi.testclient import TestClient

from altanalyze3.components.cellHarmony.scalable_discover import app as discover


def seed_job(root, job_id, *, status="completed", age_hours=9, differential="idle"):
    directory = root / job_id
    directory.mkdir()
    for name in ("uploads", "outputs", "logs"):
        (directory / name).mkdir()
        (directory / name / "payload").write_text("fixture")
    timestamp = (datetime.now(timezone.utc) - timedelta(hours=age_hours)).isoformat()
    (directory / "job.json").write_text(json.dumps({
        "job_id": job_id, "status": status, "updated_at": timestamp,
        "created_at": timestamp, "differential": {"status": differential},
    }))
    return directory


def make_app(root):
    return discover.create_discover_app({"JOB_STORAGE": str(root), "ISOLATE_JOBS": False})


def test_startup_removes_only_expired_terminal_jobs(tmp_path):
    expired = [seed_job(tmp_path, status, status=status)
               for status in ("completed", "failed", "cancelled", "canceled")]
    retained = [seed_job(tmp_path, "recent", age_hours=7)]
    retained += [seed_job(tmp_path, status, status=status, age_hours=100)
                 for status in ("uploaded", "queued", "processing", "running")]
    retained += [seed_job(tmp_path, "differential_" + status, differential=status)
                 for status in ("queued", "processing")]
    malformed = seed_job(tmp_path, "malformed")
    (malformed / "job.json").write_text("{invalid")
    retained.append(malformed)
    missing = tmp_path / "no_metadata"
    missing.mkdir()
    retained.append(missing)
    app = make_app(tmp_path)
    with TestClient(app) as client:
        assert client.get("/api/meta/species").status_code == 200
        assert all(not path.exists() for path in expired)
        assert all(path.exists() for path in retained)
        task = app.state.job_cleanup_task
        assert not task.done()
    assert task.done() and task.cancelled()


def test_periodic_cleanup_removes_job_without_an_upload(tmp_path, monkeypatch):
    monkeypatch.setattr(discover, "JOB_CLEANUP_INTERVAL_SECONDS", 0.01)
    app = make_app(tmp_path)
    removed = threading.Event()
    original = app.state.job_store.purge_old_jobs

    def purge():
        count = original()
        if count:
            removed.set()
        return count

    monkeypatch.setattr(app.state.job_store, "purge_old_jobs", purge)
    with TestClient(app) as client:
        expired = seed_job(tmp_path, "expired_after_startup")
        assert removed.wait(timeout=5)
        assert not expired.exists()
        assert client.get("/api/meta/species").status_code == 200


def test_cleanup_failure_does_not_stop_server_or_retry(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(discover, "JOB_CLEANUP_INTERVAL_SECONDS", 0.01)
    app = make_app(tmp_path)
    retried = threading.Event()
    original = app.state.job_store.purge_old_jobs

    def run():
        if purge.call_count == 1:
            raise OSError("temporary filesystem error")
        result = original()
        retried.set()
        return result

    purge = Mock(side_effect=run)
    monkeypatch.setattr(app.state.job_store, "purge_old_jobs", purge)
    expired = seed_job(tmp_path, "expired")
    with TestClient(app) as client:
        assert retried.wait(timeout=5)
        assert not expired.exists()
        assert client.get("/api/meta/species").status_code == 200
    assert "job cleanup failed" in caplog.text
