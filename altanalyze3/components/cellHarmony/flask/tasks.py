from __future__ import annotations

import contextlib
import io
import threading
import time
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path
from typing import Dict

from .job_manager import JobStore
from .pipeline import run_cellharmony_pipeline


class _JobLogStream(io.TextIOBase):
    def __init__(self, store: JobStore, job_id: str):
        self.store = store
        self.job_id = job_id
        self._buffer = ""

    def write(self, data: str) -> int:
        if not data:
            return 0
        self._buffer += data
        while "\n" in self._buffer:
            line, self._buffer = self._buffer.split("\n", 1)
            line = line.rstrip()
            if line:
                self.store.append_log(self.job_id, line)
        return len(data)

    def flush(self) -> None:
        if self._buffer.strip():
            self.store.append_log(self.job_id, self._buffer.strip())
        self._buffer = ""


class JobRunner:
    """Fire-and-forget background executor for processing jobs."""

    def __init__(self, store: JobStore, registry_path: Path, max_workers: int = 2):
        self.store = store
        self.registry_path = Path(registry_path)
        self.executor = ThreadPoolExecutor(max_workers=max_workers)
        self._futures: Dict[str, Future] = {}
        self._lock = threading.Lock()

    def submit(self, job_id: str) -> None:
        with self._lock:
            existing = self._futures.get(job_id)
            if existing and not existing.done():
                return
            future = self.executor.submit(self._run_pipeline, job_id)
            self._futures[job_id] = future

    def _run_pipeline(self, job_id: str) -> None:
        try:
            self.store.update_job(job_id, status="processing", message="Preparing inputs…", progress=10)
            self.store.append_log(job_id, "Job accepted by worker.")
            time.sleep(0.1)
            log_stream = _JobLogStream(self.store, job_id)
            with contextlib.redirect_stdout(log_stream), contextlib.redirect_stderr(log_stream):
                run_cellharmony_pipeline(job_id, self.store, self.registry_path)
            log_stream.flush()
            self.store.update_job(job_id, status="completed", message="Job finished successfully.", progress=100)
            self.store.append_log(job_id, "Job completed.")
        except Exception as exc:  # pragma: no cover - defensive
            self.store.append_log(job_id, f"Job failed: {exc}")
            self.store.update_job(job_id, status="failed", message=str(exc), progress=100)
