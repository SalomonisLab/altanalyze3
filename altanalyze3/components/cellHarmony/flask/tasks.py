from __future__ import annotations

import contextlib
import io
import threading
import time
import traceback
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path
from typing import Dict

from .job_manager import JobStore
from .pipeline import run_cellharmony_differential, run_cellharmony_pipeline


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

    def __init__(
        self,
        store: JobStore,
        registry_path: Path,
        max_workers: int = 2,
        *,
        export_approx_pdfs: bool = True,
        h5ad_compression: str = "lzf",
    ):
        self.store = store
        self.registry_path = Path(registry_path)
        self.executor = ThreadPoolExecutor(max_workers=max_workers)
        self._futures: Dict[str, Future] = {}
        self._lock = threading.Lock()
        self.export_approx_pdfs = bool(export_approx_pdfs)
        self.h5ad_compression = str(h5ad_compression or "lzf")

    def submit(self, job_id: str) -> None:
        self.submit_pipeline(job_id)

    def submit_pipeline(self, job_id: str) -> None:
        self._submit(job_id, "pipeline", self._run_pipeline)

    def submit_differential(self, job_id: str) -> None:
        self._submit(job_id, "differential", self._run_differential)

    def _submit(self, job_id: str, task_name: str, target) -> None:
        key = f"{task_name}:{job_id}"
        with self._lock:
            existing = self._futures.get(key)
            if existing and not existing.done():
                return
            future = self.executor.submit(target, job_id)
            self._futures[key] = future

    def _update_differential(self, job_id: str, **changes) -> None:
        meta = self.store.get_job(job_id)
        differential = dict(meta.get("differential") or {})
        differential.update(changes)
        self.store.update_job(job_id, differential=differential)

    @staticmethod
    def _exception_summary(exc: BaseException) -> str:
        name = type(exc).__name__
        detail = str(exc).strip()
        return f"{name}: {detail}" if detail else name

    def _log_failure(self, job_id: str, prefix: str, exc: BaseException) -> str:
        summary = self._exception_summary(exc)
        self.store.append_log(job_id, f"{prefix}: {summary}")
        for line in traceback.format_exc().strip().splitlines():
            if line:
                self.store.append_log(job_id, line)
        return summary

    def _run_pipeline(self, job_id: str) -> None:
        try:
            self.store.update_job(job_id, status="processing", message="Preparing inputs…", progress=10)
            self.store.append_log(job_id, "Job accepted by worker.")
            time.sleep(0.1)
            log_stream = _JobLogStream(self.store, job_id)
            with contextlib.redirect_stdout(log_stream), contextlib.redirect_stderr(log_stream):
                run_cellharmony_pipeline(
                    job_id,
                    self.store,
                    self.registry_path,
                    export_approx_pdfs=self.export_approx_pdfs,
                    h5ad_compression=self.h5ad_compression,
                )
            log_stream.flush()
            self.store.update_job(job_id, status="completed", message="Job finished successfully.", progress=100)
            self.store.append_log(job_id, "Job completed.")
        except Exception as exc:  # pragma: no cover - defensive
            summary = self._log_failure(job_id, "Job failed", exc)
            self.store.update_job(job_id, status="failed", message=summary, progress=100)

    def _run_differential(self, job_id: str) -> None:
        try:
            self._update_differential(job_id, status="processing", message="Preparing differential analysis.", progress=10)
            self.store.append_log(job_id, "Differential analysis accepted by worker.")
            time.sleep(0.1)
            log_stream = _JobLogStream(self.store, job_id)
            with contextlib.redirect_stdout(log_stream), contextlib.redirect_stderr(log_stream):
                run_cellharmony_differential(job_id, self.store)
            log_stream.flush()
            self._update_differential(job_id, status="completed", message="Differential analysis finished.", progress=100)
            self.store.append_log(job_id, "Differential analysis completed.")
        except Exception as exc:  # pragma: no cover - defensive
            summary = self._log_failure(job_id, "Differential analysis failed", exc)
            self._update_differential(job_id, status="failed", message=summary, progress=100)
