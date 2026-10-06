from __future__ import annotations

import contextlib
import io
import os
import subprocess
import sys
import threading
import time
import traceback
from concurrent.futures import Future, ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

from .job_manager import JobStore
from .worker_memory import process_memory, tree_rss, container_memory
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

    # The module each isolated analysis runs in; a subclass with its own pipeline names its own.
    WORKER_MODULE = "altanalyze3.components.cellHarmony.flask.worker"

    def __init__(
        self,
        store: JobStore,
        registry_path: Path,
        max_workers: int = 1,
        *,
        export_approx_pdfs: bool = True,
        h5ad_compression: str = "lzf",
        isolate_jobs: bool = False,
        worker_memory_limit_gib: float = 15,
        total_memory_limit_gib: float = 27,
    ):
        self.store = store
        self.registry_path = Path(registry_path)
        self.executor = ThreadPoolExecutor(max_workers=max_workers)
        self.max_workers = max_workers
        self._futures: Dict[str, Future] = {}
        self._lock = threading.Lock()
        self.export_approx_pdfs = bool(export_approx_pdfs)
        self.h5ad_compression = str(h5ad_compression or "lzf")
        self.isolate_jobs = bool(isolate_jobs)
        self.worker_memory_limit = int(worker_memory_limit_gib * 1024**3)
        self.total_memory_limit = int(total_memory_limit_gib * 1024**3)
        if self.worker_memory_limit <= 0 or self.total_memory_limit <= 0:
            raise ValueError("Memory admission limits must be positive.")
        self._admission = threading.Condition()
        self._active_workers = {}

    def _memory_wait_reason(self):
        try:
            processes = process_memory()
            if os.getpid() not in processes or any(pid not in processes for pid in self._active_workers):
                return "Worker memory usage is temporarily unavailable."
            if any(tree_rss(processes, pid) >= self.worker_memory_limit for pid in self._active_workers):
                return "A running analysis has reached the per-worker memory threshold."
            container = container_memory()
            used, limit = container if container else (tree_rss(processes, os.getpid()), self.total_memory_limit)
            ceiling = min(self.total_memory_limit, int(limit * 0.9)) if container else limit
            if used >= ceiling:
                return "The server is currently using its available analysis memory."
        except (OSError, ValueError, subprocess.SubprocessError):
            return "Worker memory usage is temporarily unavailable."
        return None

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
            busy = sum(not future.done() for future in self._futures.values()) >= self.max_workers
            queue_message = (
                "Your analysis is queued. Other analyses are currently running or waiting; "
                "yours will start automatically when a worker becomes available."
                if busy else
                "Your analysis is queued and will start automatically when a worker becomes available."
            )
            if task_name == "pipeline":
                # Persist before starting: a fast worker must not have its completed
                # or processing state overwritten by the request handler.
                self.store.update_job(
                    job_id, status="queued", progress=15,
                    message=queue_message, worker_pid=os.getpid(),
                )
                self.store.append_log(job_id, "Job queued by user request.")
            else:
                differential = dict(self.store.get_job(job_id).get("differential") or {})
                differential.update(status="queued", message=queue_message, worker_pid=os.getpid())
                self.store.update_job(job_id, differential=differential)
            future = (self.executor.submit(self._run_isolated, job_id, task_name)
                      if self.isolate_jobs else self.executor.submit(target, job_id))
            self._futures[key] = future

    def _run_isolated(self, job_id: str, task_name: str) -> None:
        """One fresh interpreter per analysis; only the supervisor stays in the web app."""
        log_path = self.store.logs_dir(job_id) / f"{task_name}_worker.log"
        try:
            env = dict(os.environ)
            root = str(Path(__file__).resolve().parents[4])
            env["PYTHONPATH"] = root + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
            cmd = [sys.executable, "-m", self.WORKER_MODULE,
                   str(self.store.root.resolve()), str(self.registry_path.resolve()), job_id, task_name,
                   self.h5ad_compression, str(int(self.export_approx_pdfs))]
            with log_path.open("w") as log:
                # Check and launch under one lock so two admissions cannot race.
                with self._admission:
                    last_reason = None
                    while True:
                        reason = self._memory_wait_reason()
                        if reason is None:
                            break
                        if reason != last_reason:
                            message = f"Your analysis is queued. {reason} It will start automatically when memory becomes available."
                            if task_name == "pipeline":
                                self.store.update_job(job_id, message=message)
                            else:
                                self._update_differential(job_id, message=message)
                            last_reason = reason
                        self._admission.wait(timeout=1)
                    process = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
                    self._active_workers[process.pid] = process
                try:
                    returncode = process.wait()
                finally:
                    with self._admission:
                        self._active_workers.pop(process.pid, None)
                        self._admission.notify_all()
            meta = self.store.get_job(job_id)
            state = meta if task_name == "pipeline" else meta.get("differential", {})
            if returncode or state.get("status") not in {"completed", "failed"}:
                raise RuntimeError(f"Analysis worker exited with code {returncode}; see {log_path.name}.")
        except Exception as exc:
            summary = self._log_failure(job_id, "Analysis worker failed", exc)
            if task_name == "pipeline":
                self.store.update_job(job_id, status="failed", progress=100, message=summary, worker_pid=None)
            else:
                self._update_differential(job_id, status="failed", message=summary, worker_pid=None)

    def recover_interrupted_pipeline(self, job_id: str) -> Dict:
        """Turn an orphaned pipeline into a terminal failure when status is polled."""
        with self._lock:
            meta = self.store.get_job(job_id)
            if meta.get("status") not in {"queued", "processing"}:
                return meta
            future = self._futures.get(f"pipeline:{job_id}")
            if future is not None and not future.done():
                return meta
            owner = meta.get("worker_pid")
            if owner and int(owner) != os.getpid():
                try:
                    os.kill(int(owner), 0)
                except ProcessLookupError:
                    pass
                except PermissionError:
                    return meta
                else:
                    return meta
            # A Docker restart can reuse the same PID. A fresh runner has no
            # matching Future, even when the persisted PID happens to match.
            message = "Job was interrupted by a worker restart or shutdown. Run it again to finish."
            meta = self.store.update_job(
                job_id, status="failed", message=message, progress=100,
                worker_pid=None,
            )
            self.store.append_log(job_id, message)
            return meta

    def recover_interrupted_differential(self, job_id: str) -> Dict:
        """Reconcile persisted progress with the worker that actually owns the task."""
        with self._lock:
            meta = self.store.get_job(job_id)
            differential = dict(meta.get("differential") or {})
            if differential.get("status") not in {"queued", "processing"}:
                return meta
            future = self._futures.get(f"differential:{job_id}")
            if future is not None and not future.done():
                return meta
            owner = differential.get("worker_pid")
            if owner and int(owner) != os.getpid():
                try:
                    os.kill(int(owner), 0)
                except ProcessLookupError:
                    pass
                except PermissionError:
                    return meta  # Another live process owns this task.
                else:
                    return meta
            message = "Differential analysis was interrupted. Run it again to finish the saved comparison."
            differential.update(status="failed", message=message, worker_pid=None)
            meta = self.store.update_job(job_id, differential=differential)
            self.store.append_log(job_id, message)
            return meta

    def _update_differential(self, job_id: str, **changes) -> None:
        meta = self.store.get_job(job_id)
        differential = dict(meta.get("differential") or {})
        differential.update(changes)
        updates = {"differential": differential}
        if differential.get("status") == "completed" and differential.get("run_id"):
            history = dict(meta.get("differential_history") or {})
            history[differential["run_id"]] = differential
            updates["differential_history"] = history
        self.store.update_job(job_id, **updates)

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
        started = time.perf_counter()
        try:
            self.store.update_job(job_id, status="processing", message="Preparing inputs…", progress=10,
                                  worker_pid=os.getpid(), analysis_started_at=datetime.now(timezone.utc).isoformat(),
                                  analysis_completed_at=None, analysis_duration_seconds=None)
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
            self.store.update_job(job_id, status="completed", message="Job finished successfully.", progress=100, worker_pid=None,
                                  analysis_completed_at=datetime.now(timezone.utc).isoformat(),
                                  analysis_duration_seconds=time.perf_counter() - started)
            self.store.append_log(job_id, "Job completed.")
        except Exception as exc:  # pragma: no cover - defensive
            summary = self._log_failure(job_id, "Job failed", exc)
            self.store.update_job(job_id, status="failed", message=summary, progress=100, worker_pid=None)

    def _run_differential(self, job_id: str) -> None:
        try:
            self._update_differential(job_id, status="processing", message="Preparing differential analysis.", progress=10, worker_pid=os.getpid())
            self.store.append_log(job_id, "Differential analysis accepted by worker.")
            time.sleep(0.1)
            log_stream = _JobLogStream(self.store, job_id)
            with contextlib.redirect_stdout(log_stream), contextlib.redirect_stderr(log_stream):
                run_cellharmony_differential(job_id, self.store)
            log_stream.flush()
            self._update_differential(job_id, status="completed", message="Differential analysis finished.", progress=100, worker_pid=None)
            self.store.append_log(job_id, "Differential analysis completed.")
        except Exception as exc:  # pragma: no cover - defensive
            summary = self._log_failure(job_id, "Differential analysis failed", exc)
            self._update_differential(job_id, status="failed", message=summary, progress=100, worker_pid=None)
