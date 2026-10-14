"""Short HTTP polls for expensive plots; one builder and bounded disk retention."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory
import threading
import time
import logging

from fastapi import HTTPException
from fastapi.responses import JSONResponse, Response


class PlotBuildQueue:
    def __init__(self, max_entries=8, ttl_seconds=300, clock=time.monotonic,
                 pending_message="Preparing CombPlot…", cache_errors=True):
        self.executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="plot-build")
        self.directory = TemporaryDirectory(prefix="scalable-plots-")
        self.entries = {}
        self.progress = {}
        self.lock = threading.Lock()
        self.max_entries, self.ttl, self.clock = max_entries, ttl_seconds, clock
        self.serial = 0
        self.pending_message = pending_message
        self.cache_errors = cache_errors

    def _discard(self, key):
        _, _, path = self.entries.pop(key)
        path.unlink(missing_ok=True)
        self.progress.pop(key, None)

    def report_progress(self, key, *, stage, detail):
        with self.lock:
            if key in self.progress:
                self.progress[key].update(stage=stage, detail=detail)

    def request(self, key, build):
        with self.lock:
            now = self.clock()
            for name, (future, ready_at, _) in list(self.entries.items()):
                if future.done() and ready_at is not None and now - ready_at >= self.ttl:
                    self._discard(name)
            entry = self.entries.get(key)
            if entry is None:
                if len(self.entries) >= self.max_entries:
                    finished = next((k for k, (f, _, _) in self.entries.items() if f.done()), None)
                    if finished is not None:
                        self._discard(finished)
                    else:
                        return self._pending("Waiting for a background worker.")
                self.serial += 1
                path = Path(self.directory.name) / f"{self.serial}.json"
                self.progress[key] = {'submitted_at': now}

                def run():
                    try:
                        response = build()
                    except HTTPException as exc:
                        response = JSONResponse({"detail": exc.detail}, status_code=exc.status_code)
                    except Exception as exc:
                        logging.exception("Plot build failed")
                        response = JSONResponse({"detail": f"{type(exc).__name__}: {exc}"}, status_code=500)
                    path.write_bytes(response.body)
                    with self.lock:
                        future, _, _ = self.entries[key]
                        self.entries[key] = (future, self.clock(), path)
                    return response.status_code

                future = self.executor.submit(run)
                self.entries[key] = (future, None, path)
                entry = self.entries[key]
            future, ready_at, path = entry
            if not future.done():
                info = self.progress.get(key, {})
                waiting = not future.running()
                return self._pending("Waiting for another background request to finish." if waiting else info.get('detail', self.pending_message),
                                     stage='queued' if waiting else info.get('stage', 'starting'),
                                     elapsed_seconds=round(now - info.get('submitted_at', now), 1))
            if ready_at is None:
                self.entries[key] = (future, now, path)
            # Completed results live on disk, including plots too large for a
            # RAM cache. Failures propagate through the normal API error handler.
            status = future.result()
            body = path.read_bytes()
            if status >= 400 and not self.cache_errors:
                self._discard(key)
            return Response(body, status_code=status, media_type="application/json")

    @staticmethod
    def _pending(message, **progress):
        return JSONResponse({"status": "preparing", "detail": message, **progress}, status_code=202,
                            headers={"Retry-After": "1", "Cache-Control": "no-store"})

    def close(self):
        self.executor.shutdown(wait=True, cancel_futures=True)
        self.directory.cleanup()
