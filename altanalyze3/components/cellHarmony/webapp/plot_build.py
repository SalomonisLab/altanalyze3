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
    def __init__(self, max_entries=8, ttl_seconds=300, clock=time.monotonic):
        self.executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="plot-build")
        self.directory = TemporaryDirectory(prefix="scalable-plots-")
        self.entries = {}
        self.lock = threading.Lock()
        self.max_entries, self.ttl, self.clock = max_entries, ttl_seconds, clock
        self.serial = 0

    def _discard(self, key):
        _, _, path = self.entries.pop(key)
        path.unlink(missing_ok=True)

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
                        return self._pending("Waiting for a plot worker.")
                self.serial += 1
                path = Path(self.directory.name) / f"{self.serial}.json"

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
                return self._pending("Preparing CombPlot…")
            if ready_at is None:
                self.entries[key] = (future, now, path)
            # Completed results live on disk, including plots too large for a
            # RAM cache. Failures propagate through the normal API error handler.
            status = future.result()
            return Response(path.read_bytes(), status_code=status, media_type="application/json")

    @staticmethod
    def _pending(message):
        return JSONResponse({"status": "preparing", "detail": message}, status_code=202,
                            headers={"Retry-After": "1", "Cache-Control": "no-store"})

    def close(self):
        self.executor.shutdown(wait=True, cancel_futures=True)
        self.directory.cleanup()
