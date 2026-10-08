"""The scALABLE job runner, pointed at the scALABLE-discover pipeline and worker."""
from __future__ import annotations

import contextlib
import os
import re
import time
from datetime import datetime, timezone

from altanalyze3.components.cellHarmony.flask.tasks import JobRunner, _JobLogStream

from .pipeline import run_discover_pipeline

# Log line fragment -> (percent complete, stage shown on the Run tab). QC lines come from
# cellHarmony_lite; ICGS3 lines are the ones each major ICGS3 step writes (ICGS.py `_log`).
# Progress only moves forward, so a fragment that repeats (PageRank logs one line per chunk)
# never sends the bar back. ICGS3's default path runs no PCA: downsampling builds its k-NN
# graph on dispersion-selected variable genes.
STAGES = (
    ("Running ambient RNA correction", 23, "QC: ambient RNA correction"),
    ("Cells remaining after min_genes", 27, "QC: filtering cells"),
    ("Normalization steps", 31, "QC: normalizing"),
    ("[ICGS3] loading ", 36, "ICGS3 step 1 of 10: loading the QC-retained counts"),
    ("[ICGS3] normalization mode", 38, "ICGS3 step 2 of 10: normalizing (log1p CP10K)"),
    ("RNA unsupervised gene filter before", 40, "ICGS3 step 3 of 10: protein-coding gene filter"),
    ("ICGS2 community sampling: running Louvain", 42,
     "ICGS3 step 4 of 10: downsampling, variable genes and Louvain communities"),
    ("ICGS2 PageRank sampling", 44, "ICGS3 step 4 of 10: downsampling, variable genes and PageRank"),
    ("downsampling summary", 46, "ICGS3 step 5 of 10: variable guide-gene selection and sNMF rank"),
    ("running UDON NMF", 50, "ICGS3 step 6 of 10: sNMF clustering"),
    ("pre-SVM NMF clusters", 58, "ICGS3 step 7 of 10: MarkerFinder on the sNMF clusters"),
    ("SVM reclassification target", 61, "ICGS3 step 8 of 10: SVM reclassification of every cell"),
    ("post-SVM MarkerFinder gene pool", 63, "ICGS3 step 8 of 10: MarkerFinder on the SVM cell states"),
    ("marker-robust SVM clusters", 66, "ICGS3 step 9 of 10: GO-Elite BioMarkers cell-state predictions"),
    ("running final UMAP", 68, "ICGS3 step 10 of 10: UMAP"),
    ("UMAP using graph/PCA fallback", 68, "ICGS3 step 10 of 10: UMAP"),
    ("UMAP fit mode=", 69, "ICGS3 step 10 of 10: fitting UMAP"),
    ("UMAP mapping ", 72, "ICGS3 step 10 of 10: mapping remaining cells into UMAP"),
    ("UMAP fit completed in", 75, "ICGS3 step 10 of 10: writing UMAP"),
    ("UMAP output completed", 76, "ICGS3: MarkerFinder marker set and networks"),
    ("MarkerFinder heatmap completed", 79, "ICGS3: writing results"),
    ("ICGS3 complete in", 80, "ICGS3 complete"),
)


class _StageLogStream(_JobLogStream):
    """The job log, plus a progress and stage update whenever a stage's log line appears."""

    def __init__(self, store, job_id):
        super().__init__(store, job_id)
        self._progress = 0
        self._stage_message = ""

    def _stage(self, line: str) -> None:
        mapped = re.search(r"UMAP transformed [\d,]+ remaining cells \(([\d,]+)/([\d,]+) completed;", line)
        if mapped:
            completed, total = (int(value.replace(",", "")) for value in mapped.groups())
            if 0 < completed <= total:
                progress = 72 + int(3 * completed / total)
                message = f"ICGS3 step 10 of 10: UMAP mapped {completed:,} of {total:,} remaining cells"
                if progress >= self._progress and message != self._stage_message:
                    self._progress = progress
                    self._stage_message = message
                    self.store.update_job(self.job_id, progress=progress, message=message)
            return
        for fragment, progress, message in STAGES:
            if fragment in line:
                if progress > self._progress or (progress == self._progress and message != self._stage_message):
                    self._progress = progress
                    self._stage_message = message
                    self.store.update_job(self.job_id, progress=progress, message=message)
                return

    def write(self, data: str) -> int:
        if not data:
            return 0
        self._buffer += data
        while "\n" in self._buffer:
            line, self._buffer = self._buffer.split("\n", 1)
            line = line.rstrip()
            if line:
                self.store.append_log(self.job_id, line)
                self._stage(line)
        return len(data)

    def flush(self) -> None:
        # ICGS3's Tee flushes after every write (ICGS.py Tee.write), so its lines arrive here,
        # not through the newline split in write(). Watching only write() missed every one.
        if self._buffer.strip():
            line = self._buffer.strip()
            self.store.append_log(self.job_id, line)
            self._stage(line)
        self._buffer = ""


class DiscoverJobRunner(JobRunner):
    """Same queue, memory admission and isolated workers as scALABLE-web; ICGS3 pipeline."""

    WORKER_MODULE = "altanalyze3.components.cellHarmony.scalable_discover.worker"

    def submit_differential(self, job_id: str) -> None:
        raise RuntimeError("scALABLE-discover runs no differential analysis.")

    def _run_pipeline(self, job_id: str) -> None:
        started = time.perf_counter()
        try:
            self.store.update_job(job_id, status="processing", message="Preparing inputs…", progress=10,
                                  worker_pid=os.getpid(), analysis_started_at=datetime.now(timezone.utc).isoformat(),
                                  analysis_completed_at=None, analysis_duration_seconds=None)
            self.store.append_log(job_id, "Job accepted by worker.")
            log_stream = _StageLogStream(self.store, job_id)
            with contextlib.redirect_stdout(log_stream), contextlib.redirect_stderr(log_stream):
                run_discover_pipeline(job_id, self.store, h5ad_compression=self.h5ad_compression)
            log_stream.flush()
            self.store.update_job(job_id, status="completed", message="Job finished successfully.", progress=100,
                                  worker_pid=None, analysis_completed_at=datetime.now(timezone.utc).isoformat(),
                                  analysis_duration_seconds=time.perf_counter() - started)
            self.store.append_log(job_id, "Job completed.")
        except Exception as exc:  # recorded in the job: status, message and traceback in the log
            summary = self._log_failure(job_id, "Job failed", exc)
            self.store.update_job(job_id, status="failed", message=summary, progress=100, worker_pid=None)
