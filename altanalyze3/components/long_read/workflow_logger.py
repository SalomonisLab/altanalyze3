#!/usr/bin/env python3
"""Per-step / per-sample logging for the long-read BAM workflow.

Records, for every step (and per sample where applicable):
  - wall-clock time,
  - peak resident memory during the step (process max RSS) and the delta vs. the step start,
  - whether the step used multiprocessing and, if so, the worker count and start method.

Reuses the package's cross-platform RSS helpers (``isoform_automate._get_rss_stats``) so the numbers
match the rest of the pipeline. Writes both a human-readable running log and a structured TSV
(``<logfile>.tsv``) so timing/memory can be tabulated across many samples on a cluster.

Usage:
    wl = WorkflowLogger("workflow_run.log")
    with wl.step("BAM extraction", sample="251H", multiprocessing=True, workers=8):
        ...                       # do the work
    wl.summary()                  # final per-step table
"""

from __future__ import annotations

import multiprocessing
import os
import time
from contextlib import contextmanager

from . import isoform_automate as _ia


def _rss_now_mb():
    """Current and process-peak RSS in MB (cross-platform, via isoform_automate)."""
    try:
        cur, peak = _ia._get_rss_stats()
    except Exception:
        cur, peak = None, None
    return cur, peak


def _mp_start_method():
    try:
        return multiprocessing.get_start_method(allow_none=True) or "default"
    except Exception:
        return "unknown"


class WorkflowLogger:
    def __init__(self, logfile, also_print=True, logs_dir=None):
        """Write logs under a ``logs/`` folder.

        ``logfile`` may be a bare filename or a path. If ``logs_dir`` is given the file is placed
        there; otherwise, when ``logfile`` has no directory component it is placed under ``./logs/``,
        and when it does, a ``logs/`` subfolder is created alongside it. Both the human-readable log
        and the structured ``.tsv`` live in that logs folder.
        """
        logfile = str(logfile)
        parent, name = os.path.split(logfile)
        if logs_dir is None:
            logs_dir = os.path.join(parent or ".", "logs")
        os.makedirs(logs_dir, exist_ok=True)
        self.logfile = os.path.join(logs_dir, name)
        self.tsv = self.logfile + ".tsv"
        self.also_print = also_print
        self.records = []
        self._t0 = time.time()
        with open(self.tsv, "w") as f:
            f.write("step\tsample\tstatus\twall_s\trss_peak_mb\trss_delta_mb\t"
                    "multiprocessing\tworkers\tmp_start_method\tcpu_count\tnote\n")
        self.log(f"=== workflow log started ({time.strftime('%Y-%m-%d %H:%M:%S')}) -> "
                 f"{self.logfile} ===")

    def log(self, message):
        line = f"{time.strftime('%H:%M:%S')}  {message}"
        with open(self.logfile, "a") as f:
            f.write(line + "\n")
        if self.also_print:
            print(line, flush=True)

    @contextmanager
    def step(self, name, sample=None, multiprocessing=False, workers=None, note=""):
        """Time + memory-profile a step. ``multiprocessing``/``workers`` document parallelism."""
        tag = f"{name}" + (f" [{sample}]" if sample else "")
        cur0, peak0 = _rss_now_mb()
        # ``multiprocessing`` may be a bool or, for convenience, the worker count itself.
        if workers is None and isinstance(multiprocessing, int) and not isinstance(multiprocessing, bool):
            workers = multiprocessing
        mp = bool(multiprocessing)
        nworkers = workers if mp else None
        wk = f", workers={nworkers}" if (mp and nworkers) else ""
        self.log(f"START  {tag}  (multiprocessing={mp}{wk}, start={_mp_start_method() if mp else 'n/a'}, "
                 f"rss={cur0:.0f}MB)" if cur0 is not None else f"START  {tag}  (multiprocessing={mp}{wk})")
        t = time.time()
        status = "ok"
        try:
            yield self
        except Exception as exc:
            status = f"FAILED:{type(exc).__name__}"
            raise
        finally:
            wall = time.time() - t
            cur1, peak1 = _rss_now_mb()
            peak = peak1 if peak1 is not None else (peak0 or 0.0)
            delta = (cur1 - cur0) if (cur1 is not None and cur0 is not None) else None
            rec = {
                "step": name, "sample": sample or "", "status": status, "wall_s": wall,
                "rss_peak_mb": peak, "rss_delta_mb": delta,
                "multiprocessing": mp, "workers": nworkers,
                "mp_start_method": _mp_start_method() if mp else "n/a",
                "cpu_count": multiprocessing_count(), "note": note,
            }
            self.records.append(rec)
            with open(self.tsv, "a") as f:
                f.write("\t".join(str(x) for x in [
                    rec["step"], rec["sample"], rec["status"], f"{wall:.1f}",
                    f"{peak:.0f}" if peak else "", f"{delta:.0f}" if delta is not None else "",
                    rec["multiprocessing"], rec["workers"] if rec["workers"] is not None else "",
                    rec["mp_start_method"], rec["cpu_count"], rec["note"],
                ]) + "\n")
            dtxt = f", ΔRSS={delta:+.0f}MB" if delta is not None else ""
            self.log(f"END    {tag}  status={status}  wall={wall:.1f}s  peakRSS={peak:.0f}MB{dtxt}")

    def summary(self):
        self.log("=== STEP SUMMARY (step | sample | status | wall_s | peakRSS_MB | MP | workers) ===")
        for r in self.records:
            self.log(f"  {r['step']:<28} {r['sample']:<22} {r['status']:<14} "
                     f"{r['wall_s']:>7.1f}s  {r['rss_peak_mb']:>8.0f}MB  "
                     f"MP={str(r['multiprocessing']):<5} workers={r['workers']}")
        total = time.time() - self._t0
        self.log(f"=== TOTAL wall={total:.1f}s | structured log: {self.tsv} ===")


def multiprocessing_count():
    try:
        return multiprocessing.cpu_count()
    except Exception:
        return None
