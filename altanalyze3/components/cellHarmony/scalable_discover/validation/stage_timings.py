"""Stage durations of one scALABLE-discover job, read from its logs/pipeline.log.

  python stage_timings.py <job dir> [<job dir> ...]

Each stage runs from one log marker to the next, so the stages sum to the job's duration.
"""
import re, sys
from datetime import datetime
from pathlib import Path

MARKERS = [  # (stage that ENDS at this marker, regex of the marker line)
    (None, r"Job accepted by worker\."),
    ("scALABLE load, ambient RNA, QC, ICGS3 input write", r"Running ICGS3 unsupervised clustering\."),
    ("ICGS3 load, normalize, gene filter", r"downsampling summary:"),
    ("ICGS3 feature selection, NMF rank estimate", r"running UDON NMF rank="),
    ("ICGS3 NMF", r"pre-SVM NMF clusters"),
    ("ICGS3 MarkerFinder, SVM, MarkerFinder", r"marker-robust SVM clusters"),
    ("ICGS3 GO-Elite BioMarkers", r"running final UMAP"),
    ("ICGS3 UMAP fit and plots", r"UMAP output completed"),
    ("ICGS3 marker heatmap (static render of every cell)", r"MarkerFinder heatmap completed"),
    ("ICGS3 h5ad write (gzip)", r"h5ad output completed"),
    ("combined h5ad, fastComm, downloads, bundle", r"\] Job completed\."),
]


def stage_times(job_dir):
    lines = (Path(job_dir) / "logs" / "pipeline.log").read_text().splitlines()
    stamps = []
    for _, pattern in MARKERS:
        hit = next((l for l in lines if re.search(pattern, l)), None)
        if hit is None:
            raise ValueError(f"marker {pattern!r} missing from {job_dir}")
        stamps.append(datetime.fromisoformat(hit[1:hit.index("]")]))
    render = re.search(r"render_heatmap_pdf=([0-9.]+)s", "\n".join(lines))
    rows = [(MARKERS[i][0], (stamps[i] - stamps[i - 1]).total_seconds()) for i in range(1, len(MARKERS))]
    return rows, (stamps[-1] - stamps[0]).total_seconds(), float(render.group(1)) if render else None


for job in sys.argv[1:]:
    rows, total, render = stage_times(job)
    print(f"== {Path(job).name}  total {total:.1f} s; static heatmap render alone {render} s")
    for name, seconds in rows:
        print(f"{seconds:8.1f} s  {100 * seconds / total:5.1f}%  {name}")
    print(f"{sum(s for _, s in rows):8.1f} s  sum of stages")
