"""Read-only deployment check; can be piped to Python inside an existing container.

Run with the same interpreter/environment as the server. This resolves source
locations without importing the analytical pipeline or opening any dataset.
Fresh-process source hashes do not prove what an already-running web process
has imported; compare them with its image/mounts and startup time as well.
"""
import hashlib
from importlib import machinery, metadata
import json
import os
from pathlib import Path
import sys


MODULES = {
    "pipeline": "altanalyze3.components.cellHarmony.flask.pipeline",
    "alignment": "altanalyze3.components.cellHarmony.cellHarmony_lite",
    "mapped_h5ad": "altanalyze3.components.cellHarmony.mapped_h5ad",
    "supervisor": "altanalyze3.components.cellHarmony.flask.tasks",
    "memory_accounting": "altanalyze3.components.cellHarmony.flask.worker_memory",
}


def source_path(name):
    """Resolve Python's normal package paths without executing package __init__."""
    search = None
    parts = name.split(".")
    for index in range(len(parts)):
        spec = machinery.PathFinder.find_spec(parts[index], search)
        if spec is None:
            return None
        search = list(spec.submodule_search_locations) if spec.submodule_search_locations is not None else None
    return Path(spec.origin) if spec.origin and spec.origin not in {"built-in", "frozen"} else None


def source_record(name):
    path = source_path(name)
    if path is None:
        return {"error": "module source not found"}, ""
    try:
        data = path.read_bytes()
    except OSError as exc:
        return {"path": str(path), "error": str(exc)}, ""
    return {"path": str(path), "sha256": hashlib.sha256(data).hexdigest()}, data.decode("utf-8")


def cgroup_snapshot():
    result = {}
    memberships = Path("/proc/self/cgroup")
    if not memberships.exists():
        return result
    for line in memberships.read_text().splitlines():
        _, controllers, relative = line.split(":", 2)
        if controllers == "":
            root = Path("/sys/fs/cgroup")
            names = ("memory.max", "memory.current", "memory.peak", "memory.events",
                     "memory.swap.max", "cpu.max", "memory.stat")
            marker = "memory.current"
        elif "memory" in controllers.split(","):
            root = Path("/sys/fs/cgroup/memory")
            names = ("memory.limit_in_bytes", "memory.usage_in_bytes", "memory.failcnt",
                     "memory.oom_control", "memory.stat")
            marker = "memory.usage_in_bytes"
        else:
            continue
        for base in (root / relative.lstrip("/"), root):
            if not (base / marker).exists():
                continue
            for name in names:
                try:
                    raw = (base / name).read_text().strip()
                except OSError:
                    continue
                if name == "memory.stat":
                    wanted = {"anon", "file", "inactive_file", "active_file", "rss",
                              "total_rss", "total_inactive_file", "total_active_file"}
                    result[name] = {k: int(v) for k, v in (row.split() for row in raw.splitlines())
                                    if k in wanted}
                else:
                    result[name] = raw
            result["path"] = str(base)
            return result
    return result


def report():
    sources, text = {}, {}
    for label, name in MODULES.items():
        sources[label], text[label] = source_record(name)
    versions = {}
    for package in ("anndata", "numpy", "scipy", "h5py", "scanpy"):
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = None
    return {
        "python": sys.executable,
        "sources": sources,
        "source_markers": {
            "automatic_large_h5ad_selection": "needs_disk_backed_import([header])" in text["pipeline"],
            "disk_backed_import": "Workspace(os.path.join(output_dir, '.h5ad_work')).load(h5ad_file)" in text["alignment"],
            "disk_backed_aligned_selection": "workspace.subset(adata_combined, adata_combined.obs_names.get_indexer(match_df.CellBarcode))" in text["alignment"],
            "linux_scratch_cache_release": "os.posix_fadvise" in text["mapped_h5ad"],
            "sigkill_diagnostic": "returncode == -9" in text["supervisor"],
        },
        "packages": versions,
        "settings": {key: value for key, value in os.environ.items()
                     if key in {"CELLHARMONY_JOB_WORKERS", "CELLHARMONY_ISOLATE_JOBS",
                                "CELLHARMONY_WORKER_MEMORY_LIMIT_GIB", "CELLHARMONY_TOTAL_MEMORY_LIMIT_GIB",
                                "CELLHARMONY_CACHE_MAX_GIB"}},
        "cgroup": cgroup_snapshot(),
        "scope": "Fresh interpreter source resolution; marker checks indicate source text, not proof that a job used a branch. Memory events are cumulative, not attributed to a particular job.",
    }


if __name__ == "__main__":
    print(json.dumps(report(), indent=2, sort_keys=True))
