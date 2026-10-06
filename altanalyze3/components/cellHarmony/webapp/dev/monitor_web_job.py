"""Sample a real web server and its analysis descendants; never execute analysis.

Submit uploads/QC/run through the normal web API or browser. RSS measurements on
macOS are not cgroup measurements and cannot establish a Docker memory guarantee.
The raw samples include independent server/worker RSS and CPU time deltas.
"""
import argparse
import json
import runpy
import time
from pathlib import Path

import psutil

# Importing flask's package also loads its analysis pipeline. Read this small
# production accounting module directly so the observer stays lightweight.
container_memory = runpy.run_path(str(Path(__file__).resolve().parents[2]
                                    / 'flask/worker_memory.py'))['container_memory']


def cgroup_sample():
    """Keep hard-limit events and page-cache accounting alongside process RSS."""
    usage = container_memory()
    if usage is None:
        return None
    result = dict(working_set_bytes=usage[0], limit_bytes=usage[1])
    base = Path('/sys/fs/cgroup')
    for name in ('memory.current', 'memory.peak'):
        try:
            result[name] = int((base / name).read_text())
        except (OSError, ValueError):
            pass
    for name in ('memory.events', 'memory.stat'):
        try:
            result[name] = {key: int(value) for key, value in
                            (line.split() for line in (base / name).read_text().splitlines())}
        except (OSError, ValueError):
            pass
    return result


TRANSITIONS = (
    ('Running cellHarmony_lite pipeline.', 'import'),
    ('[INFO] Running ambient RNA correction', 'ambient_correction'),
    ('[INFO] Ambient RNA correction complete.', 'qc_normalization'),
    ('Aligning cells to reference...', 'centroid_alignment'),
    ('Identifying cell-state marker genes.', 'rna_markerfinder'),
    ('Running approximate UMAP placement.', 'approximate_umap'),
    ('Running rna2lipid lipid imputation.', 'lipid_imputation'),
    ('rna2lipid lipid imputation complete.', 'post_imputation'),
    ('Running fastComm receptor-ligand', 'communication'),
    ('[bundle] building', 'serving_bundle'),
    ('cellHarmony pipeline finished.', 'finished'),
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--server-pid', type=int, required=True)
    parser.add_argument('--job-dir', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--interval', type=float, default=0.5)
    args = parser.parse_args()
    args.report.parent.mkdir(parents=True, exist_ok=True)
    server = psutil.Process(args.server_pid)
    started = time.monotonic()
    previous_cpu = {}
    phases = {}
    stage, offset, pending = 'awaiting_run', 0, ''
    terminal_since = None
    previous_elapsed = 0
    peak_cgroup_working_set = 0
    raw_path = args.report.with_suffix('.jsonl')
    with raw_path.open('w') as raw:
        while server.is_running():
            elapsed = time.monotonic() - started
            path = args.job_dir / 'logs/pipeline.log'
            if path.exists():
                with path.open() as stream:
                    stream.seek(offset)
                    incoming = stream.read()
                    offset = stream.tell()
                lines = (pending + incoming).split('\n')
                pending = lines.pop()
                for line in lines:
                    for marker, next_stage in TRANSITIONS:
                        if marker in line:
                            stage = next_stage
                            break
            detail, rss, worker_rss, cpu_seconds, threads = [], 0, 0, 0., 0
            for process in [server, *server.children(recursive=True)]:
                try:
                    with process.oneshot():
                        resident = process.memory_info().rss
                        cpu = process.cpu_times()
                        used = cpu.user + cpu.system
                        count = process.num_threads()
                        command = process.cmdline()
                    delta = max(0, used - previous_cpu.get(process.pid, used))
                    previous_cpu[process.pid] = used
                    cpu_seconds += delta
                    rss += resident
                    threads += count
                    if process.pid != server.pid:
                        worker_rss += resident
                    detail.append(dict(pid=process.pid, rss_bytes=resident,
                                       cpu_seconds=used, threads=count, command=command))
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    continue
            duration = elapsed - previous_elapsed
            cores = cpu_seconds / duration if duration > 0 else 0
            previous_elapsed = elapsed
            meta = json.loads((args.job_dir / 'job.json').read_text())
            sample = dict(elapsed_seconds=elapsed, stage=stage, status=meta.get('status'),
                          progress=meta.get('progress'), tree_rss_bytes=rss,
                          worker_rss_bytes=worker_rss, cpu_seconds_delta=cpu_seconds,
                          active_cpu_cores=cores, threads=threads, processes=detail)
            cgroup = cgroup_sample()
            if cgroup:
                sample['cgroup'] = cgroup
                peak_cgroup_working_set = max(peak_cgroup_working_set, cgroup['working_set_bytes'])
            raw.write(json.dumps(sample) + '\n')
            raw.flush()
            phase = phases.setdefault(stage, dict(seconds=0., cpu_seconds=0., peak_tree_rss_bytes=0,
                                                 peak_worker_rss_bytes=0, peak_active_cpu_cores=0.,
                                                 peak_threads=0, samples=0))
            phase['seconds'] += duration
            phase['cpu_seconds'] += cpu_seconds
            phase['peak_tree_rss_bytes'] = max(phase['peak_tree_rss_bytes'], rss)
            phase['peak_worker_rss_bytes'] = max(phase['peak_worker_rss_bytes'], worker_rss)
            phase['peak_active_cpu_cores'] = max(phase['peak_active_cpu_cores'], cores)
            phase['peak_threads'] = max(phase['peak_threads'], threads)
            phase['samples'] += 1
            for info in phases.values():
                info['average_active_cpu_cores'] = info['cpu_seconds'] / info['seconds'] if info['seconds'] else 0
            report = dict(job_id=meta['job_id'], server_pid=server.pid, elapsed_seconds=elapsed,
                          status=meta.get('status'), message=meta.get('message'), phases=phases,
                          peak_tree_rss_gib=max(v['peak_tree_rss_bytes'] for v in phases.values()) / 1024**3,
                          raw_samples=str(raw_path), platform='macOS' if psutil.MACOS else 'Linux',
                          measurement='sampled process-tree RSS; shared pages may be counted more than once',
                          cgroup_limit_measured=bool(cgroup),
                          cgroup_capacity_verified=bool(cgroup and meta.get('status') == 'completed'
                              and cgroup.get('memory.events', {}).get('oom_kill') == 0))
            if cgroup:
                report.update(cgroup=cgroup, peak_cgroup_working_set_gib=peak_cgroup_working_set / 1024**3)
            args.report.write_text(json.dumps(report, indent=2) + '\n')
            if meta.get('status') in {'completed', 'failed'}:
                terminal_since = terminal_since or elapsed
                if elapsed - terminal_since >= 5:
                    break
            time.sleep(args.interval)
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
