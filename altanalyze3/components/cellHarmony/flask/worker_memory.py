"""Read live worker-tree RSS and the container's actual memory accounting."""
from pathlib import Path
import subprocess


def host_memory():
    """Linux allocation-domain available/total bytes, independent of cgroup cap.

    MemAvailable includes reclaimable memory; MemFree alone does not. Outside
    Linux return None instead of treating absent telemetry as zero free memory.
    """
    try:
        lines = Path('/proc/meminfo').read_text().splitlines()
    except FileNotFoundError:
        return None
    entries = {key.rstrip(':'): int(value) * 1024
               for key, value, *_ in (line.split() for line in lines)}
    if 'MemAvailable' not in entries or 'MemTotal' not in entries:
        return None
    available, total = entries['MemAvailable'], entries['MemTotal']
    if not 0 <= available <= total or total <= 0:
        raise ValueError('Host memory accounting is invalid.')
    return available, total


def process_memory():
    """Return pid -> (parent pid, RSS bytes), including worker descendants."""
    if Path('/proc/self/statm').exists():
        import os
        page_size = os.sysconf('SC_PAGE_SIZE')
        result = {}
        for directory in Path('/proc').iterdir():
            if not directory.name.isdigit():
                continue
            try:
                stat = (directory / 'stat').read_text().rsplit(')', 1)[1].split()
                resident = int((directory / 'statm').read_text().split()[1])
                result[int(directory.name)] = (int(stat[1]), resident * page_size)
            except (OSError, ValueError, IndexError):
                continue  # Processes can exit while being sampled.
        return result
    output = subprocess.check_output(['ps', '-axo', 'pid=,ppid=,rss='], text=True)
    return {int(pid): (int(parent), int(rss) * 1024)
            for pid, parent, rss in (line.split() for line in output.splitlines())}


def tree_rss(processes, root):
    children = {}
    for pid, (parent, _) in processes.items():
        children.setdefault(parent, []).append(pid)
    pending, seen, total = [root], set(), 0
    while pending:
        pid = pending.pop()
        if pid in seen:
            continue
        seen.add(pid)
        total += processes.get(pid, (0, 0))[1]
        pending.extend(children.get(pid, ()))
    return total


def _working_set(base, used):
    """Charge allocations and protected cache, not clean reclaimable file pages.

    Repeated h5ad reads promote cache onto the active LRU. Subtracting only
    inactive_file makes identical jobs fail depending on earlier disk traffic.
    Keep dirty/writeback, mapped, shared-memory and unevictable pages charged.
    Bound the deduction by both file type and LRU counters, which can differ.
    The kernel cgroup cap remains in force; this only changes our early guard.
    See docs.kernel.org/admin-guide/cgroup-v2.html#memory-interface-files.
    """
    try:
        stat = dict(line.split() for line in (base / 'memory.stat').read_text().splitlines())
        # v1 total_* counters include descendants, matching usage_in_bytes.
        def count(v2, v1=None):
            return max(0, int(stat.get('total_' + (v1 or v2), stat.get(v2, stat.get(v1 or v2, 0)))))
        inactive = count('inactive_file')
        active = count('active_file')
        protected = (count('shmem') + count('file_dirty', 'dirty')
                     + count('file_writeback', 'writeback')
                     + count('file_mapped', 'mapped_file') + count('unevictable'))
        if 'file' in stat or 'cache' in stat or 'total_cache' in stat:
            clean = max(0, count('file', 'cache') - protected)
            reclaimable = min(active + inactive, clean)
        else:
            # Older/partial telemetry cannot establish that active pages are safe.
            reclaimable = max(0, inactive - protected)
    except (OSError, ValueError):
        return used
    return max(used - reclaimable, 0)


def container_memory_details():
    """Read raw and guard usage plus cache counters for diagnostic job logs."""
    try:
        memberships = Path('/proc/self/cgroup').read_text().splitlines()
    except OSError:
        return None
    for entry in memberships:
        _, controllers, relative = entry.split(':', 2)
        if controllers == '':
            bases = [Path('/sys/fs/cgroup') / relative.lstrip('/'), Path('/sys/fs/cgroup')]
            names = ('memory.current', 'memory.max')
        elif 'memory' in controllers.split(','):
            bases = [Path('/sys/fs/cgroup/memory') / relative.lstrip('/'), Path('/sys/fs/cgroup/memory')]
            names = ('memory.usage_in_bytes', 'memory.limit_in_bytes')
        else:
            continue
        for base in bases:
            try:
                used, limit = ((base / name).read_text().strip() for name in names)
                if limit != 'max':
                    result = {'raw_bytes': int(used), 'guard_bytes': _working_set(base, int(used)),
                              'limit_bytes': int(limit)}
                    try:
                        stat = dict(line.split() for line in (base / 'memory.stat').read_text().splitlines())
                        for key in ('active_file', 'inactive_file', 'file', 'cache', 'shmem',
                                    'file_dirty', 'dirty', 'file_writeback', 'writeback', 'anon', 'rss'):
                            result[key] = int(stat.get('total_' + key, stat.get(key, 0)))
                    except (OSError, ValueError):
                        pass
                    return result
            except (OSError, ValueError):
                continue
    return None


def container_memory():
    """Cgroup guard usage and hard limit; None outside a limited cgroup."""
    details = container_memory_details()
    return (details['guard_bytes'], details['limit_bytes']) if details else None
