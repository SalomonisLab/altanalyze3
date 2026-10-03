"""Read live worker-tree RSS and the container's actual memory accounting."""
from pathlib import Path
import subprocess


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
    """Usage minus inactive file cache, as docker stats and the kubelet count it.

    memory.current includes page cache from reading and writing h5ad files. The
    kernel frees it only under pressure, so after a few large jobs an idle
    container could sit above the admission ceiling and hold every new job
    (2026-10-03: 0.17 GiB anon + 3.4 GiB cache held admission indefinitely).
    """
    try:
        stat = dict(line.split() for line in (base / 'memory.stat').read_text().splitlines())
        inactive = int(stat.get('inactive_file', stat.get('total_inactive_file', 0)))
    except (OSError, ValueError):
        return used
    return max(used - inactive, 0)


def container_memory():
    """Read our cgroup v2/v1 usage and limit; return None outside a cgroup."""
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
                    return _working_set(base, int(used)), int(limit)
            except (OSError, ValueError):
                continue
    return None
