import os
import sys
import threading
import time
from pathlib import Path

import pytest

from altanalyze3.components.cellHarmony.flask import tasks, worker_memory
from altanalyze3.components.cellHarmony.flask.job_manager import JobStore


def test_tree_rss_includes_descendants_once():
    assert worker_memory.tree_rss({1: (0, 10), 2: (1, 20), 3: (2, 30), 4: (0, 99)}, 1) == 60


def test_live_process_memory_includes_current_process():
    assert worker_memory.tree_rss(worker_memory.process_memory(), os.getpid()) > 0


@pytest.mark.parametrize('membership,base,names', [
    ('0::/service', '/sys/fs/cgroup/service', ('memory.current', 'memory.max')),
    ('0::/', '/sys/fs/cgroup', ('memory.current', 'memory.max')),
    ('4:memory:/service', '/sys/fs/cgroup/memory/service',
     ('memory.usage_in_bytes', 'memory.limit_in_bytes')),
])
def test_cgroup_accounting(tmp_path, monkeypatch, membership, base, names):
    files = {'/proc/self/cgroup': membership,
             f'{base}/{names[0]}': '123', f'{base}/{names[1]}': '456'}
    real_path = Path
    class VirtualPath:
        def __init__(self, path):
            self.path = str(path)
        def __truediv__(self, suffix):
            return VirtualPath(real_path(self.path) / suffix)
        def read_text(self):
            if self.path not in files:
                raise FileNotFoundError(self.path)
            return files[self.path]
    monkeypatch.setattr(worker_memory, 'Path', VirtualPath)
    assert worker_memory.container_memory() == (123, 456)


def test_default_web_policy(monkeypatch):
    from altanalyze3.components.cellHarmony.webapp.config import load_config
    for name in ('JOB_WORKERS', 'WORKER_MEMORY_LIMIT_GIB', 'TOTAL_MEMORY_LIMIT_GIB'):
        monkeypatch.delenv('CELLHARMONY_' + name, raising=False)
    config = load_config()
    assert config['JOB_WORKERS'] == 2
    assert config['WORKER_MEMORY_LIMIT_GIB'] == 15
    assert config['TOTAL_MEMORY_LIMIT_GIB'] == 27


@pytest.mark.parametrize('worker_gib,container_gib,limit_gib,blocked', [
    (14, 26, 30, False), (15, 20, 30, True), (16, 20, 30, True),
    (1, 27, 30, True), (1, 18, 20, True),
])
def test_admission_thresholds(tmp_path, monkeypatch, worker_gib, container_gib, limit_gib, blocked):
    runner = tasks.JobRunner(JobStore(tmp_path), tmp_path / 'registry.json')
    runner._active_workers[999] = object()
    gib = 1024**3
    monkeypatch.setattr(tasks, 'process_memory', lambda: {
        os.getpid(): (0, gib), 999: (os.getpid(), gib),
        1000: (999, int((worker_gib - 1) * gib)),
    })
    monkeypatch.setattr(tasks, 'container_memory', lambda: (int(container_gib * gib), int(limit_gib * gib)))
    try:
        assert bool(runner._memory_wait_reason()) is blocked
    finally:
        runner.executor.shutdown()


def test_memory_read_failure_holds_admission(tmp_path, monkeypatch):
    runner = tasks.JobRunner(JobStore(tmp_path), tmp_path / 'registry.json')
    monkeypatch.setattr(tasks, 'process_memory', lambda: {})
    try:
        assert 'unavailable' in runner._memory_wait_reason()
    finally:
        runner.executor.shutdown()


def wait_for(predicate):
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline:
        if predicate():
            return
        time.sleep(0.02)
    pytest.fail('Worker did not reach expected state')


@pytest.mark.parametrize('memory_block', [False, True])
def test_real_workers_overlap_or_wait_for_memory(tmp_path, monkeypatch, memory_block):
    """Exercise actual child processes, shared pipeline/DE slots, and wakeup."""
    store = JobStore(tmp_path / 'jobs')
    first, second = [store.create_job('human', 'test', None, [])['job_id'] for _ in range(2)]
    store.update_job(second, status='completed', differential={'status': 'queued'})
    runner = tasks.JobRunner(store, tmp_path / 'registry.json', max_workers=2, isolate_jobs=True)
    original_popen = tasks.subprocess.Popen
    child_code = '''
import sys, time
from pathlib import Path
from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
store = JobStore(Path(sys.argv[1]))
job, task = sys.argv[2:4]
Path(sys.argv[4]).write_text('started')
while not Path(sys.argv[5]).exists():
    time.sleep(.02)
if task == 'pipeline':
    store.update_job(job, status='completed', message='Done.')
else:
    store.update_job(job, differential={'status':'completed', 'message':'Done.'})
'''
    def start_child(cmd, **kwargs):
        if len(cmd) < 3 or cmd[1:3] != ['-m', 'altanalyze3.components.cellHarmony.flask.worker']:
            return original_popen(cmd, **kwargs)
        job, task = cmd[5:7]
        return original_popen([sys.executable, '-c', child_code, str(store.root), job, task,
                               str(tmp_path / job), str(tmp_path / 'release')], **kwargs)
    monkeypatch.setattr(tasks.subprocess, 'Popen', start_child)
    blocked = threading.Event()
    original_reason = runner._memory_wait_reason
    runner._memory_wait_reason = lambda: ('A running analysis has reached the per-worker memory threshold.'
                                         if blocked.is_set() else original_reason())
    try:
        runner.submit_pipeline(first)
        wait_for(lambda: (tmp_path / first).exists())
        if memory_block:
            blocked.set()
        runner.submit_differential(second)
        if memory_block:
            wait_for(lambda: 'threshold' in store.get_job(second)['differential']['message'])
            assert not (tmp_path / second).exists()
            assert store.get_job(second)['differential']['status'] == 'queued'
            blocked.clear()
            with runner._admission:
                runner._admission.notify_all()
        wait_for(lambda: (tmp_path / second).exists())
        assert len(runner._active_workers) == 2
        (tmp_path / 'release').touch()
        for future in runner._futures.values():
            future.result(timeout=10)
        assert not runner._active_workers
        assert store.get_job(first)['status'] == 'completed'
        assert store.get_job(second)['differential']['status'] == 'completed'
    finally:
        blocked.clear()
        (tmp_path / 'release').touch()
        with runner._admission:
            runner._admission.notify_all()
        runner.executor.shutdown(wait=True)


def test_linux_host_available_memory_is_independent_of_container_limit(monkeypatch):
    class VirtualPath:
        def __init__(self, value):
            assert value == '/proc/meminfo'
        def read_text(self):
            return 'MemTotal: 33554432 kB\nMemFree: 100 kB\nMemAvailable: 20971520 kB\nHugePages_Total: 0\n'
    monkeypatch.setattr(worker_memory, 'Path', VirtualPath)
    assert worker_memory.host_memory() == (20 * 1024**3, 32 * 1024**3)


@pytest.mark.parametrize('available_gib,blocked', [(20, False), (2, True), (1, True)])
def test_host_headroom_can_queue_despite_room_in_cgroup(tmp_path, monkeypatch, available_gib, blocked):
    runner = tasks.JobRunner(JobStore(tmp_path), tmp_path / 'registry.json')
    gib = 1024**3
    monkeypatch.setattr(tasks, 'process_memory', lambda: {os.getpid(): (0, gib)})
    monkeypatch.setattr(tasks, 'container_memory', lambda: (11 * gib, 30 * gib))
    monkeypatch.setattr(tasks, 'host_memory', lambda: (available_gib * gib, 32 * gib))
    try:
        reason = runner._memory_wait_reason()
        assert bool(reason) is blocked
        if blocked:
            assert 'host' in reason
    finally:
        runner.executor.shutdown()
