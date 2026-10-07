"""Cold plot preparation never occupies a proxy request or duplicates a build."""
import json
import threading
import time
from test_grn_web_release import uploaded

from fastapi import HTTPException
from fastapi.responses import JSONResponse

from altanalyze3.components.cellHarmony.webapp.plot_build import PlotBuildQueue


def ready(queue, key, build):
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        response = queue.request(key, build)
        if response.status_code != 202:
            return response
        time.sleep(.01)
    raise AssertionError("Plot worker did not finish")


def test_pending_deduplicates_and_disk_result_expires():
    clock = [0.]
    queue = PlotBuildQueue(max_entries=2, ttl_seconds=10, clock=lambda: clock[0])
    release = threading.Event()
    calls = []

    def build():
        calls.append(1)
        assert release.wait(5)
        return JSONResponse({"values": [1., 2.], "cells": ["α", "b"]})

    try:
        start = time.monotonic()
        assert queue.request("job:filters", build).status_code == 202
        assert queue.request("job:filters", build).status_code == 202
        assert time.monotonic() - start < .2
        release.set()
        payload = ready(queue, "job:filters", build)
        assert json.loads(payload.body) == {"values": [1., 2.], "cells": ["α", "b"]}
        assert len(calls) == 1
        assert queue.request("job:filters", build).body == payload.body
        clock[0] = 11
        assert ready(queue, "job:filters", build).body == payload.body
        assert len(calls) == 2
    finally:
        release.set()
        queue.close()


def test_queue_capacity_and_errors():
    queue = PlotBuildQueue(max_entries=1)
    release = threading.Event()
    calls = []

    def blocked():
        calls.append("first")
        assert release.wait(5)
        return JSONResponse({"ok": True})

    def missing():
        calls.append("second")
        raise HTTPException(404, "No cells meet these display filters")

    try:
        assert queue.request("first", blocked).status_code == 202
        assert queue.request("second", missing).status_code == 202
        assert "second" not in calls
        release.set()
        ready(queue, "first", blocked)
        error = ready(queue, "second", missing)
        assert error.status_code == 404
        assert json.loads(error.body)["detail"] == "No cells meet these display filters"
        assert queue.request("second", missing).status_code == 404
        assert calls == ["first", "second"]
    finally:
        release.set()
        queue.close()


def test_route_defers_cache_loading_and_serializes_once(uploaded, monkeypatch):
    from importlib import import_module
    w = import_module("altanalyze3.components.cellHarmony.webapp.app")
    original = w._get_expression_cache
    release = threading.Event()
    called = []

    def slow(*args, **kwargs):
        called.append(1)
        assert release.wait(5)
        return original(*args, **kwargs)

    monkeypatch.setattr(w, "_get_expression_cache", slow)
    url = f"/api/jobs/{uploaded.job}/combplot"
    params = {"genes": "TF1", "deferred": "true", "cells_per_sample": 0}
    try:
        start = time.monotonic()
        assert uploaded.client.get(url, params=params).status_code == 202
        assert time.monotonic() - start < .5
        release.set()
        for _ in range(100):
            result = uploaded.client.get(url, params=params)
            if result.status_code != 202:
                break
            time.sleep(.01)
        assert result.status_code == 200, result.text
        assert len(called) == 1
        legacy = uploaded.client.get(url, params={**params, "deferred": "false"})
        assert result.json() == legacy.json()
        assert uploaded.client.get(url, params=params).content == result.content
        assert len(called) == 2
    finally:
        release.set()
        uploaded.app.state.plot_build_queue.close()
