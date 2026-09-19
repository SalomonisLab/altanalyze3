"""A precomputed bundle must serve its DEG tables from wherever it is mounted.

The manifest records each comparison twice: `file`, relative to
`<bundle_dir>/<prefix>_deg/`, and `path`, the absolute path of the machine that built
the bundle. Reading `path` made a bundle built on a laptop unusable on a server, where
the viewer answered every differential request with
"Differential DEG detail table is unavailable."

Run:

    /usr/bin/python3 tests/test_scalable_bundle_portable_paths.py
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from altanalyze3.components.visualization.scalable_viewer.bundle_meta import (  # noqa: E402
    _resolve_deg_table,
)

failures: list[str] = []


class _Paths:
    def __init__(self, bundle_dir: str, prefix: str) -> None:
        self.bundle_dir = bundle_dir
        self.prefix = prefix


class _DS:
    def __init__(self, bundle_dir: str, prefix: str) -> None:
        self.paths = _Paths(bundle_dir, prefix)


def check(name: str, got, want) -> None:
    if got != want:
        failures.append(f"{name}: got {got!r}, want {want!r}")


with tempfile.TemporaryDirectory() as tmp:
    prefix = "Hs-Lung-COPD-metacells"
    bundle = os.path.join(tmp, "bundles_integrated", "COPD-metacells")
    deg = os.path.join(bundle, f"{prefix}_deg", "rna")
    os.makedirs(deg)
    table = os.path.join(deg, "cancer_vs_no_cancer.tsv")
    Path(table).write_text("gene\tlog2fc\nSFTPC\t1.0\n")

    ds = _DS(bundle, prefix)

    # 1. The server case: the recorded absolute path belongs to another machine.
    absent = "/Users/someone/Dropbox/LungMAP/does-not-exist/cancer_vs_no_cancer.tsv"
    check("relocated bundle resolves",
          _resolve_deg_table(ds, {"file": "rna/cancer_vs_no_cancer.tsv", "path": absent}),
          table)

    # 2. An older manifest with no `file` still uses the recorded path.
    check("legacy manifest falls back",
          _resolve_deg_table(ds, {"path": absent}), absent)

    # 3. A `file` that names nothing falls back rather than inventing a path.
    check("missing relative falls back",
          _resolve_deg_table(ds, {"file": "rna/not_here.tsv", "path": absent}), absent)

    # 4. Neither field present is an empty answer, which the caller reports as unavailable.
    check("no location at all", _resolve_deg_table(ds, {}), "")

    # 5. The resolved path must sit inside the bundle that is mounted now.
    got = _resolve_deg_table(ds, {"file": "rna/cancer_vs_no_cancer.tsv", "path": absent})
    check("resolved path is inside the mounted bundle", got.startswith(bundle), True)

if failures:
    print(f"FAIL: {len(failures)} check(s)")
    for line in failures:
        print("  " + line)
    sys.exit(1)
print("PASS: DEG tables resolve against the mounted bundle, with the recorded path as fallback")
