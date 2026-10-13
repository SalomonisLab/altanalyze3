"""struct_rows_geometry: the bottom isoform-structure track's row positions.

Two contracts:
  1. labels_above=False reproduces the fixed ISO_ROW grid the renderer used before the figure preset
     existed, centre for centre, so every non-figure-mode PDF is byte-identical.
  2. labels_above=True (figure preset) puts each label CLOSER to the structure it names than the
     structures are to each other, and the last row of a multi-row track labels below its structure.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from altanalyze3.components.visualization.isv_web.isv_web_remote import (  # noqa: E402
    ISO_EXH, ISO_LBL_GAP, ISO_LBL_H, ISO_ROW, ISO_STRUCT_GAP, struct_rows_geometry)


@pytest.mark.parametrize("n", [0, 1, 2, 3, 7])
@pytest.mark.parametrize("fs", [1.0, 1.9, 2.6])
def test_legacy_grid_unchanged(n, fs):
    rows, h = struct_rows_geometry(n, fs, labels_above=False, y0=100.0)
    assert len(rows) == n
    assert h == pytest.approx(n * ISO_ROW * fs)
    for i, (rc, ly, va) in enumerate(rows):
        assert rc == pytest.approx(100.0 + i * ISO_ROW * fs + ISO_ROW * fs / 2.0)
        assert ly is None and va is None


def test_empty_track():
    assert struct_rows_geometry(0, 2.6, labels_above=True) == ([], 0.0)


@pytest.mark.parametrize("fs", [1.0, 1.9, 2.6])
def test_single_row_labels_above(fs):
    rows, h = struct_rows_geometry(1, fs, labels_above=True, y0=0.0)
    (rc, ly, va), = rows
    assert va == "bottom"                      # one row keeps its label above
    assert ly < rc
    assert h == pytest.approx(ISO_LBL_H * fs + ISO_LBL_GAP * fs + ISO_EXH * fs)


@pytest.mark.parametrize("n", [2, 3, 5])
@pytest.mark.parametrize("fs", [1.0, 1.9, 2.6])
def test_label_is_nearer_its_own_structure(n, fs):
    rows, h = struct_rows_geometry(n, fs, labels_above=True, y0=50.0)
    exh = ISO_EXH * fs
    assert len(rows) == n
    assert [va for _, _, va in rows] == ["bottom"] * (n - 1) + ["top"]     # last row labels below
    assert all(b[0] > a[0] for a, b in zip(rows, rows[1:]))               # rows run top to bottom
    gap_struct = min(b[0] - a[0] for a, b in zip(rows, rows[1:])) - exh   # between exon boxes
    for rc, ly, va in rows:
        edge = rc - exh / 2.0 if va == "bottom" else rc + exh / 2.0
        gap_lbl = abs(ly - edge)
        assert gap_lbl == pytest.approx(ISO_LBL_GAP * fs)
        assert gap_lbl < gap_struct                                       # the pairing the reader sees
    assert gap_struct == pytest.approx(ISO_STRUCT_GAP * fs)
    assert h == pytest.approx(n * (ISO_LBL_H + ISO_LBL_GAP + ISO_EXH) * fs + (n - 1) * ISO_STRUCT_GAP * fs)


@pytest.mark.parametrize("n", [1, 2, 4])
def test_height_covers_every_drawn_element(n):
    fs, y0 = 2.6, 17.0
    rows, h = struct_rows_geometry(n, fs, labels_above=True, y0=y0)
    exh, lbl = ISO_EXH * fs, ISO_LBL_H * fs
    top = min(min(rc - exh / 2.0, ly - lbl if va == "bottom" else ly) for rc, ly, va in rows)
    bot = max(max(rc + exh / 2.0, ly if va == "bottom" else ly + lbl) for rc, ly, va in rows)
    assert top >= y0 - 1e-9
    assert bot <= y0 + h + 1e-9
