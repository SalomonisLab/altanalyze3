"""evidence_rings: counts drawn equal the input, reflected order aligns shared evidence, PDF is vector with live text."""
import os
import pandas as pd
import pytest

from altanalyze3.components.visualization import evidence_rings as er
from altanalyze3.components.visualization.figure_style import pdf_vector_report


def _calls():
    rows = [("p1", "A", 1, "", 100), ("p1", "A", 0, "PV", 30), ("p1", "A", 0, "PV+CNV", 20), ("p1", "A", 0, "CNV", 10),
            ("p1", "B", 1, "", 40), ("p1", "C", 1, "", 5), ("p1", "C", 0, "DEG", 5),
            ("p2", "D", 1, "", 7)]
    return pd.DataFrame(rows, columns=["panel", "item", "direct", "evidence", "count"])


@pytest.mark.parametrize("ring1,band", [("aligned", "curved"), ("centered", "sweep"), ("full", "sweep"), ("full", "straight")])
def test_counts_and_vector_pdf(tmp_path, ring1, band):
    out = str(tmp_path / "rings.pdf")
    stats = er.plot_evidence_rings(_calls(), out, ring1=ring1, band=band,
                                   shells=[("PV", "#d1352b", "passenger"), ("CNV", "#e08a1e", "copy number"), ("DEG", "#7b4ea3", "expression")])
    assert stats["p1"] == dict(direct=145, total=210, shells={"PV": 50, "CNV": 30, "DEG": 5})
    assert stats["p2"] == dict(direct=7, total=7, shells={"PV": 0, "CNV": 0, "DEG": 0})
    rep = pdf_vector_report(out)
    assert os.path.getsize(out) > 2000 and rep["raster_images"] == 0 and rep["truetype_fonts"] >= 1


def test_reflected_order_keeps_each_shell_in_one_run():
    combos = [(True, False), (True, True), (False, True), (False, False)]
    order = er.reflected_order(combos, 2)
    for k in range(2):
        flags = [m[k] for m in order]
        assert sum(1 for i in range(1, len(flags)) if flags[i] != flags[i - 1]) <= 2
    assert order.index((True, True)) + 1 == order.index((False, True))


def test_rejects_imputed_call_without_evidence_and_unknown_shell(tmp_path):
    bad = pd.DataFrame([("p", "A", 0, "", 3)], columns=["panel", "item", "direct", "evidence", "count"])
    with pytest.raises(ValueError):
        er.prepare_calls(bad)
    with pytest.raises(ValueError):
        er.plot_evidence_rings(_calls(), str(tmp_path / "x.pdf"), shells=[("PV", "#d1352b", "passenger")])
