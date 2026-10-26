import pandas as pd
import pytest

fitz = pytest.importorskip("fitz")

from altanalyze3.components.goelite.plotting import write_goelite_scatter_pdf


def test_y_ticks_are_plain_live_text(tmp_path):
    fr = pd.DataFrame({"term_id": [f"GO:{i:07d}" for i in range(6)],
                       "term_name": ["a", "b", "c", "d", "e", "f"],
                       "z_score": [1.0, 2.5, 3.0, 5.0, 0.5, -1.0],
                       "p_value": [0.5, 0.01, 1e-4, 1e-8, 0.7, 0.9],
                       "fdr": [0.6, 0.02, 2e-4, 1e-7, 0.8, 0.95],
                       "selected": [False, True, True, True, False, False]})
    out = write_goelite_scatter_pdf(fr, tmp_path / "g.pdf", title_prefix="t")
    page = fitz.open(out)[0]
    words = {w[4] for w in page.get_text("words")}
    assert {"1e-7", "1e-5", "1e-3"} & words or {"1e-7", "1e-6", "1e-5"} & words
    assert "0.1" in words or "0.01" in words
    assert all(f[3].split("+")[-1].startswith("Arial") for f in page.get_fonts())
    assert not any(w.startswith("10") and len(w) <= 3 for w in words - {"10"})
