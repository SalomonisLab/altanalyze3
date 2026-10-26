"""Unit checks for the Figure 2 visualization modules: subcluster_age_distance, umap_donor_groups,
cell_frequency_boxplot and figure_style. Synthetic data only; the real-data reproduction logs live
beside the Figure 2 outputs (see the Figure2 notebooks)."""
import numpy as np
import pandas as pd
import pytest
from scipy.stats import chi2_contingency, fisher_exact

from altanalyze3.components.visualization import figure_style as fs
from altanalyze3.components.visualization import subcluster_age_distance as sad
from altanalyze3.components.visualization import umap_donor_groups as udg
from altanalyze3.components.visualization import cell_frequency_boxplot as cfb


def test_two_by_two_p_switches_between_fisher_and_chi2():
    small = [[2, 50], [30, 40]]
    p, name = sad.two_by_two_p(small)
    assert name == "fisher_exact" and p == pytest.approx(fisher_exact(small)[1])
    big = [[20, 50], [30, 40]]
    p, name = sad.two_by_two_p(big)
    assert name == "chi2_yates" and p == pytest.approx(chi2_contingency(big)[1])


def _synthetic_cells():
    rows = []
    for sub, donors in (("0", ["A", "B", "C", "D", "E"]), ("1", ["F", "G", "H", "I", "J"])):
        for d in donors:
            grp = "young" if d in "ABCDE" else "old"
            for i in range(20):
                rows.append({"uid": f"{d}_{sub}_{i}", "donor": d, "group": grp, "subcluster": sub})
    rows.append({"uid": "A_1_x", "donor": "A", "group": "young", "subcluster": "1"})
    return pd.DataFrame(rows)


def test_valid_subclusters_and_counts():
    data = _synthetic_cells()
    assert sad.valid_subclusters(data, 4, 0.9) == ["0", "1"]
    assert sad.valid_subclusters(data, 6, 0.9) == ["1"]
    t = sad.test_subclusters(data, ["0", "1"], "young", "old", 0.025)
    r0 = t.set_index("subcluster").loc["0"]
    assert r0["cells_young"] == 100 and r0["cells_old"] == 0
    assert r0["donors_young"] == 5 and r0["donors_old"] == 0
    assert r0["donors_total_young"] == 5 and r0["donors_total_old"] == 5
    r1 = t.set_index("subcluster").loc["1"]
    # donor A has 1 of 101 cells in subcluster 1: below the 2.5% donor floor, so not counted
    assert r1["donors_young"] == 0 and r1["donors_old"] == 5
    assert r1["cells_young"] == 1


def test_centroid_correlation_two_subclusters():
    data = _synthetic_cells()
    rng = np.random.default_rng(0)
    genes = [f"g{i}" for i in range(30)]
    expr = pd.DataFrame(rng.random((30, len(data))), index=genes, columns=data["uid"])
    expr.loc[:, data.loc[data["subcluster"] == "1", "uid"]] += np.linspace(0, 3, 30)[:, None]
    rho, corr = sad.centroid_correlation(expr, data)
    c0 = expr[data.loc[data["subcluster"] == "0", "uid"]].mean(axis=1)
    c1 = expr[data.loc[data["subcluster"] == "1", "uid"]].mean(axis=1)
    assert rho == pytest.approx(np.corrcoef(c0, c1)[0, 1])
    assert np.isnan(corr.values[0, 0])


def test_benjamini_hochberg_matches_reference():
    # ranks: 0.01 -> 1, 0.03 -> 2, 0.04 -> 3, 0.5 -> 4; q = min over higher ranks of p*n/rank
    p = np.array([0.01, 0.04, 0.03, 0.5])
    q = cfb.benjamini_hochberg(p)
    assert q == pytest.approx([0.04, 0.04 * 4 / 3, 0.04 * 4 / 3, 0.5])


def test_frequency_test_table_counts_zero_variance():
    rows = []
    for grp, samples in (("young", ["y1", "y2", "y3"]), ("aged", ["a1", "a2", "a3"])):
        for s in samples:
            rows.append({"Grp": grp, "Sample": s, "CellState": "flat", "pct": 1.0})
            rows.append({"Grp": grp, "Sample": s, "CellState": "var",
                         "pct": 5.0 if grp == "aged" else 1.0 + 0.1 * len(s)})
    T, flat = cfb.test_table(pd.DataFrame(rows), "CellState", "X", "young", "aged", 3)
    assert flat == 1 and len(T) == 1
    assert T.iloc[0]["higher_in"] == "aged" and T.iloc[0]["rank_biserial"] == pytest.approx(1.0)


def test_order_states_keeps_hematopoietic_order_and_appends_unknown():
    out = fs.order_states(["MPP-1", "zzz", "HSC-1", "Mac"])
    assert out == ["HSC-1", "MPP-1", "Mac", "zzz"]


def test_shade_ramp_and_short_donor():
    cols = udg.shade_ramp(udg.BLUE_RAMP, 4)
    assert len(cols) == 4 and all(c.startswith("#") and len(c) == 7 for c in cols)
    assert len(set(cols)) == 4
    assert udg.short_donor("WM72-N-1130") == "WM72" and udg.short_donor("BF21") == "BF21"


def test_read_umap_scores_and_donor_colours(tmp_path):
    p = tmp_path / "x-UMAP_scores.txt"
    p.write_text("c1\t0.5\t1.5\nc2\t-1\t2\n")
    df = udg.read_umap_scores(p)
    assert list(df.columns) == ["uid", "x", "y"] and df["x"].tolist() == [0.5, -1.0]
    cells = pd.DataFrame({"label": ["d1", "d2", "d3"], "group": ["young", "old", "old"]})
    col = udg.donor_colours(cells, ["young", "old"], ["d3", "d2"])
    assert set(col) == {"d1", "d2", "d3"}
    assert col["d3"] != col["d2"]


def test_vector_pdf_has_no_raster_and_live_text(tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fs.apply_rcparams(7)
    fig, ax = plt.subplots(figsize=(2, 2))
    fs.vector_circles(ax, [0, 1, 2], [0, 1, 2], 3.0, ["#FF0000", "#00FF00", "#0000FF"])
    fs.vector_triangles(ax, [0, 1], [1, 0], 3.0, "#000000")
    ax.set_title("t")
    fs.strip_masks(fig)
    out = tmp_path / "v.pdf"
    fig.savefig(out, format="pdf")
    rep = fs.pdf_vector_report(out)
    assert rep["raster_images"] == 0 and rep["truetype_fonts"] >= 1
