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


def _state_dir(tmp_path, name="HSC-1"):
    """A sub-cluster folder carrying only the NMF round-1 labels, as the 2024 run wrote them."""
    d = tmp_path / name / "NMF-SVM" / "SubtypeAnalyses"
    d.mkdir(parents=True)
    (d / "round1Metadata2.txt").write_text("uid\tV0\tV1\nc1\t1\t0\nc2\t0\t1\nc3\t0\t1\n")
    return tmp_path / name


def test_annotation_dir_replaces_the_run_labels(tmp_path):
    """Nathan, 2026-11-02: the panels drew NMF round-1 labels; the SVM assignment must be usable."""
    state = _state_dir(tmp_path)
    assert list(sad.read_annotation(state, None, "nmf")["subcluster"]) == ["0", "1", "1"]
    alt = tmp_path / "svm_labels"
    alt.mkdir()
    (alt / "HSC-1.txt").write_text("uid\tV0\tV1\nc1\t1\t0\nc2\t1\t0\nc3\t0\t1\n")
    assert list(sad.read_annotation(state, alt, "nmf")["subcluster"]) == ["0", "0", "1"]
    path, has_header = sad.find_annotation_file(state, alt)
    assert path == alt / "HSC-1.txt" and has_header


def test_annotation_dir_missing_file_is_an_error_not_a_fall_back(tmp_path):
    state = _state_dir(tmp_path)
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError, match="is missing"):
        sad.read_annotation(state, empty, "nmf")


def _k2_state(tmp_path, name="HSC-1", decisions=(-2.0, 2.0, 0.05), nmf=("V0", "V1", "V1")):
    """A k=2 sub-cluster folder: NMF round-1 labels, SVM decision values, marker expression."""
    from altanalyze3.components.visualization import subcluster_labels as sl  # noqa: F401
    d = tmp_path / name
    (d / "NMF-SVM" / "SubtypeAnalyses").mkdir(parents=True)
    (d / "NMF-SVM" / "SVMOutputs").mkdir(parents=True)
    (d / "NMF-SVM" / "ExpressionInput").mkdir(parents=True)
    (d / "NMF-SVM" / "centroids").mkdir(parents=True)
    uids = [f"c{i}" for i in range(1, len(nmf) + 1)]
    rows = ["uid\tV0\tV1"] + [f"{u}\t{int(v == 'V0')}\t{int(v == 'V1')}" for u, v in zip(uids, nmf)]
    (d / "NMF-SVM" / "SubtypeAnalyses" / "round1Metadata2.txt").write_text("\n".join(rows) + "\n")
    dec = ["uid\tR1-V0_vs_Others.txt\tR1-V1_vs_Others.txt"]
    dec += [f"{u}\t{v}\t{v}" for u, v in zip(uids, decisions)]
    (d / "NMF-SVM" / "SVMOutputs" / "round1SVC_decision_func.txt").write_text("\n".join(dec) + "\n")
    # c1 and c3 share a profile; c2 is the opposite one, so the centroid pass moves c3 to V0
    expr = pd.DataFrame({"c1": [10.0, 0.0, 9.0], "c2": [0.0, 10.0, 1.0], "c3": [9.0, 1.0, 10.0]},
                        index=["g1", "g2", "g3"])
    expr.index.name = "uid"
    expr.to_csv(d / "NMF-SVM" / "ExpressionInput" / "exp.NMF-MarkerFinder.txt", sep="\t")
    (d / "NMF-SVM" / "centroids" / "exp.MF.txt").write_text("g1\ng2\ng3\n")
    return d


def test_svm_labels_recover_the_assignment_altanalyze2_never_wrote(tmp_path):
    """Nathan, 2026-11-02: the k=2 run saves only the decision values, so labels must come from them."""
    from altanalyze3.components.visualization import subcluster_labels as sl
    state = _k2_state(tmp_path)
    nmf = sl.read_nmf_labels(state)
    assert list(nmf["label"]) == ["V0", "V1", "V1"]
    labels, margin = sl.svm_labels(nmf, sl.read_svm_decision(state))
    assert list(labels) == ["V0", "V1", "V1"] and margin == 3


def test_svm_orientation_follows_the_overlap_not_the_sign(tmp_path):
    from altanalyze3.components.visualization import subcluster_labels as sl
    state = _k2_state(tmp_path, decisions=(2.0, -2.0, -0.05))
    labels, margin = sl.svm_labels(sl.read_nmf_labels(state), sl.read_svm_decision(state))
    assert list(labels) == ["V0", "V1", "V1"] and margin == 3


def test_centroid_passes_move_the_small_margin_cell_and_settle(tmp_path):
    from altanalyze3.components.visualization import subcluster_labels as sl
    state = _k2_state(tmp_path)
    table, row = sl.labels_for_state(state, "svm+centroid", lambda *a, **k: None)
    assert list(table["uid"]) == ["c1", "c2", "c3"]
    assert list(table["V0"]) == [1, 0, 1] and list(table["V1"]) == [0, 1, 0]
    assert row["n_moved_centroid"] == 1 and row["centroid_converged"]


def test_nmf_method_returns_the_run_labels_unchanged(tmp_path):
    from altanalyze3.components.visualization import subcluster_labels as sl
    state = _k2_state(tmp_path)
    table, row = sl.labels_for_state(state, "nmf", lambda *a, **k: None)
    assert list(table["V1"]) == [0, 1, 1]
    assert row["n_moved_svm"] == 0 and row["n_moved_centroid"] == 0


def test_rank_other_than_two_and_missing_decision_file_are_errors(tmp_path):
    from altanalyze3.components.visualization import subcluster_labels as sl
    state = _k2_state(tmp_path, name="S1")
    (state / "NMF-SVM" / "SubtypeAnalyses" / "round1Metadata2.txt").write_text(
        "uid\tV0\tV1\tV2\nc1\t1\t0\t0\nc2\t0\t1\t0\nc3\t0\t0\t1\n")
    with pytest.raises(ValueError, match="k=2 only"):
        sl.read_nmf_labels(state)
    bare = tmp_path / "S2" / "NMF-SVM" / "SubtypeAnalyses"
    bare.mkdir(parents=True)
    (bare / "round1Metadata2.txt").write_text("uid\tV0\tV1\nc1\t1\t0\n")
    with pytest.raises(FileNotFoundError, match="no SVM decision values"):
        sl.read_svm_decision(tmp_path / "S2")


def test_zero_variance_cells_keep_their_label_and_are_counted(tmp_path):
    """A cell with no variance over the MarkerFinder genes correlates with nothing; keep and count it."""
    from altanalyze3.components.visualization import subcluster_labels as sl
    expr = pd.DataFrame({"c1": [10.0, 0.0, 9.0], "c2": [0.0, 10.0, 1.0], "flat": [5.0, 5.0, 5.0]},
                        index=["g1", "g2", "g3"])
    labels = pd.Series({"c1": "V0", "c2": "V1", "flat": "V1"})
    out, passes, converged, undecidable = sl.centroid_passes(expr, labels)
    assert undecidable == 1 and converged
    assert out["flat"] == "V1" and len(out) == 3
    assert out.notna().all()


def test_panel_b_default_labels_are_the_svm_assignment_not_the_nmf_file(tmp_path):
    """Nathan, 2026-11-02: for k=2 the panels must default to the supervised assignment."""
    state = _k2_state(tmp_path)
    assert list(sad.read_annotation(state, None, "nmf")["subcluster"]) == ["0", "1", "1"]
    assert list(sad.read_annotation(state)["subcluster"]) == ["0", "1", "1"]
    refined = sad.read_annotation(state, None, "svm+centroid")["subcluster"]
    assert list(refined) == ["0", "1", "0"]
    alt = tmp_path / "given"
    alt.mkdir()
    (alt / "HSC-1.txt").write_text("uid\tV0\tV1\nc1\t0\t1\nc2\t0\t1\nc3\t0\t1\n")
    assert list(sad.read_annotation(state, alt, "svm")["subcluster"]) == ["1", "1", "1"]


def test_panel_b_and_c_cli_defaults_name_the_k2_method():
    b = sad.build_parser().parse_args(["--outdir", "x", "--plot-only"])
    assert b.labels == "svm"
    c = udg.build_parser().parse_args(
        ["plot", "--coordinates", "c", "--groups", "g", "--out", "o"])
    assert c.labels == "svm+centroid" and c.subcluster_dir is None


def test_panel_names_follow_relative_enrichment_not_the_majority():
    """Multilin-2: young cells dominate both sub-clusters, so 'majority' would label both Young."""
    cells = pd.DataFrame({
        "x": [0, 1, 2, 10, 11, 12], "y": [0, 1, 2, 10, 11, 12],
        "subcluster": ["0", "0", "0", "1", "1", "1"],
        "group": ["young", "young", "young", "young", "young", "old"]})
    labs = udg._panel_labels(cells, ["young", "old"], ["Young", "Aged"])
    by_sub = {d["subcluster"]: d for d in labs}
    assert by_sub["0"]["name"] == "Young" and by_sub["1"]["name"] == "Aged"
    assert by_sub["1"]["minor"] == "young (n=2)"
    assert by_sub["0"]["minor"] == ""


def test_panel_labels_without_subclusters_fall_back_to_the_age_groups():
    cells = pd.DataFrame({"x": [0, 1], "y": [0, 1], "group": ["young", "old"]})
    labs = udg._panel_labels(cells, ["young", "old"], ["Young", "Aged"])
    assert [d["name"] for d in labs] == ["Young", "Aged"]
    assert all(d["minor"] == "" for d in labs)


def test_plot_only_records_what_it_drew(tmp_path):
    """Nathan, 2026-11-02: a redraw changed the page height while parameters.json still said 4.0."""
    import json
    import types
    table = pd.DataFrame({"cell_type": ["HSC-1", "HSC-2"], "distance": [0.8, 0.3],
                          "neg_log10_p": [2.5, 1.2], "status": ["ok", "ok"]})
    table.to_csv(tmp_path / "age_induced_distance.tsv", sep="\t", index=False)
    args = types.SimpleNamespace(outdir=str(tmp_path), plot_only=True, hline=0.3, cmap="coolwarm",
                                 size_scale=20.0, alpha=0.7, fig_width=12.0, fig_height=2.8,
                                 group_a="young", group_b="old", order=None)
    sad.run(args)
    rec = json.loads((tmp_path / "plot_parameters.json").read_text())
    assert rec["fig_height"] == 2.8 and rec["pdf"]["raster_images"] == 0


def test_lmpp_states_sit_together_before_lmpp_2():
    """Nathan, 2026-11-02: LMPP-1 and LMPP-1-cycling belong right before LMPP-2."""
    order = fs.HEMATOPOIETIC_STATE_ORDER
    assert len(order) == len(set(order))
    i = order.index("LMPP-2")
    assert order[i - 2:i + 1] == ["LMPP-1", "LMPP-1-cycling", "LMPP-2"]
    drawn = fs.order_states(["LMPP-2", "MPP-MEP", "LMPP-1-cycling", "HSC-1", "LMPP-1", "Multilin-1"])
    assert drawn == ["HSC-1", "MPP-MEP", "Multilin-1", "LMPP-1", "LMPP-1-cycling", "LMPP-2"]


def test_combined_panel_labels_stay_live_text(tmp_path):
    """Nathan, 2026-11-02: a white halo turned the panel labels into outlines. Text must stay text."""
    import types
    state = _k2_state(tmp_path, name="S1", decisions=(-2.0, 2.0, 1.5), nmf=("V0", "V1", "V1"))
    coords = tmp_path / "coords.txt"
    coords.write_text("c1\t0.0\t0.0\nc2\t5.0\t5.0\nc3\t5.1\t5.2\n")
    groups = tmp_path / "groups.txt"
    groups.write_text("uid\tdonor\tage-group\nc1\tAM72-N\told\nc2\tWM22-N\tyoung\nc3\tWM34-N\tyoung\n")
    panels = tmp_path / "panels.tsv"
    panels.write_text(f"title\tcoordinates\tsubcluster_dir\nS1\t{coords}\t{state}\n")
    out = tmp_path / "fig.pdf"
    udg.combined(types.SimpleNamespace(
        panels=str(panels), groups=str(groups), out=str(out), uid_col="uid", donor_col="donor",
        group_col="age-group", group_order="young,old", group_labels="Young,Aged",
        legend_group_order=None, donor_order=None, short_donor=True, labels="svm",
        aspect="auto", point_size=1.1, panel_size=2.2, seed=0))
    report = fs.pdf_vector_report(out)
    assert report["raster_images"] == 0
    # title + two axis labels + two sub-cluster names + one minority count + two legend entries
    assert report["text_chars"] >= 40, report
    rows = pd.read_csv(f"{out.with_suffix('')}_panel_labels.tsv", sep="\t")
    assert set(rows["name"]) == {"Young", "Aged"}
