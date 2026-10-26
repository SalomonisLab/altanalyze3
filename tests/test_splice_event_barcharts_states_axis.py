"""splice_event_barcharts --states-axis: the x layout draws the same counts and writes the same significant
results as the default y layout."""
import os

import pandas as pd

from altanalyze3.components.visualization import splice_event_barcharts as seb

HEADER = ["Gene_Symbol", "Isoform_1|Length", "Isoform_2|Length", "Coord_1", "Coord_2", "Event_Type", "Event_Direction",
          "Trans_Splicing_Status", "Trans_Splicing_Evidence", "Feature", "N_aged", "N_young", "Median_aged",
          "Median_young", "Mean_aged", "Mean_young", "DiffMeans_aged_m_young", "ebayesStat", "ebayesSign",
          "ebayesPval", "ebayesAdjPval"]


def _write(path, rows):
    out = []
    for k, (etype, direction, dpsi, p, n) in enumerate(rows):
        out.append(["G", f"ENST{k}|100", f"ENST{k + 50}|90(NMD)", "c1", "c2", etype, direction, "", "", f"G:E{k}|G:E{k}b",
                    n, n, 0.5, 0.5, 0.5 + dpsi, 0.5, dpsi, 1, 1, p, 1])
    pd.DataFrame(out, columns=HEADER).to_csv(path, sep="\t", index=False)


def test_states_axis_x_matches_y(tmp_path):
    d = tmp_path / "events"; d.mkdir()
    _write(d / "young-aged-HSC-1_stats-annotated.txt",
           [("Cassette Exon", "inclusion", 0.3, 0.01, 3), ("Intron Retention", "exclusion", -0.2, 0.02, 4),
            ("Cassette Exon", "inclusion", 0.05, 0.01, 3), ("AltPromoter", "proximal", 0.4, 0.2, 3)])
    _write(d / "young-aged-MPP-1_stats-annotated.txt",
           [("Cassette Exon", "inclusion", 0.3, 0.01, 2), ("Trans-Splicing", "distal", 0.5, 0.001, 5)])
    files = sorted(str(x) for x in d.glob("*_stats-annotated.txt"))
    order = ["HSC-1", "MPP-1"]
    outs = {}
    for axis in ("y", "x"):
        o = tmp_path / axis
        seb.analyze_splicing_events(files, output_dir=str(o), pval_threshold=0.05, diff_means_threshold=0.1,
                                    sample_min=3, cluster_order=order, states_axis=axis)
        assert (o / "splicing_events_all.pdf").stat().st_size > 1000
        outs[axis] = {f: open(o / f).read() for f in sorted(os.listdir(o)) if f.endswith("_significant-results.tsv")}
    assert outs["x"] == outs["y"]
    hsc = pd.read_csv(tmp_path / "x" / "young-aged-HSC-1_significant-results.tsv", sep="\t")
    assert sorted(hsc.Event_Type) == ["Cassette Exon", "Intron Retention"]
    mpp = pd.read_csv(tmp_path / "x" / "young-aged-MPP-1_significant-results.tsv", sep="\t")
    assert list(mpp.Event_Type) == ["Trans-Splicing"]
