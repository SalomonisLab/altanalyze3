"""protein-outcome mode of long_read.isoform_state_switches: per-state counts of age-differential isoforms by
protein outcome, from covariate tables or from called switches."""
import json

import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.long_read import isoform_state_switches as iss


def _classes(tmp_path):
    rows = [("G1", "ENST1", 300, "Not-NMD", 300),   # same length as the gene's longest protein
            ("G1", "nov1", 120, "Not-NMD", 300),    # shorter coding
            ("G1", "nov2", 60, "NMD", 300),         # NMD
            ("G2", "ENST2", 0, "", 200),            # no protein
            ("G2", "nov3", 150, "Potential-NMD", 200),
            ("G2", "ENST3", 200, "Not-NMD", 200)]
    p = pd.DataFrame(rows, columns=["gene", "isoform", "protein_length", "nmd_status", "longest_gene_protein"])
    p["class"] = "x"; p["known_novel"] = "x"; p["prediction_rows"] = 1
    path = tmp_path / "classes.tsv"
    p.to_csv(path, sep="\t", index=False)
    return path


def _table(path, rows, mean_scale="log2"):
    cols = ["Feature", "Gene_Symbol", "N_aged", "N_young", "Mean_aged", "Mean_young", "DiffMeans_aged_m_young",
            "ebayesPval"]
    out = []
    for feat, na, ny, ma, my, p in rows:
        out.append([feat, feat.split(":")[0], na, ny, ma, my, ma - my, p])
    pd.DataFrame(out, columns=cols).to_csv(path, sep="\t", index=False)


def test_classes_outcome(tmp_path):
    c = iss.protein_outcome_classes(_classes(tmp_path))
    assert c.loc["G1:ENST1", "outcome"] == "same length"
    assert c.loc["G1:nov1", "outcome"] == "protein coding"
    assert c.loc["G1:nov2", "outcome"] == "NMD"
    assert c.loc["G2:ENST2", "outcome"] == "no protein"
    assert c.loc["G2:nov3", "outcome"] == "NMD"          # Potential-NMD names NMD and is not Not-NMD
    assert c.loc["G2:ENST3", "outcome"] == "same length"


def test_from_tables_rules_and_counts(tmp_path):
    cls = iss.protein_outcome_classes(_classes(tmp_path))
    d = tmp_path / "diff"; d.mkdir()
    # S1: nov1 up 2-fold (log2 1.0), nov2 down, ENST1 up but P 0.2, nov3 up but only 2 young donors
    _table(d / "young-aged-S1_stats-annotated.txt",
           [("G1:nov1", 5, 3, 2.0, 1.0, 0.01), ("G1:nov2", 5, 3, 0.0, 1.0, 0.02), ("G1:ENST1", 5, 3, 3.0, 1.0, 0.2),
            ("G2:nov3", 5, 2, 3.0, 1.0, 0.001), ("G2:ENST3", 5, 3, 1.0, 0.1, 0.001)])
    _table(d / "young-aged-S2_stats-annotated.txt", [("G1:nov1", 4, 4, 1.0, 1.4, 0.01)])   # |0.4| < log2(1.5)
    t, man = iss.protein_outcome_from_tables(str(d), "aged", "young", cls, fold=1.5, rawp=0.05, min_donors=3,
                                             log=lambda *a: None)
    got = {(r.state, r.feature): (r.direction, r.outcome) for r in t.itertuples()}
    assert got == {("S1", "G1:nov1"): ("higher in aged", "protein coding"),
                   ("S1", "G1:nov2"): ("higher in young", "NMD"),
                   ("S1", "G2:ENST3"): ("higher in aged", "same length")}
    c = iss.protein_outcome_counts(t, "aged", "young", cluster_order=["S2", "S1"])
    assert int(c.isoforms.sum()) == 3
    assert list(dict.fromkeys(c.state)) == ["S1"]
    out = tmp_path / "p.pdf"
    iss.plot_protein_outcome(c, str(out), "aged", "young")
    assert out.stat().st_size > 1000


def test_ratio_scale_and_unknown_feature(tmp_path):
    cls = iss.protein_outcome_classes(_classes(tmp_path))
    d = tmp_path / "ratio"; d.mkdir()
    # usage 0.30 vs 0.15 is 2-fold (passes ratio scale); 0.30 vs 0.25 is 1.2-fold (fails)
    _table(d / "young-aged-S1_stats-annotated.txt", [("G1:nov1", 5, 3, 0.30, 0.15, 0.01), ("G1:nov2", 5, 3, 0.30, 0.25, 0.01)])
    t, _ = iss.protein_outcome_from_tables(str(d), "aged", "young", cls, fold=1.5, rawp=0.05, min_donors=3,
                                           fold_scale="ratio", log=lambda *a: None)
    assert list(t.feature) == ["G1:nov1"]
    _table(d / "young-aged-S1_stats-annotated.txt", [("G9:missing", 5, 3, 2.0, 1.0, 0.01)])
    with pytest.raises(ValueError, match="absent from the protein classes"):
        iss.protein_outcome_from_tables(str(d), "aged", "young", cls, log=lambda *a: None)


def test_from_switches_and_cli(tmp_path):
    cp = _classes(tmp_path)
    sw = pd.DataFrame({"state": ["S1"], "gene_symbol": ["G1"], "isoform_higher_in_A": ["G1:nov1"],
                       "isoform_higher_in_B": ["G1:ENST1"]})
    sp = tmp_path / "isoform_state_switches.tsv"
    sw.to_csv(sp, sep="\t", index=False)
    out = tmp_path / "out"
    iss.main(["protein-outcome", "--compare_groups", "aged,young", "--protein_classes", str(cp), "--switches", str(sp),
              "--outdir", str(out)])
    c = pd.read_csv(out / "protein_outcome_counts.tsv", sep="\t")
    pick = c.set_index(["direction", "outcome"]).isoforms
    assert pick[("higher in aged", "protein coding")] == 1 and pick[("higher in young", "same length")] == 1
    assert int(c.isoforms.sum()) == 2
    par = json.load(open(out / "protein_outcome_parameters.json"))
    assert par["isoforms"] == 2 and (out / "protein_outcome_by_state.pdf").exists()


def test_absolute_usage_shift_is_distinct_from_fold_ratio(tmp_path):
    cls=iss.protein_outcome_classes(_classes(tmp_path));d=tmp_path/'usage';d.mkdir()
    _table(d/'young-aged-S1_stats-annotated.txt',[
        ('G1:nov1',3,3,.8,.6,.01),      # delta .2 but ratio <1.5: requested shift passes
        ('G1:nov2',3,3,.03,.01,.01),    # 3-fold but delta .02: does not pass
        ('G2:nov3',3,3,.125,.025,.01),  # exact .1 boundary: does not pass
        ('G2:ENST3',3,2,.7,.4,.01)])    # insufficient biological observations
    table,_=iss.protein_outcome_from_tables(d,'aged','young',cls,fold_scale='ratio',
                                          min_usage_shift=.1,log=lambda *a:None)
    assert table.feature.tolist()==['G1:nov1']
    with pytest.raises(ValueError,match='requires ratio'):
        iss.protein_outcome_from_tables(d,'aged','young',cls,min_usage_shift=.1)
