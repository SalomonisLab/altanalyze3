"""cellHarmony_differential splicing-heatmap: the fixed-order heatmap built from long-read dPSI covariate tables."""
import json
import sys

import anndata as ad
import numpy as np
import pandas as pd

from altanalyze3.components.cellHarmony import cellHarmony_differential as chd

STATES = ["S1", "S2", "S3", "S4"]
AGED = ["A1", "A2", "A3", "A4", "A5"]
YOUNG = ["Y1", "Y2", "Y3"]


def _cohort(tmp_path):
    rng = np.random.default_rng(7)
    events = [f"G{i}:E1.1-E2.1|G{i}:E1.1-E3.1" for i in range(120)]
    psi_rows, names = [], []
    for st in STATES:
        for u in AGED + YOUNG:
            aged = u in AGED
            x = np.clip(rng.normal(0.5, 0.03, len(events)), 0, 1)
            if aged:
                x[0:5] += 0.3                          # every state: global
                if st in ("S1", "S2"):
                    x[5:10] -= 0.3                     # S1 and S2 only: co-regulated or multi-state
                if st == "S3":
                    x[10:15] += 0.3                    # S3 only: local
            # A1 and Y1 columns carry the library name, as the real PSI matrix does for some samples
            sid = f"{u}_lib" if u in ("A1", "Y1") else u
            psi_rows.append(np.clip(x, 0, 1)); names.append(f"{st}.{sid}-junction")
    psi = ad.AnnData(X=np.asarray(psi_rows, dtype=np.float32), obs=pd.DataFrame(index=names),
                     var=pd.DataFrame(index=events))
    run = tmp_path / "run"; (run / "dPSI-covariate-events").mkdir(parents=True)
    psi.write_h5ad(run / "psi_combined_pseudo_cluster.h5ad")
    from altanalyze3.components.oncosplice import metadataAnalysis as ma
    P = pd.DataFrame(psi.X.T, index=events, columns=names)
    for st in STATES:
        cols = [c for c in names if c.startswith(st + ".")]
        g = pd.DataFrame({"grp": ["aged" if c.split(".")[1][0] == "A" else "young" for c in cols]}, index=cols)
        r = ma.limmaCompute(P[cols], g, grpvar="grp", min_group_size=2)
        r = r.rename(columns={"mwuStat": "ebayesStat", "mwuSign": "ebayesSign", "mwuPval": "ebayesPval",
                              "mwuAdjPval": "ebayesAdjPval"})
        r.insert(0, "Feature", r.index); r.insert(0, "Gene_Symbol", [e.split(":")[0] for e in r.index])
        r.to_csv(run / "dPSI-covariate-events" / f"young-aged-{st}_stats-annotated.txt", sep="\t", index=False)
    meta = pd.DataFrame({"uid": AGED + YOUNG, "bam": [f"/x/{u}.bam" for u in AGED + YOUNG],
                         "library": [f"{u}_lib" for u in AGED + YOUNG], "reverse": "TRUE",
                         "groups": ["aged"] * 5 + ["young"] * 3})
    meta.to_csv(tmp_path / "groups.txt", sep="\t", index=False)
    return run


def test_splicing_heatmap_groups_and_outputs(tmp_path, monkeypatch):
    run = _cohort(tmp_path)
    out = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", ["chd", "splicing-heatmap", "--splicing_dir", str(run / "dPSI-covariate-events"),
                                      "--case", "aged", "--control", "young", "--groups", str(tmp_path / "groups.txt"),
                                      "--use_rawp", "--dpsi", "0.1", "--outdir", str(out)])
    store = chd.main()
    a = store["assigned_groups"].set_index("gene")
    glob_lab = chd._splicing_event_label("G0", "G0:E1.1-E2.1|G0:E1.1-E3.1")
    loc_lab = chd._splicing_event_label("G10", "G10:E1.1-E2.1|G10:E1.1-E3.1")
    assert a.loc[glob_lab, "group"] == "global"
    assert a.loc[loc_lab, "group"] == "local" and a.loc[loc_lab, "population_or_pattern"] == "S3"
    assert set(a.index) <= set(store["detailed_deg"].gene)
    assert (store["fold_matrix"].loc[glob_lab] > 0.2).all()
    heat = out / "heatmaps" / "heatmap_aged_vs_young_splicing_by_cell_state.pdf"
    assert heat.exists() and heat.stat().st_size > 1000
    par = json.load(open(out / "parameters_aged_vs_young_splicing.json"))
    assert par["populations"] == STATES


def test_splicing_heatmap_subset_keeps_run_grouping(tmp_path, monkeypatch):
    run = _cohort(tmp_path)
    out = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", ["chd", "splicing-heatmap", "--splicing_dir", str(run / "dPSI-covariate-events"),
                                      "--case", "aged", "--control", "young", "--groups", str(tmp_path / "groups.txt"),
                                      "--use_rawp", "--dpsi", "0.1", "--outdir", str(out)])
    full = chd.main()
    monkeypatch.setattr(sys, "argv", ["chd", "splicing-heatmap-subset", "--run_dir", str(out), "--features", "G0,G10",
                                      "--name", "two_genes", "--outdir", str(tmp_path / "sub")])
    sub = chd.main()
    a_full = full["assigned_groups"].set_index("gene")
    a_sub = sub["assigned_groups"].set_index("gene")
    assert set(a_sub.index) == {g for g in a_full.index if g.split(" ")[0] in ("G0", "G10")}
    assert (a_sub.group == a_full.loc[a_sub.index, "group"]).all()
    assert (tmp_path / "sub" / "heatmap_aged_vs_young_splicing_two_genes_by_cell_state.pdf").exists()
