"""protein-outcome predominance annotation and --states filter (long_read.isoform_state_switches): the predominant
isoform per state, group and gene, predominance switches of differential isoforms, and the longer / shorter call."""
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from altanalyze3.components.long_read import isoform_state_switches as iss

YOUNG, AGED = ["Y1", "Y2", "Y3"], ["O1", "O2", "O3"]
FEATS = ["G1:ENST1", "G1:nov1", "G1:nov2", "G2:ENST2", "G2:nov3", "G2:ENST3", "GF:f"]
# reads per (state, group): G1 switches in S1 (aged nov1, young ENST1); G2 keeps ENST3 predominant
READS = {("S1", "aged"): [20, 80, 5, 0, 10, 60, 1_000_000], ("S1", "young"): [80, 15, 10, 0, 10, 60, 1_000_000],
         ("S2", "aged"): [70, 30, 5, 0, 10, 60, 1_000_000], ("S2", "young"): [75, 25, 5, 0, 10, 60, 1_000_000]}


def _classes(tmp_path):
    rows = [("G1", "ENST1", 300, "Not-NMD", 300), ("G1", "nov1", 120, "Not-NMD", 300), ("G1", "nov2", 60, "NMD", 300),
            ("G2", "ENST2", 0, "", 200), ("G2", "nov3", 150, "Potential-NMD", 200), ("G2", "ENST3", 200, "Not-NMD", 200)]
    p = pd.DataFrame(rows, columns=["gene", "isoform", "protein_length", "nmd_status", "longest_gene_protein"])
    p["class"] = "x"; p["known_novel"] = "x"; p["prediction_rows"] = 1
    path = tmp_path / "classes.tsv"
    p.to_csv(path, sep="\t", index=False)
    return path


def _cohort(tmp_path):
    rows, names = [], []
    for k, uid in enumerate(YOUNG + AGED):
        g = "young" if uid in YOUNG else "aged"
        for s in ("S1", "S2"):
            x = np.array(READS[(s, g)], dtype=float)
            x[:-1] = x[:-1] + k % 2          # donor-to-donor variation that keeps the predominant isoforms
            rows.append(x); names.append(f"{s}.{uid}-isoform")
    a = ad.AnnData(sparse.csr_matrix(np.array(rows)), obs=pd.DataFrame(index=names), var=pd.DataFrame(index=FEATS))
    h5 = tmp_path / "counts.h5ad"
    a.write_h5ad(h5)
    meta = tmp_path / "groups.txt"
    pd.DataFrame({"uid": YOUNG + AGED, "groups": ["young"] * 3 + ["aged"] * 3}).to_csv(meta, sep="\t", index=False)
    return h5, meta, np.array(rows), names


def _by_hand(rows, names):
    """Mean donor usage (isoform / gene reads per donor) and the arg-max isoform, per state, group and gene."""
    out = {}
    for s in ("S1", "S2"):
        for g, donors in (("young", YOUNG), ("aged", AGED)):
            X = np.array([rows[names.index(f"{s}.{u}-isoform")] for u in donors])
            for gene, idx in (("G1", [0, 1, 2]), ("G2", [3, 4, 5])):
                U = X[:, idx] / X[:, idx].sum(axis=1, keepdims=True)
                out[(s, g, gene)] = FEATS[idx[int(np.argmax(U.mean(axis=0)))]]
    return out


def test_predominant_isoforms_match_by_hand(tmp_path):
    h5, meta, rows, names = _cohort(tmp_path)
    pred = iss.predominant_isoforms(str(h5), str(meta), ["aged", "young"], genes=["G1", "G2"], log=lambda *a: None)
    got = {(r.state, r.group, r.gene): r.predominant for r in pred.itertuples()}
    assert got == _by_hand(rows, names)
    assert got[("S1", "aged", "G1")] == "G1:nov1" and got[("S1", "young", "G1")] == "G1:ENST1"
    assert set(pred.donors) == {3}


def test_annotation_counts_and_cli(tmp_path):
    h5, meta, rows, names = _cohort(tmp_path)
    cp = _classes(tmp_path)
    d = tmp_path / "diff"; d.mkdir()
    cols = ["Feature", "Gene_Symbol", "N_aged", "N_young", "Mean_aged", "Mean_young", "DiffMeans_aged_m_young", "ebayesPval"]
    # S1: nov1 up in aged (switch, 120 aa < ENST1 300 aa -> shorter), ENST1 up in young (switch partner, same length as
    # the longest -> not plotted), nov2 up in young (not predominant -> no switch), ENST3 up in aged (no switch)
    s1 = [("G1:nov1", 3, 3, 2.0, 1.0), ("G1:ENST1", 3, 3, 1.0, 2.0), ("G1:nov2", 3, 3, 0.0, 1.0), ("G2:ENST3", 3, 3, 2.0, 1.0)]
    s2 = [("G1:nov1", 3, 3, 2.0, 1.0)]                    # up in aged but ENST1 is predominant in both groups
    s3 = [("G1:nov1", 3, 3, 2.0, 1.0)]                    # dropped by --states
    for st, rr in (("S1", s1), ("S2", s2), ("S3", s3)):
        pd.DataFrame([[f, f.split(":")[0], na, ny, ma, my, ma - my, 0.01] for f, na, ny, ma, my in rr],
                     columns=cols).to_csv(d / f"young-aged-{st}_stats-annotated.txt", sep="\t", index=False)
    states = tmp_path / "states.txt"; states.write_text("S1\nS2\nS4\n")
    out = tmp_path / "out"
    iss.main(["protein-outcome", "--compare_groups", "aged,young", "--protein_classes", str(cp), "--diff_dir", str(d),
              "--states", str(states), "--counts", str(h5), "--groups", str(meta), "--outdir", str(out)])
    t = pd.read_csv(out / "protein_outcome_isoforms.tsv", sep="\t", keep_default_na=False)
    assert set(t.state) == {"S1", "S2"}
    call = {(r.state, r.feature): (r.predominance_switch, r.switch_length) for r in t.itertuples()}
    assert call == {("S1", "G1:nov1"): (True, "shorter"), ("S1", "G1:ENST1"): (True, "longer"),
                    ("S1", "G1:nov2"): (False, ""), ("S1", "G2:ENST3"): (False, ""), ("S2", "G1:nov1"): (False, "")}
    c = pd.read_csv(out / "protein_outcome_counts.tsv", sep="\t")
    assert list(dict.fromkeys(c.state)) == ["S1", "S2", "S4"]          # every --states state on the x axis
    sw = pd.read_csv(out / "protein_outcome_predominance_counts.tsv", sep="\t", keep_default_na=False)
    hit = sw[sw.isoforms > 0]
    assert len(hit) == 1 and int(hit.isoforms.sum()) == 1                         # ENST1 is 'same length': not plotted
    assert hit.iloc[0][["state", "direction", "outcome", "switch_length", "longer_protein_in"]].tolist() == \
        ["S1", "higher in aged", "protein coding", "shorter", "young"]
    par = json.load(open(out / "protein_outcome_parameters.json"))
    assert par["predominance_switch_totals"]["higher in aged|shorter"] == 1
    assert (out / "protein_outcome_by_state.pdf").stat().st_size > 1000
    # gain/loss of every plotted isoform against the other group's predominant isoform (ENST1 in young, S1 and S2)
    lv = {(r.state, r.feature): r.length_vs_other_group for r in t.itertuples()}
    assert lv[("S1", "G1:nov1")] == "shorter" and lv[("S2", "G1:nov1")] == "shorter"
    gl = pd.read_csv(out / "protein_outcome_gain_loss_counts.tsv", sep="\t").set_index(["state", "direction", "category"])
    assert gl.isoforms[("S1", "higher in aged", "loss")] == 1 and gl.isoforms[("S2", "higher in aged", "loss")] == 1
    assert gl.isoforms[("S1", "higher in young", "NMD")] == 1 and int(gl.isoforms.sum()) == 3
    assert (out / "protein_outcome_gain_loss_by_state.pdf").stat().st_size > 1000
