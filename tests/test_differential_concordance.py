"""cellHarmony differential_concordance: population concordance, grouping and unique hits from a heatmap TSV."""
import numpy as np
import pandas as pd

from altanalyze3.components.cellHarmony import differential_concordance as dc


def _inputs(tmp_path):
    rng = np.random.default_rng(1)
    pops = ["H1", "H2", "M1", "M2", "L1"]
    feats = [f"G{i} E{i}" for i in range(400)]
    base_h = rng.normal(0, 0.3, 400); base_m = rng.normal(0, 0.3, 400)
    vals = {"H1": base_h + rng.normal(0, 0.05, 400), "H2": base_h + rng.normal(0, 0.05, 400),
            "M1": base_m + rng.normal(0, 0.05, 400), "M2": base_m + rng.normal(0, 0.05, 400),
            "L1": rng.normal(0, 0.3, 400)}
    h = pd.DataFrame(vals, index=feats)
    h.iloc[390:, 2:] = np.nan                                  # last 10 features untested outside the H group
    tsv = h.copy(); tsv.index = [f"grp:{f}" for f in feats]; tsv.columns = [f"{p}:{p}" for p in pops]
    tsv.to_csv(tmp_path / "heat.tsv", sep="\t")
    rows = [dict(gene=f, population=p, log2fc=h.loc[f, p], sig_metric=0.01) for p in pops for f in feats
            if np.isfinite(h.loc[f, p]) and abs(h.loc[f, p]) > 0.25]
    pd.DataFrame(rows).to_csv(tmp_path / "det.tsv", sep="\t", index=False)
    return h


def test_concordance_groups_and_unique(tmp_path):
    h = _inputs(tmp_path)
    rep = dc.main(["--heatmap_tsv", str(tmp_path / "heat.tsv"), "--deg_detailed", str(tmp_path / "det.tsv"),
                   "--k", "2", "3", "--min_hits", "10", "--outdir", str(tmp_path / "out")])
    g3 = rep["groups"]["3"]["populations"]
    sets = sorted(sorted(v) for v in g3.values())
    assert ["H1", "H2"] in sets and ["M1", "M2"] in sets and ["L1"] in sets
    R = pd.read_csv(tmp_path / "out" / "concordance_pearson_r.tsv", sep="\t", index_col=0)
    assert R.loc["H1", "H2"] > 0.9 and abs(R.loc["H1", "L1"]) < 0.3
    s = pd.read_csv(tmp_path / "out" / "k3" / "group_unique_summary.tsv", sep="\t")
    assert (s.unique_tested_in_every_other_group <= s.unique).all() and (s.unique <= s.hits).all()
    assert (tmp_path / "out" / "k3" / "venn_unique_hits.pdf").exists()
