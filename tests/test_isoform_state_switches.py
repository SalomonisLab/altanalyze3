"""Donor-paired isoform switches between cell states (altanalyze3.components.long_read.isoform_state_switches).

The moderated statistics are checked against their definitions (limma fitFDist / squeezeVar / eBayes, Smyth 2004):
prior recovery from simulated variances, the trigamma inverse, the moderated t by hand and BH against statsmodels.
The end-to-end test builds a cohort with one planted switch, decoys for every gate and 300 null genes, then checks
the calls, the paired deltas by hand, the gene-CPM missing rule, group selection and null calibration.
"""
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, special, stats
from statsmodels.stats.multitest import fdrcorrection

from altanalyze3.components.long_read import isoform_state_switches as iss


def test_trigamma_inverse_roundtrip():
    x = np.logspace(-5, 6, 200)
    y = iss.trigamma_inverse(x)
    assert np.allclose(special.polygamma(1, y), x, rtol=1e-6)
    ext = iss.trigamma_inverse(np.array([1e-7, 1e8]))
    assert np.allclose(special.polygamma(1, ext), [1e-7, 1e8], rtol=1e-3)


def test_fit_f_dist_recovers_prior():
    rng = np.random.default_rng(1)
    G, d0, s0 = 50_000, 6.0, 0.02
    df = rng.integers(2, 9, G).astype(float)
    sigma2 = d0 * s0 / rng.chisquare(d0, G)
    s2 = sigma2 * rng.chisquare(df) / df
    est_s0, est_d0, info = iss.fit_f_dist(s2, df)
    assert abs(est_d0 - d0) / d0 < 0.1, est_d0
    assert abs(est_s0 - s0) / s0 < 0.03, est_s0
    assert info["prior_features"] == G and info["zero_variance"] == 0


def test_fit_f_dist_without_prior_variation():
    rng = np.random.default_rng(2)
    df = np.full(20_000, 4.0)
    s2 = 0.05 * rng.chisquare(df) / df
    s0, d0, _ = iss.fit_f_dist(s2, df)
    assert d0 > 100 and abs(s0 - 0.05) / 0.05 < 0.03


def test_moderated_t_by_hand():
    mean = np.array([0.30, -0.05, 0.0])
    s2 = np.array([0.01, 0.02, 0.0])
    n = np.array([9.0, 5.0, 4.0])
    s0, d0 = 0.015, 4.0
    t, dft, p, post = iss.moderated_one_sample_t(mean, s2, n, s0, d0)
    exp_post = ((n - 1) * s2 + d0 * s0) / (n - 1 + d0)
    exp_t = mean / np.sqrt(exp_post / n)
    exp_df = np.minimum(n - 1 + d0, (n - 1).sum())
    assert np.allclose(post, exp_post) and np.allclose(t, exp_t) and np.allclose(dft, exp_df)
    assert np.allclose(p, 2 * stats.t.sf(np.abs(exp_t), exp_df))
    t2, dft2, _, _ = iss.moderated_one_sample_t(mean, s2, n, np.nan, np.nan)
    assert np.allclose(dft2, n - 1) and np.isnan(t2[2])


def test_bh_matches_statsmodels():
    rng = np.random.default_rng(3)
    p = rng.uniform(size=5000) ** 3
    p[::97] = np.nan
    ok = np.isfinite(p)
    assert np.allclose(iss.bh_fdr(p)[ok], fdrcorrection(p[ok])[1])
    assert np.isnan(iss.bh_fdr(p)[~ok]).all()


# ---------------------------------------------------------------- synthetic cohort
DONORS = [f"D{i}" for i in range(1, 10)]
OTHER = ["A1", "A2"]
STATES = ["S1", "S2", "S3"]
FILL = 3_000_000


def _cohort(tmp_path):
    rng = np.random.default_rng(7)
    genes = {"GSW": ["a", "b", "c"], "GLOW": ["a", "b"], "GONE": ["a", "b"], "GBORDER": ["b1", "b2", "b3"],
             "GTRUNC": ["t1", "t2"], "GMINOR": ["m1", "m2"], "GFILL": ["f"]}
    for g in range(300):
        genes[f"GN{g:03d}"] = ["x", "y"] if g % 2 else ["x", "y", "z"]
    feats = [f"{g}:{i}" for g, isos in genes.items() for i in isos]
    fi = {f: k for k, f in enumerate(feats)}
    null_mu = {g: rng.dirichlet(np.full(len(isos), 3.0)) for g, isos in genes.items() if g.startswith("GN")}
    rows, names, truth = [], [], {}
    for uid in DONORS + OTHER:
        donor_mu = {g: rng.dirichlet(150 * mu) for g, mu in null_mu.items()}
        for s in STATES + (["S4"] if uid == "D1" else []):
            x = np.zeros(len(feats))
            x[fi["GFILL:f"]] = FILL
            sw = np.array([0.7, 0.2, 0.1]) if s in ("S1", "S3") else np.array([0.2, 0.7, 0.1])
            if uid in OTHER:
                sw = sw[[1, 0, 2]]
            depth = 2 if (uid == "D9" and s == "S3") else int(rng.integers(300, 700))
            x[[fi[f"GSW:{i}"] for i in "abc"]] = rng.multinomial(depth, rng.dirichlet(300 * sw))
            x[[fi["GLOW:a"], fi["GLOW:b"]]] = rng.multinomial(15, [0.9, 0.1] if s == "S1" else [0.1, 0.9])
            x[[fi["GONE:a"], fi["GONE:b"]]] = rng.multinomial(400, [0.5, 0.5])
            x[[fi[f"GBORDER:{i}"] for i in ("b1", "b2", "b3")]] = rng.multinomial(400, [0.3, 0.3, 0.4])
            x[[fi["GTRUNC:t1"], fi["GTRUNC:t2"]]] = rng.multinomial(400, [0.5, 0.5])
            mi = np.array([0.95, 0.05]) if s in ("S1", "S3") else np.array([0.6, 0.4])
            x[[fi["GMINOR:m1"], fi["GMINOR:m2"]]] = rng.multinomial(40, rng.dirichlet(300 * mi))
            for g, mu in donor_mu.items():
                cols = [fi[f"{g}:{i}"] for i in genes[g]]
                x[cols] = rng.multinomial(int(rng.integers(300, 700)), rng.dirichlet(400 * mu))
            rows.append(x)
            names.append(f"{s}.{uid}-isoform")
            truth[(s, uid)] = x
    X = sparse.csr_matrix(np.vstack(rows).astype(np.float32))
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=names), var=pd.DataFrame(index=feats))
    counts = tmp_path / "counts.h5ad"
    a.write_h5ad(counts)

    meta = pd.DataFrame({"uid": DONORS + OTHER, "bam": "x.bam", "library": [f"{u}_lib" for u in DONORS + OTHER],
                         "reverse": "TRUE", "groups": ["healthy"] * len(DONORS) + ["AML"] * len(OTHER)})
    groups = tmp_path / "groups.txt"
    meta.to_csv(groups, sep="\t", index=False)

    de = {f: (0.0, 0.5) for f in feats}
    de.update({"GSW:a": (1.0, 0.001), "GSW:b": (-1.0, 0.001), "GSW:c": (1.0, 0.001), "GLOW:a": (1.0, 0.001),
               "GLOW:b": (-1.0, 0.001), "GONE:a": (1.0, 0.001), "GBORDER:b1": (0.5, 0.001),
               "GBORDER:b2": (1.0, 0.06), "GBORDER:b3": (-0.7, 0.01), "GONE:b": (1.0, np.nan),
               "GTRUNC:t1": (0.7, 0.01), "GTRUNC:t2": (0.7, 0.01), "GMINOR:m1": (0.7, 0.01),
               "GMINOR:m2": (0.7, 0.01)})
    for g in genes:
        if g.startswith("GN"):
            for i in genes[g]:
                de[f"{g}:{i}"] = (0.7, 0.01)
    diff = tmp_path / "diff-cluster-isoform"
    diff.mkdir()
    for s in STATES:
        fc = {f: (v if s == "S1" else (0.0, 0.5)) for f, v in de.items()}
        if s == "S3":
            fc["GBORDER:b3"] = (0.8, 0.01)
        state_minus_others = np.array([fc[f][0] for f in feats])
        df = pd.DataFrame({"Gene_Symbol": [f.split(":")[0].lower() for f in feats], "Junctions": "E1.1",
                           "Protein Length": "100/100", "Predicted NMD": "Not-NMD", "Feature": feats,
                           f"N_{s}": 9, "N_Others": 18, f"Mean_{s}": 1.0 + state_minus_others, "Mean_Others": 1.0})
        if s == "S3":
            df[f"DiffMeans_Others_m_{s}"] = -state_minus_others
        else:
            df[f"DiffMeans_{s}_m_Others"] = state_minus_others
        df["ebayesStat"], df["ebayesSign"], df["ebayesPval"] = 1.0, 1.0, 0.01
        df["ebayesAdjPval"] = [fc[f][1] for f in feats]
        df.to_csv(diff / f"healthy-{s}-Others-{s}_stats-annotated.txt", sep="\t", index=False)

    pc = pd.DataFrame({"gene": [f.split(":")[0] for f in feats], "isoform": [f.split(":")[1] for f in feats],
                       "known_novel": "known", "class": "protein_coding", "protein_length": "100",
                       "nmd_status": "Not-NMD", "longest_gene_protein": "100", "prediction_rows": "1"})
    pc.loc[pc.gene.eq("GSW") & pc.isoform.eq("c"), ["class", "nmd_status"]] = ["likely_NMD", "NMD"]
    pc.loc[pc.gene.eq("GTRUNC") & pc.isoform.eq("t2"), ["class", "protein_length"]] = \
        ["truncated_lt_25pct_longest", "20"]
    protein = tmp_path / "protein_classes.tsv"
    pc.to_csv(protein, sep="\t", index=False)
    order = tmp_path / "order.txt"
    order.write_text("S1\nS2\nS3\nS4\n")
    return dict(counts=counts, groups=groups, diff=diff, protein=protein, order=order, truth=truth, fi=fi,
                genes=genes)


def test_whitelist_thresholds(tmp_path):
    c = _cohort(tmp_path)
    wl, manifest = iss.load_differential_whitelist(str(c["diff"]), 1.5, 0.05)
    got = set(wl.feature)
    assert {"GSW:a", "GSW:b", "GSW:c", "GLOW:a", "GLOW:b", "GONE:a", "GBORDER:b3"} <= got
    assert not ({"GBORDER:b1", "GBORDER:b2", "GONE:b", "GFILL:f"} & got)
    w = wl.set_index("feature")
    assert w.loc["GSW:b", "de_down_states"] == "S1"
    assert w.loc["GBORDER:b3", "de_up_states"] == "S3" and w.loc["GBORDER:b3", "de_down_states"] == "S1"
    assert [m["state"] for m in manifest] == STATES
    assert [m["negated"] for m in manifest] == [False, False, True]

    bad = tmp_path / "bad"
    bad.mkdir()
    d = pd.read_csv(c["diff"] / "healthy-S1-Others-S1_stats-annotated.txt", sep="\t")
    d["Mean_S1"] = 1.0 - d["DiffMeans_S1_m_Others"]
    d.to_csv(bad / "healthy-S1-Others-S1_stats-annotated.txt", sep="\t", index=False)
    try:
        iss.load_differential_whitelist(str(bad), 1.5, 0.05)
        raise AssertionError("an orientation that contradicts the Mean columns must raise")
    except ValueError as e:
        assert "Mean columns" in str(e)


def test_end_to_end_switch_calls(tmp_path):
    c = _cohort(tmp_path)
    out = tmp_path / "out"
    iss.main(["--counts", str(c["counts"]), "--groups", str(c["groups"]), "--group", "healthy",
              "--diff_dir", str(c["diff"]), "--protein_classes", str(c["protein"]), "--outdir", str(out),
              "--cluster_order", str(c["order"]), "--min_reads", "200"])
    for f in ["isoform_state_switches.tsv", "isoform_state_usage_tests.tsv.gz", "isoform_state_switch_priors.tsv",
              "isoform_state_switch_whitelist.tsv",
              "isoform_state_switch_candidates.tsv", "isoform_state_switch_summary.tsv", "logs/parameters.json",
              "logs/run_log.txt"]:
        assert (out / f).exists(), f
    params = json.load(open(out / "logs" / "parameters.json"))
    assert params["donors"] == DONORS and params["states"] == STATES + ["S4"]
    assert params["donors_per_state"]["S4"] == 1

    cand = pd.read_csv(out / "isoform_state_switch_candidates.tsv", sep="\t")
    assert "GSW:c" not in set(cand.feature) and "GONE:a" not in set(cand.feature)
    assert {"GTRUNC:t1", "GTRUNC:t2"} <= set(cand.feature)
    wlo = pd.read_csv(out / "isoform_state_switch_whitelist.tsv", sep="\t").set_index("feature")
    assert wlo.loc["GSW:c", "excluded_by"].startswith("not protein coding") and not wlo.loc["GSW:c", "candidate"]
    assert wlo.loc["GONE:a", "excluded_by"].startswith("only protein-coding") and wlo.loc["GSW:a", "candidate"]
    assert set(wlo.index[wlo.candidate]) == set(cand.feature)
    assert not cand.gene_id.isin(["GBORDER", "GFILL"]).any()

    sw = pd.read_csv(out / "isoform_state_switches.tsv", sep="\t")
    got = {(r.state_A, r.state_B, r.isoform_higher_in_A, r.isoform_higher_in_B) for r in sw.itertuples()}
    assert got == {("S1", "S2", "GSW:a", "GSW:b"), ("S2", "S3", "GSW:b", "GSW:a")}, got

    allrows = pd.read_csv(out / "isoform_state_usage_tests.tsv.gz", sep="\t")
    s4 = allrows[allrows.state_B.eq("S4")]
    assert len(s4) > 0 and s4.n_donors.eq(1).all() and (~s4.tested).all() and s4.p.isna().all()
    assert s4.fdr.isna().all() and s4.higher_in.isna().all()
    assert allrows.tested.eq(allrows.n_donors >= 2).all()
    tests = allrows[allrows.tested]
    assert not tests.gene_id.eq("GLOW").any()
    row = tests[(tests.state_A == "S1") & (tests.state_B == "S2") & (tests.feature == "GSW:a")].iloc[0]
    fi, truth = c["fi"], c["truth"]
    sw_cols = [fi[f"GSW:{i}"] for i in "abc"]
    d = [truth[("S1", u)][fi["GSW:a"]] / truth[("S1", u)][sw_cols].sum()
         - truth[("S2", u)][fi["GSW:a"]] / truth[("S2", u)][sw_cols].sum() for u in DONORS]
    assert row.n_donors == 9
    assert np.isclose(row.delta_usage, np.mean(d), atol=1e-6) and np.isclose(row.sd_delta, np.std(d, ddof=1),
                                                                             atol=1e-6)
    assert row.reads_A == sum(truth[("S1", u)][fi["GSW:a"]] for u in DONORS)
    row13 = tests[(tests.state_A == "S1") & (tests.state_B == "S3") & (tests.feature == "GSW:a")].iloc[0]
    assert row13.n_donors == 8

    pri = pd.read_csv(out / "isoform_state_switch_priors.tsv", sep="\t").set_index(["state_A", "state_B"])
    sub = tests[(tests.state_A == "S1") & (tests.state_B == "S2")]
    s0, d0 = pri.loc[("S1", "S2"), "s0_squared"], pri.loc[("S1", "S2"), "d0"]
    t, dft, p, _ = iss.moderated_one_sample_t(sub.delta_usage.to_numpy(), sub.sd_delta.to_numpy() ** 2,
                                              sub.n_donors.to_numpy(float), s0, d0)
    assert np.allclose(t, sub.t_moderated, rtol=1e-4) and np.allclose(p, sub.p, rtol=1e-4, atol=1e-300)

    null = tests[tests.gene_id.str.startswith("GN")]
    assert len(null) > 1500
    frac = (null.p < 0.05).mean()
    assert 0.02 < frac < 0.09, frac
    assert (null.fdr < 0.05).sum() <= 0.05 * len(null)
    assert stats.kstest(null.drop_duplicates(["state_A", "state_B", "gene_id"]).p, "uniform").pvalue > 1e-3
    assert not sw.gene_id.str.startswith("GN").any()


def test_gene_read_rule_keeps_low_read_partner(tmp_path):
    c = _cohort(tmp_path)
    calls = {}
    for rule in ("isoform", "gene"):
        out = tmp_path / rule
        iss.main(["--counts", str(c["counts"]), "--groups", str(c["groups"]), "--group", "healthy",
                  "--diff_dir", str(c["diff"]), "--protein_classes", str(c["protein"]), "--outdir", str(out),
                  "--cluster_order", str(c["order"]), "--min_reads", "200", "--read_rule", rule])
        sw = pd.read_csv(out / "isoform_state_switches.tsv", sep="\t")
        calls[rule] = {(r.state_A, r.state_B, r.isoform_higher_in_A, r.isoform_higher_in_B) for r in sw.itertuples()}
        rows = pd.read_csv(out / "isoform_state_usage_tests.tsv.gz", sep="\t")
        m2 = rows[(rows.feature == "GMINOR:m2") & (rows.state_A == "S1") & (rows.state_B == "S2")]
        if rule == "isoform":
            assert m2.empty or (m2.reads_B <= 200).all()
        else:
            assert len(m2) == 1 and m2.reads_B.iloc[0] <= 200 and m2.gene_reads_B.iloc[0] > 200
    minor = {("S1", "S2", "GMINOR:m1", "GMINOR:m2"), ("S2", "S3", "GMINOR:m2", "GMINOR:m1")}
    gsw = {("S1", "S2", "GSW:a", "GSW:b"), ("S2", "S3", "GSW:b", "GSW:a")}
    assert calls["isoform"] == gsw, calls["isoform"]
    assert calls["gene"] == gsw | minor, calls["gene"]


def test_all_isoforms_keeps_noncoding_partners(tmp_path):
    c = _cohort(tmp_path)
    pc = pd.read_csv(c["protein"], sep="\t", dtype=str)
    pc = pc[~(pc.gene.eq("GSW") & pc.isoform.eq("a"))]
    pc.loc[pc.gene.eq("GSW") & pc.isoform.eq("b"), ["class", "nmd_status"]] = ["likely_NMD", "NMD"]
    protein = tmp_path / "protein_noncoding.tsv"
    pc.to_csv(protein, sep="\t", index=False)
    calls, cands, wls = {}, {}, {}
    for mode, extra in (("coding", []), ("all", ["--all_isoforms"])):
        out = tmp_path / mode
        iss.main(["--counts", str(c["counts"]), "--groups", str(c["groups"]), "--group", "healthy",
                  "--diff_dir", str(c["diff"]), "--protein_classes", str(protein), "--outdir", str(out),
                  "--cluster_order", str(c["order"]), "--min_reads", "200"] + extra)
        sw = pd.read_csv(out / "isoform_state_switches.tsv", sep="\t")
        calls[mode] = {(r.state_A, r.state_B, r.isoform_higher_in_A, r.isoform_higher_in_B) for r in sw.itertuples()}
        cands[mode] = pd.read_csv(out / "isoform_state_switch_candidates.tsv", sep="\t")
        wls[mode] = pd.read_csv(out / "isoform_state_switch_whitelist.tsv", sep="\t").set_index("feature")
        assert json.load(open(out / "logs" / "parameters.json"))["coding_only"] == (mode == "coding")
    assert calls["coding"] == set(), calls["coding"]
    assert calls["all"] == {("S1", "S2", "GSW:a", "GSW:b"), ("S2", "S3", "GSW:b", "GSW:a")}, calls["all"]
    assert not cands["coding"].gene_id.eq("GSW").any()
    assert {"GSW:a", "GSW:b", "GSW:c"} <= set(cands["all"].feature)
    assert cands["all"].set_index("feature").protein_length.isna()["GSW:a"]
    wl_all = set(wls["all"].index)
    multi = {f for f in wl_all if sum(g.split(":")[0] == f.split(":")[0] for g in wl_all) >= 2}
    assert set(cands["all"].feature) == multi
    assert wls["all"].loc["GONE:a", "excluded_by"] == "only whitelisted isoform of its gene"
    assert wls["coding"].loc["GSW:b", "excluded_by"].startswith("not protein coding")
    assert (wls["all"].excluded_by.fillna("") != "not protein coding (length 0 or NMD)").all()


def test_read_scope_cohort_pools_every_group(tmp_path):
    c = _cohort(tmp_path)
    a = ad.read_h5ad(c["counts"])
    X = a.X.toarray()
    other = np.array([n.split(".", 1)[1].replace("-isoform", "") in OTHER for n in a.obs_names])
    cols = [c["fi"]["GMINOR:m1"], c["fi"]["GMINOR:m2"]]
    X[np.ix_(other, cols)] *= 10
    a.X = sparse.csr_matrix(X.astype(np.float32))
    counts = tmp_path / "counts_deep_other.h5ad"
    a.write_h5ad(counts)
    calls = {}
    for scope in ("group", "cohort"):
        out = tmp_path / scope
        iss.main(["--counts", str(counts), "--groups", str(c["groups"]), "--group", "healthy",
                  "--diff_dir", str(c["diff"]), "--protein_classes", str(c["protein"]), "--outdir", str(out),
                  "--cluster_order", str(c["order"]), "--min_reads", "200", "--read_scope", scope])
        sw = pd.read_csv(out / "isoform_state_switches.tsv", sep="\t")
        calls[scope] = {(r.state_A, r.state_B, r.isoform_higher_in_A, r.isoform_higher_in_B) for r in sw.itertuples()}
        rows = pd.read_csv(out / "isoform_state_usage_tests.tsv.gz", sep="\t")
        r = rows[(rows.feature == "GMINOR:m2") & (rows.state_A == "S1") & (rows.state_B == "S2")]
        if scope == "cohort":
            st = np.array([n.split(".", 1)[0] for n in a.obs_names])
            assert len(r) == 1 and r.reads_B.iloc[0] <= 200 < r.cohort_reads_B.iloc[0]
            assert r.cohort_reads_B.iloc[0] == X[st == "S2", c["fi"]["GMINOR:m2"]].sum()
            assert r.reads_B.iloc[0] == X[(st == "S2") & ~other, c["fi"]["GMINOR:m2"]].sum()
            assert (rows.cohort_reads_A >= rows.reads_A).all() and (rows.cohort_gene_reads_B >= rows.gene_reads_B).all()
        assert json.load(open(out / "logs" / "parameters.json"))["read_scope"] == scope
    gsw = {("S1", "S2", "GSW:a", "GSW:b"), ("S2", "S3", "GSW:b", "GSW:a")}
    minor = {("S1", "S2", "GMINOR:m1", "GMINOR:m2"), ("S2", "S3", "GMINOR:m2", "GMINOR:m1")}
    assert calls["group"] == gsw, calls["group"]
    assert calls["cohort"] == gsw | minor, calls["cohort"]
