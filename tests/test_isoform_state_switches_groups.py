"""Group mode of long_read.isoform_state_switches: two donor groups compared within each cell state."""
import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse, stats

from altanalyze3.components.long_read import isoform_state_switches as iss

YOUNG = ["Y1", "Y2", "Y3", "Y4"]
AGED = ["O1", "O2", "O3", "O4", "O5"]
STATES = ["S1", "S2", "S3"]
FILL = 3_000_000


def test_two_sample_t_by_hand_and_scipy():
    rng = np.random.default_rng(11)
    A = rng.normal(0.3, 0.1, size=(5, 40))
    B = rng.normal(0.25, 0.1, size=(4, 40))
    nA, nB = np.full(40, 5.0), np.full(40, 4.0)
    delta = A.mean(0) - B.mean(0)
    s2 = (4 * A.var(0, ddof=1) + 3 * B.var(0, ddof=1)) / 7
    t, dft, p, post = iss.moderated_two_sample_t(delta, s2, nA, nB, np.nan, np.nan)
    ref = stats.ttest_ind(A, B, axis=0, equal_var=True)
    assert np.allclose(t, ref.statistic) and np.allclose(p, ref.pvalue) and np.allclose(dft, 7)
    s0, d0 = 0.012, 6.0
    t, dft, p, post = iss.moderated_two_sample_t(delta, s2, nA, nB, s0, d0)
    exp_post = (7 * s2 + d0 * s0) / (7 + d0)
    exp_t = delta / np.sqrt(exp_post * (1 / 5 + 1 / 4))
    assert np.allclose(post, exp_post) and np.allclose(t, exp_t) and np.allclose(dft, 13)
    assert np.allclose(p, 2 * stats.t.sf(np.abs(exp_t), 13))
    t0, _, p0, _ = iss.moderated_two_sample_t(np.array([0.2]), np.array([0.0]), np.array([3.0]), np.array([3.0]),
                                              np.nan, np.nan)
    assert np.isnan(t0[0]) and np.isnan(p0[0])


def _cohort(tmp_path, reverse_fold_state="S2"):
    rng = np.random.default_rng(5)
    genes = {"GSW": ["a", "b", "c"], "GEV": ["e1", "e2"], "GLOW": ["a", "b"], "GFILL": ["f"]}
    for g in range(300):
        genes[f"GN{g:03d}"] = ["x", "y"] if g % 2 else ["x", "y", "z"]
    feats = [f"{g}:{i}" for g, isos in genes.items() for i in isos]
    fi = {f: k for k, f in enumerate(feats)}
    null_mu = {g: rng.dirichlet(np.full(len(isos), 3.0)) for g, isos in genes.items() if g.startswith("GN")}
    rows, names, truth = [], [], {}
    for uid in YOUNG + AGED:
        aged = uid in AGED
        donor_mu = {g: rng.dirichlet(150 * mu) for g, mu in null_mu.items()}
        for s in STATES:
            if s == "S3" and uid in ("Y3", "Y4"):
                continue
            x = np.zeros(len(feats))
            x[fi["GFILL:f"]] = FILL
            sw = np.array([0.7, 0.2, 0.1]) if (aged and s == "S1") else np.array([0.2, 0.7, 0.1])
            x[[fi[f"GSW:{i}"] for i in "abc"]] = rng.multinomial(int(rng.integers(300, 700)),
                                                                 rng.dirichlet(300 * sw))
            ev = np.array([0.75, 0.25]) if (aged and s == "S2") else np.array([0.25, 0.75])
            x[[fi["GEV:e1"], fi["GEV:e2"]]] = rng.multinomial(int(rng.integers(300, 700)), rng.dirichlet(300 * ev))
            x[[fi["GLOW:a"], fi["GLOW:b"]]] = rng.multinomial(3, [0.9, 0.1] if aged else [0.1, 0.9])
            for g, mu in donor_mu.items():
                cols = [fi[f"{g}:{i}"] for i in genes[g]]
                x[cols] = rng.multinomial(int(rng.integers(300, 700)), rng.dirichlet(400 * mu))
            rows.append(x)
            names.append(f"{s}.{uid}-isoform")
            truth[(s, uid)] = x
    X = sparse.csr_matrix(np.vstack(rows).astype(np.float32))
    counts = tmp_path / "counts.h5ad"
    ad.AnnData(X=X, obs=pd.DataFrame(index=names), var=pd.DataFrame(index=feats)).write_h5ad(counts)
    meta = pd.DataFrame({"uid": YOUNG + AGED + ["Y1"], "bam": "x.bam",
                         "library": [f"{u}_lib" for u in YOUNG + AGED + ["Y1"]], "reverse": "TRUE",
                         "groups": ["young"] * 4 + ["aged"] * 5 + ["young"]})
    groups = tmp_path / "groups.txt"
    meta.to_csv(groups, sep="\t", index=False)

    diff = tmp_path / "diff-covariate-isoform"
    diff.mkdir()
    for s in STATES:
        aged_minus_young = np.zeros(len(feats))
        q = np.full(len(feats), 0.5)
        p = np.full(len(feats), 0.5)
        if s == "S1":
            for f, v in (("GSW:a", 1.2), ("GSW:b", -1.2), ("GSW:c", 0.0), ("GLOW:a", 1.0), ("GLOW:b", -1.0)):
                aged_minus_young[fi[f]], p[fi[f]], q[fi[f]] = v, 0.001, 0.6
        for g in genes:
            if g.startswith("GN") and s == "S1":
                for i in genes[g]:
                    aged_minus_young[fi[f"{g}:{i}"]], p[fi[f"{g}:{i}"]] = 0.7, 0.01
        df = pd.DataFrame({"Gene_Symbol": [f.split(":")[0].lower() for f in feats], "Junctions": "E1.1",
                           "Protein Length": "100/100", "Predicted NMD": "Not-NMD", "Feature": feats,
                           "N_aged": 5, "N_young": 4, "Mean_aged": 1.0 + aged_minus_young, "Mean_young": 1.0})
        if s == reverse_fold_state:
            df["DiffMeans_young_m_aged"] = -aged_minus_young
        else:
            df["DiffMeans_aged_m_young"] = aged_minus_young
        df["ebayesStat"], df["ebayesSign"], df["ebayesPval"], df["ebayesAdjPval"] = 1.0, 1.0, p, q
        df.to_csv(diff / f"young-aged-{s}_stats-annotated.txt", sep="\t", index=False)
    (diff / "young-HSC-Others-HSC_stats-annotated.txt").write_text("ignored\tfile\n")

    events = tmp_path / "dPSI-covariate-events"
    events.mkdir()
    for s in STATES:
        ev = pd.DataFrame({"Gene_Symbol": ["gev", "gev", "gsw", "gn001"],
                           "Isoform_1|Length": ["e1|120", "e1|120", "a|100", "x|90"],
                           "Isoform_2|Length": ["e2|80(NMD)", "", "b|100", "y|90"],
                           "Feature": ["GEV:E1.1-E2.1|GEV:E1.1-E3.1", "GEV:E1.1-E4.1|GEV:E1.1-E5.1",
                                       "GSW:E2.1-E3.1|GSW:E2.1-E4.1", "GN001:E1.1-E2.1|GN001:E1.1-E3.1"],
                           "Mean_aged": [0.8, 0.8, 0.55, 0.5], "Mean_young": [0.3, 0.3, 0.5, 0.2],
                           "ebayesPval": [0.001, 0.001, 0.001, 0.2]})
        if s != "S2":
            ev["ebayesPval"] = 0.5
        ev["DiffMeans_aged_m_young"] = ev.Mean_aged - ev.Mean_young
        ev.to_csv(events / f"young-aged-{s}_stats-annotated.txt", sep="\t", index=False)

    pc = pd.DataFrame({"gene": [f.split(":")[0] for f in feats], "isoform": [f.split(":")[1] for f in feats],
                       "known_novel": "known", "class": "protein_coding", "protein_length": "100",
                       "nmd_status": "Not-NMD", "longest_gene_protein": "100", "prediction_rows": "1"})
    pc.loc[pc.gene.eq("GSW") & pc.isoform.eq("c"), ["class", "nmd_status"]] = ["likely_NMD", "NMD"]
    protein = tmp_path / "protein_classes.tsv"
    pc.to_csv(protein, sep="\t", index=False)
    order = tmp_path / "order.txt"
    order.write_text("S3\nS1\nS2\n")
    return dict(counts=counts, groups=groups, diff=diff, events=events, protein=protein, order=order, truth=truth,
                fi=fi, genes=genes)


def _args(c, out, *extra):
    return ["--counts", str(c["counts"]), "--groups", str(c["groups"]), "--compare_groups", "aged,young",
            "--diff_dir", str(c["diff"]), "--protein_classes", str(c["protein"]), "--outdir", str(out),
            "--cluster_order", str(c["order"]), "--min_donors", "3", "--whitelist_rawp", "0.05"] + list(extra)


def test_covariate_whitelist_orientation_and_rawp(tmp_path):
    c = _cohort(tmp_path)
    wl, manifest = iss.load_covariate_whitelist(str(c["diff"]), "aged", "young", 1.5, 0.05, rawp=0.05)
    got = set(wl.feature)
    assert {"GSW:a", "GSW:b", "GLOW:a", "GLOW:b"} <= got and "GSW:c" not in got
    w = wl.set_index("feature")
    assert w.loc["GSW:a", "de_up_states"] == "S1" and w.loc["GSW:b", "de_down_states"] == "S1"
    assert [m["state"] for m in manifest] == STATES
    assert [m["fold_column"] for m in manifest] == ["DiffMeans_aged_m_young", "DiffMeans_young_m_aged",
                                                    "DiffMeans_aged_m_young"]
    wl_fdr, _ = iss.load_covariate_whitelist(str(c["diff"]), "aged", "young", 1.5, 0.05)
    assert wl_fdr.empty
    d = pd.read_csv(c["diff"] / "young-aged-S1_stats-annotated.txt", sep="\t")
    d["Mean_aged"] = 1.0 - d.DiffMeans_aged_m_young
    d.to_csv(c["diff"] / "young-aged-S1_stats-annotated.txt", sep="\t", index=False)
    with pytest.raises(ValueError, match="Mean columns"):
        iss.load_covariate_whitelist(str(c["diff"]), "aged", "young", 1.5, 0.05, rawp=0.05)


def test_event_whitelist_takes_both_sides(tmp_path):
    c = _cohort(tmp_path)
    ev, manifest = iss.load_event_whitelist(str(c["events"]), "aged", "young", 0.2, 0.05)
    assert set(ev.feature) == {"GEV:e1", "GEV:e2"}
    assert ev.set_index("feature").loc["GEV:e1", "event_states"] == "S2"
    s2 = [m for m in manifest if m["state"] == "S2"][0]
    assert s2["significant_events"] == 2 and s2["event_sides_without_isoform"] == 1


def test_event_isoform_keyed_to_registered_gene(tmp_path):
    ev_dir = tmp_path / "events"
    ev_dir.mkdir()
    pd.DataFrame({"Gene_Symbol": ["ga", "ga"], "Isoform_1|Length": ["iso9|0", "iso7|50"],
                  "Isoform_2|Length": ["ENST1|88", "iso8|40"],
                  "Feature": ["GA:E9.1-GB:I7.2_100|GA:E9.1-I9.1", "GA:E1.1-E2.1|GA:E1.1-E3.1"],
                  "Mean_aged": [0.8, 0.9], "Mean_young": [0.2, 0.1], "DiffMeans_aged_m_young": [0.6, 0.8],
                  "ebayesPval": [0.001, 0.001]}).to_csv(ev_dir / "young-aged-S1_stats-annotated.txt", sep="\t",
                                                         index=False)
    naive, _ = iss.load_event_whitelist(str(ev_dir), "aged", "young", 0.2, 0.05)
    assert set(naive.feature) == {"GA:iso9", "GA:ENST1", "GA:iso7", "GA:iso8"}
    reg = {"iso9": {"GB"}, "ENST1": {"GA"}, "iso7": {"GA", "GC"}, "iso8": {"GA"}}
    ev, _ = iss.load_event_whitelist(str(ev_dir), "aged", "young", 0.2, 0.05, iso_to_genes=reg)
    assert set(ev.feature) == {"GB:iso9", "GA:ENST1", "GA:iso7", "GA:iso8"}
    with pytest.raises(ValueError, match="cannot place"):
        iss.load_event_whitelist(str(ev_dir), "aged", "young", 0.2, 0.05, iso_to_genes={"iso9": {"GB"}})
    with pytest.raises(ValueError, match="several genes"):
        iss.load_event_whitelist(str(ev_dir), "aged", "young", 0.2, 0.05,
                                 iso_to_genes={**reg, "iso7": {"GC", "GD"}})


def test_end_to_end_group_switches(tmp_path):
    c = _cohort(tmp_path)
    out = tmp_path / "out"
    iss.main(_args(c, out, "--event_dir", str(c["events"]), "--required_switch", "GSW:a,GSW:b",
                   "--required_switch", "GEV:e1,GEV:e2"))
    for f in ["isoform_state_switches.tsv", "isoform_state_usage_tests.tsv.gz", "isoform_state_switch_priors.tsv",
              "isoform_state_switch_whitelist.tsv", "isoform_state_switch_candidates.tsv",
              "isoform_state_switch_summary.tsv", "positive_controls.tsv", "logs/parameters.json", "logs/run_log.txt"]:
        assert (out / f).exists(), f
    assert not (out / "RUN_INVALID_missing_positive_controls.txt").exists()
    params = json.load(open(out / "logs" / "parameters.json"))
    assert params["mode"] == "compare_groups" and params["donors_A"] == AGED and params["donors_B"] == YOUNG
    assert params["states"] == ["S3", "S1", "S2"] and params["donors_per_state_B"]["S3"] == 2
    assert params["controls_missing"] == []

    wl = pd.read_csv(out / "isoform_state_switch_whitelist.tsv", sep="\t", keep_default_na=False).set_index("feature")
    assert wl.loc["GEV:e1", "whitelist_source"] == "event" and wl.loc["GSW:a", "whitelist_source"] == "isoform"

    sw = pd.read_csv(out / "isoform_state_switches.tsv", sep="\t")
    got = {(r.state, r.isoform_higher_in_A, r.isoform_higher_in_B) for r in sw.itertuples()}
    assert got == {("S1", "GSW:a", "GSW:b"), ("S2", "GEV:e1", "GEV:e2")}, got

    rows = pd.read_csv(out / "isoform_state_usage_tests.tsv.gz", sep="\t")
    s3 = rows[rows.state.eq("S3")]
    assert len(s3) > 0 and (~s3.tested).all() and s3.p.isna().all() and s3.n_donors_B.le(2).all()
    assert not rows.gene_id.eq("GLOW").any()
    r = rows[(rows.state == "S1") & (rows.feature == "GSW:a")].iloc[0]
    fi, truth = c["fi"], c["truth"]
    cols = [fi[f"GSW:{i}"] for i in "abc"]
    ua = [truth[("S1", u)][fi["GSW:a"]] / truth[("S1", u)][cols].sum() for u in AGED]
    uy = [truth[("S1", u)][fi["GSW:a"]] / truth[("S1", u)][cols].sum() for u in YOUNG]
    assert r.n_donors_A == 5 and r.n_donors_B == 4
    assert np.isclose(r.delta_usage, np.mean(ua) - np.mean(uy), atol=1e-6)
    pooled = (4 * np.var(ua, ddof=1) + 3 * np.var(uy, ddof=1)) / 7
    assert np.isclose(r.sd_pooled, np.sqrt(pooled), atol=1e-6)
    assert r.reads_A == sum(truth[("S1", u)][fi["GSW:a"]] for u in AGED)

    pri = pd.read_csv(out / "isoform_state_switch_priors.tsv", sep="\t").set_index("state")
    sub = rows[(rows.state == "S1") & rows.tested]
    t, _, p, _ = iss.moderated_two_sample_t(sub.delta_usage.to_numpy(), sub.sd_pooled.to_numpy() ** 2,
                                            sub.n_donors_A.to_numpy(float), sub.n_donors_B.to_numpy(float),
                                            pri.loc["S1", "s0_squared"], pri.loc["S1", "d0"])
    assert np.allclose(t, sub.t_moderated, rtol=1e-4) and np.allclose(p, sub.p, rtol=1e-4, atol=1e-300)

    null = rows[rows.tested & rows.gene_id.str.startswith("GN")]
    assert len(null) > 500
    assert 0.01 < (null.p < 0.05).mean() < 0.10
    assert stats.kstest(null.drop_duplicates(["state", "gene_id"]).p, "uniform").pvalue > 1e-3
    assert not sw.gene_id.str.startswith("GN").any()


def test_missing_control_fails_the_run(tmp_path):
    c = _cohort(tmp_path)
    out = tmp_path / "out"
    with pytest.raises(RuntimeError, match="positive control"):
        iss.main(_args(c, out, "--required_switch", "GSW:b,GSW:a"))
    assert (out / "RUN_INVALID_missing_positive_controls.txt").exists()
    ctl = pd.read_csv(out / "positive_controls.tsv", sep="\t")
    assert not ctl.called.any()
    assert json.load(open(out / "logs" / "parameters.json"))["controls_missing"] == ["GSW:b > GSW:a"]


def test_without_events_the_event_only_switch_is_absent(tmp_path):
    c = _cohort(tmp_path)
    out = tmp_path / "out"
    iss.main(_args(c, out))
    sw = pd.read_csv(out / "isoform_state_switches.tsv", sep="\t")
    assert {(r.state, r.isoform_higher_in_A, r.isoform_higher_in_B) for r in sw.itertuples()} == \
        {("S1", "GSW:a", "GSW:b")}


def test_paired_mode_rejects_group_only_flags(tmp_path):
    c = _cohort(tmp_path)
    with pytest.raises(SystemExit):
        iss.main(["--counts", str(c["counts"]), "--groups", str(c["groups"]), "--group", "young",
                  "--diff_dir", str(c["diff"]), "--protein_classes", str(c["protein"]),
                  "--outdir", str(tmp_path / "o"), "--whitelist_rawp", "0.05"])
