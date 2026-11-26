"""altanalyze3.components.snp.allelic_tests: the array statistics must equal statsmodels StratifiedTable,
BH must equal statsmodels multipletests, and run_family must equal a per-test StratifiedTable over the same
strata (the PILRA rs1859788 `stratified` function)."""
import numpy as np
import pandas as pd
import pytest
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests

from altanalyze3.components.snp import allelic_tests as at


def _pilra_stratified(tables, label):
    """Copy of rs1859788_association.stratified (PILRA_rs1859788_20261119/code), the reference."""
    inf = [t for t in tables if t[0].sum() > 0 and t[1].sum() > 0 and t[:, 0].sum() > 0 and t[:, 1].sum() > 0]
    out = dict(strata_total=len(tables), strata_informative=len(inf), reads=int(sum(t.sum() for t in tables)))
    if not inf:
        return dict(out, MH_OR=np.nan, CI_low=np.nan, CI_high=np.nan, CMH_p=np.nan, BreslowDay_p=np.nan)
    st = StratifiedTable([t.astype(float) for t in inf])
    lo, hi = st.oddsratio_pooled_confint()
    bd = st.test_equal_odds().pvalue if (len(inf) > 1 and not label.startswith("cell")) else np.nan
    return dict(out, MH_OR=st.oddsratio_pooled, CI_low=lo, CI_high=hi,
                CMH_p=st.test_null_odds(correction=True).pvalue, BreslowDay_p=bd)


def test_stratified_arrays_equal_statsmodels():
    rng = np.random.default_rng(7)
    rows, ref = [], []
    for t in range(300):
        k = rng.integers(1, 9)
        tabs = []
        for _ in range(k):
            tab = rng.integers(0, 30, size=(2, 2))
            if rng.random() < 0.15:
                tab[rng.integers(0, 2), :] = 0          # uninformative stratum
            tabs.append(tab)
            rows.append((t, *tab.ravel()))
        ref.append(_pilra_stratified(tabs, "donor"))
    df = pd.DataFrame(rows, columns=["test", "a", "b", "c", "d"])
    got = at.stratified_arrays(df.a, df.b, df.c, df.d, df.test, 300)
    for i, r in enumerate(ref):
        assert got["strata_informative"][i] == r["strata_informative"]
        for k in ("MH_OR", "CI_low", "CI_high"):
            if np.isnan(r[k]) or np.isinf(r[k]):
                assert np.isnan(got[k][i]) or np.isinf(got[k][i]) or np.isnan(r[k])
            else:
                assert got[k][i] == pytest.approx(r[k], rel=1e-10)
        for k in ("CMH_p", "BreslowDay_p"):
            if np.isnan(r[k]):
                assert np.isnan(got[k][i])
            elif r[k] > 1e-12:
                assert got[k][i] == pytest.approx(r[k], rel=1e-8, abs=1e-14)


def test_bh_equals_multipletests():
    rng = np.random.default_rng(3)
    p = rng.random(500) ** 3
    p[::37] = np.nan
    q = at.bh(p)
    ok = ~np.isnan(p)
    assert np.allclose(q[ok], multipletests(p[ok], method="fdr_bh")[1], rtol=1e-12)
    assert np.isnan(q[~ok]).all()


def test_run_family_equals_per_test_stratified_table():
    rng = np.random.default_rng(11)
    recs = []
    for site in ("chr7:10:A:G", "chr7:20:C:T"):
        for donor in ("D1", "D2", "D3"):
            for cell in range(25):
                state = "S{}".format(cell % 3)
                for _ in range(rng.integers(1, 6)):
                    iso = rng.choice(["I1", "I2", "I3"], p=[0.5, 0.3, 0.2])
                    A = int(rng.random() < (0.65 if (site.endswith("G") and iso == "I1") else 0.5))
                    recs.append((site, donor, "{}:{}".format(donor, cell), state, "G1", iso, A))
    df = pd.DataFrame(recs, columns=["site_id", "donor", "cell", "cell_state", "gene", "final", "A"])
    for c in ("site_id", "donor", "cell", "cell_state", "gene", "final"):
        df[c] = df[c].astype("category")
    for strata, label in ((["donor"], "donor"), (["donor", "cell_state"], "donor_x_cell_state"),
                          (["cell"], "cell__biallelic_cells_only")):
        got = at.run_family(df, "final", strata, label)
        for r in got.itertuples():
            sub = df[(df.site_id == r.site_id) & (df.gene == r.gene)]
            tabs = []
            for _, g in sub.groupby(strata, observed=True):
                x = (g.final == r.final).values
                A = g.A.values == 1
                tabs.append(np.array([[np.sum(A & x), np.sum(A & ~x)], [np.sum(~A & x), np.sum(~A & ~x)]]))
            ref = _pilra_stratified(tabs, label)
            assert (r.strata_total, r.strata_informative, r.reads) == (ref["strata_total"], ref["strata_informative"],
                                                                       ref["reads"])
            assert r.MH_OR == pytest.approx(ref["MH_OR"], rel=1e-10)
            assert r.CMH_p == pytest.approx(ref["CMH_p"], rel=1e-8)
            if not np.isnan(ref["BreslowDay_p"]):
                assert r.BreslowDay_p == pytest.approx(ref["BreslowDay_p"], rel=1e-8)


def test_covariate_map():
    assert at.covariate_of("young") == "young" and at.covariate_of("aged") == "aged"
    assert at.covariate_of("MDS-pre-SRSF2") == "MDS" and at.covariate_of("MDS-post-SRSF2-MISTRG") == "xenograft"
    assert at.covariate_of("AML-NPM1") == "AML" and at.covariate_of("pediatric-AML-TET2-pre") == "AML"
    assert at.covariate_of("MDS-pre-SRSF2-iPSC") == "iPSC" and at.covariate_of("Anemia") == "anemia"


def test_unassigned_rows_dropped_but_counted():
    import os, tempfile
    rng = np.random.default_rng(2)
    recs = []
    for donor in ("D1", "D2"):
        for cell in range(30):
            for _ in range(3):
                fin = rng.choice(["I1", "I2", "unassigned"], p=[0.5, 0.4, 0.1])
                recs.append(("chr1:5:A:G", donor, donor + ":" + str(cell), "S1", "G", fin, int(rng.random() < 0.5),
                             donor, donor + "|young", "young", "{}{}m{}".format(donor, cell, _)))
    df = pd.DataFrame(recs, columns=["site_id", "donor", "cell", "cell_state", "gene", "final", "A", "uid", "unit",
                                     "covariate", "molecule_id"])
    for c in ("site_id", "donor", "cell", "cell_state", "gene", "final", "uid", "unit", "covariate"):
        df[c] = df[c].astype("category")
    log = {}
    with tempfile.TemporaryDirectory() as d:
        at.isoform_suite(df, at.biallelic_subset(df), d, "", log)
        iso = pd.read_csv(os.path.join(d, "isoform_tests.tsv.gz"), sep="\t")
    assert "unassigned" not in set(iso.final) and set(iso.final) == {"I1", "I2"}
    assert log["unassigned_feature_rows_dropped"][":donor"] == 1
    # the unassigned reads remain in the I1 test's 'other' counts
    r = iso[(iso.final == "I1") & (iso.strata == "donor")].iloc[0]
    assert r.A_other + r.G_other == int((df.final != "I1").sum())


def test_combine_tests_recomputes_bh_and_refuses_duplicates(tmp_path):
    import os
    from altanalyze3.components.snp import combine_tests
    rng = np.random.default_rng(1)
    for chunk, start in (("chr7_1_100", 1), ("chr7_101_200", 101)):
        d = tmp_path / "stats" / chunk
        d.mkdir(parents=True)
        n = 40
        pd.DataFrame(dict(site_id=["chr7:{}:A:G".format(start + i) for i in range(n)], gene="G", final="I1",
                          strata=["donor"] * (n // 2) + ["cell__biallelic_cells_only"] * (n // 2),
                          CMH_p=rng.random(n) ** 2)).to_csv(d / "isoform_tests.tsv.gz", sep="\t", index=False)
    log = combine_tests.combine(str(tmp_path / "stats"), str(tmp_path / "combined"))
    out = pd.read_csv(tmp_path / "combined" / "isoform_tests.tsv.gz", sep="\t")
    assert len(out) == 80 and log["tables"]["isoform_tests"]["chunks"] == 2
    for lab, g in out.groupby("strata"):
        assert np.allclose(g.BH_q_combined.values, multipletests(g.CMH_p.values, method="fdr_bh")[1])
    # an overlapping chunk duplicates tests and must stop the run
    dup = tmp_path / "stats" / "chr7_50_150"
    dup.mkdir()
    pd.read_csv(tmp_path / "stats" / "chr7_1_100" / "isoform_tests.tsv.gz", sep="\t").to_csv(
        dup / "isoform_tests.tsv.gz", sep="\t", index=False)
    with pytest.raises(SystemExit):
        combine_tests.combine(str(tmp_path / "stats"), str(tmp_path / "combined2"))


def test_donor_genotypes_keeps_tier_per_donor(tmp_path):
    """The same site carries the tier of each donor's own panel (bug fixed 2026-11-25: the tier was taken from
    the first donor listed)."""
    import os
    cols = ["site_id", "chrom", "pos", "ref", "alt", "vtype", "tier", "source"] + \
           ["all_" + c for c in ("REF", "ALT", "OTHER", "SPLICE", "DEL")] + \
           ["primary_" + c for c in ("REF", "ALT", "OTHER", "SPLICE", "DEL")]
    lib = pd.DataFrame(dict(uid=["L1", "L2"], donor=["D1", "D2"]))
    for uid, tier in (("L1", "high"), ("L2", "low")):
        d = tmp_path / uid
        d.mkdir()
        pd.DataFrame([["chr1:5:A:G", "chr1", 5, "A", "G", "SNV", tier, "panel", 10, 9, 0, 0, 0, 10, 9, 0, 0, 0]],
                     columns=cols).to_csv(d / "chr1.site_counts.tsv.gz", sep="\t", index=False)
    dg = at.donor_genotypes(str(tmp_path), "chr1", lib, 0.01)
    assert dict(zip(dg.donor, dg.tier)) == {"D1": "high", "D2": "low"}
