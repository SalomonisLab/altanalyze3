#!/usr/bin/env python3
"""Within-donor allelic tests on one chromosome: allele-specific expression (ASE), allele x isoform usage
(isoform QTL) and allele x junction usage (splicing QTL), pooled across donors, per cell state and per
covariate group.

Inputs are the read_sites outputs of every library (``<reads>/<uid>/<chrom>.reads.tsv.gz`` and
``.site_counts.tsv.gz``), the donor map and the ISV run metadata.

1. Donor genotype per site: all-alignment REF and ALT counts summed over the donor's libraries, binomial
   genotype likelihood with ``--error_rate`` (read_alleles.genotype_call). Only 0/1 donors enter steps 2-5.
2. Records: primary alignments (not secondary, not supplementary) of heterozygous donors with a REF or ALT
   call at the site and a row in the molecule index (gene, cell state, final isoform).
3. Tests per (site, gene, feature), feature = a final isoform or an intron observed in the records.
   The 2x2 table of a stratum is [[REF & feature, REF & not], [ALT & feature, ALT & not]]; for an intron the
   table holds only records whose read spans both intron flanks. Strata families follow the PILRA
   rs1859788 protocol (rs1859788_association.py, 2026-11-19): donor; donor x cell state; donor, cells
   carrying both alleles; cell, cells carrying both alleles. Mantel-Haenszel odds ratio, Robins-Breslow-
   Greenland 95% interval, CMH chi-square with continuity correction, Breslow-Day (not for cell strata),
   per-donor Fisher exact. A table with an empty row or column carries no information and is not used.
4. Cell-state-specific: per (site, gene, feature, cell state), strata = donor. Cochran's Q across states of
   the per-state MH log odds ratios (inverse RGB variance weights).
5. Covariate-specific: library groups map to covariates (young, aged, MDS, AML, xenograft, iPSC, anemia).
   Per (site, gene, feature, covariate), strata = donor x cell state; per (donor, covariate) unit the MH log
   odds ratio over cell-state strata; Cochran's Q across covariates.
6. ASE: per (site, donor) and per (site, donor, cell state), two-sided binomial test of the REF fraction
   against 0.5 on primary alignments.

The array implementation of the stratified statistics reproduces statsmodels StratifiedTable (tested in
tests/test_snp_allelic_tests.py). P-values use chi2.sf; statsmodels uses 1 - chi2.cdf, which equals 0
below about 1e-16. BH correction runs within each test family over the tests that have a p-value.

Run: python -m altanalyze3.components.snp.allelic_tests --chrom chr7 --reads READS --donor_map donor_map.tsv \\
        --isv_metadata metadata.tsv --out STATS
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
import os
import sys
import time
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy import stats

try:
    from altanalyze3.components.bam.read_alleles import genotype_call
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bam"))
    from read_alleles import genotype_call

FAMILIES = ("donor", "donor_x_cell_state", "donor__biallelic_cells_only", "cell__biallelic_cells_only")
Z975 = stats.norm.ppf(0.975)


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def covariate_of(group):
    """Coarse covariate of an ISV library group label."""
    g = str(group)
    if "MISTRG" in g:
        return "xenograft"
    if "iPSC" in g:
        return "iPSC"
    if g == "young":
        return "young"
    if g == "aged":
        return "aged"
    if g.lower().startswith("anemia"):
        return "anemia"
    if g.startswith("MDS"):
        return "MDS"
    if "AML" in g:
        return "AML"
    return "other:" + g


# ---------------------------------------------------------------------------------------------------------
# stratified 2x2 statistics, vectorised over many tests
# ---------------------------------------------------------------------------------------------------------
def stratified_arrays(a, b, c, d, test, n_tests, breslow_day=True):
    """MH OR, RGB 95% CI, CMH p (continuity corrected), Breslow-Day p for every test over its strata.
    a, b, c, d: per-stratum counts; test: per-stratum test index. Uninformative strata (an empty row or
    column) are left out, as in the PILRA protocol."""
    a, b, c, d = (np.asarray(x, dtype=float) for x in (a, b, c, d))
    test = np.asarray(test, dtype=np.int64)
    out = dict(strata_informative=np.zeros(n_tests, dtype=np.int64))
    inf = ((a + b) > 0) & ((c + d) > 0) & ((a + c) > 0) & ((b + d) > 0)
    a, b, c, d, t = a[inf], b[inf], c[inf], d[inf], test[inf]
    n = a + b + c + d
    bc = lambda w: np.bincount(t, weights=w, minlength=n_tests)
    out["strata_informative"] = np.bincount(t, minlength=n_tests)
    has = out["strata_informative"] > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        s_ad, s_bc = bc(a * d / n), bc(b * c / n)
        orr = np.where(has, s_ad / s_bc, np.nan)
        apd = a + d
        va = bc(apd * a * d / n ** 2) / s_ad ** 2
        va += bc(apd * b * c / n ** 2 + (1 - apd / n) * a * d / n) / (s_ad * s_bc)
        va += bc((1 - apd / n) * b * c / n) / s_bc ** 2
        se = np.sqrt(va / 2)
        lor = np.log(orr)
        out["MH_OR"] = orr
        out["log_MH_OR"] = np.where(has, lor, np.nan)
        out["log_MH_OR_se"] = np.where(has, se, np.nan)
        out["CI_low"] = np.where(has, np.exp(lor - Z975 * se), np.nan)
        out["CI_high"] = np.where(has, np.exp(lor + Z975 * se), np.nan)
        num = np.abs(bc(a - (a + b) * (a + c) / n)) - 0.5
        den = bc((a + b) * (a + c) * (b + d) * (c + d) / (n ** 2 * (n - 1)))
        stat = num ** 2 / den
        out["CMH_stat"] = np.where(has, stat, np.nan)
        out["CMH_p"] = np.where(has, stats.chi2.sf(stat, 1), np.nan)
        if breslow_day:
            r = orr[t]
            qa = 1 - r
            qb = r * ((a + b) + (a + c)) + (d - a)
            qc = -r * (a + b) * (a + c)
            e11 = (-qb + np.sqrt(qb ** 2 - 4 * qa * qc)) / (2 * qa)
            v11 = 1 / (1 / e11 + 1 / ((a + c) - e11) + 1 / ((a + b) - e11) + 1 / ((d - a) + e11))
            bd = bc((a - e11) ** 2 / v11)
            df = out["strata_informative"] - 1
            out["BreslowDay_p"] = np.where(df >= 1, stats.chi2.sf(bd, np.maximum(df, 1)), np.nan)
        else:
            out["BreslowDay_p"] = np.full(n_tests, np.nan)
    return out


def bh(p):
    """Benjamini-Hochberg over the finite p-values; NaN stays NaN."""
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    m = int(ok.sum())
    if m:
        order = np.argsort(p[ok])
        ranked = p[ok][order] * m / np.arange(1, m + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        qq = np.empty(m)
        qq[order] = np.minimum(ranked, 1.0)
        q[ok] = qq
    return q


def cochran_q(lor, se, key):
    """Cochran's Q heterogeneity of log odds ratios within each key (finite lor and se > 0 only)."""
    df = pd.DataFrame(dict(key=key, lor=lor, se=se))
    df = df[np.isfinite(df.lor) & np.isfinite(df.se) & (df.se > 0)]
    if df.empty:
        return pd.DataFrame(columns=["key", "k", "Q", "Q_df", "Q_p", "I2", "pooled_log_OR_fixed"])
    df["w"] = 1 / df.se ** 2
    g = df.groupby("key")
    sw = g.w.sum()
    mu = (df.w * df.lor).groupby(df.key).sum() / sw
    df["mu"] = df.key.map(mu)
    q = (df.w * (df.lor - df.mu) ** 2).groupby(df.key).sum()
    k = g.size()
    res = pd.DataFrame(dict(k=k, Q=q, pooled_log_OR_fixed=mu)).reset_index().rename(columns={"index": "key"})
    res["Q_df"] = res.k - 1
    with np.errstate(divide="ignore", invalid="ignore"):
        res["Q_p"] = np.where(res.Q_df >= 1, stats.chi2.sf(res.Q, np.maximum(res.Q_df, 1)), np.nan)
        res["I2"] = np.where((res.Q_df >= 1) & (res.Q > 0), np.clip((res.Q - res.Q_df) / res.Q, 0, 1), np.nan)
    return res


# ---------------------------------------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------------------------------------
def load_libraries(donor_map, isv_metadata):
    dm = pd.read_csv(donor_map, sep="\t", dtype=str)
    md = pd.read_csv(isv_metadata, sep="\t", dtype=str)[["uid", "groups"]]
    lib = md.merge(dm[["uid", "donor"]], on="uid", how="left", validate="1:1")
    if lib.donor.isna().any():
        raise SystemExit("libraries without a donor: {}".format(lib.uid[lib.donor.isna()].tolist()))
    lib["covariate"] = lib.groups.map(covariate_of)
    lib["unit"] = lib.donor + "|" + lib.covariate
    return lib


def donor_genotypes(reads_dir, chrom, lib, error_rate):
    rows = []
    missing = []
    for r in lib.itertuples():
        p = os.path.join(reads_dir, r.uid, "{}.site_counts.tsv.gz".format(chrom))
        if not os.path.exists(p):
            missing.append(r.uid)
            continue
        sc = pd.read_csv(p, sep="\t", dtype={"site_id": str})
        sc["uid"], sc["donor"] = r.uid, r.donor
        rows.append(sc)
    if missing:
        raise SystemExit("site_counts missing for {} libraries: {}".format(len(missing), missing))
    sc = pd.concat(rows, ignore_index=True)
    cols = ["all_REF", "all_ALT", "all_OTHER", "all_SPLICE", "all_DEL",
            "primary_REF", "primary_ALT", "primary_OTHER"]
    g = sc.groupby(["site_id", "donor"], sort=False)
    dg = g[cols].sum().reset_index()
    dg["n_libraries"] = g.size().values
    lib_calls = sc.assign(call=sc.uid + "=" + sc.all_REF.astype(str) + "/" + sc.all_ALT.astype(str))
    dg["library_ref_alt"] = lib_calls.groupby(["site_id", "donor"], sort=False).call.agg(";".join).values
    calls = [genotype_call(int(x), int(y), error_rate) for x, y in zip(dg.all_REF, dg.all_ALT)]
    dg["genotype"] = [c[0] for c in calls]
    dg["GQ"] = [c[1] for c in calls]
    with np.errstate(invalid="ignore", divide="ignore"):
        dg["alt_fraction"] = dg.all_ALT / (dg.all_REF + dg.all_ALT)
    # tier and source belong to the donor panel, so they are taken per (site, donor), never per site
    meta = sc.drop_duplicates(["site_id", "donor"]).set_index(["site_id", "donor"])[
        ["pos", "ref", "alt", "vtype", "tier", "source"]]
    tiers_per_donor_site = sc.groupby(["site_id", "donor"], sort=False).tier.nunique()
    if (tiers_per_donor_site > 1).any():
        raise AssertionError("libraries of one donor disagree on the tier of {} sites".format(
            int((tiers_per_donor_site > 1).sum())))
    dg = dg.join(meta, on=["site_id", "donor"])
    return dg


def load_records(reads_dir, chrom, lib, het, region=None):
    """Primary-alignment REF/ALT records of heterozygous donors at their heterozygous sites, one per
    (read, site, molecule-index row). Records without an index row are counted, not returned.
    read_sites writes reads in BAM order (sorted by start), so with `region` a library stops at the first read
    that starts after the region and skips reads that end before it; an unsorted table stops the run."""
    recs = []
    counts = Counter()
    for r in lib.itertuples():
        p = os.path.join(reads_dir, r.uid, "{}.reads.tsv.gz".format(chrom))
        het_sites = het.get(r.donor)
        if not het_sites:
            continue
        with gzip.open(p, "rt") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            ix = {k: head.index(k) for k in ("secondary", "supplementary", "calls", "cb_rc", "read_name",
                                             "index_rows", "molecule_id", "ref_start", "ref_end", "introns")}
            last_start = -1
            for line in fh:
                row = line.rstrip("\n").split("\t")
                start, end = int(row[ix["ref_start"]]), int(row[ix["ref_end"]])
                if start < last_start:
                    raise SystemExit("{}: reads table is not sorted by start ({} after {})".format(p, start, last_start))
                last_start = start
                if region:
                    if start > region[1]:
                        counts["libraries_stopped_after_region"] += 1
                        break
                    if end < region[0]:
                        continue
                if row[ix["secondary"]] != "0" or row[ix["supplementary"]] != "0":
                    counts["skipped_secondary_or_supplementary"] += 1
                    continue
                calls = []
                for c in row[ix["calls"]].split(";"):
                    sid, cls = c.rsplit("=", 1)
                    if cls in ("REF", "ALT"):
                        full = chrom + ":" + sid
                        if full in het_sites:
                            calls.append((full, cls == "REF"))
                if not calls:
                    continue
                cb = row[ix["cb_rc"]]
                if not cb:
                    raise SystemExit("{}: heterozygote read {} without a cell barcode".format(r.uid, row[ix["read_name"]]))
                ir = row[ix["index_rows"]]
                idx = [x.split("|") for x in ir.split(";")] if ir else []
                counts["het_site_calls"] += len(calls)
                if not idx:
                    counts["het_site_calls_without_index_row"] += len(calls)
                    continue
                for full, is_ref in calls:
                    for gene, ct, fin in idx:
                        recs.append((full, r.donor, r.uid, r.unit, r.covariate, r.uid + ":" + cb, ct,
                                     gene, fin or "unassigned", int(is_ref), row[ix["molecule_id"]],
                                     start, end, row[ix["introns"]]))
    cols = ["site_id", "donor", "uid", "unit", "covariate", "cell", "cell_state", "gene", "final", "A",
            "molecule_id", "ref_start", "ref_end", "introns"]
    df = pd.DataFrame.from_records(recs, columns=cols)
    for c in ("site_id", "donor", "uid", "unit", "covariate", "cell", "cell_state", "gene", "final"):
        df[c] = df[c].astype("category")
    return df, counts


# ---------------------------------------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------------------------------------
def feature_tables(df, feature_col, stratum_cols, test_extra=()):
    """Per-stratum 2x2 counts for every (site, gene, [extra], feature) test. df has one row per record with
    columns A (1 = REF) and `feature_col`; for junctions df carries one row per (record, spanned intron)."""
    group = ["site_id", "gene"] + list(test_extra)
    skeys = group + list(stratum_cols)
    st = df.groupby(skeys, observed=True).A.agg(["sum", "count"]).reset_index()
    st = st.rename(columns={"sum": "nA", "count": "n"})
    st["nG"] = st.n - st.nA
    ft = df.groupby(skeys + [feature_col], observed=True).A.agg(["sum", "count"]).reset_index()
    ft = ft.rename(columns={"sum": "nAX", "count": "nX"})
    ft = ft.merge(st, on=skeys, how="left", validate="m:1")
    ft["a"], ft["c"] = ft.nAX, ft.nX - ft.nAX
    ft["b"], ft["d"] = ft.nA - ft.a, ft.nG - ft.c
    tests = ft[group + [feature_col]].drop_duplicates().reset_index(drop=True)
    tests["test"] = np.arange(len(tests))
    ft = ft.merge(tests, on=group + [feature_col], how="left", validate="m:1")
    # group-level totals (all strata of the group, with or without the feature)
    gt = st.groupby(group, observed=True).agg(strata_total=("n", "size"), reads=("n", "sum"),
                                               A_total=("nA", "sum"), G_total=("nG", "sum")).reset_index()
    tests = tests.merge(gt, on=group, how="left", validate="m:1")
    fx = ft.groupby("test").agg(A_feature=("a", "sum"), G_feature=("c", "sum"))
    tests = tests.join(fx, on="test")
    tests["A_other"] = tests.A_total - tests.A_feature
    tests["G_other"] = tests.G_total - tests.G_feature
    return ft, tests


def run_family(df, feature_col, stratum_cols, label, extra=(), breslow_day=True):
    ft, tests = feature_tables(df, feature_col, stratum_cols, extra)
    res = stratified_arrays(ft.a.values, ft.b.values, ft.c.values, ft.d.values, ft.test.values, len(tests),
                            breslow_day=breslow_day and not label.startswith("cell"))
    for k, v in res.items():
        tests[k] = v
    tests["strata"] = label
    with np.errstate(invalid="ignore", divide="ignore"):
        tests["REF_fraction_in_feature"] = tests.A_feature / (tests.A_feature + tests.G_feature)
        tests["REF_fraction_in_other"] = tests.A_other / (tests.A_other + tests.G_other)
    return tests.drop(columns=["test"])


def biallelic_subset(df):
    """Records of cells that carry both alleles of the site among the (site, gene) records."""
    g = df.groupby(["site_id", "gene", "cell"], observed=True).A.agg(["sum", "count"])
    bi = g[(g["sum"] > 0) & (g["sum"] < g["count"])].reset_index()[["site_id", "gene", "cell"]]
    return df.merge(bi, on=["site_id", "gene", "cell"], how="inner")


def junction_records(df):
    """One row per (record, intron) for every intron seen in the (site, gene) records whose flanks the read
    spans; `final` is replaced by the intron and A carried over."""
    out = []
    for (site, gene), g in df.groupby(["site_id", "gene"], observed=True):
        intr = g.introns.values
        sets = [set(x.split(";")) if x else set() for x in intr]
        uniq = sorted(set().union(*sets))
        if not uniq:
            continue
        s = np.array([int(x.split("-")[0]) for x in uniq])
        e = np.array([int(x.split("-")[1]) for x in uniq])
        st, en = g.ref_start.values[:, None], g.ref_end.values[:, None]
        span = (st <= s[None, :] - 1) & (en >= e[None, :] + 1)
        ri, ji = np.nonzero(span)
        inc = np.fromiter((uniq[j] in sets[i] for i, j in zip(ri, ji)), dtype=bool, count=len(ri))
        sub = g.iloc[ri].copy()
        sub["junction"] = [uniq[j] for j in ji]
        sub["junction_included"] = inc
        out.append(sub)
    if not out:
        return pd.DataFrame(columns=list(df.columns) + ["junction", "junction_included"])
    j = pd.concat(out, ignore_index=True)
    return j


def junction_family(jr, stratum_cols, label, extra=()):
    """Allele x inclusion of one intron among the reads spanning it: [[REF&in, REF&out],[ALT&in, ALT&out]]."""
    group = ["site_id", "gene", "junction"] + list(extra)
    keys = group + list(stratum_cols)
    jr = jr.assign(inA=(jr.A == 1) & jr.junction_included, inG=(jr.A == 0) & jr.junction_included,
                   outA=(jr.A == 1) & ~jr.junction_included, outG=(jr.A == 0) & ~jr.junction_included)
    t = jr.groupby(keys, observed=True)[["inA", "outA", "inG", "outG"]].sum().reset_index()
    tests = t[group].drop_duplicates().reset_index(drop=True)
    tests["test"] = np.arange(len(tests))
    t = t.merge(tests, on=group, how="left", validate="m:1")
    res = stratified_arrays(t.inA.values, t.outA.values, t.inG.values, t.outG.values, t.test.values, len(tests),
                            breslow_day=not label.startswith("cell"))
    agg = t.groupby("test").agg(strata_total=("inA", "size"), A_included=("inA", "sum"), A_excluded=("outA", "sum"),
                                G_included=("inG", "sum"), G_excluded=("outG", "sum"))
    tests = tests.join(agg, on="test")
    for k, v in res.items():
        tests[k] = v
    tests["strata"] = label
    tests["reads"] = tests[["A_included", "A_excluded", "G_included", "G_excluded"]].sum(axis=1)
    return tests.drop(columns=["test"])


def per_donor_fisher(df, feature_col):
    """PILRA per-donor Fisher exact for every (site, gene, feature, donor)."""
    st = df.groupby(["site_id", "gene", "donor"], observed=True).A.agg(["sum", "count"]).reset_index()
    st.columns = ["site_id", "gene", "donor", "nA", "n"]
    ft = df.groupby(["site_id", "gene", "donor", feature_col], observed=True).A.agg(["sum", "count"]).reset_index()
    ft.columns = ["site_id", "gene", "donor", feature_col, "nAX", "nX"]
    feats = ft[["site_id", "gene", feature_col]].drop_duplicates()
    full = feats.merge(st, on=["site_id", "gene"], how="inner")
    full = full.merge(ft, on=["site_id", "gene", "donor", feature_col], how="left").fillna({"nAX": 0, "nX": 0})
    full["A_feature"] = full.nAX.astype(int)
    full["G_feature"] = (full.nX - full.nAX).astype(int)
    full["A_other"] = (full.nA - full.nAX).astype(int)
    full["G_other"] = (full.n - full.nA - (full.nX - full.nAX)).astype(int)
    ors, ps = [], []
    for a, b, c, d in zip(full.A_feature, full.A_other, full.G_feature, full.G_other):
        o, p = stats.fisher_exact([[a, b], [c, d]])
        ors.append(o); ps.append(p)
    full["fisher_OR"], full["fisher_p"] = ors, ps
    return full.drop(columns=["nAX", "nX", "nA", "n"])


def ase_tables(dg_het, records):
    """Binomial REF fraction vs 0.5 per (site, donor) over all primary reads at the site, and per
    (site, donor, cell state) over index-joined records."""
    a = dg_het[["site_id", "donor", "primary_REF", "primary_ALT"]].copy()
    a["n"] = a.primary_REF + a.primary_ALT
    a = a[a.n > 0]
    a["REF_fraction"] = a.primary_REF / a.n
    a["binom_p"] = [stats.binomtest(int(k), int(n), 0.5).pvalue for k, n in zip(a.primary_REF, a.n)]
    # one record per (read, site) for the state table: drop the extra rows of multi-gene molecules
    r = records.drop_duplicates(["site_id", "molecule_id"])
    s = r.groupby(["site_id", "donor", "cell_state"], observed=True).A.agg(["sum", "count"]).reset_index()
    s.columns = ["site_id", "donor", "cell_state", "REF", "n"]
    s["ALT"] = s.n - s.REF
    s["REF_fraction"] = s.REF / s.n
    s["binom_p"] = [stats.binomtest(int(k), int(n), 0.5).pvalue for k, n in zip(s.REF, s.n)]
    return a, s


def _drop_unassigned(df, log, key):
    """Molecules without a final isoform (NULL in the molecule index) stay in every table's 'other' counts,
    as in the PILRA protocol, but 'unassigned versus the rest' is not an isoform test: drop those rows."""
    m = df["final"].astype(str) == "unassigned"
    log.setdefault("unassigned_feature_rows_dropped", {})[key] = int(m.sum())
    return df[~m].reset_index(drop=True)


def isoform_suite(rec, bi, out, suffix, log):
    """Allele x isoform tests: the four PILRA strata families, cell-state-specific and covariate-specific."""
    rf = run_family
    run_family_ = lambda *x, **k: _drop_unassigned(rf(*x, **k), log, suffix + ":" + x[3])
    iso = [run_family_(rec, "final", ["donor"], "donor"),
           run_family_(rec, "final", ["donor", "cell_state"], "donor_x_cell_state"),
           run_family_(bi, "final", ["donor"], "donor__biallelic_cells_only"),
           run_family_(bi, "final", ["cell"], "cell__biallelic_cells_only")]
    iso = pd.concat(iso, ignore_index=True)
    iso["BH_q_within_strata_family_chrom"] = np.nan
    for lab, g in iso.groupby("strata"):
        iso.loc[g.index, "BH_q_within_strata_family_chrom"] = bh(g.CMH_p.values)
    # undefined statistics are recorded, never filled: 0 or infinite MH OR (a zero margin across all strata)
    # and tests with no informative stratum
    log["isoform_tests_by_family" + suffix] = {lab: dict(
        tests=int(len(g)), with_p=int(g.CMH_p.notna().sum()), no_informative_stratum=int((g.strata_informative == 0).sum()),
        MH_OR_zero=int((g.MH_OR == 0).sum()), MH_OR_infinite=int(np.isinf(g.MH_OR).sum()),
        CI_undefined=int((g.strata_informative > 0).sum() - np.isfinite(g.CI_low).sum()))
        for lab, g in iso.groupby("strata")}
    log["isoform_tests" + suffix] = write(iso, os.path.join(out, "isoform_tests{}.tsv.gz".format(suffix)))
    log["isoform_per_donor_fisher" + suffix] = write(_drop_unassigned(per_donor_fisher(rec, "final"), log, suffix + ":fisher"),
                                                     os.path.join(out, "isoform_per_donor_fisher{}.tsv.gz".format(suffix)))
    st = run_family_(rec, "final", ["donor"], "state|donor", extra=("cell_state",))
    st["BH_q_within_family_chrom"] = bh(st.CMH_p.values)
    log["isoform_state_tests" + suffix] = write(st, os.path.join(out, "isoform_state_tests{}.tsv.gz".format(suffix)))
    key = st.site_id.astype(str) + "|" + st.gene.astype(str) + "|" + st.final.astype(str)
    qs = cochran_q(st.log_MH_OR.values, st.log_MH_OR_se.values, key.values)
    qs["Q_BH_q_chrom"] = bh(qs.Q_p.values)
    log["isoform_state_heterogeneity" + suffix] = write(
        qs, os.path.join(out, "isoform_state_heterogeneity{}.tsv.gz".format(suffix)))
    cv = run_family_(rec, "final", ["donor", "cell_state"], "covariate|donor_x_cell_state", extra=("covariate",))
    cv["BH_q_within_family_chrom"] = bh(cv.CMH_p.values)
    log["isoform_covariate_tests" + suffix] = write(cv, os.path.join(out, "isoform_covariate_tests{}.tsv.gz".format(suffix)))
    un = run_family_(rec, "final", ["cell_state"], "unit|cell_state", extra=("unit", "covariate"))
    log["isoform_unit_effects" + suffix] = write(un, os.path.join(out, "isoform_unit_effects{}.tsv.gz".format(suffix)))
    key = cv.site_id.astype(str) + "|" + cv.gene.astype(str) + "|" + cv.final.astype(str)
    qc = cochran_q(cv.log_MH_OR.values, cv.log_MH_OR_se.values, key.values)
    qc["Q_BH_q_chrom"] = bh(qc.Q_p.values)
    log["isoform_covariate_heterogeneity" + suffix] = write(
        qc, os.path.join(out, "isoform_covariate_heterogeneity{}.tsv.gz".format(suffix)))


def load_haplotype_records(phase_dir, chrom, lib, het):
    """Records of the haplotype-extended test: every primary read phased into the block of a heterozygous
    site gets that site's allele from its haplotype (A = 1 when the read's haplotype carries REF), whether
    or not the read covers the site. `covers_site` marks the reads that carry a call at the site."""
    donors = lib.drop_duplicates("donor")
    unit = dict(zip(lib.uid, lib.unit)); cov = dict(zip(lib.uid, lib.covariate))
    recs, counts = [], Counter()
    for d in donors.donor:
        sp = os.path.join(phase_dir, chrom, "{}.sites.tsv.gz".format(d))
        rp = os.path.join(phase_dir, chrom, "{}.reads.tsv.gz".format(d))
        if not het.get(d):
            continue
        if not (os.path.exists(sp) and os.path.exists(rp)):
            raise SystemExit("phase outputs missing for donor {} on {}".format(d, chrom))
        sites = pd.read_csv(sp, sep="\t", dtype={"site_id": str, "block": str})
        sites = sites[sites.site_id.isin(het[d]) & (sites.block_sites > 1)]
        by_block = defaultdict(list)
        for s, b, bit in zip(sites.site_id, sites.block, sites.phase_bit):
            by_block[b].append((s, int(bit)))
        reads = pd.read_csv(rp, sep="\t", dtype=str, keep_default_na=False)
        for r in reads.itertuples():
            if r.haplotype == "-1" or r.block not in by_block or not r.index_rows:
                counts["reads_ambiguous_or_single_site_block_or_unindexed"] += 1
                continue
            h = int(r.haplotype)
            idx = [x.split("|") for x in r.index_rows.split(";")]
            for s, bit in by_block[r.block]:
                for gene, ct, fin in idx:
                    recs.append((s, d, r.uid, unit[r.uid], cov[r.uid], r.uid + ":" + r.cb_rc, ct, gene,
                                 fin or "unassigned", int(h == bit), r.molecule_id))
    cols = ["site_id", "donor", "uid", "unit", "covariate", "cell", "cell_state", "gene", "final", "A", "molecule_id"]
    df = pd.DataFrame.from_records(recs, columns=cols)
    for c in ("site_id", "donor", "uid", "unit", "covariate", "cell", "cell_state", "gene", "final"):
        df[c] = df[c].astype("category")
    counts["haplotype_records"] = len(df)
    return df, counts


def write(df, path):
    df.to_csv(path + ".tmp", sep="\t", index=False, compression="gzip")
    os.replace(path + ".tmp", path)
    return len(df)


def build_arg_parser():
    ap = argparse.ArgumentParser(prog="python -m altanalyze3.components.snp.allelic_tests",
                                 description="Within-donor allelic expression, isoform and junction tests.")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--reads", required=True, help="read_sites output folder (<uid>/<chrom>.*.tsv.gz)")
    ap.add_argument("--donor_map", required=True)
    ap.add_argument("--isv_metadata", required=True)
    ap.add_argument("--error_rate", type=float, default=0.01, help="genotype likelihood error (PILRA: 0.01)")
    ap.add_argument("--sites", default=None, help="restrict to these site ids (comma list), for checks")
    ap.add_argument("--tiers", default=None,
                    help="comma list of panel tiers to test (high, medium, low, none_pass, named); default all")
    ap.add_argument("--region", default=None, help="START-END (1-based) on --chrom; output folder <chrom>_START_END")
    ap.add_argument("--genotypes", default=None,
                    help="reuse a donor_site_genotypes.tsv.gz written by --genotypes_only instead of recomputing")
    ap.add_argument("--no_junctions", action="store_true")
    ap.add_argument("--genotypes_only", action="store_true",
                    help="write donor_site_genotypes.tsv.gz and stop (input of phase_haplotypes)")
    ap.add_argument("--phase_dir", default=None,
                    help="phase_haplotypes output: adds the haplotype-extended isoform tests (*_haplotype files)")
    ap.add_argument("--out", required=True)
    return ap


def main(argv=None):
    a = build_arg_parser().parse_args(argv)
    t0 = time.time()
    region = tuple(int(x) for x in a.region.split("-")) if a.region else None
    out = os.path.join(a.out, a.chrom + ("_{}_{}".format(*region) if region else ""))
    os.makedirs(os.path.join(out, "logs"), exist_ok=True)
    lib = load_libraries(a.donor_map, a.isv_metadata)
    write(lib, os.path.join(out, "libraries_covariates.tsv.gz"))
    keep = set(a.sites.split(",")) if a.sites else None
    if a.genotypes:
        dg = pd.read_csv(a.genotypes, sep="\t", dtype={"site_id": str, "tier": str, "genotype": str, "donor": str})
        n_dg = len(dg)
    else:
        dg = donor_genotypes(a.reads, a.chrom, lib, a.error_rate)
        if keep:
            dg = dg[dg.site_id.isin(keep)]
        n_dg = write(dg, os.path.join(out, "donor_site_genotypes.tsv.gz"))
    # heterozygous sites under test: genotype 0/1, then the requested tiers, region and sites (selection, not
    # removal: every tier is tested in its own run)
    sel = dg.genotype == "0/1"
    if a.tiers:
        sel &= dg.tier.isin(a.tiers.split(","))
    if region:
        sel &= dg.pos.between(region[0], region[1])
    if keep:
        sel &= dg.site_id.isin(keep)
    dg_het = dg[sel]
    het = defaultdict(set)
    for s, d in zip(dg_het.site_id, dg_het.donor):
        het[d].add(s)
    log = dict(module="altanalyze3.components.snp.allelic_tests", argv=sys.argv if argv is None else list(argv),
               chrom=a.chrom, region=a.region, tiers=a.tiers, genotypes_input=a.genotypes, error_rate=a.error_rate,
               donor_site_rows=n_dg, genotypes=dict(Counter(dg.genotype)),
               het_donor_sites_all_tiers=int((dg.genotype == "0/1").sum()),
               het_donor_sites_selected=int(len(dg_het)), het_selected_by_tier=dict(Counter(dg_het.tier)))
    print("[allelic_tests] {}: {:,} donor-site rows, genotypes {}".format(a.chrom, n_dg, log["genotypes"]), flush=True)
    if a.genotypes_only:
        log.update(seconds=round(time.time() - t0, 1), stage="genotypes_only",
                   module_sha256=_sha256(os.path.abspath(__file__)))
        with open(os.path.join(out, "logs", "allelic_tests.genotypes_only.parameters.json"), "w") as fh:
            json.dump(log, fh, indent=2, default=str)
        return 0
    rec, counts = load_records(a.reads, a.chrom, lib, het, region)
    log.update(records=len(rec), **counts)
    print("[allelic_tests] records {:,}; {}".format(len(rec), dict(counts)), flush=True)
    if rec.empty:
        raise SystemExit("no heterozygote records on {}".format(a.chrom))
    bi = biallelic_subset(rec)
    log["records_biallelic_cells"] = len(bi)

    isoform_suite(rec, bi, out, "", log)

    if not a.no_junctions:
        jr = junction_records(rec)
        jb = biallelic_subset(rec)
        jbr = junction_records(jb) if len(jb) else jb
        log["junction_records"] = len(jr)
        jt = [junction_family(jr, ["donor"], "donor"),
              junction_family(jr, ["donor", "cell_state"], "donor_x_cell_state")]
        if len(jbr):
            jt += [junction_family(jbr, ["donor"], "donor__biallelic_cells_only"),
                   junction_family(jbr, ["cell"], "cell__biallelic_cells_only")]
        jt = pd.concat(jt, ignore_index=True)
        jt["BH_q_within_strata_family_chrom"] = np.nan
        for lab, g in jt.groupby("strata"):
            jt.loc[g.index, "BH_q_within_strata_family_chrom"] = bh(g.CMH_p.values)
        log["junction_tests"] = write(jt, os.path.join(out, "junction_tests.tsv.gz"))
        js = junction_family(jr, ["donor"], "state|donor", extra=("cell_state",))
        js["BH_q_within_family_chrom"] = bh(js.CMH_p.values)
        log["junction_state_tests"] = write(js, os.path.join(out, "junction_state_tests.tsv.gz"))
        jc = junction_family(jr, ["donor", "cell_state"], "covariate|donor_x_cell_state", extra=("covariate",))
        jc["BH_q_within_family_chrom"] = bh(jc.CMH_p.values)
        log["junction_covariate_tests"] = write(jc, os.path.join(out, "junction_covariate_tests.tsv.gz"))

    if a.phase_dir:
        hrec, hcounts = load_haplotype_records(a.phase_dir, a.chrom, lib, het)
        log.update({"haplotype_" + k: v for k, v in hcounts.items()})
        if len(hrec):
            # covering reads: the (site, molecule) pairs that also appear in the site-level records
            cover = set(zip(rec.site_id.astype(str), rec.molecule_id.astype(str)))
            hrec_cov = np.fromiter(((s, m) in cover for s, m in zip(hrec.site_id.astype(str), hrec.molecule_id.astype(str))),
                                   dtype=bool, count=len(hrec))
            log["haplotype_records_covering_site"] = int(hrec_cov.sum())
            log["haplotype_records_not_covering_site"] = int((~hrec_cov).sum())
            isoform_suite(hrec, biallelic_subset(hrec), out, "_haplotype", log)
    ase_site, ase_state = ase_tables(dg_het, rec)
    ase_site["BH_q_chrom"] = bh(ase_site.binom_p.values)
    ase_state["BH_q_chrom"] = bh(ase_state.binom_p.values)
    log["ase_site_donor"] = write(ase_site, os.path.join(out, "ase_site_donor.tsv.gz"))
    log["ase_site_donor_state"] = write(ase_state, os.path.join(out, "ase_site_donor_state.tsv.gz"))
    # reference bias check: median REF fraction over heterozygous sites per donor (expected 0.5)
    rb = ase_site.groupby("donor").REF_fraction.agg(["median", "mean", "size"]).reset_index()
    log["ref_fraction_by_donor"] = rb.round(4).to_dict(orient="records")
    log.update(seconds=round(time.time() - t0, 1), finished=datetime.datetime.now().isoformat(timespec="seconds"),
               module_sha256=_sha256(os.path.abspath(__file__)), numpy=np.__version__, pandas=pd.__version__,
               python=sys.version.split()[0])
    with open(os.path.join(out, "logs", "allelic_tests.parameters.json"), "w") as fh:
        json.dump(log, fh, indent=2, default=str)
    print("[allelic_tests] {} done in {:.0f}s: {}".format(a.chrom, time.time() - t0,
                                                          {k: v for k, v in log.items() if k.endswith("tests")}), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
