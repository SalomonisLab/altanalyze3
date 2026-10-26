"""Isoform switches between pairs of cell states, tested within donors.

A switch between states A and B is a pair of isoforms of one gene: isoform i has higher usage in A and
isoform j has higher usage in B. Usage is isoform reads / gene reads in one donor's pseudobulk, as computed by
isoform_matrix.cpm_and_ratio_from_counts (missing where the gene CPM < 1 in that pseudobulk).

Candidates are isoforms that are differential in at least one state-vs-Others comparison
(|DiffMeans| > log2(--fold) and ebayesAdjPval < --adjp), are protein coding (predicted protein length > 0 and Not-NMD
in --protein_classes, with no length-relative-to-longest rule; --all_isoforms drops this filter), and belong to a gene
with >= 2 such isoforms.
Every state of the group is kept. For every candidate and every pair of states, an isoform is reported when the read
depth exceeds --min_reads in A or in B and >= 1 donor has a usage in both states, and tested when >= --min_donors
(>= 2) donors do. --read_rule sets the depth: isoform (pooled isoform reads) or gene (pooled reads of its gene).
--read_scope sets the pseudobulks pooled for the depth: group (the --group pseudobulks only; default) or cohort (every
pseudobulk of the state, all groups). Usage and the test always use the --group donors:
  delta     usage in A - usage in B, paired within each donor
  test      moderated one-sample t on the paired deltas; the variance prior (s0^2, d0) is fitted over every isoform
            tested for that state pair (Smyth 2004, limma fitFDist and squeezeVar)
  FDR       Benjamini-Hochberg over every isoform x state-pair test of the run
Rows with fewer than --min_donors paired donors keep their delta and reads, with no P (no across-donor variance).
An isoform is higher in A when FDR < --alpha and mean delta >= --min_delta (B: mean delta <= -min_delta). With
--read_rule isoform it also needs pooled isoform reads > --min_reads in its higher state; with --read_rule gene the
gene needs pooled reads > --min_reads in A or in B. A switch is one isoform higher in A plus another isoform of the
same gene higher in B.

python -m altanalyze3.components.long_read.isoform_state_switches \\
    --counts isoform_combined_pseudo_cluster_counts.h5ad --groups mds_metadata_bam_healthy_fixed.txt \\
    --group healthy --diff_dir analysis/healthy_diff/diff-cluster-isoform \\
    --protein_classes ALL_catalog_protein_classes.tsv --outdir isoform_state_switches

Group mode (--compare_groups A,B, e.g. aged,young) tests two donor groups WITHIN each cell state instead of two states
within donors. Candidates come from the covariate tables <A>-<B>-<state>_stats-annotated.txt (either name order) in
--diff_dir: |DiffMeans| > log2(--fold) with ebayesPval < --whitelist_rawp (or ebayesAdjPval < --adjp), plus, with
--event_dir, both isoforms named by a dPSI covariate event with |dPSI| > --event_dpsi and ebayesPval < --event_rawp.
  delta     mean usage of the A donors - mean usage of the B donors in the state (unpaired)
  test      moderated two-sample t: pooled variance on nA + nB - 2 df, prior (s0^2, d0) fitted over every isoform
            tested in that state (limma fitFDist and squeezeVar)
  FDR       Benjamini-Hochberg over every isoform x state test of the run
A state is tested when both groups have >= --min_donors donors with a usage. The read gate pools each group's
pseudobulks of the state separately (--read_scope group). A switch is one isoform higher in A plus another isoform
of the same gene higher in B, in the same state. --required_switch FEATURE_A,FEATURE_B names a positive control that
must be called in >= 1 state; the run writes its outputs and then fails when a control is missing.
"""
import argparse
import glob
import hashlib
import io
import json
import os
import platform
import re
import sys
import time
from collections import defaultdict
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import sparse, special, stats

STATS_RE = re.compile(r"^(?P<prefix>.*?)-(?P<state>.+)-Others-(?P=state)_stats-annotated\.txt$")


class _Log:
    def __init__(self, path=None):
        self.fh = open(path, "w") if path else None

    def __call__(self, *msg):
        line = " ".join(str(m) for m in msg)
        print(line, flush=True)
        if self.fh:
            self.fh.write(line + "\n")
            self.fh.flush()


def _sha256(path, chunk=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------- moderated statistics (limma, Smyth 2004)
def trigamma_inverse(x):
    """y with trigamma(y) = x (limma trigammaInverse, Newton iteration)."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    y = np.full(x.shape, np.nan)
    big, small = x > 1e7, (x >= 0) & (x < 1e-6)
    mid = (x >= 1e-6) & (x <= 1e7)
    y[big] = 1.0 / np.sqrt(x[big])
    y[small] = 1.0 / x[small]
    if mid.any():
        v = x[mid]
        z = 0.5 + 1.0 / v
        for _ in range(50):
            tri = special.polygamma(1, z)
            dif = tri * (1.0 - tri / v) / special.polygamma(2, z)
            z = z + dif
            if np.max(-dif / z) < 1e-8:
                break
        else:
            raise RuntimeError("trigamma_inverse did not converge in 50 iterations")
        y[mid] = z
    return y


def fit_f_dist(s2, df):
    """Prior variance s0^2 and prior df d0 from sample variances s2 with df (limma fitFDist, no covariate).
    Zero variances are offset to 1e-5 * median for the fit only, as limma does. Returns (s0^2, d0, info)."""
    s2 = np.asarray(s2, dtype=float)
    df = np.asarray(df, dtype=float)
    ok = np.isfinite(s2) & np.isfinite(df) & (s2 > -1e-15) & (df > 1e-15)
    x, d = np.maximum(s2[ok], 0.0), df[ok]
    info = dict(prior_features=int(x.size), zero_variance=int((x == 0).sum()), median_variance_zero=False)
    if x.size < 2:
        return np.nan, np.nan, info
    m = np.median(x)
    if m == 0:
        info["median_variance_zero"] = True
        m = 1.0
    x = np.maximum(x, 1e-5 * m)
    e = np.log(x) - special.digamma(d / 2) + np.log(d / 2)
    emean = e.mean()
    evar = np.sum((e - emean) ** 2) / (x.size - 1) - np.mean(special.polygamma(1, d / 2))
    if evar > 0:
        d0 = float(2 * trigamma_inverse(evar)[0])
        s0 = float(np.exp(emean + special.digamma(d0 / 2) - np.log(d0 / 2)))
    else:
        d0, s0 = np.inf, float(np.exp(emean))
    return s0, d0, info


def moderated_one_sample_t(mean, s2, n, s0, d0):
    """limma squeezeVar + eBayes for a one-sample mean of n paired differences. Without a prior (s0 NaN) the
    ordinary t is returned. Returns (t, df_total, two-sided P, posterior variance)."""
    df = n - 1.0
    if not np.isfinite(s0):
        s2post, df_total = s2.astype(float), df
    elif np.isinf(d0):
        s2post, df_total = np.full(s2.shape, s0), np.full(s2.shape, df.sum())
    else:
        s2post = (df * s2 + d0 * s0) / (df + d0)
        df_total = np.minimum(df + d0, df.sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        t = mean / np.sqrt(s2post / n)
    p = 2.0 * stats.t.sf(np.abs(t), df_total)
    return t, df_total, p, s2post


def moderated_two_sample_t(delta, s2, nA, nB, s0, d0):
    """limma squeezeVar + eBayes for a difference of two independent group means. s2 is the pooled variance on
    nA + nB - 2 df. Without a prior (s0 NaN) the ordinary pooled t is returned. A posterior variance of 0 leaves t
    undefined, so its t and P are NaN. Returns (t, df_total, two-sided P, posterior variance)."""
    df = nA + nB - 2.0
    if not np.isfinite(s0):
        s2post, df_total = s2.astype(float), df
    elif np.isinf(d0):
        s2post, df_total = np.full(s2.shape, s0), np.full(s2.shape, df.sum())
    else:
        s2post = (df * s2 + d0 * s0) / (df + d0)
        df_total = np.minimum(df + d0, df.sum())
    se = np.sqrt(s2post * (1.0 / nA + 1.0 / nB))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(se > 0, delta / se, np.nan)
    p = 2.0 * stats.t.sf(np.abs(t), df_total)
    return t, df_total, p, s2post


def bh_fdr(p):
    """Benjamini-Hochberg adjusted P for finite p; NaN stays NaN."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    m = int(ok.sum())
    if m == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv, kind="mergesort")
    ranked = pv[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(m)
    res[order] = np.minimum(adj, 1.0)
    out[ok] = res
    return out


# ---------------------------------------------------------------- inputs
def load_differential_whitelist(diff_dir, fold, adjp, pattern="*_stats-annotated.txt", log=print):
    """Isoforms with |DiffMeans| > log2(fold) and ebayesAdjPval < adjp in >= 1 state-vs-Others stats file.
    oncosplice.metadataAnalysis sorts the two group labels, so a state named after 'Others' is written as
    DiffMeans_Others_m_<state>; that orientation is checked against the Mean columns and negated to state - Others."""
    paths = sorted(glob.glob(os.path.join(diff_dir, pattern)))
    if not paths:
        raise FileNotFoundError(f"no files match {os.path.join(diff_dir, pattern)}")
    cut = float(np.log2(fold))
    up, down, symbols, manifest = defaultdict(list), defaultdict(list), {}, []
    for path in paths:
        m = STATS_RE.match(os.path.basename(path))
        if m is None:
            raise ValueError(f"{path}: name is not <prefix>-<state>-Others-<state>_stats-annotated.txt")
        state = m.group("state")
        raw = open(path, "rb").read()
        header = raw[:raw.index(b"\n")].decode().rstrip("\r").split("\t")
        fwd, rev = f"DiffMeans_{state}_m_Others", f"DiffMeans_Others_m_{state}"
        fc = [c for c in (fwd, rev) if c in header]
        if len(fc) != 1:
            raise ValueError(f"{path}: expected one of {fwd} or {rev} in the header")
        fc = fc[0]
        d = pd.read_csv(io.BytesIO(raw), sep="\t",
                        usecols=["Feature", "Gene_Symbol", fc, "ebayesAdjPval", f"Mean_{state}", "Mean_Others"],
                        dtype={"Feature": str, "Gene_Symbol": str}, low_memory=False)
        if not d.Feature.is_unique:
            raise ValueError(f"{path}: duplicate Feature rows")
        raw_fc = pd.to_numeric(d[fc], errors="coerce").to_numpy(float)
        ms = pd.to_numeric(d[f"Mean_{state}"], errors="coerce").to_numpy(float)
        mo = pd.to_numeric(d["Mean_Others"], errors="coerce").to_numpy(float)
        expect = ms - mo if fc == fwd else mo - ms
        both = np.isfinite(raw_fc) & np.isfinite(expect)
        if not np.allclose(raw_fc[both], expect[both], rtol=1e-4, atol=1e-5):
            raise ValueError(f"{path}: {fc} does not equal the difference of its Mean columns")
        f = raw_fc if fc == fwd else -raw_fc
        q = pd.to_numeric(d["ebayesAdjPval"], errors="coerce").to_numpy(float)
        tested = np.isfinite(f) & np.isfinite(q)
        sig = tested & (np.abs(f) > cut) & (q < adjp)
        feats = d.Feature.to_numpy()
        for x in feats[sig & (f > 0)]:
            up[x].append(state)
        for x in feats[sig & (f < 0)]:
            down[x].append(state)
        sym = d.loc[sig, ["Feature", "Gene_Symbol"]].dropna()
        symbols.update(zip(sym.Feature, sym.Gene_Symbol))
        manifest.append(dict(state=state, path=os.path.abspath(path), sha256=hashlib.sha256(raw).hexdigest(),
                             fold_column=fc, negated=fc == rev, mean_check_rows=int(both.sum()), rows=len(d),
                             tested=int(tested.sum()), differential=int(sig.sum()),
                             up=int((sig & (f > 0)).sum()), down=int((sig & (f < 0)).sum())))
        log(f"[whitelist] {state}: {int(sig.sum()):,} differential of {int(tested.sum()):,} tested "
            f"({len(d):,} rows, {len(d) - int(tested.sum()):,} untested; {fc}"
            f"{', negated to state - Others' if fc == rev else ''}; Mean check {int(both.sum()):,} rows)")
    feats = sorted(set(up) | set(down))
    wl = pd.DataFrame({"feature": feats,
                       "gene_symbol": [symbols.get(x, "") for x in feats],
                       "de_up_states": [";".join(up.get(x, [])) for x in feats],
                       "de_down_states": [";".join(down.get(x, [])) for x in feats]})
    wl["n_de_states"] = [len(up.get(x, [])) + len(down.get(x, [])) for x in feats]
    return wl, manifest


def _covariate_files(diff_dir, group_A, group_B, pattern):
    """(state, path) for <A>-<B>-<state> or <B>-<A>-<state> stats files, in sorted path order."""
    paths = sorted(glob.glob(os.path.join(diff_dir, pattern)))
    if not paths:
        raise FileNotFoundError(f"no files match {os.path.join(diff_dir, pattern)}")
    out = []
    for path in paths:
        base = os.path.basename(path)
        if not base.endswith("_stats-annotated.txt"):
            continue
        stem = base[:-len("_stats-annotated.txt")]
        for x, y in ((group_A, group_B), (group_B, group_A)):
            if stem.startswith(f"{x}-{y}-"):
                out.append((stem[len(f"{x}-{y}-"):], path))
                break
    if not out:
        raise ValueError(f"no {group_A}-{group_B}-<state> or {group_B}-{group_A}-<state> stats files in {diff_dir}")
    states = [s for s, _ in out]
    if len(set(states)) != len(states):
        raise ValueError(f"{diff_dir}: more than one stats file per state")
    return out


def _oriented_fold(path, d, header, group_A, group_B, mean_A, mean_B):
    """DiffMeans oriented to A - B, checked against the Mean columns (the stats writer sorts the group labels)."""
    fwd, rev = f"DiffMeans_{group_A}_m_{group_B}", f"DiffMeans_{group_B}_m_{group_A}"
    fc = [c for c in (fwd, rev) if c in header]
    if len(fc) != 1:
        raise ValueError(f"{path}: expected one of {fwd} or {rev} in the header")
    fc = fc[0]
    raw_fc = pd.to_numeric(d[fc], errors="coerce").to_numpy(float)
    expect = mean_A - mean_B if fc == fwd else mean_B - mean_A
    both = np.isfinite(raw_fc) & np.isfinite(expect)
    if not np.allclose(raw_fc[both], expect[both], rtol=1e-4, atol=1e-5):
        raise ValueError(f"{path}: {fc} does not equal the difference of its Mean columns")
    return (raw_fc if fc == fwd else -raw_fc), fc, int(both.sum())


def load_covariate_whitelist(diff_dir, group_A, group_B, fold, adjp, rawp=None, pattern="*_stats-annotated.txt",
                             log=print):
    """Isoforms with |DiffMeans A - B| > log2(fold) and P < threshold in >= 1 state's covariate stats file.
    The P is ebayesPval < rawp when rawp is given, else ebayesAdjPval < adjp."""
    cut = float(np.log2(fold))
    pcol, pthr = ("ebayesPval", rawp) if rawp is not None else ("ebayesAdjPval", adjp)
    up, down, symbols, manifest = defaultdict(list), defaultdict(list), {}, []
    for state, path in _covariate_files(diff_dir, group_A, group_B, pattern):
        raw = open(path, "rb").read()
        header = raw[:raw.index(b"\n")].decode().rstrip("\r").split("\t")
        need = ["Feature", "Gene_Symbol", pcol, f"Mean_{group_A}", f"Mean_{group_B}"]
        missing = [c for c in need if c not in header]
        if missing:
            raise ValueError(f"{path}: missing columns {missing}")
        fcs = [c for c in header if c.startswith("DiffMeans_")]
        d = pd.read_csv(io.BytesIO(raw), sep="\t", usecols=need + fcs, dtype={"Feature": str, "Gene_Symbol": str},
                        low_memory=False)
        if not d.Feature.is_unique:
            raise ValueError(f"{path}: duplicate Feature rows")
        mA = pd.to_numeric(d[f"Mean_{group_A}"], errors="coerce").to_numpy(float)
        mB = pd.to_numeric(d[f"Mean_{group_B}"], errors="coerce").to_numpy(float)
        f, fc, checked = _oriented_fold(path, d, header, group_A, group_B, mA, mB)
        p = pd.to_numeric(d[pcol], errors="coerce").to_numpy(float)
        tested = np.isfinite(f) & np.isfinite(p)
        sig = tested & (np.abs(f) > cut) & (p < pthr)
        feats = d.Feature.to_numpy()
        for x in feats[sig & (f > 0)]:
            up[x].append(state)
        for x in feats[sig & (f < 0)]:
            down[x].append(state)
        sym = d.loc[sig, ["Feature", "Gene_Symbol"]].dropna()
        symbols.update(zip(sym.Feature, sym.Gene_Symbol))
        manifest.append(dict(state=state, path=os.path.abspath(path), sha256=hashlib.sha256(raw).hexdigest(),
                             fold_column=fc, p_column=pcol, p_threshold=pthr, mean_check_rows=checked, rows=len(d),
                             tested=int(tested.sum()), differential=int(sig.sum()),
                             up=int((sig & (f > 0)).sum()), down=int((sig & (f < 0)).sum())))
        log(f"[whitelist] {state}: {int(sig.sum()):,} isoforms with |{fc}| > log2({fold}) and {pcol} < {pthr} of "
            f"{int(tested.sum()):,} tested ({len(d):,} rows; Mean check {checked:,} rows)")
    feats = sorted(set(up) | set(down))
    wl = pd.DataFrame({"feature": feats, "gene_symbol": [symbols.get(x, "") for x in feats],
                       "de_up_states": [";".join(up.get(x, [])) for x in feats],
                       "de_down_states": [";".join(down.get(x, [])) for x in feats]})
    wl["n_de_states"] = [len(up.get(x, [])) + len(down.get(x, [])) for x in feats]
    return wl, manifest


def load_event_whitelist(event_dir, group_A, group_B, dpsi, rawp, pattern="*_stats-annotated.txt", log=print,
                         iso_to_genes=None):
    """Both isoforms named by a dPSI covariate event (Isoform_1|Length, Isoform_2|Length) with |dPSI A - B| > dpsi
    and ebayesPval < rawp in >= 1 state. Event rows that name no isoform on a side are counted and contribute nothing
    for that side. An isoform's gene is the gene the catalogue registers it under (iso_to_genes: bare isoform id ->
    set of genes); a junction can span two genes, so the first gene of the event side is only the fallback when no
    catalogue is given, and the tie-break when the catalogue lists the id under several genes."""
    hits, symbols, manifest = defaultdict(set), {}, []
    rekeyed, unknown, ambiguous = set(), set(), set()
    for state, path in _covariate_files(event_dir, group_A, group_B, pattern):
        raw = open(path, "rb").read()
        header = raw[:raw.index(b"\n")].decode().rstrip("\r").split("\t")
        need = ["Gene_Symbol", "Feature", "Isoform_1|Length", "Isoform_2|Length", "ebayesPval",
                f"Mean_{group_A}", f"Mean_{group_B}"]
        missing = [c for c in need if c not in header]
        if missing:
            raise ValueError(f"{path}: missing columns {missing}")
        fcs = [c for c in header if c.startswith("DiffMeans_")]
        d = pd.read_csv(io.BytesIO(raw), sep="\t", usecols=need + fcs, dtype=str, low_memory=False)
        mA = pd.to_numeric(d[f"Mean_{group_A}"], errors="coerce").to_numpy(float)
        mB = pd.to_numeric(d[f"Mean_{group_B}"], errors="coerce").to_numpy(float)
        f, fc, checked = _oriented_fold(path, d, header, group_A, group_B, mA, mB)
        p = pd.to_numeric(d["ebayesPval"], errors="coerce").to_numpy(float)
        sig = np.isfinite(f) & np.isfinite(p) & (np.abs(f) > dpsi) & (p < rawp)
        s = d[sig]
        sides = s.Feature.str.split("|")
        if not sides.str.len().eq(2).all():
            raise ValueError(f"{path}: an event Feature does not have exactly two sides")
        empty = 0
        for k in (1, 2):
            iso = s[f"Isoform_{k}|Length"].fillna("").str.split("|").str[0].str.strip()
            gene = sides.str[k - 1].str.split(":").str[0]
            ok = iso.ne("") & iso.ne("nan")
            empty += int((~ok).sum())
            for g, i, sym in zip(gene[ok], iso[ok], s.Gene_Symbol[ok]):
                if iso_to_genes is not None:
                    reg = iso_to_genes.get(i)
                    if not reg:
                        unknown.add(f"{g}:{i}")
                        continue
                    if len(reg) > 1:
                        if g not in reg:
                            ambiguous.add(f"{i} ({';'.join(sorted(reg))})")
                            continue
                    else:
                        (only,) = reg
                        if only != g:
                            rekeyed.add(f"{g}:{i} -> {only}:{i}")
                            g = only
                hits[f"{g}:{i}"].add(state)
                symbols.setdefault(f"{g}:{i}", sym)
        manifest.append(dict(state=state, path=os.path.abspath(path), sha256=hashlib.sha256(raw).hexdigest(),
                             fold_column=fc, mean_check_rows=checked, rows=len(d), significant_events=int(sig.sum()),
                             event_sides_without_isoform=empty))
        log(f"[events] {state}: {int(sig.sum()):,} events with |{fc}| > {dpsi} and ebayesPval < {rawp}; "
            f"event sides naming no isoform: {empty:,}")
    if unknown or ambiguous:
        raise ValueError(f"event isoforms the catalogue cannot place: {len(unknown)} absent "
                         f"({sorted(unknown)[:10]}), {len(ambiguous)} under several genes, none the event gene "
                         f"({sorted(ambiguous)[:10]})")
    if rekeyed:
        log(f"[events] {len(rekeyed):,} isoforms keyed to their registered gene, not the first gene of a two-gene "
            f"junction: {', '.join(sorted(rekeyed)[:20])}{' ...' if len(rekeyed) > 20 else ''}")
    feats = sorted(hits)
    ev = pd.DataFrame({"feature": feats, "gene_symbol_event": [symbols.get(x, "") for x in feats],
                       "event_states": [";".join(sorted(hits[x])) for x in feats]})
    return ev, manifest


def load_protein_classes(path):
    p = pd.read_csv(path, sep="\t", dtype=str, usecols=["gene", "isoform", "class", "protein_length", "nmd_status"])
    p["feature"] = p["gene"] + ":" + p["isoform"]
    if not p.feature.is_unique:
        raise ValueError(f"{path}: duplicate gene:isoform rows")
    return p.set_index("feature")


# ---------------------------------------------------------------- protein outcome of age-differential isoforms
OUTCOME_COLORS = {"protein coding": "#1F77B4", "NMD": "#D62728"}


def protein_outcome_classes(path):
    """Per isoform: protein length, the gene's longest predicted protein, and the outcome used by protein-outcome:
    'no protein' (length 0 or missing), 'same length' (length equals the gene's longest protein), else 'NMD'
    (nmd_status names NMD and is not Not-NMD) or 'protein coding'."""
    p = pd.read_csv(path, sep="\t", dtype=str,
                    usecols=["gene", "isoform", "protein_length", "nmd_status", "longest_gene_protein"])
    p["feature"] = p["gene"] + ":" + p["isoform"]
    if not p.feature.is_unique:
        raise ValueError(f"{path}: duplicate gene:isoform rows")
    L = pd.to_numeric(p.protein_length, errors="coerce")
    M = pd.to_numeric(p.longest_gene_protein, errors="coerce")
    nmd = p.nmd_status.fillna("").str.contains("NMD") & ~p.nmd_status.fillna("").str.startswith("Not")
    p["outcome"] = np.where(~(L > 0), "no protein", np.where(L == M, "same length",
                                                             np.where(nmd, "NMD", "protein coding")))
    p["length_vs_longest"] = np.where(L > 0, L - M, np.nan)
    return p.set_index("feature")


def protein_outcome_from_tables(diff_dir, group_A, group_B, classes, fold=1.5, rawp=0.05, min_donors=3,
                                fold_scale="log2", pattern="*_stats-annotated.txt", log=print,
                                min_usage_shift=None):
    """Every isoform differential between group_A and group_B in each state's covariate table: >= min_donors donors in
    each group, ebayesPval < rawp, and |DiffMeans| > log2(fold) (fold_scale log2, expression tables) or
    max(Mean_A / Mean_B, Mean_B / Mean_A) > fold (fold_scale ratio, usage tables). Returns one row per state and
    isoform with its direction and protein outcome, plus a per-state manifest.
    ``min_usage_shift`` explicitly replaces the ratio-fold gate with a strict
    absolute native DiffMeans threshold on fractional usage tables; no defaults change.
    """
    rows, manifest = [], []
    if min_usage_shift is not None:
        if fold_scale != 'ratio' or not np.isfinite(min_usage_shift) or not 0 <= min_usage_shift <= 1:
            raise ValueError('min_usage_shift requires ratio/usage tables and a finite threshold in [0,1]')
    cut = float(np.log2(fold))
    for state, path in _covariate_files(diff_dir, group_A, group_B, pattern):
        raw = open(path, "rb").read()
        header = raw[:raw.index(b"\n")].decode().rstrip("\r").split("\t")
        need = ["Feature", "ebayesPval", f"Mean_{group_A}", f"Mean_{group_B}", f"N_{group_A}", f"N_{group_B}"]
        missing = [c for c in need if c not in header]
        if missing:
            raise ValueError(f"{path}: missing columns {missing}")
        fcs = [c for c in header if c.startswith("DiffMeans_")]
        d = pd.read_csv(io.BytesIO(raw), sep="\t", usecols=need + fcs, dtype={"Feature": str}, low_memory=False)
        if not d.Feature.is_unique:
            raise ValueError(f"{path}: duplicate Feature rows")
        mA = pd.to_numeric(d[f"Mean_{group_A}"], errors="coerce").to_numpy(float)
        mB = pd.to_numeric(d[f"Mean_{group_B}"], errors="coerce").to_numpy(float)
        f, fc, checked = _oriented_fold(path, d, header, group_A, group_B, mA, mB)
        p = pd.to_numeric(d.ebayesPval, errors="coerce").to_numpy(float)
        nA = pd.to_numeric(d[f"N_{group_A}"], errors="coerce").to_numpy(float)
        nB = pd.to_numeric(d[f"N_{group_B}"], errors="coerce").to_numpy(float)
        tested = np.isfinite(f) & np.isfinite(p) & (nA >= min_donors) & (nB >= min_donors)
        if min_usage_shift is not None:
            for means in (mA,mB):
                observed=means[np.isfinite(means)]
                if np.any((observed<0)|(observed>1)):
                    raise ValueError(f'{path}: usage means outside [0,1]; not fractional isoform usage')
            big = np.abs(f) > min_usage_shift
        elif fold_scale == "log2":
            big = np.abs(f) > cut
        else:
            with np.errstate(divide="ignore", invalid="ignore"):
                r = np.where(mB > 0, mA / mB, np.inf)
                big = (r > fold) | (r < 1.0 / fold)
            big &= (mA > 0) | (mB > 0)
        sig = tested & big & (p < rawp)
        feats = d.Feature.to_numpy()[sig]
        oc = classes.reindex(feats)
        unknown = int(oc.outcome.isna().sum())
        if unknown:
            raise ValueError(f"{path}: {unknown} differential isoforms are absent from the protein classes table")
        rows.append(pd.DataFrame({"state": state, "feature": feats,
                                  "direction": np.where(f[sig] > 0, f"higher in {group_A}", f"higher in {group_B}"),
                                  "diff": f[sig], "p": p[sig], "outcome": oc.outcome.to_numpy(),
                                  "protein_length": oc.protein_length.to_numpy(),
                                  "length_vs_longest": oc.length_vs_longest.to_numpy()}))
        manifest.append(dict(state=state, path=os.path.abspath(path), sha256=hashlib.sha256(raw).hexdigest(),
                             fold_column=fc, mean_check_rows=checked, rows=len(d), tested=int(tested.sum()),
                             differential=int(sig.sum())))
        log(f"[protein-outcome] {state}: {int(sig.sum()):,} differential of {int(tested.sum()):,} tested isoforms")
    return pd.concat(rows, ignore_index=True), manifest


def protein_outcome_from_switches(switches, classes, group_A, group_B):
    """Both isoforms of every called switch, with direction (isoform_higher_in_A -> higher in A) and outcome."""
    sw = pd.read_csv(switches, sep="\t") if isinstance(switches, str) else switches
    parts = []
    for col, lab in (("isoform_higher_in_A", f"higher in {group_A}"), ("isoform_higher_in_B", f"higher in {group_B}")):
        oc = classes.reindex(sw[col].to_numpy())
        if oc.outcome.isna().any():
            raise ValueError(f"{int(oc.outcome.isna().sum())} switch isoforms are absent from the protein classes table")
        parts.append(pd.DataFrame({"state": sw.state.to_numpy(), "feature": sw[col].to_numpy(), "direction": lab,
                                   "gene_symbol": sw.gene_symbol.to_numpy(), "outcome": oc.outcome.to_numpy(),
                                   "protein_length": oc.protein_length.to_numpy(),
                                   "length_vs_longest": oc.length_vs_longest.to_numpy()}))
    return pd.concat(parts, ignore_index=True)


def predominant_isoforms(counts_h5ad, groups_file, groups, genes=None, log=print, cpm_threshold=1,
                         return_usage=False, state_pools=None):
    """Predominant isoform of each gene in each state and group: the isoform with the highest mean donor usage, over
    every isoform of the gene in the counts. Usage is as in load_group_usage (isoform reads / gene reads in one donor's
    pseudobulk, missing where the gene CPM < cpm_threshold). The donors averaged are the ones where the gene is
    expressed, the same for every isoform of the gene, so the predominant isoform has the largest summed usage; ties go
    to the first feature id. Returns one row per state, group and gene: predominant, mean_usage, donors.
    With return_usage=True also return all nonzero context means and native pooled reads, before selecting
    the maximum. All features remain in the gene denominator; absent sparse entries are zero only when
    that gene has expressing donors in the context. Optional state_pools append explicitly named
    donor-level count sums before the same native CPM/usage calculation; every original context
    remains. Pools may overlap and are descriptive views, never independent replicate additions.
    These options do not fit or alter a statistical test."""
    import anndata as ad
    from .comparisons import _resolve_uid
    from .isoform_matrix import cpm_and_ratio_from_counts

    meta = pd.read_csv(groups_file, sep="\t", dtype=str)
    pairs = meta[["uid", "groups"]].apply(lambda s: s.str.strip()).drop_duplicates()
    if pairs.uid.duplicated().any():
        raise ValueError(f"{groups_file}: a uid maps to more than one group")
    uid_to_group = dict(zip(pairs.uid, pairs.groups))
    a = ad.read_h5ad(counts_h5ad)
    split = [n.split(".", 1) for n in np.asarray(a.obs_names, dtype=str)]
    if any(len(s) != 2 for s in split):
        raise ValueError("counts obs names must be <state>.<sample>")
    state = np.array([s[0] for s in split])
    uid = np.array([_resolve_uid(s[1], uid_to_group, "isoform") for s in split])
    keep = np.array([uid_to_group.get(u) in set(groups) for u in uid])
    if not keep.any():
        raise ValueError(f"no pseudobulk matches groups {groups}")
    var = a.var_names.to_numpy(dtype=object)
    X = sparse.csr_matrix(a.X[np.where(keep)[0]])
    del a
    codes, uniq = pd.factorize(pd.Series(list(zip(state[keep], uid[keep]))))
    if len(uniq) != X.shape[0]:
        log(f"[predominance] summing {X.shape[0] - len(uniq)} duplicate state x donor pseudobulks")
        X = (sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
             @ X).tocsr()
    if state_pools:
        original_uniq = list(uniq)
        known_states = set(state)
        extra_matrices, extra_identities = [], []
        for name, members in state_pools.items():
            if (not isinstance(name, str) or not name or name in known_states
                    or not isinstance(members, list) or not members
                    or len(set(members)) != len(members) or not set(members) <= known_states):
                raise ValueError(f"invalid or unresolved state pool {name!r}: {members!r}")
            chosen = [i for i, (s, _) in enumerate(original_uniq) if s in set(members)]
            if not chosen:
                log('[predominance-pool] no requested-group pseudobulks for', name, members)
                continue
            pool_codes, pool_donors = pd.factorize(pd.Series([original_uniq[i][1] for i in chosen]))
            pooling = sparse.csr_matrix((np.ones(len(chosen)), (pool_codes, chosen)),
                                        shape=(len(pool_donors), len(original_uniq)))
            extra_matrices.append(pooling @ X)
            extra_identities.extend((name, donor) for donor in pool_donors)
            log('[predominance-pool]', name, members, 'donors', len(pool_donors))
        if extra_matrices:
            X = sparse.vstack([X, *extra_matrices], format='csr')
            uniq = original_uniq + extra_identities
    pb_group = np.array([uid_to_group[u[1]] for u in uniq], dtype=object)
    _, ratio, E, gene_ids = cpm_and_ratio_from_counts(X, var, cpm_threshold=cpm_threshold)
    sgc, sgu = pd.factorize(pd.Series(list(zip([u[0] for u in uniq], pb_group))))
    Agg = sparse.csr_matrix((np.ones(len(sgc)), (sgc, np.arange(len(sgc)))), shape=(len(sgu), len(sgc)))
    U = (Agg @ ratio).tocoo()
    N = (Agg @ E.astype(np.float64)).tocsr()
    fgene = np.array([str(v).split(":", 1)[0] for v in var], dtype=object)
    d = pd.DataFrame({"sg": U.row, "f": U.col, "usage_sum": U.data})
    d["gene"] = fgene[d.f.to_numpy()]
    if genes is not None:
        d = d[d.gene.isin(set(genes))]
    d["predominant"] = var[d.f.to_numpy()]
    d = d.sort_values(["sg", "gene", "usage_sum", "predominant"], ascending=[True, True, False, True])
    all_usage = None
    if return_usage:
        # Preserve actual counts even where the gene fails native CPM conditioning.
        # A winner from another context may have reads only in such pseudobulks;
        # an absent sparse ratio entry must never turn those reads into zero.
        reads = (Agg @ X).tocoo()
        rgene = fgene[reads.col]
        nd = np.asarray(N[reads.row, pd.Index(gene_ids).get_indexer(rgene)]).ravel()
        sums = np.asarray(U.tocsr()[reads.row, reads.col]).ravel()
        means = np.divide(sums, nd, out=np.full(len(nd), np.nan), where=nd > 0)
        all_usage = pd.DataFrame({"state": [sgu[i][0] for i in reads.row], "group": [sgu[i][1] for i in reads.row],
                                  "gene": rgene, "feature": var[reads.col],
                                  "mean_usage": means, "donors": nd.astype(int),
                                  "assigned_reads": reads.data})
    d = d.drop_duplicates(["sg", "gene"])
    donors = np.asarray(N[d.sg.to_numpy(), pd.Index(gene_ids).get_indexer(d.gene)]).ravel()
    if (donors <= 0).any():
        raise AssertionError("a gene with usage has no expressing donor")
    log(f"[predominance] {counts_h5ad}: {len(uniq):,} state x donor pseudobulks of {'+'.join(groups)}; predominant "
        f"isoform of {d.gene.nunique():,} genes in {len(sgu)} state x group cells")
    pred = pd.DataFrame({"state": [sgu[i][0] for i in d.sg], "group": [sgu[i][1] for i in d.sg], "gene": d.gene.values,
                         "predominant": d.predominant.values, "mean_usage": d.usage_sum.values / donors,
                         "donors": donors.astype(int)})
    return (pred, all_usage) if return_usage else pred


def predominant_context_partitions(pred, usage):
    """Lossless native context partitions for repeated pair queries.

    Original row indices/values are retained. This changes lookup cost only;
    every pair receives exactly the rows selected by its native predicates.
    """
    return ({key: frame for key, frame in pred.groupby(['state', 'group'], sort=False)},
            {key: frame for key, frame in usage.groupby(['state', 'group'], sort=False)})


def predominant_pair_context(pred, usage, state_A, group_A, state_B, group_B, min_margin=.2,
                             context_partitions=None):
    """Whole-gene leader changes, with old/new leader margins in BOTH contexts.

    Preserve the outer gene roster, unavailable contexts, unchanged leaders and ties.
    The gate requires different unique maxima and a >=min_margin usage gap between
    those two leaders in at least one context. This is descriptive, not an FDR call.
    """
    if not 0 <= min_margin <= 1:
        raise ValueError("predominance margin must be a fraction between0and1")
    keys = ["gene", "feature"]
    if context_partitions is None:
        left = pred[(pred.state == state_A) & (pred.group == group_A)].set_index("gene")
        right = pred[(pred.state == state_B) & (pred.group == group_B)].set_index("gene")
        ua = usage[(usage.state == state_A) & (usage.group == group_A)]
        ub = usage[(usage.state == state_B) & (usage.group == group_B)]
    else:
        pp, uu = context_partitions
        left = pp.get((state_A, group_A), pred.iloc[:0]).set_index('gene')
        right = pp.get((state_B, group_B), pred.iloc[:0]).set_index('gene')
        ua = uu.get((state_A, group_A), usage.iloc[:0])
        ub = uu.get((state_B, group_B), usage.iloc[:0])
    if not left.index.is_unique or not right.index.is_unique or ua.duplicated(keys).any() or ub.duplicated(keys).any():
        raise ValueError("duplicate native context identity")
    r = left[['predominant','donors']].rename(columns={'predominant':'predominant_A','donors':'donors_A'}).join(
        right[['predominant','donors']].rename(columns={'predominant':'predominant_B','donors':'donors_B'}), how='outer').sort_index()
    if r.empty:
        return pd.DataFrame()
    has_a, has_b = r.predominant_A.notna(), r.predominant_B.notna()
    r['both_contexts_observed'] = has_a & has_b
    r['state_A'], r['group_A'], r['state_B'], r['group_B'] = state_A, group_A, state_B, group_B
    for label, values, has_context in [('A',ua,has_a),('B',ub,has_b)]:
        maxima = values.groupby('gene').mean_usage.transform('max')
        ties = values[values.mean_usage == maxima].groupby('gene').size()
        r['top_ties_'+label] = ties.reindex(r.index, fill_value=0).astype(int)
        index = values.set_index(keys)
        for source_label, has_model in [('A',has_a),('B',has_b)]:
            ids = pd.MultiIndex.from_arrays([r.index, r['predominant_'+source_label]])
            found = index.reindex(ids)
            available = (has_context & has_model).to_numpy()
            for field, column in [('usage','mean_usage'),('reads','assigned_reads')]:
                # Only a gene with expressing donors can supply a sparse zero.
                r[field+'_'+source_label+'_in_'+label] = np.where(available, found[column].fillna(0).to_numpy(float), np.nan)
    r['leader_changed'] = r.both_contexts_observed & r.predominant_A.ne(r.predominant_B)
    r['unique_leaders'] = r.top_ties_A.eq(1) & r.top_ties_B.eq(1)
    r['margin_A'] = r.usage_A_in_A - r.usage_B_in_A
    r['margin_B'] = r.usage_B_in_B - r.usage_A_in_B
    r['winner_A_usage_change_A_minus_B'] = r.usage_A_in_A - r.usage_A_in_B
    r['winner_B_usage_change_B_minus_A'] = r.usage_B_in_B - r.usage_B_in_A
    r['predominance_switch_at_margin'] = r.leader_changed & r.unique_leaders & r[['margin_A','margin_B']].max(axis=1).ge(min_margin)
    r[['predominant_A','predominant_B']] = r[['predominant_A','predominant_B']].fillna('')
    r[['donors_A','donors_B']] = r[['donors_A','donors_B']].fillna(0).astype(int)
    r.index.name = 'gene'
    return r.reset_index()


def main_predominance(argv):
    """Native descriptive whole-gene leader screen, independent of a differential whitelist."""
    import gzip
    ap = argparse.ArgumentParser(prog="isoform_state_switches predominance",
                                 description="Whole-gene leader identity changes and the old/new leader usage gap; no new statistical test.")
    ap.add_argument('--counts', required=True)
    ap.add_argument('--groups', required=True)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--compare_groups', help='A,B: compare the two groups within each state')
    mode.add_argument('--group', help='Compare every pair of observed states within this group')
    ap.add_argument('--min_margin', type=float, default=.2, help='Old/new leader usage gap in at least one context; absolute fraction, not relative percent')
    ap.add_argument('--min_donors', type=int, default=3, help='Coverage flag only; all rows remain in the tables')
    ap.add_argument('--min_reads', type=float, default=50, help='Native read-depth flag; each leader is measured in its own winning context')
    ap.add_argument('--cpm_threshold', type=float, default=1, help='Existing native gene-expression conditioning')
    ap.add_argument('--protein_classes', help='Also rank predicted Not-NMD coding models; gene-denominator shares remain unchanged')
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args(argv)
    groups = [g.strip() for g in a.compare_groups.split(',')] if a.compare_groups else [a.group]
    if a.compare_groups and (len(groups) != 2 or len(set(groups)) != 2 or not all(groups)):
        ap.error('--compare_groups requires two distinct nonempty groups')
    if not 0 <= a.min_margin <= 1 or a.min_donors < 1 or a.min_reads < 0:
        ap.error('invalid margin or depth parameters')
    os.makedirs(a.outdir, exist_ok=True)
    log = _Log(os.path.join(a.outdir, 'predominance_log.txt'))
    pred, usage = predominant_isoforms(a.counts, a.groups, groups, genes=None, log=log,
                                      cpm_threshold=a.cpm_threshold, return_usage=True)
    pred.to_csv(os.path.join(a.outdir, 'all_genes_predominant_isoforms.tsv'), sep='\t', index=False)
    usage.to_csv(os.path.join(a.outdir, 'all_native_observed_context_usages.tsv.gz'), sep='\t', index=False)
    coding_pred = coding_usage = None
    unresolved_genes = set()
    if a.protein_classes:
        classes = load_protein_classes(a.protein_classes)
        unknown = usage.loc[~usage.feature.isin(classes.index)]
        unknown.to_csv(os.path.join(a.outdir, 'unresolved_protein_class_contexts.tsv.gz'), sep='\t', index=False)
        unresolved_genes = set(unknown.gene)
        coding_ids = classes.index[(pd.to_numeric(classes.protein_length, errors='coerce') > 0)
                                   & classes.nmd_status.eq('Not-NMD')]
        # This is a second ranking of the original shares, not a new RNA matrix
        # or a renormalization over coding features. The all-RNA ranking stays.
        coding_usage = usage[usage.feature.isin(coding_ids)]
        # Raw support in the other context must survive even if native CPM
        # conditioning gives that model zero/undefined mean usage there.
        coding_pred = (coding_usage[coding_usage.mean_usage.gt(0)].sort_values(['state','group','gene','mean_usage','feature'],
                      ascending=[True,True,True,False,True]).drop_duplicates(['state','group','gene'])
                       .rename(columns={'feature':'predominant'})[['state','group','gene','predominant','mean_usage','donors']])
        coding_pred.to_csv(os.path.join(a.outdir, 'predicted_coding_predominant_isoforms.tsv'), sep='\t', index=False)
    states = sorted(set(pred.state))
    contexts = [(s, groups[0], s, groups[1]) for s in states] if a.compare_groups else [
        (left, a.group, right, a.group) for left, right in combinations(states, 2)]
    records = []
    partitions = predominant_context_partitions(pred, usage)
    coding_partitions = predominant_context_partitions(coding_pred, coding_usage) if coding_pred is not None else None
    folder = os.path.join(a.outdir, 'context_pairs')
    os.makedirs(folder, exist_ok=True)
    # Two partitioned tables retain every comparison without creating thousands
    # of nearly identical files. Row offsets below define the context partitions.
    path = os.path.join(folder, 'all_context_pairs.tsv.gz')
    cp = os.path.join(folder, 'predicted_coding_all_context_pairs.tsv.gz')
    combined = gzip.open(path, 'wt')
    coding_combined = gzip.open(cp, 'wt') if coding_pred is not None else None
    row_offset = coding_offset = 0
    for i, (sa, ga, sb, gb) in enumerate(contexts, 1):
        result = predominant_pair_context(pred, usage, sa, ga, sb, gb, a.min_margin,
                                          context_partitions=partitions)
        if len(result):
            result['both_donor_N_at_least_min'] = (result.donors_A >= a.min_donors) & (result.donors_B >= a.min_donors)
            result['either_winning_model_at_least_min_reads'] = (result.reads_A_in_A >= a.min_reads) | (result.reads_B_in_B >= a.min_reads)
            result['both_winning_models_at_least_min_reads'] = (result.reads_A_in_A >= a.min_reads) & (result.reads_B_in_B >= a.min_reads)
            result['well_covered_predominance_switch'] = (result.predominance_switch_at_margin
                & result.both_donor_N_at_least_min & result.both_winning_models_at_least_min_reads)
        if len(result):
            result.to_csv(combined, sep='\t', index=False, header=row_offset == 0)
        record = dict(state_A=sa, group_A=ga, state_B=sb, group_B=gb, all_gene_roles=len(result),
                      leader_changes=int(result.leader_changed.sum()) if len(result) else 0,
                      switches_at_margin=int(result.predominance_switch_at_margin.sum()) if len(result) else 0,
                      well_covered_switches=int(result.well_covered_predominance_switch.sum()) if len(result) else 0,
                      output=path, row_start=row_offset, row_count=len(result))
        row_offset += len(result)
        records.append(record)
        if coding_pred is not None:
            cr = predominant_pair_context(coding_pred, coding_usage, sa, ga, sb, gb, a.min_margin,
                                          context_partitions=coding_partitions)
            if len(cr):
                cr['coding_rank_source_coverage_complete'] = ~cr.gene.isin(unresolved_genes)
                cr['both_winning_models_at_least_min_reads'] = (cr.reads_A_in_A >= a.min_reads) & (cr.reads_B_in_B >= a.min_reads)
                cr['both_donor_N_at_least_min'] = (cr.donors_A >= a.min_donors) & (cr.donors_B >= a.min_donors)
                cr['well_covered_predominance_switch'] = (cr.predominance_switch_at_margin
                    & cr.both_winning_models_at_least_min_reads & cr.both_donor_N_at_least_min & cr.coding_rank_source_coverage_complete)
            if len(cr):
                cr.to_csv(coding_combined, sep='\t', index=False, header=coding_offset == 0)
            record.update(predicted_coding_output=cp, predicted_coding_row_start=coding_offset,
                          predicted_coding_row_count=len(cr),
                          predicted_coding_well_covered_switches=int(cr.well_covered_predominance_switch.sum()) if len(cr) else 0)
            coding_offset += len(cr)
        log('[predominance-pair]', sa, ga, sb, gb, record)
    combined.close()
    if coding_combined is not None:
        coding_combined.close()
    output_hash, coding_hash = _sha256(path), _sha256(cp) if coding_combined is not None else None
    for record in records:
        record['sha256'] = output_hash
        if coding_hash is not None:
            record['predicted_coding_sha256'] = coding_hash
    params = dict(tool='altanalyze3.components.long_read.isoform_state_switches predominance', args=vars(a),
                  inputs=dict(counts=dict(path=a.counts,sha256=_sha256(a.counts)),
                              groups=dict(path=a.groups,sha256=_sha256(a.groups)),
                              protein_classes=dict(path=a.protein_classes,sha256=_sha256(a.protein_classes)) if a.protein_classes else None),
                  native_source_sha256=_sha256(__file__),
                  full_gene_denominator_and_all_features_preserved=True, differential_whitelist_used=False,
                  native_usage_method='Mean donor isoform/gene read share, native gene CPM conditioning; equal donor weights',
                  margin_definition='max(usage(A-leader in A)-usage(B-leader in A), usage(B-leader in B)-usage(A-leader in B)) >= min_margin',
                  ties='All ties retained and flagged; a tied maximum is not a unique predominant switch',
                  statistical_test_performed=False, FDR_significance_claim=False,
                  protein_class_unresolved_genes=len(unresolved_genes),
                  source_nonzero_usage_rows=len(usage), genes=int(pred.gene.nunique()), contexts=records)
    with open(os.path.join(a.outdir, 'predominance_parameters.json'), 'w') as handle:
        json.dump(params, handle, indent=2)


def annotate_predominance(table, pred, classes, group_A, group_B):
    """Per differential isoform: its gene's predominant isoform in the group where it is higher and in the other group.
    predominance_switch = it is predominant where it is higher and another isoform is predominant in the other group;
    switch_length compares its protein length with that other isoform: 'longer', 'shorter' or 'same length' (empty when
    there is no switch)."""
    t = table.copy()
    hi = np.where(t.direction == f"higher in {group_A}", group_A, group_B)
    lo = np.where(hi == group_A, group_B, group_A)
    gene = t.feature.str.split(":", n=1).str[0]
    P = pred.set_index(["state", "group", "gene"]).predominant
    t["predominant_where_higher"] = P.reindex(pd.MultiIndex.from_arrays([t.state, hi, gene])).to_numpy()
    t["predominant_other_group"] = P.reindex(pd.MultiIndex.from_arrays([t.state, lo, gene])).to_numpy()
    feat = t.feature.to_numpy()
    sw = (t.predominant_where_higher.to_numpy() == feat) & t.predominant_other_group.notna().to_numpy() & \
         (t.predominant_other_group.to_numpy() != feat)
    oc = classes.reindex(t.predominant_other_group[sw].to_numpy())
    if oc.protein_length.isna().any():
        raise ValueError(f"{int(oc.protein_length.isna().sum())} predominant isoforms are absent from the protein "
                         f"classes table")
    Lx = pd.to_numeric(t.protein_length[sw], errors="coerce").to_numpy(float)
    Lo = pd.to_numeric(oc.protein_length, errors="coerce").to_numpy(float)
    t["predominance_switch"] = sw
    t["other_protein_length"] = np.nan
    t.loc[sw, "other_protein_length"] = Lo
    t["switch_length"] = ""
    t.loc[sw, "switch_length"] = np.where(Lx > Lo, "longer", np.where(Lx < Lo, "shorter", "same length"))
    # every isoform against the other group's predominant isoform, switch or not (gain = longer, loss = shorter)
    has = t.predominant_other_group.notna().to_numpy()
    oa = classes.reindex(t.predominant_other_group[has].to_numpy())
    if oa.protein_length.isna().any():
        raise ValueError(f"{int(oa.protein_length.isna().sum())} predominant isoforms are absent from the protein "
                         f"classes table")
    La = pd.to_numeric(t.protein_length[has], errors="coerce").to_numpy(float)
    Lb = pd.to_numeric(oa.protein_length, errors="coerce").to_numpy(float)
    t["other_group_predominant_length"] = np.nan
    t.loc[has, "other_group_predominant_length"] = Lb
    t["length_vs_other_group"] = ""
    t.loc[has, "length_vs_other_group"] = np.where(La > Lb, "longer", np.where(La < Lb, "shorter", "same length"))
    return t


GAIN_LOSS_COLORS = {"gain": "#2CA02C", "loss": "#FF7F0E", "no length change": "#1F77B4", "NMD": "#D62728"}


def gain_loss_counts(table, group_A, group_B, states):
    """Plotted isoforms (protein coding or NMD) per state, direction and category: NMD, or for coding isoforms 'gain'
    (longer protein than the gene's predominant isoform in the other group), 'loss' (shorter), 'no length change'
    (same length, or no predominant isoform of the gene in the other group)."""
    s = table[table.outcome.isin(["protein coding", "NMD"])]
    cat = np.where(s.outcome == "NMD", "NMD",
                   np.where(s.length_vs_other_group == "longer", "gain",
                            np.where(s.length_vs_other_group == "shorter", "loss", "no length change")))
    c = s.assign(category=cat).groupby(["state", "direction", "category"]).size().rename("isoforms")
    full = pd.MultiIndex.from_product([states, [f"higher in {group_A}", f"higher in {group_B}"], list(GAIN_LOSS_COLORS)],
                                      names=["state", "direction", "category"])
    c = c.reindex(full, fill_value=0).reset_index()
    if int(c.isoforms.sum()) != int(s.state.isin(states).sum()):
        raise AssertionError("gain/loss counts do not add up to the isoform table")
    return c


def plot_gain_loss(gl, outfile, group_A, group_B, title=None):
    """Mirrored bars per state (isoforms higher in group_A above zero, group_B below, same scale): protein gain and loss
    against the gene's predominant isoform in the other group stacked next to zero, then coding isoforms with no length
    change, then NMD. Written without clip groups or background boxes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "axes.linewidth": 0.5})
    states = list(dict.fromkeys(gl.state))
    n = len(states)
    fig, ax = plt.subplots(figsize=(max(5.0, 0.13 * n + 1.6), 3.2))
    xs = np.arange(n)
    top, tot = 0.0, {}
    for direction, sign in ((f"higher in {group_A}", 1), (f"higher in {group_B}", -1)):
        base = np.zeros(n)
        for cat, color in GAIN_LOSS_COLORS.items():
            v = gl[(gl.direction == direction) & (gl.category == cat)].set_index("state").isoforms
            v = v.reindex(states).fillna(0).to_numpy(float)
            tot[(direction, cat)] = int(v.sum())
            nz = v > 0
            ax.bar(xs[nz], sign * v[nz], bottom=sign * base[nz], width=0.8, linewidth=0, color=color)
            base += v
        top = max(top, base.max() if n else 0)
    lim = top * 1.08 if top > 0 else 1
    ax.set_ylim(-lim, lim)
    ax.plot([-0.6, n - 0.4], [0, 0], color="#333333", lw=0.6)
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_xticks(xs)
    ax.set_xticklabels(states, rotation=90, fontsize=min(7, max(4, int(400 / max(n, 1)))))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda y, _: f"{abs(int(round(y))):,}"))
    ax.tick_params(axis="y", labelsize=6)
    ax.set_ylabel(f"Isoforms with an altered protein\nhigher in {group_B} (below)   higher in {group_A} (above)", fontsize=7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    label = {"gain": "longer protein than the other group's predominant isoform",
             "loss": "shorter protein than the other group's predominant isoform",
             "no length change": "protein coding, same length (or no predominant isoform in the other group)",
             "NMD": "NMD"}
    handles = [matplotlib.patches.Patch(color=c, label=f"{label[k]} ({tot[(f'higher in {group_A}', k)]:,} above / "
                                                       f"{tot[(f'higher in {group_B}', k)]:,} below)")
               for k, c in GAIN_LOSS_COLORS.items()]
    ax.legend(handles=handles, fontsize=6, frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=1,
              title=title, title_fontsize=6, alignment="left")
    fig.patch.set_visible(False)
    ax.patch.set_visible(False)
    for art in ax.get_children():
        if hasattr(art, "set_clip_on"):
            art.set_clip_on(False)
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


def predominance_counts(table, group_A, group_B, states):
    """Plotted isoforms (protein coding or NMD) that are predominance switches, per state, direction, outcome and
    switch_length, with longer_protein_in naming the group whose predominant isoform has the longer protein ('' for
    same length)."""
    s = table[table.outcome.isin(["protein coding", "NMD"]) & table.predominance_switch]
    c = s.groupby(["state", "direction", "outcome", "switch_length"]).size().rename("isoforms")
    full = pd.MultiIndex.from_product([states, [f"higher in {group_A}", f"higher in {group_B}"], ["protein coding", "NMD"],
                                       ["longer", "shorter", "same length"]],
                                      names=["state", "direction", "outcome", "switch_length"])
    c = c.reindex(full, fill_value=0).reset_index()
    hi = np.where(c.direction == f"higher in {group_A}", group_A, group_B)
    lo = np.where(hi == group_A, group_B, group_A)
    c["longer_protein_in"] = np.where(c.switch_length == "longer", hi, np.where(c.switch_length == "shorter", lo, ""))
    if int(c.isoforms.sum()) != len(s[s.state.isin(states)]):
        raise AssertionError("predominance counts do not add up to the isoform table")
    return c


def protein_outcome_counts(table, group_A, group_B, cluster_order=None, states=None):
    """Isoform counts per state, direction and outcome; 'same length' and 'no protein' are counted but not plotted.
    states, when given, is the x axis (states without isoforms get zero rows)."""
    states = list(dict.fromkeys(table.state)) if states is None else list(states)
    if cluster_order:
        states = [s for s in cluster_order if s in states] + sorted(s for s in states if s not in cluster_order)
    c = table.groupby(["state", "direction", "outcome"]).size().rename("isoforms").reset_index()
    full = pd.MultiIndex.from_product([states, [f"higher in {group_A}", f"higher in {group_B}"],
                                       ["protein coding", "NMD", "same length", "no protein"]],
                                      names=["state", "direction", "outcome"])
    c = c.set_index(["state", "direction", "outcome"]).reindex(full, fill_value=0).reset_index()
    if int(c.isoforms.sum()) != int(table.state.isin(states).sum()):
        raise AssertionError("protein-outcome counts do not add up to the isoform table")
    return c


LONGER_IN_COLORS = ("#2CA02C", "#FF7F0E")   # longer protein in group A, in group B


def plot_protein_outcome(counts, outfile, group_A, group_B, title=None, switch_counts=None):
    """Mirrored bars per state: isoforms higher in group_A above zero, higher in group_B below zero on the same
    scale, each stacked protein coding (altered length) then NMD. Isoforms of the gene's longest length and isoforms
    without a protein are left out of the bars (their counts stay in the TSV). With switch_counts (predominance_counts),
    the predominance switches with a longer or shorter protein are taken out of the coding and NMD segments and stacked
    next to zero as 'longer protein in A' and 'longer protein in B' (bar heights unchanged); the PDF is then written
    without clip groups or background boxes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "axes.linewidth": 0.5})
    states = list(dict.fromkeys(counts.state))
    n = len(states)
    fig, ax = plt.subplots(figsize=(max(5.0, 0.13 * n + 1.6), 3.2))
    xs = np.arange(n)
    top = 0.0
    sw_tot = {}
    for direction, sign in ((f"higher in {group_A}", 1), (f"higher in {group_B}", -1)):
        base = np.zeros(n)
        segments = []
        taken = {o: np.zeros(n) for o in ("protein coding", "NMD")}
        if switch_counts is not None:
            sd = switch_counts[switch_counts.direction == direction]
            for g, color in zip((group_A, group_B), LONGER_IN_COLORS):
                part = sd[sd.longer_protein_in == g]
                for o in taken:
                    taken[o] += part[part.outcome == o].groupby("state").isoforms.sum().reindex(states).fillna(0).to_numpy(float)
                v = part.groupby("state").isoforms.sum().reindex(states).fillna(0).to_numpy(float)
                sw_tot[(direction, g)] = int(v.sum())
                segments.append((v, color))
        for outcome in ("protein coding", "NMD"):
            v = counts[(counts.direction == direction) & (counts.outcome == outcome)].set_index("state").isoforms
            v = v.reindex(states).to_numpy(float) - taken[outcome]
            if (v < 0).any():
                raise AssertionError("predominance switches exceed their outcome count")
            segments.append((v, OUTCOME_COLORS[outcome]))
        for v, color in segments:
            nz = v > 0
            ax.bar(xs[nz], sign * v[nz], bottom=sign * base[nz], width=0.8, linewidth=0, color=color)
            base += v
        top = max(top, base.max() if n else 0)
    lim = top * 1.08 if top > 0 else 1
    ax.set_ylim(-lim, lim)
    ax.plot([-0.6, n - 0.4], [0, 0], color="#333333", lw=0.6)
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_xticks(xs)
    ax.set_xticklabels(states, rotation=90, fontsize=min(7, max(4, int(400 / max(n, 1)))))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda y, _: f"{abs(int(round(y))):,}"))
    ax.tick_params(axis="y", labelsize=6)
    ax.set_ylabel(f"Isoforms with an altered protein\nhigher in {group_B} (below)   higher in {group_A} (above)", fontsize=7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    tot = {d: int(counts[(counts.direction == d) & counts.outcome.isin(["protein coding", "NMD"])].isoforms.sum())
           for d in (f"higher in {group_A}", f"higher in {group_B}")}
    handles = [matplotlib.patches.Patch(color=OUTCOME_COLORS[o], label=o) for o in ("protein coding", "NMD")]
    if switch_counts is not None:
        handles = [matplotlib.patches.Patch(
            color=c, label=f"predominance switch, longer protein in {g} "
                           f"({sw_tot[(f'higher in {group_A}', g)]:,} above / {sw_tot[(f'higher in {group_B}', g)]:,} below)")
            for g, c in zip((group_A, group_B), LONGER_IN_COLORS)] + handles
    ax.legend(handles=handles, fontsize=6, frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2,
              title=(title + "; " if title else "") + f"{tot[f'higher in {group_A}']:,} higher in {group_A}, "
                    f"{tot[f'higher in {group_B}']:,} higher in {group_B}", title_fontsize=6, alignment="left")
    if switch_counts is not None:
        fig.patch.set_visible(False)
        ax.patch.set_visible(False)
        for art in ax.get_children():
            if hasattr(art, "set_clip_on"):
                art.set_clip_on(False)
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


def main_protein_outcome(argv):
    """protein-outcome: plot the protein outcome of age-differential isoforms (or of called switches) per state."""
    ap = argparse.ArgumentParser(
        prog="isoform_state_switches protein-outcome",
        description="Mirrored per-state bars of isoforms with an altered protein: higher in group A above zero, higher "
                    "in group B below. An isoform counts when its predicted protein differs in length from the gene's "
                    "longest protein; it is coloured NMD or protein coding. Same-length and no-protein isoforms are "
                    "counted in the TSV and left out of the bars.")
    ap.add_argument("--compare_groups", required=True, help="A,B (e.g. aged,young); A is drawn above zero")
    ap.add_argument("--protein_classes", required=True,
                    help="ALL_catalog_protein_classes.tsv (protein_length, nmd_status, longest_gene_protein)")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--diff_dir", help="folder of <A>-<B>-<state>_stats-annotated.txt covariate tables")
    src.add_argument("--switches", help="isoform_state_switches.tsv of a --compare_groups run")
    ap.add_argument("--diff_pattern", default="*_stats-annotated.txt")
    ap.add_argument("--fold", type=float, default=1.5, help="--diff_dir: fold threshold")
    ap.add_argument("--fold_scale", choices=["log2", "ratio"], default="log2",
                    help="--diff_dir: log2 = |DiffMeans| > log2(fold) (isoform expression tables); ratio = "
                         "Mean_A / Mean_B or its inverse > fold (isoform usage tables)")
    ap.add_argument("--rawp", type=float, default=0.05, help="--diff_dir: ebayesPval < this")
    ap.add_argument("--min_usage_shift", type=float, default=None,
                    help="With --fold_scale ratio: use absolute fractional usage difference > this instead of a fold-ratio gate")
    ap.add_argument("--min_donors", type=int, default=3, help="--diff_dir: donors with a value in each group")
    ap.add_argument("--cluster_order", default=None, help="one state per line; x-axis order")
    ap.add_argument("--states", default=None,
                    help="one state per line: keep only these states (all of them on the x axis, zero when empty)")
    ap.add_argument("--counts", default=None,
                    help="cohort isoform_combined_pseudo_cluster_counts.h5ad; with --groups, marks predominance "
                         "switches (predominant isoform where higher, another isoform predominant in the other group) "
                         "and whether the isoform's protein is longer or shorter than that other isoform")
    ap.add_argument("--groups", default=None, help="metadata with uid and groups columns (with --counts)")
    ap.add_argument("--title", default=None)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    if a.min_usage_shift is not None and not a.diff_dir:
        ap.error('--min_usage_shift requires --diff_dir with fractional usage tables')
    if bool(a.counts) != bool(a.groups):
        ap.error("--counts and --groups go together")
    gA, gB = [g.strip() for g in a.compare_groups.split(",")]
    os.makedirs(a.outdir, exist_ok=True)
    log = _Log(os.path.join(a.outdir, "protein_outcome_log.txt"))
    classes = protein_outcome_classes(a.protein_classes)
    order = [s.strip() for s in open(a.cluster_order) if s.strip()] if a.cluster_order else None
    if a.diff_dir:
        table, manifest = protein_outcome_from_tables(a.diff_dir, gA, gB, classes, fold=a.fold, rawp=a.rawp,
                                                      min_donors=a.min_donors, fold_scale=a.fold_scale,
                                                      pattern=a.diff_pattern, log=log,
                                                      min_usage_shift=a.min_usage_shift)
    else:
        table, manifest = protein_outcome_from_switches(a.switches, classes, gA, gB), []
    keep_states = None
    if a.states:
        keep_states = [s.strip() for s in open(a.states) if s.strip()]
        absent = [s for s in keep_states if s not in set(table.state)]
        n0 = len(table)
        table = table[table.state.isin(keep_states)].reset_index(drop=True)
        log(f"[protein-outcome] --states: {len(keep_states)} states kept, {n0 - len(table):,} of {n0:,} isoform rows "
            f"dropped; states without differential isoforms: {absent if absent else 'none'}")
    counts = protein_outcome_counts(table, gA, gB, cluster_order=order, states=keep_states)
    sw_counts = None
    if a.counts:
        genes = sorted(set(table.feature.str.split(":", n=1).str[0]))
        pred = predominant_isoforms(a.counts, a.groups, [gA, gB], genes=genes, log=log)
        table = annotate_predominance(table, pred, classes, gA, gB)
        sw_counts = predominance_counts(table, gA, gB, list(dict.fromkeys(counts.state)))
        gl = gain_loss_counts(table, gA, gB, list(dict.fromkeys(counts.state)))
        pred.to_csv(os.path.join(a.outdir, "protein_outcome_predominant_isoforms.tsv"), sep="\t", index=False)
        sw_counts.to_csv(os.path.join(a.outdir, "protein_outcome_predominance_counts.tsv"), sep="\t", index=False)
        gl.to_csv(os.path.join(a.outdir, "protein_outcome_gain_loss_counts.tsv"), sep="\t", index=False)
        plot_gain_loss(gl, os.path.join(a.outdir, "protein_outcome_gain_loss_by_state.pdf"), gA, gB, title=a.title)
        for (d, k), v in gl.groupby(["direction", "category"]).isoforms.sum().items():
            log(f"[protein-outcome] gain/loss vs the other group's predominant isoform, {d}, {k}: {int(v):,}")
    table.to_csv(os.path.join(a.outdir, "protein_outcome_isoforms.tsv"), sep="\t", index=False)
    counts.to_csv(os.path.join(a.outdir, "protein_outcome_counts.tsv"), sep="\t", index=False)
    plot_protein_outcome(counts, os.path.join(a.outdir, "protein_outcome_by_state.pdf"), gA, gB, title=a.title,
                         switch_counts=sw_counts)
    summ = counts.groupby(["direction", "outcome"]).isoforms.sum()
    extra = {}
    if sw_counts is not None:
        extra["predominance_switch_totals"] = {f"{d}|{l}": int(v) for (d, l), v in
                                               sw_counts.groupby(["direction", "switch_length"]).isoforms.sum().items()}
    with open(os.path.join(a.outdir, "protein_outcome_parameters.json"), "w") as fh:
        json.dump(dict(tool="altanalyze3.components.long_read.isoform_state_switches protein-outcome",
                       argv=sys.argv, args=vars(a), isoforms=len(table),
                       totals={f"{d}|{o}": int(v) for (d, o), v in summ.items()}, tables=manifest, **extra), fh, indent=2)
    for (d, o), v in summ.items():
        log(f"[protein-outcome] {d}, {o}: {int(v):,}")
    for k, v in extra.get("predominance_switch_totals", {}).items():
        log(f"[protein-outcome] predominance switch, {k}: {v:,}")


def load_group_usage(counts_h5ad, groups_file, group, features, log=print, cpm_threshold=1):
    """Donor x state usage, isoform reads and gene reads for `features`, from the cohort counts h5ad.
    Usage comes from isoform_matrix.cpm_and_ratio_from_counts over ALL features of the group's pseudobulks.
    `group` is one group name, or a list of names whose pseudobulks are all kept (labelled in the returned 'group')."""
    import anndata as ad
    from .comparisons import _resolve_uid
    from .isoform_matrix import cpm_and_ratio_from_counts

    meta = pd.read_csv(groups_file, sep="\t", dtype=str)
    pairs = meta[["uid", "groups"]].apply(lambda s: s.str.strip()).drop_duplicates()
    if pairs.uid.duplicated().any():
        raise ValueError(f"{groups_file}: a uid maps to more than one group")
    uid_to_group = dict(zip(pairs.uid, pairs.groups))
    wanted = {group} if isinstance(group, str) else set(group)
    a = ad.read_h5ad(counts_h5ad)
    names = np.asarray(a.obs_names, dtype=str)
    split = [n.split(".", 1) for n in names]
    if any(len(s) != 2 for s in split):
        raise ValueError("counts obs names must be <state>.<sample>")
    state = np.array([s[0] for s in split])
    uid = np.array([_resolve_uid(s[1], uid_to_group, "isoform") for s in split])
    keep = np.array([uid_to_group.get(u) in wanted for u in uid])
    label = group if isinstance(group, str) else "+".join(group)
    log(f"[counts] {counts_h5ad}: {a.n_obs:,} pseudobulks x {a.n_vars:,} features; {label}: {int(keep.sum()):,} "
        f"pseudobulks, {len(set(uid[keep]))} donors, {len(set(state[keep]))} states")
    if not keep.any():
        raise ValueError(f"no pseudobulk matches group {group}")
    var = a.var_names.to_numpy(dtype=object)
    fpos = pd.Index(var).get_indexer(features)
    if (fpos < 0).any():
        raise ValueError(f"{int((fpos < 0).sum())} candidate features absent from {counts_h5ad}")
    fgene = np.array([f.split(":", 1)[0] for f in features])
    cand_genes, ginv = np.unique(fgene, return_inverse=True)
    in_cand = pd.Index(cand_genes).get_indexer([v.split(":", 1)[0] for v in var])
    sel = np.where(in_cand >= 0)[0]
    Gm = sparse.csr_matrix((np.ones(len(sel)), (sel, in_cand[sel])), shape=(len(var), len(cand_genes)))
    Xall = sparse.csr_matrix(a.X)
    cohort_states = np.array(sorted(set(state)))
    Sm = sparse.csr_matrix((np.ones(len(state)), (pd.Index(cohort_states).get_indexer(state), np.arange(len(state)))),
                           shape=(len(cohort_states), len(state)))
    cohort_reads = (Sm @ Xall[:, fpos]).toarray()
    cohort_gene_reads = (Sm @ (Xall @ Gm)).toarray()[:, ginv]
    del Xall
    X = sparse.csr_matrix(a.X[np.where(keep)[0]])
    del a
    key = pd.Series(list(zip(state[keep], uid[keep])))
    codes, uniq = pd.factorize(key)
    if len(uniq) != X.shape[0]:
        log(f"[counts] summing {X.shape[0] - len(uniq)} duplicate state x donor pseudobulks")
        S = sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
        X = (S @ X).tocsr()
    pb_state = np.array([u[0] for u in uniq])
    pb_uid = np.array([u[1] for u in uniq])

    cpm, ratio, E, genes = cpm_and_ratio_from_counts(X, var, cpm_threshold=cpm_threshold)
    gpos = pd.Index(genes).get_indexer(fgene)
    if (gpos < 0).any():
        raise ValueError("candidate gene missing from the ratio gene index")
    usage = ratio[:, fpos].toarray().astype(float)
    expressed = E[:, gpos].toarray().astype(bool)
    usage[~expressed] = np.nan
    reads = X[:, fpos].toarray()
    gene_reads = (X @ Gm).toarray()[:, ginv]
    return dict(state=pb_state, uid=pb_uid, group=np.array([uid_to_group[u] for u in pb_uid]), usage=usage,
                reads=reads, gene_reads=gene_reads,
                cohort_states=cohort_states, cohort_reads=cohort_reads, cohort_gene_reads=cohort_gene_reads)


# ---------------------------------------------------------------- pairwise tests
def state_arrays(data, states, donors):
    """(states x donors x features) usage with NaN where the donor has no pseudobulk, plus pooled reads."""
    F = data["usage"].shape[1]
    U = np.full((len(states), len(donors), F), np.nan)
    R = np.zeros((len(states), F))
    GR = np.zeros((len(states), F))
    si = {s: i for i, s in enumerate(states)}
    di = {d: i for i, d in enumerate(donors)}
    for r, (s, d) in enumerate(zip(data["state"], data["uid"])):
        if s not in si:
            continue
        U[si[s], di[d]] = data["usage"][r]
        R[si[s]] += data["reads"][r]
        GR[si[s]] += data["gene_reads"][r]
    return U, R, GR


def test_state_pair(UA, UB, DA, DB, min_reads, min_donors):
    """Paired moderated t for every isoform whose read depth (DA, DB) exceeds min_reads in A or B, with >= min_donors
    paired donors; isoforms with 1 to min_donors - 1 paired donors are reported without a test.
    Returns (column index, dict, prior)."""
    if min_donors < 2:
        raise ValueError("a paired test needs min_donors >= 2")
    D = UA - UB
    n = np.sum(~np.isnan(D), axis=0)
    gate = ((DA > min_reads) | (DB > min_reads)) & (n >= 1)
    cols = np.where(gate)[0]
    if cols.size == 0:
        return cols, None, None
    Dg, ng = D[:, cols], n[cols].astype(float)
    pair = ~np.isnan(Dg)
    mean = np.nanmean(Dg, axis=0)
    s2 = np.full(cols.size, np.nan)
    two = ng >= 2
    s2[two] = np.nanvar(Dg[:, two], axis=0, ddof=1)
    tested = ng >= min_donors
    t, df_total, p = (np.full(cols.size, np.nan) for _ in range(3))
    s0, d0, info = fit_f_dist(s2[tested], ng[tested] - 1)
    if tested.any():
        t[tested], df_total[tested], p[tested], _ = moderated_one_sample_t(mean[tested], s2[tested], ng[tested],
                                                                           s0, d0)
    mA = np.nanmean(np.where(pair, UA[:, cols], np.nan), axis=0)
    mB = np.nanmean(np.where(pair, UB[:, cols], np.nan), axis=0)
    res = dict(n_donors=ng.astype(int), tested=tested, mean_usage_A=mA, mean_usage_B=mB, delta_usage=mean,
               sd_delta=np.sqrt(s2), t_moderated=t, df_total=df_total, p=p)
    info.update(s0_squared=s0, d0=d0, reported=int(cols.size), tested=int(tested.sum()),
                nonfinite_t=int((~np.isfinite(t[tested])).sum()), zero_variance_tested=int((s2[tested] == 0).sum()))
    return cols, res, info


def test_groups_in_state(UA, UB, DA, DB, min_reads, min_donors):
    """Moderated two-sample t (group A donors vs group B donors) for every isoform whose read depth (DA, DB) exceeds
    min_reads in A or B and that has >= 1 donor with a usage in each group; tested when both groups have
    >= min_donors such donors. Returns (column index, dict, prior)."""
    if min_donors < 2:
        raise ValueError("a two-sample test needs min_donors >= 2")
    nA = np.sum(~np.isnan(UA), axis=0)
    nB = np.sum(~np.isnan(UB), axis=0)
    gate = ((DA > min_reads) | (DB > min_reads)) & (nA >= 1) & (nB >= 1)
    cols = np.where(gate)[0]
    if cols.size == 0:
        return cols, None, None
    A, B = UA[:, cols], UB[:, cols]
    na, nb = nA[cols].astype(float), nB[cols].astype(float)
    mA, mB = np.nanmean(A, axis=0), np.nanmean(B, axis=0)
    vA, vB = np.full(cols.size, np.nan), np.full(cols.size, np.nan)
    vA[na >= 2] = np.nanvar(A[:, na >= 2], axis=0, ddof=1)
    vB[nb >= 2] = np.nanvar(B[:, nb >= 2], axis=0, ddof=1)
    tested = (na >= min_donors) & (nb >= min_donors)
    s2 = np.full(cols.size, np.nan)
    s2[tested] = ((na[tested] - 1) * vA[tested] + (nb[tested] - 1) * vB[tested]) / (na[tested] + nb[tested] - 2)
    t, df_total, p = (np.full(cols.size, np.nan) for _ in range(3))
    s0, d0, info = fit_f_dist(s2[tested], na[tested] + nb[tested] - 2)
    if tested.any():
        t[tested], df_total[tested], p[tested], _ = moderated_two_sample_t(mA[tested] - mB[tested], s2[tested],
                                                                           na[tested], nb[tested], s0, d0)
    res = dict(n_donors_A=na.astype(int), n_donors_B=nb.astype(int), tested=tested, mean_usage_A=mA,
               mean_usage_B=mB, delta_usage=mA - mB, sd_A=np.sqrt(vA), sd_B=np.sqrt(vB), sd_pooled=np.sqrt(s2),
               t_moderated=t, df_total=df_total, p=p)
    info.update(s0_squared=s0, d0=d0, reported=int(cols.size), tested=int(tested.sum()),
                undefined_t=int((~np.isfinite(t[tested])).sum()), zero_variance_tested=int((s2[tested] == 0).sum()))
    return cols, res, info


def check_required_switches(sw, required, outdir, log=print):
    """Each (feature higher in A, feature higher in B) must be called in >= 1 row of sw. Writes
    positive_controls.tsv and returns the names of missing controls."""
    rows, missing = [], []
    for fa, fb in required:
        hit = sw[(sw.isoform_higher_in_A == fa) & (sw.isoform_higher_in_B == fb)]
        where = hit.state.tolist() if "state" in hit else [f"{a}|{b}" for a, b in zip(hit.state_A, hit.state_B)]
        rows.append(dict(feature_higher_in_A=fa, feature_higher_in_B=fb, called=len(where) > 0, n_called=len(where),
                         called_in=";".join(where)))
        if not where:
            missing.append(f"{fa} > {fb}")
        log(f"[controls] {fa} higher in A, {fb} higher in B: {'PASS' if where else 'FAIL'} "
            f"({len(where)} call(s){': ' + ', '.join(where) if where else ''})")
    pd.DataFrame(rows).to_csv(os.path.join(outdir, "positive_controls.tsv"), sep="\t", index=False)
    return missing


def run_groups(counts, groups, group_A, group_B, diff_dir, protein_classes, outdir, fold=1.5, adjp=0.05,
               rawp=None, event_dir=None, event_dpsi=0.2, event_rawp=0.05, min_reads=50, min_delta=0.2, alpha=0.05,
               min_donors=3, states=None, cluster_order=None, diff_pattern="*_stats-annotated.txt",
               read_rule="isoform", coding_only=True, required=None):
    """Isoform switches between two donor groups within each cell state (see the module docstring)."""
    if read_rule not in ("isoform", "gene"):
        raise ValueError("read_rule must be isoform or gene")
    if group_A == group_B:
        raise ValueError("--compare_groups needs two different groups")
    t0 = time.time()
    os.makedirs(os.path.join(outdir, "logs"), exist_ok=True)
    log = _Log(os.path.join(outdir, "logs", "run_log.txt"))
    params = dict(tool="altanalyze3.components.long_read.isoform_state_switches", mode="compare_groups",
                  argv=sys.argv, counts=os.path.abspath(counts), groups=os.path.abspath(groups),
                  group_A=group_A, group_B=group_B, delta="mean usage group_A - mean usage group_B",
                  diff_dir=os.path.abspath(diff_dir), diff_pattern=diff_pattern,
                  whitelist_rule=(f"|DiffMeans| > log2({fold}) and "
                                  + (f"ebayesPval < {rawp}" if rawp is not None else f"ebayesAdjPval < {adjp}")),
                  event_dir=os.path.abspath(event_dir) if event_dir else None,
                  event_rule=f"|dPSI| > {event_dpsi} and ebayesPval < {event_rawp}" if event_dir else None,
                  protein_classes=os.path.abspath(protein_classes),
                  protein_coding_rule="protein_length > 0 and nmd_status == Not-NMD", coding_only=coding_only,
                  outdir=os.path.abspath(outdir), fold=fold, log2_fold_cut=float(np.log2(fold)), adjp=adjp, rawp=rawp,
                  event_dpsi=event_dpsi, event_rawp=event_rawp, min_reads=min_reads, read_rule=read_rule,
                  read_scope="group", min_delta=min_delta, alpha=alpha, min_donors=min_donors, states=states,
                  cluster_order=os.path.abspath(cluster_order) if cluster_order else None,
                  required_switches=[list(r) for r in (required or [])],
                  module_sha256=_sha256(os.path.abspath(__file__)), python=platform.python_version(),
                  numpy=np.__version__, pandas=pd.__version__, scipy=sys.modules["scipy"].__version__,
                  host=platform.node(), started=time.strftime("%Y-%m-%d %H:%M:%S"))
    log(json.dumps(params, indent=2))

    wl, manifest = load_covariate_whitelist(diff_dir, group_A, group_B, fold, adjp, rawp, diff_pattern, log=log)
    log(f"[whitelist] {len(wl):,} isoforms pass the isoform rule in >= 1 of {len(manifest)} states "
        f"[{time.time() - t0:.0f}s]")
    pc = load_protein_classes(protein_classes)
    ev_manifest = []
    if event_dir:
        iso_to_genes = defaultdict(set)
        for g, i in zip(pc["gene"], pc["isoform"]):
            iso_to_genes[i].add(g)
        ev, ev_manifest = load_event_whitelist(event_dir, group_A, group_B, event_dpsi, event_rawp, diff_pattern,
                                               log=log, iso_to_genes=iso_to_genes)
        del iso_to_genes
        n_iso = len(wl)
        wl = wl.merge(ev, on="feature", how="outer", indicator=True)
        wl["whitelist_source"] = wl["_merge"].map({"left_only": "isoform", "right_only": "event",
                                                   "both": "isoform+event"}).astype(str)
        wl = wl.drop(columns="_merge")
        wl["gene_symbol"] = wl.gene_symbol.fillna(wl.gene_symbol_event).fillna("")
        wl = wl.drop(columns="gene_symbol_event")
        for c in ("de_up_states", "de_down_states", "event_states"):
            wl[c] = wl[c].fillna("")
        wl["n_de_states"] = wl.n_de_states.fillna(0).astype(int)
        wl = wl.sort_values("feature").reset_index(drop=True)
        log(f"[whitelist] events add {int(wl.whitelist_source.eq('event').sum()):,} isoforms; union {len(wl):,} "
            f"(isoform rule {n_iso:,}; both {int(wl.whitelist_source.eq('isoform+event').sum()):,})")
    else:
        wl["event_states"] = ""
        wl["whitelist_source"] = "isoform"

    wl["protein_class"] = wl.feature.map(pc["class"])
    wl["protein_length"] = wl.feature.map(pc["protein_length"])
    wl["nmd_status"] = wl.feature.map(pc["nmd_status"])
    miss = int(wl.protein_class.isna().sum())
    log(f"[protein] whitelist isoforms without a protein class row: {miss:,} of {len(wl):,}")
    log(f"[protein] class counts in whitelist: {wl.protein_class.value_counts(dropna=False).to_dict()}")
    is_coding = (pd.to_numeric(wl.protein_length, errors="coerce") > 0) & (wl.nmd_status == "Not-NMD")
    pool = wl[is_coding].copy() if coding_only else wl.copy()
    pool["gene_id"] = pool.feature.str.split(":", n=1).str[0]
    per_gene = pool.groupby("gene_id").feature.transform("size")
    cand = pool[per_gene >= 2].sort_values(["gene_id", "feature"]).reset_index(drop=True)
    log(f"[candidates] protein length > 0 and Not-NMD: {int(is_coding.sum()):,} of {len(wl):,} whitelist isoforms; "
        f"pool {'protein coding' if coding_only else 'every whitelisted isoform (--all_isoforms)'}: {len(pool):,}; "
        f"in genes with >= 2: {len(cand):,} isoforms, {cand.gene_id.nunique():,} genes")
    if cand.empty:
        raise ValueError("no candidate isoforms")
    wl["protein_coding"] = is_coding.to_numpy()
    wl["candidate"] = wl.feature.isin(cand.feature)
    wl["excluded_by"] = np.select([~wl.protein_coding & coding_only, ~wl.candidate],
                                  ["not protein coding (length 0 or NMD)",
                                   f"only {'protein-coding ' if coding_only else ''}whitelisted isoform of its gene"],
                                  "")
    for fa, fb in (required or []):
        for f in (fa, fb):
            log(f"[controls] {f}: whitelist {f in set(wl.feature)}, candidate {f in set(cand.feature)}")

    data = load_group_usage(counts, groups, [group_A, group_B], cand.feature.tolist(), log=log)
    present = sorted(set(data["state"]))
    if states is not None:
        unknown = sorted(set(states) - set(present))
        if unknown:
            raise ValueError(f"requested states absent from the {group_A}/{group_B} pseudobulks: {unknown}")
    order = [s.strip() for s in open(cluster_order) if s.strip()] if cluster_order else []
    ranked = [s for s in order if s in present] + sorted(s for s in present if s not in order)
    keep_states = [s for s in ranked if states is None or s in states]

    def side(g):
        m = data["group"] == g
        sub = dict(state=data["state"][m], uid=data["uid"][m], usage=data["usage"][m], reads=data["reads"][m],
                   gene_reads=data["gene_reads"][m])
        donors = sorted(set(sub["uid"]))
        n = pd.Series(sub["state"]).value_counts()
        U, R, GR = state_arrays(sub, keep_states, donors)
        return donors, n, U, R, GR

    donors_A, n_A, UA, RA, GRA = side(group_A)
    donors_B, n_B, UB, RB, GRB = side(group_B)
    nA_state = {s: int(n_A.get(s, 0)) for s in keep_states}
    nB_state = {s: int(n_B.get(s, 0)) for s in keep_states}
    few = [f"{s}({nA_state[s]}/{nB_state[s]})" for s in keep_states
           if nA_state[s] < min_donors or nB_state[s] < min_donors]
    log(f"[states] {len(keep_states)} states kept; {group_A} donors {len(donors_A)}, {group_B} donors "
        f"{len(donors_B)}; states with < {min_donors} donors in a group (reported, never tested; "
        f"{group_A}/{group_B}): {', '.join(few) or 'none'}")
    ci = pd.Index(data["cohort_states"]).get_indexer(keep_states)
    CR, CGR = data["cohort_reads"][ci], data["cohort_gene_reads"][ci]

    blocks, priors = [], []
    for k, s in enumerate(keep_states):
        DA, DB = (RA[k], RB[k]) if read_rule == "isoform" else (GRA[k], GRB[k])
        cols, res, info = test_groups_in_state(UA[k], UB[k], DA, DB, min_reads, min_donors)
        prior = dict(state=s, donors_A=nA_state[s], donors_B=nB_state[s], reported=0, tested=0)
        if res is None:
            priors.append(prior)
            continue
        prior.update(info)
        priors.append(prior)
        with np.errstate(divide="ignore", invalid="ignore"):
            res.update(ik=np.full(cols.size, k), col=cols, reads_A=RA[k, cols], reads_B=RB[k, cols],
                       gene_reads_A=GRA[k, cols], gene_reads_B=GRB[k, cols], cohort_reads=CR[k, cols],
                       cohort_gene_reads=CGR[k, cols], pooled_usage_A=RA[k, cols] / GRA[k, cols],
                       pooled_usage_B=RB[k, cols] / GRB[k, cols])
        blocks.append(res)
    priors = pd.DataFrame(priors)
    log(f"[tests] {len(keep_states):,} states; {int((priors.reported > 0).sum()):,} with >= 1 reported isoform; "
        f"{int((priors.tested > 0).sum()):,} with >= 1 tested isoform [{time.time() - t0:.0f}s]")
    if not blocks:
        raise ValueError("no isoform passed the read gate in any state")
    T = {k: np.concatenate([b[k] for b in blocks]) for k in blocks[0]}
    undefined = int((T["tested"] & ~np.isfinite(T["p"])).sum())
    T["fdr"] = bh_fdr(np.where(T["tested"], T["p"], np.nan))
    log(f"[tests] {len(T['p']):,} isoform x state rows; {int(T['tested'].sum()):,} tested (>= {min_donors} donors "
        f"per group); {int((~T['tested']).sum()):,} reported without a test; tested with undefined t (posterior "
        f"variance 0; P left NaN, excluded from BH): {undefined:,}; states without a fitted prior: "
        f"{int(priors.get('s0_squared', pd.Series(dtype=float)).isna().sum()):,}; states with d0 = Inf: "
        f"{int(np.isinf(priors.get('d0', pd.Series(dtype=float))).sum()):,}")

    sig = T["fdr"] < alpha
    if read_rule == "isoform":
        deep_A, deep_B = T["reads_A"] > min_reads, T["reads_B"] > min_reads
    else:
        deep_A = deep_B = (T["gene_reads_A"] > min_reads) | (T["gene_reads_B"] > min_reads)
    up_A = sig & (T["delta_usage"] >= min_delta) & deep_A
    up_B = sig & (T["delta_usage"] <= -min_delta) & deep_B
    st = np.array(keep_states, dtype=object)
    tests = pd.DataFrame({"state": st[T["ik"]], "group_A": group_A, "group_B": group_B,
                          "gene_id": cand.gene_id.to_numpy()[T["col"]],
                          "gene_symbol": cand.gene_symbol.to_numpy()[T["col"]],
                          "feature": cand.feature.to_numpy()[T["col"]],
                          "n_donors_A": T["n_donors_A"], "n_donors_B": T["n_donors_B"], "tested": T["tested"],
                          "reads_A": T["reads_A"].astype(np.int64), "reads_B": T["reads_B"].astype(np.int64),
                          "gene_reads_A": T["gene_reads_A"].astype(np.int64),
                          "gene_reads_B": T["gene_reads_B"].astype(np.int64),
                          "cohort_reads": T["cohort_reads"].astype(np.int64),
                          "cohort_gene_reads": T["cohort_gene_reads"].astype(np.int64),
                          "pooled_usage_A": T["pooled_usage_A"], "pooled_usage_B": T["pooled_usage_B"],
                          "mean_usage_A": T["mean_usage_A"], "mean_usage_B": T["mean_usage_B"],
                          "delta_usage": T["delta_usage"], "sd_A": T["sd_A"], "sd_B": T["sd_B"],
                          "sd_pooled": T["sd_pooled"], "t_moderated": T["t_moderated"], "df_total": T["df_total"],
                          "p": T["p"], "fdr": T["fdr"],
                          "higher_in": np.where(up_A, "A", np.where(up_B, "B", ""))})
    log(f"[calls] FDR < {alpha}: {int(sig.sum()):,}; higher in {group_A}: {int(up_A.sum()):,}; higher in "
        f"{group_B}: {int(up_B.sum()):,} of {int(tests.tested.sum()):,} tests")

    hits = tests[tests.higher_in != ""]
    keys = ["state", "gene_id"]
    stat_cols = ["feature", "reads_A", "reads_B", "pooled_usage_A", "pooled_usage_B", "mean_usage_A",
                 "mean_usage_B", "delta_usage", "t_moderated", "p", "fdr"]
    sw = hits[hits.higher_in == "A"][keys + ["gene_symbol", "n_donors_A", "n_donors_B"] + stat_cols].merge(
        hits[hits.higher_in == "B"][keys + stat_cols], on=keys, suffixes=("_iso_A", "_iso_B"))
    sw = sw.rename(columns={"feature_iso_A": "isoform_higher_in_A", "feature_iso_B": "isoform_higher_in_B"})
    sw.insert(1, "group_A", group_A)
    sw.insert(2, "group_B", group_B)
    ref = cand.set_index("feature")
    for s_ in ("A", "B"):
        f = sw[f"isoform_higher_in_{s_}"]
        for c in ("protein_length", "nmd_status", "de_up_states", "de_down_states", "event_states",
                  "whitelist_source"):
            sw[f"{c}_iso_{s_}"] = f.map(ref[c]).to_numpy()
    rank = {s: i for i, s in enumerate(keep_states)}
    sw = sw.assign(_r=sw.state.map(rank)).sort_values(["_r", "gene_symbol", "fdr_iso_A"]).drop(columns="_r")
    sw = sw.reset_index(drop=True)
    log(f"[switches] {len(sw):,} isoform pairs; {sw[keys].drop_duplicates().shape[0]:,} gene x state switches; "
        f"{sw.gene_id.nunique():,} genes; {sw.state.nunique():,} states [{time.time() - t0:.0f}s]")

    per_state = []
    for s in keep_states:
        m = sw.state == s
        tm = tests.state == s
        per_state.append(dict(state=s, donors_A=nA_state[s], donors_B=nB_state[s],
                              tested_isoforms=int((tm & tests.tested).sum()),
                              higher_in_A=int((tm & (tests.higher_in == "A")).sum()),
                              higher_in_B=int((tm & (tests.higher_in == "B")).sum()),
                              isoform_pairs=int(m.sum()), genes=int(sw.loc[m, "gene_id"].nunique())))

    wl.to_csv(os.path.join(outdir, "isoform_state_switch_whitelist.tsv"), sep="\t", index=False)
    cand.to_csv(os.path.join(outdir, "isoform_state_switch_candidates.tsv"), sep="\t", index=False)
    priors.to_csv(os.path.join(outdir, "isoform_state_switch_priors.tsv"), sep="\t", index=False)
    pd.DataFrame(per_state).to_csv(os.path.join(outdir, "isoform_state_switch_summary.tsv"), sep="\t", index=False)
    sw.to_csv(os.path.join(outdir, "isoform_state_switches.tsv"), sep="\t", index=False)
    tests.to_csv(os.path.join(outdir, "isoform_state_usage_tests.tsv.gz"), sep="\t", index=False,
                 float_format="%.9g", compression="gzip")
    missing = check_required_switches(sw, required, outdir, log=log) if required else []
    params.update(diff_files=manifest, event_files=ev_manifest, counts_sha256=_sha256(counts),
                  groups_sha256=_sha256(groups), protein_classes_sha256=_sha256(protein_classes), states=keep_states,
                  donors_A=donors_A, donors_B=donors_B, donors_per_state_A=nA_state, donors_per_state_B=nB_state,
                  whitelist_isoforms=len(wl), coding_isoforms=int(is_coding.sum()), pool_isoforms=len(pool),
                  candidate_isoforms=len(cand), candidate_genes=int(cand.gene_id.nunique()), rows=len(tests),
                  tests=int(tests.tested.sum()), undefined_t=undefined, fdr_significant=int(sig.sum()),
                  isoform_pair_switches=len(sw), controls_missing=missing,
                  finished=time.strftime("%Y-%m-%d %H:%M:%S"), seconds=round(time.time() - t0, 1))
    with open(os.path.join(outdir, "logs", "parameters.json"), "w") as f:
        json.dump(params, f, indent=2)
    if missing:
        with open(os.path.join(outdir, "RUN_INVALID_missing_positive_controls.txt"), "w") as f:
            f.write("\n".join(missing) + "\n")
        raise RuntimeError(f"required positive control(s) not called: {missing}; see {outdir}/positive_controls.tsv")
    log(f"[done] wrote {outdir} [{time.time() - t0:.0f}s]")
    return sw, tests, priors


def run(counts, groups, group, diff_dir, protein_classes, outdir, fold=1.5, adjp=0.05, min_reads=50,
        min_delta=0.2, alpha=0.05, min_donors=2, states=None, cluster_order=None,
        diff_pattern="*_stats-annotated.txt", read_rule="isoform", coding_only=True, read_scope="group"):
    if read_rule not in ("isoform", "gene"):
        raise ValueError("read_rule must be isoform or gene")
    if read_scope not in ("cohort", "group"):
        raise ValueError("read_scope must be cohort or group")
    t0 = time.time()
    os.makedirs(os.path.join(outdir, "logs"), exist_ok=True)
    log = _Log(os.path.join(outdir, "logs", "run_log.txt"))
    params = dict(tool="altanalyze3.components.long_read.isoform_state_switches", argv=sys.argv,
                  counts=os.path.abspath(counts), groups=os.path.abspath(groups), group=group,
                  diff_dir=os.path.abspath(diff_dir), diff_pattern=diff_pattern,
                  protein_classes=os.path.abspath(protein_classes),
                  protein_coding_rule="protein_length > 0 and nmd_status == Not-NMD", coding_only=coding_only,
                  outdir=os.path.abspath(outdir), fold=fold, log2_fold_cut=float(np.log2(fold)), adjp=adjp,
                  min_reads=min_reads, read_rule=read_rule, read_scope=read_scope, min_delta=min_delta, alpha=alpha, min_donors=min_donors,
                  states=states,
                  cluster_order=os.path.abspath(cluster_order) if cluster_order else None,
                  module_sha256=_sha256(os.path.abspath(__file__)), python=platform.python_version(),
                  numpy=np.__version__, pandas=pd.__version__, scipy=sys.modules["scipy"].__version__,
                  host=platform.node(), started=time.strftime("%Y-%m-%d %H:%M:%S"))
    log(json.dumps(params, indent=2))

    wl, manifest = load_differential_whitelist(diff_dir, fold, adjp, diff_pattern, log=log)
    log(f"[whitelist] {len(wl):,} isoforms differential in >= 1 of {len(manifest)} state comparisons "
        f"[{time.time() - t0:.0f}s]")
    pc = load_protein_classes(protein_classes)
    wl["protein_class"] = wl.feature.map(pc["class"])
    wl["protein_length"] = wl.feature.map(pc["protein_length"])
    wl["nmd_status"] = wl.feature.map(pc["nmd_status"])
    miss = int(wl.protein_class.isna().sum())
    log(f"[protein] whitelist isoforms without a protein class row: {miss:,} of {len(wl):,}")
    log(f"[protein] class counts in whitelist: {wl.protein_class.value_counts(dropna=False).to_dict()}")
    is_coding = (pd.to_numeric(wl.protein_length, errors="coerce") > 0) & (wl.nmd_status == "Not-NMD")
    pool = wl[is_coding].copy() if coding_only else wl.copy()
    pool["gene_id"] = pool.feature.str.split(":", n=1).str[0]
    per_gene = pool.groupby("gene_id").feature.transform("size")
    cand = pool[per_gene >= 2].sort_values(["gene_id", "feature"]).reset_index(drop=True)
    log(f"[candidates] protein length > 0 and Not-NMD: {int(is_coding.sum()):,} of {len(wl):,} whitelist isoforms "
        f"(classes {wl[is_coding].protein_class.value_counts().to_dict()}); pool "
        f"{'protein coding' if coding_only else 'every whitelisted isoform (--all_isoforms)'}: {len(pool):,}; "
        f"in genes with >= 2: {len(cand):,} isoforms, {cand.gene_id.nunique():,} genes")
    if cand.empty:
        raise ValueError("no candidate isoforms")
    wl["protein_coding"] = is_coding.to_numpy()
    wl["candidate"] = wl.feature.isin(cand.feature)
    wl["excluded_by"] = np.select([~wl.protein_coding & coding_only, ~wl.candidate],
                                  ["not protein coding (length 0 or NMD)",
                                   f"only {'protein-coding ' if coding_only else ''}whitelisted isoform of its gene"],
                                  "")

    data = load_group_usage(counts, groups, group, cand.feature.tolist(), log=log)
    present = sorted(set(data["state"]))
    donors = sorted(set(data["uid"]))
    n_by_state = pd.Series(data["state"]).value_counts()
    if states is not None:
        unknown = sorted(set(states) - set(present))
        if unknown:
            raise ValueError(f"requested states absent from the {group} pseudobulks: {unknown}")
    order = [s.strip() for s in open(cluster_order) if s.strip()] if cluster_order else []
    ranked = [s for s in order if s in present] + sorted(s for s in present if s not in order)
    keep_states = [s for s in ranked if states is None or s in states]
    few = [f"{s}({n_by_state[s]})" for s in keep_states if n_by_state[s] < min_donors]
    log(f"[states] {len(keep_states)} states kept; fewer than {min_donors} donors (reported, never tested): "
        f"{', '.join(few) or 'none'}")
    U, R, GR = state_arrays(data, keep_states, donors)
    ci = pd.Index(data["cohort_states"]).get_indexer(keep_states)
    if (ci < 0).any():
        raise ValueError("a kept state is absent from the cohort pseudobulks")
    CR, CGR = data["cohort_reads"][ci], data["cohort_gene_reads"][ci]
    if (CR < R).any() or (CGR < GR).any():
        raise ValueError("cohort reads below group reads: the cohort sum must include the group")
    gate_R, gate_GR = (CR, CGR) if read_scope == "cohort" else (R, GR)
    log(f"[reads] read gate pools {'every pseudobulk of the state (all groups)' if read_scope == 'cohort' else group + ' pseudobulks only'}")

    blocks, priors = [], []
    for ia, ib in combinations(range(len(keep_states)), 2):
        depth = gate_R if read_rule == "isoform" else gate_GR
        cols, res, info = test_state_pair(U[ia], U[ib], depth[ia], depth[ib], min_reads, min_donors)
        prior = dict(state_A=keep_states[ia], state_B=keep_states[ib], reported=0, tested=0)
        if res is None:
            priors.append(prior)
            continue
        prior.update(info)
        priors.append(prior)
        with np.errstate(divide="ignore", invalid="ignore"):
            res.update(ia=np.full(cols.size, ia), ib=np.full(cols.size, ib), col=cols, reads_A=R[ia, cols],
                       reads_B=R[ib, cols], gene_reads_A=GR[ia, cols], gene_reads_B=GR[ib, cols],
                       cohort_reads_A=CR[ia, cols], cohort_reads_B=CR[ib, cols],
                       cohort_gene_reads_A=CGR[ia, cols], cohort_gene_reads_B=CGR[ib, cols],
                       pooled_usage_A=R[ia, cols] / GR[ia, cols], pooled_usage_B=R[ib, cols] / GR[ib, cols])
        blocks.append(res)
    priors = pd.DataFrame(priors)
    n_pairs = len(keep_states) * (len(keep_states) - 1) // 2
    log(f"[tests] {n_pairs:,} state pairs; {int((priors.reported > 0).sum()):,} with >= 1 reported isoform; "
        f"{int((priors.tested > 0).sum()):,} with >= 1 tested isoform [{time.time() - t0:.0f}s]")
    if not blocks:
        raise ValueError("no isoform passed the read gate in any state pair")
    T = {k: np.concatenate([b[k] for b in blocks]) for k in blocks[0]}
    T["fdr"] = bh_fdr(np.where(T["tested"], T["p"], np.nan))
    nonfinite = int((T["tested"] & ~np.isfinite(T["p"])).sum())
    log(f"[tests] {len(T['p']):,} isoform x state-pair rows; {int(T['tested'].sum()):,} tested "
        f"(>= {min_donors} paired donors); {int((~T['tested']).sum()):,} reported without a test; "
        f"tested with non-finite P (excluded from BH): {nonfinite:,}; "
        f"pairs without a fitted prior: {int(priors.get('s0_squared', pd.Series(dtype=float)).isna().sum()):,}; "
        f"pairs with d0 = Inf: {int(np.isinf(priors.get('d0', pd.Series(dtype=float))).sum()):,}")

    sig = T["fdr"] < alpha
    pre = "cohort_" if read_scope == "cohort" else ""
    if read_rule == "isoform":
        deep_A, deep_B = T[pre + "reads_A"] > min_reads, T[pre + "reads_B"] > min_reads
    else:
        deep_A = deep_B = (T[pre + "gene_reads_A"] > min_reads) | (T[pre + "gene_reads_B"] > min_reads)
    up_A = sig & (T["delta_usage"] >= min_delta) & deep_A
    up_B = sig & (T["delta_usage"] <= -min_delta) & deep_B
    higher = np.where(up_A, "A", np.where(up_B, "B", ""))
    st = np.array(keep_states, dtype=object)
    tests = pd.DataFrame({"state_A": st[T["ia"]], "state_B": st[T["ib"]],
                          "gene_id": cand.gene_id.to_numpy()[T["col"]],
                          "gene_symbol": cand.gene_symbol.to_numpy()[T["col"]],
                          "feature": cand.feature.to_numpy()[T["col"]],
                          "n_donors": T["n_donors"], "tested": T["tested"],
                          "reads_A": T["reads_A"].astype(np.int64),
                          "reads_B": T["reads_B"].astype(np.int64),
                          "gene_reads_A": T["gene_reads_A"].astype(np.int64),
                          "gene_reads_B": T["gene_reads_B"].astype(np.int64),
                          "cohort_reads_A": T["cohort_reads_A"].astype(np.int64),
                          "cohort_reads_B": T["cohort_reads_B"].astype(np.int64),
                          "cohort_gene_reads_A": T["cohort_gene_reads_A"].astype(np.int64),
                          "cohort_gene_reads_B": T["cohort_gene_reads_B"].astype(np.int64),
                          "pooled_usage_A": T["pooled_usage_A"], "pooled_usage_B": T["pooled_usage_B"],
                          "mean_usage_A": T["mean_usage_A"], "mean_usage_B": T["mean_usage_B"],
                          "delta_usage": T["delta_usage"], "sd_delta": T["sd_delta"],
                          "t_moderated": T["t_moderated"], "df_total": T["df_total"], "p": T["p"],
                          "fdr": T["fdr"], "higher_in": higher})
    log(f"[calls] FDR < {alpha}: {int(sig.sum()):,}; higher in A: {int(up_A.sum()):,}; higher in B: "
        f"{int(up_B.sum()):,} of {int(tests.tested.sum()):,} tests")

    hits = tests[tests.higher_in != ""]
    keys = ["state_A", "state_B", "gene_id"]
    a_side = hits[hits.higher_in == "A"]
    b_side = hits[hits.higher_in == "B"]
    stat_cols = ["feature", "reads_A", "reads_B", "cohort_reads_A", "cohort_reads_B", "pooled_usage_A", "pooled_usage_B", "mean_usage_A", "mean_usage_B",
                 "delta_usage", "t_moderated", "fdr"]
    sw = a_side[keys + ["gene_symbol", "n_donors"] + stat_cols].merge(
        b_side[keys + stat_cols], on=keys, suffixes=("_iso_A", "_iso_B"))
    sw = sw.rename(columns={"feature_iso_A": "isoform_higher_in_A", "feature_iso_B": "isoform_higher_in_B"})
    for side in ("A", "B"):
        ref = cand.set_index("feature")
        f = sw[f"isoform_higher_in_{side}"]
        sw[f"protein_length_iso_{side}"] = f.map(ref.protein_length).to_numpy()
        sw[f"de_up_states_iso_{side}"] = f.map(ref.de_up_states).to_numpy()
        sw[f"de_down_states_iso_{side}"] = f.map(ref.de_down_states).to_numpy()
    sw = sw.sort_values(["state_A", "state_B", "gene_symbol", "fdr_iso_A"]).reset_index(drop=True)
    log(f"[switches] {len(sw):,} isoform pairs; {sw[keys].drop_duplicates().shape[0]:,} gene x state-pair "
        f"switches; {sw.gene_id.nunique():,} genes; {sw[['state_A', 'state_B']].drop_duplicates().shape[0]:,} "
        f"state pairs [{time.time() - t0:.0f}s]")

    per_state = []
    for s in keep_states:
        m = (sw.state_A == s) | (sw.state_B == s)
        per_state.append(dict(state=s, donors=int(n_by_state[s]), isoform_pairs=int(m.sum()),
                              gene_state_pair_switches=int(sw.loc[m, keys].drop_duplicates().shape[0]),
                              genes=int(sw.loc[m, "gene_id"].nunique()),
                              partner_states=int(pd.concat([sw.loc[m & (sw.state_A == s), "state_B"],
                                                            sw.loc[m & (sw.state_B == s), "state_A"]]).nunique())))

    wl.to_csv(os.path.join(outdir, "isoform_state_switch_whitelist.tsv"), sep="\t", index=False)
    cand.to_csv(os.path.join(outdir, "isoform_state_switch_candidates.tsv"), sep="\t", index=False)
    priors.to_csv(os.path.join(outdir, "isoform_state_switch_priors.tsv"), sep="\t", index=False)
    pd.DataFrame(per_state).to_csv(os.path.join(outdir, "isoform_state_switch_summary.tsv"), sep="\t", index=False)
    sw.to_csv(os.path.join(outdir, "isoform_state_switches.tsv"), sep="\t", index=False)
    tests.to_csv(os.path.join(outdir, "isoform_state_usage_tests.tsv.gz"), sep="\t", index=False,
                 float_format="%.9g", compression="gzip")
    params.update(diff_files=manifest, counts_sha256=_sha256(counts), groups_sha256=_sha256(groups),
                  protein_classes_sha256=_sha256(protein_classes), states=keep_states,
                  donors_per_state={s: int(n_by_state[s]) for s in keep_states}, donors=donors,
                  whitelist_isoforms=len(wl), coding_isoforms=int(is_coding.sum()), pool_isoforms=len(pool),
                  candidate_isoforms=len(cand),
                  candidate_genes=int(cand.gene_id.nunique()), state_pairs=n_pairs, rows=len(tests),
                  tests=int(tests.tested.sum()),
                  nonfinite_p=nonfinite, fdr_significant=int(sig.sum()), isoform_pair_switches=len(sw),
                  finished=time.strftime("%Y-%m-%d %H:%M:%S"), seconds=round(time.time() - t0, 1))
    with open(os.path.join(outdir, "logs", "parameters.json"), "w") as f:
        json.dump(params, f, indent=2)
    log(f"[done] wrote {outdir} [{time.time() - t0:.0f}s]")
    return sw, tests, priors


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] == "predominance":
        return main_predominance(argv[1:])
    if argv and argv[0] == "protein-outcome":
        return main_protein_outcome(argv[1:])
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True, help="cohort isoform_combined_pseudo_cluster_counts.h5ad")
    ap.add_argument("--groups", required=True, help="metadata with uid and groups columns")
    ap.add_argument("--group", default=None, help="group whose donors are paired (e.g. healthy); paired mode")
    ap.add_argument("--compare_groups", default=None,
                    help="A,B (e.g. aged,young): test group A donors against group B donors within each state")
    ap.add_argument("--diff_dir", required=True,
                    help="paired mode: folder of <group>-<state>-Others-<state>_stats-annotated.txt; "
                         "--compare_groups: folder of <A>-<B>-<state>_stats-annotated.txt (diff-covariate-isoform)")
    ap.add_argument("--whitelist_rawp", type=float, default=None,
                    help="--compare_groups: whitelist on ebayesPval < this instead of ebayesAdjPval < --adjp")
    ap.add_argument("--event_dir", default=None,
                    help="--compare_groups: dPSI-covariate-events folder; adds both isoforms of each event passing "
                         "--event_dpsi and --event_rawp to the whitelist")
    # user directive 2026-10-21: isoforms of splicing events with dPSI > 0.2 and raw P < 0.05
    ap.add_argument("--event_dpsi", type=float, default=0.2, help="--event_dir: |dPSI| > this")
    ap.add_argument("--event_rawp", type=float, default=0.05, help="--event_dir: ebayesPval < this")
    ap.add_argument("--required_switch", action="append", default=None, metavar="FEATURE_A,FEATURE_B",
                    help="--compare_groups: positive control, FEATURE_A higher in A and FEATURE_B higher in B in "
                         ">= 1 state (repeatable); the run fails when one is not called")
    ap.add_argument("--diff_pattern", default="*_stats-annotated.txt")
    ap.add_argument("--protein_classes", required=True,
                    help="ALL_catalog_protein_classes.tsv; coding = protein_length > 0 and nmd_status Not-NMD")
    ap.add_argument("--all_isoforms", action="store_true",
                    help="keep every whitelisted isoform, coding or not (default: protein coding only)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--fold", type=float, default=1.5, help="whitelist: |DiffMeans| > log2(fold)")
    ap.add_argument("--adjp", type=float, default=0.05, help="whitelist: ebayesAdjPval < adjp")
    # user directive 2026-10-20: minimum 50 healthy reads, so the TNFAIP8 positive control is kept
    ap.add_argument("--min_reads", type=float, default=50, help="pooled read depth > this (see --read_rule)")
    ap.add_argument("--read_rule", choices=["isoform", "gene"], default="isoform",
                    help="isoform: isoform reads in its higher state; gene: gene reads in A or in B")
    ap.add_argument("--read_scope", choices=["cohort", "group"], default="group",
                    help="pseudobulks pooled for the read depth: cohort = every pseudobulk of the state (all groups); "
                         "group = --group pseudobulks only")
    ap.add_argument("--min_delta", type=float, default=0.2, help="mean paired usage difference >= this")
    ap.add_argument("--alpha", type=float, default=0.05, help="BH FDR < alpha")
    ap.add_argument("--min_donors", type=int, default=2,
                    help="paired mode: paired donors needed for a test; --compare_groups: donors needed in EACH "
                         "group (>= 2)")
    ap.add_argument("--states", default=None, help="comma-separated subset of states")
    ap.add_argument("--cluster_order", default=None, help="one state per line; orients state_A before state_B")
    a = ap.parse_args(argv)
    if a.compare_groups:
        gs = [g.strip() for g in a.compare_groups.split(",")]
        if len(gs) != 2 or not all(gs):
            ap.error("--compare_groups takes exactly two group names, A,B")
        if a.group:
            ap.error("--group (paired mode) and --compare_groups cannot be combined")
        if a.read_scope != "group":
            ap.error("--compare_groups pools each group's pseudobulks separately; only --read_scope group applies")
        required = []
        for r in a.required_switch or []:
            pair = [x.strip() for x in r.split(",")]
            if len(pair) != 2 or not all(pair):
                ap.error(f"--required_switch {r}: expected FEATURE_A,FEATURE_B")
            required.append(tuple(pair))
        run_groups(a.counts, a.groups, gs[0], gs[1], a.diff_dir, a.protein_classes, a.outdir, fold=a.fold,
                   adjp=a.adjp, rawp=a.whitelist_rawp, event_dir=a.event_dir, event_dpsi=a.event_dpsi,
                   event_rawp=a.event_rawp, min_reads=a.min_reads, min_delta=a.min_delta, alpha=a.alpha,
                   min_donors=a.min_donors, states=[s.strip() for s in a.states.split(",")] if a.states else None,
                   cluster_order=a.cluster_order, diff_pattern=a.diff_pattern, read_rule=a.read_rule,
                   coding_only=not a.all_isoforms, required=required)
        return
    if not a.group:
        ap.error("paired mode needs --group (or use --compare_groups A,B)")
    for flag, val in (("--whitelist_rawp", a.whitelist_rawp), ("--event_dir", a.event_dir),
                      ("--required_switch", a.required_switch)):
        if val is not None:
            ap.error(f"{flag} applies to --compare_groups only")
    run(a.counts, a.groups, a.group, a.diff_dir, a.protein_classes, a.outdir, fold=a.fold, adjp=a.adjp,
        min_reads=a.min_reads, min_delta=a.min_delta, alpha=a.alpha, min_donors=a.min_donors,
        states=[s.strip() for s in a.states.split(",")] if a.states else None, cluster_order=a.cluster_order,
        diff_pattern=a.diff_pattern, read_rule=a.read_rule, coding_only=not a.all_isoforms,
        read_scope=a.read_scope)


if __name__ == "__main__":
    main()
