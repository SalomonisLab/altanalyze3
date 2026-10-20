"""Isoform switches between pairs of cell states, tested within donors.

A switch between states A and B is a pair of isoforms of one gene: isoform i has higher usage in A and
isoform j has higher usage in B. Usage is isoform reads / gene reads in one donor's pseudobulk, as computed by
isoform_matrix.cpm_and_ratio_from_counts (missing where the gene CPM < 1 in that pseudobulk).

Candidates are isoforms that are differential in at least one state-vs-Others comparison
(|DiffMeans| > log2(--fold) and ebayesAdjPval < --adjp), are protein coding (predicted protein length > 0 and Not-NMD
in --protein_classes, with no length-relative-to-longest rule), and belong to a gene with >= 2 such isoforms.
Every state of the group is kept. For every candidate and every pair of states, an isoform is reported when the read
depth exceeds --min_reads in A or in B and >= 1 donor has a usage in both states, and tested when >= --min_donors
(>= 2) donors do. --read_rule sets the depth: isoform (pooled isoform reads) or gene (pooled reads of its gene):
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


def load_protein_classes(path):
    p = pd.read_csv(path, sep="\t", dtype=str, usecols=["gene", "isoform", "class", "protein_length", "nmd_status"])
    p["feature"] = p["gene"] + ":" + p["isoform"]
    if not p.feature.is_unique:
        raise ValueError(f"{path}: duplicate gene:isoform rows")
    return p.set_index("feature")


def load_group_usage(counts_h5ad, groups_file, group, features, log=print, cpm_threshold=1):
    """Donor x state usage, isoform reads and gene reads for `features`, from the cohort counts h5ad.
    Usage comes from isoform_matrix.cpm_and_ratio_from_counts over ALL features of the group's pseudobulks."""
    import anndata as ad
    from .comparisons import _resolve_uid
    from .isoform_matrix import cpm_and_ratio_from_counts

    meta = pd.read_csv(groups_file, sep="\t", dtype=str)
    pairs = meta[["uid", "groups"]].apply(lambda s: s.str.strip()).drop_duplicates()
    if pairs.uid.duplicated().any():
        raise ValueError(f"{groups_file}: a uid maps to more than one group")
    uid_to_group = dict(zip(pairs.uid, pairs.groups))
    a = ad.read_h5ad(counts_h5ad)
    names = np.asarray(a.obs_names, dtype=str)
    split = [n.split(".", 1) for n in names]
    if any(len(s) != 2 for s in split):
        raise ValueError("counts obs names must be <state>.<sample>")
    state = np.array([s[0] for s in split])
    uid = np.array([_resolve_uid(s[1], uid_to_group, "isoform") for s in split])
    keep = np.array([uid_to_group.get(u) == group for u in uid])
    log(f"[counts] {counts_h5ad}: {a.n_obs:,} pseudobulks x {a.n_vars:,} features; {group}: {int(keep.sum()):,} "
        f"pseudobulks, {len(set(uid[keep]))} donors, {len(set(state[keep]))} states")
    if not keep.any():
        raise ValueError(f"no pseudobulk matches group {group}")
    X = sparse.csr_matrix(a.X[np.where(keep)[0]])
    var = a.var_names.to_numpy(dtype=object)
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
    fpos = pd.Index(var).get_indexer(features)
    if (fpos < 0).any():
        raise ValueError(f"{int((fpos < 0).sum())} candidate features absent from {counts_h5ad}")
    fgene = np.array([f.split(":", 1)[0] for f in features])
    gpos = pd.Index(genes).get_indexer(fgene)
    if (gpos < 0).any():
        raise ValueError("candidate gene missing from the ratio gene index")
    usage = ratio[:, fpos].toarray().astype(float)
    expressed = E[:, gpos].toarray().astype(bool)
    usage[~expressed] = np.nan
    reads = X[:, fpos].toarray()
    cand_genes, ginv = np.unique(fgene, return_inverse=True)
    in_cand = pd.Index(cand_genes).get_indexer([v.split(":", 1)[0] for v in var])
    sel = np.where(in_cand >= 0)[0]
    Gm = sparse.csr_matrix((np.ones(len(sel)), (sel, in_cand[sel])), shape=(len(var), len(cand_genes)))
    gene_reads = (X @ Gm).toarray()[:, ginv]
    return dict(state=pb_state, uid=pb_uid, usage=usage, reads=reads, gene_reads=gene_reads)


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


def run(counts, groups, group, diff_dir, protein_classes, outdir, fold=1.5, adjp=0.05, min_reads=200,
        min_delta=0.2, alpha=0.05, min_donors=2, states=None, cluster_order=None,
        diff_pattern="*_stats-annotated.txt", read_rule="isoform"):
    if read_rule not in ("isoform", "gene"):
        raise ValueError("read_rule must be isoform or gene")
    t0 = time.time()
    os.makedirs(os.path.join(outdir, "logs"), exist_ok=True)
    log = _Log(os.path.join(outdir, "logs", "run_log.txt"))
    params = dict(tool="altanalyze3.components.long_read.isoform_state_switches", argv=sys.argv,
                  counts=os.path.abspath(counts), groups=os.path.abspath(groups), group=group,
                  diff_dir=os.path.abspath(diff_dir), diff_pattern=diff_pattern,
                  protein_classes=os.path.abspath(protein_classes),
                  protein_coding_rule="protein_length > 0 and nmd_status == Not-NMD",
                  outdir=os.path.abspath(outdir), fold=fold, log2_fold_cut=float(np.log2(fold)), adjp=adjp,
                  min_reads=min_reads, read_rule=read_rule, min_delta=min_delta, alpha=alpha, min_donors=min_donors,
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
    coding = wl[is_coding].copy()
    coding["gene_id"] = coding.feature.str.split(":", n=1).str[0]
    per_gene = coding.groupby("gene_id").feature.transform("size")
    cand = coding[per_gene >= 2].sort_values(["gene_id", "feature"]).reset_index(drop=True)
    log(f"[candidates] protein length > 0 and Not-NMD: {len(coding):,} of {len(wl):,} whitelist isoforms "
        f"(classes {coding.protein_class.value_counts().to_dict()}); in genes with >= 2: "
        f"{len(cand):,} isoforms, {cand.gene_id.nunique():,} genes")
    if cand.empty:
        raise ValueError("no candidate isoforms")
    wl["protein_coding"] = is_coding.to_numpy()
    wl["candidate"] = wl.feature.isin(cand.feature)
    wl["excluded_by"] = np.select([~wl.protein_coding, ~wl.candidate],
                                  ["not protein coding (length 0 or NMD)", "only protein-coding whitelisted isoform "
                                   "of its gene"], "")

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

    blocks, priors = [], []
    for ia, ib in combinations(range(len(keep_states)), 2):
        depth = R if read_rule == "isoform" else GR
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
    if read_rule == "isoform":
        deep_A, deep_B = T["reads_A"] > min_reads, T["reads_B"] > min_reads
    else:
        deep_A = deep_B = (T["gene_reads_A"] > min_reads) | (T["gene_reads_B"] > min_reads)
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
    stat_cols = ["feature", "reads_A", "reads_B", "pooled_usage_A", "pooled_usage_B", "mean_usage_A", "mean_usage_B",
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
                  whitelist_isoforms=len(wl), coding_isoforms=len(coding), candidate_isoforms=len(cand),
                  candidate_genes=int(cand.gene_id.nunique()), state_pairs=n_pairs, rows=len(tests),
                  tests=int(tests.tested.sum()),
                  nonfinite_p=nonfinite, fdr_significant=int(sig.sum()), isoform_pair_switches=len(sw),
                  finished=time.strftime("%Y-%m-%d %H:%M:%S"), seconds=round(time.time() - t0, 1))
    with open(os.path.join(outdir, "logs", "parameters.json"), "w") as f:
        json.dump(params, f, indent=2)
    log(f"[done] wrote {outdir} [{time.time() - t0:.0f}s]")
    return sw, tests, priors


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True, help="cohort isoform_combined_pseudo_cluster_counts.h5ad")
    ap.add_argument("--groups", required=True, help="metadata with uid and groups columns")
    ap.add_argument("--group", required=True, help="group whose donors are paired (e.g. healthy)")
    ap.add_argument("--diff_dir", required=True, help="folder of <group>-<state>-Others-<state>_stats-annotated.txt")
    ap.add_argument("--diff_pattern", default="*_stats-annotated.txt")
    ap.add_argument("--protein_classes", required=True,
                    help="ALL_catalog_protein_classes.tsv; coding = protein_length > 0 and nmd_status Not-NMD")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--fold", type=float, default=1.5, help="whitelist: |DiffMeans| > log2(fold)")
    ap.add_argument("--adjp", type=float, default=0.05, help="whitelist: ebayesAdjPval < adjp")
    ap.add_argument("--min_reads", type=float, default=200, help="pooled read depth > this (see --read_rule)")
    ap.add_argument("--read_rule", choices=["isoform", "gene"], default="isoform",
                    help="isoform: isoform reads in its higher state; gene: gene reads in A or in B")
    ap.add_argument("--min_delta", type=float, default=0.2, help="mean paired usage difference >= this")
    ap.add_argument("--alpha", type=float, default=0.05, help="BH FDR < alpha")
    ap.add_argument("--min_donors", type=int, default=2, help="paired donors needed for a test (>= 2)")
    ap.add_argument("--states", default=None, help="comma-separated subset of states")
    ap.add_argument("--cluster_order", default=None, help="one state per line; orients state_A before state_B")
    a = ap.parse_args(argv)
    run(a.counts, a.groups, a.group, a.diff_dir, a.protein_classes, a.outdir, fold=a.fold, adjp=a.adjp,
        min_reads=a.min_reads, min_delta=a.min_delta, alpha=a.alpha, min_donors=a.min_donors,
        states=[s.strip() for s in a.states.split(",")] if a.states else None, cluster_order=a.cluster_order,
        diff_pattern=a.diff_pattern, read_rule=a.read_rule)


if __name__ == "__main__":
    main()
