#!/usr/bin/env python3.11
"""Rigorous directional validation of imputed differential TF activity.

Design (per the lab spec, generalized to K-fold for statistical power):
  - Exclude AML-7 entirely (reserved as the real-world raw-h5 speed input).
  - Control = Multiome_WT, HELD OUT OF TRAINING in every fold (RC2_TEA stays in
    training). Each fold also holds out a group of AML samples ("multiple AML +
    one control" held out together). Every AML pseudobulk is tested while held out.
  - For each held-out AML pseudobulk (sample S, cell state c) vs the cell-state-
    matched control Multiome_WT(c):
        true_diff[tf]  = trueTFactivity(S,c)   - trueTFactivity(MW,c)
        imp_diff[tf]   = impTFactivity(S,c)    - impTFactivity(MW,c)   (both imputed)
        tf_expr_diff[tf] = expr(tf gene; S,c)  - expr(tf gene; MW,c)   (RNA, CP10k+log1p)
  - Among TFs that are TRULY differential (|true_diff| large), report the % of
    examples where the imputed differential has the CORRECT sign (direction).
  - Stratify by whether the TF's own expression is also differential, to isolate
    the regulon-driven (expression-stable) case. Baselines: chance, balanced
    accuracy, TF-expression-sign, and kNN retrieval.

Outputs a long-format example table + a metrics JSON for the validation report.
"""
import json
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import binomtest

OUT = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn"
RIDGE = 1.0
N_FOLDS = 5
EXPR_DIFF_THR = 0.25     # |log1p CP10k| change to call a TF's own expression differential


def log(*a): print(*a, flush=True)


d = np.load(f"{OUT}/matched/dataset.npz", allow_pickle=True)
X = d["X"].astype(np.float64); Y = d["Y"].astype(np.float64)
genes = d["genes"].astype(str)
edge_tf = d["edge_tf"].astype(str); edge_gene = d["edge_gene"].astype(str)
sample = d["sample"].astype(str); cs = d["cell_state"].astype(str); group = d["group"].astype(str)
gp = {g: i for i, g in enumerate(genes)}; E = len(edge_tf)
tfs = np.array(sorted(set(edge_tf))); tfpos = {t: i for i, t in enumerate(tfs)}
edge_tf_idx = np.array([tfpos[t] for t in edge_tf])
tf_gene_col = {t: (gp[t] if t in gp else -1) for t in tfs}

# feature space (TF u target genes)
feat_genes = [g for g in sorted(set(edge_tf) | set(edge_gene)) if g in gp]
fpos = {g: i for i, g in enumerate(feat_genes)}
feat_cols = np.array([gp[g] for g in feat_genes])
Xf_all = X[:, feat_cols]
tgt_col = np.array([fpos[g] for g in edge_gene])
tfc_col = np.array([fpos[g] for g in edge_tf])
tf_of_edge = edge_tf_idx
# regulon averaging matrix (n_tf x n_feat)
rows_, cols_, data_ = [], [], []
for t in tfs:
    cols_t = [tgt_col[e] for e in range(E) if edge_tf[e] == t]
    w = 1.0 / len(cols_t)
    for cc in cols_t:
        rows_.append(tfpos[t]); cols_.append(cc); data_.append(w)
REGMAT = sp.csr_matrix((data_, (rows_, cols_)), shape=(len(tfs), len(feat_genes)))


def fit_peredge(tr):
    Xt = Xf_all[tr]; Yt = Y[tr]; n = tr.sum()
    t = Xt[:, tgt_col]; f = Xt[:, tfc_col]; r = np.asarray(Xt @ REGMAT.T)[:, tf_of_edge]
    mu = np.zeros((E, 3)); sd = np.ones((E, 3)); cols = [np.ones((n, E))]
    for j, A in enumerate([t, f, r]):
        m = A.mean(0); s = A.std(0); s[s < 1e-8] = 1; mu[:, j] = m; sd[:, j] = s; cols.append((A - m) / s)
    D = np.stack(cols, 2)
    AtA = np.einsum("nef,neg->efg", D, D); Aty = np.einsum("nef,ne->ef", D, Yt)
    reg = RIDGE * np.eye(4)[None].repeat(E, 0); reg[:, 0, 0] = 0
    coefs = np.linalg.solve(AtA + reg, Aty)
    return coefs, mu, sd


def pred_peredge(co, mu, sd, rows):
    Xt = Xf_all[rows]; t = Xt[:, tgt_col]; f = Xt[:, tfc_col]; r = np.asarray(Xt @ REGMAT.T)[:, tf_of_edge]
    sdc = sd.copy(); sdc[sdc < 1e-8] = 1
    ts = (t - mu[:, 0]) / sdc[:, 0]; fs = (f - mu[:, 1]) / sdc[:, 1]; rs = (r - mu[:, 2]) / sdc[:, 2]
    p = co[:, 0] + co[:, 1] * ts + co[:, 2] * fs + co[:, 3] * rs
    return np.clip(p, 0, None)


def fit_knn(tr):
    Xs = Xf_all[tr]; mu = Xs.mean(0); sd = Xs.std(0); sd[sd < 1e-8] = 1
    return mu, sd, np.where(tr)[0]


def pred_knn(state, rows, k=15):
    mu, sd, tridx = state
    A = (Xf_all[tridx] - mu) / sd; A /= (np.linalg.norm(A, axis=1, keepdims=True) + 1e-12)
    Q = (Xf_all[rows] - mu) / sd; Q /= (np.linalg.norm(Q, axis=1, keepdims=True) + 1e-12)
    sim = Q @ A.T; idx = np.argpartition(-sim, kth=k - 1, axis=1)[:, :k]
    rr = np.arange(sim.shape[0])[:, None]; w = 1 / ((1 - sim[rr, idx]) + 1e-6); w /= w.sum(1, keepdims=True)
    return np.einsum("ij,ijk->ik", w, Y[tridx][idx])


def tf_activity(edge_mat):           # (n, E) -> (n, n_tf) mean over each TF's edges
    out = np.zeros((edge_mat.shape[0], len(tfs)))
    for ti in range(len(tfs)):
        out[:, ti] = edge_mat[:, tf_of_edge == ti].mean(1)
    return out


# ---- folds over AML samples (exclude AML-7) ----
aml = np.array([s for s in pd.unique(sample) if not group[sample == s][0].endswith("control") and s != "AML-7_CCHMC"])
aml_sorted = sorted(aml)
folds = [aml_sorted[i::N_FOLDS] for i in range(N_FOLDS)]
ctrl_rows = {cs[r]: r for r in np.where(group == "Multiome_control")[0]}

recs = []
for fi, fold in enumerate(folds):
    test_aml = set(fold)
    tr = ~(np.isin(sample, list(test_aml)) | (sample == "Multiome_WT"))
    co, mu, sd = fit_peredge(tr)
    kstate = fit_knn(tr)
    # held-out AML pseudobulks
    test_rows = [r for r in np.where(np.isin(sample, list(test_aml)))[0] if cs[r] in ctrl_rows]
    crows = [ctrl_rows[cs[r]] for r in test_rows]
    if not test_rows:
        continue
    Pa = tf_activity(pred_peredge(co, mu, sd, np.array(test_rows)))
    Pc = tf_activity(pred_peredge(co, mu, sd, np.array(crows)))
    Ka = tf_activity(pred_knn(kstate, np.array(test_rows)))
    Kc = tf_activity(pred_knn(kstate, np.array(crows)))
    Ta = tf_activity(Y[np.array(test_rows)]); Tc = tf_activity(Y[np.array(crows)])
    for k, (r, cr) in enumerate(zip(test_rows, crows)):
        true_d = Ta[k] - Tc[k]; imp_d = Pa[k] - Pc[k]; knn_d = Ka[k] - Kc[k]
        expr_d = np.array([(Xf_all[r, fpos[t]] - Xf_all[cr, fpos[t]]) if t in fpos else 0.0 for t in tfs])
        mut = group[r]
        for ti, t in enumerate(tfs):
            recs.append((str(sample[r]), str(cs[r]), str(t), float(true_d[ti]),
                         float(imp_d[ti]), float(knn_d[ti]), float(expr_d[ti])))

df = pd.DataFrame(recs, columns=["sample", "cell_state", "TF", "true_diff", "imp_diff", "knn_diff", "tf_expr_diff"])
df.to_csv(f"{OUT}/benchmark/directional_examples.csv", index=False)
log(f"examples: {len(df)} over {df[['sample','cell_state']].drop_duplicates().shape[0]} held-out pseudobulks, "
    f"{df['sample'].nunique()} AML samples")

# ---- metrics ----
def dir_acc(sub, pred_col):
    s = sub[(sub[pred_col] != 0)]
    correct = np.sign(s["true_diff"]) == np.sign(s[pred_col])
    return float(correct.mean()), int(len(s))


def balanced_dir_acc(sub, pred_col):
    up = sub[sub["true_diff"] > 0]; dn = sub[sub["true_diff"] < 0]
    a_up = (np.sign(up[pred_col]) > 0).mean() if len(up) else np.nan
    a_dn = (np.sign(dn[pred_col]) < 0).mean() if len(dn) else np.nan
    return float(np.nanmean([a_up, a_dn])), float(a_up), float(a_dn), int(len(up)), int(len(dn))


metrics = {}
# magnitude strata: define "differential" by |true_diff| percentile across all examples
absd = df["true_diff"].abs()
for label, thr in [("top50pct", absd.quantile(0.50)), ("top25pct", absd.quantile(0.75)),
                   ("top10pct", absd.quantile(0.90))]:
    sub = df[absd >= thr]
    acc, n = dir_acc(sub, "imp_diff")
    bacc, aup, adn, nup, ndn = balanced_dir_acc(sub, "imp_diff")
    kacc, _ = dir_acc(sub, "knn_diff")
    # expr-sign baseline
    eacc, ne = dir_acc(sub.assign(es=np.sign(sub["tf_expr_diff"])), "es")
    bt = binomtest(int(round(acc * n)), n, 0.5, alternative="greater")
    metrics[label] = dict(thr=float(thr), n=n, peredge_dir_acc=acc, peredge_balanced_acc=bacc,
                          acc_up=aup, acc_down=adn, n_up=nup, n_down=ndn,
                          knn_dir_acc=kacc, expr_sign_acc=eacc, p_value=bt.pvalue)
    log(f"[{label}] n={n} peredge_dir={acc:.3f} balanced={bacc:.3f} (up {aup:.2f}/dn {adn:.2f}) "
        f"knn={kacc:.3f} expr_sign={eacc:.3f} p={bt.pvalue:.1e}")

# stratify the differential (top25pct) by whether TF's own expression is differential
sub = df[absd >= absd.quantile(0.75)]
stable = sub[sub["tf_expr_diff"].abs() < EXPR_DIFF_THR]
moved = sub[sub["tf_expr_diff"].abs() >= EXPR_DIFF_THR]
for nm, s in [("expr_STABLE(regulon-only)", stable), ("expr_MOVED", moved)]:
    acc, n = dir_acc(s, "imp_diff"); bacc, aup, adn, _, _ = balanced_dir_acc(s, "imp_diff")
    eacc, _ = dir_acc(s.assign(es=np.sign(s["tf_expr_diff"])), "es")
    kacc, _ = dir_acc(s, "knn_diff")
    metrics[f"strat_{nm}"] = dict(n=n, peredge_dir_acc=acc, peredge_balanced_acc=bacc,
                                  expr_sign_acc=eacc, knn_dir_acc=kacc)
    log(f"[{nm}] n={n} peredge_dir={acc:.3f} balanced={bacc:.3f} expr_sign={eacc:.3f} knn={kacc:.3f}")

# per-pseudobulk distribution (top25pct differential TFs)
per_pb = []
for (s, c), gsub in sub.groupby(["sample", "cell_state"]):
    acc, n = dir_acc(gsub, "imp_diff")
    per_pb.append(dict(sample=s, cell_state=c, n_diff_tf=n, dir_acc=acc))
ppb = pd.DataFrame(per_pb)
ppb.to_csv(f"{OUT}/benchmark/directional_per_pseudobulk.csv", index=False)
metrics["per_pseudobulk"] = dict(n_pseudobulks=len(ppb), dir_acc_median=float(ppb.dir_acc.median()),
                                 dir_acc_min=float(ppb.dir_acc.min()), dir_acc_max=float(ppb.dir_acc.max()),
                                 dir_acc_q25=float(ppb.dir_acc.quantile(.25)), dir_acc_q75=float(ppb.dir_acc.quantile(.75)))
log(f"\nper-pseudobulk dir-acc (top25pct diff TFs): median={ppb.dir_acc.median():.3f} "
    f"range [{ppb.dir_acc.min():.2f}, {ppb.dir_acc.max():.2f}]")
json.dump(metrics, open(f"{OUT}/benchmark/directional_metrics.json", "w"), indent=2)
log("\nwrote directional_metrics.json, directional_examples.csv, directional_per_pseudobulk.csv")
