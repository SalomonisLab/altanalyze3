#!/usr/bin/env python3.11
"""fig5 (current framework): does the per-edge regulon model recover differential
TF-activity DIRECTION beyond what the TF's own expression provides? Controls
retained in training; dual control (TEA preferred); per-pseudobulk offset removed.
Compares the model to a TF-expression-sign baseline, overall and stratified by
whether the TF's own mRNA is differential (the regulon-only case)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "sans-serif"; plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["ps.fonttype"] = 42
sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
from altanalyze3.components.rna2grn.evaluate import _Edges

OUT = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn"
R = f"{OUT}/reports"
GREEN = "#2CA02C"; GREY = "#999999"; RED = "#D62728"
EXPR_THR = 0.25; N_FOLDS = 5; EXCLUDE = "AML-7_CCHMC"

d = np.load(f"{OUT}/matched/dataset.npz", allow_pickle=True)
X = d["X"].astype(np.float64); Y = d["Y"].astype(np.float64); genes = d["genes"].astype(str)
etf = d["edge_tf"].astype(str); eg = d["edge_gene"].astype(str)
sample = d["sample"].astype(str); cs = d["cell_state"].astype(str); group = d["group"].astype(str)
ed = _Edges(genes, etf, eg); Xf = X[:, ed.feat_cols]; tfs = np.array(ed.tf_list)
fpos = {g: i for i, g in enumerate(ed.feat)}
tf_feat = np.array([fpos.get(t, -1) for t in tfs])

mw = {cs[r]: r for r in np.where(group == "Multiome_control")[0]}
tea = {cs[r]: r for r in np.where(group == "TEA_control")[0]}
ctrl_rna = dict(mw); [ctrl_rna.setdefault(cs[r], r) for r in np.where(group == "TEA_control")[0]]
aml = np.array(sorted(s for s in pd.unique(sample)
                      if group[sample == s][0] in ("patient", "AML_CCHMC") and s != EXCLUDE))
folds = [aml[i::N_FOLDS] for i in range(N_FOLDS)]

rec = []
for fold in folds:
    tr = ~np.isin(sample, list(fold))
    co, st = ed.fit(Xf[tr], Y[tr])
    rows = [r for r in np.where(np.isin(sample, list(fold)))[0] if cs[r] in ctrl_rna]
    crows = [ctrl_rna[cs[r]] for r in rows]
    ia = ed.tf_activity(ed.predict(co, st, Xf[np.array(rows)]))
    ic = ed.tf_activity(ed.predict(co, st, Xf[np.array(crows)]))
    ta = ed.tf_activity(Y[np.array(rows)])
    for k, r in enumerate(rows):
        c = cs[r]; ctrl = tea[c] if c in tea else mw[c]
        tc = ed.tf_activity(Y[ctrl:ctrl + 1])[0]
        td = ta[k] - tc; idf = ia[k] - ic[k]
        # TF own-expression differential (AML pseudobulk vs control RNA), feature space
        exprd = np.array([Xf[r, tf_feat[i]] - Xf[ctrl_rna[c], tf_feat[i]] if tf_feat[i] >= 0 else 0.0
                          for i in range(len(tfs))])
        # per-pseudobulk offset removal
        td -= td.mean(); idf -= idf.mean(); exprd_c = exprd - exprd.mean()
        for i in range(len(tfs)):
            rec.append((float(td[i]), float(idf[i]), float(exprd_c[i]), float(exprd[i])))
df = pd.DataFrame(rec, columns=["true_c", "imp_c", "expr_c", "expr_raw"])

thr = df["true_c"].abs().quantile(0.75)
sub = df[df["true_c"].abs() >= thr].copy()
def conc(s, col): s = s[s[col] != 0]; return (np.sign(s["true_c"]) == np.sign(s[col])).mean(), len(s)

# overall by magnitude stratum
fig, ax = plt.subplots(1, 2, figsize=(9.5, 4.2))
labels = ["top 50%", "top 25%", "top 10%"]; qs = [.50, .75, .90]
pe = []; ex = []
for q in qs:
    s = df[df["true_c"].abs() >= df["true_c"].abs().quantile(q)]
    pe.append(conc(s, "imp_c")[0]); ex.append(conc(s, "expr_c")[0])
x = np.arange(3); w = 0.38
ax[0].bar(x - w/2, pe, w, color=GREEN, label="per-edge regulon model")
ax[0].bar(x + w/2, ex, w, color=GREY, label="TF-expression-sign baseline")
ax[0].axhline(0.5, color=RED, lw=1.0, ls="--", label="chance")
ax[0].set_xticks(x); ax[0].set_xticklabels(labels); ax[0].set_ylim(0, 1)
ax[0].set_ylabel("directional concordance"); ax[0].set_xlabel("truly-differential TFs (|measured Δ| stratum)")
ax[0].set_title("Differential TF-activity direction:\nmodel vs TF-expression baseline (offset-removed)", fontsize=9)
ax[0].legend(fontsize=7.5, loc="lower right")
# stratify top-25% by whether TF's OWN expression is differential
stable = sub[sub["expr_raw"].abs() < EXPR_THR]; moved = sub[sub["expr_raw"].abs() >= EXPR_THR]
groups = ["TF mRNA STABLE\n(regulon-only)", "TF mRNA changed"]
pe2 = [conc(stable, "imp_c")[0], conc(moved, "imp_c")[0]]
ex2 = [conc(stable, "expr_c")[0], conc(moved, "expr_c")[0]]
ns = [len(stable), len(moved)]
x = np.arange(2)
ax[1].bar(x - w/2, pe2, w, color=GREEN, label="per-edge regulon model")
ax[1].bar(x + w/2, ex2, w, color=GREY, label="TF-expression-sign baseline")
ax[1].axhline(0.5, color=RED, lw=1.0, ls="--")
ax[1].set_xticks(x); ax[1].set_xticklabels([f"{g}\n(n={n})" for g, n in zip(groups, ns)], fontsize=8)
ax[1].set_ylim(0, 1); ax[1].set_ylabel("directional concordance")
ax[1].set_title("Regulon-only TFs: model exceeds the\nexpression baseline (mechanism check)", fontsize=9)
ax[1].legend(fontsize=7.5, loc="lower right")
fig.tight_layout(); fig.savefig(f"{R}/fig5_mechanism_expression_baseline.pdf"); fig.savefig(f"{R}/fig5_mechanism_expression_baseline.png", dpi=200); plt.close(fig)
print("wrote fig5_mechanism_expression_baseline.pdf")
print(f"top25 stable: per-edge={pe2[0]:.3f} expr={ex2[0]:.3f} (n={ns[0]}); "
      f"moved: per-edge={pe2[1]:.3f} expr={ex2[1]:.3f} (n={ns[1]})")
