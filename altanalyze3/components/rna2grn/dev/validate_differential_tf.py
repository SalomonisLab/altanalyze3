#!/usr/bin/env python3.11
"""User's validation: for a HELD-OUT AML sample, the TF activities truly increased
vs control (from the real GRN) must also be recovered in the IMPUTED GRN.

We compute, per cell state, the true and imputed differential TF-activity vectors
(AML - Multiome_WT control), then measure: Spearman(true,imputed) over TFs, and
recovery of the truly top-increased TFs (precision@5/@10). Held out = AML-7 (and
all AML samples via LOSO bundles)."""
import sys, time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
from altanalyze3.components.rna2grn import training, Rna2GrnBundle

DATA = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn"
npz = f"{DATA}/matched/dataset.npz"
d = np.load(npz, allow_pickle=True)
X = d["X"]; Y = d["Y"]; genes = d["genes"].astype(str)
edge_tf = d["edge_tf"].astype(str)
sample = d["sample"].astype(str); cs = d["cell_state"].astype(str); group = d["group"].astype(str)
tfs = pd.unique(edge_tf)
ctrl_rows = {cs[r]: r for r in np.where(group == "Multiome_control")[0]}


def true_tf_activity(yvec):
    return np.array([yvec[edge_tf == t].mean() for t in tfs])


def imputed_tf_activity(pred_row):
    return np.array([pred_row[edge_tf == t].mean() for t in tfs])


aml_samples = [s for s in pd.unique(sample) if not group[sample == s][0].endswith("control")]
recs = []
t0 = time.time()
for s in aml_samples:
    # build a bundle excluding this sample (true held-out)
    b = training.build_bundle(npz, "/tmp/_loso_bundle.pkl.gz", exclude_samples=[s],
                              include_controls=True)
    bundle = Rna2GrnBundle.load("/tmp/_loso_bundle.pkl.gz")
    rows = [r for r in np.where(sample == s)[0] if cs[r] in ctrl_rows]
    if not rows:
        continue
    crows = [ctrl_rows[cs[r]] for r in rows]
    Q = pd.DataFrame(X[rows], index=[cs[r] for r in rows], columns=genes)
    C = pd.DataFrame(X[crows], index=[cs[r] for r in rows], columns=genes)
    qp = bundle.predict_from_dataframe(Q, normalized=True).predictions.to_numpy()
    cp = bundle.predict_from_dataframe(C, normalized=True).predictions.to_numpy()
    for k, (r, cr) in enumerate(zip(rows, crows)):
        dt = true_tf_activity(Y[r]) - true_tf_activity(Y[cr])      # true differential TF activity
        di = imputed_tf_activity(qp[k]) - imputed_tf_activity(cp[k])
        rho = spearmanr(dt, di).correlation
        # recovery of truly top-increased TFs
        true_top = set(np.array(tfs)[np.argsort(-dt)[:10]])
        imp_top = set(np.array(tfs)[np.argsort(-di)[:10]])
        p10 = len(true_top & imp_top) / 10
        recs.append(dict(sample=s, cell_state=cs[r], tf_spearman=rho, precision_at10=p10))
res = pd.DataFrame(recs)
res.to_csv(f"{DATA}/benchmark/differential_tf_validation.csv", index=False)
print(f"=== imputed vs true DIFFERENTIAL TF activity (held-out AML vs control), {time.time()-t0:.0f}s ===")
print(f"  pseudobulks evaluated: {len(res)}")
print(f"  Spearman(true,imputed) over 217 TFs: median={res.tf_spearman.median():.3f} mean={res.tf_spearman.mean():.3f}")
print(f"  precision@10 (top-increased TFs):    median={res.precision_at10.median():.2f} mean={res.precision_at10.mean():.2f}")
print("\nper held-out sample (median across its cell states):")
g = res.groupby("sample").agg(n=("cell_state", "size"), tf_spearman=("tf_spearman", "median"),
                              p10=("precision_at10", "median")).sort_values("tf_spearman", ascending=False)
print(g.to_string())
