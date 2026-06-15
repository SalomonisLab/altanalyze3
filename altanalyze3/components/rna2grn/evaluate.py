"""rna2grn.evaluate — canonical, reproducible evaluation of the per-edge regulon
regression for differential TF-activity inference.

This is the single entry point that produces every number quoted in
``VALIDATION.md``. It is intentionally self-contained (numpy/scipy/sklearn) and
re-implements the shipped model's fit/predict inline (matching
``model.RegulonEdgeRegressor``) so a reviewer can read the whole evaluation in one
file.

Design (matches the lab spec):
  * sample = donor; pseudobulk = donor x cell-state (the prediction unit).
  * AML-7 is EXCLUDED (reserved as the raw-.h5 real-world input).
  * Grouped 5-fold over the remaining 27 AML donors; each fold holds out several
    donors. Every AML pseudobulk is tested once while its donor is held out.
  * The merged control aggregates (Multiome_WT, RC2_TEA) are UNIQUE and CANNOT be
    held out, so they are RETAINED in training every fold (they are the control
    baseline, not test donors).
  * The control RNA is the summed counts of the 12 CCHMC healthy-control donors
    (Annotation==Control & Dataset==CCHMC), per cell state, CP10k+log1p — paired
    with BOTH the Multiome_WT and RC2_TEA control GRNs.
  * For each held-out AML pseudobulk (S,c), differential TF activity is scored
    against the cell-state-matched control, using TEA where it covers the cell
    state (closer magnitude match to the 3' AML data) and Multiome otherwise.

Run:
    python -m altanalyze3.components.rna2grn.evaluate            # default paths
    python -m altanalyze3.components.rna2grn.evaluate --dataset <dataset.npz> --out <dir>

Outputs (in <out>/):
    directional_dual_control.csv     per (pseudobulk, TF, control_source) example
    directional_metrics_dual.json    summary: dir-acc by control x stratum x raw/centered
    per_pseudobulk_topworst_tf.csv    top/worst inferred TFs per pseudobulk
    detail_per_sample.csv / detail_per_cellstate.csv / detail_per_tf.csv
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

DEF_DATASET = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn/matched/dataset.npz"
DEF_OUT = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn/benchmark"
RIDGE_LAMBDA = 1.0
N_FOLDS = 5
EXCLUDE = "AML-7_CCHMC"


# ----------------------------------------------------------------- model core
class _Edges:
    """Edge bookkeeping + per-edge ridge fit/predict (== model.RegulonEdgeRegressor)."""

    def __init__(self, genes, edge_tf, edge_gene):
        gp = {g: i for i, g in enumerate(genes)}
        self.feat = [g for g in sorted(set(edge_tf) | set(edge_gene)) if g in gp]
        fpos = {g: i for i, g in enumerate(self.feat)}
        self.feat_cols = np.array([gp[g] for g in self.feat])
        self.E = len(edge_tf)
        self.tgt = np.array([fpos[g] for g in edge_gene])
        self.tfc = np.array([fpos[g] for g in edge_tf])
        self.tf_list = sorted(set(edge_tf))
        tp = {t: i for i, t in enumerate(self.tf_list)}
        self.tf_of_edge = np.array([tp[t] for t in edge_tf])
        rows, cols, data = [], [], []
        for t in self.tf_list:
            cc = [self.tgt[e] for e in range(self.E) if edge_tf[e] == t]
            w = 1.0 / len(cc)
            for c in cc:
                rows.append(tp[t]); cols.append(c); data.append(w)
        self.regmat = sp.csr_matrix((data, (rows, cols)), shape=(len(self.tf_list), len(self.feat)))

    def features(self, Xf):  # (n, E, 3): target, TF, regulon-mean
        return np.stack([Xf[:, self.tgt], Xf[:, self.tfc],
                         np.asarray(Xf @ self.regmat.T)[:, self.tf_of_edge]], axis=2)

    def fit(self, Xf, Y, lam=RIDGE_LAMBDA):
        F = self.features(Xf); n = Xf.shape[0]; cols = [np.ones((n, self.E))]; st = []
        for j in range(3):
            A = F[:, :, j]; m = A.mean(0); sd = A.std(0); sd[sd < 1e-8] = 1
            st.append((m, sd)); cols.append((A - m) / sd)
        D = np.stack(cols, 2)
        AtA = np.einsum("nef,neg->efg", D, D); Aty = np.einsum("nef,ne->ef", D, Y)
        rg = lam * np.eye(4)[None].repeat(self.E, 0); rg[:, 0, 0] = 0
        return np.linalg.solve(AtA + rg, Aty), st

    def predict(self, coefs, st, Xf):
        F = self.features(Xf); cols = [np.ones((Xf.shape[0], self.E))]
        for j, (m, sd) in enumerate(st):
            cols.append((F[:, :, j] - m) / sd)
        return np.clip(np.einsum("nef,ef->ne", np.stack(cols, 2), coefs), 0, None)

    def tf_activity(self, edge_mat):  # (n, E) -> (n, n_tf) mean of each TF's edges
        out = np.zeros((edge_mat.shape[0], len(self.tf_list)))
        for ti in range(len(self.tf_list)):
            out[:, ti] = edge_mat[:, self.tf_of_edge == ti].mean(1)
        return out


_H5AD = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk/pseudobulk_counts_hashed.h5ad"


def _mutation_map():
    """donor -> driver-mutation annotation (from the RNA h5ad obs); {} if unavailable."""
    try:
        import h5py
        f = h5py.File(_H5AD, "r")
        g = f["obs"]["Annotation"]; cats = [x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]
        codes = g["codes"][:]
        pb = [x.decode() if isinstance(x, bytes) else x for x in f["obs"]["pseudobulk"][:]]
        f.close()
        m = {}
        for i, p in enumerate(pb):
            m.setdefault(p.split("|")[0], cats[codes[i]])
        # map merged CCHMC names (e.g. AML-7_CCHMC) by stripping suffix when needed
        return {**m, **{f"{k}_CCHMC": v for k, v in m.items()}}
    except Exception:
        return {}


# ------------------------------------------------------------------ evaluate
def run(dataset=DEF_DATASET, out=DEF_OUT):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    d = np.load(dataset, allow_pickle=True)
    X = d["X"].astype(np.float64); Y = d["Y"].astype(np.float64)
    genes = d["genes"].astype(str)
    edge_tf = d["edge_tf"].astype(str); edge_gene = d["edge_gene"].astype(str)
    sample = d["sample"].astype(str); cs = d["cell_state"].astype(str); group = d["group"].astype(str)

    ed = _Edges(genes, edge_tf, edge_gene)
    Xf = X[:, ed.feat_cols]
    tfs = np.array(ed.tf_list)

    is_mw = group == "Multiome_control"
    is_tea = group == "TEA_control"
    mw_row = {cs[r]: r for r in np.where(is_mw)[0]}
    tea_row = {cs[r]: r for r in np.where(is_tea)[0]}
    # control RNA per cell state = the (shared) CCHMC aggregate stored on the control rows
    ctrl_rna_row = {cs[r]: r for r in np.where(is_mw)[0]}
    for r in np.where(is_tea)[0]:
        ctrl_rna_row.setdefault(cs[r], r)

    aml = np.array(sorted(s for s in pd.unique(sample)
                          if group[sample == s][0] in ("patient", "AML_CCHMC") and s != EXCLUDE))
    folds = [aml[i::N_FOLDS] for i in range(N_FOLDS)]

    recs = []
    for fold in folds:
        test = set(fold)
        # controls (Multiome + TEA) RETAINED in training; only AML test donors held out
        tr = ~np.isin(sample, list(test))
        coefs, st = ed.fit(Xf[tr], Y[tr])
        rows = [r for r in np.where(np.isin(sample, list(test)))[0] if cs[r] in ctrl_rna_row]
        if not rows:
            continue
        crows = [ctrl_rna_row[cs[r]] for r in rows]
        imp_aml = ed.tf_activity(ed.predict(coefs, st, Xf[np.array(rows)]))
        imp_ctrl = ed.tf_activity(ed.predict(coefs, st, Xf[np.array(crows)]))
        true_aml = ed.tf_activity(Y[np.array(rows)])
        for k, r in enumerate(rows):
            c = cs[r]
            for src, rowmap in (("TEA", tea_row), ("Multiome", mw_row)):
                if c not in rowmap:
                    continue
                true_ctrl = ed.tf_activity(Y[rowmap[c]:rowmap[c] + 1])[0]
                td = true_aml[k] - true_ctrl
                idf = imp_aml[k] - imp_ctrl[k]
                for ti, tf in enumerate(tfs):
                    recs.append((str(sample[r]), c, src, str(tf),
                                 float(td[ti]), float(idf[ti])))
    df = pd.DataFrame(recs, columns=["sample", "cell_state", "control", "TF", "true_diff", "imp_diff"])
    df.to_csv(out / "directional_dual_control.csv", index=False)

    # ---- per-pseudobulk centering (remove global assay offset) ----
    key = df["sample"] + "|" + df["cell_state"] + "|" + df["control"]
    df["true_c"] = df["true_diff"] - df.groupby(key)["true_diff"].transform("mean")
    df["imp_c"] = df["imp_diff"] - df.groupby(key)["imp_diff"].transform("mean")

    def dir_acc(sub, tcol, pcol):
        s = sub[sub[pcol] != 0]
        return float((np.sign(s[tcol]) == np.sign(s[pcol])).mean()), int(len(s))

    metrics = {}
    for ctrl in ["TEA", "Multiome"]:
        sub_all = df[df["control"] == ctrl]
        for mode, tcol, pcol in [("raw", "true_diff", "imp_diff"), ("centered", "true_c", "imp_c")]:
            absd = sub_all[tcol].abs()
            for lab, q in [("top50", .50), ("top25", .75), ("top10", .90)]:
                s = sub_all[absd >= absd.quantile(q)]
                a, n = dir_acc(s, tcol, pcol)
                metrics[f"{ctrl}|{mode}|{lab}"] = dict(dir_acc=a, n=n)
    json.dump(metrics, open(out / "directional_metrics_dual.json", "w"), indent=2)

    # ---- per (sample, cell_state) TOP / WORST inferred TFs (TEA preferred) ----
    pref = df.sort_values("control", key=lambda s: s.map({"TEA": 0, "Multiome": 1})) \
             .drop_duplicates(["sample", "cell_state", "TF"], keep="first")
    tw = []
    for (s, c), g in pref.groupby(["sample", "cell_state"]):
        ctrl_src = g["control"].iloc[0]
        g = g.reindex(g["true_diff"].abs().sort_values(ascending=False).index)  # most differential first
        diffq = g.head(25).copy()  # truly-differential TFs for this pseudobulk
        diffq["correct"] = np.sign(diffq["true_diff"]) == np.sign(diffq["imp_diff"])
        diffq["abs_true"] = diffq["true_diff"].abs()
        best = diffq[diffq["correct"]].sort_values("abs_true", ascending=False).head(5)
        worst = diffq[~diffq["correct"]].sort_values("abs_true", ascending=False).head(5)
        for kind, gg in (("top", best), ("worst", worst)):
            for _, r in gg.iterrows():
                tw.append(dict(sample=s, cell_state=c, control=ctrl_src, kind=kind, TF=r["TF"],
                               true_diff=round(r["true_diff"], 4), imp_diff=round(r["imp_diff"], 4),
                               direction_correct=bool(r["correct"])))
    twdf = pd.DataFrame(tw)
    twdf.to_csv(out / "per_pseudobulk_topworst_tf.csv", index=False)

    # ---- detail tables (per sample / cell-state / TF), TEA-preferred, top-25% diff ----
    pref = pref.copy()
    pref["abs_true"] = pref["true_diff"].abs()
    thr = pref["abs_true"].quantile(0.75)
    dd = pref[pref["abs_true"] >= thr].copy()
    dd["correct"] = np.sign(dd["true_diff"]) == np.sign(dd["imp_diff"])
    mut = _mutation_map()
    per_s = dd.groupby("sample").agg(n_cellstates=("cell_state", "nunique"),
                                     n_diff_TF=("TF", "size"), dir_acc=("correct", "mean")).reset_index()
    per_s["mutation"] = per_s["sample"].map(mut)
    per_s.sort_values("dir_acc", ascending=False).to_csv(out / "detail_per_sample.csv", index=False)
    per_c = dd.groupby("cell_state").agg(n_samples=("sample", "nunique"),
                                         n_diff_TF=("TF", "size"), dir_acc=("correct", "mean")).reset_index()
    per_c.sort_values("dir_acc", ascending=False).to_csv(out / "detail_per_cellstate.csv", index=False)
    per_t = dd.groupby("TF").agg(n_diff=("TF", "size"), frac_up_in_AML=("true_diff", lambda x: (x > 0).mean()),
                                 dir_acc=("correct", "mean")).reset_index()
    per_t[per_t["n_diff"] >= 20].sort_values("dir_acc", ascending=False).to_csv(out / "detail_per_tf.csv", index=False)

    print("dual-control directional accuracy (correct-direction fraction):")
    for k, v in metrics.items():
        print(f"  {k:30s} dir_acc={v['dir_acc']:.3f}  n={v['n']}")
    print(f"\nexamples: {len(df)};  pseudobulks scored: {pref[['sample','cell_state']].drop_duplicates().shape[0]}")
    print(f"wrote -> {out}/directional_dual_control.csv, directional_metrics_dual.json, "
          f"per_pseudobulk_topworst_tf.csv")
    return metrics, twdf


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default=DEF_DATASET)
    ap.add_argument("--out", default=DEF_OUT)
    a = ap.parse_args()
    run(a.dataset, a.out)
