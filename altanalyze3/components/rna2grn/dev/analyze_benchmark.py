#!/usr/bin/env python3.11
"""Post-process LOSO predictions with HONEST baselines + within-cell-state skill.

Reads pred_*.npy saved by benchmark_loso.py and recomputes:
  - skill vs the honest LOSO cell-state-mean baseline
  - within-cell-state residual metrics (the hard, sample-specific signal):
    residualize truth and prediction by the LOSO cell-state mean, then correlate
"""
import glob
import os
import numpy as np
import pandas as pd

OUT = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn"
B = f"{OUT}/benchmark"
d = np.load(f"{OUT}/matched/dataset.npz", allow_pickle=True)
group = d["group"].astype(str)
direct = ~np.isin(group, ["Multiome_control", "TEA_control"])
cell_state = d["cell_state"].astype(str)[direct]

Y = np.load(f"{B}/_y_true.npy")
preds = {os.path.basename(p)[5:-4]: np.load(p) for p in sorted(glob.glob(f"{B}/pred_*.npy"))}
CSm = preds["cellstate_mean"].astype(np.float64)  # honest LOSO cell-state mean


def col_pearson(A, B_):
    Az = A - A.mean(0, keepdims=True); Bz = B_ - B_.mean(0, keepdims=True)
    num = (Az * Bz).sum(0); den = np.sqrt((Az ** 2).sum(0) * (Bz ** 2).sum(0))
    out = np.full(A.shape[1], np.nan); ok = den > 1e-12; out[ok] = num[ok] / den[ok]; return out


def row_pearson(A, B_):
    Az = A - A.mean(1, keepdims=True); Bz = B_ - B_.mean(1, keepdims=True)
    num = (Az * Bz).sum(1); den = np.sqrt((Az ** 2).sum(1) * (Bz ** 2).sum(1))
    out = np.full(A.shape[0], np.nan); ok = den > 1e-12; out[ok] = num[ok] / den[ok]; return out


Yr = Y - CSm  # truth residual after removing honest cell-state mean
ss_cs = ((Y - CSm) ** 2).sum()
rows = []
for name, P in preds.items():
    P = P.astype(np.float64)
    ss_res = ((Y - P) ** 2).sum()
    Pr = P - CSm
    # within-cell-state per-edge correlation (only edges/rows with residual variance)
    ed_w = col_pearson(Yr, Pr)
    rp_w = row_pearson(Yr, Pr)
    rows.append(dict(
        model=name,
        profile_r=np.nanmedian(row_pearson(Y, P)),
        edge_r=np.nanmedian(col_pearson(Y, P)),
        skill_vs_CSmean=1 - ss_res / ss_cs,         # >0 beats honest cell-state mean
        within_cs_edge_r=np.nanmedian(ed_w),         # hard signal, per-edge
        within_cs_edge_rpos=float(np.nanmean(ed_w > 0)),
        within_cs_profile_r=np.nanmedian(rp_w),
    ))
res = pd.DataFrame(rows).sort_values("skill_vs_CSmean", ascending=False)
res.to_csv(f"{B}/loso_honest_summary.csv", index=False)
pd.set_option("display.width", 200)
print(res.to_string(index=False))
print(f"\nwrote {B}/loso_honest_summary.csv")
print("\nNOTE: skill_vs_CSmean>0 => beats the honest leave-one-sample-out cell-state-mean baseline.")
print("within_cs_* = predictive power for the HARD sample-specific (within-cell-state) signal.")
