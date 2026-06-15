#!/usr/bin/env python3.11
"""rna2grn diagnostics: gene coverage + variance decomposition by level."""
import numpy as np
import pandas as pd

OUT = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn"
d = np.load(f"{OUT}/matched/dataset.npz", allow_pickle=True)
X, Y = d["X"], d["Y"]
genes = d["genes"].astype(str)
edge_tf = d["edge_tf"].astype(str)
edge_gene = d["edge_gene"].astype(str)
sample = d["sample"].astype(str)
cell_state = d["cell_state"].astype(str)
group = d["group"].astype(str)

gene_pos = {g: i for i, g in enumerate(genes)}
tfs = pd.unique(edge_tf)
tgts = pd.unique(edge_gene)
tf_in = [t for t in tfs if t in gene_pos]
tgt_in = [t for t in tgts if t in gene_pos]
print("=== GENE COVERAGE in RNA (35702 genes) ===")
print(f"unique TFs: {len(tfs)}  in RNA: {len(tf_in)} ({100*len(tf_in)/len(tfs):.0f}%)")
print(f"unique targets: {len(tgts)}  in RNA: {len(tgt_in)} ({100*len(tgt_in)/len(tgts):.0f}%)")
miss_tf = [t for t in tfs if t not in gene_pos]
miss_tg = [t for t in tgts if t not in gene_pos]
print("missing TFs:", miss_tf[:20], "..." if len(miss_tf) > 20 else "")
print("n edges with both TF and Gene in RNA:",
      int(sum((t in gene_pos and g in gene_pos) for t, g in zip(edge_tf, edge_gene))), "/", len(edge_tf))

# --- variance decomposition by level (use direct rows only: real per-sample data) ---
direct = group != "Multiome_control"
direct &= group != "TEA_control"
Yd = Y[direct]
cs = cell_state[direct]
sm = sample[direct]
print(f"\n=== VARIANCE DECOMPOSITION (direct per-sample rows: n={Yd.shape[0]}) ===")
grand = Yd.mean(axis=0, keepdims=True)
total_ss = ((Yd - grand) ** 2).sum()

# cell-state means
cs_levels = pd.unique(cs)
cs_mean = {c: Yd[cs == c].mean(axis=0) for c in cs_levels}
CSmean = np.vstack([cs_mean[c] for c in cs])
between_cs_ss = ((CSmean - grand) ** 2).sum()
within_cs_ss = ((Yd - CSmean) ** 2).sum()
print(f"total SS                    : {total_ss:.3e}")
print(f"between-cell-state SS       : {between_cs_ss:.3e}  ({100*between_cs_ss/total_ss:.1f}%)  <- RNA-decodable (cell identity)")
print(f"within-cell-state SS        : {within_cs_ss:.3e}  ({100*within_cs_ss/total_ss:.1f}%)  <- sample-specific (hard)")

# how many edges have meaningful within-cell-state signal?
# per-edge: fraction of variance that is within-cellstate
edge_within = ((Yd - CSmean) ** 2).sum(axis=0)
edge_total = ((Yd - grand) ** 2).sum(axis=0) + 1e-12
frac_within = edge_within / edge_total
print(f"\nper-edge fraction-within-cellstate: median={np.median(frac_within):.2f}")
print(f"edges that are mostly cell-state-driven (<30% within): {int((frac_within<0.3).sum())}")
print(f"edges with substantial within-cs signal (>50% within): {int((frac_within>0.5).sum())}")

# profile-correlation floors
gm = grand.ravel()
prof_corr_gm = np.array([np.corrcoef(Yd[i], gm)[0, 1] for i in range(Yd.shape[0])])
prof_corr_cs = np.array([np.corrcoef(Yd[i], CSmean[i])[0, 1] for i in range(Yd.shape[0])])
print(f"\nprofile Pearson floor vs GLOBAL mean   : median={np.nanmedian(prof_corr_gm):.3f}")
print(f"profile Pearson floor vs CELLSTATE mean: median={np.nanmedian(prof_corr_cs):.3f}  <- baseline to beat")
