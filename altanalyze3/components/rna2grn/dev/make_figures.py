#!/usr/bin/env python3.11
"""rna2grn evaluation figures for the per-edge regulon-activity model (editable PDF, Arial)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
from altanalyze3.components.rna2grn import Rna2GrnBundle

OUT = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn"
R = f"{OUT}/reports"; B = f"{OUT}/benchmark"
BLUE = "#1F77B4"; RED = "#D62728"; GREEN = "#2CA02C"; ORANGE = "#FF7F0E"; GREY = "#999999"

# ---------- Fig 1: PERTURBATION test (per-edge vs kNN) ----------
pt = pd.read_csv(f"{B}/perturbation_tf_ranks.csv")
fig, ax = plt.subplots(1, 2, figsize=(9, 4))
# ECDF of perturbed-TF rank
for col, c, lab in [("pe_rank_med", GREEN, "per-edge regulon regression"),
                    ("knn_rank_med", RED, "kNN retrieval")]:
    v = np.sort(pt[col].values); y = np.arange(1, len(v) + 1) / len(v)
    ax[0].plot(v, y, color=c, lw=1.8, label=lab)
ax[0].set_xlabel("rank of the PERTURBED TF (1 = perfect, of 217)")
ax[0].set_ylabel("fraction of perturbations ≤ rank")
ax[0].set_title("In-silico regulon perturbation:\nonly the regression localizes the TF", fontsize=9)
ax[0].legend(fontsize=8, loc="lower right"); ax[0].set_xlim(0, 120)
top = pd.DataFrame({
    "per-edge": [(pt.pe_rank_med == 1).mean(), (pt.pe_rank_med <= 5).mean()],
    "kNN": [(pt.knn_rank_med == 1).mean(), (pt.knn_rank_med <= 5).mean()]}, index=["top-1", "top-5"])
x = np.arange(2); w = 0.35
ax[1].bar(x - w/2, top["per-edge"], w, color=GREEN, label="per-edge")
ax[1].bar(x + w/2, top["kNN"], w, color=RED, label="kNN")
ax[1].set_xticks(x); ax[1].set_xticklabels(top.index); ax[1].set_ylim(0, 1)
ax[1].set_ylabel("fraction of perturbed TFs detected")
ax[1].set_title("Perturbed-TF detection rate", fontsize=9); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{R}/fig1_perturbation_detection.pdf"); fig.savefig(f"{R}/fig1_perturbation_detection.png", dpi=200); plt.close(fig)
print("wrote fig1_perturbation_detection.pdf")

# ---------- Fig 2: variance decomposition ----------
fig, ax = plt.subplots(figsize=(4.5, 4))
ax.bar([0], [65.5], color=BLUE, label="between cell-state (65.5%)\nRNA-decodable")
ax.bar([0], [34.5], bottom=[65.5], color=ORANGE, label="within cell-state (34.5%)\nsample-specific (hard)")
ax.set_xlim(-1, 1); ax.set_xticks([]); ax.set_ylabel("% of GRN variance")
ax.set_title("GRN variance decomposition", fontsize=9)
ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5))
fig.tight_layout(); fig.savefig(f"{R}/fig2_variance_decomposition.pdf"); fig.savefig(f"{R}/fig2_variance_decomposition.png", dpi=200); plt.close(fig)
print("wrote fig2_variance_decomposition.pdf")

# ---------- Fig 3: protocol similarity ----------
ps = pd.read_csv(f"{R}/protocol_similarity_multiome_vs_tea.csv")
fig, ax = plt.subplots(figsize=(6.5, 4))
xx = np.arange(len(ps)); ax.bar(xx, ps["pearson"], color=BLUE)
ax.axhline(0.399, color=RED, lw=1.2, ls="--", label="within-Multiome across DIFFERENT cell states (0.40)")
ax.axhline(ps["pearson"].mean(), color=GREEN, lw=1.2, label=f"mean cross-protocol same cell state ({ps['pearson'].mean():.2f})")
ax.set_xticks(xx); ax.set_xticklabels(ps["cell_state"], rotation=45, ha="right", fontsize=7)
ax.set_ylabel("Pearson r (Multiome vs TEA GRN)"); ax.set_ylim(0, 1)
ax.set_title("Cross-protocol GRN similarity", fontsize=9); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f"{R}/fig3_protocol_similarity.pdf"); fig.savefig(f"{R}/fig3_protocol_similarity.png", dpi=200); plt.close(fig)
print("wrote fig3_protocol_similarity.pdf")

# ---------- Fig 4: directional concordance of differential TF activity (CANONICAL) ----------
# Built from evaluate.py dual-control outputs. AML-7 EXCLUDED; controls RETAINED in training;
# directional concordance reported raw and after per-pseudobulk offset removal.
import json
m = json.load(open(f"{B}/directional_metrics_dual.json"))
per_s = pd.read_csv(f"{B}/detail_per_sample.csv")
ex = pd.read_csv(f"{B}/directional_dual_control.csv")
strata = ["top50", "top25", "top10"]; xlab = ["top 50%", "top 25%", "top 10%"]
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
# Panel A: directional concordance by stratum, TEA vs Multiome, raw vs offset-removed
x = np.arange(3); w = 0.2
bars = [("TEA|raw", GREEN, -1.5), ("TEA|centered", "#7FCF7F", -0.5),
        ("Multiome|raw", BLUE, 0.5), ("Multiome|centered", "#9EC9E8", 1.5)]
for key, col, off in bars:
    ax[0].bar(x + off * w, [m[f"{key.split('|')[0]}|{key.split('|')[1]}|{s}"]["dir_acc"] for s in strata],
              w, color=col, label=key.replace("|", " "))
ax[0].axhline(0.5, color=RED, lw=1.0, ls="--", label="chance")
ax[0].set_xticks(x); ax[0].set_xticklabels(xlab); ax[0].set_ylim(0, 1)
ax[0].set_ylabel("directional concordance"); ax[0].set_xlabel("truly-differential TFs (|measured Δ| stratum)")
ax[0].set_title("Imputed vs measured differential TF activity\n(held-out AML vs matched control; AML-7 excluded)", fontsize=9)
ax[0].legend(fontsize=6.5, ncol=1, loc="lower right")
# Panel B: per-donor concordance (top-25%, TEA-preferred), sorted
ps = per_s.sort_values("dir_acc")
ax[1].barh(np.arange(len(ps)), ps["dir_acc"], color=GREEN)
ax[1].axvline(ps["dir_acc"].median(), color=RED, lw=1.2, label=f"median {ps['dir_acc'].median():.2f}")
ax[1].set_yticks(np.arange(len(ps))); ax[1].set_yticklabels(ps["sample"], fontsize=5)
ax[1].set_xlabel("directional concordance"); ax[1].set_xlim(0, 1)
ax[1].set_title(f"Per held-out donor (n={len(ps)})\nrange {ps['dir_acc'].min():.2f}–{ps['dir_acc'].max():.2f}", fontsize=9)
ax[1].legend(fontsize=7, loc="lower right")
# Panel C: representative held-out pseudobulk (TEA control), offset-removed differential
tea = ex[ex["control"] == "TEA"].copy()
counts = tea.groupby(["sample", "cell_state"]).size()
s_sel, c_sel = counts.idxmax()  # most-populated TEA pseudobulk
sub = tea[(tea["sample"] == s_sel) & (tea["cell_state"] == c_sel)].copy()
sub["tc"] = sub["true_diff"] - sub["true_diff"].mean()
sub["ic"] = sub["imp_diff"] - sub["imp_diff"].mean()
conc = np.sign(sub["tc"]) == np.sign(sub["ic"])
ax[2].scatter(sub["tc"][conc], sub["ic"][conc], s=16, color=GREEN, alpha=0.7, edgecolors="none", label="concordant sign")
ax[2].scatter(sub["tc"][~conc], sub["ic"][~conc], s=16, color=RED, alpha=0.7, edgecolors="none", label="discordant")
ax[2].axhline(0, color="#000000", lw=0.6); ax[2].axvline(0, color="#000000", lw=0.6)
topi = sub["tc"].abs().sort_values(ascending=False).index[:8]
for i in topi:
    ax[2].annotate(sub.loc[i, "TF"], (sub.loc[i, "tc"], sub.loc[i, "ic"]), fontsize=6.5, xytext=(2, 2), textcoords="offset points")
ax[2].set_xlabel("measured Δ TF activity (offset-removed)"); ax[2].set_ylabel("imputed Δ (offset-removed)")
ax[2].set_title(f"Representative pseudobulk\n{s_sel} | {c_sel} (TEA control)", fontsize=9)
ax[2].legend(fontsize=7, loc="lower right")
fig.tight_layout(); fig.savefig(f"{R}/fig4_directional_concordance.pdf"); fig.savefig(f"{R}/fig4_directional_concordance.png", dpi=200); plt.close(fig)
print("wrote fig4_directional_concordance.pdf")
print("done")
