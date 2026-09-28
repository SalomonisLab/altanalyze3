#!/usr/bin/env python3
"""Build a compact, auditable LungMAP hypothesis and validation package.

This script deliberately distinguishes independent measurements (Xenium expression,
spatial enrichment and donor composition) from CellRef2 RNA-derived predictions.
"""
from __future__ import annotations

import json
import math
import re
import shutil
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact, mannwhitneyu, norm


OUT = Path("/Users/saljh8/Dropbox/LungMAP/Discovery/codex_nature_predictions_20260927")
SRC = Path("/Users/saljh8/Dropbox/LungMAP/Discovery/hypotheses_20260927")
TFDIR = Path("/Users/saljh8/Dropbox/LungMAP/Discovery/runs/mechanism-20260926/tracks/M00-tf-outliers/steps/s01_primary_tf_outliers")
COMP_PATH = Path("/Users/saljh8/Dropbox/LungMAP/Discovery/evaluation/composition_atlas/hypotheses.csv")
XDIR = Path("/Users/saljh8/Dropbox/LungMAP/Xenium-IPF/inputs")
SPATIAL_PATH = Path("/Users/saljh8/Dropbox/LungMAP/Discovery/codex_spatial_pairs.tsv")

TABS, FIGS = OUT / "tables", OUT / "figures"
for d in (OUT, TABS, FIGS):
    d.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({"font.family": "Arial", "font.size": 9, "axes.titlesize": 11,
                     "axes.labelsize": 9, "figure.dpi": 130, "pdf.fonttype": 42,
                     "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
sns.set_palette("colorblind")
COL = {"Unaffected": "#6b7280", "Less Affected": "#e69f00", "More Affected": "#b2182b"}


def bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    x = p[ok]
    if not len(x):
        return out
    order = np.argsort(x)
    ranked = x[order]
    q = np.minimum.accumulate((ranked * len(x) / np.arange(1, len(x) + 1))[::-1])[::-1]
    z = np.empty_like(q)
    z[order] = np.minimum(q, 1)
    out[np.where(ok)[0]] = z
    return out


def hedges_g(case, ctrl):
    a, b = np.asarray(case, float), np.asarray(ctrl, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    sd = math.sqrt(((len(a)-1)*a.var(ddof=1) + (len(b)-1)*b.var(ddof=1)) / (len(a)+len(b)-2))
    if not sd:
        return np.nan
    d = (a.mean()-b.mean())/sd
    return d * (1 - 3/(4*(len(a)+len(b))-9))


def savefig(fig, stem):
    fig.savefig(FIGS / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGS / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def condition_rows(name):
    return pd.read_csv(SRC / f"rows_{name}.csv.gz", low_memory=False)


def selected_rows(df, population, features, modalities=None):
    x = df[(df.population == population) & df.gene.isin(features)].copy()
    if modalities is not None:
        x = x[x.modality.isin(modalities)]
    return x


# ---------- CellRef2 replicated evidence ----------
ipf = condition_rows("IPF")
inf = condition_rows("infection")
copd = condition_rows("COPD")

ipf_features = {
    "grn_tf": ["TEAD2", "ZNF322", "CREB3L1"],
    "rna": ["MMP7", "ABCA3", "SFTPC", "COL1A1", "FN1", "TEAD2", "ITGA3", "ITGB6"],
    "fastcomm": ["Alveolar macrophage|FN1->ITGAV+ITGB6",
                 "Alveolar macrophage (lipid homeostatic)|FN1->ITGAV+ITGB8",
                 "Alveolar macrophage (lipid homeostatic)|FN1->ITGAV+ITGB1",
                 "Metallothionein+ AM|SPP1->ITGAV"],
}
ev = []
for mod, genes in ipf_features.items():
    ev.append(selected_rows(ipf, "Aberrant basal__vs__AT2", genes, [mod]))
ipf_ev = pd.concat(ev, ignore_index=True)
ipf_ev["evidence_class"] = np.where(ipf_ev.modality.eq("rna"), "measured RNA differential",
                                      "RNA-derived prediction")
ipf_ev.to_csv(TABS / "ipf_multimodal_evidence.tsv", sep="\t", index=False)

# Which predicted TF targets are also independently replicated RNA changes?
grn_path = Path("/Users/saljh8/Dropbox/LungMAP/GRN/TF_to_Gene_connection_scores_log10-NOT_ordered_clusters_ALL_GENES-hybrid-r0.33-f0.05-cap50.txt")
grn = pd.read_csv(grn_path, sep="\t")
scr = pd.read_csv(SRC / "screen_IPF.csv.gz")
rep_up = set(scr[(scr.modality.eq("rna")) & (scr.population.eq("Aberrant basal__vs__AT2")) &
                 scr.replicated & (scr.direction.eq(1))].gene)
target_rows = []
for tf in ["TEAD2", "ZNF322"]:
    g = grn[grn.TF.eq(tf)].copy()
    score_cols = [c for c in g.columns if c not in ["TF", "Gene"]]
    g["max_connection_score"] = g[score_cols].max(axis=1)
    g = g[g[score_cols].gt(0).any(axis=1)]
    for r in g[g.Gene.isin(rep_up)].itertuples():
        target_rows.append({"tf": tf, "target": r.Gene, "max_connection_score": r.max_connection_score,
                            "target_is_replicated_up_RNA": True,
                            "mechanistic_module": "adhesion/cytoskeleton" if r.Gene in {"ITGB4","CDH3","JUP","ITGA3","HSPG2","LAMA5","CDC42","CTNND1","PTK7","FARP1","AFAP1","CD2AP","GIPC1","VEZT"} else "other"})
target_overlap = pd.DataFrame(target_rows).sort_values(["tf", "mechanistic_module", "max_connection_score"], ascending=[True,True,False])
target_overlap.to_csv(TABS / "ipf_tf_replicated_target_overlap.tsv", sep="\t", index=False)
target_counts = target_overlap.groupby("tf").size().to_dict()

inf_features = {
    "grn_tf": ["ETV1", "ETV5", "CEBPA", "CREB3L1", "XBP1", "MLXIPL", "CREB1", "FOS", "ELF1", "CAMTA1", "CREB3"],
    "rna": ["SFTPA1", "SFTPC", "SFTPA2", "HHIP", "PTGFR", "CA2", "HLF", "MFSD2A"],
    "adt": ["CD66c", "CD55", "Podoplanin", "CD49f"],
    "lipid": ["PC(20:4/22:6)", "PC(20:4/20:4)", "PE(P-16:0/20:4)"],
    "fastcomm": ["AT1|HMGB1->AGER", "AT2|HMGB1->AGER", "AT2-AT1 int.|CADM1->CADM1"],
}
parts = []
for mod, genes in inf_features.items():
    # Lipids define the intermediate relative to AT1; regulatory/surfactant features define it relative to AT2.
    pop = "AT2-AT1 int.__vs__AT1" if mod == "lipid" else "AT2-AT1 int.__vs__AT2"
    parts.append(selected_rows(inf, pop, genes, [mod]))
inf_ev = pd.concat(parts, ignore_index=True)
inf_ev["reference_anchor"] = np.where(inf_ev.population.str.endswith("__vs__AT1"), "AT1", "AT2")
inf_ev["evidence_class"] = np.where(inf_ev.modality.eq("rna"), "measured RNA differential",
                                     "RNA-derived prediction")
inf_ev.to_csv(TABS / "infection_multimodal_evidence.tsv", sep="\t", index=False)


# ---------- donor-pseudobulk TF effects ----------
cohort = pd.read_csv(TFDIR / "cohort_effects.tsv", sep="\t")
pooled = pd.read_csv(TFDIR / "pooled.tsv", sep="\t")
tf_keys = [("IPF", "AT2", "ARNT2"), ("COPD", "AT2", "TP63"), ("COPD", "AT2", "FOXQ1")]
tf_c = pd.concat([cohort[(cohort.condition == c) & (cohort.population == p) & (cohort.tf == t)]
                  for c, p, t in tf_keys], ignore_index=True)
tf_p = pd.concat([pooled[(pooled.condition == c) & (pooled.population == p) & (pooled.tf == t)]
                  for c, p, t in tf_keys], ignore_index=True)
tf_c.to_csv(TABS / "tf_candidate_cohort_effects.tsv", sep="\t", index=False)
tf_p.to_csv(TABS / "tf_candidate_pooled_effects.tsv", sep="\t", index=False)

# COPD pulmonary venous endothelial state, replicated stored program loss.
pv_genes = ["NFIC", "KLF5", "FOXA2", "ELF3", "SREBF2", "TEAD3", "TEAD2", "NFIX", "ATF7"]
pv = selected_rows(copd, "PVEC", pv_genes, ["grn_tf"])
pv.to_csv(TABS / "copd_pvec_regulatory_loss.tsv", sep="\t", index=False)


# ---------- independent Xenium validation ----------
meta = pd.read_csv(XDIR / "pseudobulk_metadata_final_CT_by_patient.tsv", sep="\t")
expr0 = pd.read_csv(XDIR / "pseudobulk_avgExpr_final_CT_by_patient-clean.txt", sep="\t", index_col=0)
expr = expr0.T
common = meta.pseudobulk_id.astype(str).isin(expr.index)
meta = meta.loc[common].copy()
expr = expr.loc[meta.pseudobulk_id.astype(str)].copy()
expr.index = meta.index

tests = {
    "KRT5-/KRT17+": ["ITGAV", "ITGB6", "ITGA3", "MMP7", "KRT17", "SFTPC", "TP63", "XBP1", "FN1"],
    "SPP1+ Macrophages": ["FN1", "SPP1"],
    "Activated Fibrotic FBs": ["FN1", "CTHRC1"],
    "Transitional AT2": ["MMP7", "ITGB6", "SFTPC"],
    "AT2": ["MMP7", "KRT17", "SFTPC"],
}
xrows = []
for ct, genes in tests.items():
    ix = meta.final_CT.eq(ct)
    for gene in genes:
        if gene not in expr:
            continue
        for grp in ["Unaffected", "Less Affected", "More Affected"]:
            vals = expr.loc[ix & meta.sample_affect.eq(grp), gene].astype(float)
            xrows.append({"cell_type": ct, "gene": gene, "group": grp, "n_pseudobulks": len(vals),
                          "median_log2_cp10k1": vals.median(), "mean_log2_cp10k1": vals.mean(),
                          "q25": vals.quantile(.25), "q75": vals.quantile(.75)})
        ctrl = expr.loc[ix & meta.sample_affect.eq("Unaffected"), gene].astype(float)
        dis = expr.loc[ix & meta.sample_affect.ne("Unaffected"), gene].astype(float)
        p = mannwhitneyu(dis, ctrl, alternative="two-sided").pvalue if len(dis) and len(ctrl) else np.nan
        xrows.append({"cell_type": ct, "gene": gene, "group": "Affected_vs_Unaffected_test",
                      "n_pseudobulks": len(dis), "median_log2_cp10k1": dis.median(),
                      "mean_log2_cp10k1": dis.mean(), "q25": hedges_g(dis, ctrl), "q75": p})
xstats = pd.DataFrame(xrows)
testmask = xstats.group.eq("Affected_vs_Unaffected_test")
xstats["fdr_within_selected_tests"] = np.nan
xstats.loc[testmask, "fdr_within_selected_tests"] = bh(xstats.loc[testmask, "q75"])
xstats = xstats.rename(columns={"q25": "q25_or_hedges_g_for_test", "q75": "q75_or_p_for_test"})
xstats.to_csv(TABS / "xenium_expression_stats.tsv", sep="\t", index=False)

# Composition directly from all 1.63M Xenium cell annotations, donor as analysis unit.
h5 = ad.read_h5ad(XDIR / "GSE250346_Xenium_counts_final_CT.h5ad", backed="r")
obs = h5.obs[["patient", "sample_affect", "final_CT"]].copy()
h5.file.close()
counts = obs.groupby(["patient", "sample_affect", "final_CT"], observed=True).size().rename("n").reset_index()
tot = obs.groupby(["patient", "sample_affect"], observed=True).size().rename("total").reset_index()
ct_interest = ["KRT5-/KRT17+", "SPP1+ Macrophages", "Activated Fibrotic FBs", "Transitional AT2", "AT2"]
# Region-resolved values preserve less/more-affected tissue for visualization.
grid = tot.assign(_key=1).merge(pd.DataFrame({"final_CT": ct_interest, "_key": 1}), on="_key").drop(columns="_key")
grid = grid.merge(counts[["patient", "sample_affect", "final_CT", "n"]],
                  on=["patient", "sample_affect", "final_CT"], how="left")
grid["n"] = grid.n.fillna(0)
grid["fraction"] = grid.n / grid.total
grid["logit_fraction"] = np.log((grid.n + .5) / (grid.total - grid.n + .5))
grid.to_csv(TABS / "xenium_donor_region_composition.tsv", sep="\t", index=False)

# Disease-vs-control inference collapses all sampled regions to one value per donor.
dc = obs.groupby(["patient", "final_CT"], observed=True).size().rename("n").reset_index()
dt = obs.groupby("patient", observed=True).size().rename("total").reset_index()
disease = obs.groupby("patient", observed=True).sample_affect.apply(
    lambda x: "Unaffected" if set(x.astype(str)) == {"Unaffected"} else "Affected").rename("group").reset_index()
dgrid = dt.assign(_key=1).merge(pd.DataFrame({"final_CT": ct_interest, "_key": 1}), on="_key").drop(columns="_key")
dgrid = dgrid.merge(dc, on=["patient", "final_CT"], how="left").merge(disease, on="patient")
dgrid["n"] = dgrid.n.fillna(0); dgrid["fraction"] = dgrid.n/dgrid.total
dgrid["logit_fraction"] = np.log((dgrid.n+.5)/(dgrid.total-dgrid.n+.5))
dgrid.to_csv(TABS / "xenium_donor_composition.tsv", sep="\t", index=False)
cstats = []
for ct in ct_interest:
    z = dgrid[dgrid.final_CT.eq(ct)]
    ctrl = z[z.group.eq("Unaffected")].logit_fraction
    dis = z[z.group.eq("Affected")].logit_fraction
    cstats.append({"cell_type": ct, "n_affected_donors": len(dis), "n_control_donors": len(ctrl),
                   "median_fraction_affected": z[z.group.eq("Affected")].fraction.median(),
                   "median_fraction_control": z[z.group.eq("Unaffected")].fraction.median(),
                   "hedges_g_logit": hedges_g(dis, ctrl),
                   "mannwhitney_p": mannwhitneyu(dis, ctrl).pvalue})
cstats = pd.DataFrame(cstats)
cstats["fdr_selected_states"] = bh(cstats.mannwhitney_p)
cstats.to_csv(TABS / "xenium_composition_stats.tsv", sep="\t", index=False)

# Squidpy pair z scores; section level and donor-region aggregate sensitivity.
sp = pd.read_csv(SPATIAL_PATH, sep="\t")
sp["group"] = sp.group.str.replace("_", " ", regex=False)
def donor_from_sample(s):
    return re.sub(r"(?:LA\d*|MA\d*|[AB])$", "", str(s))
sp["donor_id"] = sp["sample"].map(donor_from_sample)
sp.to_csv(TABS / "xenium_spatial_pairs_section_level.tsv", sep="\t", index=False)
spagg = sp.groupby(["group", "donor_id", "cell_a", "cell_b"], as_index=False).squidpy_z.mean()
spagg.to_csv(TABS / "xenium_spatial_pairs_donor_region.tsv", sep="\t", index=False)
presence = sp.assign(present=sp.squidpy_z.notna()).groupby(["group", "cell_a", "cell_b"]).present.agg(["sum", "count"]).reset_index()
presence.to_csv(TABS / "xenium_spatial_pair_presence.tsv", sep="\t", index=False)


# ---------- Figures ----------
# Figure 1: convergent IPF circuit.
fig, ax = plt.subplots(2, 2, figsize=(12, 9))
fig.suptitle("IPF: independent cohorts converge on an FN1–integrin–TEAD epithelial circuit", fontsize=15, fontweight="bold")

# A: composition atlas + Xenium donor composition.
comp = pd.read_csv(COMP_PATH)
c = comp[(comp.disease == "IPF") & (comp.population == "Aberrant basal")].iloc[0]
z = grid[grid.final_CT.eq("KRT5-/KRT17+")].copy()
sns.stripplot(data=z, x="sample_affect", y="fraction", order=list(COL), palette=COL, size=6, jitter=.18, ax=ax[0,0])
sns.boxplot(data=z, x="sample_affect", y="fraction", order=list(COL), color="white", width=.5, showfliers=False, ax=ax[0,0])
ax[0,0].set_yscale("symlog", linthresh=1e-5)
ax[0,0].set_title("A  Independent Xenium donor composition")
ax[0,0].set_xlabel(""); ax[0,0].set_ylabel("KRT5−/KRT17+ fraction")
ax[0,0].text(.02, .98, f"CellRef2: {c.cohorts_enriched:.0f} cohorts, {c.fold_of_mean_fraction:.1f}× mean fraction\nXenium: each point = donor",
             transform=ax[0,0].transAxes, va="top", fontsize=8)

# B: replicated CellRef2 cohorts.
plot = ipf_ev[ipf_ev.gene.isin(["TEAD2", "ZNF322", "CREB3L1", "MMP7", "ABCA3",
                                "Alveolar macrophage|FN1->ITGAV+ITGB6",
                                "Alveolar macrophage (lipid homeostatic)|FN1->ITGAV+ITGB8",
                                "Metallothionein+ AM|SPP1->ITGAV"])].copy()
plot["feature"] = plot.gene.replace({"Alveolar macrophage|FN1->ITGAV+ITGB6":"AM FN1→αvβ6",
                                     "Alveolar macrophage (lipid homeostatic)|FN1->ITGAV+ITGB8":"AM-lipid FN1→αvβ8",
                                     "Metallothionein+ AM|SPP1->ITGAV":"MT+AM SPP1→αv"})
plot["cohort"] = plot.contrast_id.str.split("__").str[0].replace({"Jaiswal2026_ILD":"Jaiswal 2026", "Adams2020":"Adams 2020"})
keep_order = ["TEAD2","ZNF322","CREB3L1","MMP7","ABCA3","AM FN1→αvβ6","AM-lipid FN1→αvβ8","MT+AM SPP1→αv"]
sns.scatterplot(data=plot, x="log2fc", y="feature", hue="cohort", style="modality", s=80, ax=ax[0,1])
ax[0,1].axvline(0, color="black", lw=.7); ax[0,1].set_title("B  Replicated CellRef2 differential evidence")
ax[0,1].set_xlabel("stored effect (log2FC or score difference)"); ax[0,1].set_ylabel("")
ax[0,1].legend(fontsize=7, loc="lower right")

# C: independent expression in predicted senders/receiver.
long = []
for ct, gene in [("SPP1+ Macrophages","FN1"),("Activated Fibrotic FBs","FN1"),
                 ("KRT5-/KRT17+","ITGB6"),("KRT5-/KRT17+","MMP7")]:
    ix = meta.final_CT.eq(ct)
    q = pd.DataFrame({"value": expr.loc[ix, gene], "group": meta.loc[ix, "sample_affect"]})
    q["feature"] = {("SPP1+ Macrophages","FN1"):"FN1 | SPP1+ macro",
                    ("Activated Fibrotic FBs","FN1"):"FN1 | fibrotic FB",
                    ("KRT5-/KRT17+","ITGB6"):"ITGB6 | KRT state",
                    ("KRT5-/KRT17+","MMP7"):"MMP7 | KRT state"}[(ct,gene)]
    long.append(q)
long = pd.concat(long)
sns.boxplot(data=long, x="feature", y="value", hue="group", hue_order=list(COL), palette=COL, showfliers=False, ax=ax[1,0])
sns.stripplot(data=long, x="feature", y="value", hue="group", hue_order=list(COL), palette=COL, dodge=True, size=3, alpha=.65, ax=ax[1,0])
handles, labels = ax[1,0].get_legend_handles_labels(); ax[1,0].legend(handles[:3], labels[:3], fontsize=7)
ax[1,0].set_title("C  Xenium expression validates sender/receiver molecules")
ax[1,0].set_xlabel(""); ax[1,0].set_ylabel("log2(CP10k + 1)"); ax[1,0].tick_params(axis="x", rotation=20)

# D: independent spatial organization.
spa = spagg[(spagg.cell_a == "KRT5-_KRT17+") & spagg.cell_b.isin(["Activated_Fibrotic_FBs","SPP1+_Macrophages","Transitional_AT2"])].copy()
spa["pair"] = spa.cell_b.replace({"Activated_Fibrotic_FBs":"KRT state ↔ fibrotic FB",
                                   "SPP1+_Macrophages":"KRT state ↔ SPP1+ macro",
                                   "Transitional_AT2":"KRT state ↔ transitional AT2"})
sns.boxplot(data=spa, x="pair", y="squidpy_z", hue="group", hue_order=list(COL), palette=COL, showfliers=False, ax=ax[1,1])
sns.stripplot(data=spa, x="pair", y="squidpy_z", hue="group", hue_order=list(COL), palette=COL, dodge=True, size=3.5, alpha=.7, ax=ax[1,1])
handles, labels = ax[1,1].get_legend_handles_labels(); ax[1,1].legend(handles[:3], labels[:3], fontsize=7)
ax[1,1].axhline(0, color="black", lw=.6); ax[1,1].set_title("D  Xenium spatial enrichment (donor-region aggregate)")
ax[1,1].set_xlabel(""); ax[1,1].set_ylabel("Squidpy neighborhood z score"); ax[1,1].tick_params(axis="x", rotation=20)
fig.tight_layout(rect=[0,0,1,.96]); savefig(fig, "Figure1_IPF_FN1_integrin_TEAD_niche")

# Figure 2: shared acute-injury checkpoint.
use = inf_ev.copy()
use = use[use.gene.isin(["ETV5","XBP1","CEBPA","MLXIPL","SFTPA1","SFTPC","MFSD2A","CD66c","Podoplanin",
                         "PC(20:4/22:6)","PC(20:4/20:4)","PE(P-16:0/20:4)"])]
use["cohort"] = use.contrast_id.str.extract(r"__(COVID19|Pneumonia)_")[0]
use["feature_label"] = use.modality + "-" + use.gene + " | vs " + use.reference_anchor
mat = use.pivot_table(index="feature_label", columns="cohort", values="log2fc", aggfunc="mean")
fig, ax = plt.subplots(figsize=(7, 7))
sns.heatmap(mat, cmap="RdBu_r", center=0, annot=True, fmt=".2f", linewidths=.4, cbar_kws={"label":"stored effect"}, ax=ax)
ax.set_title("Acute lung injury: a two-anchor transitional-state checkpoint\nregulatory/surfactant loss vs AT2; PUFA-phospholipid loss vs AT1")
ax.set_xlabel(""); ax.set_ylabel("modality-feature | reference anchor")
fig.tight_layout(); savefig(fig, "Figure2_infection_surfactant_lipid_checkpoint")

# Figure 3: TF candidates, with uncertainty.
fig, ax = plt.subplots(1, 2, figsize=(11, 4.8))
order = ["ARNT2 | IPF AT2", "TP63 | COPD AT2", "FOXQ1 | COPD AT2"]
for i, ((cond,pop,tf), label) in enumerate(zip(tf_keys, order)):
    p = tf_p[(tf_p.condition==cond)&(tf_p.population==pop)&(tf_p.tf==tf)].iloc[0]
    ax[0].errorbar(p.pooled_g, i, xerr=[[p.pooled_g-p.ci_low],[p.ci_high-p.pooled_g]], fmt="o", color="#0072b2", capsize=3)
    cs = tf_c[(tf_c.condition==cond)&(tf_c.population==pop)&(tf_c.tf==tf)]
    ax[0].scatter(cs.hedges_g, np.repeat(i,len(cs)) + np.linspace(-.12,.12,len(cs)), color="#777777", s=18, zorder=3)
ax[0].axvline(0,color="black",lw=.7); ax[0].set_yticks(range(3), order); ax[0].invert_yaxis()
ax[0].set_xlabel("Hedges g (case − control)"); ax[0].set_title("A  Donor-pseudobulk TF program effects\npooled 95% CI; gray = cohorts")
pvplot = pv.groupby("gene", as_index=False).agg(mean_log2fc=("log2fc","mean"), max_fdr=("fdr","max"), cohorts=("contrast_id","nunique")).sort_values("mean_log2fc")
sns.barplot(data=pvplot, x="mean_log2fc", y="gene", color="#56b4e9", ax=ax[1])
ax[1].axvline(0,color="black",lw=.7); ax[1].set_xlabel("mean stored GRN score difference"); ax[1].set_ylabel("")
ax[1].set_title("B  COPD pulmonary venous endothelial\nreplicated regulatory loss")
fig.tight_layout(); savefig(fig, "Figure3_TF_candidates_and_COPD_endothelium")

# Figure 4: disease-specific enriched states, only replicated positive associations.
cc = comp[(comp.cohorts_enriched >= 2) & (comp.max_g > 0) & comp.disease.isin(["IPF","COPD"])].copy()
cc = cc.sort_values(["disease","max_g"], ascending=[True,False])
fig, axes = plt.subplots(1, 2, figsize=(11, max(5, .28*max(cc.groupby('disease').size()))), sharex=True)
for a, disease in zip(axes, ["IPF","COPD"]):
    d = cc[cc.disease.eq(disease)].head(18).sort_values("max_g")
    a.hlines(d.population, 0, d.max_g, color="#cccccc")
    a.scatter(d.max_g, d.population, s=25+18*d.cohorts_enriched, c=np.log10(d.fold_of_mean_fraction.clip(lower=1)), cmap="viridis")
    a.axvline(0,color="black",lw=.5); a.set_title(f"{disease}: replicated enriched states")
    a.set_xlabel("maximum donor-level Hedges g")
axes[0].set_ylabel(""); axes[1].set_ylabel("")
fig.suptitle("Disease-enriched cell states: positive, independently replicated signals only\npoint size = cohorts; color = log10 fold of mean fraction", fontsize=13, fontweight="bold")
fig.tight_layout(rect=[0,0,1,.92]); savefig(fig, "Figure4_replicated_disease_enriched_states")


# ---------- ranked hypotheses and narrative ----------
hyp = pd.DataFrame([
 {"rank":1,"hypothesis":"Macrophage/fibroblast FN1–αv integrin signaling stabilizes a TEAD2/ZNF322 adhesion program in IPF KRT5−/KRT17+ aberrant epithelium.",
  "confidence":"high-priority prediction","independent_validation":"Xenium expression, donor composition, and spatial neighborhood enrichment",
  "decisive_experiment":"Healthy/IPF AT2 organoids + SPP1 macrophages or FN1 matrix ± fibrotic fibroblasts; factorial CRISPRi TEAD2/ZNF322 and αv/β6 blockade; quantify KRT17/MMP7/ITGB6 versus SFTPC/ABCA3 and repair."},
 {"rank":2,"hypothesis":"Loss of an ETV5–XBP1 alveolar checkpoint marks an AT2-to-AT1 intermediate with surfactant loss versus AT2 and PUFA-phospholipid loss versus AT1 across acute viral and bacterial injury.",
  "confidence":"replicated multimodal prediction; no independent lipid assay here","independent_validation":"COVID-19 and pneumonia arms agree; modalities are partially RNA-derived",
  "decisive_experiment":"In infected human alveolar organoids, CRISPRa ETV5 or inducible XBP1s; targeted lipidomics for PC(20:4/22:6), PC(20:4/20:4), PE(P-16:0/20:4), surfactant biophysics, and AT1 differentiation."},
 {"rank":3,"hypothesis":"COPD AT2 cells acquire TP63/FOXQ1 basal-secretory regulatory activity before full state conversion; pulmonary venous endothelium independently loses FOXA2/SREBF2/KLF5 identity programs.",
  "confidence":"lower-tier; TF activity is a regulon score and target support is incomplete","independent_validation":"three donor cohorts for AT2 TF scores; two stored cohorts for PVEC losses",
  "decisive_experiment":"COPD/control distal epithelial organoids and endothelial chips; perturb TP63/FOXQ1 or FOXA2/SREBF2, then test lineage markers, barrier function and smoke-response persistence."},
])
hyp.to_csv(TABS / "ranked_hypotheses.tsv", sep="\t", index=False)

def _xt(ct, gene):
    return xstats[(xstats.cell_type.eq(ct)) & (xstats.gene.eq(gene)) &
                  (xstats.group.eq("Affected_vs_Unaffected_test"))].iloc[0]

fn1_macro, fn1_fb = _xt("SPP1+ Macrophages", "FN1"), _xt("Activated Fibrotic FBs", "FN1")
itgb6_krt = _xt("KRT5-/KRT17+", "ITGB6")
krt_comp = cstats[cstats.cell_type.eq("KRT5-/KRT17+")].iloc[0]
spatial_medians = (spagg[(spagg.cell_a.eq("KRT5-_KRT17+")) &
                         (spagg.cell_b.isin(["Activated_Fibrotic_FBs", "SPP1+_Macrophages"]))]
                   .groupby(["group", "cell_b"]).squidpy_z.median().unstack())

report = f"""# LungMAP mechanistic predictions with independent Xenium validation

## Executive result

The lead result is a **testable IPF circuit**, not merely another cell-frequency observation: macrophage/fibroblast **FN1–αv-integrin signaling** is predicted to stabilize a **TEAD2/ZNF322 adhesion program** in KRT5−/KRT17+ aberrant epithelium. TEAD2 and ZNF322 activity, epithelial injury markers, and multiple FN1/integrin communication scores agree across independent CellRef2 cohorts. An independent 1.63-million-cell Xenium cohort then places the corresponding molecules in the predicted senders and receiver and shows disease-region spatial organization of the cells.

This is a directional mechanistic prediction, not causal proof. The broad aberrant-basaloid/fibroblast/SPP1-macrophage niche has prior literature support; the specific **FN1–αv → TEAD2/ZNF322** circuit and its factorial perturbation are the new prediction nominated here.

## 1. Lead hypothesis: FN1–integrin–TEAD/ZNF322 IPF circuit

### Convergent evidence

- CellRef2 identifies aberrant basal cells as IPF-enriched in two cohorts (~{c.fold_of_mean_fraction:.1f}× mean case/control fraction; minimum q={c.min_q:.2g}).
- In `Aberrant basal__vs__AT2`, TEAD2 and ZNF322 GRN activities rise and CREB3L1 falls in Adams 2020 and Jaiswal 2026. MMP7 rises while ABCA3 falls. This contrast mixes cell identity and disease and is therefore interpreted as a state signature, not a pure IPF-within-cell-state effect.
- Predicted macrophage FN1→αvβ6/αvβ8/αvβ1 and SPP1→αv inputs rise in the same state comparison. These communication and GRN modalities are RNA-derived predictions.
- The GRN-to-RNA coherence check finds **{target_counts.get('TEAD2', 0)} TEAD2** and **{target_counts.get('ZNF322', 0)} ZNF322** predicted targets that are themselves replicated-up RNA features. The TEAD2 overlap includes an adhesion/basement-membrane module (ITGB4, CDH3, JUP, ITGA3, HSPG2, LAMA5, CTNND1, PTK7 and CDC42), providing molecular endpoints for the perturbation rather than relying on a one-dimensional TF score.
- Independent Xenium pseudobulks directly measure disease-associated FN1 in SPP1+ macrophages (Hedges g={fn1_macro.q25_or_hedges_g_for_test:.2f}, selected-test FDR={fn1_macro.fdr_within_selected_tests:.2g}) and activated fibrotic fibroblasts (g={fn1_fb.q25_or_hedges_g_for_test:.2f}, FDR={fn1_fb.fdr_within_selected_tests:.2g}). KRT-state ITGB6 is also higher (g={itgb6_krt.q25_or_hedges_g_for_test:.2f}, FDR={itgb6_krt.fdr_within_selected_tests:.2g}; only three unaffected pseudobulks).
- Across all 1.63 million Xenium cells, the KRT state is independently disease-enriched at the donor level (26 affected versus 9 control donors; logit-fraction Hedges g={krt_comp.hedges_g_logit:.2f}, FDR={krt_comp.fdr_selected_states:.2g}). Donor-region median neighborhood z scores for KRT↔fibrotic-FB are {spatial_medians.loc['Unaffected','Activated_Fibrotic_FBs']:.2f} unaffected, {spatial_medians.loc['Less Affected','Activated_Fibrotic_FBs']:.2f} less affected and {spatial_medians.loc['More Affected','Activated_Fibrotic_FBs']:.2f} more affected; KRT↔SPP1-macrophage scores are {spatial_medians.loc['Unaffected','SPP1+_Macrophages']:.2f}, {spatial_medians.loc['Less Affected','SPP1+_Macrophages']:.2f} and {spatial_medians.loc['More Affected','SPP1+_Macrophages']:.2f}, respectively.

### Falsifiable experiment

Use healthy and IPF distal-lung AT2 organoids in a factorial design: control matrix versus FN1-rich matrix; control versus SPP1 macrophages; ± fibrotic fibroblasts; ± αv/β6 blockade; CRISPRi TEAD2, ZNF322, both, or non-targeting guide. Primary endpoints are KRT17, MMP7, ITGB6/ITGA3, SFTPC/ABCA3, clonogenic repair, epithelial polarity and matrix invasion. The model is supported only if FN1/macrophage exposure induces the state and integrin blockade or TEAD2/ZNF322 knockdown prevents it; epistasis should place TEAD2/ZNF322 downstream of αv-integrin input.

## 2. Acute-injury hypothesis: ETV5/XBP1 checkpoint links transition to lipid failure

COVID-19 and non-COVID pneumonia define the same AT2→AT1 intermediate from two reference anchors. Relative to AT2, ETV5, XBP1, CEBPA and MLXIPL activity falls with SFTPA1/SFTPA2/SFTPC and MFSD2A. Relative to AT1, the intermediate loses PC(20:4/22:6), PC(20:4/20:4) and PE(P-16:0/20:4). These are separate contrasts and are not treated as one effect estimate. The lipid and surface-protein values are predicted from RNA, so the experiment must use targeted lipidomics and surfactant function as truly orthogonal readouts. CRISPRa ETV5 or inducible XBP1s should rescue the exact nominated lipids and surfactant biophysics without preventing appropriate AT1 differentiation.

## 3. Lower-tier COPD hypotheses

TP63 and FOXQ1 AT2 regulon scores are positive in three COPD cohorts, but regulons were learned in other epithelial states and replicated AT2 target-RNA support is incomplete. Treat this as a basal/secretory priming reporter until perturbation demonstrates causality. In parallel, pulmonary venous endothelial cells reproducibly lose FOXA2, SREBF2, KLF5, ELF3, TEAD2/3 and related regulatory scores in two cohorts, nominating an endothelial identity/barrier experiment.

## What was excluded or downgraded

- Imputed ADT, lipid, communication and TF activities are not independent omics; all are labeled RNA-derived.
- Missing stored differential rows mean raw p≥0.05 or untested, not evidence of no change.
- The aberrant-basal-versus-AT2 comparison is identity-mixed.
- The ionocyte positive control did **not** pass the replicated positive disease-enrichment screen and is not presented as robust.
- ARNT2 is a highly reproducible IPF AT2 program score, but its state-specific regulon provenance prevents calling ARNT2 itself a validated driver.
- Spatial z scores are association, not ligand directionality. Section-level values and donor-region aggregates are both supplied.

## Primary outputs

- `figures/Figure1_IPF_FN1_integrin_TEAD_niche.pdf`: lead convergent result.
- `tables/ipf_multimodal_evidence.tsv`: exact CellRef2 rows behind the lead model.
- `tables/xenium_expression_stats.tsv`, `xenium_composition_stats.tsv`, and `xenium_spatial_pairs_*`: independent validation and analysis units.
- `tables/ranked_hypotheses.tsv`: concise experimental handoff.

## Literature boundary

Prior studies establish aberrant basaloid cells, SPP1 macrophages, CTHRC1 fibroblasts and their spatial niches in fibrosis (PMIDs 32832599, 31221805, 34969962, 36807143, 39121212, 40759629). A targeted literature check did not identify the exact FN1–αv→TEAD2/ZNF322 circuit; that absence is not proof of novelty and should be confirmed in a formal systematic search before manuscript claims.
"""
(OUT / "REPORT.md").write_text(report)

methods = """# Methods and audit notes

CellRef2 stored donor-pseudobulk differentials were screened at FDR ≤0.10. RNA, ADT and lipid features additionally required |log2FC| ≥ log2(1.2); GRN-TF and communication scores had no fold cutoff. Replication required at least two donor-independent cohorts with concordant direction and no significant disagreement. Only stored rows were available; absent rows are not zeros.

TF activity is the sum of inferred target-edge activities in donor pseudobulks. Random-effects pooled Hedges g and cohort effects are inherited from the audited M00 TF-outlier analysis. These are program scores, not TF protein measurements.

Xenium expression used the provided patient/cell-state pseudobulk matrix, log2(CP10k+1). Selected affected-versus-unaffected comparisons use two-sided Mann–Whitney tests and Hedges g; BH FDR applies only across the selected tests in the output table and is descriptive, not a discovery-wide FDR. The KRT state has only three unaffected pseudobulks, so its expression statistics are validation signals with wide uncertainty.

Xenium cell composition was recomputed from all 1,630,319 cell annotations. Patient is the disease-versus-control analysis unit: less- and more-affected samples from the same patient are collapsed before testing. Fractions use logit((n+0.5)/(N−n+0.5)) for effect sizes. Region-resolved values are retained for visualization. Spatial results are precomputed Squidpy neighborhood z scores. Figures use averages across sections within donor and affected-region category; raw section-level scores are also retained.

No causal language is warranted without perturbation. ADT, lipid, GRN and communication layers in CellRef2 are RNA-derived; only Xenium transcript expression, cell counts and spatial organization are independent measurements used here.
"""
(OUT / "METHODS.md").write_text(methods)

shutil.copy2(Path(__file__), OUT / "build_nature_predictions.py")
manifest = {"created":"2026-09-27", "output":str(OUT), "input_files":[str(SRC/"rows_IPF.csv.gz"), str(SRC/"rows_infection.csv.gz"), str(SRC/"rows_COPD.csv.gz"), str(TFDIR/"cohort_effects.tsv"), str(TFDIR/"pooled.tsv"), str(COMP_PATH), str(XDIR/"GSE250346_Xenium_counts_final_CT.h5ad"), str(XDIR/"pseudobulk_avgExpr_final_CT_by_patient-clean.txt"), str(SPATIAL_PATH)], "notes":"All derived tables and figures reproducible with build_nature_predictions.py"}
(OUT / "manifest.json").write_text(json.dumps(manifest, indent=2))
print(json.dumps({"out":str(OUT), "tables":len(list(TABS.glob('*'))), "figures":len(list(FIGS.glob('*')))}, indent=2))
