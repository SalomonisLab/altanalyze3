#!/usr/bin/env python3
"""Contrast the gene-filtered study-aware run (study_aware_protein_coding) against the original
(study_aware, no gene filter): final programs, batch dominance, ribosomal content, donor
specificity, and SATAY-UDON covariate enrichments."""
import os, sys, numpy as np, pandas as pd

PB = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk"
OLD = os.path.join(PB, "UDON", "study_aware")
NEW = os.path.join(PB, "UDON", "study_aware_protein_coding")


def programs(d, label):
    fa = pd.read_csv(os.path.join(d, "final_program_assignments.tsv"), sep="\t")
    doms = [g["Dataset"].value_counts(normalize=True).iloc[0] for _, g in fa.groupby("final_program")]
    print(f"{label:18s}: {fa['final_program'].nunique():2d} programs | {len(fa):5d} pseudobulks | "
          f"mean Dataset-dominance {np.mean(doms):.2f} (lower=less batch)")
    return fa


def ribo(d, label):
    p = os.path.join(d, "final_program_markers.txt")
    if not os.path.exists(p):
        return
    mk = pd.read_csv(p, sep="\t")
    n = mk["marker"].astype(str).str.upper().str.startswith(("RPL", "RPS")).sum()
    xist = mk["marker"].astype(str).str.upper().isin(["XIST", "TSIX"]).sum()
    print(f"{label:18s}: {n} ribosomal (RPL/RPS) + {xist} XIST/TSIX marker genes among final programs")


def covars(d, label, fdr=0.05):
    f = os.path.join(d, "SATAY-UDON", "satay_celltype_resolved.tsv")
    if not os.path.exists(f):
        print(f"{label:18s}: no SATAY-UDON table"); return None
    df = pd.read_csv(f, sep="\t")
    sig = df[df["fdr"] < fdr]
    print(f"{label:18s}: {len(sig)} sig celltype-resolved tests (FDR<{fdr}) | "
          f"{sig['covariate'].nunique()} covariates | {sig['celltype'].nunique()} cell types")
    return sig


print("=" * 78)
print("STUDY-AWARE: gene-filtered (protein-coding, no RPL/RPS/XIST/TSIX, Y kept) vs original")
print("=" * 78)
print("\n--- final programs + batch dominance ---")
old_fa = programs(OLD, "ORIGINAL"); new_fa = programs(NEW, "GENE-FILTERED")
print("\n--- ribosomal / XIST content of the final programs ---")
ribo(OLD, "ORIGINAL"); ribo(NEW, "GENE-FILTERED")

print("\n--- SATAY-UDON covariate enrichments ---")
old_sig = covars(OLD, "ORIGINAL"); new_sig = covars(NEW, "GENE-FILTERED")
if old_sig is not None and new_sig is not None:
    oc, nc = set(old_sig["covariate"]), set(new_sig["covariate"])
    print(f"\ncovariates gained (NEW only): {sorted(nc - oc)}")
    print(f"covariates lost   (OLD only): {sorted(oc - nc)}")
    print(f"covariates shared: {len(oc & nc)}")
    # top mutation associations in each
    for lab, sig in [("ORIGINAL", old_sig), ("GENE-FILTERED", new_sig)]:
        m = sig[sig["covariate"].astype(str).str.contains("|".join(
            ["NPM1", "FLT3", "inv(16)", "t(8;21)", "RUNX1", "CEBPA", "TP53", "NRAS", "KIT"]), regex=True)]
        print(f"\n{lab} top mutation-enriched (cluster:celltype:covariate, FDR):")
        for r in m.sort_values("fdr").head(8).itertuples():
            print(f"   {r.cluster}:{r.celltype}:{r.covariate}  FDR={r.fdr:.1e}")
