#!/usr/bin/env python3
"""
Integrate per-study UDON clusters into final, batch-robust UDON programs.

Pipeline (continues from study_aware_udon.py outputs):
  1. Strategy-2 marker integration: connect per-study clusters whose top-50 marker
     Jaccard >= THR; connected components = candidate conserved programs.
  2. Exclude batch-specific clusters: components present in only ONE study are
     dropped by default (EXCLUDE_SINGLE_STUDY) -> surviving clusters.
  3. SVA on the surviving clusters' centroids (combined gene space): estimate
     surrogate variables (residual components associated with STUDY after modelling
     the program) and remove them -> batch-corrected centroids. (R `sva` is not
     installed here, so this is a numpy implementation of the SVA idea.)
  4. Final clusters = cluster the SVA-corrected surviving centroids.
  5. Re-SVA + re-classify: SVA-correct the FULL combined pseudobulk fold matrix
     (study as batch) and assign EVERY pseudobulk to its nearest final program,
     then iterate once (re-SVA with the new labels) so all pseudobulks land in the
     final program set.
  6. Evaluate batch: Dataset dominance of the final programs (target: << 0.64).
"""
import os, sys, numpy as np, pandas as pd, anndata as ad
import scipy.sparse.csgraph as csgraph
from sklearn.cluster import AgglomerativeClustering

UDON_DIR = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, UDON_DIR); os.chdir(UDON_DIR)
import pseudobulk_protocol as P

PB = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk"
_CELLTYPE, _GENEFILT, _SPECIES, _SFX = P.udon_restriction()   # honours --cell-type via UDON_CELL_TYPE
SA = os.path.join(PB, "UDON", "study_aware" + _SFX)
S, CT, AN = "Sample", "Hs-BM-titrated-reference-centroid", "Annotation"
JACCARD_THR = 0.15
EXCLUDE_SINGLE_STUDY = True
SCORE_FLOOR = 0.0          # exclude pseudobulks whose best-program classification score < this (low confidence)


def regress_r2(y, X):
    X1 = np.hstack([np.ones((len(y), 1)), X])
    beta = np.linalg.lstsq(X1, y, rcond=None)[0]
    yhat = X1 @ beta
    ss = ((y - y.mean()) ** 2).sum() + 1e-12
    return 1 - ((y - yhat) ** 2).sum() / ss


def sva_remove(Y, bio_onehot, batch_labels, r2_thr=0.3, log=print):
    """Y: features x samples. Remove surrogate variables (sample-space residual
    components associated with batch) while keeping the bio model. Returns cleaned Y."""
    Xb = bio_onehot.astype(float)                                   # samples x p
    B = np.linalg.lstsq(Xb, Y.T, rcond=None)[0]
    R = Y.T - Xb @ B                                                # residuals (samples x features)
    U, Sd, Vt = np.linalg.svd(R, full_matrices=False)              # U: samples x k
    batch_oh = pd.get_dummies(pd.Series(batch_labels)).values.astype(float)
    keep = [c for c in range(U.shape[1]) if regress_r2(U[:, c], batch_oh) > r2_thr]
    if not keep:
        keep = list(range(min(len(set(batch_labels)) - 1, U.shape[1])))
    SV = U[:, keep]                                                # samples x n_sv
    XX = np.hstack([Xb, SV])
    coef = np.linalg.lstsq(XX, Y.T, rcond=None)[0]
    Yclean = (Y.T - SV @ coef[Xb.shape[1]:]).T                     # features x samples
    log(f"    SVA: removed {len(keep)} surrogate variable(s) associated with study")
    return Yclean


def main():
    rep = open(os.path.join(SA, "integration_report.txt"), "w")
    def log(m): print(m); rep.write(str(m) + "\n"); rep.flush()

    C = pd.read_csv(os.path.join(SA, "per_study_cluster_centroids.tsv"), sep="\t", index_col=0)
    # conserved vs unique are determined by direct top-50 gene-LIST overlap (study_aware_genelist.py),
    # shared-core floor; no GO/ontology.
    cvu = pd.read_csv(os.path.join(SA, "conserved_vs_unique.tsv"), sep="\t").set_index("cluster_key").reindex(C.index)
    keys = list(C.index); study = cvu["study"].astype(str).values
    genes = list(C.columns)
    log(f"[1] {len(keys)} per-study clusters; {len(genes)} centroid genes")

    # 2. surviving = gene-list CONSERVED clusters (shared core, >= 2 studies);
    #    drop 'unique' (single-study) clusters by default.
    surv = (cvu["status"] == "conserved").values if EXCLUDE_SINGLE_STUDY else np.ones(len(keys), bool)
    log(f"[2] gene-list conserved/unique: dropped {int((~surv).sum())} unique (single-study) clusters "
        f"-> {int(surv.sum())} surviving conserved clusters")
    n_groups = len(set(cvu["genelist_group"].values[surv])) if surv.any() else 0
    if int(surv.sum()) < 2 or n_groups < 2:
        log(f"[2] ABORT: only {int(surv.sum())} conserved cross-study cluster(s) in {n_groups} program(s). "
            f"Study-aware integration needs >=2 conserved programs -- the per-study data is too thin for "
            f"cross-study conservation here (typical when restricted to a single cell type). "
            f"Use the naive run instead (UDON/celltype_<cell type>/).")
        rep.close(); return
    Cs = C.iloc[np.where(surv)[0]]
    comp_s = cvu["genelist_group"].values[surv]; study_s = study[surv]
    prelim = pd.get_dummies(pd.Series(comp_s)).values             # bio model = gene-list conserved group

    # 3. SVA on surviving centroids (genes x clusters)
    log("[3] SVA on surviving-cluster centroids (study = batch)...")
    Yc = sva_remove(Cs.values.T, prelim, study_s, log=log)        # genes x clusters
    Cs_clean = pd.DataFrame(Yc.T, index=Cs.index, columns=genes)

    # 4. final clusters = agglomerate SVA-cleaned surviving centroids
    n_final = int(len(set(comp_s)))                                # one per conserved program
    fl = AgglomerativeClustering(n_clusters=n_final, metric="cosine", linkage="average"
                                 ).fit_predict(Cs_clean.values)
    log(f"[4] final UDON programs: {n_final} (from SVA-corrected surviving centroids)")
    final_centroids = pd.DataFrame({f"P{p}": Cs_clean.values[fl == p].mean(0) for p in range(n_final)},
                                   index=genes).T                  # programs x genes

    # 5. re-SVA + classify ALL pseudobulks
    log("[5] reconstructing combined fold matrix + re-SVA classifying all pseudobulks...")
    sel = P.read_control_annotation(os.path.join(PB, "UDON", "matched_sex_platform", "control_annotation.tsv"))
    a = ad.read_h5ad(os.path.join(PB, "pseudobulk_scaled_log2_hashed.h5ad"))
    import re
    sus = np.array([bool(re.search("TotalSeq|empty", str(s), re.I)) or str(x) == "0"
                    for s, x in zip(a.obs[S], a.obs[AN])]); a = a[~sus].copy()
    folds, fobs = P.build_fold_matrix(a, S, CT, AN, mode="matched", selected_controls=sel, logger=lambda m: None)
    if _CELLTYPE:                                                  # restrict reclassification to the cell type
        m = fobs[CT].astype(str) == _CELLTYPE
        folds = folds.loc[:, fobs.index[m]]; fobs = fobs[m]
        log(f"[restriction] reclassifying only '{_CELLTYPE}': {folds.shape[1]} pseudobulks")
    folds = folds.loc[genes]                                       # common (gene-filtered) gene space
    SM = pd.read_csv(os.path.join(PB, "evaluation/outputs/sample_metadata.tsv"), sep="\t", index_col=0)
    pb_study = SM.reindex(fobs[S].astype(str))["Dataset"].astype(str).values

    def assign(Fmat):
        # nearest final-program centroid by Pearson correlation
        Fz = Fmat - Fmat.mean(0, keepdims=True)
        Cz = final_centroids.values - final_centroids.values.mean(1, keepdims=True)
        Fz /= (np.linalg.norm(Fz, axis=0, keepdims=True) + 1e-12)
        Czn = Cz / (np.linalg.norm(Cz, axis=1, keepdims=True) + 1e-12)
        cor = Czn @ Fz                                             # programs x pseudobulks
        return np.argmax(cor, axis=0)

    Fmat = folds.values                                           # genes x pseudobulks
    lab = assign(Fmat); score = np.zeros(Fmat.shape[1])
    for it in range(2):                                           # re-SVA iterations
        bio = pd.get_dummies(pd.Series(lab)).values
        Fclean = sva_remove(Fmat, bio, pb_study, log=log)
        # refresh final centroids in cleaned space
        global_centroids = np.vstack([Fclean[:, lab == p].mean(1) for p in range(n_final)])
        Fz = Fclean - Fclean.mean(0, keepdims=True); Fz /= (np.linalg.norm(Fz, axis=0, keepdims=True) + 1e-12)
        Cz = global_centroids - global_centroids.mean(1, keepdims=True)
        Cz /= (np.linalg.norm(Cz, axis=1, keepdims=True) + 1e-12)
        cor = Cz @ Fz                                             # programs x pseudobulks
        lab = np.argmax(cor, axis=0); score = np.max(cor, axis=0)   # best-program correlation = score
        log(f"    re-SVA iter {it+1}: classified {len(lab)} pseudobulks into {len(set(lab))} programs")

    # exclude low-confidence classifications: best-program score < SCORE_FLOOR (analogous to SVM score<0)
    keep = score >= SCORE_FLOOR
    log(f"[5b] excluded {int((~keep).sum())} pseudobulks with classification score < {SCORE_FLOOR} "
        f"(low confidence); {int(keep.sum())} retained")
    out = pd.DataFrame({"pseudobulk": np.asarray(folds.columns)[keep], "Sample": fobs[S].values[keep],
                        "Dataset": np.asarray(pb_study)[keep], "final_program": [f"P{x}" for x in lab[keep]],
                        "score": np.round(score[keep], 3)})
    out.to_csv(os.path.join(SA, "final_program_assignments.tsv"), sep="\t", index=False)

    # 6. evaluate batch: Dataset dominance of final programs
    doms = []
    log("\n[6] FINAL programs -- Dataset dominance:")
    for p, g in out.groupby("final_program"):
        vd = g["Dataset"].value_counts(normalize=True); doms.append(vd.iloc[0])
        log(f"   {p:4s} n={len(g):5d}  top={vd.index[0]:9s} {vd.iloc[0]:.2f}  n_studies={g['Dataset'].nunique()}")
    log(f"\nMEAN dominant-Dataset fraction (final, study-aware+SVA): {np.mean(doms):.2f}")
    log("  vs combined matched-only: 0.64 ; null baseline ~0.49")
    rep.close()


if __name__ == "__main__":
    main()
