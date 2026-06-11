#!/usr/bin/env python3
"""Rigorous speed-vs-fidelity comparison: original validated UDON workflow vs
ultra-fast alternatives, on the SAME fold matrix. Assesses
  (1) selected FEATURES  (Jaccard / exact match, per stage),
  (2) NMF RANK / cluster count,
  (3) pseudobulk->cluster ASSIGNMENT (Adjusted Rand / AMI),
and runtime for each step. Reproduces the fold matrix deterministically from the
pseudobulk + saved control_annotation.tsv so it matches the reference run.

Fast alternatives:
  - feature selection: vectorized variance/intercorrelation/PCA (fast_feature_selection)
  - NMF rank: top eigenvalues via randomized SVD instead of full LA.eig on a
    samples x samples covariance (O(n^3) -> truncated)
  - clustering: sklearn NMF and KMeans-on-PCA vs nimfa SNMF (original)
The original nimfa NMF is stochastic, so we run it twice to get a stochastic-
baseline ARI and judge fast alternatives relative to that ceiling.
"""
import os, sys, time, argparse, numpy as np, pandas as pd, anndata as ad
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from sklearn.decomposition import NMF as SKNMF, PCA, TruncatedSVD
from sklearn.cluster import KMeans
from sklearn.preprocessing import scale

UDON_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, UDON_DIR)
import pseudobulk_protocol as P
import fast_feature_selection as FFS

PB = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk"
SAMPLE, CT, ANNOT = "Sample", "Hs-BM-titrated-reference-centroid", "Annotation"


def timed(fn, *a, **k):
    t = time.time(); r = fn(*a, **k); return r, time.time() - t


def jacc(a, b):
    a, b = set(a), set(b)
    return len(a & b) / max(len(a | b), 1)


def fast_rank(fold_hvg, oversample=60):
    """determine_nmf_ranks via randomized SVD: eigenvalues of scale(X)^T scale(X)
    are the squared singular values of scale(X). Count > Tracy-Widom boundary."""
    X = scale(np.asarray(fold_hvg))
    g, c = X.shape
    muTW = (np.sqrt(g - 1) + np.sqrt(c)) ** 2.0
    sigmaTW = (np.sqrt(g - 1) + np.sqrt(c)) * (1.0 / np.sqrt(g - 1) + 1.0 / np.sqrt(c)) ** (1.0 / 3.0)
    boundary = 3.273 * sigmaTW + muTW
    ncomp = min(c - 1, oversample)
    sv = TruncatedSVD(n_components=ncomp, random_state=0).fit(X).singular_values_
    k = int((sv ** 2 > boundary).sum())
    return 2 * k, boundary


def sknmf_clusters(fold_hvg, rank):
    V = fold_hvg.to_numpy()                       # genes x samples (>=0)
    model = SKNMF(n_components=int(rank), init="nndsvda", max_iter=300, random_state=0)
    model.fit(V.T)                                # samples x genes -> W (samples x rank)
    W = model.transform(V.T)
    return pd.Series(np.argmax(W, axis=1), index=fold_hvg.columns)


def kmeans_clusters(fold_hvg, rank):
    F = PCA(n_components=min(30, fold_hvg.shape[1] - 1), random_state=0).fit_transform(fold_hvg.T.to_numpy())
    km = KMeans(n_clusters=int(rank), n_init=10, random_state=0).fit(F)
    return pd.Series(km.labels_, index=fold_hvg.columns)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subset-celltypes", default="HSC-1,HSC-2,MPP-1,Classical-Mono")
    ap.add_argument("--full", action="store_true", help="use all cell types (slow original)")
    ap.add_argument("--out", default=os.path.join(PB, "UDON", "speed_comparison"))
    ap.add_argument("--skip-original-nmf", action="store_true",
                    help="skip the slow nimfa NMF; use reference udon_clusters.txt instead")
    ap.add_argument("--reference-clusters", default=None)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rep = open(os.path.join(args.out, "comparison_report.txt"), "w")
    def out(m): print(m); rep.write(m + "\n"); rep.flush()

    # ---- reconstruct fold matrix (matched, deterministic) ----
    adata = ad.read_h5ad(os.path.join(PB, "pseudobulk_scaled_log2_hashed.h5ad"))
    import re
    susp = np.array([bool(re.search("TotalSeq|empty", str(s), re.I)) or str(a) == "0"
                     for s, a in zip(adata.obs[SAMPLE], adata.obs[ANNOT])])
    adata = adata[~susp].copy()
    if not args.full:
        cts = set(args.subset_celltypes.split(","))
        adata = adata[adata.obs[CT].astype(str).isin(cts)].copy()
    sel = P.read_control_annotation(os.path.join(PB, "UDON", "control_annotation.tsv"))
    folds, fold_obs = P.build_fold_matrix(adata, SAMPLE, CT, ANNOT, mode="matched",
                                          selected_controls=sel, logger=lambda m: None)
    out(f"fold matrix: {folds.shape[0]} genes x {folds.shape[1]} pseudobulks "
        f"({'FULL' if args.full else 'subset'})")

    pc_genes = set(pd.read_csv(os.path.join(UDON_DIR, "ProteinCoding-Hs-Mm.txt"),
                               sep="\t", names=["g"])["g"])
    pc_genes = {g for g in pc_genes if str(g).isupper()}

    # ================= FEATURE SELECTION =================
    out("\n== FEATURE SELECTION ==")
    udon = P.assemble_udon_adata(folds, fold_obs)
    prev = os.getcwd(); os.chdir(UDON_DIR)
    from feature_selection import (identify_protein_coding_genes, filter_out_non_coding_genes,
                                   variance_based_feature_selection,
                                   intercorrelation_based_feature_selection, pca_feature_selection)
    (uo, t_o) = timed(lambda: _orig_fs(udon.copy(), identify_protein_coding_genes,
                                       filter_out_non_coding_genes, variance_based_feature_selection,
                                       intercorrelation_based_feature_selection, pca_feature_selection))
    os.chdir(prev)
    orig_corr = uo.var.index[uo.var["correlated_genes"]]
    orig_pca = uo.var.index[uo.var.get("pca_selected_genes", pd.Series(False, index=uo.var.index))]

    (fast, t_f) = timed(FFS.fast_feature_selection, folds, pc_genes)
    out(f"original feature selection: {t_o:.1f}s -> correlated={len(orig_corr)} pca={len(orig_pca)}")
    out(f"fast     feature selection: {t_f:.1f}s -> correlated={len(fast['correlated_genes'])} "
        f"pca={len(fast['pca_selected_genes'])}  (speedup {t_o/max(t_f,1e-6):.1f}x)")
    out(f"  correlated_genes Jaccard(original,fast) = {jacc(orig_corr, fast['correlated_genes']):.4f}")
    out(f"  pca_selected_genes Jaccard(original,fast) = {jacc(orig_pca, fast['pca_selected_genes']):.4f}")

    feat = list(orig_corr)
    fold_hvg = folds.loc[feat]

    # ================= RANK =================
    out("\n== NMF RANK ==")
    from nmf import determine_nmf_ranks, run_nmf
    (ro, tro) = timed(determine_nmf_ranks, fold_hvg)
    (rf, trf) = timed(fast_rank, fold_hvg)
    out(f"original rank (full eig): {ro} in {tro:.1f}s")
    out(f"fast rank (randomized SVD): {rf[0]} in {trf:.1f}s  (speedup {tro/max(trf,1e-6):.1f}x)")
    rank = ro if ro > 1 else 2

    # ================= CLUSTERING =================
    out("\n== CLUSTERING (rank=%d) ==" % rank)
    labels = {}
    if not args.skip_original_nmf:
        (_, c1), t1 = timed(run_nmf, fold_hvg, rank)
        (_, c2), t2 = timed(run_nmf, fold_hvg, rank)
        labels["nimfa_run1"] = c1["cluster"]; labels["nimfa_run2"] = c2["cluster"]
        out(f"nimfa SNMF run1: {t1:.1f}s, {c1['cluster'].nunique()} clusters")
        out(f"nimfa SNMF run2: {t2:.1f}s, {c2['cluster'].nunique()} clusters")
        ref = c1["cluster"]; ref_name = "nimfa_run1"
    else:
        rc = pd.read_csv(args.reference_clusters, sep="\t", index_col=0)
        ref = rc.iloc[:, 0]; ref_name = "reference"; labels[ref_name] = ref

    (sk, tsk) = timed(sknmf_clusters, fold_hvg, rank)
    (km, tkm) = timed(kmeans_clusters, fold_hvg, rank)
    labels["sklearn_NMF"] = sk; labels["kmeans_pca"] = km
    out(f"sklearn NMF: {tsk:.1f}s, {sk.nunique()} clusters")
    out(f"kmeans/PCA : {tkm:.1f}s, {km.nunique()} clusters")

    out("\n== ASSIGNMENT CONCORDANCE (ARI / AMI vs %s) ==" % ref_name)
    common = ref.index
    for name, lab in labels.items():
        if name == ref_name:
            continue
        l = lab.reindex(common).dropna()
        r = ref.reindex(l.index)
        ari = adjusted_rand_score(r, l); ami = adjusted_mutual_info_score(r, l)
        out(f"  {name:14s} vs {ref_name}: ARI={ari:.3f}  AMI={ami:.3f}  (n_clusters={lab.nunique()})")
    rep.close()
    print("\nwrote", os.path.join(args.out, "comparison_report.txt"))


def _orig_fs(udon, ipc, fonc, vbf, icf, pcf):
    udon = ipc(udon, species="Hs")
    udon = fonc(udon)
    udon = vbf(udon, fold_threshold=1, samples_differing=3)
    udon = icf(udon, corr_threshold=0.4, corr_n_events=5)
    udon = pcf(udon, corr_threshold=0.4, n_components=30)
    return udon


if __name__ == "__main__":
    main()
