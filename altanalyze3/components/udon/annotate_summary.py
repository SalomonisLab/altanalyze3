#!/usr/bin/env python3
"""Annotate the study-aware SATAY summary with, per program P#:
  - UDON_cluster_merged     : dominant merged matched_sex_platform U#
  - UDON_clusters_perstudy  : constituent per-study UDON clusters (top, >=40% share)
  - Top10_markers           : top-10 DIFFERENTIAL markers (every program gets markers,
                              including catch-alls -- ranked by program-vs-rest specificity)
  - TopEnrichedPathway      : top GO-Elite pathway on the program's top-50 markers (annotation only)
Pure Python. No R."""
import os, sys, re, numpy as np, pandas as pd, anndata as ad
UDON_DIR = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, UDON_DIR); os.chdir(UDON_DIR)
sys.path.insert(0, os.path.abspath(os.path.join(UDON_DIR, "..", "..", "..")))
import pseudobulk_protocol as P
from study_aware_integrate import sva_remove
from satay_heatmap import compute_pmap

PB = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk"
SA = os.path.join(PB, "UDON", "study_aware")
SAT = os.path.join(SA, "SATAY-UDON")          # consolidated with the other SATAY-UDON outputs
S, CT, AN = "Sample", "Hs-BM-titrated-reference-centroid", "Annotation"


def main():
    fa = pd.read_csv(os.path.join(SA, "final_program_assignments.tsv"), sep="\t")
    common = list(pd.read_csv(os.path.join(SA, "per_study_cluster_centroids.tsv"), sep="\t", index_col=0, nrows=1).columns)

    # --- reconstruct SVA-corrected fold matrix (same as the integration) ---
    print("reconstructing SVA-corrected expression for program markers ...")
    sel = P.read_control_annotation(os.path.join(PB, "UDON/matched_sex_platform/control_annotation.tsv"))
    a = ad.read_h5ad(os.path.join(PB, "pseudobulk_scaled_log2_hashed.h5ad"))
    sus = np.array([bool(re.search("TotalSeq|empty", str(s), re.I)) or str(x) == "0"
                    for s, x in zip(a.obs[S], a.obs[AN])]); a = a[~sus].copy()
    folds, fobs = P.build_fold_matrix(a, S, CT, AN, mode="matched", selected_controls=sel, logger=lambda m: None)
    folds = folds.loc[common]
    fam = fa[fa["pseudobulk"].isin(folds.columns)].copy()
    folds = folds[fam["pseudobulk"].tolist()]
    prog = fam["final_program"].values; study = fam["Dataset"].astype(str).values
    Fc = pd.DataFrame(sva_remove(folds.values, pd.get_dummies(pd.Series(prog)).values, study, log=lambda m: None),
                      index=common, columns=folds.columns)

    # --- per-program top differential markers (EVERY program, no winner-take-all threshold) ---
    progs = sorted(set(prog), key=lambda x: int(x[1:]))
    # markers = genes whose expression CORRELATES with program membership (markerFinder-style),
    # computed per program so EVERY program (incl. catch-alls) gets its top genes, and the
    # program-SPECIFIC identity genes (RUNX1T1, HOXB, ...) are recovered.
    X = Fc.values.astype(float); Xc = X - X.mean(axis=1, keepdims=True)
    Xnorm = np.linalg.norm(Xc, axis=1) + 1e-9
    top10, top50 = {}, {}
    for p in progs:
        ind = (prog == p).astype(float); ind = ind - ind.mean()
        corr = (Xc @ ind) / (Xnorm * (np.linalg.norm(ind) + 1e-9))
        gs = [common[i] for i in np.argsort(corr)[::-1]]
        top10[p] = ",".join(gs[:10]); top50[p] = gs[:50]
    pd.DataFrame({"program": progs, "top10_markers": [top10[p] for p in progs]}).to_csv(
        os.path.join(SA, "final_program_markers_all.txt"), sep="\t", index=False)

    # --- merged U# + per-study constituents ---
    pmap = compute_pmap()
    import glob
    ps = {}
    for d in glob.glob(os.path.join(SA, "per_study/*/udon_clusters.txt")):
        st = os.path.basename(os.path.dirname(d)); c = pd.read_csv(d, sep="\t", index_col=0); c.columns = ["cl"]
        for pb, cl in c["cl"].items(): ps[pb] = f"{st}|U{cl}"
    fa["psc"] = [ps.get(pb, "NA") for pb in fa["pseudobulk"]]
    psc_tot = fa["psc"].value_counts()
    constituents = {}
    for p, g in fa.groupby("final_program"):
        vc = g["psc"].value_counts(); vc = vc[vc.index != "NA"]
        keep = [k for k in vc.index if psc_tot.get(k, 0) and vc[k] / psc_tot[k] >= 0.4]
        constituents[p] = ",".join(keep[:8])

    # --- GO-Elite per program: STANDARD UDON output, written to study_aware/goelite/ ---
    print("GO-Elite per program (standard UDON output -> study_aware/goelite/) ...")
    import anndata as adata_mod
    from goelite_enrichment import run_goelite_on_udon
    mk_long = pd.DataFrame([{"marker": g, "top_cluster": p} for p in progs for g in top50[p]])
    aa = adata_mod.AnnData(X=np.zeros((1, 1), dtype="float32"))
    aa.uns["udon_marker_genes_top_n"] = mk_long
    godir = os.path.join(SA, "goelite")
    run_goelite_on_udon(aa, godir, species="Hs", background_genes=common, logger=lambda m: None)
    # top GO term per program (for the heatmap callout), read from the standard selected results
    pathway = {p: "(no enriched term)" for p in progs}
    selp = os.path.join(godir, "GOElite_UDON_selected.tsv")
    if os.path.exists(selp):
        sel = pd.read_csv(selp, sep="\t")
        if len(sel):
            for pc, grp in sel.groupby("cluster"):
                pathway[str(pc)] = grp.sort_values("fdr").iloc[0]["term_name"]
    pd.DataFrame({"program": progs, "top_gene": [top10[p].split(",")[0] for p in progs],
                  "top_go_term": [pathway.get(p, "(no enriched term)") for p in progs]}
                 ).to_csv(os.path.join(SAT, "program_callouts.tsv"), sep="\t", index=False)

    # --- build SUMMARY (p<0.01) from the canonical program-level enrichment + annotations ---
    enr = pd.read_csv(os.path.join(SAT, "satay_cluster_enrichment.tsv"), sep="\t")
    d = enr[enr["p"] < 0.01].rename(columns={"cluster": "Program", "category": "Category",
            "feature": "Feature", "overlap": "nDonors", "odds_ratio": "OddsRatio", "fdr": "FDR"})
    d = d[["Program", "Category", "Feature", "nDonors", "OddsRatio", "p", "FDR"]].sort_values(["Category", "p"])
    d.insert(1, "UDON_cluster_merged", d["Program"].map(pmap))
    d.insert(2, "UDON_clusters_perstudy", d["Program"].map(constituents))
    d.insert(3, "Top10_markers", d["Program"].map(top10))
    d.insert(4, "TopEnrichedPathway", d["Program"].map(pathway))
    d.to_csv(os.path.join(SAT, "SUMMARY_p01.tsv"), sep="\t", index=False)
    print("columns:", list(d.columns))
    print("\nper-program annotation (every program has markers):")
    for p in progs:
        print(f"  {p:4s} merged={str(pmap.get(p,'?')):4s} | pathway={str(pathway[p])[:40]:40s} | top10={top10[p]}")
        print(f"        perstudy_clusters: {constituents.get(p,'')}")


if __name__ == "__main__":
    main()
