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
SAT = os.path.join(PB, "UDON", "satay_metadata")
SA = os.path.join(PB, "UDON", "study_aware")
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
    cent = pd.DataFrame({p: Fc.iloc[:, (prog == p)].mean(axis=1).values for p in progs}, index=common)
    top10, top50 = {}, {}
    for p in progs:
        spec = (cent[p] - cent.drop(columns=p).max(axis=1)).sort_values(ascending=False)  # program vs best-other
        top10[p] = ",".join(spec.index[:10]); top50[p] = list(spec.index[:50])
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

    # --- GO-Elite top pathway per program (annotation) ---
    print("GO-Elite top pathway per program ...")
    from altanalyze3.components.goelite import GOEliteRunner, EnrichmentSettings, prepare_species_resources
    parsed = prepare_species_resources("human")
    runner = GOEliteRunner(parsed, settings=EnrichmentSettings(min_term_size=5, max_term_size=2000))
    prepared = runner.prepare_background([str(g) for g in common])
    pathway = {}
    for p in progs:
        try:
            res = runner.run_prepared(top50[p], prepared, apply_prioritization=True)
            sel2 = sorted([r for r in res if getattr(r, "selected", False)], key=lambda r: getattr(r, "fdr", 1.0))
            if sel2:
                node = runner.go_tree.get(sel2[0].term_id) if hasattr(runner, "go_tree") else None
                pathway[p] = getattr(node, "name", "") if node else sel2[0].term_id
            else:
                pathway[p] = "(no enriched term)"
        except Exception:
            pathway[p] = "(no enriched term)"

    # --- merge into SUMMARY tsv ---
    d = pd.read_csv(os.path.join(SAT, "SUMMARY_study_aware_final_p01.tsv"), sep="\t")
    for c in ["TopMarkers", "markers"]:
        if c in d.columns: d = d.drop(columns=c)
    d.insert(1, "UDON_cluster_merged", d["Program"].map(pmap))
    d.insert(2, "UDON_clusters_perstudy", d["Program"].map(constituents))
    d.insert(3, "Top10_markers", d["Program"].map(top10))
    d.insert(4, "TopEnrichedPathway", d["Program"].map(pathway))
    d.to_csv(os.path.join(SAT, "SUMMARY_study_aware_final_p01.tsv"), sep="\t", index=False)
    print("columns:", list(d.columns))
    print("\nper-program annotation (every program has markers):")
    for p in progs:
        print(f"  {p:4s} merged={str(pmap.get(p,'?')):4s} | pathway={str(pathway[p])[:40]:40s} | top10={top10[p]}")
        print(f"        perstudy_clusters: {constituents.get(p,'')}")


if __name__ == "__main__":
    main()
