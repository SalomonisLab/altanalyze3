#!/usr/bin/env python3
"""Final UDON heatmap for the batch-largely-excluded result (study-aware + SVA).

Uses the final program assignment of every pseudobulk (final_program_assignments.tsv),
recomputes the SVA-corrected (batch-removed) combined fold matrix the same way the
integration did, runs MarkerFinder to get program markers, and renders the canonical
MarkerFinder heatmap (markers x pseudobulks, ordered + color-barred by final program).
Pure Python (numpy/scipy/sklearn). No R.
"""
import os, sys, re, numpy as np, pandas as pd, anndata as ad
UDON_DIR = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, UDON_DIR); os.chdir(UDON_DIR)
import pseudobulk_protocol as P
from study_aware_integrate import sva_remove
from markerFinder import marker_finder_wrapper
from visualizations import plot_markers_df

PB = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk"
SA = os.path.join(PB, "UDON", "study_aware")
S, CT, AN = "Sample", "Hs-BM-titrated-reference-centroid", "Annotation"
TOP_N = 40


def main():
    fa = pd.read_csv(os.path.join(SA, "final_program_assignments.tsv"), sep="\t")
    common = list(pd.read_csv(os.path.join(SA, "per_study_cluster_centroids.tsv"),
                              sep="\t", index_col=0).columns)
    print(f"final programs: {fa['final_program'].nunique()} ; pseudobulks: {len(fa)} ; genes: {len(common)}")

    # reconstruct combined matched folds (common genes) -- same as the integration
    sel = P.read_control_annotation(os.path.join(PB, "UDON", "matched_sex_platform", "control_annotation.tsv"))
    a = ad.read_h5ad(os.path.join(PB, "pseudobulk_scaled_log2_hashed.h5ad"))
    sus = np.array([bool(re.search("TotalSeq|empty", str(s), re.I)) or str(x) == "0"
                    for s, x in zip(a.obs[S], a.obs[AN])]); a = a[~sus].copy()
    folds, fobs = P.build_fold_matrix(a, S, CT, AN, mode="matched", selected_controls=sel, logger=lambda m: None)
    folds = folds.loc[common]

    fa = fa[fa["pseudobulk"].isin(folds.columns)].copy()
    folds = folds[fa["pseudobulk"].tolist()]
    prog = fa["final_program"].values
    study = fa["Dataset"].astype(str).values

    # SVA-correct (final program = bio model, study = batch) -> batch-removed expression
    bio = pd.get_dummies(pd.Series(prog)).values
    Fclean = sva_remove(folds.values, bio, study, log=print)        # genes x pseudobulks
    Fc = pd.DataFrame(Fclean, index=common, columns=folds.columns)

    # MarkerFinder on the batch-removed matrix
    groups = pd.DataFrame({"cluster": prog}, index=folds.columns)
    _, markers, heat = marker_finder_wrapper(input_df=Fc.T, groups=groups, top_n=TOP_N,
                                             rho_threshold=0.2, marker_finder_rho=0.2)
    print(f"markers: {len(markers)} across {markers['top_cluster'].nunique()} programs")
    markers.to_csv(os.path.join(SA, "final_program_markers.txt"), sep="\t", index=False)
    heat.to_csv(os.path.join(SA, "final_program_heatmap.txt"), sep="\t")
    plot_markers_df(heat, markers, groups, os.path.join(SA, "final_program_heatmap.pdf"))
    print("wrote final_program_heatmap.png + .pdf + final_program_markers.txt + final_program_heatmap.txt")
    # program sizes (for the color bar legend)
    print(groups["cluster"].value_counts().to_string())


if __name__ == "__main__":
    main()
