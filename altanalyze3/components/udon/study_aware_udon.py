#!/usr/bin/env python3
"""
Study-aware UDON with GO-Elite-based integration (no black box).

1. Combined matched-control fold matrix (sex+platform), with the non-negative
   correction applied ONCE on the combined matrix (shared per-gene shift), then
   split by study.
2. Per study: UDON (feature selection -> NMF -> MarkerFinder). SAVE per-study
   outputs so each study is inspectable:
      per_study/<Study>/udon_clusters.txt, marker_genes.txt, marker_heatmap.txt,
      marker_heatmap.png/.pdf (canonical MarkerFinder layout).
3. GO-Elite on the top-50 marker genes of EVERY (study,cluster); cluster-cluster
   similarity = Jaccard of significant GO terms. Conserved cluster = its GO group
   spans >= 2 studies; unique cluster = GO group in a single study.
4. Save centroids + top-50 markers + GO terms + conserved/unique calls for the
   downstream SVA integration (study_aware_integrate.py).

Pure Python (numpy/scipy/sklearn/anndata + altanalyze3 GO-Elite). No R.
"""
import os, sys, numpy as np, pandas as pd, anndata as ad
import scipy.sparse.csgraph as csgraph

UDON_DIR = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, UDON_DIR); os.chdir(UDON_DIR)
import pseudobulk_protocol as P
import fast_feature_selection as FFS
from nmf import run_nmf
from markerFinder import marker_finder_wrapper
from visualizations import plot_markers_df

PB = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk"
OUT = os.path.join(PB, "UDON", "study_aware"); os.makedirs(OUT, exist_ok=True)
PERSTUDY = os.path.join(OUT, "per_study"); os.makedirs(PERSTUDY, exist_ok=True)
S, CT, AN = "Sample", "Hs-BM-titrated-reference-centroid", "Annotation"
MIN_STUDY_PB = 200
TOP_MARKERS = 50
GO_TERM_JACCARD_THR = 0.10     # conserved if GO-term Jaccard >= this links clusters across studies


def build_combined_folds(log):
    SM = pd.read_csv(os.path.join(PB, "evaluation/outputs/sample_metadata.tsv"), sep="\t", index_col=0)
    chem = SM["chemistry"].astype(str).to_dict(); study = SM["Dataset"].astype(str).to_dict()
    a = ad.read_h5ad(os.path.join(PB, "pseudobulk_scaled_log2_hashed.h5ad"))
    import re
    sus = np.array([bool(re.search("TotalSeq|empty", str(s), re.I)) or str(x) == "0"
                    for s, x in zip(a.obs[S], a.obs[AN])]); a = a[~sus].copy()
    cad = ad.read_h5ad(os.path.join(PB, "pseudobulk_counts_hashed.h5ad"))
    sex = P.predict_sample_sex(a, S, counts_adata=cad)
    sel, _ = P.select_matched_controls(a, S, CT, AN, "n_cells", sex, chemistry_of=chem, logger=lambda m: None)
    folds, fobs = P.build_fold_matrix(a, S, CT, AN, mode="matched", selected_controls=sel, logger=log)
    fobs["Dataset"] = fobs[S].map(study)
    log(f"[combined] global non-negative folds {folds.shape}; min={float(folds.to_numpy().min()):.3f}")
    return folds, fobs


def udon_one_study(fst, log):
    fs = FFS.fast_feature_selection(fst, apply_gene_name_filter=False)
    corr = list(fs["correlated_genes"])
    if len(corr) < 50 or fst.shape[1] < 12:
        return None, None, None
    fhvg = fst.loc[corr]; rank = max(2, FFS.fast_rank(fhvg))
    _, clusters = run_nmf(df=fhvg, rank=rank); clusters.index = fst.columns
    keep = clusters["cluster"].value_counts(); keep = keep[keep >= 5].index
    clusters = clusters[clusters["cluster"].isin(keep)]
    if clusters["cluster"].nunique() < 2:
        return None, None, None
    mdf_og, markers, heat = marker_finder_wrapper(input_df=fst[clusters.index].transpose(), groups=clusters,
                                                  top_n=TOP_MARKERS, rho_threshold=0.2, marker_finder_rho=0.2)
    return clusters, markers, heat


def main():
    rep = open(os.path.join(OUT, "study_aware_report.txt"), "w")
    def log(m): print(m); rep.write(str(m) + "\n"); rep.flush()

    folds, fobs = build_combined_folds(log)
    common = list(FFS.fast_feature_selection(folds, apply_gene_name_filter=False)["correlated_genes"])
    log(f"[combined] common centroid gene space: {len(common)} genes")
    studies = fobs["Dataset"].value_counts(); studies = studies[studies >= MIN_STUDY_PB].index.tolist()
    log(f"[studies] {dict(fobs['Dataset'].value_counts().loc[studies])}")

    centroids = {}; markers_by = {}
    for st in studies:
        cols = fobs.index[fobs["Dataset"] == st]; fst = folds[cols]
        clusters, markers, heat = udon_one_study(fst, log)
        if clusters is None:
            log(f"  {st}: skipped"); continue
        sd = os.path.join(PERSTUDY, st.replace(" ", "_")); os.makedirs(sd, exist_ok=True)
        clusters.to_csv(os.path.join(sd, "udon_clusters.txt"), sep="\t")
        markers.to_csv(os.path.join(sd, "marker_genes.txt"), sep="\t", index=False)
        heat.to_csv(os.path.join(sd, "marker_heatmap.txt"), sep="\t")
        try:
            plot_markers_df(heat, markers, clusters, os.path.join(sd, "marker_heatmap.pdf"))
        except Exception as e:
            log(f"    {st}: heatmap skipped: {e}")
        log(f"  {st}: {len(cols)} pseudobulks -> {clusters['cluster'].nunique()} clusters "
            f"(saved per_study/{st.replace(' ','_')}/)")
        for cl, g in clusters.groupby("cluster"):
            key = f"{st}|U{cl}"
            centroids[key] = fst.loc[common, g.index].mean(axis=1)
            markers_by[key] = markers[markers["top_cluster"] == cl]["marker"].tolist()[:TOP_MARKERS]

    pd.DataFrame(centroids).T.to_csv(os.path.join(OUT, "per_study_cluster_centroids.tsv"), sep="\t")
    pd.DataFrame([{"cluster_key": k, "study": k.split("|")[0], "markers": ",".join(v)}
                  for k, v in markers_by.items()]).to_csv(os.path.join(OUT, "per_cluster_markers.tsv"),
                                                          sep="\t", index=False)
    keys = list(centroids); log(f"\n[clusters] {len(keys)} per-study UDON clusters total")

    # ---- GO-Elite per cluster (top-50 markers) ----
    log("[go-elite] running GO-Elite on each cluster's top-50 markers ...")
    sys.path.insert(0, os.path.abspath(os.path.join(UDON_DIR, "..", "..", "..")))
    from altanalyze3.components.goelite import GOEliteRunner, EnrichmentSettings, prepare_species_resources
    parsed = prepare_species_resources("human")
    runner = GOEliteRunner(parsed, settings=EnrichmentSettings(min_term_size=5, max_term_size=2000))
    background = [str(g) for g in pd.Index(common).dropna().unique()]
    prepared = runner.prepare_background(background)
    go_terms = {}
    rows = []
    for k in keys:
        try:
            res = runner.run_prepared(markers_by[k], prepared, apply_prioritization=True)
        except Exception:
            go_terms[k] = set(); continue
        sel_terms = set()
        for r in res:
            if bool(getattr(r, "selected", False)):
                sel_terms.add(r.term_id)
                node = runner.go_tree.get(r.term_id) if hasattr(runner, "go_tree") else None
                rows.append({"cluster_key": k, "study": k.split("|")[0], "term_id": r.term_id,
                             "term_name": getattr(node, "name", "") if node else "",
                             "z_score": getattr(r, "z_score", None), "fdr": getattr(r, "fdr", None)})
        go_terms[k] = sel_terms
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "per_cluster_goelite_selected.tsv"), sep="\t", index=False)
    nwith = sum(1 for k in keys if go_terms[k])
    log(f"[go-elite] {nwith}/{len(keys)} clusters have >=1 prioritized GO term")

    # ---- GO-term Jaccard similarity + conserved/unique ----
    n = len(keys); G = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            a_, b_ = go_terms[keys[i]], go_terms[keys[j]]
            G[i, j] = len(a_ & b_) / max(len(a_ | b_), 1) if (a_ or b_) else 0.0
    pd.DataFrame(G, index=keys, columns=keys).to_csv(os.path.join(OUT, "goterm_jaccard.tsv"), sep="\t")
    A = (G >= GO_TERM_JACCARD_THR).astype(int); np.fill_diagonal(A, 0)
    ncomp, comp = csgraph.connected_components(A, directed=False)
    study_of = np.array([k.split("|")[0] for k in keys])
    comp_studies = {c: set(study_of[comp == c]) for c in range(ncomp)}
    out = pd.DataFrame({"cluster_key": keys, "study": study_of, "go_group": comp,
                        "n_studies_in_group": [len(comp_studies[c]) for c in comp],
                        "status": ["conserved" if len(comp_studies[comp[i]]) >= 2 else "unique"
                                   for i in range(n)],
                        "n_go_terms": [len(go_terms[k]) for k in keys]})
    out.to_csv(os.path.join(OUT, "conserved_vs_unique.tsv"), sep="\t", index=False)
    ncons = int((out["status"] == "conserved").sum()); nuniq = int((out["status"] == "unique").sum())
    ngrp_cons = int(sum(1 for c in range(ncomp) if len(comp_studies[c]) >= 2))
    log(f"\n[integration:GO-Elite] {ncomp} GO groups | conserved clusters={ncons} (in {ngrp_cons} multi-study groups) "
        f"| unique clusters={nuniq}")
    log("[integration] example conserved GO groups (>=3 studies):")
    for c in range(ncomp):
        if len(comp_studies[c]) >= 3:
            mem = [keys[i] for i in range(n) if comp[i] == c]
            shared = set.intersection(*[go_terms[k] for k in mem]) if len(mem) > 1 else set()
            names = [r["term_name"] for r in rows if r["term_id"] in shared][:4]
            log(f"   group {c} ({len(comp_studies[c])} studies, {len(mem)} clusters): shared GO={list(set(names))[:4]}")
    rep.close()
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
