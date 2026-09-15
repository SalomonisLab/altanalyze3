#!/usr/bin/env python3
"""Transcription-factor enrichment for UDON clusters built on GRN TF-target edges.

Gene Ontology asks which biological process a set of GENES belongs to. A GRN UDON run has no
gene set. Its features are regulatory edges named `TF|target`, so the question that matches the
data is which TF drives a cluster's marker edges. This module answers that and replaces GO-Elite
whenever the input modality is `grn`.

Method, and why it is the one the counts justify:

  Each cluster's marker edges are drawn from the full edge universe without replacement, so the
  count of a TF's edges among them follows the hypergeometric distribution. Fisher's exact test,
  one sided in the `greater` direction, IS that hypergeometric tail. It is exact at every count,
  which matters here because most TFs own only a handful of edges and a chi-square would violate
  its expected-count-of-5 requirement on them. Benjamini-Hochberg then runs WITHIN each cluster,
  because each cluster is a separate question, which is how SATAY-UDON corrects per covariate.

  For every (cluster, TF) the 2x2 table is

      a  edges of this TF among this cluster's markers
      b  edges of other TFs among this cluster's markers
      c  edges of this TF elsewhere in the universe
      d  every remaining edge

  Effect size is reported two ways: the odds ratio and the fold enrichment, observed over the
  count expected if the cluster drew edges at random.

Output mirrors `goelite_enrichment.run_goelite_on_udon` so the marker-heatmap callout code reads
either source through the same columns (`cluster`, `term_name`, `fdr`, `selected`):

    TF_enrichment_UDON.tsv           every tested (cluster, TF)
    TF_enrichment_UDON_selected.tsv  the rows passing the FDR cut
    cluster_edge_lists/cluster_<c>_marker_edges.txt

Pure Python (numpy/scipy/pandas/statsmodels). No R.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

EDGE_SEP = "|"          # feature names arrive as TF|target
DISPLAY_SEP = "-"       # and are shown as TF-target


def split_edge(name, sep=EDGE_SEP):
    """`SPI1|SIGLEC9` -> ('SPI1', 'SIGLEC9'). Returns (None, None) when the name is not an edge."""
    s = str(name)
    if sep not in s:
        return None, None
    tf, target = s.split(sep, 1)
    return tf.strip(), target.strip()


def display_feature(name, sep=EDGE_SEP, display_sep=DISPLAY_SEP):
    """Feature name as it should be PRINTED: a GRN edge reads `TF-target`, not `TF|target`.
    A plain gene symbol carries no separator and comes back unchanged."""
    s = str(name)
    return s.replace(sep, display_sep) if sep in s else s


def run_tf_enrichment_on_udon(adata, outdir, marker_key="udon_marker_genes_top_n",
                              background_features=None, fdr_alpha=0.05, min_edges=2,
                              logger=print):
    """Fisher's exact TF enrichment per UDON cluster over `TF|target` features.
    Returns the combined DataFrame, or None when the input carries no edge-shaped features."""
    markers = adata.uns.get(marker_key)
    if markers is None or len(markers) == 0:
        logger("[tf-enrich] SKIPPED — no marker features found")
        return None
    markers = pd.DataFrame(markers)

    if background_features is None:
        background_features = (list(adata.varm["pseudobulk_folds"].index)
                               if "pseudobulk_folds" in adata.varm else list(adata.var_names))
    universe = pd.Index([str(f) for f in background_features]).dropna().unique()
    uni_tf = np.array([split_edge(f)[0] for f in universe], dtype=object)
    is_edge = uni_tf != None                                             # noqa: E711
    n_bad = int((~is_edge).sum())
    if not is_edge.any():
        logger(f"[tf-enrich] SKIPPED — no feature carries the '{EDGE_SEP}' separator; "
               f"this input is not a TF-target matrix")
        return None
    if n_bad:
        logger(f"[tf-enrich] {n_bad} of {len(universe)} features carry no '{EDGE_SEP}' and are "
               f"excluded from the universe")
    universe, uni_tf = universe[is_edge], uni_tf[is_edge]
    N = len(universe)
    tf_total = pd.Series(uni_tf).value_counts()                          # edges per TF, universe
    logger(f"[tf-enrich] universe: {N} edges over {len(tf_total)} TFs "
           f"(median {int(tf_total.median())} edges/TF)")

    os.makedirs(outdir, exist_ok=True)
    edges_dir = os.path.join(outdir, "cluster_edge_lists")
    os.makedirs(edges_dir, exist_ok=True)

    rows = []
    for cluster, grp in markers.groupby("top_cluster"):
        feats = [str(f) for f in grp["marker"].tolist()]
        with open(os.path.join(edges_dir, f"cluster_{cluster}_marker_edges.txt"), "w") as fh:
            fh.write("\n".join(display_feature(f) for f in feats) + "\n")
        tfs = [split_edge(f)[0] for f in feats]
        tfs = [t for t in tfs if t is not None]
        K = len(tfs)                                                     # marker edges in cluster
        if K == 0:
            logger(f"[tf-enrich] cluster {cluster}: no edge-shaped marker; skipped")
            continue
        obs = pd.Series(tfs).value_counts()
        recs = []
        for tf, a in obs.items():
            if a < min_edges:                                            # too few to test
                continue
            n_tf = int(tf_total.get(tf, 0))
            b = K - int(a)
            c = n_tf - int(a)
            d = N - K - c
            if c < 0 or d < 0:
                raise ValueError(f"[tf-enrich] cluster {cluster} TF {tf}: marker edges are not a "
                                 f"subset of the universe (a={a}, n_tf={n_tf}, K={K}, N={N})")
            orr, p = fisher_exact([[int(a), int(b)], [int(c), int(d)]], alternative="greater")
            expected = K * n_tf / float(N)
            recs.append({"cluster": cluster, "term_id": tf, "term_name": tf,
                         "namespace": "TF_regulon",
                         "p_value": float(p), "odds_ratio": float(orr),
                         "fold_enrichment": (float(a) / expected) if expected > 0 else np.nan,
                         "overlap": int(a), "term_size": n_tf, "query_size": int(K),
                         "expected": round(expected, 3),
                         "overlap_genes": ",".join(
                             display_feature(f) for f in feats
                             if split_edge(f)[0] == tf)})
        if not recs:
            logger(f"[tf-enrich] cluster {cluster}: {K} marker edges, no TF with >= {min_edges}")
            continue
        sub = pd.DataFrame(recs)
        sub["fdr"] = multipletests(sub["p_value"].values, alpha=fdr_alpha, method="fdr_bh")[1]
        sub["selected"] = sub["fdr"] < fdr_alpha
        rows.append(sub)
        top = sub.sort_values(["fdr", "p_value"]).iloc[0]
        logger(f"[tf-enrich] cluster {cluster}: {K} marker edges, {len(sub)} TFs tested, "
               f"{int(sub['selected'].sum())} at FDR<{fdr_alpha}; top {top['term_name']} "
               f"({int(top['overlap'])}/{int(top['term_size'])} edges, "
               f"{top['fold_enrichment']:.1f}x, FDR {top['fdr']:.2e})")

    if not rows:
        logger("[tf-enrich] no enrichment results")
        return None
    df = pd.concat(rows, ignore_index=True)
    out_tsv = os.path.join(outdir, "TF_enrichment_UDON.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)
    sel_tsv = os.path.join(outdir, "TF_enrichment_UDON_selected.tsv")
    df[df["selected"]].to_csv(sel_tsv, sep="\t", index=False)
    logger(f"[tf-enrich] wrote {out_tsv} ({len(df)} rows, {int(df['selected'].sum())} selected)")
    logger(f"[tf-enrich] wrote {sel_tsv}")
    return df
