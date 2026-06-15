#!/usr/bin/env python3.11
"""rna2grn — Phase 0 data inspection.

Inspect the GRN connection-score matrix and the RNA pseudobulk h5ads to
understand structure, naming, protocol grouping, and (critically) whether there
is across-sample signal to predict beyond the static mean GRN.
"""
import sys
import json
import numpy as np
import pandas as pd

GRN = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/tf_to_gene_connection_scores_log10-NOT_ordered_clusters_COMBINED_AML_Multiome_RC2.txt"
RNA_COUNTS = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk/pseudobulk_counts_hashed.h5ad"
RNA_SCALED = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk/pseudobulk_scaled_log2_hashed.h5ad"
OUT = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn/reports"


def sep(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


def inspect_grn():
    sep("GRN MATRIX")
    df = pd.read_csv(GRN, sep="\t")
    print("shape:", df.shape)
    print("first two cols:", df.columns[:2].tolist())
    sample_cols = df.columns[2:].tolist()
    print("n sample-cellstate columns:", len(sample_cols))
    tf = df.iloc[:, 0].astype(str)
    gene = df.iloc[:, 1].astype(str)
    print("n unique TF_CON:", tf.nunique())
    print("n unique Gene:", gene.nunique())
    print("n rows (edges):", len(df))
    # are rows unique (TF,Gene) pairs?
    pair = tf + "||" + gene
    print("n unique (TF,Gene) pairs:", pair.nunique(), "(== n rows?)", pair.nunique() == len(df))
    # self edges
    print("n self edges (TF==Gene):", int((tf == gene).sum()))
    # how many target genes per TF
    per_tf = df.groupby(tf).size()
    print("targets per TF: min/median/max =", per_tf.min(), int(per_tf.median()), per_tf.max())
    print("\nfirst 20 TF_CON unique:", sorted(tf.unique())[:20])
    print("first 20 Gene unique:", sorted(gene.unique())[:20])

    # value distribution
    V = df[sample_cols].to_numpy(dtype=float)
    print("\nvalue matrix shape:", V.shape, "dtype float")
    print("global min/median/max:", np.nanmin(V), np.nanmedian(V), np.nanmax(V))
    print("n NaN:", int(np.isnan(V).sum()), " n negative:", int((V < 0).sum()), " n zero:", int((V == 0).sum()))
    print("(name says log10 but values look positive small -> likely abs connection scores)")

    # ---- column name parsing: identify protocol/sample groups ----
    sep("GRN COLUMN NAME PARSING")
    print("ALL", len(sample_cols), "columns:")
    for c in sample_cols:
        print("  ", c)

    # group by leading token to see structure
    import collections
    lead = collections.Counter(c.split("_")[0] for c in sample_cols)
    sep("leading token frequency")
    for k, v in lead.most_common():
        print(f"  {k:20s} {v}")

    # ---- across-sample signal vs static mean (THE key diagnostic) ----
    sep("ACROSS-SAMPLE SIGNAL STRUCTURE (signal vs static mean)")
    Vc = np.nan_to_num(V, nan=0.0)
    edge_mean = Vc.mean(axis=1)
    edge_std = Vc.std(axis=1)
    cv = edge_std / (edge_mean + 1e-12)
    print("per-edge mean across samples: min/median/max =",
          np.min(edge_mean), np.median(edge_mean), np.max(edge_mean))
    print("per-edge std across samples : min/median/max =",
          np.min(edge_std), np.median(edge_std), np.max(edge_std))
    print("per-edge CV (std/mean)      : median =", np.median(cv))
    # how much of total variance is sample-specific vs edge-baseline?
    grand = Vc.mean()
    total_ss = ((Vc - grand) ** 2).sum()
    edge_baseline_ss = (((edge_mean - grand) ** 2)[:, None] * np.ones((1, Vc.shape[1]))).sum()
    resid_ss = ((Vc - edge_mean[:, None]) ** 2).sum()
    print(f"\nVariance decomposition of GRN matrix:")
    print(f"  total SS                : {total_ss:.3e}")
    print(f"  edge-baseline SS (mean) : {edge_baseline_ss:.3e}  ({100*edge_baseline_ss/total_ss:.1f}% of total)")
    print(f"  residual SS (sample-spec): {resid_ss:.3e}  ({100*resid_ss/total_ss:.1f}% of total)")
    print("  -> residual % is the 'predictable beyond mean' headroom; if tiny, task is trivial")

    # correlation of each sample's profile with the global mean profile (inflation check)
    prof_corr = []
    for j in range(Vc.shape[1]):
        prof_corr.append(np.corrcoef(Vc[:, j], edge_mean)[0, 1])
    prof_corr = np.array(prof_corr)
    print(f"\nper-sample profile corr with GLOBAL MEAN profile: median={np.median(prof_corr):.3f} "
          f"min={prof_corr.min():.3f}  (high => profile-corr metric is inflated by static structure)")

    # save edge index + parsed columns
    edge_df = pd.DataFrame({"TF": tf, "Gene": gene, "edge_mean": edge_mean, "edge_std": edge_std})
    edge_df.to_csv(f"{OUT}/grn_edge_stats.csv", index=False)
    pd.Series(sample_cols, name="grn_column").to_csv(f"{OUT}/grn_columns.csv", index=False)
    print(f"\nwrote {OUT}/grn_edge_stats.csv and grn_columns.csv")
    return sample_cols, tf.unique().tolist(), gene.unique().tolist()


def inspect_h5ad(path, label):
    sep(f"RNA pseudobulk h5ad: {label}")
    import anndata as adata_mod
    A = adata_mod.read_h5ad(path)
    print("shape (obs x var):", A.shape)
    print("obs columns:", list(A.obs.columns))
    print("var columns:", list(A.var.columns))
    print("layers:", list(A.layers.keys()))
    print("obsm:", list(A.obsm.keys()))
    print("\nfirst 25 obs_names:")
    for n in A.obs_names[:25].tolist():
        print("   ", n)
    print("\nfirst 15 var_names:", A.var_names[:15].tolist())
    # X stats
    X = A.X
    import scipy.sparse as sp
    if sp.issparse(X):
        xs = X[:50].toarray() if X.shape[0] > 50 else X.toarray()
        print("\nX is sparse", X.dtype, "density~", X.nnz / (X.shape[0] * X.shape[1]))
    else:
        xs = np.asarray(X[:50])
        print("\nX is dense", X.dtype)
    print("X[:50] min/median/max:", float(np.min(xs)), float(np.median(xs)), float(np.max(xs)))
    print("X[:50] max int-ness (looks like counts?):", np.allclose(xs, np.round(xs)))
    # show a few obs rows
    print("\nobs head:")
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(A.obs.head(8).to_string())
    # unique values of small-cardinality obs columns
    for c in A.obs.columns:
        try:
            nu = A.obs[c].nunique()
            if nu <= 60:
                print(f"\nobs['{c}'] ({nu} unique): {sorted(map(str, A.obs[c].unique()))[:60]}")
            else:
                print(f"\nobs['{c}'] has {nu} unique values (sample): {list(map(str, A.obs[c].unique()[:10]))}")
        except Exception as e:
            print(f"obs['{c}'] err {e}")
    return A


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("grn", "all"):
        inspect_grn()
    if which in ("counts", "all"):
        inspect_h5ad(RNA_COUNTS, "counts_hashed")
    if which in ("scaled", "all"):
        inspect_h5ad(RNA_SCALED, "scaled_log2_hashed")
