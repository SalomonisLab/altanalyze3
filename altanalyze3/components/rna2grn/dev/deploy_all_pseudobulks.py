#!/usr/bin/env python3.11
"""Deploy the final rna2grn model on every pseudobulk in the RNA atlas and write
an imputed-GRN AnnData.

Output h5ad:
  X    : imputed TF->target edge activity scores (n_pseudobulks x 7,486), float32
  var  : one row per TF-target edge; var_names = 'TF|Gene', columns TF and Gene
  obs  : all original obs columns from the RNA pseudobulk h5ad (carried verbatim)
"""
import sys
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
from altanalyze3.components.rna2grn import Rna2GrnBundle

BUNDLE = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2grn/rna2grn_bundle.pkl.gz"
RNA = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk/pseudobulk_counts_hashed.h5ad"
OUT = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn/imputed_grn_all_pseudobulks.h5ad"
CHUNK = 2000


def main():
    b = Rna2GrnBundle.load(BUNDLE)
    edge_ids = list(b.output_edges)
    var = pd.DataFrame({"TF": list(b.metadata["edge_tf"]), "Gene": list(b.metadata["edge_gene"])},
                       index=pd.Index(edge_ids, name="edge"))

    A = ad.read_h5ad(RNA, backed="r")
    obs = A.obs.copy()
    genes = A.var_names.tolist()
    n = A.n_obs
    print(f"atlas: {n} pseudobulks x {A.n_vars} genes; imputing {len(edge_ids)} edges", flush=True)

    scores = np.zeros((n, len(edge_ids)), dtype=np.float32)
    t0 = time.time()
    matched_genes = None
    for start in range(0, n, CHUNK):
        stop = min(start + CHUNK, n)
        sub = A[start:stop].to_memory().X
        sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
        aligned, matched, totals = b._align_matrix(sub, genes)
        matched_genes = matched
        scores[start:stop] = b._predict_aligned_counts(aligned, totals)
        print(f"  {stop}/{n}  ({time.time()-t0:.0f}s)", flush=True)

    out = ad.AnnData(X=scores, obs=obs, var=var)
    out.uns["rna2grn"] = {
        "model": "per_edge_regulon_regression",
        "normalization": "cp10k_log1p",
        "matched_feature_genes": int(matched_genes),
        "n_feature_genes": len(b.input_genes),
        "source_rna_h5ad": RNA,
        "bundle": BUNDLE,
        "value": "imputed TF->target edge connection/activity score (>=0)",
    }
    out.write_h5ad(OUT)
    import os
    print(f"\nwrote {OUT}  ({os.path.getsize(OUT)/1e6:.0f} MB)")
    print(f"shape {out.shape}; obs columns carried: {list(obs.columns)}")
    print(f"var head:\n{var.head(3).to_string()}")
    print(f"X range: {float(scores.min()):.4f} .. {float(scores.max()):.4f}; "
          f"matched genes {matched_genes}/{len(b.input_genes)}")


if __name__ == "__main__":
    main()
