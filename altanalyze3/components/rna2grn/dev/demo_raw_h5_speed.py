#!/usr/bin/env python3.11
"""Raw 10x .h5 pipeline + pure-inference speed benchmark."""
import sys, time
import numpy as np
import pandas as pd
sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
from altanalyze3.components.rna2grn import Rna2GrnBundle

REPO = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2grn"
H5 = "/Users/saljh8/Dropbox/Collaborations/Grimes/RC2/AML-CITE-Seq/AML-7/filtered_feature_bc_matrix.h5"

t0 = time.time()
b = Rna2GrnBundle.load(f"{REPO}/rna2grn_bundle.pkl.gz")
print(f"bundle load: {time.time()-t0:.2f}s")

# ---- pure inference speed (warm) ----
DATA = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn"
d = np.load(f"{DATA}/matched/dataset.npz", allow_pickle=True)
genes = d["genes"].astype(str)
# build a 100-pseudobulk normalized query by tiling reference rows
Xq = np.tile(d["X"][:50], (2, 1))
dfq = pd.DataFrame(Xq, columns=genes)
_ = b.predict_from_dataframe(dfq.iloc[:2], normalized=True)  # warmup
t0 = time.time()
res = b.predict_from_dataframe(dfq, normalized=True)
dt = time.time() - t0
print(f"pure inference: {len(dfq)} pseudobulks -> GRN in {dt*1000:.0f} ms "
      f"({1000*dt/len(dfq):.2f} ms/pseudobulk)")

# ---- raw 10x h5 end-to-end ----
import scanpy as sc
t0 = time.time()
adata = sc.read_10x_h5(H5)
adata.var_names_make_unique()
read_s = time.time() - t0
print(f"\nraw 10x h5: {adata.shape[0]} cells x {adata.shape[1]} features, read {read_s:.1f}s")

# whole-sample pseudobulk GRN (mechanics + speed; production pairs with cellHarmony labels)
t0 = time.time()
res = b.predict_from_adata(adata, groupby=None, return_neighbors=True)
# predict_from_adata with groupby=None treats each cell as a row -> too many; instead pseudobulk:
print("note: whole-sample pseudobulk path below")

# proper whole-sample pseudobulk: sum all cells
adata.obs["_all"] = "AML7_whole"
t0 = time.time()
res = b.predict_from_adata(adata, groupby="_all", return_neighbors=True)
dt = time.time() - t0
print(f"whole-sample pseudobulk GRN: predicted in {dt:.2f}s "
      f"(matched {res.summary['matched_genes']}/{res.summary['model_gene_count']} genes, "
      f"{res.summary['n_cells_per_group']} cells)")
print("predicted GRN shape:", res.predictions.shape)
print("top reference neighbors:", res.summary["neighbors"].iloc[0, 0][:200])
print("top-10 strongest predicted edges:")
top = res.predictions.iloc[0].sort_values(ascending=False).head(10)
for e, v in top.items():
    print(f"   {e:22s} {v:.3f}")
