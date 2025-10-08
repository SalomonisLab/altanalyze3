import anndata
import pandas as pd

adata = anndata.read_h5ad("/Users/saljh8/Dropbox/Revio/BoneMarrow/test_env/altanalyze3/tests/BEDs/aggregated_counts.h5ad")

# Safely construct UID by casting all components to string
uids = (
    adata.var["chr"].astype(str)
    + ":" + adata.var["start"].astype(str)
    + "-" + adata.var["end"].astype(str)
    + ":" + adata.var["strand"].astype(str)
)
uids.name = "uid"

# Convert to dense DataFrame (features × samples)
df_dense = pd.DataFrame(
    adata.X.toarray().T,
    index=uids,
    columns=adata.obs_names
)

# Add metadata if desired
df_dense["strand"] = adata.var["strand"].astype(str).values
df_dense["type"] = adata.var["type"].astype(str).values

# Reorder columns
cols = list(adata.obs_names) + ["strand", "type"]
df_dense = df_dense[cols]

# Save as TSV
df_dense.to_csv("aggregated_counts_dense.tsv", sep="\t")
