import pandas as pd
import numpy as np
import pickle
from itertools import combinations
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib.pyplot as plt


with open("ElasticNet_model.pkl", "rb") as f:
    bundle = pickle.load(f)

scaler = bundle["scaler"]
model = bundle["enet_fs"]
X_cols = bundle["X_columns"]
Y_cols = bundle["Y_columns"]

# with open("elastic_multitask_try.pkl", "rb") as f:
#     bundle = pickle.load(f)

# model = bundle["model"]
# scaler_x = bundle["scaler_x"]
# scaler_y = bundle.get("scaler_y", None)
# X_cols = bundle["X_columns"]
# Y_cols = bundle["Y_columns"]

X1 = pd.read_csv("data/feature_blankreduiction.csv", index_col=0)
X3 = pd.read_csv("data/newrna_cell_clair_filtered_symbol.csv", index_col=0).T
Y1 = pd.read_csv("data/Bulk_lipids_cleaned_normalized_median_527.csv", index_col=0)
Y3 = pd.read_csv("data/cell_lipids_cleaned_norm_median_286.csv", index_col=0)

# collapse duplicated genes in X3
X3 = X3.groupby(X3.columns, axis=1).mean()

for df in [X1, X3, Y1, Y3]:
    df.index = df.index.astype(str).str.strip()

X_df = pd.concat([X1, X3], axis=0, join="inner")
Y_df = pd.concat([Y1, Y3], axis=0, join="inner")

# keep matched samples only
common_idx = X_df.index.intersection(Y_df.index)
X_df = X_df.loc[common_idx].copy()
Y_df = Y_df.loc[common_idx].copy()

# numeric cleanup
X_df = X_df.apply(pd.to_numeric, errors="coerce")
X_df = X_df.loc[:, X_df.notna().sum() > 0].fillna(0)

# remove duplicate samples
X_df = X_df[~X_df.index.duplicated(keep="first")]
Y_df = Y_df[~Y_df.index.duplicated(keep="first")]

print("X_df shape:", X_df.shape)
print("Y_df shape:", Y_df.shape)

#test inside train
# test_samples = {
#     'D041_EPI', 'D022_MIC', 'D024_EPI', 'D043_MES', 'D024_END',
#     'D043_EPI', 'D071_EPI', 'D044_MIC', 'D044_MES', 'D019_MIC',
#     'D071_END', 'D038_MIC', 'D043_END', 'D071_PMX', 'D044_PMX',
#     'D022_MES', 'D019_END', 'D024_PMX', 'D038_PMX', 'D018_END',
#     'D022_PMX'
# }

# sample_series = Y_df.index
# train_mask = ~sample_series.isin(test_samples)
# test_mask = sample_series.isin(test_samples)

#Lets exclude 4 patient identifiers from all training sets if the index starts with  "D018","D022","D024","D036" for example "D022_END" and "D022_EPI", "D024_END","D036_MES" 


# X_df = X_df[~X_df.index.duplicated(keep="first")]  # remove duplicate samples

for df in [X1, X3, Y1, Y3]:
    df.index = df.index.astype(str).str.strip()

Y_df = pd.concat([Y1, Y3], axis=0, join="inner")
X_df = pd.concat([X1, X3], axis=0, join = "inner")
print(Y_df.shape)  
print(X_df.shape)

common_idx = X_df.index.intersection(Y_df.index)
print("Matched samples:", len(common_idx))

X_df = X_df.loc[common_idx].copy()
Y_df = Y_df.loc[common_idx].copy()

X_df = X_df.apply(pd.to_numeric, errors="coerce")
X_df = X_df.loc[:, X_df.notna().sum() > 0]

X_df = X_df.fillna(0)

common_samples = Y_df.index.intersection(X_df.index)
X_df = X_df.loc[common_samples]
Y_df = Y_df.loc[common_samples]

X_df = X_df[~X_df.index.duplicated(keep="first")]
Y_df = Y_df[~Y_df.index.duplicated(keep="first")]

# # HELD-OUT INTERNAL VALIDATION
# test_samples = {
#     'D041_EPI', 'D022_MIC', 'D024_EPI', 'D043_MES', 'D024_END',
#     'D043_EPI', 'D071_EPI', 'D044_MIC', 'D044_MES', 'D019_MIC',
#     'D071_END', 'D038_MIC', 'D043_END', 'D071_PMX', 'D044_PMX',
#     'D022_MES', 'D019_END', 'D024_PMX', 'D038_PMX', 'D018_END',
#     'D022_PMX'
# }

# sample_series = Y_df.index
# train_mask = ~sample_series.isin(test_samples)
# test_mask  = sample_series.isin(test_samples)

# X_train = X_df.loc[train_mask].copy()
# Y_train = Y_df.loc[train_mask].copy()

# X_holdout = X_df.loc[test_mask].copy()
# Y_holdout = Y_df.loc[test_mask].copy()

X_train = X_df
Y_train = Y_df
X_holdout = X_df
Y_holdout = Y_df


common_samples = X_holdout.index.intersection(Y_holdout.index)
common_lipids = Y_holdout.columns.intersection(Y_df.columns)

X_holdout_aligned = X_holdout.loc[common_samples].copy()
Y_holdout_aligned = Y_holdout.loc[common_samples, common_lipids].copy()

nan_counts = Y_holdout_aligned.isna().sum()
Y_holdout_aligned = Y_holdout_aligned.loc[:, nan_counts <= 16].fillna(0)

# align genes to model input
X_holdout_aligned = X_holdout_aligned.reindex(columns=X_cols, fill_value=0.0)

X_holdout_scaled = scaler.transform(X_holdout_aligned)
Y_holdout_pred = model.predict(X_holdout_scaled)

pred_holdout_df = pd.DataFrame(
    Y_holdout_pred,
    index=X_holdout_aligned.index,
    columns=Y_cols
)

# restrict predicted lipids to those in truth
common_lipids = Y_holdout_aligned.columns.intersection(pred_holdout_df.columns)
Y_holdout_aligned = Y_holdout_aligned.loc[:, common_lipids]
pred_holdout_df = pred_holdout_df.loc[Y_holdout_aligned.index, common_lipids]

meta_holdout = pd.DataFrame(index=Y_holdout_aligned.index)
meta_holdout["Organ"] = meta_holdout.index.str.split("_").str[0]
meta_holdout["CellType"] = meta_holdout.index.str.split("_").str[-1]


#real validationish
counts = pd.read_csv(
    "data/exp.GSE161382_counts_matrix_CPTT-sample-revised-gene-pediatric.txt",
    sep="\t"
)

groups = pd.read_csv(
    "data/groups.GSE161382_counts_matrix_CPTT-sample-revised-prediatric.txt",
    sep="\t"
)

counts = counts.set_index(counts.columns[0])
groups.columns = ["sample", "cell_type"]

groups = groups[groups["sample"].isin(counts.columns)].copy()
counts = counts[groups["sample"]]

groups["donor"] = groups["sample"].str.split("__").str[0]
groups["celltype"] = groups["cell_type"]

pseudobulk = []
for (donor, ct), rows in groups.groupby(["donor", "celltype"]):
    cols = rows["sample"].tolist()
    expr = counts[cols].mean(axis=1)
    name = f"{donor}_{ct.replace(' ', '_').replace('/', '_').replace('-', '_')}"
    pseudobulk.append(expr.rename(name))

pediatric_X_test = pd.DataFrame(pseudobulk)
pediatric_X_test = pediatric_X_test.apply(pd.to_numeric, errors="coerce").fillna(0.0)

# align to training/model genes
pediatric_X_test_aligned = pediatric_X_test.reindex(columns=X_cols, fill_value=0.0)

print("Pediatric X shape:", pediatric_X_test_aligned.shape)
print("Same feature count?", pediatric_X_test_aligned.shape[1] == len(X_cols))
print("Same order?", (pediatric_X_test_aligned.columns == pd.Index(X_cols)).all())

pediatric_X_scaled = scaler.transform(pediatric_X_test_aligned)
pediatric_X_scaled = np.nan_to_num(pediatric_X_scaled, nan=0.0, posinf=0.0, neginf=0.0)

pediatric_Y_pred = model.predict(pediatric_X_scaled)

pediatric_pred_df = pd.DataFrame(
    pediatric_Y_pred,
    index=pediatric_X_test_aligned.index,
    columns=Y_cols
)


meta_pediatric = pd.DataFrame(index=pediatric_pred_df.index)
meta_pediatric["Organ"] = meta_pediatric.index.str.split("_").str[0]
meta_pediatric["CellType_raw"] = meta_pediatric.index.str.split("_", n=1).str[1]

def map_to_coarse(ct):
    ct = str(ct).lower()

    if ("pmn" in ct) or ("polymorph" in ct):
        return "PMX"
    if ("macrophage" in ct) or ("mono" in ct) or ("dc" in ct) or ("dendritic" in ct):
        return "MIC"
    if ("neutro" in ct):
        return "PMX"   # put neutrophils here if that is your intent
    if ("endothelial" in ct) or ("cap" in ct) or ("vein" in ct) or ("arter" in ct):
        return "END"
    if ("fibro" in ct) or ("myofibro" in ct) or ("pericyte" in ct) or ("smooth_muscle" in ct) or ("stromal" in ct):
        return "MES"
    if ("at1" in ct) or ("at2" in ct) or ("club" in ct) or ("ciliated" in ct) or ("epithelial" in ct):
        return "EPI"
    return "OTHER"

meta_pediatric["CellType"] = meta_pediatric["CellType_raw"].map(map_to_coarse)

print(meta_pediatric["CellType"].value_counts(dropna=False))

# keep only desired coarse groups
meta_pediatric5 = meta_pediatric[
    meta_pediatric["CellType"].isin(["END", "EPI", "MES", "MIC", "PMX"])
].copy()

pediatric_pred5 = pediatric_pred_df.loc[meta_pediatric5.index].copy()

print(meta_pediatric5["CellType"].value_counts())

celltype_pairs = [
    ('END', 'MES'),
    ('END', 'MIC'),
    ('END', 'PMX'),
    ('EPI', 'MES'),
    ('EPI', 'MIC'),
    ('EPI', 'PMX'),
    ('MES', 'MIC'),
    ('MES', 'PMX'),
    ('MIC', 'PMX')
]


from statsmodels.stats.multitest import multipletests
def compute_pairwise_lipids(df, metadata, group_col, group1, group2, pseudocount=1e-6):
    common_idx = df.index.intersection(metadata.index)
    df = df.loc[common_idx].copy()
    meta = metadata.loc[common_idx].copy()

    mask = meta[group_col].isin([group1, group2])
    df_sub = df.loc[mask].apply(pd.to_numeric, errors="coerce")
    meta_sub = meta.loc[mask]

    g1_idx = meta_sub[group_col] == group1
    g2_idx = meta_sub[group_col] == group2

    n1, n2 = int(g1_idx.sum()), int(g2_idx.sum())
    print(f"{group1} vs {group2} | n={n1} vs n={n2}")

    required_cols = ["Lipid","Group1","Group2","log2FC","t_stat","pval","FDR","EffectSize","Direction"]

    if n1 < 2 or n2 < 2:
        return pd.DataFrame(columns=required_cols)

    results = []
    for lipid in df_sub.columns:
        vals1 = df_sub.loc[g1_idx, lipid].astype(float).dropna().values
        vals2 = df_sub.loc[g2_idx, lipid].astype(float).dropna().values

        if len(vals1) < 2 or len(vals2) < 2:
            continue

        mean1, mean2 = np.mean(vals1), np.mean(vals2)
        logfc = np.log2((mean1 + pseudocount) / (mean2 + pseudocount))
        stat, pval = ttest_ind(vals1, vals2, equal_var=False)

        results.append({
            "Lipid": lipid,
            "Group1": group1,
            "Group2": group2,
            "log2FC": logfc,
            "t_stat": stat,
            "pval": pval,
            "EffectSize": abs(mean1 - mean2),
            "Direction": f"{group1}_up" if logfc > 0 else f"{group2}_up"
        })

    out = pd.DataFrame(results)
    if out.empty:
        return pd.DataFrame(columns=required_cols)

    out["FDR"] = multipletests(out["pval"].values, method="fdr_bh")[1]
    out = out.sort_values(["FDR","EffectSize"], ascending=[True, False]).reset_index(drop=True)
    return out[required_cols]


def build_pred_top_table(pred_df, top_n=202):
    # pred_df already has FDR and EffectSize
    if pred_df is None or pred_df.empty:
        return pd.DataFrame(columns=["Lipid","logFC_pred","FDR_pred","Effect_pred","Direction"])

    pred_top = pred_df.sort_values(["FDR","EffectSize"], ascending=[True, False]).head(top_n)

    return pred_top[["Lipid","log2FC","FDR","EffectSize","Direction"]].rename(
        columns={"log2FC":"logFC_pred","FDR":"FDR_pred","EffectSize":"Effect_pred"}
    ).reset_index(drop=True)

from itertools import combinations

celltypes_present = sorted(meta5["CellType"].unique())
celltype_pairs = list(combinations(celltypes_present, 2))

def build_comparison_table(true_stats, pred_stats, top_n=202):
    """
    true_stats: output of compute_pairwise_lipids_true(...) with columns Lipid, log2FC, FDR/pval, EffectSize
    pred_stats: output of compute_pairwise_lipids_pred(...) with columns Lipid, log2FC, FDR, EffectSize
    """
    if true_stats.empty and pred_stats.empty:
        return pd.DataFrame(columns=[
            "Lipid","logFC_true","FDR_true","Effect_true",
            "logFC_pred","FDR_pred","Effect_pred","Category"
        ])

    # take top_n from each (ranked already)
    true_top = true_stats.head(top_n).copy() if not true_stats.empty else pd.DataFrame(columns=["Lipid"])
    pred_top = pred_stats.head(top_n).copy() if not pred_stats.empty else pd.DataFrame(columns=["Lipid"])

    true_top = true_top.rename(columns={"log2FC":"logFC_true","FDR":"FDR_true","EffectSize":"Effect_true"})
    pred_top = pred_top.rename(columns={"log2FC":"logFC_pred","FDR":"FDR_pred","EffectSize":"Effect_pred"})

    merged = pd.merge(
        true_top[["Lipid","logFC_true","FDR_true","Effect_true"]],
        pred_top[["Lipid","logFC_pred","FDR_pred","Effect_pred"]],
        on="Lipid",
        how="outer"
    )

    def cat(row):
        t = pd.notna(row["logFC_true"])
        p = pd.notna(row["logFC_pred"])
        if t and p: return "Both"
        if t: return "True_only"
        return "Pred_only"

    merged["Category"] = merged.apply(cat, axis=1)

    # nice ordering: Both first, then True_only, then Pred_only
    order = pd.Categorical(merged["Category"], ["Both","True_only","Pred_only"], ordered=True)
    merged = merged.assign(Category=order).sort_values(["Category"])

    return merged.reset_index(drop=True)
comparison_tables = {}

for g1, g2 in celltype_pairs:
    key = f"{g1}_vs_{g2}"

    pred_stats = compute_pairwise_lipids_predonly(pred5,meta5, "CellType", g1, g2)
    true_stats = compute_pairwise_lipids_predonly(pred5, meta5, "CellType", g1, g2)  # <-- replace with your TRUE function

    if pred_stats.empty or true_stats.empty:
        print(f"Skipping {key}: not enough samples or empty result")
        continue

    comparison_tables[key] = build_comparison_table(true_stats, pred_stats, top_n=202)


# summary = []
# for k, tbl in comparison_tables.items():
#     if tbl.empty:
#         summary.append([k, 0, 0])
#     else:
#         n_sig = (tbl["FDR_pred"] < 0.05).sum()
#         summary.append([k, len(tbl), int(n_sig)])

# summary_df = pd.DataFrame(summary, columns=["Comparison","TopN_returned","N_FDR<0.05_inTopN"])
# summary_df.sort_values("N_FDR<0.05_inTopN", ascending=False)


def build_comparison_table(true_df, pred_df, top_n=25):

    # Guard against empties
    if true_df is None or true_df.empty:
        true_top = pd.DataFrame(columns=["Lipid","log2FC","pval"])
    else:
        if "EffectSize" not in true_df.columns:
            true_df = true_df.assign(EffectSize=true_df["log2FC"].abs())
        true_top = true_df.sort_values(["pval","EffectSize"], ascending=[True,False]).head(top_n)

    if pred_df is None or pred_df.empty:
        pred_top = pd.DataFrame(columns=["Lipid","log2FC","pval"])
    else:
        if "EffectSize" not in pred_df.columns:
            pred_df = pred_df.assign(EffectSize=pred_df["log2FC"].abs())
        pred_top = pred_df.sort_values(["pval","EffectSize"], ascending=[True,False]).head(top_n)

    true_set = set(true_top["Lipid"].dropna())
    pred_set = set(pred_top["Lipid"].dropna())
    all_lipids = sorted(true_set | pred_set)

    combined = pd.DataFrame({"Lipid": all_lipids}).set_index("Lipid")

    combined = combined.join(
        true_top.set_index("Lipid")[["log2FC","pval"]]
        .rename(columns={"log2FC":"logFC_true","pval":"p_true"}),
        how="left"
    )

    combined = combined.join(
        pred_top.set_index("Lipid")[["log2FC","pval"]]
        .rename(columns={"log2FC":"logFC_pred","pval":"p_pred"}),
        how="left"
    )

    combined["Category"] = np.where(
        combined["logFC_true"].notna() & combined["logFC_pred"].notna(), "Both",
        np.where(combined["logFC_true"].notna(), "True only",
                 np.where(combined["logFC_pred"].notna(), "Pred only", "Neither"))
    )

    return combined.reset_index()


import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# If pred_test_df has same samples but different order/format, force align to meta_test
# pred_test_df = pred_test_df.copy()
# pred_test_df.index = pred_test_df.index.astype(str)
# meta_test.index = meta_test.index.astype(str)

# pred_test_df = pred_test_df.loc[pred_test_df.index.intersection(meta_test.index)]
# meta_test_aligned = meta_test.loc[pred_test_df.index]
def build_pred_top_table(pred_df):
    if pred_df.empty:
        return pd.DataFrame(columns=["Lipid","logFC_pred","FDR_pred","Effect_pred","Direction"])
    top = pred_df.head(top_n).copy()
    return top.rename(columns={
        "log2FC":"logFC_pred",
        "FDR":"FDR_pred",
        "EffectSize":"Effect_pred"
    })[["Lipid","logFC_pred","FDR_pred","Effect_pred","Direction"]]


from itertools import combinations

celltypes_present = sorted(meta5["CellType"].unique())
celltype_pairs = list(combinations(celltypes_present, 2))
pred_valid_list = []
def build_valid_df(df, metadata, group_col="CellType", top_n=202, pairs=None):
    if pairs is None:
        celltypes_present = sorted(metadata[group_col].dropna().unique())
        pairs = list(combinations(celltypes_present, 2))

    valid_list = []

    for g1, g2 in pairs:
        stats = compute_pairwise_lipids(df, metadata, group_col, g1, g2)

        if stats is None or stats.empty:
            print(f"Skipping {g1}_vs_{g2}: empty")
            continue

        top = stats.head(top_n).copy()
        top["Comparison"] = f"{g1}_vs_{g2}"
        valid_list.append(top)

    if len(valid_list) == 0:
        print("No valid pairwise comparisons found.")
        return pd.DataFrame(columns=[
            "Lipid","Group1","Group2","log2FC","t_stat","pval","FDR",
            "EffectSize","Direction","Comparison"
        ])

    return pd.concat(valid_list, ignore_index=True)
# Combine all comparison

celltype_pairs = [
    ("END", "EPI"),
    ("END", "MES"),
    ("END", "MIC"),
    ("END", "PMX"),
    ("EPI", "MES"),
    ("EPI", "MIC"),
    ("EPI", "PMX"),
    ("MES", "MIC"),
    ("MES", "PMX"),
    ("MIC", "PMX"),
]

# internal held-out truth
true_valid_df = build_valid_df(
    df=Y_holdout_aligned,
    metadata=meta_holdout,
    group_col="CellType",
    top_n=202,
    pairs=celltype_pairs
)

# internal held-out prediction
pred_valid_df = build_valid_df(
    df=pred_holdout_df,
    metadata=meta_holdout,
    group_col="CellType",
    top_n=202,
    pairs=celltype_pairs
)

# pediatric external prediction
pediatric_pred_valid_df = build_valid_df(
    df=pediatric_pred5,
    metadata=meta_pediatric5,
    group_col="CellType",
    top_n=202,
    pairs=celltype_pairs
)

print("true_valid_df:", true_valid_df.shape)
print("pred_valid_df:", pred_valid_df.shape)
print("pediatric_pred_valid_df:", pediatric_pred_valid_df.shape)


#OPTIONAL: HEATMAPS
def make_logfc_heatmap(valid_df, title, outname=None):
    if valid_df.empty:
        print(f"No data for {title}")
        return

    heat_df = valid_df.pivot_table(
        index="Lipid",
        columns="Comparison",
        values="log2FC",
        aggfunc="mean"
    )

    plt.figure(figsize=(10, 8))
    sns.heatmap(heat_df, cmap="RdBu_r", center=0)
    plt.title(title)
    plt.tight_layout()

    if outname is not None:
        plt.savefig(outname, dpi=300, bbox_inches="tight")

    plt.show()


make_logfc_heatmap(
    true_valid_df,
    title="True lipid log2FC across celltype comparisons",
    outname="true_lipid_heatmap.png"
)

make_logfc_heatmap(
    pred_valid_df,
    title="Predicted lipid log2FC across celltype comparisons",
    outname="pred_lipid_heatmap.png"
)

#verification OF the working or not



true_summary = (
    df[df["Source"] == "True"]
    .groupby("Lipid")
    .agg(
        mean_pval_true=("pval", "mean")
    )
)

pred_summary = (
    df[df["Source"] == "Pred"]
    .groupby("Lipid")
    .agg(
        mean_pval_pred=("pval", "mean")
    )
)



from upsetplot import UpSet, from_contents
import matplotlib.pyplot as plt

alpha = 0.05

# significant rows
true_sig = true_valid_df[true_valid_df["pval"] < alpha].copy()
pred_sig = pediatric_pred_valid_df[pediatric_pred_valid_df["pval"] < alpha].copy()

# clean lipid names (important for matching)
true_sig["Lipid"] = true_sig["Lipid"].astype(str).str.strip()
pred_sig["Lipid"] = pred_sig["Lipid"].astype(str).str.strip()

# find all comparisons
comparisons = sorted(set(true_sig["Comparison"]) | set(pred_sig["Comparison"]))

print("Comparisons found:", comparisons)

for comp in comparisons:

    true_set = set(true_sig.loc[true_sig["Comparison"] == comp, "Lipid"])
    pred_set = set(pred_sig.loc[pred_sig["Comparison"] == comp, "Lipid"])

    print(f"\n{comp}")
    print("True:", len(true_set))
    print("Pred:", len(pred_set))
    print("Overlap:", len(true_set & pred_set))

    sets = {
        "True_DE": true_set,
        "Pred_DE": pred_set
    }

    upset_data = from_contents(sets)

    plt.figure(figsize=(6,4))

    UpSet(
        upset_data,
        subset_size="count",
        show_counts=True
    ).plot()

    plt.suptitle(f"{comp} differential lipids (p < {alpha})")

    plt.tight_layout()
    plt.savefig(f"clea_{comp}_diff_lipids.png", dpi=300)

    plt.show()










from upsetplot import UpSet, from_contents
import matplotlib.pyplot as plt

comparisons = [
    "END_vs_EPI",
    "END_vs_MES",
    "END_vs_MIC",
    "EPI_vs_MES",
    "EPI_vs_MIC",
    "MES_vs_MIC"
]

agreement_sets = {}

for comp in comparisons:
    true_lipids = set(
        combined_df[
            (combined_df["Source"] == "True") &
            (combined_df["Comparison"] == comp)
        ]["Lipid"].astype(str).str.strip()
    )

    pred_lipids = set(
        combined_df[
            (combined_df["Source"] == "Pred") &
            (combined_df["Comparison"] == comp)
        ]["Lipid"].astype(str).str.strip()
    )

    # only matched-pair overlap
    shared_lipids = true_lipids & pred_lipids
    agreement_sets[comp] = shared_lipids

# drop empty ones like END_vs_EPI if needed
agreement_sets = {k: v for k, v in agreement_sets.items() if len(v) > 0}

print({k: len(v) for k, v in agreement_sets.items()})

upset_data = from_contents(agreement_sets)

plt.figure(figsize=(14, 7))
UpSet(
    upset_data,
    subset_size="count",
    show_counts=True,
    sort_categories_by="input",
    sort_by="cardinality",
    min_degree=1,
    max_subset_rank=25
).plot()

plt.suptitle("UpSet: Shared lipids between True and Pred within matching cell-type pairs")
plt.tight_layout()
plt.savefig("matched_true_pred_pairs.png", dpi=300)
plt.show()





from upsetplot import UpSet, from_contents
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

alpha = 0.05

comparisons = [
    "END_vs_EPI",
    "END_vs_MES",
    "END_vs_MIC",
    "EPI_vs_MES",
    "EPI_vs_MIC",
    "MES_vs_MIC"
]

# start from original dfs
true_df = true_valid_df.copy()
pred_df = pred_valid_df.copy()

# clean lipid names
true_df["Lipid"] = true_df["Lipid"].astype(str).str.strip()
pred_df["Lipid"] = pred_df["Lipid"].astype(str).str.strip()

# add source labels
true_df["Source"] = "True"
pred_df["Source"] = "Pred"

# combine
combined_df = pd.concat([true_df, pred_df], ignore_index=True)

# build matched-pair significant sets, keeping True/Pred labels
lipid_sets = {}

for comp in comparisons:
    true_lipids = set(
        combined_df[
            (combined_df["Source"] == "True") &
            (combined_df["Comparison"] == comp) &
            (combined_df["pval"] < alpha)
        ]["Lipid"]
    )

    pred_lipids = set(
        combined_df[
            (combined_df["Source"] == "Pred") &
            (combined_df["Comparison"] == comp) &
            (combined_df["pval"] < alpha)
        ]["Lipid"]
    )

    shared = true_lipids & pred_lipids

    # keep both labels, but only shared lipids for the matching 

    if len(shared) > 0:
        lipid_sets[f"True_{comp}"] = shared
        lipid_sets[f"Pred_{comp}"] = shared

print({k: len(v) for k, v in lipid_sets.items()})

upset_data = from_contents(lipid_sets)

plt.figure(figsize=(16, 8))
UpSet(
    upset_data,
    subset_size="count",
    show_counts=True,
    sort_categories_by="input",
    sort_by="cardinality",
    min_degree=2,
    max_subset_rank=25
).plot()

plt.suptitle(
    f"UpSet: Significant shared lipids between True and Pred for matching comparisons (p < {alpha})"
)

plt.tight_layout()
plt.savefig("upset_true_pred_significant_matched_pairs.png", dpi=300)
plt.show()




import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators

alpha = 0.05

# 1. keep only significant lipids
sig_df = combined_df[combined_df["pval"] < alpha].copy()

# 2. keep only True and Pred rows
sig_df = sig_df[sig_df["Source"].isin(["True", "Pred"])].copy()

# 3. make paired labels
sig_df["Comparison_Source"] = sig_df["Source"] + "_" + sig_df["Comparison"]

# 4. only keep comparisons that have BOTH True and Pred
valid_comparisons = []
for comp, sub in sig_df.groupby("Comparison"):
    sources = set(sub["Source"])
    if {"True", "Pred"}.issubset(sources):
        valid_comparisons.append(comp)

sig_df = sig_df[sig_df["Comparison"].isin(valid_comparisons)].copy()

# 5. build membership matrix
membership = (
    sig_df.assign(value=True)
    .pivot_table(
        index="Lipid",
        columns="Comparison_Source",
        values="value",
        aggfunc="any",
        fill_value=False
    )
    .astype(bool)
)

# 6. explicitly order columns as True/Pred pairs
ordered_cols = []
for comp in sorted(valid_comparisons):
    true_col = f"True_{comp}"
    pred_col = f"Pred_{comp}"
    if true_col in membership.columns and pred_col in membership.columns:
        ordered_cols.extend([true_col, pred_col])

membership = membership[ordered_cols]

# optional: remove lipids that are in none of the columns
membership = membership.loc[membership.any(axis=1)]

# 7. convert for upsetplot
upset_data = from_indicators(membership.columns, membership)

# 8. plot one large UpSet
up = UpSet(
    upset_data,
    subset_size="count",
    show_counts=True,
    sort_by="cardinality",
    sort_categories_by=None
)

up.plot()
plt.suptitle("Significant Lipid Overlap Across Matched True vs Pred Cell-Type Comparisons")
plt.show()
plt.savefig("giantcomparison_corresponding.png")

summary = []

for comp in comparisons:
    sub = sig_df[sig_df["Comparison"] == comp].copy()
    
    true_set = set(sub.loc[sub["Source"] == "True", "Lipid"])
    pred_set = set(sub.loc[sub["Source"] == "Pred", "Lipid"])
    
    both = true_set & pred_set
    true_only = true_set - pred_set
    pred_only = pred_set - true_set
    
    summary.append({
        "Comparison": comp,
        "True_only": len(true_only),
        "Pred_only": len(pred_only),
        "Both": len(both),
        "True_total": len(true_set),
        "Pred_total": len(pred_set),
        "Jaccard": len(both) / len(true_set | pred_set) if len(true_set | pred_set) > 0 else 0
    })

summary_df = pd.DataFrame(summary).sort_values("Jaccard", ascending=False)
print(summary_df)


import math
import pandas as pd
import matplotlib.pyplot as plt

alpha = 0.05

sig_df = combined_df.loc[
    (combined_df["pval"] < alpha) &
    (combined_df["Source"].isin(["True", "Pred"]))
].copy()

valid_comparisons = []
for comp, sub in sig_df.groupby("Comparison"):
    if {"True", "Pred"}.issubset(set(sub["Source"])):
        valid_comparisons.append(comp)

valid_comparisons = sorted(valid_comparisons)

summary = []

for comp in valid_comparisons:
    sub = sig_df[sig_df["Comparison"] == comp]
    true_set = set(sub.loc[sub["Source"] == "True", "Lipid"])
    pred_set = set(sub.loc[sub["Source"] == "Pred", "Lipid"])

    true_only = len(true_set - pred_set)
    pred_only = len(pred_set - true_set)
    both = len(true_set & pred_set)

    summary.append({
        "Comparison": comp,
        "True only": true_only,
        "Both": both,
        "Pred only": pred_only
    })

summary_df = pd.DataFrame(summary)

n = len(summary_df)
ncols = 2
nrows = math.ceil(n / ncols)

fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4 * nrows))
axes = axes.flatten()

for ax, (_, row) in zip(axes, summary_df.iterrows()):
    vals = [row["True only"], row["Both"], row["Pred only"]]
    labels = ["True only", "Both", "Pred only"]

    ax.bar(labels, vals)
    ax.set_title(row["Comparison"])
    ax.set_ylabel("Number of significant lipids")
    ax.tick_params(axis="x", rotation=20)

for ax in axes[len(summary_df):]:
    ax.axis("off")

fig.suptitle("Significant Lipid Overlap by Cell-Type Pair", fontsize=14)
plt.tight_layout()
plt.savefig("Significant Lipid Overlap by Cell-Type Pair")
