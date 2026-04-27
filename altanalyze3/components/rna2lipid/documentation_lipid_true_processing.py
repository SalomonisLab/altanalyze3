
import numpy as np
import pandas as pd
import numpy as np
import pandas as pd

# =========================================================
# Shared helpers
# =========================================================
def make_unique(index):
    counts = {}
    new_index = []
    for name in index:
        if name not in counts:
            counts[name] = 1
            new_index.append(name)
        else:
            new_index.append(f"{name}.{counts[name]}")
            counts[name] += 1
    return new_index


def clean_lipid_table(df):
    df = df.copy()

    # Clean lipid names
    df.index = df.index.astype(str).str.strip()
    df = df[df.index != ""]

    # Make duplicate lipid names unique
    df.index = make_unique(df.index.tolist())

    return df


def summarize_missingness(df, name="dataset"):
    neg_count = (df < 0).sum().sum()
    na_count_before = df.isna().sum().sum()
    total_values = df.size
    row_missing = df.isna().sum(axis=1)

    print(f"\n--- Missingness summary: {name} ---")
    print("Total negative values:", neg_count)
    print("Total missing (before):", na_count_before)
    print("Total values:", total_values)
    print("Rows with >50% missing:", (row_missing > (df.shape[1] / 2)).sum())
    print("Rows fully missing:", (row_missing == df.shape[1]).sum())


def impute_row_median(df, name="dataset"):
    df = df.copy()
    impute_counter = {"filled": 0, "all_missing_rows": 0}

    def fill_row_median(row):
        row_new = row.copy()

        # Treat negative values as missing
        neg_mask = row_new < 0
        row_new.loc[neg_mask] = np.nan

        n_missing = row_new.isna().sum()
        if n_missing > 0:
            impute_counter["filled"] += n_missing

        # If whole row is missing, fill with 0
        if row_new.isna().all():
            impute_counter["all_missing_rows"] += 1
            return row_new.fillna(0)

        # Otherwise fill with row median
        return row_new.fillna(row_new.median())

    df = df.apply(fill_row_median, axis=1)

    print(f"\n--- Imputation summary: {name} ---")
    print("Total values imputed:", impute_counter["filled"])
    print("Rows fully missing and set to 0:", impute_counter["all_missing_rows"])

    return df


# =========================================================
# Transform helper
# =========================================================
def log2_1p_10x(df):
    """
    Apply log2(1 + 10*x)
    """
    return np.log2(1 + 10 * df)


# =========================================================
# Y1 pipeline
# Assumption: Y1 is already log2-scaled
# Goal: reverse log2, then apply log2(1 + 10*x)
# =========================================================
Y1 = pd.read_csv("/users/ramkd9/Lipid_Predict/Clair/bulk_unduplicated_487.csv", index_col=0).T
Y1 = clean_lipid_table(Y1)
summarize_missingness(Y1, name="Y1 raw (assumed log2 scale)")

Y1_imputed = impute_row_median(Y1, name="Y1 raw (assumed log2 scale)")

# Reverse assumed log2 scale back to linear
Y1_linear = np.power(2, Y1_imputed)

# Apply same final transform
Y1_norm_old = log2_1p_10x(Y1_linear)

print("\nY1 summary:")
print("Raw mean:", Y1.mean().mean())
print("Linear mean after 2**Y1:", Y1_linear.mean().mean())
print("Final normalized mean after log2(1 + 10*x):", Y1_norm.mean().mean())
# Shared samples (donors)
common_samples = Y1.index.intersection(Y3.index)

# Shared lipids
common_lipids = Y1.columns.intersection(Y3.columns)

print("Shared samples:", len(common_samples))
print("Shared lipids:", len(common_lipids))

# =========================================================
# Y3 pipeline
# Assumption: Y3 is already raw / linear
# Goal: directly apply log2(1 + 10*x)
# =========================================================
Y3 = pd.read_csv("/users/ramkd9/Lipid_Predict/Clair/justcombined_ppstNeg.csv", index_col=0).T
Y3 = clean_lipid_table(Y3)
summarize_missingness(Y3, name="Y3 raw (linear scale)")

Y3_imputed = impute_row_median(Y3, name="Y3 raw (linear scale)")

# Apply transform directly
Y3_normold = log2_1p_10x(Y3_imputed)

print("\nY3 summary:")
print("Raw mean:", Y3.mean().mean())
print("Final normalized mean after log2(1 + 10*x):", Y3_normold.mean().mean())

print("\nShapes:")
print("Y1 columns:", Y1_norm_old.shape[1])
print("Y3 columns:", Y3_normold.shape[1])


Y1_common = Y1.loc[Y1.index, common_lipids].copy()
Y3_common = Y3.loc[Y3.index, common_lipids].copy()
# bulk = reference, do not transform
Y1_norm = Y1_common.copy()

# combined = transform candidate
Y3_norm = np.log2(1 + 10 * Y3_common)


# print("Common columns:", len(common_cols))
# common = Y1_imputed.columns.intersection(Y3_norm.columns)
# common_cols = Y1_imputed.columns.intersection(Y3_norm.columns)

# # -----------------------------------
# # Subset both datasets
# # -----------------------------------
# Y1_norm = Y1_imputed[common_cols].copy()
# Y3_norm = Y3_norm[common_cols].copy()

# # Extract only D018 sorted cells
# Y3_D018 = Y3_common[Y3_common.index.str.startswith("D019")].copy()

# # Extract cell type labels
# Y3_D018["CellType"] = Y3_D018.index.str.split("_").str[-1]

# # Use these directly (no averaging)
# Y3_D018_mat = Y3_D018.drop(columns=["CellType"])

# # Align lipids

# # -----------------------------------
# # Find common lipid columns
# # -----------------------------------

# y = Y1_common.loc["D019", common].values
# X = Y3_D018_mat[common].values.T  # lipids × celltypes

# from scipy.optimize import nnls

# w, _ = nnls(X, y)
# w = w / w.sum()

# pred = X @ w
# corr = safe_corr(y, pred)

# import matplotlib.pyplot as plt
# import numpy as np

# sample = "D019"  # adjust if needed

# y = Y1_common.loc[sample].values
# w = weights_df.loc[sample].values
# pred = X @ w

# # Clean values
# mask = np.isfinite(y) & np.isfinite(pred)
# y = y[mask]
# pred = pred[mask]

# # Safe correlation
# def safe_corr(a, b):
#     if len(a) < 2 or np.std(a) == 0 or np.std(b) == 0:
#         return np.nan
#     return np.corrcoef(a, b)[0, 1]

# corr = safe_corr(y, pred)

# # Plot
# plt.figure(figsize=(5, 5))
# plt.scatter(y, pred, alpha=0.7)

# lims = [min(y.min(), pred.min()), max(y.max(), pred.max())]
# plt.plot(lims, lims, "--")

# plt.xlabel("Observed bulk (Y1)")
# plt.ylabel("Reconstructed (from Y3)")
# plt.title(f"{sample}\nPearson r = {corr:.3f}")

# plt.tight_layout()
# plt.show()
# # =========================================================
# 5) Interpret the missingness stats
# =========================================================
# For Y1 (bulk), printed results were:
# Total negative values: 0
# Total missing (before): 11
# Total values: 14229
# Rows with >50% missing: 0
# Rows fully missing: 0
# Interpretation:
# - No negative values, so nothing is being flagged as invalid for being below zero.
# - Only 11 missing values out of 14,229 total values.
#       11 / 14229 = ~0.000773 = 0.077%
# - No rows are mostly missing.
# - No rows are fully missing.
# Median imputation is applied to the 11 missing values

# For Y3 (cell), printed results were:
# Total negative values: 0
# Total missing (before): 0
# Total values: 14300
# Rows with >50% missing: 0
# Rows fully missing: 0
#
# Interpretation:
# - No missingness at all before filling.
# - No negative values.
# - The imputation step does nothing
# =========================================================
# 6) Set up counters to track imputation
# =========================================================
# - how many individual values were filled?
# - how many rows were completely missing and had to be filled with zero?
impute_counter = {"filled": 0, "all_missing_rows": 0}

# =========================================================
# 7) Fill missing or negative values with row median
# =========================================================
def fill_row_median(row):
    row_new = row.copy()
    # Identify negative values in this row.
    # These are treated as invalid and converted to missing.
    neg_mask = row < 0
    row_new.loc[neg_mask] = pd.NA
    # Count how many missing values this row has after marking negatives as missing.
    n_missing = row_new.isna().sum()
    # Add that number to the running total of imputed values.
    if n_missing > 0:
        impute_counter["filled"] += n_missing
    # fill missing values with the median of that lipid across samples.
        return row_new.fillna(row_new.median())
# Print a summary of how much imputation actually happened.
print("Total values imputed:", impute_counter["filled"])
# Y3 ONLY Y1 =0
# Total values imputed: 11
# Rows fully missing: 0

lipids
# =========================================================
# 9) Save output
# =========================================================
# Toggle whichever output file to write.
# lipids_norm.to_csv("Bulk_lipids_cleaned_normalized_median_527.csv")
# lipids_norm.to_csv("cell_lipids_cleaned_norm_median_286.csv")



#model train
# ===============================
# Load datasets
# ===============================
Y1 = pd.read_csv("/users/ramkd9/rna_lipid_clean/data/Bulk_lipids_cleaned_normalized_median_527.csv", index_col=0)
# Load bulk lipid data → rows = samples, columns = lipids

Y3 = pd.read_csv("/users/ramkd9/rna_lipid_clean/data/cell_lipids_cleaned_norm_median_286.csv", index_col=0)
# Load cell lipid data → rows = samples, columns = lipids
# X1: (19, 1348)
# X3 (after initial T): (150, 1849)
# Y1: (27, 527)
# Y3: (50, 286)
# ===============================
# Clean index and column names
# ===============================

for df in (X1, X3, Y1, Y3):
    df.index = df.index.astype(str).str.strip()
    # Ensure sample IDs are strings and remove leading/trailing whitespace

    df.columns = df.columns.astype(str).str.strip()
    # Ensure feature names (genes/lipids) are strings and remove whitespace


# ===============================
# Combine lipid datasets
# ===============================

Y_train = pd.concat([Y1, Y3], axis=0, join="inner")
# Stack bulk + cell lipid data by rows (samples)
# join="inner" → keep only lipids present in BOTH Y1 and Y3


# ===============================
# Align samples between RNA and lipid data
# ===============================

common_train = X_train.index.intersection(Y_train.index)
# Find sample IDs that exist in BOTH RNA (X_train) and lipid (Y_train)

# ===============================
# Subset to matched samples
# ==============================+

Y_train = Y_train.loc[common_train].copy()
# filtering for lipid data alignment with X_train

