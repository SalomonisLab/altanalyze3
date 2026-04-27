import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNet
from sklearn.metrics import r2_score
from sklearn.multioutput import MultiOutputRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import MultiTaskElasticNetCV
from tqdm import tqdm
import time


X1 = pd.read_csv("data/feature_blankreduiction.csv", index_col=0)
X3 = pd.read_csv("data/newrna_cell_clair_filtered_symbol.csv", index_col=0).T
Y1 = pd.read_csv("data/Bulk_lipids_cleaned_normalized_median_527.csv", index_col=0)
Y3 = pd.read_csv("data/cell_lipids_cleaned_norm_median_286.csv", index_col=0)

# collapse duplicate gene symbols
X3 = X3.T.groupby(level=0).mean().T

for df in (X1, X3, Y1, Y3):
    df.index = df.index.astype(str).str.strip()

# align samples
X_df = pd.concat([X1, X3], axis=0, join="inner")
Y_df = pd.concat([Y1, Y3], axis=0, join="inner")

X_df = X_df[~X_df.index.duplicated(keep="first")]
Y_df = Y_df[~Y_df.index.duplicated(keep="first")]

common_idx = X_df.index.intersection(Y_df.index)
X_df = X_df.loc[common_idx].copy()
Y_df = Y_df.loc[common_idx].copy()


# clean matrices
X_df = X_df.apply(pd.to_numeric, errors="coerce").fillna(0)
X_df = X_df.loc[:, (X_df != 0).any(axis=0)]

Y_df = Y_df.apply(pd.to_numeric, errors="coerce")
Y_df_clean = Y_df.dropna().copy()
X_df_clean = X_df.loc[Y_df_clean.index].copy()

print("Final X shape:", X_df_clean.shape)
print("Final Y shape:", Y_df_clean.shape)


# scale X and Y

scaler_x = StandardScaler()
X_scaled = scaler_x.fit_transform(X_df_clean)

scaler_y = StandardScaler()
Y_scaled = scaler_y.fit_transform(Y_df_clean)

# fit model
model = MultiTaskElasticNetCV(
    l1_ratio=[0.2, 0.5, 0.8],
    alphas=np.logspace(-3, 0, 10),
    cv=3,
    max_iter=30000,
    n_jobs=-1
)
print("Fitting MultiTaskElasticNetCV...")

t0 = time.time()
model.fit(X_scaled, Y_scaled)
print(f"Training finished in {time.time() - t0:.2f} seconds")

# enet_fs = MultiOutputRegressor(
#     ElasticNet(alpha=0.005, l1_ratio=0.7, max_iter=10000)
# )

# enet_fs.fit(X_scaled, Y_df_clean)

# Y_pred = enet_fs.predict(X_scaled)

# predict and inverse transform
Y_pred_scaled = model.predict(X_scaled)
Y_pred = scaler_y.inverse_transform(Y_pred_scaled)

pred_df = pd.DataFrame(
    Y_pred,
    index=Y_df_clean.index,
    columns=Y_df_clean.columns
)

print(pred_df.head())


# fast global metrics
true_vals = Y_df_clean.to_numpy().ravel()
pred_vals = pred_df.to_numpy().ravel()

safe_true = np.where(true_vals == 0, np.nan, true_vals)

rmse = np.sqrt(np.mean((true_vals - pred_vals) ** 2))
norm_rmse = rmse / np.mean(true_vals)
pearson_corr = np.corrcoef(true_vals, pred_vals)[0, 1]
mape = np.nanmean(np.abs((true_vals - pred_vals) / safe_true))

denom = np.sum((true_vals - np.mean(true_vals)) ** 2)
nse = 1 - (np.sum((true_vals - pred_vals) ** 2) / denom)

r = pearson_corr
alpha_kge = np.std(pred_vals) / np.std(true_vals)
beta_kge = np.mean(pred_vals) / np.mean(true_vals)
kge = 1 - np.sqrt((1 - r) ** 2 + (1 - alpha_kge) ** 2 + (1 - beta_kge) ** 2)

metrics = {
    "RMSE": rmse,
    "Normalized_RMSE": norm_rmse,
    "Pearson_r": pearson_corr,
    "MAPE": mape,
    "NSE": nse,
    "KGE": kge,
}

print(metrics)


pred_long = (
    pred_df.reset_index(names="sample")
    .melt(id_vars="sample", var_name="Lipid", value_name="Predicted")
)

true_long = (
    Y_df_clean.reset_index(names="sample")
    .melt(id_vars="sample", var_name="Lipid", value_name="True")
)

eval_df = pred_long.merge(true_long, on=["sample", "Lipid"], how="inner")
eval_df["error"] = eval_df["True"] - eval_df["Predicted"]

eval_df["Organ"] = np.where(eval_df["sample"].str.startswith("D"), "lung", "serum")

suffix_to_celltype = {
    "_EPI": "EPI",
    "_END": "END",
    "_MIC": "MIC",
    "_PMX": "PMX",
    "_MES": "MES",
}

eval_df["CellType"] = "non_cell_type"
for suffix, label in suffix_to_celltype.items():
    eval_df.loc[eval_df["sample"].str.endswith(suffix), "CellType"] = label

r2_per_lipid = (
    eval_df.groupby("Lipid")
    .apply(lambda df: r2_score(df["True"], df["Predicted"]))
    .rename("R2")
)


import pickle

bundle = {
    "model": model,
    "scaler_x": scaler_x,
    "scaler_y": scaler_y,
    "X_columns": X_df_clean.columns.tolist(),
    "Y_columns": Y_df_clean.columns.tolist(),
}

with open("otherelastic_multitask_try.pkl", "wb") as f:
    pickle.dump(bundle, f)


# bundle = {
#     "scaler": scaler,
#     "enet_fs": enet_fs,
#     "X_columns": X_df_clean.columns.tolist(),
#     "Y_columns": Y_df_clean.columns.tolist(),
# }   
# with open("elasticnet_model.pkl", "wb") as f:
#     pickle.dump(bundle, f)