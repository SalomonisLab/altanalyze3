

# SPARSE LIPID-BY-LIPID ELASTIC NET RNA-TO-LIPID MODEL
#
#   1. Scale RNA and lipid data using training samples only
#   2. Rank RNA genes separately for each lipid
#   3. Test multiple top-gene set sizes
#   4. Fit ElasticNetCV for each candidate gene set
#   5. Select the model balancing CV performance and sparsity
#   6. Evaluate on an untouched holdout set
#   7. Refit the selected architecture on all samples
#   8. Save a prediction-compatible model bundle


import os
import time
import pickle
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tqdm import tqdm
from scipy.stats import pearsonr, spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import (
    mean_squared_error,
    mean_absolute_error,
    r2_score
)

warnings.filterwarnings("ignore")

RANDOM_SEED = 1

TOP_GENE_OPTIONS = [
    25,
    50,
    100,
    200
]

L1_RATIO_GRID = [
    0.70,
    0.90,
    0.95,
    1.00
]

ALPHA_GRID = np.logspace(
    -3,
    1,
    30
)

CV_FOLDS = 3
MAX_ITER = 50000
N_JOBS = -1

# Penalizes retaining too many nonzero RNA coefficients.
SPARSITY_PENALTY = 0.002

OUTPUT_DIR = (
    "sparse_lipidwise_elasticnet_outputs"
)

os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)



# 2. LOAD DATA


X1 = pd.read_csv(
    "/users/ramkd9/Lipid_Predict/"
    "feature_blankreduiction.csv",
    index_col=0
)

X3 = pd.read_csv(
    "/users/ramkd9/Lipid_Predict/"
    "newrna_cell_clair_filtered_symbol.csv",
    index_col=0
).T

Y1 = pd.read_csv(
    "/users/ramkd9/Lipid_Predict/Clair/"
    "Bulk_lipids_cleaned_normalized_median_527_log2.csv",
    index_col=0
)

Y3 = pd.read_csv(
    "/users/ramkd9/Lipid_Predict/Clair/"
    "cell_lipids_cleaned_norm_median_286_log2.csv",
    index_col=0
)



# 3. CLEAN AND ALIGN MATRICES


# Collapse duplicate RNA gene symbols.
X3 = (
    X3.T
    .groupby(level=0)
    .mean()
    .T
)

for dataframe in [
    X1,
    X3,
    Y1,
    Y3
]:
    dataframe.index = (
        dataframe.index
        .astype(str)
        .str.strip()
    )

    dataframe.columns = (
        dataframe.columns
        .astype(str)
        .str.strip()
    )


# Combine bulk and sorted-cell RNA matrices.
X_df = pd.concat(
    [
        X1,
        X3
    ],
    axis=0,
    join="inner"
)

# Combine corresponding lipid matrices.
Y_df = pd.concat(
    [
        Y1,
        Y3
    ],
    axis=0,
    join="inner"
)


# Remove duplicated sample IDs.
X_df = X_df.loc[
    ~X_df.index.duplicated(
        keep="first"
    )
].copy()

Y_df = Y_df.loc[
    ~Y_df.index.duplicated(
        keep="first"
    )
].copy()


# Retain matched RNA and lipid samples.
common_samples = (
    X_df.index
    .intersection(
        Y_df.index
    )
)

X_df = X_df.loc[
    common_samples
].copy()

Y_df = Y_df.loc[
    common_samples
].copy()


# Convert all values to numeric.
X_df = X_df.apply(
    pd.to_numeric,
    errors="coerce"
)

Y_df = Y_df.apply(
    pd.to_numeric,
    errors="coerce"
)


# Fill missing values using feature medians.
X_df = X_df.fillna(
    X_df.median(
        numeric_only=True
    )
)

Y_df = Y_df.fillna(
    Y_df.median(
        numeric_only=True
    )
)


# Remove RNA genes containing no information.
X_df = X_df.loc[
    :,
    X_df.notna().sum(axis=0) > 0
]

X_df = X_df.loc[
    :,
    (X_df != 0).any(axis=0)
]


# Final cleaned matrices.
X_df_clean = X_df.copy()
Y_df_clean = Y_df.copy()


print("=" * 70)
print("CLEANED DATA")
print("=" * 70)

print(
    "Final RNA matrix:",
    X_df_clean.shape
)

print(
    "Final lipid matrix:",
    Y_df_clean.shape
)

print(
    "Matched samples:",
    len(
        X_df_clean.index
    )
)



# 4. BUILD SAMPLE METADATA


meta = pd.DataFrame(
    index=Y_df_clean.index.copy()
)

meta.index = (
    meta.index
    .astype(str)
    .str.strip()
)

meta["Donor"] = (
    meta.index
    .str.split("_")
    .str[0]
)

meta["Label"] = (
    meta.index
    .str.split(
        "_",
        n=1
    )
    .str[1]
    .fillna("")
)


def parse_celltype(label):
    label = (
        str(label)
        .strip()
        .lower()
    )

    if (
        label in {
            "bulk",
            "lung",
            "wholelung",
            "whole_lung"
        }
        or "bulk" in label
    ):
        return "BULK"

    if "epi" in label:
        return "EPI"

    if "mes" in label:
        return "MES"

    if "mic" in label:
        return "MIC"

    if "end" in label:
        return "END"

    if (
        "pmx" in label
        or "pmn" in label
    ):
        return "PMX"

    return "UNKNOWN"


meta["CellType"] = (
    meta["Label"]
    .apply(
        parse_celltype
    )
)

required_celltypes = [
    "EPI",
    "MES",
    "MIC",
    "END",
    "PMX"
]

meta = meta.loc[
    meta["CellType"].isin(
        required_celltypes
    )
].copy()


# Restrict matrices to usable samples.
X_df_clean = X_df_clean.loc[
    meta.index
].copy()

Y_df_clean = Y_df_clean.loc[
    meta.index
].copy()


print("\nUsable sample types:")
print(
    meta["CellType"]
    .value_counts()
)

print(
    "\nUsable samples:",
    len(meta)
)



# 5. CREATE CONSTRAINED HOLDOUT SPLIT


def generate_holdout_constrained(
    metadata,
    seed=42,
    total_holdout=28,
    min_per_celltype=5,
    max_per_donor=3,
    required_celltypes=(
        "EPI",
        "MES",
        "MIC",
        "END",
        "PMX"
    ),
    max_tries=10000
):
    rng = np.random.default_rng(
        seed
    )

    sample_ids = (
        metadata.index
        .tolist()
    )

    for _ in range(max_tries):

        donor_counts = {}
        holdout_samples = []
        valid_split = True

        # First guarantee minimum representation
        # from each required cell type.
        for celltype in required_celltypes:

            celltype_pool = (
                metadata.index[
                    metadata["CellType"]
                    == celltype
                ]
                .tolist()
            )

            rng.shuffle(
                celltype_pool
            )

            selected = []

            for sample_id in celltype_pool:

                donor = metadata.at[
                    sample_id,
                    "Donor"
                ]

                donor_count = (
                    donor_counts.get(
                        donor,
                        0
                    )
                )

                if donor_count < max_per_donor:

                    selected.append(
                        sample_id
                    )

                    donor_counts[donor] = (
                        donor_count + 1
                    )

                if (
                    len(selected)
                    == min_per_celltype
                ):
                    break

            if (
                len(selected)
                < min_per_celltype
            ):
                valid_split = False
                break

            holdout_samples.extend(
                selected
            )

        if not valid_split:
            continue

        remaining_needed = (
            total_holdout
            - len(holdout_samples)
        )

        remaining_pool = [
            sample_id
            for sample_id in sample_ids
            if sample_id
            not in holdout_samples
        ]

        rng.shuffle(
            remaining_pool
        )

        for sample_id in remaining_pool:

            if remaining_needed == 0:
                break

            donor = metadata.at[
                sample_id,
                "Donor"
            ]

            donor_count = (
                donor_counts.get(
                    donor,
                    0
                )
            )

            if donor_count < max_per_donor:

                holdout_samples.append(
                    sample_id
                )

                donor_counts[donor] = (
                    donor_count + 1
                )

                remaining_needed -= 1

        if (
            len(holdout_samples)
            != total_holdout
        ):
            continue

        holdout_metadata = metadata.loc[
            holdout_samples
        ].copy()

        # Confirm donor constraint.
        if (
            holdout_metadata[
                "Donor"
            ]
            .value_counts()
            .max()
            > max_per_donor
        ):
            continue

        # Confirm cell-type constraint.
        counts = (
            holdout_metadata[
                "CellType"
            ]
            .value_counts()
        )

        valid_celltypes = all(
            counts.get(
                celltype,
                0
            )
            >= min_per_celltype

            for celltype
            in required_celltypes
        )

        if not valid_celltypes:
            continue

        train_metadata = metadata.drop(
            index=holdout_samples
        ).copy()

        return (
            holdout_samples,
            train_metadata,
            holdout_metadata
        )

    raise RuntimeError(
        "Could not generate a valid "
        "constrained holdout split."
    )


(
    holdout_samples,
    meta_train,
    meta_holdout
) = generate_holdout_constrained(
    metadata=meta,
    seed=42,
    total_holdout=28,
    min_per_celltype=5,
    max_per_donor=3
)


meta_train.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "training_metadata.csv"
    )
)

meta_holdout.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "holdout_metadata.csv"
    )
)


X_train = X_df_clean.loc[
    meta_train.index
].copy()

Y_train = Y_df_clean.loc[
    meta_train.index
].copy()

X_holdout = X_df_clean.loc[
    meta_holdout.index
].copy()

Y_holdout = Y_df_clean.loc[
    meta_holdout.index
].copy()


print("\n" + "=" * 70)
print("TRAIN/HOLDOUT SPLIT")
print("=" * 70)

print(
    "Training RNA:",
    X_train.shape
)

print(
    "Training lipids:",
    Y_train.shape
)

print(
    "Holdout RNA:",
    X_holdout.shape
)

print(
    "Holdout lipids:",
    Y_holdout.shape
)

print("\nHoldout cell types:")
print(
    meta_holdout[
        "CellType"
    ]
    .value_counts()
)



# 6. HELPER FUNCTIONS


def safe_pearson(
    observed,
    predicted
):
    observed = np.asarray(
        observed,
        dtype=float
    )

    predicted = np.asarray(
        predicted,
        dtype=float
    )

    valid = (
        np.isfinite(observed)
        & np.isfinite(predicted)
    )

    observed = observed[
        valid
    ]

    predicted = predicted[
        valid
    ]

    if len(observed) < 3:
        return np.nan

    if (
        np.std(observed) == 0
        or np.std(predicted) == 0
    ):
        return np.nan

    return pearsonr(
        observed,
        predicted
    )[0]


def safe_spearman(
    observed,
    predicted
):
    observed = np.asarray(
        observed,
        dtype=float
    )

    predicted = np.asarray(
        predicted,
        dtype=float
    )

    valid = (
        np.isfinite(observed)
        & np.isfinite(predicted)
    )

    observed = observed[
        valid
    ]

    predicted = predicted[
        valid
    ]

    if len(observed) < 3:
        return np.nan

    if (
        np.std(observed) == 0
        or np.std(predicted) == 0
    ):
        return np.nan

    return spearmanr(
        observed,
        predicted
    )[0]


def calculate_global_metrics(
    y_true,
    y_pred,
    training_seconds
):
    observed = (
        y_true
        .to_numpy()
        .ravel()
    )

    predicted = (
        y_pred
        .to_numpy()
        .ravel()
    )

    valid = (
        np.isfinite(observed)
        & np.isfinite(predicted)
    )

    observed = observed[
        valid
    ]

    predicted = predicted[
        valid
    ]

    true_sd = np.std(
        observed
    )

    predicted_sd = np.std(
        predicted
    )

    return {
        "Model":
            "SparseLipidwiseElasticNetCV",

        "N_samples":
            y_true.shape[0],

        "N_lipids":
            y_true.shape[1],

        "N_values":
            len(observed),

        "RMSE":
            np.sqrt(
                mean_squared_error(
                    observed,
                    predicted
                )
            ),

        "MAE":
            mean_absolute_error(
                observed,
                predicted
            ),

        "Pearson_r":
            safe_pearson(
                observed,
                predicted
            ),

        "Spearman_r":
            safe_spearman(
                observed,
                predicted
            ),

        "R2":
            r2_score(
                observed,
                predicted
            ),

        "True_SD":
            true_sd,

        "Predicted_SD":
            predicted_sd,

        "SD_Ratio":
            (
                predicted_sd
                / true_sd
                if true_sd != 0
                else np.nan
            ),

        "Training_seconds":
            training_seconds
    }



# FIT SPARSE LIPID-BY-LIPID ELASTIC NET


def fit_sparse_lipidwise_elasticnet(
    x_train,
    y_train,
    x_test,
    top_gene_options,
    l1_ratio_grid,
    alpha_grid,
    cv_folds,
    max_iter,
    sparsity_penalty,
    random_seed
):
    # Scale RNA data using training samples only.

    scaler_x = StandardScaler()

    x_train_scaled = pd.DataFrame(
        scaler_x.fit_transform(
            x_train
        ),
        index=x_train.index,
        columns=x_train.columns
    )

    x_test_scaled = pd.DataFrame(
        scaler_x.transform(
            x_test
        ),
        index=x_test.index,
        columns=x_test.columns
    )

    # Scale lipid data using training samples only
    scaler_y = StandardScaler()

    y_train_scaled = pd.DataFrame(
        scaler_y.fit_transform(
            y_train
        ),
        index=y_train.index,
        columns=y_train.columns
    )

    all_models = {}
    predictions_scaled = pd.DataFrame(
        index=x_test.index
    )

    coefficient_records = []
    summary_records = []

    start_time = time.perf_counter()

    for lipid in tqdm(
        y_train_scaled.columns,
        desc=(
            "Training sparse "
            "lipid-by-lipid models"
        )
    ):
        y = y_train_scaled[
            lipid
        ]

        # Rank RNA genes by training-set correlation.

        correlations = (
            x_train_scaled
            .apply(
                lambda feature:
                    feature.corr(y),
                axis=0
            )
            .replace(
                [
                    np.inf,
                    -np.inf
                ],
                np.nan
            )
            .fillna(0)
        )

        ranked_genes = (
            correlations
            .abs()
            .sort_values(
                ascending=False
            )
        )

        best_model = None
        best_genes = None
        best_test_prediction = None
        best_training_prediction = None
        best_score = -np.inf
        best_top_n = None

        candidate_records = []


        for top_n in top_gene_options:

            actual_top_n = min(
                top_n,
                x_train_scaled.shape[1]
            )

            selected_genes = (
                ranked_genes
                .head(actual_top_n)
                .index
                .tolist()
            )

            x_train_subset = (
                x_train_scaled[
                    selected_genes
                ]
            )

            x_test_subset = (
                x_test_scaled[
                    selected_genes
                ]
            )

            model = ElasticNetCV(
                l1_ratio=
                    l1_ratio_grid,

                alphas=
                    alpha_grid,

                cv=
                    cv_folds,

                max_iter=
                    max_iter,

                n_jobs=
                    N_JOBS,

                random_state=
                    random_seed,

                selection=
                    "cyclic"
            )

            model.fit(
                x_train_subset,
                y
            )

            training_prediction = (
                model.predict(
                    x_train_subset
                )
            )

            test_prediction = (
                model.predict(
                    x_test_subset
                )
            )

            training_r2 = r2_score(
                y,
                training_prediction
            )

            number_nonzero = int(
                np.sum(
                    model.coef_ != 0
                )
            )

            # Preserve your original model-selection rule.
            sparsity_adjusted_score = (
                training_r2
                - sparsity_penalty
                * number_nonzero
            )

            candidate_records.append({
                "Lipid":
                    lipid,

                "Top_N":
                    actual_top_n,

                "Alpha":
                    float(
                        model.alpha_
                    ),

                "L1_ratio":
                    float(
                        model.l1_ratio_
                    ),

                "Nonzero_Coefficients":
                    number_nonzero,

                "Train_R2_scaled":
                    training_r2,

                "Sparsity_Adjusted_Score":
                    sparsity_adjusted_score
            })

            if (
                sparsity_adjusted_score
                > best_score
            ):
                best_score = (
                    sparsity_adjusted_score
                )

                best_model = model

                best_genes = (
                    selected_genes
                )

                best_training_prediction = (
                    training_prediction
                )

                best_test_prediction = (
                    test_prediction
                )

                best_top_n = (
                    actual_top_n
                )


        all_models[lipid] = {
            "model":
                best_model,

            "genes":
                best_genes,

            "top_n":
                best_top_n,

            "selected_alpha":
                float(
                    best_model.alpha_
                ),

            "selected_l1_ratio":
                float(
                    best_model.l1_ratio_
                ),

            "sparsity_adjusted_score":
                float(
                    best_score
                )
        }

        predictions_scaled[
            lipid
        ] = best_test_prediction

        nonzero_mask = (
            best_model.coef_
            != 0
        )

        nonzero_genes = (
            np.asarray(
                best_genes
            )[
                nonzero_mask
            ]
        )

        nonzero_coefficients = (
            best_model.coef_[
                nonzero_mask
            ]
        )

        for (
            gene,
            coefficient
        ) in zip(
            nonzero_genes,
            nonzero_coefficients
        ):
            coefficient_records.append({
                "Lipid":
                    lipid,

                "Gene":
                    gene,

                "Coefficient":
                    coefficient,

                "Abs_Coefficient":
                    abs(coefficient),

                "Correlation_with_lipid":
                    correlations.loc[
                        gene
                    ]
            })

        training_pearson = (
            safe_pearson(
                y,
                best_training_prediction
            )
        )

        summary_records.append({
            "Lipid":
                lipid,

            "Top_N_Correlation_Filter":
                best_top_n,

            "Candidate_Top_N_Options":
                ",".join(
                    map(
                        str,
                        top_gene_options
                    )
                ),

            "Nonzero_Coefficients":
                int(
                    len(
                        nonzero_genes
                    )
                ),

            "Alpha":
                float(
                    best_model.alpha_
                ),

            "L1_ratio":
                float(
                    best_model.l1_ratio_
                ),

            "Train_R2_scaled":
                r2_score(
                    y,
                    best_training_prediction
                ),

            "Train_Pearson_scaled":
                training_pearson,

            "Sparsity_Adjusted_Score":
                best_score
        })

    training_seconds = (
        time.perf_counter()
        - start_time
    )

    # Ensure exact lipid order.
    predictions_scaled = (
        predictions_scaled[
            y_train.columns
        ]
    )

    # Return predictions to original lipid scale.
    predictions = (
        scaler_y.inverse_transform(
            predictions_scaled
        )
    )

    prediction_df = pd.DataFrame(
        predictions,
        index=x_test.index,
        columns=y_train.columns
    )

    summary_df = pd.DataFrame(
        summary_records
    )

    coefficient_df = pd.DataFrame(
        coefficient_records
    )

    candidate_df = pd.DataFrame(
        candidate_records
    )

    bundle = {
        "model_name":
            "SparseLipidwiseElasticNetCV",

        "architecture":
            (
                "Separate ElasticNetCV model "
                "for each lipid"
            ),

        "models":
            all_models,

        "scaler_x":
            scaler_x,

        "scaler_y":
            scaler_y,

        "X_columns":
            list(
                x_train.columns
            ),

        "Y_columns":
            list(
                y_train.columns
            ),

        "top_gene_options":
            list(
                top_gene_options
            ),

        "l1_ratio_grid":
            list(
                l1_ratio_grid
            ),

        "alpha_grid":
            list(
                alpha_grid
            ),

        "cv_folds":
            cv_folds,

        "max_iter":
            max_iter,

        "sparsity_penalty":
            sparsity_penalty,

        "random_seed":
            random_seed,

        "summary":
            summary_df,

        "coefficients":
            coefficient_df,

        "candidate_models":
            candidate_df
    }

    return (
        prediction_df,
        bundle,
        training_seconds
    )



# 8. TRAIN AND PREDICT HOLDOUT


(
    pred_holdout_df,
    holdout_bundle,
    training_seconds
) = fit_sparse_lipidwise_elasticnet(
    x_train=X_train,
    y_train=Y_train,
    x_test=X_holdout,
    top_gene_options=
        TOP_GENE_OPTIONS,
    l1_ratio_grid=
        L1_RATIO_GRID,
    alpha_grid=
        ALPHA_GRID,
    cv_folds=
        CV_FOLDS,
    max_iter=
        MAX_ITER,
    sparsity_penalty=
        SPARSITY_PENALTY,
    random_seed=
        RANDOM_SEED
)


print("\nTraining completed in:")
print(
    f"{training_seconds:.2f} seconds"
)



# 9. HOLDOUT METRICS


global_metrics = (
    calculate_global_metrics(
        y_true=Y_holdout,
        y_pred=pred_holdout_df,
        training_seconds=
            training_seconds
    )
)

global_metrics_df = pd.DataFrame(
    [
        global_metrics
    ]
)


per_lipid_records = []

for lipid in Y_holdout.columns:

    observed = (
        Y_holdout[
            lipid
        ]
        .to_numpy()
    )

    predicted = (
        pred_holdout_df[
            lipid
        ]
        .to_numpy()
    )

    valid = (
        np.isfinite(observed)
        & np.isfinite(predicted)
    )

    observed = observed[
        valid
    ]

    predicted = predicted[
        valid
    ]

    if len(observed) < 3:
        continue

    model_information = (
        holdout_bundle[
            "models"
        ][lipid]
    )

    per_lipid_records.append({
        "Lipid":
            lipid,

        "N":
            len(observed),

        "RMSE":
            np.sqrt(
                mean_squared_error(
                    observed,
                    predicted
                )
            ),

        "MAE":
            mean_absolute_error(
                observed,
                predicted
            ),

        "Pearson_r":
            safe_pearson(
                observed,
                predicted
            ),

        "Spearman_r":
            safe_spearman(
                observed,
                predicted
            ),

        "R2":
            r2_score(
                observed,
                predicted
            ),

        "True_SD":
            np.std(
                observed
            ),

        "Predicted_SD":
            np.std(
                predicted
            ),

        "Top_N":
            model_information[
                "top_n"
            ],

        "Alpha":
            model_information[
                "selected_alpha"
            ],

        "L1_ratio":
            model_information[
                "selected_l1_ratio"
            ],

        "Sparsity_Adjusted_Score":
            model_information[
                "sparsity_adjusted_score"
            ]
    })


per_lipid_metrics_df = pd.DataFrame(
    per_lipid_records
)


print("\n" + "=" * 70)
print("HOLDOUT PERFORMANCE")
print("=" * 70)

print(
    global_metrics_df
    .T
)



# 10. SAVE HOLDOUT OUTPUTS


pred_holdout_df.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "holdout_predictions.csv"
    )
)

Y_holdout.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "holdout_true_values.csv"
    )
)

global_metrics_df.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "holdout_global_metrics.csv"
    ),
    index=False
)

per_lipid_metrics_df.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "holdout_per_lipid_metrics.csv"
    ),
    index=False
)

holdout_bundle[
    "summary"
].to_csv(
    os.path.join(
        OUTPUT_DIR,
        "holdout_model_summary.csv"
    ),
    index=False
)

holdout_bundle[
    "coefficients"
].to_csv(
    os.path.join(
        OUTPUT_DIR,
        "holdout_nonzero_coefficients.csv"
    ),
    index=False
)

holdout_bundle[
    "candidate_models"
].to_csv(
    os.path.join(
        OUTPUT_DIR,
        "holdout_candidate_top_gene_models.csv"
    ),
    index=False
)

with open(
    os.path.join(
        OUTPUT_DIR,
        "holdout_model_bundle.pkl"
    ),
    "wb"
) as file:
    pickle.dump(
        holdout_bundle,
        file,
        protocol=
            pickle.HIGHEST_PROTOCOL
    )



#  HOLDOUT SCATTERPLOT


scatter_true = (
    Y_holdout
    .to_numpy()
    .ravel()
)

scatter_predicted = (
    pred_holdout_df
    .to_numpy()
    .ravel()
)

valid = (
    np.isfinite(
        scatter_true
    )
    & np.isfinite(
        scatter_predicted
    )
)

scatter_true = scatter_true[
    valid
]

scatter_predicted = (
    scatter_predicted[
        valid
    ]
)


fig, ax = plt.subplots(
    figsize=(
        7,
        7
    )
)

ax.scatter(
    scatter_true,
    scatter_predicted,
    alpha=0.35,
    s=15
)

minimum_value = min(
    scatter_true.min(),
    scatter_predicted.min()
)

maximum_value = max(
    scatter_true.max(),
    scatter_predicted.max()
)

ax.plot(
    [
        minimum_value,
        maximum_value
    ],
    [
        minimum_value,
        maximum_value
    ],
    linestyle="--",
    linewidth=1
)

ax.set_xlabel(
    "True lipid abundance"
)

ax.set_ylabel(
    "Predicted lipid abundance"
)

ax.set_title(
    "Sparse Lipid-by-Lipid Elastic Net\n"
    f"Pearson r = "
    f"{global_metrics['Pearson_r']:.3f}, "
    f"RMSE = "
    f"{global_metrics['RMSE']:.3f}, "
    f"MAE = "
    f"{global_metrics['MAE']:.3f}"
)

plt.tight_layout()

plt.savefig(
    os.path.join(
        OUTPUT_DIR,
        "holdout_true_vs_predicted_scatter.png"
    ),
    dpi=300,
    bbox_inches="tight"
)

plt.show()
plt.close()



# 12. SUMMARY OF CHOSEN TOP-GENE COUNTS


print("\nSelected top-gene counts:")

print(
    holdout_bundle[
        "summary"
    ][
        "Top_N_Correlation_Filter"
    ]
    .value_counts()
    .sort_index()
)

print("\nSelected L1 ratios:")

print(
    holdout_bundle[
        "summary"
    ][
        "L1_ratio"
    ]
    .value_counts()
    .sort_index()
)

print("\nNonzero coefficients per lipid:")

print(
    holdout_bundle[
        "summary"
    ][
        "Nonzero_Coefficients"
    ]
    .describe()
)



# The holdout model above measures performance.
#
# The model below is the final production model. It uses every
# available matched sample after the evaluation is completed.


(
    final_training_predictions,
    final_bundle,
    final_training_seconds
) = fit_sparse_lipidwise_elasticnet(
    x_train=X_df_clean,
    y_train=Y_df_clean,

    # Predicting the training samples here is only used so the
    # helper function can return a correctly shaped dataframe.
    # These values are not used as unbiased validation.
    x_test=X_df_clean,

    top_gene_options=
        TOP_GENE_OPTIONS,

    l1_ratio_grid=
        L1_RATIO_GRID,

    alpha_grid=
        ALPHA_GRID,

    cv_folds=
        CV_FOLDS,

    max_iter=
        MAX_ITER,

    sparsity_penalty=
        SPARSITY_PENALTY,

    random_seed=
        RANDOM_SEED
)


final_bundle[
    "training_samples"
] = list(
    X_df_clean.index
)

final_bundle[
    "training_metadata"
] = meta.copy()

final_bundle[
    "final_training_seconds"
] = (
    final_training_seconds
)

final_bundle[
    "validation_global_metrics"
] = (
    global_metrics
)


final_model_path = os.path.join(
    OUTPUT_DIR,
    "rna2lipid_sparse_lipidwise_"
    "ElasticNetCV_final_bundle.pkl"
)

with open(
    final_model_path,
    "wb"
) as file:
    pickle.dump(
        final_bundle,
        file,
        protocol=
            pickle.HIGHEST_PROTOCOL
    )


final_bundle[
    "summary"
].to_csv(
    os.path.join(
        OUTPUT_DIR,
        "final_model_summary.csv"
    ),
    index=False
)

final_bundle[
    "coefficients"
].to_csv(
    os.path.join(
        OUTPUT_DIR,
        "final_nonzero_coefficients.csv"
    ),
    index=False
)

final_bundle[
    "candidate_models"
].to_csv(
    os.path.join(
        OUTPUT_DIR,
        "final_candidate_top_gene_models.csv"
    ),
    index=False
)


print("\n" + "=" * 70)
print("FINAL MODEL SAVED")
print("=" * 70)

print(
    final_model_path
)

print(
    "Final model training time:",
    f"{final_training_seconds:.2f} seconds"
)

print(
    "\nFinal top-gene selections:"
)

print(
    final_bundle[
        "summary"
    ][
        "Top_N_Correlation_Filter"
    ]
    .value_counts()
    .sort_index()
)

print(
    "\nFinal nonzero coefficients:"
)

print(
    final_bundle[
        "summary"
    ][
        "Nonzero_Coefficients"
    ]
    .describe()
)













# import time
# import pickle
# from pathlib import Path

# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt

# from scipy.stats import pearsonr, spearmanr
# from sklearn.preprocessing import StandardScaler
# from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
# from sklearn.linear_model import MultiTaskLassoCV

# from sklearn.linear_model import MultiTaskElasticNetCV
# from sklearn.preprocessing import StandardScaler

# import time
# import numpy as np
# import pandas as pd

# ALPHAS = np.logspace(-3, 1, 30)

# L1_RATIOS = [
#     0.70,
#     0.90,
#     0.95,
#     1.00,
# ]

# CV_FOLDS = 3
# MAX_ITER = 50000
# N_JOBS = -1
# RANDOM_SEED = 1
# 
# # SETTINGS
# 

# X_BULK_PATH = Path(
#     "/users/ramkd9/Lipid_Predict/feature_blankreduiction.csv"
# )
# X_CELL_PATH = Path(
#     "/users/ramkd9/Lipid_Predict/newrna_cell_clair_filtered_symbol.csv"
# )
# Y_BULK_PATH = Path(
#     "/users/ramkd9/Lipid_Predict/Clair/"
#     "Bulk_lipids_cleaned_normalized_median_527_log2.csv"
# )
# Y_CELL_PATH = Path(
#     "/users/ramkd9/Lipid_Predict/Clair/"
#     "cell_lipids_cleaned_norm_median_286_log2.csv"
# )

# OUTPUT_DIR = Path("multitask_lasso_outputs")
# OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# RANDOM_SEED = 42
# TOTAL_HOLDOUT = 28
# MIN_PER_CELL_TYPE = 5
# MAX_PER_DONOR = 3
# REQUIRED_CELL_TYPES = ("EPI", "MES", "MIC", "END", "PMX")

# ALPHAS = np.logspace(-3, 1, 30)
# CV_FOLDS = 3
# MAX_ITER = 50000
# N_JOBS = -1


# 
# # HELPERS
# 

# def parse_celltype(label):
#     lab = str(label).strip().lower()

#     if lab in {"bulk", "lung", "wholelung", "whole_lung"} or "bulk" in lab:
#         return "BULK"
#     if "epi" in lab:
#         return "EPI"
#     if "mes" in lab:
#         return "MES"
#     if "mic" in lab:
#         return "MIC"
#     if "end" in lab:
#         return "END"
#     if "pmx" in lab or "pmn" in lab:
#         return "PMX"

#     return "UNKNOWN"


# def generate_holdout_constrained(
#     meta,
#     seed=RANDOM_SEED,
#     total_holdout=TOTAL_HOLDOUT,
#     min_per_ct=MIN_PER_CELL_TYPE,
#     max_per_donor=MAX_PER_DONOR,
#     required_cts=REQUIRED_CELL_TYPES,
#     max_tries=10000,
# ):
#     """Create a reproducible holdout with cell-type and donor constraints."""

#     rng = np.random.default_rng(seed)
#     sample_ids = meta.index.tolist()

#     for _ in range(max_tries):
#         donor_counts = {}
#         holdout = []
#         success = True

#         for cell_type in required_cts:
#             cell_pool = meta.index[
#                 meta["CellType"] == cell_type
#             ].tolist()
#             rng.shuffle(cell_pool)

#             chosen = []

#             for sample_id in cell_pool:
#                 donor = meta.at[sample_id, "Donor"]

#                 if donor_counts.get(donor, 0) < max_per_donor:
#                     chosen.append(sample_id)
#                     donor_counts[donor] = donor_counts.get(donor, 0) + 1

#                 if len(chosen) == min_per_ct:
#                     break

#             if len(chosen) < min_per_ct:
#                 success = False
#                 break

#             holdout.extend(chosen)

#         if not success:
#             continue

#         remaining_needed = total_holdout - len(holdout)
#         remaining_pool = [
#             sample_id
#             for sample_id in sample_ids
#             if sample_id not in holdout
#         ]
#         rng.shuffle(remaining_pool)

#         extras = []

#         for sample_id in remaining_pool:
#             donor = meta.at[sample_id, "Donor"]

#             if donor_counts.get(donor, 0) < max_per_donor:
#                 extras.append(sample_id)
#                 donor_counts[donor] = donor_counts.get(donor, 0) + 1

#             if len(extras) == remaining_needed:
#                 break

#         holdout_final = holdout + extras

#         if len(holdout_final) != total_holdout:
#             continue

#         holdout_meta = meta.loc[holdout_final].copy()

#         if holdout_meta["Donor"].value_counts().max() > max_per_donor:
#             continue

#         if any(
#             holdout_meta["CellType"].value_counts().get(cell_type, 0)
#             < min_per_ct
#             for cell_type in required_cts
#         ):
#             continue

#         train_meta = meta.drop(index=holdout_final).copy()
#         return holdout_final, train_meta, holdout_meta

#     raise RuntimeError("Could not generate a valid constrained holdout split.")


# def safe_correlation(x, y, method="pearson"):
#     x = np.asarray(x, dtype=float)
#     y = np.asarray(y, dtype=float)

#     mask = np.isfinite(x) & np.isfinite(y)
#     x = x[mask]
#     y = y[mask]

#     if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
#         return np.nan

#     if method == "spearman":
#         return spearmanr(x, y).statistic

#     return pearsonr(x, y).statistic


# def compute_global_metrics(y_true, y_pred, train_seconds):
#     true_values = y_true.to_numpy(dtype=float).ravel()
#     pred_values = y_pred.to_numpy(dtype=float).ravel()

#     mask = np.isfinite(true_values) & np.isfinite(pred_values)
#     true_values = true_values[mask]
#     pred_values = pred_values[mask]

#     true_sd = np.std(true_values)
#     pred_sd = np.std(pred_values)

#     return pd.DataFrame([{
#         "Model": "MultiTaskLassoCV",
#         "N_samples": y_true.shape[0],
#         "N_lipids": y_true.shape[1],
#         "N_values": len(true_values),
#         "RMSE": np.sqrt(mean_squared_error(true_values, pred_values)),
#         "MAE": mean_absolute_error(true_values, pred_values),
#         "Pearson_r": safe_correlation(true_values, pred_values, "pearson"),
#         "Spearman_r": safe_correlation(true_values, pred_values, "spearman"),
#         "R2": r2_score(true_values, pred_values),
#         "True_SD": true_sd,
#         "Pred_SD": pred_sd,
#         "SD_ratio_pred_true": (
#             pred_sd / true_sd if true_sd != 0 else np.nan
#         ),
#         "Train_seconds": train_seconds,
#     }])


# def compute_per_lipid_metrics(y_true, y_pred):
#     records = []

#     for lipid in y_true.columns:
#         true_values = y_true[lipid].to_numpy(dtype=float)
#         pred_values = y_pred[lipid].to_numpy(dtype=float)

#         mask = np.isfinite(true_values) & np.isfinite(pred_values)
#         true_values = true_values[mask]
#         pred_values = pred_values[mask]

#         if len(true_values) < 3:
#             continue

#         records.append({
#             "Model": "MultiTaskLassoCV",
#             "Lipid": lipid,
#             "N": len(true_values),
#             "RMSE": np.sqrt(mean_squared_error(true_values, pred_values)),
#             "MAE": mean_absolute_error(true_values, pred_values),
#             "Pearson_r": safe_correlation(true_values, pred_values, "pearson"),
#             "Spearman_r": safe_correlation(true_values, pred_values, "spearman"),
#             "R2": r2_score(true_values, pred_values),
#             "True_SD": np.std(true_values),
#             "Pred_SD": np.std(pred_values),
#         })

#     return pd.DataFrame(records)


# # def fit_multitask_elasticnet_sparse(x_train, y_train, x_test):
# #     scaler_x = StandardScaler()
# #     scaler_y = StandardScaler()

# #     x_train_scaled = scaler_x.fit_transform(x_train)
# #     x_test_scaled = scaler_x.transform(x_test)
# #     y_train_scaled = scaler_y.fit_transform(y_train)

# #     model = MultiTaskLassoCV(
# #         alphas=ALPHAS,
# #         cv=CV_FOLDS,
# #         max_iter=MAX_ITER,
# #         n_jobs=N_JOBS,
# #         random_state=RANDOM_SEED,
# #     )

# #     start = time.perf_counter()
# #     model.fit(x_train_scaled, y_train_scaled)
# #     train_seconds = time.perf_counter() - start

# #     pred_scaled = model.predict(x_test_scaled)
# #     pred = scaler_y.inverse_transform(pred_scaled)

# #     pred_df = pd.DataFrame(
# #         pred,
# #         index=x_test.index,
# #         columns=y_train.columns,
# #     )

# #     bundle = {
# #         "model_name": "MultiTaskLassoCV",
# #         "model": model,
# #         "scaler_x": scaler_x,
# #         "scaler_y": scaler_y,
# #         "X_columns": list(x_train.columns),
# #         "Y_columns": list(y_train.columns),
# #         "alpha_grid": list(ALPHAS),
# #         "selected_alpha": float(model.alpha_),
# #         "cv_folds": CV_FOLDS,
# #         "max_iter": MAX_ITER,
# #     }

# #     return pred_df, bundle, train_seconds

# def fit_multitask_elasticnet_sparse(
#     x_train,
#     y_train,
#     x_test,
# ):
#     """
#     Train a sparse MultiTaskElasticNetCV model.

#     The model predicts all lipid outcomes simultaneously and
#     encourages shared RNA feature selection across lipids.

#     L1_RATIOS controls the mixture of:
#         1.0  = pure multitask Lasso
#         <1.0 = multitask Lasso plus Ridge stabilization
#     """

#     # ========================================================
#     # SCALE RNA FEATURES
#     # ========================================================

#     scaler_x = StandardScaler()

#     x_train_scaled = scaler_x.fit_transform(
#         x_train
#     )

#     x_test_scaled = scaler_x.transform(
#         x_test
#     )

#     # ========================================================
#     # SCALE LIPID TARGETS
#     # ========================================================

#     scaler_y = StandardScaler()

#     y_train_scaled = scaler_y.fit_transform(
#         y_train
#     )

#     # ========================================================
#     # DEFINE SPARSE MULTITASK ELASTIC NET
#     # ========================================================

#     model = MultiTaskElasticNetCV(
#         l1_ratio=L1_RATIOS,
#         alphas=ALPHAS,
#         cv=CV_FOLDS,
#         max_iter=MAX_ITER,
#         n_jobs=N_JOBS,
#         random_state=RANDOM_SEED,
#         selection="cyclic",
#     )

#     # ========================================================
#     # TRAIN MODEL
#     # ========================================================

#     start = time.perf_counter()

#     model.fit(
#         x_train_scaled,
#         y_train_scaled
#     )

#     train_seconds = (
#         time.perf_counter() - start
#     )

#     # ========================================================
#     # PREDICT HOLDOUT DATA
#     # ========================================================

#     pred_scaled = model.predict(
#         x_test_scaled
#     )

#     pred = scaler_y.inverse_transform(
#         pred_scaled
#     )

#     pred_df = pd.DataFrame(
#         pred,
#         index=x_test.index,
#         columns=y_train.columns,
#     )

#     # ========================================================
#     # CALCULATE SPARSITY INFORMATION
#     # ========================================================

#     # coef_ shape:
#     # number of lipids × number of RNA features
#     coefficient_matrix = model.coef_

#     # A gene is retained if it has a nonzero coefficient
#     # for at least one lipid.
#     retained_gene_mask = np.any(
#         coefficient_matrix != 0,
#         axis=0
#     )

#     retained_genes = (
#         np.asarray(x_train.columns)[
#             retained_gene_mask
#         ]
#         .astype(str)
#         .tolist()
#     )

#     removed_genes = (
#         np.asarray(x_train.columns)[
#             ~retained_gene_mask
#         ]
#         .astype(str)
#         .tolist()
#     )

#     number_nonzero_coefficients = int(
#         np.sum(coefficient_matrix != 0)
#     )

#     number_total_coefficients = int(
#         coefficient_matrix.size
#     )

#     coefficient_sparsity = (
#         1
#         - (
#             number_nonzero_coefficients
#             / number_total_coefficients
#         )
#     )

#     # ========================================================
#     # SAVE MODEL BUNDLE
#     # ========================================================

#     bundle = {
#         "model_name": "MultiTaskElasticNetCV_sparse",
#         "model": model,
#         "scaler_x": scaler_x,
#         "scaler_y": scaler_y,

#         "X_columns": list(
#             x_train.columns
#         ),

#         "Y_columns": list(
#             y_train.columns
#         ),

#         "alpha_grid": list(
#             ALPHAS
#         ),

#         "l1_ratio_grid": list(
#             L1_RATIOS
#         ),

#         "selected_alpha": float(
#             model.alpha_
#         ),

#         "selected_l1_ratio": float(
#             model.l1_ratio_
#         ),

#         "cv_folds": CV_FOLDS,
#         "max_iter": MAX_ITER,

#         "retained_genes": retained_genes,
#         "removed_genes": removed_genes,

#         "number_input_genes": int(
#             x_train.shape[1]
#         ),

#         "number_retained_genes": int(
#             len(retained_genes)
#         ),

#         "number_removed_genes": int(
#             len(removed_genes)
#         ),

#         "number_nonzero_coefficients":
#             number_nonzero_coefficients,

#         "number_total_coefficients":
#             number_total_coefficients,

#         "coefficient_sparsity": float(
#             coefficient_sparsity
#         ),
#     }

#     return (
#         pred_df,
#         bundle,
#         train_seconds,
#     )
# 
# # 1. LOAD AND ALIGN DATA
# 

# x_bulk = pd.read_csv(X_BULK_PATH, index_col=0)
# x_cell = pd.read_csv(X_CELL_PATH, index_col=0).T
# y_bulk = pd.read_csv(Y_BULK_PATH, index_col=0)
# y_cell = pd.read_csv(Y_CELL_PATH, index_col=0)

# # Collapse duplicate gene symbols in the cell-level RNA matrix.
# x_cell = x_cell.T.groupby(level=0).mean().T

# for dataframe in (x_bulk, x_cell, y_bulk, y_cell):
#     dataframe.index = dataframe.index.astype(str).str.strip()
#     dataframe.columns = dataframe.columns.astype(str).str.strip()

# x_df = pd.concat([x_bulk, x_cell], axis=0, join="inner")
# y_df = pd.concat([y_bulk, y_cell], axis=0, join="inner")

# x_df = x_df.loc[~x_df.index.duplicated(keep="first")].copy()
# y_df = y_df.loc[~y_df.index.duplicated(keep="first")].copy()

# common_samples = x_df.index.intersection(y_df.index)
# x_df = x_df.loc[common_samples].copy()
# y_df = y_df.loc[common_samples].copy()

# x_df = x_df.apply(pd.to_numeric, errors="coerce")
# y_df = y_df.apply(pd.to_numeric, errors="coerce")

# # Impute each feature with its training-table median.
# x_df = x_df.fillna(x_df.median(numeric_only=True))
# y_df = y_df.fillna(y_df.median(numeric_only=True))

# # Remove invalid or all-zero RNA features.
# x_df = x_df.loc[:, x_df.notna().sum(axis=0) > 0]
# x_df = x_df.loc[:, (x_df != 0).any(axis=0)]

# # Remove invalid lipid columns if any remain.
# y_df = y_df.loc[:, y_df.notna().sum(axis=0) > 0]

# print("Aligned RNA shape:", x_df.shape)
# print("Aligned lipid shape:", y_df.shape)


# 
# # 2. BUILD METADATA AND HOLDOUT SPLIT
# 

# meta = pd.DataFrame(index=y_df.index.copy())
# meta.index = meta.index.astype(str).str.strip()
# meta["Donor"] = meta.index.str.split("_").str[0]
# meta["Label"] = meta.index.str.split("_", n=1).str[1].fillna("")
# meta["CellType"] = meta["Label"].apply(parse_celltype)

# # The selected benchmark evaluates the five sorted cell populations.
# meta = meta.loc[
#     meta["CellType"].isin(REQUIRED_CELL_TYPES)
# ].copy()

# x_df = x_df.loc[meta.index].copy()
# y_df = y_df.loc[meta.index].copy()

# print("Usable samples:", len(meta))
# print(meta["CellType"].value_counts())

# holdout_samples, meta_train, meta_holdout = generate_holdout_constrained(meta)

# meta_train.to_csv(OUTPUT_DIR / "train_metadata.csv")
# meta_holdout.to_csv(OUTPUT_DIR / "holdout_metadata.csv")
# pd.Series(
#     holdout_samples,
#     name="SampleID",
# ).to_csv(
#     OUTPUT_DIR / "holdout_samples.csv",
#     index=False,
# )

# x_train = x_df.loc[meta_train.index].copy()
# y_train = y_df.loc[meta_train.index].copy()
# x_holdout = x_df.loc[meta_holdout.index].copy()
# y_holdout = y_df.loc[meta_holdout.index].copy()

# print("Training shapes:", x_train.shape, y_train.shape)
# print("Holdout shapes:", x_holdout.shape, y_holdout.shape)


# 
# # 3. FIT AND EVALUATE MULTITASK LASSO
# 

# holdout_pred, holdout_bundle, train_seconds = fit_multitask_elasticnet_sparse(
#     x_train,
#     y_train,
#     x_holdout,
# )

# holdout_pred.to_csv(
#     OUTPUT_DIR / "MultiTaskLassoCV_holdout_predictions.csv"
# )

# global_metrics = compute_global_metrics(
#     y_holdout,
#     holdout_pred,
#     train_seconds,
# )
# per_lipid_metrics = compute_per_lipid_metrics(
#     y_holdout,
#     holdout_pred,
# )

# global_metrics.to_csv(
#     OUTPUT_DIR / "MultiTaskLassoCV_global_metrics.csv",
#     index=False,
# )
# per_lipid_metrics.to_csv(
#     OUTPUT_DIR / "MultiTaskLassoCV_per_lipid_metrics.csv",
#     index=False,
# )

# with open(
#     OUTPUT_DIR / "MultiTaskLassoCV_holdout_bundle.pkl",
#     "wb",
# ) as handle:
#     pickle.dump(
#         holdout_bundle,
#         handle,
#         protocol=pickle.HIGHEST_PROTOCOL,
#     )

# print("\nHoldout metrics:")
# print(global_metrics.T)
# print("Selected alpha:", holdout_bundle["selected_alpha"])


# 
# # 4. SAVE MULTITASK COEFFICIENTS
# 

# # MultiTaskLassoCV coef_ shape is: n_lipids x n_genes.
# coef_df = pd.DataFrame(
#     holdout_bundle["model"].coef_,
#     index=y_train.columns,
#     columns=x_train.columns,
# )
# coef_df.index.name = "Lipid"
# coef_df.to_csv(
#     OUTPUT_DIR / "MultiTaskLassoCV_coefficients.csv"
# )

# nonzero_counts = (coef_df != 0).sum(axis=1)
# nonzero_summary = pd.DataFrame({
#     "Lipid": nonzero_counts.index,
#     "Nonzero_RNA_Features": nonzero_counts.values,
# }).sort_values(
#     "Nonzero_RNA_Features",
#     ascending=False,
# )
# nonzero_summary.to_csv(
#     OUTPUT_DIR / "MultiTaskLassoCV_nonzero_feature_counts.csv",
#     index=False,
# )

# print("\nNonzero RNA features per lipid:")
# print(nonzero_counts.describe())


# 
# # 5. HOLDOUT SCATTERPLOT
# 

# true_flat = y_holdout.to_numpy(dtype=float).ravel()
# pred_flat = holdout_pred.to_numpy(dtype=float).ravel()
# mask = np.isfinite(true_flat) & np.isfinite(pred_flat)
# true_flat = true_flat[mask]
# pred_flat = pred_flat[mask]

# if len(true_flat) > 10000:
#     rng = np.random.default_rng(RANDOM_SEED)
#     selected = rng.choice(len(true_flat), size=10000, replace=False)
#     true_plot = true_flat[selected]
#     pred_plot = pred_flat[selected]
# else:
#     true_plot = true_flat
#     pred_plot = pred_flat

# fig, ax = plt.subplots(figsize=(7, 7))
# ax.scatter(true_plot, pred_plot, alpha=0.35, s=12)

# minimum = min(true_plot.min(), pred_plot.min())
# maximum = max(true_plot.max(), pred_plot.max())
# ax.plot(
#     [minimum, maximum],
#     [minimum, maximum],
#     linestyle="--",
#     linewidth=1,
#     color="black",
# )

# metrics_row = global_metrics.iloc[0]
# ax.set_title(
#     "MultiTaskLassoCV Holdout Performance\n"
#     f"Pearson r={metrics_row['Pearson_r']:.3f} | "
#     f"RMSE={metrics_row['RMSE']:.3f} | "
#     f"MAE={metrics_row['MAE']:.3f}"
# )
# ax.set_xlabel("True lipid abundance")
# ax.set_ylabel("Predicted lipid abundance")
# plt.tight_layout()
# plt.savefig(
#     OUTPUT_DIR / "MultiTaskLassoCV_holdout_scatter.png",
#     dpi=300,
#     bbox_inches="tight",
# )
# plt.show()
# plt.close()


# 
# # 6. REFIT FINAL PRODUCTION MODEL ON ALL USABLE SAMPLES
# 

# final_scaler_x = StandardScaler()
# final_scaler_y = StandardScaler()

# x_all_scaled = final_scaler_x.fit_transform(x_df)
# y_all_scaled = final_scaler_y.fit_transform(y_df)

# final_model = MultiTaskLassoCV(
#     alphas=ALPHAS,
#     cv=CV_FOLDS,
#     max_iter=MAX_ITER,
#     n_jobs=N_JOBS,
#     random_state=RANDOM_SEED,
# )

# final_start = time.perf_counter()
# final_model.fit(x_all_scaled, y_all_scaled)
# final_train_seconds = time.perf_counter() - final_start

# final_bundle = {
#     "model_name": "MultiTaskLassoCV",
#     "model": final_model,
#     "scaler_x": final_scaler_x,
#     "scaler_y": final_scaler_y,
#     "X_columns": list(x_df.columns),
#     "Y_columns": list(y_df.columns),
#     "alpha_grid": list(ALPHAS),
#     "selected_alpha": float(final_model.alpha_),
#     "cv_folds": CV_FOLDS,
#     "max_iter": MAX_ITER,
#     "training_samples": list(x_df.index),
#     "training_metadata": meta.copy(),
#     "target_scale": "Input Y files are log2-scaled lipid abundance tables",
# }

# with open(
#     OUTPUT_DIR / "rna2lipid_MultiTaskLassoCV_final_bundle.pkl",
#     "wb",
# ) as handle:
#     pickle.dump(
#         final_bundle,
#         handle,
#         protocol=pickle.HIGHEST_PROTOCOL,
#     )

# final_coef_df = pd.DataFrame(
#     final_model.coef_,
#     index=y_df.columns,
#     columns=x_df.columns,
# )
# final_coef_df.index.name = "Lipid"
# final_coef_df.to_csv(
#     OUTPUT_DIR / "rna2lipid_MultiTaskLassoCV_final_coefficients.csv"
# )

# print("\nFinal model saved:")
# print(OUTPUT_DIR / "rna2lipid_MultiTaskLassoCV_final_bundle.pkl")
# print("Final selected alpha:", final_model.alpha_)
# print(f"Final full-data fit time: {final_train_seconds:.2f} seconds")




