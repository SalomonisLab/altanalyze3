"""Prove the ported sparse lipid-by-lipid trainer is numerically identical to
the delivered training script.

The delivered script is a top-level program that reads cluster paths at import
time, so it cannot be imported. This module slices the two functions it needs
out of provenance/train_sparse_lipidwise_delivered_2026-08-19.py by line range
and executes only those, which removes any transcription risk.

Both implementations are then fitted on the same real lung matrices and every
selected gene, every coefficient and every prediction is compared.

Run:
  /usr/bin/python3 components/rna2lipid/validation/validate_trainer_equivalence.py \
      [--max-lipids N]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import types
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

MODULE_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = MODULE_DIR.parent.parent.parent
DELIVERED = MODULE_DIR / "provenance" / "train_sparse_lipidwise_delivered_2026-08-19.py"

# 1-indexed, inclusive line ranges inside DELIVERED.
SAFE_PEARSON_LINES = (648, 689)
FIT_FUNCTION_LINES = (838, 1313)

TOP_GENE_OPTIONS = [25, 50, 100, 200]
L1_RATIO_GRID = [0.70, 0.90, 0.95, 1.00]
ALPHA_GRID = np.logspace(-3, 1, 30)
CV_FOLDS = 3
MAX_ITER = 50000
SPARSITY_PENALTY = 0.002
RANDOM_SEED = 1
N_JOBS = -1


def load_delivered_reference() -> types.ModuleType:
    """Execute only the delivered safe_pearson and fit function."""
    lines = DELIVERED.read_text().splitlines()

    def slice_lines(bounds):
        start, stop = bounds
        return "\n".join(lines[start - 1 : stop])

    header = "\n".join([
        "import time",
        "import numpy as np",
        "import pandas as pd",
        "from scipy.stats import pearsonr",
        "from sklearn.preprocessing import StandardScaler",
        "from sklearn.linear_model import ElasticNetCV",
        "from sklearn.metrics import r2_score",
        "def tqdm(iterable, **kwargs):",
        "    return iterable",
        f"N_JOBS = {N_JOBS}",
    ])
    source = "\n\n".join([header, slice_lines(SAFE_PEARSON_LINES), slice_lines(FIT_FUNCTION_LINES)])

    module = types.ModuleType("rna2lipid_delivered_reference")
    exec(compile(source, str(DELIVERED), "exec"), module.__dict__)
    if not hasattr(module, "fit_sparse_lipidwise_elasticnet"):
        raise RuntimeError("The extracted line range does not define the fit function")
    return module


def build_matrices():
    """Rebuild the delivered preprocessing on the tables present in this repo.

    The delivered run used log2 lipid tables that are not in this repository.
    Equivalence does not depend on which lipid table is used, so this test uses
    the two lipid tables the repository does carry.
    """
    x_bulk = pd.read_csv(MODULE_DIR / "data" / "feature_blankreduiction.csv", index_col=0)
    x_cell = pd.read_csv(MODULE_DIR / "data" / "newrna_cell_clair_filtered_symbol.csv", index_col=0).T
    y_bulk = pd.read_csv(MODULE_DIR / "data" / "Bulk_lipids_cleaned_normalized_median_527.csv", index_col=0)
    y_cell = pd.read_csv(MODULE_DIR / "data" / "cell_lipids_cleaned_norm_median_286.csv", index_col=0)

    x_cell = x_cell.T.groupby(level=0).mean().T
    for frame in (x_bulk, x_cell, y_bulk, y_cell):
        frame.index = frame.index.astype(str).str.strip()
        frame.columns = frame.columns.astype(str).str.strip()

    X = pd.concat([x_bulk, x_cell], axis=0, join="inner")
    Y = pd.concat([y_bulk, y_cell], axis=0, join="inner")
    X = X.loc[~X.index.duplicated(keep="first")].copy()
    Y = Y.loc[~Y.index.duplicated(keep="first")].copy()
    common = X.index.intersection(Y.index)
    X = X.loc[common].copy()
    Y = Y.loc[common].copy()
    X = X.apply(pd.to_numeric, errors="coerce")
    Y = Y.apply(pd.to_numeric, errors="coerce")
    X = X.fillna(X.median(numeric_only=True))
    Y = Y.fillna(Y.median(numeric_only=True))
    X = X.loc[:, X.notna().sum(axis=0) > 0]
    X = X.loc[:, (X != 0).any(axis=0)]
    return X, Y


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-lipids", type=int, default=None)
    args = parser.parse_args()

    sys.path.insert(0, str(REPO_ROOT))
    from altanalyze3.components.rna2lipid.training import fit_sparse_lipidwise_elasticnet as ported

    delivered_module = load_delivered_reference()
    delivered = delivered_module.fit_sparse_lipidwise_elasticnet

    X, Y = build_matrices()
    if args.max_lipids is not None:
        Y = Y.iloc[:, : args.max_lipids]

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "delivered_script": str(DELIVERED),
        "delivered_line_ranges": {
            "safe_pearson": list(SAFE_PEARSON_LINES),
            "fit_sparse_lipidwise_elasticnet": list(FIT_FUNCTION_LINES),
        },
        "matrix_shapes": {"X": list(X.shape), "Y": list(Y.shape)},
    }

    print(f"Fitting delivered reference on X{X.shape} Y{Y.shape} ...", flush=True)
    delivered_pred, delivered_bundle, _ = delivered(
        x_train=X,
        y_train=Y,
        x_test=X,
        top_gene_options=TOP_GENE_OPTIONS,
        l1_ratio_grid=L1_RATIO_GRID,
        alpha_grid=ALPHA_GRID,
        cv_folds=CV_FOLDS,
        max_iter=MAX_ITER,
        sparsity_penalty=SPARSITY_PENALTY,
        random_seed=RANDOM_SEED,
    )

    print("Fitting ported altanalyze3 trainer ...", flush=True)
    ported_pred, ported_bundle, _ = ported(
        X,
        Y,
        X,
        top_gene_options=TOP_GENE_OPTIONS,
        l1_ratio_grid=L1_RATIO_GRID,
        alpha_grid=ALPHA_GRID,
        cv_folds=CV_FOLDS,
        max_iter=MAX_ITER,
        sparsity_penalty=SPARSITY_PENALTY,
        random_seed=RANDOM_SEED,
        n_jobs=N_JOBS,
        verbose=False,
    )

    mismatches = []
    coefficient_max_diff = 0.0
    intercept_max_diff = 0.0
    for lipid in Y.columns:
        a = delivered_bundle["models"][lipid]
        b = ported_bundle["models"][lipid]
        if list(a["genes"]) != list(b["genes"]):
            mismatches.append({"lipid": lipid, "field": "genes"})
        if int(a["top_n"]) != int(b["top_n"]):
            mismatches.append({"lipid": lipid, "field": "top_n"})
        if float(a["selected_alpha"]) != float(b["selected_alpha"]):
            mismatches.append({"lipid": lipid, "field": "selected_alpha"})
        if float(a["selected_l1_ratio"]) != float(b["selected_l1_ratio"]):
            mismatches.append({"lipid": lipid, "field": "selected_l1_ratio"})
        coefficient_max_diff = max(
            coefficient_max_diff,
            float(np.max(np.abs(a["model"].coef_ - b["model"].coef_))),
        )
        intercept_max_diff = max(
            intercept_max_diff,
            float(abs(float(a["model"].intercept_) - float(b["model"].intercept_))),
        )

    prediction_max_diff = float(
        np.max(np.abs(delivered_pred.to_numpy(float) - ported_pred.to_numpy(float)))
    )
    scaler_x_identical = bool(
        np.array_equal(delivered_bundle["scaler_x"].mean_, ported_bundle["scaler_x"].mean_)
        and np.array_equal(delivered_bundle["scaler_x"].scale_, ported_bundle["scaler_x"].scale_)
    )
    scaler_y_identical = bool(
        np.array_equal(delivered_bundle["scaler_y"].mean_, ported_bundle["scaler_y"].mean_)
        and np.array_equal(delivered_bundle["scaler_y"].scale_, ported_bundle["scaler_y"].scale_)
    )

    report["check_E_trainer_equivalence"] = {
        "lipids_compared": int(Y.shape[1]),
        "models_compared": int(len(ported_bundle["models"])),
        "coefficients_compared": int(
            sum(len(v["model"].coef_) for v in ported_bundle["models"].values())
        ),
        "predictions_compared": int(delivered_pred.size),
        "selection_field_mismatches": mismatches,
        "max_abs_coefficient_difference": coefficient_max_diff,
        "max_abs_intercept_difference": intercept_max_diff,
        "max_abs_prediction_difference": prediction_max_diff,
        "scaler_x_identical": scaler_x_identical,
        "scaler_y_identical": scaler_y_identical,
        "x_columns_identical": list(delivered_bundle["X_columns"]) == list(ported_bundle["X_columns"]),
        "y_columns_identical": list(delivered_bundle["Y_columns"]) == list(ported_bundle["Y_columns"]),
        "summary_table_equal": bool(
            delivered_bundle["summary"].equals(ported_bundle["summary"])
        ),
        "pass": bool(
            not mismatches
            and coefficient_max_diff == 0.0
            and intercept_max_diff == 0.0
            and prediction_max_diff == 0.0
            and scaler_x_identical
            and scaler_y_identical
        ),
    }
    report["all_pass"] = report["check_E_trainer_equivalence"]["pass"]

    suffix = "" if args.max_lipids is None else f"_first{args.max_lipids}"
    out_path = MODULE_DIR / "validation" / f"trainer_equivalence{suffix}.json"
    out_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print(f"\nWrote {out_path}")
    return 0 if report["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
