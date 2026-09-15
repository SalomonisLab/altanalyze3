"""Validate the lipid-wise rna2lipid api against the prior api and against the
per-lipid coefficients stored in the bundle.

Check A  backward compatibility. The new api must reproduce the prior api on the
         prior MultiTaskElasticNetCV bundle, value for value.
Check B  prediction math. The new api on the lipid-wise bundle must equal an
         independent reimplementation built from each lipid's stored coef_,
         intercept_ and gene list, plus the scaler_y inverse transform.
Check C  h5ad path. Per-cell prediction through predict_from_adata must equal
         the dataframe path on the same values.

Run:
  /usr/bin/python3 components/rna2lipid/validation/validate_api_lipidwise.py
"""
from __future__ import annotations

import importlib.util
import json
import pickle
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

MODULE_DIR = Path(__file__).resolve().parent.parent
VALIDATION_DIR = MODULE_DIR / "validation"
NEW_BUNDLE = MODULE_DIR / "rna2lipid_hs_lung_lipidwise_bundle.pkl"
OLD_BUNDLE = MODULE_DIR / "newnormelastic_multitask_try.pkl"
RNA_CSV = MODULE_DIR / "data" / "feature_blankreduiction.csv"
CELLTYPE_RNA_CSV = MODULE_DIR / "data" / "newrna_cell_clair_filtered_symbol.csv"
PREVIOUS_API = VALIDATION_DIR / "_api_previous_f57f8b3.py"


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _max_abs_diff(a: pd.DataFrame, b: pd.DataFrame) -> float:
    if list(a.columns) != list(b.columns):
        raise AssertionError("column order differs")
    if list(a.index) != list(b.index):
        raise AssertionError("row order differs")
    return float(np.max(np.abs(a.to_numpy(float) - b.to_numpy(float))))


def main() -> int:
    report: dict = {"generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds")}

    bulk_rna = pd.read_csv(RNA_CSV, index_col=0)
    celltype_rna = pd.read_csv(CELLTYPE_RNA_CSV, index_col=0).T
    celltype_rna = celltype_rna.T.groupby(level=0).mean().T
    expression = pd.concat([bulk_rna, celltype_rna], axis=0, join="inner")
    expression = expression.loc[~pd.Index(expression.index).duplicated(keep="first")]
    report["input"] = {
        "bulk_rna_path": str(RNA_CSV),
        "bulk_rna_rows": int(bulk_rna.shape[0]),
        "celltype_rna_path": str(CELLTYPE_RNA_CSV),
        "celltype_rna_rows": int(celltype_rna.shape[0]),
        "rows": int(expression.shape[0]),
        "columns": int(expression.shape[1]),
    }

    new_api = _load_module(MODULE_DIR / "api.py", "rna2lipid_api_new")
    old_api = _load_module(PREVIOUS_API, "rna2lipid_api_previous")

    # ---------------- Check A: backward compatibility on the prior bundle
    old_on_old = old_api.load_bundle(OLD_BUNDLE).predict_from_dataframe(expression)
    new_on_old = new_api.load_bundle(OLD_BUNDLE).predict_from_dataframe(expression)
    diff_a = _max_abs_diff(old_on_old.predictions, new_on_old.predictions)
    report["check_A_backward_compatibility"] = {
        "bundle": str(OLD_BUNDLE),
        "rows_compared": int(old_on_old.predictions.shape[0]),
        "lipids_compared": int(old_on_old.predictions.shape[1]),
        "values_compared": int(old_on_old.predictions.size),
        "max_abs_difference": diff_a,
        "bitwise_identical": bool(
            np.array_equal(
                old_on_old.predictions.to_numpy(float),
                new_on_old.predictions.to_numpy(float),
            )
        ),
        "target_scaling_mode": new_api.load_bundle(OLD_BUNDLE).target_scaling_mode,
        "pass": bool(diff_a == 0.0),
    }

    # ---------------- Check B: prediction math on the lipid-wise bundle
    bundle_obj = new_api.load_bundle(NEW_BUNDLE)
    result = bundle_obj.predict_from_dataframe(expression)

    with NEW_BUNDLE.open("rb") as handle:
        raw = pickle.load(handle)

    genes = [str(g).strip() for g in raw["X_columns"]]
    lipids = [str(l).strip() for l in raw["Y_columns"]]
    aligned = expression.copy()
    aligned.index = [str(i).strip() for i in aligned.index]
    aligned.columns = [str(c).strip() for c in aligned.columns]
    aligned = aligned.T.groupby(level=0).mean().T
    aligned = aligned.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    aligned = aligned.reindex(columns=genes, fill_value=0.0)

    scaled = pd.DataFrame(
        raw["scaler_x"].transform(aligned), index=aligned.index, columns=genes
    )

    manual_scaled = np.zeros((scaled.shape[0], len(lipids)), dtype=float)
    nonzero_counts = []
    for j, lipid in enumerate(lipids):
        entry = raw["models"][lipid]
        estimator = entry["model"]
        selected = [str(g).strip() for g in entry["genes"]]
        design = scaled.loc[:, selected].to_numpy(float)
        manual_scaled[:, j] = design @ estimator.coef_ + float(estimator.intercept_)
        nonzero_counts.append(int(np.count_nonzero(estimator.coef_)))

    manual = pd.DataFrame(
        raw["scaler_y"].inverse_transform(manual_scaled),
        index=aligned.index,
        columns=lipids,
    )
    diff_b = _max_abs_diff(result.predictions, manual)
    report["check_B_prediction_math"] = {
        "bundle": str(NEW_BUNDLE),
        "architecture": bundle_obj.architecture,
        "target_scaling_mode": bundle_obj.target_scaling_mode,
        "rows_compared": int(result.predictions.shape[0]),
        "lipids_compared": int(result.predictions.shape[1]),
        "values_compared": int(result.predictions.size),
        "max_abs_difference": diff_b,
        "per_lipid_nonzero_genes_min": int(min(nonzero_counts)),
        "per_lipid_nonzero_genes_median": int(np.median(nonzero_counts)),
        "per_lipid_nonzero_genes_max": int(max(nonzero_counts)),
        "predicted_min": float(result.predictions.to_numpy(float).min()),
        "predicted_median": float(np.median(result.predictions.to_numpy(float))),
        "predicted_max": float(result.predictions.to_numpy(float).max()),
        "pass": bool(diff_b < 1e-9),
    }

    # ---------------- Check C: adata path equals dataframe path
    import anndata as ad

    adata = ad.AnnData(
        X=aligned.to_numpy(np.float32),
        obs=pd.DataFrame(index=pd.Index(aligned.index, dtype=str)),
        var=pd.DataFrame(index=pd.Index(genes, dtype=str)),
    )
    adata_result = bundle_obj.predict_from_adata(adata, chunk_size=7)
    diff_c = float(
        np.max(
            np.abs(
                adata_result.predictions.loc[result.predictions.index].to_numpy(float)
                - result.predictions.to_numpy(float)
            )
        )
    )
    report["check_C_adata_path"] = {
        "chunk_size": 7,
        "rows_compared": int(adata_result.predictions.shape[0]),
        "values_compared": int(adata_result.predictions.size),
        "max_abs_difference": diff_c,
        "pass": bool(diff_c < 1e-4),
        "note": "adata path casts the design matrix to float32, so the tolerance is 1e-4, not 0.",
    }

    # ---------------- Check D: scope of the shipped lipid-wise bundle
    training_samples = [str(s) for s in raw["training_samples"]]
    training_metadata = raw["training_metadata"]
    zero_model_lipids = [
        lipid
        for lipid in lipids
        if int(np.count_nonzero(raw["models"][lipid]["model"].coef_)) == 0
    ]
    y_mean = np.asarray(raw["scaler_y"].mean_, dtype=float)
    y_sd = np.asarray(raw["scaler_y"].scale_, dtype=float)
    lo = y_mean - 3.0 * y_sd
    hi = y_mean + 3.0 * y_sd
    values = result.predictions.to_numpy(float)
    outside = (values < lo[None, :]) | (values > hi[None, :])
    bulk_mask = np.array([("_" not in idx) for idx in result.predictions.index])
    report["check_D_bundle_scope"] = {
        "training_profiles": len(training_samples),
        "training_donors": int(training_metadata["Donor"].nunique()),
        "training_cell_types": sorted(training_metadata["CellType"].unique().tolist()),
        "training_bulk_profiles": sum(1 for s in training_samples if "_" not in s),
        "holdout_reported_metrics": {
            k: v for k, v in raw["validation_global_metrics"].items()
        },
        "lipids_with_all_zero_coefficients": len(zero_model_lipids),
        "lipids_with_all_zero_coefficients_names": zero_model_lipids,
        "predictions_outside_train_mean_pm_3sd": {
            "all_rows": float(outside.mean()),
            "bulk_rows": float(outside[bulk_mask].mean()) if bulk_mask.any() else None,
            "celltype_rows": float(outside[~bulk_mask].mean()) if (~bulk_mask).any() else None,
        },
        "pass": True,
        "note": (
            "Descriptive, not a gate. The shipped bundle was fit on sorted "
            "cell-type profiles only, so bulk RNA input is extrapolation."
        ),
    }

    report["all_pass"] = all(
        report[key]["pass"] for key in report if key.startswith("check_")
    )

    out_path = VALIDATION_DIR / "api_lipidwise_validation.json"
    out_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print(f"\nWrote {out_path}")
    return 0 if report["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
