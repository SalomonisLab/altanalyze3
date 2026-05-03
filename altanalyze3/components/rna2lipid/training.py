from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import pickle
from typing import Any, Dict, Mapping

import numpy as np
import pandas as pd
from sklearn.linear_model import MultiTaskElasticNetCV
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler

from .pipeline import PreparedTrainingData, build_training_dataset, dump_json, json_ready, load_json, utc_timestamp


@dataclass(frozen=True)
class TrainedBundleArtifacts:
    bundle: Dict[str, Any]
    training_summary: Dict[str, Any]
    train_predictions: pd.DataFrame
    per_lipid_metrics: pd.DataFrame


def _safe_float(value) -> float | None:
    if value is None:
        return None
    value = float(value)
    if np.isnan(value) or np.isinf(value):
        return None
    return value


def normalize_target_scaling(value: Any) -> str:
    raw = str("standard" if value is None else value).strip().lower()
    if raw in {"standard", "standardize", "standardized", "zscore", "z-score"}:
        return "standard"
    if raw in {"none", "off", "false", "identity", "raw", "unscaled"}:
        return "none"
    raise ValueError(f"Unsupported target scaling mode: {value!r}")


def fit_target_scaler(y_values: pd.DataFrame, *, target_scaling: str):
    mode = normalize_target_scaling(target_scaling)
    if mode == "standard":
        scaler_y = StandardScaler()
        transformed = scaler_y.fit_transform(y_values)
        return scaler_y, transformed
    return None, y_values.to_numpy(dtype=float)


def inverse_target_scaler(predicted, *, scaler_y, target_scaling: str) -> np.ndarray:
    mode = normalize_target_scaling(target_scaling)
    if mode == "standard":
        if scaler_y is None:
            raise ValueError("Expected scaler_y for standardized target outputs.")
        return scaler_y.inverse_transform(predicted)
    return np.asarray(predicted, dtype=float)


def regression_metrics(y_true: pd.DataFrame, y_pred: pd.DataFrame) -> Dict[str, Any]:
    true_vals = y_true.to_numpy(dtype=float).ravel()
    pred_vals = y_pred.to_numpy(dtype=float).ravel()

    rmse = float(np.sqrt(np.mean((true_vals - pred_vals) ** 2)))
    denom = float(np.mean(np.abs(true_vals)))
    normalized_rmse = rmse / denom if denom else None

    nonzero_mask = true_vals != 0
    mape = None
    if nonzero_mask.any():
        mape = float(np.mean(np.abs((true_vals[nonzero_mask] - pred_vals[nonzero_mask]) / true_vals[nonzero_mask])))

    pearson_r = None
    if true_vals.size > 1 and np.std(true_vals) > 0 and np.std(pred_vals) > 0:
        pearson_r = float(np.corrcoef(true_vals, pred_vals)[0, 1])

    r2_global = None
    if true_vals.size > 1 and np.std(true_vals) > 0:
        r2_global = float(r2_score(true_vals, pred_vals))

    per_lipid_rows = []
    for lipid in y_true.columns:
        truth = y_true[lipid].to_numpy(dtype=float)
        pred = y_pred[lipid].to_numpy(dtype=float)
        lipid_r2 = None
        if truth.size > 1 and np.std(truth) > 0:
            lipid_r2 = float(r2_score(truth, pred))
        per_lipid_rows.append({
            "lipid": lipid,
            "r2": _safe_float(lipid_r2),
            "rmse": float(np.sqrt(np.mean((truth - pred) ** 2))),
            "mean_true": float(np.mean(truth)),
            "mean_pred": float(np.mean(pred)),
        })

    per_lipid_metrics = pd.DataFrame(per_lipid_rows).set_index("lipid")
    valid_r2 = per_lipid_metrics["r2"].dropna()

    return {
        "global": {
            "sample_count": int(y_true.shape[0]),
            "output_count": int(y_true.shape[1]),
            "rmse": rmse,
            "normalized_rmse": _safe_float(normalized_rmse),
            "mape": _safe_float(mape),
            "pearson_r": _safe_float(pearson_r),
            "r2_global": _safe_float(r2_global),
            "mean_output_r2": _safe_float(valid_r2.mean()) if not valid_r2.empty else None,
            "median_output_r2": _safe_float(valid_r2.median()) if not valid_r2.empty else None,
            "positive_output_r2_count": int((valid_r2 > 0).sum()) if not valid_r2.empty else 0,
        },
        "per_lipid_metrics": per_lipid_metrics,
    }


def train_multitask_bundle(
    data: PreparedTrainingData,
    *,
    model_config: Mapping[str, Any],
    training_label: str,
    source_config_path: str,
) -> TrainedBundleArtifacts:
    target_scaling = normalize_target_scaling(model_config.get("target_scaling"))
    scaler_x = StandardScaler()
    X_scaled = scaler_x.fit_transform(data.X)

    scaler_y, y_model_values = fit_target_scaler(data.Y, target_scaling=target_scaling)

    model = MultiTaskElasticNetCV(
        l1_ratio=model_config["l1_ratio"],
        alphas=np.asarray(model_config["alphas"], dtype=float),
        cv=int(model_config["cv"]),
        max_iter=int(model_config["max_iter"]),
        n_jobs=int(model_config["n_jobs"]),
    )
    model.fit(X_scaled, y_model_values)

    predicted_model_values = model.predict(X_scaled)
    predicted = inverse_target_scaler(predicted_model_values, scaler_y=scaler_y, target_scaling=target_scaling)
    pred_df = pd.DataFrame(predicted, index=data.X.index, columns=data.Y.columns)

    metrics = regression_metrics(data.Y, pred_df)
    per_lipid_metrics = metrics.pop("per_lipid_metrics")

    metadata = {
        "bundle_format_version": 2,
        "training_label": training_label,
        "created_at": utc_timestamp(),
        "source_config_path": source_config_path,
        "training_summary": {
            "sample_count": int(data.X.shape[0]),
            "input_gene_count": int(data.X.shape[1]),
            "output_lipid_count": int(data.Y.shape[1]),
            "donor_count": int(data.sample_metadata["donor_id"].nunique()),
            "profile_kind_counts": {str(key): int(value) for key, value in data.sample_metadata["profile_kind"].value_counts().sort_index().items()},
        },
        "dataset_manifest": data.manifest,
        "model_hyperparameters": json_ready(dict(model_config)),
        "target_scaling": {
            "mode": target_scaling,
            "enabled": bool(target_scaling == "standard"),
        },
        "fit_metrics": json_ready(metrics["global"]),
    }

    bundle = {
        "model": model,
        "scaler_x": scaler_x,
        "scaler_y": scaler_y,
        "X_columns": data.X.columns.tolist(),
        "Y_columns": data.Y.columns.tolist(),
        "metadata": metadata,
    }

    training_summary = {
        "metadata": metadata,
        "fit_metrics": json_ready(metrics["global"]),
    }

    return TrainedBundleArtifacts(
        bundle=bundle,
        training_summary=training_summary,
        train_predictions=pred_df,
        per_lipid_metrics=per_lipid_metrics,
    )


def save_bundle_artifacts(
    artifacts: TrainedBundleArtifacts,
    *,
    bundle_path: str | Path,
    output_dir: str | Path,
) -> Dict[str, str]:
    bundle_path = Path(bundle_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    bundle_path.parent.mkdir(parents=True, exist_ok=True)

    with bundle_path.open("wb") as handle:
        pickle.dump(artifacts.bundle, handle)

    summary_path = output_dir / "training_summary.json"
    predictions_path = output_dir / "training_predictions.tsv"
    per_lipid_path = output_dir / "training_per_lipid_metrics.tsv"
    manifest_path = output_dir / "bundle_manifest.json"

    dump_json(summary_path, json_ready(artifacts.training_summary))
    artifacts.train_predictions.to_csv(predictions_path, sep="\t")
    artifacts.per_lipid_metrics.to_csv(per_lipid_path, sep="\t")
    dump_json(manifest_path, json_ready(artifacts.bundle["metadata"]))

    return {
        "bundle_path": str(bundle_path),
        "summary_path": str(summary_path),
        "predictions_path": str(predictions_path),
        "per_lipid_metrics_path": str(per_lipid_path),
        "manifest_path": str(manifest_path),
    }


def run_training(
    *,
    config_path: str | Path,
    bundle_path: str | Path | None = None,
    output_dir: str | Path | None = None,
) -> Dict[str, Any]:
    config = load_json(config_path)
    config_file = Path(config["_config_path"])
    base_dir = config_file.parent
    train_cfg = config["training"]
    data = build_training_dataset(config_file)

    resolved_bundle_path = Path(bundle_path) if bundle_path else (base_dir / train_cfg["bundle_output"]).resolve()
    resolved_output_dir = Path(output_dir) if output_dir else (base_dir / train_cfg["report_output_dir"]).resolve()

    artifacts = train_multitask_bundle(
        data,
        model_config=train_cfg["model"],
        training_label=str(config.get("name", "rna2lipid-training")),
        source_config_path=str(config_file),
    )
    written = save_bundle_artifacts(
        artifacts,
        bundle_path=resolved_bundle_path,
        output_dir=resolved_output_dir,
    )

    return {
        "bundle_metadata": json_ready(artifacts.bundle["metadata"]),
        "written_files": written,
    }
