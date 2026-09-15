from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import pickle
import time
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import ElasticNetCV, MultiTaskElasticNetCV
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler

from .pipeline import (
    PreparedTrainingData,
    build_dataset,
    build_training_dataset,
    dump_json,
    json_ready,
    load_json,
    utc_timestamp,
)


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

    written = {
        "bundle_path": str(bundle_path),
        "summary_path": str(summary_path),
        "predictions_path": str(predictions_path),
        "per_lipid_metrics_path": str(per_lipid_path),
        "manifest_path": str(manifest_path),
    }

    # Lipid-wise bundles carry three extra tables. Write them beside the report
    # so the selected gene sets and coefficients are inspectable without
    # unpickling the bundle.
    for key, filename in (
        ("summary", "model_summary.tsv"),
        ("coefficients", "nonzero_coefficients.tsv"),
        ("candidate_models", "candidate_top_gene_models.tsv"),
    ):
        table = artifacts.bundle.get(key)
        if isinstance(table, pd.DataFrame):
            table_path = output_dir / filename
            table.to_csv(table_path, sep="\t", index=False)
            written[f"{key}_path"] = str(table_path)

    return written


# ---------------------------------------------------------------------------
# Sparse lipid-by-lipid ElasticNetCV
#
# One ElasticNetCV per lipid. For each lipid the genes are ranked by their
# training-set Pearson correlation with that lipid, several top-N gene sets are
# fitted, and the set with the best sparsity-adjusted training R2 is kept.
#
# This is a port of the delivered training script, kept at
# provenance/train_sparse_lipidwise_delivered_2026-08-19.py. The numerics are
# unchanged: validation/validate_trainer_equivalence.py fits both on the same
# data and compares every coefficient.
# ---------------------------------------------------------------------------


def resolve_alpha_grid(value: Any) -> np.ndarray:
    """Accept an explicit list of alphas or a {'logspace': [start, stop, num]} spec."""
    if isinstance(value, Mapping):
        if "logspace" in value:
            start, stop, num = value["logspace"]
            return np.logspace(float(start), float(stop), int(num))
        if "linspace" in value:
            start, stop, num = value["linspace"]
            return np.linspace(float(start), float(stop), int(num))
        raise ValueError(f"Unsupported alpha grid specification: {dict(value)!r}")
    return np.asarray(list(value), dtype=float)


def fit_sparse_lipidwise_elasticnet(
    x_train: pd.DataFrame,
    y_train: pd.DataFrame,
    x_test: pd.DataFrame,
    *,
    top_gene_options: Sequence[int],
    l1_ratio_grid: Sequence[float],
    alpha_grid: Sequence[float],
    cv_folds: int,
    max_iter: int,
    sparsity_penalty: float,
    random_seed: int,
    n_jobs: int = -1,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, Any], float]:
    """Fit one ElasticNetCV per lipid and return (predictions, bundle, seconds)."""
    scaler_x = StandardScaler()
    x_train_scaled = pd.DataFrame(
        scaler_x.fit_transform(x_train), index=x_train.index, columns=x_train.columns
    )
    x_test_scaled = pd.DataFrame(
        scaler_x.transform(x_test), index=x_test.index, columns=x_test.columns
    )

    scaler_y = StandardScaler()
    y_train_scaled = pd.DataFrame(
        scaler_y.fit_transform(y_train), index=y_train.index, columns=y_train.columns
    )

    all_models: Dict[str, Any] = {}
    predictions_scaled = pd.DataFrame(index=x_test.index)
    coefficient_records: List[Dict[str, Any]] = []
    summary_records: List[Dict[str, Any]] = []
    candidate_records: List[Dict[str, Any]] = []

    start_time = time.perf_counter()
    total_lipids = int(y_train_scaled.shape[1])

    for position, lipid in enumerate(y_train_scaled.columns, start=1):
        if verbose and (position == 1 or position % 25 == 0 or position == total_lipids):
            print(f"  lipid {position}/{total_lipids}: {lipid}", flush=True)

        y = y_train_scaled[lipid]

        correlations = (
            x_train_scaled.apply(lambda feature: feature.corr(y), axis=0)
            .replace([np.inf, -np.inf], np.nan)
            .fillna(0)
        )
        ranked_genes = correlations.abs().sort_values(ascending=False)

        best_model = None
        best_genes = None
        best_test_prediction = None
        best_training_prediction = None
        best_score = -np.inf
        best_top_n = None

        for top_n in top_gene_options:
            actual_top_n = min(int(top_n), x_train_scaled.shape[1])
            selected_genes = ranked_genes.head(actual_top_n).index.tolist()
            x_train_subset = x_train_scaled[selected_genes]
            x_test_subset = x_test_scaled[selected_genes]

            model = ElasticNetCV(
                l1_ratio=list(l1_ratio_grid),
                alphas=np.asarray(alpha_grid, dtype=float),
                cv=int(cv_folds),
                max_iter=int(max_iter),
                n_jobs=int(n_jobs),
                random_state=int(random_seed),
                selection="cyclic",
            )
            model.fit(x_train_subset, y)

            training_prediction = model.predict(x_train_subset)
            test_prediction = model.predict(x_test_subset)
            training_r2 = r2_score(y, training_prediction)
            number_nonzero = int(np.sum(model.coef_ != 0))
            sparsity_adjusted_score = training_r2 - sparsity_penalty * number_nonzero

            candidate_records.append({
                "Lipid": lipid,
                "Top_N": actual_top_n,
                "Alpha": float(model.alpha_),
                "L1_ratio": float(model.l1_ratio_),
                "Nonzero_Coefficients": number_nonzero,
                "Train_R2_scaled": training_r2,
                "Sparsity_Adjusted_Score": sparsity_adjusted_score,
            })

            if sparsity_adjusted_score > best_score:
                best_score = sparsity_adjusted_score
                best_model = model
                best_genes = selected_genes
                best_training_prediction = training_prediction
                best_test_prediction = test_prediction
                best_top_n = actual_top_n

        if best_model is None:
            raise RuntimeError(f"No candidate model was fitted for lipid {lipid!r}")

        all_models[lipid] = {
            "model": best_model,
            "genes": best_genes,
            "top_n": best_top_n,
            "selected_alpha": float(best_model.alpha_),
            "selected_l1_ratio": float(best_model.l1_ratio_),
            "sparsity_adjusted_score": float(best_score),
        }

        predictions_scaled[lipid] = best_test_prediction

        nonzero_mask = best_model.coef_ != 0
        nonzero_genes = np.asarray(best_genes)[nonzero_mask]
        nonzero_coefficients = best_model.coef_[nonzero_mask]
        for gene, coefficient in zip(nonzero_genes, nonzero_coefficients):
            coefficient_records.append({
                "Lipid": lipid,
                "Gene": gene,
                "Coefficient": coefficient,
                "Abs_Coefficient": abs(coefficient),
                "Correlation_with_lipid": correlations.loc[gene],
            })

        summary_records.append({
            "Lipid": lipid,
            "Top_N_Correlation_Filter": best_top_n,
            "Candidate_Top_N_Options": ",".join(map(str, top_gene_options)),
            "Nonzero_Coefficients": int(len(nonzero_genes)),
            "Alpha": float(best_model.alpha_),
            "L1_ratio": float(best_model.l1_ratio_),
            "Train_R2_scaled": r2_score(y, best_training_prediction),
            "Train_Pearson_scaled": _safe_pearson(y, best_training_prediction),
            "Sparsity_Adjusted_Score": best_score,
        })

    training_seconds = time.perf_counter() - start_time

    predictions_scaled = predictions_scaled[y_train.columns]
    predictions = scaler_y.inverse_transform(predictions_scaled)
    prediction_df = pd.DataFrame(predictions, index=x_test.index, columns=y_train.columns)

    bundle = {
        "model_name": "SparseLipidwiseElasticNetCV",
        "architecture": "Separate ElasticNetCV model for each lipid",
        "models": all_models,
        "scaler_x": scaler_x,
        "scaler_y": scaler_y,
        "X_columns": list(x_train.columns),
        "Y_columns": list(y_train.columns),
        "top_gene_options": list(top_gene_options),
        "l1_ratio_grid": list(l1_ratio_grid),
        "alpha_grid": list(np.asarray(alpha_grid, dtype=float)),
        "cv_folds": int(cv_folds),
        "max_iter": int(max_iter),
        "sparsity_penalty": float(sparsity_penalty),
        "random_seed": int(random_seed),
        "summary": pd.DataFrame(summary_records),
        "coefficients": pd.DataFrame(coefficient_records),
        "candidate_models": pd.DataFrame(candidate_records),
    }
    return prediction_df, bundle, training_seconds


def _safe_pearson(observed, predicted) -> float | None:
    observed = np.asarray(observed, dtype=float).ravel()
    predicted = np.asarray(predicted, dtype=float).ravel()
    if observed.size < 2:
        return None
    if np.std(observed) == 0 or np.std(predicted) == 0:
        return None
    return float(pearsonr(observed, predicted)[0])


# ---------------------------------------------------------------------------
# Holdout construction
# ---------------------------------------------------------------------------


def generate_holdout_constrained(
    metadata: pd.DataFrame,
    *,
    seed: int = 42,
    total_holdout: int = 28,
    min_per_group: int = 5,
    max_per_donor: int = 3,
    required_groups: Sequence[str],
    max_tries: int = 10000,
) -> Tuple[List[str], pd.DataFrame, pd.DataFrame]:
    """The delivered holdout rule: every group is represented, no donor over a cap.

    NOT donor-disjoint. A donor with more profiles than ``max_per_donor``
    appears in both the training and the holdout set. Use holdout mode
    ``donor_disjoint`` when a donor-held-out estimate is required.
    """
    rng = np.random.default_rng(seed)
    sample_ids = metadata.index.tolist()
    required_groups = [str(value) for value in required_groups]

    for _ in range(int(max_tries)):
        donor_counts: Dict[str, int] = {}
        holdout_samples: List[str] = []
        valid_split = True

        for group in required_groups:
            group_pool = metadata.index[metadata["group"] == group].tolist()
            rng.shuffle(group_pool)
            selected: List[str] = []
            for sample_id in group_pool:
                donor = metadata.at[sample_id, "donor_id"]
                if donor_counts.get(donor, 0) < max_per_donor:
                    selected.append(sample_id)
                    donor_counts[donor] = donor_counts.get(donor, 0) + 1
                if len(selected) == min_per_group:
                    break
            if len(selected) < min_per_group:
                valid_split = False
                break
            holdout_samples.extend(selected)

        if not valid_split:
            continue

        remaining_needed = total_holdout - len(holdout_samples)
        remaining_pool = [s for s in sample_ids if s not in holdout_samples]
        rng.shuffle(remaining_pool)
        for sample_id in remaining_pool:
            if remaining_needed == 0:
                break
            donor = metadata.at[sample_id, "donor_id"]
            if donor_counts.get(donor, 0) < max_per_donor:
                holdout_samples.append(sample_id)
                donor_counts[donor] = donor_counts.get(donor, 0) + 1
                remaining_needed -= 1

        if len(holdout_samples) != total_holdout:
            continue

        holdout_metadata = metadata.loc[holdout_samples].copy()
        if holdout_metadata["donor_id"].value_counts().max() > max_per_donor:
            continue
        counts = holdout_metadata["group"].value_counts()
        if not all(counts.get(group, 0) >= min_per_group for group in required_groups):
            continue

        train_metadata = metadata.drop(index=holdout_samples).copy()
        return holdout_samples, train_metadata, holdout_metadata

    raise RuntimeError(
        "Could not generate a valid constrained holdout split. Relax "
        "total_holdout, min_per_group or max_per_donor for this dataset."
    )


def generate_holdout_donor_disjoint(
    metadata: pd.DataFrame,
    *,
    seed: int = 42,
    holdout_donor_fraction: float = 0.3,
    min_holdout_donors: int = 2,
) -> Tuple[List[str], pd.DataFrame, pd.DataFrame]:
    """ADDITION, not the delivered rule. Hold out whole donors.

    Every profile of a held-out donor leaves the training set, so the resulting
    metrics carry no donor leakage.
    """
    rng = np.random.default_rng(seed)
    donors = sorted(metadata["donor_id"].unique().tolist())
    if len(donors) < 2:
        raise ValueError("Donor-disjoint holdout needs at least 2 donors")
    n_holdout = max(int(min_holdout_donors), int(round(len(donors) * float(holdout_donor_fraction))))
    n_holdout = min(n_holdout, len(donors) - 1)
    order = list(donors)
    rng.shuffle(order)
    holdout_donors = set(order[:n_holdout])

    holdout_metadata = metadata.loc[metadata["donor_id"].isin(holdout_donors)].copy()
    train_metadata = metadata.loc[~metadata["donor_id"].isin(holdout_donors)].copy()
    return holdout_metadata.index.tolist(), train_metadata, holdout_metadata


def build_holdout(metadata: pd.DataFrame, config: Mapping[str, Any]) -> Tuple[List[str], pd.DataFrame, pd.DataFrame]:
    mode = str(config.get("mode", "constrained")).strip().lower()
    if mode == "constrained":
        required = config.get("required_groups")
        if not required:
            required = sorted(metadata["group"].unique().tolist())
        return generate_holdout_constrained(
            metadata,
            seed=int(config.get("seed", 42)),
            total_holdout=int(config.get("total_holdout", 28)),
            min_per_group=int(config.get("min_per_group", 5)),
            max_per_donor=int(config.get("max_per_donor", 3)),
            required_groups=required,
            max_tries=int(config.get("max_tries", 10000)),
        )
    if mode == "donor_disjoint":
        return generate_holdout_donor_disjoint(
            metadata,
            seed=int(config.get("seed", 42)),
            holdout_donor_fraction=float(config.get("holdout_donor_fraction", 0.3)),
            min_holdout_donors=int(config.get("min_holdout_donors", 2)),
        )
    if mode == "none":
        return [], metadata.copy(), metadata.iloc[0:0].copy()
    raise ValueError(f"Unsupported holdout mode: {mode!r}")


def holdout_global_metrics(y_true: pd.DataFrame, y_pred: pd.DataFrame, *, model_name: str) -> Dict[str, Any]:
    true_vals = y_true.to_numpy(dtype=float).ravel()
    pred_vals = y_pred.to_numpy(dtype=float).ravel()
    finite = np.isfinite(true_vals) & np.isfinite(pred_vals)
    true_vals = true_vals[finite]
    pred_vals = pred_vals[finite]
    metrics: Dict[str, Any] = {
        "Model": model_name,
        "N_samples": int(y_true.shape[0]),
        "N_lipids": int(y_true.shape[1]),
        "N_values": int(true_vals.size),
        "N_values_dropped_non_finite": int(int(y_true.size) - int(true_vals.size)),
        "RMSE": float(np.sqrt(np.mean((true_vals - pred_vals) ** 2))),
        "MAE": float(np.mean(np.abs(true_vals - pred_vals))),
        "True_SD": float(np.std(true_vals)),
        "Predicted_SD": float(np.std(pred_vals)),
    }
    metrics["SD_Ratio"] = (
        metrics["Predicted_SD"] / metrics["True_SD"] if metrics["True_SD"] else None
    )
    metrics["Pearson_r"] = _safe_pearson(true_vals, pred_vals)
    if true_vals.size > 1 and np.std(true_vals) > 0 and np.std(pred_vals) > 0:
        metrics["Spearman_r"] = float(spearmanr(true_vals, pred_vals)[0])
        metrics["R2"] = float(r2_score(true_vals, pred_vals))
    else:
        metrics["Spearman_r"] = None
        metrics["R2"] = None
    return metrics


def train_lipidwise_bundle(
    data: PreparedTrainingData,
    *,
    model_config: Mapping[str, Any],
    holdout_config: Mapping[str, Any],
    training_label: str,
    tissue: str | None,
    source_config_path: str,
) -> TrainedBundleArtifacts:
    """Train the sparse lipid-by-lipid model, evaluate on a holdout, then refit on all samples."""
    top_gene_options = [int(v) for v in model_config.get("top_gene_options", [25, 50, 100, 200])]
    l1_ratio_grid = [float(v) for v in model_config.get("l1_ratio", [0.70, 0.90, 0.95, 1.00])]
    alpha_grid = resolve_alpha_grid(model_config.get("alphas", {"logspace": [-3, 1, 30]}))
    cv_folds = int(model_config.get("cv", 3))
    max_iter = int(model_config.get("max_iter", 50000))
    n_jobs = int(model_config.get("n_jobs", -1))
    sparsity_penalty = float(model_config.get("sparsity_penalty", 0.002))
    random_seed = int(model_config.get("random_seed", 1))

    fit_kwargs = dict(
        top_gene_options=top_gene_options,
        l1_ratio_grid=l1_ratio_grid,
        alpha_grid=alpha_grid,
        cv_folds=cv_folds,
        max_iter=max_iter,
        sparsity_penalty=sparsity_penalty,
        random_seed=random_seed,
        n_jobs=n_jobs,
    )

    holdout_metrics: Dict[str, Any] | None = None
    holdout_report: Dict[str, Any] = {"mode": str(holdout_config.get("mode", "constrained"))}

    holdout_samples, meta_train, meta_holdout = build_holdout(data.sample_metadata, holdout_config)
    if holdout_samples:
        x_train = data.X.loc[meta_train.index]
        y_train = data.Y.loc[meta_train.index]
        x_holdout = data.X.loc[meta_holdout.index]
        y_holdout = data.Y.loc[meta_holdout.index]
        print(
            f"Holdout evaluation: train {x_train.shape[0]} profiles, "
            f"holdout {x_holdout.shape[0]} profiles.",
            flush=True,
        )
        holdout_predictions, _, holdout_seconds = fit_sparse_lipidwise_elasticnet(
            x_train, y_train, x_holdout, **fit_kwargs
        )
        holdout_metrics = holdout_global_metrics(
            y_holdout, holdout_predictions, model_name="SparseLipidwiseElasticNetCV"
        )
        holdout_metrics["Training_seconds"] = holdout_seconds
        shared_donors = sorted(
            set(meta_train["donor_id"]).intersection(set(meta_holdout["donor_id"]))
        )
        holdout_report.update({
            "train_profiles": int(x_train.shape[0]),
            "holdout_profiles": int(x_holdout.shape[0]),
            "train_donors": int(meta_train["donor_id"].nunique()),
            "holdout_donors": int(meta_holdout["donor_id"].nunique()),
            "donors_in_both_sets": shared_donors,
            "donor_disjoint": bool(not shared_donors),
            "holdout_group_counts": {
                str(k): int(v) for k, v in meta_holdout["group"].value_counts().sort_index().items()
            },
            "holdout_samples": list(holdout_samples),
        })
    else:
        holdout_report["note"] = "No holdout was constructed; no held-out metric exists."

    print(f"Final model: refitting on all {data.X.shape[0]} profiles.", flush=True)
    final_predictions, bundle, final_seconds = fit_sparse_lipidwise_elasticnet(
        data.X, data.Y, data.X, **fit_kwargs
    )

    metrics = regression_metrics(data.Y, final_predictions)
    per_lipid_metrics = metrics.pop("per_lipid_metrics")

    zero_coefficient_lipids = [
        str(lipid)
        for lipid, entry in bundle["models"].items()
        if int(np.count_nonzero(entry["model"].coef_)) == 0
    ]

    metadata = {
        "bundle_format_version": 3,
        "training_label": training_label,
        "tissue": tissue,
        "created_at": utc_timestamp(),
        "source_config_path": source_config_path,
        "model_architecture": "sparse_lipidwise_elasticnet_cv",
        "training_summary": {
            "sample_count": int(data.X.shape[0]),
            "input_gene_count": int(data.X.shape[1]),
            "output_lipid_count": int(data.Y.shape[1]),
            "donor_count": int(data.sample_metadata["donor_id"].nunique()),
            "group_counts": {
                str(k): int(v)
                for k, v in data.sample_metadata["group"].value_counts().sort_index().items()
            },
            "profile_kind_counts": {
                str(k): int(v)
                for k, v in data.sample_metadata["profile_kind"].value_counts().sort_index().items()
            },
            "final_training_seconds": final_seconds,
            "lipids_with_all_zero_coefficients": len(zero_coefficient_lipids),
            "lipids_with_all_zero_coefficients_names": zero_coefficient_lipids,
        },
        "dataset_manifest": data.manifest,
        "model_hyperparameters": json_ready({
            "top_gene_options": top_gene_options,
            "l1_ratio": l1_ratio_grid,
            "alphas": alpha_grid.tolist(),
            "cv": cv_folds,
            "max_iter": max_iter,
            "n_jobs": n_jobs,
            "sparsity_penalty": sparsity_penalty,
            "random_seed": random_seed,
        }),
        "target_scaling": {"mode": "standard", "enabled": True},
        "holdout": json_ready(holdout_report),
        "holdout_metrics": json_ready(holdout_metrics),
        "resubstitution_metrics": json_ready(metrics["global"]),
    }

    bundle["metadata"] = metadata
    bundle["training_samples"] = list(data.X.index)
    bundle["training_metadata"] = data.sample_metadata.copy()
    bundle["final_training_seconds"] = final_seconds
    bundle["validation_global_metrics"] = holdout_metrics

    training_summary = {
        "metadata": metadata,
        "holdout_metrics": json_ready(holdout_metrics),
        "resubstitution_metrics": json_ready(metrics["global"]),
    }

    return TrainedBundleArtifacts(
        bundle=bundle,
        training_summary=training_summary,
        train_predictions=final_predictions,
        per_lipid_metrics=per_lipid_metrics,
    )


ARCHITECTURES = {"multitask", "sparse_lipidwise"}


def normalize_architecture(value: Any) -> str:
    raw = str("multitask" if value is None else value).strip().lower()
    if raw in {"multitask", "multi_task", "multitask_elasticnet", "multitaskelasticnetcv"}:
        return "multitask"
    if raw in {"sparse_lipidwise", "lipidwise", "sparse_lipid_wise", "per_lipid"}:
        return "sparse_lipidwise"
    raise ValueError(
        f"Unsupported architecture: {value!r}. Choose one of {sorted(ARCHITECTURES)}."
    )


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
    architecture = normalize_architecture(train_cfg.get("architecture"))
    data = build_dataset(config_file)

    resolved_bundle_path = Path(bundle_path) if bundle_path else (base_dir / train_cfg["bundle_output"]).resolve()
    resolved_output_dir = Path(output_dir) if output_dir else (base_dir / train_cfg["report_output_dir"]).resolve()

    print(
        f"Architecture: {architecture}. "
        f"X {data.X.shape[0]} x {data.X.shape[1]}, "
        f"Y {data.Y.shape[0]} x {data.Y.shape[1]}, "
        f"{data.sample_metadata['donor_id'].nunique()} donors.",
        flush=True,
    )

    if architecture == "sparse_lipidwise":
        artifacts = train_lipidwise_bundle(
            data,
            model_config=train_cfg["model"],
            holdout_config=train_cfg.get("holdout", {}),
            training_label=str(config.get("name", "rna2lipid-training")),
            tissue=config.get("tissue"),
            source_config_path=str(config_file),
        )
    else:
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
        "architecture": architecture,
        "bundle_metadata": json_ready(artifacts.bundle["metadata"]),
        "written_files": written,
    }
