from __future__ import annotations

from itertools import combinations
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import MultiTaskElasticNetCV, RidgeCV
from sklearn.metrics import silhouette_score
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.multioutput import MultiOutputRegressor
from sklearn.preprocessing import StandardScaler

from .api import load_bundle
from .pipeline import build_training_dataset, dump_json, json_ready, load_json, resolve_path
from .training import fit_target_scaler, inverse_target_scaler, normalize_target_scaling, regression_metrics


class MeanLipidBaseline:
    def __init__(self) -> None:
        self._mean_vector: np.ndarray | None = None
        self.columns: List[str] = []

    def fit(self, y_train: pd.DataFrame) -> "MeanLipidBaseline":
        self._mean_vector = y_train.mean(axis=0).to_numpy(dtype=float)
        self.columns = list(y_train.columns)
        return self

    def predict(self, n_rows: int, index: Iterable[str]) -> pd.DataFrame:
        if self._mean_vector is None:
            raise RuntimeError("MeanLipidBaseline.fit must be called first.")
        preds = np.repeat(self._mean_vector[None, :], n_rows, axis=0)
        return pd.DataFrame(preds, index=list(index), columns=self.columns)


class PcaRidgeBaseline:
    def __init__(self, *, max_components: int, ridge_alphas: Iterable[float]) -> None:
        self.max_components = int(max_components)
        self.ridge_alphas = np.asarray(list(ridge_alphas), dtype=float)
        self.scaler_x = StandardScaler()
        self.pca: PCA | None = None
        self.model = MultiOutputRegressor(RidgeCV(alphas=self.ridge_alphas))
        self.columns: List[str] = []

    def fit(self, x_train: pd.DataFrame, y_train: pd.DataFrame) -> "PcaRidgeBaseline":
        self.columns = list(y_train.columns)
        x_scaled = self.scaler_x.fit_transform(x_train)
        n_components = min(self.max_components, x_train.shape[0] - 1, x_train.shape[1])
        n_components = max(1, n_components)
        self.pca = PCA(n_components=n_components)
        x_pca = self.pca.fit_transform(x_scaled)
        self.model.fit(x_pca, y_train.to_numpy(dtype=float))
        return self

    def predict(self, x_test: pd.DataFrame) -> pd.DataFrame:
        if self.pca is None:
            raise RuntimeError("PcaRidgeBaseline.fit must be called first.")
        x_scaled = self.scaler_x.transform(x_test)
        x_pca = self.pca.transform(x_scaled)
        preds = self.model.predict(x_pca)
        return pd.DataFrame(preds, index=x_test.index, columns=self.columns)


def _fit_main_model(x_train: pd.DataFrame, y_train: pd.DataFrame, model_config: Mapping[str, Any]):
    target_scaling = normalize_target_scaling(model_config.get("target_scaling"))
    scaler_x = StandardScaler()
    x_scaled = scaler_x.fit_transform(x_train)
    scaler_y, y_model_values = fit_target_scaler(y_train, target_scaling=target_scaling)
    model = MultiTaskElasticNetCV(
        l1_ratio=model_config["l1_ratio"],
        alphas=np.asarray(model_config["alphas"], dtype=float),
        cv=int(model_config["cv"]),
        max_iter=int(model_config["max_iter"]),
        n_jobs=int(model_config["n_jobs"]),
    )
    model.fit(x_scaled, y_model_values)
    return model, scaler_x, scaler_y, target_scaling


def _predict_main_model(model, scaler_x, scaler_y, target_scaling: str, x_test: pd.DataFrame, y_columns: List[str]) -> pd.DataFrame:
    x_scaled = scaler_x.transform(x_test)
    y_model_values = model.predict(x_scaled)
    y_pred = inverse_target_scaler(y_model_values, scaler_y=scaler_y, target_scaling=target_scaling)
    return pd.DataFrame(y_pred, index=x_test.index, columns=y_columns)


def evaluate_internal_holdout(
    *,
    config_path: str | Path,
    output_dir: str | Path,
    max_folds: int | None = None,
) -> Dict[str, Any]:
    config = load_json(config_path)
    data = build_training_dataset(config_path)
    model_config = config["training"]["model"]
    baseline_config = config["evaluation"]["baselines"]
    groups = data.sample_metadata["donor_id"].to_numpy()
    logo = LeaveOneGroupOut()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fold_records: List[Dict[str, Any]] = []
    pred_frames: Dict[str, List[pd.DataFrame]] = {
        "multitask_elastic_net": [],
        "mean_lipid_profile": [],
        "pca_ridge": [],
    }
    truth_frames: List[pd.DataFrame] = []

    for fold_index, (train_idx, test_idx) in enumerate(logo.split(data.X, data.Y, groups=groups), start=1):
        if max_folds is not None and fold_index > max_folds:
            break
        x_train = data.X.iloc[train_idx]
        y_train = data.Y.iloc[train_idx]
        x_test = data.X.iloc[test_idx]
        y_test = data.Y.iloc[test_idx]
        test_donors = sorted(set(data.sample_metadata.iloc[test_idx]["donor_id"].tolist()))

        model, scaler_x, scaler_y, target_scaling = _fit_main_model(x_train, y_train, model_config)
        main_pred = _predict_main_model(model, scaler_x, scaler_y, target_scaling, x_test, list(y_train.columns))

        mean_baseline = MeanLipidBaseline().fit(y_train)
        mean_pred = mean_baseline.predict(len(x_test), x_test.index)

        pca_ridge = PcaRidgeBaseline(
            max_components=int(baseline_config["pca_ridge"]["max_components"]),
            ridge_alphas=baseline_config["pca_ridge"]["ridge_alphas"],
        ).fit(x_train, y_train)
        pca_pred = pca_ridge.predict(x_test)

        truth_frames.append(y_test)
        pred_frames["multitask_elastic_net"].append(main_pred)
        pred_frames["mean_lipid_profile"].append(mean_pred)
        pred_frames["pca_ridge"].append(pca_pred)

        for model_name, frame in [
            ("multitask_elastic_net", main_pred),
            ("mean_lipid_profile", mean_pred),
            ("pca_ridge", pca_pred),
        ]:
            metrics = regression_metrics(y_test, frame)["global"]
            fold_records.append({
                "fold": fold_index,
                "model": model_name,
                "held_out_donors": ",".join(test_donors),
                "n_test_samples": int(len(test_idx)),
                **json_ready(metrics),
            })

    if not truth_frames:
        raise ValueError("No folds were evaluated. Check that donor groups are present.")

    truth_all = pd.concat(truth_frames, axis=0)
    aggregated_metrics = {}
    for model_name, frames in pred_frames.items():
        pred_all = pd.concat(frames, axis=0).loc[truth_all.index]
        aggregated_metrics[model_name] = json_ready(regression_metrics(truth_all, pred_all)["global"])
        pred_all.to_csv(output_dir / f"{model_name}.heldout_predictions.tsv", sep="\t")

    fold_df = pd.DataFrame(fold_records)
    fold_df.to_csv(output_dir / "internal_holdout_fold_metrics.tsv", sep="\t", index=False)
    truth_all.to_csv(output_dir / "internal_holdout_truth.tsv", sep="\t")

    summary = {
        "grouping": "donor_id",
        "fold_count": int(fold_df["fold"].nunique()),
        "models": aggregated_metrics,
        "truth_shape": [int(truth_all.shape[0]), int(truth_all.shape[1])],
    }
    dump_json(output_dir / "internal_holdout_summary.json", summary)
    return summary


def _map_external_celltype(value: str) -> str:
    value = str(value).lower()
    if ("pmn" in value) or ("polymorph" in value) or ("neutro" in value):
        return "PMX"
    if ("macrophage" in value) or ("mono" in value) or ("dc" in value) or ("dendritic" in value):
        return "MIC"
    if ("endothelial" in value) or ("cap" in value) or ("vein" in value) or ("arter" in value):
        return "END"
    if ("fibro" in value) or ("myofibro" in value) or ("pericyte" in value) or ("smooth_muscle" in value) or ("stromal" in value):
        return "MES"
    if ("at1" in value) or ("at2" in value) or ("club" in value) or ("ciliated" in value) or ("epithelial" in value):
        return "EPI"
    return "OTHER"


def _external_pairwise_contrasts(predictions: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    grouped_types = sorted([
        celltype
        for celltype, count in metadata["coarse_cell_type"].value_counts().items()
        if count >= 2 and celltype != "OTHER"
    ])
    for group1, group2 in combinations(grouped_types, 2):
        idx1 = metadata.index[metadata["coarse_cell_type"] == group1]
        idx2 = metadata.index[metadata["coarse_cell_type"] == group2]
        for lipid in predictions.columns:
            vals1 = predictions.loc[idx1, lipid].to_numpy(dtype=float)
            vals2 = predictions.loc[idx2, lipid].to_numpy(dtype=float)
            mean1 = float(np.mean(vals1))
            mean2 = float(np.mean(vals2))
            # The lipid values used here are already on a log-transformed scale,
            # so log2 fold-change is the direct difference between group means.
            log2_fc = float(mean1 - mean2)
            rows.append({
                "group1": group1,
                "group2": group2,
                "lipid": lipid,
                "log2_fc": log2_fc,
                "effect_size": float(abs(mean1 - mean2)),
                "direction": f"{group1}_up" if log2_fc >= 0 else f"{group2}_up",
            })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["group1", "group2", "effect_size"], ascending=[True, True, False]).reset_index(drop=True)


def evaluate_external_pediatric(
    *,
    config_path: str | Path,
    bundle_path: str | Path,
    output_dir: str | Path,
) -> Dict[str, Any]:
    config = load_json(config_path)
    config_file = Path(config["_config_path"])
    base_dir = config_file.parent
    external_cfg = config["evaluation"]["external_pediatric"]
    counts_path = resolve_path(external_cfg["counts_path"], relative_to=base_dir)
    groups_path = resolve_path(external_cfg["groups_path"], relative_to=base_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(counts_path, sep="\t")
    groups = pd.read_csv(groups_path, sep="\t")
    counts = counts.set_index(counts.columns[0])
    counts.index = counts.index.astype(str).str.strip()
    groups.columns = ["sample", "cell_type"]
    groups["sample"] = groups["sample"].astype(str).str.strip()
    groups["cell_type"] = groups["cell_type"].astype(str).str.strip()
    groups = groups.loc[groups["sample"].isin(counts.columns)].copy()
    counts = counts.loc[:, groups["sample"]]

    groups["donor_id"] = groups["sample"].str.split("__").str[0]
    pseudobulk = []
    for (donor_id, cell_type), rows in groups.groupby(["donor_id", "cell_type"]):
        expr = counts.loc[:, rows["sample"]].mean(axis=1)
        pseudobulk.append(expr.rename(f"{donor_id}_{cell_type.replace(' ', '_').replace('/', '_').replace('-', '_')}"))
    pseudobulk_df = pd.DataFrame(pseudobulk).apply(pd.to_numeric, errors="coerce").fillna(0.0)

    bundle = load_bundle(bundle_path)
    prediction_result = bundle.predict_from_dataframe(pseudobulk_df)
    predictions = prediction_result.predictions

    metadata = pd.DataFrame(index=predictions.index)
    metadata["donor_id"] = metadata.index.to_series().str.split("_").str[0]
    metadata["raw_cell_type"] = metadata.index.to_series().str.split("_", n=1).str[1]
    metadata["coarse_cell_type"] = metadata["raw_cell_type"].map(_map_external_celltype)

    valid_mask = metadata["coarse_cell_type"] != "OTHER"
    silhouette = None
    if valid_mask.sum() > 2 and metadata.loc[valid_mask, "coarse_cell_type"].nunique() > 1:
        silhouette = float(
            silhouette_score(
                predictions.loc[valid_mask].to_numpy(dtype=float),
                metadata.loc[valid_mask, "coarse_cell_type"].to_numpy(),
                metric="euclidean",
            )
        )

    contrasts = _external_pairwise_contrasts(predictions, metadata)
    if contrasts.empty:
        top_contrasts = pd.DataFrame(columns=["group1", "group2", "lipid", "log2_fc", "effect_size", "direction"])
    else:
        top_contrasts = contrasts.groupby(["group1", "group2"], as_index=False).head(1).reset_index(drop=True)

    predictions.to_csv(output_dir / "external_pediatric_predictions.tsv", sep="\t")
    metadata.to_csv(output_dir / "external_pediatric_metadata.tsv", sep="\t")
    contrasts.to_csv(output_dir / "external_pediatric_pairwise_contrasts.tsv", sep="\t", index=False)
    top_contrasts.to_csv(output_dir / "external_pediatric_top_contrasts.tsv", sep="\t", index=False)

    summary = {
        "bundle_path": str(Path(bundle_path).resolve()),
        "prediction_summary": prediction_result.summary,
        "profile_count": int(predictions.shape[0]),
        "coarse_cell_type_counts": {str(key): int(value) for key, value in metadata["coarse_cell_type"].value_counts().sort_index().items()},
        "silhouette_score_coarse_cell_types": silhouette,
        "top_pairwise_markers": top_contrasts.to_dict(orient="records"),
        "note": "This external pediatric evaluation is a structure-preservation report on predicted lipids, not a direct lipid-truth benchmark.",
    }
    dump_json(output_dir / "external_pediatric_summary.json", json_ready(summary))
    return summary


def run_evaluation(
    *,
    config_path: str | Path,
    bundle_path: str | Path,
    output_dir: str | Path,
    max_folds: int | None = None,
    skip_external: bool = False,
) -> Dict[str, Any]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    internal_summary = evaluate_internal_holdout(
        config_path=config_path,
        output_dir=output_dir / "internal_holdout",
        max_folds=max_folds,
    )
    external_summary = None
    if not skip_external:
        external_summary = evaluate_external_pediatric(
            config_path=config_path,
            bundle_path=bundle_path,
            output_dir=output_dir / "external_pediatric",
        )

    summary = {
        "internal_holdout": internal_summary,
        "external_pediatric": external_summary,
    }
    dump_json(output_dir / "evaluation_summary.json", json_ready(summary))
    return summary
