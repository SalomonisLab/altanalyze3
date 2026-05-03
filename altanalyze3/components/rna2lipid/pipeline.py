from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping

import numpy as np
import pandas as pd
import sklearn

PACKAGE_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH = PACKAGE_DIR / "configs" / "default_training.json"


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def resolve_path(value: str | Path, *, relative_to: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (relative_to / path).resolve()


def load_json(path: str | Path) -> Dict[str, Any]:
    config_path = Path(path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Expected a JSON object in {config_path}")
    data["_config_path"] = str(config_path)
    return data


def dump_json(path: str | Path, payload: Mapping[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def clean_labels(values: Iterable[object]) -> List[str]:
    return [str(value).strip() for value in values]


def _standardize_frame(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.index = clean_labels(out.index)
    out.columns = clean_labels(out.columns)
    out = out.loc[~pd.Index(out.index).duplicated(keep="first")]
    return out


def _load_expression_table(path: Path, *, transpose: bool = False, collapse_duplicate_genes: bool = False) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    if transpose:
        df = df.T
    df = _standardize_frame(df)
    if collapse_duplicate_genes:
        df = df.T.groupby(level=0).mean().T
    return df


def _load_lipid_table(path: Path) -> pd.DataFrame:
    return _standardize_frame(pd.read_csv(path, index_col=0))


def _extract_donor(sample_id: str) -> str:
    sample_id = str(sample_id).strip()
    if not sample_id:
        return ""
    return sample_id.split("_", 1)[0]


def _extract_compartment(sample_id: str) -> str:
    sample_id = str(sample_id).strip()
    if "_" not in sample_id:
        return "bulk"
    return sample_id.rsplit("_", 1)[-1]


@dataclass(frozen=True)
class PreparedTrainingData:
    X: pd.DataFrame
    Y: pd.DataFrame
    sample_metadata: pd.DataFrame
    manifest: Dict[str, Any]


def build_training_dataset(config_path: str | Path = DEFAULT_CONFIG_PATH) -> PreparedTrainingData:
    config = load_json(config_path)
    config_file = Path(config["_config_path"])
    base_dir = config_file.parent
    data_cfg = config["data"]
    preprocessing_cfg = config["preprocessing"]

    bulk_rna_path = resolve_path(data_cfg["bulk_rna"], relative_to=base_dir)
    celltype_rna_path = resolve_path(data_cfg["celltype_rna"], relative_to=base_dir)
    bulk_lipid_path = resolve_path(data_cfg["bulk_lipids"], relative_to=base_dir)
    celltype_lipid_path = resolve_path(data_cfg["celltype_lipids"], relative_to=base_dir)

    bulk_rna = _load_expression_table(bulk_rna_path)
    celltype_rna = _load_expression_table(
        celltype_rna_path,
        transpose=True,
        collapse_duplicate_genes=bool(preprocessing_cfg.get("collapse_duplicate_genes", True)),
    )
    bulk_lipids = _load_lipid_table(bulk_lipid_path)
    celltype_lipids = _load_lipid_table(celltype_lipid_path)

    bulk_common_samples = bulk_rna.index.intersection(bulk_lipids.index)
    celltype_common_samples = celltype_rna.index.intersection(celltype_lipids.index)

    bulk_rna = bulk_rna.loc[bulk_common_samples].copy()
    bulk_lipids = bulk_lipids.loc[bulk_common_samples].copy()
    celltype_rna = celltype_rna.loc[celltype_common_samples].copy()
    celltype_lipids = celltype_lipids.loc[celltype_common_samples].copy()

    x_join_mode = str(preprocessing_cfg.get("x_join_mode", "inner")).strip().lower()
    y_join_mode = str(preprocessing_cfg.get("y_join_mode", "inner")).strip().lower()
    if x_join_mode not in {"inner", "outer"}:
        raise ValueError("x_join_mode must be 'inner' or 'outer'")
    if y_join_mode not in {"inner", "outer"}:
        raise ValueError("y_join_mode must be 'inner' or 'outer'")

    X = pd.concat([bulk_rna, celltype_rna], axis=0, join=x_join_mode)
    Y = pd.concat([bulk_lipids, celltype_lipids], axis=0, join=y_join_mode)

    X = X.loc[~pd.Index(X.index).duplicated(keep="first")].copy()
    Y = Y.loc[~pd.Index(Y.index).duplicated(keep="first")].copy()

    common_samples = X.index.intersection(Y.index)
    X = X.loc[common_samples].copy()
    Y = Y.loc[common_samples].copy()

    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    if preprocessing_cfg.get("remove_zero_only_genes", True):
        X = X.loc[:, (X != 0).any(axis=0)]

    Y = Y.apply(pd.to_numeric, errors="coerce")
    if preprocessing_cfg.get("drop_rows_with_missing_lipids", True):
        Y = Y.dropna(axis=0, how="any")
        X = X.loc[Y.index].copy()
    else:
        Y = Y.fillna(0.0)

    sample_ids = clean_labels(X.index)
    sample_metadata = pd.DataFrame(index=sample_ids)
    sample_metadata["sample_id"] = sample_ids
    sample_metadata["donor_id"] = sample_metadata["sample_id"].map(_extract_donor)
    sample_metadata["profile_kind"] = np.where(sample_metadata["sample_id"].str.contains("_"), "celltype", "bulk")
    sample_metadata["cell_compartment"] = sample_metadata["sample_id"].map(_extract_compartment)

    manifest = {
        "created_at": utc_timestamp(),
        "config_path": str(config_file),
        "data_paths": {
            "bulk_rna": str(bulk_rna_path),
            "celltype_rna": str(celltype_rna_path),
            "bulk_lipids": str(bulk_lipid_path),
            "celltype_lipids": str(celltype_lipid_path),
        },
        "preprocessing": {
            "collapse_duplicate_genes": bool(preprocessing_cfg.get("collapse_duplicate_genes", True)),
            "remove_zero_only_genes": bool(preprocessing_cfg.get("remove_zero_only_genes", True)),
            "drop_rows_with_missing_lipids": bool(preprocessing_cfg.get("drop_rows_with_missing_lipids", True)),
            "x_join_mode": x_join_mode,
            "y_join_mode": y_join_mode,
        },
        "source_shapes": {
            "bulk_rna": [int(bulk_rna.shape[0]), int(bulk_rna.shape[1])],
            "celltype_rna": [int(celltype_rna.shape[0]), int(celltype_rna.shape[1])],
            "bulk_lipids": [int(bulk_lipids.shape[0]), int(bulk_lipids.shape[1])],
            "celltype_lipids": [int(celltype_lipids.shape[0]), int(celltype_lipids.shape[1])],
        },
        "combined_shapes": {
            "X": [int(X.shape[0]), int(X.shape[1])],
            "Y": [int(Y.shape[0]), int(Y.shape[1])],
        },
        "counts": {
            "bulk_matched_samples": int(len(bulk_common_samples)),
            "celltype_matched_samples": int(len(celltype_common_samples)),
            "shared_lipids_between_bulk_and_celltype": int(len(bulk_lipids.columns.intersection(celltype_lipids.columns))),
            "unique_donors": int(sample_metadata["donor_id"].nunique()),
        },
        "sample_summary": {
            "profile_kind_counts": {str(key): int(value) for key, value in sample_metadata["profile_kind"].value_counts().sort_index().items()},
            "cell_compartment_counts": {str(key): int(value) for key, value in sample_metadata["cell_compartment"].value_counts().sort_index().items()},
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "sklearn": sklearn.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }

    return PreparedTrainingData(
        X=X,
        Y=Y,
        sample_metadata=sample_metadata,
        manifest=manifest,
    )


def json_ready(value: Any) -> Any:
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): json_ready(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def merge_json(base: MutableMapping[str, Any], updates: Mapping[str, Any]) -> MutableMapping[str, Any]:
    for key, value in updates.items():
        if isinstance(value, Mapping) and isinstance(base.get(key), MutableMapping):
            merge_json(base[key], value)
        else:
            base[key] = value
    return base
