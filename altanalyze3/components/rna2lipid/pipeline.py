from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import re
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


# ---------------------------------------------------------------------------
# Tissue-agnostic paired-table dataset builder
#
# The legacy builder above hard-codes the four lung tables. A new tissue
# declares its own tables under ``data.tables``. Each entry pairs one RNA
# matrix with one lipid matrix that share sample identifiers.
# ---------------------------------------------------------------------------


ORIENTATIONS = {"samples_by_features", "features_by_samples"}


def _read_table(path: Path, *, orientation: str) -> pd.DataFrame:
    orientation = str(orientation).strip().lower()
    if orientation not in ORIENTATIONS:
        raise ValueError(
            f"orientation must be one of {sorted(ORIENTATIONS)}, got {orientation!r}"
        )
    df = pd.read_csv(path, index_col=0)
    if orientation == "features_by_samples":
        df = df.T
    return _standardize_frame(df)


def _derive_field(sample_ids: pd.Index, spec: Mapping[str, Any]) -> pd.Series:
    """Derive one metadata column from the sample identifiers."""
    strategy = str(spec.get("strategy", "split_field")).strip().lower()
    ids = pd.Series(list(sample_ids), index=list(sample_ids), dtype=str)

    if strategy == "constant":
        return ids.map(lambda _: str(spec["value"]))

    if strategy == "split_field":
        delimiter = str(spec.get("delimiter", "_"))
        field = int(spec.get("field", 0))
        maxsplit = spec.get("maxsplit")
        def _split(value: str) -> str:
            parts = value.split(delimiter) if maxsplit is None else value.split(delimiter, int(maxsplit))
            try:
                return parts[field]
            except IndexError:
                return str(spec.get("missing", ""))
        return ids.map(_split)

    if strategy == "regex":
        pattern = re.compile(str(spec["pattern"]))
        group = spec.get("group", 1)
        def _match(value: str) -> str:
            found = pattern.search(value)
            if found is None:
                return str(spec.get("missing", ""))
            return str(found.group(group))
        return ids.map(_match)

    raise ValueError(f"Unsupported metadata strategy: {strategy!r}")


def _derive_group(labels: pd.Series, spec: Mapping[str, Any]) -> pd.Series:
    """Map a per-sample label onto the group used for stratification."""
    strategy = str(spec.get("strategy", "label")).strip().lower()

    if strategy == "label":
        return labels.astype(str).str.strip()

    if strategy == "keyword_map":
        rules = list(spec.get("rules", []))
        default = str(spec.get("default", "UNKNOWN"))
        if not rules:
            raise ValueError("keyword_map requires a non-empty 'rules' list")
        def _assign(label: str) -> str:
            text = str(label).strip().lower()
            for rule in rules:
                group = str(rule["group"])
                exact = {str(item).strip().lower() for item in rule.get("equals", [])}
                if text in exact:
                    return group
                for token in rule.get("any_of", []):
                    if str(token).strip().lower() in text:
                        return group
            return default
        return labels.map(_assign)

    if strategy == "regex":
        return _derive_field(pd.Index(labels.index), spec)

    raise ValueError(f"Unsupported group strategy: {strategy!r}")


def build_sample_metadata(sample_ids: Iterable[str], config: Mapping[str, Any]) -> pd.DataFrame:
    """Build donor / label / group metadata from sample identifiers.

    ``config`` is the ``sample_metadata`` block. Defaults reproduce the lung
    convention ``<donor>_<compartment>``.
    """
    ids = pd.Index(clean_labels(sample_ids))
    metadata = pd.DataFrame(index=ids)
    metadata["sample_id"] = list(ids)

    donor_spec = dict(config.get("donor") or {"strategy": "split_field", "delimiter": "_", "field": 0})
    metadata["donor_id"] = _derive_field(ids, donor_spec).to_numpy()

    label_spec = dict(
        config.get("label")
        or {"strategy": "split_field", "delimiter": "_", "field": 1, "maxsplit": 1, "missing": ""}
    )
    metadata["label"] = _derive_field(ids, label_spec).to_numpy()

    group_spec = dict(config.get("group") or {"strategy": "label"})
    metadata["group"] = _derive_group(metadata["label"], group_spec).to_numpy()

    metadata["profile_kind"] = np.where(metadata["sample_id"].str.contains("_"), "celltype", "bulk")
    return metadata


def build_paired_training_dataset(config_path: str | Path) -> PreparedTrainingData:
    """Build X and Y from an arbitrary list of paired RNA / lipid tables.

    Every table entry contributes the samples present in BOTH its RNA and its
    lipid matrix. Tables are then stacked and intersected on features.
    """
    config = load_json(config_path)
    config_file = Path(config["_config_path"])
    base_dir = config_file.parent
    data_cfg = config["data"]
    preprocessing_cfg = config.get("preprocessing", {})
    metadata_cfg = config.get("sample_metadata", {})

    table_specs = list(data_cfg.get("tables") or [])
    if not table_specs:
        raise ValueError("data.tables must list at least one paired RNA/lipid table")

    collapse_duplicate_genes = bool(preprocessing_cfg.get("collapse_duplicate_genes", True))

    rna_frames: List[pd.DataFrame] = []
    lipid_frames: List[pd.DataFrame] = []
    table_report: List[Dict[str, Any]] = []

    for spec in table_specs:
        table_id = str(spec.get("id") or f"table_{len(table_report) + 1}")
        rna_path = resolve_path(spec["rna"], relative_to=base_dir)
        lipid_path = resolve_path(spec["lipids"], relative_to=base_dir)
        rna = _read_table(rna_path, orientation=spec.get("rna_orientation", "samples_by_features"))
        if collapse_duplicate_genes:
            rna = rna.T.groupby(level=0).mean().T
        lipids = _read_table(lipid_path, orientation=spec.get("lipid_orientation", "samples_by_features"))

        matched = rna.index.intersection(lipids.index)
        table_report.append({
            "id": table_id,
            "rna_path": str(rna_path),
            "lipid_path": str(lipid_path),
            "rna_shape": [int(rna.shape[0]), int(rna.shape[1])],
            "lipid_shape": [int(lipids.shape[0]), int(lipids.shape[1])],
            "matched_samples": int(len(matched)),
            "rna_samples_dropped_unmatched": int(rna.shape[0] - len(matched)),
            "lipid_samples_dropped_unmatched": int(lipids.shape[0] - len(matched)),
        })
        if len(matched) == 0:
            raise ValueError(
                f"Table {table_id!r} shares no sample identifiers between "
                f"{rna_path} and {lipid_path}"
            )
        rna_frames.append(rna.loc[matched].copy())
        lipid_frames.append(lipids.loc[matched].copy())

    x_join_mode = str(preprocessing_cfg.get("x_join_mode", "inner")).strip().lower()
    y_join_mode = str(preprocessing_cfg.get("y_join_mode", "inner")).strip().lower()
    if x_join_mode not in {"inner", "outer"} or y_join_mode not in {"inner", "outer"}:
        raise ValueError("x_join_mode and y_join_mode must be 'inner' or 'outer'")

    X = pd.concat(rna_frames, axis=0, join=x_join_mode)
    Y = pd.concat(lipid_frames, axis=0, join=y_join_mode)

    X = X.loc[~pd.Index(X.index).duplicated(keep="first")].copy()
    Y = Y.loc[~pd.Index(Y.index).duplicated(keep="first")].copy()

    common_samples = X.index.intersection(Y.index)
    X = X.loc[common_samples].copy()
    Y = Y.loc[common_samples].copy()

    X = X.apply(pd.to_numeric, errors="coerce")
    Y = Y.apply(pd.to_numeric, errors="coerce")

    missing_x_before = int(X.isna().to_numpy().sum())
    missing_y_before = int(Y.isna().to_numpy().sum())

    if preprocessing_cfg.get("fill_missing_with_feature_median", True):
        X = X.fillna(X.median(numeric_only=True))
        Y = Y.fillna(Y.median(numeric_only=True))
    if preprocessing_cfg.get("drop_rows_with_missing_lipids", False):
        Y = Y.dropna(axis=0, how="any")
        X = X.loc[Y.index].copy()

    genes_before = int(X.shape[1])
    X = X.loc[:, X.notna().sum(axis=0) > 0]
    if preprocessing_cfg.get("remove_zero_only_genes", True):
        X = X.loc[:, (X != 0).any(axis=0)]
    genes_after = int(X.shape[1])

    sample_metadata = build_sample_metadata(Y.index, metadata_cfg)

    keep_groups = metadata_cfg.get("keep_groups")
    samples_before_group_filter = int(sample_metadata.shape[0])
    dropped_group_counts: Dict[str, int] = {}
    if keep_groups:
        keep = [str(value) for value in keep_groups]
        dropped = sample_metadata.loc[~sample_metadata["group"].isin(keep)]
        dropped_group_counts = {
            str(key): int(value) for key, value in dropped["group"].value_counts().items()
        }
        sample_metadata = sample_metadata.loc[sample_metadata["group"].isin(keep)].copy()
        X = X.loc[sample_metadata.index].copy()
        Y = Y.loc[sample_metadata.index].copy()

    retained_fraction = (
        float(sample_metadata.shape[0]) / float(samples_before_group_filter)
        if samples_before_group_filter
        else 0.0
    )

    manifest = {
        "created_at": utc_timestamp(),
        "config_path": str(config_file),
        "builder": "build_paired_training_dataset",
        "tissue": config.get("tissue"),
        "tables": table_report,
        "preprocessing": {
            "collapse_duplicate_genes": collapse_duplicate_genes,
            "fill_missing_with_feature_median": bool(
                preprocessing_cfg.get("fill_missing_with_feature_median", True)
            ),
            "drop_rows_with_missing_lipids": bool(
                preprocessing_cfg.get("drop_rows_with_missing_lipids", False)
            ),
            "remove_zero_only_genes": bool(preprocessing_cfg.get("remove_zero_only_genes", True)),
            "x_join_mode": x_join_mode,
            "y_join_mode": y_join_mode,
        },
        "combined_shapes": {
            "X": [int(X.shape[0]), int(X.shape[1])],
            "Y": [int(Y.shape[0]), int(Y.shape[1])],
        },
        "counts": {
            "matched_samples_before_group_filter": samples_before_group_filter,
            "matched_samples_after_group_filter": int(sample_metadata.shape[0]),
            "sample_retained_fraction": retained_fraction,
            "samples_dropped_by_group_filter": dropped_group_counts,
            "genes_before_information_filter": genes_before,
            "genes_after_information_filter": genes_after,
            "missing_rna_values_before_fill": missing_x_before,
            "missing_lipid_values_before_fill": missing_y_before,
            "unique_donors": int(sample_metadata["donor_id"].nunique()),
        },
        "sample_summary": {
            "group_counts": {
                str(key): int(value)
                for key, value in sample_metadata["group"].value_counts().sort_index().items()
            },
            "profile_kind_counts": {
                str(key): int(value)
                for key, value in sample_metadata["profile_kind"].value_counts().sort_index().items()
            },
            "donor_counts": {
                str(key): int(value)
                for key, value in sample_metadata["donor_id"].value_counts().sort_index().items()
            },
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


def build_dataset(config_path: str | Path = DEFAULT_CONFIG_PATH) -> PreparedTrainingData:
    """Dispatch on the config shape.

    A config with ``data.tables`` uses the tissue-agnostic paired-table builder.
    A config with the four named lung keys uses the legacy builder, unchanged.
    """
    config = load_json(config_path)
    if config.get("data", {}).get("tables"):
        return build_paired_training_dataset(config_path)
    return build_training_dataset(config_path)
