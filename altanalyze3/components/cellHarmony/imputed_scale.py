"""Explicit prediction encodings shared by export, aggregation and differentials.

The source transform's pseudocount is reversed before computing means. It is
never added to those means to manufacture a finite abundance fold change.
"""
from __future__ import annotations

import numpy as np

IMPUTED_MODALITIES = {"lipids", "adt", "metabolite", "lipid", "grn", "grn_tf"}


def uses_imputed_scale(adata):
    return (adata.uns.get("modality") in IMPUTED_MODALITIES
            or adata.uns.get("expression_scale") == "native_relative_log2")


def normalize_scale(value):
    raw = str(value or "").strip().lower()
    return {"natural_log1p": "log1p",
            "log_2": "log2", "raw": "linear", "none": "linear"}.get(raw, raw)


def log_base_value(value, *, default):
    if value is None or str(value).strip().lower() in {"", "none", "null"}:
        base = float(default)
    elif str(value).strip().lower() in {"e", "ln", "natural", "loge"}:
        base = float(np.e)
    else:
        base = float(value)
    if not np.isfinite(base) or base <= 1:
        raise ValueError("The declared log base must be finite and greater than one.")
    return base


def prediction_encoding(adata):
    """Require explicit log2 offset; never infer an inverse from value ranges."""
    uns = adata.uns
    scale = normalize_scale(uns.get("expression_scale"))
    if scale == "linear":
        return scale, None, 0.0
    if scale in {"native_relative_log2", "log2p1"}:
        base, offset = 2.0, (0.0 if scale == "native_relative_log2" else 1.0)
    elif scale == "log1p":
        base = log_base_value(uns.get("log_base", uns.get("log1p", {}).get("base")), default=np.e)
        offset = 1.0
    elif scale == "log2" and "log_pseudocount" in uns:
        base, offset = 2.0, float(uns["log_pseudocount"])
    else:
        raise ValueError(f"{uns.get('modality', 'Imputed modality')}: target encoding is not fully declared. "
                         "Provide the target preprocessing protocol (log base and pseudocount); do not guess its inverse.")
    if not np.isfinite(offset) or offset < 0:
        raise ValueError("The declared source pseudocount must be finite and nonnegative.")
    if "log_base" in uns and not np.isclose(log_base_value(uns["log_base"], default=base), base):
        raise ValueError("Conflicting declared log bases; resolve the source metadata.")
    if "log1p" in uns and not np.isclose(log_base_value(uns["log1p"].get("base"), default=np.e), base):
        raise ValueError("Conflicting declared log bases; resolve the source metadata.")
    if "log_pseudocount" in uns and float(uns["log_pseudocount"]) != offset:
        raise ValueError("Conflicting declared source pseudocounts; resolve the source metadata.")
    return scale, base, offset


def inverse_predictions(values, encoding):
    scale, base, offset = encoding
    values = np.asarray(values, dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite predictions; do not fill values or drop features/samples.")
    if scale == "linear":
        return values
    with np.errstate(over="ignore", under="ignore"):
        linear = (np.expm1(values * np.log(base)) if offset == 1
                  else np.power(base, values) - offset)
    if not np.isfinite(linear).all() or np.any(linear < 0) or (offset == 0 and np.any(linear <= 0)):
        raise ValueError("Predictions cannot be inverted to finite nonnegative abundances in the declared scale; "
                         "do not clip values or add a pseudocount.")
    return linear


def encode_predictions(linear, encoding):
    scale, base, offset = encoding
    linear = np.asarray(linear, dtype=np.float64)
    if scale == "linear":
        return linear
    shifted = linear + offset
    if not np.isfinite(shifted).all() or np.any(shifted <= 0):
        raise ValueError("A mean cannot be represented in the declared log scale; do not add a new pseudocount.")
    return np.log(shifted) / np.log(base)


def fold_convention(adata):
    scale, base, offset = prediction_encoding(adata)
    return {"expression_scale": scale, "log_base": base, "source_pseudocount": offset,
            "fold_statistic": "ratio_of_arithmetic_mean_linear_predictions",
            "fold_pseudocount": 0.0, "feature_total_normalized": False,
            "BH_policy": "complete_imputed_feature_panel_no_abundance_filter",
            "zero_mean_policy": "one_zero: infinite_log2_fold; both_zero: undefined_log2_fold"}
