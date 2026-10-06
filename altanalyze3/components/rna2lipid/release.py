"""Fail-closed identity and estimator checks for the promoted lung bundle.

This promotes an already fitted, validated model; it does not fit a candidate
or relax the historical training/source gates in candidate_integrity.py.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
RELEASE_BUNDLE_PATH = HERE / "rna2lipid_hs_lung_native_log2_20261005_bundle.pkl"
RELEASE_MANIFEST_PATH = HERE / "native_log2_release.json"
MANIFEST_SHA256 = "4d8adc116d7de0e2abbfc10aa2d3ca1d6dac8b8f5efe544bf6af2506b6468dd7"


def file_sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def release_manifest():
    if file_sha256(RELEASE_MANIFEST_PATH) != MANIFEST_SHA256:
        raise ValueError("Lipid release manifest changed after review; stop and reconcile it.")
    return json.loads(RELEASE_MANIFEST_PATH.read_text())


def exact_identifiers(expected, actual, label):
    if list(actual) != list(expected) or len(set(actual)) != len(actual):
        raise ValueError(f"{label} differs from the complete approved release; do not discard or intersect identities.")


def verify_release_bundle(path, bundle):
    manifest = release_manifest()
    if file_sha256(path) != manifest["bundle_sha256"]:
        raise ValueError("The default lipid bundle is not the exact validated release.")
    for key, label in [("X_columns", "RNA inputs"), ("Y_columns", "Lipid outputs"),
                       ("training_samples", "Training profiles")]:
        exact_identifiers(manifest[key], bundle[key], label)
    exact_identifiers(manifest["Y_columns"], bundle["models"], "Lipid estimators")
    for key, columns in [("scaler_x", "X_columns"), ("scaler_y", "Y_columns")]:
        scaler = bundle[key]
        if type(scaler) is not StandardScaler or scaler.n_features_in_ != len(manifest[columns]):
            raise ValueError(f"{key} is not the validated StandardScaler.")
        exact_identifiers(manifest[columns], scaler.feature_names_in_, key + " feature order")
        if not np.all(np.asarray(scaler.n_samples_seen_) == 45):
            raise ValueError(f"{key} has an unapproved training roster size.")
        if not all(np.isfinite(getattr(scaler, field)).all() for field in ["mean_", "scale_", "var_"]):
            raise ValueError(f"{key} contains nonfinite values.")
    settings = manifest["method"]["fit_settings"]
    for lipid, entry in bundle["models"].items():
        model = entry["model"]
        if type(model) is not ElasticNetCV:
            raise ValueError(f"{lipid} is not the established ElasticNetCV estimator.")
        genes = entry["genes"]
        if not genes or len(set(genes)) != len(genes) or not set(genes).issubset(manifest["X_columns"]):
            raise ValueError(f"{lipid} has unresolved selected RNA identities.")
        exact_identifiers(genes, model.feature_names_in_, lipid + " selected RNA order")
        expected = {"cv": settings["cv_folds"], "max_iter": settings["max_iter"],
                    "random_state": settings["random_seed"], "n_jobs": settings["n_jobs"],
                    **manifest["method"]["estimator_defaults"], "positive": False}
        for key, value in expected.items():
            if getattr(model, key) != value:
                raise ValueError(f"{lipid}: ElasticNetCV {key} differs from the baseline.")
        np.testing.assert_array_equal(model.alphas, settings["alpha_grid"])
        np.testing.assert_array_equal(model.l1_ratio, settings["l1_ratio_grid"])
        if model.coef_.shape != (len(genes),) or not np.isfinite(model.coef_).all() or not np.isfinite(model.intercept_):
            raise ValueError(f"{lipid}: incomplete/nonfinite fitted coefficients.")
    return manifest
