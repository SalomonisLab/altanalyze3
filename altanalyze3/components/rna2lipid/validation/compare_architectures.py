"""Head-to-head: the new sparse lipid-by-lipid ElasticNet against the prior
MultiTaskElasticNetCV, on identical data and an identical donor-disjoint split.

Each architecture runs with the hyperparameter grid it ships with:
  sparse_lipidwise -> configs/lung_lipidwise.json
  multitask        -> configs/default_training.json

A mean-lipid baseline gives the floor. Any model that cannot beat the per-lipid
training mean has R2 <= 0 on held-out data.

The lipid tables here are the two this repository carries, NOT the log2 tables
the shipped bundle was trained on. This measures which ARCHITECTURE generalizes
better on the same inputs. It is not the shipped model's absolute performance.

Run:
  /usr/bin/python3 components/rna2lipid/validation/compare_architectures.py
"""
from __future__ import annotations

import json
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

MODULE_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = MODULE_DIR.parent.parent.parent
CONFIG = MODULE_DIR / "validation" / "config_smoke_generic_tissue.json"

# Grid shipped with configs/lung_lipidwise.json
LIPIDWISE_GRID = dict(
    top_gene_options=[25, 50, 100, 200],
    l1_ratio_grid=[0.70, 0.90, 0.95, 1.00],
    alpha_grid=np.logspace(-3, 1, 30),
    cv_folds=3,
    max_iter=50000,
    sparsity_penalty=0.002,
    random_seed=1,
    n_jobs=-1,
    verbose=False,
)

# Grid shipped with configs/default_training.json
MULTITASK_GRID = {
    "target_scaling": "none",
    "l1_ratio": [0.2, 0.5, 0.8],
    "alphas": [0.001, 0.00215443469, 0.00464158883, 0.01, 0.0215443469,
               0.0464158883, 0.1, 0.215443469, 0.464158883, 1.0],
    "cv": 3,
    "max_iter": 30000,
    "n_jobs": -1,
}


def main() -> int:
    sys.path.insert(0, str(REPO_ROOT))
    from altanalyze3.components.rna2lipid.evaluation import _fit_main_model, _predict_main_model
    from altanalyze3.components.rna2lipid.pipeline import build_dataset, load_json
    from altanalyze3.components.rna2lipid.training import (
        build_holdout,
        fit_sparse_lipidwise_elasticnet,
        holdout_global_metrics,
    )

    # Keep every matched sample, bulk included, so both architectures see the
    # data the PRIOR model was trained on.
    config = load_json(CONFIG)
    config["sample_metadata"].pop("keep_groups", None)
    unfiltered = MODULE_DIR / "validation" / "_config_all_groups.json"
    unfiltered.write_text(json.dumps({k: v for k, v in config.items() if k != "_config_path"}, indent=2))

    data = build_dataset(unfiltered)
    holdout_samples, meta_train, meta_holdout = build_holdout(
        data.sample_metadata,
        {"mode": "donor_disjoint", "seed": 42, "holdout_donor_fraction": 0.3, "min_holdout_donors": 2},
    )
    x_train, y_train = data.X.loc[meta_train.index], data.Y.loc[meta_train.index]
    x_test, y_test = data.X.loc[meta_holdout.index], data.Y.loc[meta_holdout.index]
    shared = sorted(set(meta_train["donor_id"]) & set(meta_holdout["donor_id"]))

    print(f"X {list(data.X.shape)}  Y {list(data.Y.shape)}", flush=True)
    print(f"train {len(x_train)} / holdout {len(x_test)} / donors in both {len(shared)}", flush=True)

    results = {}

    print("Fitting sparse_lipidwise ...", flush=True)
    t0 = time.perf_counter()
    lipidwise_pred, _, _ = fit_sparse_lipidwise_elasticnet(x_train, y_train, x_test, **LIPIDWISE_GRID)
    results["sparse_lipidwise_NEW"] = holdout_global_metrics(
        y_test, lipidwise_pred, model_name="SparseLipidwiseElasticNetCV")
    results["sparse_lipidwise_NEW"]["fit_seconds"] = time.perf_counter() - t0

    print("Fitting multitask ...", flush=True)
    t0 = time.perf_counter()
    model, sx, sy, scaling = _fit_main_model(x_train, y_train, MULTITASK_GRID)
    multitask_pred = _predict_main_model(model, sx, sy, scaling, x_test, list(y_train.columns))
    results["multitask_PRIOR"] = holdout_global_metrics(
        y_test, multitask_pred, model_name="MultiTaskElasticNetCV")
    results["multitask_PRIOR"]["fit_seconds"] = time.perf_counter() - t0

    mean_pred = pd.DataFrame(
        np.repeat(y_train.mean(axis=0).to_numpy(float)[None, :], len(x_test), axis=0),
        index=x_test.index, columns=y_train.columns)
    results["mean_lipid_baseline"] = holdout_global_metrics(
        y_test, mean_pred, model_name="MeanLipidBaseline")

    # Per-lipid win counts, so one aggregate number does not hide the spread.
    def per_lipid_r2(pred):
        out = {}
        for lipid in y_test.columns:
            truth = y_test[lipid].to_numpy(float)
            if np.std(truth) == 0:
                continue
            ss_res = float(np.sum((truth - pred[lipid].to_numpy(float)) ** 2))
            ss_tot = float(np.sum((truth - truth.mean()) ** 2))
            out[lipid] = 1.0 - ss_res / ss_tot if ss_tot else np.nan
        return pd.Series(out)

    r2_new, r2_old = per_lipid_r2(lipidwise_pred), per_lipid_r2(multitask_pred)
    both = pd.concat([r2_new.rename("new"), r2_old.rename("old")], axis=1).dropna()

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "design": {
            "matrices": {"X": list(data.X.shape), "Y": list(data.Y.shape)},
            "train_profiles": int(len(x_train)),
            "holdout_profiles": int(len(x_test)),
            "train_donors": int(meta_train["donor_id"].nunique()),
            "holdout_donors": int(meta_holdout["donor_id"].nunique()),
            "donors_in_both_sets": shared,
            "donor_disjoint": not shared,
            "holdout_group_counts": {str(k): int(v) for k, v in meta_holdout["group"].value_counts().items()},
            "note": "Both architectures share the split, the matrices and the samples.",
        },
        "holdout_metrics": results,
        "per_lipid": {
            "lipids_scored": int(len(both)),
            "new_beats_old": int((both["new"] > both["old"]).sum()),
            "old_beats_new": int((both["old"] > both["new"]).sum()),
            "new_median_r2": float(both["new"].median()),
            "old_median_r2": float(both["old"].median()),
            "new_lipids_with_positive_r2": int((both["new"] > 0).sum()),
            "old_lipids_with_positive_r2": int((both["old"] > 0).sum()),
        },
    }
    out_path = MODULE_DIR / "validation" / "architecture_comparison.json"
    out_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print(f"\nWrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
