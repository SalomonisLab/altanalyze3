"""ADDITION, not the delivered protocol.

The delivered holdout rule caps each donor at ``max_per_donor`` profiles in the
holdout set. A donor with more profiles than that cap therefore appears in both
the training and the holdout set, so the delivered Pearson r and R2 are not
donor-held-out numbers.

This script measures the size of that effect. It runs the SAME trainer, on the
SAME matrices, with the SAME hyperparameter grid, and changes only the holdout
rule: the delivered constrained split versus a donor-disjoint split.

The lipid tables used here are the two tables this repository carries, not the
log2 tables the shipped bundle was trained on, so these numbers describe the
size of the leakage effect and are NOT the shipped model's performance.

Run:
  /usr/bin/python3 components/rna2lipid/validation/compare_holdout_leakage.py
"""
from __future__ import annotations

import json
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

MODULE_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = MODULE_DIR.parent.parent.parent
CONFIG = MODULE_DIR / "validation" / "config_smoke_generic_tissue.json"

FIT_KWARGS = dict(
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

HOLDOUT_MODES = {
    "delivered_constrained": {
        "mode": "constrained",
        "seed": 42,
        "total_holdout": 28,
        "min_per_group": 5,
        "max_per_donor": 3,
        "required_groups": ["EPI", "MES", "MIC", "END", "PMX"],
    },
    "addition_donor_disjoint": {
        "mode": "donor_disjoint",
        "seed": 42,
        "holdout_donor_fraction": 0.3,
        "min_holdout_donors": 2,
    },
}


def main() -> int:
    sys.path.insert(0, str(REPO_ROOT))
    from altanalyze3.components.rna2lipid.pipeline import build_dataset
    from altanalyze3.components.rna2lipid.training import (
        build_holdout,
        fit_sparse_lipidwise_elasticnet,
        holdout_global_metrics,
    )

    data = build_dataset(CONFIG)
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "config": str(CONFIG),
        "label": "ADDITION. The delivered protocol is delivered_constrained.",
        "matrices": {
            "X": list(data.X.shape),
            "Y": list(data.Y.shape),
            "donors": int(data.sample_metadata["donor_id"].nunique()),
        },
        "hyperparameters": {
            k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in FIT_KWARGS.items()
        },
        "modes": {},
    }

    for name, holdout_config in HOLDOUT_MODES.items():
        holdout_samples, meta_train, meta_holdout = build_holdout(data.sample_metadata, holdout_config)
        shared = sorted(set(meta_train["donor_id"]).intersection(set(meta_holdout["donor_id"])))
        print(
            f"{name}: train {len(meta_train)} profiles / {meta_train['donor_id'].nunique()} donors, "
            f"holdout {len(meta_holdout)} profiles / {meta_holdout['donor_id'].nunique()} donors, "
            f"donors in both = {len(shared)}",
            flush=True,
        )
        predictions, _, seconds = fit_sparse_lipidwise_elasticnet(
            data.X.loc[meta_train.index],
            data.Y.loc[meta_train.index],
            data.X.loc[meta_holdout.index],
            **FIT_KWARGS,
        )
        metrics = holdout_global_metrics(
            data.Y.loc[meta_holdout.index], predictions, model_name="SparseLipidwiseElasticNetCV"
        )
        report["modes"][name] = {
            "holdout_rule": holdout_config,
            "train_profiles": int(len(meta_train)),
            "holdout_profiles": int(len(meta_holdout)),
            "train_donors": int(meta_train["donor_id"].nunique()),
            "holdout_donors": int(meta_holdout["donor_id"].nunique()),
            "donors_in_both_sets": shared,
            "donor_disjoint": bool(not shared),
            "metrics": metrics,
            "fit_seconds": seconds,
        }

    a = report["modes"]["delivered_constrained"]["metrics"]
    b = report["modes"]["addition_donor_disjoint"]["metrics"]
    report["difference"] = {
        "Pearson_r_constrained_minus_donor_disjoint": float(a["Pearson_r"] - b["Pearson_r"]),
        "R2_constrained_minus_donor_disjoint": float(a["R2"] - b["R2"]),
        "note": (
            "Both runs share the trainer, the matrices and the hyperparameter grid. "
            "The holdout profile counts differ, so the gap mixes donor leakage with "
            "training-set size. Direction, not exact magnitude, is the claim."
        ),
    }

    out_path = MODULE_DIR / "validation" / "holdout_leakage_comparison.json"
    out_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["modes"]["delivered_constrained"]["metrics"], indent=2))
    print(json.dumps(report["modes"]["addition_donor_disjoint"]["metrics"], indent=2))
    print(json.dumps(report["difference"], indent=2))
    print(f"\nWrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
