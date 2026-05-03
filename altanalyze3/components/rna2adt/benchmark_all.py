"""End-to-end benchmark for all rna2adt candidates including the aberrant
expression test. Trains each candidate on the bone marrow atlas and reports
side-by-side cell-holdout, LODO[BM27], and synthetic-aberrant lifts.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import anndata as ad
import numpy as np
import pandas as pd

from .aberrant_test import report_to_dataframe, run_aberrant_test
from .training import (
    SplitData,
    _evaluate,
    _fit_centroid,
    _fit_generalizable,
    _fit_ridge,
    _fit_rna2lipid_arch,
    _fit_rna2lipid_arch_per_protein,
    _fit_targeted,
    _predict,
    load_split_arrays,
    save_bundle,
)


def _make_candidate_recipes() -> List[Dict[str, object]]:
    return [
        {"label": "centroid_only", "fn": _fit_centroid,
         "kwargs": dict(use_residual=False, residual_alpha=10.0)},
        {"label": "centroid_global_ridge_residual", "fn": _fit_centroid,
         "kwargs": dict(use_residual=True, residual_alpha=10.0, residual_kind="global_ridge")},
        {"label": "centroid_targeted_residual", "fn": _fit_centroid,
         "kwargs": dict(use_residual=True, residual_alpha=10.0, residual_kind="targeted")},
        {"label": "global_ridge", "fn": _fit_ridge, "kwargs": dict(alpha=1.0)},
        {"label": "targeted_ridge", "fn": _fit_targeted,
         "kwargs": dict(method="ridge", alpha=1.0, l1_ratio=0.5)},
        {"label": "targeted_elastic", "fn": _fit_targeted,
         "kwargs": dict(method="elastic", alpha=0.001, l1_ratio=0.5)},
        {"label": "genz_ridge", "fn": _fit_generalizable,
         "kwargs": dict(head_kind="ridge", alpha=1.0, n_pca=50, n_neighbors=25)},
        {"label": "genz_elastic", "fn": _fit_generalizable,
         "kwargs": dict(head_kind="elastic", alpha=0.001, l1_ratio=0.5,
                        n_pca=50, n_neighbors=25)},
        {"label": "genz_elastic_rank", "fn": _fit_generalizable,
         "kwargs": dict(head_kind="elastic_rank", alpha=0.001, l1_ratio=0.5,
                        n_pca=50, n_neighbors=25)},
        {"label": "genz_elastic_spline", "fn": _fit_generalizable,
         "kwargs": dict(head_kind="elastic_spline", alpha=0.001, l1_ratio=0.5,
                        n_pca=50, n_neighbors=25)},
        {"label": "genz_spline", "fn": _fit_generalizable,
         "kwargs": dict(head_kind="spline", alpha=1.0, n_pca=50, n_neighbors=25)},
        {"label": "genz_xgboost", "fn": _fit_generalizable,
         "kwargs": dict(head_kind="xgboost", alpha=1.0, n_pca=50, n_neighbors=25,
                        xgb_n_estimators=200, xgb_max_depth=4, xgb_lr=0.08)},
        {"label": "rna2lipid_arch_full", "fn": _fit_rna2lipid_arch,
         "kwargs": dict(feature_subset="full", cv=3, max_iter=1500, n_jobs=4,
                        l1_ratio=0.5, alphas=[0.001, 0.01, 0.1])},
        {"label": "rna2lipid_arch_panel", "fn": _fit_rna2lipid_arch,
         "kwargs": dict(feature_subset="panel", cv=3, max_iter=1500, n_jobs=4,
                        l1_ratio=0.5, alphas=[0.001, 0.01, 0.1])},
        {"label": "rna2lipid_arch_panel_loose", "fn": _fit_rna2lipid_arch,
         "kwargs": dict(feature_subset="panel", use_cv=False, fixed_alpha=0.001,
                        l1_ratio=0.3, max_iter=2000, label_suffix="loose")},
        {"label": "rna2lipid_arch_panel_minshrink", "fn": _fit_rna2lipid_arch,
         "kwargs": dict(feature_subset="panel", use_cv=False, fixed_alpha=0.0001,
                        l1_ratio=0.1, max_iter=3000, label_suffix="minshrink")},
        {"label": "rna2lipid_arch_panel_per_protein", "fn": _fit_rna2lipid_arch_per_protein,
         "kwargs": dict(alpha=0.01, l1_ratio=0.5, max_iter=2000, feature_source="curated")},
        {"label": "rna2lipid_arch_panel_per_protein_whitelist", "fn": _fit_rna2lipid_arch_per_protein,
         "kwargs": dict(alpha=0.01, l1_ratio=0.5, max_iter=2000, feature_source="whitelist")},
        {"label": "rna2lipid_arch_panel_per_protein_whitelist_union", "fn": _fit_rna2lipid_arch_per_protein,
         "kwargs": dict(alpha=0.01, l1_ratio=0.5, max_iter=2000, feature_source="whitelist_union")},
    ]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--train-h5ad", required=True, type=Path)
    p.add_argument("--heldout-h5ad", required=True, type=Path)
    p.add_argument("--lodo-donor", default="BM27")
    p.add_argument("--cluster-col", default="Level 3 Multimodal")
    p.add_argument("--max-train-cells", type=int, default=30000)
    p.add_argument("--max-test-cells", type=int, default=15000)
    p.add_argument("--rna-top-n", type=int, default=3000)
    p.add_argument("--output-dir", type=Path, default=Path("artifacts/all_benchmark"))
    p.add_argument("--bundle-dir", type=Path, default=None,
                   help="If set, save each trained candidate's bundle into this directory.")
    p.add_argument("--only", default="",
                   help="Optional comma-separated subset of candidate labels to run.")
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    from .adt_rna_map import load_curated_adt_rna_map
    from .training import _surface_marker_companions
    curated_partners = sorted({g for genes in load_curated_adt_rna_map().values() for g in genes})

    # Also pull in every gene listed in the empirical whitelist so the variance
    # filter doesn't drop them. Without this, whitelist genes that happen to be
    # below the top-N variance cutoff would silently disappear from rna_genes
    # and the per-protein head would have no features to fit.
    whitelist_extra: set = set()
    whitelist_path = Path(__file__).parent / "configs" / "empirical_whitelist.tsv"
    if whitelist_path.exists():
        import pandas as _pd
        wl = _pd.read_csv(whitelist_path, sep="\t")
        for entry in wl["feature_genes"].astype(str):
            for g in entry.split(","):
                g = g.strip()
                if g:
                    whitelist_extra.add(g)
        print(f"[bench] empirical whitelist contributes {len(whitelist_extra)} genes to extra_rna_genes", flush=True)

    extra_genes = sorted(set(_surface_marker_companions()) | set(curated_partners) | whitelist_extra)

    print(f"[bench] loading cell-holdout split", flush=True)
    cell_split = load_split_arrays(
        args.train_h5ad, test_path=args.heldout_h5ad,
        max_train_cells=args.max_train_cells,
        max_test_cells=args.max_test_cells,
        rna_top_n=args.rna_top_n,
        extra_rna_genes=extra_genes,
    )
    print(f"[bench] loading LODO split (held-out donor={args.lodo_donor})", flush=True)
    lodo_split = load_split_arrays(
        args.train_h5ad, leave_out_donor=args.lodo_donor,
        max_train_cells=args.max_train_cells,
        max_test_cells=args.max_test_cells,
        rna_top_n=args.rna_top_n,
        extra_rna_genes=extra_genes,
    )

    summary_rows: List[Dict[str, object]] = []
    aberrant_rows: List[Dict[str, object]] = []

    # We train each recipe twice: once on cell-holdout primary, once on LODO primary.
    # For aberrant test we use the LODO-trained model so the host cells are fully
    # held out (their donor was BM27, never seen in training).
    holdout_adata: ad.AnnData = ad.read_h5ad(args.heldout_h5ad)

    recipes = _make_candidate_recipes()
    if args.only:
        wanted = {k.strip() for k in args.only.split(",") if k.strip()}
        recipes = [r for r in recipes if r["label"] in wanted]
        if not recipes:
            raise SystemExit(f"--only={args.only!r} matched no candidates.")
    for recipe in recipes:
        label = recipe["label"]
        print(f"\n[bench] === {label} ===", flush=True)
        # Train on each split
        t0 = time.time()
        cand_cell = recipe["fn"](cell_split, **recipe["kwargs"])
        secs_cell = time.time() - t0
        t0 = time.time()
        cand_lodo = recipe["fn"](lodo_split, **recipe["kwargs"])
        secs_lodo = time.time() - t0

        pred_cell, _ = _predict(cand_cell, cell_split.X_test, cell_split.Y_test,
                                test_clusters=cell_split.test_clusters)
        pred_lodo, _ = _predict(cand_lodo, lodo_split.X_test, lodo_split.Y_test,
                                test_clusters=lodo_split.test_clusters)

        rep_cell = _evaluate(label, "cell_holdout", cell_split.Y_test, pred_cell,
                             cell_split.adt_names, fit_seconds=secs_cell)
        rep_lodo = _evaluate(label, f"LODO[{args.lodo_donor}]", lodo_split.Y_test, pred_lodo,
                             lodo_split.adt_names, fit_seconds=secs_lodo)
        print(f"[bench] {label} cell_holdout: mean_pearson={rep_cell.mean_pearson:.3f} "
              f"median={rep_cell.median_pearson:.3f}", flush=True)
        print(f"[bench] {label} LODO:         mean_pearson={rep_lodo.mean_pearson:.3f} "
              f"median={rep_lodo.median_pearson:.3f}", flush=True)
        summary_rows.append({
            "model": label, "split": "cell_holdout",
            "mean_pearson": round(rep_cell.mean_pearson, 4),
            "median_pearson": round(rep_cell.median_pearson, 4),
            "rmse": round(rep_cell.rmse, 4),
            "fit_s": round(secs_cell, 1),
        })
        summary_rows.append({
            "model": label, "split": f"LODO[{args.lodo_donor}]",
            "mean_pearson": round(rep_lodo.mean_pearson, 4),
            "median_pearson": round(rep_lodo.median_pearson, 4),
            "rmse": round(rep_lodo.rmse, 4),
            "fit_s": round(secs_lodo, 1),
        })

        # --- Aberrant expression test ---
        aberrant_reports = run_aberrant_test(
            model=cand_lodo.model,
            adata=holdout_adata,
            rna_genes=lodo_split.rna_genes,
            adt_names=lodo_split.adt_names,
            cluster_col=args.cluster_col,
        )
        df = report_to_dataframe(aberrant_reports)
        df.insert(0, "model", label)
        aberrant_rows.append(df)
        if not df.empty:
            print(f"[bench] {label} aberrant lift_fraction: "
                  f"mean={df['lift_fraction'].mean():.2f} probes={len(df)}", flush=True)

        # Optional bundle save
        if args.bundle_dir is not None:
            args.bundle_dir.mkdir(parents=True, exist_ok=True)
            save_bundle(
                cand_lodo,
                bundle_path=args.bundle_dir / f"{label}.pkl",
                rna_genes=lodo_split.rna_genes,
                adt_names=lodo_split.adt_names,
                metadata={"label": label, "split": "lodo_trained"},
            )

    # --- Persist combined reports ---
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    print("\n=== Summary across splits ===")
    print(summary_df.to_string(index=False))

    if aberrant_rows:
        aberrant_df = pd.concat(aberrant_rows, axis=0, ignore_index=True)
        aberrant_df.to_csv(args.output_dir / "aberrant_lift.tsv", sep="\t", index=False)
        print("\n=== Aberrant lift (LODO-trained models, boosted partner-gene RNA) ===")
        print(aberrant_df.to_string(index=False))


if __name__ == "__main__":
    main()
