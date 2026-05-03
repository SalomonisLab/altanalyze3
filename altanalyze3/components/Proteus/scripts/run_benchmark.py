"""
run_benchmark.py: Comprehensive Proteus benchmark against baselines.

Evaluates the trained Proteus model on held-out test data across all tasks,
computes per-task AUROC/AUPRC/F1/ECE with bootstrap 95% confidence intervals,
and compares against three baselines:
  1. SequenceDeltaBaseline (ESM2 difference + logistic regression)
  2. LengthRatioBaseline (alt/ref length ratio)
  3. StructuralAnnotationBaseline (rule-based from structural features)

Generates two output files:
  - benchmark_results.tsv    : flat metrics table (one row per task × model × stratum)
  - benchmark_summary.txt    : formatted human-readable summary

Usage
-----
python scripts/run_benchmark.py \
    --checkpoint data/processed/finetune/best_model.pt \
    --test_tsv   data/processed/proteus_test.tsv \
    --train_tsv  data/processed/proteus_train.tsv \
    --cache_dir  data/interim \
    --output_dir data/processed/benchmark \
    --device     cuda

# CPU-only (slower)
python scripts/run_benchmark.py \
    --checkpoint data/processed/finetune/best_model.pt \
    --test_tsv   data/processed/proteus_test.tsv \
    --cache_dir  data/interim \
    --device     cpu
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Run Proteus benchmark")
    p.add_argument("--checkpoint", required=True, help="Trained Proteus model checkpoint (.pt)")
    p.add_argument("--test_tsv", required=True, help="Test split TSV")
    p.add_argument("--train_tsv", default=None,
                   help="Train split TSV (needed to fit parametric baselines). "
                        "If omitted, only StructuralAnnotationBaseline is run.")
    p.add_argument("--cache_dir", required=True, help="Root cache dir (interim/)")
    p.add_argument("--output_dir", default="data/processed/benchmark")
    p.add_argument("--device", default="cpu", choices=["cpu", "cuda", "mps"])
    p.add_argument("--batch_size", type=int, default=64)
    p.add_argument("--n_bootstrap", type=int, default=500,
                   help="Bootstrap resamples for CI (default 500, use 100 for quick run)")
    p.add_argument("--skip_baselines", action="store_true",
                   help="Skip baseline evaluation (faster, Proteus only)")
    return p.parse_args()


def main():
    args = parse_args()

    script_dir = Path(__file__).parent
    proteus_root = script_dir.parent
    if str(proteus_root) not in sys.path:
        sys.path.insert(0, str(proteus_root))

    import torch
    from torch.utils.data import DataLoader

    from proteus.model import ProteusModel
    from proteus.data.dataset import ProteusDataset
    from proteus.data.collate import proteus_collate_fn
    from proteus.benchmarking import ProteusEvaluator
    from proteus.benchmarking.baselines import (
        SequenceDeltaBaseline, LengthRatioBaseline, StructuralAnnotationBaseline,
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    device = args.device
    if device == "cuda" and not torch.cuda.is_available():
        print("WARNING: CUDA requested but not available, falling back to CPU.")
        device = "cpu"

    # ---- Load test dataset ----
    print(f"Loading test data from {args.test_tsv} ...")
    test_ds = ProteusDataset(
        tsv_path=Path(args.test_tsv),
        cache_dir=Path(args.cache_dir),
        split="test",
    )
    test_loader = DataLoader(
        test_ds,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=0,
        collate_fn=proteus_collate_fn,
    )

    # ---- Evaluate Proteus ----
    print(f"\n=== Evaluating Proteus ===")
    evaluator = ProteusEvaluator(device=device, n_bootstrap=args.n_bootstrap)
    proteus_results = evaluator.run(
        model_checkpoint=Path(args.checkpoint),
        test_tsv=Path(args.test_tsv),
        cache_dir=Path(args.cache_dir),
    )
    evaluator.print_summary(proteus_results)

    all_dfs = [evaluator.to_dataframe(proteus_results)]

    # ---- Baselines ----
    if not args.skip_baselines:
        # Optional: fit parametric baselines on training data
        if args.train_tsv:
            print(f"\n=== Fitting parametric baselines on {args.train_tsv} ===")
            train_ds = ProteusDataset(
                tsv_path=Path(args.train_tsv),
                cache_dir=Path(args.cache_dir),
                split="train",
            )
            train_loader = DataLoader(
                train_ds,
                batch_size=args.batch_size,
                shuffle=False,
                num_workers=0,
                collate_fn=proteus_collate_fn,
            )

            seq_baseline = SequenceDeltaBaseline(C=0.1)
            seq_baseline.fit(train_loader)

            len_baseline = LengthRatioBaseline()
            len_baseline.fit(train_loader)
        else:
            seq_baseline = None
            len_baseline = None

        struct_baseline = StructuralAnnotationBaseline()

        # Evaluate each baseline
        baselines = {"StructuralAnnotation": struct_baseline}
        if seq_baseline:
            baselines["SequenceDelta"] = seq_baseline
        if len_baseline:
            baselines["LengthRatio"] = len_baseline

        for bname, baseline in baselines.items():
            print(f"\n=== Evaluating {bname} baseline ===")
            all_preds = baseline.predict(test_loader)

            # Build labels + masks from test dataset
            import numpy as np
            from proteus.heads.task_heads import BINARY_TASKS

            all_labels: dict = {}
            all_masks: dict = {}
            for batch in test_loader:
                for task in BINARY_TASKS:
                    lk = "label_preservation" if task == "global_preservation" else f"label_{task}"
                    lab = batch.get(lk, batch.get("label_preservation"))
                    if lab is not None:
                        all_labels.setdefault(task, []).append(lab.numpy().flatten())
                    mk = batch.get("task_masks", {}).get(task)
                    if mk is not None:
                        all_masks.setdefault(task, []).append(mk.numpy().flatten())

            concat_labels = {k: np.concatenate(v) for k, v in all_labels.items() if v}
            concat_masks = {k: np.concatenate(v) for k, v in all_masks.items() if v}

            import pandas as pd
            rows = []
            for task in BINARY_TASKS:
                preds = all_preds.get(task)
                labels = concat_labels.get(task)
                masks = concat_masks.get(task)
                if preds is None or labels is None:
                    continue
                tm = evaluator._compute_task_metrics(
                    task, preds, labels, masks, "all", bname
                )
                rows.append(vars(tm))
            if rows:
                baseline_df = pd.DataFrame(rows)
                all_dfs.append(baseline_df)
                print(f"[{bname}] Tasks evaluated: {len(rows)}")
                for row in rows:
                    auroc = row.get("auroc", float("nan"))
                    auprc = row.get("auprc", float("nan"))
                    if not (auroc != auroc):  # not NaN
                        print(f"  {row['task']:<30} AUROC={auroc:.4f}  AUPRC={auprc:.4f}")

    # ---- Write outputs ----
    import pandas as pd
    combined = pd.concat(all_dfs, ignore_index=True)
    results_path = output_dir / "benchmark_results.tsv"
    combined.to_csv(results_path, sep="\t", index=False)
    print(f"\nBenchmark results written to: {results_path}")

    # Summary text
    summary_path = output_dir / "benchmark_summary.txt"
    with open(summary_path, "w") as f:
        f.write(f"Proteus Benchmark\n")
        f.write(f"Test pairs: {proteus_results.n_total_pairs:,}\n")
        f.write(f"Checkpoint: {args.checkpoint}\n\n")

        # Per-model AUROC table
        models = combined["model_name"].unique()
        tasks = combined[combined["stratum"] == "all"]["task"].unique()

        header = f"{'Task':<30} " + " ".join(f"{m:<14}" for m in models)
        f.write(header + "\n")
        f.write("-" * len(header) + "\n")

        for task in tasks:
            row_str = f"{task:<30} "
            for model in models:
                sub = combined[(combined["task"] == task) &
                               (combined["model_name"] == model) &
                               (combined["stratum"] == "all")]
                if not sub.empty and not (sub["auroc"].values[0] != sub["auroc"].values[0]):
                    row_str += f"{sub['auroc'].values[0]:<14.4f} "
                else:
                    row_str += f"{'N/A':<14} "
            f.write(row_str + "\n")

        if not (proteus_results.deviation_class_accuracy != proteus_results.deviation_class_accuracy):
            f.write(f"\nDeviation class accuracy: {proteus_results.deviation_class_accuracy:.4f}\n")
            f.write(f"Deviation class macro-F1: {proteus_results.deviation_class_f1_macro:.4f}\n")

    print(f"Summary written to: {summary_path}")


if __name__ == "__main__":
    main()
