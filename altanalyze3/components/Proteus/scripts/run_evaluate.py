"""
run_evaluate.py: Evaluate a trained ProteusModel checkpoint on the test set.

Reports per-task AUROC, AUPRC, Brier score, balanced accuracy,
and within-gene ranking accuracy. Saves a markdown evaluation report.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import torch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate a trained ProteusModel checkpoint on the Proteus test set."
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        required=True,
        help="Path to trained model checkpoint (.pt file).",
    )
    parser.add_argument(
        "--data_processed_dir",
        type=Path,
        default=Path("data/processed"),
        help="Directory containing proteus_test.tsv.",
    )
    parser.add_argument(
        "--cache_dir",
        type=Path,
        default=Path("data/interim"),
        help="Root cache directory for embeddings.",
    )
    parser.add_argument(
        "--output_report",
        type=Path,
        default=Path("data/processed/evaluation_report.md"),
        help="Path to save markdown evaluation report.",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=512,
        help="Batch size for evaluation. Default: 512",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="auto",
        help="Device (auto/cuda/mps/cpu). Default: auto",
    )
    parser.add_argument(
        "--split",
        type=str,
        default="test",
        choices=["train", "val", "test"],
        help="Which split to evaluate. Default: test",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sys.path.insert(0, str(Path(__file__).parent.parent))

    from proteus.model import ProteusModel
    from proteus.data.dataset import ProteusDataset
    from proteus.data.collate import proteus_collate_fn
    from proteus.training.metrics import (
        compute_task_metrics,
        within_gene_ranking_accuracy,
    )
    from proteus.heads.task_heads import BINARY_TASKS, DEVIATION_CLASSES
    from proteus.utils.device import get_device

    device = get_device(args.device)
    print(f"Device: {device}")

    # Load model
    print(f"Loading checkpoint: {args.checkpoint}")
    model = ProteusModel.load_checkpoint(args.checkpoint)
    model.to(device)
    model.eval()

    # Dataset
    tsv_path = args.data_processed_dir / f"proteus_{args.split}.tsv"
    if not tsv_path.exists():
        print(f"ERROR: {args.split} TSV not found: {tsv_path}")
        sys.exit(1)

    print(f"Loading {args.split} dataset from: {tsv_path}")
    dataset = ProteusDataset(tsv_path=tsv_path, cache_dir=args.cache_dir)

    import torch.utils.data as tud
    loader = tud.DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=2,
        collate_fn=proteus_collate_fn,
    )

    # Collect predictions
    all_outputs: dict = {}
    all_labels_preservation = []
    all_deviation_labels = []
    all_gene_ids = []
    all_task_masks: dict = {}

    print(f"Running evaluation on {len(dataset)} samples...")
    with torch.no_grad():
        for batch in loader:
            # Move to device
            for k, v in batch.items():
                if isinstance(v, torch.Tensor):
                    batch[k] = v.to(device)
                elif isinstance(v, dict):
                    batch[k] = {
                        kk: vv.to(device) if isinstance(vv, torch.Tensor) else vv
                        for kk, vv in v.items()
                    }

            outputs = model(batch)

            for task in BINARY_TASKS:
                if task in outputs:
                    all_outputs.setdefault(task, []).append(
                        torch.sigmoid(outputs[task]).detach().cpu()
                    )

            if "deviation_class" in outputs:
                all_outputs.setdefault("deviation_class", []).append(
                    torch.softmax(outputs["deviation_class"], dim=-1).detach().cpu()
                )

            lbl = batch.get("label_preservation")
            if lbl is not None:
                all_labels_preservation.append(lbl.cpu())

            dev = batch.get("deviation_class")
            if dev is not None:
                all_deviation_labels.append(dev.cpu())

            all_gene_ids.extend(batch.get("gene_id", []))

            task_masks = batch.get("task_masks", {})
            for task, mask in task_masks.items():
                all_task_masks.setdefault(task, []).append(mask.cpu())

    # Concatenate
    cat_outputs = {k: torch.cat(vs, dim=0) for k, vs in all_outputs.items()}
    y_true = (
        torch.cat(all_labels_preservation, dim=0).numpy()
        if all_labels_preservation else np.array([])
    )
    dev_true = (
        torch.cat(all_deviation_labels, dim=0).numpy()
        if all_deviation_labels else np.array([])
    )
    cat_masks = {
        k: torch.cat(vs, dim=0).numpy()
        for k, vs in all_task_masks.items()
    }

    # Compute metrics
    results: dict = {}
    print("\n--- Per-Task Metrics ---")

    for task in BINARY_TASKS:
        if task not in cat_outputs:
            continue
        probs = cat_outputs[task].numpy().squeeze(-1)
        mask = cat_masks.get(task)
        y = y_true.copy()
        if mask is not None:
            y[mask < 0.5] = -1

        m = compute_task_metrics(y, probs, task_name=task)
        results[task] = m

        auroc = f"{m['auroc']:.4f}" if m['auroc'] is not None else "N/A"
        auprc = f"{m['auprc']:.4f}" if m['auprc'] is not None else "N/A"
        brier = f"{m['brier']:.4f}" if m['brier'] is not None else "N/A"
        bal_acc = f"{m['balanced_accuracy']:.4f}" if m['balanced_accuracy'] is not None else "N/A"

        print(
            f"  {task:25s}  n={m['n']:5d}  "
            f"AUROC={auroc}  AUPRC={auprc}  "
            f"Brier={brier}  BalAcc={bal_acc}"
        )

    # Within-gene ranking accuracy
    wgra = None
    if all_gene_ids and "global_preservation" in cat_outputs:
        import pandas as pd
        gp_probs = cat_outputs["global_preservation"].numpy().squeeze(-1)
        pair_df = pd.DataFrame({
            "gene_id": all_gene_ids,
            "weak_label": y_true,
        })
        wgra = within_gene_ranking_accuracy(pair_df, gp_probs)
        print(f"\n  Within-gene ranking accuracy: {wgra:.4f}")

    # Deviation class metrics
    dev_metrics = None
    if "deviation_class" in cat_outputs and len(dev_true) > 0:
        dev_probs = cat_outputs["deviation_class"].numpy()
        dev_preds = dev_probs.argmax(axis=-1)
        valid = dev_true != -1
        if valid.sum() > 0:
            from sklearn.metrics import balanced_accuracy_score, accuracy_score
            try:
                bal_acc = balanced_accuracy_score(dev_true[valid], dev_preds[valid])
                acc = accuracy_score(dev_true[valid], dev_preds[valid])
                dev_metrics = {"balanced_accuracy": bal_acc, "accuracy": acc, "n": int(valid.sum())}
                print(f"\n  deviation_class: n={int(valid.sum())}  "
                      f"BalAcc={bal_acc:.4f}  Acc={acc:.4f}")
            except Exception:
                pass

    # Write markdown report
    report_lines = [
        f"# Proteus Evaluation Report",
        f"",
        f"**Checkpoint:** `{args.checkpoint}`",
        f"**Split:** {args.split}",
        f"**Total samples:** {len(dataset)}",
        f"",
        f"## Per-Task Metrics",
        f"",
        f"| Task | N | AUROC | AUPRC | Brier | Balanced Acc |",
        f"|------|---|-------|-------|-------|--------------|",
    ]

    for task, m in results.items():
        auroc = f"{m['auroc']:.4f}" if m['auroc'] is not None else "N/A"
        auprc = f"{m['auprc']:.4f}" if m['auprc'] is not None else "N/A"
        brier = f"{m['brier']:.4f}" if m['brier'] is not None else "N/A"
        bal_acc = f"{m['balanced_accuracy']:.4f}" if m['balanced_accuracy'] is not None else "N/A"
        report_lines.append(
            f"| {task} | {m['n']} | {auroc} | {auprc} | {brier} | {bal_acc} |"
        )

    if wgra is not None:
        report_lines += [
            f"",
            f"## Within-Gene Ranking Accuracy",
            f"",
            f"**{wgra:.4f}** — fraction of multi-isoform genes where the reference-"
            f"preferred isoform is ranked highest by global_preservation probability.",
        ]

    if dev_metrics is not None:
        report_lines += [
            f"",
            f"## Deviation Class Classification",
            f"",
            f"| Metric | Value |",
            f"|--------|-------|",
            f"| N | {dev_metrics['n']} |",
            f"| Balanced Accuracy | {dev_metrics['balanced_accuracy']:.4f} |",
            f"| Accuracy | {dev_metrics['accuracy']:.4f} |",
            f"",
            f"**Classes:**",
        ]
        for idx, name in DEVIATION_CLASSES.items():
            report_lines.append(f"- {idx}: `{name}`")

    report_lines += [
        f"",
        f"---",
        f"*Generated by run_evaluate.py*",
    ]

    report_text = "\n".join(report_lines)
    args.output_report.parent.mkdir(parents=True, exist_ok=True)
    args.output_report.write_text(report_text)
    print(f"\nEvaluation report saved to: {args.output_report}")


if __name__ == "__main__":
    main()
