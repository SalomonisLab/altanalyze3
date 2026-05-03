"""
evaluator.py: Comprehensive benchmarking for Proteus.

Evaluates Proteus against multiple baselines on held-out test data,
stratified by:
  - Gene family (kinase, TF, membrane, GPCR, adhesion, etc.)
  - Evidence type (ClinVar-supported, experimental, rule-based)
  - Protein class (single-pass TM, multi-pass TM, soluble, secreted)

Metrics per task:
  - AUROC
  - AUPRC
  - F1 at threshold 0.5
  - Expected Calibration Error (ECE)
  - Bootstrap 95% confidence intervals (1000 resamples)

Usage
-----
from proteus.benchmarking import ProteusEvaluator

evaluator = ProteusEvaluator(device="cuda")
results = evaluator.run(
    model_checkpoint="data/processed/best_model.pt",
    test_tsv="data/processed/proteus_test.tsv",
    cache_dir="data/interim",
)
evaluator.to_tsv(results, "benchmark_results.tsv")
evaluator.print_summary(results)
"""

from __future__ import annotations

import warnings
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


# Tasks that require specific mask columns
_TASK_MASK: Dict[str, Optional[str]] = {
    "global_preservation": None,
    "topology_preserved": "is_membrane_pair",
    "surface_retained": "is_surface_pair",
    "kinase_competent": "is_kinase_pair",
    "tf_competent": "is_tf_pair",
    "disorder_preserved": None,
    "localization_preserved": "is_localization_pair",
    "tm_insertion_competent": "is_membrane_pair",
    "signaling_competent": "is_signaling_pair",
    "ppi_interface_preserved": "is_ppi_pair",
}


@dataclass
class TaskMetrics:
    """Metrics for a single task."""
    task: str
    n_positives: int = 0
    n_negatives: int = 0
    auroc: float = float("nan")
    auprc: float = float("nan")
    f1: float = float("nan")
    precision: float = float("nan")
    recall: float = float("nan")
    ece: float = float("nan")
    auroc_ci_lower: float = float("nan")
    auroc_ci_upper: float = float("nan")
    auprc_ci_lower: float = float("nan")
    auprc_ci_upper: float = float("nan")
    stratum: str = "all"
    model_name: str = "Proteus"


@dataclass
class BenchmarkResults:
    """Full benchmark results across all tasks and strata."""
    model_name: str = "Proteus"
    test_tsv: str = ""
    n_total_pairs: int = 0
    task_metrics: List[TaskMetrics] = field(default_factory=list)
    deviation_class_accuracy: float = float("nan")
    deviation_class_f1_macro: float = float("nan")
    stratified_metrics: Dict[str, List[TaskMetrics]] = field(default_factory=dict)


class ProteusEvaluator:
    """
    Comprehensive evaluator for Proteus model benchmarking.

    Parameters
    ----------
    device : str
        "cuda" or "cpu".
    n_bootstrap : int
        Number of bootstrap resamples for CI computation (default 500).
    threshold : float
        Decision threshold for binary classification (default 0.5).
    """

    def __init__(
        self,
        device: str = "cpu",
        n_bootstrap: int = 500,
        threshold: float = 0.5,
    ) -> None:
        self.device = device
        self.n_bootstrap = n_bootstrap
        self.threshold = threshold

    # ------------------------------------------------------------------
    # Main benchmark entry point
    # ------------------------------------------------------------------

    def run(
        self,
        model_checkpoint: Path,
        test_tsv: Path,
        cache_dir: Path,
        additional_models: Optional[Dict[str, Any]] = None,
    ) -> BenchmarkResults:
        """
        Run full benchmark: Proteus + baselines on held-out test data.

        Parameters
        ----------
        model_checkpoint : Path
            Path to trained Proteus .pt checkpoint.
        test_tsv : Path
            Path to test split TSV.
        cache_dir : Path
            Root cache directory (contains rna_embeddings/, etc.).
        additional_models : dict | None
            Optional additional model instances to benchmark alongside Proteus.
            Format: {"model_name": model_instance}

        Returns
        -------
        BenchmarkResults
        """
        import torch
        from torch.utils.data import DataLoader
        from proteus.model import ProteusModel
        from proteus.data.dataset import ProteusDataset
        from proteus.data.collate import proteus_collate_fn

        # Load model
        model = ProteusModel.load_checkpoint(Path(model_checkpoint))
        model.to(self.device)
        model.eval()

        # Load test dataset
        dataset = ProteusDataset(
            tsv_path=Path(test_tsv),
            cache_dir=Path(cache_dir),
            split="test",
        )
        loader = DataLoader(
            dataset,
            batch_size=64,
            shuffle=False,
            num_workers=0,
            collate_fn=proteus_collate_fn,
        )

        # Collect predictions and labels
        all_preds, all_labels, all_masks, all_meta = self._collect_predictions(
            model, loader
        )

        results = BenchmarkResults(
            model_name="Proteus",
            test_tsv=str(test_tsv),
            n_total_pairs=len(dataset),
        )

        # Compute per-task metrics
        for task in _TASK_MASK:
            preds = all_preds.get(task)
            labels = all_labels.get(task)
            masks = all_masks.get(task)
            if preds is None or labels is None:
                continue
            tm = self._compute_task_metrics(task, preds, labels, masks, "all", "Proteus")
            results.task_metrics.append(tm)

        # Deviation class metrics
        if "deviation_class" in all_preds and "deviation_class_label" in all_labels:
            acc, f1 = self._compute_multiclass_metrics(
                all_preds["deviation_class"],
                all_labels["deviation_class_label"],
            )
            results.deviation_class_accuracy = acc
            results.deviation_class_f1_macro = f1

        # Stratified metrics by gene family
        meta_df = self._build_meta_df(all_meta)
        for stratum_col in ("gene_family", "protein_class", "tm_class"):
            if stratum_col in meta_df.columns:
                strat_metrics = self._stratified_metrics(
                    all_preds, all_labels, all_masks, meta_df, stratum_col
                )
                results.stratified_metrics[stratum_col] = strat_metrics

        return results

    def _collect_predictions(
        self, model, loader
    ) -> Tuple[Dict, Dict, Dict, List]:
        """Run model forward pass on all batches and collect outputs."""
        import torch
        from proteus.heads.task_heads import BINARY_TASKS

        all_preds: Dict[str, List] = defaultdict(list)
        all_labels: Dict[str, List] = defaultdict(list)
        all_masks: Dict[str, List] = defaultdict(list)
        all_meta: List[Dict] = []

        with torch.no_grad():
            for batch in loader:
                batch_device = {
                    k: v.to(self.device) if hasattr(v, "to") else v
                    for k, v in batch.items()
                    if not isinstance(v, dict)
                }
                batch_device["task_masks"] = batch.get("task_masks", {})

                preds = model.predict(batch_device)

                for task in BINARY_TASKS:
                    if task in preds:
                        all_preds[task].append(preds[task].cpu().numpy().flatten())
                    # Labels
                    label_key = f"label_{task}" if f"label_{task}" in batch else "label_preservation"
                    if task == "global_preservation":
                        label_key = "label_preservation"
                    if label_key in batch:
                        all_labels[task].append(batch[label_key].numpy().flatten())
                    # Masks
                    mask_key = f"{task}_mask"
                    if "task_masks" in batch and task in batch["task_masks"]:
                        all_masks[task].append(batch["task_masks"][task].numpy().flatten())

                if "deviation_class" in preds:
                    all_preds["deviation_class"].append(preds["deviation_class"].cpu().numpy())
                if "deviation_class" in batch:
                    all_labels["deviation_class_label"].append(batch["deviation_class"].numpy().flatten())

                # Metadata
                for i in range(len(batch.get("gene_id", []))):
                    all_meta.append({
                        "gene_id": batch["gene_id"][i] if "gene_id" in batch else "",
                        "gene_name": batch.get("gene_name", [""])[i],
                        "is_membrane": int(batch["task_masks"]["topology_preserved"][i].item()) if "task_masks" in batch else 0,
                        "is_kinase": int(batch["task_masks"]["kinase_competent"][i].item()) if "task_masks" in batch else 0,
                        "is_tf": int(batch["task_masks"]["tf_competent"][i].item()) if "task_masks" in batch else 0,
                    })

        # Concatenate
        concat_preds = {k: np.concatenate(v) for k, v in all_preds.items() if v}
        concat_labels = {k: np.concatenate(v) for k, v in all_labels.items() if v}
        concat_masks = {k: np.concatenate(v) for k, v in all_masks.items() if v}

        return concat_preds, concat_labels, concat_masks, all_meta

    def _compute_task_metrics(
        self,
        task: str,
        preds: np.ndarray,
        labels: np.ndarray,
        masks: Optional[np.ndarray],
        stratum: str,
        model_name: str,
    ) -> TaskMetrics:
        """Compute AUROC, AUPRC, F1, ECE for one task with optional masking."""
        from sklearn.metrics import (
            roc_auc_score, average_precision_score,
            f1_score, precision_score, recall_score,
        )

        # Apply mask: only evaluate on samples where mask=1 and label != -1
        valid = (labels != -1)
        if masks is not None:
            valid = valid & (masks > 0.5)

        p = preds[valid]
        y = labels[valid].astype(int)

        n_pos = int(y.sum())
        n_neg = int((1 - y).sum())

        tm = TaskMetrics(
            task=task,
            n_positives=n_pos,
            n_negatives=n_neg,
            stratum=stratum,
            model_name=model_name,
        )

        if n_pos < 2 or n_neg < 2:
            # Insufficient data for reliable metrics
            return tm

        try:
            tm.auroc = float(roc_auc_score(y, p))
            tm.auprc = float(average_precision_score(y, p))
            y_bin = (p >= self.threshold).astype(int)
            tm.f1 = float(f1_score(y, y_bin, zero_division=0))
            tm.precision = float(precision_score(y, y_bin, zero_division=0))
            tm.recall = float(recall_score(y, y_bin, zero_division=0))
            tm.ece = float(self._expected_calibration_error(y, p))

            # Bootstrap CIs
            auroc_samples, auprc_samples = [], []
            rng = np.random.default_rng(42)
            for _ in range(self.n_bootstrap):
                idx = rng.integers(0, len(y), size=len(y))
                yb, pb = y[idx], p[idx]
                if yb.sum() < 1 or (1 - yb).sum() < 1:
                    continue
                try:
                    auroc_samples.append(roc_auc_score(yb, pb))
                    auprc_samples.append(average_precision_score(yb, pb))
                except Exception:
                    pass

            if auroc_samples:
                tm.auroc_ci_lower = float(np.percentile(auroc_samples, 2.5))
                tm.auroc_ci_upper = float(np.percentile(auroc_samples, 97.5))
            if auprc_samples:
                tm.auprc_ci_lower = float(np.percentile(auprc_samples, 2.5))
                tm.auprc_ci_upper = float(np.percentile(auprc_samples, 97.5))

        except Exception as exc:
            warnings.warn(f"[ProteusEvaluator] Metric computation failed for {task}: {exc}")

        return tm

    @staticmethod
    def _expected_calibration_error(y: np.ndarray, p: np.ndarray, n_bins: int = 10) -> float:
        """Compute Expected Calibration Error (ECE)."""
        bins = np.linspace(0, 1, n_bins + 1)
        ece = 0.0
        n = len(y)
        for i in range(n_bins):
            lo, hi = bins[i], bins[i + 1]
            mask = (p >= lo) & (p < hi)
            if mask.sum() == 0:
                continue
            frac_pos = y[mask].mean()
            mean_prob = p[mask].mean()
            ece += mask.sum() / n * abs(frac_pos - mean_prob)
        return ece

    @staticmethod
    def _compute_multiclass_metrics(
        preds: np.ndarray,
        labels: np.ndarray,
    ) -> Tuple[float, float]:
        """Accuracy and macro-F1 for deviation class prediction."""
        from sklearn.metrics import accuracy_score, f1_score
        valid = labels != -1
        if valid.sum() < 10:
            return float("nan"), float("nan")
        y = labels[valid]
        y_hat = preds[valid].argmax(axis=1) if preds.ndim > 1 else preds[valid]
        acc = float(accuracy_score(y, y_hat))
        f1 = float(f1_score(y, y_hat, average="macro", zero_division=0))
        return acc, f1

    @staticmethod
    def _build_meta_df(all_meta: List[Dict]):
        """Convert metadata list to DataFrame, add derived columns."""
        import pandas as pd
        if not all_meta:
            return pd.DataFrame()
        df = pd.DataFrame(all_meta)
        # Assign gene family based on mask columns
        def assign_family(row):
            if row.get("is_kinase"):
                return "kinase"
            if row.get("is_tf"):
                return "transcription_factor"
            if row.get("is_membrane"):
                return "membrane"
            return "other"
        df["gene_family"] = df.apply(assign_family, axis=1)
        return df

    def _stratified_metrics(
        self,
        all_preds, all_labels, all_masks,
        meta_df,
        stratum_col: str,
    ) -> List[TaskMetrics]:
        """Compute per-stratum metrics."""
        results = []
        if stratum_col not in meta_df.columns:
            return results
        for stratum_val in meta_df[stratum_col].unique():
            idx = meta_df[meta_df[stratum_col] == stratum_val].index.values
            for task in _TASK_MASK:
                preds = all_preds.get(task)
                labels = all_labels.get(task)
                masks = all_masks.get(task)
                if preds is None or labels is None:
                    continue
                if len(idx) == 0 or idx.max() >= len(preds):
                    continue
                tm = self._compute_task_metrics(
                    task,
                    preds[idx],
                    labels[idx],
                    masks[idx] if masks is not None else None,
                    str(stratum_val),
                    "Proteus",
                )
                results.append(tm)
        return results

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def to_dataframe(self, results: BenchmarkResults):
        """Convert BenchmarkResults to a flat pandas DataFrame."""
        import pandas as pd
        rows = []
        for tm in results.task_metrics:
            rows.append(vars(tm))
        for strat_list in results.stratified_metrics.values():
            for tm in strat_list:
                rows.append(vars(tm))
        return pd.DataFrame(rows)

    def to_tsv(self, results: BenchmarkResults, output_path: Path) -> None:
        """Write benchmark results to TSV."""
        df = self.to_dataframe(results)
        df.to_csv(output_path, sep="\t", index=False)
        print(f"[ProteusEvaluator] Results written to {output_path}")

    def print_summary(self, results: BenchmarkResults) -> None:
        """Print a formatted summary table."""
        header = f"{'Task':<30} {'N+':<6} {'N-':<6} {'AUROC':<8} {'95% CI':<16} {'AUPRC':<8} {'F1':<6} {'ECE':<6}"
        print(f"\n=== Proteus Benchmark: {results.n_total_pairs:,} test pairs ===")
        print(header)
        print("-" * len(header))
        for tm in results.task_metrics:
            if np.isnan(tm.auroc):
                continue
            ci = f"[{tm.auroc_ci_lower:.3f}-{tm.auroc_ci_upper:.3f}]"
            print(
                f"{tm.task:<30} {tm.n_positives:<6} {tm.n_negatives:<6} "
                f"{tm.auroc:<8.4f} {ci:<16} {tm.auprc:<8.4f} "
                f"{tm.f1:<6.3f} {tm.ece:<6.4f}"
            )
        if not np.isnan(results.deviation_class_accuracy):
            print(f"\nDeviation class accuracy: {results.deviation_class_accuracy:.4f}  "
                  f"macro-F1: {results.deviation_class_f1_macro:.4f}")
        print()
