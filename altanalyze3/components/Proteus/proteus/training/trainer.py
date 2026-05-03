"""
ProteusTrainer: Plain PyTorch training loop for Proteus.

No PyTorch Lightning. Includes cosine LR with linear warmup,
mixed precision, early stopping, and comprehensive logging.
"""

from __future__ import annotations

import math
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from ..model import ProteusModel
from .losses import ProteusLoss
from .metrics import compute_task_metrics, within_gene_ranking_accuracy


class WarmupCosineScheduler:
    """Linear warmup followed by cosine annealing."""

    def __init__(
        self,
        optimizer: torch.optim.Optimizer,
        warmup_steps: int,
        total_steps: int,
        min_lr_ratio: float = 0.0,
    ) -> None:
        self.optimizer = optimizer
        self.warmup_steps = warmup_steps
        self.total_steps = total_steps
        self.min_lr_ratio = min_lr_ratio
        self._step = 0
        self._base_lrs = [pg["lr"] for pg in optimizer.param_groups]

    def step(self) -> None:
        self._step += 1
        lr_scale = self._get_lr_scale()
        for pg, base_lr in zip(self.optimizer.param_groups, self._base_lrs):
            pg["lr"] = base_lr * lr_scale

    def _get_lr_scale(self) -> float:
        if self._step <= self.warmup_steps:
            return float(self._step) / max(1, self.warmup_steps)
        progress = float(self._step - self.warmup_steps) / max(
            1, self.total_steps - self.warmup_steps
        )
        cosine_scale = 0.5 * (1.0 + math.cos(math.pi * progress))
        return self.min_lr_ratio + (1.0 - self.min_lr_ratio) * cosine_scale

    def state_dict(self) -> dict:
        return {"_step": self._step, "_base_lrs": self._base_lrs}

    def load_state_dict(self, sd: dict) -> None:
        self._step = sd["_step"]
        self._base_lrs = sd["_base_lrs"]


class ProteusTrainer:
    """
    Training orchestrator for ProteusModel.

    Handles the full training loop including:
    - AdamW optimizer
    - Cosine LR scheduler with linear warmup
    - Mixed precision (torch.autocast)
    - Gradient clipping
    - Early stopping
    - Per-task metric logging
    - Checkpoint saving

    Parameters
    ----------
    model : ProteusModel
    train_loader : DataLoader
    val_loader : DataLoader
    cfg : dict
        Training configuration (from configs/training.yaml).
    device : str
        Device string ("cuda", "cpu", "mps").
    """

    def __init__(
        self,
        model: ProteusModel,
        train_loader: DataLoader,
        val_loader: DataLoader,
        cfg: Dict[str, Any],
        device: str = "cpu",
    ) -> None:
        self.model = model.to(device)
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.cfg = cfg
        self.device = device

        # Loss
        task_weights = {}
        for task_name, task_cfg in cfg.get("heads", {}).items():
            if isinstance(task_cfg, dict):
                task_weights[task_name] = task_cfg.get("loss_weight", 1.0)
        self.loss_fn = ProteusLoss(task_weights=task_weights or None)

        # Optimizer
        self.optimizer = torch.optim.AdamW(
            model.parameters(),
            lr=cfg.get("lr", 3e-4),
            weight_decay=cfg.get("weight_decay", 1e-4),
        )

        # LR scheduler
        steps_per_epoch = len(train_loader)
        total_steps = steps_per_epoch * cfg.get("max_epochs", 100)
        warmup_steps = cfg.get("warmup_steps", 500)
        self.scheduler = WarmupCosineScheduler(
            self.optimizer,
            warmup_steps=warmup_steps,
            total_steps=total_steps,
        )

        # Mixed precision
        self.use_amp = cfg.get("mixed_precision", True) and device == "cuda"
        self.scaler = torch.cuda.amp.GradScaler(enabled=self.use_amp)

        # Early stopping
        self.patience = cfg.get("early_stopping_patience", 10)
        self.es_metric = cfg.get("early_stopping_metric", "val/global_preservation/auroc")
        self.best_metric = -float("inf")
        self.patience_counter = 0

        # History
        self.history: List[Dict[str, Any]] = []

        # wandb
        self._wandb = None
        try:
            import wandb
            self._wandb = wandb
        except ImportError:
            pass

    def _to_device(self, batch: Dict[str, Any]) -> Dict[str, Any]:
        """Move batch tensors to device."""
        result = {}
        for k, v in batch.items():
            if isinstance(v, torch.Tensor):
                result[k] = v.to(self.device)
            elif isinstance(v, dict):
                result[k] = {
                    kk: vv.to(self.device) if isinstance(vv, torch.Tensor) else vv
                    for kk, vv in v.items()
                }
            else:
                result[k] = v
        return result

    def train_epoch(self) -> Dict[str, float]:
        """Run one training epoch."""
        self.model.train()
        epoch_losses: Dict[str, List[float]] = {}
        n_batches = 0

        for batch in self.train_loader:
            batch = self._to_device(batch)

            self.optimizer.zero_grad()

            with torch.autocast(
                device_type=self.device if self.device != "mps" else "cpu",
                enabled=self.use_amp,
            ):
                outputs = self.model(batch)
                losses = self.loss_fn(outputs, batch)

            total_loss = losses["total"]

            if self.use_amp:
                self.scaler.scale(total_loss).backward()
                self.scaler.unscale_(self.optimizer)
                torch.nn.utils.clip_grad_norm_(
                    self.model.parameters(),
                    max_norm=self.cfg.get("grad_clip", 1.0),
                )
                self.scaler.step(self.optimizer)
                self.scaler.update()
            else:
                total_loss.backward()
                torch.nn.utils.clip_grad_norm_(
                    self.model.parameters(),
                    max_norm=self.cfg.get("grad_clip", 1.0),
                )
                self.optimizer.step()

            self.scheduler.step()
            n_batches += 1

            for k, v in losses.items():
                if isinstance(v, torch.Tensor):
                    epoch_losses.setdefault(k, []).append(v.item())

        return {k: float(np.mean(vs)) for k, vs in epoch_losses.items()}

    def val_epoch(self) -> Dict[str, Any]:
        """Run one validation epoch, computing full metrics."""
        self.model.eval()

        all_outputs: Dict[str, List] = {}
        all_labels_preservation: List = []
        all_deviation_labels: List = []
        all_gene_ids: List = []
        all_task_masks: Dict[str, List] = {}
        val_losses: Dict[str, List[float]] = {}

        with torch.no_grad():
            for batch in self.val_loader:
                batch = self._to_device(batch)
                outputs = self.model(batch)
                losses = self.loss_fn(outputs, batch)

                for k, v in losses.items():
                    if isinstance(v, torch.Tensor):
                        val_losses.setdefault(k, []).append(v.item())

                # Collect predictions
                for task in list(outputs.keys()):
                    if not task.endswith("_mask") and isinstance(outputs[task], torch.Tensor):
                        all_outputs.setdefault(task, []).append(
                            outputs[task].detach().cpu()
                        )

                lbl = batch.get("label_preservation")
                if lbl is not None:
                    all_labels_preservation.append(lbl.cpu())

                dev_lbl = batch.get("deviation_class")
                if dev_lbl is not None:
                    all_deviation_labels.append(dev_lbl.cpu())

                gids = batch.get("gene_id", [])
                all_gene_ids.extend(gids if isinstance(gids, list) else [])

                task_masks = batch.get("task_masks", {})
                for task, mask in task_masks.items():
                    all_task_masks.setdefault(task, []).append(mask.cpu())

        # Concatenate
        concatenated_outputs = {
            k: torch.cat(vs, dim=0) for k, vs in all_outputs.items()
            if vs and isinstance(vs[0], torch.Tensor)
        }
        y_true = (
            torch.cat(all_labels_preservation, dim=0).numpy()
            if all_labels_preservation else np.array([])
        )
        concatenated_masks = {
            k: torch.cat(vs, dim=0).numpy()
            for k, vs in all_task_masks.items()
            if vs
        }

        # Compute metrics
        from ..heads.task_heads import BINARY_TASKS

        metrics: Dict[str, Any] = {
            k: float(np.mean(vs)) for k, vs in val_losses.items()
        }

        for task in BINARY_TASKS:
            if task not in concatenated_outputs:
                continue
            probs = torch.sigmoid(concatenated_outputs[task]).numpy().squeeze(-1)
            task_mask = concatenated_masks.get(task)
            y_true_task = y_true.copy()
            if task_mask is not None:
                y_true_task[task_mask < 0.5] = -1

            task_metrics = compute_task_metrics(y_true_task, probs, task_name=task)
            for k, v in task_metrics.items():
                if isinstance(v, (int, float)) and v is not None:
                    metrics[f"{task}/{k}"] = v

        # Within-gene ranking accuracy
        if all_gene_ids and "global_preservation" in concatenated_outputs:
            import pandas as pd
            gp_probs = (
                torch.sigmoid(concatenated_outputs["global_preservation"])
                .numpy()
                .squeeze(-1)
            )
            pair_df = pd.DataFrame({
                "gene_id": all_gene_ids,
                "weak_label": y_true,
            })
            wgra = within_gene_ranking_accuracy(pair_df, gp_probs)
            metrics["within_gene_ranking_accuracy"] = wgra

        return metrics

    def train(self, n_epochs: Optional[int] = None) -> List[Dict[str, Any]]:
        """
        Main training loop.

        Parameters
        ----------
        n_epochs : int | None
            Number of epochs. Defaults to cfg["max_epochs"].

        Returns
        -------
        list[dict]
            Training history (one dict per epoch).
        """
        n_epochs = n_epochs or self.cfg.get("max_epochs", 100)
        checkpoint_every = self.cfg.get("checkpoint_every_n_epochs", 5)

        best_ckpt_path = Path("data/processed/best_model.pt")

        for epoch in range(1, n_epochs + 1):
            t0 = time.time()
            train_metrics = self.train_epoch()
            val_metrics = self.val_epoch()
            elapsed = time.time() - t0

            epoch_record = {
                "epoch": epoch,
                "elapsed": elapsed,
                **{f"train/{k}": v for k, v in train_metrics.items()},
                **{f"val/{k}": v for k, v in val_metrics.items()},
            }
            self.history.append(epoch_record)

            # Log
            print(
                f"Epoch {epoch:3d}/{n_epochs} [{elapsed:.1f}s] "
                f"train_loss={train_metrics.get('total', 0):.4f} "
                f"val_loss={val_metrics.get('total', 0):.4f} "
                f"gp_auroc={val_metrics.get('global_preservation/auroc', 0):.4f}"
            )

            if self._wandb:
                self._wandb.log(epoch_record, step=epoch)

            # Early stopping check
            es_val = val_metrics.get(
                self.es_metric.replace("val/", ""),
                val_metrics.get("global_preservation/auroc", 0),
            )
            if es_val is None:
                es_val = 0.0

            if es_val > self.best_metric:
                self.best_metric = es_val
                self.patience_counter = 0
                self.save_checkpoint(best_ckpt_path, epoch, val_metrics)
                print(f"  New best model saved (metric={es_val:.4f})")
            else:
                self.patience_counter += 1
                if self.patience_counter >= self.patience:
                    print(f"Early stopping at epoch {epoch} (patience={self.patience})")
                    break

            # Periodic checkpoint
            if epoch % checkpoint_every == 0:
                ckpt_path = Path(f"data/processed/checkpoint_epoch{epoch:04d}.pt")
                self.save_checkpoint(ckpt_path, epoch, val_metrics)

        return self.history

    def save_checkpoint(
        self,
        path: Path,
        epoch: int,
        metrics: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Save full training checkpoint (model + optimizer + scheduler)."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        torch.save(
            {
                "model_state_dict": self.model.state_dict(),
                "optimizer_state_dict": self.optimizer.state_dict(),
                "scheduler_state_dict": self.scheduler.state_dict(),
                "cfg": self.cfg,
                "epoch": epoch,
                "metrics": metrics or {},
                "best_metric": self.best_metric,
                "model_class": "ProteusModel",
            },
            path,
        )

    def load_checkpoint(self, path: Path) -> None:
        """Restore model + optimizer + scheduler from a checkpoint."""
        path = Path(path)
        checkpoint = torch.load(path, map_location=self.device, weights_only=False)
        self.model.load_state_dict(checkpoint["model_state_dict"])
        self.optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
        if "scheduler_state_dict" in checkpoint:
            self.scheduler.load_state_dict(checkpoint["scheduler_state_dict"])
        self.best_metric = checkpoint.get("best_metric", -float("inf"))
        print(
            f"[ProteusTrainer] Restored checkpoint: {path} "
            f"(epoch {checkpoint.get('epoch', '?')})"
        )
