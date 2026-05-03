"""
run_pretrain.py: Self-supervised pretraining of ProteusPretrainModel on
synthetic exon-skip events.

Loads configs/pretrain.yaml and configs/model.yaml, instantiates model,
runs training loop, saves best checkpoint to data/processed/pretrain_best.pt.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import torch
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Proteus SSL pretraining on synthetic exon-skip pairs."
    )
    parser.add_argument(
        "--pretrain_config",
        type=Path,
        default=Path("configs/pretrain.yaml"),
        help="Path to pretrain config. Default: configs/pretrain.yaml",
    )
    parser.add_argument(
        "--model_config",
        type=Path,
        default=Path("configs/model.yaml"),
        help="Path to model config. Default: configs/model.yaml",
    )
    parser.add_argument(
        "--synthetic_tsv",
        type=Path,
        default=Path("data/interim/pretrain_events/synthetic_exon_skips.tsv"),
        help="Path to synthetic exon-skip TSV.",
    )
    parser.add_argument(
        "--cache_dir",
        type=Path,
        default=Path("data/interim"),
        help="Root cache directory with rna_embeddings/, protein_embeddings/, etc.",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=Path("data/processed"),
        help="Directory for saving checkpoints.",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="auto",
        help="Device (auto/cuda/mps/cpu). Default: auto",
    )
    parser.add_argument(
        "--max_samples",
        type=int,
        default=None,
        help="Limit dataset to N samples (for debugging).",
    )
    return parser.parse_args()


def _load_yaml(path: Path) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


def main() -> None:
    args = parse_args()
    sys.path.insert(0, str(Path(__file__).parent.parent))

    from proteus.pretrain_model import ProteusPretrainModel
    from proteus.data.pretrain_dataset import SyntheticExonSkipDataset
    from proteus.data.collate import proteus_collate_fn
    from proteus.pretraining.objectives import PretrainObjective
    from proteus.utils.device import get_device

    # Load configs
    pretrain_cfg = _load_yaml(args.pretrain_config) if args.pretrain_config.exists() else {}
    model_cfg = _load_yaml(args.model_config) if args.model_config.exists() else {}

    device = get_device(args.device)
    print(f"Device: {device}")

    # Check synthetic TSV exists
    if not args.synthetic_tsv.exists():
        print(f"ERROR: Synthetic TSV not found: {args.synthetic_tsv}")
        print("Run generate_pretrain_events.py first.")
        sys.exit(1)

    # Dataset
    dataset = SyntheticExonSkipDataset(
        tsv_path=args.synthetic_tsv,
        cache_dir=args.cache_dir,
        max_samples=args.max_samples,
    )

    batch_size = pretrain_cfg.get("batch_size", 512)
    loader = torch.utils.data.DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=2,
        collate_fn=proteus_collate_fn,
        pin_memory=(str(device) == "cuda"),
    )

    # Model
    model = ProteusPretrainModel(cfg=model_cfg)
    model.to(device)
    print(f"Model parameters: {model.count_parameters():,}")

    # Loss
    task_weights = {}
    for task, tcfg in pretrain_cfg.get("tasks", {}).items():
        if isinstance(tcfg, dict):
            task_weights[task] = tcfg.get("loss_weight", 1.0)
    objective = PretrainObjective(loss_weights=task_weights or None)

    # Optimizer
    lr = pretrain_cfg.get("lr", 1e-3)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)

    n_epochs = pretrain_cfg.get("epochs", 30)
    total_steps = n_epochs * len(loader)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=total_steps)

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    best_checkpoint = output_dir / "pretrain_best.pt"

    use_amp = str(device) == "cuda"
    scaler = torch.cuda.amp.GradScaler(enabled=use_amp)

    best_loss = float("inf")

    print(f"\nStarting pretraining for {n_epochs} epochs ({len(loader)} batches/epoch)...")

    for epoch in range(1, n_epochs + 1):
        model.train()
        epoch_losses: dict = {}
        t0 = time.time()

        for batch in loader:
            # Move to device
            for k, v in batch.items():
                if isinstance(v, torch.Tensor):
                    batch[k] = v.to(device)

            optimizer.zero_grad()

            with torch.autocast(device_type="cuda" if str(device) == "cuda" else "cpu",
                                enabled=use_amp):
                outputs = model(batch)
                ssl_labels = {
                    "domain_disrupted": batch.get("domain_disrupted"),
                    "nmd_triggered": batch.get("nmd_triggered"),
                    "tm_lost": batch.get("tm_lost"),
                    "disorder_delta": batch.get("disorder_delta"),
                }
                ssl_labels = {k: v for k, v in ssl_labels.items() if v is not None}
                losses = objective.compute(outputs, ssl_labels)

            total_loss = losses["total"]
            if use_amp:
                scaler.scale(total_loss).backward()
                scaler.step(optimizer)
                scaler.update()
            else:
                total_loss.backward()
                optimizer.step()
            scheduler.step()

            for k, v in losses.items():
                if isinstance(v, torch.Tensor):
                    epoch_losses.setdefault(k, []).append(v.item())

        elapsed = time.time() - t0
        mean_losses = {k: float(np.mean(vs)) for k, vs in epoch_losses.items()}
        total_mean = mean_losses.get("total", 0)

        print(
            f"Epoch {epoch:3d}/{n_epochs} [{elapsed:.1f}s] "
            + " ".join(f"{k}={v:.4f}" for k, v in mean_losses.items())
        )

        if total_mean < best_loss:
            best_loss = total_mean
            model.save_checkpoint(best_checkpoint, epoch=epoch, metrics=mean_losses)
            print(f"  New best checkpoint saved (loss={best_loss:.4f})")

    print(f"\nPretraining complete. Best checkpoint: {best_checkpoint}")
    print(f"Transfer weights to ProteusModel with:")
    print(f"  model = ProteusModel.from_pretrain('{best_checkpoint}')")


if __name__ == "__main__":
    main()
