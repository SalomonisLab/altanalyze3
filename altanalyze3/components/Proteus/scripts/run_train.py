"""
run_train.py: Fine-tune ProteusModel on Daedalus Phase A isoform pairs.

Loads all configs, creates datasets and DataLoaders, instantiates ProteusTrainer,
runs training, and reports final metrics.
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import numpy as np
import torch
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fine-tune ProteusModel on Daedalus Phase A isoform pairs."
    )
    parser.add_argument(
        "--model_config",
        type=Path,
        default=Path("configs/model.yaml"),
        help="Path to model config. Default: configs/model.yaml",
    )
    parser.add_argument(
        "--training_config",
        type=Path,
        default=Path("configs/training.yaml"),
        help="Path to training config. Default: configs/training.yaml",
    )
    parser.add_argument(
        "--data_processed_dir",
        type=Path,
        default=Path("data/processed"),
        help="Directory containing proteus_train/val/test.tsv.",
    )
    parser.add_argument(
        "--cache_dir",
        type=Path,
        default=Path("data/interim"),
        help="Root cache directory for embeddings.",
    )
    parser.add_argument(
        "--pretrain_checkpoint",
        type=Path,
        default=None,
        help="Path to pretrain checkpoint for weight transfer. Optional.",
    )
    parser.add_argument(
        "--resume_checkpoint",
        type=Path,
        default=None,
        help="Path to a previous training checkpoint to resume from.",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=Path("data/processed"),
        help="Directory for saving checkpoints and logs.",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="auto",
        help="Device (auto/cuda/mps/cpu). Default: auto",
    )
    parser.add_argument(
        "--n_epochs",
        type=int,
        default=None,
        help="Override max_epochs from config.",
    )
    parser.add_argument(
        "--use_wandb",
        action="store_true",
        help="Enable Weights & Biases logging.",
    )
    return parser.parse_args()


def _load_yaml(path: Path) -> dict:
    if path.exists():
        with open(path, "r") as f:
            return yaml.safe_load(f) or {}
    return {}


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def main() -> None:
    args = parse_args()
    sys.path.insert(0, str(Path(__file__).parent.parent))

    from proteus.model import ProteusModel
    from proteus.data.dataset import ProteusDataset
    from proteus.data.collate import proteus_collate_fn
    from proteus.training.trainer import ProteusTrainer
    from proteus.utils.device import get_device, get_device_info

    # Load configs
    model_cfg = _load_yaml(args.model_config)
    train_cfg = _load_yaml(args.training_config)

    # Merge configs for trainer
    combined_cfg = {**model_cfg, **train_cfg}

    seed = train_cfg.get("seed", 42)
    set_seed(seed)

    device = str(get_device(args.device))
    print(f"Device: {device}")
    print(f"Device info: {get_device_info()}")

    # WandB
    if args.use_wandb:
        try:
            import wandb
            wandb.init(
                project="proteus",
                config={**model_cfg, **train_cfg},
            )
        except ImportError:
            print("wandb not installed. Skipping.")

    # Datasets
    train_tsv = args.data_processed_dir / "proteus_train.tsv"
    val_tsv = args.data_processed_dir / "proteus_val.tsv"

    if not train_tsv.exists():
        print(f"ERROR: Training TSV not found: {train_tsv}")
        print("Run init_proteus.py first.")
        sys.exit(1)

    print(f"Loading datasets...")
    train_dataset = ProteusDataset(tsv_path=train_tsv, cache_dir=args.cache_dir)
    val_dataset = ProteusDataset(tsv_path=val_tsv, cache_dir=args.cache_dir)

    batch_size = train_cfg.get("batch_size", 256)
    train_loader = torch.utils.data.DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=4,
        collate_fn=proteus_collate_fn,
        pin_memory=(device == "cuda"),
        drop_last=True,
    )
    val_loader = torch.utils.data.DataLoader(
        val_dataset,
        batch_size=batch_size * 2,
        shuffle=False,
        num_workers=2,
        collate_fn=proteus_collate_fn,
        pin_memory=(device == "cuda"),
    )

    print(f"Train batches: {len(train_loader)}, Val batches: {len(val_loader)}")

    # Model
    if args.pretrain_checkpoint is not None and args.pretrain_checkpoint.exists():
        print(f"Loading pretrained weights from: {args.pretrain_checkpoint}")
        model = ProteusModel.from_pretrain(
            args.pretrain_checkpoint,
            cfg=model_cfg,
        )
    else:
        model = ProteusModel(cfg=model_cfg)

    print(f"Model parameters: {model.count_parameters():,}")

    # Trainer
    trainer = ProteusTrainer(
        model=model,
        train_loader=train_loader,
        val_loader=val_loader,
        cfg=combined_cfg,
        device=device,
    )

    # Optionally resume from checkpoint
    if args.resume_checkpoint is not None and args.resume_checkpoint.exists():
        print(f"Resuming from checkpoint: {args.resume_checkpoint}")
        trainer.load_checkpoint(args.resume_checkpoint)

    # Train
    n_epochs = args.n_epochs or train_cfg.get("max_epochs", 100)
    print(f"\nStarting training for {n_epochs} epochs...")
    history = trainer.train(n_epochs=n_epochs)

    # Final metrics
    if history:
        last = history[-1]
        print("\nFinal epoch metrics:")
        for k, v in sorted(last.items()):
            if isinstance(v, float):
                print(f"  {k}: {v:.4f}")

    best_ckpt = args.output_dir / "best_model.pt"
    print(f"\nTraining complete. Best model saved to: {best_ckpt}")
    print(f"Run evaluation with:")
    print(f"  python scripts/run_evaluate.py --checkpoint {best_ckpt}")


if __name__ == "__main__":
    main()
