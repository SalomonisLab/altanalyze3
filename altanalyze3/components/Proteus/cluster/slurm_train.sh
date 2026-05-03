#!/bin/bash
#SBATCH --job-name=proteus_train
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/train_%j.out
#SBATCH --error=logs/train_%j.err

# ---------------------------------------------------------------
# Proteus Stage 2: Supervised fine-tuning on Daedalus pairs
#
# Usage:
#   sbatch cluster/slurm_train.sh
#
# Requires pretrain checkpoint from slurm_pretrain.sh.
# ---------------------------------------------------------------

DATA_DIR="/scratch/${USER}/proteus_data"
SIF_PATH="/scratch/${USER}/proteus.sif"
PRETRAIN_CKPT="data/processed/pretrain/best_model.pt"

mkdir -p "${DATA_DIR}/logs"

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
nvidia-smi

singularity exec --nv \
    --bind "${DATA_DIR}:/workspace/Proteus/data" \
    "${SIF_PATH}" \
    python /workspace/Proteus/scripts/run_train.py \
        --pretrain_checkpoint "${PRETRAIN_CKPT}" \
        --config /workspace/Proteus/configs/training.yaml \
        --output_dir data/processed/finetune

echo "Fine-tuning complete."
