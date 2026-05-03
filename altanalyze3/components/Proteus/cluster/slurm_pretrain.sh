#!/bin/bash
#SBATCH --job-name=proteus_pretrain
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/pretrain_%j.out
#SBATCH --error=logs/pretrain_%j.err

# ---------------------------------------------------------------
# Proteus Stage 1: Self-supervised pretraining (synthetic exon skips)
#
# Usage:
#   sbatch cluster/slurm_pretrain.sh
#
# Edit DATA_DIR and SIF_PATH before submitting.
# ---------------------------------------------------------------

DATA_DIR="/scratch/${USER}/proteus_data"
SIF_PATH="/scratch/${USER}/proteus.sif"
PROTEUS_DIR="/scratch/${USER}/Proteus"

mkdir -p "${DATA_DIR}/logs"

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "GPUs:   ${CUDA_VISIBLE_DEVICES}"
nvidia-smi

singularity exec --nv \
    --bind "${DATA_DIR}:/workspace/Proteus/data" \
    "${SIF_PATH}" \
    python /workspace/Proteus/scripts/run_pretrain.py \
        --config /workspace/Proteus/configs/pretrain.yaml \
        --output_dir data/processed/pretrain

echo "Pretraining complete."
