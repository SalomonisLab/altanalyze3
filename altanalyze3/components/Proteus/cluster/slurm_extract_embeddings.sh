#!/bin/bash
#SBATCH --job-name=proteus_extract
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --output=logs/extract_%j.out
#SBATCH --error=logs/extract_%j.err

# ---------------------------------------------------------------
# Proteus: GPU embedding extraction (RNA via Orthrus, Protein via ESM2)
#
# CPU-only steps (disorder, structural, topology) do NOT need this
# job — run them directly on a login/compute node without a GPU.
#
# Usage:
#   sbatch cluster/slurm_extract_embeddings.sh
# ---------------------------------------------------------------

DATA_DIR="/scratch/${USER}/proteus_data"
SIF_PATH="/scratch/${USER}/proteus.sif"
GENCODE_DIR="${DATA_DIR}/raw/gencode"

mkdir -p "${DATA_DIR}/logs"

echo "Job ID: ${SLURM_JOB_ID}"
nvidia-smi

# RNA embeddings (Orthrus, ~3-4 hours for full GENCODE v49)
echo "=== Extracting RNA embeddings ==="
singularity exec --nv \
    --bind "${DATA_DIR}:/workspace/Proteus/data" \
    "${SIF_PATH}" \
    python /workspace/Proteus/scripts/extract_rna_embeddings.py \
        --fasta data/raw/gencode/gencode.v49.transcripts.fa.gz \
        --output_dir data/interim/rna_embeddings \
        --batch_size 32

# Protein embeddings (ESM2, ~2-3 hours for full GENCODE v49)
echo "=== Extracting protein embeddings ==="
singularity exec --nv \
    --bind "${DATA_DIR}:/workspace/Proteus/data" \
    "${SIF_PATH}" \
    python /workspace/Proteus/scripts/extract_protein_embeddings.py \
        --fasta data/raw/gencode/gencode.v49.pc_translations.fa.gz \
        --output_dir data/interim/protein_embeddings \
        --batch_size 16

echo "Embedding extraction complete."
