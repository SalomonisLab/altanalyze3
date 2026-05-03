# Proteus Docker image
#
# Base: PyTorch 2.1.2 with CUDA 12.1 (cudnn8 devel — needed to compile mamba_ssm)
#
# Build:
#   docker build -t proteus:latest .
#
# Convert to Singularity SIF (for HPC):
#   singularity build proteus.sif docker-daemon://proteus:latest
#
# Run on GPU:
#   docker run --gpus all -v /host/data:/workspace/data proteus:latest \
#       python scripts/run_train.py
#
# On SLURM with Singularity:
#   singularity exec --nv \
#       --bind /scratch/$USER/data:/workspace/data \
#       proteus.sif python scripts/run_train.py

FROM pytorch/pytorch:2.1.2-cuda12.1-cudnn8-devel

LABEL maintainer="AltAnalyze3/Proteus"
LABEL description="Proteus isoform function prediction model"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        git \
        build-essential \
        wget \
        && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

# Copy requirements first (layer cache — only invalidated if requirements change)
COPY requirements.txt /workspace/requirements.txt

# Install Python dependencies
# mamba-ssm requires CUDA headers (cudnn8-devel base) and --no-build-isolation
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir mamba-ssm --no-build-isolation

# Optional: STARLING disorder predictor (heavier, comment out if not needed)
# RUN pip install --no-cache-dir starling

# Copy Proteus source code
COPY . /workspace/Proteus

# Install Proteus as editable package
RUN pip install --no-cache-dir -e /workspace/Proteus

# Set PYTHONPATH so imports work from any working directory
ENV PYTHONPATH="/workspace/Proteus:${PYTHONPATH}"

# Default working directory for scripts
WORKDIR /workspace/Proteus

# Default command — print available scripts
CMD ["python", "-c", "import proteus; print('Proteus ready. Run: python scripts/run_train.py')"]
