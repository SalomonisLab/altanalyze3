# Proteus: Reference-Conditioned Isoform Function Prediction

> *Proteus, the shape-shifting Greek god, retained essential knowledge even as his form changed. Alternative isoforms may change shape — Proteus determines if they retain function.*

## What is Proteus?

Proteus is a **reference-conditioned delta model** for predicting whether an alternative RNA isoform preserves the core biological function of its canonical reference isoform. Given a pair (reference isoform, alternative isoform), Proteus outputs multi-task predictions covering seven functional dimensions.

It also provides a **surfaceome validation pipeline** that takes novel junction-spanning peptides from AltAnalyze3 long-read output and determines whether they are surface-accessible — the key input for CAR-T epitope discovery and cell-surface proteomics prioritization.

## What Proteus is NOT

Proteus is **not a foundation model trained from scratch**. It does not learn RNA or protein representations from raw sequences. Instead, it:

1. Uses **frozen, pre-trained foundation model encoders** (Orthrus for RNA, ESM2 for protein, STARLING/metapredict for disorder) to extract rich embeddings
2. Learns a **lightweight delta architecture** that reasons about the *difference* between reference and alternative isoform embeddings
3. Is trained on isoform pairs from the **Daedalus Phase A** dataset with weak supervision from biological priors

This means Proteus trains in hours (not months), requires modest GPU resources, and is interpretable at the level of what each modality contributes.

---

## Architecture

```
                    Reference Isoform                   Alternative Isoform
                         |                                       |
          +--------------+----------------+      +---------------+---------------+
          |              |                |      |               |               |
     [Orthrus]       [ESM2]        [metapredict] [Orthrus]   [ESM2]      [metapredict]
    RNA emb [512]  Prot emb [1280]  Dis [50]   RNA emb [512] Prot emb [1280]  Dis [50]
     Struct [96]                                 Struct [96]

          |              |                |          |               |               |
          +-------- ReferenceDelta -------+   +------- ReferenceDelta --------+
          |   (Asymmetric: loss gate,     |   |  + Structural features [96]   |
          |    gain gate, cross-attn)     |   |                               |
          |              |                |   |              |                |
       RNA delta [256]  Prot delta [256]  |   Disorder delta [256] Struct [256]
                                          |
                            CrossModalFusion
                    (cross-attention: RNA <-> Protein, MLP)
                                    |
                             Fused repr [256]
                                    |
                 +------------------+------------------+
                 |         |        |        |         |
           global_    topology_ surface_ kinase_  deviation_
         preservation preserved retained competent  class
            [2]        [2]      [2]      [2]        [6]
              + tf_competent [2] + disorder_preserved [2]
```

### Frozen Encoders (NOT trained)

| Encoder | Model | Output dim | Purpose |
|---------|-------|-----------|---------|
| **Orthrus** | `quietflamingo/orthrus-large-6-track` | 512 | RNA transcript embedding (Mamba SSM, 6-track) |
| **ESM2** | `facebook/esm2_t33_650M_UR50D` | 1280 | Protein language model embedding |
| **metapredict / STARLING** | — | 50 | Per-residue disorder, aggregated to 50-dim features |
| **StructuralFeatureEncoder** | Daedalus tables | 96 | UniProt topology + BioGRID PPI + HPA surfaceome |

### Trainable Architecture (5.06M parameters)

- **ReferenceDelta** (×4): Asymmetric cross-attention delta. Captures what the reference had that the alternative **lost** (loss gate), what the alternative gained (gain gate), and how the alternative attends over the reference. Swapping ref↔alt gives a different output.
- **CrossModalFusion**: Cross-attention between RNA and protein deltas, then 3-layer MLP fusion of all four modalities → 256-dim fused representation.
- **TaskHeads**: Seven prediction heads.

### Task Heads

| Task | Type | Description |
|------|------|-------------|
| `global_preservation` | Binary | Does the alt isoform preserve any known function of the reference? |
| `topology_preserved` | Binary | For membrane proteins: is topological orientation maintained? |
| `surface_retained` | Binary | Are surface-exposed functional residues retained? |
| `kinase_competent` | Binary | For kinases: is the activation loop/catalytic triad intact? |
| `tf_competent` | Binary | For TFs: are DNA-binding domains retained? |
| `disorder_preserved` | Binary | Is the intrinsic disorder profile functionally maintained? |
| `deviation_class` | 6-class | Type of functional divergence |

**Deviation classes:** 0=preserved, 1=truncation_partial, 2=truncation_nmd, 3=topology_loss, 4=domain_loss, 5=novel_fusion

### Structural Feature Vector (96 dimensions)

The 96-dim structural feature vector per protein is assembled from four Daedalus tables:

| Features | Dims | Source |
|----------|------|--------|
| UniProt domain / PTM annotations | [0–18] | `uniprot_features.tsv` (long-format, aggregated per protein) |
| Sequence-predicted TM topology | [19–27] | `gencode_protein_sequence_features.tsv` (Kyte-Doolittle) |
| BioGRID PPI + UniProt region priors | [28–37] | `biogrid_gene_interactions.tsv` + `uniprot_region_priors.tsv` |
| Phospho/DNA-binding/surfaceome | [38–61] | HPA + UniProt phospho/DNA features |
| Reserved | [62–95] | zeros (future expansion) |

---

## Surfaceome Validation Pipeline

This is the primary use case for novel AltAnalyze3 isoforms.

### Conceptual workflow

```
AltAnalyze3 long-read GFF
         |
   Novel junction BED/TSV
         |
   Junction-spanning peptide       ← unique AA sequence across the novel exon junction,
         |                            absent from reference protein
   validate_junction_topology.py
         |
   Predict TM topology of full alternative protein (Kyte-Doolittle)
         |
   Map junction peptide position onto topology
         |
   ┌──────────────────────────────────────────────┐
   │  transmembrane_helix  → novel TM topology    │
   │  extracellular        → surface epitope ✓    │
   │  cytoplasmic          → not accessible       │
   │  signal_peptide       → secreted/GPI anchor  │
   │  nmd_absent           → protein not made     │
   └──────────────────────────────────────────────┘
         |
   junction_topology_validation.tsv
   (surface-accessible junctions → CAR-T / antibody candidates)
```

The **junction peptide** is the short amino acid sequence that *spans the novel exon junction and is unique to the alternative isoform* — it does not appear in the reference protein sequence or any other GENCODE isoform. This is what AltAnalyze3's long-read pipeline identifies.

### Run surfaceome validation

```bash
# AltAnalyze3 output → junction topology validation
python scripts/validate_junction_topology.py \
    --input  altanalyze3_novel_junctions.tsv \
    --output junction_topology_validation.tsv

# Only surface-accessible results
python scripts/validate_junction_topology.py \
    --input  altanalyze3_novel_junctions.tsv \
    --output surface_accessible_junctions.tsv \
    --surface_only

# If protein sequence is in a FASTA instead of the TSV
python scripts/validate_junction_topology.py \
    --input  altanalyze3_novel_junctions.tsv \
    --fasta  novel_isoform_proteins.fa \
    --output junction_topology_validation.tsv
```

**Required input TSV columns:**

| Column | Description |
|--------|-------------|
| `pair_id` | Unique identifier for this reference/alternative pair |
| `alternative_protein_seq` | Full translated AA sequence of the alternative isoform |
| `junction_peptide` | 6–30 AA sequence spanning the novel junction (unique to alt) |
| `has_alt_protein` | 1 or 0 — 0 if this isoform is predicted to undergo NMD |

**Output columns added:**

| Column | Description |
|--------|-------------|
| `topology_of_junction` | transmembrane_helix / extracellular / cytoplasmic / signal_peptide / soluble / nmd_absent / unknown |
| `is_surface_accessible` | 1 if junction is in TM helix, extracellular, or signal peptide region |
| `n_tm_helices` | Number of TM helices predicted in the full alt protein |
| `tm_segments` | TM helix positions, e.g. `23-45;67-89` |
| `junction_start_aa` | 0-based start of junction peptide in alt protein |
| `confidence` | 0–1 confidence in topology call |
| `evidence_summary` | Human-readable description |

---

## Full Training Workflow

### Prerequisites

```bash
cd /path/to/Proteus
pip install -r requirements.txt

# For Orthrus RNA encoder (requires CUDA for mamba_ssm)
pip install mamba-ssm

# For STARLING disorder prediction (optional, metapredict is the fallback)
pip install starling

# For experiment tracking (optional)
pip install wandb
```

Daedalus Phase A must be run first; its `data/interim/` directory contains
the source tables used by all extraction scripts.

### Step 1: Initialize from Daedalus pairs

```bash
python scripts/init_proteus.py \
    --daedalus_interim /path/to/daedalus/phase_a/data/interim
```

This copies the isoform pair TSV, splits into train/val/test, and symlinks
the Daedalus interim directory into `data/interim/`.

### Step 2: Extract features (CPU steps — parallelizable)

```bash
# Disorder features (metapredict/STARLING, 50-dim per protein)
python scripts/extract_disorder_features.py \
    --fasta data/raw/gencode/gencode.v49.pc_translations.fa.gz \
    --output_dir data/interim/disorder_features \
    --n_workers 16

# Structural features (96-dim, from Daedalus tables)
python scripts/extract_structural_features.py \
    --daedalus_interim data/interim \
    --output_dir data/interim/structural_features

# Sequence topology features for novel isoforms (9-dim, Kyte-Doolittle)
python scripts/extract_sequence_topology_features.py \
    --fasta novel_isoform_proteins.fa \
    --output_dir data/interim/sequence_topology_features \
    --n_workers 8
```

### Step 3: Extract embeddings (GPU steps)

```bash
# RNA embeddings via Orthrus (GPU required)
python scripts/extract_rna_embeddings.py \
    --fasta data/raw/gencode/gencode.v49.transcripts.fa.gz \
    --output_dir data/interim/rna_embeddings \
    --batch_size 32

# Protein embeddings via ESM2 (GPU recommended)
python scripts/extract_protein_embeddings.py \
    --fasta data/raw/gencode/gencode.v49.pc_translations.fa.gz \
    --output_dir data/interim/protein_embeddings \
    --batch_size 16
```

### Step 4: Generate synthetic pretraining events

```bash
python scripts/generate_pretrain_events.py \
    --gtf data/raw/gencode/gencode.v49.annotation.gtf.gz \
    --daedalus_interim data/interim \
    --output_dir data/interim/pretrain_events \
    --n_workers 8
```

Generates ~2M synthetic exon-skip events with rule-based labels:
`domain_disrupted`, `nmd_triggered`, `tm_lost`, `disorder_delta`.

### Step 5: Pretrain (GPU)

```bash
python scripts/run_pretrain.py \
    --config configs/pretrain.yaml \
    --output_dir data/processed/pretrain
```

### Step 6: Fine-tune (GPU)

```bash
python scripts/run_train.py \
    --pretrain_checkpoint data/processed/pretrain/best_model.pt \
    --config configs/training.yaml \
    --output_dir data/processed/finetune
```

### Step 7: Evaluate

```bash
python scripts/run_evaluate.py \
    --checkpoint data/processed/finetune/best_model.pt \
    --split test
```

---

## GPU Cluster Deployment (Singularity/Docker)

Yes — if you are running on an HPC cluster, you will need a container.
Most HPC systems use **Singularity** (now Apptainer), not Docker directly.

### Why containerization is required

- `mamba_ssm` (required by Orthrus) requires a specific CUDA version and compiled CUDA kernels
- ESM2 and STARLING also have CUDA dependencies that conflict with cluster system libraries
- Container images ensure reproducible CUDA environments across nodes

### Recommended approach

**Option A: Build a Docker image, convert to Singularity SIF**

```dockerfile
# Dockerfile (base: CUDA 12.1 + PyTorch 2.1)
FROM pytorch/pytorch:2.1.2-cuda12.1-cudnn8-devel

RUN pip install --no-cache-dir \
    transformers>=4.36 \
    einops \
    pandas \
    numpy \
    scikit-learn \
    tqdm \
    pyyaml \
    h5py

# mamba_ssm requires pre-compiled CUDA kernels
RUN pip install mamba-ssm --no-build-isolation

WORKDIR /workspace
COPY . /workspace/Proteus
RUN pip install -e /workspace/Proteus
```

```bash
# Build Docker image
docker build -t proteus:latest .

# Convert to Singularity SIF (on a machine with both Docker and Singularity)
singularity build proteus.sif docker-daemon://proteus:latest

# Or pull from Docker Hub if you push the image
singularity pull docker://yourorg/proteus:latest
```

**Option B: Singularity definition file directly**

```singularity
Bootstrap: docker
From: pytorch/pytorch:2.1.2-cuda12.1-cudnn8-devel

%post
    pip install transformers einops pandas numpy scikit-learn tqdm pyyaml
    pip install mamba-ssm --no-build-isolation

%environment
    export PYTHONPATH=/workspace/Proteus:$PYTHONPATH
```

```bash
sudo singularity build proteus.sif proteus.def
```

### Running on SLURM

Ready-to-submit SLURM scripts are in [cluster/](cluster/):

```bash
# Step 1: GPU embedding extraction (Orthrus + ESM2)
sbatch cluster/slurm_extract_embeddings.sh

# Step 2: Pretraining
sbatch cluster/slurm_pretrain.sh

# Step 3: Fine-tuning (after pretraining completes)
sbatch cluster/slurm_train.sh
```

Edit `DATA_DIR` and `SIF_PATH` at the top of each script to match your cluster paths. The `--nv` flag passes NVIDIA GPU access into the Singularity container.

### CPU-only steps (no container required)

Feature extraction steps (disorder, structural, topology) can run without a container on any Linux/macOS machine with Python ≥3.9:

```bash
pip install metapredict pandas numpy tqdm
python scripts/extract_disorder_features.py ...
python scripts/extract_structural_features.py ...
python scripts/validate_junction_topology.py ...
```

---

## Inference: Novel Isoforms

### Predict function preservation for a novel isoform pair

```python
from proteus import ProteusModel
import torch

model = ProteusModel.load_checkpoint("data/processed/best_model.pt")
model.eval()

batch = {
    "ref_rna_emb":        torch.zeros(1, 512),   # from Orthrus
    "alt_rna_emb":        torch.zeros(1, 512),
    "ref_protein_emb":    torch.zeros(1, 1280),  # from ESM2
    "alt_protein_emb":    torch.zeros(1, 1280),
    "ref_disorder_raw":   torch.zeros(1, 50),
    "alt_disorder_raw":   torch.zeros(1, 50),
    "ref_structural_raw": torch.zeros(1, 96),    # 96-dim structural features
    "alt_structural_raw": torch.zeros(1, 96),
    "has_alt_protein":    torch.ones(1, dtype=torch.bool),
}

with torch.no_grad():
    preds = model.predict(batch)

print(f"Global preservation:    {preds['global_preservation'].item():.3f}")
print(f"Surface retained:       {preds['surface_retained'].item():.3f}")
print(f"Topology preserved:     {preds['topology_preserved'].item():.3f}")
print(f"Deviation class:        {preds['deviation_class'].argmax().item()}")
```

### Predict TM topology for a novel protein sequence

```python
from proteus.encoders.topology import predict_tm_segments, predict_signal_peptide, predict_topology_type

seq = "MKTIIALSYIFCLVFA..."  # novel alternative isoform

tm_segs = predict_tm_segments(seq)
sig_candidate, sig_score = predict_signal_peptide(seq)
topology_type = predict_topology_type(seq)

print(f"TM helices: {len(tm_segs)} at positions {tm_segs}")
print(f"Signal peptide: {sig_candidate} (score={sig_score:.2f})")
print(f"Topology class: {topology_type}")  # single_pass, multi_pass, soluble, etc.
```

---

## File Structure

```
Proteus/
├── README.md                   # This file: workflow and run instructions
├── PROTEUS_METHODS.md          # Architecture, databases, benchmarking details
├── requirements.txt
├── configs/
│   ├── model.yaml              # Architecture hyperparameters
│   ├── encoders.yaml           # Encoder settings and cache dirs
│   ├── training.yaml           # Fine-tuning hyperparameters
│   └── pretrain.yaml           # Pretraining hyperparameters
├── proteus/
│   ├── model.py                # ProteusModel (main trainable model)
│   ├── pretrain_model.py       # ProteusPretrainModel (SSL pretraining)
│   ├── encoders/
│   │   ├── rna.py              # OrthrusRNAEncoder (frozen)
│   │   ├── protein.py          # ESM2ProteinEncoder (frozen)
│   │   ├── disorder.py         # DisorderEncoder (50-dim)
│   │   ├── structural.py       # StructuralFeatureEncoder (96-dim)
│   │   ├── topology.py         # SequenceTopologyEncoder (novel isoforms, 9-dim)
│   │   └── ppi_interface.py    # PPIInterfaceEncoder (29-dim standalone)
│   ├── delta/
│   │   ├── reference_delta.py  # ReferenceDelta (asymmetric cross-attention)
│   │   └── cross_modal.py      # CrossModalFusion
│   ├── heads/
│   │   └── task_heads.py       # Multi-task prediction heads
│   ├── data/
│   │   ├── dataset.py          # ProteusDataset
│   │   ├── pretrain_dataset.py # SyntheticExonSkipDataset
│   │   └── collate.py
│   ├── training/
│   │   ├── losses.py           # ProteusLoss, PretrainLoss
│   │   ├── metrics.py          # AUROC, AUPRC, ranking metrics
│   │   └── trainer.py          # ProteusTrainer
│   ├── pretraining/
│   │   ├── synthetic.py        # SyntheticExonSkipGenerator
│   │   └── objectives.py       # PretrainObjective
│   ├── validation/
│   │   ├── junction_topology.py  # JunctionTopologyValidator (surfaceome)
│   │   └── surfaceome.py         # SurfaceomeValidator (multi-evidence)
│   └── utils/
│       ├── device.py
│       └── sequences.py
├── scripts/
│   ├── init_proteus.py
│   ├── extract_rna_embeddings.py
│   ├── extract_protein_embeddings.py
│   ├── extract_disorder_features.py
│   ├── extract_structural_features.py
│   ├── extract_sequence_topology_features.py  # Novel isoform topology (9-dim)
│   ├── validate_junction_topology.py           # Surfaceome validation
│   ├── generate_pretrain_events.py
│   ├── run_pretrain.py
│   ├── run_train.py
│   └── run_evaluate.py
├── tests/
│   ├── test_model.py
│   └── test_reference_delta.py
├── cluster/
│   ├── slurm_extract_embeddings.sh  # GPU: RNA + protein embedding extraction
│   ├── slurm_pretrain.sh            # GPU: Stage 1 pretraining
│   └── slurm_train.sh               # GPU: Stage 2 fine-tuning
├── Dockerfile                       # CUDA 12.1 + PyTorch 2.1 image
├── singularity.def                  # Singularity definition (convert from Docker or build directly)
└── data/
    ├── raw/                    # GENCODE/UniProt source files (not tracked in git)
    ├── interim/                # Cached embeddings and features (not tracked)
    └── processed/              # Train/val/test TSVs, model checkpoints
```

---

## Citation

If you use Proteus in your research, please cite the Daedalus/AltAnalyze3 framework.

## License

See the AltAnalyze3 repository license.
