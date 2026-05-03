# Daedalus

## Purpose

`Daedalus` is the planning and implementation workspace for a reference-conditioned
isoform function prediction model built on top of AltAnalyze3 long-read outputs.
The target problem is not generic transcript embedding. The target is:

- given a canonical reference isoform and an observed alternative isoform,
- predict whether core function is preserved or deviates,
- with specialized outputs for:
  - transmembrane / surface proteins
  - kinases
  - transcription factors

This document focuses on feasibility under a practical hardware constraint:

- `1 x H100`

## Feasibility Assessment

### Bottom line

A credible first-generation methods model is feasible on `1 x H100` if we do **not**
attempt Orthrus-scale or ESM-scale pretraining from scratch.

What is feasible:

- build a large paired isoform dataset
- reuse established pretrained protein representations
- train a medium-size reference-conditioned multimodal predictor
- run strong baselines and ablations
- obtain a publishable first model if the benchmark is designed correctly

What is not realistic on `1 x H100`:

- foundation-scale pretraining from raw transcript corpora
- ESM3-scale generative multimodal protein training
- Orthrus-style pretraining across tens of millions of transcripts as a first step

So the correct strategy is:

1. use existing pretrained backbones where useful
2. train a task-specific paired model
3. make the benchmark and supervision the real innovation

## How To Run

Use the project-local interpreter:

```bash
Daedalus/.venv/bin/python
```

### Phase A: build the benchmark substrate

Download resources:

```bash
python3 Daedalus/phase_a/scripts/download_resources.py
python3 Daedalus/phase_a/scripts/download_uniprot_reviewed.py
python3 Daedalus/phase_a/scripts/download_appris_references.py
python3 Daedalus/phase_a/scripts/download_gencode_references.py
```

Build the Phase A tables:

```bash
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_hpa_localization_table.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_clinvar_splice_table.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_mane_table.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_gencode_reference_tables.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_uniprot_isoform_tables.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_uniprot_region_priors.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_uniprot_functional_priors.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_gencode_protein_sequence_features.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_biogrid_gene_interactions.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_uniprot_gencode_map.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_gene_supervision_catalog.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_transcript_supervision_catalog.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_isoform_pair_candidates.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_phase_a_splits.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_priority_pair_subsets.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_baseline_feature_matrix.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/build_seed_task_summary.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/run_sklearn_seed_baselines.py
Daedalus/.venv/bin/python Daedalus/phase_a/scripts/run_torch_seed_baselines.py
```

### Phase B: build the reference-conditioned matrix and baselines

```bash
Daedalus/.venv/bin/python Daedalus/phase_b/scripts/build_reference_delta_matrix.py
Daedalus/.venv/bin/python Daedalus/phase_b/scripts/build_multitask_instances.py
Daedalus/.venv/bin/python Daedalus/phase_b/scripts/run_reference_delta_sklearn_baseline.py
Daedalus/.venv/bin/python Daedalus/phase_b/scripts/train_reference_delta_multitask.py
```

### Query a novel isoform

```bash
Daedalus/.venv/bin/python Daedalus/phase_b/scripts/predict_query_isoform.py \
  --species human \
  --gene-name EGFR \
  --alt-protein-fasta novel_isoform.fa
```

This path now supports automatic `MANE/APPRIS` reference selection. It returns:

- reference-selected query context
- current Phase B task probabilities
- explicit evidence channels for:
  - PPI impact
  - PDI / DPI impact
  - localization shift
  - TF DNA binding
  - kinase signaling
  - signal peptide retention
  - TM insertion / fold support

### Main outputs

- Phase A feasibility: `Daedalus/phase_a/FEASIBILITY.md`
- Phase B status: `Daedalus/phase_b/STATUS.md`
- Phase B validation plan: `Daedalus/phase_b/VALIDATION.md`
- Phase B query workflow: `Daedalus/phase_b/QUERY_WORKFLOW.md`
- Phase C status: `Daedalus/phase_c/STATUS.md`
- Phase D objectives and proxy benchmark status: `Daedalus/phase_d/README.md`

## Documentation Map

- `Daedalus/METHODOLOGY.md`
  end-to-end construction of the approach, external databases used, workflow
  steps, and current benchmark/validation state
- `Daedalus/COMPARISON_WITH_PROTEUS.md`
  grounded comparison of Daedalus versus Proteus, including pros/cons,
  strengths/weaknesses, and the current benchmark state of each project
- `Daedalus/phase_b/QUERY_WORKFLOW.md`
  query-time workflow for novel isoforms, including automatic `MANE/APPRIS`
  reference selection and structured channel evidence outputs
- `Daedalus/phase_c/STATUS.md`
  Phase C segment-aware biochemical benchmark status and current model ranking
- `Daedalus/phase_a/FEASIBILITY.md`
  Phase A benchmark substrate and initial feasibility results
- `Daedalus/phase_b/STATUS.md`
  current Phase B baseline results
- `Daedalus/phase_b/VALIDATION.md`
  external TF validation plan

## Container Images

Container specs are in:

- `Daedalus/container/Dockerfile`
- `Daedalus/container/apptainer.def`
- `Daedalus/container/requirements.txt`
- `Daedalus/container/run_apptainer.sh`

### Build the Docker/OCI image

Run from `altanalyze3/altanalyze3/components`:

```bash
docker build -f Daedalus/container/Dockerfile -t daedalus:latest .
```

### Build the Apptainer/Singularity image

```bash
apptainer build daedalus.sif Daedalus/container/apptainer.def
```

### Run on a GPU node with Apptainer

```bash
apptainer exec --nv \
  --bind /path/to/altanalyze3/altanalyze3/components:/workspace \
  daedalus.sif \
  bash
```

Inside the container:

```bash
cd /workspace
python Daedalus/phase_b/scripts/train_reference_delta_multitask.py
```

Or use the helper wrapper:

```bash
Daedalus/container/run_apptainer.sh /path/to/daedalus.sif
Daedalus/container/run_apptainer.sh /path/to/daedalus.sif \
  python Daedalus/phase_b/scripts/train_reference_delta_multitask.py
```

### Practical recommendation

For cluster jobs, the clean path is:

1. build a Docker/OCI image locally or on a build node
2. build or pull an Apptainer/Singularity image from that environment
3. run with `--nv` on GPU nodes

This is not strictly required for GPU execution, but it is the safest path for
cluster portability and reproducibility.

## Phase C

Phase C is now implemented as a new path for:

- segment-aware biochemical features
- class-specific autoencoder modeling
- explicit latent dimensional reduction
- task-specific heads

The current Phase C subset benchmark indicates:

- `XGBoost` on the rich biochemical feature matrix is currently the strongest
  classifier
- the deep autoencoder is implemented and benchmarked, but it is not yet the
  default recommendation
- `segment_change` features are consistently useful
- several current typed-overlap and broad physicochemical blocks are noisy and
  need tightening
- a separate DeepImmuno-style benchmark path is now implemented with:
  - PCA reduction on physicochemical blocks
  - systematic ML/DL comparison
  - CNN and CNN+autoencoder classifiers
- on the completed five-task DeepImmuno-style benchmark:
  - `hist_gradient_boosting` ranks first on all five tasks
  - `random_forest` ranks second on all five tasks
  - strongest deep model by task:
    - `global`: `deepimmuno_cnn_autoencoder`
    - `membrane`: `deepimmuno_cnn`
    - `surface`: `deepimmuno_cnn`
    - `kinase`: `deepimmuno_cnn_autoencoder`
    - `transcription_factor`: `deepimmuno_cnn_autoencoder`
  - consolidated results are in:
    - `Daedalus/phase_c/checkpoints/deepimmuno_style_benchmark_all_tasks.md`

See:

- `Daedalus/phase_c/README.md`
- `Daedalus/phase_c/STATUS.md`

## Phase D

Phase D refactors the broad membrane and TF outputs into explicit task heads.

Membrane heads:

- `secretory_pathway_compatibility`
- `membrane_insertion_topology`
- `extracellular_altered_segment_accessibility`
- `folding_stability_qc_escape`
- `cell_surface_localization`
- `antibody_targetability`
- aggregate: `tm_stable_functional`, `tm_non_stable_non_functional`

TF heads:

- `nuclear_localization`
- `dbd_integrity_specificity_shift`
- `cofactor_ppi_rewiring`
- `activation_repression_competence`
- `dominant_negative_neomorphic_behavior`
- aggregate: `overall_regulatory_functionality`, `tf_non_stable_non_functional`

Phase D now has:

- objective scoring
- benchmark-instance generation
- family-aware multitask training
- model comparison against logistic / boosted trees / xgboost

Current boundary:

- the Phase D proxy benchmark runs end-to-end
- many thresholded proxy tasks are extremely one-sided
- so the current Phase D numbers validate the pipeline, not the biology

See:

- `Daedalus/phase_d/README.md`
- `Daedalus/phase_d/TRAINING_AND_BENCHMARK_SPEC.md`

## Recommended Framework Stack

### 1. Protein representation backbone

Recommended:

- `ESM2` or a comparable open protein encoder as the initial protein sequence backbone

Why:

- mature and practical
- strong residue-level and sequence-level protein features
- realistic on `1 x H100`
- suitable for fine-tuning or frozen-encoder feature extraction

Not recommended as the initial backbone:

- `ESM3`

Reason:

- conceptually attractive, but too large and too ambitious for first-pass model
- more useful as inspiration for multimodal design than as the first engineering target

Practical use:

- encode reference and alternative translated ORFs
- extract whole-protein and local edited-region embeddings
- optionally freeze for v1, then selectively unfreeze later

### 2. Transcript structure encoder

Recommended:

- a custom lightweight splice-graph / transcript-delta encoder

Implementation options:

- small Transformer over exon/junction tokens
- Graph Neural Network over exon-junction graph
- compact Mamba/SSM encoder over transcript feature tokens

Best choice for v1:

- a compact Transformer or graph encoder, not a large Mamba stack

Reason:

- the transcript input is structured and sparse
- the problem is pairwise delta reasoning, not raw long-context transcript language modeling
- easier to debug and benchmark than a custom large SSM stack

### 3. Fusion / comparison module

Recommended:

- reference-conditioned cross-attention or siamese-delta fusion

The model should not embed each isoform independently and compare by cosine only.
It should explicitly model:

- what regions changed
- where the ORF changed
- whether critical segments were deleted, truncated, inserted, or frame-shifted

Best choice:

- paired encoder with cross-attention focused on altered segments

### 4. Family-specific heads

Recommended:

- multitask heads for:
  - global preservation
  - membrane topology / surface retention
  - kinase catalytic competence
  - TF regulatory competence

This is necessary. A single scalar “function score” is too weak.

### 5. Teacher models and auxiliary annotation engines

Recommended:

- `DeepTMHMM` for TM topology pseudo-labels
- `SignalP 6.0` for signal peptide pseudo-labels
- `InterPro` / `Pfam` / `UniProt` domain mappings
- optional AlphaFold DB feature summaries where available

These should be used as:

- supervision sources
- auxiliary targets
- benchmark strata

Not as the final method by themselves.

## Proposed Main Algorithm

## Name

`Daedalus`: a reference-conditioned multimodal isoform delta network

## Inputs

For each pair `(reference_isoform, alternative_isoform)`:

- transcript structure
  - exon coordinates
  - junction chain
  - CDS coordinates
  - UTR/CDS segmentation
  - frame / phase transitions
- translated ORF / protein sequence
- pairwise edit map
  - deleted exons
  - inserted exons
  - altered splice junctions
  - frame-preserving vs frame-disrupting edits
- domain and topology priors
  - InterPro/Pfam/UniProt
  - TM helices
  - signal peptide
  - extracellular/intracellular segmentation
- family flags
  - membrane protein
  - kinase
  - transcription factor

## Architecture

### Branch A: transcript delta encoder

Encodes the structure difference between the reference and alternative transcript.

Suggested tokenization:

- exon tokens
- junction tokens
- CDS segment tokens
- edit-operation tokens

Suggested backbone:

- 4-8 layer compact Transformer or graph encoder

### Branch B: protein encoder

Encodes both reference and alternative proteins using a pretrained protein model.

Outputs:

- global reference embedding
- global alternative embedding
- local embeddings around altered protein regions

### Branch C: region delta module

Cross-attention over:

- reference altered neighborhood
- alternative altered neighborhood
- aligned domain/topology landmarks

This branch is load-bearing. It is where the model learns that not all changes are
equal.

### Fusion

Fuse:

- transcript delta representation
- protein sequence representation
- local altered-region representation
- family metadata

### Heads

- `P_preserved`
- `P_topology_preserved`
- `P_surface_retained`
- `P_kinase_competent`
- `P_tf_competent`
- `deviation_class`
- attribution map over altered regions

## Why this is more innovative than Orthrus

Orthrus learns transcript embeddings and then infers functional relatedness from
distance in embedding space.

`Daedalus` would instead:

- reason relative to a canonical reference isoform
- directly model the structural delta between isoforms
- combine transcript structure with protein representation
- use task-specific functional supervision
- output calibrated probabilities for defined functional competencies

That is a different methodological contribution.

## Training Data Strategy

### Hard labels

Primary curated supervision sources:

- `ClinVar` splice/pathogenic transcript effects
- `UniProt` reviewed isoforms and functional annotations
- `APPRIS` principal isoforms
- `TRIFID` as a weak functional-priority signal
- `Human Protein Atlas` localization labels
- curated kinase and TF resources

### Teacher / pseudo-label sources

- `DeepTMHMM`
- `SignalP`
- domain architecture from `InterPro` / `Pfam`
- AlphaFold-derived region summaries where available

### Pair construction

Positive-ish / preserved pairs:

- canonical vs functionally supported alternative isoforms
- MANE/APPRIS-supported close alternatives
- isoforms with proteomics support and preserved localization/domain evidence

Negative / deviated pairs:

- NMD-like truncations
- splice-disrupted pathogenic isoforms
- membrane-to-secreted or membrane-topology-breaking isoforms
- kinase catalytic motif disruptions
- TF DNA-binding or NLS disruptions

### Split design

Must split by:

- gene
- optionally family / orthogroup for hardest generalization test

Never allow isoforms from the same gene in both train and test.

## Training Plan for 1 x H100

### Phase 0: dataset and baselines

Time:

- `2-4 weeks`

Compute:

- mostly CPU
- light GPU optional

Work:

- build pair dataset
- compute domain/topology annotations
- build frozen train/val/test splits
- train logistic regression / XGBoost / small MLP baselines

This phase is mandatory.

### Phase 1: medium paired model

Time:

- `1-2 weeks` initial training
- `1-3 weeks` ablations and fixes

Compute:

- feasible on `1 x H100`

Suggested model size:

- `50M-300M` trainable parameters

Practical approach:

- freeze most of the protein backbone initially
- train fusion and task heads
- then optionally unfreeze top protein layers

### Phase 2: refinement

Time:

- `2-4 weeks`

Work:

- family-specific ablations
- calibration
- error analysis
- external validation

## Reliable Evaluation

This project only works if the benchmark is harder than the model.

Required metrics:

- AUROC / AUPRC for binary heads
- macro-F1 for deviation class
- calibration error / Brier score
- within-gene ranking accuracy

Required baselines:

- ORF length / NMD only
- domain-overlap count
- TRIFID
- APPRIS principal-vs-alternative rule
- protein embedding cosine
- Orthrus-style transcript embedding distance if available

Required external tests:

- held-out membrane proteins
- held-out kinases
- held-out TF families
- held-out genes with known functional isoform differences

## Recommended Software Frameworks

### Training

- `PyTorch`
- `PyTorch Lightning` or `Lightning Fabric`
- `Hydra` for configs

### Protein embeddings

- `fair-esm` / open ESM2 implementation

### Structure/topology teachers

- `DeepTMHMM`
- `SignalP`

### Annotation

- `InterPro`
- `Pfam`
- `UniProt`
- `APPRIS`
- `TRIFID`
- `Human Protein Atlas`

### Experiment tracking

- `Weights & Biases` or `MLflow`

## Recommendation

For a `1 x H100` budget, the best main framework is:

- a **paired multimodal predictor**
- with a **frozen or lightly fine-tuned protein backbone**
- plus a **custom transcript delta encoder**
- plus **multitask family-specific heads**

The recommended starting backbone choice is:

- `ESM2`-class protein encoder

The recommended transcript-side framework is:

- compact Transformer or graph encoder

The recommended overall strategy is:

- not Orthrus-scale transcript pretraining
- not ESM3-scale generative modeling
- build the strongest benchmark first
- then train the smallest model that beats all simple baselines

## Immediate Next Steps

1. define exact supervision targets
2. build the isoform-pair training table
3. create frozen splits
4. run simple baselines
5. only then train the paired model
