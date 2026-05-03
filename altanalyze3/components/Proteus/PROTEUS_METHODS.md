# Proteus: Methods, Databases, and Benchmarking

This document describes the construction of Proteus in detail — the databases used, how each feature type was derived, the model architecture rationale, the training strategy, and how validation results should be interpreted.

---

## 1. Biological Motivation

Alternative splicing affects >90% of human multi-exon genes. The functional consequence of a given isoform — whether it preserves the reference protein's activity, disrupts a domain, alters membrane topology, or is degraded by NMD — determines whether it is a therapeutic target, a biomarker, or irrelevant noise.

Existing tools (SpliceAI, MMSplice, etc.) predict *splicing patterns* from sequence. Proteus answers a different question: **given that this alternative isoform is expressed, does it maintain the function of the canonical protein?**

For surfaceome biology specifically, a key question is: *Does the novel exon junction create an epitope that is surface-accessible?* This determines antibody accessibility, CAR-T target potential, and mass-spectrometry detectability.

---

## 2. Training Data: Daedalus Phase A

Proteus is trained and evaluated on isoform pairs derived from the **Daedalus Phase A** dataset, which was constructed from:

### 2.1 Transcript Sources

| Source | Version | Species | Transcripts |
|--------|---------|---------|-------------|
| GENCODE | v49 | Human (GRCh38) | ~252,000 |
| GENCODE | vM38 | Mouse (GRCm39) | ~149,000 |

Protein-coding transcripts with annotated CDS were retained. Pairs were constructed by grouping transcripts by gene (ENSG), designating the **MANE SELECT** or **APPRIS principal** isoform as the reference, and pairing it with each alternative isoform.

- **Total isoform pairs**: ~257,851
- **Human pairs**: ~83%
- **Mouse pairs**: ~17%

### 2.2 Weak Labels

Because no large-scale experimental dataset of isoform function exists, Proteus uses **biological priors** as weak labels. These are derived from:

| Label | Source | Logic |
|-------|--------|-------|
| `domain_disrupted` | UniProt features | Exon removal overlaps a UniProt domain or active site |
| `nmd_triggered` | CDS structure | PTC is >55 nt upstream of last exon junction complex (EJC) |
| `tm_lost` | UniProt + sequence | Exon removal removes a transmembrane helix |
| `kinase_competent` | UniProt domains | Kinase domain and activation loop retained in alt |
| `tf_competent` | UniProt domains | DNA-binding domain retained in alt |
| `surface_retained` | HPA + UniProt extracellular | Extracellular-annotated residues retained |
| `topology_preserved` | UniProt topology | TM helix count and order preserved |

These are noisy labels, not ground truth. The model is trained with uncertainty-aware masking — pairs without evidence for a given task are masked out of that task's loss.

### 2.3 Data Split

| Split | Pairs | Notes |
|-------|-------|-------|
| Train | ~206,280 (80%) | |
| Validation | ~25,785 (10%) | Held out per-gene (no gene appears in both train and val) |
| Test | ~25,786 (10%) | Held out per-gene |

Mouse pairs are included in training only (insufficient human-equivalent annotations for fair benchmarking).

---

## 3. Feature Encoders

### 3.1 RNA Encoder: Orthrus

**Model**: `quietflamingo/orthrus-large-6-track` (Mamba SSM, 6-track)
**Output**: 512-dimensional mean-pooled transcript embedding
**What it captures**: Splicing-aware RNA sequence representation. The 6-track input encodes nucleotide identity, splice signals, and reading frame simultaneously.
**Usage**: Embeddings are extracted once and cached as `.pt` files. The Orthrus model itself is never fine-tuned.

Orthrus was chosen over a simple k-mer embedding because it has been pretrained to distinguish functionally distinct isoforms and captures long-range exon-exon dependencies that affect splicing.

### 3.2 Protein Encoder: ESM2

**Model**: `facebook/esm2_t33_650M_UR50D` (650M parameter protein LM)
**Output**: 1280-dimensional mean-pooled residue embedding
**What it captures**: Protein sequence evolutionary context, secondary structure propensity, domain composition.
**Usage**: Embeddings extracted once and cached. ESM2 never fine-tuned.

For NMD-predicted alternative isoforms (where no protein is produced), the alternative protein embedding is zeroed out, forcing the protein delta module to represent "complete loss."

### 3.3 Disorder Encoder: metapredict / STARLING

**Output**: 50-dimensional feature vector
**What it captures**: Intrinsic disorder profile of the protein.

The 50 features summarize the per-residue disorder score distribution: mean, std, percentiles (5th, 10th, ..., 95th), IDR count, IDR fraction, longest IDR length, mean IDR length, and N/C-terminal disorder flags.

Both metapredict (CPU, fast) and STARLING (GPU, more accurate) are supported. metapredict uses a biLSTM trained on MoRFpred/disprot annotations. STARLING is an ensemble approach.

Disorder is relevant because IDRs are functionally significant — they mediate phase separation, hub interactions, and are enriched at regulatory protein interfaces. An isoform that loses or gains a large IDR may have altered interaction capacity.

### 3.4 Structural Feature Encoder (96 dimensions)

The structural encoder assembles a 96-dimensional feature vector per protein from four Daedalus data tables. This encoder represents the structural and functional *context* of the protein — topology, PTMs, PPIs — rather than deriving information from the raw protein sequence.

#### Feature layout

**[0–18] UniProt domain and PTM features** (from `uniprot_features.tsv`, 1.1M rows, long-format)

| Idx | Feature | Description |
|-----|---------|-------------|
| 0 | domain_count | log1p count of annotated UniProt domains |
| 1 | active_site_count | catalytic active sites |
| 2 | binding_site_count | substrate/cofactor binding sites |
| 3 | disulfide_count | disulfide bonds |
| 4 | glycosylation_count | N-/O-/C-glycosylation sites |
| 5 | phosphorylation_count | phosphoserine/threonine/tyrosine sites |
| 6 | signal_peptide_present | binary: signal peptide annotated |
| 7 | tm_helix_count | transmembrane helix count |
| 8 | tm_total_span | total TM span in amino acids |
| 9 | coiled_coil_count | coiled-coil region count |
| 10 | repeat_count | repeat region count |
| 11 | zinc_finger_count | zinc finger domain count |
| 12 | is_membrane_protein | binary: ≥1 TM helix |
| 13 | is_surface_protein | binary: HPA extracellular or UniProt extracellular topology |
| 14 | protein_length_norm | protein length / 1000 |
| 15 | n_terminal_tm | TM helix in first 50 aa |
| 16 | motif_count | linear motif count |
| 17 | modified_residue_count | other modified residues |
| 18 | nuclear_topology_count | nuclear localization regions |

**[19–27] Sequence-predicted TM topology** (from `gencode_protein_sequence_features.tsv`)

These are predicted by Daedalus using the same Kyte-Doolittle sliding-window algorithm as `SequenceTopologyEncoder` in `proteus/encoders/topology.py`:

| Idx | Feature |
|-----|---------|
| 19 | predicted_tm_count (log1p) |
| 20 | predicted_signal_candidate (binary) |
| 21 | predicted_signal_score |
| 22 | seq_tm_span_fraction |
| 23 | max_hydropathy_19aa |
| 24 | n_terminal_hydropathy |
| 25 | glyco_motif_count (log1p) |
| 26 | cysteine_fraction |
| 27 | predicted_topology_type / 5 |

**[28–37] BioGRID PPI + UniProt region priors**

| Idx | Feature | Source |
|-----|---------|--------|
| 28 | biogrid_partner_count (log1p) | BioGRID TAB3 |
| 29 | biogrid_physical_partner_count (log1p) | BioGRID TAB3 |
| 30 | uniprot_ppi_binding_feature_count (log1p) | uniprot_region_priors.tsv |
| 31 | uniprot_phospho_interaction_count (log1p) | uniprot_region_priors.tsv |
| 32 | extracellular_topology_count (log1p) | uniprot_region_priors.tsv |
| 33 | extracellular_aa_fraction | uniprot_region_priors.tsv |
| 34 | cytoplasmic_topology_count (log1p) | uniprot_region_priors.tsv |
| 35 | signal_count (log1p) | uniprot_region_priors.tsv |
| 36 | interaction_partner_count (log1p) | uniprot_region_priors.tsv |
| 37 | tm_total_span_log1p | uniprot_region_priors.tsv |

**[38–61] Phospho / DNA-binding / surfaceome**

| Idx | Feature | Source |
|-----|---------|--------|
| 38 | phospho_feature_count (log1p) | uniprot_region_priors.tsv |
| 39 | dna_binding_feature_count (log1p) | uniprot_region_priors.tsv |
| 40 | zinc_finger_feature_count (log1p) | uniprot_region_priors.tsv |
| 41 | dna_region_feature_count (log1p) | uniprot_region_priors.tsv |
| 42 | coiled_coil_count_rp (log1p) | uniprot_region_priors.tsv |
| 43 | glycosylation_count_rp (log1p) | uniprot_region_priors.tsv |
| 44 | disulfide_count_rp (log1p) | uniprot_region_priors.tsv |
| 45 | motif_count_rp (log1p) | uniprot_region_priors.tsv |
| 46 | domain_count_rp (log1p) | uniprot_region_priors.tsv |
| 47 | has_hpa_extracellular (binary) | hpa_localization.tsv |
| 48 | is_hpa_membrane_or_surface (binary) | hpa_localization.tsv |
| 49–55 | (extended modification features) | uniprot_region_priors.tsv |
| 56–61 | topology_type one-hot (6 classes) | gencode_protein_sequence_features.tsv |

**[62–95] Reserved** (zeros, future expansion)

---

## 4. Model Architecture

### 4.1 ReferenceDelta

The core novel contribution of Proteus. Unlike a simple vector difference (alt − ref), `ReferenceDelta` is **asymmetric and directional**.

```python
class ReferenceDelta(nn.Module):
    # Inputs: ref [B, D], alt [B, D]
    # Output: delta [B, out_dim]

    # Loss gate: what did the alternative LOSE relative to reference?
    loss_context = cross_attention(query=alt, key=ref, value=ref)
    loss_gate = sigmoid(mlp(concat(alt, ref, alt - ref, loss_context)))

    # Gain gate: what did the alternative GAIN not present in reference?
    gain_context = cross_attention(query=ref, key=alt, value=alt)
    gain_gate = sigmoid(mlp(concat(ref, alt, alt - ref, gain_context)))

    delta = concat(loss_gate * loss_context, gain_gate * gain_context)
    delta = layer_norm(projection(delta))
```

Swapping reference and alternative produces a different output — this directional asymmetry is critical because "the alternative lost a TM helix" is biologically different from "the reference lacks a TM helix that the alternative gained."

### 4.2 CrossModalFusion

Four delta vectors (RNA, protein, disorder, structural) are fused:

1. **Cross-modal attention**: RNA delta attends over protein delta and vice versa (RNA-protein co-variation)
2. **Concatenation**: all four deltas concatenated → 1024-dim
3. **3-layer MLP with residual connections**: 1024 → 512 → 256-dim fused representation

Modality dropout (10% probability during training) randomly zeroes one modality's delta, preventing the model from over-relying on any single input type and improving robustness to missing data at inference.

### 4.3 TaskHeads

Seven prediction heads share the 256-dim fused representation:

```python
class TaskHeads(nn.Module):
    # Binary heads (2-class logits)
    global_preservation:  Linear(256 → 2) + LayerNorm
    topology_preserved:   Linear(256 → 2)
    surface_retained:     Linear(256 → 2)
    kinase_competent:     Linear(256 → 2)
    tf_competent:         Linear(256 → 2)
    disorder_preserved:   Linear(256 → 2)

    # Multi-class head
    deviation_class:      Linear(256 → 6)
```

Task masks are applied at loss computation: kinase_competent is only computed for known kinases, tf_competent only for TFs, topology_preserved only for membrane proteins.

### 4.4 Parameter Count

| Component | Parameters |
|-----------|-----------|
| disorder_proj (50→64) | 3,264 |
| structural_proj (96→64) | 6,208 |
| rna_delta | ~1.31M |
| protein_delta | ~1.31M |
| disorder_delta | ~330K |
| structural_delta | ~330K |
| fusion (CrossModalFusion) | ~1.58M |
| task_heads | ~200K |
| **Total trainable** | **~5.06M** |

---

## 5. Training

### 5.1 Stage 1: Self-Supervised Pretraining

`SyntheticExonSkipGenerator` enumerates all single-exon skip variants of GENCODE multi-exon transcripts. For each synthetic skip:

- **domain_disrupted**: skipped exon overlaps a UniProt domain or active site (from `uniprot_features.tsv`)
- **nmd_triggered**: skip introduces a PTC >55 nt before the last exon junction (standard NMD rule)
- **tm_lost**: skipped exon overlaps a UniProt TM helix (from `uniprot_region_priors.tsv`)
- **disorder_delta**: signed fractional change in mean per-residue disorder score

Scale: ~1,960,432 synthetic events across human + mouse GENCODE.

**Pretraining objective** (from `proteus/pretraining/objectives.py`):
- BCE on domain_disrupted, nmd_triggered, tm_lost (binary labels)
- MSE on disorder_delta (continuous)
- Contrastive loss: embed (ref, skip-with-domain-loss) closer together than (ref, random-alt)

**Hyperparameters**: AdamW lr=3e-4, warmup 5%, cosine decay, batch=256, 30 epochs.

### 5.2 Stage 2: Supervised Fine-Tuning

Fine-tuning on Daedalus Phase A isoform pairs with task-masked BCE losses:

```
total_loss = sum_t [ w_t * mask_t * BCE(pred_t, label_t) ]
           + lambda_dev * CrossEntropy(deviation_class_pred, deviation_class)
```

Task weights: global_preservation=1.0, topology_preserved=1.5 (rarer, upweighted), surface_retained=1.5, kinase_competent=2.0, tf_competent=2.0, disorder_preserved=0.8.

**Hyperparameters**: AdamW lr=5e-5, warmup 10%, cosine decay, batch=128, 50 epochs, early stopping (patience=8 on val AUROC).

---

## 6. Databases

### 6.1 UniProt Protein Knowledgebase

**URL**: https://www.uniprot.org/
**Version**: 2024-04 (Swiss-Prot reviewed entries, human + mouse)
**Files used**:
- `uniprot_features.tsv`: Long-format feature annotations. ~1.1M rows. Each row is one feature (domain, active site, TM helix, phospho site, etc.) for one protein. Columns: `primary_accession`, `feature_type`, `description`, `start`, `end`.
- `uniprot_region_priors.tsv`: Wide-format aggregated per-protein statistics. Columns include `tm_total_span`, `extracellular_topology_count`, `extracellular_topology_aa`, `ppi_binding_feature_count`, `phospho_interaction_feature_count`, `dna_binding_feature_count`, `zinc_finger_feature_count`, etc.
- `uniprot_gencode_map.tsv`: Mapping table between UniProt primary accessions and GENCODE protein/transcript IDs.

**Key bridging issue**: UniProt accessions (e.g., `P00533` for EGFR) must be mapped to GENCODE ENSP IDs via `uniprot_gencode_map.tsv`. UniProt covers ~80% of GENCODE proteins; the remainder get zero-vectors for UniProt-derived features.

### 6.2 BioGRID

**URL**: https://thebiogrid.org/
**Version**: BioGRID 4.4.x (TAB3 format, human + mouse)
**File used**: `biogrid_gene_interactions.tsv` (aggregated from raw TAB3 by Daedalus)

**Columns used**:
- `gene_name`: HGNC/MGI symbol
- `biogrid_partner_count`: unique interaction partners (physical + genetic)
- `biogrid_physical_partner_count`: physical interaction partners only
- `biogrid_interaction_count`: total interaction records
- `biogrid_physical_interaction_count`: physical interaction records

**Why physical PPIs matter for surfaceome**: Membrane receptor complexes (e.g., EGFR + ERBB2, CD3ζ + LCK) require specific domain interfaces. An isoform that removes a domain involved in physical PPI may disrupt complex assembly and surface expression.

**Bridging to ENSP**: BioGRID uses gene symbols; Daedalus bridges via `gencode_protein_reference.tsv` → `gene_supervision_catalog.tsv` → HGNC symbol → ENSP.

### 6.3 Human Protein Atlas (HPA)

**URL**: https://www.proteinatlas.org/
**Version**: HPA v23
**File used**: `hpa_localization.tsv`

**Columns used**:
- `ensembl_gene_id`: ENSG (unversioned, e.g., `ENSG00000146648`)
- `has_extracellular_annotation`: binary, HPA detected extracellular localization
- `is_membrane_or_surface_annotated`: binary, any membrane/surface category

**Key implementation note**: HPA uses *unversioned* ENSG IDs (`ENSG00000146648`) while GENCODE uses *versioned* IDs (`ENSG00000146648.7`). The ID must be stripped to the base ENSG before lookup:
```python
gid_base = gid_versioned.split(".")[0]
hpa_row = hpa_by_gid.get(gid_versioned) or hpa_by_gid.get(gid_base)
```
Without this fix, HPA returns 0 entries for all proteins.

### 6.4 GENCODE

**Human**: GENCODE v49 (GRCh38.p14)
**Mouse**: GENCODE vM38 (GRCm39)
**Files used**: annotation GTF, transcript FASTA, protein translation FASTA, protein sequence feature TSV (generated by Daedalus).

**Protein sequence features** (`gencode_protein_sequence_features.tsv`): Generated by Daedalus `build_gencode_protein_sequence_features.py` using Kyte-Doolittle hydropathy sliding window (window=19, threshold=1.6 for TM calling, window=9 for signal peptide). This is the same algorithm implemented in `proteus/encoders/topology.py` for use on novel sequences not in the reference annotation.

---

## 7. Novel Isoform Topology Prediction

For novel alternative isoforms from AltAnalyze3 long-read sequencing (which by definition are not in GENCODE v49), sequence-derived topology must be predicted de novo.

### 7.1 SequenceTopologyEncoder

`proteus/encoders/topology.py` implements the same Kyte-Doolittle algorithm as Daedalus, applied to any amino acid sequence:

**Kyte-Doolittle hydropathy values** for TM helix detection:
- Window size: 19 residues (canonical TM helix length)
- Threshold: 1.6 (default) — segments where mean hydropathy ≥ 1.6 over ≥15 consecutive residues are called as TM helices
- Signal peptide detection: window=9, check first 35 aa, threshold=2.0

**Output features** (9-dim):
1. TM helix count (log1p)
2. Signal peptide candidate (binary)
3. Signal peptide score (normalized)
4. TM span fraction (TM aa / total aa)
5. Max mean hydropathy in any 19-aa window (normalized)
6. Mean hydropathy of first 19 aa (N-terminal)
7. N-glycosylation motif count (NxS/T, log1p)
8. Cysteine fraction (disulfide potential)
9. Topology type (0=soluble, 1=single-pass type I, 2=single-pass type II, 3=multi-pass, 4=GPI-anchored estimate, 5=other) / 5

### 7.2 JunctionTopologyValidator

`proteus/validation/junction_topology.py` addresses the specific surfaceome validation question:

> *Is the unique junction-spanning peptide that distinguishes this alternative isoform located in a surface-accessible region?*

**Algorithm**:
1. Locate the junction peptide (string search) in the full alternative protein sequence
2. Predict TM helices and signal peptide in the full alternative protein
3. Count what fraction of junction peptide residues overlap TM helices
   - If ≥50%: classify as `transmembrane_helix`
4. If not TM, use alternating topology model to assign inter-TM loops:
   - Loop 0 (before first TM): extracellular
   - Loop 1 (between TM1–TM2): cytoplasmic
   - Loop N: alternates (even = extracellular, odd = cytoplasmic)
   - C-terminal tail: extracellular for type I (single pass), cytoplasmic for multi-pass
5. If no TM helices: classify as `soluble` (or `signal_peptide` if signal detected and peptide is N-terminal)

**Surface accessibility definition**: topology ∈ {transmembrane_helix, extracellular, signal_peptide}

**Confidence scores**:
- TM helix overlap: 0.5 + overlap_fraction (max 1.0)
- Multi-pass extracellular/cytoplasmic: 0.55
- Single-pass: 0.45 (ambiguous orientation for novel proteins)
- NMD: 1.0 (high confidence — not accessible because protein absent)

---

## 8. Surfaceome Validation: Multi-Evidence Composite Scoring

`proteus/validation/surfaceome.py` provides `SurfaceomeValidator`, which combines five independent lines of evidence for surface-membrane status:

| Evidence | Weight | Source |
|----------|--------|--------|
| Proteus `surface_retained` logit → probability | 0.35 | ProteusModel prediction |
| Sequence-predicted TM helices (min(count/4, 1)) | 0.20 | SequenceTopologyEncoder |
| UniProt TM annotation span (min(span/120, 1)) | 0.15 | uniprot_region_priors.tsv |
| UniProt extracellular topology count (min(count/3, 1)) | 0.10 | uniprot_region_priors.tsv |
| HPA extracellular annotation (binary) | 0.12 | hpa_localization.tsv |
| BioGRID physical partner count (min(count/10, 1)) | 0.08 | biogrid_gene_interactions.tsv |

**Validation threshold**: composite score ≥ 0.55.

**Example output for EGFR (ENSP00000275493)**:
```
confidence_score: 0.586
is_validated: True
evidence_flags:
  model_surface_prob=0.924
  uniprot_tm_span=23aa
  HPA_extracellular
  biogrid_physical_partners=1734
```

---

## 9. Benchmarking and Validation

*Note: Proteus v1.0 has not yet completed supervised fine-tuning. The benchmarks below describe the evaluation protocol and baselines.*

### 9.1 Evaluation Protocol

**Primary metric**: Area Under ROC Curve (AUROC) per task on held-out test set.
**Secondary metrics**: Area Under Precision-Recall Curve (AUPRC), calibration (ECE).

Tasks are evaluated only on proteins with available evidence for that task (task masks).

### 9.2 Baselines

| Baseline | Description |
|----------|-------------|
| **Majority class** | Always predict the majority label |
| **Sequence difference** | ESM2(alt) − ESM2(ref) linear classifier |
| **Protein length ratio** | alt_length / ref_length linear classifier |
| **SpliceAI + logistic** | SpliceAI junction scores → logistic regression |

### 9.3 Expected Benchmark Structure

```
Task                    Proteus AUROC    Sequence-diff AUROC
global_preservation        TBD               TBD
topology_preserved         TBD               TBD
surface_retained           TBD               TBD
kinase_competent           TBD               TBD
tf_competent               TBD               TBD
disorder_preserved         TBD               TBD
```

### 9.4 Pretraining Validation (completed)

Synthetic pretraining dataset statistics:
- Total synthetic events: 1,960,432
- `domain_disrupted` = 1 (positive): 214,122 (10.9%)
- `nmd_triggered` = 1: 1,143,796 (58.3%)
- `tm_lost` = 1: 34,142 (1.7%)

These rates are biologically plausible:
- NMD: ~58% of skips introduce a PTC (consistent with published ~55-60% estimates)
- Domain disruption: ~11% (consistent with ~10-15% domain coverage of exons)
- TM loss: ~1.7% (consistent with ~3-5% of proteins being membrane proteins × ~50% TM exon fraction)

### 9.5 Feature Quality (structural encoder)

After HPA versioning fix:
- Structural features non-zero: 290,255 / 290,255 (100%)
- HPA entries mapped: 185,472
- Disorder features non-zero: 289,697 / 290,255 (99.8%)

EGFR (ENSP00000275493.3) structural feature vector sum: 5.124 (non-trivial, HPA-annotated).

---

## 10. Limitations and Known Issues

1. **Alternating topology model is simplified**: The `_assign_topology()` method in `JunctionTopologyValidator` uses a canonical even/odd loop assignment. Type I vs type II orientation is not explicitly determined from sequence alone without full experimental topology data. For single-pass proteins, confidence is reduced to 0.45.

2. **Novel sequences lack UniProt context**: For AltAnalyze3 novel isoforms not in GENCODE, structural features default to sequence-predicted topology only. UniProt TM, BioGRID, and HPA features are zero-padded. The model was trained primarily on reference-annotated proteins, so out-of-distribution novel sequences may have reduced prediction quality.

3. **Weak label noise**: Training labels are derived from biological rules, not experimental assays. Estimated label noise: ~15-20% for `global_preservation`, ~10% for `tm_lost`, ~5% for `nmd_triggered`.

4. **Mouse model coverage**: Mouse proteins have lower UniProt/BioGRID coverage than human. Mouse pairs are included in pretraining for data augmentation but should not be used as the primary validation benchmark.

5. **ESM2 length limit**: ESM2 truncates at 1022 residues. Proteins longer than this have their C-terminal features truncated. This affects ~3% of human proteins.

---

## 11. Software Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| PyTorch | ≥2.1.0 | Core deep learning |
| transformers | ≥4.36 | ESM2 model and tokenizer |
| numpy | ≥1.24 | Feature computation |
| pandas | ≥2.0 | TSV loading |
| scikit-learn | ≥1.3 | Metrics, calibration |
| tqdm | ≥4.65 | Progress bars |
| pyyaml | ≥6.0 | Config loading |
| metapredict | ≥2.6 | Disorder prediction (CPU) |
| mamba-ssm | ≥1.2 | Orthrus RNA encoder (GPU/CUDA) |
| h5py | optional | HDF5 cache format |
| wandb | optional | Experiment tracking |

**Python**: ≥3.9
**CUDA**: ≥12.1 for mamba-ssm; ≥11.8 for ESM2

---

## 12. Reproducibility

All random seeds are set in the config files (`configs/training.yaml`, `configs/pretrain.yaml`). Feature extraction is deterministic (no randomness). Embedding extraction with Orthrus/ESM2 is deterministic given the same model checkpoint and tokenization.

To reproduce the full pipeline from scratch:

```bash
# 1. Run Daedalus Phase A to generate interim tables
# 2. python scripts/init_proteus.py --daedalus_interim .../data/interim
# 3. python scripts/extract_disorder_features.py ...
# 4. python scripts/extract_structural_features.py ...
# 5. python scripts/extract_rna_embeddings.py ...
# 6. python scripts/extract_protein_embeddings.py ...
# 7. python scripts/generate_pretrain_events.py ...
# 8. python scripts/run_pretrain.py
# 9. python scripts/run_train.py --pretrain_checkpoint ...
# 10. python scripts/run_evaluate.py --checkpoint ...
```

Expected total runtime: ~2-4 hours feature extraction (CPU, 16 cores), ~6-8 hours GPU pretraining (A100 40GB), ~4-6 hours GPU fine-tuning.
