# Daedalus Methodology

## Scope

`Daedalus` is a reference-conditioned isoform function prediction workflow built
on top of AltAnalyze3 long-read outputs. The model target is not generic
transcript similarity. The target is functional preservation or deviation of an
alternative isoform relative to a reference isoform.

Main modeled branches:

- global isoform preservation/deviation
- membrane/surface retention
- kinase competence
- transcription-factor competence

The query-time workflow also exposes structured evidence channels for:

- PPI impact
- PDI / DPI impact
- localization shift
- TF DNA-binding retention
- kinase signaling retention
- signal peptide retention
- TM insertion support
- TM fold support

## Workflow construction

### Phase A

Phase A constructs the training and benchmark substrate.

Steps:

1. download external public resources
2. normalize each source into stable interim tables
3. assemble gene-level and transcript-level supervision catalogs
4. construct reference-vs-alternative isoform pair tables
5. freeze gene-disjoint train/val/test splits
6. run leakage-resistant baseline models

### Phase B

Phase B converts the Phase A pair tables into a reference-conditioned modeling
matrix and trains baseline predictive models.

Steps:

1. build the `reference_delta_matrix`
2. expand it into multitask instances
3. train sklearn and torch baseline models
4. evaluate on frozen held-out gene-disjoint splits
5. compare weak-label benchmark behavior across tasks

## Databases and resources used

Transcript/reference anchors:

- GENCODE human `v49`
- GENCODE mouse `vM38`
- MANE `GRCh38`
- APPRIS

Protein/function/topology annotations:

- UniProt reviewed entries
- UniProt isoform comments and feature annotations
- Human Protein Atlas subcellular localization
- BioGRID interactions

Clinical weak-label support:

- ClinVar splice/pathogenic summary

Sequence-derived features:

- in-house lightweight signal/TM candidate predictions from translated GENCODE
  protein sequences

Planned but not yet fully integrated:

- InterPro / protein2ipr for larger domain coverage

## Current modeled variables

The current Phase B matrix includes:

- transcript and protein length context
- reference-known membrane/signal/kinase/TF annotations
- reference topology and region counts
- glycosylation, disulfide, phospho, PPI-binding, DNA-binding, zinc-finger,
  domain, and motif context
- BioGRID partner counts
- alternative sequence-derived TM/signal features
- alternative glyco-motif and cysteine context

At query time, the workflow additionally derives:

- alternative TM segment coordinates
- predicted signal-sequence support
- TM hydropathy context
- known BioGRID partner lists for the reference gene

To reduce weak-label leakage:

- direct alternative UniProt region annotations are excluded from the Phase B
  modeled matrix
- reference-known annotations are retained
- alternative isoforms contribute sequence-derived rather than directly curated
  region features

## Initial validation and benchmarking

### Phase A baseline feasibility

The initial non-foundation baselines showed that:

- `global`, `membrane`, `kinase`, and `transcription_factor` tasks are viable
- the original `surface` task was too weak before topology-aware expansion

See:

- `Daedalus/phase_a/FEASIBILITY.md`

### Phase B current held-out benchmark state

Current cleaned Phase B results on gene-disjoint held-out test data:

Torch multitask baseline:

- global `AP 0.840682`, `AUROC 0.882048`
- membrane `AP 0.839828`, `AUROC 0.875073`
- kinase `AP 0.866310`, `AUROC 0.903348`
- transcription_factor `AP 0.858565`, `AUROC 0.884078`
- surface `AP 0.837000`, `AUROC 0.871491`

sklearn logistic baseline:

- global `AP 0.810107`, `AUROC 0.845348`
- membrane `AP 0.809486`, `AUROC 0.843945`
- kinase `AP 0.853772`, `AUROC 0.874458`
- transcription_factor `AP 0.827434`, `AUROC 0.842662`
- surface `AP 0.802420`, `AUROC 0.835832`

Interpretation:

- all five current tasks are learnable on the weak-label benchmark
- the torch baseline is stronger than the sklearn baseline overall
- `surface` is improved but still needs stronger external biological validation

See:

- `Daedalus/phase_b/STATUS.md`

## External validation plan

The TF branch will not rely only on weak labels. It will be validated against
precomputed hESC experimental outputs:

- TF isoform PPI
- TF isoform PDI
- TF isoform transcriptional activation
- TF isoform localization
- TF isoform perturbation / bootstrap DEGs

See:

- `Daedalus/phase_b/VALIDATION.md`

## Query-time workflow

`Daedalus` now includes a query workflow for novel isoforms:

- `Daedalus/phase_b/scripts/predict_query_isoform.py`
- `Daedalus/phase_b/QUERY_WORKFLOW.md`

The reference isoform is selected automatically if not supplied by the user
using:

1. `MANE Select`
2. `APPRIS PRINCIPAL:1`
3. other `APPRIS PRINCIPAL:*`
4. `UNIPROT_SUPPORTED_LONGEST`
5. `LONGEST_PROTEIN`
6. `LONGEST_TRANSCRIPT`

This keeps the query path consistent with the Phase A pair-construction logic.

The current query output combines:

- learned Phase B task probabilities
- structured evidence for PPI / PDI / localization / DNA binding / kinase
  signaling / signal peptide / TM insertion / TM fold support

Only the five Phase B task probabilities are currently learned heads. The
additional channel outputs are structured evidence layers built from the learned
probabilities plus reference annotations and alternative sequence-derived
features.

## Phase C: segment-aware biochemical autoencoder

Phase C has now been added as the next model family beyond the Phase B
count-based multitask MLP.

The core changes are:

- explicit changed-segment detection between reference and alternative isoforms
- typed biochemical feature extraction
- class-specific expert branches for membrane, kinase, TF, and other proteins
- an autoencoder bottleneck for latent dimensional reduction
- task-specific supervised heads instead of a single shared scalar head

Implemented files:

- `Daedalus/phase_c/models/class_specific_biochem_autoencoder.py`
- `Daedalus/phase_c/scripts/build_biochem_delta_matrix.py`
- `Daedalus/phase_c/scripts/build_multitask_instances.py`
- `Daedalus/phase_c/scripts/train_biochem_multitask.py`
- `Daedalus/phase_c/scripts/run_model_comparison.py`
- `Daedalus/phase_c/scripts/run_feature_ablation.py`

### Current Phase C benchmark

The current completed Phase C benchmark uses the richer biochemical subset:

- `40,000` biochemical pair rows
- `75,337` multitask instances

This benchmark is sufficient to rank currently implemented model families and to
run feature-block ablations.

Benchmark commands:

```bash
Daedalus/.venv/bin/python Daedalus/phase_c/scripts/build_multitask_instances.py
Daedalus/.venv/bin/python Daedalus/phase_c/scripts/train_biochem_multitask.py \
  --max-epochs 40 \
  --patience 8 \
  --batch-size 1024
Daedalus/.venv/bin/python Daedalus/phase_c/scripts/run_model_comparison.py
Daedalus/.venv/bin/python Daedalus/phase_c/scripts/run_feature_ablation.py
```

Current classifier ranking on the Phase C subset:

- `phase_c_xgboost`
  - `global` AP `0.878873`, AUROC `0.914486`
  - `membrane` AP `0.886331`, AUROC `0.903270`
  - `surface` AP `0.878424`, AUROC `0.895676`
  - `kinase` AP `0.881222`, AUROC `0.919501`
  - `transcription_factor` AP `0.874922`, AUROC `0.902396`

- `phase_c_logistic`
  - `global` AP `0.821485`, AUROC `0.869970`
  - `membrane` AP `0.812549`, AUROC `0.850795`
  - `surface` AP `0.806671`, AUROC `0.845223`
  - `kinase` AP `0.761381`, AUROC `0.815748`
  - `transcription_factor` AP `0.827372`, AUROC `0.865278`

- `phase_c_autoencoder`
  - `global` AP `0.799552`, AUROC `0.853325`
  - `membrane` AP `0.808293`, AUROC `0.836624`
  - `surface` AP `0.803444`, AUROC `0.831492`
  - `kinase` AP `0.783835`, AUROC `0.847317`
  - `transcription_factor` AP `0.845512`, AUROC `0.881434`

Interpretation:

- the richer biochemical feature layer clearly carries useful signal
- `XGBoost` is currently the best efficient classifier with available tools
- the autoencoder is not yet the best general classifier on this matrix
- the autoencoder remains worth keeping as a deep-learning comparator and latent
  bottleneck path, especially for TF-focused follow-up work

### DeepImmuno-style benchmark extension

To explicitly mirror the strongest methodological elements from DeepImmuno, a
second Phase C benchmark path is now implemented with:

- systematic ML/DL comparison
- PCA dimensional reduction on the high-dimensional physicochemical blocks
- CNN and CNN+autoencoder classifiers

Implemented files:

- `Daedalus/phase_c/models/deepimmuno_style.py`
- `Daedalus/phase_c/scripts/run_deepimmuno_style_benchmark.py`

Tested models:

- `elasticnet`
- `knn`
- `svm_rbf`
- `random_forest`
- `adaboost`
- `hist_gradient_boosting`
- `residual_mlp_autoencoder`
- `deepimmuno_cnn`
- `deepimmuno_cnn_autoencoder`

Completed benchmark:

- `global`
- `membrane`
- `surface`
- `kinase`
- `transcription_factor`

Best model by task:

- `global`: `hist_gradient_boosting` AP `0.878693`, AUROC `0.916296`
- `membrane`: `hist_gradient_boosting` AP `0.883712`, AUROC `0.906119`
- `surface`: `hist_gradient_boosting` AP `0.886910`, AUROC `0.903786`
- `kinase`: `hist_gradient_boosting` AP `0.874902`, AUROC `0.909127`
- `transcription_factor`: `hist_gradient_boosting` AP `0.870329`,
  AUROC `0.899536`

Strongest deep model by task:

- `global`: `deepimmuno_cnn_autoencoder` AP `0.826830`, AUROC `0.878698`
- `membrane`: `deepimmuno_cnn` AP `0.832216`, AUROC `0.862649`
- `surface`: `deepimmuno_cnn` AP `0.814196`, AUROC `0.852984`
- `kinase`: `deepimmuno_cnn_autoencoder` AP `0.796933`, AUROC `0.857060`
- `transcription_factor`: `deepimmuno_cnn_autoencoder` AP `0.850442`,
  AUROC `0.891157`

PCA reduction across tasks:

- `full_sequence_physchem`: `74 -> 23-24`, explained variance `0.9455-0.9509`
- `terminal_sequence_physchem`: `120 -> 24`, explained variance `0.8852-0.8958`
- `typed_overlap`: `51 -> 16`, explained variance `0.7188-0.7747`

Interpretation:

- the DeepImmuno-style PCA reduction step is implemented and operational
- the deep CNN models are competitive but not currently the best classifiers on
  this problem
- the strongest current deep-learning option overall in this benchmark family
  is the CNN+autoencoder, although the plain CNN is stronger on the membrane
  and surface tasks
- the strongest overall model in the completed DeepImmuno-style benchmark is a
  boosted/tree ensemble rather than a deep model
- the stable aggregate benchmark artifact is:
  - `Daedalus/phase_c/checkpoints/deepimmuno_style_benchmark_all_tasks.md`

### Phase C feature-block ablation

The first ablation sweep gives a clear feature ranking.

Useful blocks:

- `base_reference_delta`
- `segment_change`

Weak or currently noisy blocks:

- `typed_overlap`
- `terminal_sequence_physchem`
- `full_sequence_physchem`
- `motifs` (mostly neutral to weak)

Interpretation:

- the model is learning most of its signal from structured reference-vs-alt
  context plus explicit changed-segment features
- broad full-sequence physicochemical descriptors are too diffuse in the current
  implementation
- the current typed-overlap feature map is biologically plausible but too noisy
  and should be tightened around higher-value functional grammars rather than
  broad feature buckets

## Promotion to supervised task claims

To convert those channel outputs into real model claims, Phase B now includes a
supervised benchmark scaffold:

- `Daedalus/phase_b/config/task_registry.json`
- `Daedalus/phase_b/schemas/supervised_benchmark_record.schema.json`
- `Daedalus/phase_b/scripts/build_supervised_benchmark_table.py`

The supervised task registry currently defines:

- `ppi_impact`
- `pdi_impact`
- `localization_shift`
- `tf_dna_binding`
- `tf_transcriptional_activity`
- `tf_perturbation_response`
- `kinase_signaling`
- `signal_retention`
- `tm_insertion`
- `tm_fold`

This is the correct next step scientifically:

1. map external assays onto frozen reference-vs-alternative pairs
2. train these as explicit heads rather than treating them as evidence-only
   channels
3. benchmark each head on external held-out datasets
4. only then elevate them to manuscript-level functional claims

The holdout mechanism is now explicit:

- external assay rows inherit the frozen Phase B `train` / `val` / `test`
  assignment from the matched reference-delta pairs
- supervised baseline evaluation is run through:
  - `Daedalus/phase_b/scripts/run_supervised_holdout_baselines.py`

This avoids the common failure mode where external validation is reported on
matched assay rows without preserving the original gene-disjoint holdout split.

### Current promoted-task proxy benchmark state

The promoted-task scaffold has now been exercised using a proxy supervised
benchmark assembled from the frozen Phase A pair table. This is not the same as
true external validation, but it is enough to establish which promoted tasks
look learnable on the current reference-delta matrix.

Current holdout test metrics:

- `signal_retention`: AP `0.459793`, AUROC `0.805449`
- `tm_insertion`: AP `0.433830`, AUROC `0.838737`
- `tm_fold`: AP `0.357609`, AUROC `0.823234`
- `ppi_impact`: AP `0.168298`, AUROC `0.636388`
- `pdi_impact`: AP `0.161066`, AUROC `0.521164`
- `tf_dna_binding`: AP `0.161066`, AUROC `0.521164`
- `tf_transcriptional_activity`: AP `0.161066`, AUROC `0.521164`
- `tf_perturbation_response`: AP `0.161066`, AUROC `0.521164`
- `localization_shift`: AP `0.159922`, AUROC `0.525888`

Current interpretation:

- membrane-focused promoted tasks are the most credible immediate extension of
  the current model family
- `signal_retention`, `tm_insertion`, and `tm_fold` are already in a range that
  justifies keeping them as active supervised heads
- `ppi_impact` is feasible but not yet robust enough to claim as a strong
  biological task without external assays
- the TF-related promoted tasks are not yet convincing under proxy supervision
  and need the real hESC assay tables to become meaningful
- `localization_shift` is unstable across splits and should be treated as
  unresolved until external localization supervision is loaded
- `kinase_signaling` is not yet benchmarkable because the current proxy label
  construction collapsed to a single class

The concrete output files are:

- `Daedalus/phase_b/data/processed/supervised_holdout_metrics.tsv`
- `Daedalus/phase_b/data/processed/supervised_holdout_report.md`

## Phase D: explicit TM and TF objective heads

Phase D replaces the broad membrane / TF readout with explicit task heads.

Membrane heads:

- `secretory_pathway_compatibility`
- `membrane_insertion_topology`
- `extracellular_altered_segment_accessibility`
- `folding_stability_qc_escape`
- `cell_surface_localization`
- `antibody_targetability`
- aggregate:
  - `tm_stable_functional`
  - `tm_non_stable_non_functional`

TF heads:

- `nuclear_localization`
- `dbd_integrity_specificity_shift`
- `cofactor_ppi_rewiring`
- `activation_repression_competence`
- `dominant_negative_neomorphic_behavior`
- aggregate:
  - `overall_regulatory_functionality`
  - `tf_non_stable_non_functional`

Phase D is now implemented as:

- objective scoring:
  - `Daedalus/phase_d/scripts/build_phase_d_scores.py`
- benchmark-instance generation:
  - `Daedalus/phase_d/scripts/build_phase_d_benchmark_instances.py`
- family-aware multitask deep model:
  - `Daedalus/phase_d/scripts/train_phase_d_multitask.py`
- model comparison:
  - `Daedalus/phase_d/scripts/run_phase_d_model_comparison.py`

Current benchmark state:

- the proxy benchmark runs end-to-end
- the family-aware multitask model trains successfully
- logistic / hist-gradient-boosting / xgboost / multitask outputs are all
  generated

Current limitation:

- many Phase D proxy tasks are extremely one-sided after thresholding
- several test tasks have zero negatives or zero positives
- many current proxy metrics are therefore trivially perfect and should not be
  interpreted as biological validation

So the current Phase D benchmark is valid as:

- pipeline validation
- model-family plumbing validation

It is not yet valid as:

- final biological performance evidence

One important exception now exists for the membrane branch:

- `cell_surface_localization` includes explicit curated TM-negative held-out
  examples derived from UniProt isoforms that retain at least one TM domain but
  lose `Cell membrane` localization

Those labels are built in:

- `Daedalus/phase_a/scripts/build_tm_negative_pair_labels.py`

and saved to:

- `Daedalus/phase_a/data/interim/tm_negatives/tm_negative_pair_labels.tsv`

They are injected into the Phase D benchmark as explicit negative overrides for
`cell_surface_localization`, giving the current test set:

- `533` rows total
- `520` positives
- `13` negatives
- `7` curated held-out negatives

Current task-level benchmark after those injected negatives:

- `phase_d_logistic`
  - AP `0.999370`
  - AUROC `0.975148`
  - balanced accuracy `0.790385`
- `phase_d_hist_gradient_boosting`
  - AP `0.999830`
  - AUROC `0.993195`
  - balanced accuracy `0.692308`
- `phase_d_xgboost`
  - AP `0.999802`
  - AUROC `0.992160`
  - balanced accuracy `0.692308`
- `phase_d_family_multitask`
  - AP `0.999322`
  - AUROC `0.973669`
  - balanced accuracy `0.729808`

This is still a highly imbalanced test, but it is materially better than the
pure objective-threshold proxy because it now includes external curated
non-surface TM isoforms.

Primary Phase D deliverables now on disk:

- `Daedalus/phase_d/checkpoints/phase_d_multitask_metrics.json`
- `Daedalus/phase_d/checkpoints/phase_d_multitask_report.md`
- `Daedalus/phase_d/checkpoints/phase_d_model_comparison_metrics.tsv`
- `Daedalus/phase_d/checkpoints/phase_d_model_comparison_report.md`
- `Daedalus/phase_d/TRAINING_AND_BENCHMARK_SPEC.md`

## Cluster / GPU execution

This project does not have to be containerized first in order to run on a GPU
cluster. If the cluster allows direct Python environments on GPU nodes, you can
run it from a venv or conda environment.

Containerization becomes useful when:

- the cluster requires Apptainer/Singularity for reproducible jobs
- dependency portability is a problem
- you want frozen CUDA / Python / package versions

For Singularity/Apptainer, the usual path is:

1. create a Docker/OCI image or definition
2. convert or pull it into an Apptainer/Singularity image
3. run the training jobs against that image on the cluster

So:

- `required first`: no
- `recommended for reproducibility and cluster portability`: yes

## Container artifacts

This project now includes:

- `Daedalus/container/Dockerfile`
- `Daedalus/container/apptainer.def`
- `Daedalus/container/requirements.txt`
- `Daedalus/container/run_apptainer.sh`

Recommended cluster path:

```bash
docker build -f Daedalus/container/Dockerfile -t daedalus:latest .
apptainer build daedalus.sif Daedalus/container/apptainer.def
apptainer exec --nv --bind /path/to/components:/workspace daedalus.sif bash
```

Convenience wrapper:

```bash
Daedalus/container/run_apptainer.sh /path/to/daedalus.sif
Daedalus/container/run_apptainer.sh /path/to/daedalus.sif \
  python Daedalus/phase_b/scripts/train_reference_delta_multitask.py
```

## Comparison with Proteus

See:

- `Daedalus/COMPARISON_WITH_PROTEUS.md`

The short version is:

- `Daedalus` is the stronger implemented benchmark substrate right now, with
  completed held-out baseline metrics
- `Proteus` is the richer architectural proposal right now, with pretrained
  multimodal encoders and a clearer end-state model design
- `Proteus` documentation explicitly states supervised fine-tuning is not yet
  completed, so it currently lacks comparable held-out benchmark outputs
