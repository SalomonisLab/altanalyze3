# Daedalus Phase B Status

Phase B is now in a usable baseline state.

## Completed

- built the reference-conditioned feature table
- expanded it into a multitask instance table on the frozen Phase A splits
- trained a shared torch multitask baseline
- added a leakage-resistant sklearn baseline for cross-checking
- verified that train/val/test splits are disjoint at the gene level

## Main outputs

- `data/processed/reference_delta_matrix.parquet`
- `data/processed/multitask_instances.parquet`
- `data/processed/reference_delta_sklearn_metrics.tsv`
- `checkpoints/reference_delta_multitask_metrics.json`
- `checkpoints/reference_delta_multitask_report.md`
- `checkpoints/reference_delta_sklearn_report.md`

## Current interpretation

- `global`, `membrane`, `kinase`, and `transcription_factor` are all learnable on the current corpus
- `surface` is now materially improved after adding topology, signal, PPI, PTM, and cysteine/glycosylation context
- the current Phase B matrix is de-leaked:
  - reference-known UniProt/HPA/BioGRID annotations are retained
  - alternative isoforms contribute sequence-derived features and coarse length deltas
  - direct alternative UniProt region annotations are excluded from the modeled matrix
- the Phase B torch baseline remains stronger than the sklearn baseline on the same cleaned weak-label regime

## Test summary

Torch multitask baseline:

- `global` AP `0.840682`, AUROC `0.882048`
- `membrane` AP `0.839828`, AUROC `0.875073`
- `kinase` AP `0.866310`, AUROC `0.903348`
- `transcription_factor` AP `0.858565`, AUROC `0.884078`
- `surface` AP `0.837000`, AUROC `0.871491`

sklearn logistic baseline:

- `global` AP `0.810107`, AUROC `0.845348`
- `membrane` AP `0.809486`, AUROC `0.843945`
- `kinase` AP `0.853772`, AUROC `0.874458`
- `transcription_factor` AP `0.827434`, AUROC `0.842662`
- `surface` AP `0.802420`, AUROC `0.835832`

## Next validation layer

The TF branch will be validated against external hESC assay outputs:

- TF isoform PPI
- TF isoform PDI
- TF isoform transcriptional activation
- TF isoform localization
- TF isoform perturbation / bootstrap DEGs

See `VALIDATION.md`.
