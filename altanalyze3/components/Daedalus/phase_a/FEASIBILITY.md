# Daedalus Phase A Feasibility

## Baseline setup

Phase A now includes:

- leakage-safe gene-level train/validation/test splits
- weakly supervised preservation seeds
  - positive: shared UniProt accession support across reference and alternative isoforms
  - negative: reference-preferred APPRIS/MANE contrast seeds
- family-specific benchmark subsets for membrane, surface, kinase, and transcription-factor genes

Baselines were run with pure Python only. No external ML stack was required.

## Important benchmark caveat

`logistic_all_features` reaches perfect performance because it still uses UniProt-derived alternative-function features that are too close to the weak labels. It is useful as an upper bound, not as a realistic feasibility estimate.

The relevant feasibility estimate is `logistic_structural`, which excludes direct UniProt-derived alternative function-retention features and uses structural/reference-context proxies instead.

## Non-leaky feasibility results

Best non-leaky test-set results:

- `task_global_seed`
  - `n=12,858`
  - positives `3,467`, negatives `9,391`
  - best model: `logistic_structural`
  - `AP=0.773746`, `AUROC=0.790579`

- `task_membrane_seed`
  - `n=4,246`
  - positives `1,216`, negatives `3,030`
  - best model: `logistic_structural`
  - `AP=0.799868`, `AUROC=0.828088`

- `task_kinase_seed`
  - `n=1,682`
  - positives `467`, negatives `1,215`
  - best model: `logistic_structural`
  - `AP=0.849596`, `AUROC=0.865343`

- `task_tf_seed`
  - `n=2,255`
  - positives `645`, negatives `1,610`
  - best model: `logistic_structural`
  - `AP=0.825644`, `AUROC=0.843121`

- `task_surface_seed`
  - `n=226`
  - positives `32`, negatives `194`
  - best model: `logistic_structural`
  - `AP=0.393975`, `AUROC=0.515303`

## Interpretation

The Phase A corpus is already strong enough to support meaningful baseline modeling for:

- global preservation/deviation
- membrane-associated isoform preservation
- kinase-focused preservation
- transcription-factor-focused preservation

The surface/topology subset is not yet strong enough. It likely needs:

- more direct extracellular/topology supervision
- better transmembrane and signal-peptide labels
- possibly InterPro / topology-specific resources in addition to HPA

## Practical conclusion

With one `H100`, the next stage is feasible.

The recommended next step is:

1. keep the current Phase A corpus frozen
2. move to a stronger reference-conditioned baseline model
3. prioritize `global`, `membrane`, `kinase`, and `TF` tasks first
4. defer strong claims about `surface` predictions until topology supervision improves
