# Daedalus Query Workflow

This workflow scores a novel alternative isoform relative to an automatically
selected canonical reference transcript for the same gene.

## Automatic reference selection

If you do not provide a reference transcript, Daedalus selects one
automatically using:

1. `MANE Select`
2. `APPRIS PRINCIPAL:1`
3. other `APPRIS PRINCIPAL:*`
4. `UNIPROT_SUPPORTED_LONGEST`
5. `LONGEST_PROTEIN`
6. `LONGEST_TRANSCRIPT`

This matches the ranking logic used in Phase A pair construction.

## Query inputs

Required:

- `--species`
- `--gene-id` or `--gene-name`
- `--alt-protein-seq` or `--alt-protein-fasta`

Optional:

- `--reference-transcript-id`
- `--alt-transcript-length`
- `--alternative-transcript-id`
- `--include-biogrid-partners` to enumerate specific BioGRID partners for the
  reference gene

## Current outputs

The query workflow returns:

- selected reference transcript and selection reason
- current trained task probabilities:
  - `global`
  - `membrane`
  - `surface`
  - `kinase`
  - `transcription_factor`
- explicit channel evidence for:
  - `ppi_impact_risk_score`
  - `pdi_impact_risk_score`
  - `dpi_impact_risk_score`
  - `localization_shift_risk_score`
  - `tf_dna_binding_retained_score`
  - `kinase_signaling_retained_score`
  - `signal_peptide_retained_score`
  - `tm_insertion_support_score`
  - `tm_fold_support_score`
- known BioGRID partners for the reference gene
- predicted TM segments in the query protein

Specific BioGRID partner enumeration is optional because it requires scanning
the raw BioGRID archive. The default fast path leaves that list empty but still
uses gene-level BioGRID counts in the risk evidence.

## Interpretation boundary

Only the five Phase B task probabilities are currently learned from the trained
baseline model.

The additional PPI/PDI/localization/signal/TM channel outputs are currently
structured evidence scores derived from:

- the trained task probabilities
- reference UniProt/HPA/BioGRID context
- alternative sequence-derived TM/signal/glyco/cysteine features

So the workflow now supports these outputs explicitly, but they are not yet
separately trained heads with external benchmark reports.

## Example

```bash
Daedalus/.venv/bin/python Daedalus/phase_b/scripts/predict_query_isoform.py \
  --species human \
  --gene-name EGFR \
  --alt-protein-fasta novel_egfr_isoform.fa
```
