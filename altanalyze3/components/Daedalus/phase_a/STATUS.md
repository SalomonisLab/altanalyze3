# Daedalus Phase A Status

Current Phase A assets are sufficient to start leakage-safe baseline modeling on a single `H100`.

## Core resources acquired

- Human Protein Atlas subcellular localization
- ClinVar variant summary
- UniProt reviewed human and mouse entries
- GENCODE human `v49` and mouse `vM38`
- APPRIS human `GRCh38` and mouse `GRCm39`
- MANE human `GRCh38 v1.5`

## Current interim tables

- `hpa_localization.tsv`
- `clinvar_splice_pathogenic.tsv`
- `uniprot_entries.tsv`
- `uniprot_isoforms.tsv`
- `uniprot_features.tsv`
- `uniprot_functional_priors.tsv`
- `gencode_transcript_reference.tsv`
- `gencode_protein_reference.tsv`
- `uniprot_gencode_map.tsv`
- `appris_principal_isoforms.tsv`
- `mane_transcripts.tsv`
- `gene_supervision_catalog.tsv`
- `transcript_supervision_catalog.tsv`
- `isoform_pair_candidates.tsv`
- `isoform_pair_candidates.with_splits.tsv`
- `priority_pair_subsets.tsv`

## Current scale

- Transcript supervision rows: `785,691`
- Human MANE Select transcripts recovered: `19,297`
- Isoform pair candidates: `257,851`
- Benchmark genes with pair candidates: `32,404`

## Pair family coverage

- Membrane pairs: `76,403`
- Surface-oriented membrane pairs: `2,290`
- Kinase pairs: `33,994`
- Transcription-factor pairs: `42,091`
- ClinVar splice-sensitive gene pairs: `59,899`
- Shared-UniProt preserved seeds: `35,214`
- APPRIS/MANE contrast seeds: `93,223`
- High-confidence pairs with UniProt support on both isoforms: `17,417`

## Current split sizes

- Train pairs: `205,485`
- Validation pairs: `26,352`
- Test pairs: `26,014`

Validation and test counts remain large enough for global baselines and family-specific baselines.

## Current reference-anchor policy

- Human: prefer `MANE Select`, then `APPRIS`, then UniProt-supported longest protein, then longest protein
- Mouse: prefer `APPRIS`, then UniProt-supported longest protein, then longest protein

Transcript and gene joins ignore Ensembl version suffixes where needed for APPRIS, MANE, and HPA integration.
