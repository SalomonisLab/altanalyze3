# data_ingest

UniProt-only data preparation: downloads the reviewed human + mouse UniProt
release, reconstructs every alternative-isoform protein sequence from the
canonical sequence plus VSP (Variant Sequence Patch) features, and emits the
per-isoform feature tables that `benchmark/` consumes.

## Inputs

A network connection on first run; subsequent runs read cached files in
`data/raw/`.

## Outputs

| File | Producer | Description |
|---|---|---|
| `data/raw/uniprot_reviewed_human_mouse.jsonl.gz` | `download_uniprot.py` | Reviewed UniProt records (human + mouse) |
| `data/interim/uniprot_entries.tsv` | `download_uniprot.py` | Per-accession metadata (organism, gene name, subcellular_locations, interaction partners) |
| `data/interim/uniprot_features.tsv` | `download_uniprot.py` | Per-accession feature annotations (Transmembrane, Signal, Glycosylation, Disulfide bond, Topological domain, Region, Domain, Motif, etc.) |
| `data/interim/tm_negatives/uniprot_isoform_catalog_full.tsv` | `download_uniprot.py` | Per-isoform Cell Membrane / Transmembrane catalog used by the benchmark builder |
| `data/interim/uniprot_isoform_proteins.tsv` | `build_isoform_proteins.py` | One row per catalog isoform with reconstructed full protein sequence (canonical sequence + VSP patches applied) |
| `data/interim/uniprot_isoform_features.tsv` | `build_isoform_features.py` | 74 sequence-derived + UniProt-feature counts per isoform |
| `data/interim/uniprot_isoform_trafficking_features.tsv` | `build_trafficking_features.py` | 48 trafficking-biology descriptors (signal-peptide integrity, TM topology parity, ER motif gain/loss, disulfide-partner orphaning, PDZ/GPI/polybasic C-terminal status) |
| `data/interim/uniprot_isoform_esm2_t30_embeddings.parquet` | `build_isoform_esm_embeddings.py` | 640-dim ESM2-150M (t30) mean-pooled embeddings for every reconstructed isoform |
| `data/interim/uniprot_isoform_esm2_t33_embeddings.parquet` | (optional) re-run with `--model-name facebook/esm2_t33_650M_UR50D` | 1,280-dim ESM2-650M embeddings; benchmarked but not used in the production model |
| `data/interim/uniprot_canonical_esm2_t30_embeddings.parquet` | `build_canonical_esm_embeddings.py` | ESM2 t30 embeddings for all 37,683 canonical UniProt entries (training set for the localization classifier) |
| `data/interim/uniprot_canonical_pm_evidence.tsv` | one-time helper (see archive `intermediate_data.tar.gz`) | Per-accession booleans for the strict-canonical-PM label filter (whether the canonical UniProt entry has unambiguous Cell Membrane vs ambiguous multi-compartment annotation) |
| `data/interim/uniprot_localization_classifier.pt` | `build_localization_classifier.py` | DeepLoc-style 10-class subcellular-localization head trained on the 37,683 canonical embeddings |
| `data/interim/uniprot_isoform_localization_probs.tsv` | `build_localization_classifier.py` | Per-isoform 10-class softmax probabilities (PM, Cytoplasm, Nucleus, ER, Golgi, Mitochondrion, Lysosome, Secreted, Other-membrane, Other) |

## Reproduction

```bash
cd components/Daedalus
python data_ingest/scripts/download_uniprot.py
python data_ingest/scripts/build_isoform_proteins.py
python data_ingest/scripts/build_isoform_features.py
python data_ingest/scripts/build_trafficking_features.py
python data_ingest/scripts/build_isoform_esm_embeddings.py
python data_ingest/scripts/build_canonical_esm_embeddings.py
python data_ingest/scripts/build_localization_classifier.py
```

ESM2 extraction times on CPU: t30 ~27 minutes for the 5,882 isoforms and
~3 hours for the 37,683 canonical sequences (overnight); t33 (optional)
~9 hours for the isoforms.

## Compute

CPU-only by design. `transformers` and `huggingface_hub` are required for
the ESM2 extraction scripts.
