# benchmark

Assembles the (canonical, alternative)-isoform-pair feature parquet that
the training stage consumes. Pair construction: every non-canonical
isoform in the catalog is paired against its gene's canonical reference
isoform, with reference, alternative, and delta-from-reference values
computed for each numeric feature.

## Inputs

Reads files in `../data_ingest/data/interim/` produced by the data-ingest
stage. Required:
- `uniprot_isoform_features.tsv`
- `uniprot_isoform_trafficking_features.tsv`
- `uniprot_isoform_esm2_t30_embeddings.parquet`
- `uniprot_isoform_localization_probs.tsv`
- `uniprot_canonical_pm_evidence.tsv`
- `tm_negatives/uniprot_isoform_catalog_full.tsv`

## Outputs

| File | Description |
|---|---|
| `data/uniprot_isoform_benchmark_instances.parquet` | The benchmark itself: one row per (canonical, alt) pair with reference/alternative/delta features (1,486 numeric features in the production strict-canonical configuration) |
| `data/uniprot_isoform_benchmark_instances.tsv` | Same content as the parquet, TSV form |
| `data/uniprot_isoform_benchmark_schema.json` | Feature-column inventory + label/group column names + category encoding metadata |
| `data/uniprot_isoform_benchmark_summary.tsv` | Per-split label counts |

## Production filter

`--strict-canonical-pm` excludes ~419 genes whose canonical UniProt entry
is annotated at both Cell Membrane and another major compartment
(Cytoplasm / Nucleus / Mitochondrion / Lysosome). The filter eliminates
the dominant source of label noise and lifts AUROC from 0.81 to 0.91
relative to the unfiltered version.

```bash
cd components/Daedalus
python benchmark/scripts/build_benchmark.py --strict-canonical-pm
```

Without that flag the benchmark emits all 3,982 pairs (3,703 PM-positive
plus 279 cleaned non-surface negatives across 1,815 genes); the
production stratification reduces this to 3,022 pairs (2,839 + 183) over
1,396 genes with cleaner labels.

## Other CLI flags

| Flag | Effect |
|---|---|
| `--keep-organellar-negatives` | Include alternative isoforms on non-PM organellar membranes as negatives. Off by default — these isoforms are TM proteins on intracellular membranes and are indistinguishable from PM proteins by sequence-derived features. |
| `--keep-missing-negatives` | Include alternative isoforms with no `isoform_locations` annotation as negatives. Off by default. |
| `--keep-pm-like-negatives` | Include alternative isoforms tagged Cell-surface or Cell-junction as negatives. Off by default — almost certainly mislabelled. |

## Pair semantics

Each non-canonical isoform appears in exactly one pair with its gene's
single UniProt-displayed canonical reference. A gene with five
alternatives produces five rows. The canonical itself is never the
"alternative" side. Reference features are constant within a gene;
delta features (alternative − reference) carry the discriminating
signal.
