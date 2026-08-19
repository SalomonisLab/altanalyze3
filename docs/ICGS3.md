# ICGS3

## Overview

ICGS3 clusters single-cell data with sparse graph construction, intelligent downsampling, UDON
sparse NMF, MarkerFinder, linear-SVM reclassification, UMAP and GO-Elite summaries. The module
lives at `altanalyze3/components/clustering/ICGS.py`.

Run it through its CLI:

```bash
python3.11 -m altanalyze3.components.clustering.ICGS --input <matrix> --output-dir <folder>
```

The `--modality` flag accepts `rna`, `adt`, `grn`, `metabolite`, `lipid` and `psi`. RNA is the
default.

## Recommended ADT workflow

We recommend this command for antibody-derived tag (ADT) panels:

```bash
python3.11 -m altanalyze3.components.clustering.ICGS \
  --input <adt_matrix.h5ad> \
  --output-dir <folder> \
  --modality adt \
  --species <Hs|Mm> \
  --umap-covariates <obs columns> \
  --heatmap-covariates <obs columns> \
  --heatmap-goelite-terms
```

Two points carry the recommendation:

1. **Pass `--modality adt` and let ICGS3 pick the marker gate.** ICGS3 then resolves the gate to
   `--marker-rho 0.15` and `--marker-min-per-cluster 1`. An explicit flag still wins.
2. **Omit `--input-normalized` when the matrix holds linear values.** ICGS3 then applies
   `normalize_total(1e4)` and `log1p`. Pass the flag only when the matrix already holds
   log-scale values.

Point 2 matters for a TotalVI-denoised panel. TotalVI writes linear values, so the log step still
applies. Check `uns['expression_scale']` when the exporter records it.

## Why we recommend it

Internal evaluation, 2026-08-17. Two mouse CITE-seq panels from the Grimes Thymus study supplied
the test data: ADT112 (1,592 cells x 112 antibodies) and ADT195 (2,140 cells x 195 antibodies).
Both matrices hold linear TotalVI-denoised values. We crossed two scale options against two marker
gates and ran all four cells of the design on both panels.

| Panel | Scale | Marker gate | Final clusters | Cells retained |
|---|---|---|---|---|
| ADT112 | log | rho 0.3, 2 markers | 3 | 1,478/1,592 = 92.8% |
| ADT112 | log | rho 0.15, 1 marker | **4** | **1,472/1,592 = 92.5%** |
| ADT112 | linear | rho 0.3, 2 markers | 2 | 1,592/1,592 = 100% |
| ADT112 | linear | rho 0.15, 1 marker | 3 | 873/1,592 = 54.8% |
| ADT195 | log | rho 0.3, 2 markers | 2 | 2,140/2,140 = 100% |
| ADT195 | log | rho 0.15, 1 marker | 2 | 2,140/2,140 = 100% |
| ADT195 | linear | rho 0.3, 2 markers | 2 | 2,140/2,140 = 100% |
| ADT195 | linear | rho 0.15, 1 marker | 4 | 1,844/2,140 = 86.2% |

The bold row won. The reviewer picked it on marker-heatmap quality. The log scale also protects
cell retention: the linear ADT112 run at the same gate dropped 45.2% of the cells.

The log step tames a long right tail. ADT112 holds a median value of 2.45 and a maximum of 489.2.
Pearson MarkerFinder and NMF then answer to a few bright antibodies. Per-cell totals also spread
280-fold across ADT112 cells, from 11.9 to 3,329.3, so the `normalize_total` step equalizes
staining depth.

## Resolved defaults by modality

`apply_modality_defaults` (`ICGS.py:424`) resolves these before a run starts. An explicit CLI flag
always overrides the resolved value.

| Setting | rna | adt | grn, metabolite, lipid, psi |
|---|---|---|---|
| `--marker-rho` | 0.3 | 0.15 | 0.3 |
| `--marker-min-per-cluster` | 2 | 1 | 2 |
| `--min-genes` | 500 | 0 | 0 |
| `--min-counts` | 1000 | 0 | 0 |
| `--mito-percent` | 30.0 | none | none |
| `--n-top-features` | 3000 | 0 | 0 |

The marker gate keeps a cluster only when that cluster owns at least `--marker-min-per-cluster`
markers at Pearson rho of `--marker-rho` or above. The gate runs at `ICGS.py:1517` and
`ICGS.py:1566`.

`--marker-rho` and `--marker-min-per-cluster` carry no parser default. They stay `None` until
`apply_modality_defaults` reads the modality. An earlier version detected an unset value by
comparing it against the literal default of the day. Those defaults later changed, the ADT branch
became dead code, and every ADT run then used the RNA gate. The test
`test_marker_gate_defaults_resolve_by_modality` in `altanalyze3/components/clustering/test_ICGS.py`
guards the current behaviour.

## Open limitation: the NMF rank estimator on small panels

`ICGS.py:1340` calls `determine_nmf_ranks(expr_sampled, small_feature=False)` for every modality.
The comment at `altanalyze3/components/udon/nmf.py:62-65` names the risk: few features and many
cells inflate the Tracy-Widom boundary, so only the dominant global component clears it, and the
estimator returns rank 2.

Measured on the same two panels, 2026-08-17:

| Panel | Scale | `small_feature=False` (what runs today) | `small_feature=True` |
|---|---|---|---|
| ADT112 | log | 4 | 8 |
| ADT112 | linear | 6 | 4 |
| ADT195 | log | 2 | 8 |
| ADT195 | linear | 4 | 3 |

The ADT195 log run returned rank 2, so it produced 2 clusters before the marker gate ever ran. The
rank sets a hard ceiling on the cluster count.

`altanalyze3/components/udon/clustering_wrapper.py:31` passes the flag through, and ICGS3 does not
expose it. Status: OPEN. Override the rank with `--nmf-k` when a panel looks under-clustered.

## Validation of the ADT gate change

Command with no marker flags, ADT112, log scale: the run reproduced the approved run at
1,472/1,472 identical cluster labels (100.00%), in identical barcode order, with identical cluster
sizes. The files `icgs3_marker_heatmap_markers.tsv`, `icgs3_marker_heatmap_fold_matrix.tsv` and
`icgs3_markers.tsv` matched by MD5. ADT195 matched at 2,140/2,140 labels.

RNA regression, same input and flags as the run before the change: 3,638/3,638 identical cluster
labels (100.00%), identical sizes C1=1,252, C2=514, C3=1,413, C4=459.

Unit tests: 21 of 21 passed in `altanalyze3/components/clustering/test_ICGS.py`.

Evaluation artifacts sit under
`/Users/saljh8/Dropbox/Collaborations/Grimes/Thymus/ICGS3/`, with one `run_icgs3.sh`,
`icgs3_config.json` and `logs/` per run, plus
`/Users/saljh8/Dropbox/Collaborations/Grimes/Thymus/ICGS3/README_ICGS3_20260817.md`.

## Provenance and maintenance

| Fact | Re-verification command | Still true when |
|---|---|---|
| Gate resolves by modality | `pytest altanalyze3/components/clustering/test_ICGS.py -k marker_gate` | The test passes |
| Resolved default table | `sed -n '424,450p' altanalyze3/components/clustering/ICGS.py` | The values match the table |
| Rank estimator limitation | `grep -n "small_feature" altanalyze3/components/clustering/ICGS.py` | The call still reads `small_feature=False` |

Evaluation date: 2026-08-17. Panels evaluated: 2 of 2 available ADT cocktails. Species evaluated:
mouse only. We did not evaluate the `grn`, `metabolite`, `lipid` or `psi` modalities, so their gate
stays at the RNA values.
