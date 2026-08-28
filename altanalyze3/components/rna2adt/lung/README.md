# Human lung rna2adt — per-protein ElasticNet on a protein-coding empirical whitelist

This package predicts CITE-seq ADT values from RNA for human lung datasets. It
uses the architecture and the training protocol of the bone marrow model in
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/` and
the murine port in
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/mouse/`.
Only the panel, the curated map, the training atlas and the bundle differ.

## Files

All paths below start at `/`.

| Path | What it holds |
|---|---|
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/lung/adt_rna_map.py` | Curated ADT to HGNC map for the 56-antibody lung panel. |
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/lung/configs/adt_rna_map.tsv` | The curated map as a TSV. |
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/lung/data.py` | Loads the ADT table, loads the RNA `.h5ad`, and joins the two cell-name conventions. |
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/lung/build_whitelist.py` | Generates the empirical feature whitelist. |
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/lung/configs/empirical_whitelist.tsv` | The whitelist the shipped bundle uses. 56 ADTs, 81 RNA genes. |
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/lung/model.py` | The `HumanLungRna2AdtModel` class. |
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/lung/train_lung.py` | Trainer. |
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/lung/rna2adt_hs_lung_bundle.pkl` | The trained bundle. |
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/lung/rna2adt_hs_lung_per_adt_metrics_cells.tsv` | Per-ADT scores on the cell holdout. |
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/lung/validate_port.py` | Proves the ported whitelist code reproduces the mouse whitelist. |
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/lung/validate_bundle.py` | Loads the bundle through the inference API and drives the cellHarmony-web builder. |
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/lung/compare_target_scale.py` | Compares log1p targets against linear targets. |

Run logs and run artifacts live under
`/Users/saljh8/Dropbox/Transfer/COVID-TotalVI/rna2adt_lung/`.

## Source data

The RNA lives in
`/Users/saljh8/Dropbox/Transfer/COVID_with_umap_and_markers.h5ad`. The file
holds 302,922 cells and 35,545 genes. Its `X` layer carries CP10k+log1p values
on a natural-log scale; 20 of 20 sampled cells sum to 10000.0 after `expm1`.
Its `layers['counts']` holds raw counts. Its `obs` carries `HLCA` (58 cell
states), `Library` (52), `group` (52) and UMAP coordinates.

The ADT lives in
`/Users/saljh8/Dropbox/Transfer/COVID-TotalVI/denoised_ADT_HTC_covid_all_samples.txt`.
The file holds 400,870 cells and 56 TotalVI-denoised antibodies on a linear
scale. The panel carries no isotype control.

`data.py` joins the two files by this rule:

```
Sample_<Library>_<BC>-1-<k>   ->   <BC>-1.<Library>
```

The join keeps **302,876 cells**. That number is 302,876 of 302,922 RNA cells
(99.98%) and 302,876 of 400,870 ADT rows (75.55%). The 97,994 unmatched ADT
rows are cells the RNA file dropped at quality control. 0 of 400,870 ADT row
names failed to parse.

## Target scale

`data.py` applies `np.log1p` to the denoised ADT values. The bone marrow atlas
stores its ADT on the same log scale, and
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/flask/reference_config.json`
declares `expression_scale: log1p, log_base: e` for the bone marrow bundle. The
lung values after `log1p` run 0.000 to 9.496 with mean 1.872; the bone marrow
values run 0.000 to 10.201 with mean 2.343.

`compare_target_scale.py` measured the choice against linear targets on the same
cells and the same features. log1p wins on both splits:

| Target | Split | Train cells | Mean Pearson | Median Pearson |
|---|---|---|---|---|
| log1p | cells | 40,000 | **0.442** | **0.395** |
| linear | cells | 40,000 | 0.310 | 0.298 |
| log1p | donor | 40,000 | **0.331** | **0.268** |
| linear | donor | 40,000 | 0.279 | 0.235 |
| log1p | cells | 150,000 | 0.447 | 0.399 |
| log1p | donor | 150,000 | 0.340 | 0.273 |

Raising the training set from 40,000 to 150,000 cells adds 0.005 mean Pearson on
the cell split and 0.009 on the donor split. The shipped bundle therefore keeps
the 40,000-cell protocol the mouse README documents.

## Curated map

`adt_rna_map.py` maps all 56 panel antibodies to HGNC symbols. An audit against
the 35,545 atlas genes matched **56 of 56 ADTs** with no missing partner gene.
32 of the 56 entries copy the partner set already curated for the bone marrow
panel in
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/adt_rna_map.py`.

`clean_adt_name` strips the BioLegend catalogue suffix and the punctuation, so
`CD274 (PD-L1)_A0007` becomes `CD274_PD-L1`. `atlas_var_name` then adds the
`Hu.` prefix the bone marrow bundle uses.

## Empirical whitelist

`build_whitelist.py` ranks all 35,545 genes against each ADT by Spearman
correlation over a 30,000-cell subsample, then applies the bone marrow rule:

1. A curated partner inside the top 100 gives that partner alone.
2. No curated partner inside the top 100 gives the partner plus 2 correlates.
3. A fallback gene chosen by more than 10 ADTs is dropped and replaced.

The build produced 56 rows and **81 unique RNA genes**:

| feature_source | n ADTs |
|---|---|
| `curated_in_top100` | 13 |
| `curated+fallback` | 43 |
| `fallback_only` | 0 |

The dedup step dropped `CD163` and `FOS`.

### Deviation: protein-coding fallbacks

The bone marrow and mouse builders draw fallback genes from all genes. This
build draws them only from the 20,014 protein-coding symbols in
`/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl-annotations_simple-PC.txt`.
18,810 of the 35,545 atlas genes (52.9%) are protein-coding. The top-100 ranking
still runs over all genes, so the `curated_in_top100` decision is unchanged.

A control build with `--no-protein-coding-filter` writes
`/Users/saljh8/Dropbox/Transfer/COVID-TotalVI/rna2adt_lung/artifacts/empirical_whitelist_no_pc_filter.tsv`.
The filter changed **1 of 56 ADTs**: `Hu.CD38` took `JARID2` instead of the
lncRNA `NEAT1`. All other 55 rows are identical.

## Trained model

`HumanLungRna2AdtModel` fits one `sklearn.linear_model.ElasticNet` per ADT with
`alpha=0.01` and `l1_ratio=0.5`. Every head reads the full 81-gene panel union;
L1 inside each head selects that head's genes. A `StandardScaler` z-scores the
input, and a second `StandardScaler` maps the output back to the log1p ADT
scale. The settings match the bone marrow bundle and the mouse bundle.

* input: 81 genes
* output: 56 ADTs, named `Hu.<marker>`
* training cells: 40,000
* test cells: 15,000

## Performance

| Split | Mean Pearson | Median Pearson | Mean Spearman | Valid ADTs |
|---|---|---|---|---|
| Random cell holdout | 0.442 | 0.395 | 0.428 | 56 of 56 |
| 10 held-out donors | 0.331 | 0.268 | 0.323 | 56 of 56 |

The donor split holds out D105, D116, D239, D292, D303, D307, D334, D346, D370
and D386. No donor appears in both halves.

Ten ADTs reach Pearson above 0.55 on the cell holdout: CD206, CD45, CD31, CD326,
CD90, CD49a, CD16, CD49f, HLA-DR and CD14. Per-ADT scores sit in
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/lung/rna2adt_hs_lung_per_adt_metrics_cells.tsv`
and
`/Users/saljh8/Dropbox/Transfer/COVID-TotalVI/rna2adt_lung/artifacts/rna2adt_hs_lung_per_adt_metrics_donor.tsv`.

The mouse bundle reaches mean Pearson 0.637 on its cell holdout. The lung bundle
scores lower, and the next section gives the reason.

## Limits — read before you use a marker

**The panel, not the model, sets the ceiling.**
`/Users/saljh8/Dropbox/Transfer/COVID-TotalVI/rna2adt_lung/artifacts/diagnose_ceiling.py`
correlated each measured ADT against its own cognate transcript over 30,000
cells. **32 of 56 ADTs score below 0.1 in absolute Spearman.** Only 10 of 56
exceed 0.3. Those 32 antibodies also show a p99/p50 ratio between 1.3 and 1.9,
so the denoised value barely moves across cells. B-cell and T-cell markers
dominate that group, and lung tissue holds few of those cells. Read
`/Users/saljh8/Dropbox/Transfer/COVID-TotalVI/rna2adt_lung/artifacts/adt_ceiling_diagnosis.tsv`
before you trust any single marker.

Predictions for those 32 ADTs carry the panel's background, not a protein
measurement. The bundle still emits all 56 ADTs, because dropping a tier of
low-confidence features would hide the problem rather than report it.

**A few fallback genes carry many heads.** `THBS1` feeds 16 of 56 ADTs, `IL1R2`
feeds 11, and `ACSL1` and `FMN1` feed 10 each. The dedup rule drops a gene that
exceeds 10 ADTs in its first pass, then lets replacement picks push another gene
past the same limit. The mouse and bone marrow builders behave the same way, and
`validate_port.py` proves the port reproduces that behaviour.

**Cell-state coverage is untested.** These scores pool all cell states. No
result here reports accuracy inside a rare state.

## Deviations from the bone marrow and mouse code

1. Fallback genes come only from the EnsMart100 protein-coding list. The section
   above measures the effect at 1 of 56 ADTs.
2. `train_lung.py` draws the cell holdout from `rng.permutation`.
   `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/mouse/train_mouse.py`
   sorts a random sample and then slices the head and the tail, which makes its
   test set the highest-indexed cells rather than a random set. On this atlas the
   sorted-sample split reported mean Pearson 0.257 where the permutation split
   reports 0.442.
3. `model.py` holds the model class. A class defined inside a `python -m` entry
   point pickles as `__main__.<name>` and then fails to unpickle inside
   cellHarmony-web.

## Validation performed

* `validate_port.py` drove the ported ranking and selection code with the mouse
  atlas, the mouse curated map and seed 0, then compared every cell against
  `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/mouse/configs/empirical_whitelist.tsv`.
  Result: **all 721 cells identical across 103 ADT rows and 7 columns**,
  including the `Gm42418` dedup the mouse README records.
* `compare_target_scale.py` reproduced the shipped bundle's cell-holdout score
  (0.442 mean, 0.395 median) and the donor-holdout score (0.331, 0.268) from an
  independent code path.
* `validate_bundle.py` loaded the bundle through
  `altanalyze3.components.rna2adt.api.load_bundle`, predicted 5,000 cells
  straight from the query `.h5ad`, and ran the cellHarmony-web ADT builder on
  the `hs_lung_hlca_reference` entry. Result: **PASS**. The API matched 81 of
  81 model genes with 0 missing, and returned 5,000 rows by 56 ADTs. The
  builder produced an AnnData of 5,000 by 56 with `expression_scale=log1p`,
  base e, a linear `counts` layer running 0.000 to 2708.6, and 62 clipped
  negative predictions. Scores on those cells were 0.441 mean Pearson and
  0.390 median, which agree with the clean cell holdout; those 5,000 cells are
  not a clean holdout, because roughly 13% of them sat in the training set.

## cellHarmony-web wiring

`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/flask/reference_config.json`
now offers the ADT modality on all 4 human lung references:
`hs_lung_cellref_reference`, `hs_lung_hlca_reference`, `hs_lung_natri_reference`
and `hs_lung_bpd_sun_reference`. Each entry carries:

```json
"adt": {
  "bundle_path": "../../rna2adt/lung/rna2adt_hs_lung_bundle.pkl",
  "expression_scale": "log1p",
  "log_base": "e"
}
```

The registry resolves `bundle_path` against its own directory. The bone marrow
entry and the mouse entries are unchanged.

## Regenerating

```bash
cd /Users/saljh8/Documents/GitHub/altanalyze3
export PYTHONPATH=.
PY=/opt/homebrew/opt/python@3.11/bin/python3.11

# 1. Write and audit the curated map
$PY -m altanalyze3.components.rna2adt.lung.adt_rna_map --write-tsv \
  --audit-h5ad /Users/saljh8/Dropbox/Transfer/COVID_with_umap_and_markers.h5ad

# 2. Prove the port still matches the mouse whitelist
$PY -m altanalyze3.components.rna2adt.lung.validate_port

# 3. Rebuild the whitelist
$PY -m altanalyze3.components.rna2adt.lung.build_whitelist

# 4. Train
$PY -m altanalyze3.components.rna2adt.lung.train_lung --holdout cells

# 5. Score held-out donors
$PY -m altanalyze3.components.rna2adt.lung.train_lung --holdout donor --no-save-bundle

# 6. Check the bundle and the web wiring
$PY -m altanalyze3.components.rna2adt.lung.validate_bundle --n-cells 5000
```

Steps 3 to 5 give the same result for the same `--seed` (default 0).
