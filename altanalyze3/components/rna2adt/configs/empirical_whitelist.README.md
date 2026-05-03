# Empirical RNA → ADT feature whitelist

This file (`empirical_whitelist.tsv`) defines which RNA features each ADT
prediction head consumes when training the `rna2adt` model. It supersedes
the strict curated-only map in `adt_rna_map.tsv` for production training.

## Why an empirical whitelist

The strict curated map (one row per ADT, listing the canonical RNA partner
gene(s) — e.g. `Hu.CD19 → CD19`, `Hu.CD3 → CD3D, CD3E, CD3G, CD247`) is
biologically the right starting point but is mechanically too narrow on
single-cell data:

* Many surface proteins are stable on the cell surface long after their
  encoding RNA has dropped out, so single-cell `RNA → protein` correlation
  for the canonical gene is often weak (e.g. `Hu.CD8` vs `CD8A` at the
  global level).
* For multi-subunit complexes (CD3, CD8, integrins, FcγRs, HLA class I/II),
  any single subunit's RNA carries only part of the surface signal.
* Single-cell dropout zeros most of the canonical gene's expression on
  most cells, leaving a co-expressed partner of comparable rank-correlation
  but lower zero-fraction as the more informative feature.

The empirical whitelist preserves the curated partner whenever it carries
real signal in the training atlas, and otherwise *augments* it with the
most-correlated RNA genes empirically observed.

## Source data

* Training atlas: `adata_combined_rna_adt_annotated-titrated.h5ad` (Zhang
  2024 bone marrow CITE-seq, 4 donors, 71460 cells).
* ADT panel as released: 141 antibodies = 129 `Hu.*` (human) + 3 `HuMs.*`
  (cross-reactive Hu/Mouse) + 9 `Isotype_*` (controls).
* **The 9 `Isotype_*` controls and `Hu.IgG.Fc` are excluded** from the
  whitelist (and from the model's output set) because they have no
  biological RNA partner. Their nonspecific binding profile would force
  the per-head ElasticNet to fit on whichever genes happen to co-vary with
  panel-wide background, polluting the panel-union with non-informative
  features and hurting the predictive Pearson on real ADTs. After the
  drop, **131 ADTs** remain (129 `Hu.*` + 3 `HuMs.*` minus the
  non-existent `Hu.IgG.Fc`).
* RNA matrix: 36602 genes after removing the ADT entries from the var
  index.

## Generation algorithm

The whitelist is produced by
[`build_whitelist.py`](../build_whitelist.py). For each ADT:

1. **Compute global Spearman correlation** between the ADT and every RNA
   gene across a 30 000-cell random subsample of the training atlas. Keep
   the top-100 ranked genes per ADT.
2. **Look up curated partner gene(s)** from
   [`adt_rna_map.py`](../adt_rna_map.py) (the curated map is the source of
   truth for which RNA gene encodes the protein each antibody targets).
3. **Curated-in-top100 case:** if at least one curated partner appears
   among the ADT's top-100 correlates, the ADT's feature set is
   *exactly* the curated partner(s) that hit the top-100. No extras are
   added; the empirical evidence has confirmed the canonical relationship.
4. **Curated-not-in-top100 case:** the ADT's feature set is the curated
   partner(s) (always retained for biological priority, regardless of
   global rank) plus the top-2 ranked global correlates that aren't
   curated partners.
5. **No-curated case** (HuMs / Isotype ADTs that have no listed partner):
   the feature set is the top-3 ranked global correlates.
6. **Promiscuity dedup.** A *fallback* gene (one selected as a top-correlate
   substitute, **not** as a curated partner) that gets selected for more
   than 10 ADTs is dropped from those ADTs' lists and replaced with the
   next-highest top-100 correlate that hasn't already exceeded the limit.
   Curated partner genes are exempt from this rule (e.g. `PTPRC` is the
   curated partner of all four CD45 isoform ADTs and stays on each one's
   feature set even though that's >1 ADT pointing to the same gene).

The promiscuity rule prevents a small number of broadly co-expressed
genes (e.g. `CD230 = PRNP`, `LYZ`, `PRTN3`) from dominating dozens of
ADTs' feature sets and making the per-protein heads collapse to the same
signal. In the current build the deduplicator removed:

* `CD230` (broadly elevated in HSC/progenitor compartments)
* `LYZ` (myeloid lineage marker)
* `PRTN3` (granulocyte marker)

These were each picked as a fallback by >10 ADTs and replaced.

## Output schema

| column | type | meaning |
|---|---|---|
| `adt_clean` | str | ADT name with prefix stripped (`CD19`, not `Hu.CD19`) |
| `adt_raw` | str | original var-name from the atlas (`Hu.CD19`, `HuMs.CD44`, `Isotype_RTK4530`) |
| `has_curated_partner` | TRUE/FALSE | whether the curated map has any RNA gene listed for this ADT |
| `curated_in_top100` | TRUE/FALSE | whether at least one curated partner is in the global top-100 correlates |
| `feature_source` | enum | `curated_in_top100` / `curated+fallback` / `fallback_only` |
| `n_features` | int | number of genes in this ADT's feature set after dedup |
| `feature_genes` | comma-sep str | the actual RNA gene symbols used as features |

## Composition

From the current build (with isotype/IgG controls excluded):

| feature_source | n ADTs |
|---|---|
| `curated_in_top100` | 43 |
| `curated+fallback` | 85 |
| `fallback_only` | 3 |
| **total** | **131** |

Unique RNA genes across all feature sets: **214**. Compare this to the
strict curated-only set (`adt_rna_map.tsv`) which used 132–136 unique
RNA genes — the whitelist adds ~75 empirical genes that downstream
heads can depend on.

## How it's consumed by the model

Replacing the curated map with this whitelist requires a single change
to the rna2adt model classes (`rna2lipid_arch.py`,
`generalizable.py`): instead of building each head's input columns
from `load_curated_adt_rna_map()`, read this TSV and use the
`feature_genes` column. The downstream pipeline (cellHarmony imputation,
bundle save/load) is unaffected.

## Regenerating

```bash
python -m altanalyze3.components.rna2adt.build_whitelist \
  --adata /path/to/training_atlas.h5ad \
  --out altanalyze3/components/rna2adt/configs/empirical_whitelist.tsv \
  --top-n 100 --fallback-k 2 --promiscuity-limit 10 \
  --max-cells 30000 --rna-chunk 512
```

Re-run whenever the curated partner map (`adt_rna_map.py`) changes or
when adopting a new training atlas with a different ADT panel. The TSV is
deterministic given the same atlas and the same `--seed` (default 0).
