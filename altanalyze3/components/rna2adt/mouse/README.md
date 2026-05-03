# Mouse rna2adt — per-protein ElasticNet on empirical whitelist union

Murine adaptation of the [`rna2adt`](../README.md) bone marrow CITE-seq
ADT imputation model. Same architecture, same training protocol; only the
panel definition, the curated ADT-to-MGI gene map, and the trained bundle
differ.

## Files

| Path | What it is |
|---|---|
| [`adt_mgi_map.py`](adt_mgi_map.py) | Curated ADT (`ADT-*`) → MGI gene map (103 ADTs after dropping isotype controls). |
| [`configs/adt_mgi_map.tsv`](configs/adt_mgi_map.tsv) | TSV of the curated map written by `adt_mgi_map.py --write-tsv`. |
| [`build_whitelist.py`](build_whitelist.py) | Empirical whitelist generator (top-100 Spearman correlates per ADT, curated-priority + 2 fallback genes, promiscuity dedup). |
| [`configs/empirical_whitelist.tsv`](configs/empirical_whitelist.tsv) | The generated whitelist used for production training. 103 ADTs, 172 unique RNA features. |
| [`train_mouse.py`](train_mouse.py) | Trainer script: per-protein ElasticNet on the panel-union of whitelist genes. |
| [`rna2adt_mm_bundle.pkl`](rna2adt_mm_bundle.pkl) | Production-ready trained bundle. Loadable via `altanalyze3.components.rna2adt.api.load_bundle`. |
| [`rna2adt_mm_per_adt_metrics.tsv`](rna2adt_mm_per_adt_metrics.tsv) | Per-ADT cell-holdout Pearson and Spearman scores. |

## Source data

* Training atlas: `adata_combined_60k_rna_adt.h5ad` (Salomonis lab murine bone marrow CITE-seq, 57,636 cells, 28,804 vars = 28,692 mouse RNA + 112 ADTs).
* ADT panel: 112 antibodies prefixed `ADT-`, of which **9 are isotype controls** (`ADT-Mouse_IgG1/2a/2b___isotype_Ctrl`, `ADT-Rat_IgG1/1/2a/2b/2c___Isotype_Ctrl`, `ADT-Armenian_Hamster_IgG_Isotype_Ctrl`). These are excluded from both the whitelist and the model output set, leaving **103 ADTs**.
* Cluster annotations available in `obs`: `Level 1`–`Level 4`, `cluster`, `Capture`. Cluster column not used for prediction (the model is RNA-only).

## Curated ADT → MGI gene map

103 ADTs are mapped to MGI official mouse gene symbols. Coverage is complete after two corrections:

* `ADT-PIR_A_B`: mapped to `(Pira2, Pirb)` because `Pira` is no longer a current symbol but `Pira2` is the dominant family member present in the atlas.
* `ADT-TER_119_Erythroid_Cells`: mapped to `(Gypa,)` because `Ly76` (the legacy gene symbol for the Ter-119 antigen) has been retired from current MGI annotation. The Ter-119 epitope sits on the glycophorin-A complex.
* `ADT-CD85k_gp49_Receptor`: mapped to `(Lilr4b, Lilrb4a)` to capture both family members; the user's authoritative symbol `Lilrb4` is not in the atlas under that name.

Multi-subunit complexes list every chain (e.g. `ADT-CD3 → Cd3d, Cd3e, Cd3g, Cd247`; `ADT-CD8a → Cd8a, Cd8b1`; `ADT-I_A_I_E → H2-Aa, H2-Ab1, H2-Eb1`).

The audit script reports all 103 panel ADTs as fully matched to atlas RNA genes:

```bash
python -m altanalyze3.components.rna2adt.mouse.adt_mgi_map \
  --write-tsv \
  --audit-h5ad path/to/adata_combined_60k_rna_adt.h5ad
```

## Empirical whitelist

For each of the 103 panel ADTs we computed Spearman correlation against all 28,692 mouse RNA genes across a 30,000-cell random subsample of the atlas, then kept the top-100 ranked genes per ADT. The whitelist for each ADT is then:

1. **Curated partner(s) found in top-100** (35 ADTs): use only those genes. The empirical evidence has confirmed the canonical RNA→protein relationship.
2. **Curated partner(s) not in top-100** (68 ADTs): use the curated partner(s) plus the top-2 ranked global correlates that aren't curated partners.
3. **No curated partner** (0 ADTs after dropping isotypes): would use top-3 global correlates.

A *fallback* gene that gets selected as a substitute for >10 ADTs is dropped and replaced with the next-ranked correlate. In the current build the deduplicator removed:

* `Gm42418` — a broadly expressed mouse pseudogene/lncRNA that contaminates many cells' top correlates.

The resulting whitelist contains **172 unique mouse RNA genes** spanning the panel.

Composition:

| feature_source | n ADTs |
|---|---|
| `curated_in_top100` | 35 |
| `curated+fallback` | 68 |
| `fallback_only` | 0 |
| **total** | **103** |

## Trained model

Architecture: `MouseRna2AdtModel` in [`train_mouse.py`](train_mouse.py) — equivalent to the human `Rna2LipidArchPanelPerProtein` with `feature_source="whitelist_union"`.

* Input: per-cell RNA expression on the 172 panel-union genes, z-scored across the training set.
* Output: 103 ADTs (one head per ADT) z-scored back to the original TotalVI-denoised scale via the inverse `StandardScaler`.
* Per-protein head: `sklearn.linear_model.ElasticNet(alpha=0.01, l1_ratio=0.5, max_iter=2000)`. L1 sparsity within each head selects which of the 172 panel-union genes that head actually uses.
* Independent fits across the 103 ADTs (no joint multi-task penalty).

## Cell-holdout performance

40,000 training cells / 15,000 held-out cells:

| Metric | Mean | Median |
|---|---|---|
| Pearson | **0.637** | **0.669** |
| Spearman | 0.651 | 0.690 |

103/103 ADTs produced valid (non-degenerate) predictions.

Strongest predictions (Pearson > 0.85): `CD71 (Tfrc)`, `CD11b (Itgam)`, `CD49b (Itga2)`, `CD31 (Pecam1)`, `CD205_DEC_205 (Ly75)`, `CD55 (Cd55)`, `CD371_CLEC12A (Clec12a)`, `CD62L (Sell)`, `Ly_6G (Ly6g)`, `CD27 (Cd27)`, `Ly_6A_E_Sca_1 (Ly6a, Ly6e)`, `CD9 (Cd9)`. Weakest (Pearson < 0.35): `CD223_LAG_3`, `CD40`, `CD140a`, `CD301b`, `Siglec_H`, `CX3CR1`, `CD69`, `CD95_Fas`, `CD49a`. Per-ADT details in [`rna2adt_mm_per_adt_metrics.tsv`](rna2adt_mm_per_adt_metrics.tsv).

## Loading the bundle

```python
from altanalyze3.components.rna2adt.api import load_bundle
b = load_bundle('altanalyze3/components/rna2adt/mouse/rna2adt_mm_bundle.pkl')
print(b.metadata['approach'])    # rna2adt_mm_panel_per_protein_whitelist_union
print(len(b.input_genes))        # 172
print(len(b.output_adts))        # 103
preds = b.predict_from_adata(query_rna_adata)
```

The bundle's output ADT names are **clean marker names** (e.g. `CD19`, `Ly_6C`, `CD11b`) — the `ADT-` prefix and the `anti_<species>_<species>_` wrapper are stripped at training time so downstream consumers see consistent names. The training-atlas raw names (`ADT-CD19`, `ADT-anti_mouse_human_CD11b`) are recoverable by reading the curated map.

## Regenerating

```bash
# 1. Re-audit / refresh the curated map
python -m altanalyze3.components.rna2adt.mouse.adt_mgi_map --write-tsv \
  --audit-h5ad path/to/atlas.h5ad

# 2. Regenerate the empirical whitelist
python -m altanalyze3.components.rna2adt.mouse.build_whitelist \
  --adata path/to/atlas.h5ad \
  --out altanalyze3/components/rna2adt/mouse/configs/empirical_whitelist.tsv

# 3. Train the bundle
python -m altanalyze3.components.rna2adt.mouse.train_mouse \
  --adata path/to/atlas.h5ad
```

Steps 2 and 3 are deterministic given the same `--seed` (default 0).
