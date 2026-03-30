# cellHarmony-lite reference examples

These files are synthetic examples showing the minimum formats expected by the
web app reference registry and the downstream `cellHarmony_lite` and
`approximate_umap` loaders.

## Files

- `reference_config.example.json`
  Example registry entry pointing to the files in this directory.
- `example_states.tsv`
  Reference expression/state matrix for `cellHarmony_lite`.
- `example_reference_clusters.tsv`
  Reference barcode-to-cluster/population assignments for approximate UMAP.
- `example_reference_umap.tsv`
  Reference barcode-to-UMAP coordinates for approximate UMAP.

## Required formats

### 1. `states_tsv`

Loaded by `cellHarmony_lite.py` with:

```python
pd.read_csv(cellharmony_ref, sep="\t", index_col=0)
```

So the required format is:

- first column: gene identifier
- remaining columns: reference populations/states
- values: numeric expression/signature values

Example:

```tsv
Gene	HSC-1	MPP-1	Myeloid-1
MPO	0.1	1.4	3.2
GATA2	2.8	1.6	0.2
```

### 2. `reference_clusters_tsv`

Loaded by `approximate_umap._load_reference_from_tsv(...)`.

Required columns:

- one barcode column named like `barcode`, `cellbarcode`, `cell_barcode`, or `cell`
- one cluster label column matching the configured `cluster_key`
  or named like `cluster`, `population`, or `label`

Example:

```tsv
barcode	Population
ref_cell_1	HSC-1
ref_cell_2	MPP-1
```

### 3. `reference_coords_tsv`

Also loaded by `approximate_umap._load_reference_from_tsv(...)`.

Required columns:

- one barcode column named like `barcode`, `cellbarcode`, `cell_barcode`, or `cell`
- two UMAP coordinate columns named like `UMAP1`/`UMAP2`

Example:

```tsv
barcode	UMAP1	UMAP2
ref_cell_1	-2.1	1.4
ref_cell_2	0.3	0.2
```

## Using the example registry

Copy `reference_config.example.json` and replace the example file paths with
real ones on your system.
