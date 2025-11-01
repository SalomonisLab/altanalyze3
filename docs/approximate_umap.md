# Approximate UMAP placement for query AnnData sets

This guide documents the Python helper introduced in
`altanalyze3.components.visualization.approximate_umap`. The helper ports the
legacy AltAnalyze2 `approximateUMAP.py` script to work directly with AnnData
inputs, allowing you to reuse UMAP coordinates from a reference dataset when
visualizing related query data.

At a high level the function:

- loads query and reference AnnData objects (paths or in-memory instances),
- validates that both contain the required cluster annotations,
- samples UMAP coordinates from the reference for every query cluster,
- inserts placeholder “Null-\<cluster\>” rows when the query is missing any reference clusters,
- returns a structured object that mirrors the text artifacts written by the
  original script, and
- optionally copies cluster colour palettes onto the query AnnData copy.

## Quick start (Python API)

```python
import anndata as ad
from altanalyze3.components.visualization.approximate_umap import approximate_umap

reference_adata = ad.read_h5ad("reference_umap.h5ad")

result = approximate_umap(
    query="query_cells.h5ad",
    reference=reference_adata,
    query_cluster_key="leiden",
    reference_cluster_key="cellHarmony_cluster",
    umap_key="X_umap",
    jitter=0.05,
    num_reference_cells=2,
    random_state=42,
)

print(result.query_adata)  # AnnData copy with `.obsm["X_umap"]` newly populated
result.write_text_outputs("outputs/my_query")
result.write_comparison_pdf(
    reference_adata=reference_adata,
    reference_cluster_key="cellHarmony_cluster",
    query_cluster_key="leiden",
    umap_key="X_umap",
    output_path="outputs/my_query-comparison.pdf",
)
```

## Quick start (command line)

The module exports a lightweight CLI. Supply the query and reference `.h5ad`
files, the relevant cluster keys, and an output directory:

```bash
python -m altanalyze3.components.visualization.approximate_umap \
  --query path/to/query.h5ad \
  --reference path/to/reference.h5ad \
  --query-cluster-key leiden \
  --reference-cluster-key cellHarmony_cluster \
  --outdir outputs/ \
  --output-prefix my_query
```

The CLI writes:
- `outputs/my_query-approximate-umap.tsv` – barcode-to-UMAP table.
- `outputs/my_query-augmented.tsv` – simple barcode listing.
- `outputs/my_query-comparison.pdf` – side-by-side reference/query UMAP plot.
- optionally, `outputs/my_query.h5ad` if `--save-updated-h5ad` is provided (or if
  `--output-h5ad custom_name.h5ad` is supplied, which auto-enables the export).

Use `--save-updated-h5ad` to write the augmented AnnData (compressed with gzip),
or pass `--output-h5ad custom_name.h5ad` to both enable the export and choose a
filename within `--outdir`. `--output-pdf` accepts either a filename (written
inside `--outdir`) or an absolute path.

### Inputs

- `query`: AnnData object or `.h5ad` path for the dataset that needs UMAP
  coordinates. `query.obs` must include a cluster column.
- `reference`: AnnData object or `.h5ad` path with a matching cluster column and
  a pre-computed UMAP embedding stored in `.obsm[umap_key]`.

### Key parameters

- `query_cluster_key`: Name of the column in `query.obs` describing cluster
  assignments.
- `reference_cluster_key`: Optional column in `reference.obs` defining the
  reference clusters. When omitted the query key is reused.
- `umap_key`: Matrix entry inside `reference.obsm` that holds the embedding.
  Only the first two dimensions are used when producing placeholder tables.
- `num_reference_cells`: Number of reference barcodes sampled per query barcode.
  Values >1 compute an average if the primary sampled cell was already used.
- `jitter`: Magnitude of uniform noise applied the first time a reference cell
  seeds an approximate coordinate. Set to `0.0` to disable.
- `random_state`: Seed for the NumPy random number generator, enabling reproducible
  layouts.

## Returned data

`approximate_umap` returns an `ApproximateUMAPResult` dataclass with the
following attributes:

- `query_adata`: Copy of the query AnnData with `.obsm[umap_key]` populated and
  cluster categories ordered to match the reference.
- `coordinates`: DataFrame indexed by barcode containing `umap_0` and `umap_1`.
  Includes any placeholder rows.
- `augmented_obs`: DataFrame describing the augmented clusters with an
  `is_placeholder` flag.
- `placeholder_expression`: Wide placeholder matrix identical to the original
  AltAnalyze2 export format.
- `reference_choices`: Mapping from each barcode to the IDs of the reference
  cells that were sampled when computing its coordinates.
- `cluster_order`: Preserved ordering of clusters derived from the reference,
  or from the query `uns["lineage_order"]` when present, and reused for plotting
  and categorical alignment.

Use `write_text_outputs` to persist the TSV files (`-augmented.tsv`,
`-approximate-umap.tsv`, `-placeholder-exp.tsv`). Call `write_comparison_pdf` to
produce a side-by-side PDF of the reference UMAP and the approximated query
UMAP.

## Error messages and troubleshooting

- **Missing cluster annotations**: If either AnnData object lacks the requested
  cluster column, a `KeyError` is raised.
- **Unknown clusters in the query**: The helper stops with a `ValueError` when a
  query cluster is absent from the reference. Either harmonise cluster labels or
  relax matching requirements.
- **Reference cluster without cells**: An empty cluster in the reference leads
  to a `ValueError`. Ensure your reference AnnData contains UMAP coordinates for
  every advertised cluster.
- **High-dimensional embeddings**: When `reference.obsm[umap_key]` has more than
  two columns only the first two are exported. A debug-level log message is
  emitted to highlight this.

## Suggested workflow integration

1. Run your reference pipeline once to compute a UMAP embedding and cluster
   assignments; save the resulting AnnData.
2. For each new query dataset, run `approximate_umap` with the reference to
   place the query cells. Consider pre-filtering or relabelling clusters to
   match the reference nomenclature.
3. Persist the `query_adata` copy for downstream plotting (e.g. Scanpy,
   Seaborn) or write legacy TSV files for compatibility with existing AltAnalyze
   utilities.
4. Keep track of the `reference_choices` mapping if you need to audit which
   reference cells influenced each query coordinate.

By consolidating the UMAP approximation logic inside AltAnalyze3’s visualization
component, you can migrate old workflows to the new AnnData-based ecosystem
while retaining reproducible layout behaviour and artifact formats.
