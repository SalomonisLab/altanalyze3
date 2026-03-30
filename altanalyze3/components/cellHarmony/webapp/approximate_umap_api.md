# approximate_umap API

This FastAPI endpoint exposes the full server-side `approximate_umap` tool from
`altanalyze3.components.visualization.approximate_umap` through the cellHarmony
web app.

Endpoint:

```text
POST /api/tools/approximate-umap
```

## Request body

JSON fields mirror the CLI arguments.

Required:

- `query_cluster_key`
- `outdir`
- either `query` or `query_clusters_tsv`
- either `reference` or both `reference_coords_tsv` and `reference_clusters_tsv`

Optional:

- `reference_cluster_key`
- `umap_key`
- `jitter`
- `num_reference_cells`
- `random_state`
- `custom_colors_tsv`
- `restrict_obs_field`
- `restrict_obs_value`
- `output_prefix`
- `output_h5ad`
- `output_pdf`
- `save_updated_h5ad`
- `verbose`

## Example using a query h5ad and TSV reference

```bash
curl -X POST http://127.0.0.1:8000/api/tools/approximate-umap \
  -H "Content-Type: application/json" \
  -d '{
    "query": "/path/to/user_query.h5ad",
    "reference_coords_tsv": "/path/to/Hs-MarrowAtlas-L3M_reference_umap.tsv",
    "reference_clusters_tsv": "/path/to/Hs-MarrowAtlas-L3M_reference_clusters.tsv",
    "query_cluster_key": "cell_state",
    "reference_cluster_key": "Population",
    "outdir": "/tmp/approximate_umap_run",
    "save_updated_h5ad": true,
    "verbose": true
  }'
```

## Example using query/reference TSV inputs only

```bash
curl -X POST http://127.0.0.1:8000/api/tools/approximate-umap \
  -H "Content-Type: application/json" \
  -d '{
    "query_clusters_tsv": "/path/to/query_assignments.tsv",
    "reference_coords_tsv": "/path/to/reference_umap.tsv",
    "reference_clusters_tsv": "/path/to/reference_clusters.tsv",
    "query_cluster_key": "Population",
    "reference_cluster_key": "Population",
    "jitter": 0.05,
    "num_reference_cells": 1,
    "random_state": 0,
    "outdir": "/tmp/approximate_umap_run"
  }'
```

## Response

Successful responses return:

- `status`
- `query_cluster_key`
- `reference_cluster_key`
- `umap_key`
- `cluster_order`
- `n_query_cells`
- `outputs`

`outputs` contains the generated file paths:

- `coordinates`
- `augmented`
- `placeholder_expression`
- `comparison_pdf`
- `comparison_pdf_plain`
- `updated_h5ad` if `save_updated_h5ad=true` or `output_h5ad` is supplied

## Notes

- Paths are server-side filesystem paths.
- `output_h5ad` and `output_pdf` can be absolute paths or paths relative to
  `outdir`.
- `output_prefix` controls the shared base name for the TSV outputs.
