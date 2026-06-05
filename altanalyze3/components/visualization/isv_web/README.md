# ISV Web — Interactive Isoform Structure Viewer

A self-contained web app over the AltAnalyze3 isoform structure viewer engine
(`isoform_structure_view.py`) and its fast precompute caches (`gene_indexes_v2/`).
Independent of the cellHarmony web app.

## What it does
Two tabs over one shared, IGV-style genomic track renderer (D3/SVG): a **genomic coordinate ruler** +
**reference gene-model track** sit above the isoform structures, with **mouse-wheel zoom / drag-pan**
on the track area (hover + right-click preserved).

**Tab 1 — Cell type × covariate** (the primary view): cell-type-specific expression **heatmap columns
grouped into covariate blocks** with a labeled block + separator per covariate
(e.g. `HSC-1 HSC-2 MPP-1 │ HSC-1 HSC-2 MPP-1`, left block = covariate 1, right = covariate 2). Each row
is an isoform: heat strip on the left, genomic structure on the right. Cyan→yellow heat ramp,
row-normalized per isoform by default. Clustering thresholds live under **Advanced**.

**Tab 2 — Molecule view (read-level)**: the true ISV pileup. One **stacked read-pileup panel per
covariate** (e.g. `young`, `AML-NPM1`), each row = **one individual molecule** drawn at its exact
genomic structure (`tokens_raw`), colored by shared cluster, row height ∝ read count, with the gene
model + genomic axis at the bottom — the interactive analogue of the ISV PDF. This is produced by the
**engine itself**: the app drives `plot_isoform_structures_by_conditions` and parses the per-condition
`*_isoform_ids.tsv` it writes, so the molecules/clusters match the ISV output exactly (validated: HOPX
`young`+`HSC-1+HSC-2+MPP-1` → 300 molecules).

Both tabs:
- **Mouseover** shows exon-region id + genomic coords (and protein length on isoform exons); the row
  tooltip adds **protein length + NMD + read counts** (molecule id + sample + cluster in the read view).
- **Right-click** an isoform/molecule → export its **protein** or **mRNA** sequence (FASTA), all
  proteins in its cluster (heatmap tab), or all visible proteins.
- Search **genes** (symbol or ENSG); the column order follows your cell-type selection order.
- **Deep links / shareable URLs**: `?gene=HOPX&tab=molecule&cells=HSC-1,HSC-2,MPP-1&groups=young,AML-NPM1`
  auto-renders on load.

### Performance (read-level view)
First render of a new (gene × cell-states) selection drives the engine and may take **~10–30 s** (it
builds the per-sample isoform counts caches under each sample's `gene_indexes_v2/`). The result is
persisted to `<run_dir>/_isv_web_cache/reads/` **and** memoized in-process, so repeat views (including
after a restart) are near-instant. For uniformly fast first renders, pre-build the atomic per-cell-type
caches once with `precompute_viewer_index.precompute(sample_dict, barcode_sample_dict)`.

## Run
```bash
pip install -r requirements.txt   # fastapi, uvicorn, jinja2 (component-local)

python -m altanalyze3.components.visualization.isv_web.run \
    --metadata   /path/to/metadata.txt \
    --gene_model /path/to/Hs_Ensembl_exon.txt \
    --gene_symbol /path/to/Hs_Ensembl-annotations.txt \
    [--run_dir <dir with per-sample *-isoform.h5ad + gff-output>] \
    [--port 8050]
# then open http://127.0.0.1:8050
```
`run_dir` defaults to the metadata's directory. Files are assumed local.

## Fast queries
The backend loads the **combined isoform pseudobulk h5ad**
(`isoform_combined_pseudo_cluster_counts.h5ad`, obs = `<cellType>.<library>-isoform` columns, var =
`<ENSG>:<isoform>`) **once at startup** into a CSC matrix. Every query is a pure in-memory slice:
take the gene's isoform columns, select the obs rows matching the chosen samples/groups/cell-types,
and sum per output column. No full-h5ad read at query time. Measured on the local 4-sample run:
startup ~3 s, first gene query ~1.5 s, **warm query ~0 ms** (per-query memoization).

Isoform token structures come from `transcript_associations.txt` read **directly** (not via
`transcripts.db`, which indexes only ENST references and drops novel `molecule.sample` isoforms).
Protein length/NMD come from `protein_summary.txt` (indexed at startup); protein sequences from
`protein_sequences.fasta` (byte-offset index; seek-read on demand for export).

## Data sources (per run dir)
| menu / field | source |
|---|---|
| samples / groups | metadata (`library`, `groups`) |
| cell types | `cellHarmony/*/_barcode_clusters.txt` (`obs['cluster']`) |
| genes | per-sample `var_names` sidecars + `Hs_Ensembl-annotations.txt` |
| isoform structures | `transcript_associations.txt` (read directly; known + novel) |
| reference track + exon coords | `Hs_Ensembl_exon.txt` (`gene_segments` for the IGV ref track + ruler; `exon_lookup` for per-isoform segments / mouseover) |
| expression / read counts (heatmap) | `isoform_combined_pseudo_cluster_counts.h5ad` (column slice + sum) |
| individual molecules (read-level) | engine `plot_isoform_structures_by_conditions` over the per-sample `<library>.h5ad` + per-sample `gff-output/transcript_associations.txt`; parsed from its `*_isoform_ids.tsv` (cached in `_isv_web_cache/reads/`) |
| protein length / NMD | `gff-output/protein_summary.txt` |
| protein / mRNA sequence (export) | `gff-output/protein_sequences.fasta`, `transcript_sequences.fasta` |

## API
- `POST /api/isoforms` — heatmap tab: clustered isoforms + per-(cell-type×covariate) expression + `gene_model`.
- `POST /api/reads` — read-level molecule tab: per-covariate panels of individual molecules
  (`{gene, cell_types, conditions, max_isoforms}`) parsed from the engine's `*_isoform_ids.tsv`.
- `POST /api/molecules` — (legacy) pseudobulk molecule rollup; superseded by `/api/reads`.
- `GET  /api/isoform/{id}/protein|mrna`, `POST /api/proteins` — FASTA export (single / batch).
- `GET  /api/catalog`, `/api/genes`, `/api/junctions` — menus + autocomplete.

## Files
- `data_api.py` — data layer (RunContext + `query_isoforms` heatmap / `query_reads` read-level / `gene_model_track`); wraps the engine.
- `server.py` — FastAPI endpoints (responses memoized per signature).
- `run.py` — CLI launcher.
- `templates/index.html`, `static/app.js`, `static/styles.css` — frontend (two tabs + shared IGV-style renderer).
- `_tests/` — engine-vs-API parity + speed + e2e smoke.
