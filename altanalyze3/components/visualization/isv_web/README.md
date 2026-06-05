# ISV Web — Interactive Isoform Structure Viewer

A self-contained web app over the AltAnalyze3 isoform structure viewer engine
(`isoform_structure_view.py`) and its fast precompute caches (`gene_indexes_v2/`).
Independent of the cellHarmony web app.

## What it does
- Select **samples**, **groups** (covariate), and **cell types**; combine columns **by group or by cell type**.
- Search **genes** (symbol or ENSG), **restrict by junctions**, and adjust **clustering thresholds** live.
- Interactive **D3/SVG isoform tracks** with an expression heat strip per column.
- **Mouseover** shows exon-region id + genomic coords, and (per isoform) **protein length + NMD**.
- **Right-click** an isoform → export its **protein sequence** (FASTA), or export all isoforms in its cluster.

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
| exon coords (mouseover) | `Hs_Ensembl_exon.txt` (`exon_lookup`) |
| expression | `isoform_combined_pseudo_cluster_counts.h5ad` (column slice + sum) |
| protein length / NMD | `gff-output/protein_summary.txt` |
| protein sequence (export) | `gff-output/protein_sequences.fasta` |

## Files
- `data_api.py` — data layer (RunContext + query); wraps the engine, no drawing.
- `server.py` — FastAPI endpoints.
- `run.py` — CLI launcher.
- `templates/index.html`, `static/app.js`, `static/styles.css` — frontend.
- `_tests/` — engine-vs-API parity + speed + e2e smoke.
