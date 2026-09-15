# scALABLE and scALABLE-viewer

scALABLE is the web application of AltAnalyze3 for single-cell analysis. scALABLE-viewer serves a
finished atlas through the same interface without a Run tab. This page names the two programs,
their components, and the document that owns each fact, as of 2026-09-14.

| Program | Location | Purpose |
| --- | --- | --- |
| scALABLE | `altanalyze3/components/cellHarmony/webapp/` | upload, QC, cellHarmony alignment, MarkerFinder, approximate UMAP, fastComm, multimodal imputation, group differentials, chat |
| scALABLE-viewer | `altanalyze3/components/visualization/scalable_viewer/` | precompute a bundle from an h5ad and serve it; catalog of datasets; Study tab; chat over precomputed statistics |

## Documents of record

| Fact | Owner |
| --- | --- |
| how to run the server, environment variables, Docker, job retention | `altanalyze3/components/cellHarmony/webapp/README.md` |
| the reference registry and how to add a reference | the same file, section "The reference registry" |
| every pipeline step and its parameters | the same file, section "The pipeline" |
| modality ids, imputation models, output files, value scales | the same file, section "Modalities" |
| differential tests, thresholds, GO-Elite, networks, output files | the same file, section "The differential engine" |
| chat architecture, protocols, statuses | the same file, section "Chat" |
| the interface, control by control | `altanalyze3/components/cellHarmony/webapp/HOW_TO_USE.md` |
| bundle format, `precompute.py`, `prepare_assets.py`, `validate.py`, `run.py` flags | `altanalyze3/components/visualization/scalable_viewer/README.md` |
| what the viewer adds or changes for a reader | `altanalyze3/components/visualization/scalable_viewer/HOW_TO_USE.md` |
| measured chat routing accuracy | `altanalyze3/components/visualization/scalable_viewer/VALIDATION.md` |
| the alignment CLI `cellHarmony_lite.py` | `docs/cellHarmony.md` |
| the differential CLI `cellHarmony_differential.py` | `docs/cellHarmony_differential.md` |
| the approximate UMAP method | `docs/approximate_umap.md` |

## Components scALABLE calls

| Component | Path | Role |
| --- | --- | --- |
| cellHarmony_lite | `altanalyze3/components/cellHarmony/cellHarmony_lite.py` | cosine alignment of cells to reference centroids |
| MarkerFinder heatmap | `altanalyze3/components/visualization/marker_heatmap_h5ad.py` | top 50 markers per state and the heatmap cache |
| NetPerspective | `altanalyze3/components/visualization/NetPerspective.py` | marker and differential interaction networks |
| approximate UMAP | `altanalyze3/components/visualization/approximate_umap.py` | placement of query cells on the reference map |
| fastComm | `altanalyze3/components/fastComm/` | receptor-ligand communication scores |
| cellHarmony_differential | `altanalyze3/components/cellHarmony/cellHarmony_differential.py` | pseudobulk moderated t-test, Wilcoxon on cells, GO-Elite, fold matrices |
| GO-Elite | `altanalyze3/components/goelite/` | hypergeometric enrichment with DAG prioritisation |
| rna2adt | `altanalyze3/components/rna2adt/` | ADT imputation per cell |
| rna2lipid | `altanalyze3/components/rna2lipid/` | lung lipid imputation per cell; AML lipid imputation per pseudobulk under `aml/` |
| rna2metabolite | `altanalyze3/components/rna2metabolite/` | AML metabolite imputation per pseudobulk |
| rna2grn | `altanalyze3/components/rna2grn/` | TF-to-target edge scores per pseudobulk and per-cell TF activity |
| fastCNV | `altanalyze3/components/fastCNV/` | optional clone analysis behind `CELLHARMONY_ENABLE_FASTCNV` |

## Entry points

```bash
# scALABLE
cd /path/to/altanalyze3
python3.11 -m uvicorn altanalyze3.components.cellHarmony.webapp.app:app --host 127.0.0.1 --port 8000

# scALABLE-viewer
PYTHONPATH=/path/to/altanalyze3 python3.11 -m altanalyze3.components.visualization.scalable_viewer.run \
  --root /path/to/bundles --assets /path/to/assets --port 8062
```

Both chat tabs call an intent router on port 8001 that runs inside the LungMAP site process;
neither program holds a language model. The router's address comes from
`CELLHARMONY_ASSISTANT_URL` and `SCALABLE_ASSISTANT_URL`.
