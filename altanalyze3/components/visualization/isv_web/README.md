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
- **Right-click** an isoform/molecule → export its **protein** (aa), **ORF/CDS** (nt), or **mRNA** (nt)
  sequence (FASTA), all proteins in its cluster (heatmap tab), or all visible proteins.
- Search **genes** (symbol or ENSG); the column order follows your cell-type selection order.
- **Deep links / shareable URLs**: `?gene=HOPX&tab=molecule&cells=HSC-1,HSC-2,MPP-1&groups=young,AML-NPM1`
  auto-renders on load.

### Performance (read-level view)
The read-level view has a **fast retrieval index**: a per-sample SQLite (`_isv_web_cache/mol_index/<lib>.reads.db`)
mapping `gene → (molecule_id, cell_type, count)`, built once at indexing time. When present, a gene
fetches its molecules in **~ms** (no whole-matrix counts load), then clusters **dynamically at viz time**
— measured: HOPX read view served in **~20–60 ms** (`source: "index"` in the response), reproducing the
engine's molecule counts (young+HSC = 300). Build it with `precompute_reads_index.build_reads_index(...)`,
or pass `--build_reads_index` to `run.py` (idempotent; ~3 min / ~4 GB for the 4-sample test set).

If the index is **absent** for any selected sample, the view falls back to driving the engine
(`plot_isoform_structures_by_conditions`), which is correct but slower (first render ~10–30 s, then
cached to `_isv_web_cache/reads/` + memoized).

### Molecule-view grouping (final collapsed isoform vs live clustering)
The read-level view groups reads two ways (`ctx.molecule_grouping`):
- **`final_isoform`** (default) — each read is grouped/coloured/labelled by the **final collapsed isoform**
  it maps to (the isoform-collapse protocol's read→isoform assignment), and **all reads in a group share one
  RNA + one protein FASTA** (the final isoform's). The final isoforms are *ordered* via `group_structures`
  so similar ones sit together (clustered, not random). Requires the **mol→final index**:
  `_isv_web_cache/mol_index/<lib>.mol2final.db`, built by `precompute_final_isoform_index` (joins each
  sample's `tier1/<lib>.mol2struct.tsv` with `gff-output/FINAL_structure_to_exemplar.tsv`; ~98–99 % of reads
  resolve, the rest fall back to their own id). It is built automatically alongside the reads index.
- **`structure`** (legacy) — live structure-clustering (`group_structures`) over the reads, kept intact.

Switch at runtime without restart: `curl -X POST 'http://HOST:PORT/api/grouping?mode=structure'` (or
`?mode=final_isoform`); `GET /api/grouping` reports the current mode. Env `ISV_GROUPING` sets the startup default.
Thresholds/strategy stay dynamic in `structure` mode.

## Optional filter for reads extending into model introns

In Molecule view → Advanced, enable **Hide reads extending into model introns**.
The option is off by default. It examines the full eligible read pool for the selected
isoforms, samples and populations **before** drawing up to **200 valid reads per isoform**.
It therefore refills each pileup from passing reads that the previous sample did not include.
If fewer than 200 pass, all passing reads are shown. The read cap is adjustable.

A read is hidden when one of its continuous aligned blocks crosses an assigned isoform's
splice-site boundary into an intron removed in that isoform. The default tolerance is 0 nt.
The comparison uses original alignment coordinates and the registered model's exon coordinates,
before any display compression. Whole read rows are hidden; exons are never clipped to make a
read appear compatible. Free transcript-end extensions are allowed. This is specifically a
splice-edge-extension filter, not a test of every possible difference between a read and model.

The response records examined, hidden, eligible and drawn counts separately, with every hidden
read identity and conflicting coordinates in `hidden_bleeding_reads`. Source read assignments,
count matrices, model definitions and other analyses are unchanged. Unresolved required alignment
or model coordinates produce an error rather than silently excluding the read.

Headless export supports `--hide-bleeding-reads --bleed-tolerance 0
--max-reads-per-isoform 200`. Omit `--max_reads` when requesting this per-isoform cap:
that older option imposes an additional cap across the whole panel. API equivalents are
`hide_bleeding_reads`, `bleed_tolerance` and `max_reads_per_isoform` on `/api/reads`.

## Trans-splicing candidates and other multi-gene transcripts

Transcripts containing gene-qualified exon tokens from multiple genes automatically use **separate
genomic axes for each locus** in the heatmap, molecule and usage views, and in PDF export. When a trans-splicing candidate is displayed, each locus has its own vertically stacked reference row, naming its gene, chromosome and strand. Ordinary views retain one reference row; multiple or overlapping gene labels alone do not add rows. Equal panel widths and the gaps between panels do not encode
genomic distance. Hover labels retain the original genomic coordinates, RNA support and read/final
isoform identities. Counts, collapse assignments, clustering and sequence exports are unchanged.

Read pileups use solid thin connectors to avoid aligned-dash stripe artifacts. Isoform model views retain dashed inter-gene connectors.
For these stacked reference rows, each genomic backbone continues through every displayed window on
the same chromosome. Real annotated exons are clipped to each window; the backbone does not imply
that the gene annotation extends across that whole span. Other chromosomes remain separate.

An orange dashed connector denotes a **candidate connecting non-overlapping reference loci**.
Reference gene extents must be pairwise disjoint on the same chromosome (different chromosomes are
also distinct). Overlapping gene annotations use gray connectors and are **not** called distinct-locus
trans-splicing candidates. Unresolved reference loci are labeled unresolved rather than classified
as candidates. Neither category establishes a trans-splicing mechanism: read-through,
rearrangements and sequencing/alignment artifacts remain possible.
This display flag does not replace the catalogue's physical-evidence audit or recalculate its
trans-splicing statistics.

Multi-gene exon coordinates are reconstructed separately for each gene from the annotated token
sequence. The legacy read-coordinate store contains untagged coordinate pairs, so it cannot safely
identify their chromosomes/genes for a multi-locus read. The view explicitly labels this coordinate
source. Saved coding-region annotations now distinguish thick CDS portions from thin UTR portions
in multi-gene models: the original CDS is intersected with each gene's exons before independent
axis projection. This requires every participating gene to match the saved CDS chromosome.
A single-chromosome CDS record is not extrapolated onto another chromosome; unsupported CDS/UTR
coverage remains unknown. Single-gene CDS/UTR rendering is unchanged.
The API adds `gene_models`, `has_multiple_loci` and per-row `gene_path`, `locus_status`,
`locus_joins` and `trans_splicing_candidate`; exon segments retain `gene`, `chrom` and `strand`.
Existing single-gene queries retain their previous coordinate and rendering path. No index rebuild,
recollapse or requantification is needed.

Headless PDFs extend the existing `isv_web_remote.render()` pileup, reference exon track and
registered isoform structure track, including its two-isoform gold/blue palette and predicted protein
length labels. There is no separate multi-locus plot renderer. Every RNA retains its individual
indexed exon annotation; registered structures identify the final isoforms independently of those
reads. Multi-locus registered models require gene-qualified segments; the CLI obtains these from
final transcript associations and retains registered GFF exons for single-gene structures. Shared
registered lists are projected once. Single-gene rendering follows its existing path.

For headless multi-locus PDFs, omit the ambiguous single
`--region` option; `--intron_scale` applies independently within each locus.

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
| individual molecules (read-level) | per-gene retrieval index `_isv_web_cache/mol_index/<lib>.reads.db` (counts) + per-sample `transcripts.db` (structures); clustered live. Fallback: engine `plot_isoform_structures_by_conditions` → `*_isoform_ids.tsv` (cached in `_isv_web_cache/reads/`) |
| protein length / NMD | `gff-output/protein_summary.txt` |
| protein / ORF / mRNA sequence (export) | `gff-output/protein_sequences.fasta`, `orf_sequences.fasta` (CDS), `transcript_sequences.fasta` |

## API
- `POST /api/isoforms` — heatmap tab: clustered isoforms + per-(cell-type×covariate) expression + `gene_model`.
- `POST /api/reads` — read-level molecule tab: per-covariate panels of individual molecules
  (`{gene, cell_types, conditions, max_isoforms}`) parsed from the engine's `*_isoform_ids.tsv`.
- `POST /api/molecules` — (legacy) pseudobulk molecule rollup; superseded by `/api/reads`.
- `GET  /api/isoform/{id}/protein|orf|mrna`, `POST /api/proteins` — FASTA export (protein aa / ORF-CDS nt / full-transcript nt; batch proteins).
- `GET  /api/catalog`, `/api/genes`, `/api/junctions` — menus + autocomplete.

## Files
- `data_api.py` — data layer (RunContext + `query_isoforms` heatmap / `query_reads` read-level / `gene_model_track`); wraps the engine.
- `server.py` — FastAPI endpoints (responses memoized per signature).
- `precompute_reads_index.py` — builds the per-gene molecule retrieval index (`build_reads_index`); run at indexing time or via `run.py --build_reads_index`.
- `run.py` — CLI launcher (`--build_reads_index` builds/refreshes the retrieval index at startup).
- `templates/index.html`, `static/app.js`, `static/styles.css` — frontend (two tabs + shared IGV-style renderer).
- `_tests/` — engine-vs-API parity + speed + e2e smoke.

## Optional protein-domain architecture tracks

The headless ISV exporter accepts `--protein-features FIGURE_DIRECTORY`. The directory
must be a completed standard `altanalyze3 sclr-iso2func protein-impact-plots` output
created with `--protein-plot-style compact --coding-exon-map EXON_MAP`. ISV verifies
the saved manifest hashes and reads only the small domain-placement and coding-exon
tables. It does not realign sequences, download annotations, or scan FASTAs/GFFs.

```bash
python -m altanalyze3.components.visualization.isv_web.isv_web_remote \
  --run_dir RUN --gene NR4A2 \
  --cell_groups 'HSC-MPP=HSC-1,HSC-2,MPP-1,MPP-2' --conditions young,aged \
  --isoforms ENST00000339562 NATIVE_ALTERNATIVE_ID \
  --figure-mode --intron-scale 0.1 \
  --protein-features ISO2FUNCTION_FIGURES \
  --save-render-input NR4A2.render_input.json.gz --out NR4A2_with_protein.pdf
```

This adds native N-to-C protein tracks under the genomic read/reference/transcript
tracks. Protein tracks share their own amino-acid scale; they do **not** use the
genomic x-coordinate scale, even on the minus strand. Exon colors and AltAnalyze exon
IDs come from the saved coding-exon map. Domain bars use the existing iso2function
placements and conventional labels. `partial` denotes sequence-alignment coverage,
not measured activity. Native IDs appear beneath each protein in the same color as
its RNA pileup/transcript. The optional view uses regular panel/transcript labels,
an italic gene title, and no narrative footer. Existing caption arguments remain
available. Every input feature and its source remain in the original complete
iso2function annotation audit. Unplaced annotations are explicitly listed.

Missing, conflicting, or ambiguous native IDs, exon labels, reference projections,
and unexpected protein lengths stop export; ISV never drops those isoforms. An
explicit source-table/native-FASTA length discrepancy remains labeled on the protein
track. All selected isoforms require registered transcript models. No residues are
highlighted as substitutions. The default ISV layout is unchanged when the option
is omitted. A `.protein_tracks.json` sidecar records all native exon/feature
coordinates, source-manifest hashes, display colors, RNA counts, and PDF hash.

A previously saved response can be rendered without reopening the dataset:

```bash
python -m altanalyze3.components.visualization.isv_web.isv_web_remote \
  --run_dir RUN --gene NR4A2 --render-input NR4A2.render_input.json.gz \
  --figure-mode --intron-scale 0.1 --protein-features ISO2FUNCTION_FIGURES \
  --out NR4A2_with_protein.pdf
```

The saved response already defines the read/isoform selection. Query/filter arguments
are rejected with `--render-input`; rendering options remain available. Save each
selection's source provenance with the input. This path preserves its read counts,
original alignment coordinates, sampling and transcript models.

For the interactive viewer, add the same `--protein-features FIGURE_DIRECTORY` to
`python -m altanalyze3.components.visualization.isv_web.run ...`, then enable
**Molecule view → Advanced → Show protein domains and coding exons**. The `/api/reads`
request uses `protein_features: true`; its response contains the same vector
primitives as the PDF renderer. Protein tracks stay on their amino-acid scale during
genomic zoom/pan and are included in the existing vector SVG/PDF export.
