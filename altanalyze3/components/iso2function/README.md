# iso2function

Isoform-resolved molecular-interaction & functional annotation for AltAnalyze3.

Attaches **definitively measured, isoform-resolved** interaction/function data — protein–protein
interactions (PPI), protein–DNA interactions (PDI), transcriptional activation (M1H), and condensate /
localization behavior — to long-read isoforms, anchored on the canonical **Ensembl91 structure
string**. From there it builds Cytoscape-ready interaction networks and runs an isoform-focused
gene-set enrichment, so that an isoform switch or expression difference can be interpreted as a
concrete change in a molecular-interaction network.

See **[GOALS.md](GOALS.md)** for the full design, identifier model, phased plan, risks, and constraints.

> **Scope:** *measured* interactions from the TFIso atlas (paper2, TF) **and** Yang 2016 (paper1,
> incl. non-TF) — plus a **UniProt domain-contact layer** (`evidence=domain_inferred`, kept strictly
> separate from measured). Every crosswalked isoform carries BOTH its Ensembl91 structure string and
> portable genomic splice-site coordinates.

## Data sources & crosswalk

| Source | What | Isoform key | -> structure |
|---|---|---|---|
| paper2 (TFIso atlas) | PPI/PDI/M1H/switch categories (measured, TF) | clone_id `SYMBOL-N` | gene+index -> bridge structure |
| paper1 (Yang 2016) | Y2H PPI + domain-domain + linear motifs (measured, incl. non-TF) | `SYMBOL_N` + ORF seq | ORF translate -> GENCODE ENST -> Ensembl91 structure |
| UniProt | isoform domain/interface features + curated PPIs (domain_inferred) | accession/isoform_id + seq | seq -> GENCODE ENST -> Ensembl91 structure |

ENST -> Ensembl91 structure via `ENST_reference_structures.tsv` (244k); structure -> genomic splice
coords via `structure_coords.sqlite`. All reused from existing Daedalus/MDS-AML outputs (no download).

## Status

| Layer | Module | Status |
|---|---|---|
| 1. Ingest | `ingest/paper2.py` | ✅ runnable + QC'd (10/10 tables) + tested; M1H call + Data_S8 categories + Data_S10 GSEA added |
| 2. Crosswalk | `crosswalk/structure_key.py`, `crosswalk/clone_map.py` | ✅ Leg A resolved by gene+isoform-index join → **724/745 atlas clones keyed to a structure** (520 known ENST); 21 unresolved reported, never forced |
| 2b. Associate (product) | `associate.py` | ✅ `isoform_function.tsv` (724 structure-keyed isoforms + PPI/PDI/M1H) and `switch_function.tsv` (447 ref→alt switches, 434 fully mapped, authors' Data_S8 categories) |
| 3. Network | `network/build.py`, `network/cytoscape.py` | ✅ binary-only graph (Y2H + eY1H), nodes carry the **structure key** + M1H call; `switch_consequence()` (Data_S8) + Cytoscape.js/SIF. N2H excluded (no cutoff) |
| 4. Enrichment | `enrichment/isoform_gsea.py`, `enrichment/switch_enrichment.py` | ✅ per-switch differential interactions (gained/lost PPI/PDI) + authors' GSEA re-keyed to structures; hypergeometric + BH over assayed background. **Validated**: my binary gained/lost matches Data_S8 PPI categories by direction. GO-Elite gene sets pending |
| 5. CLI / integration | `cli.py` + `utilities/parser.py` | ✅ **registered** as `altanalyze3 sclr-iso2func <action>`; also runnable standalone via `python -m ...cli` |

### Calling policy (binary-only, paper-defined)

Interaction calls use the authors' own thresholds (verified verbatim from `mmc12.pdf`): **Y2H** Boolean
(growth score ≥2 vs auto-activator) for PPI, **eY1H** Boolean (manual 3-researcher consensus) for PDI,
and **M1H** signal **≥1 = activator / ≤−1 = repressor**. **N2H** (`log2 NLR`) and **luciferase
validation** (`Log2 FC`) have **no numeric cutoff in the metadata** → they are *not* thresholded (N2H is
excluded from called edges; scores retained for reference). Switch consequences come from the authors'
`Data_S8` categories. No cutoff is ever invented.

## Quick start

```bash
# from the repo root, with the altanalyze3 python env
python -m altanalyze3.components.iso2function.cli all          # ingest + crosswalk + associate + network
python -m altanalyze3.components.iso2function.cli ingest       # just parse the supp tables
python -m altanalyze3.components.iso2function.cli associate    # build structure-keyed function tables
python -m altanalyze3.components.iso2function.cli network --clone TCF4-1 --clone TCF4-6
altanalyze3 sclr-iso2func query --gene MAX --interaction-type all --artifacts-dir results/MAX
```

Outputs land in the (gitignored) `data/` (tidy tables, crosswalk) and `artifacts/` (graph + Cytoscape).

## Standard measured-interaction queries

`altanalyze3 sclr-iso2func query` and
`python -m altanalyze3.components.iso2function.cli query` use **the same argument registration and
handler**. This is a reusable component action, with no dependency on dated analysis scripts.

By default it examines **all within-gene clone pairs**, including pairs absent from the authors'
classified-switch table. `--gene`, repeatable `--clone`, or `--pairs` can restrict the requested
universe. Pair TSVs accept `gene,A_clone,B_clone` (the existing switch-plot schema), or
`gene,clone_A,clone_B`. Every pair decision and target call is retained in audit outputs.

```bash
altanalyze3 sclr-iso2func query \
  --pseudobulk RUN/isv/isoform_combined_pseudo_cluster_counts.h5ad \
  --protein-fasta RUN/isv/gff-output/protein_sequences.fasta.gz \
  --metadata RUN/isv/metadata.tsv \
  --state HSC-1 --state HSC-2 --state MPP-1 --state MPP-2 \
  --condition young --condition aged \
  --min-isoform-reads 30 \
  --gene-cpm DATASET/gene_combined_pseudo_cluster_cpm.h5ad \
  --gene-counts DATASET/gene_combined_pseudo_cluster_counts.h5ad \
  --min-partner-cpm 1 --min-differential-partners 1 \
  --interaction-type all --combined-network --artifacts-dir results/HSPC
```

Use `--out-dir` to select the ingested atlas table directory. `--clonelist`, `--gene-annotation`,
`--mane-summary` and additive `--symbol-map` override reference inputs. Symbol-map TSV columns are
`symbol,ensg,source`, tab separated; record a verifiable source for retired symbols. Default references
come from `config.py`. Metadata requires `library,groups,original_library`; optional `original_uid`
records donor identity. States and conditions are repeatable and never hard-coded.

Rules enforced by this query:

- Native support requires **literal full AA sequence equality within the verified gene locus**.
  FASTA wrapping is removed, but stops, case and residues are not normalized. MANE supplies a gene
  locus only; its selected transcript never substitutes for an alternative protein. Native IDs,
  including novel IDs with periods, remain intact. All exact aliases are retained independently;
  each side needs at least one individual native isoform meeting the inclusive read threshold.
- Differences require **explicit True versus False within the same assay modality**. Missing,
  blank, untested and contradictory calls do not become negatives. Contradictory calls fail the query.
  Unknown targets remain in `all_pair_target_calls.tsv`. N2H scores are not thresholded.
- Protein partners pass when the **arithmetic donor mean of native gene CPM is strictly greater**
  than the threshold in at least one selected state. Zero-expression donors remain in the mean.
  Source-qualified duplicate libraries remain separate in isoform read support but represent one
  native gene library in donor means. All required state/library combinations must be present.
  `--gene-counts` verifies native CPM across the entire gene panel without renormalizing it.
- **PDI DNA baits are retained as DNA baits**, including mutant baits. Gene CPM cannot be assigned to
  a DNA sequence. Consequently `--min-partner-cpm` filters PPI partners, while measured PDI differences
  remain eligible based on native TF support. `DNA_bait_audit.tsv` preserves available bait sequences.
  Binding to a bait does not establish target-gene regulation or its direction in HSPCs.
- With `--interaction-type all` (default), a pair is retained if its supported PPI differences and/or
  measured PDI differences meet `--min-differential-partners`. Use `ppi` to reproduce a PPI-only screen,
  or `pdi` for the PDI screen. M1H measurements annotate pairs without an invented difference cutoff.
- Unresolved symbol/locus mappings and missing required observations fail with an incomplete manifest
  and diagnostic audits. A mapped gene absent from the quantified panel is **expression unavailable**,
  never zero, and cannot satisfy a CPM threshold. Missing native assay genes remain in coverage audits.

Outputs include `retained_pairs.tsv`, `divergent_interactions.tsv`, `all_pair_decisions.tsv`,
`all_pair_target_calls.tsv`, complete native/metadata/mapping/expression audits, and
`query_manifest.json` with source SHA256 hashes, exact roster hashes, parameters and coverage.
Interpret only a manifest with `status=complete`. Use a fresh output directory for each query.

`--combined-network` writes `combined_divergent_interactions.pdf` (combined overview plus complete
per-gene detail pages), node/edge TSVs and Cytoscape JSON. It reuses **ISV's `_combined_layout` /
igraph Fruchterman–Reingold** layout and the existing AltAnalyze3 edge renderer. Both positive and
explicit negative calls are shown for divergent targets; shared/unknown/unchanged comparisons remain
in the audit tables. This is an assay-call visualization, not inferred interaction activity.
It requires `python-igraph`, `matplotlib` and `adjustText` (label separation). From the repository root,
install the optional network dependencies with `pip install '.[iso2function]'`. The overview labels TFs;
all target labels appear on the detail pages. No generative graphics are used.

Export completed query results, including native isoform aliases/read support and raw calls, through
the existing export action (no recomputation or second scientific filter):

```bash
altanalyze3 sclr-iso2func export --query-results results/HSPC \
  --interaction-type all --out HSPC_PPI_PDI_divergent_interactions.txt
altanalyze3 sclr-iso2func export --query-results results/HSPC \
  --interaction-type ppi --out HSPC_expressed_PPI_divergent_interactions.txt
altanalyze3 sclr-iso2func export --query-results results/HSPC \
  --interaction-type pdi --out HSPC_TF_supported_PDI_divergent_interactions.txt
```

Each text file has a provenance sidecar with source and output hashes. PDI expression support refers
to the TF isoforms; a DNA bait has no gene CPM. The `expression_status` column retains this distinction.

## Inputs (read-only; never modified)

- **paper2 atlas** — `~/Downloads/TFIso/paper2/*.zip` (Mol Cell 2025 TFIso atlas). Members are streamed
  out of the zips; nothing is extracted into the download folder.
- **bridge** — `~/Claude_Project/TFiso_AltAnalyze3_annotation/tfiso_to_final_assignment.protein.tsv`
  (GFF clone → structure → observed isoform + protein consequences).
- **reference** — Ensembl91 exon table, genome FASTA, MDS-AML-KINNEX-4 collapse outputs. Shared-drive
  paths auto-resolve between local `/Volumes/...` and cluster `/data/...` (see `config.py`).

All paths are centralized in `config.py`.

## Outputs

`data/` (regenerable, gitignored):
- `ppi_y2h.tsv`, `ppi_n2h.tsv`, `pdi_ey1h.tsv`, `pdi_dna_baits.tsv`, `pdi_validation.tsv`,
  `activation_m1h.tsv`, `paralog_divergence.tsv`, `condensate.tsv` — tidy long-form tables.
- `ingest_manifest.tsv` — per-table QC (expected vs observed rows).
- `crosswalk_structure.tsv` — one structure-keyed row per bridge clone (+ observed isoform, ENST, NMD).
- `clone_to_structure.tsv` — **resolved** atlas clone_id → structure key (+ ENST/observed/NMD).
- `clone_unresolved.tsv`, `crosswalk_coverage.tsv`.
- `isoform_function.tsv` — structure-keyed per-isoform PPI partners / PDI targets / M1H call (product).
- `switch_function.tsv` — structure-keyed ref→alt switches with the authors' Data_S8 categories (product).
- `switch_differential_interactions.tsv` — per switch: PPI/PDI partners gained/lost (binary, co-tested).
- `isoform_gsea_by_structure.tsv` — authors' alt-vs-ref isoform GSEA (MSigDB) re-keyed to structures.

`artifacts/` (gitignored): `iso2function_graph.{nodes,edges}.tsv`, `*.cyjs.json` (Cytoscape.js),
`*.sif` (Cytoscape Desktop).

## Design notes (honesty / critical-thinking)

- **Structure string is the canonical key.** All other IDs (clone, ENST, ENSP, UniProt, the per-cohort
  `final_isoform_id` hash) are crosswalked *to* it, and every hop is verified — never guessed.
- **Tested ≠ untested.** A Y2H/Y1H `False` is *assayed-and-not-detected* and is kept as a real
  negative edge (`detected=False, tested=True`). Blank/not-determined graph calls remain unknown
  (`detected=None, tested=False`), never negative. The isoform-switch differential calls a partner "gained/lost" only within the
  *co-tested* space of both isoforms.
- **Leg A is not yet verified.** Interaction data is keyed on paper2 `clone_id`; the link to the GFF
  clone (and thus to the structure string) is emitted as candidates pending protein-sequence identity.
  Isoform nodes carry NO asserted ENST/structure until that verification runs.

## Testing

```bash
python -m pytest altanalyze3/components/iso2function/tests/ -q
```
Pure-logic tests need no data; the ingest integration test self-skips if the supp zips are absent.

## CLI integration — DONE (registered in `altanalyze3/utilities/parser.py`)

Available as `altanalyze3 sclr-iso2func <action>`. The global parser delegates all component flags and
action choices to `cli.add_arguments`; do not maintain a second manually copied list.

```python
# top of file, with the other component imports:
from altanalyze3.components.iso2function.cli import run_iso2function, add_arguments

# inside ArgsParser.get_parser(), with the other sclr-* subparsers:
iso2func = subparsers.add_parser(
    "sclr-iso2func", parents=[parent_parser],
    help="Annotate isoform PPI/PDI/function and build interaction networks (iso2function)")
iso2func.set_defaults(func=run_iso2function)
add_arguments(iso2func, include_loglevel=False)
self.add_common_arguments(iso2func)   # provides --loglevel/--threads etc.
```
(Do not also add `--loglevel` here — `add_common_arguments` provides it.)

## Constraints

No git/gh commands (the user owns all git operations and pushes). No mass deletions. Never modify
input files. Analysis artifacts (`data/`, `external/`, `artifacts/`) are gitignored. Python over
numpy/scipy/scikit-learn only (no R). No guessing identifiers/coordinates.

## Putative protein impacts without a direct native-isoform assay

### Protein domain and motif figures

`protein-impact-plots` renders completed, hash-verified `protein-impact` results through one
gene-independent plotting implementation. It performs no online annotation lookup or new
functional prediction. First run the `protein-impact` analysis described below, then:

```bash
python -m altanalyze3 sclr-iso2func protein-impact-plots \
  --impact-results protein_impacts \
  --genes TOX,CITED2,HOXB5,NR4A2,LYL1 \
  --artifacts-dir new_empty_figure_directory
```

To label those blocks with **AltAnalyze-assigned exon numbers**, enrich the saved map
using the same gene-model index as the ISV run. This reads only the required genes from
the existing SQLite index and does not scan the GFF or CDS table again:

```bash
python -m altanalyze3 sclr-iso2func protein-exon-map \
  --impact-results protein_impacts --coding-exon-map coding_exon_map \
  --gene-model-index REFERENCE/gene_indexes_v2/gene_model.db \
  --artifacts-dir coding_exon_map_labeled

python -m altanalyze3 sclr-iso2func protein-impact-plots \
  --impact-results protein_impacts --coding-exon-map coding_exon_map_labeled \
  --genes TOX,CITED2,HOXB5,NR4A2,LYL1 --protein-plot-style compact \
  --artifacts-dir new_empty_labeled_figures
```

Inside-block labels preserve reference exon-group IDs: `E4.1` and `E4.2` are regions
of `E4`. They are not renumbered per isoform. `I3` identifies retained sequence assigned
to reference intron 3; a block containing both is labeled `E3/I3`. Full region IDs and
genomic coordinates remain in the tables. Color-group ordinals are never used as exon
numbers. Unresolved reference coverage stays in the complete map and prevents plotting
the affected model with assigned labels. Small blocks use rotated inside labels when
necessary; the renderer rejects labels that cannot fit legibly. Queried reference rows,
their hash, and the parent map hash are exported for reproduction.

The default compact view displays reviewed UniProt functional annotations plus explicitly
identified literature regions. The CITED2 CBP/p300 segment is literature-derived. EnsMart100
protein-feature annotations remain in the complete audit and can be selected explicitly
with `--annotation-source`; Ensembl116 is used for protein catalog matching, not the displayed
feature boundaries. Exon numbering follows the explicitly supplied ISV gene-model index.

To combine all selected comparisons into **one figure on a single PDF page**, add
`--protein-figure-layout stacked` to `protein-impact-plots`. The reusable
`render_compact_figure()` function stacks the complete standard panels without rescaling
them or rasterizing their contents. The combined `protein_feature_comparisons.pdf` has
one page; matching combined SVG/PNG, panel geometry audit, individual panels and all
annotation/coordinate tables are also exported. The default `pages` layout retains the
existing multipage PDF. Both layouts use the same cached inputs and feature selection.

Use `--protein-figure-width 16` when a small coding exon needs more physical space for
its inside label; widths remain proportional to the original peptide coordinates.
The compact key explicitly defines partial marks as sequence alignment only, without
asserting retention of an intact functional domain. Conflicting reported protein lengths
are shown alongside the verified FASTA length as `169 aa (table: 144)`, rather than
silently replacing either source value.

For a gene with annotations from several reference proteins, select the desired reference
explicitly with `--feature-reference CRHBP=ENSP00000274368`. This selects the figure view
for that gene only; projections and all other source records remain unchanged in the
completed impact audit and in the figure's display audit (`reference_not_in_view`).
Unresolved gene/reference requests fail explicitly. Genes without this option retain
their existing annotation selection.

The default `--protein-plot-style compact` uses two proportional protein backbones with
colored domain/motif marks and a short feature key. Each row uses its own native amino-acid
coordinates; both rows have the same scale. The complete mapped domain interval includes
small insertions and substitutions without residue-by-residue difference highlighting.
Annotations retain their identities/colors across the two proteins; partial truncations are
labeled. There is no narrative footer. Native IDs, lengths, annotation names and literature
tags remain visible; read support, sequence differences and full provenance remain in the
unchanged audit tables and manifest. This compact renderer was added to AltAnalyze3; it is
not a pre-existing ISV plotting implementation.

Use `--protein-plot-style detailed` for the earlier audit-oriented layout, including its
alignment tracks, residue-difference colors and full caption. This is an optional diagnostic
view. Both styles use exactly the same verified proteins, feature projections and selection.

#### Coding-exon-colored protein backbones and ISV index reuse

The exon-colored protein-bar convention follows the visual reference in
`Gene._get_exon_colors()`, `Gene.exon_diagram()` and `Gene.protein_diagram()` from the TF isoform
paper code (Lambourne, Mattioli, Santoso et al., 2025). That package is not a runtime dependency.
The implementation retains the existing AltAnalyze3 annotations and coordinate projections;
it does not adopt the reference's sequence-change palette overrides or its practice of hiding
most annotations on lower protein tracks.

Build one reusable coding-exon map for the **complete audited native roster**, then pass it
to any subsequent plot selection:

```bash
python -m altanalyze3 sclr-iso2func protein-exon-map \
  --impact-results protein_impacts \
  --isv-structure-index RUN/isv/_isv_web_cache/final_structures.db \
  --artifacts-dir coding_exon_map

python -m altanalyze3 sclr-iso2func protein-impact-plots \
  --impact-results protein_impacts --coding-exon-map coding_exon_map \
  --genes TOX,CITED2,HOXB5,NR4A2,LYL1 --protein-plot-style compact \
  --artifacts-dir new_empty_figure_directory
```

`protein-exon-map` reuses ISV's existing registered `StructureIndex` with exact full native
identities; missing or ambiguous rows fail explicitly. It performs one source attestation
against the GFF/CDS hashes pinned in the completed impact audit, verifying every required
indexed exon coordinate and the existing genomic audit. The current ISV structure index
does not store its original source checksum, and the CDS loader otherwise scans a large text
file, so this registration is performed **once**, not once per figure or gene. Reuse the saved
map for reruns. Plotting reads only the small map and its manifest: no GFF, CDS, full SQLite
index, protein FASTA catalog or expression-matrix scan is performed.

Outputs are `coding_exon_map.json`, `coding_exon_coordinates.tsv` and
`exon_map_manifest.json`. The map preserves all native identities, original metadata
discrepancies and source protein hashes. Shared overlapping genomic coding-exon groups use
the same color across isoforms. Peptide blocks use exact CDS nucleotide lengths divided by
three, preserving fractional boundaries for codons split across exons. Only the terminal
stop codon, when present, is removed from the drawn protein extent. Unexplained CDS/peptide
length differences fail instead of truncating or dropping models. Colors never encode
amino-acid substitutions or sequence-change categories. Native coordinate scales and the
existing domain/motif annotations are unchanged. The figure export includes the selected
coding-exon coordinates and pins the complete map manifest.

Omit `--genes` to select all genes with audited native pairs. An explicitly requested gene or
pair without a verified comparison fails; no requested pair is silently skipped. Optional
`--pairs` accepts a TSV with `gene_symbol,isoform_A,isoform_B`. There is no read-count exclusion
in the renderer: the original counts and NMD calls are shown. A is the longer native protein
as defined by the analysis, not necessarily the canonical protein.

The default view displays `Domain`, `DNA binding`, `Zinc finger`, `Motif`, `Coiled coil` and
`Region` annotations from pinned reviewed UniProt and literature supplements. Repeat
`--feature-type` or `--annotation-source` to select another declared view. Provider names are
`reviewed_uniprot_snapshot`, `literature_construct_snapshot`, `uniprot_features` and
`protein_features`, subject to availability in the pinned audit. Every source feature row for
the selected pairs remains in `annotation_display_audit.tsv` with its display decision.
Identical displayed annotations may share one row; all original feature IDs are retained.
One annotation repeated on overlapping intervals from several reference proteins (the
EnsMart100 Ensembl table lists a family domain once per reference protein) is drawn once, on its
longest source interval; the other rows stay in the audit as
`coalesced_overlapping_same_annotation`.

**Essential feature classes are always drawn**, whatever `--feature-type` and
`--annotation-source` select: transmembrane helices (labelled `TM1`, `TM2`, ...), signal peptides,
kinase and phosphatase domains, transcription-factor DNA-binding domains, nuclear-receptor
ligand-binding domains, RNA-binding, helicase, GTPase, E3-ligase, protease, death-fold and
chromatin reader/writer domains. They take the lane nearest the protein backbone, are drawn above
any annotation they overlap, and carry their class as the label and key entry. Overlapping
rows of one class from several providers become the bars of the provider that resolves the most
separate intervals (ties to curated UniProt coordinates), so two zinc fingers stay two bars under
a family row that spans both; the other rows are `coalesced_same_essential_class` in the audit.
Every bar is labelled with its source name (UniProt rows without the `DOMAIN-` style key prefix;
an explicit ellipsis past 40 characters); the running id `F<n>` identifies a row only in the key
and the audit. Annotations absent from a native protein are listed by name with counts
("Not in this protein: Kinase x3, TM x2"). The rules are the
data file `plotting/essential_feature_whitelist.tsv` (class, label, priority, colour, provider and
regexes on the source name, accession and description), matched against the pinned vocabulary:
EnsMart100 UniProt feature keys (`TRANSMEM-`, `SIGNAL-`, `DNA_BIND-`, `ZN_FING-C2H2`,
`DOMAIN-Protein kinas`, ...), Pfam and InterPro accessions, and UniProt REST feature types.
`--essential-whitelist FILE` substitutes a reviewed copy. The figure manifest records the
whitelist path and SHA-256 and the counts `essential_rows_A`, `essential_rows_A_shown` and
`essential_rows_A_hidden_placeable`; the last must be 0. Rule set 2026-10-23 after a CLDND1 panel
drew a Pfam family hit eleven times and none of its four UniProt transmembrane helices.

Named UniProt `MOTIF` annotations are also retained by essential selection, including when
the requested feature types are Pfam accessions. The label preserves the literal source
description; the renderer does not infer a signaling motif from sequence or repair a
truncated source name. Distinct named motifs remain separate essential classes.

In compact and integrated ISV/domain plots, a coding segment too narrow for its exon label
keeps its true amino-acid width. Its literal label and amino-acid offsets appear in the
`Small coding segments` key instead. The layout gate still rejects a key that cannot fit;
the renderer neither stretches the segment nor silently removes its label.

Literature annotations are explicitly labeled. CITED2's CBP/p300-binding region is such an
annotation, not a built-in EnsMart100 domain.

Detailed plots use native A's 1-based inclusive amino-acid axis. B's actual residue coordinates
are also printed and exported. Missing sequence is gray, identical B sequence is green,
substitutions/insertions are orange, and uncertain placement is hatched. Partial domains
retain their true remaining widths; N-terminal shortening does not shift downstream domains
onto incorrect A coordinates. Repeated motifs and conflicting projections remain uncertain.
Features without an A coordinate remain visible as unplaced rows. The AA alignment and
feature projection use the existing `protein_isoform_audit` methods without modification.

Outputs include a combined multipage PDF; separate PDF, editable SVG and PNG files for each
pair/page; `plotted_feature_coordinates.tsv`; `alignment_events.tsv`; `plot_data.json`;
`annotation_display_audit.tsv`; and `label_layout_audit.json`. More than 12 features are
paginated with the same geometry. Label clipping/overlap checks fail execution. The final
`figure_manifest.json` records input/output/code hashes, annotation policy, RNA-support units,
software, font and FreeType versions, and projection discrepancies. A failed run retains an
incomplete manifest. Fixed export metadata and SVG IDs permit byte-identical reruns with the
same frozen source, inputs and software environment; cross-version rendering parity is not
assumed. These are sequence-based predictions, not direct measurements of native activity.

### Reproducible protein feature analysis

`protein-impact` is a separate, deterministic component action. It annotates native protein
sequences; it does **not** create experimental PPI/PDI calls or establish protein production in vivo.
It uses the existing `annotation/domain_impact.py` genomic coverage engine, the AltAnalyze2
`FeatureAlignment.importEnsemblFTdata` / `ExonAnalyze.compareProteinFeatures` sequence baseline,
and the existing `protein_isoform_audit` AA projection. No generative model or runtime database
query is involved.

```bash
python -m altanalyze3 sclr-iso2func protein-impact \
  --native-models native_models.tsv --required-roster required_roster.json \
  --reference-bundle references/bundle.json \
  --gff RUN/isv/gff-output --coding-regions RUN/isv/gff-output/coding_regions.txt \
  --pseudobulk RUN/isv/isoform_combined_pseudo_cluster_counts.h5ad \
  --metadata RUN/isv/metadata.tsv \
  --state HSC-1 --state HSC-2 --state MPP-1 --state MPP-2 \
  --condition young --condition aged --min-isoform-reads 30 \
  --prior-feature-audit prior/reference_mapped_features.tsv \
  --artifacts-dir new_empty_output_directory
```

Minimal native columns are `gene_symbol,ensg,isoform_id,protein_length,protein_sequence,
cds_sequence,cdna_sequence,annotated,nmd_status`. All identifiers, protein/CDS/cDNA sequences,
reported lengths and NMD flags are retained. The independent roster JSON requires
`source_table_sha256` and an ordered `required_roster` list of `{gene_symbol,isoform_id}` objects.
Required identities cannot be dropped because they fail a lookup. Source records outside the
independently defined roster are explicitly listed in `source_roster_audit.tsv`. CDS translation
must exactly reproduce the native protein (one terminal stop is removed from the *translation*,
never internal stops or arbitrary sequence characters). If `--gff` is provided, all required
proteins must also match its literal `protein_sequences.fasta.gz`, GFF exons and coding records.
The read threshold is an annotation of **every pair**, not an exclusion.

A reference bundle has `schema: "altanalyze-native-protein-impact-v1"`, `release_decision` and
`resources`. Each resource requires `path` (relative to bundle or absolute), `sha256`, `source`
and `release`. Explicit release choices are persisted; changed reference bytes fail execution.
Mandatory resource roles correspond to these built-in AltAnalyze database files:

| Resource role | EnsMart100 human archive member |
|---|---|
| `protein_relationships` | `ensembl/Hs/Hs_Ensembl_Protein__100_38.tab` |
| `protein_features` | `ensembl/Hs/Hs_ProteinFeatures_build_100_38.tab` |
| `protein_coordinates` | `ensembl/Hs/Hs_ProteinCoordinates_build_100_38.tab` |
| `uniprot_features` | `uniprot/Hs/Hs_FeatureCoordinate.txt` |
| `protein_sequences` | `Hs/RNASeq/SEQUENCE-protein-dbase_exoncomp.txt` |
| `ensembl_peptides` | Explicitly selected, pinned Ensembl peptide FASTA for exact catalog checks |

Official archive: https://altanalyze.org/archiveDBs/AltDatabase/updated/EnsMart100/Hs.zip .
The feature-reference release and novelty-search release are separate, declared roles; selecting
another release requires an explicit user decision. Optional `reviewed_uniprot` (REST JSON),
`uniprot_isoforms` (FASTA), `all_uniprot` (REST JSON) and `refseq_proteins` (FASTA) are frozen,
hashed supplements, **not** silently refreshed databases. Optional `curated_regions` uses the
existing literature-construct JSON schema (`gene_symbol,reference_accession,type,description,
location,evidences`) and is labeled `literature_construct_snapshot`, separately from built-ins.
Missing a region in the database cannot be interpreted as evidence that the protein retains it.

Outputs:

- `native_model_audit.tsv`, `native_genomic_audit.tsv`: complete native identity, translation,
  read-support, source-length and NMD checks, including discrepancies.
- `feature_projections.tsv`: **every source feature instance**, its reference protein and AA
  bounds, descriptive aligned/identical residues, exact-sequence presence, historical substring
  presence, and native genomic coding coverage. Invalid/missing/uncertain annotations remain
  unresolved rows. Source feature descriptions are preserved verbatim.
- `pair_feature_impacts.tsv`: B relative to A (A is longer; ID breaks ties), with the AA, historical
  substring and genomic comparisons separately exposed. Full features absent from both native
  proteins are `absent_from_both`, not a pair-specific loss. Partial sequence retention remains
  `altered_or_partial`; repeated exact hits retain sequence presence but ambiguous placement.
- `pair_summary.tsv`: all within-gene pairs, read eligibility, NMD and terminal changes. Feature
  counts are source-row counts, **not independent domains or interactions**.
- `protein_catalog_audit.tsv`, `exact_reference_matches.tsv`: literal peptide matching, kept
  separate from native transcript annotation. No-match is scoped to the pinned references.
- `prior_feature_regression.tsv`: re-derivation of every earlier feature projection when supplied.
- `impact_manifest.json`: complete source/reference/code/output hashes, method, software versions,
  unresolved annotation counts and interpretation limits. Failures leave it `incomplete`.

The historical feature slice uses `[start:end]` and expands short fragments. Its coordinate
convention is preserved **only as a labeled baseline**; verified feature intervals use 1-based
inclusive coordinates. Empty fragments are rejected. Genomic overlap cannot prove an intact AA
feature in a changed reading frame. Exact short-site matches cannot prove position-specific
modification or activity. A novel tail is never assigned a new domain by this workflow.

### Assay support when native proteins differ by one or two residues

Small AA differences do not make an annotated isoform a novel isoform. The default `query`
policy remains exact-only. An explicitly requested substitution allowance is reproducible:

```bash
altanalyze3 sclr-iso2func query ... --max-aa-substitutions 2 --artifacts-dir new_results
```

This allows **same-gene, same-length proteins with at most two substitutions**, with no indels
or stop normalization. An available exact assay match takes precedence. The allowance is a
reporting policy, not evidence that the substitutions are benign. Every mismatch is recorded
in `native_isoform_audit.tsv:match_details` and the pair's native-alias JSON. Nonidentical
matches are `putative_transfer_with_AA_substitutions`. Two different native protein sequences,
each independently reaching the read threshold, are required; one native protein cannot
supply both sides of a newly supported comparison. All 803 (or explicitly requested) assay
pairs remain in the audit, regardless of eligibility. Called interactions are still the
original measured Y2H/eY1H values; their transfer to a differing native sequence is predicted.

For explicit native/protein pairs, including pairs that do not meet the interaction filter:

```bash
altanalyze3 sclr-iso2func assay-match \
  --pairs native_assay_pairs.tsv --protein-fasta RUN/isv/gff-output/protein_sequences.fasta.gz \
  --clonelist SuppTable_CloneList.txt --out-dir parsed_atlas_tables \
  --pseudobulk RUN/isv/isoform_combined_pseudo_cluster_counts.h5ad --metadata RUN/isv/metadata.tsv \
  --state HSC-1 --state HSC-2 --state MPP-1 --state MPP-2 --condition young --condition aged \
  --artifacts-dir new_assay_audit
```

Pair columns: `gene,ensg,native_A,native_B,clone_A,clone_B`. This action calls the existing
`query` implementation for every assay comparison and retains all proposed native pairs.
`amino_acid_differences.tsv` reports exact residues and 1-based positions; edit blocks use
`difflib.SequenceMatcher(autojunk=False)`, the earlier mismatch audit's baseline. Zero-length
insertions/deletions have blank 1-based coordinates and explicit half-open coordinates.
Optional partner-expression flags are the same as `query`; unfiltered assay differences and
expression-qualified differences have separate columns. Reporter means and author-defined
activator/repressor/neither categories are retained even when there is no PPI/PDI difference.

Generate a **single component-produced report** containing all supplied evidence:

```bash
altanalyze3 sclr-iso2func report \
  --impact-results protein_impacts --query-results interaction_query \
  --assay-results native_assay_audit --artifacts-dir new_combined_report
```

This produces `iso2function_evidence.html`, `gene_evidence_summary.tsv` and a source/hash
manifest. Domain predictions, exact-protein assay support, assay support with AA differences,
and reporter results appear under each gene. No generative model writes its contents.

### Embed protein architectures in ISV

Completed compact `protein-impact-plots` outputs that include a coding-exon map can
be passed directly to ISV with `--protein-features FIGURE_DIRECTORY`. The shared
`plotting/isv_architecture.py` component verifies the cached placement/exon hashes,
requires complete coverage of every selected native isoform, and supplies one vector
layout to the ISV PDF and web SVG renderers. It introduces no new functional
predictions or sequence matching. See the [ISV interface](../visualization/isv_web/README.md#optional-protein-domain-architecture-tracks)
for CLI examples, the optional web control, coordinate conventions and provenance
sidecars.
