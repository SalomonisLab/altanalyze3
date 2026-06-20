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
```

Outputs land in the (gitignored) `data/` (tidy tables, crosswalk) and `artifacts/` (graph + Cytoscape).

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
  negative edge (`detected=False, tested=True`). Genuinely blank/not-determined cells are dropped, not
  rendered as negatives. The isoform-switch differential calls a partner "gained/lost" only within the
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

Available as `altanalyze3 sclr-iso2func {ingest|crosswalk|associate|network|enrich|all}`. For reference,
the registration that was applied (top-of-file import + subparser in `ArgsParser.get_parser()`):

```python
# top of file, with the other component imports:
from altanalyze3.components.iso2function.cli import run_iso2function

# inside ArgsParser.get_parser(), with the other sclr-* subparsers:
iso2func = subparsers.add_parser(
    "sclr-iso2func", parents=[parent_parser],
    help="Annotate isoform PPI/PDI/function and build interaction networks (iso2function)")
iso2func.set_defaults(func=run_iso2function)
iso2func.add_argument("action", choices=["ingest", "crosswalk", "network", "enrich", "all"])
iso2func.add_argument("--out-dir", default=None)
iso2func.add_argument("--artifacts-dir", default=None)
iso2func.add_argument("--clone", action="append", default=None)
iso2func.add_argument("--detected-only", action="store_true")
self.add_common_arguments(iso2func)   # provides --loglevel/--threads etc.
```
(Do not also add `--loglevel` here — `add_common_arguments` provides it.)

## Constraints

No git/gh commands (the user owns all git operations and pushes). No mass deletions. Never modify
input files. Analysis artifacts (`data/`, `external/`, `artifacts/`) are gitignored. Python over
numpy/scipy/scikit-learn only (no R). No guessing identifiers/coordinates.
