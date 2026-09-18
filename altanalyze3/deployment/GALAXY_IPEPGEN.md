# Opt-in SNAF outputs for OneClick iPepGen

The exporter implements input contracts extracted from the supplied
`Galaxy-Workflow-OneClick-iPepGen-Neoantigen-workflow.ga` (SHA-256
`bd42a760092e26ad387550fdea6281064e9d8543ca504077d26905480a8b40d3`).
The small [bundled contract](../components/neoantigen/data/ipepgen_contract.json)
records tool versions, nested step locations, length restrictions and file formats.
The original workflow file is unchanged. pyNeoQuant is not needed for this export.

## Enable during a SNAF run

```bash
altanalyze3 snaf --juncounts counts.tsv --hla sample_hla.tsv \
  --db_dir reference --genome_fasta hg38.fa --output snaf_run \
  --galaxy_integration
```

Normal runs do not create these additional files. The independent
`--export_proteomics` switch continues to control the pyNeoQuant exchange bundle.
Both switches can be used together. Combined candidate reports are selected
once; their per-sample copies are not imported a second time.

`--galaxy_workflow /path/to/workflow.ga` reads the actual FragPipe, PepQuery2 and
IEDB settings from a compatible replacement workflow, including edited lengths.
The bundled OneClick settings are used when this option is absent. Unsupported
input modes or ambiguous duplicate prediction tools fail explicitly.

## Export existing results without rerunning prediction

```bash
altanalyze3-neo galaxy-export \
  --candidates snaf_run/T_candidates/T_antigen_candidates_all.txt \
  --hla sample_hla.tsv --predictor MHCflurry --outdir galaxy_export
```

Use `--sample S1` to select one sample. Otherwise all supplied HLA samples have
separate output files, including samples with no reported candidates. Every
candidate's sample and allele must be present in the supplied HLA table.
Use a fresh export destination; existing nonempty destinations are rejected to
prevent stale sample files from being collected into a later Galaxy run.

## Files and workflow connections

The SNAF flag writes `snaf_run/galaxy_export/`:

| Output | Format and connection |
|---|---|
| `snaf.database.fasta` | Cohort candidate sequences, headers `>generic\|SNAF_<source hash>\|SNAF_splice_candidate`. Add as another FASTA at **11/4**, preserving the human/reference, nonreference and fusion inputs. |
| `by_sample/*.database.fasta` | Same format, restricted to one sample; preferred for independent sample searches. |
| `by_sample/*.pepquery.txt` | Unique unmodified peptides, **one per line without a header**. Connect directly to **12 input 0**, `Neoantigen_Peptide_Candidates_for_PepQuery`. Do not run the upstream FragPipe header-removal step on these files. |
| `by_sample/*.iedb.fasta` | Unique peptide FASTA; connect to **14 input 1**, `FASTA-for-IEDB`. |
| `by_sample/*.alleles.txt` | One normalized class-I allele per line without a header; connect to **14 input 0**, `IEDB-optityle-seq2hla-alleles`. |
| `samples.tsv` | Original sample IDs, safe collection identifiers, candidate counts, and paths for matching the four sample collections. |
| `candidate_map.tsv` | Search accession and full FASTA identifier → stable source/candidate IDs, peptide, sample, event, gene, isoform, HLA and original SNAF report fields, including junction coordinates when present. |
| `annotations.tsv` | The annotation workflow's eight named columns: `Peptide, Chromosome, Start, End, Strand, Annotation, IGV_Genome_Coordinate, UCSC_Genome_Browser`. Unavailable peptide coordinates remain blank and are explicitly marked. |
| `peptides.bed` | Validated supplied complete-peptide BED12 mappings; empty when mappings were not supplied. |
| `length_eligibility.tsv` | Per-candidate length and inclusion flags for each destination. No peptide is trimmed or padded. |
| `hla_alleles.tsv` | Sample/allele provenance table; use the separate headerless allele files for IEDB. |
| `integration.json` | Source workflow checksum/settings, input and output checksums, counts, coordinate semantics and wiring notes. |

Cohort `snaf.pepquery.txt` and `snaf.iedb.fasta` are also provided for inspection
or intentional pooled analyses. Use the matched **sample collections** with
sample-specific raw MS data and HLA calls. An empty sample file indicates no
eligible candidates; skip the corresponding downstream search rather than
submitting an empty input to an external service.

The supplied workflow accepts **7–25 aa in FragPipe, 9–11 aa in PepQuery, and
8–12 aa in IEDB**. For example, an 8-mer is retained for FragPipe/IEDB and excluded
from the default PepQuery list. These bounds come from tool state in the `.ga`
file and are recorded in every export. IEDB's configured multiple lengths may
also generate shorter subsequences of an input peptide; SNAF has not assigned
scores to those new subsequences.

## Annotation branch and coordinate semantics

The existing **step 13** branch parses StringTie, SAV and INDEL accessions using
positional columns and regular expressions. SNAF accessions require a separate
join to `candidate_map.tsv`; sending them through that branch does not establish
valid SNAF annotations. The provided eight-column table gives the downstream
schema while retaining this distinction in its `Annotation` field.

SNAF report `coord` fields describe junctions, not the translated peptide's genomic
footprint. The exporter does not convert a junction span into a peptide BED row.
To include independently mapped peptides, pass `--galaxy_peptide_bed` during a
SNAF run, or `--peptide-bed` to `galaxy-export`. BED12 names must match the stable
`SNAF_<source hash>` accessions in the export. Coordinates must be zero-based,
half-open, with nonoverlapping coding blocks totaling three bases per amino
acid, and `thickStart/chromStart` and `thickEnd/chromEnd` equal. All mappings for
an accession are preserved. Browser links convert the start to one-based
coordinates; annotation Start/End remain BED coordinates. Set
`--galaxy_assembly` / `--assembly` if the supplied mappings use a different build.

The FASTA contains **reported candidate peptides**, not reconstructed full-length
proteins or all possible translated peptides. Keep the reference proteome and
FragPipe decoy/contaminant generation in the search. The workflow's sequence-based
FASTA deduplication can collapse accessions for a shared peptide: retain the full
candidate map and all compatible sources when interpreting returned evidence.
An exported prediction is not a FragPipe match or a PepQuery-validated peptide.

## Galaxy and Nextflow

Rebuild the SNAF container from this checkout before installing the updated Galaxy
wrappers; an older image with the same local build tag will not contain the exporter.

- The SNAF prediction wrapper has a **Create iPepGen integration outputs** switch,
  default off. Additional datasets and four matching sample collections appear
  only when enabled.
- The standalone [SNAF iPepGen export tool](galaxy/snaf_ipepgen_export.xml) produces
  these outputs from an existing combined report and HLA table.
- [SNAF-iPepGen-export.ga](galaxy/workflows/SNAF-iPepGen-export.ga) is a small
  importable export workflow for an installation containing that local tool.
  Its structure is checked locally; server import/execution remains to be verified.
- Nextflow `--mode full` or `--mode bam` accepts `--galaxy_integration`, with optional
  `--galaxy_workflow`, `--galaxy_peptide_bed`, and `--galaxy_assembly`. The SNAF task
  stages optional inputs and publishes `galaxy_export/` with the cohort results.

## Validation

Tests cover extracted workflow settings, sample/allele isolation, length filtering,
shared peptides, empty results, default-off behavior, combined-report selection,
path-safe collection IDs, and complete split-block peptide mappings. Local XML
lint, rendered command/collection tests, and Nextflow stub checks complement the
Python tests. These checks do not execute FragPipe, PepQuery, IEDB or a Galaxy server.
