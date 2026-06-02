# Canonical Long-Read Isoform Collapse Design

## Goal

Create a scalable, evidence-aware isoform collapse strategy for single-cell long-read data that can process millions of reads per sample and up to roughly 200 samples without loading all BAM reads or all raw molecules together. The pipeline should identify known and truly novel isoforms, quantify final standardized isoform abundance per sample, and preserve rare high-confidence isoforms with unique exon structures. Isoform structures observed only once across all samples may be ignored.

The design preserves the existing AltAnalyze3 long-read annotation behavior in `components/long_read/gff_process.py`. In particular, the validated exon/intron annotation rules, Ensembl-prioritized structure annotation, splicing event annotation, and terminal-position-insensitive isoform annotation semantics should not be changed.

## Current Problem

The current single-cell long-read flow can extract raw reads from BAM files using `components/bam/isoform_structure_extract.py`, write a raw-read GFF, and write an h5ad matrix keyed by molecule IDs. This works well for splicing and junction analyses because those analyses quantify validated splice junctions directly.

It is not sufficient for differential isoform abundance because final isoforms need to be collapsed full-length or near-full-length transcript models with stable IDs across samples. Pigeon-derived isoform assignments were too variable between samples, and the current `components/long_read/isoform-collapse/collapse_transcript_associations.py` is only a preliminary postprocessor over transcript association text. It does not fully solve Ensembl prioritization, expression-aware representative selection, terminal handling, or h5ad remapping.

## Constraints

- Do not load all BAM reads from all samples together.
- Do not rely on EM-style isoform reference quantification.
- Do not discard rare but structurally unique isoforms solely because they are sample-specific.
- Do discard isoform structures observed only once across all samples, unless explicitly configured otherwise.
- Do not let unreliable 5-prime starts or 3-prime ends define isoform identity.
- Do not over-cluster biologically distinct isoforms with similar internal structure.
- Keep outputs compatible with existing AltAnalyze3 splicing, isoform ratio, and differential comparison workflows.

## High-Level Architecture

Use a three-stage pipeline:

1. Per-sample candidate reduction.
2. Cross-sample canonical isoform resolution.
3. Per-sample molecule and h5ad remapping.

The key distinction is that read annotation and final isoform identity are separate. Existing `gff_process` logic should continue to annotate structures. A new canonical resolver should decide which annotated structures represent the same final isoform across the cohort.

## Stage 1: Per-Sample Candidate Reduction

### Inputs

- Raw-read GFF or GFF.GZ from `isoform_structure_extract.py`.
- Raw molecule h5ad from `isoform_structure_extract.py`.
- Ensembl exon model file, e.g. `Hs_Ensembl_exon.txt`.
- Optional reference GFF/GTF for known transcript models.

### Behavior

For each sample independently:

- Parse raw-read GFF streamingly.
- Use the existing AltAnalyze3 annotation logic to assign each read to:
  - gene
  - strand
  - AltAnalyze exon/intron structure string
  - raw full structure including terminal exons
  - known Ensembl transcript match, where applicable
- Build a terminal-insensitive structure key by trimming uncertain terminal coordinate differences.
- Collapse exact duplicate structures and safe within-sample containments.
- Track the best exemplar molecule for each candidate.
- Preserve terminal coordinate ranges separately from the collapse key.
- Write molecule-to-candidate mappings for later remapping.

### Candidate Identity

The preferred structure key should be the AltAnalyze token string already used elsewhere in AltAnalyze3, for example:

```text
E1.1|E2.2|E2.4
```

The key should not include unreliable transcript start and end positions. Terminal start and end observations should be retained as metadata.

### Per-Sample Candidate Output

Write compact sample-local outputs next to the source GFF/h5ad:

```text
<sample>.isoform_candidates.tsv.gz
<sample>.molecule_to_candidate.tsv.gz
<sample>.candidate_counts.h5ad
```

Recommended `isoform_candidates.tsv.gz` columns:

```text
sample
gene
strand
candidate_id
structure_key
raw_structure_key
first_exon_start_min
first_exon_start_max
last_exon_end_min
last_exon_end_max
junction_count
exon_token_count
read_count
cell_count
barcode_count
is_ensembl
ensembl_transcript_id
best_exemplar_molecule
best_exemplar_score
```

Recommended `molecule_to_candidate.tsv.gz` columns:

```text
sample
molecule_id
barcode
gene
candidate_id
structure_key
read_count
```

For a single-cell BAM extraction where each molecule appears once, `read_count` will usually be `1`, but the column allows future support for molecule-level consolidation.

### Within-Sample Representative Selection

Candidate representative scoring should prefer:

1. Ensembl transcript match.
2. Higher read count.
3. More cell/barcode support.
4. More internal junctions or exon tokens.
5. Longer observed transcript span.

Example scoring tuple:

```python
(
    is_ensembl,
    read_count,
    cell_count,
    junction_count,
    observed_span,
    -terminal_uncertainty,
)
```

Within-sample collapse should be conservative. A shorter structure can collapse into a longer candidate only if:

- the shorter internal token sequence is a contiguous subsequence of the longer sequence;
- all shorter junctions are present in the longer candidate;
- the longer candidate has equal or stronger evidence;
- the collapse does not remove an internal alternative exon, intron retention, alternative donor, or alternative acceptor event.

## Stage 2: Cross-Sample Canonical Isoform Resolution

### Inputs

- All per-sample `*.isoform_candidates.tsv.gz`.
- Reference Ensembl/Gencode transcript structures annotated into the same AltAnalyze token space.
- Optional sample metadata.

### Behavior

Resolve final isoform IDs gene by gene. The resolver should never need all genes in memory at once.

For each gene:

1. Load all observed sample candidates for that gene.
2. Add reference Ensembl/Gencode isoform models for that gene.
3. Group exact terminal-insensitive structure matches.
4. Assign Ensembl IDs to exact or equivalent known structures.
5. Select observed representatives for novel structures.
6. Collapse safe lower-evidence contained structures into stronger parent structures.
7. Retain unique novel structures that have sufficient evidence.
8. Drop structures observed only once across all samples, unless they are protected by configuration.

### Final Isoform Priority

Reference isoforms should be preferred when structure-equivalent to observed candidates. For novel observed candidates, representative selection should consider both structure and expression evidence.

Recommended score:

```text
final_score =
    ensembl_priority
  + log1p(total_reads)
  + log1p(sample_count) * sample_weight
  + log1p(cell_count) * cell_weight
  + internal_junction_count * structure_weight
  + full_length_bonus
```

The score is used for representative selection, not for loose clustering. Two abundant structures that differ by a biologically meaningful internal event should remain separate even if their scores are similar.

### Cross-Sample Collapse Rules

Collapse candidate A into final candidate B only when:

- A and B are from the same gene and strand.
- A is an exact structure match to B, or A is a safe internal contiguous subsequence of B.
- Every junction in A exists in B.
- B has stronger aggregate support or is an Ensembl-prioritized model.
- A is low evidence relative to B, or A is terminal-truncated relative to B.

Do not collapse when:

- A contains a junction absent from B.
- A skips an internal exon present in B.
- A represents intron retention not present in B.
- A and B differ by internal donor or acceptor choice.
- A is independently abundant, even if shorter.

### Singleton Filtering

Default filtering:

```text
drop if total_reads_across_all_samples <= 1
```

Optional stricter filters can be supported:

```text
--min-total-reads 2
--min-samples 1
--min-cells 1
```

Do not use a default threshold that removes rare but replicated structures.

### Global Outputs

Write global resolver outputs under the existing `gff-output` directory:

```text
gff-output/final_isoforms.tsv.gz
gff-output/candidate_to_final.tsv.gz
gff-output/final_isoform_links.txt
gff-output/final_isoform_annotations.tsv
gff-output/final_isoform_resolver_stats.tsv
```

Recommended `final_isoforms.tsv.gz` columns:

```text
final_isoform_id
gene
strand
structure_key
source
ensembl_transcript_id
total_reads
sample_count
cell_count
candidate_count
is_novel
is_reference
representative_sample
representative_candidate_id
first_exon_start_min
first_exon_start_max
last_exon_end_min
last_exon_end_max
```

Recommended `candidate_to_final.tsv.gz` columns:

```text
sample
gene
strand
candidate_id
candidate_structure_key
final_isoform_id
final_structure_key
collapse_reason
collapse_score
candidate_total_reads
final_total_reads
```

`final_isoform_links.txt` should be compatible with current downstream expectations where possible. If compatibility is not clean, update downstream functions to consume `candidate_to_final.tsv.gz` directly.

## Stage 3: Per-Sample h5ad Remapping

### Inputs

- Original raw molecule h5ad or candidate-count h5ad.
- Per-sample `molecule_to_candidate.tsv.gz`.
- Global `candidate_to_final.tsv.gz`.

### Behavior

For each sample independently:

- Load one h5ad at a time.
- Build a sparse column mapping from molecule or candidate IDs to final isoform IDs.
- Aggregate counts into final isoform columns.
- Preserve barcode rows and sample metadata.
- Write a standardized final isoform h5ad.

### Sparse Matrix Strategy

The remapping should avoid dense matrices.

For raw molecule matrices:

1. Read h5ad as CSR.
2. Build `old_col_index -> final_col_index`.
3. Convert to COO or iterate CSR rows.
4. Replace column indices with final column indices.
5. Sum duplicates by converting back to CSR.

For candidate-count matrices:

1. Read candidate columns.
2. Map candidate columns to final isoform columns.
3. Aggregate sparse columns.

### Per-Sample Final Outputs

```text
<sample>-final_isoform.h5ad
<sample>-final_molecule_map.tsv.gz
<sample>-final_isoform_counts.tsv.gz
```

The h5ad `var` table should include:

```text
gene
structure_key
source
ensembl_transcript_id
is_novel
total_reads_global
sample_count_global
```

The h5ad `obs` table should preserve existing barcode/sample metadata.

## Suggested Command-Line Interface

Replace the current `collapse_transcript_associations.py` with a new resolver-oriented module. Suggested file:

```text
components/long_read/isoform-collapse/canonical_isoform_resolver.py
```

Suggested commands:

```bash
python canonical_isoform_resolver.py reduce-sample \
  --gff sample.gff.gz \
  --matrix sample.h5ad \
  --ensembl-exons Hs_Ensembl_exon.txt \
  --outdir sample_dir
```

```bash
python canonical_isoform_resolver.py resolve-global \
  --candidate-files candidates.list.txt \
  --reference-gff gencode.v45.annotation.gff3 \
  --ensembl-exons Hs_Ensembl_exon.txt \
  --outdir gff-output
```

```bash
python canonical_isoform_resolver.py remap-h5ad \
  --matrix sample.h5ad \
  --molecule-map sample.molecule_to_candidate.tsv.gz \
  --candidate-to-final gff-output/candidate_to_final.tsv.gz \
  --final-isoforms gff-output/final_isoforms.tsv.gz \
  --out sample-final_isoform.h5ad
```

An orchestration command can wrap all three:

```bash
python canonical_isoform_resolver.py run \
  --metadata mds_metadata_bam.txt \
  --ensembl-exons Hs_Ensembl_exon.txt \
  --reference-gff gencode.v45.annotation.gff3 \
  --genome-fasta genome.fa \
  --outdir gff-output
```

## Integration With Existing Workflow

Current user-facing workflow:

```python
sample_dict = isoa.import_metadata(
    metadata_file,
    include_hashed_samples=True,
    extract_from_bams=False,
    reference_model=ensembl_exon_dir,
)
isoa.pre_process_samples(metadata_file, barcode_cluster_dirs, ensembl_exon_dir)
isoa.combine_processed_samples(
    metadata_file,
    barcode_cluster_dirs,
    ensembl_exon_dir,
    gencode_gff,
    genome_fasta,
)
```

Recommended integration:

- `pre_process_samples` should continue to produce raw GFF/h5ad files and junction outputs.
- After raw GFF generation, add per-sample candidate reduction.
- `combine_processed_samples` should run global canonical resolution after all per-sample reductions exist.
- `export_isoform_h5ad` should use final standardized isoform mappings instead of the old `isoform_links.txt` behavior when final resolver outputs are available.
- Differential isoform comparisons should consume `<sample>-final_isoform.h5ad`.

This keeps splicing analysis untouched while replacing only the isoform abundance model.

## Data Backend

Use sorted gzipped TSV files for primary outputs because they are portable and easy to inspect. For global resolution, either in-memory gene chunks or SQLite can be used.

SQLite is recommended if candidate counts approach millions:

```sql
CREATE TABLE candidate (
    sample TEXT,
    gene TEXT,
    strand TEXT,
    candidate_id TEXT,
    structure_key TEXT,
    read_count INTEGER,
    cell_count INTEGER,
    is_ensembl INTEGER,
    ensembl_transcript_id TEXT
);

CREATE TABLE candidate_to_final (
    sample TEXT,
    gene TEXT,
    candidate_id TEXT,
    final_isoform_id TEXT,
    collapse_reason TEXT
);

CREATE INDEX idx_candidate_gene ON candidate(gene);
CREATE INDEX idx_candidate_structure ON candidate(gene, strand, structure_key);
CREATE INDEX idx_map_sample_candidate ON candidate_to_final(sample, candidate_id);
```

DuckDB would also work, but SQLite avoids introducing a heavier dependency.

## Algorithmic Notes

### Token Handling

Avoid ad hoc substring matching on raw strings. Tokenize structures first:

```python
tokens = [token for token in structure_key.split("|") if token]
```

Containment should be sequence-based:

```python
def is_contiguous_subsequence(short_tokens, long_tokens):
    if len(short_tokens) > len(long_tokens):
        return False
    for i in range(len(long_tokens) - len(short_tokens) + 1):
        if long_tokens[i:i + len(short_tokens)] == short_tokens:
            return True
    return False
```

This prevents accidental matches where one token string is a substring of another token.

### Expression-Aware Parent Choice

When multiple longer parents contain a shorter candidate:

1. Prefer exact Ensembl parent.
2. Prefer parent with all candidate junctions and highest total read support.
3. Prefer parent with broader sample support.
4. Prefer parent with smaller structural distance to the candidate.
5. If unresolved, keep the shorter candidate separate.

### Terminal Coordinates

Final isoform IDs should not change because one sample has a longer 5-prime or 3-prime terminal extent. However, terminal ranges should be retained for QC and potential future export:

```text
first_exon_start_min
first_exon_start_max
last_exon_end_min
last_exon_end_max
```

### Known Versus Novel

An isoform is known if its terminal-insensitive internal structure matches an Ensembl/Gencode model under the existing AltAnalyze annotation rules. It is novel if it has a supported internal structure not assigned to a known transcript.

Novel IDs should be deterministic:

```text
<gene>|NOVEL.<rank>
```

Rank should be stable by sorting:

```text
gene, strand, structure_key, total_reads descending, representative candidate
```

## Operational Efficiency

- Process samples independently and in parallel where possible.
- Process global resolution by gene.
- Store molecule maps as gzip TSV and stream them.
- Use sparse matrix aggregation for h5ad remapping.
- Do not parse reference GFF repeatedly for every sample if a precomputed reference structure cache exists.
- Add resumability by skipping outputs that already exist unless `--force` is supplied.
- Record per-stage stats so large runs can be audited without opening h5ad files.

Recommended stats:

```text
sample
raw_molecules
annotated_molecules
candidate_count
singleton_candidate_count
final_isoform_count
known_final_isoform_count
novel_final_isoform_count
dropped_singletons
remapped_matrix_nnz
```

## Validation Plan

Create small synthetic GFF/h5ad fixtures with known expected outputs.

Required cases:

1. Exact Ensembl match wins over observed candidate ID.
2. Terminal-truncated read collapses to stronger full candidate.
3. Alternative internal exon isoforms remain separate.
4. Alternative donor or acceptor isoforms remain separate.
5. Intron retention isoform remains separate from spliced isoform.
6. Rare novel structure with two or more total reads is retained.
7. Structure observed once across all samples is dropped.
8. Isoform expressed in only one sample is retained if sufficiently supported.
9. Molecule-to-final remapping preserves total counts per sample.
10. Sparse h5ad remapping preserves total counts per barcode.
11. Final isoform IDs are deterministic across repeated runs.
12. Existing junction/splicing outputs are unchanged.

## Implementation Milestones

### Milestone 1: Candidate Reducer

- Add structure-token utilities.
- Add per-sample candidate reducer using existing `gff_process` annotation.
- Write candidate and molecule map outputs.
- Add synthetic tests for token containment and singleton handling.

### Milestone 2: Global Resolver

- Load all sample candidate files by gene.
- Add reference isoform structures.
- Implement Ensembl-prioritized exact matching.
- Implement conservative expression-aware containment collapse.
- Write global final isoform and candidate-to-final tables.

### Milestone 3: h5ad Remapper

- Implement sparse molecule/candidate-to-final matrix remapping.
- Write standardized final isoform h5ad.
- Verify barcode and total-count preservation.

### Milestone 4: Workflow Integration

- Hook candidate reduction into `isoform_automate.pre_process_samples`.
- Hook global resolution and h5ad remapping into `isoform_automate.combine_processed_samples`.
- Preserve backward compatibility when final resolver outputs are absent.

### Milestone 5: Scale Testing

- Test on one large sample.
- Test on 5-10 samples.
- Test on full cohort with resumability.
- Record memory, runtime, candidate reduction ratio, and final isoform count.

## Open Design Decisions

- Whether the first implementation should use pure sorted TSV processing or SQLite for global resolution.
- Whether candidate-count h5ad should be mandatory or whether raw molecule h5ad remapping is sufficient.
- Whether to expose a default `--min-total-reads 2` singleton filter or leave it configurable with default `2`.
- Whether final known isoform IDs should use Ensembl transcript IDs directly or AltAnalyze-compatible structure IDs with Ensembl annotation metadata.
- Whether final novel isoforms should be named by stable rank or by a hash of `gene|strand|structure_key`.

## Recommended Defaults

```text
min_total_reads = 2
min_samples = 1
collapse_terminal_truncations = true
collapse_internal_alternatives = false
prefer_ensembl = true
use_sqlite_global_index = true
force = false
```

These defaults keep rare sample-specific isoforms when they have replicated read evidence, remove only true singleton structures, and avoid over-collapsing biologically distinct isoforms.
