# Synthetic isoforms from sample-level junction evidence

SNAF-B now has an evidence-driven prediction path in `surface/evidence_isoform.py`.
It reconstructs candidate transcripts from reference annotation and measured short-read
junction combinations, then ranks their junction-spanning ORFs. An optional small linear
ranker learns candidate preferences from long-read associations in other genes. The learned
method is now the default for `altanalyze3 snaf-b --mode short_read` and the standalone CLI.
The explicit evidence rule is available as a second option.

The inference path does not read a long-read catalog. The existing reference-only API
retains its behavior when neither `evidence` nor `first_exon=True` is supplied.

## Construction

1. Start from each reference transcript of the gene.
2. Search combinations of measured junctions with a nonempty **joint** sample intersection.
   Pairwise co-occurrence alone is insufficient for a three-event hypothesis. A sample-specific
   run requires that sample to support every added novel junction.
3. Reject overlapping or competing splice gaps and verify that every requested junction
   remains present after all edits. Adjacent exon/intron coverage markers are not splice gaps.
4. Identify complementary novel-exon arms by their inner coordinates within the same local
   intron. Their outer splice sites may skip annotated exons. Default exon length is at most
   500 bases; longer exons require joint detection in at least three sample columns.
5. Reject extensions across intervening annotated exons. Also reject newly inserted blocks
   that cover whole annotated introns without explicit retention-boundary evidence, unless
   that continuous exon is already annotated. This check includes AltAnalyze's historical
   mRNA exon annotations, which can be absent from Ensembl transcript models.
6. Limit an exon extension without a second measured splice boundary to 500 bases. This is
   a configurable construction prior, not a biological assertion that longer extensions
   cannot exist. Selected first exons have the separate fallback described below.
7. Translate the resulting transcripts, retain junction-spanning ORFs, and rank hypotheses.
   NMD is a risk annotation in evidence mode; it is not a measured transcript fate.

Reference transcript starts and ends remain priors. Junction counts do not determine complete
5′/3′ ends, translation initiation, or molecular phasing. Several isoforms may coexist in one
sample. Preserve alternatives for downstream inspection.

### First exons with an unknown 5′ boundary

Set `first_exon=True` in the prediction API, or `first_exon=true` on a target row in
the CLI input, to select that junction's donor as the end of the first exon. This
selection removes upstream transcript structure while preserving compatible downstream
structure. A junction alone does not establish that its donor exon is a first exon.

If the donor falls in an annotated first exon from a reference transcript on the same
strand, preserve that exon's annotated 5′ boundary. Otherwise assume the **entire first
exon is 250 nt long**, using 1-based inclusive coordinates:

- Plus strand: `[donor − 249, donor]`.
- Minus strand: `[donor, donor + 249]`.

The output field `first_exon_boundary_source` records `annotated`, `assumed_250nt`, or
`reference` (ordinary backbone inheritance). The GTF also records this provenance.
The assumed boundary is not an inferred transcription start site; translation initiation
is still selected from the reconstructed sequence. Candidates that cannot contain all
250 bases at a chromosome boundary remain unresolved. Splice compatibility and exon
integrity checks still apply. The SLC24A4 annotated 356-nt first exon stays 356 nt.

This fallback was added after the September benchmark; its accuracy has not yet been
measured in that benchmark.

### Search limits

Defaults: three reads for detection, up to four simultaneous edits, 16 candidate neighbors,
beam width 24, 160 exon chains per search depth, and the three longest distinct ORFs per chain
plus a non-NMD alternative when needed. ORFs must have at least 90 nucleotides. The existing
200,000-nucleotide transcript bound also applies. Hypotheses from smaller search depths are
preserved when a larger search reaches its cap. This is bounded inference, not exhaustive assembly.

## Python API

For the complete surface-antigen workflow, call after `snaf.initialize(...)` and
`surface.initialize(...)`:

```python
from altanalyze3.components.snaf import surface

surface.run(
    membrane_tuples, outdir="results", prediction_mode="short_read",
    junction_counts=complete_junction_dataframe,  # all measured junctions, all samples
    isoform_method="learned",                    # default when counts are supplied
    genome_fasta="genome.fa",                    # optional; Alt91 gene FASTA otherwise
    first_exon_junctions=selected_first_exon_uids, # optional UID list
)
surface.generate_full_results(
    outdir="results", freq_path=frequency_table,
    mode="short_read", validation_gtf=None,
)
```

Use `isoform_method="evidence"` for the explicit two-edit rule. Existing Python calls that
provide no count matrix retain legacy reconstruction. Standalone lower-level functions
continue to use the explicit rule unless a model is supplied:

```python
from altanalyze3.components.snaf.surface.predict_isoform import predict_isoform
from altanalyze3.components.snaf.surface.evidence_isoform import (
    Evidence, load_exon_annotation, generate_hypotheses, rank_candidates, load_ranker,
)

# Coordinates: 1-based inclusive exon boundaries, ascending on both strands.
# Every count vector uses the same sample order; use all measured junctions of the gene.
evidence = Evidence(junction_counts, min_reads=3, sample_index=None)
annotation = load_exon_annotation(exon_file, {gene})[gene]
fitted_model = load_ranker()  # bundled portable linear model; no long-read input

prediction = predict_isoform(
    reference_models, lo, hi, fetch_forward_sequence,
    evidence=evidence, exon_annotation=annotation,
    ranking_model=fitted_model,  # omit for the explicit sample-support rule
    junction_label=junction_id,
)

# Retain scores, construction labels, and supporting sample indices:
candidates = generate_hypotheses(
    reference_models, (lo, hi), fetch_forward_sequence, evidence,
    exon_annotation=annotation, junction_label=junction_id,
)
ranked = rank_candidates(candidates, fitted_model)
```

`predict_isoforms` also accepts `evidence_by_gene={gene: Evidence(...)}` and
`exon_annotation_by_gene`. A missing evidence row is not silently converted into a measured
zero row. A cohort-level reference-only hypothesis may still be returned for such a target;
the CLI labels it explicitly. A sample-specific run abstains without target support.

## CLI and output

```bash
SNAF_OFFLINE=1 python -m altanalyze3.components.snaf.surface.predict_evidence \
  --junctions targets.tsv \
  --models reference_models.pkl \
  --evidence cohort_counts.pkl \
  --exon-annotation Hs_Ensembl_exon.txt \
  --genome genome.fa \
  --method learned \
  --top-k 5 \
  --out predictions
```

`targets.tsv` columns: `junction_id`, `gene`, `chrom`, `junction_start`, `junction_end`.
Reference model pickle: `{gene: {transcript: (strand, [(start, end), ...])}}`.
The optional `first_exon` column selects first-exon targets (`true`/`false` or `1`/`0`).
The bundled learned ranker is used by default; `--ranker` overrides it with a JSON model or
a fitted pickle bundle. Use `--method evidence` for the explicit rule. Add `--sample
SAMPLE_COLUMN` for sample-specific inference. Learned defaults to `--max-edits 4 --search-mode
union`; evidence defaults to `--max-edits 2 --search-mode depth`. Explicit search flags
override those defaults. The learned candidate pool preserves shorter searches and supports
more complex multi-junction hypotheses.

Outputs:

- `predictions.tsv`: ranked alternatives, full mRNA/CDS/protein, NMD annotation, parent
  transcript, applied junctions, construction category, and supporting sample columns.
- `top1.gtf`, `top1.mrna.fa`, `top1.protein.fa`: top hypothesis per target. Stable artifact
  identifiers incorporate the strand, exon chain, and CDS boundaries so alternative ORFs
  cannot overwrite one another.
- `manifest.json`: input paths, file stamps, configuration, and feature schema.

These are predicted reference artifacts for review and downstream SNAF-B processing.
They do not carry claims of long-read confirmation or experimental protein validation.

### Deployment checks

The installed local `altanalyze3` command loads the revised workflow. Integration tests
exercise both ranking options, first-exon reconstruction, sample restriction, protein
filtering, and final candidate reporting. The real SLC7A5, MPZL1, and SLC24A4 inputs
reproduce the revised structures with both deployed options. The portable bundled model
matches the fitted model's scores to floating-point precision, and package-data checks
include the model and these guides. Nextflow inputs and outputs have been updated; its
runtime was not tested in this environment because Nextflow is not installed.

## Learning and evaluation

### Why the learned method is the default

In the September 15 internal evaluation of **191 verified long-read cases**, the fixed
linear ranker optimized for combined protein/structure agreement had the highest mean
protein identity among the two deployed options:

| Method | Mean protein identity | Mean intron-chain F1 | Exact splice chains |
|---|---:|---:|---:|
| Legacy single/pair baseline | 0.6239 | 0.5856 | 20/191 |
| Explicit evidence, two edits | 0.6879 | 0.6598 | **38/191** |
| **Learned linear, joint objective (default)** | **0.7009** | **0.6662** | 34/191 |

Nested selection among the evaluated model/objective combinations reached 0.7012 mean
protein identity, 0.6672 F1, and 32/191 exact chains. Those nested-selection figures are
distinct from the fixed linear model chosen for deployment. Use the explicit option when
exact splice-chain recovery is the preferred criterion. The bundled final model was fitted
on 111 training genes, excluding all three exemplar genes. Its portable JSON contains feature
names, coefficients, scaling, training gene IDs, and the source model's SHA-256 fingerprint.
No per-gene identifiers or long-read observations are needed at inference.

The first-exon fallback was added after this benchmark. These numbers describe the earlier
internal evaluation, not a new accuracy estimate for the deployed first-exon extension.

### Three exemplar structures

| Example | Original construction | Revised construction |
|---|---|---|
| SLC7A5 | 11 exons, including a 177-nt novel exon at 87,843,643–87,843,819; 566 aa | Same structure and protein, now supported by joint co-detection of both splice arms |
| MPZL1 | A 23,038-nt first exon spanning 167,721,950–167,744,987 swallowed intervening annotated structure; 6 exons, 283 aa | First exon 167,721,950–167,722,242 plus an 85-nt novel exon 167,744,903–167,744,987; 7 exons, 265 aa. The measured upstream junction skips E4–E6 |
| SLC24A4 | First exon 92,322,581–92,323,960 merged E1–E3; 5,883-nt mRNA | Annotated alternative first exon 92,323,605–92,323,960 (356 nt), plus a longer annotated 3′ end; 10,091-nt mRNA. Same 141-nt novel exon, CDS, and 669-aa protein |

The 250-nt fallback does not replace SLC24A4's known first-exon boundary. Sample support
for the complementary arms was 68 columns for SLC7A5, 63 for MPZL1, and 31 for SLC24A4 at
the three-read threshold. These are sample columns, not necessarily distinct patients.

The evaluator compares legacy single/pair editing, local construction at one through four
edits, explicit evidence ranking, ridge regression, tree ensembles, gradient boosting, and
a linear ranker trained on within-junction candidate differences. Twenty features describe
sequence length, ORF position, NMD risk, reference structure, sample support, added junctions,
paired-exon geometry, and possible intron retention. Gene identifiers, long-read sequences,
and long-read exon chains are not predictor features.

Genes, including all their junctions, remain together in five outer folds. An inner gene split
chooses model and training objective. The objective is either intron-chain F1 or the mean of
that F1 and full-protein identity. Invalid target-to-long-read associations are excluded from
supervised fitting and exported in `truth_audit.tsv`. Unpredicted events count as zero in the
historical full-denominator comparison; a verified-truth subset is also reported.

Evaluation entry points in `components/snaf/dev/`:

- `prepare_isoform_evidence.py`: extract gene subsets from TSV counts or sample-by-junction
  H5AD. Coordinate aliases merge by maximum count, not sum.
- `benchmark_evidence_isoforms.py`: novel isoform evaluation, nested gene splits, per-event
  predictions, bootstrap intervals, fitted ranker, and manifest.
- `benchmark_known_isoforms.py`: separate structural control using annotated intron chains
  also observed in long reads. Models exclude the control gene from training. No protein
  accuracy is claimed from this chain-only control.
- `report_isoform_exemplars.py`: SLC7A5, MPZL1, and SLC24A4, including historical models,
  alternative hypotheses, sample support at multiple read thresholds, and GTF/FASTA output.

The September 15, 2026 evaluation is exploratory and internal. Gene-held-out results reduce
memorization but do not replace an independent cohort test after method selection. Sample
columns can include technical replicates; recurrence is not automatically patient recurrence.

## Biological interpretation

Sample co-detection supports a candidate combination but does not phase distant junctions
onto one molecule. Splice-graph assembly and transcript-end ambiguity are established concerns
in short-read reconstruction ([StringTie](https://www.nature.com/articles/nbt.3122)). NMD rules
have experimentally observed exceptions, so predicted NMD should not be described as certainty
of destruction ([boundary-dependent NMD experiments](https://pmc.ncbi.nlm.nih.gov/articles/PMC1084009/)).
