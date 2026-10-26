Healthy hematopoietic splice and promoter survey
================================================

This workflow queries existing differential-splicing and transcript-usage tests.
It does not fit a new statistical model or substitute new reference annotations.
Run the modules from the repository parent, or set ``PYTHONPATH`` to that parent.
Use fresh output directories; completed results are immutable provenance inputs.

1. Link the complete healthy differential-splicing catalog to native isoforms::

    python -m altanalyze3.components.annotation.splice_isoform_survey \
      --counts RUN/isoform_combined_pseudo_cluster_counts.h5ad \
      --metadata RUN/mds_metadata_bam_healthy_fixed.txt --group healthy \
      --events RUN/analysis/healthy_diff/dPSI-cluster-events \
      --transcript-index RUN/isv/gff-output/gene_indexes_v2/transcripts.db \
      --min-reads 100 --q 0.05 --outdir SURVEY

The threshold is strictly greater than 100 pooled assigned reads for **each**
isoform in at least one healthy state; the states may differ. The complete count
feature/sample rosters are retained. All opposite exclusive junction carriers
are included, without a coding, NMD, protein-length, top-N or effect-size filter.
Carriers of both competing junctions remain explicit and do not create exclusive
pairs. Literal native IDs, including molecule/source suffixes, are preserved.
Unknown samples or required native structures reject execution.

2. Query existing reciprocal usage statistics and extract literal proteins::

    python -m altanalyze3.components.annotation.splice_survey_review \
      --survey SURVEY --run RUN \
      --usage RUN/analysis/healthy_diff/diff-cluster-ratio \
      --transcript-index RUN/isv/gff-output/gene_indexes_v2/transcripts.db \
      --literature LITERATURE_CATALOG.json --outdir REVIEW

The versioned primary-study catalog is supplied explicitly. The repository
provides ``components/annotation/resources/healthy_hematopoietic_splice_literature_20261020.json``.
Study reproduction must use the frozen catalog listed in the run script.
Both isoforms must have original usage FDR below the threshold and opposite
usage changes in the same state for a reciprocal-usage call. No new tests are
computed. The full literature-gene screen includes failed read criteria and
non-significant/untested cases. Protein metadata/FASTA discrepancies remain
visible; peptide sequence is never reconstructed from an annotation length.

3. Audit coding and transcript differences::

    python -m altanalyze3.components.annotation.splice_survey_correspondence \
      --survey SURVEY --review REVIEW \
      --donor-switches EXISTING_DONOR_PAIRED_SWITCHES.tsv \
      --outdir CORRESPONDENCE

``--donor-switches`` is optional. Its source is an existing
``isoform_state_switches.tsv`` and preserves that analysis's original selection
limits. It supplements the global event survey; a coding-only, DE-whitelisted
source is not a complete test of all isoforms. The protein alignment reuses
``iso2function.protein_isoform_audit.align_mapping`` (global BLOSUM62, gap-open
-10, gap-extension -0.5). Alignment gaps, terminal replacement, and substitutions
are evidence for review, not automatic literature identity or function calls.

4. Validate explicit literature correspondence and export the report::

    python -m altanalyze3.components.annotation.splice_survey_report \
      --survey SURVEY --review REVIEW --correspondence CORRESPONDENCE \
      --curation NATIVE_LITERATURE_CORRESPONDENCE.json \
      --gene-symbols EXISTING_Hs_Ensembl-annotations.txt.gz --outdir REPORT

Curation is a hashed input with native identities, cited primary mechanisms,
qualifications and executable sequence expectations. A changed expected coding
signature fails the report. This is not a generative literature classifier.
The report additionally flags the stricter rule requiring both partners to exceed
100 reads together in at least one state; broader pairs remain in the complete
table. Outputs include the HTML report, verified correspondence with original
statistics, all screened pairs, all screened genes, every native count feature
in those genes, and the complete genome-wide qualifying gene roster.

Interpretation
--------------

* Original event FDR is per state-versus-Others table, not a new global correction
  over all 80 tables. Event significance links junction carriers; it is not a new
  pair-specific significance test.
* Pooled gene-usage dominance reversal is descriptive. Original event/usage
  tests and optional donor-paired tests have different comparison units.
* Protein predictions and read support do not demonstrate translation, binding,
  or function in the measured cells. Published assays are identified separately.
* AltAnalyze exon-group numbers are not assumed to equal publication exon numbers.
* A new molecule/transcript ID does not establish a novel protein.
* The data survey is genome-wide within the existing test catalogs. The explicit
  literature catalog defines the scope of literature verification; it is not an
  exhaustive index of every published hematopoietic event.

The 2026-10-20 cluster run, scripts, frozen code, source hashes and curation live
under ``isoform_final_20261016/analysis/healthy_splice_literature_survey_20261020``.
