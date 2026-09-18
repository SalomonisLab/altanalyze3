# Portable SNAF workflow and optional pyNeoQuant integration

## Implementation delivered

1. **HLA evidence matching:** pyNeoQuant normalizes allele names, scopes predictions
   to samples, retains multiple allele/predictor results, and rejects conflicting
   duplicate predictions. Explicit allele mismatches remain unmatched. Predictor,
   version and score meaning survive evidence-table import/export.
2. **SNAF/proteomics interface:** AltAnalyze3 exports stable source/candidate IDs,
   FASTA and manifest files, invokes the independent pyNeoQuant CLI optionally,
   and joins MS evidence back by sample, peptide and source. Explicit MS HLA
   annotations must also match. Decoys, entrapments and invalid/missing q-values
   do not count as passing evidence; competing sources remain visible.
3. **HLA genotype adapter:** supplied sample alleles take precedence over BAM
   inference. Inferred calls use the existing bam2hla engine, with no-calls and
   per-gene QC retained. Class-I only; this implementation does not establish new
   genotype accuracy claims.
4. **Nextflow orchestration:** reference indexing, per-BAM junction/intron counts,
   aggregation, annotation, ≥20-read filtering, optional HLA inference, SNAF and
   optional pyNeoQuant stages. Resource-related failures retry twice with increased
   memory/time; other failures stop the run and preserve resumable completed work.
5. **Distribution:** separate wheels and container recipes, CPU dependency locks
   for Linux x86_64/aarch64, model/reference manifest tools, Galaxy wrappers/tests,
   a draft-07 parameter schema, and CI definitions. pyNeoQuant remains a separate
   library and image, absent from AltAnalyze3's required dependencies.

## Install

Python 3.11 is the deployment baseline. From the AltAnalyze3 repository root:

```bash
pip install -e '.[snaf]'
altanalyze3-neo --help
```

For the optional proteomics stage, install pyNeoQuant **separately** from its
checkout, wheel, or future published distribution. Version 0.1.0a1 adds the
`fasta_record` manifest selector required by this integration:

```bash
pip install /path/to/pyneoquant-0.1.0a1-py3-none-any.whl
pyneoquant schema-version
```

The two packages may live in different environments. Use `--executable` with the
pyNeoQuant environment's executable; AltAnalyze3 communicates through files and
subprocesses, not a Python import. The Nextflow stage uses a separate container.

## Export, analyze and return MS evidence

```bash
altanalyze3 snaf --juncounts counts.tsv --hla hla.tsv --db_dir reference \
  --genome_fasta hg38.fa --output snaf --export_proteomics

# Or export an existing SNAF report:
altanalyze3-neo export --candidates snaf/T_candidates/T_antigen_candidates_*.txt \
  --predictor MHCflurry --canonical-fasta human_proteome.fasta --outdir bundle

altanalyze3-neo proteomics --bundle bundle --psm-table sample_psms.tsv \
  --sample S1 --outdir ms/S1 --executable /path/to/pyneoquant
```

Outputs include `source_manifest.tsv`, `candidates.fasta`, `candidates.tsv`,
provenance JSON, the pyNeoQuant analysis directory and `candidates_with_ms.tsv`.
`fasta_record` selects exactly one record from a shared FASTA, so metadata from
one junction cannot be assigned to every other sequence in that file.

The SNAF export contains **reported candidate peptides**, not reconstructed
full-length proteins. The source/event association is retained; junction geometry
is not inferred from a peptide sequence. A matching spectrum supports the
peptide/source association, not a particular presenting HLA allele. Including the
canonical reference in the manifest is recommended to expose competing sources.

Inputs must agree on sample identifiers. MS tables may include multiple samples;
joins preserve sample identity. Candidate rows without qualifying evidence are
retained. If a search engine used different sample labels, map those labels before
running this step. `--q-threshold` defaults to 0.01.

## Optional Galaxy iPepGen exports

Use `altanalyze3 snaf ... --galaxy_integration` or `altanalyze3-neo galaxy-export`
to create the additional input datasets required by the supplied OneClick iPepGen
workflow. The export is off by default and independent of pyNeoQuant. It includes
search FASTA, matched per-sample PepQuery/IEDB collections, accession mappings,
length eligibility, and annotation schemas. See [GALAXY_IPEPGEN.md](GALAXY_IPEPGEN.md)
for exact files, workflow connections, coordinate handling and verification.

## HLA genotype and binding interfaces

```bash
altanalyze3-neo hla --sample S1 --bam S1.bam --supplied supplied_hla.tsv \
  --output S1.hla.tsv --qc S1.hla_qc.tsv
altanalyze3-neo combine-hla --inputs S1.hla.tsv S2.hla.tsv --output cohort.hla.tsv

# Predict peptide binding for the external pyNeoQuant adapter contract:
altanalyze3-neo bind --evidence peptides.tsv --hla cohort.hla.tsv \
  --method MHCflurry --output binding.tsv
pyneoquant merge-hla --evidence evidence.tsv --predictions binding.tsv --out annotated.tsv
```

Supplied HLA tables use a header followed by `sample<TAB>allele1,allele2,...`.
Missing loci appear as no-calls in QC. `--require-all` requires A/B/C calls; a
sample with no usable alleles always fails. Expression-status suffixes require
explicit review instead of being silently collapsed into an ordinary allele.
Binding output records whether ranks are MHCflurry presentation percentiles or
NetMHCpan EL percentiles. Missing models fail in offline mode, and external
NetMHCpan process failures propagate to the caller.

## Nextflow

The pipeline is `altanalyze3/components/snaf/nextflow/main.nf` from the repository
root. It requires Nextflow ≥25.04.8. Samplesheet paths are resolved relative to
the samplesheet. Sample/cohort IDs use letters, digits, `_` and `-`; BAM filenames
may retain their original names. Inputs must be coordinate-sorted and indexed.

BAM samplesheet:

```csv
id,bam,bai,psm
S1,S1.bam,S1.bam.bai,S1.psms.tsv
S2,S2.bam,S2.bam.bai,
```

```bash
nextflow run altanalyze3/components/snaf/nextflow/main.nf -profile local \
  --mode bam --input samples.csv --gtf reference.gtf \
  --db_dir snaf_reference --genome_fasta hg38.fa --hla supplied_hla.tsv \
  --infer_hla --outdir results

# Optional MS stage; pyNeoQuant must be installed in the local environment:
# append --with_pyneoquant

# Self-contained real BAM/counting test; no scientific model downloads:
nextflow run altanalyze3/components/snaf/nextflow/main.nf -profile test,local
```

`--stop_after_counts` runs only the BAM-to-filtered-matrix portion. Matrix mode
`--mode full` uses `id,juncounts,hla`; with proteomics, provide one cohort and a
separate `--ms_input` CSV with `id,psm`. Existing `ts` and `surface` modes remain.
The surface samplesheet still includes `freq`; this does not alter the newer
standalone SNAF-B CLI's ability to derive its own frequency table.

GTF/GFF3 indexing supports quoted/custom exon identifiers and exon records lacking
an exon ID. References and SNAF's Alt91/healthy-tissue databases must use compatible
assembly and junction naming. This work does not regenerate the healthy-tissue
reference or establish equivalence for an arbitrary replacement annotation.
The corrected versus legacy annotation mode is explicit in workflow parameters.

## Containers and model provenance

From the repository root:

```bash
pip install build
bash altanalyze3/deployment/containers/build.sh
```

Build pyNeoQuant using its own `deployment/containers/build.sh`. The tags
`altanalyze3-snaf:0.1.3` and `pyneoquant:0.1.0a1` are **local build tags**. After
publication, pass immutable registry digests through `--container` and
`--pyneoquant_container`; replace the Galaxy runtime image references too.

Dependency lock regeneration requires `uv`:

```bash
bash altanalyze3/deployment/containers/lock.sh
```

Locks pin transitive versions and hashes for Linux x86_64/aarch64 and use CPU
PyTorch. Both platform resolutions passed during development. Container builds
also run `pip check`. Model downloads are separate from image creation. Record
and verify prepared reference/MHCflurry bundles before using them:

```bash
python -m altanalyze3.components.neoantigen.resources record \
  --reference reference --mhcflurry mhcflurry_models --manifest resources.json
python -m altanalyze3.components.neoantigen.resources verify \
  --reference reference --mhcflurry mhcflurry_models --manifest resources.json
```

For the CLI, set `MHCFLURRY_DOWNLOADS_DIR` to the prepared model directory. In
Nextflow, pass `--mhcflurry_models /path/to/models`; the workflow stages the
directory and sets the environment variable inside the prediction task. The
bundled DeepImmuno and bam2hla artifacts are checked in the wheel and recorded in
`containers/bundled_models.json`. NetMHCpan remains a user-provisioned binary.

Apptainer/Singularity can consume a published Docker image or a locally built SIF;
pass a SIF path through `--container`. A conversion recipe is included. The CI
container job builds a SIF from a Docker archive and exercises the same BAM fixture.

## Galaxy and CI

AltAnalyze3 wrappers live in `deployment/galaxy`; pyNeoQuant's wrapper lives in
its own repository. The SNAF prediction wrapper accepts a data-only tar archive
with `Alt91_db/` and `controls/` at its root. Archive traversal/links are rejected.
The pyNeoQuant wrapper accepts one FASTA and a manifest pointing at
`candidates.fasta`; CLI/Nextflow support richer manifests including canonical data.

```bash
planemo lint altanalyze3/deployment/galaxy/*.xml
planemo test altanalyze3/deployment/galaxy/*.xml
```

The checked-in tests cover the lightweight tools plus invalid-reference rejection
for full prediction. The main prediction acceptance test needs the real reference
and provisioned models. `deployment/tests/galaxy_commands.py` also runs the actual
Cheetah templates and fixture assertions without requiring a Galaxy server.

CI definitions cover interface/regression tests, wheel resources, Galaxy lint,
Nextflow 25.04.8/25.10.0, and an explicitly dispatched Docker/Apptainer runtime job.
They have not been submitted to or executed by GitHub Actions in this session.
This is an nf-core-style implementation, not an accepted nf-core release.

## Verification and current limitation

See [VALIDATION.md](VALIDATION.md) for commands, results and the container-storage
blocker. No community submission, registry push, or production cluster deployment
was performed as part of this local implementation.
