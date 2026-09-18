# Local validation — September 15, 2026

## September 17 follow-up: OneClick iPepGen exports

- The extracted bundled contract matches the user-supplied `.ga` checksum and settings.
- Neoantigen/interface suite: **15 passed, 1 skipped**. The skipped test requires
  a separately installed pyNeoQuant executable; the Galaxy path has no such dependency.
- All **six AltAnalyze3 Galaxy wrappers** pass Planemo lint without warnings and
  rendered-command tests, including two-sample collection discovery for the new exporter.
- Nextflow 25.10.0 stub runs pass with Galaxy exports enabled in BAM mode and
  disabled in matrix mode; expected presence/absence of the export directory was checked.
- Main SNAF CLI help exposes the new opt-in flags; the importable export workflow's
  tool/version/output/connection references are checked against the local wrapper.
- Wheel building and the bundled contract/resource check pass.

See [GALAXY_IPEPGEN.md](GALAXY_IPEPGEN.md) for outputs and connections. These are
local export/orchestration checks; no Galaxy server, FragPipe, PepQuery or IEDB
prediction was executed. This follow-up does not establish the older container
or scientific cross-platform acceptance milestones below.

## Results

| Check | Observed result |
|---|---|
| pyNeoQuant complete test suite | 161 passed; independent Python environment |
| AltAnalyze3 neoantigen interfaces | 9 passed, including external pyNeoQuant roundtrip, sample/allele binding, no-calls, and filter boundaries |
| Existing SNAF modernization tests | 8 passed; plus 2 new tests passed for NetMHCpan failure propagation and missing offline models |
| Real BAM pipeline | Passed on Nextflow 25.04.8 and 25.10.0: GTF indexing, junction/intron counts, aggregation and filtering |
| BAM expected counts | S1/S2 junction counts 20/19; S1 has 3 reads at each intron boundary; exactly one junction passes the any-sample ≥20 filter |
| Retry recovery | Both BAM junction tasks intentionally exited 137 on attempt 1, then succeeded on retry; downstream workflow completed |
| Resume | All 7 completed upstream tasks were reused from cache with `-resume` |
| Prediction orchestration | BAM/HLA/SNAF and matrix/SNAF/optional-pyNeoQuant/merge stub runs passed; surface stub passed, including identical input basenames |
| Optional Nextflow proteomics | Real independent pyNeoQuant CLI and real AltAnalyze3 return merge passed using a pre-exported candidate bundle; expected MS-supported candidate retained |
| Galaxy XML | All 6 wrappers passed Planemo lint without warnings; all 6 rendered Cheetah command tests passed |
| Packaging | Both wheels built; required model/reference assets and independent-package boundary checked |
| Dependency resolution | Hash locks resolved for Linux x86_64 and aarch64, Python 3.11 |

The existing modernization tests and two new failure tests were executed in separate
runs. The interface test uses the separately installed pyNeoQuant CLI when
`PYNEOQUANT_EXECUTABLE` is set; otherwise that one test is explicitly skipped.
Selected test logs are retained in [validation_logs](validation_logs).

## Reproduce

Commands below start at the AltAnalyze3 repository root, with the documented
Python dependencies, Nextflow and (for the optional checks) pyNeoQuant installed.

```bash
PYNEOQUANT_EXECUTABLE=/path/to/pyneoquant python -m pytest \
  altanalyze3/components/neoantigen/tests \
  altanalyze3/components/snaf/tests/test_snaf_modernization.py -q

NXF_VER=25.04.8 bash altanalyze3/components/snaf/nextflow/tests/smoke.sh
NXF_VER=25.10.0 bash altanalyze3/components/snaf/nextflow/tests/smoke.sh
bash altanalyze3/components/snaf/nextflow/tests/smoke.sh \
  -c "$PWD/altanalyze3/components/snaf/nextflow/tests/retry.config"

python -m build --wheel
python altanalyze3/deployment/tests/check_artifacts.py \
  --wheel dist/altanalyze3-0.1.3-py3-none-any.whl
planemo lint altanalyze3/deployment/galaxy/*.xml
python altanalyze3/deployment/tests/galaxy_commands.py \
  altanalyze3/deployment/galaxy/*.xml /path/to/pyNeoQuant/deployment/galaxy/*.xml
```

The self-contained BAM fixture occupies approximately 32 KB on this filesystem.
The fixture's SNAF database directory is intentionally a stub; it cannot support
scientific prediction. `smoke.sh` runs real upstream counting followed by a stub
prediction workflow. To exercise the optional MS stages with real execution:

```bash
altanalyze3-neo export \
  --candidates altanalyze3/deployment/galaxy/test-data/candidates.tsv \
  --predictor MHCflurry --outdir /tmp/snaf-bundle
nextflow run altanalyze3/components/snaf/nextflow/tests/proteomics_smoke.nf \
  -c altanalyze3/components/snaf/nextflow/nextflow.config \
  --bundle /tmp/snaf-bundle \
  --psm altanalyze3/deployment/galaxy/test-data/evidence.tsv \
  --sample S1 --outdir /tmp/snaf-ms-test
```

Check `merge_ms/S1.candidates_with_ms.tsv` for `peptide_source_supported` with
one distinct spectrum and q-value 0.001. This harness may report an unused
SNAF_TS configuration selector because it runs only the two proteomics processes.

## Fixes established by executable checks

- Standard quoted/custom GTF exon identifiers and missing exon IDs are handled.
- Intron counting retains the caller-owned indexed reference, allowing concurrent BAM jobs.
- The stage adapter supplies the intron counter's strand enum and numeric logger level.
- Aggregation resolves the annotation companion for both `.bed` and `.bed.gz` inputs.
- Filtering stages its input separately and rejects direct in-place file overwrites.
- Surface inputs with identical basenames are staged in separate directories.
- HLA merging does not assign another allele's rank or another sample's predictions.
- A shared FASTA manifest selects each source record explicitly, preventing metadata reassignment.

## Unverified acceptance work and container blocker

Fresh Docker builds were attempted, including the small independent pyNeoQuant
image, but Colima's internal Linux storage failed with **`no space left on device`**.
The host filesystem had free space; the Docker VM was full. No user images,
volumes, or cache were pruned, and the VM was not restarted. Consequently the new
images, their runtime `pip check`, and Docker/Apptainer execution are unverified.
Freeing or expanding that VM's storage is required before rerunning the build
scripts and the container runtime tests. Registry tags remain local build targets.

The following were not established by these local tests:

- A full biological BAM-to-SNAF prediction run with compatible Alt91/healthy-tissue
  references and provisioned MHCflurry/DeepImmuno models.
- New HLA genotype accuracy or peptide-binding calibration on biological cohorts.
- Galaxy server/runner acceptance (`planemo test`), MSI deployment, or IUC review.
  Rendering XML commands locally does not validate the Galaxy job runner.
- Hosted GitHub Actions execution, nf-core zero-warning lint/community acceptance,
  image publication, Dockstore/Zenodo registration, or CCHMC cluster deployment.
- Galaxy-versus-Nextflow scientific equivalence, raw MS search benchmarks, or the
  separate iPepGen/FragPipe/PepQuery deliverables.

pyNeoQuant remains optional, separately packaged, and absent from AltAnalyze3's
required dependencies. Its checkout was not initialized as a Git repository or
published during this work.
