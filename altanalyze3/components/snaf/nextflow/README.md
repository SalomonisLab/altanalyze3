# SNAF — nf-core-style Nextflow workflow

DSL2 modules wrapping the `altanalyze3 snaf-ts` and `altanalyze3 snaf` CLI subcommands.
Cross-platform (macOS / Windows / Linux) and offline-capable.

## Layout
```
nextflow/
  main.nf                              # workflow (modes: ts | full)
  nextflow.config                      # manifest, resources, profiles
  conf/test.config                     # minimal test profile
  assets/samplesheet_test.csv          # example samplesheet
  modules/local/
    snaf_ts/{main.nf,meta.yml,environment.yml}   # tumor-specificity only
    snaf/main.nf                                 # full T-antigen pipeline
```

## Samplesheet
```csv
id,juncounts,hla
cohort1,/abs/path/junction_counts.txt,/abs/path/hla.txt
```
`hla` is required only for `--mode full`.

## Run — tumor specificity (DB-free, offline; adds mean + MLE + BayesTS)
```bash
nextflow run main.nf -profile local \
    --mode ts \
    --input samplesheet.csv \
    --control_h5ad /abs/path/GTEx_junction_counts.h5ad \
    --outdir results
```

## Run — full MHC-bound T-antigen pipeline
```bash
nextflow run main.nf -profile docker \
    --mode full \
    --input samplesheet.csv \
    --db_dir /abs/path/snaf_reference \      # contains Alt91_db/ and controls/
    --genome_fasta /abs/path/hg38.fa \       # optional: offline UTR sequence retrieval
    --outdir results
```

## Software provisioning (conda NOT required)
Nextflow can supply `altanalyze3` three ways — pick one; conda is optional:
- **Native (`-profile local`)** — simplest. Uses the `altanalyze3` already on your
  `PATH`. Install the current GitHub build first (the PyPI release lags the refactor):
  ```bash
  pip install -e /path/to/altanalyze3        # or: pip install git+https://github.com/SalomonisLab/altanalyze3.git
  nextflow run main.nf -profile local --mode ts ...
  ```
  For `--mode full` also install `mhcflurry tensorflow tf-keras` into that same env.
- **Conda (`-profile conda`)** — uses each module's `environment.yml` (both now install
  altanalyze3 from GitHub, not PyPI).
- **Container (`-profile docker|singularity`)** — requires a published `altanalyze3`
  image; the `community.wave.seqera.io/...` reference is a placeholder, so build/point
  it at your own image before using these profiles.

## Notes
- The full pipeline defaults to MHCflurry (cross-platform). netMHCpan (Linux only)
  is opt-in via `ext.args = '--binding_method netMHCpan --software_path ...'`.
- Set `SNAF_OFFLINE=1` (or pass `--genome_fasta`) to forbid any network access.
- `-profile test` runs the TS module on the example samplesheet (supply your own
  small `--control_h5ad`).
