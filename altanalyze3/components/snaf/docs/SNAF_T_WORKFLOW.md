# SNAF-T (T-antigen) workflow

Predicts MHC-I–presented, immunogenic neoantigens from a splice-junction count matrix.
Cross-platform pure-Python (macOS/Linux/Windows), no compiled external binaries required.

## CLI

```bash
altanalyze3 snaf \
  --juncounts counts.txt --db_dir <SNAF_reference> --hla hla.txt \
  --binding_method MHCflurry --genome_fasta genome.fa \
  --min_samples 1 --cpus 8 --output out/
```

`snaf-ts` is the long-read (isoform) variant (same engine, adds the isoform reference).

| Flag | Meaning | Default |
|------|---------|---------|
| `--juncounts` | junction × sample count matrix (AltAnalyze format; `uid=extra` index is split, duplicates dropped) | required |
| `--db_dir` | SNAF reference root (`Alt91_db/`, `controls/GTEx_junction_counts.h5ad`, …) | required |
| `--hla` | per-sample HLA alleles | required |
| `--genome_fasta` | local hg38 FASTA → offline in-silico translation (replaces per-junction UCSC DAS HTTP) | recommended |
| `--gtex_db` | override control h5ad (honored; was silently ignored before) | Ensembl91 GTEx |
| `--download_ref` | fetch the reference bundle if missing, else error with the download command | off |
| `--min_samples` | keep only junctions expressed (per the sifting condition) in ≥ N tumor samples | 1 |
| `--cpus` | worker processes | all |

Environment toggles: `SNAF_FAST_NN=1` (pure-NumPy NN inference, see below), `MHCFLURRY_DOWNLOADS_DIR` (offline models).

## Pipeline stages

1. **initialize** — load reference; control matrix is bulk-read then subset in-memory (see perf #4).
2. **sift** (`multiple_crude_sifting_maxmin_fast`) — keep junctions tumor-specific vs GTEx (max−mean>t_min & mean<n_max), optionally drop Ensembl-documented junctions, apply `--min_samples` recurrence filter.
3. **translate** (kind=1) — in-silico ORF translation → candidate 9/10-mer peptides (offline via `--genome_fasta`).
4. **bind + immunogenicity** (kind=3) — MHCflurry presentation percentile per (peptide, allele) + DeepImmuno CNN immunogenicity.
5. **generate_results** — burden/frequency tables (stage0–3), gene-symbol + coordinate + tumor-specificity annotation, per-sample `T_antigen_candidates`.

## Speed optimizations (all validated — see `benchmarks/BENCHMARKS.md`)

| # | Optimization | Speedup | Identity |
|---|---|---|---|
| 1 | Antigen-processing **dedup**: MHCflurry recomputes the allele-independent processing score per pair; compute it once per unique peptide + reuse mhcflurry's own affinity/logistic/percentile | 4.3× (grows with #HLA) | bit-identical (0.0) |
| 2 | DeepImmuno: cache model + aaindex tables + memoize scores (was reloaded every call) | 2192× | bit-identical (0.0) |
| 3 | Tumor-specificity: vectorized/chunked sparse (was per-uid loop) | 100× | 1.11e-16 |
| 4 | Control-subset load: bulk read + in-memory subset (anndata-0.12 backed fancy-index was pathological) | 26 min → 48 s | same rows |
| 5 | Offline hg38 translation via `--genome_fasta` (was per-junction UCSC DAS HTTP) | ~10⁵× | identical peptides |
| 6 | **Sift recurrence-filter O(n²) fix**: `set(keep)` was rebuilt per iteration of `--min_samples>1` (latent at default 1). Hoisted out. | 15–20 min stall → 15.7 s | byte-identical (valid unchanged) |
| 7 | generate_results: mygene dedup + persistent cache, drop per-candidate deepcopy | 2.46×/2.88× | byte-identical |

### Optional: `SNAF_FAST_NN=1` — pure-NumPy NN inference (no TensorFlow/PyTorch)
MHCflurry's processing + affinity nets and DeepImmuno's CNN are re-implemented in BLAS-backed NumPy
(`fast_infer_processing.py`, `fast_infer_affinity.py`, `deepimmuno/fast_cnn.py`). Matches stock to
float32 epsilon (~1e-7); the **called candidate set is byte-identical** (verified on the 55k×4 subset:
frequency_stage3 0 differing rows, T_candidates identical; only `binding_affinity` differs ≤0.9%
from Accelerate float32 rounding). Processing goes 11× faster (the 96.5%-of-runtime bottleneck).
Default is OFF (bit-identical stock path); enable for large cohorts where 3 h wall-clock matters.

## Cross-platform notes
- macOS forking after TF/MHCflurry load: `OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES` (set automatically for darwin).
- netMHCpan/MixMHCpred are academic-licensed and cannot ship; MHCflurry (Apache) is the distributable default. netMHCpan remains supported if the user installs it (`--binding_method netMHCpan --software_path`).
