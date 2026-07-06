# SNAF modernization (AltAnalyze3)

This documents the changes that make the SNAF neoantigen workflow cross-platform,
faster, offline-capable, scalable past 100k features, BayesTS-integrated, and
nf-core compatible. Every change is covered by `tests/test_snaf_modernization.py`
(self-contained, synthetic data; also validated against a real 6,000-junction
GTEx slice during development).

## What changed

### 1. Scale (>100k junctions) & speed — `gtex.py`
- `multiple_crude_sifting_maxmin_fast` / `_prevalance_fast`: replace the
  `pd.concat([col]*n_samples)` dense tiling (3–7 full float64 copies of the
  junctions×samples matrix) with **numpy broadcasting** + a single boolean
  `cond_df`. Proven identical to the legacy output; runs at 250k junctions.
- Control DB read in **`backed='r'`** mode and subset on disk — only the
  intersecting rows enter RAM (not the whole multi-GB h5ad).
- `add_control` combine now uses **scipy.sparse** row-union + `hstack`
  (`_combine_two_sparse`) instead of `adata.to_df().join()` densification.
- GTEx prevalence counts computed **directly on the sparse matrix** (no
  `.toarray()`).
- Exon "in Ensembl db" test uses a **precomputed per-gene pair set** instead of
  3 regex compiles per junction.
- Half-normal MLE tumor-specificity is now the **closed form** `sqrt(mean(y²))`
  (clipped to the original (0,1) bound) instead of a per-junction optimizer.
- Set `SNAF_LEGACY_SIFTING=1` to force the original dense path (kept for
  equivalence testing).

### 2. Direct BayesTS integration — `bayests.py` (new)
- Faithful, importable adaptation of BayesTS (pyro, RNA modes `XY`/`Y`) — a joint
  model over all queried junctions returning per-junction sigma in [0,1].
- Runs on **CPU** (works on macOS/Windows), is seeded/deterministic, offline, and
  has no file-roundtrip/plot side effects.
- Wired into `gtex.add_tumor_specificity_frequency_table(df, method='bayesian')`,
  replacing the retired per-junction **pymc3/theano** path (which did not install
  on macOS/Windows). `torch`/`pyro` are optional (lazy import).

### 3. Cross-platform portability
- **Multiprocessing** (`snaf.py`): all 5 `mp.Pool` sites route through
  `_parallel_apply`, which uses a **fork context** on posix (so macOS works — its
  default `spawn` loses the runtime-set globals the workers need) and a **serial
  fallback** on Windows. `.iteritems()` → `.items()` (pandas 2).
- **MHC binding** (`binding.py`): MHCflurry is the cross-platform default; the
  loader is cached, fixes the post-download `NameError`, and honors
  `MHCFLURRY_DOWNLOADS_DIR` for offline models. netMHCpan output is parsed in
  **pure Python** (`_parse_netMHCpan_stdout`) — no `awk`/`grep`/`shell=True`.
- **Offline sequence retrieval** (`snaf.py`): UTR / novel-exon sequences read from
  a local indexed genome FASTA (`set_genome_fasta` / `SNAF_GENOME_FASTA` /
  `initialize(..., genome_fasta=...)`) instead of the UCSC DAS web service.
  `SNAF_OFFLINE=1` forbids the network entirely.
- Optional heavy deps (**tensorflow** via deepimmuno, **dash**, **torch/pyro**)
  import lazily so `import snaf` works without them.

### 4. CLI + nf-core — `cli.py`, `nextflow/`
- `altanalyze3 snaf-ts` — tumor-specificity only (DB-free, offline; mean + MLE +
  BayesTS).
- `altanalyze3 snaf` — full MHC-bound T-antigen pipeline.
- `altanalyze3 snaf-b` — surface / B-antigen pipeline (see §6).
- nf-core-style DSL2 modules/workflow under `nextflow/` (`SNAF_TS`, `SNAF`,
  `SNAF_B`; `--mode ts|full|surface`). All modules install altanalyze3 from GitHub
  (not PyPI) and run under `-profile local` with no conda/container.
- Reference auto-resolution: `reference.py` (`ensure_reference`) locates the bundle
  (tolerating the tarball `data/` nesting), validates every file up front, and
  either errors with the exact `curl` command or downloads it (`--download_ref`).
- Optional install extras in `pyproject.toml`: `snaf-bayes`, `snaf-immuno`,
  `snaf-surface`, `snaf`.

### 5. Speed / memory (results-preserving)
- **Tumor specificity vectorized** (`gtex.tumor_specificity_batch`): the per-junction
  `tumor_specificity(uid, ...)` loop (405k iterations on the 2.6M-junction cohort) is
  replaced by batched sparse ops — `mean` from `obs['mean']`, `MLE = min(sqrt(mean(cpm^2)),1)`
  over the sparse control matrix with no densification. Output is identical
  (mean diff 0, MLE diff < 1e-15); `snaf-ts` CPU on the small cohort dropped 3.8x.
- **macOS fork safety**: `snaf._mp_context` sets `OBJC_DISABLE_INITIALIZE_FORK_SAFETY`
  on darwin so fork workers don't abort once tensorflow/mhcflurry are loaded.

### 6. Surface / B-antigen — now ported & CLI-driven (`surface/`)
Previously Linux-only; now pure-Python and offline-first, with graceful degradation
(any missing optional dep is logged and the remaining steps still run):
- **Transmembrane topology**: TMHMM 2.0c binary → **pure-Python `tmhmm.py`**
  (`surface/alignment.py`); missing → the TM gate is skipped with a note. A legacy
  binary is still used if an explicit `--tmhmm_path` is given.
- **Pairwise alignment**: EMBOSS-needle EBI REST → **Biopython `PairwiseAligner`**
  (`surface/api.py`); the `emboss.py` REST client is retired.
- **UTR-event sequence** (batch path): `surface`'s `retrieveSeqFromUCSCapi` delegates
  to the T-antigen offline-first resolver (local genome FASTA / `SNAF_OFFLINE`).
- **Gene symbols**: `mygene` wrapped to degrade to Ensembl IDs if absent/offline.
- **Redundant work removed** (results-preserving): the validation GTF is parsed once
  (`_GTF_PARSE_CACHE`), the extracellular-region tables built once, and
  `is_support_by_est_or_long_read` memoized by `(uid, ORF, strict)` — cutting the 6x
  recompute across the style×overlap×stringency passes.

## Ported since (see §6): TMHMM 2.0c, EMBOSS-needle, Ensembl REST, MyGene
The surface/B-antigen pipeline is now pure-Python and CLI/nf-core driven. The
remaining Linux-only / network pieces are genuinely optional:
- **netMHCpan** binary (Linux only) — MHCflurry is the cross-platform default.
- **singularity + ggsashimi + samtools** sashimi plots (`downstream.py`, visualization).
- **Ensembl REST** in the interactive B-antigen Dash viewer only (batch path is offline);
  degrades gracefully. **MyGene** is optional (falls back to Ensembl IDs).

## Not addressed (data-format migration, separate effort)
SNAF still consumes the AltAnalyze2 `Alt91_db` reference (Ensembl-v91 exon table,
`mRNA-ExonIDs`, 2000-flank FASTA, v91 GTF). This is a **data** dependency, not
Python-2 code. The AltAnalyze3 replacements were mapped (long_read/gff_process,
annotation/junction_isoform, long_read/isoform_translation, iso2function/uniprot),
but swapping them in is a larger, separate migration.
