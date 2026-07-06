# SNAF workflow documentation

Per-workflow docs for the modernized, cross-platform SNAF (pure-Python; no external binaries required).
See `../MODERNIZATION.md` for the overall change log and `benchmarks/` (in the test directory) for
before/after timings with identity validation.

| Workflow | CLI | Doc |
|---|---|---|
| T-antigen (MHC-I neoantigens) | `altanalyze3 snaf` / `snaf-ts` | [SNAF_T_WORKFLOW.md](SNAF_T_WORKFLOW.md) |
| B-antigen (cell-surface antigens) | `altanalyze3 snaf-b` | [SNAF_B_WORKFLOW.md](SNAF_B_WORKFLOW.md) |
| Control-reference precompute (BayesTS + GTEx stats) | `altanalyze3 snaf-precompute-control` | [CONTROL_PRECOMPUTE_WORKFLOW.md](CONTROL_PRECOMPUTE_WORKFLOW.md) |
| nf-core pipeline (all of the above) | `nextflow run … --mode {tantigen,surface}` | [../nextflow/README.md](../nextflow/README.md) |

## Key runtime toggles
- `SNAF_FAST_NN=1` — pure-NumPy NN inference (no TF/PyTorch); candidate-set byte-identical, ~1e-7 score rounding, big speedup on large cohorts.
- `--min_samples N` — recurrence feature filter (default 1 = analyze every sifted junction = original behavior).
- `--genome_fasta` — offline in-silico translation.
- `--max_bayests_percentile` — control-side tumor-specificity filter (requires a precomputed control-stats table).
