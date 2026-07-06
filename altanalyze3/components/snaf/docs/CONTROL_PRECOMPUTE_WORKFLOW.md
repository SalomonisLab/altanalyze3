# Control-reference precompute workflow (`snaf-precompute-control`)

BayesTS and the GTEx control summaries depend **only on the control reference**, not on the tumor
cohort. They should be computed **once per control** (shipped for the default Ensembl91 control, or
computed on first use of a new control and cached) and reused at runtime as a small lookup table —
instead of loading the multi-GB control h5ad and re-running BayesTS on every analysis.

## CLI

```bash
altanalyze3 snaf-precompute-control \
  --control_h5ad controls/GTEx_junction_counts.h5ad \
  --output controls/GTEx_junction_counts.snaf_stats.tsv.gz \
  --bayes_mode XY --bayes_epoch 2000 --bayes_batch 50000 --normal_cutoff 5
# --no_bayes to write only the fast stats (mean/std/mle/normal_prevalence)
```

## Output: one small per-junction table (`<control>.snaf_stats.tsv.gz`)

| Column | Meaning |
|---|---|
| `mean`, `std` | per-junction count mean/std across control samples |
| `mle` | closed-form half-normal tumor-specificity MLE (control-side) |
| `normal_prevalence` | fraction of control samples with count > cutoff |
| `bayests_sigma` | BayesTS per-junction sigma ∈ [0,1]; **lower = more tumor-specific** |
| `bayests_percentile` | global rank of `bayests_sigma` in (0,1] — the SNAF "p-value" |

## Runtime use (default tunable filter)

At analysis time SNAF filters sifted junctions on a default `bayests_percentile` threshold
(`--max_bayests_percentile`, user-tunable). Junctions **absent** from the control are maximally
tumor-specific and always kept. This drops the large fraction of junctions expressed in normal
tissue **before** the expensive translate/bind steps — a principled feature-space reduction
(not a crude "never-in-normal" binary filter).

## Cost & shipping model
- BayesTS is a pyro/torch SVI joint model (CPU, cross-platform). One-time cost over the full
  Ensembl91 GTEx control ≈ 2.6 h (epoch 300) – 9 h (epoch 1000+); measured 3.75–13.5 ms/junction.
- **Default Ensembl91**: ship the precomputed `.snaf_stats.tsv.gz` (~tens of MB) so users never run BayesTS.
- **New control**: compute on first use and cache next to the h5ad; reuse thereafter.

## Status / TODO
- `control_stats.py` (`precompute_control_stats`, `load_control_stats`, `default_stats_path`) + the
  `snaf-precompute-control` CLI + the runtime `--max_bayests_percentile` filter are implemented.
- Pending: run the one-time BayesTS precompute for the shipped Ensembl91 table, and speed up BayesTS
  (vectorize SVI / reduce epochs at equal sigma) before distribution.
