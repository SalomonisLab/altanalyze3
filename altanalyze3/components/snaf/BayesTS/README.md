# Bundled BayesTS references

Precomputed per-junction BayesTS tumor-specificity scores that ship with SNAF. **SNAF never
recomputes BayesTS.** It reads these tables.

BayesTS is a pyro SVI model costing hours over millions of junctions, and its answer depends
only on the **reference** cohort — not on the tumor cohort being analysed. So the finished
scores are stored once and reused by every run.

## Files

| File | Reference | Junctions | Source cohort |
|---|---|---|---|
| `GTEx_BayesTS.tsv.gz` | GTEx | 2,476,734 | 2,629 samples across 51 tissues |
| `TabulaSapiens_BayesTS.tsv.gz` | Tabula Sapiens | 2,852,713 | 762 pseudobulk samples across 112 cell types |

Each is a gzipped TSV, ~50 MB, in the shape `control_stats.load_control_stats` already reads:

```
uid <tab> bayests_sigma <tab> bayests_percentile
```

* UIDs carry **no** `=chr...` coordinate suffix — they are the key SNAF sifts on.
* `bayests_sigma` ∈ [0,1]; **lower = more tumor-specific** (less expressed across normal).
* `bayests_percentile` is sigma's rank within that reference, in [0,1].
* Values are stored to 6 significant digits, far finer than any threshold decision.

## Provenance

Built by `build_bundled_bayests.py` from:

* **GTEx** — `GTEx_junction_counts.snaf_stats.tsv.gz`, produced 2026-07-06 by
  `altanalyze3 snaf-precompute-control`.
* **Tabula Sapiens** — `Run01b_full_results.zip`, identity confirmed by a 100.0% junction
  match against `counts.TabulaSapiens.h5ad` and a modal `X_mean` of 1/112. See that
  archive's README.

Rebuild with `python build_bundled_bayests.py` (the source paths are recorded in its
docstring). Neither source is modified.

## How SNAF uses them

`gtex_configuration()` attaches `bayests_sigma` and `bayests_percentile` to the control
`AnnData` whenever they are not already present, so `--max_bayests_percentile` works in
**every** configuration.

This matters: the summary-backed precompute path is skipped whenever `--add_control` is
given (`gtex.py`, "Skipped when add_control is given"), which previously left the control
without a `bayests_percentile` column and made `--max_bayests_percentile` a **silent
no-op**. It is no longer silent — if no BayesTS score can be attached, SNAF now prints that
the filter is not being applied.

The scores are also written into every frequency table by
`snaf.add_bayests_frequency_table()`, so `bayests_sigma` and `bayests_percentile` appear in
the standard SNAF-T and SNAF-B outputs alongside `tumor_specificity_mean` and
`tumor_specificity_mle`.

### Combining the two references

`load_bundled_bayests()` uses both by default and combines them by taking the **maximum
percentile**: a junction must look tumor-specific in *every* reference to survive, so broad
expression in either normal cohort disqualifies it. That direction is conservative — it can
reject candidates, never invent them.

### Junctions with no score

A junction absent from a reference gets `NaN`, meaning *never observed in that normal
cohort* — i.e. maximally tumor-specific. It is **kept**, and must never be read as sigma 0.
Measured on the 3,325 pediatric-AML SNAF-B neojunctions: 2,449 scored by GTEx, 1,256 by
Tabula Sapiens, 2,578 by their union, 747 absent from both.

## API

```python
from altanalyze3.components.snaf.control_stats import (
    bundled_bayests_path, load_bundled_bayests, DEFAULT_BAYESTS_REFERENCES)

bundled_bayests_path('tabulasapiens')            # -> absolute path, or None if absent
load_bundled_bayests(uids=my_junction_uids)      # -> DataFrame, both references combined
load_bundled_bayests(references=('gtex',))       # -> one reference only
```

## Adding another reference

1. Produce `uid <tab> bayests_sigma <tab> bayests_percentile`, gzipped, UIDs without the
   `=` suffix.
2. Drop it in this directory.
3. Register it in `BUNDLED_BAYESTS` in `control_stats.py`, and add it to
   `DEFAULT_BAYESTS_REFERENCES` if it should apply by default.

## Size note

These two files add ~99 MB to the repository. They are data, not code, and change only when
a reference cohort is rebuilt. If that is unwanted in git, move them to a shared reference
directory and point `BAYESTS_DIR` at it — the loader resolves the path at import time and
degrades with a warning when a table is absent.
