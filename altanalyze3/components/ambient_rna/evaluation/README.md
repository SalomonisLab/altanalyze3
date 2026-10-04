# Ambient RNA evaluation

Benchmark and simulation drivers are kept here, separately from production
correction code. Results remain in the sibling `../benchmarking/` directory.

Run modules from the repository root containing the `altanalyze3` package:

```bash
.venv/bin/python -m altanalyze3.components.ambient_rna.evaluation.benchmark_nd20_167
.venv/bin/python -m altanalyze3.components.ambient_rna.evaluation.simulate_nd20_167_release
.venv/bin/python -m altanalyze3.components.ambient_rna.evaluation.summarize_nd20_167_release
.venv/bin/python -m altanalyze3.components.ambient_rna.evaluation.write_nd20_167_report
```

`simulate_nd20_167_release.py` implements the corrected-baseline RNA release model:
ambient gene probabilities are proportional to population mean transcript counts
weighted by the number of contributing cells. Its loading models and fractions
are saved before execution. `STRUCTURED_CASE` optionally selects a single case,
for example `heterogeneous_rho0.20_seed101`. The summarizer requires all cases.

`simulation_engine.py` executes the release simulations; `simulation_analysis.py`
provides shared cellHarmony tests, embeddings and sparse count sampling.
`simulation_summary.py` exports the comparison tables and figures.
`assess_simulation_restoration.py` calculates cell-type gene proportions and
introduced-DEG resolution.

`simulation_soupx.R` runs the unmodified official SoupX R functions from the
pinned source checkout recorded with the results. Python correction uses the
production scALABLE/cellHarmony interface. Reports document dependency versions,
input provenance, statistical thresholds, completed conditions and failures.

Production `ambient_subtract.py` and `soupx_correct.py` remain in the parent
directory. Moving these evaluation drivers does not change correction defaults.
