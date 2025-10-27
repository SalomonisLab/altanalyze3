## GO-Elite: first-time human setup and analysis

**Assumed working directory:** `/Users/saljh8/Dropbox/altanalyze3/analyses/goelite`

These instructions use the modern GO-Elite component bundled with AltAnalyze3 (`altanalyze3/components/goelite`).  
Prerequisites: Python 3.11+, `pytest`, and the standard AltAnalyze3 dependencies.

---

### 1. Create a virtual environment (recommended)

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e /Users/saljh8/Documents/GitHub/altanalyze3
```

*(The editable install makes AltAnalyze3 modules available while you work from this analysis folder.)*

---

### 2. Inspect the GO-Elite CLI help

```bash
python3 -m altanalyze3.components.goelite.main --help
```

Key options:

* `--species`: `human` or `mouse` (enables cached downloads)
* `--cache-dir`: optional custom cache location (default: `~/.altanalyze3/goelite`)
* `--force-refresh`: rebuild cached resources even if already present
* `--query` / `--background`: newline-delimited gene lists
* Output files are written to `--outdir` (TSV + JSON)

---

### 3. Build (or refresh) the human GO cache

The first run downloads `go-basic.obo` and `goa_human.gaf`, attaches gene sets, and stores everything in the cache.

```bash
python3 -m altanalyze3.components.goelite.main \
    --species human \
    --force-refresh \
    --query /path/to/query_genes.txt \
    --background /path/to/background_genes.txt \
    --outdir ./results/human_test
```

*Use your actual gene lists. The `--force-refresh` flag is optional after the initial build.*

After this step, the cache (default `~/.altanalyze3/goelite`) contains:

```
human/<version>/downloads/   # raw OBO + GAF used to build the cache
human/<version>/data/        # go_tree.json and term_gene parquet/TSV
human/latest.json            # points to the most recent version
human/versions.json          # list of cached builds
```

---

### 4. Run analyses using the cached resources

Once the cache exists you can omit `--force-refresh`.  Use any query/background lists.

```bash
python3 -m altanalyze3.components.goelite.main \
    --species human \
    --query /path/to/my_query_genes.txt \
    --background /path/to/my_background_genes.txt \
    --outdir ./results/human_run1
```

Outputs for each run:

* `goelite_results.tsv` – full table with z-score, p-value, FDR, overlap, selection flag, etc.
* `goelite_results.json` – same content as records.

---

### 5. Optional: run with explicit GO files (no caching)

If you want to bypass caching and point to specific GO versions:

```bash
python3 -m altanalyze3.components.goelite.main \
    --obo /path/to/go-basic.obo \
    --gaf /path/to/goa_human.gaf \
    --query /path/to/query_genes.txt \
    --background /path/to/background_genes.txt \
    --outdir ./results/custom_go
```

---

### 6. (Optional) Verify tests

To confirm the component operates as expected:

```bash
python3 -m pytest \
  /Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/tests/test_goelite_basic.py
```

This exercises prioritisation logic, enrichment scoring, and cache round trips.

---

### 7. Subsequent analyses

1. Activate the environment (`source .venv/bin/activate`).
2. Reuse `--species human` runs (fast, since cache is already built).
3. Adjust enrichment settings with CLI flags if needed (`--min-term-size`, `--max-fdr`, etc).

Keep per-project runs organised under subdirectories of this folder (e.g.
`results/human_runN`) to make it easy to compare outputs across analyses.
