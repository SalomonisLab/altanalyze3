# Landmark UMAP comparison — October 6, 2026

## Authorization and current status

Nathan explicitly authorized benchmarking landmark UMAP against the current full fit
(“Yes you can benchmark that”). This authorizes an isolated comparison with all cells
receiving coordinates. It does not authorize changing production defaults.

The real-job comparison is **pending an original-source answer**. Job
`d5fee6ea361e4ee2bdfce51b9c6ee312` used 149,745 final cells and 1,307 final post-SVM
MarkerFinder genes, with correlation distance, 50 neighbors, min_dist=0.75 and seed=0.
The 700.3-second UMAP phase is recorded in its original log. Its `minimal_outputs=true`
configuration skipped the original ordered marker/UMAP-feature export. The persisted
canonical heatmap markers come from a separate MarkerFinder pass; they are not a
verified substitute for the original UMAP feature list.

A source question was presented on October 6: provide the original
`icgs3_markers.tsv`/`icgs3_umap_features.tsv`, or authorize regeneration using the
original post-SVM MarkerFinder procedure. No source reconstruction or real-job fit
may proceed while that answer is pending. An engineering smoke test uses a fully
specified synthetic fixture and makes no claim about this job's performance.

## Prepared benchmark

`benchmark_landmark_umap.py` accepts a verified prepared matrix:

- `X.npy`: the exact C-contiguous float32 input the baseline UMAP consumes.
- `cells.txt` and `features.txt`: complete ordered identities, one per line.
- `states.npy`: original cell-state labels aligned to those cell IDs.
- `manifest.json`: hashes for every preceding file, pinned UMAP parameters,
  recorded baseline source/method verification, and explicit benchmark authorization.

The manifest may only record baseline verification after the missing source has been
resolved and original expression preprocessing has been checked. No bypass flag or
feature intersection is provided. Hash mismatches, duplicate identities, roster
mismatches and nonfinite inputs reject fitting.

Each fit runs in its own process. The candidate retains every feature and returns
coordinates for every cell, fitting only selected landmarks and transforming the
remaining cells in bounded batches. Fitting, transformation, verification, peak RSS,
constructor parameters and library versions are reported separately. The default
sampling candidate now reserves at least 200 cells per existing state, every cell in
smaller states, and allocates remaining places proportionally within larger states.
The budget expands when the minimum coverage requires more places. Historical
random-selection comparisons remain explicitly selectable for reproduction only. This changes
which cells are fitted, and is reported as a methodological difference. It does not
remove cells or change upstream clustering, MarkerFinder, imputation or differentials.

Comparison measures input-neighbor recall, coordinate-neighbor retention against the
full fit, and state purity for each state. Evaluation samples up to 50 cells per
state, including every cell in smaller states. State purity evaluates separation of
existing cluster labels; it is not an independent biological truth or proof that a
layout is accurate. The saved full-fit neighbor graph provides a common reference;
small engineering fixtures use brute-force input neighbors when UMAP did not retain
that graph.

## Commands after baseline verification

From the repository root with the original job's Python/library environment:

```bash
OPENBLAS_NUM_THREADS=2 NUMBA_NUM_THREADS=4 python \
  altanalyze3/components/clustering/benchmarking/benchmark_landmark_umap.py \
  fit /path/to/verified-input /path/to/full --mode full

OPENBLAS_NUM_THREADS=2 NUMBA_NUM_THREADS=4 python \
  altanalyze3/components/clustering/benchmarking/benchmark_landmark_umap.py \
  fit /path/to/verified-input /path/to/landmark-30k \
  --mode landmark --landmarks 30000 --strategy state-stratified

python altanalyze3/components/clustering/benchmarking/benchmark_landmark_umap.py \
  compare /path/to/verified-input /path/to/full /path/to/landmark-30k \
  /path/to/comparison.json
```

The proposed real comparison includes 10,000 and 30,000 landmarks, with all cells
transformed and full feature coverage verified. Automatic fitting and transformation
epochs depend on input/batch size in UMAP 0.5.7; these defaults are retained and must be
reported, rather than described as identical fitting work. Full preprocessing and
pipeline timing are outside this embedding-only benchmark.

## Safeguard verification

The runner completed a 1,000-cell, 32-feature synthetic full-fit and 300-landmark
fit/transform smoke test, including all-cell coordinate and per-state neighborhood
checks. This fixture uses 30 fitting epochs to exercise the harness; it is not the
real-job benchmark. Raw smoke measurements are saved in
`landmark_umap_smoke_results.json`; JIT/setup costs dominate this fixture and its
timings must not be extrapolated to the real job.

Future minimal-export ICGS3 marker-feature UMAP runs now persist their small ordered
`UMAPs/icgs3_umap_features.tsv` provenance file. A regression test checks exact order
under minimal exports. This does not recover the missing panel for the existing job
or authorize reconstruction while its source question is pending.

## Subsequent Discover option authorization

Nathan subsequently requested an accelerated UMAP option for scALABLE-discover, while
explicitly retaining ICGS3's default. Discover now selects a state-stratified 30,000-cell
landmark fit by default, with a Full fit selector. ICGS3 defaults to full fitting.
`umap_fit.py` loads the fitting matrix and remaining cells in bounded blocks and preserves
the baseline feature panel. Normalized expression is passed with `input_normalized=True`.

Verification uses a new run of the 65,662-cell upload from failed job
`3e0bce435c54416d831997a3f74d1772`. Its feature panel is captured during the original
pipeline invocation, so this work is independent of the unresolved source question for
the older 149,745-cell job. That older comparison remains paused.

## Verified Discover comparison

Saved evidence: `discover_umap_results.json`. On the new run's complete ordered 60,164
clustered cells and 1,273 original UMAP features, fresh processes measured:

| Embedding | Fit | Transform | Total | Peak RSS |
| --- | ---: | ---: | ---: | ---: |
| Full fit | 252.72 s | — | 252.72 s | 3.09 GiB |
| 30,000 landmarks, maximum 50,000-cell transform blocks | 122.08 s | 60.89 s | 182.97 s | 2.65 GiB |

The production landmark helper is used by the paired candidate. Its selected barcodes
exactly match the completed pipeline's saved landmark roster, and every state's quota
passes. Both outputs contain finite coordinates for every cell. The isolated embedding is
28% faster and uses 14% less peak RSS; these are not whole-pipeline reductions.

Landmark fitting changes local geometry. On 3,255 state-stratified evaluation cells, recall
of 15 input neighbors was 0.1714 full versus 0.1439 landmark. Mean per-state neighborhood
purity was 0.7166 versus 0.7013; C4 (122 cells) had the largest purity drop, 0.248. The
interface offers full fitting and explicitly notes that layouts can differ. These cluster
separation checks are not independent biological validation.

The first complete pipeline, job `6b5c827ff6234dd5b53890753e339e74`, finished in 379.45 s
at 13.28 GiB peak RSS. It used the initially tested 10,000-cell transform blocks, which took
74.55 s; the final default is 50,000, verified above. Smaller blocks trigger 100 transform
epochs in UMAP 0.5.7, versus 30 above 10,000 cells. Clustering, genes and all assigned cells
were preserved. The serving bundle built successfully, and browser checks covered UMAP,
RNA expression, violin, filters and both cell-state layers.

That upload's separate cell-communication step failed to match ligand/receptor genes, and
its predicted cell-state labels are UNK. Diagnosis of identifier mapping is paused pending
the user's answer to the required source/diagnostic question. It was not repaired by
guessing mappings or changing feature identities.
