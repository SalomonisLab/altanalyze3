# fastComm

`fastComm` is an early receptor-ligand communication component for scoring
sender-state to receiver-state interactions from one transcriptome capture.

The design keeps ligand-receptor expression as a required signal, then raises
confidence when the receiver state shows a downstream transcriptional response
consistent with that ligand, pathway, or ligand-receptor axis.

## Current Status

This first scaffold includes:

- fast state pseudobulk construction from cell-by-gene matrices or `.h5ad`
- complex-aware ligand and receptor expression scoring
- optional receiver-response scoring from a ligand/pathway/LR response matrix
- a source registry for autonomous collection targets
- ignored working directories for downloaded data, trained bundles, caches, and benchmarks

The bundled `resources/seed_ligand_receptor.tsv` is only a tiny smoke-test
resource. Production use should build a versioned LR table from connectomeDB,
CellPhoneDB, CellChatDB, LIANA/OmniPath, ICELLNET, and CellTalkDB.

## Ignored Artifact Layout

Large files are intentionally excluded from git:

- `components/fastComm/training_data/`
- `components/fastComm/models/`
- `components/fastComm/benchmarks/`
- `components/fastComm/cache/`
- `components/fastComm/artifacts/`

The source registry lives at:

`components/fastComm/configs/source_registry.json`

Prototype evaluation notes and benchmark commands are tracked in:

`components/fastComm/EVALUATION.md`

## Score Command

From a cell-by-gene expression matrix and metadata:

```bash
python -m altanalyze3.components.fastComm.cli score \
  --expression expression.tsv \
  --metadata metadata.tsv \
  --state-key cell_state \
  --output /tmp/fastcomm_scores.tsv
```

From AnnData:

```bash
python -m altanalyze3.components.fastComm.cli score \
  --h5ad sample.h5ad \
  --state-key cell_state \
  --output /tmp/fastcomm_scores.tsv
```

With receiver-response evidence:

```bash
python -m altanalyze3.components.fastComm.cli score \
  --h5ad sample.h5ad \
  --state-key cell_state \
  --response-matrix ligand_response.tsv \
  --output /tmp/fastcomm_scores.tsv
```

The response matrix should have rows keyed by `lr_key`, ligand, or pathway and
columns as genes. Positive values indicate expected induction in the receiver;
negative values indicate expected repression.

## Human And Mouse Gene Symbols

`fastComm` matches required LR/response genes case-insensitively and renames
matched input columns to the canonical symbols used by the LR and response
tables. This allows mouse title-case symbols such as `Cxcl12`, `Cxcr4`,
`Tgfb1`, and `Tgfbr1` to work with the bundled canonical seed resources.

The run summary reports:

- `n_required_genes`
- `n_matched_required_genes`
- `n_missing_required_genes`
- `missing_required_genes`

For production mouse resources, a species-specific LR/response matrix is still
preferred when interactions or response targets are not one-to-one orthologs.

## Scoring Model

The current initial score is:

```text
fastcomm_score =
  0.55 * scaled_lr_expression_score +
  0.30 * receiver_response_score +
  0.15 * positive_state_promotion_score
```

`lr_expression_score` combines ligand expression in the sender, receptor complex
expression in the receiver, expression specificity, complex completeness, and
database evidence weight.

`receiver_response_score` is the positive cosine agreement between the receiver
state delta and the expected ligand/pathway/LR response vector.

`state_promotion_score` keeps the signed cosine so antagonistic or mismatched
responses remain visible in the output.

`lr_empirical_p` and `fastcomm_empirical_p` are permutation-free upper-tail
rank statistics within the scored run. The companion percentile columns are
often easier to read for prioritization; they are empirical ranks, not
calibrated probabilities.

The CLI defaults filter weak raw evidence with:

- `--min-ligand-expr 0.01`
- `--min-receptor-expr 0.01`
- `--min-lr-expression-score 0.001`

These are deliberately permissive for normalized single-cell matrices. Raise
them for cleaner discovery reports.

## Benchmark A Dataset

Run split stability on one `.h5ad`, for example by donor:

```bash
python -m altanalyze3.components.fastComm.cli benchmark-h5ad \
  --h5ad sample.h5ad \
  --state-key "Level 2" \
  --split-key Donor \
  --output-dir components/fastComm/benchmarks/sample_level2_donor
```

Outputs:

- `full_scores.tsv`
- `split_<value>_scores.tsv`
- `split_stability.tsv`
- `benchmark_summary.json`

`split_stability.tsv` reports top-edge Jaccard and Spearman score correlation
against the full run.

## Build Learned Response Signatures

Convert long perturbation signatures into a scorer-ready matrix:

```bash
python -m altanalyze3.components.fastComm.cli build-response-matrix \
  --input cytokine_response_long.tsv \
  --output components/fastComm/models/cytosig_response_matrix.tsv \
  --manifest components/fastComm/models/cytosig_response_matrix.manifest.json \
  --id-col signature \
  --gene-col gene \
  --score-col score \
  --top-genes 200 \
  --min-abs-score 0.1
```

The output should be passed to `score --response-matrix`. Model outputs belong
under `components/fastComm/models/`, which is ignored by git.

## Development Roadmap

1. Add source collectors that materialize versioned LR, complex, signaling,
   TF-target, and perturbation-response tables under ignored paths.
2. Build a normalized evidence graph with signed signaling and provenance.
3. Train ligand-to-response priors from perturbation datasets and distill them
   into sparse matrices for fast inference.
4. Add calibration models for probability-like confidence scores.
5. Add benchmarks against perturbation recovery, spatial adjacency, receptor
   protein support, and runtime at 10k, 100k, and 1M cells.
