# Daedalus

UniProt-only isoform-aware classifier that predicts whether a transmembrane
alternative protein isoform reaches the plasma membrane. The deployable
model is a stacked meta-learner over twelve base classifiers (six
gradient-boosted-tree / linear baselines and six residual-MLP autoencoder
variants) trained on a curated UniProt benchmark of 2,839 plasma-membrane
positive isoform pairs and 183 stringently filtered non-surface negative
pairs across 1,396 genes. Under gene-grouped five-fold cross-validation
the stacked model achieves **95.1% sensitivity, 71.0% specificity, 98.1%
positive predictive value, and AUROC 0.89** at its operating threshold.

See [training/checkpoints/BENCHMARK_REPORT.md](training/checkpoints/BENCHMARK_REPORT.md)
for the full algorithm comparison and confusion matrices.

## Pipeline overview

```
data_ingest/   →   benchmark/   →   training/   →   model artifacts
  (UniProt          (build the         (kfold +         (per-fold OOF
   ingest +          per-pair           stacking +        probabilities,
   features +        feature             consensus)        stacked model
   embeddings)       matrix)                              outputs)
```

Each stage is a directory with its own `scripts/`, `data/` (or
`checkpoints/`), and `README.md` documenting inputs, outputs, and
reproduction commands.

## Directory layout

| Directory | Purpose |
|---|---|
| `data_ingest/` | UniProt fetch, isoform-sequence reconstruction, UniProt-feature counts, trafficking-biology descriptors, ESM2 embeddings, DeepLoc-style localization classifier |
| `benchmark/` | Assemble the (canonical, alternative)-isoform-pair feature parquet from data_ingest outputs; apply strict-canonical-PM label filter |
| `models/` | Model class definitions (residual-MLP autoencoder variants, DeepImmuno-style CNNs, family multitask network, training objectives) |
| `training/` | k-fold runner, stacking meta-learner, consensus ensemble, OOF prediction artifacts |
| `archive/` | Compressed tarballs of obsolete development material (gencode pipeline, proxy tasks, biochem benchmarks, intermediate kfold rounds). Excluded from git via `.gitignore`. |

## Reproduction

End-to-end from a fresh clone (CPU-only, ~30 minutes excluding optional
ESM2-650M extraction):

```bash
source .venv/bin/activate
cd components/Daedalus

# 1. Data ingest — one-time UniProt fetch and feature extraction
python data_ingest/scripts/download_uniprot.py
python data_ingest/scripts/build_isoform_proteins.py
python data_ingest/scripts/build_isoform_features.py
python data_ingest/scripts/build_trafficking_features.py
python data_ingest/scripts/build_isoform_esm_embeddings.py
python data_ingest/scripts/build_canonical_esm_embeddings.py
python data_ingest/scripts/build_localization_classifier.py

# 2. Benchmark assembly with the production strict-canonical-PM filter
python benchmark/scripts/build_benchmark.py --strict-canonical-pm

# 3. k-fold training (classical and residual-MLP base models)
python -m Daedalus.training.scripts.run_kfold_benchmark \
    --models phase_d_logistic,phase_c_elasticnet,phase_c_adaboost,phase_d_hist_gradient_boosting,phase_d_xgboost,phase_c_random_forest \
    --tune-threshold --sensitivity-floor 0.85 \
    --out-dir training/checkpoints/kfold_classical

python -m Daedalus.training.scripts.run_kfold_benchmark \
    --models phase_c_residual_mlp_autoencoder,residual_mlp_ae_wide,residual_mlp_ae_deep,residual_mlp_ae_highdrop,residual_mlp_ae_strongrecon,residual_mlp_ae_narrow \
    --tune-threshold --sensitivity-floor 0.85 --max-epochs 30 --patience 6 \
    --out-dir training/checkpoints/kfold_residual_mlp

# 4. Stacking + consensus on combined OOF
python -m Daedalus.training.scripts.stacking_meta_learner
python -m Daedalus.training.scripts.ensemble_consensus \
    --out training/checkpoints/ensemble_consensus.tsv

# 5. Inspect performance
cat training/checkpoints/BENCHMARK_REPORT.md
```

## Notes on legacy identifiers in source code

Model name strings like `phase_d_xgboost` and `phase_c_adaboost` are
**identifier labels**, not import paths. They serve as keys in the kfold
runner and as column-prefix names in the OOF prediction parquets
(`prob_phase_d_xgboost`, etc.). Renaming them would break OOF-checkpoint
backward compatibility, so they are retained verbatim. A handful of
historical `phase_a/...` and `phase_d/...` references remain inside
docstrings of inherited files; they are accurate-at-the-time comments
and have been preserved.

## Status

- **Production**: yes. Outputs at the operating threshold are used
  upstream as the localization gate for ExNeoEpitope target
  prioritization.
- **Hardware**: CPU-only.
- **Active development**: AlphaFold3 structural features, expanded
  literature-mined negative-control reference (~300 isoforms), and
  prospective MOLM13 surfaceome validation are in progress under
  the parent grant proposal.
