# Archive manifest

Archived material from earlier development phases of Daedalus. None of these
files are required by the production pipeline. Tarballs are excluded from git
via `.gitignore`. To inspect contents without extracting, use
`tar -tzf <archive>.tar.gz | head`.

## gencode_pipeline.tar.gz (~42 KB)

Original gencode-keyed feature pipeline. Replaced by the UniProt-only flow.

Includes: `phase_a/scripts/build_gencode_*`, `build_isoform_pair_candidates.py`,
`build_baseline_feature_matrix.py`, `build_appris_*`, `build_clinvar_*`,
`build_hpa_*`, `build_mane_*`, `build_priority_*`, `build_seed_*`,
`build_tm_negative_*`, `build_transcript_supervision_*`,
`build_uniprot_functional_priors.py`, `build_uniprot_gencode_map.py`,
`build_uniprot_isoform_tables.py` (single-isoform-per-entry, superseded by
`build_uniprot_isoform_proteins.py`), `build_uniprot_region_priors.py`,
`build_biogrid_gene_interactions.py`, `download_appris_references.py`,
`download_gencode_references.py`, `download_resources.py`, `init_phase_a.py`,
`resolve_uniprot_to_gencode_isoforms.py`, `run_seed_baselines.py`,
`run_sklearn_seed_baselines.py`, `run_torch_seed_baselines.py`,
`train_structural_task.py`, plus the `phase_a/models/` package and the
`phase_a/STAGE2.md`, `STATUS.md`, `FEASIBILITY.md`, and `README.md` notes.

Reason for retirement: the entire feature pipeline was rebuilt as
UniProt-only because gencode-keyed alternative_is_membrane / alternative_tm_feature_count
fields turned out to be parent-entry-level rather than isoform-resolved, which
caused systematic label-noise.

## proxy_tasks.tar.gz (~50 MB)

Entire phase_b: proxy supervised tasks built from gencode pair_candidates,
reference-delta multitask training, query inference. Replaced by direct
UniProt-isoform binary classification.

## biochem_benchmarks.tar.gz (~536 MB)

Phase_c biochem-delta-matrix and DeepImmuno-style multitask benchmarks.
Includes `phase_c/scripts/`, `phase_c/data/`, `phase_c/checkpoints/`,
`phase_c/models/class_specific_biochem_autoencoder.py`. The other model
classes (`deepimmuno_style.py`, `residual_mlp_variants.py`) were retained
in `models/` because they are still imported by the production kfold runner.

Reason for retirement: superseded by `training/run_kfold_benchmark.py` on
the UniProt-only benchmark.

## intermediate_rounds.tar.gz (~14 MB)

All Round 0 through Round 11 kfold/stacking checkpoint directories from the
hyperparameter and label-set iteration history:

- `kfold_classical{,_v2,_v3,_v4,_v5}/`, `kfold_residual_mlp{,_v2,_v3,_v4,_v5}/`,
  `stacking{,_v2,_v3,_v4,_v5}/` (Rounds 0–4: gencode-keyed)
- `uniprot_kfold_classical{,_clean,_final,_r1,_r3,_tier1,_tuned}/` and
  `uniprot_kfold_residual{_mlp,_mlp_tuned,_clean,_esm2,_esm2t30,_final,_r1,_r2,_tier1,_traffic}/`
  (Rounds 5–11: UniProt iteration)
- `uniprot_stacking{,_clean,_final,_v3,_tier1}/`
- `ensemble_consensus_binary*.tsv`, `uniprot_ensemble_*.tsv` (all variants
  except the production strict-canonical artifact)
- Older `phase_d_*` checkpoints, `tm_isoform_model_comparison_*` artifacts

The single deployable result (Round 10, strict-canonical Tier-1.1) is
preserved live under `training/checkpoints/`.

## intermediate_data.tar.gz (~1.2 GB)

All gencode-derived TSVs, the v1 ESM2-t12 embeddings (superseded by t30),
the entire `phase_d/data/processed/phase_d_*` and `phase_d/data/processed/tm_isoform_*`
tree (superseded by `benchmark/data/uniprot_isoform_benchmark_*`), plus
the multi-task instance and objective-score files from earlier rounds.

## To restore an archive

```bash
cd <repo-root>/components/Daedalus
tar -xzf archive/<archive>.tar.gz
```

Tarballs preserve the original `phase_a/`, `phase_b/`, `phase_c/`, `phase_d/`
relative paths.
