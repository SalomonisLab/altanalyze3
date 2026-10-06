# Enforced analysis requirements

The LungMAP scientific analysis is **blocked**. Safeguard implementation and tests are authorized; they do not authorize reconstruction, retraining, inference or a differential rerun. The outstanding D071 question is recorded in `decision_state.json`.

The user's full instruction is in the global `/Users/saljh8/.codex/AGENTS.md`. The reusable skill is `/Users/saljh8/.codex/skills/preserve-analysis-integrity/SKILL.md`, also discoverable through the symlink `/Users/saljh8/.agents/skills/preserve-analysis-integrity`. Repository `AGENTS.md` requires it for covered work. Global and project guidance discovery follows the [official AGENTS.md documentation](https://learn.chatgpt.com/docs/agent-configuration/agents-md); the symlink supports the [documented skills discovery locations](https://learn.chatgpt.com/docs/build-skills).

## What the checks enforce

- The independent baseline contract contains every ordered production identifier: 202 lipid outputs, 50 training profiles and 1,303 RNA genes. Counts alone are insufficient. Missing, extra, duplicated, substituted and reordered identifiers fail.
- The production bundle, trainer, API, delivered trainer and original RNA source are SHA-256 pinned. The guard also pins the contract itself so redefining its roster cannot silently make a reduced run pass.
- The method record covers the estimator, fitted scalers, per-lipid gene selection, ElasticNetCV hyperparameters, CV/seed, original RNA handling, training scope, inference and differential policy. A candidate label does not authorize differences.
- Every required sample and lipid needs a verified source record and source identifier. Sources must exist and retain their reviewed hashes. Target provenance must bind to the exact target file, enumerate imputation, describe transformations and have no unresolved dependencies. A legacy reconstruction cannot be labeled as corrected data merely because it reproduces saved model statistics.
- The reviewed run manifest must match the contract and exact input hashes. An approval must reference that manifest's SHA-256, the actual user message, and its conversation reference. Pending questions block execution even if a document is labeled approved. There is no force/skip or approve CLI.
- `fit_verified_candidate` checks actual tables before estimator fitting and uses the pinned baseline trainer. It does not automatically fill target values or write/deploy a bundle.
- The ten older narrowed/alternative candidate drivers in `protected_entrypoints.json` reject both CLI and callable orchestration before scientific work or output writes. Further source reconstruction diagnostics are separately blocked while the missing-source question is pending.

## Before any future run

Ask the user for the outstanding source, then wait. After their answer, perform only the investigation they authorize. Keep all expected identities in the audit; do not reinterpret poor QC or absence from an intermediate table as permission to drop them.

Prepare a concrete full-panel proposal with the baseline comparison, exact sources, mappings and transformations. Discuss every material difference. Only after actual user approval may its literal answer and manifest hash be recorded; never create an approval because implementation or tests finished. A resolved question needs its answer and supporting evidence. If the user explicitly changes scope, document and review that change rather than silently changing the contract.

An approved future manifest has schema version 1, `baseline_contract_sha256`, the reviewed `method`, and `inputs` containing `RNA`, `corrected_targets`, and `provenance` file records (`path` and `sha256`). Corrected targets additionally require `data_role: corrected_training_targets`. The provenance JSON binds `target_sha256`, lists all `sample_sources` and `feature_sources` in baseline order, includes verified status, source records, source identifiers and evidence, and explicitly lists `unresolved` and `imputed_entries`. See `candidate_integrity.py` for the checked fields. No approved manifest exists now.

## Verification

```bash
/tmp/rna2lipid_elasticnet161/bin/python components/rna2lipid/validation/test_candidate_integrity.py
python3 components/rna2lipid/candidate_integrity.py --status
```

The status command intentionally exits 2 while blocked. Test approval records are synthetic, confined to temporary directories, and never modify the real decision state. Tests exercise a valid synthetic preflight and invalid coverage/method/provenance/approval cases, actual reduced-candidate rejection, CLI refusal without writes, and callable driver refusal. No real model is fitted by the safeguard tests.

## Enforcement boundary

These are tested local workflow protections and persistent instructions, not a repair or audit of OpenAI's internal harness. They do not authenticate conversation authorship, prove scientific truth from metadata, or prevent a process with unrestricted filesystem access from modifying code or bypassing a library wrapper. The skill explicitly prohibits those bypasses. The unchanged general training library remains usable for unrelated authorized work; using it to evade this candidate's guards violates the user's instruction. Structural skill validation does not prove future model compliance. Report these limits rather than claiming an unbreakable guarantee.

Production models, prediction matrices and statistical results were not regenerated while implementing these safeguards.
