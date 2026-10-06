# Model and dataset integrity

These requirements are explicit user instructions.

- Never exclude model features or dataset samples because an initial name, identifier, annotation, or upstream mapping does not match. Treat the mismatch as an unresolved failure. Immediately ask the user for the source or mapping and wait before further dependent investigation or reconstruction. Do not assume the records are unavailable.
- Preserve the full instructed baseline feature panel and sample roster. Reconcile every feature and sample against the baseline before training or reporting a replacement or comparison candidate.
- Always verify the proposed training and inference algorithms against the instructed baseline model, including preprocessing, feature selection, hyperparameters, training samples, target panel, transformations, and missing-input handling. Report any differences explicitly.
- Do not proceed with a major algorithm change, exclusion of features or samples, or a reduced-coverage substitute until the discrepancy, evidence, and proposed change have been fully discussed with the user and the user explicitly authorizes it. General authorization to complete an analysis does not authorize these changes.
- If baseline coverage cannot yet be reconciled, mark the requested comparison as incomplete. Do not present a reduced-panel or reduced-sample run as satisfying it. Diagnostic investigation requires explicit user authorization when the missing-source question is pending. Safeguard implementation does not authorize resuming the scientific analysis.
- For the current rna2lipid comparison, the required baseline is the deployed lipid-wise ElasticNetCV procedure, all 202 production lipid outputs, all 50 original training profiles, and all 1,303 production RNA input genes. The 47-lipid/45-profile candidate does not satisfy this requirement.

## Required integrity skill and execution gate

Read `/Users/saljh8/.codex/skills/preserve-analysis-integrity/SKILL.md` for covered work. The user's verbatim prohibition is also persisted in `/Users/saljh8/.codex/AGENTS.md`.

The LungMAP analysis remains stopped pending source answers and review. `components/rna2lipid/integrity/decision_state.json` is authoritative for pending decisions; do not mark it approved without an actual user response. Existing reduced/alternative candidate drivers must reject execution. Future candidate runs must use `components/rna2lipid/candidate_integrity.py`, retain every baseline identity, and pass method/source/authorization checks before fitting or writing model outputs. Do not bypass these checks through another script, change the contract to fit a reduced candidate, or add a skip/force option.

## Authorized scope change — 2026-10-05

The user explicitly instructed removal of D071 because it is a 60-year-old donor. The candidate must therefore contain the other 45 original training profiles, all 202 lipids and all 1,303 RNA inputs, using the deployed ElasticNetCV procedure. The original production baseline remains 50 profiles for provenance. D071 source provenance is no longer a prerequisite for this candidate; its exclusion is authorized, not inferred from a failed match. Other missing-source and methodological decisions still require discussion. This decision is recorded in `components/rna2lipid/integrity/decision_state.json`.
