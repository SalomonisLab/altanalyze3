# AltAnalyze3 / scALABLE model repository

This repository records model discovery, community proposals, and reviewed default history. It contains metadata, not executable community payloads.
The runtime catalog is packaged in `altanalyze3.components.model_registry/catalog.json`.
Repository: https://github.com/SalomonisLab/scalable-models (private).

## IDs recorded in results

- `model_version_id`: component plus SHA-256 of the named artifact hashes. Renaming or
  relocating an artifact preserves its identity; replacing its bytes changes the ID.
- `inference_version_id`: SHA-256 of the explicitly listed inference source files.
- `analysis_model_version_id`: SHA-256 of both artifacts and inference source hashes.
- `runtime_versions`: package versions, recorded separately from model identity.

Model IDs identify saved artifacts, including learned parameters and embedded metadata.
They do not certify training provenance. Inference IDs cover the explicit file inventory
in `registry.py`, not arbitrary dynamically imported plugins. Custom model implementations
must extend that inventory when integrated. scALABLE's postprocessing and run settings
remain in its existing summaries/configuration; these IDs alone are not a complete run ID.

RNA prediction CLI exports always write `<output>.model_provenance.json`. fastComm does
likewise for score, state-pair and expression tables. scALABLE and scALABLE-discover expose
`model_provenance.json` as a downloadable artifact and store `model_versions` in job
metadata. Imputed H5AD prediction summaries also contain the IDs. The standalone
scALABLE viewer preserves supplied prediction IDs in its bundle metadata and manifest. Old results without
these records have unknown versions; do not assign the current default retrospectively.

## Current defaults and time frames

`catalog.json` records the actual bytes selected by the current working tree as observed
on **2026-10-06**, including API defaults and each scALABLE tissue/reference selection.
`observed-default` is an inventory status, not a new scientific approval. Unknown training
and deployment dates remain null. No trained model or selection was changed by this work.

`default_history.json` records committed AltAnalyze3 API and scALABLE default intervals, commit IDs and artifact Git
blob IDs. Git commit times are repository evidence, not training dates or deployment
start/end dates. Historical artifact identifiers are explicitly Git blob IDs, not current
SHA-256 IDs. An artifact absent from a commit remains null; history is not a claim that
it was deployed. Registry ingestion never loads historical model pickles.

Refresh from the parent directory of the `altanalyze3` Python package:

```bash
python -m altanalyze3.components.model_registry.snapshot
python -m altanalyze3.components.model_registry.history
```

The snapshot command is a maintenance operation: inspect its diff, then release matching
catalog copies together. It never changes a default. Default history should also record
actual deployment events when maintainers can supply that evidence.

## Methods and model cards

The inventory distinguishes lung lipid models, AML lipid models, human bone marrow,
human lung and mouse ADT models, GRN reference variants, and human/mouse fastComm resources.
Existing component documentation remains the method source:

| Family | Method documentation |
| --- | --- |
| Lung rna2lipid | `components/rna2lipid/README.md`, `api.py` (current lipid-wise ElasticNetCV; legacy multitask retained separately) |
| AML rna2lipid | `components/rna2lipid/aml/api.py`, `_impute.py` (per-target ridge) |
| rna2metabolite | `components/rna2metabolite/api.py`, `_impute.py` (per-target ridge) |
| rna2adt | `components/rna2adt/README.md`, `training.py`, `lung/model.py` and mouse documentation |
| rna2grn | `components/rna2grn/README.md`, `model.py` and reference-specific bundle metadata |
| fastComm | `components/fastComm/README.md`, `scoring.py`, upstream resource manifests |

These references document methods; the new catalog does not claim to have independently
reconciled historical training rosters. The registry does not approve LungMAP source/method decisions; its existing integrity
gates and recorded user decisions remain authoritative. Separate working-tree release
edits are inventoried as observed selectors, without certifying their scientific review.

## Community submissions

1. Fork the model repository and add `proposals/<unique-name>.json`, using
   `proposal.example.json` as a template. Host the payload and ordered input, output and
   training-sample identifier files in a stable HTTPS release/archive. Include hashes,
   license, training date, species/tissue, full method/code revision, sources, held-out
   evaluation and comparison against a named baseline version.
2. Run the metadata validator below. The validator checks fields and content IDs without
   fetching URLs or unpickling files. Uploaded payload verification and roster/method
   comparison are separate, mandatory maintainer review steps.
3. Open a pull request. Maintainers review provenance, artifact hashes, executable code,
   complete feature/sample identities, preprocessing, missing-input handling, validation
   and applicability. Document every difference relative to the baseline.
4. A merged proposal remains `proposed` until reviewed. Neither upload nor a passing
   metadata check changes inference or makes the proposal a default.

```bash
python -m altanalyze3.components.model_registry.proposals model_repository/proposals/*.json
```

Use proposals with a different panel or cohort as explicitly reviewed models; do not
silently replace a baseline with a reduced match. Existing candidate integrity gates and
pending decisions continue to apply.

## Promotion and rollback

Promotion requires an explicit maintainer approval record linking the PR, reviewers,
complete roster/method comparison, evaluation report and applicable tissue/species.
In one reviewed application PR: add an immutable catalog record, install verified
artifacts, update the appropriate API/reference default, and append the event to
`deployment_events.json` with application, context, old/new IDs, approved_by,
approval_url, effective_from, effective_until, release/version and commit. Close the
prior interval rather than editing its identity. Community PRs cannot write approvals.
Default selectors currently remain the existing API/reference configuration; registry
entries alone never drive automatic downloads or selection. Rollback uses a new event
pointing to a previous immutable ID. Never rewrite old run manifests or version IDs.

The CI template at `ci/validate-proposals.yml` validates proposal metadata. To activate
it, copy it to `.github/workflows/validate-proposals.yml` and push with credentials
that include GitHub’s `workflow` scope. The publishing token lacked that scope, so
GitHub Actions is not yet enabled; the standalone validator runs locally. No community model is automatically executed in CI.

## Integrating another model family

Add the family and its inference source inventory to `INFERENCE_FILES`, call
`describe_model(component, artifacts)` at load time, and retain that record on the loaded
object and every prediction summary. Use stable artifact roles such as `bundle`,
`ligand_receptor` and `response_matrix`; hash all resources that affect the result.
Use `write_provenance` for table exports and preserve summaries in H5AD `uns`.
Extend the snapshot selector inventory and add loader/export tests. Runtime integration
is required before a new community model family can become an application default.
