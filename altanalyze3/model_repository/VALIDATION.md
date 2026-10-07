# Validation on 2026-10-06

- Registry, five synthetic RNA model loaders/prediction values, H5AD round trips,
  resource provenance, immutable IDs, proposal checks, static selectors and standalone
  viewer metadata: covered by `components/model_registry/tests/test_registry.py`.
- Combined registry, existing dense modality builder and scALABLE-discover app checks:
  **32 passed**. The synthetic ADT fixture emits a harmless sklearn feature-name warning.
- Existing fastComm/scALABLE-discover regressions: **14 passed, 4 failed**. All four failures
  also reproduce against the unchanged HEAD fastComm API. AnnData 0.10.9 lacks the
  `anndata.io` imports used by those existing tests; auto species detection consequently
  fails in one of them. No analytical method or installed dependency was changed to
  accommodate these unrelated failures.
- Python compilation, JavaScript syntax, proposal example validation and whitespace checks
  passed.
- Snapshot records 11 observed model/resource versions and 27 default selections.
  History records 11 committed default intervals. A separate working-tree lipid release
  selector changed during implementation and is inventoried without being altered or
  approved by the registry work.

No model was trained, no scientific dataset was analyzed, and no live website or public
model repository was deployed for these checks. An end-to-end production dataset run
has not been performed.


## Latest registry refresh: registry-2026.10.06.1

- 19 registry/versioning tests passed, including retention of historical artifacts,
  inference variants, release evidence and download locations.
- Combined registry, dense modality builder and scALABLE-discover checks: 32 passed,
  4 failed. The four existing count-conversion fixtures exercise scale/negative-value
  expectations that differ from the recently committed scALABLE behavior. The pipeline,
  imputed-scale implementation and fixture file are byte-identical to AltAnalyze3 HEAD
  b550c6b for this registry update. No numerical algorithm was changed to satisfy them.
- Both model bundle hashes match their existing release manifests. No model payload was
  executed for this verification and no biological data were rerun.
- Catalog checks verify 12 distinct saved artifact IDs (10 current, 2 historical),
  all 27 unique default contexts, 7 changed default selections, and preserved before/after
  catalog snapshots. Git-backed default history now contains 12 repository intervals.
- The release manifests contain recorded training rosters and authorizations. Their
  statements are reported as source evidence, not as a new scientific approval.


## Modality/version organization

- 23 registry/versioning tests passed, including readable result labels, version folder
  generation and rejection of conflicting version assignments.
- All 12 pre-existing saved-artifact identities and hashes, all 27 default contexts,
  and archived catalogs were retained. No model, sample roster or feature panel changed.
- Models are grouped under 5 modalities with v1.0/v1.1 directories and tissue/species
  variants. Seven matching modality/version releases distribute the registered assets.
- Proposal metadata now requires explicit modality, readable version and variant.
