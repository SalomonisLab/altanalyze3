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
