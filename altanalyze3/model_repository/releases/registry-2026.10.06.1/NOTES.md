# Registry update — 2026-10-06

Release: `registry-2026.10.06.1`

- Lung lipid default: the five configured scALABLE lung references now select the
  native-relative-log2 release, matching the standalone API. The saved artifact ID
  was already inventoried; its inference-code ID has changed. The release manifest
  declares 202 lipids, 1,303 RNA inputs and 45 training profiles, with original
  ElasticNetCV and an explicitly authorized donor exclusion.
- AML metabolite default: API and bone-marrow scALABLE select the NA30 release.
  The manifest records the authorized >30% target-missingness rule: 2,023 of 2,533
  targets retained; all 84 training cases and 12,416 RNA inputs preserved. Retained
  weights are unchanged; no refitting or new CV was performed by that release.
- Shared scALABLE source-scale conversion, arithmetic-linear pseudobulk aggregation,
  fold and BH code now has an explicit application-method ID, separate from model
  weights. Registry changes do not alter these numerical algorithms.
- Unchanged artifact IDs remain unchanged. Every previously registered artifact and
  complete prior catalog snapshot is retained. Full model IDs, source hashes and
  seven changed default selections are in `update.json`.

The two model bundles and SHA256SUMS are attached to the private GitHub release.
Source release manifests are stored here. Metadata records release authors' method
and authorization statements; registry publication does not independently certify
scientific performance. Lung release activation remains subject to the recorded
historical-API integrity decision. No service restart or biological rerun was performed.
Repository/default observation times are separate from actual deployment dates.
