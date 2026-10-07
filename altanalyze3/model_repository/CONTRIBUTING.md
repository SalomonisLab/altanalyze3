# Contributing

See README.md for proposal requirements and the promotion procedure. Add metadata under
`proposals/`; host binary models in immutable, licensed HTTPS release assets or archives.
Do not submit patient-identifying sample metadata. A stable deidentified ordered roster
is sufficient to support identity and coverage checks.

Contributors must not edit the default catalog, deployment history or approval records
as part of a proposal. Maintainer review is required before loading executable pickles,
changing application selectors or recording a default promotion.


Include `modality`, `version` (for example `v1.1`) and `variant` (for example `lung`)
in proposal metadata. Maintainers verify the version against existing entries before
merging. A revision increments the affected modality/variant series; unchanged tissue
models keep their versions. Never reuse a version for different artifact hashes.
