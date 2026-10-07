# scALABLE imputation models

Models are organized by **modality → version → tissue/species**. Open a modality to
compare versions, then select the variant appropriate for your analysis.

| Modality | Available versions | Current defaults |
| --- | --- | --- |
| [rna2lipid](models/rna2lipid/README.md) | v1.0, v1.1 | Lung **v1.1**; AML **v1.0** |
| [rna2metabolite](models/rna2metabolite/README.md) | v1.0, v1.1 | AML **v1.1** |
| [rna2adt](models/rna2adt/README.md) | v1.0 | Human bone marrow, human lung, mouse **v1.0** |
| [rna2grn](models/rna2grn/README.md) | v1.0 | Leukemia and lung **v1.0** |
| [fastComm](models/fastComm/README.md) | v1.0 | Human and mouse **v1.0** |

```text
models/
  rna2lipid/
    v1.0/
      lung/           # previous lung model
      aml/            # current AML model
    v1.1/
      lung/           # current lung model
  rna2metabolite/
    v1.0/aml/         # previous AML model
    v1.1/aml/         # current NA30 AML model
  rna2adt/v1.0/
    human-bone-marrow/
    human-lung/
    mouse/
  rna2grn/v1.0/
    leukemia/
    lung/
  fastComm/v1.0/
    human/
    mouse/
```

Each model directory contains a readable `README.md` with downloads and a `model.json`
with the full identity, file hashes, method records and provenance. Model artifacts are
attached to modality releases such as **rna2lipid-v1.1** and **rna2metabolite-v1.1**.
[Current model defaults](models/current.json) also have a machine-readable index.

## Version rules

Versions start at **v1.0** for the first registered model in each modality/variant
series and progress through **v1.1**, **v1.2**, etc. when that model is revised.
Different tissues/species are distinct model variants; changing the lung model does
not renumber the unchanged AML model. These registry labels are not original training
dates or claims that different tissue models are interchangeable.

Analysis results record readable names such as **rna2lipid/lung/v1.1** together with
immutable artifact and inference-code hashes. Changing a file cannot reuse its old
hash identity. Archived identities and old model versions remain available.
Version labels already assigned to an artifact must not be changed or reused.
Maintainers assign a new version after review; a community proposal never makes itself
a default. Shared scALABLE downstream code has a separate method identity.

## Updates and contributions

[Recent update](releases/registry-2026.10.06.1/NOTES.md): native-log2 lung lipids and the
AML metabolite NA30 panel. Default entries reflect the checked-in configuration;
actual running-service deployment dates require separate evidence.

Submit proposed models through [CONTRIBUTING.md](CONTRIBUTING.md). Specify the modality,
proposed version and tissue/species, with complete source, method and evaluation records.
See [detailed provenance and default history](docs/PROVENANCE.md),
[validation](VALIDATION.md), and [default history](default_history.json).

This repository and its model release assets are private to SalomonisLab.
