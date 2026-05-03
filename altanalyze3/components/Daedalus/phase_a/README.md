# Daedalus Phase A

Phase A establishes the data and benchmark substrate for the `Daedalus` model.

## Goals

1. Acquire core external resources
2. Normalize them into a shared schema
3. Build isoform-pair supervision tables
4. Train strong non-foundation baselines before any large paired model

## Deliverables

- `resources/manifest.json`
- `schemas/isoform_pair_record.schema.json`
- `scripts/download_resources.py`
- `scripts/init_phase_a.py`
- `data/raw/`
- `data/interim/`
- `data/processed/`
- `benchmarks/`

## Resource classes

- transcript reference sets
- isoform reference labels
- protein function / domain annotations
- topology / localization annotations
- curated transcript anchors (`APPRIS`, `MANE`)
- curated family membership annotations
- pathogenic / consequence annotations

## Initial download policy

The initial automated download step only fetches resources that are:

- public
- reasonably sized
- stable enough for scripted retrieval

Large resources such as full `InterPro protein2ipr.dat.gz` are listed in the
manifest but disabled by default because they are large enough to warrant a
deliberate fetch.
