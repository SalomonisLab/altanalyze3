# AML metabolite identity audit — 17 September 2026

The unidentified labels originate in the study data. They were not substituted
for compound names by scALABLE. No verified replacement identities were recovered.
This audit restores the source annotations and coordinates needed to investigate
them without changing model IDs, predictions, or saved differential results.

## Verified source chain

1. Chu et al., [Nature Cancer, DOI 10.1038/s43018-026-01175-6](https://www.nature.com/articles/s43018-026-01175-6).
2. [Supplementary Tables 1–40](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs43018-026-01175-6/MediaObjects/43018_2026_1175_MOESM3_ESM.xlsx):
   Table21, positive-mode HILIC; Table22, positive-mode reversed-phase (RP).
3. Local training extraction in `Human-MS-impute/code/extract_ms_matrices.py`
   retains the literal `Metabolite` names. `build_unique_ms_tables.py` selects
   one assay per repeated name using detection fraction, then median signal.
   Its selected-feature annotation file kept assay/detection/signal but omitted
   the richer mass-spectrometry annotation columns. The names themselves survived.
4. `rna2metabolite_aml_bundle.pkl.gz` contains those 2,533 selected targets.
5. The metabolite matrix in saved session `1b0faa928a9848fbb3cd225ca26a0e42`
   has exactly the same target IDs and order as the model.
6. Independent public-deposit check: [PDC000561, CPTAC AML Study — Metabolome](https://pdc.cancer.gov/pdc/study/PDC000561),
   `cptac.metabo.table.csv`, file UUID `0286f9e2-69bc-11f1-b328-0a200668d9ef`.
   Downloaded MD5 `3f51f1f144b3ec8672a65bd84c738ab9` matches the PDC manifest.
   All 2,655 assay-qualified row IDs match the workbook, including the unknowns.

## Counts and recovered metadata

| Scope | Named | Unknown | Total |
|---|---:|---:|---:|
| Original HILIC + RP measurement rows / PDC table | 540 | 2,115 | 2,655 |
| Unique model targets / saved scALABLE session | 418 | 2,115 | 2,533 |
| Paper's filtered ICA analysis, Supplementary Table24 | 433 | 1,347 | 1,780 |

The article reports 322 unique named metabolites in its filtered analysis and
states that the 1,347 remaining features could not be identified using the
spectral libraries tested. The model uses the broader source set, not just that
filtered subset. The 122-row reduction from 2,655 to 2,533 is the training
procedure's collapse of repeated names across assays. Of the assay rows selected
for the model, 1,675 occur in Table24. These different feature sets explain the
different totals; they do not indicate lost names.

The recovered catalog includes all 2,533 target IDs, assay, sheet and exact Excel
row, m/z, retention time, source tag, synonyms/isomers, formula, KEGG identifier,
Table24 membership, DOI, URL and verbatim source annotation fields. Among the
unknowns, 33 have formula text, one has candidate text, and zero have a KEGG ID.
Among the 418 named targets, 205 have a source KEGG annotation. A source name
does not establish a uniform identification-confidence level.

| Example | Exact source location | m/z | RT, minutes | Result |
|---|---|---:|---:|---|
| Unknown 0235 | Table22!A464 | 234.09706 | 1.221 | Unidentified in workbook and PDC |
| Unknown 0309 | Table22!A538 | 262.07476 | 6.054 | Unidentified in workbook and PDC |
| Unknown0 1621 | Table21!A1622 | 92.07052 | 5.364 | Candidate annotation needs clarification |

`Unknown0 1621` has `Sarcosine;` in `Synonyms.or.Isomers`, but the source formula
is `C3 H9 N O2`. [PubChem's sarcosine record](https://pubchem.ncbi.nlm.nih.gov/compound/1088)
gives C3H7NO2. Its neutral formula and expected protonated mass do not match this
row. We therefore retain the candidate as a provenance note, not as a name mapping.
Furthermore, Table21 columns H:K have apparent header/value inconsistencies
(e.g. an adduct string under `Annot.DeltaMass.ppm`). The catalog preserves them
literally and does not silently shift columns or infer confidence from `MS2`.

## Additional sources checked

- The PDC API lists 210 raw mass-spectrometry files and one processed CSV for
  PDC000561. That CSV provides measurements, not an updated compound crosswalk.
  Raw files were not downloaded; public raw data are available if reannotation
  is pursued.
- [Author analysis repository](https://github.com/Nesvilab/CPTAC_AML), commit
  `07ffd50fb0c5e4fe57d8147d1ca9359482469c0d`: the metabolomics script references
  the original PNNL HILIC/RP spreadsheets and explicitly analyzes unknown
  features. No additional identity crosswalk is included.
- [FunMap repository](https://github.com/bzhanglab/funmap_aml), commit
  `3abd736a355150796bb859a5899f42e91315fda4`: repository inventory contains
  gene-network display data; no metabolite annotation file was located.
- Supplementary Table24 and Source Data Extended Data Fig.4 retain 1,347
  unknowns. Source Data Fig.3 supplies selected biological analyses, not an
  unknown-to-compound lookup. Other locally available figure source workbooks
  (MOESM4,7,9,11) yielded no such lookup.
- The currently downloadable Supplementary Information PDF (MOESM1) contains
  consortium membership, not the detailed metabolomics identification workflow
  referenced in the article. The Reporting Summary (MOESM2) does not supply the
  missing feature identities.

These checks establish the status of the inspected public files, not that no
newer private annotation exists.

## What would enable further mapping

No additional paper download is needed. The useful next input is an updated
feature-to-compound annotation export from the study's PNNL metabolomics team
(Jennifer E. Kyle / corresponding author Tao Liu), if one exists. A ready-to-send
request is in [AUTHOR_REQUEST.md](AUTHOR_REQUEST.md); it has not been sent.

The export should link the exact HILIC/RP unknown IDs to compound candidates,
adduct/charge, MS/MS library accessions and match scores, identification confidence,
and stable chemical identifiers. An accompanying feature-to-spectrum/scan mapping
and MGF/MSP or mzML export would permit new spectral searches if the authors
have no newer assignments. Formula or precursor mass alone is insufficient to
choose among isomers or adduct explanations. New spectral hits would be retained
as candidates until supported by the evidence appropriate to their confidence.

Updated annotations can join to this catalog by `(assay, source_name)` and be
checked against mass and retention time. Keep `model_feature_id` stable; store
any accepted display name alongside its source, version, match evidence and
confidence. Do not infer a metabolite's identity from the genes that predict it.

## Files and validation

- [aml_metabolite_provenance.tsv](aml_metabolite_provenance.tsv): full catalog.
- [unidentified_features_for_authors.tsv](unidentified_features_for_authors.tsv):
  unknown rows ready for annotation review.
- [audit.json](audit.json): source checksums, counts and exact-ID comparisons.
- [session_validation.json](session_validation.json): model/session validation.
- [source_manifest.json](source_manifest.json): downloaded-source hashes and PDC metadata.

Reproduce using `components/rna2metabolite/provenance.py --workbook <MOESM3.xlsx>
--selection <metabolomics_unique_annotation.tsv> --pdc-table <cptac.metabo.table.csv>
--output <audit-directory>`. The command reads source files without modifying them.
The source CSV can be obtained through PDC's documented `filesPerStudy` API using
PDC000561 and its filename. Two regression tests verify assay-specific joins,
literal unknown retention, candidate separation, source coordinates, zero-valued
metadata and rejection of missing deposited IDs. Both pass, and all 2,533 real
targets pass the model/session checks. No app behavior or analytical values were
changed by this audit.
