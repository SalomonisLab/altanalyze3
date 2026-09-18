# Draft request — not sent

We are tracing metabolite identities from the CPTAC AML study
(10.1038/s43018-026-01175-6; PDC000561) into an RNA-based imputation model.
Supplementary Tables 21/22 and the deposited `cptac.metabo.table.csv` agree on
the original feature names, including 2,115 unknown entries in the broader
measurement tables. We recognize that the paper reports 1,347 unidentified
features after filtering.

Do you have a newer HILIC/RP feature annotation export with assignments for any
of these unknowns? We can provide the attached `unidentified_features_for_authors.tsv`,
which preserves the exact original IDs, source rows, m/z and retention times.

The most useful fields would be:

- Original feature ID and assay (HILIC versus RP), m/z and retention time.
- Proposed name, formula, adduct/charge, and KEGG/HMDB/PubChem/InChIKey identifiers.
- Identification-confidence level, reference-standard support, library name,
  library accession and MS/MS match score.
- Feature-to-spectrum or scan identifiers; MGF/MSP or mzML exports if available,
  particularly if no newer compound assignments exist.

Your analysis script references
`Metabolites_HILIC_Pos_log2_GlobalMedian_Backtrans_20230426.xlsx` and
`Metabolites_RP_Pos_log2_GlobalMedian_Backtrans_20230502.xlsx`. Do their original
annotation sheets or later versions contain information omitted from the published
tables? Could you also provide the detailed metabolomics identification methods
and library search settings referenced in the article?

Two specific annotation questions:

1. Table21!A1622 (`Unknown0 1621`, m/z 92.07052, RT 5.364 min) lists `Sarcosine;`
   as a synonym/isomer but formula `C3 H9 N O2`. Is that candidate intentional,
   an isotope/adduct annotation, or a spreadsheet artifact?
2. Table21 H:K labels appear misaligned with their values (for example,
   `[M+H]+1` occurs under `Annot.DeltaMass.ppm`). What are the correct column
   meanings?

We will retain the original feature IDs and distinguish tentative annotations
from confirmed identities rather than assigning names from precursor mass alone.
