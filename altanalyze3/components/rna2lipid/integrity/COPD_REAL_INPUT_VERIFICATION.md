# COPD Zhang local imputation input verification

Verified 2026-10-05 by read-only inspection. The local original imputation and prior candidate replay use the real per-sample COPD inputs.

| Input | Rows | Samples | Donors |
|---|---:|---:|---:|
| copd_pb_ln1pcp10k.h5ad | 4,781 | 141 | 141 |
| copd_mc_ln1pcp10k.h5ad | 84,646 | 141 | 141 |

Both are under `/Volumes/salomonis2/LungMAP/CellRef2/inputs`. Their study acronym is `Zhang 2026 NatGenet`.

The complete COPD metacell membership file was joined through CellBarcode -> Library -> Sample/Donor using the uncensored harmonized metadata. All 1,678,054 member cells mapped. There were zero mixed-sample metacells, zero mixed-donor metacells, zero sample/donor label mismatches and zero membership-count mismatches. All 84,646 input metacell IDs matched the membership roster.

The pseudobulk builder starts from the library-by-cell-state UNION counts, maps 148 libraries to 141 real samples, and sums libraries within each sample and cell state before normalization. The source metacell builder groups by container, study and real Sample. It explicitly distinguishes this v8 build from pooled COPD release metacells.

The older `CellRef2.0_metacells_v7_with_COPD.h5ad` release contains 83,416 COPD metacells and uses meta-samples. It is not the source named in these local imputation inputs.

Provenance and numeric checks: [JSON audit](COPD_real_input_audit_20261005.json). Original builders and imputation driver are recorded there. The previously saved candidate inference audit also names these inputs and reports production replay errors below 0.000001 for COPD.

Separate scale observation: these historical COPD RNA inputs use natural log1p(CP10k), whereas the other atlas replay inputs use log2(1+CP10k). No normalization or model change was made in this verification. Donor separation does not establish that both scales are appropriate for the model.

The user authorized per-lipid median NA filling on log2 training targets, retaining negative observed values, and reaffirmed the deployed lipidwise ElasticNetCV. This decision is recorded in decision_state.json. This audit did not retrain a model or rerun differentials.
