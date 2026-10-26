# Lung rna2metabolite bundles, 2026-10-21

Two lung bundles sit in `artifacts/` beside the AML default. `api.DEFAULT_BUNDLE_PATH` still names the
AML NA30 release; load a lung bundle by path. Both bundles come from `training.py` with the AML recipe
(1000 genes per metabolite, ridge alpha 100, protein-coding genes with TPM >= 5 in >= 3 samples,
CP10k+log1p, Recon3D prior). `validation/reproduce_aml_bundle.py` shows that this code rebuilds the AML
model exactly.

| Bundle | SHA-256 | Training data | Held-out | Null (50 permutations) |
|---|---|---|---|---|
| `rna2metabolite_lung_monash_A_aml_recipe_bundle.pkl.gz` | `09cba19c83282568af01ba235f737de783b532b5abd1a8f93f2471cfb9659731` | Monash IPF explants: 22 samples, 11 donors (6 IPF, 5 non-diseased), 465 MS-DIAL metabolites | leave one donor out: median Spearman 0.215, 191 of 465 > 0.3 | null 95th 0.076; p = 0.0196 (floor) |
| `rna2metabolite_lung_lund_A_aml_recipe_bundle.pkl.gz` | `0f76836b32ce7762647ab5d6f831dcdbf1933b4b34c54686a55d3c2944c6de88` | Lund LUAD: 37 tumors, 139 targeted HILIC metabolites | KFold 5: median Spearman 0.418, 83 of 139 > 0.3 | null max 0.050; p = 0.0196 (floor) |

Input scale for both: raw gene counts by symbol; the bundle applies CP10k over all input genes. The Lund
bundle was trained on counts back-computed from StringTie FPKM (FPKM x Ensembl 103 union-exon length). On
37 Monash libraries that back-computation kept prediction ranks (per-metabolite r median 0.987) and shifted
levels by a median 0.88 SD against true counts.

Cross-cohort agreement is weak: Monash A on Lund tumors median Spearman 0.057 (15 of 92 matched
metabolites > 0.3); Lund A on Monash samples 0.082 (30 of 92 > 0.3). Use each bundle for its own tissue.

## Where everything else lives

The source data, processing code, logs, the nested-sweep results and the matched-donor rna2lipid test are
outside this repository, in `/data/salomonis-archive/LabFiles/Nathan/Other/rna2metabolite`. Its
`README.md` holds the methods, every result table and the file list.

| Source | Location |
|---|---|
| Monash RNA (re-aligned, bulk-rnaseq-splicing counts mode) | `cohorts/Monash_IPF_GSE213001/featurecounts/counts_matrix.tsv` |
| Monash MS | Mendeley 10.17632/4zr8jjmzhy.1 (`cohorts/Monash_IPF_GSE213001/raw_downloads/mendeley_4zr8jjmzhy_v1/`) |
| Lund RNA | ArrayExpress E-MTAB-16082 `SUH_RNAseq_FPKM.txt` |
| Lund MS | Arbajian et al. 2025 Table S2 sheet A; Metabolomics Workbench ST003599 |
| Training drivers | `code/train_monash_rna2metabolite.py`, `code/train_lund_and_cross_cohort.py` |
