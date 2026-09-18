# Modality imputation in scALABLE — reviewer statement

Drafted 2026-09-18. Word count of the statement below: 189 of a 200-word limit.

## Statement

scALABLE imputes four modalities from RNA. Every model is a per-target penalized linear
regression fit on paired multi-omic data, then applied to aligned query cells.

**ADT (129 proteins):** per-protein ElasticNet (alpha 0.01, l1_ratio 0.5) over a 3,132-gene
whitelist, trained on a CITE-seq bone-marrow atlas. Leave-one-donor-out mean Pearson r = 0.737
(median 0.807, RMSE 0.494) matches the cell holdout (0.726 / 0.782), so accuracy does not depend
on donor identity.

**Metabolites (2,533) and lipids (1,009):** per-molecule ridge (alpha 100) over Recon3D
mechanistic prior genes, trained on CPTAC-3 AML (84 and 87 cases). Out-of-fold 5-fold CV median
Spearman is 0.267 and 0.355; 1,084 and 625 molecules exceed Spearman 0.3, giving classification
accuracy 0.663 and 0.678 (AUROC 0.718 and 0.723).

**TF activity (217 TFs, 7,486 edges):** per-edge ridge on three covariates, target expression,
TF expression and regulon mean. Leave-one-donor-out R2 = 0.773. Among differential TFs the
predicted direction is correct for 84.5% (top quartile, TEA control, 71 pseudobulks), median
0.83 per donor across 27 held-out donors.

**Specificity:** raising one regulon while holding TF mRNA constant ranks the perturbed TF
median 2 of 217 (top-5 83%); kNN retrieval ranks it 74 (top-5 5%).

## Source of every number

| Claim | File |
|---|---|
| ADT algorithm, 3,132 genes, 129 ADTs | `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/rna2adt_bm_bundle.pkl` (`metadata` key) |
| ADT Pearson r, RMSE | `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/artifacts/all_benchmark_v10_no_isotype/summary.tsv` |
| Metabolite and lipid CV | `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2metabolite/VALIDATION.md`, `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2metabolite/README.md` |
| TF activity, perturbation test | `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2grn/VALIDATION.md` sections 4, 5 and 9 |

## Scope and exclusions

- The statement covers the human bone-marrow bundles that scALABLE loads by default.
- The statement excludes RNA. scALABLE measures RNA and does not impute it.
- The statement excludes the lung bundles. Those are separate models. The lung lipid bundle
  reports Pearson r 0.871, and that value is not donor-held-out
  (`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2lipid/VALIDATION.md`
  lines 159 to 183).

## Two numbers a reader must not take from elsewhere

1. `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/artifacts/summary.tsv`
   reports mean Pearson 0.8975. That row is model `centroid_resid_a10.0`, which scALABLE does not
   ship. The same model falls to 0.5275 under leave-one-donor-out. The shipped bundle is
   `rna2lipid_arch_panel_per_protein_whitelist_union`.
2. The ADT counterpart of the TF specificity test is
   `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/aberrant_test.py`.
   For the shipped bundle it recovers CD8 0.818, CD14 0.212, CD19 0.196 and CD34 0.059 of the true
   gap, over 4 probes at 1,000 cells each. One probe of 4 supports the claim, so the statement
   omits ADT from the specificity sentence.

## Not yet computed

The reviewer also asked for ARI against Azimuth. No cellHarmony-versus-Azimuth label comparison
exists in this repository or under `/Users/saljh8/Dropbox/Manuscripts`. That statement needs the
atlas samples and a run.
