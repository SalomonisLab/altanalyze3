# rna2psi — validation

The **current default model (v2, 2026-07-09)** — event set, statistical thresholds,
features, training/assessment datasets, and full BEAT-AML performance — is documented in
**`METHODS.md`** (§5–6). Headline: trained on all 437 Leucegene AML samples; assessed
cross-cohort on BEAT-AML (671 clinically-annotated samples) at a 98%-specificity operating
point, where splicing-factor (SF3B1 0.57, ZRSR2 0.50, SRSF2 0.39) and NPM1 (0.25, precision
0.77) signatures transfer with usable precision.

The **prior v1 model** (protein-coding genes, ElasticNet/Ridge CV-pick, 18,168 events;
leakage-free honest OOF median Spearman 0.646) and its validation notes are archived in
`legacy/v1_2026-07-07_protein_coding_marginal/`.
