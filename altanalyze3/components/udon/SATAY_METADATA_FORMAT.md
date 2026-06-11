# SATAY-UDON metadata file — format specification

SATAY-UDON tests each UDON cluster for enrichment of sample covariates. The covariates are
supplied as **one standardized, tab-delimited text file** passed with `--metadata`. This is the
**only** metadata input the tool accepts — there is no built-in/default file and no
dataset-specific logic in the software.

## Requirements (enforced; the tool errors if violated)

1. **Tab-delimited** plain text (`.tsv` or `.txt`), UTF-8, with a single header row.
2. Exactly one **`Sample`** column (required, case-sensitive). Values are the sample identifiers
   and must match the `Sample` portion of the pseudobulk IDs (`celltype__Sample`). Duplicate
   `Sample` values are an error.
3. Optional **`Donor_ID`** column. If present, enrichment is collapsed to **distinct donors**
   (≥3 positive donors required) to avoid pseudoreplication. If absent, testing is per sample.
4. Optional **`Study`** column (batch/cohort label; carried for context).
5. **Every other column is a covariate** and must be cleanly one of:
   - **binary** — values are a subset of `{0, 1}` (blanks allowed) → a single *presence* test;
   - **categorical** — string levels → one test *per level* that has ≥ `MIN_LEVEL` (5) samples.
   A column that is neither (free text, or a continuous numeric such as raw age or VAF) is an
   error — **continuous variables must be pre-binned** by the user (e.g. `age_ge_60` as 0/1).
6. Missing values: empty, `NA`, `nan`, `Unknown`, `None`, `n.a.`, `na` → undefined for that
   sample (excluded from that covariate's test).
7. There must be ≥1 covariate column and ≥1 `Sample` overlapping the pseudobulks.

## Example (`metadata_template.tsv`)

```
Sample      Donor_ID   Study    NPM1   FLT3_ITD   sex      ELN_risk        age_ge_60
PT01        D01        NYU-1    1      0          M        Adverse         1
PT02        D02        NYU-1    0      1          F        Intermediate    0
PT03        D03        WashU    0      0          F        Favorable       0
PT04        D03        WashU                      M        Intermediate    1
```

- `NPM1`, `FLT3_ITD`, `age_ge_60` → binary (presence tests).
- `sex`, `ELN_risk` → categorical (one test per level ≥5 samples).
- `PT04` has blank mutation cells → undefined for `NPM1`/`FLT3_ITD`, still tested for `sex`/`ELN_risk`.
- `PT03`/`PT04` share `Donor_ID=D03` → counted as one donor.

## Usage

```
python satay_udon.py --metadata /path/to/metadata.tsv
```

## Converting an existing (non-standard) source

Dataset-specific cleaning/harmonization (e.g. collapsing redundant ELN_risk labels, or merging an
xlsx's separate mutation/clinical tabs) is **not** part of the tool. Write a small prep script that
emits this standard file, and keep it **in the analysis/output directory with your data** (not in
the software). See `UDON/prepare_satay_metadata.py` for the AML example.
