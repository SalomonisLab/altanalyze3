# Variant Discovery from Long-Read Single-Cell BAM

Per-cell genotyping (SNV / indel → cell barcode → WT/MUT) directly from a long-read
single-cell RNA-seq BAM, **given a file of known variant positions**. This is the "discovery /
read-genotyping" half of the workflow; the downstream expression-marker imputation that expands
these calls to undetected cells is documented in **[VARIANT_IMPACT.md](VARIANT_IMPACT.md)**.

All code lives in `altanalyze3/components/bam/`. Nothing here was invented for this document —
every command, path, and number below is copied from the source or from the production run
directory `/Volumes/salomonis-archive/LabFiles/Nathan/Revio/MDS-AML-KINNEX-4/variant_impact_RUNX1_prototype/`
(cluster path: `/data/salomonis-archive/...`).

---

## 1. Modules

There are three genotypers, ordered by how directly they answer *"call these variants in this BAM"*:

| Module | Role | Reference FASTA? | Barcode tags | CLI |
|---|---|---|---|---|
| `variant_impact.py` → `genotype_from_bam()` | Reference-free per-cell caller from the `=`/`X` extended CIGAR. Core of the `variant-impact` subcommand. | **No** (BAM encodes match/mismatch) | `CB,BC,XC,UB` (first found) | in central CLI (`altanalyze3 variant-impact`) |
| `variant_extraction.py` | Standalone supervised SNV **and** indel caller; forensic per-read table + per-cell genotype matrix. | **Yes** (`--reference genome.fa`) | `CB,BC,XC,UB` (first found) | standalone `python3 variant_extraction.py ...` |
| `global_snv.py` (+ `combineVariants.py`) | Older `pysam.pileup` path; de-novo genome scan **or** supplied `chr:pos` list. | Yes (positional) | `CB` only | standalone |

For long-read KINNEX/Revio BAMs that carry the extended CIGAR, `genotype_from_bam()` is the
primary path (no FASTA needed). `variant_extraction.py` is used when you want the full per-read
forensic table or a plain-`M` BAM (which needs the FASTA). The two agree in practice: on
sample `5801M_pre` both returned **256 MUT / 435 WT** for RUNX1 W279* (see §5).

---

## 2. Variant-position file format

Tab-delimited, one variant per line. This is the file consumed by both `variant-impact`
(`--variants <file>`) and `variant_extraction.py` (`--mutations <file>`):

```
chrom   start   end        label          type
chr21   34799432  34799432  RUNX1p.W279*   SNV
chr17   76736877  76736877  SRSF2p.P95R    SNV
```

- `type` ∈ `SNV` (default) | `Indel`.
- `label` carries the protein change and is propagated to output paths / genotype columns.
- `chr` prefix is auto-added if absent.

Real files used (verbatim):
`.../variant_impact_RUNX1_prototype/RUNX1_mutation.txt` and `.../SRSF2_mutation.txt`.

> **Parser note (verified in code).** `variant_impact._parse_variants` reads columns
> `[0]=chrom [1]=pos [3]=label [4]=type`. `variant_extraction.parse_mutation_file` reads them
> positionally as `chrom pos ref alt mode gene`, so for the shared 5-column file only column 1
> (`pos`) and column 5 (`type`) drive the calling logic — columns 2–4 (`end`, `label`) are ignored
> or mislabeled internally by `variant_extraction.py`. This is harmless for SNV/indel calling
> (the reference base comes from the FASTA, the observed base from the reads) but is a real
> docstring-vs-parser inconsistency to be aware of.

`variant-impact` also accepts an **inline** spec instead of a file:
`--variants 'chr21:34799432:RUNX1_W279*;chr17:76736877:SRSF2_P95R'` (`chr:pos:label[:type]`,
semicolon-separated).

---

## 3. `genotype_from_bam()` — the reference-free caller

`altanalyze3/components/bam/variant_impact.py`

```python
genotype_from_bam(bam_path, chrom, pos1, vtype="SNV",
                  barcode_tags=("CB","BC","XC","UB"), min_mapq=20, indel_window=5)
    -> DataFrame[core_barcode] -> (mut, wt);  df.attrs["ref_allele"], df.attrs["alt_allele"]
```

Logic (only the **two major alleles** are counted; minor/error alleles are dropped):

- **SNV**: `reference allele` = consensus base of `=` (match) reads (= the genome reference);
  `variant allele` = major base of `X` (mismatch) reads. `mut` = variant-allele reads,
  `wt` = reference reads.
- **Indel**: `variant allele` = the major indel `(type,length)` carried by the most reads within
  `indel_window` bp (absorbs left/right-alignment); `wt` = spanning reads with no indel. Fetches
  a `±indel_window` window so shifted alignments of the same event are seen.
- Filters: `is_unmapped`, `mapping_quality < min_mapq`, missing barcode.
- Barcode: first tag found among `CB,BC,XC,UB`, reduced to its `_core` (strips `-1`/`.sample`).

**Hard limitation (emitted as a warning in code):** if the BAM uses plain `M` CIGAR (no `=`/`X`),
`ref_allele` is `None` and the locus **cannot be called** — a reference FASTA would be required
(use `variant_extraction.py` for that case):

```
genotype_from_bam chr:pos uses standard 'M' CIGAR (no =/X) -- a reference FASTA would be needed to call this locus.
```

---

## 4. `variant_extraction.py` — standalone forensic caller

`altanalyze3/components/bam/variant_extraction.py`

### CLI (verified argparse)

```
python3 variant_extraction.py \
  --sample  <name>          # --s
  --bam     <scisoseq.mapped.bam>   # --b   (BAM index required)
  --mutations <variants.txt>  # --m   (format in §2)
  --reference <genome.fa>     # --r   (FASTA index .fai required)
  --output-dir <dir>          # --o   (docstring says --output; the flag is --output-dir)
```

### Modes (from the module docstring)

- **SNV mode**: examine the exact variant position only → `V` (variant) or `R` (reference) per read.
- **Indel mode**: check for an indel within ±50 bp; only a 0/+1-distance indel is called `I`,
  otherwise `R`.

### Key functions

- `analyze_comprehensive(bam_file, mutations_list, ref_genome, context_bases=5, min_qual=20)` —
  opens BAM (`pysam.AlignmentFile(...,"rb")`) + FASTA (`pysam.FastaFile`); emits per-read/per-cell
  rows with allele status, base/mapping quality, indel details, homopolymer length, tandem-repeat
  flag, strand.
- `create_mutation_matrix(input_file, output_file, sample_name)` — collapses the per-read TSV to a
  per-cell genotype matrix. **Reverse-complements every barcode** before matrixing (orientation
  assumption baked in — relevant when joining to an h5ad; see VARIANT_IMPACT.md §orientation).

### Outputs (per sample)

- `<sample>_complete_analysis.tsv` — ~32-col per-read forensic table
  (`cell_barcode, chr, position, ref, alt, gene, ..., observed_base, allele_status (R/V/I),
  base_quality, mapping_quality, num_indels, ..., position_total_reads, position_alt_reads`).
- `<sample>_mutation_matrix.csv` —
  `cell_barcode, genotype, mutations_present, <label>_WT, <label>_MUT, <label>_GENOTYPE`.
- `<sample>_mutation_matrix_summary.txt` — WT/MUT cell counts and percentages.

### Known caveats (visible in the source — do not treat these columns as trustworthy)

- Strand-bias Fisher test is **degenerate**: `fisher_exact([[f,r],[f,r]])` compares identical
  rows → p ≈ 1.0 always.
- Several PCR-artifact columns are hard-coded placeholders (`"0.0","0.0","1.0000"`; `context_changes = 0`).
- Indel calling only counts indels at distance ≤ 1 as `I`, despite fetching a 200 bp window.

---

## 5. Exact commands and inputs actually run (sample `5801M_pre`)

**BAM** (long-read, PacBio Revio KINNEX, 10x barcodes):
```
/Volumes/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/2-16-23PabioRevio/SamplePre34_outputs/BAM-clusters/scisoseq.mapped.bam
```
(library `WM71-0121-5801_CD34`, group `MDS-pre-SRSF2`; index required.)

**Reference** (only for `variant_extraction.py`; symlinked in the run dir):
`genome.fa -> /Volumes/salomonis-archive/LabFiles/Nathan/Revio/Hs/genome.fa`
(only `chr17.fa`/`chr21.fa` staged locally with `.fai`).

### 5a. Standalone forensic call (`variant_extraction.py`)

```bash
/opt/homebrew/opt/python@3.11/bin/python3.11 variant_extraction.py \
  --bam       /Volumes/.../SamplePre34_outputs/BAM-clusters/scisoseq.mapped.bam \
  --sample    WM71-0121-5801_CD34 \
  --mutations RUNX1_mutation.txt \
  --reference genome.fa \
  --output-dir variants
```
Produced `variants/WM71-0121-5801_CD34_*` (matrix + per-read + summary). The SRSF2 run wrote to
`variants_SRSF2/`.

### 5b. Reference-free call inside `variant-impact`

The `variants` file is genotyped straight from the BAM as step 1 of the `variant-impact`
subcommand — no separate command needed (see VARIANT_IMPACT.md §Commands).

### Genotype results (from `*_mutation_matrix_summary.txt`)

| Variant | Total genotyped cells | WT | MUT |
|---|---|---|---|
| RUNX1 W279* (chr21:34799432) | 691 | 435 (63.0%) | 256 (37.0%) |
| SRSF2 P95R (chr17:76736877) | 3339 | 1286 (38.5%) | 2053 (61.5%) |

The `variant-impact` reference-free genotyper reproduced RUNX1 exactly (combined confident set
`256 MUT / 435 WT` at `--mut-min 1 --wt-min 1`), a useful cross-method sanity check.

---

## 6. Limitations of variant discovery (the real ones)

1. **Long-read coverage is ~1 read/cell.** WT calls need ≥2 reference reads by default (`--wt-min 2`).
   For a largely-mutant clone (e.g. RUNX1 here) this collapsed WT to ~20 cells total → the runs were
   redone at `--wt-min 1` (435 WT). But `≥1`-read WT admits **allelic-dropout false-WT** — a cell
   with the mutation but only the reference allele sequenced is miscalled WT. This is the central
   accuracy tension of read-genotyping at this depth.
2. **`genotype_from_bam` needs `=`/`X` CIGAR.** Plain-`M` BAMs are refused (warning above); use
   `variant_extraction.py` + FASTA instead.
3. **Two-major-allele model.** Multiallelic sites and minor/error alleles are ignored by design.
4. **Barcode orientation.** Long-read isoform BAM barcodes are the **reverse-complement** of the
   matched short-read/junction barcodes; `genotype_from_bam` handles this via `_pick_orientation`,
   but `variant_extraction.create_mutation_matrix` unconditionally RCs — mismatched orientation is
   the first thing to check if a join yields 0 overlapping cells.
5. **`variant_extraction.py` QC columns are placeholders** (§4) — do not report its strand-bias or
   PCR-artifact fields as real.
6. **`global_snv.py`** reads only the `CB` tag and has a hard-coded input path in `combineVariants.py`;
   treat it as legacy.

---

## 7. Reproduce

```bash
PY=/opt/homebrew/opt/python@3.11/bin/python3.11    # numpy/pysam-equipped interpreter
cd /Volumes/salomonis-archive/LabFiles/Nathan/Revio/MDS-AML-KINNEX-4/variant_impact_RUNX1_prototype
# ensure BAM and reference .fai indices exist (pysam requires them)
$PY /Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/bam/variant_extraction.py \
    --bam .../scisoseq.mapped.bam --sample WM71-0121-5801_CD34 \
    --mutations RUNX1_mutation.txt --reference genome.fa --output-dir variants
```

**Reproducibility gaps to close (observed in the run dir):** no `run_parameters.json`, no pinned
altanalyze3 commit recorded, reference genome is a symlink (not self-contained). For a publishable
run, pin the commit and stage/record the exact reference build alongside the outputs.

---

*Cross-reference: `variant-impact-deg-concordance` memory; `bam2hla` shares the pysam-pileup
infrastructure style. Downstream expansion → [VARIANT_IMPACT.md](VARIANT_IMPACT.md).*

---

## 8. Code changes — 2026-07-19 (MDS-AML-KINNEX-5 run of record)

Made so the whole supervised scan runs from the installed package rather than ad hoc scripts, and so
indel calling detects the variant it was handed. Sections 4–7 above describe the pre-change module;
the CLI additions and behaviour changes below supersede them.

### 8.1 `variant_extraction.py`

| Change | Why |
|---|---|
| `load_mutations()`; `--mutations` takes **multiple files** | Panels can be passed directly, no pre-concatenation step. |
| `parse_mutation_file()` auto-detects a **headered panel TSV** (`chrom pos vtype gene label source expected_uids ref_aa codon notes`) alongside the historical whitespace format; parses `end=` / `edges=` out of `notes` | The module previously crashed on the panel header (`ValueError: invalid literal for int(): 'pos'`) and could not read the ITD region or exon edges. |
| `load_exon_edges()` / `_symbol_to_ensg()` — exon boundaries loaded by default from the bundled `components/long_read/resources/Hs_Ensembl_exon.txt.gz`; `--exon_annot` overrides, `--gff` is a fallback | The FLT3-ITD soft-clip caller silently ran with an **empty** splice blacklist unless `--gff` was passed, removing its only specificity control. Deleted the hard-coded hg38 literal `28033500 <= b <= 28034700`. Verified: the bundled file reproduces the panel-derived FLT3 edges exactly (E15.1 28033887–28033991, E14.1 28034082–28034214, E13.1 28034301–28034407). |
| Exon edges resolved **per gene** (`gene_edges`), and per-locus `edges` from the panel take precedence | FLT3's edges were previously applied to every ITD locus regardless of gene. |
| `run_from_metadata()` + `--metadata` / `--samples`; writes `input_manifest.tsv` (`uid`, `bam`, `n_alt_bams`, `seconds`) | Whole-cohort scan without a wrapper script, and a record of which BAM produced which output. |
| Writes `<sample>_variant_readcounts.tsv` (18 columns: `uid chrom pos gene label vtype source expected depth ref_base ref_reads alt_base alt_reads other_reads vaf n_cells_mut n_cells_wt note`) | Only a human-readable `.txt` summary existed; nothing machine-readable for cross-sample aggregation. Column set is a contract — do not reorder. |
| Restored `other_reads` and the `note` evidence field (SNV FASTA-vs-reads mismatch, indel allele spectrum, INS dominant-length fraction, ITD mechanism counts) | These QC fields were computed by the primitives and discarded. |
| Each distinct locus scanned **once**; one summary row per input row; per-cell rows emitted once per label | A panel lists a variant once per expecting donor (153 rows / 146 labels / 139 loci). Cells were being double- and triple-counted (ALK 70 → 210, ATRX 1487 → 2974). |
| Indel `note` rendered from a sorted dict | Dict order followed hash-seeded set iteration, so the field was not byte-reproducible. |

### 8.2 `variant_scan.py`

| Change | Why |
|---|---|
| New `expected_indel_class(label)` → `ins` / `del` / `any` / `None`, parsed from the panel HGVS (`c.487_488ins`→ins, `c.377delA`→del, `delins`/`p.…fs`→any) | A supervised scan already knows the variant class. |
| `scan_indel(..., expected=...)`: a read is mutant if it carries an indel **of the expected class** within the window. Cell = MUT if it has one; WT if it has only clean spanning reads; **not reported** if it carries only some other indel class | The previous "two-major-allele" model held a popularity contest and called a cell mutant only if it carried the winning allele. At `RUNX1_c.487_488ins` with 7 insertion vs 7 deletion reads, cells carrying the **insertion** — the variant under test — were dropped, and some were labelled wild-type. Example: barcode `CGCTTGGGAAGTAACA` (`0018_Af_N1c`) had `mut_reads=1 wt_reads=1` → MUT; under the vote it became `mut_reads=0` → **WT**, while carrying an insertion. |
| Major-allele vote retained only as the fallback for labels with no indel class, with an explicit `(-count, type, length)` tie-break | `inds` is a `set` of `(str, int)` tuples; string hashing is seeded per process → set iteration order varies → `Counter` insertion order varies → `most_common()` broke ties differently **between runs of the same code on the same BAM**. |

### 8.3 Testing performed

| Test | Result |
|---|---|
| Barcode-level validation, KINNEX (39/39) | MUT **142,056/142,056**, WT **3,852,135/3,852,135** confirmed from the BAM, **0 miscalls** |
| Barcode-level validation, BEAT-AML (8/8, `--bulk`) | MUT **510/510**, WT **568/568** confirmed, **0 miscalls** |
| Determinism | two independent runs, all 9 readcount fields incl. `note`: **0 differences** |
| Reference genotypes reproduced | `5801M_pre` SRSF2 P95R **2053 MUT / 1286 WT**; RUNX1 W279\* **256 / 435** |
| FLT3-ITD with no flags | 83 FLT3 exon boundaries auto-loaded; identical call to the panel-supplied edges (`ITD9bp`, depth 1299, alt 3, vaf 0.0023) |
| Divergence from the superseded `results/` | 229 barcodes dropped — **0** carry the expected class, **219 confirmed prior-run errors**, 10 unverifiable; 10 read-count changes, all increases, all matching the BAM |

Validation is **independent**: `scripts/validate_barcode_calls.py` (in the analysis folder) re-opens
each BAM and parses CIGARs and base composition itself, importing neither `variant_scan` nor
`variant_extraction`. Comparing against the old `results/` would prove nothing, since it was produced
by the major-allele rule this change rejects.

### 8.4 Not covered

* `STAG2_c.516A>C` (29,856 barcodes) and `ASXL1_c.1888` (484) are typed `Indel` in the panel but their
  labels state no indel class, so no expected allele can be derived. They fall back to the major-allele
  vote and are reported as `skipped_no_expected_class` — **not** as validated. `STAG2_c.516A>C` is
  written as a substitution; the panel entry looks wrong at source.
* Panel construction, de-novo CDS scanning and cross-sample calling (`build_panels.py`,
  `denovo_gene_scan.py`, `aggregate.py`, `aggregate_barcodes.py`) remain outside the installed package.

Run of record, outputs and full command lines:
`/Volumes/salomonis-archive/LabFiles/Nathan/Revio/MDS-AML-KINNEX-5/analysis/genomic_variants/variant_discovery/KINNEX_supervised_variant_detection/README_variant_extraction.md`
