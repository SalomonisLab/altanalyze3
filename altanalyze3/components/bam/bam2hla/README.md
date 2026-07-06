# bam2hla — fast class-I HLA genotyping from a genome BAM (pysam pileup)

Calls **HLA-A / -B / -C at 2-field (4-digit) resolution** directly from a
genome-aligned BAM, with **no read re-alignment and no OptiType**. Reads that map
to the HLA loci on chr6 already carry the sample's alleles; we read them off with
`pysam` pileup at a precomputed set of allele-defining positions and solve for the
diploid genotype. Works for **bulk RNA-seq, WGS/WES, and single-cell long-read
(PacBio Iso-Seq / KINNEX)** BAMs, in **hg19 or hg38** (auto-detected). Written in
the pileup style of `components/bam/variant_extraction.py`.

## Why this instead of OptiType
OptiType/arcasHLA extract reads and *re-align* them to the whole IMGT/HLA allele
set (razers3/kallisto) then solve an ILP — heavy, and needs FASTQ + a container.
Here the reads are already aligned; we only need to know **what each IMGT allele
looks like at each genomic position**, which we precompute once. Runtime is a few
seconds of pileup + a vectorized solve.

## Method

### 1. Reference build (`build_reference.py`, one-time, self-validating)
For each gene we project every IMGT/HLA allele's **exon-2/3/4** sequence (the
antigen-recognition domain that defines the 2-field type) onto chr6 coordinates:

1. Reconstruct the reference transcript CDS from a GTF + chr6 FASTA (strand-aware),
   giving a cDNA-index → genomic-coordinate map.
2. Find the single IMGT allele **identical** to the genome's exon-2/3/4 sequence —
   this is the allele the reference genome carries. An exact match *proves* the
   coordinate mapping; recovered genome alleles are **A\*03:01:01:01 / B\*07:02:01:01
   / C\*07:02:01:01** for both GRCh38 and GRCh37, the known reference alleles.
3. Walk the IMGT alignment columns; each column where the genome allele has a base
   → one genomic position with its reference base and every allele's expected
   forward-strand base. A hard gate requires the projected reference base to equal
   the FASTA base at every position.

Output: `build/signatures_{hg38,hg19}.json.gz` (822 exon-2/3/4 positions/gene;
~680–693 polymorphic; ~5.5–7k 2-field alleles/gene).

### 2. Runtime typing (`bam2hla.py`)
1. **Auto-detect build** from the chr6 contig length (hg38 170,805,979 vs hg19
   171,115,067), confirmed by requiring read coverage at the HLA loci for that
   build (the two builds place the genes ~300 kb apart, so only the right one has
   coverage). Handles both `chr6` and `6` naming.
2. **Pileup** the signature positions (spliced/`N`-CIGAR aware, so RNA works).
3. **Denoise** each position to its *real base-set* = major base + any minor base
   ≥ 15% (and ≥3 reads). This removes sequencing-error reads so zygosity is decided
   by real signal, not noise.
4. **Fit**: choose the allele pair whose predicted base-set per position
   `{e_X, e_Y}` matches the observed real base-sets at the most positions
   (numpy-bitmask vectorized). Homozygous vs heterozygous falls out naturally
   (a hom position wants one base, a het position two).
5. **Tie-break** equally-fitting solutions by the 21-population allele-frequency
   prior, then completeness, then prefer homozygous — this stops rare partially
   typed alleles from beating the true common one.

## Usage
```bash
PY=...GitHub/altanalyze3/.venv/bin/python   # needs pysam + numpy

# one-time reference build (already run; artifacts in build/)
$PY build_reference.py

# type a BAM (build auto-detected)
$PY bam2hla.py --bam sample.bam --out sample_hla.tsv
# -> [hla] HLA-A*02:01,HLA-A*02:05,HLA-B*50:01,HLA-B*18:01,HLA-C*06:02,HLA-C*05:01

# benchmark vs OptiType (fast local disk, e.g. cluster)
$PY benchmark.py run --build hg38 --bam-dir <dir> --truth <optitype.txt>
```

## Files
| file | purpose |
|---|---|
| `bam2hla.py` | runtime typer: build detection, pileup, set-mismatch diploid solver, CLI |
| `imgt_parser.py` | parse IPD-IMGT/HLA `*_nuc.txt` alignments → per-allele exon sequences |
| `build_reference.py` | anchor IMGT alleles to chr6 coords → `build/signatures_*.json.gz` |
| `benchmark.py` | OptiType concordance: `cache` (resilient count caching) / `eval` / `run` |
| `validate_subset.py` | ready-to-run check on 7 diverse samples (embedded OptiType truth) |
| `run_cluster_benchmark.lsf` | bsub script for the full 362-sample benchmark on the cluster |
| `build/signatures_*.json.gz` | genomic allele-signature DBs (hg38, hg19) |
| `data/imgt/…`, `data/genome/…`, `data/HLA_Allele_frequency_21_populations.csv` | inputs (see Provenance) |

## Provenance / versions
- IPD-IMGT/HLA **3.64.0** (`ANHIG/IMGTHLA/Latest`): `A_nuc.txt`, `B_nuc.txt`, `C_nuc.txt`, `hla_nuc.fasta`
- Exon coords: Ensembl **104** (GRCh38, CCDS transcripts ENST00000396634 / ENST00000412585 / ENST00000376228) and Ensembl **87** (GRCh37)
- chr6 FASTA: UCSC hg38 / hg19 (chr6 sequence is identical across GRCh38/GRCh37 sources, so Ensembl coords apply)
- Frequency prior: `HLA_Allele_frequency_21_populations.csv` (from `components/snaf`)

## Validation status
- **293T (HEK293T), hg38** — bam2hla: `A*02:01/A*02:01, B*07:02/B*07:02, C*07:02/C*07:02`
  → **exact match** to OptiType (all 6/6 alleles).
- **MK-Platelet (iPSC, long-read Iso-Seq, hg38)** — `A*02:01/A*02:05, B*50:01/B*18:01,
  C*06:02/C*05:01`. No independent truth, but all calls are common alleles and the
  B\*18:01~C\*05:01 and B\*50:01~C\*06:02 pairings are known linked haplotypes
  (internal consistency).

## Notes / limitations
- **Class I only** (HLA-A/-B/-C), matching OptiType. Class II (DRB1/DQB1/DPB1) is a
  straightforward extension (add those `*_nuc.txt` + gene coords) but not built here.
- Resolution is **2-field**, matching the OptiType aggregate.
- Requires expressed/covered loci; RNA is fine because exon 2/3 (the ARD) is
  covered. Strong allelic-expression imbalance in RNA can mask a true heterozygote
  (a limitation shared with any RNA-based caller, incl. OptiType).
