#!/bin/bash
# =====================================================================
# LSF parallel launcher for the AltAnalyze3 long-read single-cell workflow
# (sclr). An ALTERNATIVE to running MDS.py-style drivers in one process:
# fans out the per-sample step (phase 1) as one bsub per sample, then runs
# the cross-sample integration steps (phases 2-4) once. Same analyses, same
# outputs as the single-process driver.
#
#   phase 1  sclr            BAM -> molecule/junction h5ad (+ cellHarmony)   [PARALLEL, 1 job/sample]
#   phase 2  sclr-junctions  junction aggregation + PSI + diff splicing      [1 integration job]
#   phase 3  sclr-isoforms   two-tier isoform collapse + isoform h5ads        [1 integration job]
#   phase 4  sclr-diff       differential isoform/junction between groups      [1 integration job]
#
# No install needed: all deps already exist on the cluster, so we set
# PYTHONPATH to the source checkout and invoke `python3 -m altanalyze3`.
# Per-sample BAM directory and reverse-complement are read from each metadata
# row, so nothing about either appears on the command line.
#
#   bash submit_longread_cluster.sh
# =====================================================================

# ---- HARD-CODED PATHS (edit these) ----
METADATA="/data/salomonis-archive/LabFiles/Nathan/Revio/MDS-AML-KINNEX-2/mds_metadata_bam_revised.txt"
SPECIES="human"                       # human | mouse
ALTANALYZE3_SRC="/data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3"

# Reference databases (from MDS.py). exon_annot/gene_symbol are OPTIONAL: omit to use the
# species defaults bundled with the package. ref_gff + genome_fasta are required for the collapse.
EXON_ANNOT="/data/salomonis-archive/LabFiles/Nathan/Revio/Hs/Ensembl91/Hs_Ensembl_exon.txt"
GENE_SYMBOL="/data/salomonis-archive/LabFiles/Nathan/Revio/Hs/Ensembl91/Hs_Ensembl-annotations.txt"
GENCODE_GFF="/data/salomonis-archive/LabFiles/Nathan/Revio/Hs/gencode.v45.annotation.gff3"
GENOME_FASTA="/data/salomonis-archive/LabFiles/Nathan/Revio/Hs/genome.fa"

# Cell cluster labels -- TWO mutually exclusive options (set ONE, leave the other empty):
#   CELL_ANNOT      : a pre-existing barcode->cluster file (cellHarmony format, "barcode-1.sample<TAB>cluster").
#                     Used as-is (no alignment). This is the combined reference built for this metadata.
#   CELLHARMONY_REF : a registry id (e.g. hs_bm_reference) or centroid .txt -- RE-RUNS cellHarmony alignment.
CELL_ANNOT="/data/salomonis-archive/LabFiles/Nathan/Revio/MDS-AML-KINNEX-2/combined_cellHarmony_unique_barcodes.txt"
CELLHARMONY_REF=""
# Isoform collapse method: wta (winner-takes-all, default) or em.
COLLAPSE_METHOD="wta"
# Differential comparisons: "groupA,groupB" pairs (groups from the metadata 'groups' column),
# semicolon-separated for multiple. Leave empty to skip phase 4.
CONDITIONS="young,aged"
ANALYSES="junction,isoform"
DIFF_METHOD="limma"                     # mwu | limma

PYTHON_MODULE="python3/3.11.3"          # `module load` target on the cluster

# ---- DIRECTORIES ----
ROOT=$PWD
LOG_DIR="${ROOT}/logs"
mkdir -p "$LOG_DIR"

# Invocation: no install, run the package from source.
AA3="env PYTHONPATH=${ALTANALYZE3_SRC} python3 -m altanalyze3"

# ---- VALIDATION ----
if [[ ! -f "$METADATA" ]]; then
    echo "ERROR: metadata file not found:"; echo "       $METADATA"; exit 1
fi
if [[ ! -d "$ALTANALYZE3_SRC/altanalyze3" ]]; then
    echo "ERROR: altanalyze3 source dir not found (expected an 'altanalyze3' package inside):"
    echo "       $ALTANALYZE3_SRC"; exit 1
fi
if [[ ! -f "$GENCODE_GFF" ]]; then
    echo "ERROR: reference GFF not found:"; echo "       $GENCODE_GFF"; exit 1
fi
if [[ ! -f "$GENOME_FASTA" ]]; then
    echo "ERROR: genome FASTA not found:"; echo "       $GENOME_FASTA"; exit 1
fi
# Cluster labels: exactly ONE of CELL_ANNOT / CELLHARMONY_REF must be set.
if [[ -n "$CELL_ANNOT" && -n "$CELLHARMONY_REF" ]]; then
    echo "ERROR: set only ONE of CELL_ANNOT or CELLHARMONY_REF (they are mutually exclusive)."; exit 1
fi
if [[ -z "$CELL_ANNOT" && -z "$CELLHARMONY_REF" ]]; then
    echo "ERROR: set either CELL_ANNOT (existing labels) or CELLHARMONY_REF (re-run alignment)."; exit 1
fi
if [[ -n "$CELL_ANNOT" && ! -f "$CELL_ANNOT" ]]; then
    echo "ERROR: CELL_ANNOT file not found:"; echo "       $CELL_ANNOT"; exit 1
fi

# Resolve the cluster-labeling flags.
#  - CLUSTER_OPT      : for phase 1 (sclr accepts --cell_annot OR --cellHarmony_ref).
#  - INT_CLUSTER_OPT  : for phases 2-4 (integration accepts ONLY --cell_annot). When re-aligning
#                       (CELLHARMONY_REF), the integration steps auto-discover the per-sample
#                       barcode->cluster TSVs that phase 1 wrote, so no flag is passed.
if [[ -n "$CELL_ANNOT" ]]; then
    CLUSTER_OPT="--cell_annot $CELL_ANNOT"
    INT_CLUSTER_OPT="--cell_annot $CELL_ANNOT"
    CLUSTER_DESC="$CELL_ANNOT (pre-existing labels)"
else
    CLUSTER_OPT="--cellHarmony_ref $CELLHARMONY_REF"
    INT_CLUSTER_OPT=""
    CLUSTER_DESC="$CELLHARMONY_REF (re-run alignment; integration auto-discovers per-sample TSVs)"
fi

echo "Metadata:          $METADATA"
echo "Species:           $SPECIES"
echo "AltAnalyze3 src:   $ALTANALYZE3_SRC"
echo "Reference GFF:     $GENCODE_GFF"
echo "Genome FASTA:      $GENOME_FASTA"
echo "Cell labels:       $CLUSTER_DESC"
echo "Collapse method:   $COLLAPSE_METHOD"
echo "Logs:              $LOG_DIR"

# Optional exon/gene-symbol overrides (only added to the commands if the files exist).
EXON_OPT=""; [[ -f "$EXON_ANNOT" ]] && EXON_OPT="--exon_annot $EXON_ANNOT"
GENE_OPT=""; [[ -f "$GENE_SYMBOL" ]] && GENE_OPT="--gene_symbol $GENE_SYMBOL"

# =====================================================================
# PHASE 1: per-sample, one bsub per uid (PARALLEL).
# A uid with multiple BAM rows (hashed/merged libraries) is processed as a
# UNIT by a single job via --sample <uid>; the metadata row(s) carry each
# BAM path and the per-sample reverse-complement flag.
# =====================================================================
JOB_ID=0
DEP=""
# Unique uids = column 1, skip header, strip \r, drop blanks.
# NOTE: do NOT name the loop variable UID -- it is a read-only bash builtin (the user id), which
# silently makes the loop run zero times. Use SAMPLE_UID.
for SAMPLE_UID in $(awk -F'\t' 'NR>1 && $1!="" {gsub(/\r/,"",$1); print $1}' "$METADATA" | sort -u); do
    JOB_ID=$((JOB_ID+1))
    JNAME="sclr_${SAMPLE_UID}"
    bsub -W 24:00 -n 12 -M 64000 \
         -R "span[hosts=1]" \
         -J "$JNAME" \
         -o "${LOG_DIR}/${SAMPLE_UID}.phase1.out" \
         -e "${LOG_DIR}/${SAMPLE_UID}.phase1.err" \
         "module load ${PYTHON_MODULE}; \
          ${AA3} sclr \
            --metadata \"$METADATA\" --sample \"$SAMPLE_UID\" \
            --ref_gff \"$GENCODE_GFF\" --species \"$SPECIES\" \
            ${CLUSTER_OPT} \
            ${EXON_OPT} ${GENE_OPT} --cpus 12"
    if [[ -z "$DEP" ]]; then DEP="done($JNAME)"; else DEP="$DEP && done($JNAME)"; fi
done
echo "Submitted $JOB_ID per-sample (phase 1) jobs."

# Guard: if no phase-1 jobs were submitted, DEP is empty and the integration jobs would run with
# `-w ''` against non-existent per-sample outputs. Abort instead of silently launching them.
if [[ "$JOB_ID" -eq 0 || -z "$DEP" ]]; then
    echo "ERROR: 0 phase-1 jobs submitted (no uids parsed from $METADATA). Not launching integration."
    exit 1
fi

# =====================================================================
# PHASE 2-4: integration jobs, each gated on the previous finishing.
# Phase 2 waits for ALL phase-1 jobs (the -w dependency below). For a MANUAL
# two-phase run, delete the `-w "..."` flags and submit these after the
# sclr_* jobs clear.
# =====================================================================
bsub -W 24:00 -n 16 -M 196000 -R "span[hosts=1]" \
     -J sclr_junctions -w "$DEP" \
     -o "${LOG_DIR}/phase2_junctions.out" -e "${LOG_DIR}/phase2_junctions.err" \
     "module load ${PYTHON_MODULE}; \
      ${AA3} sclr-junctions --metadata \"$METADATA\" --species \"$SPECIES\" ${EXON_OPT} ${INT_CLUSTER_OPT}"

bsub -W 24:00 -n 16 -M 196000 -R "span[hosts=1]" \
     -J sclr_isoforms -w "done(sclr_junctions)" \
     -o "${LOG_DIR}/phase3_isoforms.out" -e "${LOG_DIR}/phase3_isoforms.err" \
     "module load ${PYTHON_MODULE}; \
      ${AA3} sclr-isoforms --metadata \"$METADATA\" --ref_gff \"$GENCODE_GFF\" \
        --genome_fasta \"$GENOME_FASTA\" --collapse_method \"$COLLAPSE_METHOD\" \
        --species \"$SPECIES\" ${EXON_OPT} ${INT_CLUSTER_OPT}"

if [[ -n "$CONDITIONS" ]]; then
    bsub -W 12:00 -n 8 -M 128000 -R "span[hosts=1]" \
         -J sclr_diff -w "done(sclr_isoforms)" \
         -o "${LOG_DIR}/phase4_diff.out" -e "${LOG_DIR}/phase4_diff.err" \
         "module load ${PYTHON_MODULE}; \
          ${AA3} sclr-diff --metadata \"$METADATA\" --conditions \"$CONDITIONS\" \
            --analyses \"$ANALYSES\" --method \"$DIFF_METHOD\" --species \"$SPECIES\" ${GENE_OPT} ${INT_CLUSTER_OPT}"
fi

echo "Submitted integration chain (phase 2 -> 3 -> 4)."
echo "Phase 1 runs in parallel; phases 2-4 start automatically as each dependency completes."

# NOTE: per-sample reverse-complement is NOT a command-line flag -- the 'reverse' (TRUE/FALSE)
# column of the metadata is read per row and applied at junction/isoform matrix build, so it is
# honored automatically for every sample that indicates it.

# To debug submission (echo each command + variable as it runs, e.g. when a job fails to launch):
#   bash -x submit_longread_cluster.sh
