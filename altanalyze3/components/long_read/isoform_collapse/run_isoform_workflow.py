#!/usr/bin/env python3
"""End-to-end long-read isoform-collapse workflow (as a user would run it).

From the GENCODE reference + per-sample long-read outputs, this produces:
  1. an ENST reference cache (from gff_process output -- NOT re-annotated),
  2. an integrated cross-sample isoform collapse with ENST injection (combined, not per-sample),
  3. quantified per-sample isoform h5ads + the FINAL catalog/report (reads-desc within each gene),
  4. protein prediction (if a genome FASTA is available).

Edit the CONFIG block for your data, then:  python run_isoform_workflow.py
"""

import os
import sys
import time

# Make the altanalyze3 package importable when run as a plain script.
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from altanalyze3.components.long_read import gff_process
from altanalyze3.components.long_read.isoform_collapse import reference as R
from altanalyze3.components.long_read.isoform_collapse import pipeline as P


# ----------------------------------------------------------------- CONFIG ----
# Ensembl exon reference DB (defines the fixed exon-token namespace shared by ALL gff_process runs).
ENSEMBL_EXON_DB = "/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt"
GENE_SYMBOLS   = "/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl-annotations.txt"

# GENCODE/Ensembl reference annotation (run separately through gff_process, like any other GFF).
REFERENCE_GFF  = "/Users/saljh8/Dropbox/Revio/BAMs/iPSC/gencode.v49.annotation.gtf.gz"

# Per-sample inputs: (name, raw_per_read_h5ad, transcript_associations.txt, raw_gff.gz)
#   raw_per_read_h5ad: the UNCOLLAPSED molecule-level h5ad (var_names = "gene:molecule",
#   millions of vars) -- NOT the "-isoform.h5ad" output of a prior collapse.
SAMPLES = [
    ("EP_SRSF2_Ctrl",
     "/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Ctrl/EP_SRSF2_Ctrl.h5ad",
     "/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Ctrl/gff-output/transcript_associations.txt",
     "/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Ctrl/EP_SRSF2_Ctrl.gff.gz"),
    ("EP_SRSF2_Aza",
     "/Users/saljh8/Dropbox/Revio/BAMs/iPSC/AZA/EP_SFSF2_AZA.h5ad",
     "/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Aza/gff-output/transcript_associations.txt",
     "/Users/saljh8/Dropbox/Revio/BAMs/iPSC/AZA/EP_SFSF2_AZA.gff.gz"),
]

OUTDIR       = "/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Ctrl/gff-output"
GENOME_FASTA = None   # set to the genome .fa for protein prediction; None => protein step is skipped
NPROC        = 4
MIN_TOTAL    = 3      # drop integrated isoforms with <= 2 total reads (keep >= 3)


def log(msg):
    print(f"{time.strftime('%H:%M:%S')}  {msg}", flush=True)


def write_report(outdir, symbols_file, log=log):
    """FINAL_isoform_report.txt: one row per final isoform, isoforms sorted reads-desc within gene."""
    import collections
    cat = os.path.join(outdir, "FINAL_isoform_catalog.tsv")
    mp  = os.path.join(outdir, "FINAL_structure_to_exemplar.tsv")
    out = os.path.join(outdir, "FINAL_isoform_report.txt")

    sym = {}
    if os.path.exists(symbols_file):
        with open(symbols_file) as f:
            for line in f:
                t = line.rstrip("\n").split("\t")
                if len(t) >= 2:
                    sym[t[0]] = t[1]

    # exemplar -> its longest mapped structure (the representative full-length form)
    exemplar_struct = {}
    with open(mp) as f:
        next(f, None)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            _g, struct, ex = p
            cur = exemplar_struct.get(ex)
            if cur is None or len(struct) > len(cur):
                exemplar_struct[ex] = struct

    by_gene = collections.defaultdict(list)
    with open(cat) as f:
        next(f, None)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            g, ex, _blk, reads, _bin, known = p
            by_gene[g].append((ex, int(reads), known == "known"))

    nk = nn = rk = rn = 0
    with open(out, "w") as o:
        o.write("isoform_id\tgene\tgene_symbol\ttotal_reads\tsource\tstructure\n")
        for g in sorted(by_gene):
            for ex, reads, kn in sorted(by_gene[g], key=lambda x: -x[1]):
                src = "Ensembl" if kn else "novel"
                o.write(f"{ex}\t{g}\t{sym.get(g,'NA')}\t{reads}\t{src}\t{exemplar_struct.get(ex,'')}\n")
                if kn: nk += 1; rk += reads
                else:  nn += 1; rn += reads
    total = nk + nn
    log(f"[report] {out}")
    log(f"[report] total={total:,}  Ensembl={nk:,} ({100*nk/total:.1f}%, {rk:,} reads)  "
        f"novel={nn:,} ({100*nn/total:.1f}%, {rn:,} reads)")
    return out


def main():
    t0 = time.time()
    os.makedirs(OUTDIR, exist_ok=True)

    # --- STEP 1: GENCODE reference -> gff_process transcript_associations (run separately) -------
    ref_ta = os.path.join(os.path.dirname(REFERENCE_GFF), "gff-output", "transcript_associations.txt")
    if not os.path.exists(ref_ta):
        log(f"[gff_process] annotating reference {os.path.basename(REFERENCE_GFF)} ...")
        ref_ta = gff_process.consolidateLongReadGFFs(REFERENCE_GFF, ENSEMBL_EXON_DB, mode="Ensembl")
    log(f"[gff_process] reference transcript_associations: {ref_ta}")

    # --- STEP 2: ENST cache from gff_process output (filter ENST rows, strip version) ------------
    enst_cache = os.path.join(OUTDIR, "ENST_reference_structures.tsv")
    ref = R.annotate_reference(ref_gff=None, ref_exons=None, cache_path=enst_cache,
                               force=True, transcript_associations=ref_ta, log=log)
    log(f"[reference] ENST cache: {sum(len(v) for v in ref.values()):,} structures, "
        f"{len(ref):,} genes -> {enst_cache}")

    # --- STEP 3: integrated collapse (ENST injected on COMBINED set) -> catalog + quant h5ads ----
    #            + STEP 4: protein prediction (only if GENOME_FASTA is provided) -----------------
    result = P.run_pipeline(
        samples=SAMPLES, outdir=OUTDIR, nproc=NPROC, min_total=MIN_TOTAL,
        write_h5ad=True, genome_fasta=GENOME_FASTA, ref_gff=REFERENCE_GFF,
        enst_cache=enst_cache, log=log,
    )

    # --- STEP 5: final isoform report (reads-desc within each gene) ------------------------------
    write_report(OUTDIR, GENE_SYMBOLS, log=log)

    if result.get("protein"):
        log(f"[protein] protein outputs: {result['protein']}")
    elif GENOME_FASTA is None:
        log("[protein] SKIPPED (GENOME_FASTA is None -- set it to enable protein prediction)")

    log(f"[done] full workflow in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
