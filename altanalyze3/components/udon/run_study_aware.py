#!/usr/bin/env python3
"""Hard-coded study-aware UDON workflow -> integrated, batch-removed programs + SATAY-UDON.

Runs, in order, on matched-control (sex+platform) pseudobulk folds:
  1. study_aware_udon       per-study UDON (feature selection -> NMF -> MarkerFinder) on the
                            gene-filtered combined folds; saves per-study clusters/markers/heatmaps
  2. study_aware_genelist   conserved vs unique programs by top-50 marker gene-list overlap
  3. study_aware_integrate  SVA integration of conserved centroids -> final programs; re-SVA
                            classifies EVERY pseudobulk -> final_program_assignments.tsv
  4. final_heatmap          MarkerFinder heatmap of the final batch-removed programs (+ GO callouts,
                            donor/cell-type composition heatmaps)
  5. satay_udon             SATAY-UDON donor-level covariate enrichment of the final programs

Env config (pseudobulk_protocol.udon_restriction):
  UDON_GENE_FILTER  default on: protein-coding only (species annotation), drop RPL/RPS ribosomal
                    + XIST/TSIX, KEEP Y-chromosome genes (loss-of-Y signal)
  UDON_CELL_TYPE    restrict the whole workflow to one cell type
  UDON_TAG          output-folder suffix (keep a gene-filtered run separate from an old one)
  UDON_SPECIES      default Hs
Outputs -> UDON/study_aware[<_celltype>][<_tag>]/.
"""
import os, sys, subprocess

UDON_DIR = os.path.dirname(os.path.abspath(__file__))
STEPS = ["study_aware_udon.py", "study_aware_genelist.py", "study_aware_integrate.py",
         "final_heatmap.py", "satay_udon.py"]


def main():
    import pseudobulk_protocol as P
    ct, gf, sp, sfx = P.udon_restriction()
    print(f"=== study-aware workflow | cell_type={ct or 'ALL'} | gene_filter={gf} | "
          f"species={sp} | output suffix='{sfx}' ===", flush=True)
    for s in STEPS:
        print(f"\n{'='*72}\n=== {s} ===\n{'='*72}", flush=True)
        r = subprocess.run([sys.executable, os.path.join(UDON_DIR, s)], env=os.environ)
        if r.returncode != 0:
            print(f"!!! {s} exited with {r.returncode}; stopping workflow.")
            sys.exit(r.returncode)
    print("\n=== study-aware workflow complete ===")


if __name__ == "__main__":
    sys.path.insert(0, UDON_DIR)
    main()
