"""Standard cutoffs for single-cell variant evidence (known variants and their inferred carrier cells).

One place for every threshold that decides whether a variant, a marker or an inferred cell is reported.
`variant_impact.py` and `mt_variants.py` take their argparse defaults from this file, so a run with no threshold
flags reproduces the run of record (MDS-AML-KINNEX-5, 2026-07 to 2026-11). Each entry names its source.

Print the table:  python -m altanalyze3.components.bam.variant_evidence_defaults [--tsv OUT.tsv]
"""
from __future__ import annotations

# name: (value, applies to, source)
CUTOFFS = {
    # ---- passenger-variant markers (variant_impact.py --passenger-association, mt_variants.py --nominate-fdr)
    "PASSENGER_MAX_Q": (0.05, "nuclear and mtDNA passenger markers: Benjamini-Hochberg q below this",
                        "module default; July 2026 runs"),
    "PASSENGER_MIN_MUT_CELLS": (10, "nuclear and mtDNA passenger markers: minimum driver-mutant cells for a unit to be tested",
                                "module default; July 2026 runs"),
    "NUCLEAR_PASSENGER_MAX_BACKGROUND": (0.15, "nuclear passenger marker: maximum carrier fraction among driver-wild-type cells",
                                         "July 2026 results of record: 50 of 50 markers are at or below 0.1496. The July command line says 0.5, "
                                         "but the July software applied 0.15; a default-0.5 rerun gives 93 markers (step 38, 2026-11-14)"),
    "MT_PASSENGER_MAX_BACKGROUND": (0.15, "mtDNA passenger marker (--nominate-fdr): maximum carrier fraction among driver-wild-type cells",
                                    "July 2026 run of record (scripts/nominate_bam_matrices.sh); the CLI default was 0.10"),
    "PASSENGER_MIN_MAPQ": (1, "reads counted when the passenger matrices are built", "module behaviour (None resolves to 1)"),
    "SELF_MARKER_MAX_DISTANCE_BP": (10, "a nuclear marker on the driver's chromosome within this distance of the driver position "
                                        "is the driver itself and is not a passenger",
                                    "MDS-AML-KINNEX-5 2026-11-13: the NPM1c insertion was nominated as its own marker in 8 hit rows"),
    # ---- expression classifiers (variant_impact.py gene and isoform level)
    "CLASSIFIER_MIN_CELLS_PER_CLASS": (15, "minimum mutant and minimum wild-type cells in a cell state to train a model",
                                       "July 2026 gene run and November 2026 isoform run; the CLI default was 20"),
    "CLASSIFIER_TOP_MARKERS": (50, "marker features per class", "module default; both runs"),
    "CLASSIFIER_MIN_ACCURACY": (0.85, "a per-cell-state model deploys calls at or above this resubstitution accuracy",
                                "variant_impact.ACCURACY_MIN"),
    "DEG_CONCORDANCE_MIN_PEARSON": (0.6, "gene level: Pearson r over the DEG-in-both overlap, above this",
                                    "Variant-Impact benchmark (deg_concordance), July 2026"),
    # ---- LOH classifier
    "LOH_MIN_PROBABILITY": (0.5, "a cell is called LOH at or above this classifier probability", "MDS-AML-KINNEX-5 step 04"),
    "LOH_RETENTION_MAX_P": (0.05, "LOH cells count as evidence for a variant only when predicted-LOH cells carry the variant more "
                                  "often than predicted-heterozygous cells (one-sided Fisher exact p below this)",
                            "MDS-AML-KINNEX-5 2026-11-14: 5801 chr21 p 1e-225 passes, AM75 chr4 p 0.43 fails"),
    # ---- CNV clone association (CNV-mutation/method_comparison/code/cnv_method_comparison.py)
    "CNV_MIN_MUT_CELLS": (15, "minimum mutant cells among the clone cells for a variant to be tested", "July 2026 MIN_MUT"),
    "CNV_MIN_MUT_CELLS_IN_CLONE": (3, "minimum mutant cells in a clone for the clone to be a candidate", "July 2026 MIN_MUT_IN_CLONE"),
    "CNV_MAX_FDR": (0.05, "best-clone one-sided Fisher exact, Bonferroni over clones, then BH; below this", "July 2026"),
    # ---- germline-like variants (excluded from the somatic variant list)
    "GERMLINE_MIN_PATIENT_VAF": (0.95, "patient-wide alt read fraction at or above this", "MDS-AML-KINNEX-5 2026-11-14"),
    "GERMLINE_MIN_PATIENT_READS": (100, "minimum patient-wide reads over the position", "MDS-AML-KINNEX-5 2026-11-14"),
    "GERMLINE_MIN_LYMPHOID_VAF": (0.95, "alt read fraction in mature lymphoid cells at or above this", "MDS-AML-KINNEX-5 2026-11-14"),
    "GERMLINE_LYMPHOID_NULL_FRACTION": (0.5, "binomial null: at most this alt fraction in lymphoid cells", "MDS-AML-KINNEX-5 2026-11-14"),
    "GERMLINE_LYMPHOID_MAX_P": (0.05, "one-sided binomial p below this rejects the null", "MDS-AML-KINNEX-5 2026-11-14"),
    "GERMLINE_CLONE_ABSENT_MAX_FRACTION": (0.05, "a sample is clone-absent when every other variant of the patient with enough reads "
                                                 "is below this alt read fraction", "MDS-AML-KINNEX-5 2026-11-14 (AM75 on therapy)"),
    "GERMLINE_CLONE_ABSENT_MIN_READS": (50, "minimum reads over another variant for it to be assessed in a sample", "MDS-AML-KINNEX-5 2026-11-14"),
    # ---- clonal flag
    "CLONAL_MIN_FOLD": (4.0, "share of a variant among the alt reads of a sample changes at least this fold between two samples",
                        "Nathan Salomonis, 2026-11-13"),
    "CLONAL_MIN_VARIANT_ALT_READS": (50, "alt reads of the variant in at least one sample", "Nathan Salomonis, 2026-11-13"),
    "CLONAL_MIN_SAMPLE_ALT_READS": (100, "alt reads over the patient's variants for a sample to enter", "Nathan Salomonis, 2026-11-13"),
}

for _k, _v in CUTOFFS.items():
    globals()[_k] = _v[0]

# Mature lymphoid cell states of the marrow reference (Hs-MarrowAtlas-L3M). T, NK, mature B and plasma cells sit outside
# the myeloid clone in most MDS and AML patients, so a somatic clonal variant leaves reference reads in them.
MATURE_LYMPHOID_STATES = frozenset({
    "B Memory-1", "B Memory-2", "CD4 TCM", "CD8 TEM", "MAIT", "NK-Mature", "Plasma Cell", "T CD4 Activated",
    "T CD4 Naive-1", "T CD4 Naive-2", "T CD8 Naive", "Transitional-B-1", "Transitional-B-2"})


def binomial_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p)."""
    from scipy.stats import binom
    return float(binom.sf(k - 1, n, p))


def germline_call(alt: int, reads: int, lymphoid_alt: int, lymphoid_reads: int) -> str:
    """'germline_like', 'undetermined' or 'somatic_pattern' for one variant in one patient.

    germline_like  patient VAF and lymphoid VAF both at or above their minimum, and the lymphoid alt count rejects the
                   null of at most GERMLINE_LYMPHOID_NULL_FRACTION alt reads (one-sided binomial).
    undetermined   patient VAF qualifies but the lymphoid reads cannot reject or confirm (too few reads).
    somatic_pattern  everything else: patient VAF below the minimum, too few patient reads, or lymphoid cells with
                   reference reads.
    """
    if reads < GERMLINE_MIN_PATIENT_READS or alt / reads < GERMLINE_MIN_PATIENT_VAF:  # noqa: F821
        return "somatic_pattern"
    if lymphoid_reads == 0:
        return "undetermined"
    p = binomial_sf(lymphoid_alt, lymphoid_reads, GERMLINE_LYMPHOID_NULL_FRACTION)  # noqa: F821
    if lymphoid_alt / lymphoid_reads >= GERMLINE_MIN_LYMPHOID_VAF and p < GERMLINE_LYMPHOID_MAX_P:  # noqa: F821
        return "germline_like"
    if lymphoid_alt / lymphoid_reads >= GERMLINE_MIN_LYMPHOID_VAF:  # noqa: F821
        return "undetermined"
    return "somatic_pattern"


def inherited_by_clone_absent_sample(alt: int, reads: int, others) -> bool:
    """True when a variant stays in nearly every read of a sample from which the patient's other variants are gone.

    alt, reads   the candidate variant in ONE sample.
    others       [(alt, reads), ...] of the patient's other variants in the same sample; candidates that are themselves
                 near-fixed in the patient (alt fraction >= GERMLINE_MIN_PATIENT_VAF) must be left out by the caller.
    The sample is clone-absent when at least one other variant has >= GERMLINE_CLONE_ABSENT_MIN_READS reads and every such
    variant is below GERMLINE_CLONE_ABSENT_MAX_FRACTION. A variant at >= GERMLINE_MIN_LYMPHOID_VAF there (one-sided binomial
    p < GERMLINE_LYMPHOID_MAX_P against GERMLINE_LYMPHOID_NULL_FRACTION) is carried by the non-malignant cells: inherited.
    """
    seen = [(a, t) for a, t in others if t >= GERMLINE_CLONE_ABSENT_MIN_READS]  # noqa: F821
    if not seen or any(a / t >= GERMLINE_CLONE_ABSENT_MAX_FRACTION for a, t in seen) or reads == 0:  # noqa: F821
        return False
    return (alt / reads >= GERMLINE_MIN_LYMPHOID_VAF  # noqa: F821
            and binomial_sf(alt, reads, GERMLINE_LYMPHOID_NULL_FRACTION) < GERMLINE_LYMPHOID_MAX_P)  # noqa: F821


def _locus(name: str):
    """('5', 171410542) from 'NPM1_chr5:171410542', 'NPM1_chr5:171410539_C>CTCTG' or 'DDX10::SKA3_chr11:108715038'."""
    loc = name.rsplit("_chr", 1)[1]
    chrom, rest = loc.split(":", 1)
    return chrom, int(rest.split("_")[0])


def is_self_marker(driver: str, marker: str) -> bool:
    """True when a nuclear marker is the driver itself: same chromosome, within SELF_MARKER_MAX_DISTANCE_BP."""
    try:
        dc, dp = _locus(driver)
        mc, mp = _locus(marker)
    except (IndexError, ValueError):
        return False
    return dc == mc and abs(dp - mp) <= SELF_MARKER_MAX_DISTANCE_BP  # noqa: F821


def loh_retention_p(mut_in_loh: int, wt_in_loh: int, mut_in_het: int, wt_in_het: int) -> float:
    """One-sided Fisher exact p that predicted-LOH cells carry the variant more often than predicted-heterozygous cells."""
    from scipy.stats import fisher_exact
    return float(fisher_exact([[mut_in_loh, wt_in_loh], [mut_in_het, wt_in_het]], alternative="greater")[1])


def loh_is_evidence(mut_in_loh: int, wt_in_loh: int, mut_in_het: int, wt_in_het: int) -> bool:
    return loh_retention_p(mut_in_loh, wt_in_loh, mut_in_het, wt_in_het) < LOH_RETENTION_MAX_P  # noqa: F821


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="Print the standard variant-evidence cutoffs.")
    ap.add_argument("--tsv", default=None, help="also write the table to this path")
    a = ap.parse_args(argv)
    lines = ["cutoff\tvalue\tapplies_to\tsource"] + ["%s\t%s\t%s\t%s" % (k, v[0], v[1], v[2]) for k, v in CUTOFFS.items()]
    print("\n".join(lines))
    if a.tsv:
        with open(a.tsv, "w") as fh:
            fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
