#!/usr/bin/env python3
"""altanalyze3 LOH scan -- cohort loss-of-heterozygosity (LOH) from bulk allelic imbalance of germline
heterozygous SNPs, per KINNEX sample.

Consumes each sample's unfiltered germline VCF (annotation/gatk/germline.vcf.gz, the LOH-panel source)
+ the BAM, re-measures the true REF/ALT allele fraction at each germline het SNP directly from the BAM
(GATK AD is unusable on KINNEX), and calls LOH regions where het-SNP allelic fraction collapses away
from 0.5 coherently; a cohort aggregate compares mutant vs wild-type.

EFFICIENCY: the per-sample scan does ONE pileup pass per chromosome and tallies only at the het-SNP
positions (every read visited once). The earlier design called pileup() ONCE PER SNP (~66k per sample),
re-fetching the overlapping reads tens of thousands of times -- hours per sample; this is ~50-100x
faster (minutes).

Modes:
  (A) per-sample scan:  --germline VCF --bam BAM --sample S --out DIR
        -> S_het_snp_vaf.tsv, S_loh_segments.bed, S_arm_imbalance.tsv
  (B) cohort aggregate: --aggregate --manifest M --out DIR
        -> cohort_arm_imbalance.tsv (samples x chr-arm mean |VAF-0.5|), cohort_loh_heatmap.pdf

Run: python -m altanalyze3.components.bam.loh_scan ...
"""
import argparse
import csv
import os
from collections import defaultdict, OrderedDict

# chromosome-arm centromere table (GRCh38, approx Mb)
CHROM_ARMS = {
    'chr1': 123.4e6, 'chr2': 93.9e6, 'chr3': 90.9e6, 'chr4': 50.4e6, 'chr5': 48.4e6, 'chr6': 61.0e6,
    'chr7': 59.9e6, 'chr8': 45.6e6, 'chr9': 49.0e6, 'chr10': 40.2e6, 'chr11': 53.7e6, 'chr12': 35.8e6,
    'chr13': 17.9e6, 'chr14': 17.6e6, 'chr15': 19.0e6, 'chr16': 36.6e6, 'chr17': 25.1e6, 'chr18': 18.5e6,
    'chr19': 26.2e6, 'chr20': 28.1e6, 'chr21': 12.0e6, 'chr22': 15.0e6, 'chrX': 61.0e6,
}
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def arm_label(chrom, pos):
    c = CHROM_ARMS.get(chrom, 0)
    return chrom.replace("chr", "") + ("p" if pos < c else "q")


def segment_loh(snps, win=25, dev_thresh=0.20, min_run=15):
    """snps: sorted list of (pos, dev) where dev=|vaf-0.5| for covered het SNPs on ONE chrom. Rolling
    median deviation; call LOH runs of >= min_run SNPs whose windowed median dev stays >= dev_thresh.
    Returns list of (start_pos, end_pos, n_snp, mean_dev)."""
    n = len(snps)
    if n < min_run:
        return []
    devs = [d for _, d in snps]
    flag = [False] * n
    for i in range(n):
        lo = max(0, i - win // 2); hi = min(n, i + win // 2 + 1)
        w = sorted(devs[lo:hi]); med = w[len(w) // 2]
        flag[i] = med >= dev_thresh
    segs = []
    i = 0
    while i < n:
        if not flag[i]:
            i += 1
            continue
        j = i
        while j < n and flag[j]:
            j += 1
        if j - i >= min_run:
            run = snps[i:j]
            segs.append((run[0][0], run[-1][0], j - i, sum(d for _, d in run) / (j - i)))
        i = j
    return segs


def per_sample_scan(args):
    """Re-measure het-SNP allele fractions from the BAM and call LOH. ONE read-pass per chromosome."""
    import bisect
    import pysam
    vf = pysam.VariantFile(args.germline)
    bam = pysam.AlignmentFile(args.bam, "rb")
    os.makedirs(args.out, exist_ok=True)
    chroms = args.chroms.split(",") if args.chroms else AUTOSOMES
    reflen = dict(zip(bam.references, bam.lengths))
    by_chrom = defaultdict(list)          # chrom -> [(pos, vaf)]
    vaf_rows = []
    for chrom in chroms:
        if chrom not in reflen:
            continue
        # het-SNP positions on this chrom (ONE germline-VCF pass) -> pos0 -> (ref, alt)
        want = {}
        try:
            it = vf.fetch(chrom)
        except (ValueError, KeyError):
            continue
        for rec in it:
            if not rec.alts or len(rec.ref) != 1 or any(len(a) != 1 for a in rec.alts):
                continue
            s = rec.samples[0]
            if s.get("GT") not in ((0, 1), (1, 0)):
                continue
            want[rec.pos - 1] = (rec.ref.upper(), rec.alts[0].upper())
        if not want:
            continue
        # ONE fetch over the chromosome's reads; walk each read's CIGAR once and tally the base ONLY at
        # het-SNP positions (binary-searched within each aligned block). No per-column pileup is built --
        # cost is O(reads), not O(covered columns) or O(SNPs x pileup-setup). This is the mt_variants
        # single-pass pattern and is the fast primitive for scattered positions.
        want_pos = sorted(want)
        counts = {}                                     # pos0 -> [ref_c, alt_c]
        for read in bam.fetch(chrom):
            if read.is_unmapped or read.cigartuples is None:
                continue
            seq = read.query_sequence
            if not seq:
                continue
            refpos = read.reference_start
            qpos = 0
            for op, ln in read.cigartuples:
                if op == 0 or op == 7 or op == 8:       # M / = / X : consume ref + query
                    i0 = bisect.bisect_left(want_pos, refpos)
                    i1 = bisect.bisect_left(want_pos, refpos + ln)
                    for ii in range(i0, i1):
                        p0 = want_pos[ii]
                        b = seq[qpos + (p0 - refpos)].upper()
                        ref, alt = want[p0]
                        c = counts.get(p0)
                        if c is None:
                            c = counts[p0] = [0, 0]
                        if b == ref:
                            c[0] += 1
                        elif b == alt:
                            c[1] += 1
                    refpos += ln; qpos += ln
                elif op == 1 or op == 4:                # I / S : consume query only
                    qpos += ln
                elif op == 2 or op == 3:                # D / N : consume ref only
                    refpos += ln
        for p0 in sorted(counts):
            rc, ac = counts[p0]
            dp = rc + ac
            if dp < args.min_cov:
                continue
            ref, alt = want[p0]
            vaf = ac / dp
            by_chrom[chrom].append((p0 + 1, vaf))
            vaf_rows.append((chrom, p0 + 1, ref, alt, rc, ac, round(vaf, 4)))
    # write per-SNP VAF
    with open(os.path.join(args.out, f"{args.sample}_het_snp_vaf.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom", "pos", "ref", "alt", "ref_c", "alt_c", "vaf"])
        for row in vaf_rows:
            w.writerow(row)
    # segment LOH per chrom + arm imbalance
    seg_rows = []
    arm_dev = defaultdict(list)
    for chrom, snps in by_chrom.items():
        snps.sort()
        dev_snps = [(p, abs(v - 0.5)) for p, v in snps]
        for (a0, a1, nsnp, md) in segment_loh(dev_snps, args.win, args.dev_thresh, args.min_run):
            seg_rows.append((chrom, a0, a1, nsnp, round(md, 4)))
        for p, v in snps:
            arm_dev[arm_label(chrom, p)].append(abs(v - 0.5))
    with open(os.path.join(args.out, f"{args.sample}_loh_segments.bed"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom", "start", "end", "n_het_snp", "mean_abs_dev"])
        for row in sorted(seg_rows):
            w.writerow(row)
    with open(os.path.join(args.out, f"{args.sample}_arm_imbalance.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["arm", "n_het_snp", "mean_abs_dev_from_0.5"])
        for arm in sorted(arm_dev):
            vals = arm_dev[arm]
            w.writerow([arm, len(vals), round(sum(vals) / len(vals), 4)])
    print(f"[loh_scan] {args.sample}: {len(vaf_rows)} het SNPs, {len(seg_rows)} LOH segments "
          f"-> {args.out}/{args.sample}_loh_segments.bed")


def cohort_aggregate(args):
    man = OrderedDict()
    with open(args.manifest) as f:
        for d in csv.DictReader(f, delimiter="\t"):
            man[d["sample"]] = d
    arms_order = [f"{i}{a}" for i in list(range(1, 23)) for a in ("p", "q")]
    mat = {}
    present = []
    for sample in man:
        p = os.path.join(args.out, f"{sample}_arm_imbalance.tsv")
        if not os.path.exists(p):
            continue
        present.append(sample); mat[sample] = {}
        with open(p) as f:
            for d in csv.DictReader(f, delimiter="\t"):
                mat[sample][d["arm"]] = float(d["mean_abs_dev_from_0.5"])
    out_tsv = os.path.join(args.out, "cohort_arm_imbalance.tsv")
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "group"] + arms_order)
        for s in present:
            w.writerow([s, man[s].get("group", "")] + [mat[s].get(a, "") for a in arms_order])
    made = draw_heatmap(present, arms_order, man, mat, os.path.join(args.out, "cohort_loh_heatmap.pdf"))
    print(f"[loh_scan] cohort arm-imbalance for {len(present)} samples -> {out_tsv}"
          + ("" if made else "  (heatmap skipped)"))


def draw_heatmap(samples, arms, man, mat, out_pdf):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
    except Exception:
        return False
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    cmap = LinearSegmentedColormap.from_list("cyan_yellow", ["#00FFFF", "#FFFF00"])
    order = sorted(samples, key=lambda s: (man.get(s, {}).get("group", ""), s))
    import numpy as np
    M = np.array([[mat[s].get(a, np.nan) for a in arms] for s in order], dtype=float)
    fig, ax = plt.subplots(figsize=(max(8, len(arms) * 0.28), max(4, len(order) * 0.24)))
    im = ax.imshow(M, aspect="auto", cmap=cmap, vmin=0, vmax=0.5)
    ax.set_xticks(range(len(arms))); ax.set_xticklabels(arms, rotation=90, fontsize=6)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{s} [{man.get(s, {}).get('group', '')}]" for s in order], fontsize=5)
    ax.set_title("Bulk het-SNP allelic imbalance (|VAF-0.5|) -- LOH landscape", fontsize=10)
    cb = fig.colorbar(im, ax=ax, fraction=0.02); cb.set_label("mean |VAF-0.5| (0.5 = full LOH)")
    fig.tight_layout(); fig.savefig(out_pdf); plt.close(fig)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aggregate", action="store_true", help="cohort aggregate mode")
    ap.add_argument("--manifest")
    ap.add_argument("--out", required=True)
    ap.add_argument("--germline")
    ap.add_argument("--bam")
    ap.add_argument("--sample")
    ap.add_argument("--chroms", default="", help="comma list (default chr1..chr22)")
    ap.add_argument("--min-cov", dest="min_cov", type=int, default=20)
    ap.add_argument("--win", type=int, default=25)
    ap.add_argument("--dev-thresh", dest="dev_thresh", type=float, default=0.20)
    ap.add_argument("--min-run", dest="min_run", type=int, default=15)
    args = ap.parse_args()
    if args.aggregate:
        if not args.manifest:
            ap.error("--aggregate needs --manifest")
        cohort_aggregate(args)
    else:
        for r in ("germline", "bam", "sample"):
            if not getattr(args, r):
                ap.error(f"per-sample scan needs --{r}")
        per_sample_scan(args)


if __name__ == "__main__":
    main()
