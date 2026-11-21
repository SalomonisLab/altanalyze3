#!/usr/bin/env python3
"""Per-read allele at known SNVs in long-read BAMs, keyed by the ISV molecule id.

Each read that aligns across a listed site gets one allele call. The call is the base the read carries
there, or the reason it carries none: the read splices over the site, deletes it, or stops before it.
Library counts and a binomial genotype call are written beside the per-read table.

Base calls reuse the validated reference-free genotyper ``variant_scan.base_and_class_at``. The read
filter is the one ``variant_scan.scan_snv`` applies: skip unmapped reads, reads under ``--min_mapq`` and
reads without a CIGAR. Secondary and supplementary alignments are flagged, not dropped.

Molecule id. The isoform workflow names a read ``molecule/<N>`` and the ISV read index names the same
molecule ``<N>.<uid>``. This module writes that id, so the table joins to ``isv_web_remote
--variant_alleles`` and to ``--export_molecules`` without any renaming. A read name that does not start
with ``molecule/`` keeps its full name before the ``.<uid>`` suffix.

Samples come from one of two inputs:
  --samples TSV          columns uid, bam
  --isv_metadata TSV     the ISV run metadata (uid, original_metadata_ordinal, original_library), joined
  --bam_metadata TSV     to the workflow BAM metadata (bam, library) by row ordinal; the library names
                         must agree on every row or the run stops.

Run:
  python -m altanalyze3.components.bam.read_alleles \\
      --variant rs1859788:chr7:100374211:A:G \\
      --isv_metadata <run>/isv/metadata.tsv --bam_metadata mds_metadata_bam_revised.txt \\
      --fasta genome.fa --region chr7:100367530-100400000 --out <outdir>

Outputs in --out:
  read_alleles.tsv            one row per (variant, read)
  library_allele_counts.tsv   one row per (variant, uid) with counts and the genotype call
  logs/read_alleles.parameters.json
"""
import argparse
import csv
import datetime
import hashlib
import json
import math
import os
import sys
from collections import Counter
from multiprocessing import Pool

import pysam

try:
    from altanalyze3.components.bam.variant_scan import base_and_class_at, core_bc, get_bc, rc
except ImportError:  # run as a plain script beside variant_scan.py
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from variant_scan import base_and_class_at, core_bc, get_bc, rc

ALLELE_CLASSES = ("REF", "ALT", "OTHER", "SPLICE", "DEL", "NOT_SPANNED")
READ_FIELDS = ["variant_id", "chrom", "pos", "ref", "alt", "uid", "molecule_id", "read_name", "cb_bam",
               "cb_rc", "flag", "secondary", "supplementary", "mapq", "base", "allele_class", "cigar_op"]
COUNT_FIELDS = ["variant_id", "chrom", "pos", "ref", "alt", "uid", "bam", "reads", "REF", "ALT", "OTHER",
                "SPLICE", "DEL", "NOT_SPANNED", "secondary_or_supplementary", "duplicate_molecule_ids",
                "alt_fraction", "genotype", "GQ", "log10_L_hom_ref", "log10_L_het", "log10_L_hom_alt"]


def parse_variant(text):
    """'rs1859788:chr7:100374211:A:G' -> dict(id, chrom, pos, ref, alt)."""
    parts = str(text).strip().split(":")
    if len(parts) != 5:
        raise ValueError("--variant expects ID:CHROM:POS:REF:ALT, got {!r}".format(text))
    vid, chrom, pos, ref, alt = parts
    ref, alt = ref.upper(), alt.upper()
    if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
        raise ValueError("only single-base REF/ALT SNVs are supported, got {!r}".format(text))
    return dict(id=vid, chrom=chrom, pos=int(pos), ref=ref, alt=alt)


def load_variants(args):
    variants = [parse_variant(v) for v in (args.variant or [])]
    if args.variants:
        with open(args.variants) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                variants.append(parse_variant(":".join([row["id"], row["chrom"], row["pos"], row["ref"], row["alt"]])))
    if not variants:
        raise ValueError("give at least one --variant or a --variants TSV")
    ids = [v["id"] for v in variants]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate variant ids: {}".format(sorted(k for k, n in Counter(ids).items() if n > 1)))
    return variants


def parse_region(text):
    if not text:
        return None
    chrom, span = str(text).replace(",", "").split(":", 1)
    a, b = (int(x) for x in span.split("-", 1))
    if b <= a:
        raise ValueError("--region END must exceed START, got {!r}".format(text))
    return chrom, a, b


def load_samples(args):
    """[(uid, bam)] from --samples, or from --isv_metadata joined to --bam_metadata by ordinal."""
    if args.samples:
        with open(args.samples) as fh:
            rows = [(r["uid"], r["bam"]) for r in csv.DictReader(fh, delimiter="\t")]
    else:
        if not (args.isv_metadata and args.bam_metadata):
            raise ValueError("give --samples, or both --isv_metadata and --bam_metadata")
        with open(args.bam_metadata) as fh:
            bam_rows = list(csv.DictReader(fh, delimiter="\t"))
        rows = []
        with open(args.isv_metadata) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                k = int(r["original_metadata_ordinal"])
                if not 0 <= k < len(bam_rows):
                    raise ValueError("ordinal {} of {} is outside the {} BAM metadata rows".format(k, r["uid"], len(bam_rows)))
                b = bam_rows[k]
                if b["library"] != r["original_library"]:
                    raise ValueError("library mismatch at ordinal {}: BAM metadata {!r} vs ISV {!r}".format(
                        k, b["library"], r["original_library"]))
                rows.append((r["uid"], b["bam"]))
    uids = [u for u, _ in rows]
    if len(set(uids)) != len(uids):
        raise ValueError("duplicate uids in the sample list")
    if not rows:
        raise ValueError("no samples")
    return rows


def molecule_id(read_name, uid):
    num = read_name.split("/", 1)[1] if read_name.startswith("molecule/") else read_name
    return "{}.{}".format(num, uid)


def site_state(read, pos0):
    """Where pos0 falls in the read's alignment: 'SPLICE' (N), 'DEL' (D) or 'NOT_SPANNED'."""
    ref = read.reference_start
    for op, ln in read.cigartuples:
        if op in (0, 2, 3, 7, 8):
            if ref <= pos0 < ref + ln:
                if op == 3:
                    return "SPLICE"
                if op == 2:
                    return "DEL"
                return "ALIGNED"
            ref += ln
        if ref > pos0:
            break
    return "NOT_SPANNED"


def classify(read, variant):
    """(base, allele_class, cigar_op) of one read at one variant."""
    pos0 = variant["pos"] - 1
    base, op = base_and_class_at(read, pos0)
    if base is not None:
        if base == variant["ref"]:
            return base, "REF", op
        if base == variant["alt"]:
            return base, "ALT", op
        return base, "OTHER", op
    state = site_state(read, pos0)
    if state == "ALIGNED":  # base_and_class_at found no base on an aligned op: query sequence missing
        return "", "OTHER", ""
    return "", state, ""


def genotype_call(n_ref, n_alt, error_rate):
    """Binomial genotype likelihoods over REF and ALT reads, flat prior.
    Returns (genotype, GQ, log10 L(0/0), log10 L(0/1), log10 L(1/1)). GQ is capped at 99."""
    n = n_ref + n_alt
    if n == 0:
        return "./.", 0, float("nan"), float("nan"), float("nan")
    e = float(error_rate)
    ll = [n_alt * math.log10(e) + n_ref * math.log10(1 - e),
          n * math.log10(0.5),
          n_alt * math.log10(1 - e) + n_ref * math.log10(e)]
    best = max(range(3), key=lambda i: ll[i])
    mx = ll[best]
    post = [10 ** (x - mx) for x in ll]
    p_best = post[best] / sum(post)
    gq = 99 if p_best >= 1.0 else min(99, int(round(-10 * math.log10(max(1e-300, 1 - p_best)))))
    return ("0/0", "0/1", "1/1")[best], gq, ll[0], ll[1], ll[2]


def scan_sample(job):
    """Every read at every variant in one BAM. Returns (uid, rows)."""
    uid, bam_path, variants, region, min_mapq = job
    rows = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for v in variants:
            pos0 = v["pos"] - 1
            if region and region[0] == v["chrom"] and region[1] <= v["pos"] <= region[2]:
                fetch = (v["chrom"], region[1] - 1, region[2])
            else:
                fetch = (v["chrom"], pos0, pos0 + 1)
            if fetch[0] not in bam.references:
                raise ValueError("{} has no contig {}".format(bam_path, fetch[0]))
            for read in bam.fetch(*fetch):
                if read.is_unmapped or read.mapping_quality < min_mapq or read.cigartuples is None:
                    continue
                if read.reference_end is not None and read.reference_start <= pos0 < read.reference_end:
                    base, cls, op = classify(read, v)
                else:
                    base, cls, op = "", "NOT_SPANNED", ""
                bc = get_bc(read) or ""
                cb = core_bc(bc) if bc else ""
                rows.append(dict(variant_id=v["id"], chrom=v["chrom"], pos=v["pos"], ref=v["ref"], alt=v["alt"],
                                 uid=uid, molecule_id=molecule_id(read.query_name, uid), read_name=read.query_name,
                                 cb_bam=cb, cb_rc=rc(cb) if cb else "", flag=read.flag,
                                 secondary=int(read.is_secondary), supplementary=int(read.is_supplementary),
                                 mapq=read.mapping_quality, base=base, allele_class=cls, cigar_op=op))
    return uid, rows


def library_counts(uid, bam, variant, rows, error_rate):
    c = Counter(r["allele_class"] for r in rows)
    dup = sum(n - 1 for n in Counter(r["molecule_id"] for r in rows).values() if n > 1)
    n_ref, n_alt = c.get("REF", 0), c.get("ALT", 0)
    gt, gq, l0, l1, l2 = genotype_call(n_ref, n_alt, error_rate)
    af = n_alt / (n_ref + n_alt) if (n_ref + n_alt) else float("nan")
    out = dict(variant_id=variant["id"], chrom=variant["chrom"], pos=variant["pos"], ref=variant["ref"],
               alt=variant["alt"], uid=uid, bam=bam, reads=len(rows),
               secondary_or_supplementary=sum(1 for r in rows if r["secondary"] or r["supplementary"]),
               duplicate_molecule_ids=dup, alt_fraction=round(af, 4) if af == af else "NA",
               genotype=gt, GQ=gq,
               log10_L_hom_ref=round(l0, 4) if l0 == l0 else "NA",
               log10_L_het=round(l1, 4) if l1 == l1 else "NA",
               log10_L_hom_alt=round(l2, 4) if l2 == l2 else "NA")
    for k in ALLELE_CLASSES:
        out[k] = c.get(k, 0)
    if sum(c.values()) != len(rows):
        raise AssertionError("allele classes do not partition the reads of {}".format(uid))
    return out


def _sha256(path, limit=None):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read() if limit is None else fh.read(limit))
    return h.hexdigest()


def build_arg_parser():
    ap = argparse.ArgumentParser(prog="python -m altanalyze3.components.bam.read_alleles",
                                 description="Per-read allele at known SNVs in long-read BAMs, keyed by ISV molecule id.")
    ap.add_argument("--variant", action="append", help="ID:CHROM:POS:REF:ALT, 1-based POS (repeatable)")
    ap.add_argument("--variants", help="TSV with columns id, chrom, pos, ref, alt")
    ap.add_argument("--samples", help="TSV with columns uid, bam")
    ap.add_argument("--isv_metadata", help="ISV run metadata.tsv (uid, original_metadata_ordinal, original_library)")
    ap.add_argument("--bam_metadata", help="workflow BAM metadata (bam, library), row ordinal = original_metadata_ordinal")
    ap.add_argument("--fasta", required=True, help="reference FASTA; every REF must match it")
    ap.add_argument("--region", default=None,
                    help="CHROM:START-END, usually the gene. Reads in it that stop before the site are written as "
                         "NOT_SPANNED, so a molecule that never reaches the site is told apart from a missing one")
    ap.add_argument("--min_mapq", type=int, default=20, help="minimum MAPQ (default 20, as variant_scan)")
    ap.add_argument("--error_rate", type=float, default=0.01, help="per-read base error for the genotype likelihood")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--isv_molecules", default=None,
                    help="isv_web_remote --export_molecules TSV for a gene inside --region. Checks the BAM pairing: "
                         "every library's ISV molecules are looked up among its own BAM reads, the counts go to "
                         "isv_molecule_join.tsv, and the run stops when a library with ISV molecules matches none")
    ap.add_argument("--out", required=True)
    return ap


def check_isv_pairing(isv_molecules, results, out_dir):
    """Per library: ISV molecules found among that library's BAM reads. A library whose ISV molecules are all
    absent was read from the wrong BAM (two KINNEX-5 metadata rows once swapped their BAM paths)."""
    exported = {}
    with open(isv_molecules) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            exported.setdefault(r["sample"], set()).add(r["molecule_id"])
    rows, broken = [], []
    for uid, rr in results.items():
        found = {r["molecule_id"] for r in rr}
        ids = exported.get(uid, set())
        hit = len(ids & found)
        rows.append((uid, len(ids), hit))
        if ids and not hit:
            broken.append(uid)
    unknown = sorted(set(exported) - set(results))
    with open(os.path.join(out_dir, "isv_molecule_join.tsv"), "w") as fh:
        fh.write("uid\tisv_molecules\tfound_in_bam\tmissing\n")
        for uid, n, hit in rows:
            fh.write("{}\t{}\t{}\t{}\n".format(uid, n, hit, n - hit))
    tot, hit = sum(n for _, n, _ in rows), sum(h for _, _, h in rows)
    print("[read_alleles] ISV molecules found in their own library BAM: {:,} of {:,}".format(hit, tot), flush=True)
    if unknown:
        raise SystemExit("ISV export libraries absent from the sample list: {}".format(unknown))
    if broken:
        raise SystemExit("BAM pairing failure: none of the ISV molecules of {} are in the BAM given for them".format(broken))
    return tot, hit


def main(argv=None):
    a = build_arg_parser().parse_args(argv)
    started = datetime.datetime.now().isoformat(timespec="seconds")
    variants = load_variants(a)
    region = parse_region(a.region)
    samples = load_samples(a)
    fa = pysam.FastaFile(a.fasta)
    for v in variants:
        contig = v["chrom"] if v["chrom"] in fa.references else v["chrom"].replace("chr", "", 1)
        if contig not in fa.references:
            raise ValueError("FASTA has no contig {} for {}".format(v["chrom"], v["id"]))
        base = fa.fetch(contig, v["pos"] - 1, v["pos"]).upper()
        if base != v["ref"]:
            raise ValueError("{}: FASTA base {} at {}:{} differs from REF {}".format(v["id"], base, v["chrom"], v["pos"], v["ref"]))
        print("[read_alleles] {} {}:{} REF {} ALT {} (FASTA base {})".format(v["id"], v["chrom"], v["pos"], v["ref"], v["alt"], base), flush=True)
    os.makedirs(os.path.join(a.out, "logs"), exist_ok=True)
    jobs = [(uid, bam, variants, region, a.min_mapq) for uid, bam in samples]
    if a.threads > 1:
        with Pool(a.threads) as pool:
            results = dict(pool.map(scan_sample, jobs))
    else:
        results = dict(scan_sample(j) for j in jobs)

    isv_join = check_isv_pairing(a.isv_molecules, results, a.out) if a.isv_molecules else None
    read_path = os.path.join(a.out, "read_alleles.tsv")
    count_path = os.path.join(a.out, "library_allele_counts.tsv")
    n_reads = 0
    count_rows = []
    with open(read_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=READ_FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for uid, bam in samples:
            rows = results[uid]
            for r in rows:
                w.writerow(r)
            n_reads += len(rows)
            for v in variants:
                vrows = [r for r in rows if r["variant_id"] == v["id"]]
                count_rows.append(library_counts(uid, bam, v, vrows, a.error_rate))
    with open(count_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COUNT_FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(count_rows)
    if sum(r["reads"] for r in count_rows) != n_reads:
        raise AssertionError("library counts do not sum to the per-read table")
    for v in variants:
        vr = [r for r in count_rows if r["variant_id"] == v["id"]]
        tot = Counter()
        for r in vr:
            for k in ALLELE_CLASSES:
                tot[k] += r[k]
        print("[read_alleles] {}: {} libraries, {} reads: {}; genotypes {}".format(
            v["id"], len(vr), sum(tot.values()), dict(tot), dict(Counter(r["genotype"] for r in vr))), flush=True)
    params = dict(module="altanalyze3.components.bam.read_alleles", argv=sys.argv if argv is None else list(argv),
                  started=started, finished=datetime.datetime.now().isoformat(timespec="seconds"),
                  variants=variants, region=a.region, min_mapq=a.min_mapq, error_rate=a.error_rate,
                  samples=[dict(uid=u, bam=b) for u, b in samples], fasta=a.fasta,
                  module_sha256=_sha256(os.path.abspath(__file__)),
                  variant_scan_sha256=_sha256(sys.modules[base_and_class_at.__module__].__file__),
                  pysam=pysam.__version__, python=sys.version.split()[0],
                  isv_molecules=a.isv_molecules,
                  isv_molecules_found=None if isv_join is None else dict(total=isv_join[0], found=isv_join[1]),
                  outputs=dict(read_alleles=read_path, read_rows=n_reads,
                               library_allele_counts=count_path, library_rows=len(count_rows)))
    with open(os.path.join(a.out, "logs", "read_alleles.parameters.json"), "w") as fh:
        json.dump(params, fh, indent=2)
    print("[read_alleles] wrote {} read rows -> {}".format(n_reads, read_path), flush=True)
    print("[read_alleles] wrote {} library rows -> {}".format(len(count_rows), count_path), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
