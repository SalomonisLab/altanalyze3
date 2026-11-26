#!/usr/bin/env python3
"""Per-read alleles at every candidate germline site on one chromosome of one library BAM, joined to the
molecule's gene, cell state and final isoform.

One pass streams the chromosome's alignments. Each read gets a call at every panel site inside its aligned
span: REF, ALT, OTHER, SPLICE (the site sits in an N gap of the read), DEL (inside a deletion) for SNVs, and
REF, ALT, OTHER, SPLICE for indels. The read table keeps the base calls (REF, ALT, OTHER); SPLICE and DEL are
counted in the site table only, because a read across a long intron would otherwise carry one SPLICE call for
every candidate site inside the intron. SNV calls reuse ``read_alleles.classify`` (and through it the validated
``variant_scan.base_and_class_at``); indel calls reuse ``variant_scan.indels_in_window`` with an exact
type-and-length match inside +/- ``--indel_window`` bases. The read filter is read_alleles' filter: skip
unmapped reads, reads under ``--min_mapq`` and reads without a CIGAR; flag secondary and supplementary
alignments without dropping them.

Molecule join. The isoform workflow names a read ``molecule/<N>``. The ISV molecule index
``<library>.reads_final.db`` (table rf: gene, mol, ct, cnt, final) gives each molecule its gene, cell state
and final isoform. When ``--composed_index_bindings`` is given, the module installs the approved
cleanup-viewer composition (``cleanup_composed_index.install``) before opening the index, so final IDs pass
through the direct redirect table exactly as the ISV viewer serves them.

Run:
  python -m altanalyze3.components.snp.read_sites --uid BM27_CD34 --chrom chr7 \\
      --isv_metadata isv/metadata.tsv --bam_metadata mds_metadata_bam_fixed.txt --donor_map donor_map.tsv \\
      --panel_dir PANEL --molecule_index MOL_INDEX_DIR --gene_model Hs_Ensembl_exon.txt --fasta genome.fa \\
      [--composed_index_bindings CLEANUP_VIEWER_INPUT_STAGING.json --composed_index_helper cleanup_composed_index.py] \\
      [--extra_sites controls.tsv] --out OUT
Outputs in OUT/<uid>/:
  <chrom>.reads.tsv.gz        one row per alignment with at least one site call inside its span
  <chrom>.site_counts.tsv.gz  per site: allele classes over all alignments and over primary alignments
  logs/<chrom>.parameters.json
"""
import argparse
import bisect
import csv
import datetime
import gzip
import hashlib
import importlib.util
import json
import os
import sqlite3
import sys
import time
from collections import Counter, defaultdict

import pysam

try:
    from altanalyze3.components.bam.read_alleles import classify, load_samples, molecule_id, site_state
    from altanalyze3.components.bam.variant_scan import core_bc, get_bc, indels_in_window, rc
except ImportError:  # run beside a staged copy
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bam"))
    from read_alleles import classify, load_samples, molecule_id, site_state
    from variant_scan import core_bc, get_bc, indels_in_window, rc

SNV_CLASSES = ("REF", "ALT", "OTHER", "SPLICE", "DEL")
READ_TABLE_CLASSES = ("REF", "ALT", "OTHER")   # SPLICE and DEL are counted per site, not written per read
READ_FIELDS = ["molecule_id", "read_name", "cb_rc", "flag", "secondary", "supplementary", "mapq", "strand",
               "ref_start", "ref_end", "introns", "n_ref_alt", "calls", "index_rows"]
COUNT_FIELDS = ["site_id", "chrom", "pos", "ref", "alt", "vtype", "tier", "source"] + \
               ["all_" + c for c in SNV_CLASSES] + ["primary_" + c for c in SNV_CLASSES]


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def site_id(s):
    return "{}:{}:{}:{}".format(s["chrom"], s["pos"], s["ref"], s["alt"])


def load_sites(panel_path, chrom, extra_path=None):
    """Panel rows on `chrom` (SNV, INS, DEL) plus named extra sites; sorted by position."""
    sites = {}
    with gzip.open(panel_path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["chrom"] != chrom or r["vtype"] not in ("SNV", "INS", "DEL"):
                continue
            s = dict(chrom=r["chrom"], pos=int(r["pos"]), ref=r["ref"], alt=r["alt"], vtype=r["vtype"],
                     tier=r["tier"], source="panel")
            sites[site_id(s)] = s
    if extra_path:
        with open(extra_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["chrom"] != chrom:
                    continue
                ref, alt = r["ref"].upper(), r["alt"].upper()
                vt = "SNV" if len(ref) == len(alt) == 1 else ("INS" if len(alt) > len(ref) else "DEL")
                s = dict(chrom=chrom, pos=int(r["pos"]), ref=ref, alt=alt, vtype=vt, tier="named", source=r["id"])
                k = site_id(s)
                if k in sites:
                    sites[k]["source"] = "panel+" + r["id"]
                else:
                    sites[k] = s
    out = sorted(sites.values(), key=lambda s: (s["pos"], s["ref"], s["alt"]))
    for i, s in enumerate(out):
        s["sid"] = site_id(s)
    return out


def check_fasta(sites, fasta):
    """Every REF must equal the FASTA bases. Returns the number checked."""
    fa = pysam.FastaFile(fasta)
    bad = []
    for s in sites:
        contig = s["chrom"] if s["chrom"] in fa.references else s["chrom"].replace("chr", "", 1)
        if contig not in fa.references:
            raise SystemExit("FASTA has no contig {}".format(s["chrom"]))
        seq = fa.fetch(contig, s["pos"] - 1, s["pos"] - 1 + len(s["ref"])).upper()
        if seq != s["ref"]:
            bad.append((s["sid"], seq))
    if bad:
        raise SystemExit("{} sites have REF different from the FASTA, e.g. {}".format(len(bad), bad[:5]))
    return len(sites)


def classify_indel(read, s, window):
    """REF / ALT / OTHER / SPLICE / NOT_SPANNED for one INS or DEL site."""
    length = abs(len(s["alt"]) - len(s["ref"]))
    kind = "ins" if s["vtype"] == "INS" else "del"
    target0 = s["pos"]                   # 0-based position of the first base after the VCF anchor base
    state = site_state(read, target0 if kind == "del" else s["pos"] - 1)
    if state == "SPLICE":
        return "SPLICE"
    near = indels_in_window(read, target0, window)
    if any(t == kind and ln == length for t, ln, _ in near):
        return "ALT"
    if near:
        return "OTHER"
    # REF needs aligned bases on both sides of the event
    left, right = s["pos"] - 1, s["pos"] + (length if kind == "del" else 0)
    if site_state(read, left) == "ALIGNED" and site_state(read, right) == "ALIGNED":
        return "REF"
    return "NOT_SPANNED"


def introns_of(read):
    out = []
    ref = read.reference_start
    for op, ln in read.cigartuples:
        if op == 3:
            out.append("{}-{}".format(ref + 1, ref + ln))
        if op in (0, 2, 3, 7, 8):
            ref += ln
    return ";".join(out)


def load_genes(gene_model, chrom):
    """Genes on `chrom` from an Ensembl exon table (gene in column 1, chromosome in column 3)."""
    genes = set()
    bare = chrom.replace("chr", "", 1)
    with open(gene_model) as fh:
        for line in fh:
            p = line.split("\t", 4)
            if len(p) > 3 and p[2].replace("chr", "", 1) == bare:
                genes.add(p[0])
    genes.discard("gene")
    return sorted(genes)


def install_composed(bindings_json, helper_path):
    q = json.load(open(bindings_json))
    spec = importlib.util.spec_from_file_location("cleanup_composed_index", helper_path)
    ci = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ci)
    if q.get("helper_sha256") and _sha256(helper_path) != q["helper_sha256"]:
        raise SystemExit("composed-index helper sha256 differs from the bindings record")
    ci.install(q["descriptor_bindings"])
    by_name = defaultdict(list)
    for p in q["descriptor_bindings"]:
        by_name[os.path.basename(p)].append(p)
    return by_name, q.get("redirect_database")


def molecule_rows(db_path, genes, wanted):
    """{bare molecule id: [(gene, ct, final), ...]} for molecules in `wanted`, over the given genes."""
    con = sqlite3.connect("file:{}?mode=ro".format(db_path), uri=True)
    out = defaultdict(list)
    n = 0
    for i in range(0, len(genes), 500):
        chunk = genes[i:i + 500]
        q = "SELECT gene, mol, ct, final FROM rf WHERE gene IN ({})".format(",".join("?" * len(chunk)))
        for gene, mol, ct, fin in con.execute(q, chunk):
            n += 1
            if mol in wanted:
                out[mol].append((gene, ct, fin if fin is not None else ""))
    con.close()
    return out, n


def scan(bam_path, chrom, sites, min_mapq, window):
    """Yield one dict per alignment with at least one site call inside its span; fill per-site counts."""
    pos = [s["pos"] for s in sites]
    counts = {s["sid"]: Counter() for s in sites}
    n_seen = n_kept = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if chrom not in bam.references:
            raise SystemExit("{} has no contig {}".format(bam_path, chrom))
        for read in bam.fetch(chrom):
            if read.is_unmapped or read.mapping_quality < min_mapq or read.cigartuples is None:
                continue
            n_seen += 1
            start, end = read.reference_start, read.reference_end
            lo = bisect.bisect_left(pos, start)           # VCF pos is 1-based: pos-1 >= start  <=> pos > start
            hi = bisect.bisect_right(pos, end)
            if lo >= hi:
                continue
            primary = not (read.is_secondary or read.is_supplementary)
            calls = []
            n_ra = 0
            for s in sites[lo:hi]:
                if not (start <= s["pos"] - 1 < end):
                    continue
                if s["vtype"] == "SNV":
                    _base, cls, _op = classify(read, dict(pos=s["pos"], ref=s["ref"], alt=s["alt"]))
                else:
                    cls = classify_indel(read, s, window)
                if cls == "NOT_SPANNED":
                    continue
                if cls in READ_TABLE_CLASSES:
                    calls.append("{}={}".format(s["sid"].split(":", 1)[1], cls))
                counts[s["sid"]]["all_" + cls] += 1
                if primary:
                    counts[s["sid"]]["primary_" + cls] += 1
                n_ra += cls in ("REF", "ALT")
            if not calls:
                continue
            n_kept += 1
            bc = get_bc(read) or ""
            cb = core_bc(bc) if bc else ""
            yield dict(read_name=read.query_name, cb_rc=rc(cb) if cb else "", flag=read.flag,
                       secondary=int(read.is_secondary), supplementary=int(read.is_supplementary),
                       mapq=read.mapping_quality, strand="-" if read.is_reverse else "+",
                       ref_start=start + 1, ref_end=end, introns=introns_of(read), n_ref_alt=n_ra,
                       calls=";".join(calls))
    scan.stats = dict(alignments_passing_filter=n_seen, alignments_with_calls=n_kept)
    scan.counts = counts


def build_arg_parser():
    ap = argparse.ArgumentParser(prog="python -m altanalyze3.components.snp.read_sites",
                                 description="Per-read alleles at candidate germline sites, joined to isoforms.")
    ap.add_argument("--uid", required=True, help="ISV library uid (a row of --isv_metadata)")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--isv_metadata", required=True)
    ap.add_argument("--bam_metadata", required=True)
    ap.add_argument("--donor_map", required=True, help="TSV with uid (ISV uid) and donor columns")
    ap.add_argument("--panel_dir", required=True, help="het_panel output folder (donors/<donor>.sites.tsv.gz)")
    ap.add_argument("--extra_sites", default=None, help="TSV id, chrom, pos, ref, alt: always scanned")
    ap.add_argument("--molecule_index", required=True, help="folder with <uid>.reads_final.db")
    ap.add_argument("--composed_index_bindings", default=None)
    ap.add_argument("--composed_index_helper", default=None)
    ap.add_argument("--gene_model", required=True, help="Ensembl exon table: gene, exon, chromosome, ...")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--min_mapq", type=int, default=20, help="default 20, as read_alleles and variant_scan")
    ap.add_argument("--indel_window", type=int, default=5, help="default 5, as variant_scan.scan_indel")
    ap.add_argument("--out", required=True)
    return ap


def main(argv=None):
    a = build_arg_parser().parse_args(argv)
    t0 = time.time()
    started = datetime.datetime.now().isoformat(timespec="seconds")
    samples = dict(load_samples(argparse.Namespace(samples=None, isv_metadata=a.isv_metadata,
                                                   bam_metadata=a.bam_metadata)))
    if a.uid not in samples:
        raise SystemExit("uid {} is not in {}".format(a.uid, a.isv_metadata))
    bam = samples[a.uid]
    donors = {}
    with open(a.donor_map) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            donors[r["uid"]] = r["donor"]
    donor = donors.get(a.uid)
    if not donor:
        raise SystemExit("uid {} has no donor in {}".format(a.uid, a.donor_map))
    panel = os.path.join(a.panel_dir, "donors", donor + ".sites.tsv.gz")
    if not os.path.exists(panel):
        raise SystemExit("no donor panel {}".format(panel))
    sites = load_sites(panel, a.chrom, a.extra_sites)
    n_fa = check_fasta(sites, a.fasta)
    print("[read_sites] {} ({}) {}: {:,} sites ({}), REF checked against FASTA for {:,}".format(
        a.uid, donor, a.chrom, len(sites), dict(Counter(s["vtype"] for s in sites)), n_fa), flush=True)

    db_name = "{}.reads_final.db".format(a.uid)
    db_path = os.path.join(a.molecule_index, db_name)
    composed = None
    if a.composed_index_bindings:
        if not a.composed_index_helper:
            raise SystemExit("--composed_index_bindings needs --composed_index_helper")
        by_name, redirect = install_composed(a.composed_index_bindings, a.composed_index_helper)
        match = [p for p in by_name.get(db_name, []) if os.path.realpath(p) == os.path.realpath(db_path)]
        if not match:
            raise SystemExit("no composed binding resolves to {} (bindings named {}: {})".format(
                db_path, db_name, by_name.get(db_name, [])))
        db_path = match[0]       # the binding key itself: the installed connect() matches on this path
        composed = dict(bindings=a.composed_index_bindings, helper=a.composed_index_helper, redirect=redirect)
    if not os.path.exists(db_path):
        raise SystemExit("no molecule index {}".format(db_path))

    out_dir = os.path.join(a.out, a.uid)
    os.makedirs(os.path.join(out_dir, "logs"), exist_ok=True)
    read_path = os.path.join(out_dir, "{}.reads.tsv.gz".format(a.chrom))
    stage1 = read_path + ".stage1.tmp"
    wanted = set()
    n_reads = n_ra_reads = n_calls = 0
    with gzip.open(stage1, "wt") as fh:
        w = csv.DictWriter(fh, fieldnames=READ_FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in scan(bam, a.chrom, sites, a.min_mapq, a.indel_window):
            r["molecule_id"] = molecule_id(r["read_name"], a.uid)
            r["index_rows"] = ""
            if r["read_name"].startswith("molecule/"):
                wanted.add(r["read_name"].split("/", 1)[1])
            n_reads += 1
            n_ra_reads += r["n_ref_alt"] > 0
            n_calls += r["calls"].count("=")
            w.writerow(r)
    stats = dict(scan.stats)
    t_scan = time.time() - t0
    genes = load_genes(a.gene_model, a.chrom)
    mrows, n_rf = molecule_rows(db_path, genes, wanted)
    joined = n_rows_out = 0
    with gzip.open(stage1, "rt") as fin, gzip.open(read_path + ".tmp", "wt") as fout:
        rd = csv.DictReader(fin, delimiter="\t")
        w = csv.DictWriter(fout, fieldnames=READ_FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rd:
            bare = r["read_name"].split("/", 1)[1] if r["read_name"].startswith("molecule/") else None
            rows = mrows.get(bare, []) if bare else []
            r["index_rows"] = ";".join("{}|{}|{}".format(*x) for x in rows)
            joined += bool(rows)
            n_rows_out += 1
            w.writerow(r)
    if n_rows_out != n_reads:
        raise AssertionError("stage-2 rows {} != stage-1 rows {}".format(n_rows_out, n_reads))
    os.replace(read_path + ".tmp", read_path)
    os.remove(stage1)
    count_path = os.path.join(out_dir, "{}.site_counts.tsv.gz".format(a.chrom))
    n_cov = 0
    with gzip.open(count_path + ".tmp", "wt") as fh:
        w = csv.DictWriter(fh, fieldnames=COUNT_FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for s in sites:
            c = scan.counts[s["sid"]]
            if not c:
                continue
            n_cov += 1
            row = dict(site_id=s["sid"], chrom=s["chrom"], pos=s["pos"], ref=s["ref"], alt=s["alt"],
                       vtype=s["vtype"], tier=s["tier"], source=s["source"])
            for k in COUNT_FIELDS[8:]:
                row[k] = c.get(k, 0)
            w.writerow(row)
    os.replace(count_path + ".tmp", count_path)
    # partition check: every call in the read table is counted once in the site table
    n_counted = sum(sum(c.get("all_" + k, 0) for k in READ_TABLE_CLASSES) for c in scan.counts.values())
    if n_calls != n_counted:
        raise AssertionError("read calls {} != site counts {}".format(n_calls, n_counted))
    log = dict(module="altanalyze3.components.snp.read_sites", argv=sys.argv if argv is None else list(argv),
               started=started, finished=datetime.datetime.now().isoformat(timespec="seconds"),
               seconds_scan=round(t_scan, 1), seconds_total=round(time.time() - t0, 1),
               uid=a.uid, donor=donor, bam=bam, chrom=a.chrom, panel=panel, extra_sites=a.extra_sites,
               sites=len(sites), sites_by_type=dict(Counter(s["vtype"] for s in sites)), sites_covered=n_cov,
               min_mapq=a.min_mapq, indel_window=a.indel_window, fasta=a.fasta, molecule_index=db_path,
               composed_index=composed, genes_queried=len(genes), rf_rows_read=n_rf,
               reads_with_calls=n_reads, reads_with_ref_or_alt=n_ra_reads,
               reads_joined_to_index=joined, molecules_queried=len(wanted), calls=n_calls, **stats,
               pysam=pysam.__version__, python=sys.version.split()[0],
               module_sha256=_sha256(os.path.abspath(__file__)),
               outputs=dict(reads=read_path, site_counts=count_path))
    with open(os.path.join(out_dir, "logs", "{}.parameters.json".format(a.chrom)), "w") as fh:
        json.dump(log, fh, indent=2)
    print("[read_sites] {} {}: {:,} alignments with calls ({:,} with REF/ALT), {:,} joined to the index; "
          "{:,} calls; scan {:.0f}s, total {:.0f}s".format(a.uid, a.chrom, n_reads, log["reads_with_ref_or_alt"],
                                                          joined, n_calls, t_scan, time.time() - t0), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
