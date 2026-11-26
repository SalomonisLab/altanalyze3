#!/usr/bin/env python3
"""Candidate germline variant panel per donor from the four-caller genome variant discovery outputs.

The discovery pipeline (KINNEX_unsupervised_variant_detection/kinnex_variant_pipeline.sh) writes four caller
VCFs beside each BAM, in ``<bam dir>/genome_variant_detection/results/``. Its merged confidence tiers
(``results/merged/*_confidence.vcf.gz``) leave out the GERMLINE bucket by design (merge_callers.py,
SOMATIC_BUCKETS), so a germline panel must come from the per-caller VCFs. This module reads them directly.

Stage ``--library``: one BAM. Every record whose genotype carries an ALT allele in at least one caller
becomes one row per (chrom, pos, ref, alt), with each caller's genotype and FILTER, the number of callers
that call ALT and heterozygous with FILTER PASS or '.', the longcallR phase set and RDS label, and VEP
fields from the annotated VCFs.

Stage ``--donor``: the union of one donor's library tables. The tier follows merge_callers.py:
high = ALT in >= 3 callers, medium = 2, low = 1. No row is dropped; every tier is kept and labelled.
The donor genotype itself is called later from the pooled reads (read_sites + allelic_tests), because
RNA callers see one library at a time.

Library -> donor goes through the BAM metadata (uid, bam, library) and a donor map (original_uid, donor).
The genome_variant_detection folder is taken from the BAM path in the BAM metadata, never from the discovery
manifest: the manifest carries a swapped WF86/WM63 BAM pairing.

Run:
  python -m altanalyze3.components.snp.het_panel --library --bam_metadata mds_metadata_bam_fixed.txt \\
      --bam_row 7 --donor_map donor_map.tsv --out PANEL_DIR
  python -m altanalyze3.components.snp.het_panel --donor AM75 --bam_metadata ... --donor_map ... --out PANEL_DIR
Outputs:
  PANEL_DIR/libraries/<row>__<library>.sites.tsv.gz
  PANEL_DIR/donors/<donor>.sites.tsv.gz
  PANEL_DIR/logs/*.parameters.json
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
import os
import sys
from collections import Counter, OrderedDict

CALLERS = OrderedDict([
    ("gatk", "results/haplotypecaller/variants.vcf.gz"),
    ("deepvariant", "results/deepvariant/variants.vcf.gz"),
    ("clair3", "results/clair3_rna/variants.vcf.gz"),
    ("longcallr", "results/longcallr/variants.vcf.gz"),
])
ANNOTATED = ("annotation/gatk/annotated.vcf.gz", "annotation/clair3/annotated.vcf.gz")
CSQ_FIELDS = ("SYMBOL", "Gene", "Consequence", "Existing_variation", "MAX_AF")
PASS_FILTERS = ("PASS", ".")
LIB_FIELDS = (["chrom", "pos", "ref", "alt", "vtype"]
              + ["{}_{}".format(c, k) for c in CALLERS for k in ("gt", "filter")]
              + ["n_alt_callers", "n_het_callers", "longcallr_ps", "longcallr_rds",
                 "symbol", "gene", "consequence", "existing_variation", "max_af", "somatic_class"])
DONOR_FIELDS = ["chrom", "pos", "ref", "alt", "vtype", "tier", "max_alt_callers", "max_het_callers",
                "n_libraries", "n_libraries_alt", "n_libraries_het", "libraries_alt", "longcallr_ps_by_library",
                "longcallr_rds", "symbol", "gene", "consequence", "existing_variation", "max_af", "somatic_class"]


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def normalize(pos, ref, alt):
    """Trim the shared suffix, then the shared prefix, keeping one anchor base (bcftools norm without
    left-alignment). Returns (pos, ref, alt)."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def vtype(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == 1 and len(alt) > 1 and alt[0] == ref[0]:
        return "INS"
    if len(alt) == 1 and len(ref) > 1 and ref[0] == alt[0]:
        return "DEL"
    return "COMPLEX"


def parse_gt(text):
    """'0/1', '1|0', '1/2', './.' -> (list of allele indices or None, phased)."""
    if not text or text in (".", "./.", ".|."):
        return None, False
    phased = "|" in text
    try:
        return [None if a == "." else int(a) for a in text.replace("|", "/").split("/")], phased
    except ValueError:
        return None, phased


def zygosity(alleles, index):
    """Zygosity of ALT allele `index` in a genotype: 'het' (index with REF), 'hom', 'multi' (with another ALT)."""
    if not alleles or index not in alleles:
        return None
    others = [a for a in alleles if a != index]
    if not others:
        return "hom"
    if all(a == 0 for a in others):
        return "het"
    return "multi"


def read_caller(path):
    """{(chrom,pos,ref,alt): dict(gt, filter, zyg, ps, rds)} for every ALT allele present in a caller's genotype."""
    out = {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            fmt = p[8].split(":")
            val = p[9].split(":")
            fv = dict(zip(fmt, val))
            alleles, _phased = parse_gt(fv.get("GT", ""))
            if not alleles:
                continue
            rds = ""
            if "RDS=" in p[7]:
                rds = p[7].split("RDS=", 1)[1].split(";", 1)[0]
            for i, alt in enumerate(p[4].split(","), start=1):
                z = zygosity(alleles, i)
                if z is None or alt in ("*", "<NON_REF>", "."):
                    continue
                pos, ref, alt2 = normalize(int(p[1]), p[3].upper(), alt.upper())
                out[(p[0], pos, ref, alt2)] = dict(gt=fv.get("GT", ""), filter=p[6], zyg=z,
                                                   ps=fv.get("PS", "") if fv.get("PS", ".") != "." else "", rds=rds)
    return out


def read_annotation(path, keys):
    """VEP fields and SOMATIC_CLASS for the given site keys from an annotated VCF."""
    out = {}
    idx = None
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                fmt = line.split("Format: ", 1)[1].split('"', 1)[0].split("|")
                idx = {f: fmt.index(f) for f in CSQ_FIELDS if f in fmt}
                allele_i = fmt.index("Allele")
                continue
            if line.startswith("#"):
                continue
            p = line.split("\t", 8)
            info = p[7]
            for alt in p[4].split(","):
                key = (p[0],) + normalize(int(p[1]), p[3].upper(), alt.upper())
                if key not in keys or key in out:
                    continue
                row = {"somatic_class": ""}
                if "SOMATIC_CLASS=" in info:
                    row["somatic_class"] = info.split("SOMATIC_CLASS=", 1)[1].split(";", 1)[0]
                if idx is not None and "CSQ=" in info:
                    csq = info.split("CSQ=", 1)[1].split(";", 1)[0].split(",")
                    # VEP writes the allele without the shared anchor base for indels; take the first entry
                    # whose Allele matches either form, else the first entry.
                    vep_allele = alt[1:] if len(alt) != len(p[3]) and alt[:1] == p[3][:1] else alt
                    vep_allele = vep_allele or "-"
                    pick = next((c for c in csq if c.split("|")[allele_i] in (alt, vep_allele)), csq[0]).split("|")
                    for f, k in (("SYMBOL", "symbol"), ("Gene", "gene"), ("Consequence", "consequence"),
                                 ("Existing_variation", "existing_variation"), ("MAX_AF", "max_af")):
                        row[k] = pick[idx[f]] if f in idx and idx[f] < len(pick) else ""
                out[key] = row
    return out


def load_bam_metadata(path):
    with open(path) as fh:
        rows = list(csv.DictReader((l.replace("\r", "") for l in fh), delimiter="\t"))
    for i, r in enumerate(rows):
        r["row"] = i
        r["gvd"] = os.path.join(os.path.dirname(r["bam"]), "genome_variant_detection")
    return rows


def load_donor_map(path):
    """original_uid -> donor. One original uid must map to one donor."""
    m = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            u, d = r["original_uid"], r["donor"]
            if m.get(u, d) != d:
                raise SystemExit("original_uid {} maps to two donors: {} and {}".format(u, m[u], d))
            m[u] = d
    return m


def library_tag(r):
    return "{:02d}__{}".format(r["row"], r["library"])


def run_library(a):
    meta = load_bam_metadata(a.bam_metadata)
    r = meta[a.bam_row]
    dmap = load_donor_map(a.donor_map)
    if r["uid"] not in dmap:
        raise SystemExit("BAM metadata uid {} has no donor in {}".format(r["uid"], a.donor_map))
    calls, used, missing = {}, {}, []
    for caller, sub in CALLERS.items():
        path = os.path.join(r["gvd"], sub)
        if not os.path.exists(path):
            missing.append(caller)
            continue
        used[caller] = dict(path=path, size=os.path.getsize(path), mtime=os.path.getmtime(path))
        calls[caller] = read_caller(path)
        print("[het_panel] {} {}: {:,} ALT alleles".format(library_tag(r), caller, len(calls[caller])), flush=True)
    if not calls:
        raise SystemExit("no caller VCF under {}".format(r["gvd"]))
    keys = set()
    for c in calls.values():
        keys.update(c)
    ann = {}
    if not a.no_annotation:
        for sub in ANNOTATED:
            path = os.path.join(r["gvd"], sub)
            if os.path.exists(path):
                todo = keys - set(ann)
                ann.update(read_annotation(path, todo))
                used[sub] = dict(path=path, size=os.path.getsize(path), mtime=os.path.getmtime(path))
    os.makedirs(os.path.join(a.out, "libraries"), exist_ok=True)
    out_path = os.path.join(a.out, "libraries", library_tag(r) + ".sites.tsv.gz")
    tmp = out_path + ".tmp"
    tiers = Counter()
    with gzip.open(tmp, "wt") as fh:
        w = csv.DictWriter(fh, fieldnames=LIB_FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for key in sorted(keys, key=lambda k: (k[0], k[1], k[2], k[3])):
            row = dict(chrom=key[0], pos=key[1], ref=key[2], alt=key[3], vtype=vtype(key[2], key[3]))
            n_alt = n_het = 0
            for caller in CALLERS:
                c = calls.get(caller, {}).get(key)
                row[caller + "_gt"] = c["gt"] if c else ("NA" if caller not in calls else "")
                row[caller + "_filter"] = c["filter"] if c else ""
                if c and c["filter"] in PASS_FILTERS:
                    n_alt += 1
                    n_het += c["zyg"] == "het"
            lc = calls.get("longcallr", {}).get(key)
            row.update(n_alt_callers=n_alt, n_het_callers=n_het,
                       longcallr_ps=lc["ps"] if lc else "", longcallr_rds=lc["rds"] if lc else "")
            an = ann.get(key, {})
            for k in ("symbol", "gene", "consequence", "existing_variation", "max_af", "somatic_class"):
                row[k] = an.get(k, "")
            w.writerow(row)
            tiers[n_alt] += 1
    os.replace(tmp, out_path)
    log = dict(module="altanalyze3.components.snp.het_panel", stage="library", argv=sys.argv,
               finished=datetime.datetime.now().isoformat(timespec="seconds"),
               bam_metadata=a.bam_metadata, bam_row=a.bam_row, uid=r["uid"], library=r["library"], bam=r["bam"],
               donor=dmap[r["uid"]], inputs=used, callers_missing=missing, sites=len(keys),
               sites_by_n_pass_alt_callers={str(k): v for k, v in sorted(tiers.items())},
               annotated_sites=len(ann), output=out_path, module_sha256=_sha256(os.path.abspath(__file__)))
    os.makedirs(os.path.join(a.out, "logs"), exist_ok=True)
    with open(os.path.join(a.out, "logs", library_tag(r) + ".parameters.json"), "w") as fh:
        json.dump(log, fh, indent=2)
    print("[het_panel] wrote {:,} sites -> {} ({})".format(len(keys), out_path, dict(tiers)), flush=True)


def tier_of(n):
    return "high" if n >= 3 else ("medium" if n == 2 else ("low" if n == 1 else "none_pass"))


def run_donor(a):
    meta = load_bam_metadata(a.bam_metadata)
    dmap = load_donor_map(a.donor_map)
    libs = [r for r in meta if dmap.get(r["uid"]) == a.donor]
    if not libs:
        raise SystemExit("no BAM metadata row belongs to donor {}".format(a.donor))
    merged, present, absent = {}, [], []
    for r in libs:
        path = os.path.join(a.out, "libraries", library_tag(r) + ".sites.tsv.gz")
        if not os.path.exists(path):
            absent.append(library_tag(r))
            continue
        present.append(library_tag(r))
        with gzip.open(path, "rt") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                key = (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
                m = merged.get(key)
                if m is None:
                    m = merged[key] = dict(chrom=key[0], pos=key[1], ref=key[2], alt=key[3], vtype=row["vtype"],
                                           max_alt_callers=0, max_het_callers=0, n_libraries_alt=0,
                                           n_libraries_het=0, libraries_alt=[], ps=[], rds=set(), ann={})
                na, nh = int(row["n_alt_callers"]), int(row["n_het_callers"])
                m["max_alt_callers"] = max(m["max_alt_callers"], na)
                m["max_het_callers"] = max(m["max_het_callers"], nh)
                m["n_libraries_alt"] += na > 0
                m["n_libraries_het"] += nh > 0
                m["libraries_alt"].append(library_tag(r))
                if row["longcallr_ps"]:
                    m["ps"].append("{}={}".format(library_tag(r), row["longcallr_ps"]))
                if row["longcallr_rds"]:
                    m["rds"].add(row["longcallr_rds"])
                if not m["ann"] and row["symbol"] + row["existing_variation"] + row["consequence"]:
                    m["ann"] = {k: row[k] for k in ("symbol", "gene", "consequence", "existing_variation",
                                                    "max_af", "somatic_class")}
    if not present:
        raise SystemExit("donor {}: no library panel written yet ({} expected)".format(a.donor, len(libs)))
    os.makedirs(os.path.join(a.out, "donors"), exist_ok=True)
    out_path = os.path.join(a.out, "donors", a.donor + ".sites.tsv.gz")
    tmp = out_path + ".tmp"
    tiers = Counter()
    with gzip.open(tmp, "wt") as fh:
        w = csv.DictWriter(fh, fieldnames=DONOR_FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for key in sorted(merged):
            m = merged[key]
            t = tier_of(m["max_alt_callers"])
            tiers[t] += 1
            row = dict(chrom=m["chrom"], pos=m["pos"], ref=m["ref"], alt=m["alt"], vtype=m["vtype"], tier=t,
                       max_alt_callers=m["max_alt_callers"], max_het_callers=m["max_het_callers"],
                       n_libraries=len(present), n_libraries_alt=m["n_libraries_alt"],
                       n_libraries_het=m["n_libraries_het"], libraries_alt=",".join(m["libraries_alt"]),
                       longcallr_ps_by_library=",".join(m["ps"]), longcallr_rds=",".join(sorted(m["rds"])))
            for k in ("symbol", "gene", "consequence", "existing_variation", "max_af", "somatic_class"):
                row[k] = m["ann"].get(k, "")
            w.writerow(row)
    os.replace(tmp, out_path)
    log = dict(module="altanalyze3.components.snp.het_panel", stage="donor", argv=sys.argv, donor=a.donor,
               finished=datetime.datetime.now().isoformat(timespec="seconds"),
               libraries_in_metadata=[library_tag(r) for r in libs], libraries_used=present,
               libraries_without_panel=absent, sites=len(merged), sites_by_tier=dict(tiers), output=out_path,
               module_sha256=_sha256(os.path.abspath(__file__)))
    os.makedirs(os.path.join(a.out, "logs"), exist_ok=True)
    with open(os.path.join(a.out, "logs", "donor__{}.parameters.json".format(a.donor)), "w") as fh:
        json.dump(log, fh, indent=2)
    print("[het_panel] donor {}: {:,} sites from {} of {} libraries ({}) -> {}".format(
        a.donor, len(merged), len(present), len(libs), dict(tiers), out_path), flush=True)


def run_concordance(a):
    """Donor identity by measurement over every pair of library tables, on autosomal SNVs that at least
    --concordance_callers callers call ALT. Two overlap measures: Jaccard |A&B|/|A|B| (depth-sensitive) and
    containment |A&B|/min(|A|,|B|) (robust to a shallow library). Each library's best match by each measure
    must belong to its mapped donor."""
    meta = load_bam_metadata(a.bam_metadata)
    dmap = load_donor_map(a.donor_map)
    sets = {}
    for r in meta:
        path = os.path.join(a.out, "libraries", library_tag(r) + ".sites.tsv.gz")
        if not os.path.exists(path):
            continue
        s = set()
        with gzip.open(path, "rt") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            ic, ip, ia, iv, ina = (head.index(k) for k in ("chrom", "pos", "alt", "vtype", "n_alt_callers"))
            for line in fh:
                p = line.split("\t")
                c = p[ic][3:] if p[ic].startswith("chr") else p[ic]
                if p[iv] == "SNV" and c.isdigit() and int(p[ina]) >= a.concordance_callers:
                    s.add((c, p[ip], p[ia]))
        sets[library_tag(r)] = (s, dmap.get(r["uid"], ""))
    tags = sorted(sets)
    out_dir = a.concordance_out or os.path.join(a.out, "concordance")
    os.makedirs(out_dir, exist_ok=True)
    n = len(tags)
    jac = [[float("nan")] * n for _ in range(n)]
    con = [[float("nan")] * n for _ in range(n)]
    for i in range(n):
        si = sets[tags[i]][0]
        for j in range(i, n):
            sj = sets[tags[j]][0]
            inter = len(si & sj)
            uni = len(si) + len(sj) - inter
            jac[i][j] = jac[j][i] = inter / uni if uni else float("nan")
            m = min(len(si), len(sj))
            con[i][j] = con[j][i] = inter / m if m else float("nan")
    for name, mat in (("library_jaccard.tsv", jac), ("library_containment.tsv", con)):
        with open(os.path.join(out_dir, name), "w") as fh:
            fh.write("library\tdonor\tn_sites\t" + "\t".join(tags) + "\n")
            for i, t in enumerate(tags):
                fh.write("{}\t{}\t{}\t{}\n".format(t, sets[t][1], len(sets[t][0]),
                                                  "\t".join("{:.4f}".format(x) for x in mat[i])))
    best_rows = []
    for i, t in enumerate(tags):
        d = sets[t][1]
        row = dict(library=t, donor=d, n_sites=len(sets[t][0]))
        for name, mat in (("jaccard", jac), ("containment", con)):
            others = [(mat[i][j], tags[j]) for j in range(n) if j != i and mat[i][j] == mat[i][j]]
            same = [x for x, u in others if sets[u][1] == d]
            diff = [x for x, u in others if sets[u][1] != d]
            bx, bu = max(others) if others else (float("nan"), "")
            row.update({"best_" + name: round(bx, 4), "best_" + name + "_library": bu,
                        "best_" + name + "_donor": sets[bu][1] if bu else "",
                        "max_same_donor_" + name: round(max(same), 4) if same else "",
                        "max_other_donor_" + name: round(max(diff), 4) if diff else "",
                        "status_" + name: ("no_sites" if not sets[t][0] else ("single_library_donor" if not same else
                                           ("PASS" if (not diff or max(same) > max(diff)) else "MISMATCH")))})
        best_rows.append(row)
    bpath = os.path.join(out_dir, "library_best_match.tsv")
    with open(bpath, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(best_rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(best_rows)
    for name in ("jaccard", "containment"):
        st = Counter(r["status_" + name] for r in best_rows)
        print("[het_panel] concordance ({}, >= {} callers) over {} libraries: {}".format(
            name, a.concordance_callers, len(tags), dict(st)), flush=True)
        bad = [(r["library"], r["donor"], r["best_" + name + "_library"], r["best_" + name + "_donor"])
               for r in best_rows if r["status_" + name] == "MISMATCH"]
        if bad:
            print("[het_panel] {} MISMATCH: {}".format(name, bad), flush=True)


def build_arg_parser():
    ap = argparse.ArgumentParser(prog="python -m altanalyze3.components.snp.het_panel",
                                 description="Candidate germline variant panel per library and per donor.")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--library", action="store_true", help="build one library table (needs --bam_row)")
    mode.add_argument("--donor", help="merge the library tables of this donor")
    mode.add_argument("--concordance", action="store_true",
                      help="pairwise Jaccard of library SNV sets; checks the donor map by measurement")
    ap.add_argument("--concordance_callers", type=int, default=3,
                    help="callers that must call ALT for a site to enter the concordance (default 3 = high tier)")
    ap.add_argument("--concordance_out", default=None, help="concordance folder (default <out>/concordance)")
    ap.add_argument("--bam_metadata", required=True, help="workflow BAM metadata (uid, bam, library)")
    ap.add_argument("--bam_row", type=int, help="0-based data row of --bam_metadata")
    ap.add_argument("--donor_map", required=True, help="TSV with original_uid and donor columns")
    ap.add_argument("--no_annotation", action="store_true", help="skip the VEP fields of the annotated VCFs")
    ap.add_argument("--out", required=True)
    return ap


def main(argv=None):
    a = build_arg_parser().parse_args(argv)
    if a.library:
        if a.bam_row is None:
            raise SystemExit("--library needs --bam_row")
        run_library(a)
    elif a.concordance:
        run_concordance(a)
    else:
        run_donor(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
