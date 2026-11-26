#!/usr/bin/env python3
"""Read-backed phasing of one donor's heterozygous sites on one chromosome, and a haplotype label for every
read that carries at least one of them.

Why: a read that skips the exon holding a SNP carries no allele at that SNP, so a per-site allele x isoform
test can never see the skipping reads (exon skipping driven by a SNP inside the skipped exon, as CD33
rs12459419, is invisible to it). Full-length reads usually carry several heterozygous sites; phasing them
gives the skipping read the haplotype, and so the allele, of the SNP it does not cover.

Method. For every pair of heterozygous sites observed on the same primary read, count cis reads (both REF or
both ALT) and trans reads. A pair is an edge when cis != trans; its orientation is "same phase" when cis >
trans, its weight |cis - trans|. Edges are taken by decreasing weight into a union-find with parity (Kruskal
maximum spanning forest); an edge that closes a cycle with the opposite parity is counted as a conflict, not
used. Each connected component is a phase block; each site gets a phase bit (0: REF on haplotype 1).
A read votes per site (REF at a phase-0 site -> haplotype 1, ...); it takes the majority haplotype of its
votes within a block, and the counts of agreeing and disagreeing votes are written beside it.

Inputs: read_sites reads tables, and allelic_tests donor_site_genotypes.tsv.gz for the genotype calls.
Outputs in OUT/<chrom>/: <donor>.sites.tsv.gz, <donor>.reads.tsv.gz, <donor>.edges.tsv.gz,
logs/<donor>.parameters.json
Run: python -m altanalyze3.components.snp.phase_haplotypes --chrom chr7 --donor AML7 --reads READS \\
        --genotypes STATS/chr7/donor_site_genotypes.tsv.gz --donor_map donor_map.tsv --out PHASE
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
import os
import sys
import time
from collections import Counter, defaultdict
from itertools import combinations


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class ParityUnionFind:
    """Union-find where each node carries its parity relative to the root."""

    def __init__(self):
        self.parent, self.parity, self.rank = {}, {}, {}

    def add(self, x):
        if x not in self.parent:
            self.parent[x], self.parity[x], self.rank[x] = x, 0, 0

    def find(self, x):
        path = []
        while self.parent[x] != x:
            path.append(x)
            x = self.parent[x]
        root = x
        # path compression with parity accumulation, nearest-to-root first
        acc = 0
        for node in reversed(path):
            acc ^= self.parity[node]
            self.parity[node] = acc
            self.parent[node] = root
        return root

    def parity_of(self, x):
        self.find(x)
        return self.parity[x] if self.parent[x] != x else 0

    def union(self, x, y, diff):
        """Join x and y with parity(x) xor parity(y) == diff. Returns 'joined', 'consistent' or 'conflict'."""
        rx, ry = self.find(x), self.find(y)
        px, py = self.parity_of(x), self.parity_of(y)
        if rx == ry:
            return "consistent" if (px ^ py) == diff else "conflict"
        if self.rank[rx] < self.rank[ry]:
            rx, ry, px, py = ry, rx, py, px
        self.parent[ry] = rx
        self.parity[ry] = px ^ py ^ diff
        if self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1
        return "joined"


def phase(read_calls):
    """read_calls: {read_key: [(site_id, is_ref), ...]} with >= 1 call each.
    Returns (site_phase {site: (block, bit)}, edges list, read_haps {read: [(block, hap, agree, disagree)]})."""
    pair = defaultdict(lambda: [0, 0])          # (s1, s2) -> [cis, trans]
    for calls in read_calls.values():
        if len(calls) < 2:
            continue
        seen = {}
        for s, r in calls:
            seen[s] = r                          # one call per site per read
        items = sorted(seen.items())
        for (s1, r1), (s2, r2) in combinations(items, 2):
            pair[(s1, s2)][0 if r1 == r2 else 1] += 1
    edges = []
    for (s1, s2), (cis, trans) in pair.items():
        edges.append(dict(site1=s1, site2=s2, cis=cis, trans=trans, weight=abs(cis - trans),
                          orientation="same" if cis > trans else ("opposite" if trans > cis else "tie")))
    uf = ParityUnionFind()
    for calls in read_calls.values():
        for s, _ in calls:
            uf.add(s)
    for e in sorted(edges, key=lambda e: (-e["weight"], e["site1"], e["site2"])):
        if e["orientation"] == "tie":
            e["status"] = "tie_unused"
            continue
        e["status"] = uf.union(e["site1"], e["site2"], 0 if e["orientation"] == "same" else 1)
    site_phase = {s: (uf.find(s), uf.parity_of(s)) for s in uf.parent}
    read_haps = {}
    for key, calls in read_calls.items():
        votes = defaultdict(Counter)
        for s, is_ref in calls:
            block, bit = site_phase[s]
            hap = bit if is_ref else 1 - bit     # REF at a phase-0 site -> haplotype 0
            votes[block][hap] += 1
        out = []
        for block, c in votes.items():
            h0, h1 = c.get(0, 0), c.get(1, 0)
            hap = 0 if h0 > h1 else (1 if h1 > h0 else -1)
            out.append((block, hap, max(h0, h1), min(h0, h1)))
        read_haps[key] = out
    return site_phase, edges, read_haps


def load_het(genotypes, donor, tiers=None):
    het = set()
    with gzip.open(genotypes, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["donor"] == donor and r["genotype"] == "0/1" and (tiers is None or r["tier"] in tiers):
                het.add(r["site_id"])
    return het


def build_arg_parser():
    ap = argparse.ArgumentParser(prog="python -m altanalyze3.components.snp.phase_haplotypes",
                                 description="Read-backed phase blocks and per-read haplotypes for one donor.")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--donor", required=True)
    ap.add_argument("--tiers", default=None, help="comma list of panel tiers whose heterozygous sites are phased")
    ap.add_argument("--reads", required=True, help="read_sites output folder")
    ap.add_argument("--genotypes", required=True, help="allelic_tests donor_site_genotypes.tsv.gz of the chromosome")
    ap.add_argument("--donor_map", required=True, help="TSV with uid (ISV uid) and donor")
    ap.add_argument("--panel_dir", default=None,
                    help="het_panel output: compare the phase with longcallR phase sets of the donor's libraries")
    ap.add_argument("--bam_metadata", default=None, help="BAM metadata of het_panel (needed with --panel_dir)")
    ap.add_argument("--out", required=True)
    return ap


def longcallr_check(panel_dir, bam_metadata, donor_map, donor, chrom, site_phase):
    """Concordance of the relative phase of adjacent site pairs that longcallR placed in one phase set (GT 0|1
    or 1|0) and this module placed in one block. longcallR bit: 0 for 0|1 (REF on haplotype 1), 1 for 1|0."""
    from altanalyze3.components.snp.het_panel import library_tag, load_bam_metadata, load_donor_map
    dmap = load_donor_map(donor_map)
    rows = []
    for r in load_bam_metadata(bam_metadata):
        if dmap.get(r["uid"]) != donor:
            continue
        p = os.path.join(panel_dir, "libraries", library_tag(r) + ".sites.tsv.gz")
        if not os.path.exists(p):
            continue
        by_ps = defaultdict(list)
        with gzip.open(p, "rt") as fh:
            for x in csv.DictReader(fh, delimiter="\t"):
                if x["chrom"] != chrom or not x["longcallr_ps"] or x["longcallr_gt"] not in ("0|1", "1|0"):
                    continue
                sid = "{}:{}:{}:{}".format(x["chrom"], x["pos"], x["ref"], x["alt"])
                if sid in site_phase:
                    by_ps[x["longcallr_ps"]].append((int(x["pos"]), sid, 0 if x["longcallr_gt"] == "0|1" else 1))
        conc = disc = other_block = 0
        for ps, v in by_ps.items():
            v.sort()
            for (p1, s1, l1), (p2, s2, l2) in zip(v, v[1:]):
                (b1, o1), (b2, o2) = site_phase[s1], site_phase[s2]
                if b1 != b2:
                    other_block += 1
                elif (o1 ^ o2) == (l1 ^ l2):
                    conc += 1
                else:
                    disc += 1
        rows.append(dict(library=library_tag(r), adjacent_pairs_same_block=conc + disc, concordant=conc,
                         discordant=disc, pairs_in_different_blocks=other_block,
                         concordance=round(conc / (conc + disc), 4) if conc + disc else None))
    return rows


def main(argv=None):
    a = build_arg_parser().parse_args(argv)
    t0 = time.time()
    tiers = set(a.tiers.split(",")) if a.tiers else None
    het = load_het(a.genotypes, a.donor, tiers)
    uids = []
    with open(a.donor_map) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["donor"] == a.donor:
                uids.append(r["uid"])
    if not uids:
        raise SystemExit("donor {} has no library in {}".format(a.donor, a.donor_map))
    read_calls, meta, n = {}, {}, Counter()
    for uid in uids:
        p = os.path.join(a.reads, uid, "{}.reads.tsv.gz".format(a.chrom))
        if not os.path.exists(p):
            raise SystemExit("missing {}".format(p))
        with gzip.open(p, "rt") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["secondary"] != "0" or row["supplementary"] != "0":
                    continue
                calls = []
                for c in row["calls"].split(";"):
                    sid, cls = c.rsplit("=", 1)
                    full = a.chrom + ":" + sid
                    if cls in ("REF", "ALT") and full in het:
                        calls.append((full, cls == "REF"))
                if not calls:
                    continue
                key = row["molecule_id"]
                if key in read_calls:
                    raise SystemExit("molecule id {} appears twice among primary alignments".format(key))
                read_calls[key] = calls
                meta[key] = (uid, row["cb_rc"], row["index_rows"])
                n["reads"] += 1
                n["reads_multi_site"] += len({s for s, _ in calls}) > 1
    site_phase, edges, read_haps = phase(read_calls)
    out = os.path.join(a.out, a.chrom)
    os.makedirs(os.path.join(out, "logs"), exist_ok=True)
    blocks = Counter(b for b, _ in site_phase.values())
    with gzip.open(os.path.join(out, a.donor + ".sites.tsv.gz"), "wt") as fh:
        fh.write("donor\tsite_id\tblock\tphase_bit\tblock_sites\n")
        for s in sorted(site_phase, key=lambda x: int(x.split(":")[1])):
            b, bit = site_phase[s]
            fh.write("{}\t{}\t{}\t{}\t{}\n".format(a.donor, s, b, bit, blocks[b]))
    with gzip.open(os.path.join(out, a.donor + ".edges.tsv.gz"), "wt") as fh:
        cols = ["site1", "site2", "cis", "trans", "weight", "orientation", "status"]
        fh.write("\t".join(cols) + "\n")
        for e in edges:
            fh.write("\t".join(str(e[c]) for c in cols) + "\n")
    hap_counts = Counter()
    with gzip.open(os.path.join(out, a.donor + ".reads.tsv.gz"), "wt") as fh:
        fh.write("molecule_id\tuid\tcb_rc\tblock\thaplotype\tvotes_agree\tvotes_disagree\tindex_rows\n")
        for key, hs in read_haps.items():
            uid, cb, idx = meta[key]
            for b, h, ag, dis in hs:
                hap_counts["ambiguous" if h < 0 else ("conflicting_votes" if dis else "clean")] += 1
                fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(key, uid, cb, b, h, ag, dis, idx))
    st = Counter(e["status"] for e in edges)
    lcr = None
    if a.panel_dir:
        if not a.bam_metadata:
            raise SystemExit("--panel_dir needs --bam_metadata")
        lcr = longcallr_check(a.panel_dir, a.bam_metadata, a.donor_map, a.donor, a.chrom, site_phase)
    log = dict(module="altanalyze3.components.snp.phase_haplotypes", argv=sys.argv if argv is None else list(argv),
               donor=a.donor, chrom=a.chrom, tiers=a.tiers, libraries=uids, het_sites=len(het), het_sites_with_reads=len(site_phase),
               blocks=len(blocks), blocks_with_2plus_sites=sum(1 for v in blocks.values() if v > 1),
               sites_in_multi_site_blocks=sum(v for v in blocks.values() if v > 1), edges=len(edges),
               edge_status=dict(st), read_block_assignments=dict(hap_counts), longcallr_phase_check=lcr, **n,
               seconds=round(time.time() - t0, 1), finished=datetime.datetime.now().isoformat(timespec="seconds"),
               module_sha256=_sha256(os.path.abspath(__file__)))
    with open(os.path.join(out, "logs", a.donor + ".parameters.json"), "w") as fh:
        json.dump(log, fh, indent=2)
    print("[phase_haplotypes] {} {}: {:,} het sites with reads in {:,} blocks ({:,} multi-site); edges {}; reads {}".format(
        a.donor, a.chrom, len(site_phase), len(blocks), log["blocks_with_2plus_sites"], dict(st), dict(hap_counts)),
        flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
