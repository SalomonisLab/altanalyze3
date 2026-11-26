"""Synthetic checks of altanalyze3.components.snp.read_sites and het_panel: allele classes per read, indel
matching, site-count partition, molecule-index join, multi-allelic normalisation and tiering."""
import csv
import gzip
import json
import os
import sqlite3

import pysam
import pytest

from altanalyze3.components.snp import het_panel, read_sites

SEQ = ("ACGTTGCAAGCTTGACCGTAGGCTAACGTTAGCCATGGCATTCAGGTACCGATTGCAAGTCCAGTTAGCAACG" * 4)[:300]


def _fasta(tmp):
    p = os.path.join(tmp, "g.fa")
    with open(p, "w") as fh:
        fh.write(">chr7\n" + SEQ + "\n")
    pysam.faidx(p)
    return p


def _read(name, start0, cigar, seq, cb="AAAACCCCGGGGTTTT", flag=0):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.flag = flag
    a.reference_id = 0
    a.reference_start = start0
    a.mapping_quality = 60
    a.cigarstring = cigar
    a.query_sequence = seq
    a.set_tag("CB", cb)
    return a


def _bam(tmp, reads):
    p = os.path.join(tmp, "x.bam")
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr7", "LN": len(SEQ)}]}
    with pysam.AlignmentFile(p + ".u.bam", "wb", header=hdr) as out:
        for r in reads:
            out.write(r)
    pysam.sort("-o", p, p + ".u.bam")
    pysam.index(p)
    return p


def _mutate(seq, pos1, base):
    return seq[:pos1 - 1] + base + seq[pos1:]


@pytest.fixture
def world(tmp_path):
    tmp = str(tmp_path)
    fa = _fasta(tmp)
    snv = 50                                   # 1-based SNV position
    alt_base = "T" if SEQ[snv - 1] != "T" else "G"
    ref_base = SEQ[snv - 1]
    dpos, dlen = 120, 2                        # DEL: anchor at 120, deletes 121-122
    ipos, iseq = 200, "GG"                     # INS after anchor 200
    reads = []
    # r1: REF at the SNV, REF across the DEL
    reads.append(_read("molecule/1", 30, "120M", SEQ[30:150]))
    # r2: ALT at the SNV, carries the 2-bp deletion
    s2 = _mutate(SEQ, snv, alt_base)
    reads.append(_read("molecule/2", 30, "{}M{}D{}M".format(dpos - 30, dlen, 150 - dpos - dlen),
                       s2[30:dpos] + s2[dpos + dlen:150]))
    # r3: splices over the SNV (N 40..69 0-based covers pos 50)
    reads.append(_read("molecule/3", 20, "20M30N40M", SEQ[20:40] + SEQ[70:110]))
    # r4: carries the insertion after 200
    reads.append(_read("molecule/4", 180, "20M2I30M", SEQ[180:200] + iseq + SEQ[200:230]))
    # r5: REF across the insertion site
    reads.append(_read("molecule/5", 180, "50M", SEQ[180:230]))
    # r6: secondary alignment, ALT at the SNV
    reads.append(_read("molecule/6", 40, "30M", s2[40:70], flag=256))
    # r7: a 1-bp deletion near the DEL site: OTHER
    reads.append(_read("molecule/7", 100, "{}M1D{}M".format(dpos - 100, 140 - dpos - 1), SEQ[100:dpos] + SEQ[dpos + 1:140]))
    bam = _bam(tmp, reads)
    panel_dir = os.path.join(tmp, "panel")
    os.makedirs(os.path.join(panel_dir, "donors"))
    with gzip.open(os.path.join(panel_dir, "donors", "D1.sites.tsv.gz"), "wt") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(het_panel.DONOR_FIELDS)
        for pos, ref, alt, vt in ((snv, ref_base, alt_base, "SNV"),
                                  (dpos, SEQ[dpos - 1:dpos + dlen], SEQ[dpos - 1], "DEL"),
                                  (ipos, SEQ[ipos - 1], SEQ[ipos - 1] + iseq, "INS")):
            w.writerow(["chr7", pos, ref, alt, vt, "high", 3, 3] + [""] * (len(het_panel.DONOR_FIELDS) - 8))
    bam_meta = os.path.join(tmp, "bam_meta.txt")
    with open(bam_meta, "w") as fh:
        fh.write("uid\tbam\tlibrary\treverse\tgroups\nU1\t{}\tLIB1\tTRUE\tyoung\n".format(bam))
    isv_meta = os.path.join(tmp, "isv_meta.tsv")
    with open(isv_meta, "w") as fh:
        fh.write("uid\toriginal_uid\toriginal_library\toriginal_metadata_ordinal\nLIB1\tU1\tLIB1\t0\n")
    dmap = os.path.join(tmp, "donor_map.tsv")
    with open(dmap, "w") as fh:
        fh.write("uid\toriginal_uid\tdonor\nLIB1\tU1\tD1\n")
    mol = os.path.join(tmp, "mol")
    os.makedirs(mol)
    con = sqlite3.connect(os.path.join(mol, "LIB1.reads_final.db"))
    con.execute("CREATE TABLE rf (gene TEXT NOT NULL, mol TEXT NOT NULL, ct TEXT NOT NULL, cnt REAL, final TEXT, "
                "PRIMARY KEY (gene, mol, ct)) WITHOUT ROWID")
    con.executemany("INSERT INTO rf VALUES (?,?,?,?,?)",
                    [("ENSG1", "1", "HSC-1", 1.0, "ENST_A"), ("ENSG1", "2", "HSC-1", 1.0, "ENST_B"),
                     ("ENSG1", "99", "HSC-1", 1.0, "ENST_A")])
    con.commit(); con.close()
    gm = os.path.join(tmp, "exon.txt")
    with open(gm, "w") as fh:
        fh.write("gene\texon-id\tchromosome\tstrand\texon-region-start(s)\texon-region-stop(s)\n"
                 "ENSG1\tE1.1\tchr7\t+\t1\t300\nENSG2\tE1.1\tchr8\t+\t1\t300\n")
    return dict(tmp=tmp, fa=fa, bam=bam, panel_dir=panel_dir, bam_meta=bam_meta, isv_meta=isv_meta, dmap=dmap,
                mol=mol, gm=gm, snv=snv, ref_base=ref_base, alt_base=alt_base, dpos=dpos, ipos=ipos)


def _run(w):
    out = os.path.join(w["tmp"], "out")
    read_sites.main(["--uid", "LIB1", "--chrom", "chr7", "--isv_metadata", w["isv_meta"], "--bam_metadata",
                     w["bam_meta"], "--donor_map", w["dmap"], "--panel_dir", w["panel_dir"], "--molecule_index",
                     w["mol"], "--gene_model", w["gm"], "--fasta", w["fa"], "--out", out])
    with gzip.open(os.path.join(out, "LIB1", "chr7.reads.tsv.gz"), "rt") as fh:
        reads = {r["read_name"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with gzip.open(os.path.join(out, "LIB1", "chr7.site_counts.tsv.gz"), "rt") as fh:
        counts = {r["vtype"]: r for r in csv.DictReader(fh, delimiter="\t")}
    log = json.load(open(os.path.join(out, "LIB1", "logs", "chr7.parameters.json")))
    return reads, counts, log


def test_allele_classes_and_join(world):
    reads, counts, log = _run(world)
    w = world
    snv = "{}:{}:{}=".format(w["snv"], w["ref_base"], w["alt_base"])
    assert snv + "REF" in reads["molecule/1"]["calls"]
    assert snv + "ALT" in reads["molecule/2"]["calls"]
    assert "molecule/3" not in reads or snv + "SPLICE" in reads["molecule/3"]["calls"]
    assert snv + "ALT" in reads["molecule/6"]["calls"] and reads["molecule/6"]["secondary"] == "1"
    dcall = [c for c in reads["molecule/2"]["calls"].split(";") if c.startswith(str(w["dpos"]) + ":")][0]
    assert dcall.endswith("=ALT")
    rcall = [c for c in reads["molecule/1"]["calls"].split(";") if c.startswith(str(w["dpos"]) + ":")][0]
    assert rcall.endswith("=REF")
    ocall = [c for c in reads["molecule/7"]["calls"].split(";") if c.startswith(str(w["dpos"]) + ":")][0]
    assert ocall.endswith("=OTHER")
    assert reads["molecule/4"]["calls"].endswith("=ALT")
    assert reads["molecule/5"]["calls"].endswith("=REF")
    # SNV counts: all alignments REF 1 ALT 2 SPLICE 1; primary ALT 1
    c = counts["SNV"]
    assert (c["all_REF"], c["all_ALT"], c["all_SPLICE"], c["primary_ALT"]) == ("1", "2", "1", "1")
    # join: molecules 1 and 2 are in the index; molecule id carries the uid
    assert reads["molecule/1"]["index_rows"] == "ENSG1|HSC-1|ENST_A"
    assert reads["molecule/2"]["index_rows"] == "ENSG1|HSC-1|ENST_B"
    assert reads["molecule/4"]["index_rows"] == ""
    assert reads["molecule/1"]["molecule_id"] == "1.LIB1"
    assert reads["molecule/1"]["cb_rc"] == "AAAACCCCGGGGTTTT"[::-1].translate(str.maketrans("ACGT", "TGCA"))
    assert log["reads_joined_to_index"] == 2 and log["genes_queried"] == 1


def test_fasta_mismatch_stops(world):
    p = os.path.join(world["panel_dir"], "donors", "D1.sites.tsv.gz")
    with gzip.open(p, "rt") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    wrong = [b for b in "ACGT" if b != world["ref_base"] and b != world["alt_base"]][0]
    rows[1][2] = wrong
    with gzip.open(p, "wt") as fh:
        csv.writer(fh, delimiter="\t").writerows(rows)
    with pytest.raises(SystemExit):
        _run(world)


def test_normalize_and_zygosity():
    assert het_panel.normalize(10, "AT", "ATT") == (10, "A", "AT")
    assert het_panel.normalize(10, "AT", "A") == (10, "AT", "A")
    assert het_panel.normalize(10, "CAT", "CGT") == (11, "A", "G")
    assert het_panel.zygosity([0, 1], 1) == "het"
    assert het_panel.zygosity([1, 1], 1) == "hom"
    assert het_panel.zygosity([1, 2], 2) == "multi"
    assert het_panel.zygosity([0, 0], 1) is None
    assert het_panel.vtype("A", "G") == "SNV" and het_panel.vtype("A", "AT") == "INS"
    assert het_panel.vtype("AT", "A") == "DEL" and het_panel.vtype("AT", "GC") == "COMPLEX"


def _vcf(path, records, fmt="GT"):
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
        for r in records:
            fh.write("\t".join(map(str, r)) + "\n")


def test_het_panel_library_and_donor(tmp_path):
    tmp = str(tmp_path)
    rows = []
    for k, lib in enumerate(("LIBA", "LIBB")):
        d = os.path.join(tmp, lib)
        g = os.path.join(d, "genome_variant_detection", "results")
        for sub in ("haplotypecaller", "deepvariant", "clair3_rna", "longcallr"):
            os.makedirs(os.path.join(g, sub))
        _vcf(os.path.join(g, "haplotypecaller", "variants.vcf.gz"),
             [("chr7", 50, ".", "A", "G", 50, ".", ".", "GT", "0/1"),
              ("chr7", 90, ".", "AT", "A,ATT", 50, ".", ".", "GT", "1/2")])
        _vcf(os.path.join(g, "deepvariant", "variants.vcf.gz"),
             [("chr7", 50, ".", "A", "G", 50, "PASS", ".", "GT", "0/1"),
              ("chr7", 70, ".", "C", "T", 5, "RefCall", ".", "GT", "0/1")])
        _vcf(os.path.join(g, "clair3_rna", "variants.vcf.gz"),
             [("chr7", 50, ".", "A", "G", 50, "PASS", ".", "GT", "1/1")])
        _vcf(os.path.join(g, "longcallr", "variants.vcf.gz"),
             [("chr7", 50, ".", "A", "G", 50, "PASS", "RDS=select", "GT:PS", "0|1:40")])
        rows.append("U{}\t{}\t{}\tTRUE\tyoung".format(k, os.path.join(d, "x.bam"), lib))
    meta = os.path.join(tmp, "meta.txt")
    with open(meta, "w") as fh:
        fh.write("uid\tbam\tlibrary\treverse\tgroups\n" + "\n".join(rows) + "\n")
    dmap = os.path.join(tmp, "dmap.tsv")
    with open(dmap, "w") as fh:
        fh.write("uid\toriginal_uid\tdonor\nLIBA\tU0\tD\nLIBB\tU1\tD\n")
    out = os.path.join(tmp, "panel")
    for i in (0, 1):
        het_panel.main(["--library", "--bam_metadata", meta, "--bam_row", str(i), "--donor_map", dmap,
                        "--no_annotation", "--out", out])
    het_panel.main(["--donor", "D", "--bam_metadata", meta, "--donor_map", dmap, "--out", out])
    with gzip.open(os.path.join(out, "donors", "D.sites.tsv.gz"), "rt") as fh:
        got = {(r["pos"], r["ref"], r["alt"]): r for r in csv.DictReader(fh, delimiter="\t")}
    snv = got[("50", "A", "G")]
    assert snv["tier"] == "high" and snv["max_alt_callers"] == "4" and snv["max_het_callers"] == "3"
    assert snv["n_libraries"] == "2" and snv["n_libraries_alt"] == "2"
    assert snv["longcallr_ps_by_library"] == "00__LIBA=40,01__LIBB=40" and snv["longcallr_rds"] == "select"
    assert got[("90", "AT", "A")]["tier"] == "low" and got[("90", "A", "AT")]["vtype"] == "INS"
    assert got[("70", "C", "T")]["tier"] == "none_pass"         # RefCall is not a PASS call; the row is kept
    het_panel.main(["--concordance", "--bam_metadata", meta, "--donor_map", dmap, "--out", out])
    with open(os.path.join(out, "concordance", "library_best_match.tsv")) as fh:
        best = list(csv.DictReader(fh, delimiter="\t"))
    assert [b["status_jaccard"] for b in best] == ["PASS", "PASS"] and best[0]["best_jaccard"] == "1.0"
    assert [b["status_containment"] for b in best] == ["PASS", "PASS"] and best[0]["best_containment"] == "1.0"


def test_allelic_tests_end_to_end_tiers_region(world):
    """read_sites -> allelic_tests --genotypes_only -> allelic_tests --genotypes --tiers --region on the
    synthetic BAM: every output table is written and the selection is logged."""
    from altanalyze3.components.snp import allelic_tests
    w = world
    _run(w)
    reads = os.path.join(w["tmp"], "out")
    with open(w["isv_meta"]) as fh:
        rows = fh.read().splitlines()
    isv = os.path.join(w["tmp"], "isv_meta_groups.tsv")
    with open(isv, "w") as fh:
        fh.write(rows[0] + "\tgroups\n" + rows[1] + "\tyoung\n")
    stats = os.path.join(w["tmp"], "stats")
    base = ["--chrom", "chr7", "--reads", reads, "--donor_map", w["dmap"], "--isv_metadata", isv]
    assert allelic_tests.main(base + ["--genotypes_only", "--out", stats]) == 0
    g = os.path.join(stats, "chr7", "donor_site_genotypes.tsv.gz")
    assert os.path.exists(g)
    allelic_tests.main(base + ["--genotypes", g, "--tiers", "high", "--region", "1-300", "--out", stats])
    out = os.path.join(stats, "chr7_1_300")
    for f in ("isoform_tests.tsv.gz", "ase_site_donor.tsv.gz", "isoform_state_tests.tsv.gz"):
        assert os.path.exists(os.path.join(out, f)), f
    log = json.load(open(os.path.join(out, "logs", "allelic_tests.parameters.json")))
    assert log["tiers"] == "high" and log["region"] == "1-300"
    assert log["het_donor_sites_selected"] == sum(log["het_selected_by_tier"].values())
