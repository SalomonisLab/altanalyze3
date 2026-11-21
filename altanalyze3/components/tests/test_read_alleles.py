"""read_alleles: one allele class per read at a known SNV, keyed by the ISV molecule id."""
import csv
import json

import pysam
import pytest

from altanalyze3.components.bam import read_alleles as ra

SEQ = ("ACGT" * 50)            # 200 bp contig; 1-based position 101 holds 'A'
POS = 101                      # 1-based site
REF, ALT = SEQ[POS - 1], "G"


def _read(name, start0, cigar, seq, mapq=60, flag=0, cb="AAAACCCCGGGGTTTT"):
    r = pysam.AlignedSegment()
    r.query_name = name
    r.reference_id = 0
    r.reference_start = start0
    r.cigarstring = cigar
    r.query_sequence = seq
    r.mapping_quality = mapq
    r.flag = flag
    r.set_tag("CB", cb)
    return r


def _with_base(start0, length, base):
    s = list(SEQ[start0:start0 + length])
    s[POS - 1 - start0] = base
    return "".join(s)


@pytest.fixture()
def fixture_paths(tmp_path):
    fa = tmp_path / "ref.fa"
    fa.write_text(">chrT\n" + SEQ + "\n")
    pysam.faidx(str(fa))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrT", "LN": len(SEQ)}]}
    reads = [
        _read("molecule/1", 80, "50M", _with_base(80, 50, REF)),
        _read("molecule/2", 80, "50M", _with_base(80, 50, ALT)),
        _read("molecule/3", 80, "50M", _with_base(80, 50, "T")),
        _read("molecule/4", 70, "20M30N20M", SEQ[70:90] + SEQ[120:140]),   # site 100 (0-based) in the gap
        _read("molecule/5", 80, "20M5D25M", SEQ[80:100] + SEQ[105:130]),   # site in the deletion
        _read("molecule/6", 20, "40M", SEQ[20:60]),                        # in the region, stops before the site
        _read("molecule/7", 80, "50M", _with_base(80, 50, ALT), mapq=5),   # filtered by MAPQ
        _read("molecule/8", 80, "50M", _with_base(80, 50, ALT), flag=256),  # secondary, kept and flagged
    ]
    reads.sort(key=lambda r: r.reference_start)
    bam = tmp_path / "s1.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.index(str(bam))
    samples = tmp_path / "samples.tsv"
    samples.write_text("uid\tbam\nS1\t{}\n".format(bam))
    return fa, samples, tmp_path


def _rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_every_read_gets_one_class_and_the_isv_molecule_id(fixture_paths):
    fa, samples, tmp = fixture_paths
    out = tmp / "out"
    assert ra.main(["--variant", "rsT:chrT:{}:{}:{}".format(POS, REF, ALT), "--samples", str(samples),
                    "--fasta", str(fa), "--region", "chrT:1-200", "--out", str(out)]) == 0
    rows = {r["molecule_id"]: r for r in _rows(out / "read_alleles.tsv")}
    assert {k: v["allele_class"] for k, v in rows.items()} == {
        "1.S1": "REF", "2.S1": "ALT", "3.S1": "OTHER", "4.S1": "SPLICE", "5.S1": "DEL",
        "6.S1": "NOT_SPANNED", "8.S1": "ALT"}
    assert rows["8.S1"]["secondary"] == "1"
    assert rows["3.S1"]["base"] == "T"
    assert rows["1.S1"]["cb_rc"] == "AAAACCCCGGGGTTTT"[::-1].translate(str.maketrans("ACGT", "TGCA"))
    (lib,) = _rows(out / "library_allele_counts.tsv")
    assert (lib["reads"], lib["REF"], lib["ALT"], lib["OTHER"], lib["SPLICE"], lib["DEL"], lib["NOT_SPANNED"]) == \
        ("7", "1", "2", "1", "1", "1", "1")
    assert lib["secondary_or_supplementary"] == "1"
    params = json.loads((out / "logs" / "read_alleles.parameters.json").read_text())
    assert params["outputs"]["read_rows"] == 7 and params["min_mapq"] == 20


def test_without_region_only_reads_over_the_site_are_listed(fixture_paths):
    fa, samples, tmp = fixture_paths
    out = tmp / "site_only"
    ra.main(["--variant", "rsT:chrT:{}:{}:{}".format(POS, REF, ALT), "--samples", str(samples),
             "--fasta", str(fa), "--out", str(out)])
    classes = sorted(r["allele_class"] for r in _rows(out / "read_alleles.tsv"))
    assert classes == ["ALT", "ALT", "DEL", "OTHER", "REF", "SPLICE"]


def test_ref_must_match_the_fasta(fixture_paths):
    fa, samples, tmp = fixture_paths
    wrong = "C" if REF != "C" else "T"
    with pytest.raises(ValueError, match="differs from REF"):
        ra.main(["--variant", "rsT:chrT:{}:{}:G".format(POS, wrong), "--samples", str(samples),
                 "--fasta", str(fa), "--out", str(tmp / "bad")])


def test_isv_metadata_join_refuses_a_library_mismatch(tmp_path):
    (tmp_path / "bam_meta.tsv").write_text("uid\tbam\tlibrary\tgroups\nA\t/x/a.bam\tLIB_A\tyoung\n")
    (tmp_path / "isv_meta.tsv").write_text("uid\toriginal_library\toriginal_metadata_ordinal\nLIB_A\tLIB_B\t0\n")
    args = ra.build_arg_parser().parse_args(["--fasta", "x", "--out", "y", "--isv_metadata", str(tmp_path / "isv_meta.tsv"),
                                             "--bam_metadata", str(tmp_path / "bam_meta.tsv")])
    with pytest.raises(ValueError, match="library mismatch"):
        ra.load_samples(args)


@pytest.mark.parametrize("n_ref,n_alt,expected", [(0, 29, "1/1"), (14, 19, "0/1"), (11, 1, "0/0"),
                                                  (2, 825, "1/1"), (0, 0, "./.")])
def test_genotype_call(n_ref, n_alt, expected):
    assert ra.genotype_call(n_ref, n_alt, 0.01)[0] == expected


def test_isv_pairing_check_stops_on_a_library_read_from_the_wrong_bam(fixture_paths):
    fa, samples, tmp = fixture_paths
    good = tmp / "isv_good.tsv"
    good.write_text("molecule_id\tfinal_isoform_id\tcell_state\tsample\tcount\n1.S1\tENST1\tMono\tS1\t1.0\n"
                    "99.S1\tENST1\tMono\tS1\t1.0\n")
    out = tmp / "paired"
    ra.main(["--variant", "rsT:chrT:{}:{}:{}".format(POS, REF, ALT), "--samples", str(samples), "--fasta", str(fa),
             "--region", "chrT:1-200", "--isv_molecules", str(good), "--out", str(out)])
    (row,) = _rows(out / "isv_molecule_join.tsv")
    assert (row["isv_molecules"], row["found_in_bam"], row["missing"]) == ("2", "1", "1")
    swapped = tmp / "isv_swapped.tsv"
    swapped.write_text("molecule_id\tfinal_isoform_id\tcell_state\tsample\tcount\n500.S1\tENST1\tMono\tS1\t1.0\n")
    with pytest.raises(SystemExit, match="BAM pairing failure"):
        ra.main(["--variant", "rsT:chrT:{}:{}:{}".format(POS, REF, ALT), "--samples", str(samples), "--fasta", str(fa),
                 "--region", "chrT:1-200", "--isv_molecules", str(swapped), "--out", str(tmp / "bad_pair")])
