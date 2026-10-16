"""build_isoform_segments must reproduce the source GFF's genomic blocks.

The token string is a boundary list. The old implementation grouped consecutive tokens sharing a
base and spanned min..max of their coordinates, so 'E4.1|E4.1_66859679|E4.1_66859746' became one
68 bp block over the GFF's intron and the real 43 bp exon 66,859,637-66,859,679 disappeared.

Ground truth here is the per-sample GFF record for read 35292108 of AML-7_1_CITE_GEX, LRFN4, read
off gff-output/combined.gff.gz and _isv_web_cache/final_structures.db, which agree:
    66856538-66856657, 66856776-66857296, 66857513-66857544,
    66858374-66859093, 66859637-66859679, 66859746-66860473
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from altanalyze3.components.visualization.isoform_structure_view import (  # noqa: E402
    build_isoform_segments,
)

G = "ENSG00000173621"
EXONS = {"E1.1": (66856647, 66856864, "E"), "I1.1": (66856865, 66857117, "I"),
         "E2.1": (66857118, 66857296, "E"), "I2.1": (66857297, 66857501, "I"),
         "E3.1": (66857502, 66857735, "E"), "E3.2": (66857735, 66857790, "E"),
         "E3.3": (66857790, 66858876, "E"), "E3.4": (66858876, 66859093, "E"),
         "I3.1": (66859094, 66859636, "I"), "E4.1": (66859637, 66860247, "E"),
         "E4.2": (66860247, 66860475, "E")}
LOOKUP = {(G, k): {"chrom": "chr11", "strand": "+", "start": v[0], "end": v[1], "type": v[2]}
          for k, v in EXONS.items()}
RAW = ("U1.1_66856538|E1.1_66856657|E1.1_66856776|E1.1|I1.1|E2.1|E3.1_66857513|E3.1_66857544|"
       "E3.3_66858374|E3.3|E3.4|E4.1|E4.1_66859679|E4.1_66859746|E4.2_66860473").split("|")
GFF = [(66856538, 66856657), (66856776, 66857296), (66857513, 66857544),
       (66858374, 66859093), (66859637, 66859679), (66859746, 66860473)]


def _blocks(segments):
    """Merge adjacent segments back into the genomic blocks the renderer draws."""
    out = []
    for s in sorted(segments, key=lambda s: s["start"]):
        if out and s["start"] <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], s["end"]))
        else:
            out.append((s["start"], s["end"]))
    return out


def test_blocks_equal_the_source_gff():
    assert _blocks(build_isoform_segments(RAW, LOOKUP, G)) == GFF


def test_the_43bp_exon_is_present_and_the_intron_is_not_drawn_as_one():
    blocks = _blocks(build_isoform_segments(RAW, LOOKUP, G))
    assert (66859637, 66859679) in blocks
    assert (66859679, 66859746) not in blocks


def test_retained_intron_stays_typed_I():
    segs = build_isoform_segments(RAW, LOOKUP, G)
    introns = [(s["start"], s["end"]) for s in segs if s["type"] == "I"]
    assert (66856865, 66857117) in introns        # I1.1 inside the 521 bp block
    assert all(s <= 66857117 and e >= 66856865 for s, e in introns)


def test_five_prime_extension_outside_any_exon_is_kept():
    segs = build_isoform_segments(RAW, LOOKUP, G)
    assert any(s["start"] == 66856538 for s in segs)


def test_empty_and_unknown_tokens_are_safe():
    assert build_isoform_segments([], LOOKUP, G) == []
    assert build_isoform_segments(["ZZ9.9"], LOOKUP, G) == []


# ---------------------------------------------------------------------------------------------
# All five reads of final isoform 35292108.AML-7_1_CITE_GEX, LRFN4, root index.
# Token strings come from each library's gff-output/gene_indexes_v2/transcripts.db, read through
# isv.load_transcript_associations_from_sqlite. The blocks come from that library's own sample GFF
# with gene_id required to be ENSG00000173621, and for 35292108 they also equal the aligned blocks
# of molecule/35292108 in AML7.mapped.bam.
#
# 39053794 is the case that abutment alone cannot parse: 'E3.1_66857513|E3.3_66858795' is one GFF
# block 66,857,513-66,858,795, bridged by E3.2.
READS = {
    "35292108.AML-7_1_CITE_GEX": (
        "U1.1_66856538|E1.1_66856657|E1.1_66856776|E1.1|I1.1|E2.1|E3.1_66857513|E3.1_66857544|"
        "E3.3_66858374|E3.3|E3.4|E4.1|E4.1_66859679|E4.1_66859746|E4.2_66860473",
        [(66856538, 66856657), (66856776, 66857296), (66857513, 66857544),
         (66858374, 66859093), (66859637, 66859679), (66859746, 66860473)]),
    "35800172.AML-7_1_CITE_GEX": (
        "E2.1_66857072|E2.1|E3.1_66857513|E3.3_66857821|E3.3_66858431|E3.3|E3.4|E4.1|E4.2_66860474",
        [(66857072, 66857296), (66857513, 66857821), (66858431, 66859093), (66859637, 66860474)]),
    "20304321.AML-13_CITE_GEX": (
        "E2.1_66857076|E2.1|E3.1_66857513|E3.3_66858333|E3.4_66858967|E3.4|E4.1|E4.2_66860472",
        [(66857076, 66857296), (66857513, 66858333), (66858967, 66859093), (66859637, 66860472)]),
    "32769458.Y3936": (
        "E2.1_66857092|E2.1|E3.1_66857513|E3.1_66857691|E3.3_66858822|E3.3|E3.4|E4.1|E4.2|"
        "U5.1_66860477",
        [(66857092, 66857296), (66857513, 66857691), (66858822, 66859093), (66859637, 66860477)]),
    "39053794.Y863": (
        "I1.1_66857064|I1.1|E2.1|E3.1_66857513|E3.3_66858795|E3.4_66858929|E3.4|E4.1|E4.2_66860472",
        [(66857064, 66857296), (66857513, 66858795), (66858929, 66859093), (66859637, 66860472)]),
}


def test_every_read_of_the_isoform_reproduces_its_own_gff_record():
    wrong = {}
    for read, (raw, gff) in READS.items():
        got = _blocks(build_isoform_segments(raw.split("|"), LOOKUP, G))
        if got != gff:
            wrong[read] = (got, gff)
    assert not wrong, wrong


def test_exon_bridged_block_is_not_split():
    """E3.2 bridges E3.1 and E3.3, so 'E3.1_66857513|E3.3_66858795' is ONE block."""
    raw, gff = READS["39053794.Y863"]
    blocks = _blocks(build_isoform_segments(raw.split("|"), LOOKUP, G))
    assert (66857513, 66858795) in blocks
    assert blocks == gff


def test_intron_between_exon_blocks_still_splits():
    """I2.1 is not exonic, so E2.1 -> E3.1 stays two blocks in every read that spans both."""
    for read, (raw, gff) in READS.items():
        blocks = _blocks(build_isoform_segments(raw.split("|"), LOOKUP, G))
        assert not any(s <= 66857296 and e >= 66857513 for s, e in blocks), read


def test_terminal_unknown_token_closes_the_last_block():
    """'U5.1_66860477' is not in the gene model; it is the 3' end of the read, not a new block."""
    raw, gff = READS["32769458.Y3936"]
    blocks = _blocks(build_isoform_segments(raw.split("|"), LOOKUP, G))
    assert blocks[-1] == (66859637, 66860477)


# ---------------------------------------------------------------------------------------------
# Minus strand. Exon ids are numbered in TRANSCRIPTION order, so E1.1 sits at the highest genomic
# coordinate and the token string descends. The LRFN4 read 35292108 is mirrored about 133,717,013
# (= 66,856,538 + 66,860,475) onto a minus-strand gene: the same structure, read the other way.
MIRROR = 66856538 + 66860475
G_MINUS = "ENSGMINUS000001"


def _flip(a, b):
    return (MIRROR - b, MIRROR - a)


EXONS_MINUS = {k: _flip(v[0], v[1]) + (v[2],) for k, v in EXONS.items()}
LOOKUP_MINUS = {(G_MINUS, k): {"chrom": "chr11", "strand": "-", "start": v[0], "end": v[1],
                               "type": v[2]} for k, v in EXONS_MINUS.items()}
RAW_MINUS = [t if "_" not in t else "%s_%d" % (t.split("_")[0], MIRROR - int(t.split("_")[1]))
             for t in RAW]
GFF_MINUS = sorted(_flip(a, b) for a, b in GFF)


def test_minus_strand_tokens_are_read_in_genomic_order():
    assert _blocks(build_isoform_segments(RAW_MINUS, LOOKUP_MINUS, G_MINUS)) == GFF_MINUS


def test_minus_strand_keeps_the_43bp_exon_and_the_retained_intron():
    segs = build_isoform_segments(RAW_MINUS, LOOKUP_MINUS, G_MINUS)
    blocks = _blocks(segs)
    assert _flip(66859637, 66859679) in blocks
    assert _flip(66856865, 66857117) in [(s["start"], s["end"]) for s in segs if s["type"] == "I"]
