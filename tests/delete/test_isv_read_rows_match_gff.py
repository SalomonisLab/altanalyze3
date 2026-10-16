"""Read rows must be drawn at the read's own genomic coordinates.

Rebuilding a read from its exon-ID tokens does not reproduce its GFF block boundaries: a base
carrying two coordinates, 'E4.1|E4.1_66859679|E4.1_66859746', marks where one block ends and the
next begins, and build_isoform_segments spanned min..max. On LRFN4 the drawn blocks disagreed with
the GFF for 5 of 5 reads of 35292108.AML-7_1_CITE_GEX, losing the 43 bp exon
chr11:66,859,637-66,859,679 and drawing the intron in its place.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from altanalyze3.components.visualization.isv_web.data_api import _read_segments  # noqa: E402

LRFN4 = "ENSG00000173621"
READ = "35292108.AML-7_1_CITE_GEX"
TOKENS = ["E1.1_66856657", "E1.1_66856776", "E1.1", "I1.1", "E2.1", "E3.1_66857513",
          "E3.1_66857544", "E3.3_66858374", "E3.3", "E3.4", "E4.1", "E4.1_66859679",
          "E4.1_66859746"]
GFF = [(66856538, 66856657), (66856776, 66857296), (66857513, 66857544),
       (66858374, 66859093), (66859637, 66859679), (66859746, 66860473)]


class _Ctx:
    def __init__(self, store=True):
        self.store = store
        self.exon_lookup = {
            (LRFN4, "E1.1"): {"chrom": "chr11", "strand": "+", "start": 66856647, "end": 66856864, "type": "E"},
            (LRFN4, "I1.1"): {"chrom": "chr11", "strand": "+", "start": 66856865, "end": 66857117, "type": "I"},
            (LRFN4, "E2.1"): {"chrom": "chr11", "strand": "+", "start": 66857118, "end": 66857296, "type": "E"},
            (LRFN4, "E3.1"): {"chrom": "chr11", "strand": "+", "start": 66857502, "end": 66857735, "type": "E"},
            (LRFN4, "E3.3"): {"chrom": "chr11", "strand": "+", "start": 66857790, "end": 66858876, "type": "E"},
            (LRFN4, "E3.4"): {"chrom": "chr11", "strand": "+", "start": 66858876, "end": 66859093, "type": "E"},
            (LRFN4, "E4.1"): {"chrom": "chr11", "strand": "+", "start": 66859637, "end": 66860247, "type": "E"},
        }

    def read_exons(self, gene, read_id):
        if not self.store or gene != LRFN4 or str(read_id).split(".", 1)[0] != "35292108":
            return None
        return list(GFF)


def test_read_row_equals_the_gff_blocks():
    segs = _read_segments(_Ctx(), LRFN4, READ, TOKENS)
    assert [(s["start"], s["end"]) for s in segs] == GFF


def test_the_43bp_exon_is_drawn_and_the_intron_is_not():
    spans = {(s["start"], s["end"]) for s in _read_segments(_Ctx(), LRFN4, READ, TOKENS)}
    assert (66859637, 66859679) in spans
    assert (66859679, 66859746) not in spans


def test_token_fallback_is_used_when_the_read_is_absent():
    segs = _read_segments(_Ctx(store=False), LRFN4, READ, TOKENS)
    assert segs and [(s["start"], s["end"]) for s in segs] != GFF
