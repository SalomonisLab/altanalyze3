"""Read rows must be drawn from the genomic exons the source GFF recorded, not re-derived.

build_isoform_segments spans min..max over a base that carries two coordinates, so the tokens
'E4.1|E4.1_66859679|E4.1_66859746' were drawn as one 68 bp block 66,859,679-66,859,746. Those two
numbers are where one GFF exon ends and the next begins: the real block is the 43 bp exon
chr11:66,859,637-66,859,679. gff_process.py:769 already records the true blocks in
gff-output/structure_coords.sqlite, keyed by the same token string.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from altanalyze3.components.visualization.isv_web.data_api import _read_segments  # noqa: E402

LRFN4 = "ENSG00000173621"
TOKENS = ["E1.1_66856657", "E1.1_66856776", "E1.1", "I1.1", "E2.1", "E3.1_66857513",
          "E3.1_66857544", "E3.3_66858374", "E3.3", "E3.4", "E4.1", "E4.1_66859679",
          "E4.1_66859746"]
GFF_EXONS = [(66856538, 66856657), (66856776, 66857296), (66857513, 66857544),
             (66858374, 66859093), (66859637, 66859679), (66859746, 66860473)]


class _Ctx:
    """Stands in for RunContext: a coord store that holds one structure, plus an exon lookup."""
    def __init__(self, with_store=True):
        self.with_store = with_store
        self.exon_lookup = {
            (LRFN4, "E1.1"): {"chrom": "chr11", "strand": "+", "start": 66856647, "end": 66856864, "type": "E"},
            (LRFN4, "I1.1"): {"chrom": "chr11", "strand": "+", "start": 66856865, "end": 66857117, "type": "I"},
            (LRFN4, "E2.1"): {"chrom": "chr11", "strand": "+", "start": 66857118, "end": 66857296, "type": "E"},
            (LRFN4, "E3.1"): {"chrom": "chr11", "strand": "+", "start": 66857502, "end": 66857735, "type": "E"},
            (LRFN4, "E3.3"): {"chrom": "chr11", "strand": "+", "start": 66857790, "end": 66858876, "type": "E"},
            (LRFN4, "E3.4"): {"chrom": "chr11", "strand": "+", "start": 66858876, "end": 66859093, "type": "E"},
            (LRFN4, "E4.1"): {"chrom": "chr11", "strand": "+", "start": 66859637, "end": 66860247, "type": "E"},
        }

    def structure_exons(self, gene, structure):
        if not self.with_store or gene != LRFN4 or structure != "|".join(TOKENS):
            return None
        return list(GFF_EXONS)


def test_read_segments_use_the_recorded_gff_blocks():
    segs = _read_segments(_Ctx(), LRFN4, TOKENS)
    assert [(s["start"], s["end"]) for s in segs] == GFF_EXONS


def test_the_43bp_exon_is_drawn_and_the_intron_is_not():
    segs = _read_segments(_Ctx(), LRFN4, TOKENS)
    spans = {(s["start"], s["end"]) for s in segs}
    assert (66859637, 66859679) in spans          # the 43 bp exon
    assert (66859679, 66859746) not in spans      # the intron the old path drew


def test_fallback_when_the_coord_store_has_no_row():
    segs = _read_segments(_Ctx(with_store=False), LRFN4, TOKENS)
    assert segs and [(s["start"], s["end"]) for s in segs] != GFF_EXONS
