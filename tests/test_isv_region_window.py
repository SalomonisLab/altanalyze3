"""Tests for the --region genomic window of isv_web_remote.

The window exists because the x-axis is otherwise the union of every plotted molecule, so an event inside
one exon of a long gene is drawn under a point wide. Two invariants matter:

1. parse_region accepts the coordinate spellings a user types and REFUSES anything it cannot read, rather
   than guessing a window.
2. segs_in_window / mol_in_window are a no-op whenever the window is the molecules' own extent, which is
   what the renderer uses without --region. That is what keeps every existing panel byte-identical.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from altanalyze3.components.visualization.isv_web.isv_web_remote import (  # noqa: E402
    mol_in_window,
    parse_region,
    segs_in_window,
)


def seg(a, b):
    return {"a": a, "b": b}


@pytest.mark.parametrize("text,expect", [
    ("1887200-1892300", (1887200, 1892300)),
    ("chr11:1887200-1892300", (1887200, 1892300)),
    ("chr11:1,887,200-1,892,300", (1887200, 1892300)),
    ("1887200..1892300", (1887200, 1892300)),
    ("  chr21:35,047,800-35,050,400  ", (35047800, 35050400)),
    ("chrX:1-2", (1, 2)),
])
def test_parse_region_accepts(text, expect):
    assert parse_region(text) == expect


@pytest.mark.parametrize("text", [None, "", 0])
def test_parse_region_passes_through_empty(text):
    assert parse_region(text) is None


@pytest.mark.parametrize("text", [
    "1887200",             # one bound
    "1887200-",            # missing end
    "a-b",                 # not integers
    "1892300-1887200",     # end before start
    "1887200-1887200",     # zero width
    "chr11:",              # nothing after the colon
])
def test_parse_region_refuses_rather_than_guesses(text):
    with pytest.raises(ValueError):
        parse_region(text)


def test_segs_in_window_keeps_overlapping_and_drops_outside():
    segs = [seg(100, 200), seg(250, 300), seg(400, 500), seg(900, 1000)]
    assert segs_in_window(segs, 240, 520) == [seg(250, 300), seg(400, 500)]


def test_segs_in_window_keeps_a_segment_that_straddles_an_edge():
    assert segs_in_window([seg(100, 300)], 200, 400) == [seg(100, 300)]
    assert segs_in_window([seg(300, 500)], 200, 400) == [seg(300, 500)]


def test_segs_in_window_touching_an_edge_counts_as_inside():
    assert segs_in_window([seg(100, 200)], 200, 400) == [seg(100, 200)]
    assert segs_in_window([seg(400, 500)], 200, 400) == [seg(400, 500)]


def test_segs_in_window_is_a_noop_over_the_segments_own_extent():
    """The renderer's lo/hi without --region is exactly this extent, so no existing panel can change."""
    segs = [seg(100, 200), seg(250, 300), seg(400, 500)]
    lo = min(s["a"] for s in segs)
    hi = max(s["b"] for s in segs)
    assert segs_in_window(segs, lo, hi) == segs


def test_mol_in_window_uses_exons_and_introns():
    merged = {"exons": [seg(100, 200)], "introns": [seg(200, 900)]}
    assert mol_in_window(merged, 300, 400) is True      # only the intron reaches the window
    assert mol_in_window(merged, 950, 1000) is False
    assert mol_in_window(merged, 150, 160) is True


def test_mol_in_window_handles_missing_keys():
    assert mol_in_window({}, 1, 10) is False
    assert mol_in_window({"exons": None, "introns": None}, 1, 10) is False
    assert mol_in_window({"exons": [seg(5, 6)]}, 1, 10) is True


def test_lsp1_window_separates_the_two_bins():
    """The real case: a 3'-fragment bin and a full-length bin over chr11:1,887,200-1,892,300.

    Both molecules overlap the window, so neither bin is dropped, and the full-length molecule keeps only
    its last two segments. The event is the donor at 1,889,609 against the E34.2 donor at 1,887,576.
    """
    novel = {"exons": [seg(1891773, 1892263)], "introns": [seg(1889000, 1889609)]}
    reference = {"exons": [seg(1853084, 1853197), seg(1880087, 1880224), seg(1887474, 1887576),
                           seg(1891773, 1892263)],
                 "introns": []}
    lo, hi = 1887200, 1892300
    assert mol_in_window(novel, lo, hi) is True
    assert mol_in_window(reference, lo, hi) is True
    assert segs_in_window(reference["exons"], lo, hi) == [seg(1887474, 1887576), seg(1891773, 1892263)]


def test_runx1_window_drops_reads_that_never_reach_the_event():
    """RUNX1 reference molecules reach E42 in a minority of reads; those that do not must be dropped, or
    they occupy a row and draw an edge-pinned tick at the window boundary."""
    lo, hi = 35047800, 35050400
    reaches = {"exons": [seg(35048958, 35049168), seg(34792022, 34792300)], "introns": []}
    misses = {"exons": [seg(34792022, 34792300), seg(34800000, 34800500)], "introns": []}
    assert mol_in_window(reaches, lo, hi) is True
    assert mol_in_window(misses, lo, hi) is False


from altanalyze3.components.visualization.isv_web.isv_web_remote import build_xz  # noqa: E402


def test_build_xz_fills_the_axis_when_the_window_ends_inside_an_intron():
    """The RUNX1 defect: an exon block beyond the window must not add its distance to the axis.

    Window chr21:35,048,300-35,049,600 holds one exon block (35,048,842-35,049,344). The next block of the
    gene sits at 35,086,000, far above the window. Before the clamp, that block's start became the end of a
    phantom 'intron gap' and the drawn content collapsed into the first few percent of the track.
    """
    lo, hi = 35048300, 35049600
    ivs = [(34792022, 34792300), (35048842, 35049344), (35086000, 35086400)]
    xz = build_xz(ivs, lo, hi, 1.0, 0.0, 1000.0)
    assert xz(lo) == pytest.approx(0.0, abs=1e-6)
    assert xz(hi) == pytest.approx(1000.0, abs=1e-6)
    # intron_scale 1.0 is a plain linear map, so every coordinate lands at its true fraction of the window
    for c in (35048842, 35048958, 35049168, 35049344):
        assert xz(c) == pytest.approx(1000.0 * (c - lo) / (hi - lo), abs=1e-6)
    # the event occupies its true share of the axis, not a sliver. The splice removes 35,048,959-35,049,167,
    # 209 bp, so donor and acceptor sit 210 bp apart on the axis.
    assert (xz(35049168) - xz(35048958)) / 1000.0 == pytest.approx(210 / 1300, abs=1e-6)


def test_build_xz_window_ending_inside_an_intron_with_no_exon_beyond():
    """The LSP1 window: nothing lies above it, so only a short trailing intron is added."""
    lo, hi = 1887300, 1892400
    ivs = [(1887474, 1887576), (1891773, 1892267)]
    xz = build_xz(ivs, lo, hi, 1.0, 0.0, 1000.0)
    assert xz(lo) == pytest.approx(0.0, abs=1e-6)
    assert xz(hi) == pytest.approx(1000.0, abs=1e-6)
    assert xz(1889609) == pytest.approx(1000.0 * (1889609 - lo) / (hi - lo), abs=1e-6)


def test_build_xz_unchanged_when_hi_is_the_last_exon_end():
    """Every panel rendered without --region has hi at a drawn exon's end, where the clamp is a no-op."""
    lo, hi = 100, 500
    ivs = [(100, 200), (300, 500)]
    a = build_xz(ivs, lo, hi, 0.1, 0.0, 1000.0)
    b = build_xz(ivs, lo, hi, 0.1, 0.0, 1000.0)
    for c in (100, 150, 200, 250, 300, 400, 500):
        assert a(c) == b(c)
    # exons keep true scale, the intron is shrunk 10x: 100 + 10 + 200 = 310 units
    assert a(200) == pytest.approx(1000.0 * 100 / 310, abs=1e-6)
    assert a(300) == pytest.approx(1000.0 * 110 / 310, abs=1e-6)


def test_retained_intron_keeps_true_scale():
    """A retained intron is sequence the read carries; intron_scale must not shrink it.

    Real RUNX1 novel-bin read 45456890 (Y863, MPP-1): I:35045994-35048841 (2,848 bp retained),
    E:35048842-35048958 (117 bp), E:35049168-35049260 (93 bp). With only the exons given true scale the
    read is drawn far shorter than its 3,058 bp of sequence; including the retained intron fixes it.

    The invariant: on one axis, a retained intron of length L must draw the same width as an exon of
    length L. A reference exon of the same 2,848 bp sits elsewhere on the axis as the yardstick.
    """
    lo, hi = 35040000, 35060000
    yardstick = (35052000, 35054848)                  # a 2,848 bp exon elsewhere on the same axis
    retained = (35045994, 35048841)                   # the read's 2,848 bp retained intron
    exons = [(35048842, 35048958), (35049168, 35049260), yardstick]
    W = 1000.0
    broken = build_xz(sorted(exons), lo, hi, 0.1, 0.0, W)
    fixed = build_xz(sorted(exons + [retained]), lo, hi, 0.1, 0.0, W)
    yard_b = broken(yardstick[1]) - broken(yardstick[0])
    ret_b = broken(retained[1]) - broken(retained[0])
    yard_f = fixed(yardstick[1]) - fixed(yardstick[0])
    ret_f = fixed(retained[1]) - fixed(retained[0])
    assert ret_b == pytest.approx(0.1 * yard_b, rel=1e-3)   # before: shrunk 10x against an equal exon
    assert ret_f == pytest.approx(yard_f, rel=1e-3)         # after: identical width, as it must be
    # the spliced gap between the two exons (35,048,959-35,049,167) stays compressed after the fix
    assert (fixed(35049168) - fixed(35048958)) == pytest.approx(0.1 * 210 / 2848 * yard_f, rel=1e-2)


def test_spliced_gap_between_segments_is_still_compressed():
    """Only segments the read covers get true scale. A gap with no segment stays at intron_scale."""
    lo, hi = 1000, 11000
    ivs = [(1000, 1500), (10500, 11000)]          # two covered blocks, 9,000 bp spliced out between
    xz = build_xz(ivs, lo, hi, 0.1, 0.0, 1000.0)
    total = 500 + 9000 * 0.1 + 500                 # 1,900 units
    assert (xz(1500) - xz(1000)) / 1000.0 == pytest.approx(500 / total, abs=1e-6)
    assert (xz(10500) - xz(1500)) / 1000.0 == pytest.approx(900 / total, abs=1e-6)


def molecule_backbone(merged, lo, hi):
    """The span the renderer draws the connector line across: every segment of the molecule."""
    span = segs_in_window(merged["exons"], lo, hi) + segs_in_window(merged["introns"], lo, hi)
    return (min(s["a"] for s in span), max(s["b"] for s in span)) if span else None


def test_backbone_spans_a_retained_intron_far_from_the_exons():
    """LSP1 gold read: I:1889516-1889609 then E:1891773-1892262. The connector must reach the I block.

    Taking min/max over the exons alone started the line at 1,891,773 and left the retained-intron block
    floating 2,164 bp to its left with nothing joining it to the C-terminal exon.
    """
    merged = {"exons": [seg(1891773, 1891974), seg(1891974, 1892011), seg(1892011, 1892262)],
              "introns": [seg(1889516, 1889609)]}
    lo, hi = 1887300, 1892400
    assert molecule_backbone(merged, lo, hi) == (1889516, 1892262)
    exons_only = min(s["a"] for s in merged["exons"])
    assert exons_only == 1891773                      # what the old code used
    assert molecule_backbone(merged, lo, hi)[0] < exons_only


def test_backbone_unchanged_when_every_segment_is_an_exon():
    """The reference bin starts on an exon, so its backbone is identical before and after the change."""
    merged = {"exons": [seg(1853099, 1853197), seg(1880087, 1880224), seg(1892011, 1892262)],
              "introns": []}
    lo, hi = 1850000, 1893000
    assert molecule_backbone(merged, lo, hi) == (1853099, 1892262)


def test_backbone_covers_a_leading_retained_intron_runx1():
    """RUNX1 gold read 45456890: I:35045994-35048841 then two exons. The I block is contiguous with the
    first exon, so the old code happened to look connected; the backbone must still start at the I block."""
    merged = {"exons": [seg(35048842, 35048958), seg(35049168, 35049260)],
              "introns": [seg(35045994, 35048841)]}
    lo, hi = 35040000, 35060000
    assert molecule_backbone(merged, lo, hi) == (35045994, 35049260)


from altanalyze3.components.visualization.isv_web.isv_web_remote import (  # noqa: E402
    mol_extent, row_sort_key)


def mol(ci, segs_e, segs_i=(), iso="m", count=1):
    merged = {"exons": [seg(a, b) for a, b in segs_e], "introns": [seg(a, b) for a, b in segs_i]}
    return {"cluster_index": ci, "merged": merged, "isoform_id": iso, "count": count}


def test_mol_extent_covers_exons_and_retained_introns():
    m = mol(0, [(200, 300)], [(100, 150)])
    assert mol_extent(m["merged"]) == (100, 300, 200)
    assert mol_extent({}) == (0, 0, 0)


def test_rows_order_longest_first_within_a_colour_group():
    short = mol(0, [(1000, 1100)], iso="a")
    mid = mol(0, [(500, 1100)], iso="b")
    long_ = mol(0, [(100, 1100)], iso="c")
    out = sorted([short, long_, mid], key=row_sort_key)
    assert [m["isoform_id"] for m in out] == ["c", "b", "a"]


def test_colour_group_still_dominates_the_order():
    """Gold and blue must not interleave: the colour key stays the primary sort term."""
    g_short = mol(0, [(1000, 1100)], iso="g")
    b_long = mol(1, [(100, 1100)], iso="b")
    out = sorted([b_long, g_short], key=row_sort_key)
    assert [m["cluster_index"] for m in out] == [0, 1]


def test_order_is_deterministic_for_equal_length_molecules():
    """Count is constant at 1 for long reads, so ties must break on position then id, not dict order."""
    a = mol(0, [(100, 200)], iso="zzz")
    b = mol(0, [(100, 200)], iso="aaa")
    c = mol(0, [(300, 400)], iso="mmm")
    assert [m["isoform_id"] for m in sorted([a, b, c], key=row_sort_key)] == ["aaa", "zzz", "mmm"]
    assert [m["isoform_id"] for m in sorted([c, b, a], key=row_sort_key)] == ["aaa", "zzz", "mmm"]


def test_retained_intron_counts_toward_the_length_used_for_ordering():
    """A RUNX1 gold read is mostly retained intron; it must rank by its true extent, not its exons."""
    mostly_intron = mol(0, [(35048842, 35048958)], [(35045994, 35048841)], iso="runx1_gold")
    exon_only = mol(0, [(35048842, 35049260)], iso="exon_only")
    assert mol_extent(mostly_intron["merged"])[2] == 35048958 - 35045994
    assert [m["isoform_id"] for m in sorted([exon_only, mostly_intron], key=row_sort_key)] \
        == ["runx1_gold", "exon_only"]
