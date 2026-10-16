"""Tests for the isoform ordering of isv_web_remote: LONGEST ISOFORM FIRST, by default.

Isoform order used to be the cluster index, which the read index assigns by abundance. The bottom
structure track and the colour groups above it therefore ran most-read-first, in no length order:
on FLNB chr3:57,994,127-58,172,251 the 163,830 bp ENST00000295956 was followed by the 12,979 bp
ENST00000419752 and then the 26,420 bp ENST00000481470.

Three invariants:

1. The length of an isoform is the length of the structure the figure DRAWS: the registered
   combined.gff.gz exons when the isoform is registered, the representative read otherwise.
2. One rank orders both the read rows and the structure rows, so a panel and the track below it
   list the isoforms in the same order.
3. An explicit --isoforms order (color_index) is the user's own order and is never overridden.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from altanalyze3.components.visualization.isv_web.isv_web_remote import (  # noqa: E402
    color_key,
    group_key,
    iso_order_map,
    iso_struct_extent,
    merge_segments,
    render,
    row_sort_key,
)


def mol(mid, iso, segs, cluster, registered=None, count=1.0):
    return {"isoform_id": mid, "final_isoform_id": iso, "cluster_index": cluster, "count": count,
            "known": True, "registered_exons": registered,
            "exon_segments": [{"start": a, "end": b, "type": "E", "exon_id": "E%d.1" % i}
                              for i, (a, b) in enumerate(segs, 1)]}


def test_iso_struct_extent_prefers_the_registered_structure():
    m = mol("r1", "ENST1", [(58159623, 58170783)], 0, registered=[(58008422, 58008500), (58172000, 58172251)])
    m["merged"] = merge_segments(m["exon_segments"])
    lo, hi, bp = iso_struct_extent(m)
    assert (lo, hi) == (58008422, 58172251)          # the registered span, not the read's 11,161 bp
    assert bp == 78 + 251


def test_iso_struct_extent_falls_back_to_the_read():
    m = mol("r1", "NOVEL1", [(100, 400), (900, 1000)], 0)
    m["merged"] = merge_segments(m["exon_segments"])
    assert iso_struct_extent(m) == (100, 1000, 400)


def test_iso_order_map_is_longest_isoform_first():
    mols = [
        mol("r1", "SHORT", [(58159623, 58172251)], 0),      # 12,628 bp, most reads -> cluster 0
        mol("r2", "LONG", [(58008422, 58172251)], 1),       # 163,829 bp
        mol("r3", "MID", [(58145831, 58172251)], 2),        # 26,420 bp
    ]
    for m in mols:
        m["merged"] = merge_segments(m["exon_segments"])
    rank = iso_order_map(mols)
    assert sorted(rank, key=rank.get) == ["LONG", "MID", "SHORT"]


def test_iso_order_map_ties_break_on_exonic_bp_then_start_then_id():
    a = mol("r1", "A", [(100, 200), (900, 1000)], 0)        # span 900, 200 exonic bp
    b = mol("r2", "B", [(100, 1000)], 1)                    # span 900, 900 exonic bp
    for m in (a, b):
        m["merged"] = merge_segments(m["exon_segments"])
    rank = iso_order_map([a, b])
    assert sorted(rank, key=rank.get) == ["B", "A"]


def test_group_key_uses_the_length_rank_only_when_stamped():
    m = {"cluster_index": 7}
    assert group_key(m) == color_key(m) == 7
    m["_iso_rank"] = 0
    assert group_key(m) == 0


def test_row_sort_key_groups_by_isoform_then_longest_molecule_first():
    short = {"_iso_rank": 0, "_span": (100, 300, 200), "isoform_id": "r2"}
    long_ = {"_iso_rank": 0, "_span": (100, 900, 800), "isoform_id": "r1"}
    other = {"_iso_rank": 1, "_span": (100, 999, 899), "isoform_id": "r0"}
    assert sorted([short, other, long_], key=row_sort_key) == [long_, short, other]


def _res():
    """Three isoforms whose read abundance is the REVERSE of their length."""
    gm = {"chrom": "chr3", "strand": "+", "extent": [58008422, 58172251],
          "blocks": [{"start": 58008422, "end": 58008600}, {"start": 58172000, "end": 58172251}]}
    mols = []
    for i in range(6):                                   # most reads, shortest isoform
        mols.append(mol("s%d" % i, "SHORT", [(58159623, 58172251)], 0))
    for i in range(3):
        mols.append(mol("m%d" % i, "MID", [(58145831, 58172251)], 1))
    for i in range(2):                                   # fewest reads, longest isoform
        mols.append(mol("l%d" % i, "LONG", [(58008422, 58010000), (58170000, 58172251)], 2))
    return {"gene": "G", "symbol": "FLNB", "gene_model": gm,
            "panels": [{"condition": "HSC-MPP", "n_molecules": len(mols), "molecules": mols}]}


def _drawn_order(panel):
    out = []
    for m in panel["mol"]:
        k = m["final_isoform_id"]
        if k not in out:
            out.append(k)
    return out


def test_render_orders_isoform_groups_longest_first(tmp_path):
    res = _res()
    out = str(tmp_path / "order.pdf")
    assert render(res, out, isoform_track=True) is not False
    assert os.path.getsize(out) > 0
    assert _drawn_order(res["panels"][0]) == ["LONG", "MID", "SHORT"]


def test_render_keeps_an_explicit_isoforms_order(tmp_path):
    """--isoforms SHORT,MID,LONG: main() stamps color_index and passes iso_order="given"."""
    res = _res()
    order = {"SHORT": 0, "MID": 1, "LONG": 2}
    for m in res["panels"][0]["molecules"]:
        m["color_index"] = order[m["final_isoform_id"]]
    out = str(tmp_path / "explicit.pdf")
    assert render(res, out, isoform_track=True, iso_order="given") is not False
    assert _drawn_order(res["panels"][0]) == ["SHORT", "MID", "LONG"]


def test_render_refuses_an_unknown_iso_order(tmp_path):
    with pytest.raises(ValueError):
        render(_res(), str(tmp_path / "bad.pdf"), isoform_track=True, iso_order="count")


def test_main_passes_given_only_with_isoforms():
    """The one line in main() that chooses the order, read off the source rather than guessed."""
    src = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "altanalyze3", "components", "visualization", "isv_web",
                            "isv_web_remote.py")).read()
    assert 'iso_order=("given" if a.isoforms else "length")' in src


def test_sparse_panel_rows_are_visible(tmp_path):
    """A 5-read single-isoform panel must not collapse into one band.

    The panel height floor counted colour groups, not rows, so 5 reads of one isoform were scaled to
    13 px and each row drew 2.6 px tall. LRFN4 35292108.AML-7_1_CITE_GEX rendered as a single block.
    """
    from altanalyze3.components.visualization.isv_web.isv_web_remote import ROW_VIS
    res = _res()
    res["panels"][0]["molecules"] = [m for m in res["panels"][0]["molecules"]
                                     if m["final_isoform_id"] == "LONG"]          # 2 reads, 1 isoform
    out = str(tmp_path / "sparse.pdf")
    assert render(res, out, isoform_track=True) is not False
    rows = res["panels"][0]["_rows"]
    assert len(rows) == 2
    assert min(r["h"] for r in rows) >= ROW_VIS - 1e-9


def test_dense_panel_height_is_unchanged(tmp_path):
    """The repair must not fire where every row already fits: 100 rows would need 800 px."""
    from altanalyze3.components.visualization.isv_web.isv_web_remote import MAXPANEL, ROW_VIS
    res = _res()
    base = res["panels"][0]["molecules"][0]
    res["panels"][0]["molecules"] = [dict(base, isoform_id="r%d" % i) for i in range(100)]
    out = str(tmp_path / "dense.pdf")
    assert render(res, out, isoform_track=True) is not False
    rows = res["panels"][0]["_rows"]
    assert len(rows) == 100
    assert 100 * ROW_VIS > MAXPANEL
    assert sum(r["h"] for r in rows) <= MAXPANEL + 1e-6
