import re
import xml.etree.ElementTree as ET

import pandas as pd
import pytest

from altanalyze3.components.visualization import chromosome_ideogram as ci


@pytest.fixture(scope="module")
def bands():
    return ci.load_cytobands()


def _row(**kw):
    base = {
        "chimeric_rna": "FOSB>SF1", "donor_gene": "FOSB", "acceptor_gene": "SF1",
        "junction": "chr19:45468712:+ > chr11:64776626:-",
        "overall_bestp_source": "Signature_GSE242273", "overall_bestp_p": "0.0001157296291649",
        "overall_bestp_q": "0.0983616979507364", "overall_bestp_z": "3.855", "overall_bestp_dataset": "Plasma_Bcell",
        "overall_bestp_cell_type_or_tissue": "Plasma cells (Blood)",
        "overall_bestp_donor_side_contact": "chr19:45000001-46000000",
        "overall_bestp_acceptor_side_contact": "chr11:64000001-65000000", "overall_bestp_tied_datasets": "1.0",
        "overall_bestp_significant": "p<0.05",
        "overall_bestq_source": "Signature_GSE242273", "overall_bestq_p": "0.0001639028876974",
        "overall_bestq_q": "0.0437906715191017", "overall_bestq_z": "3.769", "overall_bestq_dataset": "Islet_Lawlor",
        "overall_bestq_cell_type_or_tissue": "Islet (Pancreas)",
        "overall_bestq_donor_side_contact": "chr19:45000001-46000000",
        "overall_bestq_acceptor_side_contact": "chr11:64000001-65000000", "overall_bestq_tied_datasets": "1.0",
        "overall_bestq_significant": "q<0.05",
        "kinnex_junction_reads_sum_all_samples": "63", "kinnex_samples_positive": "10",
        "kinnex_patients_positive": "10", "ogawa_junction_fragments_sum_all_362_samples": "6",
        "ogawa_samples_positive_of_362": "6",
    }
    base.update(kw)
    return pd.Series(base)


def test_band_table_tiles_every_chromosome(bands):
    assert len(bands["chr19"]) == 19 and bands["chr19"][-1].end == 58617616
    assert bands["chr11"][-1].end == 135086622


def test_band_lookup_known_genes(bands):
    assert ci.band_at(bands["chr19"], 45468712).name == "q13.32"
    assert ci.band_at(bands["chr11"], 64776626).name == "q13.1"
    # 1-based coordinate equal to a 0-based band end belongs to that band; the next base starts the next band
    assert ci.band_at(bands["chr19"], 44700000).name == "q13.31"
    assert ci.band_at(bands["chr19"], 44700001).name == "q13.32"
    with pytest.raises(ValueError):
        ci.band_at(bands["chr19"], 58617617)


def test_gap_in_band_table_is_refused(tmp_path):
    p = tmp_path / "bad.txt"
    p.write_text("chrX\t0\t10\tp1\tgneg\nchrX\t11\t20\tq1\tgneg\n")
    with pytest.raises(ValueError):
        ci.load_cytobands(str(p))


def test_contact_rule():
    assert ci.select_contact(_row()).dataset == "Islet_Lawlor"
    r = ci.select_contact(_row(overall_bestq_significant="q>=0.05"))
    assert r.dataset == "Plasma_Bcell" and r.rule.startswith("best by p")
    empty = {k: "" for k in _row().index if k.startswith("overall_best")}
    assert ci.select_contact(_row(**empty)) is None


def test_hatch_segments_stay_inside_band(bands):
    ideo = ci.Ideogram("chr19", bands["chr19"], 50, 40, 400 / 58617616).build()
    n = 0
    for s in ideo.shapes:
        if s["band"].stain in ci.HATCHED:
            for a, b in ideo.hatch_segments(s["poly"]):
                assert ci.point_in_convex(s["poly"], a, tol=1e-3) and ci.point_in_convex(s["poly"], b, tol=1e-3)
                n += 1
    assert n > 10


def test_svg_is_plain_vector(bands):
    svg, values = ci.render_chimera(_row(), bands)
    root = ET.fromstring(svg.encode())
    tags = {el.tag.split('}')[-1] for el in root.iter()}
    assert not tags & {"clipPath", "mask", "pattern", "image", "use", "filter"}
    for el in root.iter():
        for attr in ("fill", "stroke"):
            v = el.get(attr)
            if v is not None:
                assert v == "none" or re.fullmatch(r"#[0-9A-F]{6}", v), (attr, v)
        if el.tag.endswith("text"):
            assert el.get("font-family").startswith("Arial")
    assert values["donor_band"] == "19q13.32" and values["acceptor_band"] == "11q13.1"
    assert values["contact_dataset"] == "Islet_Lawlor"


def test_same_chromosome_refused(bands):
    with pytest.raises(ValueError):
        ci.render_chimera(_row(junction="chr19:45468712:+ > chr19:1000:+"), bands)


def test_horizontal_layout_is_plain_vector_and_places_genes(bands):
    svg, values = ci.render_chimera_horizontal(_row(), bands)
    root = ET.fromstring(svg.encode())
    tags = {el.tag.split('}')[-1] for el in root.iter()}
    assert not tags & {"clipPath", "mask", "pattern", "image", "use", "filter"}
    texts = [el.text for el in root.iter() if el.tag.endswith("text")]
    assert "SF1  11q13.1" in texts and "FOSB  19q13.32" in texts and "Chr 11" in texts and "Chr 19" in texts
    assert float(root.get("width")) == 720.0
    assert values["layout"] == "horizontal" and values["contact_dataset"] == "Islet_Lawlor"


def test_horizontal_geometry_matches_vertical_bands(bands):
    v = ci.Ideogram("chr19", bands["chr19"], 50, 40, 0.00001).build()
    h = ci.Ideogram("chr19", bands["chr19"], 50, 40, 0.00001, orientation="horizontal").build()
    for sv, sh in zip(v.shapes, h.shapes):
        assert [(round(y, 6), round(x, 6)) for x, y in sv["poly"]] == [(round(x, 6), round(y, 6)) for x, y in sh["poly"]]


def test_horizontal_strip_without_title_legend_caption(bands):
    svg, _ = ci.render_chimera_horizontal(_row(), bands, title=False, legend=False, caption=False)
    root = ET.fromstring(svg.encode())
    texts = [el.text for el in root.iter() if el.tag.endswith("text")]
    assert "FOSB>SF1" not in texts and "gneg" not in texts
    assert float(root.get("height")) < 140


def _cubic(p, t):
    return tuple((1 - t) ** 3 * a + 3 * (1 - t) ** 2 * t * b + 3 * (1 - t) * t ** 2 * c + t ** 3 * d for a, b, c, d in zip(*p))


def test_horizontal_junction_line_lies_inside_contact_ribbon(bands):
    svg, _ = ci.render_chimera_horizontal(_row(), bands)
    root = ET.fromstring(svg.encode())
    paths = [el for el in root.iter() if el.tag.endswith("path")]
    ribbon = [el for el in paths if el.get("fill") == ci.CONTACT_COLOR and el.get("fill-opacity") and " C" in el.get("d")]
    line = [el for el in paths if el.get("stroke") == ci.INK and el.get("fill") == "none" and " C" in el.get("d")]
    assert len(ribbon) == 1 and len(line) == 1
    num = lambda d: [float(v) for v in re.findall(r"-?\d+\.\d+", d)]
    r = num(ribbon[0].get("d"))  # M d0,yB C .. .. a0,yT L a1,yT C .. .. d1,yB Z
    left = [(r[0], r[1]), (r[2], r[3]), (r[4], r[5]), (r[6], r[7])]
    assert len(r) == 16, len(r)
    right = [(r[14], r[15]), (r[12], r[13]), (r[10], r[11]), (r[8], r[9])]  # second curve, reversed to run donor to acceptor
    j = num(line[0].get("d"))
    junction = [(j[0], j[1]), (j[2], j[3]), (j[4], j[5]), (j[6], j[7])]
    for i in range(101):
        t = i / 100
        xl, yl = _cubic(left, t)
        xr, yr = _cubic(right, t)
        xj, yj = _cubic(junction, t)
        assert abs(yl - yj) < 1e-6 and abs(yr - yj) < 1e-6
        assert min(xl, xr) - 1e-6 <= xj <= max(xl, xr) + 1e-6, (t, xl, xj, xr)


def test_junction_option_restricts_rows(tmp_path):
    rows = [_row(), _row(junction="chr19:45468712:+ > chr11:64000000:-")]
    tab = tmp_path / "t.tsv"
    pd.DataFrame(rows).to_csv(tab, sep="\t", index=False)
    out = tmp_path / "o"
    ci.main(["--junctions", str(tab), "--chimera", "FOSB>SF1", "--junction", "chr19:45468712:+ > chr11:64776626:-",
             "--no-caption", "--outdir", str(out)])
    v = pd.read_csv(out / "IDEOGRAM_VALUES.tsv", sep="\t")
    assert len(v) == 1 and v.loc[0, "acceptor_pos"] == 64776626
    assert (out / "FOSB_to_SF1_ideogram_horizontal.svg").exists()
