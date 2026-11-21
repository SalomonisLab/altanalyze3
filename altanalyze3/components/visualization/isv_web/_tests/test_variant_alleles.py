"""--variant_alleles: one allele mark per read at a known SNV, from read_alleles.tsv."""
from copy import deepcopy
from unittest.mock import patch

import pytest
from matplotlib.patches import Rectangle

from altanalyze3.components.visualization.isv_web import isv_web_remote as remote

HEADER = ("variant_id\tchrom\tpos\tref\talt\tuid\tmolecule_id\tread_name\tcb_bam\tcb_rc\tflag\tsecondary\t"
          "supplementary\tmapq\tbase\tallele_class\tcigar_op\n")


def _row(mid, cls, vid="rsX", pos=150, secondary=0):
    base = {"REF": "A", "ALT": "G", "OTHER": "T"}.get(cls, "")
    return "{}\tchr1\t{}\tA\tG\tS1\t{}\tmolecule/{}\tAC\tGT\t0\t{}\t0\t60\t{}\t{}\t=\n".format(
        vid, pos, mid, mid.split(".")[0], secondary, base, cls)


def table(tmp_path, rows):
    p = tmp_path / "read_alleles.tsv"
    p.write_text(HEADER + "".join(rows))
    return p


def response(n=8):
    mols = []
    for i in range(n):
        mols.append(dict(isoform_id="{}.S1".format(i), final_isoform_id="A" if i < n // 2 else "B",
                         known=True, count=1, cluster_index=0 if i < n // 2 else 1,
                         exon_segments=[dict(start=100, end=200, type="E"), dict(start=300, end=400, type="E")]))
    return dict(gene="G", symbol="G", gene_model=dict(chrom="chr1", strand="+", extent=[100, 400],
                blocks=[dict(start=100, end=200, exon_id="E1"), dict(start=300, end=400, exon_id="E2")]),
                panels=[dict(condition="P1", n_molecules=n, molecules=mols)])


CLASSES = ["REF", "ALT", "REF", "SPLICE", "ALT", "OTHER", "NOT_SPANNED"]   # molecule 7 has no row


def _variant(tmp_path):
    return remote.load_variant_alleles(table(tmp_path, [_row("{}.S1".format(i), c) for i, c in enumerate(CLASSES)]))


def _marks(fig, variant):
    """Base letter of every filled mark, or HOLLOW for a spliced/deleted read."""
    fills = {c.upper(): b for b, c in variant["colors"].items()}
    out = []
    for p in fig.axes[0].patches:
        if isinstance(p, Rectangle) and p.get_zorder() == 6:
            if p.get_facecolor()[3] == 0:
                out.append("HOLLOW")
                continue
            rgb = "#%02X%02X%02X" % tuple(int(round(255 * c)) for c in p.get_facecolor()[:3])
            out.append(fills.get(rgb, rgb))
    return out


def _lab(h):
    h = h.lstrip("#")
    lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in (int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))]
    r, g, b = lin
    x, y, z = ((0.4124 * r + 0.3576 * g + 0.1805 * b) / 0.95047, 0.2126 * r + 0.7152 * g + 0.0722 * b,
               (0.0193 * r + 0.1192 * g + 0.9505 * b) / 1.08883)
    f = lambda t: t ** (1 / 3) if t > 0.008856 else 7.787 * t + 16 / 116
    return 116 * f(y) - 16, 500 * (f(x) - f(y)), 200 * (f(y) - f(z))


def _delta_e(a, b):
    return sum((p - q) ** 2 for p, q in zip(_lab(a), _lab(b))) ** 0.5


def test_one_mark_per_read_with_a_class_and_none_otherwise(tmp_path):
    v = _variant(tmp_path)
    res = response()
    stats = remote.attach_variant_alleles(res, v)
    assert dict(stats["P1"]) == {"REF": 2, "ALT": 2, "SPLICE": 1, "OTHER": 1, "NOT_SPANNED": 1, "absent": 1}
    figs = []
    with patch.object(remote.plt, "close", side_effect=lambda f: figs.append(f)):
        assert remote.render(res, str(tmp_path / "v.pdf"), show_read_labels=False, variant=v)
    marks = sorted(_marks(figs[-1], v))
    # IGV base colours: A (ref) twice, G (alt) twice, T (other base) once; the spliced read is a hollow box
    assert marks == sorted(["A", "G", "A", "HOLLOW", "G", "T"])
    texts = [t.get_text() for t in figs[-1].axes[0].texts]
    assert any(t.startswith("rsX chr1:150 A>G") for t in texts)
    assert any(t == "rsX: A 2 · G 2 · other 1 · spliced/del 1 · no site 2" for t in texts)


def test_no_variant_draws_no_mark(tmp_path):
    v = _variant(tmp_path)
    res = response()
    remote.attach_variant_alleles(res, v)
    figs = []
    with patch.object(remote.plt, "close", side_effect=lambda f: figs.append(f)):
        remote.render(res, str(tmp_path / "plain.pdf"), show_read_labels=False)
    assert _marks(figs[-1], v) == []


def test_variant_outside_the_window_is_refused(tmp_path):
    v = remote.load_variant_alleles(table(tmp_path, [_row("0.S1", "REF", pos=5000)]))
    res = response()
    remote.attach_variant_alleles(res, v)
    with pytest.raises(ValueError, match="outside the drawn window"):
        remote.render(res, str(tmp_path / "x.pdf"), show_read_labels=False, variant=v)


def test_variant_sort_puts_ref_before_alt_inside_each_isoform(tmp_path):
    v = _variant(tmp_path)
    res = response()
    remote.attach_variant_alleles(res, v, sort_by_allele=True)
    for m in res["panels"][0]["molecules"]:
        m["merged"] = remote.merge_segments(m["exon_segments"])
    order = [m["variant_class"] for m in sorted(res["panels"][0]["molecules"], key=remote.row_sort_key)]
    # isoform A (molecules 0-3) then isoform B (4-7); allele order inside each, absent last
    assert order == ["REF", "REF", "ALT", "SPLICE", "ALT", "OTHER", "NOT_SPANNED", None]
    plain = response()
    for m in plain["panels"][0]["molecules"]:
        m["merged"] = remote.merge_segments(m["exon_segments"])
    assert [m["isoform_id"] for m in sorted(plain["panels"][0]["molecules"], key=remote.row_sort_key)] == \
        sorted(m["isoform_id"] for m in plain["panels"][0]["molecules"])


def test_primary_row_wins_and_conflicting_primaries_stop(tmp_path):
    v = remote.load_variant_alleles(table(tmp_path, [_row("1.S1", "ALT", secondary=1), _row("1.S1", "REF")]))
    assert v["calls"]["1.S1"] == "REF"
    with pytest.raises(ValueError, match="two primary rows"):
        remote.load_variant_alleles(table(tmp_path, [_row("1.S1", "ALT"), _row("1.S1", "REF")]))


def test_two_variants_need_an_id(tmp_path):
    p = table(tmp_path, [_row("1.S1", "ALT", vid="rsA"), _row("1.S1", "REF", vid="rsB")])
    with pytest.raises(ValueError, match="--variant_id"):
        remote.load_variant_alleles(p)
    assert remote.load_variant_alleles(p, "rsB")["calls"] == {"1.S1": "REF"}


def test_colors_must_be_hex_and_name_a_base():
    assert remote.parse_variant_colors("A=#112233") == {"A": "#112233"}
    with pytest.raises(ValueError):
        remote.parse_variant_colors("A=green")
    with pytest.raises(ValueError):
        remote.parse_variant_colors("REF=#000000")


def test_igv_base_colours_match_igv_desktop_defaults():
    # igvteam/igv main, src/main/java/org/igv/renderer/SequenceRenderer.java and Globals.DARK_MODE_BLUE
    assert remote.IGV_BASE_COLORS == {"A": "#009600", "C": "#0096FF", "G": "#D17105", "T": "#FF0000", "N": "#808080"}


def test_isoform_fills_never_use_a_base_colour(tmp_path):
    for fill in remote.VARIANT_SAFE_ISOFORM_COLORS:
        assert min(_delta_e(fill, c) for c in remote.IGV_BASE_COLORS.values()) >= 40, fill
    v = _variant(tmp_path)
    res = response()
    for m in res["panels"][0]["molecules"]:
        m["display_color"] = "#2f8a4e"          # a caller green near A must not survive in variant mode
    remote.attach_variant_alleles(res, v)
    figs = []
    with patch.object(remote.plt, "close", side_effect=lambda f: figs.append(f)):
        remote.render(res, str(tmp_path / "safe.pdf"), show_read_labels=False, isoform_track=True, variant=v)
    safe = {c.upper() for c in remote.VARIANT_SAFE_ISOFORM_COLORS}
    exon_fills = {"#%02X%02X%02X" % tuple(int(round(255 * c)) for c in p.get_facecolor()[:3])
                  for p in figs[-1].axes[0].patches
                  if isinstance(p, Rectangle) and p.get_zorder() != 6 and p.get_facecolor()[3] > 0}
    assert "#2F8A4E" not in exon_fills
    assert {"#7A4FB5", "#BCBD22"} <= exon_fills and exon_fills & safe == {"#7A4FB5", "#BCBD22"}


def test_cli_flags_need_the_table_and_render_input_refuses_the_display_filter(tmp_path):
    a = remote.build_arg_parser().parse_args(["--run_dir", str(tmp_path), "--gene", "G", "--variant_sort"])
    with pytest.raises(SystemExit):
        remote._prepare_variant(a)
    p = table(tmp_path, [_row("0.S1", "REF")])
    with pytest.raises(SystemExit, match="variant_spanning_only"):
        remote.run_from_args(remote.build_arg_parser().parse_args(
            ["--run_dir", str(tmp_path), "--gene", "G", "--render-input", str(tmp_path / "none.json"),
             "--variant_alleles", str(p), "--variant_spanning_only"]))


def test_render_input_path_passes_the_variant_to_render(tmp_path):
    import json
    res = response()
    src = tmp_path / "in.json"
    src.write_text(json.dumps(res))
    p = table(tmp_path, [_row("{}.S1".format(i), c) for i, c in enumerate(CLASSES)])
    args = remote.build_arg_parser().parse_args(["--run_dir", str(tmp_path), "--gene", "G", "--render-input", str(src),
                                                 "--variant_alleles", str(p)])
    with patch.object(remote, "render", return_value=True) as draw, pytest.raises(SystemExit):
        remote.run_from_args(args)
    v = draw.call_args.kwargs["variant"]
    assert v["id"] == "rsX" and v["pos"] == 150
    drawn = draw.call_args.args[0]["panels"][0]["molecules"]
    assert [m["variant_class"] for m in drawn] == CLASSES + [None]
    before = deepcopy(response())
    assert [m["isoform_id"] for m in drawn] == [m["isoform_id"] for m in before["panels"][0]["molecules"]]


def test_figure_mode_text_stays_on_the_page_and_headers_do_not_collide(tmp_path):
    v = _variant(tmp_path)
    res = response()
    res["panels"][0]["condition"] = "All_cell_states | Healthy"
    remote.attach_variant_alleles(res, v)
    figs = []
    with patch.object(remote.plt, "close", side_effect=lambda f: figs.append(f)):
        remote.render(res, str(tmp_path / "fig.pdf"), show_read_labels=False, font_scale=2.6, show_ruler=False,
                      trim_to_reads=True, title_gene="G", structure_labels_above=True, show_read_counts=False,
                      show_structure_header=False, isoform_track=True, variant=v)
    fig = figs[-1]
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    page = fig.bbox
    texts = [t for t in fig.axes[0].texts if t.get_text()]
    for t in texts:
        b = t.get_window_extent(rend)
        assert b.x0 >= page.x0 - 1 and b.x1 <= page.x1 + 1, t.get_text()
    header = [t for t in texts if t.get_text().startswith("rsX:") or t.get_text() == "All_cell_states | Healthy"]
    assert len(header) == 2
    assert not header[0].get_window_extent(rend).overlaps(header[1].get_window_extent(rend))
