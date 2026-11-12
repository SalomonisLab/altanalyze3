#!/usr/bin/env python3
"""Catalog flow diagram: reads -> iso-structures -> collapsed isoforms -> UMI>=3 classes -> protein coding.

Layout supplied by Nathan Salomonis (2026-11-11). The linear branch (known and novel isoforms) and the
chimeric branch (intra- and inter-chromosomal trans-splicing, circRNA) each show their per-source
iso-structures and their collapsed isoforms or junctions. Every number comes from a two-column
table (``box``, ``value``) written by ``altanalyze3.components.long_read.noncanonical.catalog_flow``;
nothing is computed here.

    python -m altanalyze3 visualize catalog-flow --values CATALOG_FLOW_VALUES.tsv --out-prefix <dir>/catalog_flow

Output: <prefix>.pdf (vector, live Arial text) and <prefix>.png.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

REQUIRED = ('umi_total', 'linear_structures', 'known', 'novel', 'known_protein_coding', 'novel_protein_coding',
            'chimeric_structures', 'chimeric_collapsed', 'intra_trans_splice', 'inter_trans_splice', 'circrna',
            'intra_trans_splice_protein_coding', 'inter_trans_splice_protein_coding')
COLORS = dict(umi='#E5E5E5', linear_structures='#F2F5BA', chimeric_structures='#F7E6F5',
              linear_collapsed='#A6F5C3', chimeric_collapsed='#C97BEA', known='#5C5CF2', novel='#E86B6B',
              intra='#E05AB4', inter='#6A1FCC', circrna='#5CC08A', protein_fill='#F7F7F7',
              protein_edge='#E03030', edge='#6E6E6E', text='#000000', light_text='#FFFFFF')


def read_values(path):
    """Two-column table: box, value. Integers are parsed; other values stay text."""
    values = {}
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames is None or not {'box', 'value'} <= set(reader.fieldnames):
            raise ValueError(f'{path}: columns box and value are required')
        for row in reader:
            key = row['box'].strip()
            if key in values:
                raise ValueError(f'{path}: duplicate box {key}')
            raw = row['value'].strip()
            values[key] = int(raw) if raw.lstrip('-').isdigit() else raw
    return values


def collapsed_label(values):
    """One exact number, or an explicit range when only bounds are known."""
    if 'linear_collapsed' in values and values['linear_collapsed'] != '':
        return f"{values['linear_collapsed']:,}" if isinstance(values['linear_collapsed'], int) else str(values['linear_collapsed'])
    lo, hi = values.get('linear_collapsed_low'), values.get('linear_collapsed_high')
    if not isinstance(lo, int) or not isinstance(hi, int) or lo > hi:
        raise ValueError('linear_collapsed, or integer linear_collapsed_low <= linear_collapsed_high, is required')
    return f'{lo:,}' if lo == hi else f'{lo:,}–{hi:,}'


def render(values, out_prefix):
    missing = [k for k in REQUIRED if k not in values]
    if missing:
        raise ValueError('Missing catalog flow values: ' + ', '.join(missing))
    for k in REQUIRED:
        if not isinstance(values[k], int):
            raise ValueError(f'{k} must be an integer, got {values[k]!r}')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, Polygon
    from .figure_style import apply_rcparams, strip_masks

    apply_rcparams(9)
    W, H = 2000.0, 910.0                       # drawing units; 200 units per inch
    fig = plt.figure(figsize=(W / 200.0, H / 200.0))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.set_aspect('equal'); ax.axis('off')
    drawn = []

    renderer = fig.canvas.get_renderer()

    def fit(artist, width):
        """Shrink a label until it fits the box width; installed fonts differ in width (Arial vs DejaVu)."""
        limit = (ax.transData.transform((width, 0))[0] - ax.transData.transform((0, 0))[0]) * 0.92
        while artist.get_window_extent(renderer).width > limit and artist.get_fontsize() > 8:
            artist.set_fontsize(artist.get_fontsize() - 0.5)

    def box(x0, y0, x1, y1, face, title, value, light=False, tsize=19, vsize=20):
        ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle='round,pad=0,rounding_size=34',
                                    facecolor=face, edgecolor=COLORS['edge'], linewidth=1.5))
        col = COLORS['light_text'] if light else COLORS['text']
        cx = (x0 + x1) / 2.0
        text = value if isinstance(value, str) else f'{value:,}'
        fit(ax.text(cx, y0 + (y1 - y0) * 0.30, title, ha='center', va='center', fontsize=tsize, color=col), x1 - x0)
        fit(ax.text(cx, y0 + (y1 - y0) * 0.72, text, ha='center', va='center', fontsize=vsize, fontweight='bold',
                    color=col), x1 - x0)
        drawn.extend([title, text])

    def arrow(x0, y0, x1, y1):
        a = math.atan2(y1 - y0, x1 - x0); length, half = 26.0, 13.0
        xb, yb = x1 - length * math.cos(a), y1 - length * math.sin(a)
        ax.plot([x0, xb], [y0, yb], color=COLORS['text'], linewidth=2.6, solid_capstyle='butt')
        px, py = -math.sin(a) * half, math.cos(a) * half
        ax.add_patch(Polygon([(x1, y1), (xb + px, yb + py), (xb - px, yb - py)], closed=True,
                             facecolor=COLORS['text'], edgecolor='none'))

    def label(x, y, text, size=17):
        ax.text(x, y, text, ha='center', va='center', fontsize=size, color=COLORS['text'])
        drawn.append(text)

    umi_title = values.get('umi_label', 'UMI')
    support = values.get('support_label', 'UMI≥3')
    box(548, 30, 1225, 152, COLORS['umi'], umi_title, values['umi_total'])
    box(42, 203, 680, 325, COLORS['linear_structures'], 'Linear Iso-Structures', values['linear_structures'])
    box(1097, 203, 1735, 325, COLORS['chimeric_structures'], 'Chimeric Iso-Structures', values['chimeric_structures'])
    box(42, 403, 680, 525, COLORS['linear_collapsed'], 'Collapsed Isoforms', collapsed_label(values))
    box(1097, 400, 1735, 522, COLORS['chimeric_collapsed'], 'Collapsed Isoforms', values['chimeric_collapsed'])
    box(45, 613, 320, 737, COLORS['known'], 'Known', values['known'], light=True)
    box(352, 613, 647, 737, COLORS['novel'], 'Novel', values['novel'], light=True)
    box(718, 613, 1170, 737, COLORS['intra'], 'Intra-Trans-splice', values['intra_trans_splice'], light=True)
    box(1190, 613, 1640, 737, COLORS['inter'], 'Inter-Trans-splice', values['inter_trans_splice'], light=True)
    box(1662, 613, 1957, 737, COLORS['circrna'], 'circRNA', values['circrna'], light=True)

    arrow(548, 142, 482, 190); arrow(1228, 142, 1294, 190)
    arrow(360, 330, 360, 396); arrow(1415, 330, 1415, 393)
    arrow(190, 532, 190, 606); arrow(505, 532, 505, 606); label(345, 557, support)
    arrow(1205, 530, 1100, 606); arrow(1415, 530, 1415, 606); arrow(1640, 530, 1746, 606)
    label(1306, 557, support); label(1543, 557, support)

    ax.add_patch(FancyBboxPatch((45, 812), 1515, 90, boxstyle='round,pad=0,rounding_size=40',
                                facecolor=COLORS['protein_fill'], edgecolor=COLORS['protein_edge'], linewidth=2.0))
    for x, key in ((190, 'known_protein_coding'), (505, 'novel_protein_coding'),
                   (940, 'intra_trans_splice_protein_coding'), (1415, 'inter_trans_splice_protein_coding')):
        arrow(x, 744, x, 805)
        text = f'{values[key]:,}'
        ax.text(x, 858, text, ha='center', va='center', fontsize=20, fontweight='bold', color=COLORS['text'])
        drawn.append(text)

    strip_masks(fig)
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    pdf, png = out_prefix.with_suffix('.pdf'), out_prefix.with_suffix('.png')
    fig.savefig(pdf, format='pdf')
    fig.savefig(png, format='png', dpi=200, facecolor='#FFFFFF')
    plt.close(fig)
    return dict(pdf=str(pdf), png=str(png), drawn_text=drawn)


def add_arguments(parser):
    parser.add_argument('--values', required=True, help='Two-column table (box, value) from the catalog_flow stage')
    parser.add_argument('--out-prefix', required=True, help='Output path prefix; .pdf and .png are written')


def run(args):
    result = render(read_values(args.values), args.out_prefix)
    print(result['pdf'])
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    add_arguments(parser)
    return run(parser.parse_args(argv))


if __name__ == '__main__':
    main()
