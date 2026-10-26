"""Hashed iso2function protein architectures for optional ISV tracks.

Only saved projections and AltAnalyze coding-exon labels are used. Genomic and
protein coordinates remain separate; no sequence alignment is recomputed here.
The same vector primitives drive the headless PDF and interactive SVG views.
"""
from copy import deepcopy
import csv
import hashlib
import json
import math
from pathlib import Path

from .protein_features import (architecture_features, architecture_label, reference_architecture,
                               unplaced_summary)

PALETTE = ('#4477AA', '#EE6677', '#228833', '#AA3377', '#66CCEE', '#CCBB44',
           '#009988', '#BB5566', '#332288', '#AA4499', '#44AA99', '#999933')
SCHEMA = 'isv-protein-architecture-v1'


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_architectures(directory):
    """Read three small indexed exports; never scan native FASTAs/GFFs or matrices."""
    root = Path(directory)
    manifest = json.loads((root / 'figure_manifest.json').read_text())
    if manifest.get('status') != 'complete' or manifest.get('plot_style') != 'compact':
        raise ValueError('ISV protein tracks require complete compact protein-impact-plots output')
    if not manifest.get('coding_exon_map'):
        raise ValueError('ISV protein tracks require an AltAnalyze coding-exon map')
    for name in ('plot_data.json', 'coding_exon_coordinates.tsv', 'annotation_display_audit.tsv'):
        if _sha(root / name) != manifest['outputs'].get(name):
            raise ValueError('Changed or unpinned ISV architecture input: ' + name)
    plans = json.loads((root / 'plot_data.json').read_text())
    with (root / 'coding_exon_coordinates.tsv').open() as fh:
        exons = list(csv.DictReader(fh, delimiter='\t'))
    models = {}
    for row in exons:
        iid = row['isoform_id']
        meta = {k: row[k] for k in ('gene_symbol', 'protein_aa', 'reported_aa', 'protein_sha256', 'source_nmd')}
        model = models.setdefault(iid, dict(meta, coding_exons=[]))
        if any(model[k] != v for k, v in meta.items()):
            raise ValueError('Conflicting native protein metadata: ' + iid)
        if row['exon_label_status'] != 'assigned' or not row['exon_label']:
            raise ValueError('Unresolved AltAnalyze exon label: ' + iid)
        model['coding_exons'].append(dict(start=float(row['aa_start0']), end=float(row['aa_end0']),
                                          color=row['color'], label=row['exon_label']))
    for iid, model in models.items():
        model['protein_aa'] = int(model['protein_aa'])
        blocks = model['coding_exons']
        if (not blocks or abs(blocks[0]['start']) > 1e-7
                or abs(blocks[-1]['end'] - model['protein_aa']) > 1e-7
                or any(not math.isfinite(e['start'] + e['end']) or e['end'] <= e['start'] for e in blocks)
                or any(abs(a['end'] - b['start']) > 1e-7 for a, b in zip(blocks, blocks[1:]))):
            raise ValueError('Incomplete coding-exon coverage: ' + iid)
    return dict(schema=SCHEMA, root=str(root.resolve()), plans=plans, models=models,
                manifest_sha256=_sha(root / 'figure_manifest.json'),
                source_manifest_sha256=manifest['source_manifest_sha256'], selection=manifest['selection'])


def select_tracks(bundle, isoforms, gene_symbol=None):
    """Resolve every requested native ID; ambiguity/missing records are errors."""
    if not isoforms or len(isoforms) != len(set(isoforms)):
        raise ValueError('Protein tracks require unique native isoform IDs')
    tracks = []
    keys = []
    for iid in isoforms:
        if iid not in bundle['models']:
            raise ValueError('No audited protein architecture for selected ISV isoform: ' + iid)
        candidates = []
        for plan in bundle['plans']:
            pair = plan['pair']
            for side in ('A', 'B'):
                if pair['isoform_' + side] != iid:
                    continue
                placed, unplaced = architecture_features(plan['features'], side)
                features = []
                for item in placed:
                    f = item['feature']
                    key = (f['feature_type'], f['description'], f['provider'])
                    if key not in keys:
                        keys.append(key)
                    features.append(dict(key=list(key), feature_id=f['feature_id'],
                        label=architecture_label(f), description=f['description'], provider=f['provider'],
                        reference_protein=f['reference_protein'], segments=item['segments'], partial=item['partial'],
                        essential_class=f.get('essential_class', ''), class_color=f.get('class_color', '')))
                candidates.append(dict(features=features, unplaced=unplaced, gene_symbol=pair['gene_symbol'],
                                       unplaced_text=unplaced_summary(plan['features'], unplaced)))
        if not candidates or any(c != candidates[0] for c in candidates[1:]):
            raise ValueError('Missing or ambiguous pair projections for selected ISV isoform: ' + iid)
        track = dict(deepcopy(bundle['models'][iid]), isoform_id=iid, **candidates[0])
        if gene_symbol and track['gene_symbol'] != gene_symbol:
            raise ValueError('Protein architecture belongs to a different gene: ' + iid)
        for feature in track['features']:
            if any(a < 1 or b < a or b > track['protein_aa'] for a, b in feature['segments']):
                raise ValueError('Protein feature outside native protein: ' + iid)
        tracks.append(track)
    for track in tracks:
        for f in track['features']:
            # An essential class keeps its whitelist colour; other annotations cycle the palette.
            f['color'] = f['class_color'] if f.get('essential_class') else PALETTE[keys.index(tuple(f['key'])) % len(PALETTE)]
    # An essential class that NO selected isoform can place is shown on its own reference protein row,
    # never dropped to the unplaced text line (Nathan, 2026-10-25: IKZF1 drew no DNA-binding bar at all).
    seen = set()
    for plan in bundle['plans']:
        pair = plan['pair']
        if not any(pair['isoform_' + side] in isoforms for side in ('A', 'B')):
            continue
        for ref in reference_architecture(plan['features']):
            if ref['reference_protein'] in seen:
                continue
            seen.add(ref['reference_protein'])
            features = []
            for item in ref['features']:
                f = item['feature']
                features.append(dict(key=[f['feature_type'], f['description'], f['provider']],
                    feature_id=f['feature_id'], label=architecture_label(f), description=f['description'],
                    provider=f['provider'], reference_protein=f['reference_protein'], segments=item['segments'],
                    partial=False, essential_class=f.get('essential_class', ''),
                    class_color=f.get('class_color', ''), color=f.get('class_color', '') or PALETTE[0]))
            tracks.append(dict(gene_symbol=pair['gene_symbol'], isoform_id=ref['reference_protein'], reference=True,
                protein_aa=ref['extent'], reported_aa='', coding_exons=[], features=features,
                unplaced=[], unplaced_text='', display_color='#555555'))
    return tracks


def _text_width(text, size):
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    return TextPath((0, 0), text, size=size, prop=FontProperties(family='DejaVu Sans')).get_extents().width


def layout_tracks(tracks, width, font_scale=1.0):
    """Return reproducible point-space geometry, increasing height to avoid collisions."""
    if not tracks or not math.isfinite(width + font_scale) or width <= 0 or font_scale <= 0:
        raise ValueError('Protein track width and font scale must be positive and finite')
    fs = 8 * font_scale
    maximum = max(t['protein_aa'] for t in tracks)
    primitives = []
    def text(x, y, value, size=fs, color='#222222', **kwargs):
        primitives.append(dict(kind='text', x=x, y=y, text=value, size=size, color=color, **kwargs))
    def rect(x, y, w, h, color, **kwargs):
        primitives.append(dict(kind='rect', x=x, y=y, width=w, height=h, color=color, **kwargs))
    text(0, fs, 'Predicted protein features (N → C; amino-acid scale)', size=6 * font_scale)
    y = 20 * font_scale
    for track in tracks:
        lanes = []
        packed = []
        scale = width / maximum
        # Essential classes pack first (lane nearest the backbone) and are drawn last (on top).
        for f in sorted(track['features'], key=lambda f: (0 if f.get('essential_class') else 1, f['segments'][0][0], f['segments'][-1][1], f['feature_id'])):
            a, b = f['segments'][0][0] - 1, f['segments'][-1][1]
            label = f['label'] + (' (partial)' if f['partial'] else '')
            center = (a + b) * scale / 2
            half = _text_width(label, fs) / 2 + 3
            cx = max(half, min(width - half, center))
            left, right = min(a * scale, cx - half), max(b * scale, cx + half)
            lane = next((i for i, last in enumerate(lanes) if left > last + 5), len(lanes))
            if lane == len(lanes): lanes.append(right)
            else: lanes[lane] = right
            packed.append((f, label, a, b, cx, lane))
        lane_height = 15 * font_scale
        base = y + max(1, len(lanes)) * lane_height
        for f, label, a, b, cx, lane in sorted(packed, key=lambda x: 1 if x[0].get('essential_class') else 0):
            fy = base - (lane + 1) * lane_height
            tip = f['description'] + ' [' + f['provider'] + '; ' + f['reference_protein'] + ']'
            if f.get('essential_class'):
                tip = f['essential_class'] + ' (essential class): ' + tip
            rect(a * scale, fy + 2 * font_scale, (b - a) * scale, 4.5 * font_scale, f['color'],
                 role='protein_feature', feature_id=f['feature_id'], isoform_id=track['isoform_id'], tooltip=tip,
                 essential=bool(f.get('essential_class')))
            text(cx, fy - 1.5 * font_scale, label, anchor='middle', role='feature_label',
                 weight='bold' if f.get('essential_class') else 'normal')
        body_height = 12 * font_scale
        if not track['coding_exons']:
            # A reference-architecture row has no measured coding exons: draw its backbone explicitly.
            rect(0, base, track['protein_aa'] * scale, body_height, '#E8E8E8', edge='#777777',
                 role='reference_backbone', isoform_id=track['isoform_id'])
        small_exon_labels = []
        for e in track['coding_exons']:
            x = e['start'] * scale; w = (e['end'] - e['start']) * scale
            rect(x, base, w, body_height, e['color'], edge='#777777', role='coding_exon',
                 exon_label=e['label'], isoform_id=track['isoform_id'])
            fitted = None
            for rotation, size in ((0, fs), (0, fs*.875), (90, fs*.875), (90, fs*.8125)):
                tw = _text_width(e['label'], size)
                if ((tw < w - 1 and size < body_height - 1) if rotation == 0
                        else (size < w - 1 and tw < body_height - 1)):
                    fitted = rotation, size; break
            if fitted is None:
                # Keep the true coding width; retain the ID and exact AA offsets
                # in a legible key instead of shrinking text or losing the exon.
                small_exon_labels.append(e['label'] + ' [' + format(e['start'], 'g') + ',' +
                                         format(e['end'], 'g') + ')')
                continue
            rgb = [int(e['color'][i:i+2], 16)/255 for i in (1, 3, 5)]
            linear = [c/12.92 if c <= .04045 else ((c+.055)/1.055)**2.4 for c in rgb]
            luminance = sum(c*w for c, w in zip(linear, (.2126, .7152, .0722)))
            ink = 'black' if (luminance+.05)/.05 >= 1.05/(luminance+.05) else 'white'
            text(x+w/2, base+body_height/2, e['label'], size=fitted[1], color=ink,
                 anchor='middle', valign='center', rotation=fitted[0], role='exon_label')
        label = track['isoform_id'] + (' (reference)' if track.get('reference') else '')
        if track['reported_aa'] and int(track['reported_aa']) != track['protein_aa']:
            label += ' (' + str(track['protein_aa']) + ' aa; table: ' + track['reported_aa'] + ')'
        label_size = 9 * font_scale
        if _text_width(label, label_size) > width:
            raise ValueError('Native protein identifier does not fit; increase figure width: ' + label)
        half_label = _text_width(label, label_size)/2
        label_x = max(half_label, min(width-half_label, track['protein_aa'] * scale / 2))
        text(label_x, base+body_height+10*font_scale, label,
             size=label_size, color=track.get('display_color', '#222222'), anchor='middle', role='protein_id')
        y = base + body_height + 25*font_scale
        if small_exon_labels:
            line = 'Small coding segments (AA offsets):'
            for value in small_exon_labels:
                candidate = line + ' ' + value
                if _text_width(candidate, fs) > width:
                    text(0, y, line, size=fs, role='small_exon_key')
                    y += 12*font_scale
                    line = value
                else:
                    line = candidate
            text(0, y, line, size=fs, role='small_exon_key')
            y += 12*font_scale
        if track['unplaced']:
            # Wrap the list to the track width; every id stays on the figure.
            size = 6*font_scale
            head = 'Unplaced annotations:'
            lines, line = [], head
            for ident in (track.get('unplaced_text') or ', '.join(track['unplaced'])).split(', '):
                candidate = line + ' ' + ident + ','
                if _text_width(candidate, size) > width and line != head:
                    lines.append(line.rstrip(',')); line = '    ' + ident + ','
                else:
                    line = candidate
            lines.append(line.rstrip(','))
            for ln in lines:
                text(0, y, ln, size=size)
                y += 12*font_scale
    return dict(schema=SCHEMA, width=width, height=y, tracks=tracks, primitives=primitives,
                coordinate_system='Native protein N-to-C; separate amino-acid scale from genomic tracks',
                evidence='Predicted annotations; partial denotes sequence alignment coverage only')


def draw_matplotlib(ax, layout, x0, y0):
    """Draw the shared primitive layout into the native ISV point-space axes."""
    from matplotlib.patches import Rectangle
    artists = []
    for p in layout['primitives']:
        if p['kind'] == 'rect':
            artist = Rectangle((x0+p['x'], y0+p['y']), p['width'], p['height'],
                               facecolor=p['color'], edgecolor=p.get('edge', p['color']), linewidth=.5)
            ax.add_patch(artist)
        else:
            artist = ax.text(x0+p['x'], y0+p['y'], p['text'], fontsize=p['size'], color=p['color'],
                             ha={'middle':'center'}.get(p.get('anchor'), 'left'),
                             va=p.get('valign', 'bottom'), rotation=p.get('rotation', 0), family='DejaVu Sans',
                             fontweight=p.get('weight', 'normal'))
        artist.set_gid(p.get('role', 'protein_annotation'))
        artists.append(artist)
    return artists
