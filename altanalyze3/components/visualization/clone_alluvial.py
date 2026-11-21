"""Cell-state/clone alluvial; complete clinical roster, configurable state aggregation, assigned-cell denominators by default (WT included). Grafts retained in source and audited separately."""
from pathlib import Path
import pandas as pd, numpy as np, json, hashlib, re, textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.colors import to_rgb, to_hex
from matplotlib import font_manager

def render(figures, output, config=None, aggregation='configured', include_unassigned=False):
    if not isinstance(include_unassigned, bool):
        raise TypeError('include_unassigned must be a boolean')
    P = Path(figures)
    O = Path(output)
    O.mkdir(parents=True, exist_ok=True)
    config = json.loads(Path(config).read_text()) if config else {}
    patient = str(config.get('patient', pd.read_csv(P / 'cell_assignments.tsv.gz', sep='\t', dtype={'cell_id':str,'sample':str,'patient':str}).patient.iloc[0]))
    from altanalyze3.components.bam.clonal_inference.detection_plotting import require_arial
    require_arial()
    plt.rcParams.update({'font.family': 'Arial', 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    a = pd.read_csv(P / 'cell_assignments.tsv.gz', sep='\t', dtype={'cell_id':str,'sample':str,'patient':str})
    k = pd.read_csv(P / 'clone_key.tsv', sep='\t')
    key = k.set_index('node').clone
    if not (a.cell_id.is_unique and (not a.cell_state.isna().any())):
        raise ValueError('Integrity check failed: a.cell_id.is_unique and (not a.cell_state.isna().any())')
    if not a.is_xenograft.dtype == bool:
        raise ValueError('is_xenograft must be boolean')
    clinical = a.loc[~a.is_xenograft]
    meta = clinical.drop_duplicates('sample').sort_values('stage_order', kind='stable')
    samples = config.get('samples', meta['sample'].tolist())
    if not (set(samples) == set(clinical['sample']) and len(samples) == len(set(samples))):
        raise ValueError('All clinical specimens are required')
    expected = clinical.groupby('sample').size().reindex(samples).tolist()
    if not (not config.get('expected_counts') or expected == config['expected_counts']):
        raise ValueError("Integrity check failed: not config.get('expected_counts') or expected == config['expected_counts']")
    d = clinical.copy()
    d['clone'] = d.node.map(key)
    if not d.clone.notna().all():
        raise ValueError('Integrity check failed: d.clone.notna().all()')
    # Validate the complete saved baseline before applying the requested display filter.
    frozen = pd.read_csv(P / 'clone_frequencies.tsv', sep='\t')
    observed = d.groupby(['sample', 'clone']).size()
    for r in frozen[frozen['sample'].isin(samples)].itertuples():
        if observed.get((r.sample, r.clone), 0) != r.cells:
            raise ValueError('Saved clone counts do not match source assignments')
    input_counts = clinical.groupby('sample').size().reindex(samples)
    unresolved = d['clone'].eq('U') | d['node'].eq(-1)
    excluded = d.loc[unresolved & (not include_unassigned)].copy()
    excluded_counts = excluded.groupby('sample').size().reindex(samples, fill_value=0)
    if not include_unassigned:
        d = d.loc[~unresolved].copy()
        k = k.loc[~(k['clone'].eq('U') | k['node'].eq(-1))].copy()
    expected = d.groupby('sample').size().reindex(samples, fill_value=0).tolist()
    if any(n == 0 for n in expected):
        raise ValueError('No assigned cells remain in a clinical specimen; cannot normalize an empty sample')
    if not set(d.cell_id).isdisjoint(excluded.cell_id) or set(d.cell_id) | set(excluded.cell_id) != set(clinical.cell_id):
        raise ValueError('Reporting filter failed complete cell accounting')
    excluded.to_csv(O / 'excluded_unassigned_cells.tsv.gz', sep='\t', index=False)
    pd.DataFrame({'sample': samples, 'input_cells': input_counts.to_numpy(),
                  'excluded_unassigned_cells': excluded_counts.to_numpy(),
                  'plotted_cells': expected}).to_csv(O / 'cell_accounting.tsv', sep='\t', index=False)
    group_map = config.get('state_map', {})
    if aggregation == 'fine':
        group_map = {}
    d['original_cell_state'] = d.cell_state
    d['cell_state'] = d.original_cell_state.map(lambda s: group_map.get(s, s))
    if not d.cell_id.tolist() == clinical.loc[~clinical.cell_id.isin(excluded.cell_id)].cell_id.tolist():
        raise ValueError('Integrity check failed: d.cell_id.tolist() == clinical.loc[~clinical.cell_id.isin(excluded.cell_id)].cell_id.tolist()')
    mapping = d[['original_cell_state', 'cell_state']].drop_duplicates().sort_values(['cell_state', 'original_cell_state'])
    mapping.to_csv(O / 'cell_state_aggregation.tsv', sep='\t', index=False)
    order = config.get('state_order', [])
    if aggregation == 'fine':
        order = config.get('fine_state_order', [])
    states = [s for s in order if s in set(d.cell_state)] + sorted(set(d.cell_state) - set(order))
    if not (len(states) == len(set(states)) and set(states) == set(d.cell_state)):
        raise ValueError('Integrity check failed: len(states) == len(set(states)) and set(states) == set(d.cell_state)')
    clones = [x for x in k.clone if x != 'WT'] + (['WT'] if 'WT' in set(k.clone) else [])
    if not set(d.clone) <= set(clones):
        raise ValueError('Integrity check failed: set(d.clone) <= set(clones)')
    idx = pd.MultiIndex.from_product([states, clones], names=['cell_state', 'clone'])
    counts = pd.crosstab([d.cell_state, d.clone], d['sample']).reindex(index=idx, columns=samples, fill_value=0)
    freq = counts.div(expected, axis=1) * 100
    if not np.allclose(freq.sum(), 100):
        raise ValueError('Integrity check failed: np.allclose(freq.sum(), 100)')
    state_palette = {'HSC / MPP': '#677bc0', 'Myeloid intermediate': '#96775d', 'Mono / Mac lineage': '#765033', 'DC': '#b28adc', 'Lymphoid': '#1b8e91', 'Stromal': '#737d88', 'MKP': '#536ce5', 'MEP': '#ff0000', 'MEP-Eryth': '#c25d7c', 'Multilin': '#8450bf', 'MultiLin-GMP': '#a443b2', 'ERP': '#b43c66', 'LMPP': '#4ebf43', 'MDP': '#8fae94', 'preNeu / immNeu': '#ef8a00', 'ERP-1': '#b43c66', 'ERP-2': '#e58a9b', 'ERP-3': '#81244e', 'ERP-4': '#d55191', 'ERP-5': '#c56a63', 'ERP-6': '#edaa91', 'ERP-7': '#9b454a', 'ERP-8': '#cb747e', 'LMPP-1': '#168455', 'LMPP-1-cycling': '#93b83e', 'LMPP-2': '#4ebf43', 'Multilin-1': '#b28ad7', 'Multilin-2': '#8450bf', 'MultiLin-GMP-1': '#5c2a8d', 'MultiLin-GMP-2': '#a443b2', 'MultiLin-GMP-3': '#6861b8'}
    for family in [('ERP-',), ('LMPP-',), ('Multilin', 'MultiLin')]:
        members = [s for s in states if s.startswith(family)]
        if not all((s in state_palette for s in members)):
            raise ValueError('Integrity check failed: all((s in state_palette for s in members))')
        if not len({state_palette[s] for s in members}) == len(members):
            raise ValueError('Integrity check failed: len({state_palette[s] for s in members}) == len(members)')

    def base(s):
        if s in config.get('state_colors', {}):
            return config['state_colors'][s]
        if s in state_palette:
            return state_palette[s]
        if s.startswith('HSC'):
            return '#214ca4'
        if s.startswith('MPP-') and s != 'MPP-MEP':
            return '#3156ac'
        if s == 'MPP-MEP':
            return '#3895aa'
        if s.startswith(('Multilin', 'MultiLin')):
            return '#8146cb'
        if s.startswith('BMCP'):
            return '#cd65c8'
        if s.startswith('MEP'):
            return '#942149'
        if s.startswith('MKP'):
            return '#536ce5'
        if s.startswith('ERP'):
            return '#d8839b'
        if s.startswith('Erythroblast'):
            return '#d2ad85'
        if s.startswith('LMPP'):
            return '#4ebf43'
        if s == 'MDP-1':
            return '#a2bd31'
        if s.startswith('MDP'):
            return '#8fae94'
        if s.startswith(('pre-DC', 'cDC', 'pDC', 'ASDC')):
            return '#9070d1'
        if 'Mono' in s or s in ['Mac', 'cMOP']:
            return '#b8946c'
        if 'Neu' in s or s.startswith('Myeloid'):
            return '#b8a08b'
        if s.startswith(('Pro-B', 'B ', 'Plasma', 'CLP')):
            return '#53aaaa'
        return '#818b98'
    light = {'SC1': 0.76, 'SC2': 0.65, 'SC3': 0.52, 'SC4': 0.39, 'SC5': 0.27, 'SC6': 0.14, 'SC7': 0.03, 'WT': -0.19}
    for c, t in zip([x for x in clones if x not in light], np.linspace(0.7, 0.05, max(1, len([x for x in clones if x not in light])))):
        light[c] = float(t)

    def color(s, c):
        b = np.array(to_rgb(base(s)))
        v = light[c]
        return to_hex(b + (1 - b) * v if v >= 0 else b * (1 + v))
    bw = 0.96
    gap = bw * 0.5
    x = np.arange(len(samples)) * (bw + gap)
    v = freq.to_numpy()
    lower = 100 - np.cumsum(v, axis=0)
    upper = lower + v
    fig = plt.figure(figsize=(18.8, max(12.7, len(states)*.33)))
    ax = fig.add_axes([0.17, 0.25, 0.66, 0.65])
    ax.set_xlim(-0.7, x[-1] + 0.7)
    ax.set_ylim(0, 100)
    ax.axis('off')

    def curve(y, left, right, yl, yr):
        return [(left, yl), (left + (right - left) * 0.44, yl), (left + (right - left) * 0.56, yr), (right, yr)]

    def ribbon(x0, x1, lo0, hi0, lo1, hi1, col):
        verts = curve(None, x0, x1, hi0, hi1) + [(x1, lo1), (x1 - (x1 - x0) * 0.44, lo1), (x1 - (x1 - x0) * 0.56, lo0), (x0, lo0), (x0, hi0)]
        codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.LINETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.CLOSEPOLY]
        ax.add_patch(PathPatch(MPath(verts, codes), facecolor=col, edgecolor='white', lw=0.32))
    for i, (s, c) in enumerate(idx):
        col = color(s, c)
        for j in range(len(samples)):
            if v[i, j] > 0:
                ax.add_patch(Rectangle((x[j] - bw / 2, lower[i, j]), bw, v[i, j], facecolor=col, edgecolor='white', lw=0.35))
        for j in range(len(samples) - 1):
            if v[i, j] + v[i, j + 1] > 0:
                ribbon(x[j] + bw / 2, x[j + 1] - bw / 2, lower[i, j], upper[i, j], lower[i, j + 1], upper[i, j + 1], col)
    sf = freq.groupby(level=0, sort=False).sum().reindex(states)
    bottoms = 100 - sf.cumsum()
    tops = bottoms + sf
    for i, s in enumerate(states):
        y = bottoms.loc[s].values
        lw = 1.15 if sf.loc[s].max() >= 1.5 else 0.14
        for j in range(len(samples)):
            ax.plot([x[j] - bw / 2, x[j] + bw / 2], [y[j], y[j]], color='black', lw=lw)
        for j in range(len(samples) - 1):
            ax.add_patch(PathPatch(MPath(curve(None, x[j] + bw / 2, x[j + 1] - bw / 2, y[j], y[j + 1]), [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]), facecolor='none', edgecolor='black', lw=lw))
    ax.plot([-bw / 2, x[-1] + bw / 2], [100, 100], color='black', lw=1.4)
    for j in range(len(samples) - 1):
        ax.add_patch(Rectangle((x[j] + bw / 2, 0), gap, 100, facecolor='white', edgecolor='none', alpha=0.5, zorder=2.5))
    bar_width = 0.055
    bar_gap = 0.02
    clone_colors = k.set_index('clone')['color'].to_dict()
    bar_rows = []
    for state in states:
        y = float(bottoms.loc[state].iloc[0])
        h = float(sf.loc[state].iloc[0])
        if h > 0:
            ax.add_patch(Rectangle((x[0] - bw / 2 - bar_gap - bar_width, y), bar_width, h, facecolor=base(state), edgecolor='none', zorder=3))
            bar_rows.append({'side': 'left', 'sample': samples[0], 'cell_state': state, 'clone': '', 'bottom_percent': y, 'height_percent': h, 'color': base(state)})
    for i, (state, clone) in enumerate(idx):
        y = float(lower[i, -1])
        h = float(v[i, -1])
        if h > 0:
            ax.add_patch(Rectangle((x[-1] + bw / 2 + bar_gap, y), bar_width, h, facecolor=clone_colors[clone], edgecolor='none', zorder=3))
            bar_rows.append({'side': 'right', 'sample': samples[-1], 'cell_state': state, 'clone': clone, 'bottom_percent': y, 'height_percent': h, 'color': clone_colors[clone]})
    for side in ['left', 'right']:
        if not abs(sum((r['height_percent'] for r in bar_rows if r['side'] == side)) - 100) < 1e-09:
            raise ValueError("Integrity check failed: abs(sum((r['height_percent'] for r in bar_rows if r['side'] == side)) - 100) < 1e-09")
    pd.DataFrame(bar_rows).to_csv(O / 'aligned_color_bars.tsv', sep='\t', index=False)
    labels = []
    for i, (s, c) in enumerate(idx):
        j = int(v[i].argmax())
        h = v[i, j]
        if h >= 1.05:
            rgb = np.array(to_rgb(color(s, c)))
            textcolor = 'white' if rgb @ np.array([0.2126, 0.7152, 0.0722]) < 0.62 else '#222222'
            ax.text(x[j], (lower[i, j] + upper[i, j]) / 2, c, ha='center', va='center', fontsize=11 if h < 1.7 else 14, color=textcolor)
            labels.append({'cell_state': s, 'clone': c, 'sample': samples[j], 'percent_all': h})

    def side_labels(j, side):
        selected = list(states)
        targets = np.array([(bottoms.loc[s].iloc[j] + tops.loc[s].iloc[j]) / 2 for s in selected])
        positions = targets.copy()
        positions = np.clip(targets, 1, 99)
        for i in range(len(positions) - 2, -1, -1):
            positions[i] = max(positions[i], positions[i + 1] + min(2.05, 97 / max(1, len(positions)-1)))
        positions[0] = min(positions[0], 99)
        for i in range(1, len(positions)):
            positions[i] = min(positions[i], positions[i - 1] - min(2.05, 97 / max(1, len(positions)-1)))
        # Roundoff at 1.0 can be a few ulps negative after repeated spacing arithmetic.
        if not (positions[-1] >= 1-1e-9 and positions[0] <= 99+1e-9):
            raise ValueError('Integrity check failed: positions[-1] >= 1 and positions[0] <= 99')
        positions = np.clip(positions, 1, 99)
        edge = x[j] + side * (bw / 2 + bar_gap + bar_width)
        tx = edge + side * 0.13
        for s, y, t in zip(selected, positions, targets):
            ax.plot([edge + side * 0.015, edge + side * 0.075, tx - side * 0.018], [t, y, y], color=base(s), lw=0.65, clip_on=False)
            label = s.replace('cDC1', 'DC-1').replace('Myeloid intermediate ', 'Myeloid int. ').replace('Intermediate Mono-', 'Int. Mono-').replace('LMPP-1-cycling', 'LMPP-1 cycle')
            ax.text(tx, y, label, ha='left' if side > 0 else 'right', va='center', fontsize=13, color=base(s), clip_on=False)
    side_labels(0, -1)
    titles = config.get('titles', [str(meta.set_index('sample').loc[s, 'stage']) for s in samples])
    if not len(titles) == len(samples):
        raise ValueError('Integrity check failed: len(titles) == len(samples)')
    for j, t in enumerate(titles):
        ax.text(x[j], 102.8, t, ha='center', va='bottom', fontsize=19, fontweight='bold', linespacing=1.12)
        ax.text(x[j], -2.3, 'n = ' + format(expected[j], ','), ha='center', va='top', fontsize=12)
    fig.text(0.5, 0.976, config.get('title', patient + ' | Cell states and inferred clones'), ha='center', va='top', fontsize=21)
    fig.text(0.5, 0.185, 'SC labels are placed where each cell-state / clone combination is most frequent.', ha='center', fontsize=11)
    legend_y = 0.153
    for first in range(0, len(k), 3):
        lines = []
        for col, (_, r) in enumerate(k.iloc[first:first+3].iterrows()):
            label = r['clone'] + ' = ' + r.mutation_combination.replace('No direct or imputed mutant call', 'no mutant call').replace('|', ' | ')
            wrapped = textwrap.fill(label, 61)
            fig.text(0.07 + col * 0.305, legend_y, wrapped, fontsize=9.2, va='top')
            lines.append(len(wrapped.splitlines()))
        legend_y -= max(0.036, (max(lines) * 11.2 + 8) / (fig.get_figheight() * 72))
    caption_y = min(0.055, legend_y - 0.016)
    denominator = 'all profiled cells' if include_unassigned else 'assigned cells (SC and WT; U excluded)'
    exclusion_note = ('U retains unresolved mutant evidence. ' if include_unassigned else f'{len(excluded):,} U cells excluded from displayed frequencies and n; original assignments retained. ')
    caption = f'One patient; {len(samples)} clinical specimens. Band height = percentage of {denominator} at that timepoint. {exclusion_note} Hue identifies cell state; shades distinguish the listed SC clones and WT within a state. The left edge bar matches cell-state label colors; the right edge bar gives fishplot clone colors, aligned to every final-timepoint state/clone band. SC labels follow the supplied assignments; minimum supported lineages are not proven exclusive genotypes. Clone calls use direct and imputed evidence; WT denotes no mutant call. All ' + str(len(states)) + ' categories retain all ' + str(d.original_cell_state.nunique()) + ' original clinical cell states;  every resulting state/lineage band is labeled, and SC labels require ≥1.05% at their largest timepoint. Connecting regions are half a timepoint block wide and have a 50% white overlay. Curves connect category frequencies, not tracked cells or proven lineage transitions. MISTRG excluded as requested. Descriptive frequencies; no statistical test or uncertainty estimate.'
    fig.text(0.04, caption_y, textwrap.fill(caption, 220), fontsize=8.6, va='top')
    for ext in ['svg', 'pdf', 'png']:
        fig.savefig(O / (config.get('output_stem', 'state_clone_alluvial') + '.' + ext), dpi=210, bbox_inches='tight' if caption_y < 0.055 else None)
    plt.close(fig)
    counts.stack().rename('cells').reset_index().rename(columns={'level_2': 'sample'}).to_csv(O / 'cell_state_clone_counts.tsv', sep='\t', index=False)
    freq.stack().rename('percent_all_cells' if include_unassigned else 'percent_assigned_cells').reset_index().rename(columns={'level_2': 'sample'}).to_csv(O / 'cell_state_clone_frequencies.tsv', sep='\t', index=False)
    sf.T.to_csv(O / 'cell_state_percentages.tsv', sep='\t')
    pd.DataFrame(labels).to_csv(O / 'plot_labels.tsv', sep='\t', index=False)
    d.to_csv(O / 'clinical_cell_assignments.tsv.gz', sep='\t', index=False)
    k.to_csv(O / 'clone_key.tsv', sep='\t', index=False)
    pd.DataFrame([{'cell_state': s, 'base_color': base(s)} for s in states]).to_csv(O / 'cell_state_color_key.tsv', sep='\t', index=False)
    pd.DataFrame([{'cell_state': s, 'clone': c, 'color': color(s, c)} for s, c in idx]).to_csv(O / 'color_key.tsv', sep='\t', index=False)
    json.dump({'source': str(P / 'cell_assignments.tsv.gz'), 'source_sha256': hashlib.sha256((P / 'cell_assignments.tsv.gz').read_bytes()).hexdigest(), 'clinical_cells': len(d), 'clinical_states': len(states), 'original_clinical_states': int(d.original_cell_state.nunique()), 'aggregation_map': group_map, 'aggregation_authorization': config.get('authorization', 'Identity mapping; no state aggregation'), 'timepoint_cells': dict(zip(samples, expected)), 'MISTRG_separate_from_clinical_alluvial': int(a.is_xenograft.sum()), 'clone_counts_match_saved_frequencies': True, 'each_timepoint_sums_to_100_percent': True, 'all_clinical_cell_ids_preserved': len(excluded) == 0, 'all_input_cells_accounted_for': True, 'include_unassigned': include_unassigned, 'excluded_unassigned_cells': int(len(excluded)), 'input_timepoint_cells': {s:int(input_counts[s]) for s in samples}, 'excluded_by_sample': {s:int(excluded_counts[s]) for s in samples}, 'denominator': denominator, 'inference_changed': False, 'connector_to_timepoint_width_ratio': 0.5, 'connector_white_overlay_opacity': 0.5, 'right_cell_state_labels': False, 'edge_color_bars': 'left: first-timepoint cell states; right: final-timepoint state/clone bands; exact matching positions and heights', 'edge_color_bars_each_sum_to_100': True, 'ERP_LMPP_MultiLin_state_base_colors_unique_within_family': True}, (O / 'validation.json').open('w'), indent=2)
    (O / 'figure_legend.txt').write_text(caption + '\n')
    return json.loads((O / 'validation.json').read_text())

def add_arguments(p):
    p.add_argument('--figures', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--config')
    p.add_argument('--include-unassigned', action='store_true', default=False, help='Include U cells in bands and denominators; default excludes U and retains WT')
    p.add_argument('--aggregation', choices=['configured', 'fine'], default='configured')

def run(a):
    return render(a.figures, a.output, a.config, a.aggregation, a.include_unassigned)

def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    add_arguments(p)
    return run(p.parse_args(argv))
if __name__ == '__main__':
    main()
