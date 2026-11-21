"""Approved compact source heatmap: 600 dpi pixel layer, 0.15 pt marks, editable labels. No inference or cell filtering."""
from pathlib import Path
import argparse, json, hashlib, textwrap
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch, Rectangle
from altanalyze3.components.bam.clonal_inference.detection_plotting import require_arial
from altanalyze3.components.bam.clonal_inference.integrity import sha256

def render(input_dir, figures, output, patient=None, include_unassigned=False):
    I = Path(input_dir)
    F = Path(figures)
    O = Path(output)
    patient = str(patient or pd.read_csv(I / 'cells.tsv', sep='\t', dtype={'cell_id':str, 'sample':str, 'patient':str}).patient.iloc[0])
    pc = json.loads((I / 'plot_source_contract.json').read_text())
    if not all((sha256(I / n) == h for n, h in pc['sha256'].items())):
        raise ValueError("Integrity check failed: all((sha256(I / n) == h for n, h in pc['sha256'].items()))")
    z = np.load(I / 'source_evidence.npz')
    er = pd.read_csv(I / 'evidence_rows.tsv', sep='\t')
    c = pd.read_csv(I / 'cells.tsv', sep='\t', dtype={'cell_id':str, 'sample':str, 'patient':str})
    order = pd.read_csv(F / 'cell_order.tsv.gz', sep='\t', dtype={'cell_id':str,'sample':str})
    key = pd.read_csv(F / 'clone_key.tsv', sep='\t')
    widths = pd.read_csv(F / 'MAP_widths.tsv', sep='\t')
    if not list(z['cell_ids']) == list(c.cell_id):
        raise ValueError("Integrity check failed: list(z['cell_ids']) == list(c.cell_id)")
    if not list(order.cell_id.unique()).__len__() == len(c):
        raise ValueError('Integrity check failed: list(order.cell_id.unique()).__len__() == len(c)')
    if not set(order.cell_id) == set(c.cell_id):
        raise ValueError('Integrity check failed: set(order.cell_id) == set(c.cell_id)')
    if not z['matrix'].shape == (len(er), len(c)):
        raise ValueError("Integrity check failed: z['matrix'].shape == (len(er), len(c))")
    if not (len(er) > 0 and er.event.nunique() > 0):
        raise ValueError('Integrity check failed: len(er) > 0 and er.event.nunique() > 0')
    if not np.isin(z['matrix'], [0, 1, 2, 3]).all():
        raise ValueError("Integrity check failed: np.isin(z['matrix'], [0, 1, 2, 3]).all()")
    input_cell_count = len(c)
    excluded = order.loc[order['clone'].eq('U') & (not include_unassigned)].copy()
    if not include_unassigned:
        order = order.loc[~order['clone'].eq('U')].copy()
        order['column'] = order.groupby('sample', sort=False).cumcount()
        key = key.loc[~key['clone'].eq('U')].copy()
    if set(order['sample']) != set(c['sample']):
        raise ValueError('No assigned cells remain in a specimen')
    if set(order.cell_id) | set(excluded.cell_id) != set(c.cell_id):
        raise ValueError('Display filter failed complete cell accounting')
    O.mkdir(parents=True, exist_ok=True)
    excluded.to_csv(O/'excluded_unassigned_cells.tsv.gz', sep='\t', index=False)
    require_arial()
    plt.rcParams.update({'font.family': 'Arial', 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'font.size': 6.5, 'svg.hashsalt': '5801-compact-evidence-v1'})
    modalities = {'Accepted call': ('#a6a6a6', 'Accepted'), 'Direct reads': ('#303030', 'Real'), 'PV': ('#d1352b', 'PV'), 'MPV': ('#2f6fb7', 'MPV'), 'LOH': ('#2e9a62', 'LOH'), 'CNV': ('#e08a1e', 'CNV'), 'DEG': ('#7b4ea3', 'DEG'), 'DEI': ('#1f9aa6', 'DEI'), 'LOH classifier': ('#2e9a62', 'LOH')}
    if not set(er.channel) <= set(modalities):
        raise ValueError('Integrity check failed: set(er.channel) <= set(modalities)')
    palette = ['#ffffff', '#acd9a5', '#f06b48', '#f06b48']
    samples = list(order['sample'].drop_duplicates())
    meta = c[['sample', 'stage', 'is_xenograft']].drop_duplicates()
    if not (meta['sample'].is_unique and set(samples) == set(meta['sample'])):
        raise ValueError("Integrity check failed: meta['sample'].is_unique and set(samples) == set(meta['sample'])")
    stages = dict(zip(meta['sample'], meta.stage.str.replace(' MISTRG', '\nMISTRG', regex=False)))
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    labels = [e.replace('RUNX1_LOH', 'RUNX1-LOH').replace('TP53_syn', 'TP53 (syn.)') for e in er.event.unique()]
    left_inches = max(1.1775, max((TextPath((0, 0), t, prop=FontProperties(family='Arial', style='italic', size=6.8)).get_extents().width for t in labels)) / 72 + 0.19)
    clone_lines = []
    for r in key.itertuples():
        txt = r.clone + ' = ' + ('no mutant call' if r.clone == 'WT' else r.mutation_combination.replace('|', ' | '))
        clone_lines.append(textwrap.wrap(txt, 48, break_long_words=False))
    key_row_heights = [max((len(x) for x in clone_lines[j:j + 3])) * 5.2 * 1.1 / 72 + 0.11 for j in range(0, len(clone_lines), 3)]
    key_height = sum(key_row_heights)
    caption_height = 0.58
    heat_height = 2.4225
    height = max(4.75, caption_height + key_height + 0.15 + heat_height + 1.02)
    heat_bottom = caption_height + key_height + 0.15
    fig, axs = plt.subplots(1, len(samples), figsize=(7.5, height), gridspec_kw={'left': left_inches / 7.5, 'right': 0.992, 'top': (heat_bottom + heat_height) / height, 'bottom': heat_bottom / height, 'wspace': 0.035})
    axs = np.atleast_1d(axs)
    R = len(er)
    groups = []
    for event, g in er.groupby('event', sort=False):
        ix = g.index.to_numpy()
        if not np.array_equal(ix, np.arange(ix[0], ix[-1] + 1)):
            raise ValueError('Integrity check failed: np.array_equal(ix, np.arange(ix[0], ix[-1] + 1))')
        groups.append((event, ix[0], ix[-1] + 1))
    byid = {v: i for i, v in enumerate(z['cell_ids'])}
    counts = []
    orders_out = []
    clone_texts = []
    for si, (s, ax) in enumerate(zip(samples, axs)):
        sub = order[order['sample'] == s].sort_values('column')
        if not np.array_equal(sub.column, np.arange(len(sub))):
            raise ValueError('Integrity check failed: np.array_equal(sub.column, np.arange(len(sub)))')
        ix = np.array([byid[x] for x in sub.cell_id])
        mat = z['matrix'][:, ix]
        display = mat.copy()
        display[display == 3] = 2
        if not np.array_equal(display, np.where(mat == 3, 2, mat)):
            raise ValueError('Integrity check failed: np.array_equal(display, np.where(mat == 3, 2, mat))')
        verts = []
        colors = []
        recon = np.zeros_like(display)
        for row, vals in enumerate(display):
            br = np.r_[0, np.flatnonzero(vals[1:] != vals[:-1]) + 1, len(vals)]
            for lo, hi in zip(br[:-1], br[1:]):
                if vals[lo]:
                    verts.append([(lo, row), (hi, row), (hi, row + 1), (lo, row + 1)])
                    colors.append(palette[int(vals[lo])])
                    recon[row, lo:hi] = vals[lo]
        if not np.array_equal(recon, display):
            raise ValueError('Integrity check failed: np.array_equal(recon, display)')
        ax.add_collection(PolyCollection(verts, facecolors=colors, edgecolors=colors, linewidths=0.15, rasterized=True))
        ax.set_xlim(0, len(ix))
        ax.set_ylim(R, -1.8)
        ax.set_xticks([])
        ax.set_yticks([])
        for _, lo, hi in groups:
            if lo:
                ax.axhline(lo, color='#777', lw=0.3)
        for spine in ax.spines.values():
            spine.set_linewidth(0.45)
            spine.set_color('#555')
        start = 0
        for row in key.itertuples():
            n = int(sub.clone.eq(row.clone).sum())
            orig = widths[(widths['sample'] == s) & (widths.clone == row.clone)].iloc[0]
            if not n == orig.MAP_cells:
                raise ValueError('Integrity check failed: n == orig.MAP_cells')
            if n:
                ax.add_patch(Rectangle((start, -1.55), n, 1.25, facecolor=row.color, edgecolor='none'))
                ax.axvline(start, color='#333', lw=0.25)
                if n / len(ix) > 0.045:
                    tx = ax.text(start + n / 2, -0.91, row.clone, fontsize=4.5, ha='center', va='center')
                    clone_texts.append((ax, tx, n / len(ix)))
            if not sub.clone.iloc[start:start + n].eq(row.clone).all():
                raise ValueError('Clone order is not contiguous')
            start += n
        if not start == len(ix):
            raise ValueError('Integrity check failed: start == len(ix)')
        ax.set_title(stages[s] + '\nn=' + f'{len(ix):,}', fontsize=6.5, pad=5)
        counts.append(dict(sample=s, cells=len(ix), original_conflict_entries=int((mat == 3).sum()), drawn_positive_entries=int((display == 2).sum()), all_source_values_retained=True))
        orders_out.extend(sub.to_dict('records'))
    pos = axs[0].get_position()
    strip = fig.add_axes([pos.x0 - 0.0135, pos.y0, 0.009, pos.height])
    strip.set_xlim(0, 1)
    strip.set_ylim(R, -1.8)
    strip.axis('off')
    for i, row in er.iterrows():
        strip.add_patch(Rectangle((0, i), 1, 1, facecolor=modalities[row.channel][0], edgecolor='white', lw=0.12))
    for event, lo, hi in groups:
        label = event.replace('RUNX1_LOH', 'RUNX1-LOH').replace('TP53_syn', 'TP53 (syn.)')
        strip.text(-0.6, (lo + hi) / 2, label, ha='right', va='center', fontsize=6.8, fontstyle='italic', clip_on=False)
        if lo:
            strip.axhline(lo, color='#666', lw=0.3)
    fig.suptitle(patient + ' | Mutation evidence across serial specimens', fontsize=10, y=1 - 0.07 / height)
    handles = [Patch(facecolor=modalities[k][0], label=modalities[k][1]) for k in ['Direct reads', 'PV', 'MPV', 'LOH', 'CNV', 'DEG', 'DEI', 'Accepted call']]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.56, 1 - 0.29 / height), ncol=8, frameon=False, fontsize=6.5, handlelength=1, handleheight=0.75, columnspacing=1.1, handletextpad=0.35)
    fig.legend(handles=[Patch(facecolor=col, edgecolor='#aaa', linewidth=0.3, label=lab) for col, lab in zip(palette[:3], ['No call', 'WT/reference', 'MUT'])], loc='upper center', bbox_to_anchor=(0.56, 1 - 0.55 / height), ncol=3, frameon=False, fontsize=6.5, handlelength=1, handleheight=0.75)
    key_ax = fig.add_axes([0.025, caption_height / height, 0.96, key_height / height])
    key_ax.axis('off')
    for i, r in enumerate(key.itertuples()):
        col = i % 3
        row = i // 3
        x = col / 3
        y = 1 - sum(key_row_heights[:row]) / key_height
        key_ax.add_patch(Rectangle((x, y - 0.04 / key_height), 0.009, 0.04 / key_height, transform=key_ax.transAxes, facecolor=r.color, edgecolor='none'))
        txt = r.clone + ' = ' + ('no mutant call' if r.clone == 'WT' else r.mutation_combination.replace('|', ' | '))
        key_ax.text(x + 0.015, y, '\n'.join(textwrap.wrap(txt, 48, break_long_words=False)), transform=key_ax.transAxes, fontsize=5.2, va='top', fontstyle='italic' if r.clone != 'WT' else 'normal', linespacing=1.1)
    legend = f'{patient}: one patient; {int((~meta.is_xenograft).sum())} clinical and {int(meta.is_xenograft.sum())} matched MISTRG specimens; {len(order):,} displayed cells; {len(excluded):,} U cells excluded. Columns retain the selected cells, ordered by clone then barcode; source rows and clone assignments are unchanged. Small clone-strip labels are omitted when they do not fit; all clones remain in the key. Left colors identify evidence modalities; mutation names label each complete group once. White = no call, green = WT/reference, orange = MUT. Clone widths show maximum reporting assignments, not fractional frequencies. Cell marks have a 0.15 pt matching-color outline for visibility. Direct and imputed evidence are descriptive; no hypothesis test or confidence interval. Imputation and sparse RNA coverage limit genotype certainty.'
    fig.text(0.025, 0.04 / height, '\n'.join(textwrap.wrap(legend, 167)), fontsize=5.2, va='bottom', linespacing=1.1)
    fig.canvas.draw()
    for ax, tx, fraction in clone_texts:
        if tx.get_window_extent(fig.canvas.get_renderer()).width > fraction * ax.get_window_extent().width - 2:
            tx.set_visible(False)
    for t in fig.findobj(match=matplotlib.text.Text):
        if t.get_visible() and t.get_text():
            bb = t.get_window_extent(fig.canvas.get_renderer())
            if not (bb.x0 >= -1 and bb.y0 >= -1 and (bb.x1 <= fig.bbox.width + 1) and (bb.y1 <= fig.bbox.height + 1)):
                raise ValueError(t.get_text())
    O.mkdir(parents=True, exist_ok=True)
    for ext in ['svg', 'pdf', 'png']:
        fig.savefig(O / (patient + '_compact_source_heatmap_thicker.' + ext), dpi=600)
    if not '<image' in (O / (patient + '_compact_source_heatmap_thicker.svg')).read_text():
        raise ValueError("Integrity check failed: '<image' in (O / (patient + '_compact_source_heatmap_thicker.svg')).read_text()")
    plt.close(fig)
    (O / 'figure_legend.txt').write_text(legend + '\nModalities: Real = direct reads; Accepted = consolidated upstream call; PV = nuclear passenger; MPV = mitochondrial passenger; LOH = loss of heterozygosity; CNV = copy-number association; DEG/DEI = gene/isoform expression inference. PV/MPV/LOH/CNV/DEG/DEI colors match the evidence ring. Real and Accepted have neutral charcoal and gray colors, respectively. The original ring has variant-specific colors for direct calls, not one global Real color.\nPer user request, original MUT+WT source entries use MUT display color; raw codes remain unchanged.\n')
    er.assign(modality_color=er.channel.map(lambda x: modalities[x][0])).to_csv(O / 'row_modality_key.tsv', sep='\t', index=False)
    key.to_csv(O / 'clone_key.tsv', sep='\t', index=False)
    pd.DataFrame(orders_out).to_csv(O / 'cell_order.tsv.gz', sep='\t', index=False)
    (O / 'validation.json').write_text(json.dumps(dict(source_heatmap_sha256=sha256(F / 'cell_source_heatmaps.svg') if (F/'cell_source_heatmaps.svg').exists() else None, source_matrix_sha256=sha256(I / 'source_evidence.npz'), source_rows_sha256=sha256(I / 'evidence_rows.tsv'), source_cell_order_sha256=sha256(F / 'cell_order.tsv.gz'), cells=len(order), input_cells=input_cell_count, excluded_unassigned_cells=len(excluded), include_unassigned=include_unassigned, events=int(er.event.nunique()), source_rows=len(er), one_mutation_label_per_event=True, source_rows_and_order_unchanged=True, relative_cell_order_and_clone_assignments_unchanged=True, no_downsampling=True, conflict_display='MUT color; no conflict legend or distinct color', specimens=counts, canvas_inches=[7.5, height], heatmap_outline_pt=0.15, raster_dpi=600), indent=2))
    return json.loads((O / 'validation.json').read_text())

def add_arguments(p):
    p.add_argument('--input', required=True)
    p.add_argument('--figures', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--patient')
    p.add_argument('--include-unassigned', action='store_true', default=False)

def run(a):
    return render(a.input, a.figures, a.output, a.patient, a.include_unassigned)

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    add_arguments(p)
    return run(p.parse_args(argv))
if __name__ == '__main__':
    main()
