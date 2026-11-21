"""Companion clone/timepoint panels on saved coordinates. Native AltAnalyze2 remains the primary renderer."""
from pathlib import Path
import argparse, json, hashlib, textwrap
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from altanalyze3.components.bam.clonal_inference.detection_plotting import require_arial

def render(inputs, output, config=None):
    I = Path(inputs)
    O = Path(output)
    O.mkdir(parents=True, exist_ok=True)
    config = json.loads(Path(config).read_text()) if config else {}
    patient = config.get('patient', 'Patient')
    v = json.loads((I / 'mapping_validation.json').read_text())
    if not all((hashlib.sha256((I / f).read_bytes()).hexdigest() == s for f, s in v['render_inputs_sha256'].items())):
        raise ValueError("Integrity check failed: all((hashlib.sha256((I / f).read_bytes()).hexdigest() == s for f, s in v['render_inputs_sha256'].items()))")
    d = pd.read_csv(I / 'mapped_umap_cells.tsv.gz', sep='\t')
    k = pd.read_csv(I / 'clone_key.tsv', sep='\t')
    if not (len(d) == v['UMAP_cells'] and d.cell_id.notna().sum() == v['matched_cells'] and d.umap_id.is_unique):
        raise ValueError("Integrity check failed: len(d) == v['UMAP_cells'] and d.cell_id.notna().sum() == v['matched_cells'] and d.umap_id.is_unique")
    require_arial()
    plt.rcParams.update({'font.family': 'Arial', 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'font.size': 8, 'svg.hashsalt': '5801-illumina-umap-v1'})
    colors = k.set_index('clone').color.to_dict()
    colors.update({'Unmatched': '#d5d5d5', 'Aged healthy': '#555555'})
    limits = []
    for col in ['UMAP1', 'UMAP2']:
        lo, hi = (d[col].min(), d[col].max())
        pad = (hi - lo) * 0.06
        limits.append((lo - pad, hi + pad))
    matched = d[d.cell_id.notna()]
    rng = np.random.default_rng(222)
    caption = config.get('caption', f'Saved coordinates unchanged. {len(d):,} cells; {len(matched):,} exact library/barcode clone matches. Unmatched cells are gray, controls dark gray, and are never assigned WT. WT means no accepted mutation in the original clone analysis. Descriptive overlay; no UMAP refitting or hypothesis test. Matching and RNA inference limit interpretation.')
    records = []

    def base(ax):
        ax.scatter(d.UMAP1, d.UMAP2, s=1.25, c=colors['Unmatched'], linewidths=0, rasterized=True, zorder=1)
        h = d[d.plot_clone == 'Aged healthy']
        ax.scatter(h.UMAP1, h.UMAP2, s=2.1, c=colors['Aged healthy'], linewidths=0, rasterized=True, zorder=2)
        ax.set_xlim(*limits[0])
        ax.set_ylim(*limits[1])
        ax.set_aspect('equal', adjustable='box')
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ['top', 'right']:
            ax.spines[side].set_visible(False)
        for side in ['left', 'bottom']:
            ax.spines[side].set_linewidth(0.65)
        ax.set_xlabel('UMAP-1', fontsize=8)
        ax.set_ylabel('UMAP-2', fontsize=8)

    def overlay(ax, x, size=3):
        x = x.iloc[rng.permutation(len(x))]
        ax.scatter(x.UMAP1, x.UMAP2, c=x.clone.map(colors), s=size, linewidths=0, alpha=1, rasterized=True, zorder=3)

    def footer(fig, top=0.23):
        ax = fig.add_axes([0.03, 0.095, 0.94, top - 0.095])
        ax.axis('off')
        for i, r in enumerate(k.itertuples()):
            x = i % 2 * 0.5
            y = 1 - i // 2 * 0.235
            label = r.clone + ' = ' + ('no mutant call' if r.clone == 'WT' else r.mutation_combination.replace('|', ' | '))
            ax.plot([x + 0.006], [y - 0.025], marker='o', markersize=4, color=r.color, transform=ax.transAxes)
            ax.text(x + 0.02, y, '\n'.join(textwrap.wrap(label, 68, break_long_words=False)), fontsize=6.5, va='top', fontstyle='italic' if r.clone != 'WT' else 'normal', transform=ax.transAxes, linespacing=1.05)
        fig.text(0.03, 0.012, '\n'.join(textwrap.wrap(caption, int(fig.get_figwidth() * 20))), fontsize=5.6, va='bottom', linespacing=1.1)

    def export(fig, name):
        fig.canvas.draw()
        for t in fig.findobj(match=matplotlib.text.Text):
            if t.get_visible() and t.get_text():
                bb = t.get_window_extent(fig.canvas.get_renderer())
                if not (bb.x0 >= -2 and bb.y0 >= -2 and (bb.x1 <= fig.bbox.width + 2) and (bb.y1 <= fig.bbox.height + 2)):
                    raise ValueError((name, t.get_text(), bb.bounds))
        for ext in ['svg', 'pdf', 'png']:
            fig.savefig(O / (name + '.' + ext), dpi=600)
        records.append({'figure': name, 'background_cells': len(d), 'coordinates_unchanged': True, 'point_layer_raster_dpi': 600, 'labels_editable': True})
        plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.5, 7))
    fig.subplots_adjust(left=0.08, right=0.96, top=0.9, bottom=0.29)
    base(ax)
    overlay(ax, matched, 3)
    ax.set_title(config.get('title', patient + ' | Clones on saved UMAP') + f'\n{len(matched):,} matched / {len(d):,} coordinate cells', fontsize=11, pad=10)
    ax.legend(handles=[Line2D([], [], marker='o', ls='', color=colors['Unmatched'], label='Unmatched patient cells'), Line2D([], [], marker='o', ls='', color=colors['Aged healthy'], label='Aged healthy controls')], loc='upper left', frameon=False, fontsize=7, markerscale=0.7)
    footer(fig)
    export(fig, config.get('prefix', patient + '_') + 'UMAP_clones')
    stages = config.get('stages', matched.umap_timepoint.drop_duplicates().tolist())
    fig, axs = plt.subplots(1, len(stages), squeeze=False, figsize=(11, 5.7))
    fig.subplots_adjust(left=0.045, right=0.99, top=0.85, bottom=0.35, wspace=0.18)
    fig.suptitle(patient + ' | Clone assignments across clinical timepoints', fontsize=12, y=0.97)
    for ax, st in zip(axs.flat, stages):
        base(ax)
        x = matched[matched.umap_timepoint == st]
        overlay(ax, x, 3.7)
        ax.set_title(st.replace(' / ', ' /\n') + '\n' + f'{len(x):,} matched cells', fontsize=9, pad=8)
    footer(fig, 0.28)
    export(fig, config.get('prefix', patient + '_') + 'UMAP_clones_by_timepoint')
    fig, axs = plt.subplots(int(np.ceil(len(k) / 4)), 4, squeeze=False, figsize=(11, 9))
    fig.subplots_adjust(left=0.05, right=0.99, top=0.91, bottom=0.27, wspace=0.18, hspace=0.28)
    fig.suptitle(patient + ' | Individual clones on the same saved UMAP', fontsize=12, y=0.97)
    for ax, r in zip(axs.flat, k.itertuples()):
        base(ax)
        x = matched[matched.clone == r.clone]
        overlay(ax, x, 4)
        ax.set_title(r.clone + '  |  ' + f'{len(x):,} matched cells', fontsize=10, color=r.color if r.clone != 'WT' else '#795c1a', pad=6)
    footer(fig, 0.22)
    export(fig, config.get('prefix', patient + '_') + 'UMAP_individual_clones')
    (O / 'figure_legend.txt').write_text(caption + '\nPoints rasterized at 600 dpi; labels and axes remain editable vectors. Marker draw order is shuffled using seed 222 to avoid systematic clone-order occlusion. All panels retain all saved coordinate cells; timepoint/clone panels highlight only the stated subset.\n')
    (O / 'plot_validation.json').write_text(json.dumps(records, indent=2))
    return records

def add_arguments(p):
    p.add_argument('--inputs', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--config')

def run(a):
    return render(a.inputs, a.output, a.config)

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    add_arguments(p)
    return run(p.parse_args(argv))
if __name__ == '__main__':
    main()
