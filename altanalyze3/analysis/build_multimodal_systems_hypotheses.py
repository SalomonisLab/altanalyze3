#!/usr/bin/env python3
"""Build context-complete multimodal systems hypotheses from the LungMAP evidence atlas."""

from pathlib import Path
import shutil
import textwrap

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import seaborn as sns


ROOT = Path('/Users/saljh8/Dropbox/LungMAP/Discovery/codex_nature_predictions_20260927')
EXT = ROOT / 'extended'
TAB = EXT / 'tables'
OUT = EXT / 'systems_biology'
FIG = OUT / 'figures'
OTAB = OUT / 'tables'
for directory in [OUT, FIG, OTAB]:
    directory.mkdir(parents=True, exist_ok=True)

sns.set_theme(style='whitegrid', font_scale=.85)
assoc = pd.read_csv(TAB / 'power_aware_associations.tsv.gz', sep='\t', low_memory=False)
assoc = assoc[~assoc.population.str.contains('__vs__', regex=False)].copy()


def save(fig, stem):
    for suffix in ['png', 'pdf']:
        fig.savefig(FIG / f'{stem}.{suffix}', dpi=300 if suffix == 'png' else None,
                    bbox_inches='tight', facecolor='white')
    plt.close(fig)


def significance_label(q, confirmed=False, conflict=False):
    label = '*' if q < .05 else ('†' if q <= .10 else '')
    if confirmed:
        label += '‡'
    if conflict:
        label += '×'
    return label


def volcano(ax, data, labels, title, xlabel, fold_threshold=None):
    z = data.copy()
    z['mlogq'] = -np.log10(z.discovery_fdr.clip(lower=1e-300))
    pass_color = z.discovery_fdr.le(.10)
    if fold_threshold is not None:
        pass_color &= z.discovery_log2fc.abs().ge(fold_threshold)
    colors = np.where(pass_color & z.discovery_log2fc.gt(0), '#cb181d',
             np.where(pass_color & z.discovery_log2fc.lt(0), '#2171b5', '#bdbdbd'))
    ax.scatter(z.discovery_log2fc, z.mlogq, c=colors, s=18, alpha=.75, linewidth=0)
    ax.axvline(0, color='black', lw=.7)
    ax.axhline(-np.log10(.10), color='#555', ls='--', lw=.8)
    if fold_threshold is not None:
        ax.axvline(fold_threshold, color='#777', ls=':', lw=.8)
        ax.axvline(-fold_threshold, color='#777', ls=':', lw=.8)
    chosen = set(labels)
    chosen.update(z.nsmallest(5, 'discovery_fdr').gene)
    for _, r in z[z.gene.isin(chosen)].iterrows():
        ax.annotate(r.gene, (r.discovery_log2fc, r.mlogq), xytext=(3, 3),
                    textcoords='offset points', fontsize=7)
    ax.set_title(title, loc='left', weight='bold')
    ax.set_xlabel(xlabel); ax.set_ylabel('−log10 discovery FDR q')


def evidence_heatmap(ax, data, genes, populations, title):
    columns = pd.MultiIndex.from_product([genes, ['COPD', 'IPF']], names=['gene', 'disease'])
    values = data.pivot_table(index='population', columns=['gene', 'disease'],
                              values='discovery_log2fc', aggfunc='first').reindex(index=populations, columns=columns)
    q = data.pivot_table(index='population', columns=['gene', 'disease'],
                         values='discovery_fdr', aggfunc='first').reindex(index=populations, columns=columns)
    confirmed = data.pivot_table(index='population', columns=['gene', 'disease'],
                                 values='independently_confirmed', aggfunc='first').reindex(index=populations, columns=columns)
    conflict = data.pivot_table(index='population', columns=['gene', 'disease'],
                                values='n_conflicting_cohorts', aggfunc='first').reindex(index=populations, columns=columns)
    lim = np.nanquantile(np.abs(values.to_numpy()), .98)
    sns.heatmap(values, cmap='vlag', center=0, vmin=-lim, vmax=lim, mask=values.isna(),
                linewidths=.25, linecolor='#eeeeee', ax=ax,
                cbar_kws={'label': 'largest-effective-N cohort RNA log2FC\n(disease − control)'})
    ax.set_facecolor('#d9d9d9')
    ax.set_xticklabels([f'{gene}\n{disease}' for gene, disease in columns], rotation=45, ha='right', fontsize=7)
    ax.set_yticklabels(populations, rotation=0, fontsize=7)
    ax.set_xlabel('gene; diseases are paired'); ax.set_ylabel('all tested cell states')
    ax.set_title(title, loc='left', weight='bold')
    for iy, population in enumerate(populations):
        for ix, column in enumerate(columns):
            if pd.isna(values.loc[population, column]):
                continue
            symbol = significance_label(q.loc[population, column],
                                        bool(confirmed.loc[population, column]) if not pd.isna(confirmed.loc[population, column]) else False,
                                        conflict.loc[population, column] > 0)
            if symbol:
                ax.text(ix + .5, iy + .5, symbol, ha='center', va='center', fontsize=6, weight='bold')
        

def draw_experiment(ax, boxes, arrows, title):
    ax.axis('off'); ax.set_title(title, loc='left', weight='bold')
    for x, y, w, h, label, color in boxes:
        patch = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=.02', fc=color, ec='white', transform=ax.transAxes)
        ax.add_patch(patch)
        ax.text(x+w/2, y+h/2, label, ha='center', va='center', transform=ax.transAxes,
                fontsize=8, weight='bold', color='white' if color not in ['#f4d35e', '#e9ecef'] else 'black')
    for x1, y1, x2, y2 in arrows:
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1), xycoords='axes fraction',
                    arrowprops=dict(arrowstyle='->', lw=1.8, color='#444'))


# ---------------------------------------------------------------------------
# Hypothesis 1: IPF fibroblast/epithelial niche.
sp_grn = assoc[(assoc.disease.eq('IPF')) & assoc.population.eq('SPFB') & assoc.modality.eq('grn_tf')].copy()
sp_rna = assoc[(assoc.disease.eq('IPF')) & assoc.population.eq('SPFB') & assoc.modality.eq('rna')].copy()
sp_adt = assoc[(assoc.disease.eq('IPF')) & assoc.population.eq('SPFB') & assoc.modality.eq('adt')].copy()

fig = plt.figure(figsize=(18, 12.5))
gs = fig.add_gridspec(2, 2, hspace=.40, wspace=.28)
ax = fig.add_subplot(gs[0, 0])
volcano(ax, sp_grn, ['TEAD2', 'STAT3', 'NFATC2', 'ETV1', 'SREBF2'],
        'A. SPFB regulatory remodeling includes TEAD2, STAT3, NFATC2 and ETV1',
        'GRN activity difference (IPF − control)')
ax = fig.add_subplot(gs[0, 1])
volcano(ax, sp_rna, ['THY1', 'POSTN', 'CCN2', 'FN1', 'COL1A2', 'RGCC', 'SFRP2', 'PLA2G2A', 'PLTP'],
        'B. The same state expresses matrix, complement and lipid-remodeling genes',
        'RNA log2FC (IPF − control)', fold_threshold=np.log2(1.2))
ax = fig.add_subplot(gs[1, 0])
z = sp_adt.sort_values('discovery_log2fc').copy()
colors = np.where(z.discovery_fdr.le(.10), '#2a9d8f', '#bdbdbd')
ax.barh(z.gene, z.discovery_log2fc, color=colors)
ax.axvline(0, color='black', lw=.8)
for iy, (_, r) in enumerate(z.iterrows()):
    ax.text(r.discovery_log2fc + (.03 if r.discovery_log2fc >= 0 else -.03), iy,
            f'q={r.discovery_fdr:.2g}', va='center',
            ha='left' if r.discovery_log2fc >= 0 else 'right', fontsize=7)
ax.set_xlabel('RNA-imputed ADT effect (IPF − control)')
ax.set_title('C. The full predicted surface panel nominates CD90-high/CD49a+ isolation', loc='left', weight='bold')

sub = gs[1, 1].subgridspec(1, 3, wspace=.55)
palette = {'Control': '#6c757d', 'IPF': '#d95f59'}
donor_comp = pd.read_csv(ROOT / 'tables/xenium_donor_composition.tsv', sep='\t')
dc = donor_comp[donor_comp.final_CT.eq('Activated Fibrotic FBs')].copy()
dc['status'] = np.where(dc.group.eq('Unaffected'), 'Control', 'IPF')
ax = fig.add_subplot(sub[0, 0])
sns.boxplot(data=dc, x='status', y='logit_fraction', palette=palette, hue='status', legend=False,
            whis=1.5, showfliers=False, ax=ax)
sns.stripplot(data=dc, x='status', y='logit_fraction', palette=palette, hue='status', legend=False,
              jitter=.12, size=3.4, ax=ax)
ax.set_title('D1. Fibrotic FBs expand', weight='bold', fontsize=9); ax.set_xlabel(''); ax.set_ylabel('donor logit fraction')

xdir = Path('/Users/saljh8/Dropbox/LungMAP/Xenium-IPF/inputs')
meta = pd.read_csv(xdir / 'pseudobulk_metadata_final_CT_by_patient.tsv', sep='\t')
expr = pd.read_csv(xdir / 'pseudobulk_avgExpr_final_CT_by_patient-clean.txt', sep='\t', index_col=0).T
meta = meta[meta.pseudobulk_id.astype(str).isin(expr.index)].copy()
expr = expr.loc[meta.pseudobulk_id.astype(str)].copy(); expr.index = meta.index
ix = meta.final_CT.eq('Activated Fibrotic FBs')
ex = pd.DataFrame({'FN1': expr.loc[ix, 'FN1'].astype(float),
                   'status': np.where(meta.loc[ix, 'sample_affect'].eq('Unaffected'), 'Control', 'IPF')})
ax = fig.add_subplot(sub[0, 1])
sns.boxplot(data=ex, x='status', y='FN1', palette=palette, hue='status', legend=False,
            whis=1.5, showfliers=False, ax=ax)
sns.stripplot(data=ex, x='status', y='FN1', palette=palette, hue='status', legend=False,
              jitter=.12, size=3.4, ax=ax)
ax.set_title('D2. Fibrotic FB FN1 rises', weight='bold', fontsize=9); ax.set_xlabel(''); ax.set_ylabel('log2(CP10K+1)')

spatial = pd.read_csv(ROOT / 'tables/xenium_spatial_pairs_donor_region.tsv', sep='\t')
q = spatial[(spatial.cell_a.eq('KRT5-_KRT17+')) & (spatial.cell_b.eq('Activated_Fibrotic_FBs'))].copy()
q['status'] = np.where(q.group.eq('Unaffected'), 'Control', 'IPF')
q = q.groupby(['donor_id', 'status'], as_index=False).squidpy_z.median()
ax = fig.add_subplot(sub[0, 2])
sns.boxplot(data=q, x='status', y='squidpy_z', palette=palette, hue='status', legend=False,
            whis=1.5, showfliers=False, ax=ax)
sns.stripplot(data=q, x='status', y='squidpy_z', palette=palette, hue='status', legend=False,
              jitter=.12, size=3.4, ax=ax)
ax.set_title('D3. KRT–FB proximity rises', weight='bold', fontsize=9); ax.set_xlabel(''); ax.set_ylabel('median donor Squidpy z')

fig.suptitle('A TEAD2/STAT3-active, CD90-high subpleural-fibroblast program converges with an FN1-rich KRT epithelial niche in IPF',
             fontsize=16, weight='bold', y=.985)
caption = (
    'Systems hypothesis. Panels A–C show the complete SPFB feature universe for the indicated modality; gray points/bars are retained nonsignificant tests and labels include prespecified mechanistic genes plus the five lowest-q features. Discovery used the cohort with the largest effective donor N. The exploratory rule was q≤0.10; RNA additionally required |log2FC|≥log2(1.2), chosen to balance unequal cohort power with a minimum biological effect. Exact q values remain continuous. '
    'GRN activity and ADT are RNA-derived predictions, whereas RNA is measured. Panels D1–D3 are independent Xenium donor-level controls: boxes show median and IQR, whiskers 1.5×IQR, and every point is a biological donor pseudobulk. Fibrotic-fibroblast abundance and FN1 expression increase, and KRT–fibroblast proximity differs by two-sided Mann–Whitney U P=0.0022 after donor aggregation. '
    'The data nominate, but do not prove, a lineage relationship between CellRef SPFB and Xenium activated fibrotic fibroblasts. Decisive test: sort CD45−EPCAM−CD31−CD90-high/CD49a+/PDGFRα+ fibroblasts, perturb TEAD2/STAT3 or FN1, and quantify matrix assembly and induction of KRT17/ITGB6 in donor-matched alveolar organoids.'
)
fig.text(.02, .012, textwrap.fill(caption, 230), ha='left', va='bottom', fontsize=8.5)
fig.subplots_adjust(left=.07, right=.98, top=.92, bottom=.13)
save(fig, 'Systems1_IPF_fibroblast_epithelial_niche')


# ---------------------------------------------------------------------------
# Hypothesis 2: shared macrophage state switch.
switch_genes = ['CD163', 'PLIN2', 'VSIG4', 'SGK1', 'USP53', 'FN1', 'RGCC']
mac_context = assoc[(assoc.modality.eq('rna')) & assoc.gene.isin(switch_genes)].copy()
preferred = ['AT1', 'AT2', 'Aberrant basal', 'AdvFB', 'SPFB', 'PVEC',
             'AM', 'AM-lipid', 'AM-prolif', 'MT+ AM', 'IM', 'tMDM', 'iMON', 'pMON',
             'Neutrophil', 'cDC1', 'maDC', 'pDC', 'B', 'CD4 T', 'CD8 T', 'NK']
present = set(mac_context.population)
pops = [p for p in preferred if p in present] + sorted(present - set(preferred))
fig = plt.figure(figsize=(18, 20))
gs = fig.add_gridspec(2, 2, height_ratios=[2.5, .8], hspace=.28, wspace=.28)
ax = fig.add_subplot(gs[0, :])
evidence_heatmap(ax, mac_context, switch_genes, pops,
                 'A. All tested cell states show where the resident-loss and remodeling genes are—and are not—specific')

ax = fig.add_subplot(gs[1, 0])
tf_genes = ['CEBPA', 'FOXA2', 'SREBF2', 'KLF5', 'TEAD2', 'JDP2']
tf = assoc[(assoc.modality.eq('grn_tf')) & assoc.gene.isin(tf_genes) &
           assoc.population.isin(['AM', 'AM-lipid', 'AM-prolif', 'MT+ AM', 'IM', 'tMDM'])].copy()
tf['row'] = tf.population + ' | ' + tf.disease
tm = tf.pivot_table(index='row', columns='gene', values='discovery_log2fc', aggfunc='first')
tq = tf.pivot_table(index='row', columns='gene', values='discovery_fdr', aggfunc='first')
lim = np.nanquantile(np.abs(tm.to_numpy()), .98)
sns.heatmap(tm, cmap='vlag', center=0, vmin=-lim, vmax=lim, mask=tm.isna(), linewidths=.4,
            cbar_kws={'label': 'GRN activity difference'}, ax=ax)
ax.set_facecolor('#d9d9d9'); ax.set_xlabel('RNA-derived TF activity'); ax.set_ylabel('macrophage state | disease')
ax.set_title('B. Lipid/homeostatic TF loss is strongest in IPF MT+ alveolar macrophages', loc='left', weight='bold')
for iy, row in enumerate(tm.index):
    for ix, gene in enumerate(tm.columns):
        if not pd.isna(tm.loc[row, gene]) and tq.loc[row, gene] <= .10:
            ax.text(ix+.5, iy+.5, '*' if tq.loc[row, gene] < .05 else '†', ha='center', va='center', fontsize=7)

ax = fig.add_subplot(gs[1, 1])
draw_experiment(ax,
    [(0.02, .63, .25, .20, 'Resident AM/IM\nCD163 PLIN2 VSIG4', '#457b9d'),
     (0.38, .63, .25, .20, 'Disease macrophage\nSGK1 FN1 RGCC', '#d95f59'),
     (0.72, .63, .25, .20, 'Matrix/integrin\nepithelial response', '#6a994e'),
     (0.20, .18, .25, .18, 'Perturb\nSGK1 / SYK / ROCK', '#6c757d'),
     (0.57, .18, .25, .18, 'Rescue criteria\nresident lipid handling', '#f4d35e')],
    [(0.27, .73, .38, .73), (0.63, .73, .72, .73), (0.45, .62, .34, .36), (0.45, .27, .57, .27)],
    'C. A donor-matched perturbation distinguishes state conversion from macrophage replacement')
fig.suptitle('IPF and COPD share loss of macrophage lipid/scavenger identity and SGK1 remodeling, with FN1 dominance in IPF',
             fontsize=16, weight='bold', y=.985)
caption = (
    'Systems hypothesis. Panel A retains every available direct disease-versus-control RNA test for seven prespecified state-switch genes across all available epithelial, stromal, vascular, myeloid and lymphoid cell states; gray is unavailable, not null. Red/blue are largest-effective-N cohort log2FC. * denotes discovery q<0.05, † denotes 0.05≤q≤0.10, ‡ adds an independent concordant nominal two-sided P<0.05 without a significant opposing cohort, and × denotes directional conflict. RNA discovery also required |log2FC|≥log2(1.2). '
    'Panel B shows all available macrophage-state GRN tests for six prespecified homeostatic/regulatory TFs; GRN scores are inferred from RNA and are not TF-protein measurements. The shared RNA component is strongest for AM CD163/PLIN2 loss and IM SGK1/USP53 gain, while FN1 and lipid-TF changes are more IPF-weighted; therefore the hypothesis is a partially shared state axis, not identical macrophage biology in both diseases. '
    'Decisive test: lineage-resolved primary AM/IM cultures from both diseases with SGK1, SYK and ROCK perturbation; require restoration of CD163/PLIN2/VSIG4 and lipid handling without suppressing viability, plus reduced FN1 matrix and epithelial ITGB6/KRT17 induction in coculture.'
)
fig.text(.02, .012, textwrap.fill(caption, 230), ha='left', va='bottom', fontsize=8.5)
fig.subplots_adjust(left=.08, right=.98, top=.94, bottom=.09)
save(fig, 'Systems2_shared_macrophage_state_switch')


# ---------------------------------------------------------------------------
# Hypothesis 3: COPD PVEC regulatory failure and incoming signaling.
pvec = assoc[(assoc.disease.eq('COPD')) & assoc.population.eq('PVEC')].copy()
pvec_tf = pvec[pvec.modality.eq('grn_tf')].copy()
selected_tf = ['FOXA2', 'SREBF2', 'KLF5', 'TEAD2', 'TEAD3', 'NFIC', 'NFIX', 'ELF3', 'ATF7']
fig = plt.figure(figsize=(18, 11.5))
gs = fig.add_gridspec(2, 2, hspace=.42, wspace=.30)
ax = fig.add_subplot(gs[0, 0])
volcano(ax, pvec_tf, selected_tf,
        'A. COPD PVEC loses a coherent FOXA2/SREBF2/KLF5/TEAD regulatory program',
        'GRN activity difference (COPD − control)')

ax = fig.add_subplot(gs[0, 1])
rna_tf = pvec[(pvec.modality.eq('rna')) & pvec.gene.isin(selected_tf)].set_index('gene').reindex(selected_tf)
available = rna_tf.discovery_log2fc.notna()
bar_values = rna_tf.discovery_log2fc.fillna(0)
ax.barh(rna_tf.index, bar_values, color=np.where(rna_tf.discovery_fdr.le(.10), '#2a9d8f', '#bdbdbd'))
ax.axvline(0, color='black', lw=.8)
for iy, (_, r) in enumerate(rna_tf.iterrows()):
    label = f'q={r.discovery_fdr:.2g}' if pd.notna(r.discovery_fdr) else 'not available'
    ax.text(.98, iy, label, transform=ax.get_yaxis_transform(),
            va='center', ha='right', fontsize=7)
ax.set_xlabel('measured TF RNA log2FC (COPD − control)')
ax.set_title('B. TF RNA does not reproduce the inferred activity loss', loc='left', weight='bold')

ax = fig.add_subplot(gs[1, 0])
comm_all = pvec[pvec.modality.eq('fastcomm')].copy()
comm_all['mlogq'] = -np.log10(comm_all.discovery_fdr.clip(lower=1e-300))
comm_confirmed = comm_all[comm_all.independently_confirmed].copy()
ax.scatter(comm_all.discovery_log2fc, comm_all.mlogq, s=7, c='#c7c7c7', alpha=.35, linewidth=0,
           label=f'all tested edges (n={len(comm_all):,})')
ax.scatter(comm_confirmed.discovery_log2fc, comm_confirmed.mlogq, s=26,
           c=np.where(comm_confirmed.discovery_log2fc.gt(0), '#d95f59', '#457b9d'),
           alpha=.9, linewidth=.3, edgecolor='white', label=f'independently confirmed (n={len(comm_confirmed)})')
ax.axvline(0, color='black', lw=.8)
ax.axhline(-np.log10(.10), color='#555', ls='--', lw=.8)
label_edges = comm_confirmed.nsmallest(8, 'discovery_fdr')
for _, r in label_edges.iterrows():
    label = r.gene.replace(' cell|', '|').replace('Aberrant basaloid', 'Aberrant epi')
    ax.annotate(label, (r.discovery_log2fc, r.mlogq), xytext=(3, 2), textcoords='offset points', fontsize=6)
ax.set_xlabel('RNA-derived communication-score difference (COPD − control)')
ax.set_ylabel('−log10 discovery FDR q')
ax.set_title('C. Confirmed mesothelial/epithelial inputs are highlighted within the complete incoming-edge universe', loc='left', weight='bold')
ax.legend(frameon=False, fontsize=7)

ax = fig.add_subplot(gs[1, 1])
draw_experiment(ax,
    [(0.02, .65, .25, .20, 'Mesothelial ligands\nPLA2G2A / ANXA2 / NAMPT', '#6a994e'),
     (0.38, .65, .25, .20, 'COPD PVEC\nαv/β1 / PLAT / PTPRB', '#d95f59'),
     (0.72, .65, .25, .20, 'Barrier + lipid\ntransport phenotype', '#457b9d'),
     (0.18, .18, .28, .18, 'Perturb ligand/receptor\n± FOXA2/SREBF2 rescue', '#6c757d'),
     (0.58, .18, .28, .18, 'Readouts\nTEER, permeability, lipid flux', '#f4d35e')],
    [(0.27, .75, .38, .75), (0.63, .75, .72, .75), (0.50, .64, .34, .36), (0.46, .27, .58, .27)],
    'D. An endothelial-chip experiment can separate incoming signaling from intrinsic TF failure')

fig.suptitle('COPD pulmonary-venous endothelium shows inferred identity-program failure despite preserved TF RNA and altered incoming signals',
             fontsize=16, weight='bold', y=.985)
caption = (
    'Systems hypothesis. Panel A retains the complete COPD PVEC GRN feature universe; colored points pass exploratory discovery q≤0.10 and labels include the prespecified nine-factor identity module plus the five lowest-q factors. Panel B shows measured RNA for the same TFs, including nonsignificant and opposite estimates, demonstrating that RNA abundance is not evidence of the inferred activity loss. '
    'Panel C displays every tested incoming fastcomm edge as gray background and highlights only edges that passed exploratory discovery and independent concordant confirmation without directional conflict; positive and negative edges are both retained. The eight lowest-q confirmed edges are labeled to prevent text overlap. fastcomm and GRN are RNA-derived, reuse expression information, and are not orthogonal evidence of ligand binding or TF activity. COPD discovery cohorts and donor counts vary by modality/state and are retained in the evidence workbook. '
    'The model is innovative but less mature than the IPF niche because it lacks spatial/protein validation. Decisive test: donor-matched pulmonary-venous endothelial chips exposed to mesothelial conditioned medium, with PLA2G2A/ITGAV/ITGB1 blockade and FOXA2 or SREBF2 rescue; quantify barrier resistance, permeability, lipid transport, APLNR/LAMA4, and chromatin accessibility. Failure of ligand blockade and TF rescue to alter these endpoints would reject the proposed coupling.'
)
fig.text(.02, .012, textwrap.fill(caption, 230), ha='left', va='bottom', fontsize=8.5)
fig.subplots_adjust(left=.16, right=.98, top=.92, bottom=.14)
save(fig, 'Systems3_COPD_PVEC_regulatory_communication')

# SYS-4 already has a context-complete figure with separate measured and imputed modalities.
for suffix in ['png', 'pdf']:
    source = ROOT / 'figures' / f'Figure2_infection_surfactant_lipid_checkpoint.{suffix}'
    if source.exists():
        shutil.copy2(source, FIG / f'Systems4_acute_injury_surfactant_lipid_checkpoint.{suffix}')


# ---------------------------------------------------------------------------
# Evidence tables and explicit falsification criteria.
hypotheses = pd.DataFrame([
    dict(hypothesis_id='SYS-1', claim='A TEAD2/STAT3-active, CD90-high SPFB program converges with an FN1-rich KRT epithelial niche in IPF.',
         compartments='SPFB/fibrotic fibroblast; KRT5−/KRT17+ epithelium',
         key_statistics='SPFB TEAD2 GRN Δ=0.92, q=8.8e-5, independently confirmed; THY1 RNA log2FC=1.81, q=6.9e-10, independently confirmed; Xenium fibrotic-FB abundance g=3.19, q=5.4e-5; fibrotic-FB FN1 g=1.77, q=0.0066; donor-level KRT–FB proximity P=0.0022.',
         direct_evidence='SPFB RNA matrix program; Xenium donor composition, FN1 expression and KRT–fibroblast proximity',
         inferred_evidence='SPFB GRN activity; RNA-imputed CD90/CD49a surface phenotype; ligand direction',
         negative_context='SPFB and Xenium activated fibrotic fibroblasts are not proven to be the same lineage; KRT–SPP1-macrophage proximity P=0.17',
         decisive_experiment='Sort CD90-high/CD49a+/PDGFRα+ fibroblasts; TEAD2/STAT3/FN1 perturbation in donor-matched alveolar coculture',
         rejection_criterion='No change in FN1 matrix or epithelial KRT17/ITGB6 after target-engaged perturbation.'),
    dict(hypothesis_id='SYS-2', claim='IPF and COPD share macrophage resident lipid/scavenger loss and SGK1 remodeling, with an IPF-weighted FN1 arm.',
         compartments='AM, IM, MT+ AM; epithelium and matrix',
         key_statistics='COPD AM CD163 log2FC=−1.02, q=5.8e-4 and PLIN2=−0.55, q=0.017; IPF IM SGK1=1.05, q=1.2e-7 and COPD IM SGK1=0.74, q=0.0015; IPF IM FN1=2.63, q=3.8e-6. Listed associations have independent confirmation.',
         direct_evidence='Cross-disease AM CD163/PLIN2 loss and IM SGK1/USP53 gain; IPF macrophage FN1 increase',
         inferred_evidence='IPF macrophage lipid/homeostatic GRN loss; communication and drug reversal',
         negative_context='FN1 and TF changes are more IPF-weighted; this is a partially shared axis, not identical disease biology',
         decisive_experiment='Lineage-resolved AM/IM SGK1, SYK and ROCK perturbation with lipid-handling and epithelial-coculture readouts',
         rejection_criterion='No restoration of resident markers/lipid handling or reduction of FN1/epithelial ITGB6 at noncytotoxic exposure.'),
    dict(hypothesis_id='SYS-3', claim='COPD PVEC has regulatory identity failure coupled to altered mesothelial/epithelial input.',
         compartments='Mesothelium/aberrant epithelium; pulmonary venous endothelium',
         key_statistics='PVEC KLF5 GRN Δ=−0.57, q=0.0043; FOXA2=−0.38, q=0.0098; SREBF2=−0.68, q=0.0188; corresponding TF RNA q values are 0.76–1.0 or unavailable; 13 incoming fastcomm edges independently confirm.',
         direct_evidence='PVEC disease RNA endpoints including APLNR/LAMA4/IFI6',
         inferred_evidence='FOXA2/SREBF2/KLF5/TEAD GRN loss and fastcomm ligand–receptor edges',
         negative_context='Corresponding TF RNA is largely nonsignificant; no spatial or protein validation',
         decisive_experiment='PVEC chip with mesothelial conditioned medium, ligand/receptor blockade and TF rescue',
         rejection_criterion='Neither incoming-signal blockade nor TF rescue changes barrier, lipid flux or PVEC identity.'),
    dict(hypothesis_id='SYS-4', claim='Acute-injury AT2 intermediates lose an ETV5/XBP1 surfactant checkpoint with predicted PUFA-phospholipid depletion.',
         compartments='AT2→AT1 intermediate; surfactant/lipid metabolism',
         key_statistics='COVID/pneumonia ETV5 GRN Δ=−1.60/−2.31 and XBP1=−0.46/−0.44; SFTPC RNA log2FC=−4.66/−8.43; predicted PC(20:4/22:6)=−0.47/−0.58.',
         direct_evidence='Measured SFTPA1/SFTPA2/SFTPC/MFSD2A RNA loss in COVID-19 and pneumonia state contrasts',
         inferred_evidence='ETV5/XBP1/CEBPA/MLXIPL GRN loss; RNA-imputed ADT and lipid changes',
         negative_context='State-versus-state contrast is not longitudinal transition; no direct lipid or protein assay',
         decisive_experiment='ETV5 CRISPRa or XBP1s rescue in infected alveolospheres with LC–MS/MS and surfactant biophysics',
         rejection_criterion='Target engagement fails to restore surfactant RNA/function or the three prespecified lipids.')
])
hypotheses.to_csv(OTAB / 'systems_hypotheses.tsv', sep='\t', index=False)
sp_grn.to_csv(OTAB / 'SYS1_SPFB_complete_GRN.tsv', sep='\t', index=False)
sp_rna.to_csv(OTAB / 'SYS1_SPFB_complete_RNA.tsv.gz', sep='\t', index=False)
sp_adt.to_csv(OTAB / 'SYS1_SPFB_complete_ADT.tsv', sep='\t', index=False)
mac_context.to_csv(OTAB / 'SYS2_complete_cell_context.tsv.gz', sep='\t', index=False)
pvec.to_csv(OTAB / 'SYS3_PVEC_complete_evidence.tsv.gz', sep='\t', index=False)

with pd.ExcelWriter(OTAB / 'multimodal_systems_hypotheses.xlsx', engine='openpyxl') as writer:
    hypotheses.to_excel(writer, sheet_name='Hypotheses', index=False)
    sp_grn.to_excel(writer, sheet_name='SYS1 SPFB GRN', index=False)
    sp_rna.to_excel(writer, sheet_name='SYS1 SPFB RNA', index=False)
    sp_adt.to_excel(writer, sheet_name='SYS1 SPFB ADT', index=False)
    mac_context.to_excel(writer, sheet_name='SYS2 macrophage context', index=False)
    pvec.to_excel(writer, sheet_name='SYS3 PVEC evidence', index=False)

report = f'''# Multimodal systems-biology hypotheses

This expansion prioritizes mechanisms that connect at least two biological compartments and at least two modalities. It does not assign a subjective score. Each hypothesis is accompanied by complete reference distributions, explicit negative evidence, and a rejection experiment.

## Principal hypotheses

1. **SYS-1 — fibroblast/epithelial matrix niche.** A TEAD2/STAT3-active, CD90-high subpleural-fibroblast program converges with the independently measured FN1-rich activated-fibroblast/KRT epithelial niche in IPF. This is the most experimentally mature expansion, but SPFB-to-Xenium-fibroblast lineage equivalence remains unproven.
2. **SYS-2 — shared macrophage state axis.** IPF and COPD share loss of macrophage resident lipid/scavenger identity and gain of SGK1-centered remodeling, while the FN1 arm and lipid-TF loss are more IPF-weighted. This predicts a shared axis with disease-specific endpoints rather than a single universal macrophage state.
3. **SYS-3 — COPD venous endothelial control failure.** PVEC shows coherent RNA-inferred FOXA2/SREBF2/KLF5/TEAD activity loss without corresponding TF RNA loss, alongside altered predicted input from mesothelial and aberrant epithelial compartments. It is innovative but requires orthogonal spatial/protein validation.
4. **SYS-4 — acute-injury surfactant/lipid checkpoint.** ETV5/XBP1/CEBPA/MLXIPL regulatory loss accompanies measured surfactant-gene loss and predicted PUFA-phospholipid depletion in COVID-19 and pneumonia intermediates. Direct lipidomics is the gating validation.

## Statistical framing

The discovery cutoff `q≤0.10` is an exploratory rule chosen for unequal donor power and combined with an RNA effect floor of `|log2FC|≥log2(1.2)`. It is not treated as a biological truth. Exact continuous effects and q values are retained; `q<0.05` and `0.05≤q≤0.10` are displayed separately. Independent confirmation requires concordant direction and nominal two-sided `P<0.05` in another cohort without a significant opposite cohort. The figures retain nonsignificant, opposite, conflicting and unavailable tests.

## Files

- `tables/multimodal_systems_hypotheses.xlsx` — claims, complete evidence, negative context and rejection criteria.
- `figures/Systems1_IPF_fibroblast_epithelial_niche.*`
- `figures/Systems2_shared_macrophage_state_switch.*`
- `figures/Systems3_COPD_PVEC_regulatory_communication.*`
- `figures/Systems4_acute_injury_surfactant_lipid_checkpoint.*`
'''
(OUT / 'REPORT.md').write_text(report)
(OUT / 'scripts').mkdir(exist_ok=True)
shutil.copy2(Path(__file__), OUT / 'scripts' / Path(__file__).name)
print({'hypotheses': len(hypotheses), 'sys1_grn': len(sp_grn), 'sys1_rna': len(sp_rna),
       'sys2_context': len(mac_context), 'sys3_evidence': len(pvec)})
