#!/usr/bin/env python3
"""Make remaining figures self-contained and remove positive-only drug displays."""

from pathlib import Path
import shutil
import textwrap

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = Path('/Users/saljh8/Dropbox/LungMAP/Discovery/codex_nature_predictions_20260927')
EXT = ROOT / 'extended'
TAB = EXT / 'tables'
FIG = EXT / 'figures'
BASEFIG = ROOT / 'figures'
sns.set_theme(style='whitegrid', font_scale=.86)


def save(fig, directory, stem):
    for suffix in ['png', 'pdf']:
        fig.savefig(directory / f'{stem}.{suffix}', dpi=300 if suffix == 'png' else None,
                    bbox_inches='tight', facecolor='white')
    plt.close(fig)


# Figure 1: same evidence, now with the exact biological question and full caption.
ipf = pd.read_csv(ROOT / 'tables/ipf_multimodal_evidence.tsv', sep='\t')
fig = plt.figure(figsize=(20, 13.5))
gs = fig.add_gridspec(2, 6, hspace=.48, wspace=.72)
ax = fig.add_subplot(gs[0, :3]); ax.axis('off')
ax.set_title('A. The measurements motivate a testable multicellular mechanism; arrows remain untested', loc='left', weight='bold')
nodes = [(.08, .62, 'SPP1 macrophage\nFN1 / SPP1', '#457b9d'),
         (.08, .22, 'Fibrotic fibroblast\nFN1 matrix', '#6a994e'),
         (.53, .43, 'αvβ6 / β8 / β1\naberrant epithelium', '#e9c46a'),
         (.86, .43, 'TEAD2 / ZNF322\nadhesion program', '#d95f59')]
for x, y, label, color in nodes:
    ax.text(x, y, label, ha='center', va='center', transform=ax.transAxes,
            bbox=dict(boxstyle='round,pad=.65', fc=color, ec='white'),
            color='black' if color == '#e9c46a' else 'white', weight='bold')
for y in [.62, .22]:
    ax.annotate('', xy=(.405, .43), xytext=(.20, y), xycoords='axes fraction',
                arrowprops=dict(arrowstyle='->', lw=2, color='#555'))
ax.annotate('', xy=(.735, .43), xytext=(.655, .43), xycoords='axes fraction',
            arrowprops=dict(arrowstyle='->', lw=2, color='#555'))
ax.text(.5, .04, 'Observed: cell-state abundance, FN1/ITGB6 RNA, donor-level proximity\n'
                    'Inferred: ligand direction, receptor engagement, TF causality',
        ha='center', transform=ax.transAxes, fontsize=9)

ax = fig.add_subplot(gs[0, 3:])
selected = ['TEAD2', 'ZNF322', 'CREB3L1', 'ITGA3', 'MMP7', 'ABCA3']
z = ipf[ipf.gene.isin(selected) & ipf.modality.isin(['grn_tf', 'rna'])].copy()
z['cohort'] = z.contrast_id.str.split('__').str[0].replace({'Adams2020': 'Adams', 'Jaiswal2026_ILD': 'Jaiswal'})
z['row'] = z.modality.replace({'grn_tf': 'GRN activity', 'rna': 'RNA'}) + ' | ' + z.gene
order = [f'GRN activity | {g}' for g in ['TEAD2', 'ZNF322', 'CREB3L1']] + [f'RNA | {g}' for g in ['ITGA3', 'MMP7', 'ABCA3']]
mat = z.pivot_table(index='row', columns='cohort', values='log2fc', aggfunc='first').reindex(order)
sns.heatmap(np.sign(mat), cmap='vlag', center=0, vmin=-1, vmax=1,
            annot=mat.map(lambda v: '' if pd.isna(v) else f'{v:.2f}'), fmt='', cbar=False,
            linewidths=.8, ax=ax)
ax.set_title('B. The aberrant-epithelial TEAD/adhesion program reproduces across two IPF cohorts', loc='left', weight='bold')
ax.set_xlabel('IPF donor-pseudobulk cohort'); ax.set_ylabel('selected mechanistic readout')

ax = fig.add_subplot(gs[1, :2])
comp = pd.read_csv(ROOT / 'tables/xenium_composition_stats.tsv', sep='\t')
xe = pd.read_csv(ROOT / 'tables/xenium_expression_stats.tsv', sep='\t')
tests = xe[xe.group.str.contains('test')]
donor_comp = pd.read_csv(ROOT / 'tables/xenium_donor_composition.tsv', sep='\t')
comp_map = {'KRT5-/KRT17+': 'KRT state', 'Activated Fibrotic FBs': 'Fibrotic FB', 'AT2': 'AT2'}
dc = donor_comp[donor_comp.final_CT.isin(comp_map)].copy()
dc['feature'] = dc.final_CT.map(comp_map)
dc['status'] = np.where(dc.group.eq('Unaffected'), 'Control', 'IPF')
qmap = comp.set_index('cell_type').fdr_selected_states.to_dict()
dc['label'] = dc.apply(lambda r: f"{r.feature}\n{'*' if qmap[r.final_CT] < .05 else ''}q={qmap[r.final_CT]:.2g}", axis=1)
palette = {'Control': '#6c757d', 'IPF': '#d95f59'}
sns.boxplot(data=dc, x='label', y='logit_fraction', hue='status', palette=palette,
            whis=1.5, showfliers=False, width=.65, ax=ax)
sns.stripplot(data=dc, x='label', y='logit_fraction', hue='status', palette=palette,
              dodge=True, jitter=.12, size=3.4, alpha=.8, ax=ax)
handles, labels = ax.get_legend_handles_labels(); ax.legend(handles[:2], labels[:2], frameon=False, fontsize=7)
ax.set_xlabel('cell state'); ax.set_ylabel('donor logit cell fraction')
ax.set_title('C. Independent Xenium confirms KRT/fibroblast expansion\nand AT2 loss', loc='left', weight='bold', fontsize=10)
ax.tick_params(axis='x', rotation=20)

ax = fig.add_subplot(gs[1, 2:4])
xdir = Path('/Users/saljh8/Dropbox/LungMAP/Xenium-IPF/inputs')
meta_x = pd.read_csv(xdir / 'pseudobulk_metadata_final_CT_by_patient.tsv', sep='\t')
expr_x = pd.read_csv(xdir / 'pseudobulk_avgExpr_final_CT_by_patient-clean.txt', sep='\t', index_col=0).T
meta_x = meta_x[meta_x.pseudobulk_id.astype(str).isin(expr_x.index)].copy()
expr_x = expr_x.loc[meta_x.pseudobulk_id.astype(str)].copy(); expr_x.index = meta_x.index
expr_specs = [('SPP1+ Macrophages', 'FN1', 'SPP1 Mac: FN1'),
              ('Activated Fibrotic FBs', 'FN1', 'Fibrotic FB: FN1'),
              ('KRT5-/KRT17+', 'ITGB6', 'KRT state: ITGB6'),
              ('Transitional AT2', 'ITGB6', 'Trans. AT2: ITGB6')]
expr_parts = []
for cell, gene, label in expr_specs:
    ix = meta_x.final_CT.eq(cell)
    z = pd.DataFrame({'value': expr_x.loc[ix, gene].astype(float),
                      'status': np.where(meta_x.loc[ix, 'sample_affect'].eq('Unaffected'), 'Control', 'IPF')})
    stat = tests[(tests.cell_type.eq(cell)) & (tests.gene.eq(gene))].iloc[0]
    z['feature'] = f"{label}\n{'*' if stat.fdr_within_selected_tests < .05 else ''}q={stat.fdr_within_selected_tests:.2g}"
    expr_parts.append(z)
ed = pd.concat(expr_parts, ignore_index=True)
sns.boxplot(data=ed, x='feature', y='value', hue='status', palette=palette,
            whis=1.5, showfliers=False, width=.65, ax=ax)
sns.stripplot(data=ed, x='feature', y='value', hue='status', palette=palette,
              dodge=True, jitter=.12, size=3.4, alpha=.8, ax=ax)
handles, labels = ax.get_legend_handles_labels(); ax.legend(handles[:2], labels[:2], frameon=False, fontsize=7)
ax.set_xlabel('cell state and transcript'); ax.set_ylabel('donor-pseudobulk log2(CP10K+1)')
ax.set_title('D. Sender FN1 and receiver ITGB6 RNA increase\nin independent Xenium donors', loc='left', weight='bold', fontsize=10)
ax.tick_params(axis='x', rotation=25)

ax = fig.add_subplot(gs[1, 4:])
sp = pd.read_csv(ROOT / 'tables/xenium_spatial_pairs_donor_region.tsv', sep='\t')
parts = []
for a, b, label in [('KRT5-_KRT17+', 'Activated_Fibrotic_FBs', 'KRT ↔ fibrotic FB'),
                    ('KRT5-_KRT17+', 'SPP1+_Macrophages', 'KRT ↔ SPP1 macrophage')]:
    q = sp[(sp.cell_a == a) & (sp.cell_b == b)].copy()
    q['status'] = np.where(q.group.eq('Unaffected'), 'Control', 'IPF')
    q = q.groupby(['donor_id', 'status'], as_index=False).squidpy_z.median()
    disease = q[q.status.eq('IPF')].squidpy_z.dropna(); control = q[q.status.eq('Control')].squidpy_z.dropna()
    p = mannwhitneyu(disease, control, alternative='two-sided').pvalue
    q['pair'] = f'{label}\ntwo-sided MW P={p:.2g}'
    parts.append(q)
sr = pd.concat(parts, ignore_index=True)
sns.boxplot(data=sr, x='pair', y='squidpy_z', hue='status', palette=palette,
            whis=1.5, showfliers=False, width=.62, ax=ax)
sns.stripplot(data=sr, x='pair', y='squidpy_z', hue='status', dodge=True, jitter=.12,
              palette=palette, size=3.7, alpha=.85, ax=ax)
ax.axhline(0, color='black', lw=.7); ax.set_xlabel(''); ax.set_ylabel('median Squidpy neighborhood z per donor')
ax.set_title('E. Only KRT–fibroblast proximity differs significantly\nbetween IPF and control donors', loc='left', weight='bold', fontsize=10)
handles, labels = ax.get_legend_handles_labels(); ax.legend(handles[:2], labels[:2], frameon=False, fontsize=7)

fig.suptitle('An FN1-rich macrophage–fibroblast niche co-occurs with αv-integrin+ aberrant epithelium and TEAD activation in IPF',
             fontsize=15, weight='bold', y=.985)
caption = (
    'Comparison and interpretation. Panel B compares KRT5−/KRT17+ aberrant basal cells with AT2 cells within IPF donor pseudobulks; GRN values are RNA-inferred activity differences and RNA values are log2 fold changes. It tests state identity, not an IPF-versus-control effect. '
    'Panels C and D use an independent Xenium cohort: composition values are donor-level logit cell fractions (26 IPF/affected and 9 unaffected donors) and expression values are donor-level cell-type pseudobulk log2(CP10K+1). Boxes show the interquartile range, center lines the median, whiskers extend to 1.5×IQR, and every point is one biological donor pseudobulk. Two-sided Mann–Whitney U P values were Benjamini–Hochberg corrected within the prespecified composition-state or expression-test family; * denotes q<0.05 and exact q is printed. '
    'Panel E aggregates tissue regions to donors before a two-sided Mann–Whitney U test; boxes and points are defined as above. Proximity is nondirectional and abundance-sensitive. The KRT–fibroblast association is supported (P=0.0022), whereas KRT–SPP1-macrophage proximity is not (P=0.17). '
    'Together these data support co-occurrence and a testable FN1–integrin–TEAD model, but do not demonstrate ligand binding, signaling direction, or TF causality. Full cell-state and disease context for the TFs is shown in Figure 3.'
)
fig.text(.02, .012, textwrap.fill(caption, 225), ha='left', va='bottom', fontsize=8.5)
fig.subplots_adjust(left=.055, right=.985, top=.92, bottom=.24, hspace=.48, wspace=.90)
save(fig, BASEFIG, 'Figure1_IPF_FN1_integrin_TEAD_niche')


# Figure 2: explicitly define cell-state contrasts, assays, and inference limits.
inf = pd.read_csv(ROOT / 'tables/infection_multimodal_evidence.tsv', sep='\t')
inf['cohort'] = inf.contrast_id.str.extract(r'__(COVID19|Pneumonia)_')[0]
fig, axes = plt.subplots(1, 3, figsize=(15, 8))
panels = [('A. RNA-derived TF activity', 'grn_tf', ['ETV5', 'XBP1', 'CEBPA', 'MLXIPL']),
          ('B. Measured RNA', 'rna', ['SFTPA1', 'SFTPA2', 'SFTPC', 'MFSD2A']),
          ('C. RNA-imputed—not measured', None, ['Podoplanin', 'CD55', 'CD49f', 'PC(20:4/22:6)', 'PC(20:4/20:4)', 'PE(P-16:0/20:4)'])]
for ax, (title, modality, genes) in zip(axes, panels):
    z = inf[inf.gene.isin(genes) & (inf.modality.eq(modality) if modality else inf.modality.isin(['adt', 'lipid']))]
    mat = z.pivot_table(index='gene', columns='cohort', values='log2fc', aggfunc='mean').reindex(genes)
    lim = np.nanmax(np.abs(mat.to_numpy()))
    sns.heatmap(mat, cmap='vlag', center=0, vmin=-lim, vmax=lim, annot=True, fmt='.2f',
                linewidths=.6, cbar=False, ax=ax)
    ax.set_title(title, loc='left', weight='bold'); ax.set_xlabel('injury cohort'); ax.set_ylabel('')
    ax.tick_params(axis='y', rotation=0)
axes[0].text(.5, -.12, 'intermediate − AT2', ha='center', transform=axes[0].transAxes)
axes[1].text(.5, -.12, 'intermediate − AT2', ha='center', transform=axes[1].transAxes)
axes[2].text(.5, -.12, 'ADT: intermediate − AT2; lipids: intermediate − AT1', ha='center', transform=axes[2].transAxes)
fig.suptitle('Does the AT2→AT1 intermediate lose the AT2 regulatory and surfactant program in both COVID-19 and pneumonia?',
             fontsize=15, weight='bold', y=.98)
caption = (
    'Each column is a separate acute-injury cohort and each value is the indicated state contrast in native model units. Blue denotes lower and red higher values in the intermediate state relative to the stated reference. '
    'Panels A and B compare the same AT2→AT1 intermediate with mature AT2 cells: panel A is RNA-derived GRN activity and panel B is measured donor-pseudobulk RNA log2 fold change. Panel C contains predictions only: surface proteins are RNA-imputed ADT values versus AT2, whereas named phospholipids are RNA-imputed lipid values versus AT1. '
    'The repeated direction across COVID-19 and pneumonia supports a reproducible transitional-state association, not longitudinal lineage progression. ADT predictions require CITE-seq/flow validation and lipid predictions require isomer-resolved LC–MS/MS; neither constitutes orthogonal protein or lipid evidence.'
)
fig.text(.02, .015, textwrap.fill(caption, 220), ha='left', va='bottom', fontsize=8.7)
fig.subplots_adjust(left=.06, right=.98, top=.87, bottom=.19, wspace=.32)
save(fig, BASEFIG, 'Figure2_infection_surfactant_lipid_checkpoint')


# Retire the non-statistical strong/partial evidence matrix from the primary set.
for suffix in ['png', 'pdf']:
    path = FIG / f'Figure5_intuitive_evidence_matrix.{suffix}'
    if path.exists():
        retired = EXT / 'legacy_figures_before_context_rebuild' / f'RETIRED_Figure5_subjective_evidence_matrix.{suffix}'
        if not retired.exists():
            shutil.copy2(path, retired)
        path.unlink()


# Figures 8 and 10: show all power-aware query programs, including no coverage and mimicry.
queries = pd.read_csv(TAB / 'power_aware_drug_query_programs.tsv', sep='\t')
qorder = queries['query'].tolist()
qlabel = {r['query']: f"{r['disease']} | {r['state']}\n({int(r['n_genes'])} genes)" for _, r in queries.iterrows()}
era = pd.read_csv(TAB / 'existing_drug_state_scores_power_aware.tsv', sep='\t')
existing = ['nintedanib', 'pirfenidone', 'roflumilast', 'budesonide', 'fluticasone',
            'prednisone', 'dexamethasone', 'salmeterol', 'formoterol']
mat = era.pivot(index='compound', columns='query', values='median_rho').reindex(index=existing, columns=qorder)
nmat = era.pivot(index='compound', columns='query', values='n_signatures').reindex(index=existing, columns=qorder)
fig, ax = plt.subplots(figsize=(16, 8))
sns.heatmap(mat, cmap='vlag', center=0, vmin=-.2, vmax=.2, mask=mat.isna(), linewidths=.7,
            annot=np.where(mat.notna(), np.char.add(np.char.mod('%.2f', mat.fillna(0).to_numpy()),
                  np.char.add('\nn=', nmat.fillna(0).astype(int).astype(str).to_numpy())), ''), fmt='',
            xticklabels=[qlabel[q] for q in qorder], cbar_kws={'label': 'median Spearman ρ'}, ax=ax)
ax.set_facecolor('#d9d9d9')
ax.set_xlabel('power-aware disease program queried against LINCS signatures'); ax.set_ylabel('existing IPF/COPD therapy')
ax.set_title('Do existing IPF and COPD therapies reverse or mimic each replicated cell-state program in LINCS?',
             weight='bold', fontsize=14)
caption = (
    'Cells report the median Spearman correlation (ρ) across eligible LINCS L1000 signatures for the drug and query; n is the number of signatures. Blue (ρ<0) indicates transcriptomic reversal and red (ρ>0) mimicry. Gray denotes no eligible comparison after requiring at least 50 shared genes—not a null effect. '
    'Every prespecified power-aware query program is retained, including programs without drug-signature coverage. Query genes passed largest-effective-N discovery FDR q≤0.10, the modality-specific fold threshold, and independent concordant nominal P<0.05 without directional conflict. '
    'These correlations are cell-line perturbational screens; they neither establish target-level off-target binding nor predict clinical benefit or adverse events.'
)
fig.text(.02, .015, textwrap.fill(caption, 220), ha='left', va='bottom', fontsize=8.7)
fig.subplots_adjust(left=.13, right=.95, top=.90, bottom=.25)
save(fig, FIG, 'Figure8_drug_reversal_and_state_liabilities')


ranked = pd.read_csv(TAB / 'ranked_drug_predictions_power_aware.tsv', sep='\t')
candidates = ranked.sort_values('priority').compound_name.head(5).tolist()
engine = EXT / 'drug_engine_power_aware'
meta = pd.read_csv(engine / 'signature_meta.tsv', sep='\t', dtype=str, keep_default_na=False).set_index('lmd_key')
rows = []
for query in qorder:
    npz = np.load(engine / f'{query}.npz', allow_pickle=False)
    z = pd.DataFrame({'rho': npz['rho'], 'n_shared': npz['n_shared']}, index=npz['key'].astype(str)).join(meta[['compound_name']])
    z = z[(z.n_shared >= 50) & z.compound_name.isin(candidates)]
    for compound, g in z.groupby('compound_name'):
        rows.append({'compound': compound, 'query': query, 'median_rho': g.rho.median(), 'n_signatures': len(g),
                     'fraction_reversing': (g.rho < 0).mean()})
drug = pd.DataFrame(rows)
drug.to_csv(TAB / 'top_candidate_full_program_context.tsv', sep='\t', index=False)
dmat = drug.pivot(index='compound', columns='query', values='median_rho').reindex(index=candidates, columns=qorder)
dn = drug.pivot(index='compound', columns='query', values='n_signatures').reindex(index=candidates, columns=qorder)
fig, ax = plt.subplots(figsize=(16, 6.8))
sns.heatmap(dmat, cmap='vlag', center=0, vmin=-.35, vmax=.35, mask=dmat.isna(), linewidths=.7,
            annot=np.where(dmat.notna(), np.char.add(np.char.mod('%.2f', dmat.fillna(0).to_numpy()),
                  np.char.add('\nn=', dn.fillna(0).astype(int).astype(str).to_numpy())), ''), fmt='',
            xticklabels=[qlabel[q] for q in qorder], cbar_kws={'label': 'median Spearman ρ'}, ax=ax)
ax.set_facecolor('#d9d9d9')
ax.set_xlabel('all prespecified power-aware disease programs'); ax.set_ylabel('nominated perturbagen')
ax.set_title('Do nominated ROCK, SYK, and MEK perturbagens reverse all replicated disease-state programs—or only selected ones?',
             weight='bold', fontsize=14)
caption = (
    'Rows are the five prespecified mechanistic candidates from the full compound screen; columns retain every power-aware query program. Values are median Spearman ρ across LINCS L1000 signatures with at least 50 shared genes; n is the signature count. Blue indicates reversal, red mimicry, near-white weak association, and gray unavailable coverage. '
    'Candidates were nominated because they reached the top 0.5% reversal tail for both the confirmed IPF aberrant-epithelial and COPD interstitial-macrophage queries; this selection criterion makes those two columns discovery evidence rather than independent validation. The remaining columns provide an explicit specificity and liability context. '
    'Cell-line, dose, exposure, cytostasis, and tissue mismatch limit inference. Prospective dose-matched primary lung assays are required before therapeutic interpretation.'
)
fig.text(.02, .015, textwrap.fill(caption, 220), ha='left', va='bottom', fontsize=8.7)
fig.subplots_adjust(left=.13, right=.95, top=.87, bottom=.29)
save(fig, FIG, 'Figure10_power_aware_drug_reversals')

shutil.copy2(Path(__file__), EXT / 'scripts' / Path(__file__).name)
print({'rebuilt': [1, 2, 8, 10], 'retired': [5], 'drug_context_rows': len(drug)})
