#!/usr/bin/env python3
"""Replace claim-selected graphics with complete, statistically annotated context plots."""

from pathlib import Path
import shutil
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch
import seaborn as sns
from openpyxl import load_workbook


ROOT = Path('/Users/saljh8/Dropbox/LungMAP/Discovery/codex_nature_predictions_20260927')
EXT = ROOT / 'extended'
TAB = EXT / 'tables'
FIG = EXT / 'figures'
BASEFIG = ROOT / 'figures'
LEGACY = EXT / 'legacy_figures_before_context_rebuild'
LEGACY.mkdir(exist_ok=True)

for path in list(BASEFIG.glob('Figure*.png')) + list(BASEFIG.glob('Figure*.pdf')) + list(FIG.glob('Figure*.png')) + list(FIG.glob('Figure*.pdf')):
    target = LEGACY / path.name
    if not target.exists():
        shutil.copy2(path, target)

sns.set_theme(style='white', font_scale=.78)
ASSOC = TAB / 'power_aware_associations.tsv.gz'
TF = ['TEAD2', 'TEAD3', 'ZNF322', 'CREB3L1', 'ARNT2', 'FOXA2', 'SREBF2', 'KLF5', 'TP63', 'FOXQ1']
SHARED = ['RGCC', 'CXCL12', 'CFI', 'SGK1', 'USP53', 'PLPP3', 'CD163', 'PLIN2']


def load_context(genes, modalities):
    pieces = []
    cols = ['disease', 'modality', 'population', 'gene', 'discovery_cohort',
            'discovery_n_case', 'discovery_n_control', 'discovery_effective_n',
            'discovery_log2fc', 'discovery_p', 'discovery_fdr', 'discovery_pass',
            'n_confirmation_cohorts', 'confirmation_cohorts',
            'n_conflicting_cohorts', 'conflicting_cohorts', 'cohorts_tested',
            'total_case_donors', 'total_control_donors', 'independently_confirmed']
    for chunk in pd.read_csv(ASSOC, sep='\t', usecols=cols, chunksize=100000, low_memory=False):
        keep = chunk.gene.isin(genes) & chunk.modality.isin(modalities)
        keep &= ~chunk.population.str.contains('__vs__', regex=False)
        if keep.any():
            pieces.append(chunk.loc[keep].copy())
    return pd.concat(pieces, ignore_index=True)


context = load_context(set(TF + SHARED), {'grn_tf', 'rna'})
context.to_csv(TAB / 'complete_context_for_highlighted_TFs_and_genes.tsv.gz', sep='\t', index=False)

all_assoc = pd.read_csv(ASSOC, sep='\t', low_memory=False)
direct = all_assoc[~all_assoc.population.str.contains('__vs__', regex=False)].copy()


def ordered_populations(frame):
    preferred = [
        'AT1', 'AT2', 'AT2-AT1 int.', 'AT2-prolif', 'Aberrant basal', 'Basal', 'Hillock basal',
        'Hillock squamous', 'Club-Bronch', 'Club-TB', 'Club-nasal', 'Goblet', 'Goblet-nasal',
        'Ciliated-Bronch', 'Ciliated-TB', 'Ciliated-axon', 'Deuterosomal', 'RASC', 'SAEC',
        'PNEC', 'Tuft', 'Ionocyte', 'Serous-nasal', 'SMG Duct',
        'CAP1', 'CAP1-int', 'CAP2', 'CAP2-int', 'PAEC', 'PArEC', 'PCEC', 'PVEC',
        'SVEC', 'SVEC-act', 'LEC', 'LEC-cycling',
        'AdvFB', 'IntFB', 'LipoFB1', 'PBFB', 'SPFB', 'AF', 'SCMF', 'iALF', 'Mesothelial',
        'Pericyte', 'ASMC', 'LSMC', 'VSMC', 'cVSMC',
        'AM', 'AM-lipid', 'AM-prolif', 'MT+ AM', 'IM', 'tMDM', 'ADM', 'cADM', 'iMON',
        'pMON', 'Neutrophil', 'Mast', 'cDC1', 'maDC', 'pDC', 'Langerhans',
        'B', 'Plasma', 'CD4 T', 'CD8 T', 'Treg', 'T-prolif', 'NK', 'iNKT', 'ILC1',
        'HSC', 'Erythrocyte'
    ]
    present = set(frame.population.unique())
    return [x for x in preferred if x in present] + sorted(present - set(preferred))


def matrix_panel(ax, frame, genes, populations, title, modality_label):
    columns = pd.MultiIndex.from_product([genes, ['COPD', 'IPF']], names=['feature', 'disease'])
    idx = pd.Index(populations, name='cell state')
    val = frame.pivot_table(index='population', columns=['gene', 'disease'], values='discovery_log2fc', aggfunc='first').reindex(index=idx, columns=columns)
    q = frame.pivot_table(index='population', columns=['gene', 'disease'], values='discovery_fdr', aggfunc='first').reindex(index=idx, columns=columns)
    rep = frame.pivot_table(index='population', columns=['gene', 'disease'], values='independently_confirmed', aggfunc='first').reindex(index=idx, columns=columns)
    conflict = frame.pivot_table(index='population', columns=['gene', 'disease'], values='n_conflicting_cohorts', aggfunc='first').reindex(index=idx, columns=columns)

    lim = float(np.nanquantile(np.abs(val.to_numpy()), .98))
    if not np.isfinite(lim) or lim == 0:
        lim = 1
    sns.heatmap(val, ax=ax, cmap='vlag', center=0, vmin=-lim, vmax=lim,
                mask=val.isna(), linewidths=.22, linecolor='#eeeeee',
                cbar_kws={'label': f'largest-effective-N cohort effect\n({modality_label}; disease − control)'})
    ax.set_facecolor('#d9d9d9')
    labels = []
    for gene, disease in columns:
        labels.append(f'{gene}\n{disease}')
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=7)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=6.5)
    ax.set_xlabel('Candidate factor; adjacent columns show both diseases')
    ax.set_ylabel('All tested lung cell states (unfiltered)')
    ax.set_title(title, loc='left', weight='bold', fontsize=12)
    for y, pop in enumerate(populations):
        for x, col in enumerate(columns):
            if pd.isna(val.loc[pop, col]):
                continue
            symbol = ''
            if q.loc[pop, col] < .05:
                symbol = '*'
            elif q.loc[pop, col] <= .10:
                symbol = '†'
            if bool(rep.loc[pop, col]) if not pd.isna(rep.loc[pop, col]) else False:
                symbol += '‡'
            if conflict.loc[pop, col] > 0:
                symbol += '×'
            if symbol:
                ax.text(x + .5, y + .5, symbol, ha='center', va='center', fontsize=6.5,
                        color='black', weight='bold')
    for x in range(2, len(columns), 2):
        ax.axvline(x, color='black', lw=.75)


tf_context = context[context.gene.isin(TF)]
pops = ordered_populations(tf_context)
fig, axes = plt.subplots(1, 2, figsize=(25, 27), gridspec_kw={'wspace': .28})
matrix_panel(axes[0], tf_context[tf_context.modality.eq('grn_tf')], TF, pops,
             'A. Are inferred TF-activity changes restricted to the proposed cell state?',
             'GRN activity difference, native score units')
matrix_panel(axes[1], tf_context[tf_context.modality.eq('rna')], TF, pops,
             'B. Does measured TF RNA support the inferred activity change?',
             'RNA log2 fold change')
fig.suptitle('Nominated TF programs are not uniformly cell-state specific, and inferred activity often diverges from TF RNA',
             fontsize=18, weight='bold', y=.975)
legend = (
    'Question and comparison. Each row is an annotated lung cell state; every available direct disease-versus-control test is shown for COPD and IPF. '
    'Columns pair the two diseases for each nominated TF. Panel A reports RNA-derived GRN activity; panel B reports measured TF RNA. '
    'Color is the signed effect in the cohort with the largest effective donor N = n(disease)×n(control)/[n(disease)+n(control)]: red, higher in disease; blue, lower; gray, not tested/available. '
    'Symbols: *, discovery Benjamini–Hochberg FDR q<0.05; †, 0.05≤q≤0.10; ‡, additionally confirmed in an independent cohort with concordant direction, nominal two-sided P<0.05, the modality-specific fold criterion, and no significant opposing cohort; ×, at least one independent cohort had a significant opposite direction. '
    'RNA discovery also required |log2FC|≥log2(1.2); no fold threshold was imposed on GRN scores. Nonsignificant and opposite-direction estimates are deliberately retained. '
    'The GRN panel is expression-inferred and does not establish TF protein activity, binding, or causality. Exact cohorts, donor counts, P values, q values, and conflicts are supplied in the companion table.'
)
fig.text(.02, .006, textwrap.fill(legend, 235), ha='left', va='bottom', fontsize=9)
fig.subplots_adjust(left=.07, right=.96, top=.92, bottom=.08, wspace=.34)
for ext in ['png', 'pdf']:
    fig.savefig(BASEFIG / f'Figure3_TF_candidates_and_COPD_endothelium.{ext}', dpi=300 if ext == 'png' else None, bbox_inches='tight')
plt.close(fig)


# Complete disease-state RNA summary: all states, all tested genes, both directions.
rna = direct[direct.modality.eq('rna')].copy()
summary = (rna.groupby(['disease', 'population'], as_index=False)
           .agg(genes_tested=('gene', 'size'),
                discovery_up=('discovery_pass', lambda x: 0),
                median_effective_n=('discovery_effective_n', 'median')))
counts = []
for (disease, pop), z in rna.groupby(['disease', 'population']):
    discovered = z.discovery_pass
    confirmed = z.independently_confirmed
    counts.append(dict(disease=disease, population=pop, genes_tested=len(z),
                       discovery_up=int((discovered & (z.discovery_log2fc > 0)).sum()),
                       discovery_down=int((discovered & (z.discovery_log2fc < 0)).sum()),
                       confirmed_up=int((confirmed & (z.discovery_log2fc > 0)).sum()),
                       confirmed_down=int((confirmed & (z.discovery_log2fc < 0)).sum()),
                       median_effective_n=float(z.discovery_effective_n.median())))
state = pd.DataFrame(counts)
for c in ['discovery_up', 'discovery_down', 'confirmed_up', 'confirmed_down']:
    state[c + '_pct'] = 100 * state[c] / state.genes_tested
state.to_csv(TAB / 'complete_cell_state_RNA_discovery_replication_summary.tsv', sep='\t', index=False)

state_pops = ordered_populations(state)
fig, axes = plt.subplots(1, 2, figsize=(15, 24), sharey=True)
for ax, disease in zip(axes, ['COPD', 'IPF']):
    z = state[state.disease.eq(disease)].set_index('population').reindex(state_pops)
    y = np.arange(len(state_pops))
    ax.barh(y, -z.discovery_down_pct, color='#9ecae1', label='discovered lower')
    ax.barh(y, z.discovery_up_pct, color='#f4a6a1', label='discovered higher')
    ax.barh(y, -z.confirmed_down_pct, color='#2171b5', label='independently confirmed lower')
    ax.barh(y, z.confirmed_up_pct, color='#cb181d', label='independently confirmed higher')
    ax.axvline(0, color='black', lw=.8)
    ax.set_yticks(y, state_pops, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel('percentage of all tested RNA features\nlower in disease  ←  0  →  higher in disease')
    ax.set_title(f'{disease}: all tested cell states', weight='bold')
    ax.grid(axis='x', color='#dddddd', lw=.5)
axes[0].set_ylabel('Annotated lung cell state')
handles = [Patch(color='#9ecae1', label='discovered lower'), Patch(color='#2171b5', label='confirmed lower'),
           Patch(color='#f4a6a1', label='discovered higher'), Patch(color='#cb181d', label='confirmed higher')]
fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .963), ncol=4, frameon=False)
fig.suptitle('IPF shows broader cell-state RNA remodeling than COPD, while independent confirmation remains partial',
             fontsize=17, weight='bold', y=.995)
legend = (
    'Each bar uses every RNA feature tested in the indicated disease–cell-state comparison; no cell state or effect direction was selected for display. '
    'Pale bars are largest-effective-N cohort discoveries (Benjamini–Hochberg FDR q≤0.10 and |log2FC|≥log2[1.2]); dark overlays are the subset confirmed in at least one independent cohort with concordant direction and nominal two-sided P<0.05, with no significant opposing cohort. '
    'Bar length is the percentage of all tested genes, not the raw number of hits, so states with different feature coverage are comparable. Left indicates lower expression in disease and right indicates higher expression. '
    'A large bar denotes a broad association program, not a causal or lineage-specific mechanism. Exact tested denominators, hit counts, median effective donor N, and state-level values are provided in the companion table.'
)
fig.text(.02, .006, textwrap.fill(legend, 185), ha='left', va='bottom', fontsize=9)
fig.tight_layout(rect=[0, .065, 1, .955])
for ext in ['png', 'pdf']:
    fig.savefig(BASEFIG / f'Figure4_replicated_disease_enriched_states.{ext}', dpi=300 if ext == 'png' else None, bbox_inches='tight')
plt.close(fig)


# Full cell-state context for the highlighted cross-disease genes.
shared = context[(context.modality.eq('rna')) & context.gene.isin(SHARED)]
shared_pops = ordered_populations(shared)
fig, ax = plt.subplots(figsize=(14, 24))
matrix_panel(ax, shared, SHARED, shared_pops,
             'Are the highlighted shared IPF–COPD RNA associations specific to the nominated cell state?',
             'RNA log2 fold change')
fig.suptitle('Genes reproduced in IPF and COPD frequently extend beyond their nominated cell state', fontsize=17, weight='bold', y=.995)
legend = (
    'Every available direct disease-versus-control RNA test is shown across all annotated cell states; columns pair COPD and IPF for each gene. '
    'Red indicates higher and blue lower expression in disease in the largest-effective-N cohort; gray means unavailable. * marks discovery FDR q<0.05; † marks 0.05≤q≤0.10; both also require |log2FC|≥log2(1.2). ‡ marks independent concordant confirmation at nominal two-sided P<0.05 with no significant opposite cohort; × marks a significant directional conflict. '
    'The originally highlighted gene–cell-state pairs are therefore evaluated against all other tested cell states and both diseases rather than displayed as positive examples alone. '
    'The result supports cell-state specificity only when the nominated row is exceptional relative to this full reference distribution. It does not identify a shared upstream cause.'
)
fig.text(.02, .006, textwrap.fill(legend, 185), ha='left', va='bottom', fontsize=9)
fig.tight_layout(rect=[0, .055, 1, .975])
for ext in ['png', 'pdf']:
    fig.savefig(FIG / f'Figure7_shared_IPF_COPD_biology.{ext}', dpi=300 if ext == 'png' else None, bbox_inches='tight')
plt.close(fig)


# Cohort power figure: all cohorts and all state-level replication yields.
audit = pd.read_csv(TAB / 'cohort_power_summary.tsv', sep='\t')
cohort = (audit.groupby(['disease', 'cohort'], as_index=False)
          .agg(raw_total_n=('raw_total_n', 'max'), effective_n=('effective_n', 'max'),
               minimum_case_n=('n_case', 'min'), minimum_control_n=('n_control', 'min')))
cohort['label'] = cohort.disease + ' | ' + cohort.cohort
cohort = cohort.sort_values(['disease', 'effective_n'])
rep = state.copy()
rep['discovered_total'] = rep.discovery_up + rep.discovery_down
rep['confirmed_total'] = rep.confirmed_up + rep.confirmed_down
rep['confirmation_pct_of_discoveries'] = 100 * rep.confirmed_total / rep.discovered_total.replace(0, np.nan)

fig, axes = plt.subplots(1, 2, figsize=(17, 9), gridspec_kw={'width_ratios': [.85, 1.35]})
axes[0].barh(cohort.label, cohort.raw_total_n, color='#cbd5e1', label='largest raw total N')
axes[0].barh(cohort.label, cohort.effective_n, color='#264653', label='largest effective N')
axes[0].set_xlabel('donors')
axes[0].set_title('A. How much independent donor information is available?', loc='left', weight='bold')
axes[0].legend(frameon=False)

colors = {'COPD': '#d95f59', 'IPF': '#457b9d'}
for disease, z in rep.groupby('disease'):
    axes[1].scatter(z.median_effective_n, z.confirmation_pct_of_discoveries,
                    s=np.clip(z.discovered_total, 5, 250), alpha=.65,
                    color=colors[disease], edgecolor='white', linewidth=.4, label=disease)
    for _, r in z.nlargest(6, 'discovered_total').iterrows():
        axes[1].annotate(r.population, (r.median_effective_n, r.confirmation_pct_of_discoveries),
                         xytext=(3, 2), textcoords='offset points', fontsize=7)
axes[1].set_xlabel('median discovery effective donor N across tested RNA features')
axes[1].set_ylabel('independently confirmed discoveries / all discoveries (%)')
axes[1].set_title('B. Does discovery power predict independent confirmation yield?', loc='left', weight='bold')
axes[1].legend(frameon=False)
axes[1].grid(color='#dddddd', lw=.5)
fig.suptitle('Unequal case–control balance constrains effective donor information, and confirmation yield varies widely by cell state',
             fontsize=16, weight='bold', y=.985)
legend = (
    'Panel A reports each cohort’s maximum raw donor total and effective N = n(disease)×n(control)/[n(disease)+n(control)]; effective N is limited by the smaller group and is the prespecified rule used to choose the discovery cohort separately for each disease, modality, cell state, and feature. '
    'Panel B includes every direct RNA cell-state comparison. Point area is proportional to the number of discovery features; only the six largest programs per disease are labeled to prevent text overlap. Confirmation yield is the percentage of discoveries reproduced in an independent cohort with concordant direction, nominal two-sided P<0.05, the same fold criterion, and no significant directional conflict. '
    'Discovery required Benjamini–Hochberg FDR q≤0.10 and |log2FC|≥log2(1.2). States with zero discoveries have undefined yield and are omitted only from panel B.'
)
fig.text(.02, .012, textwrap.fill(legend, 215), ha='left', va='bottom', fontsize=9)
fig.tight_layout(rect=[0, .10, 1, .95])
for ext in ['png', 'pdf']:
    fig.savefig(FIG / f'Figure9_power_aware_discovery_confirmation.{ext}', dpi=300 if ext == 'png' else None, bbox_inches='tight')
plt.close(fig)


# Remove the subjective score from compact tables and both Excel deliverables.
claim_path = TAB / 'adversarial_claim_inventory.tsv'
claims = pd.read_csv(claim_path, sep='\t')
if 'rigor_score_0_10' in claims.columns:
    claims = claims.drop(columns='rigor_score_0_10')
claims = claims.sort_values(['domain', 'claim_id'])
claims.to_csv(claim_path, sep='\t', index=False)

for book_path in [TAB / 'intuitive_evidence_matrix.xlsx', TAB / 'adversarial_evidence_review.xlsx']:
    wb = load_workbook(book_path)
    for ws in wb.worksheets:
        headers = [cell.value for cell in ws[1]]
        if 'rigor_score_0_10' in headers:
            ws.delete_cols(headers.index('rigor_score_0_10') + 1)
    # Append/replace compact context sheets without rewriting the 110k-row atlas sheet.
    for name in ['TF context', 'State RNA summary']:
        if name in wb.sheetnames:
            del wb[name]
    ws = wb.create_sheet('TF context')
    ws.append(list(context.columns))
    for row in context.itertuples(index=False, name=None):
        ws.append(list(row))
    ws = wb.create_sheet('State RNA summary')
    ws.append(list(state.columns))
    for row in state.itertuples(index=False, name=None):
        ws.append(list(row))
    wb.save(book_path)

note = EXT / 'FIGURE_INTERPRETATION_GUIDE.md'
note.write_text('''# Context-complete figure revision

The revised figures use biological-question titles and retain negative, nonsignificant, and opposing estimates. Figure 3 shows every available cell state and both diseases for nominated TFs, separately for inferred GRN activity and measured RNA. Figure 4 shows the complete signed RNA discovery and confirmation burden for every disease–cell-state comparison. Figure 7 tests highlighted shared genes against all other cell states and both diseases. Figure 9 documents donor imbalance and replication yield.

Symbols and statistical rules are defined in each figure legend. The former subjective “adversarial rigor score” has been removed from figures, TSVs, and Excel workbooks.
''')

shutil.copy2(Path(__file__), EXT / 'scripts' / Path(__file__).name)
print({'context_rows': len(context), 'cell_state_summary_rows': len(state),
       'figures_replaced': ['Figure3', 'Figure4', 'Figure7', 'Figure9']})
