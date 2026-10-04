"""Render the technical report from saved benchmark tables without reanalysis."""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import re
from datetime import datetime, timezone

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'benchmarking' / 'ND20_167_20261002'
METHODS = ['uncorrected', 'python_auto', 'python_0.10', 'python_0.20', 'python_0.30',
           'SoupX_0.10', 'SoupX_0.15', 'SoupX_0.20', 'SoupX_0.25']


def markdown(df, digits=3):
    def cell(x):
        if isinstance(x, (float, np.floating)):
            return f'{x:.{digits}f}' if np.isfinite(x) else '—'
        return str(x)
    return '\n'.join(['| ' + ' | '.join(map(str, df.columns)) + ' |',
                      '| ' + ' | '.join(['---'] * len(df.columns)) + ' |'] +
                     ['| ' + ' | '.join(map(cell, row)) + ' |'
                      for row in df.itertuples(index=False, name=None)])


def main():
    read = lambda name: pd.read_csv(OUT / (name + '.tsv'), sep='\t')
    summary = read('benchmark_summary').set_index('method').loc[METHODS]
    loss = read('transcript_loss')
    distances = read('embedding_distances')
    populations = pd.read_csv(OUT / 'population_counts.tsv', sep='\t', index_col=0)
    populations.index.name = 'Reference state'
    primary = populations.index[populations[['HSC', 'MPP']].min(axis=1) >= 100]
    shared = populations.index[populations[['HSC', 'MPP']].min(axis=1) >= 20]
    d = distances[distances.population.isin(primary)]
    auto = pd.read_csv(OUT / 'python_auto/ambient/soupx_summary.tsv', sep='\t')
    qc = read('qc')
    sample = qc.groupby('Library').agg(filtered_cells=('included', 'size'),
             analysis_cells=('included', 'sum'), median_input_UMIs=('UMIs', 'median')).reset_index()
    tables = {'SAMPLES': markdown(sample), 'POPULATIONS': markdown(populations.loc[shared].reset_index())}
    columns = {'predominant_DEG_gene_state_pairs': 'DEG pairs', 'predominant_unique_DEGs': 'Unique DEG genes',
               'DEG_reduction_percent': 'DEG reduction (%)', 'UMAP_within_scaled': 'UMAP D/R',
               'UMAP_reduction_percent': 'UMAP reduction (%)', 'PCA_within_scaled': 'PCA D/R',
               'PCA_mixing': 'PCA mixing ratio'}
    tables['PRIMARY'] = markdown(summary[list(columns)].rename(columns=columns).reset_index())
    comparisons = []
    weighted_loss = loss.groupby('method')[['original_UMIs', 'retained_UMIs']].sum()
    weighted_loss['loss'] = 100 * (1 - weighted_loss.retained_UMIs / weighted_loss.original_UMIs)
    for method in METHODS[5:]:
        comparisons.append({'SoupX comparator': method,
            'Python auto / SoupX DEG pairs': f'90 / {int(summary.loc[method, "predominant_DEG_gene_state_pairs"])}',
            'DEG difference (%)': 100 * (90 / summary.loc[method, 'predominant_DEG_gene_state_pairs'] - 1),
            'UMAP difference (%)': 100 * (summary.loc['python_auto', 'UMAP_within_scaled'] / summary.loc[method, 'UMAP_within_scaled'] - 1),
            'PCA difference (%)': 100 * (summary.loc['python_auto', 'PCA_within_scaled'] / summary.loc[method, 'PCA_within_scaled'] - 1),
            'Count loss difference (percentage points)': weighted_loss.loc['python_auto', 'loss'] - weighted_loss.loc[method, 'loss']})
    tables['HEAD_TO_HEAD'] = markdown(pd.DataFrame(comparisons), 2)
    tables['AUTO'] = markdown(auto[['library', 'rho', 'ambient_cells_used', 'rho_eval_cells',
          'rho_eval_baseline_residual', 'rho_eval_selected_residual', 'rho_eval_fraction_removed',
          'rho_eval_selected_score']].rename(columns={'rho_eval_baseline_residual':'r(0)',
          'rho_eval_selected_residual':'r(selected)', 'rho_eval_fraction_removed':'f(selected)',
          'rho_eval_selected_score':'S(selected)'}), 6)
    tables['DE_STATES'] = markdown(read('DEGs_by_state'))
    tables['FOLD'] = markdown(summary[['DEG_FC_1.1', 'predominant_DEG_gene_state_pairs',
          'DEG_FC_1.5', 'DEG_FC_2.0']].rename(columns={'predominant_DEG_gene_state_pairs':'DEG_FC_1.2'}).reset_index())
    for token, name, value in [('NORMALIZED', 'normalized_effect_sensitivity', 'CP10k_DEGs'),
                               ('FIXED_FAMILY', 'fixed_BH_family_sensitivity', 'fixed_family_DEGs')]:
        tables[token] = markdown(read(name).pivot(index='method', columns='population', values=value).reindex(METHODS).reset_index())
    tables['SEEDS'] = markdown(d.groupby(['method', 'seed']).UMAP_distance_within_scaled.mean().unstack().reindex(METHODS).rename(columns={0:'Seed 0',1:'Seed 1',2:'Seed 2'}).reset_index())
    state_geometry = d.groupby(['method', 'population']).agg(
        UMAP_D=('UMAP_distance', 'mean'), UMAP_D_R=('UMAP_distance_within_scaled','mean'),
        PCA_D_R=('PCA_distance_within_scaled','mean'), mixing=('PCA_neighbor_mixing','mean'))
    tables['GEOMETRY'] = markdown(state_geometry.reindex(pd.MultiIndex.from_product([METHODS, primary], names=['method','population'])).reset_index())
    loss_table = loss[['method','capture','filtered_cells','original_UMIs','retained_UMIs','loss_percent',
                     'QC_loss_percent','median_cell_loss_percent']].copy()
    for col in ['original_UMIs','retained_UMIs']:
        loss_table[col] = loss_table[col].map(lambda x: f'{x:,.0f}')
    tables['LOSS'] = markdown(loss_table, 2)
    tables['STABILITY'] = markdown(read('assignment_stability'), 2)
    tables['STATE_STABILITY'] = markdown(read('primary_state_reassignment'), 2)
    tables['RUNTIME'] = markdown(loss[loss.capture=='HSC'][['method','correction_seconds','alignment_seconds']], 3)
    packages = ['numpy','scipy','pandas','anndata','scanpy','scikit-learn','umap-learn','numba',
                'pynndescent','h5py','matplotlib','statsmodels','tqdm','threadpoolctl','llvmlite']
    environment = {'recorded_at_utc':datetime.now(timezone.utc).isoformat(),
        'scope':'Report revision environment; supplementary to preserved original analysis provenance',
        'python':platform.python_version(), 'platform':platform.platform(),
        'packages':{p:importlib.metadata.version(p) for p in packages}}
    (OUT / 'report_environment.json').write_text(json.dumps(environment, indent=2) + '\n')
    purposes = ['Numerical arrays and vectorized subtraction','Sparse matrices, Matrix Market I/O and statistical distributions',
        'Annotations, summaries and tabular I/O','Annotated sparse count matrices and H5AD persistence',
        'Normalization, Wilcoxon ranking, PCA/neighbor/UMAP orchestration',
        'PCA and within-state nearest-neighbor diagnostics','UMAP optimization',
        'Compiled acceleration for UMAP and neighbor search','Approximate neighbor construction',
        '10x HDF5 and H5AD I/O','Scientific figure export','Installed cellHarmony dependency; primary BH uses local NumPy implementation',
        'Progress reporting','Numerical thread-pool control','Numba compiler dependency']
    dependency = pd.DataFrame({'Library':packages,'Version':[environment['packages'][p] for p in packages],'Role':purposes})
    tables['DEPENDENCIES'] = markdown(dependency)
    template = ROOT / 'benchmarking' / 'ambient_subtract_benchmarking.template.md'
    report = template.read_text()
    for token, table in tables.items():
        assert report.count('@@' + token + '@@') == 1, token
        report = report.replace('@@' + token + '@@', table)
    assert not re.search(r'@@\w+@@', report)
    simulation = OUT.parent / 'ambient_subtract_release_simulation.md'
    if simulation.exists():
        report = report.replace('## Technical conclusions and threshold assessment',
                                '## Empirical-capture conclusions and threshold assessment')
        report = report.replace('## Dependencies and reproducibility',
                                simulation.read_text() + '\n\n## Dependencies and reproducibility')
        abstract = OUT.parent / 'ND20_167_release_20261003/abstract_extension.md'
        if abstract.exists():
            report = report.replace('## Background and analytical rationale',
                                    abstract.read_text() + '\n\n## Background and analytical rationale')
        report = report.replace('Analysis completed 2 October 2026. Report revised to document the statistical estimands, algorithmic assumptions and comparator provenance. The analytical results and production scALABLE thresholds were retained.',
            'Empirical analysis completed 2 October 2026; population-weighted RNA-release simulations completed 3 October 2026. Production correction parameters were retained.')
    target = OUT.parent / 'ambient_subtract_benchmarking.md'
    empirical_legacy = OUT / 'legacy' / 'ambient_subtract_benchmarking_empirical_20261002.md'
    if simulation.exists() and target.exists() and not empirical_legacy.exists():
        empirical_legacy.parent.mkdir(exist_ok=True)
        empirical_legacy.write_bytes(target.read_bytes())
    legacy = OUT / 'legacy' / 'ambient_subtract_benchmarking_initial.md'
    if target.exists() and not legacy.exists():
        legacy.parent.mkdir(exist_ok=True)
        legacy.write_bytes(target.read_bytes())
    target.write_text(report)
    sources = [Path(__file__), template, ROOT/'evaluation/summarize_nd20_167.py']
    if simulation.exists():sources.append(simulation)
    (OUT / 'report_revision.json').write_text(json.dumps({
        'recorded_at_utc':environment['recorded_at_utc'], 'scope':'Report rendering only; saved analytical results unchanged',
        'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        'style_reference':'/Users/saljh8/Downloads/nihms-2119165.pdf'}, indent=2) + '\n')
    print(target)


if __name__ == '__main__':
    main()
