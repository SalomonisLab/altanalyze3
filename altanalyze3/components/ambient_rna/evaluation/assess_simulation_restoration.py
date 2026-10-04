"""Assess ambient gene burden and introduced DEGs against the corrected baseline.

Uses completed simulations and actual correction outputs; no molecular labels.
"""
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse

from .. import ambient_subtract as ambient
from .simulation_analysis import OUT, fingerprint
from .write_nd20_167_report import markdown


def deg_set(directory, states):
    t = pd.read_csv(directory / 'all_gene_statistics.tsv.gz', sep='\t')
    t = t[t.population.isin(states) & (t.fdr < .05) &
          (t.log2fc.abs() > np.log2(1.2))]
    return set(zip(t.population, t.gene))


def main(out=OUT, models=('heterogeneous','depth_proportional'), write_report=True):
    OUT = out
    base = ad.read_h5ad(OUT / 'baseline_full.h5ad')
    populations = pd.read_csv(OUT / 'population_counts.tsv', sep='\t', index_col=0)
    states = populations.index[populations.min(axis=1) >= 100].tolist()
    baseline_de = deg_set(OUT / 'baseline', states)
    summary = pd.read_csv(OUT / 'benchmark_summary.tsv', sep='\t')
    scenarios = summary.loc[(summary.model.isin(models)) &
                            (summary.rho >= .20), 'scenario'].unique()
    gene_rows, de_rows = [], []
    for scenario in scenarios:
        directory = OUT / scenario
        observed = ad.read_h5ad(directory / 'contaminated.h5ad')
        assert observed.obs_names.equals(base.obs_names)
        assert observed.var_names.equals(base.var_names)
        contaminated_de = deg_set(directory / 'contaminated', states)
        introduced = contaminated_de - baseline_de
        methods = summary.loc[summary.scenario == scenario, 'method'].tolist()
        for method in methods:
            if (directory / method / 'corrected_counts.npz').exists():
                corrected = sparse.load_npz(directory / method / 'corrected_counts.npz')
                assert fingerprint(corrected) == (directory / method / 'counts.sha256').read_text().strip()
            elif method == 'contaminated':
                corrected = observed.X
            elif method.startswith('python'):
                inp = observed.copy()
                result = ambient.process_anndata(inp, rho='auto' if method == 'python_auto'
                    else float(summary.loc[summary.scenario == scenario, 'rho'].iloc[0]),
                    outdir=directory / 'restoration_recheck' / method,
                    write_individual=False, write_merged=False, inplace=True,
                    store_corrected_layer=False)
                corrected = result.X
                assert fingerprint(corrected) == (directory / method / 'counts.sha256').read_text().strip()
            else:
                corrected = sparse.vstack([io.mmread(directory / f'R_{cap}' /
                    f'{method}.mtx').T.tocsr().astype(np.float32)
                    for cap in ('HSC', 'MPP')], format='csr')
                assert fingerprint(corrected) == (directory / method / 'counts.sha256').read_text().strip()
            recovered_de = deg_set(directory / method, states)
            for state in ['ALL'] + states:
                restrict = lambda pairs: pairs if state == 'ALL' else {p for p in pairs if p[0] == state}
                target, injected, remaining = map(restrict, (baseline_de, introduced, recovered_de))
                resolved = injected - remaining
                de_rows.append(dict(scenario=scenario, method=method, population=state,
                    introduced_DEG_pairs=len(injected), resolved_introduced_DEG_pairs=len(resolved),
                    resolved_percent=100 * len(resolved) / len(injected) if injected else np.nan,
                    baseline_DEG_pairs=len(target), preserved_baseline_DEG_pairs=len(target & remaining),
                    new_after_correction_not_introduced=len(remaining - target - injected)))
            for cap in ('HSC', 'MPP'):
                support = pd.read_csv(directory / f'true_profile_{cap}.tsv', sep='\t').probability.to_numpy() > 0
                for state in states:
                    mask = ((base.obs.Library.to_numpy() == cap) &
                            (base.obs.fixed_population.to_numpy() == state) &
                            base.obs.analysis_included.to_numpy())
                    b = np.asarray(base.layers['counts'][mask].sum(0), dtype=np.float64).ravel()
                    y = np.asarray(observed.X[mask].sum(0), dtype=np.float64).ravel()
                    c = np.asarray(corrected[mask].sum(0), dtype=np.float64).ravel()
                    added, removed = y - b, y - c
                    tv_before = .5 * np.abs(y / y.sum() - b / b.sum()).sum()
                    tv_after = .5 * np.abs(c / c.sum() - b / b.sum()).sum()
                    gene_rows.append(dict(scenario=scenario, method=method, capture=cap,
                        population=state, cells=int(mask.sum()), added_counts=added.sum(),
                        removed_counts=removed.sum(), observed_loss_percent=100 * removed.sum() / y.sum(),
                        added_fraction_percent=100 * added.sum() / y.sum(),
                        ambient_support_removed_percent_of_added=100 * removed[support].sum() / added.sum(),
                        off_support_removed_counts=removed[~support].sum(),
                        off_support_removed_percent_of_baseline=100 * removed[~support].sum() / b[~support].sum()
                            if b[~support].sum() else np.nan,
                        ambient_removal_profile_TV=.5 * np.abs(removed / removed.sum() - added / added.sum()).sum()
                            if removed.sum() else np.nan,
                        baseline_expression_TV_before=tv_before, baseline_expression_TV_after=tv_after,
                        expression_distortion_reduction_percent=100 * (1 - tv_after / tv_before)))
            print(scenario, method, 'restoration assessed', flush=True)
    genes, de = pd.DataFrame(gene_rows), pd.DataFrame(de_rows)
    genes.to_csv(OUT / 'cell_type_ambient_restoration.tsv', sep='\t', index=False)
    de.to_csv(OUT / 'introduced_DEG_resolution.tsv', sep='\t', index=False)
    table = de[de.population == 'ALL'].copy()
    table['model'] = table.scenario.str.split('_rho').str[0]
    table['nominal_rho'] = table.scenario.str.extract(r'rho([0-9.]+)_')[0].astype(float)
    grouped = table.groupby(['model', 'nominal_rho', 'method']).agg(
        draws=('scenario', 'size'), introduced=('introduced_DEG_pairs', 'mean'),
        resolved=('resolved_introduced_DEG_pairs', 'mean'),
        resolution_percent=('resolved_percent', 'mean'),
        baseline_preserved=('preserved_baseline_DEG_pairs', 'mean'),
        newly_created=('new_after_correction_not_introduced', 'mean')).reset_index()
    grouped.to_csv(OUT / 'introduced_DEG_resolution_summary.tsv', sep='\t', index=False)
    gt = genes.copy()
    gt['model'] = gt.scenario.str.split('_rho').str[0]
    gt['nominal_rho'] = gt.scenario.str.extract(r'rho([0-9.]+)_')[0].astype(float)
    gsummary = gt.groupby(['model', 'nominal_rho', 'method']).agg(
        gene_state_capture_observations=('cells', 'size'),
        mean_loss_percent=('observed_loss_percent', 'mean'),
        mean_ambient_support_removal_percent=('ambient_support_removed_percent_of_added', 'mean'),
        mean_off_support_baseline_depletion_percent=('off_support_removed_percent_of_baseline', 'mean'),
        mean_baseline_expression_TV_before=('baseline_expression_TV_before', 'mean'),
        mean_baseline_expression_TV_after=('baseline_expression_TV_after', 'mean'),
        mean_ambient_removal_profile_TV=('ambient_removal_profile_TV', 'mean'),
        mean_expression_distortion_reduction_percent=('expression_distortion_reduction_percent', 'mean')).reset_index()
    gsummary.to_csv(OUT / 'cell_type_ambient_restoration_summary.tsv', sep='\t', index=False)
    return genes, de, grouped, gsummary

if __name__ == '__main__':
    main()
