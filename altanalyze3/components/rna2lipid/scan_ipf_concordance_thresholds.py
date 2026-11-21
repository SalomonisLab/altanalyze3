"""Audit existing IPF results across explicitly displayed significance/fold gates.

No fitting, imputation, statistical testing, or BH recalculation is performed.
"""
from pathlib import Path
import hashlib
import itertools
import json

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

HERE = Path(__file__).resolve().parent
SOURCE = HERE / 'artifacts/LungMAP_full202_BH_corrected_20261006'
OUT = SOURCE / 'threshold_scan'
KEYS = ['contrast', 'unit', 'model', 'population']
CANDIDATE = 'candidate_native_log2_202'
PRIOR = 'current_explicit_log2'
CUTS = [.001, .005, .01, .025, .05, .1]
FOLDS = [1., 1.05, 1.1, 1.2, 1.3, 1.5, 2.]
MACROPHAGE_STATES = {'AM', 'AM-lipid', 'AM-prolif', 'MT+ AM', 'tMDM', 'IM'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def md(frame):
    return frame.to_markdown(index=False, floatfmt='.5g')


def main():
    paths = [SOURCE / n for n in ['all_measured_vs_imputed_tests.csv',
             'all_cellHarmony_tests.csv', 'lipid_matching_audit.csv',
             'population_replication_census.csv', 'BH_correction_manifest.json']]
    hashes = {str(p): sha(p) for p in paths}
    manifest = json.loads(paths[-1].read_text())
    for p in paths[:3]:
        expected_hash = manifest['sources'][p.name]['sha256'] if p.name == 'lipid_matching_audit.csv' else None
        if expected_hash and sha(p) != expected_hash:
            raise RuntimeError('Saved matching audit changed: ' + str(p))
    d = pd.read_csv(paths[0], float_precision='round_trip')
    full = pd.read_csv(paths[1], float_precision='round_trip')
    matching = pd.read_csv(paths[2])
    census = pd.read_csv(paths[3])
    contract = json.loads((HERE / 'integrity/baseline_contract.json').read_text())
    expected_genes = set(contract['lipids'])
    assert len(expected_genes) == 202
    assert set(full.model) == set(d.model) == {CANDIDATE, PRIOR}
    full_keys = set()
    bh_error = 0.
    for key, g in full.groupby(KEYS):
        assert len(g) == 202 and set(g.gene) == expected_genes, key
        assert np.isfinite(g[['pval', 'fdr', 'log2fc']].to_numpy()).all(), key
        assert g.pval.between(0, 1).all() and g.fdr.between(0, 1).all(), key
        err = np.max(np.abs(false_discovery_control(g.pval.to_numpy(), method='bh') - g.fdr.to_numpy()))
        assert err < 1e-12, (key, err)
        bh_error = max(bh_error, float(err))
        full_keys.add(key)
    assert not d.duplicated(KEYS + ['measured_feature']).any()
    assert d.BH_panel_size.eq(202).all()
    assert d.BH_policy.eq('full_panel_no_abundance_filter').all()
    assert np.isfinite(d[['pval', 'fdr', 'log2fc', 'signed_fold',
                          'supplied_raw_p', 'supplied_adjusted_p', 'supplied_log_effect']].to_numpy()).all()
    assert np.max(np.abs(np.log2(np.abs(d.signed_fold)) - np.abs(d.log2fc))) < 1e-12
    assert set(d.groupby(KEYS).groups) == full_keys
    expected_features = {}
    for model, m in matching.groupby('model'):
        assert len(m) == 202 and set(m.gene) == expected_genes
        hits = m[m.status.eq('matched')]
        assert not hits.measured_feature.duplicated().any()
        expected_features[model] = set(hits.measured_feature)
    assert expected_features[CANDIDATE] == expected_features[PRIOR]
    source_index = full.set_index(KEYS + ['gene'])
    joined_index = d.set_index(KEYS + ['gene'])
    for column in ['pval', 'fdr', 'log2fc', 'signed_fold']:
        assert np.array_equal(joined_index[column].to_numpy(),
                              source_index.loc[joined_index.index, column].to_numpy()), column
    rep = census[census.status.eq('tested')].set_index(KEYS)
    assert not rep.index.duplicated().any() and set(rep.index) == full_keys
    cross = pd.read_csv(HERE / 'artifacts/LungMAP_full202_native_log2_20261005/provenance/crosswalk.tsv', sep='\t')
    names = dict(zip(cross.cell_type, cross.cell_state))
    grid = pd.DataFrame(itertools.product(['rawp', 'BH'], CUTS, CUTS, FOLDS),
                        columns=['significance_kind', 'bulk_cutoff', 'predicted_cutoff', 'minimum_predicted_fold'])
    grid.insert(0, 'threshold_id', [f'{k}_B{b:g}_P{p:g}_F{f:g}' for k, b, p, f in grid.itertuples(index=False, name=None)])
    groups = {}
    results = []
    for key, g in d.groupby(KEYS, sort=True):
        assert set(g.measured_feature) == expected_features[key[2]], key
        g = g.sort_values('measured_feature').copy()
        groups[key] = g
        direction = np.sign(g.log2fc.to_numpy()) * np.sign(g.supplied_log_effect.to_numpy())
        magnitude = np.abs(g.signed_fold.to_numpy())
        masks = np.empty((len(grid), len(g)), dtype=bool)
        for kind, expcol, predcol in [('rawp', 'supplied_raw_p', 'pval'), ('BH', 'supplied_adjusted_p', 'fdr')]:
            ii = np.flatnonzero(grid.significance_kind.eq(kind).to_numpy())
            settings = grid.iloc[ii]
            masks[ii] = ((g[expcol].to_numpy()[None, :] < settings.bulk_cutoff.to_numpy()[:, None]) &
                         (g[predcol].to_numpy()[None, :] < settings.predicted_cutoff.to_numpy()[:, None]) &
                         (magnitude[None, :] >= settings.minimum_predicted_fold.to_numpy()[:, None]))
        c = (masks & (direction[None, :] > 0)).sum(axis=1)
        o = (masks & (direction[None, :] < 0)).sum(axis=1)
        z = (masks & (direction[None, :] == 0)).sum(axis=1)
        n = masks.sum(axis=1)
        assert np.array_equal(n, c + o + z)
        counts = grid.copy()
        for field, value in zip(KEYS, key):
            counts[field] = value
        counts['cell_state_full'] = names.get(key[3], key[3])
        counts['comparison_kind'] = 'state_override' if '__vs__' in key[3] else 'same_state'
        counts['concordant'] = c
        counts['opposite'] = o
        counts['zero_or_ambiguous'] = z
        counts['both_significant_after_fold'] = n
        counts['concordant_up'] = (masks & (direction[None, :] > 0) & (g.supplied_log_effect.to_numpy()[None, :] > 0)).sum(axis=1)
        counts['concordant_down'] = c - counts.concordant_up
        counts['agreement_percent'] = np.divide(100. * c, c + o, out=np.full(len(c), np.nan), where=(c + o) > 0)
        # Integer comparison avoids rounding a 75.0% result into strict >75%.
        counts['above_75_percent'] = 4 * c > 3 * (c + o)
        counts['directional_denominator'] = c + o
        counts['small_denominator_lt10'] = (c + o) < 10
        for column in ['n_case', 'n_control', 'n_case_donors', 'n_control_donors']:
            counts[column] = rep.loc[key, column]
        results.append(counts)
    all_rows = pd.concat(results, ignore_index=True)
    assert len(all_rows) == len(full_keys) * len(grid)
    paired_keys = ['contrast', 'unit', 'population', 'threshold_id']
    candidate = all_rows[all_rows.model.eq(CANDIDATE)].copy()
    prior = all_rows[all_rows.model.eq(PRIOR)].copy()
    assert set(candidate.set_index(paired_keys).index) == set(prior.set_index(paired_keys).index)
    paired = candidate.merge(prior[paired_keys + ['concordant', 'opposite', 'agreement_percent', 'above_75_percent']],
                             on=paired_keys, how='left', validate='one_to_one', suffixes=('', '_prior'))
    assert paired.concordant_prior.notna().all()
    paired['concordant_count_difference'] = paired.concordant - paired.concordant_prior
    paired['agreement_percentage_point_difference'] = paired.agreement_percent - paired.agreement_percent_prior
    winners = all_rows[all_rows.above_75_percent].copy()
    winners = winners.sort_values(['model', 'concordant', 'agreement_percent', 'opposite'], ascending=[True, False, False, True])
    ids_c, ids_o, ids_z = [], [], []
    for row in winners.itertuples(index=False):
        key = tuple(getattr(row, k) for k in KEYS)
        g = groups[key]
        expcol, predcol = ('supplied_raw_p', 'pval') if row.significance_kind == 'rawp' else ('supplied_adjusted_p', 'fdr')
        mask = g[expcol].lt(row.bulk_cutoff) & g[predcol].lt(row.predicted_cutoff) & g.signed_fold.abs().ge(row.minimum_predicted_fold)
        relation = np.sign(g.log2fc) * np.sign(g.supplied_log_effect)
        concordant = g.loc[mask & relation.gt(0), 'measured_feature'].tolist()
        opposite = g.loc[mask & relation.lt(0), 'measured_feature'].tolist()
        ambiguous = g.loc[mask & relation.eq(0), 'measured_feature'].tolist()
        assert len(concordant) == row.concordant and len(opposite) == row.opposite
        ids_c.append('; '.join(concordant)); ids_o.append('; '.join(opposite)); ids_z.append('; '.join(ambiguous))
    winners['concordant_lipids'] = ids_c
    winners['opposite_lipids'] = ids_o
    winners['zero_or_ambiguous_lipids'] = ids_z
    # These are display rankings only; every grid row is retained in all_rows.
    sort_columns = ['concordant', 'agreement_percent', 'opposite', 'bulk_cutoff', 'predicted_cutoff', 'minimum_predicted_fold', 'population']
    ascending = [False, False, True, True, True, False, True]
    primary = winners[winners.comparison_kind.eq('same_state')].sort_values(sort_columns, ascending=ascending)
    equal = primary[primary.bulk_cutoff.eq(primary.predicted_cutoff)]
    best = equal.drop_duplicates(['contrast', 'unit', 'model', 'significance_kind'])
    best_independent = primary.drop_duplicates(['contrast', 'unit', 'model', 'significance_kind'])
    best_state = equal.drop_duplicates(KEYS + ['significance_kind'])
    # MT+ AM is spelled "Metallothionein+ AM" in the crosswalk. Include its
    # established identity explicitly rather than relying on a substring.
    is_macrophage = equal.population.isin(MACROPHAGE_STATES) | equal.cell_state_full.str.contains('macrophage', case=False)
    macrophages = equal[is_macrophage].drop_duplicates(['contrast', 'unit', 'model', 'significance_kind'])
    # Verify the existing BH/no-fold summary at both previously reported cutoffs.
    earlier = pd.read_csv(SOURCE / 'cell_state_FDR_rankings.csv')
    verified = 0
    for cut in [.05, .1]:
        a = all_rows[all_rows.significance_kind.eq('BH') & all_rows.bulk_cutoff.eq(cut) &
                     all_rows.predicted_cutoff.eq(cut) & all_rows.minimum_predicted_fold.eq(1)].set_index(KEYS)
        e = earlier[earlier.bulk_BH_cutoff.eq(cut) & earlier.predicted_BH_cutoff.eq(cut)].set_index(KEYS)
        assert set(a.index) == set(e.index)
        e = e.loc[a.index]
        for field in ['concordant', 'opposite', 'zero_or_ambiguous']:
            assert np.array_equal(a[field].to_numpy(), e[field].to_numpy()), field
        verified += len(a)
    assert all(sha(p) == hashes[str(p)] for p in paths)
    OUT.mkdir(exist_ok=True)
    all_rows.to_csv(OUT / 'all_threshold_results.csv', index=False)
    winners.to_csv(OUT / 'all_above75_comparisons_with_lipids.csv', index=False)
    paired.to_csv(OUT / 'candidate_vs_prior_same_thresholds.csv', index=False)
    best.to_csv(OUT / 'best_per_cohort_equal_cutoffs.csv', index=False)
    best_independent.to_csv(OUT / 'best_per_cohort_independent_cutoffs.csv', index=False)
    best_state.to_csv(OUT / 'best_per_cell_state.csv', index=False)
    macrophages.to_csv(OUT / 'best_macrophage_per_cohort.csv', index=False)
    census.to_csv(OUT / 'all_comparison_replication_status.csv', index=False)
    fields = ['contrast', 'unit', 'population', 'significance_kind', 'bulk_cutoff', 'predicted_cutoff',
              'minimum_predicted_fold', 'concordant_up', 'concordant_down', 'concordant', 'opposite',
              'agreement_percent', 'directional_denominator', 'n_case_donors', 'n_control_donors']
    cn = best[best.model.eq(CANDIDATE)]
    cm = macrophages[macrophages.model.eq(CANDIDATE)]
    lines = ['# IPF concordance above 75%: fold and significance threshold scan', '',
             'The full 202-lipid native-log2 ElasticNetCV candidate and prior ElasticNet model are evaluated '
             'using their saved cellHarmony tests and the supplied lung-tissue IPF lipidomics statistics. '
             'Models, features, samples, predictions, raw tests and BH values are unchanged.', '',
             '**Threshold grid:** predicted fold magnitude ≥ ' + ', '.join(f'{v:g}' for v in FOLDS) +
             '; raw p or BH-adjusted p < ' + ', '.join(f'{v:g}' for v in CUTS) +
             '. Experimental and predicted cutoffs are varied independently within the same significance kind. '
             'There are 504 threshold settings per tested cohort/state/input/model group.', '',
             '**Definition:** agreement = concordant/(concordant + opposite) among results significant in both '
             'and passing the predicted fold threshold. The inequality is strictly >75%; exactly 75% is not included. '
             'Rounded-zero experimental effects are counted separately as directionally ambiguous. '
             'A threshold of 1 imposes no additional fold-magnitude cutoff. '
             'Signed +2 means twofold higher; -2 means twofold lower. Experimental fold magnitudes are not thresholded '
             'because the supplied effect log base is unconfirmed.', '',
             '**Ranking:** largest concordant lipid count among qualifying settings, then greater percentage agreement, '
             'fewer opposite lipids, smaller p cutoffs, greater fold cutoff, and cell-state name. '
             'No minimum lipid-count restriction is imposed. Every qualifying setting, including one-lipid results, '
             'is in the complete CSV. Main tables use equal experimental/predicted significance cutoffs; '
             'independent-cutoff results are separately retained.', '',
             'These are agreements selected by a threshold scan, not independent validation of the selected thresholds. '
             'BH remains over all 202 outputs within each original comparison; it does not correct the search across thresholds '
             'or states. Metacell p-values use metacell observations and are nominal relative to donor-level replication. '
             'Donor counts are displayed. Insufficient-replication groups remain in the census and have no fabricated p-values. '
             'Natri pooled/less/more-fibrotic contrasts share donors and are not independent cohorts.', '',
             '## Best candidate comparisons at equal p cutoffs', '', md(cn[fields].sort_values(['significance_kind', 'contrast', 'unit'])), '',
             'A missing cohort/unit in this table means no tested same-state setting exceeded 75%, '
             'or the group could not be tested; see the complete scan and replication census.', '',
             '## Best candidate macrophage comparisons at equal p cutoffs', '', md(cm[fields].sort_values(['significance_kind', 'contrast', 'unit'])), '',
             '## Identical-threshold comparison with the prior model', '']
    paired_best = cn.merge(paired[paired_keys + ['concordant_prior', 'opposite_prior', 'agreement_percent_prior']],
                           on=paired_keys, how='left', validate='one_to_one')
    lines += [md(paired_best[fields + ['concordant_prior', 'opposite_prior', 'agreement_percent_prior']]), '',
              '## Best candidate comparisons with independent experimental and predicted cutoffs', '',
              md(best_independent[best_independent.model.eq(CANDIDATE)][fields].sort_values(['significance_kind', 'contrast', 'unit'])), '',
              '## All qualifying cell states: best equal-cutoff setting per state and significance kind', '',
              md(best_state[best_state.model.eq(CANDIDATE)][fields].sort_values(['contrast', 'unit', 'population', 'significance_kind'])), '',
              '## Specific lipids for each best candidate comparison', '']
    lipid_detail = []
    focus = pd.concat([cn, cm, best_independent[best_independent.model.eq(CANDIDATE)]], ignore_index=True).drop_duplicates(KEYS + ['threshold_id'])
    for row in focus.sort_values(['contrast', 'unit', 'population', 'significance_kind']).itertuples(index=False):
        g = groups[tuple(getattr(row, k) for k in KEYS)].copy()
        expcol, predcol = ('supplied_raw_p', 'pval') if row.significance_kind == 'rawp' else ('supplied_adjusted_p', 'fdr')
        selected = g[(g[expcol] < row.bulk_cutoff) & (g[predcol] < row.predicted_cutoff) &
                     (g.signed_fold.abs() >= row.minimum_predicted_fold)].copy()
        relation = np.sign(selected.log2fc) * np.sign(selected.supplied_log_effect)
        selected['direction_relation'] = np.where(relation > 0, 'concordant', np.where(relation < 0, 'opposite', 'zero_or_ambiguous'))
        selected['threshold_id'] = row.threshold_id
        lipid_detail.append(selected)
        lines += ['', f'### {row.contrast} / {row.unit} / {row.population}: {row.significance_kind}<{row.bulk_cutoff:g} '
                  f'in both, fold≥{row.minimum_predicted_fold:g}', '',
                  md(selected[['measured_feature', 'gene', 'signed_fold', 'pval', 'fdr', 'supplied_raw_p',
                               'supplied_adjusted_p', 'direction_relation']])]
    pd.concat(lipid_detail, ignore_index=True).to_csv(OUT / 'best_comparison_specific_lipids.csv', index=False)
    lines += ['', '## Complete results', '',
              '- [Every >75% setting, both models, with concordant and opposite lipid identities](all_above75_comparisons_with_lipids.csv)',
              '- [Every scanned setting, including failures and empty intersections](all_threshold_results.csv)',
              '- [Both models at identical thresholds](candidate_vs_prior_same_thresholds.csv)',
              '- [Best cohort results allowing different experimental/predicted cutoffs](best_per_cohort_independent_cutoffs.csv)',
              '- [All tested and insufficient-replication comparison identities](all_comparison_replication_status.csv)', '',
              f'Computed {len(all_rows):,} model/comparison/threshold rows; {len(winners):,} exceed 75%. '
              f'Candidate qualifying rows: {int(winners.model.eq(CANDIDATE).sum()):,}; '
              f'prior qualifying rows: {int(winners.model.eq(PRIOR).sum()):,}. '
              'Same-state and explicit state-override rows are labeled separately in the complete outputs.']
    (OUT / 'THRESHOLD_SCAN_REPORT.md').write_text('\n'.join(lines) + '\n')
    audit = {'source_sha256': hashes, 'script_sha256': sha(Path(__file__)),
             'p_cutoffs': CUTS, 'predicted_fold_cutoffs': FOLDS, 'settings_per_group': len(grid),
             'full202_groups_verified': len(full_keys), 'matched_experimental_features_per_model': len(expected_features[CANDIDATE]),
             'BH_values_independently_verified': len(full), 'BH_max_error': bh_error,
             'existing_no_fold_count_rows_reproduced': verified, 'scan_rows': len(all_rows),
             'above75_rows': len(winners), 'above75_by_model': winners.groupby('model').size().to_dict(),
             'new_fitting_inference_tests_or_BH_recalculation': False,
             'all_original_comparison_and_model_identities_preserved': True,
             'source_hashes_unchanged': all(sha(p) == hashes[str(p)] for p in paths)}
    (OUT / 'threshold_scan_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(audit, indent=2))
    print('\nBEST CANDIDATE BH COMPARISONS:\n' + cn[cn.significance_kind.eq('BH')][fields].to_string(index=False))


if __name__ == '__main__':
    main()
