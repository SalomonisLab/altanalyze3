"""Describe corrected IPF concordance and every significant opposite-direction result."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / 'artifacts/LungMAP_full202_BH_corrected_20261006'
CANDIDATE = 'candidate_native_log2_202'
KEYS = ['contrast', 'unit', 'model', 'population']
CUTS = [(b, p) for b in [.05, .1] for p in [.05, .1]]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    source = OUT / 'all_measured_vs_imputed_tests.csv'
    d = pd.read_csv(source, float_precision='round_trip')
    full = pd.read_csv(OUT / 'all_cellHarmony_tests.csv', float_precision='round_trip')
    contract = json.loads((HERE / 'integrity/baseline_contract.json').read_text())
    expected = set(contract['lipids'])
    for key, group in full.groupby(KEYS):
        if len(group) != 202 or set(group.gene) != expected or not np.isfinite(group.fdr).all():
            raise RuntimeError('Incomplete corrected input: ' + str(key))
    if d.duplicated(KEYS + ['measured_feature']).any():
        raise RuntimeError('Duplicate measured lipid identity.')
    if not d.BH_panel_size.eq(202).all() or not d.BH_policy.eq('full_panel_no_abundance_filter').all():
        raise RuntimeError('The input does not use approved complete-panel BH.')
    census = pd.read_csv(OUT / 'population_replication_census.csv')
    measured = pd.read_csv(OUT / 'experimental_lung_IPF_lipids_all_544.csv', index_col=0)
    matched = set(d.measured_feature)
    cross = pd.read_csv(HERE / 'artifacts/LungMAP_full202_native_log2_20261005/provenance/crosswalk.tsv', sep='\t')
    names = dict(zip(cross.cell_type, cross.cell_state))
    rows, detailed, unions = [], [], []
    for bulk_cut, predicted_cut in CUTS:
        selected = d[d.supplied_adjusted_p.lt(bulk_cut) & d.fdr.lt(predicted_cut)].copy()
        selected['bulk_BH_cutoff'] = bulk_cut
        selected['predicted_BH_cutoff'] = predicted_cut
        product = selected.log2fc * selected.supplied_log_effect
        selected['direction_relation'] = np.where(product.gt(0), 'concordant', np.where(product.lt(0), 'opposite', 'zero_or_ambiguous'))
        selected['bulk_direction'] = np.where(selected.supplied_log_effect.gt(0), 'up', np.where(selected.supplied_log_effect.lt(0), 'down', 'zero'))
        selected['predicted_direction'] = np.where(selected.log2fc.gt(0), 'up', np.where(selected.log2fc.lt(0), 'down', 'zero'))
        selected['comparison_kind'] = np.where(selected.population.str.contains('__vs__'), 'state_override', 'same_state')
        detailed.append(selected)
        count = int(measured.index.isin(matched).astype(int).dot(measured.supplied_adjusted_p.lt(bulk_cut).astype(int)))
        grouped = {key: group for key, group in selected.groupby(KEYS)}
        for key, source_group in d.groupby(KEYS):
            contrast, unit, model, population = key
            g = grouped.get(key, selected.iloc[:0])
            c = g[g.direction_relation.eq('concordant')]
            opposite = g[g.direction_relation.eq('opposite')]
            n = g[g.direction_relation.eq('zero_or_ambiguous')]
            rep = census[(census.contrast == contrast) & (census.unit == unit) &
                         (census.model == model) & (census.population == population)]
            if len(rep) != 1:
                raise RuntimeError('Unresolved replication metadata: ' + str(key))
            r = rep.iloc[0]
            rows.append(dict(zip(KEYS, key)) | {
                'cell_state_full': names.get(population, population),
                'comparison_kind': 'state_override' if '__vs__' in population else 'same_state',
                'bulk_BH_cutoff': bulk_cut, 'predicted_BH_cutoff': predicted_cut,
                'matched_bulk_significant': count, 'both_significant': g.measured_feature.nunique(),
                'concordant': c.measured_feature.nunique(), 'opposite': opposite.measured_feature.nunique(),
                'zero_or_ambiguous': n.measured_feature.nunique(),
                'concordant_up': c.bulk_direction.eq('up').sum(), 'concordant_down': c.bulk_direction.eq('down').sum(),
                'agreement_percent': 100 * len(c) / (len(c) + len(opposite)) if len(c) + len(opposite) else np.nan,
                'n_case': r.n_case, 'n_control': r.n_control,
                'n_case_donors': r.n_case_donors, 'n_control_donors': r.n_control_donors})
        primary = selected[selected.comparison_kind.eq('same_state')]
        for key, g in primary.groupby(['contrast', 'unit', 'model']):
            concordant = set(g.loc[g.direction_relation.eq('concordant'), 'measured_feature'])
            opposite = set(g.loc[g.direction_relation.eq('opposite'), 'measured_feature'])
            unions.append(dict(zip(['contrast', 'unit', 'model'], key)) | {
                'bulk_BH_cutoff': bulk_cut, 'predicted_BH_cutoff': predicted_cut,
                'concordant_any_state': len(concordant), 'opposite_any_state': len(opposite),
                'both_directions_in_different_states': len(concordant & opposite),
                'opposite_without_concordant_state': len(opposite - concordant)})
    ranking = pd.DataFrame(rows)
    assert (ranking.both_significant == ranking.concordant + ranking.opposite + ranking.zero_or_ambiguous).all()
    ranking = ranking.sort_values(['concordant', 'opposite', 'population'], ascending=[False, True, True])
    grouping = ['contrast', 'unit', 'model', 'bulk_BH_cutoff', 'predicted_BH_cutoff']
    ranking['coverage_rank'] = ranking.groupby(grouping + ['comparison_kind']).cumcount() + 1
    best = ranking[ranking.comparison_kind.eq('same_state') & ranking.concordant.gt(0)].drop_duplicates(grouping)
    significant = pd.concat(detailed, ignore_index=True)
    discordant = significant[significant.direction_relation.eq('opposite')].copy()
    discordant = discordant.sort_values(grouping + ['population', 'gene'])
    unions = pd.DataFrame(unions)
    ranking.to_csv(OUT / 'cell_state_FDR_rankings.csv', index=False)
    best.to_csv(OUT / 'best_cell_states_by_FDR.csv', index=False)
    significant.to_csv(OUT / 'both_significant_cell_state_results.csv', index=False)
    discordant.to_csv(OUT / 'discordant_lipids_all_cell_states.csv', index=False)
    unions.to_csv(OUT / 'concordant_and_opposite_state_unions.csv', index=False)
    # Validate the equal-cutoff union counts against the existing corrected report.
    existing = pd.read_csv(OUT / 'cohort_concordant_union_summary.csv')
    for cutoff, gate in [(.05, 'both_BH005'), (.1, 'both_BH010')]:
        expected_counts = existing[existing.overlap_gate.eq(gate) & existing.minimum_imputed_fold.eq(1)].set_index(['contrast', 'unit', 'model'])
        for r in unions[unions.bulk_BH_cutoff.eq(cutoff) & unions.predicted_BH_cutoff.eq(cutoff)].itertuples():
            assert r.concordant_any_state == expected_counts.loc[(r.contrast, r.unit, r.model), 'concordant_union']
    focus = best[best.model.eq(CANDIDATE) & best.bulk_BH_cutoff.eq(best.predicted_BH_cutoff)]
    fields = ['contrast', 'unit', 'bulk_BH_cutoff', 'predicted_BH_cutoff', 'population', 'concordant', 'opposite',
              'agreement_percent', 'n_case_donors', 'n_control_donors']
    def table(frame):
        return frame.to_markdown(index=False, floatfmt='.4g')
    lines = ['# IPF cell-state concordance and significant opposite-direction lipids', '',
             'Reference: Dr. Clair’s supplied lung-tissue `10_results_with_statistics.csv`. '
             'A publication citation is not present in the supplied provenance; none is inferred. '
             'The RNA cohorts and experimental lipidomics are unmatched.', '',
             'Primary results use the full 202-lipid corrected candidate. BH is applied to all 202 outputs within each existing '
             'cell-state/cohort/input-unit comparison. The experimental gate uses `supplied_adjusted_p`; '
             'the prediction gate uses corrected `fdr`. The log-effect product must be positive for concordance or negative '
             'for an opposite direction; a zero effect is explicitly recorded separately. No additional fold threshold is applied.', '',
             'Best means the largest number of distinct concordant experimental features, with fewer opposite features breaking ties. '
             'Agreement is concordant/(concordant + opposite) among results significant in both. Selection across cell states and '
             'FDR cutoffs is descriptive and does not supply a new independent validation or FDR control for the union.', '',
             '## Best candidate cell states at equal FDR cutoffs', '', table(focus[fields].sort_values(['bulk_BH_cutoff', 'contrast', 'unit'])), '',
             'MC = metacells; PB = pseudobulks. AM = alveolar macrophage; AM-lipid = lipid-homeostatic alveolar macrophage; '
             'MT+ AM = metallothionein-positive alveolar macrophage; AT2 = alveolar type 2; AT2-prolif = proliferating AT2; '
             'Ciliated-Bronch = bronchiolar multiciliated; Ciliated-axon = axonemal multiciliated.', '',
             'Metacell statistical tests use metacells, rather than independent donors, as observations. '
             'The donor counts describe source coverage and do not change that test. UPenn v3 has only one IPF donor in its best state; '
             'it does not establish replicated disease effects. UPenn PB comparisons lack the required replicate coverage. '
             'Jaiswal PB comparisons were tested but have no concordant overlap at the two equal BH thresholds.', '',
             '## All FDR combinations and prior-model comparison', '',
             'Bulk BH and predicted BH are independently shown at 0.05 and 0.1. Larger cutoffs admit more results; '
             'counts with different bulk cutoffs have different experimental denominators (85 at 0.05, 97 at 0.1).', '',
             table(best[fields + ['model']].sort_values(grouping)), '',
             '## Discordant lipids in each best candidate state', '',
             'Every row below is significant in both datasets and has the opposite direction. '
             'Signed predicted folds use +2 for a twofold increase and -2 for a twofold decrease. '
             'Experimental direction comes from the supplied effect; its log base is unconfirmed, so experimental fold magnitudes are not inferred.', '']
    primary_best = focus.sort_values(['bulk_BH_cutoff', 'contrast', 'unit'])
    for r in primary_best.itertuples():
        g = discordant[(discordant.model == CANDIDATE) & (discordant.contrast == r.contrast) &
                       (discordant.unit == r.unit) & (discordant.population == r.population) &
                       (discordant.bulk_BH_cutoff == r.bulk_BH_cutoff) &
                       (discordant.predicted_BH_cutoff == r.predicted_BH_cutoff)]
        assert g.measured_feature.nunique() == r.opposite
        lines += ['', f'### {r.contrast} / {r.unit} / {r.population}; bulk BH<{r.bulk_BH_cutoff}, predicted BH<{r.predicted_BH_cutoff}', '',
                  table(g[['gene', 'measured_feature', 'bulk_direction', 'predicted_direction', 'signed_fold',
                           'pval', 'fdr', 'supplied_raw_p', 'supplied_adjusted_p']])]
    lines += ['', '## Opposite-direction lipids across any candidate cell state, by comparison', '',
              'A feature may be concordant in one state and opposite in another. These sets overlap; neither is silently excluded. '
              'Cell-state-specific biological effects remain hypotheses because the experimental reference is bulk tissue.', '',
              table(unions[unions.model.eq(CANDIDATE) & unions.bulk_BH_cutoff.eq(unions.predicted_BH_cutoff)].sort_values(['bulk_BH_cutoff','contrast','unit']))]
    for cutoff in [.05, .1]:
        primary_discordant = discordant[discordant.model.eq(CANDIDATE) & discordant.bulk_BH_cutoff.eq(cutoff) &
                                      discordant.predicted_BH_cutoff.eq(cutoff) & discordant.comparison_kind.eq('same_state')]
        primary_significant = significant[significant.model.eq(CANDIDATE) & significant.bulk_BH_cutoff.eq(cutoff) &
                                        significant.predicted_BH_cutoff.eq(cutoff) & significant.comparison_kind.eq('same_state')]
        for (contrast, unit), group in primary_discordant.groupby(['contrast', 'unit']):
            lipid_rows = []
            for feature, g in group.groupby('measured_feature'):
                agree = primary_significant[(primary_significant.contrast == contrast) & (primary_significant.unit == unit) &
                                             (primary_significant.measured_feature == feature) &
                                             primary_significant.direction_relation.eq('concordant')]
                lipid_rows.append({'lipid': g.gene.iloc[0], 'measured_feature': feature,
                                   'bulk_direction': g.bulk_direction.iloc[0], 'predicted_opposite_direction': g.predicted_direction.iloc[0],
                                   'opposite_cell_states': '; '.join(sorted(g.population.unique())),
                                   'concordant_cell_states': '; '.join(sorted(agree.population.unique())),
                                   'bulk_BH': g.supplied_adjusted_p.iloc[0], 'lowest_opposite_predicted_BH': g.fdr.min()})
            lines += ['', f'### {contrast} / {unit}; both BH<{cutoff}; every discordant lipid and its states', '', table(pd.DataFrame(lipid_rows))]
    lines += ['', '## Complete state-level rows and override comparisons', '',
              '`discordant_lipids_all_cell_states.csv` contains every opposite result for both models and all four FDR combinations, '
              'including explicit disease-state/reference-state overrides. Same-state and override comparisons are labeled separately. '
              '`cell_state_FDR_rankings.csv` contains every tested state, including states with zero concordant overlap. '
              '`population_replication_census.csv` retains the untested insufficient-replication states. '
              'Model features, samples, predictions, raw tests and the BH correction were not changed by this summary.']
    report = OUT / 'CELL_STATE_CONCORDANCE_AND_DISCORDANCE.md'
    report.write_text('\n'.join(lines) + '\n')
    audit = {'source': str(source), 'source_sha256': sha(source), 'source_rows': len(d),
             'matched_experimental_features': len(matched), 'full_model_test_groups_verified': full.groupby(KEYS).ngroups,
             'FDR_combinations': CUTS, 'ranked_state_rows': len(ranking), 'opposite_rows_all_cutoffs_models': len(discordant),
             'equal_cutoff_union_counts_match_existing_report': True,
             'new_training_inference_tests_or_BH_performed': False, 'report': str(report)}
    (OUT / 'cell_state_summary_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
