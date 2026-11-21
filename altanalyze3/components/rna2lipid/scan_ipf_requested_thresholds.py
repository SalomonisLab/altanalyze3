"""Report strict >1.1/>1.2/>1.5 folds at the user's three p-value gates."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

HERE = Path(__file__).resolve().parent
SOURCE = HERE / 'artifacts/LungMAP_full202_BH_corrected_20261006'
OUT = SOURCE / 'requested_thresholds'
KEYS = ['contrast', 'unit', 'model', 'population']
CANDIDATE = 'candidate_native_log2_202'
PRIOR = 'current_explicit_log2'
FOLDS = [1.1, 1.2, 1.5]
GATES = [('rawp005', 'supplied_raw_p', 'pval', .05),
         ('BH005', 'supplied_adjusted_p', 'fdr', .05),
         ('BH010', 'supplied_adjusted_p', 'fdr', .1)]
MACROPHAGES = {'AM', 'AM-lipid', 'AM-prolif', 'MT+ AM', 'tMDM', 'IM'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    paths = [SOURCE / n for n in ['all_measured_vs_imputed_tests.csv',
             'all_cellHarmony_tests.csv', 'population_replication_census.csv', 'lipid_matching_audit.csv']]
    hashes = {str(p): sha(p) for p in paths}
    d = pd.read_csv(paths[0], float_precision='round_trip')
    full = pd.read_csv(paths[1], float_precision='round_trip')
    rep = pd.read_csv(paths[2])
    mapping = pd.read_csv(paths[3])
    expected = set(json.loads((HERE / 'integrity/baseline_contract.json').read_text())['lipids'])
    assert len(expected) == 202 and set(full.model) == set(d.model) == {CANDIDATE, PRIOR}
    verified = 0
    for key, g in full.groupby(KEYS):
        assert len(g) == 202 and set(g.gene) == expected, key
        assert np.isfinite(g[['pval', 'fdr', 'log2fc']].to_numpy()).all()
        assert np.max(np.abs(false_discovery_control(g.pval.to_numpy(), method='bh') - g.fdr.to_numpy())) < 1e-12
        verified += 1
    assert not d.duplicated(KEYS + ['measured_feature']).any()
    assert d.BH_panel_size.eq(202).all() and d.BH_policy.eq('full_panel_no_abundance_filter').all()
    assert set(d.groupby(KEYS).groups) == set(full.groupby(KEYS).groups)
    assert np.isfinite(d[['pval', 'fdr', 'signed_fold', 'log2fc', 'supplied_log_effect',
                          'supplied_raw_p', 'supplied_adjusted_p']].to_numpy()).all()
    assert np.max(np.abs(np.log2(d.signed_fold.abs()) - d.log2fc.abs())) < 1e-12
    full_index = full.set_index(KEYS + ['gene'])
    ix = d.set_index(KEYS + ['gene'])
    for col in ['pval', 'fdr', 'log2fc', 'signed_fold']:
        assert np.array_equal(ix[col].to_numpy(), full_index.loc[ix.index, col].to_numpy())
    features = {model: set(g.loc[g.status.eq('matched'), 'measured_feature']) for model, g in mapping.groupby('model')}
    assert features[CANDIDATE] == features[PRIOR]
    census = rep[rep.status.eq('tested')].set_index(KEYS)
    assert not census.index.duplicated().any() and set(census.index) == set(full.groupby(KEYS).groups)
    names_df = pd.read_csv(HERE / 'artifacts/LungMAP_full202_native_log2_20261005/provenance/crosswalk.tsv', sep='\t')
    names = dict(zip(names_df.cell_type, names_df.cell_state))
    rows, details = [], []
    for key, g in d.groupby(KEYS, sort=True):
        assert set(g.measured_feature) == features[key[2]], key
        for gate, expcol, predcol, cutoff in GATES:
            for fold in FOLDS:
                s = g[g[expcol].lt(cutoff) & g[predcol].lt(cutoff) & g.signed_fold.abs().gt(fold)].copy()
                relation = np.sign(s.log2fc) * np.sign(s.supplied_log_effect)
                c = s[relation.gt(0)]; o = s[relation.lt(0)]; z = s[relation.eq(0)]
                count = len(c); opposite = len(o); denominator = count + opposite
                row = dict(zip(KEYS, key)) | {
                    'cell_state_full': names.get(key[3], key[3]),
                    'comparison_kind': 'state_override' if '__vs__' in key[3] else 'same_state',
                    'gate': gate, 'p_cutoff_both': cutoff, 'predicted_fold_strictly_greater_than': fold,
                    'concordant_up': int(c.supplied_log_effect.gt(0).sum()),
                    'concordant_down': int(c.supplied_log_effect.lt(0).sum()),
                    'concordant': count, 'opposite': opposite, 'zero_or_ambiguous': len(z),
                    'both_significant_after_fold': len(s), 'directional_denominator': denominator,
                    'agreement_percent': 100 * count / denominator if denominator else np.nan,
                    'above_75_percent': 4 * count > 3 * denominator,
                    'concordant_lipids': '; '.join(sorted(c.measured_feature)),
                    'opposite_lipids': '; '.join(sorted(o.measured_feature))}
                for col in ['n_case', 'n_control', 'n_case_donors', 'n_control_donors']:
                    row[col] = census.loc[key, col]
                assert row['both_significant_after_fold'] == count + opposite + len(z)
                rows.append(row)
                if row['above_75_percent']:
                    s['direction_relation'] = np.where(relation.gt(0), 'concordant', np.where(relation.lt(0), 'opposite', 'zero_or_ambiguous'))
                    s['gate'] = gate; s['p_cutoff_both'] = cutoff; s['predicted_fold_strictly_greater_than'] = fold
                    details.append(s)
    all_rows = pd.DataFrame(rows)
    assert len(all_rows) == verified * 9
    rank = ['concordant', 'agreement_percent', 'opposite', 'population']
    qualify = all_rows[all_rows.above_75_percent].sort_values(rank, ascending=[False, False, True, True])
    candidate = qualify[qualify.model.eq(CANDIDATE) & qualify.comparison_kind.eq('same_state')]
    best = candidate.drop_duplicates(['contrast', 'unit', 'gate', 'predicted_fold_strictly_greater_than'])
    strongest = candidate.drop_duplicates(['gate', 'predicted_fold_strictly_greater_than'])
    strongest_multiple_donors = candidate[candidate.n_case_donors.ge(2) & candidate.n_control_donors.ge(2)].drop_duplicates(['gate', 'predicted_fold_strictly_greater_than'])
    macro = candidate[candidate.population.isin(MACROPHAGES)].drop_duplicates(['contrast', 'unit', 'gate', 'predicted_fold_strictly_greater_than'])
    pairkeys = ['contrast', 'unit', 'population', 'gate', 'predicted_fold_strictly_greater_than']
    c = all_rows[all_rows.model.eq(CANDIDATE)]; p = all_rows[all_rows.model.eq(PRIOR)]
    assert set(c.set_index(pairkeys).index) == set(p.set_index(pairkeys).index)
    paired = c.merge(p[pairkeys + ['concordant', 'opposite', 'agreement_percent', 'above_75_percent']],
                     on=pairkeys, how='left', validate='one_to_one', suffixes=('', '_prior'))
    assert paired.concordant_prior.notna().all()
    # Independently count every report winner using scalar row comparisons.
    records = {key: g.to_dict('records') for key, g in d.groupby(KEYS)}
    manual_checks = 0
    for r in best.to_dict('records'):
        expcol, predcol = ('supplied_raw_p', 'pval') if r['gate']=='rawp005' else ('supplied_adjusted_p', 'fdr')
        nc = no = nz = 0
        for lipid in records[tuple(r[k] for k in KEYS)]:
            if lipid[expcol] >= r['p_cutoff_both'] or lipid[predcol] >= r['p_cutoff_both'] or abs(lipid['signed_fold']) <= r['predicted_fold_strictly_greater_than']:
                continue
            if lipid['log2fc'] == 0 or lipid['supplied_log_effect'] == 0: nz += 1
            elif (lipid['log2fc'] > 0) == (lipid['supplied_log_effect'] > 0): nc += 1
            else: no += 1
        assert (nc,no,nz) == (r['concordant'],r['opposite'],r['zero_or_ambiguous'])
        manual_checks += 1
    assert all(sha(p) == hashes[str(p)] for p in paths)
    OUT.mkdir(exist_ok=True)
    all_rows.to_csv(OUT / 'all_requested_threshold_results.csv', index=False)
    qualify.to_csv(OUT / 'all_above75_with_lipid_identities.csv', index=False)
    best.to_csv(OUT / 'best_candidate_by_cohort_and_threshold.csv', index=False)
    strongest_multiple_donors.to_csv(OUT / 'strongest_candidate_multiple_donor_arms.csv', index=False)
    macro.to_csv(OUT / 'best_candidate_macrophages.csv', index=False)
    paired.to_csv(OUT / 'candidate_vs_prior_identical_thresholds.csv', index=False)
    rep.to_csv(OUT / 'complete_replication_census.csv', index=False)
    fields = ['gate', 'predicted_fold_strictly_greater_than', 'contrast', 'unit', 'population',
              'concordant_up', 'concordant_down', 'concordant', 'opposite', 'agreement_percent',
              'n_case_donors', 'n_control_donors']
    counts = candidate.groupby(['gate','predicted_fold_strictly_greater_than']).size().rename('qualifying_same_state_comparisons').reset_index()
    counts.to_csv(OUT / 'qualifying_comparison_counts.csv', index=False)
    def table(frame): return frame.to_markdown(index=False, floatfmt='.5g')
    lines = ['# Exact requested IPF fold and p-value thresholds', '',
             'Strict predicted fold magnitude >1.1, >1.2, or >1.5; raw p<0.05 in both datasets, '
             'BH-adjusted p<0.05 in both, or BH-adjusted p<0.1 in both. '
             'Every result above 75% directional concordance is retained, including small overlaps. '
             'Experimental fold magnitudes are not gated because the supplied effect log base is unconfirmed.', '',
             'Agreement = concordant/(concordant + opposite), among lipids significant in both and passing the '
             'predicted fold cutoff. Exactly 75% does not qualify. Zero experimental effects are listed as ambiguous. '
             'Signed +2 represents twofold higher and -2 twofold lower; the magnitude threshold is symmetric.', '',
             'These are summaries of unchanged full-202-lipid candidate/prior ElasticNetCV results. '
             'BH remains over all 202 outputs within each original comparison. No refitting, imputation, raw testing, '
             'or BH recalculation was performed. Metacell p-values are nominal relative to donor replication; '
             'single-donor arms do not establish replicated disease validation. '
             'Selection over thresholds/states is descriptive, not independent validation. '
             'Natri pooled and fibrotic strata share donors.', '',
             '## Candidate: largest concordant overlap for each requested setting', '', table(strongest[fields].sort_values(['gate','predicted_fold_strictly_greater_than'])), '',
             '## Candidate: largest overlaps among states with multiple donors in both arms', '',
             'This display facet requires at least two source donors per arm. Every single-donor comparison remains '
             'in the complete results above and below; no model feature or sample is removed. '
             'Metacell testing still uses metacell observations, so source donor coverage does not make these donor-level tests.', '',
             table(strongest_multiple_donors[fields].sort_values(['gate','predicted_fold_strictly_greater_than'])), '',
             '## Number of candidate same-state comparisons exceeding 75%', '', table(counts), '',
             '## Candidate: best cell state for every cohort/input/threshold setting', '',
             table(best[fields].sort_values(['gate','predicted_fold_strictly_greater_than','contrast','unit'])), '',
             '## Candidate: macrophage comparisons', '', table(macro[fields].sort_values(['gate','predicted_fold_strictly_greater_than','contrast','unit'])), '',
             '## Every qualifying comparison for both models', '',
             table(qualify[['model','comparison_kind'] + fields].sort_values(['model','gate','predicted_fold_strictly_greater_than','contrast','unit','population'])), '',
             '## Every specific lipid in every qualifying comparison', '']
    for gkey, g in pd.concat(details,ignore_index=True).groupby(KEYS + ['gate','predicted_fold_strictly_greater_than']):
        lines += ['', '### ' + ' / '.join(str(x) for x in gkey), '',
                  table(g[['measured_feature','gene','signed_fold','pval','fdr','supplied_raw_p','supplied_adjusted_p','direction_relation']])]
    lines += ['', '## Downloads', '',
              '- [Every >75% comparison and lipid identity](all_above75_with_lipid_identities.csv)',
              '- [All nine settings for every tested comparison, including nonqualifying results](all_requested_threshold_results.csv)',
              '- [Candidate versus prior at identical thresholds](candidate_vs_prior_identical_thresholds.csv)',
              '- [Replication census, including untested comparisons](complete_replication_census.csv)']
    (OUT / 'REQUESTED_THRESHOLDS_REPORT.md').write_text('\n'.join(lines)+'\n')
    audit = {'source_sha256': hashes, 'script_sha256': sha(Path(__file__)), 'fold_operator': 'strictly greater than',
             'fold_cutoffs': FOLDS, 'gates': [g[0] for g in GATES], 'full202_comparison_groups':verified,
             'settings_per_group':9, 'scan_rows':len(all_rows), 'above75_rows_by_model':qualify.groupby('model').size().to_dict(),
             'independent_scalar_count_verifications':manual_checks, 'source_hashes_unchanged':all(sha(p)==hashes[str(p)] for p in paths),
             'new_training_inference_raw_tests_or_BH_recalculation':False}
    (OUT / 'requested_threshold_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))
    print('\nSTRONGEST:\n'+strongest[fields].sort_values(['gate','predicted_fold_strictly_greater_than']).to_string(index=False))
    print('\nCOUNTS:\n'+counts.to_string(index=False))


if __name__ == '__main__':
    main()
