"""Audit the corrected lung-tissue table and compare independently fitted RNA predictions."""
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from ipf_utils import canonical, fold_table, concordance


def mode(feature):
    return 'P' if str(feature).endswith(('_POS', '_P')) else 'N' if str(feature).endswith(('_NEG', '_N')) else None


def load_lung(source, out, rng, bootstraps):
    path = Path(source) / '10_results_with_statistics.csv'
    raw = pd.read_csv(path, index_col=0)
    columns = [c for c in raw if re.fullmatch(r'L\d+_\d+', c)]
    if len(columns) != 40 or raw.index.duplicated().any():
        raise ValueError('Unexpected or ambiguous lung source matrix')
    # Column assignment is checked against both provided group summaries.
    case, control = columns[:30], columns[30:]
    if max(abs(raw[case].mean(axis=1) - raw.IPF_mean).max(), abs(raw[control].mean(axis=1) - raw.Ctrl_mean).max()) > .11:
        raise ValueError('IPF/control sample assignment does not reproduce rounded source summaries')
    values = raw[columns].T.astype(float)
    md = pd.DataFrame({'donor': [s.split('_')[0] for s in columns],
                       'group': ['IPF'] * 30 + ['NDC'] * 10}, index=columns)
    if (md.groupby('donor').group.nunique() > 1).any():
        raise ValueError('Conflicting tissue diagnoses within donor')
    donor_log = values.groupby(md.donor).mean()
    donor_linear = np.exp2(values).groupby(md.donor).mean()
    groups = md.groupby('donor').group.first().loc[donor_log.index]
    stats = fold_table(donor_linear, donor_log, groups, rng, bootstraps)
    stats['original_annotation'] = raw.Lipid
    stats['supplied_log_effect'] = raw['log(IPF/Ctrl)']
    stats['supplied_raw_p'] = raw.IPF_vs_Ctrl_Ttest_p
    stats['supplied_adjusted_p'] = raw.IPF_vs_Ctrl_Ttest_padj
    stats['supplied_BH_recomputed_from_rounded_p'] = multipletests(raw.IPF_vs_Ctrl_Ttest_p, method='fdr_bh')[1]
    stats['source_row'] = np.arange(1, len(raw) + 1)
    stats['canonical'] = raw.Lipid.map(canonical)
    stats['ion_mode'] = [mode(f) for f in raw.index]
    stats['feature_key'] = [f'{c}|{m}' if c and m else None for c, m in zip(stats.canonical, stats.ion_mode)]
    values.to_csv(Path(out) / 'measured_lung_native_samples.csv')
    donor_log.to_csv(Path(out) / 'measured_lung_donor_native.csv')
    donor_linear.to_csv(Path(out) / 'measured_lung_donor_relative_linear.csv')
    md.to_csv(Path(out) / 'measured_lung_sample_metadata.csv')
    stats.to_csv(Path(out) / 'measured_lung_IPF_vs_control.csv')
    audit = {'source_file': str(path), 'tissue': 'human lung, user-confirmed', 'features': len(raw),
             'IPF_tissue_samples': 30, 'control_tissue_samples': 10,
             'IPF_donors': int((groups == 'IPF').sum()), 'control_donors': int((groups == 'NDC').sum()),
             'native_range': [float(values.min().min()), float(values.max().max())],
             'maximum_group_mean_rounding_error': float(max(abs(raw[case].mean(axis=1) - raw.IPF_mean).max(), abs(raw[control].mean(axis=1) - raw.Ctrl_mean).max())),
             'maximum_supplied_adjustment_vs_BH_of_rounded_p': float(abs(stats.supplied_BH_recomputed_from_rounded_p - stats.supplied_adjusted_p).max()),
             'scale': 'Transformed relative signal, not linear abundance; log2 interpretation provisional pending protocol confirmation. Supplied Z-score columns are not used.',
             'replication': 'L-prefix interpreted as donor; three spatial tissue measurements per IPF donor are averaged. Independent donor-level Welch/BH tests use 10 versus 10.',
             'paper': 'Dr. Clair email, IPF and BPD lipid differentials: 10_ is IPF lipidomics, statistics exported from a romics object; matching proteomics and metabolomics exist; comparison paper planned around fall. No published citation or log base supplied.'}
    audit['rounded_zero_raw_p_features'] = int(raw.IPF_vs_Ctrl_Ttest_p.eq(0).sum())
    audit['rounded_zero_adjusted_p_features'] = int(raw.IPF_vs_Ctrl_Ttest_padj.eq(0).sum())
    mode_rows = []
    for ion, group in stats.groupby('ion_mode'):
        mode_rows.append({'ion_mode': ion, 'features': len(group),
                          'supplied_increases': int(group.supplied_log_effect.gt(0).sum()),
                          'supplied_decreases': int(group.supplied_log_effect.lt(0).sum()),
                          'median_supplied_effect': float(group.supplied_log_effect.median())})
    pd.DataFrame(mode_rows).to_csv(Path(out) / 'source_ion_mode_summary.csv', index=False)
    (Path(out) / 'lung_source_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    return path, stats, stats[['original_annotation','canonical','ion_mode','feature_key','source_row']], audit


def compare(stats, annotation, predicted):
    rows, audit = [], []
    for feature, pred in predicted.iterrows():
        key = f'{canonical(feature)}|{mode(feature)}' if canonical(feature) and mode(feature) else None
        ids = annotation.index[annotation.feature_key.eq(key)] if key else []
        audit.append({'model_feature': feature, 'feature_key': key, 'source_matches': len(ids),
                      'status': 'matched' if len(ids) == 1 else 'absent_or_ambiguous_exact_mode_match'})
        if len(ids) != 1:
            continue
        m = stats.loc[ids[0]]
        rows.append({'model_feature': feature, 'measured_feature': ids[0], 'paper_row': ids[0],
                     'canonical': m.canonical, 'ion_mode': m.ion_mode,
                     'measured_log2FC': m.log2FC, 'predicted_log2FC': pred.log2FC,
                     'measured_fold': m.fold_change, 'predicted_fold': pred.fold_change,
                     'measured_geometric_log2FC': m.geometric_log2FC,
                     'predicted_geometric_log2FC': pred.geometric_log2FC,
                     'supplied_log_effect': m.supplied_log_effect, 'supplied_raw_p': m.supplied_raw_p,
                     'supplied_adjusted_p': m.supplied_adjusted_p,
                     'measured_raw_p': m.pvalue, 'measured_BH_full_544': m.FDR,
                     'predicted_raw_p': pred.pvalue, 'predicted_BH_full_panel': pred.FDR})
    frame = pd.DataFrame(rows)
    summary = {'all_arithmetic_mean_folds': concordance(frame), 'matched_features': len(frame)}
    for label, selected in [('top25_matched_by_supplied_p', frame.sort_values(['supplied_raw_p','measured_feature']).head(25)),
                            ('global_top25_by_supplied_p', frame[frame.measured_feature.isin(stats.sort_values(['supplied_raw_p','source_row']).head(25).index)]),
                            ('top25_matched_by_donor_p', frame.sort_values(['measured_raw_p','measured_feature']).head(25))]:
        primary = concordance(selected)
        exact = selected.copy()
        exact['measured_log2FC'] = exact.supplied_log_effect
        exact['predicted_log2FC'] = exact.predicted_geometric_log2FC
        summary[label] = {'arithmetic_mean_folds': primary, 'supplied_effect_vs_predicted_geometric': concordance(exact),
                         'features': selected.measured_feature.tolist()}
    up = frame.measured_fold.gt(1) & frame.predicted_fold.gt(1)
    down = frame.measured_fold.lt(1) & frame.predicted_fold.lt(1)
    both = frame.measured_BH_full_544.le(.1) & frame.predicted_BH_full_panel.le(.1)
    summary['concordance_counts'] = {'up': int(up.sum()), 'down': int(down.sum()),
                                   'both_BH_up': int((up & both).sum()), 'both_BH_down': int((down & both).sum())}
    summary['by_ion_mode'] = {ion: {'matched_features': len(g),
                                  'supplied_geometric_direction_concordant': int((g.supplied_log_effect * g.predicted_geometric_log2FC > 0).sum())}
                             for ion, g in frame.groupby('ion_mode')}
    tied = frame[frame.supplied_raw_p.eq(0)]
    summary['all_matched_zero_reported_p_ties'] = {'features': len(tied),
        'supplied_geometric_direction_concordant': int((tied.supplied_log_effect * tied.predicted_geometric_log2FC > 0).sum())}
    return frame, pd.DataFrame(audit), summary
