"""Read-only source audit and original mean/statistics tables; no model fitting."""

from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control


ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'artifacts' / 'Geremy_original_lipid_statistics_20261008'
SOURCES = {
    'IPF': Path('/Users/saljh8/Downloads/Geremy/10_results_with_statistics.xlsx'),
    'BPD': Path('/Users/saljh8/Downloads/Geremy/21_results_lipids_with_stats.xlsx'),
}
EXPECTED = {'IPF': (544, 40), 'BPD': (530, 23)}
CSV = Path('/Users/saljh8/Downloads/Lipidomics/IPF/10_results_with_statistics.csv')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def number(value):
    if pd.isna(value):
        return 'NA'
    if isinstance(value, (float, int, np.number)):
        return format(float(value), '.12g')
    return str(value).replace('|', r'\|').replace('\n', ' ')


def markdown(frame):
    return '\n'.join([
        '| ' + ' | '.join(map(str, frame.columns)) + ' |',
        '| ' + ' | '.join(['---'] * len(frame.columns)) + ' |',
        *['| ' + ' | '.join(number(x) for x in row) + ' |'
          for row in frame.itertuples(index=False, name=None)],
    ])


def prepare_bulk_source_tables(ipf):
    """Attach original workbook summaries to verified, unchanged model statistics.

    This audit needs no lipid log-base assumption and performs no threshold scan.
    """
    source = ROOT / 'artifacts' / 'bulk_IPF_three_version_comparison_20261008'
    audit = json.loads((source / 'THREE_VERSION_AUDIT.json').read_text())
    assert audit['all_checks_passed']
    contract = json.loads((ROOT / 'integrity/baseline_contract.json').read_text())
    frame = pd.read_csv(source / 'three_versions_all_202_lipids.csv',
                        float_precision='round_trip')
    lookup = ipf.set_index(ipf.columns[0])
    assert lookup.index.is_unique and len(lookup) == 544
    source_columns = {'Ctrl_mean': 'measured_control_mean',
                      'Ctrl_sd': 'measured_control_sd',
                      'IPF_mean': 'measured_IPF_mean',
                      'IPF_sd': 'measured_IPF_sd',
                      'log(IPF/Ctrl)': 'supplied_log_effect',
                      'IPF_vs_Ctrl_Ttest_p': 'supplied_raw_p',
                      'IPF_vs_Ctrl_Ttest_padj': 'supplied_adjusted_p'}
    frame.insert(0, 'source_row', np.nan)
    for column in ['supplied_log_effect', 'supplied_raw_p', 'supplied_adjusted_p']:
        frame[f'previous_CSV_{column}'] = frame[column]
    for dest in source_columns.values():
        frame[dest] = np.nan
    for version, group in frame.groupby('version', sort=False):
        assert group.model_lipid.tolist() == contract['lipids']
        assert len(group) == 202 and not group.model_lipid.duplicated().any()
        assert group['model'].nunique() == group.RNA_representation.nunique() == 1
        model, representation = group.iloc[0][['model', 'RNA_representation']]
        stats = pd.read_csv(source / 'full_panel_rerun' /
                            f'{model}_{representation}_donor_balanced_all_202_differentials.csv',
                            float_precision='round_trip').set_index('feature')
        assert stats.index.tolist() == contract['lipids']
        for dest, original in [('predicted_IPF_mean_linear', 'IPF_mean'),
                               ('predicted_control_mean_linear', 'control_mean')]:
            frame.loc[group.index, dest] = stats.loc[group.model_lipid, original].to_numpy()
        for existing, original in [('predicted_raw_p', 'pvalue'), ('predicted_BH', 'FDR'),
                                   ('predicted_signed_fold', 'signed_fold')]:
            np.testing.assert_allclose(group[existing], stats.loc[group.model_lipid, original],
                                       rtol=1e-12, atol=1e-14)
        matched = group.status.eq('matched')
        ids = group.loc[matched, 'measured_feature']
        if not ids.isin(lookup.index).all():
            raise ValueError('Required measured feature unresolved: stop and ask for mapping')
        ix = group.loc[matched].index
        frame.loc[ix, 'source_row'] = lookup.index.get_indexer(ids) + 2
        for original, dest in source_columns.items():
            frame.loc[ix, dest] = lookup.loc[ids, original].to_numpy()
        assert matched.sum() == 134
    assert len(frame) == 606
    assert frame.loc[frame.status.eq('matched'), list(source_columns.values())].notna().all().all()
    display = ['model_lipid', 'measured_feature', 'status', 'measured_control_mean',
               'measured_control_sd', 'measured_IPF_mean', 'measured_IPF_sd',
               'supplied_log_effect', 'supplied_raw_p', 'supplied_adjusted_p',
               'predicted_control_mean_linear', 'predicted_IPF_mean_linear',
               'predicted_signed_fold', 'predicted_raw_p', 'predicted_BH']
    lines = ['# Bulk IPF: original means and existing prediction statistics', '',
             'The full-precision XLSX statistics are shown next to the existing verified '
             'bulk prediction results. All 202 model outputs are retained for each of '
             'the three versions. Measured means and SDs remain on the original centered '
             'scale; imputed means are arithmetic means on each model’s linear output '
             'scale. Their numerical values are not directly interchangeable.', '',
             'RNA biological replicates: 20 IPF and 14 control donors in independent '
             'groups. These donors are not paired to the lipidomics samples. Predicted '
             'raw tests and full-202 BH values are unchanged. Original experimental '
             'raw p and full-544 BH values are unchanged; XLSX retains more precision '
             'than the previous CSV. No measured fold conversion or cutoff is applied.', '',
             'Each model has the same 134 established measured-feature matches. The '
             'remaining 68 outputs retain their established unresolved matching status '
             'and their prediction statistics. No new mapping or feature exclusion is '
             'introduced. Prior reports are preserved.', '']
    for version, group in frame.groupby('version', sort=False):
        lines += [f'## {version}', '', markdown(group[display]), '']
    return frame, '\n'.join(lines)


def main():
    # All checks run before creating any report. Original workbooks stay untouched.
    evidence = {'scope': 'Original workbook audit and tables only; no model comparison',
                'verified_log_base': None, 'hypothesis': 'Sample-median-centered log2',
                'sources': {}}
    data = {}
    summary = []
    for disease, path in SOURCES.items():
        original_hash = sha(path)
        book = pd.ExcelFile(path)
        assert len(book.sheet_names) == 1
        frame = pd.read_excel(book, sheet_name=book.sheet_names[0])
        # BPD's final Excel row is a COUNTIF summary, not a lipid feature.
        mask = frame.iloc[:, 0].notna()
        features = frame.loc[mask].copy()
        summary_rows = frame.loc[~mask].copy()
        assert not features.iloc[:, 0].duplicated().any()
        samples = (list(frame.columns[7:47]) if disease == 'IPF'
                   else [c for c in frame if c.startswith('LMap_')])
        assert (len(features), len(samples)) == EXPECTED[disease]
        assert len(summary_rows) == (0 if disease == 'IPF' else 1)
        effect = f'log({disease}/Ctrl)'
        p = f'{disease}_vs_Ctrl_Ttest_p'
        q = f'{disease}_vs_Ctrl_Ttest_padj'
        delta = features[f'{disease}_mean'] - features['Ctrl_mean']
        residual = float(np.max(np.abs(delta - features[effect])))
        assert residual < 2e-13
        assert features[[p, q, effect, 'Ctrl_mean', f'{disease}_mean']].notna().all().all()
        bh_error = float(np.max(np.abs(
            false_discovery_control(features[p].to_numpy(), method='bh') - features[q]
        )))
        assert bh_error < 1e-12
        values = features[samples].to_numpy(dtype=float)
        means = [c for c in frame if c.endswith('_mean')]
        stats = [c for c in frame if c.endswith(('_sd', '_p', '_padj'))
                 or c.startswith('log(') or c.endswith('_directionality')]
        # All original contrasts and group summaries, with original column names.
        cols = [frame.columns[0]] + means + stats
        complete_stats = features[cols].copy()
        selected = features[[frame.columns[0], 'Ctrl_mean', 'Ctrl_sd',
                             f'{disease}_mean', f'{disease}_sd', effect, p, q]].copy()
        source_rows = np.flatnonzero(mask) + 2
        selected.insert(0, 'Excel_row', source_rows)
        complete_stats.insert(0, 'Excel_row', source_rows)
        sample_medians = np.nanmedian(values, axis=0)
        entry = {
            'path': str(path), 'sha256': original_hash, 'sheet': book.sheet_names[0],
            'feature_rows': len(features), 'source_rows': len(frame),
            'sample_columns': samples, 'sample_count': len(samples),
            'missing_sample_values_retained': int(np.isnan(values).sum()),
            'sample_minimum': float(np.nanmin(values)),
            'sample_maximum': float(np.nanmax(values)),
            'max_abs_sample_median': float(np.max(np.abs(sample_medians))),
            'max_abs_effect_minus_mean_difference': residual,
            'max_abs_original_padj_minus_independent_BH': bh_error,
            'non_feature_summary_rows_retained': (np.flatnonzero(~mask) + 2).tolist(),
            'original_summary_rows': json.loads(summary_rows.to_json(orient='records')),
        }
        for label, column, cutoff in [('Raw p < 0.05', p, .05),
                                       ('BH < 0.05', q, .05), ('BH < 0.1', q, .1)]:
            passing = features[column] < cutoff
            counts = {'Disease': disease, 'Criterion': label,
                      'Increased': int((passing & (features[effect] > 0)).sum()),
                      'Decreased': int((passing & (features[effect] < 0)).sum()),
                      'Total': int(passing.sum())}
            summary.append(counts)
        if disease == 'IPF':
            csv = pd.read_csv(CSV)
            assert list(csv) == list(frame)
            assert csv.iloc[:, 0].tolist() == features.iloc[:, 0].tolist()
            assert np.allclose(np.round(features[effect], 1), csv[effect], atol=1e-14)
            entry['prior_csv_sha256'] = sha(CSV)
            entry['prior_csv_identical_features_and_order'] = True
            entry['csv_rounded_effect_zeros'] = int((csv[effect] == 0).sum())
            entry['xlsx_exact_effect_zeros'] = int((features[effect] == 0).sum())
            entry['csv_rounded_raw_p_zeros'] = int((csv[p] == 0).sum())
            entry['csv_rounded_BH_zeros'] = int((csv[q] == 0).sum())
            entry['xlsx_raw_p_zeros'] = int((features[p] == 0).sum())
            entry['xlsx_BH_zeros'] = int((features[q] == 0).sum())
            entry['threshold_call_changes_vs_csv'] = {
                f'{c}<{t}': int(((features[c] < t) != (csv[c] < t)).sum())
                for c, t in [(p, .05), (q, .05), (q, .1)]
            }
            assert not any(entry['threshold_call_changes_vs_csv'].values())
        evidence['sources'][disease] = entry
        data[disease] = (selected, complete_stats, frame)

    OUT.mkdir(parents=True, exist_ok=True)
    bulk_table, bulk_report = prepare_bulk_source_tables(data['IPF'][2])
    bulk_table.to_csv(OUT / 'three_versions_all_202_with_original_means.csv', index=False)
    (OUT / 'BULK_IPF_ORIGINAL_MEANS_AND_PREDICTIONS.md').write_text(bulk_report)
    report = [
        '# Original measured lipid means and statistics', '',
        'These tables reproduce the original workbook values. No means, SDs, tests, '
        'or adjusted p-values were refitted. All 544 IPF and 530 BPD lipid features '
        'are retained, including distinct ion-mode and isomer annotations.', '',
        'Means and SDs remain in the supplied normalized measurement space. '
        'Every sample median is zero. The reported log effect equals disease mean '
        'minus control mean to numerical precision. Sample-median-centered log2 '
        'is the current best-supported hypothesis, but the workbooks do not establish '
        'the logarithm base. Dividing these centered means is inappropriate. '
        'Under a log2 assumption, normalized geometric-mean fold = 2^(mean difference). '
        'No such assumption or fold cutoff is applied in these original-statistics tables.', '',
        'The original Ttest_padj values independently reproduce Benjamini–Hochberg '
        'adjustment across all 544 IPF features and all 530 BPD features, respectively. '
        'Counts below use the original pooled disease-versus-control comparison '
        'without a fold cutoff and without a model-match restriction.', '',
        markdown(pd.DataFrame(summary)), '',
        'The IPF CSV contains the same ordered feature panel but rounds values. '
        'It rounds 110 raw p-values and 75 BH values to zero. The XLSX contains '
        'their nonzero values. Raw p < 0.05, BH < 0.05 and BH < 0.1 calls remain '
        'identical. Nineteen small nonzero effects are rounded to zero in the CSV; '
        'none of those has BH < 0.05. Ranking tied CSV p-values can differ from '
        'ranking full-precision workbook p-values.', '',
        'The BPD workbook has 530 lipid rows and a separate final Excel row 532 '
        'containing COUNTIF summary formulas. That row is preserved in the complete '
        'worksheet CSV and audit, and is identified as a summary rather than a lipid. '
        'Four missing BPD sample measurements are retained. BPD subgroup comparisons '
        'are supplied in the full statistics CSV; pooled BPD-versus-control counts '
        'do not describe those subgroup contrasts.', '',
    ]
    for disease, (selected, complete, frame) in data.items():
        selected.to_csv(OUT / f'{disease}_original_means_and_statistics.csv', index=False)
        complete.to_csv(OUT / f'{disease}_all_original_group_means_and_statistics.csv', index=False)
        frame.to_csv(OUT / f'{disease}_complete_original_worksheet.csv', index=False)
        entry = evidence['sources'][disease]
        report.extend([
            f'## {disease}: all original lipid results', '',
            f'Source: `{entry["path"]}`, sheet `{entry["sheet"]}`. '
            f'{entry["feature_rows"]} lipid features and '
            f'{entry["sample_count"]} measured sample/profile columns. '
            'Sample columns are not assumed to be independent donors.', '',
            'Control and disease means and SDs are reproduced on the original '
            'scale. The effect is the original reported log difference. '
            'Markdown displays 12 significant digits; CSV files retain Python’s '
            'round-trip numerical precision. Excel_row links each lipid to its '
            'source row. Tables follow original source order.', '',
            markdown(selected), '',
        ])
        # Verify every exported value, including source identifiers and source order.
        reread = pd.read_csv(OUT / f'{disease}_original_means_and_statistics.csv',
                             float_precision='round_trip')
        pd.testing.assert_frame_equal(reread, selected.reset_index(drop=True),
                                      check_dtype=False, check_exact=True)
        assert sha(SOURCES[disease]) == entry['sha256']
    (OUT / 'ORIGINAL_MEANS_AND_STATISTICS.md').write_text('\n'.join(report))
    (OUT / 'workbook_audit.json').write_text(json.dumps(evidence, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'report': str(OUT / 'ORIGINAL_MEANS_AND_STATISTICS.md'),
                      'counts': summary, 'source_hashes_unchanged': True}))


if __name__ == '__main__':
    main()
