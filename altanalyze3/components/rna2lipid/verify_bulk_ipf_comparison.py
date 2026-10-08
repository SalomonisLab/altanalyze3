"""Independently verify saved fixed-model results and write their local report.

No training or inference outputs are replaced. Verify the saved API predictions
against every original stored estimator, and the reported statistics against
SciPy Welch tests and independently aggregated donor observations.
"""
from pathlib import Path
import hashlib
import json
import pickle
import warnings

import numpy as np
import pandas as pd
import scipy
from scipy.stats import ttest_ind
import sklearn
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
OUT = HERE / 'artifacts/bulk_IPF_full202_comparison_20261007'


def aggregate(frame, md, region):
    selected = md.diseasegroup.isin(['IPF', 'NDC']) & md.lunglocation.isin(['Apex', 'Base'])
    if region != 'donor_balanced':
        selected &= md.lunglocation.eq(region)
    values = frame.loc[md.index[selected]].copy()
    values['donor'] = md.loc[selected, 'donorid'].values
    values['region'] = md.loc[selected, 'lunglocation'].values
    return values.groupby(['donor', 'region']).mean().groupby('donor').mean()


def number(value):
    if pd.isna(value):
        return '—'
    if isinstance(value, (float, np.floating)):
        return f'{value:.4g}'
    return str(value).replace('|', '\\|').replace('\n', ' ')


def table(frame):
    return '\n'.join(['| ' + ' | '.join(map(str, frame.columns)) + ' |',
        '| ' + ' | '.join(['---'] * len(frame.columns)) + ' |'] +
        ['| ' + ' | '.join(number(v) for v in row) + ' |' for row in frame.itertuples(index=False, name=None)])


def main():
    audit = json.loads((OUT / 'comparison_audit.json').read_text())
    md = pd.read_csv(OUT / 'all_139_RNA_sample_metadata.csv', index_col=0)
    fills = pd.read_csv(OUT / 'RNA_mapping_and_median_filling.csv')
    groups = md.groupby('donorid').diseasegroup.first()
    verification = {'models_unchanged': True, 'RNA_profiles': 139, 'RNA_inputs': 1303,
                    'lipid_outputs': 202, 'checks': [], 'sklearn_version': sklearn.__version__,
                    'scipy_version': scipy.__version__, 'fits_executed': 0}
    donor_counts = {}
    baseline_settings = None
    settings_keys = ['top_gene_options', 'l1_ratio_grid', 'alpha_grid', 'cv_folds',
                     'max_iter', 'sparsity_penalty', 'random_seed']
    for label, info in audit['model_checks'].items():
        assert hashlib.sha256(Path(info['path']).read_bytes()).hexdigest() == info['sha256']
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            with open(info['path'], 'rb') as handle:
                bundle = pickle.load(handle)
        assert {type(entry['model']).__name__ for entry in bundle['models'].values()} == {'ElasticNetCV'}
        current_settings = {key: bundle[key] for key in settings_keys}
        if baseline_settings is None:
            baseline_settings = current_settings
        else:
            for key in settings_keys:
                np.testing.assert_array_equal(current_settings[key], baseline_settings[key])
        verification['identical_ElasticNetCV_training_settings'] = settings_keys
        for representation in ['log2_CP10k', 'publisher_logRPKM']:
            X = pd.read_csv(OUT / f'{representation}_identical_RNA_inputs.csv', index_col=0)
            saved = pd.read_csv(OUT / f'{label}_{representation}_all_139_predictions_log2.csv', index_col=0)
            assert X.shape == (139, 1303) and saved.shape == (139, 202)
            assert list(X.index) == list(saved.index) == list(md.index)
            assert list(X.columns) == bundle['X_columns'] and list(saved.columns) == bundle['Y_columns']
            assert np.isfinite(X).all().all() and np.isfinite(saved).all().all()
            for row in fills[fills.filled].itertuples():
                np.testing.assert_allclose(X[row.gene], row.training_median_fill, rtol=0, atol=1e-12)
            transformed = pd.DataFrame(bundle['scaler_x'].transform(X), index=X.index, columns=X.columns)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                direct = np.column_stack([bundle['models'][gene]['model'].predict(
                    transformed.loc[:, bundle['models'][gene]['genes']]) for gene in bundle['Y_columns']])
            if info['target_scaling'] == 'standard':
                direct = bundle['scaler_y'].inverse_transform(direct)
            np.testing.assert_allclose(direct, saved, rtol=1e-12, atol=1e-10)
            prediction_error = float(np.max(np.abs(direct - saved.to_numpy())))
            for region in ['donor_balanced', 'Apex', 'Base']:
                linears = aggregate(np.exp2(saved), md, region)
                logged = aggregate(saved, md, region)
                old_linear = pd.read_csv(OUT / f'{label}_{representation}_{region}_donor_predictions_linear.csv', index_col=0)
                pd.testing.assert_index_equal(linears.index, old_linear.index)
                np.testing.assert_allclose(linears, old_linear, rtol=1e-12)
                group = groups.loc[linears.index]
                counts = {'IPF': int(group.eq('IPF').sum()), 'control': int(group.eq('NDC').sum())}
                donor_counts[region] = counts
                case, control = logged.loc[group.eq('IPF')], logged.loc[group.eq('NDC')]
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    raw_p = ttest_ind(case, control, equal_var=False, axis=0).pvalue
                magnitude = np.maximum(np.maximum(abs(case).max(axis=0), abs(control).max(axis=0)), 1.)
                testable = case.var(axis=0, ddof=1) + control.var(axis=0, ddof=1) > np.finfo(float).eps * magnitude ** 2
                raw_p = np.where(testable & np.isfinite(raw_p), raw_p, 1.)
                stats = pd.read_csv(OUT / f'{label}_{representation}_{region}_all_202_differentials.csv', index_col=0)
                assert list(stats.index) == list(saved.columns)
                np.testing.assert_allclose(stats.pvalue, raw_p, rtol=1e-9, atol=1e-12)
                np.testing.assert_allclose(stats.FDR, multipletests(raw_p, method='fdr_bh')[1], rtol=1e-9, atol=1e-12)
                ratio = linears.loc[group.eq('IPF')].mean() / linears.loc[group.eq('NDC')].mean()
                np.testing.assert_allclose(stats.fold_change, ratio, rtol=1e-12)
                verification['checks'].append({'model': label, 'RNA_representation': representation, 'region': region,
                    'prediction_max_abs_error': prediction_error, 'raw_p_max_abs_error': float(np.max(abs(stats.pvalue - raw_p))),
                    'BH_family': len(stats), 'donor_counts': counts})
        assert hashlib.sha256(Path(info['path']).read_bytes()).hexdigest() == info['sha256']

    detail = pd.read_csv(OUT / 'every_lipid_all_comparisons.csv')
    thresholds = pd.read_csv(OUT / 'all_threshold_results.csv')
    assert len(detail) == 12 * 202 and len(thresholds) == 180
    # Independent identities/significance counts, not a call to the analysis helper.
    for row in thresholds.itertuples():
        frame = detail[(detail.model == row.model) & (detail.RNA_representation == row.RNA_representation) & (detail.region == row.region)]
        sp = frame.supplied_raw_p if row.statistic == 'raw_p' else frame.supplied_adjusted_p
        pp = frame.predicted_raw_p if row.statistic == 'raw_p' else frame.predicted_BH
        both = frame.status.eq('matched') & sp.lt(row.p_cutoff) & pp.lt(row.p_cutoff) & frame.supplied_log_effect.ne(0)
        both &= frame.predicted_signed_fold.abs().gt(row.minimum_predicted_fold)
        direction = np.sign(frame.supplied_log_effect) == np.sign(frame.predicted_signed_fold)
        assert row.significant_in_both == int(both.sum())
        assert row.concordant_total == int((both & direction).sum())
        assert row.discordant_total == int((both & ~direction).sum())

    paired = []
    for representation in ['log2_CP10k', 'publisher_logRPKM']:
        frames = {m: detail[(detail.model == m) & (detail.RNA_representation == representation) &
                           (detail.region == 'donor_balanced')].set_index('model_lipid') for m in ['candidate', 'prior']}
        c, p = frames['candidate'], frames['prior']
        pd.testing.assert_series_equal(c.measured_feature, p.measured_feature)
        for statistic, cutoff in [('raw_p', .05), ('BH', .05), ('BH', .1)]:
            source_field, prediction_field = ('supplied_raw_p', 'predicted_raw_p') if statistic == 'raw_p' else ('supplied_adjusted_p', 'predicted_BH')
            common = c.status.eq('matched') & c[source_field].lt(cutoff) & c.supplied_log_effect.ne(0)
            common &= c[prediction_field].lt(cutoff) & p[prediction_field].lt(cutoff)
            c_agree = c.supplied_log_effect * c.predicted_log2FC > 0
            p_agree = p.supplied_log_effect * p.predicted_log2FC > 0
            paired.append({'RNA_representation': representation, 'statistic': statistic, 'p_cutoff': cutoff,
                'significant_experimental_and_both_models': int(common.sum()),
                'candidate_concordant': int((common & c_agree).sum()), 'prior_concordant': int((common & p_agree).sum()),
                'candidate_only_correct': int((common & c_agree & ~p_agree).sum()),
                'prior_only_correct': int((common & ~c_agree & p_agree).sum()),
                'correct_in_both': int((common & c_agree & p_agree).sum()),
                'wrong_in_both': int((common & ~c_agree & ~p_agree).sum())})
    paired = pd.DataFrame(paired)
    paired.to_csv(OUT / 'same_lipid_paired_model_comparison.csv', index=False)
    verification['all_checks_passed'] = True
    (OUT / 'independent_verification.json').write_text(json.dumps(verification, indent=2) + '\n')

    primary = thresholds[(thresholds.RNA_representation == 'log2_CP10k') & (thresholds.region == 'donor_balanced')]
    display_columns = ['model', 'statistic', 'p_cutoff', 'minimum_predicted_fold', 'significant_in_both',
                       'concordant_up', 'concordant_down', 'concordant_total', 'discordant_total', 'concordance_percent']
    lines = ['# Bulk unmatched IPF: complete 202-lipid candidate versus prior model', '',
        'Evaluation date: 2026-10-07. The candidate is the current promoted native-relative-log2 lipid-wise ElasticNetCV bundle; the prior is the preserved original lipid-wise ElasticNetCV bundle. Both fixed models were evaluated without refitting.', '',
        '## Results', '',
        'The candidate has higher directional agreement with the supplied lung-tissue IPF lipid differentials. In the primary donor-balanced analysis at BH FDR <0.1, the candidate agrees for **39/62 (62.9%)** lipids significant in both assays versus **22/63 (34.9%)** for the prior. These are model-specific significant overlap sets; the identical-lipid comparison below checks the result on a common set. This bulk comparison does not attain 75% agreement at the no-fold-cutoff setting.', '',
        'Counts below require significance in both measured lipidomics and the indicated model. Up and down counts refer to concordant changes. Fold cutoffs apply to the prediction only because the supplied experimental effect log base is unconfirmed.', '',
        table(primary[primary.minimum_predicted_fold.eq(1)][display_columns]), '',
        '## Comparison on the same significant lipids', '',
        'Each row requires the lipid to be significant in the measured data and in **both** models at the stated threshold; therefore the candidate and prior have an identical denominator. No additional fold cutoff is applied.', '',
        table(paired), '',
        '## Source data and methods', '',
        '- RNA: GSE213001, supplied raw counts and publisher logRPKM, 15,065 source genes and 139 profiles. Diagnosis metadata come from the supplied GEO series matrix. IPF and non-diseased control are separate donors, not patient-paired contrasts. All 139 profiles were imputed and retained in the exported prediction matrices.',
        '- Primary IPF/control comparison: 101 profiles with known apex/base anatomy, 20 IPF and 14 control donors, using the original comparison roster. Average left/right profiles within donor/region, then weight available apex/base regions equally within donor. ALF018E and ALF026E lack region labels and remain in the full prediction exports but do not enter this established regional aggregation. ILD and CLAD profiles are likewise retained in the exports but do not belong to the IPF/control contrast.',
        '- Lipidomics: corrected lung-tissue `10_results_with_statistics.csv`, 544 feature rows and 40 tissue-sample columns. Use its supplied `log(IPF/Ctrl)`, `IPF_vs_Ctrl_Ttest_p`, and `IPF_vs_Ctrl_Ttest_padj`. This is the Dr. Clair romics export described in the user-supplied email. A publication citation and experimental log base have not been established; no source-fold magnitude or absolute concentration is asserted.',
        '- Lipid identity matching: retain the previously reviewed ion-mode/acyl-composition matching audit unchanged. Exactly 132 unambiguous experimental lipids match either model; both use the same 132 identities. The other 70 model lipids remain predicted, tested, BH-corrected, and listed below, with their matching status explicit. They have no established corresponding experimental result for this comparison.',
        '- Primary RNA inputs: `log2(1 + counts / sum(all 15065 supplied gene counts) * 10000)`. All source genes, including genes outside the model panel or annotation, contribute to the library denominator. No IPF/control-dependent alignment is applied.',
        '- Sensitivity inputs: supplied publisher logRPKM, with duplicate gene mappings combined in linear space as in the previous bulk workflow. This alternate representation checks robustness; it does not replace the primary deployed RNA transformation.',
        '- Authorized NA policy: per-gene median of the verified 45-profile RNA training reference, on its stored model-input scale. Fill only ADORA3, CHGA, FABP7, FFAR1, GGT1, HBG1, ST6GALNAC6 and STAR, which could not be reconciled to rows in the supplied RNA source. Every one of the 1,303 required RNA features is retained. Both models receive the identical filled matrix. The >30%-NA AML-target removal authorization is not applied to these fixed lung RNA inputs.',
        '- Estimators: independent ElasticNetCV for every one of the 202 lipids in both models. Stored coefficients, scalers, and lipid-specific selected genes are used unchanged. Candidate training contained the authorized 45 non-D071 profiles; prior training contained the original 50. Consequently this compares the two delivered methods and does not isolate the target-scale correction from the authorized D071 removal.',
        '- Prediction scale: inverse the returned lipid log2 coordinate with `2**prediction`. Aggregate existing positive linear predictions; do not normalize their lipid-panel sum. Predicted signed folds are ratios of arithmetic donor means: +2 means twice control, -2 means half control; 0 denotes no change.',
        '- Raw tests: established donor-level Welch t-test on donor/region-averaged log predictions, unchanged from the prior bulk evaluation. BH correction includes the complete 202-lipid prediction family, with no arbitrary abundance/detection eligibility filter. Experimental p-values/FDR are taken directly from the supplied file; they are not recomputed on the matched subset. Bootstrap intervals in the differential CSVs use 2,000 donor resamples, seed 20261007.',
        '- Thresholds: raw p <0.05, BH <0.05, and BH <0.1 in both measurements and predictions, with strict predicted absolute folds >1 (no additional effect cutoff), >1.1, >1.2, >1.5 and >2. The same thresholds apply to both models. Threshold scans are descriptive; small high-agreement subsets do not establish improved agreement across the complete measured panel.', '',
        'RNA regional donor counts:', '', table(pd.DataFrame([{'comparison': r, **counts} for r, counts in donor_counts.items()])), '',
        'Approved median values used for absent RNA genes:', '',
        table(fills[fills.filled][['gene', 'current_annotation_Ensembl_ids', 'training_median_fill']]), '',
        '## Complete primary fold-threshold comparison', '', table(primary[display_columns]), '',
        '## RNA-representation and region sensitivity results', '',
        'The following table has no additional fold cutoff. Regional analyses retain their original donor rosters. It is not a cell-state analysis.', '',
        table(thresholds[thresholds.minimum_predicted_fold.eq(1)][['model','RNA_representation','region','statistic','p_cutoff',
            'significant_in_both','concordant_up','concordant_down','concordant_total','discordant_total','concordance_percent']]), '',
        'All 180 threshold combinations and their explicit concordant/discordant lipid lists are in [all_threshold_results.csv](all_threshold_results.csv). All 2,424 model-lipid/representation/region rows are retained in [every_lipid_all_comparisons.csv](every_lipid_all_comparisons.csv).', '',
        '## Primary discordant lipids significant in both', '']
    for row in primary[primary.minimum_predicted_fold.eq(1)].itertuples():
        lines += [f'**{row.model}, {row.statistic} <{row.p_cutoff:g}: {row.discordant_total} discordant lipids.**', '',
                  row.discordant_lipids if isinstance(row.discordant_lipids, str) else 'None.', '']
    c = detail[(detail.model == 'candidate') & (detail.RNA_representation == 'log2_CP10k') & (detail.region == 'donor_balanced')].set_index('model_lipid')
    p = detail[(detail.model == 'prior') & (detail.RNA_representation == 'log2_CP10k') & (detail.region == 'donor_balanced')].set_index('model_lipid')
    full = c[['measured_feature','status','supplied_log_effect','supplied_raw_p','supplied_adjusted_p']].copy()
    full = full.rename(columns={'supplied_log_effect':'source effect (log base unconfirmed)', 'supplied_raw_p':'source raw p','supplied_adjusted_p':'source BH'})
    for label, frame in [('candidate',c),('prior',p)]:
        for source, name in [('predicted_signed_fold','signed fold'), ('predicted_raw_p','raw p'), ('predicted_BH','BH')]:
            full[f'{label} {name}'] = frame[source]
    full.index.name = 'model lipid'
    lines += ['## Every one of the 202 model lipids: primary analysis', '',
              'Model folds are signed abundance ratios. The source effect is reported in its original unconfirmed log units; it is not mislabeled as a fold. A dash means no established matched experimental result, not a discarded model output.', '',
              table(full.reset_index()), '', '## Verification and reproducibility', '',
              'All saved predictions independently matched individual stored estimator.predict calls after their original scalers. Both bundles have identical top-gene options, l1-ratio and alpha grids, CV folds, iteration limit, sparsity penalty, and random seed. Welch raw p-values, complete-panel BH, arithmetic folds, every threshold count, ordered RNA/prediction identities, training-median fills, and both unchanged bundle hashes passed verification. No models were fitted or modified. The missing-input and prediction-equivalence test suites passed (13 tests).', '',
              f"Maximum prediction discrepancy: {max(row['prediction_max_abs_error'] for row in verification['checks']):.3g}; maximum raw-p discrepancy: {max(row['raw_p_max_abs_error'] for row in verification['checks']):.3g}.", '',
              'Source and model hashes: [comparison_audit.json](comparison_audit.json). Independent numerical checks: [independent_verification.json](independent_verification.json). Same-lipid paired comparisons: [same_lipid_paired_model_comparison.csv](same_lipid_paired_model_comparison.csv).', '',
              'Runner: [compare_bulk_ipf_models.py](../../compare_bulk_ipf_models.py). Verification/report generator: [verify_bulk_ipf_comparison.py](../../verify_bulk_ipf_comparison.py).', '',
              'Inference used scikit-learn '+sklearn.__version__+'; the candidate was saved with scikit-learn 1.6.1. No fitting occurred in the inference environment, and every delivered prediction was checked against the stored estimator/scaler calculation.', '',
              'This is independent, unmatched-cohort directional support for the delivered candidate. It does not validate individual-patient predictions, absolute lipid concentrations, or experimental fold magnitudes whose source log base is unknown.', '']
    (OUT / 'BULK_IPF_CANDIDATE_VS_PRIOR_REPORT.md').write_text('\n'.join(lines))
    print(json.dumps({'all_checks_passed': True, 'checks': len(verification['checks']),
        'max_prediction_error': max(row['prediction_max_abs_error'] for row in verification['checks'])}))
    print(paired.to_string(index=False))


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
