"""Fresh bulk-only rerun: supplied Abby model, last CP10k, confirmed CPM.

Retain the approved fixed models, full panels, biological replicates and tests.
No reproduction claim for Abby's missing historical saved outputs.
"""
from pathlib import Path
import hashlib
import json
import zipfile

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from .compare_bulk_ipf_models import main as rerun, sha
from .verify_bulk_ipf_comparison import main as verify, table

HERE = Path(__file__).resolve().parent
OUT = HERE / 'artifacts/bulk_IPF_three_version_comparison_20261008'
RAW = OUT / 'full_panel_rerun'
PREVIOUS = HERE / 'artifacts/bulk_IPF_full202_comparison_20261007'
LAST_CPM = HERE / 'artifacts/bulk_IPF_CPM_harmonized_20261008'
VERSIONS = [('Abby supplied model / CPM', 'prior', 'log2_CPM'),
            ('Last corrected model / CP10k', 'candidate', 'log2_CP10k'),
            ('New corrected model / CPM', 'candidate', 'log2_CPM')]
REQUEST = 'Perform the IPF bulk analyses again and compare to the prior set of Abby provided, your last version and this new version (bulk only lipid and RNA with different thresholds).'


def select_versions(frame):
    rows = []
    for name, model, representation in VERSIONS:
        item = frame[frame.model.eq(model) & frame.RNA_representation.eq(representation) & frame.region.eq('donor_balanced')].copy()
        item.insert(0, 'version', name)
        rows.append(item)
    return pd.concat(rows, ignore_index=True)


def frozen_reference_comparison(detail):
    """Same measured set at every cutoff; predicted p never changes denominator."""
    rows = []
    for statistic, cutoff in [('raw_p', .05), ('BH', .05), ('BH', .1)]:
        source_field = 'supplied_raw_p' if statistic == 'raw_p' else 'supplied_adjusted_p'
        for name, frame in detail.groupby('version', sort=False):
            eligible = frame.status.eq('matched') & frame[source_field].lt(cutoff) & frame.supplied_log_effect.ne(0)
            for estimand, effect in [('mean_log_prediction', 'predicted_geometric_log2FC'), ('arithmetic_abundance', 'predicted_log2FC')]:
                direction = frame.supplied_log_effect * frame[effect]
                up = eligible & direction.gt(0) & frame.supplied_log_effect.gt(0)
                down = eligible & direction.gt(0) & frame.supplied_log_effect.lt(0)
                rows.append({'version': name, 'statistic': statistic, 'p_cutoff': cutoff, 'effect_estimand': estimand,
                    'measured_significant': int(eligible.sum()), 'concordant_up': int(up.sum()), 'concordant_down': int(down.sum()),
                    'concordant_total': int((eligible & direction.gt(0)).sum()),
                    'discordant_total': int((eligible & direction.lt(0)).sum()), 'zero_effect': int((eligible & direction.eq(0)).sum()),
                    'concordance_percent': float(100 * direction[eligible].gt(0).mean())})
    return pd.DataFrame(rows)


def common_three_comparison(detail):
    """Require the same measured lipid to pass all three versions at each setting."""
    frames = {name: frame.set_index('model_lipid') for name, frame in detail.groupby('version', sort=False)}
    names = list(frames)
    base = frames[names[0]]
    for frame in frames.values():
        assert list(frame.index) == list(base.index) and list(frame.measured_feature.fillna('')) == list(base.measured_feature.fillna(''))
    rows, records = [], []
    for statistic, cutoff in [('raw_p', .05), ('BH', .05), ('BH', .1)]:
        source = 'supplied_raw_p' if statistic == 'raw_p' else 'supplied_adjusted_p'
        pp = 'predicted_raw_p' if statistic == 'raw_p' else 'predicted_BH'
        for fold in [1., 1.1, 1.2, 1.5, 2.]:
            common = base.status.eq('matched') & base[source].lt(cutoff) & base.supplied_log_effect.ne(0)
            for frame in frames.values():
                common &= frame[pp].lt(cutoff) & frame.predicted_signed_fold.abs().gt(fold)
            for name, frame in frames.items():
                direction = frame.supplied_log_effect * frame.predicted_signed_fold
                agree = common & direction.gt(0)
                rows.append({'version': name, 'statistic': statistic, 'p_cutoff': cutoff, 'minimum_predicted_fold': fold,
                    'identical_lipid_denominator': int(common.sum()), 'concordant_up': int((agree & frame.supplied_log_effect.gt(0)).sum()),
                    'concordant_down': int((agree & frame.supplied_log_effect.lt(0)).sum()), 'concordant_total': int(agree.sum()),
                    'discordant_total': int((common & direction.lt(0)).sum()),
                    'concordance_percent': float(100 * agree.sum() / common.sum()) if common.sum() else np.nan})
                selected = frame.loc[common].copy()
                selected['statistic'], selected['p_cutoff'], selected['minimum_predicted_fold'] = statistic, cutoff, fold
                records.append(selected.reset_index())
    return pd.DataFrame(rows), pd.concat(records, ignore_index=True)


def preservation_checks():
    checks = []
    # Do not rewrite previous directories. Verify all 202 rows and all 139 profiles.
    for original, representations in [(PREVIOUS, ['log2_CP10k', 'publisher_logRPKM']),
                                      (LAST_CPM, ['log2_CPM', 'log2_CP10k', 'publisher_logRPKM'])]:
        for model in ['candidate', 'prior']:
            for representation in representations:
                for suffix in ['all_139_predictions_log2', 'donor_balanced_all_202_differentials']:
                    filename = f'{model}_{representation}_{suffix}.csv'
                    old, new = [pd.read_csv(root / filename, index_col=0) for root in [original, RAW]]
                    assert list(old.index) == list(new.index) and list(old.columns) == list(new.columns)
                    np.testing.assert_allclose(old, new, atol=1e-10, rtol=1e-10, equal_nan=True)
                    checks.append({'original': str(original), 'file': filename,
                                   'max_abs_difference': float(np.nanmax(abs(old.to_numpy() - new.to_numpy())))})
    return checks


def main():
    state = json.loads((HERE / 'integrity/decision_state.json').read_text())
    authorization = state['bulk_IPF_three_version_rerun_20261008']
    if authorization.get('granted_by') != 'user' or authorization.get('user_message') != REQUEST:
        raise ValueError('Three-version rerun requires the actual user request')
    Abby = json.loads((HERE / 'artifacts/Abby_latest_bundle_comparison_20261008/comparison_audit.json').read_text())
    if sha(Abby['archive']) != Abby['archive_sha256']:
        raise ValueError('Abby archive changed; source must be reconciled')
    with zipfile.ZipFile(Abby['archive']) as archive:
        payload = archive.read(Abby['supplied_bundle_member'])
    if hashlib.sha256(payload).hexdigest() != Abby['supplied_bundle_sha256'] or sha(HERE / 'rna2lipid_hs_lung_lipidwise_bundle.pkl') != Abby['supplied_bundle_sha256']:
        raise ValueError('Prior comparator is not exactly the supplied Abby model')
    rerun(RAW)
    verify(RAW)
    preservation = preservation_checks()
    detail = select_versions(pd.read_csv(RAW / 'every_lipid_all_comparisons.csv'))
    thresholds = select_versions(pd.read_csv(RAW / 'all_threshold_results.csv'))
    refs = select_versions(pd.read_csv(RAW / 'reference_only_concordance.csv'))
    reference_thresholds = frozen_reference_comparison(detail)
    common, common_details = common_three_comparison(detail)
    top_detail = select_versions(pd.read_csv(RAW / 'reference_only_every_lipid.csv'))
    top_detail = top_detail[top_detail.selection.eq('top15_up_top15_down')]
    assert len(top_detail) == 90
    top_detail.to_csv(OUT / 'top30_all_versions_every_lipid.csv', index=False)
    assert detail.model_lipid.nunique() == 202 and len(detail) == 606 and len(thresholds) == 45
    for label, frame in [('three_versions_all_202_lipids',detail),('three_versions_all_thresholds', thresholds),
                         ('three_versions_reference_only',refs),('reference_only_pvalue_thresholds',reference_thresholds),
                         ('same_lipids_significant_in_all_three',common),('same_lipids_all_three_details',common_details)]:
        frame.to_csv(OUT / f'{label}.csv', index=False)
    checkpoint = pd.read_csv(PREVIOUS / 'all_threshold_results.csv')
    checkpoint = checkpoint[checkpoint.RNA_representation.eq('log2_CP10k') & checkpoint.region.eq('donor_balanced')]
    checkpoint.to_csv(OUT / 'historical_CP10k_132_match_thresholds.csv', index=False)
    write_report(detail, thresholds, refs, reference_thresholds, common, checkpoint, top_detail)
    audit = {'authorization':authorization, 'Abby_supplied_model_verified_byte_identical':True,
        'model_versions':[{'version':n,'model':m,'RNA_representation':r} for n,m,r in VERSIONS],
        'full_panel_rerun_directory':str(RAW), 'prior_result_preservation':preservation,
        'all_checks_passed':True, 'raw_estimator_statistical_verification':json.loads((RAW/'independent_verification.json').read_text()),
        'RNA_IPF_donors':20,'RNA_control_donors':14,'RNA_input_genes':1303,'lipid_outputs':202,
        'harmonized_measured_matches':134,'threshold_rows':45,
        'Abby_original_historical_saved_results_available':False,'training_performed':False,'deployment_changed':False,
        'script_sha256':sha(__file__)}
    (OUT/'THREE_VERSION_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(thresholds[thresholds.minimum_predicted_fold.eq(1)][['version','statistic','p_cutoff','significant_in_both','concordant_total','discordant_total','concordance_percent']].to_string(index=False),flush=True)


def write_report(detail, thresholds, refs, reference_thresholds, common, historical, top_detail):
    display=['version','statistic','p_cutoff','minimum_predicted_fold','significant_in_both','concordant_up','concordant_down','concordant_total','discordant_total','concordance_percent']
    top=refs[refs.selection.eq('top15_up_top15_down') & refs.effect_estimand.eq('mean_log_prediction')]
    lines=['# Bulk IPF: Abby supplied model, last CODEX version, new CODEX version','','Fresh rerun: 2026-10-08. Bulk RNA and lung-tissue lipidomics only. Fixed models; no retraining or deployment.','',
        'The new CPM-input version agrees with 18/30 (60%) of the fixed top measured lipids versus 15/30 (50%) for the last CP10k-input version and 12/30 (40%) for the supplied Abby model rerun with the same donor-level CPM analysis. For lipids significant in both at BH <0.1, without an additional fold cutoff, these are 48/70 (68.6%), 40/63 (63.5%), and 23/61 (37.7%), respectively. The overlap denominators differ, so identical-lipid comparisons are also provided below.','',
        '## What is being compared','',table(pd.DataFrame([{'version':name,'fixed_model':'Original 50-profile lipid-wise ElasticNetCV' if model=='prior' else 'Corrected 45-profile native-relative-log2 ElasticNetCV','RNA':'log2(CPM + 1)' if rep=='log2_CPM' else 'log2(CP10k + 1)'} for name,model,rep in VERSIONS])),'',
        'Abby’s supplied latest pickle is byte-identical to the preserved prior model. The last/new CODEX versions use the same corrected pickle and weights; their only input difference is the RNA normalization multiplier. All other bulk analysis steps are identical. This comparison does not substitute another estimator or fit new weights.','',
        '**Abby’s original approximate 50% is author-reported, not a verified result column.** The original ZIP contains no saved result tables. Her supplied historical script used 88 profiles (52 IPF/36 controls), profile-level Mann–Whitney testing, sample-specific median RNA filling, and lowest-FDR duplicate lipid selection. The common-method supplied-model rerun above uses the approved donor-level workflow and should not be described as reproduction of that original script’s execution.','',
        '## Source-only top-30 directional concordance','',
        'Select up to 15 measured increases and 15 measured decreases at source BH <0.05, ranked only by source BH and source effect. Prediction significance is not required. All three versions use exactly the same 30 lipids. The table uses mean-log prediction direction, closest to Abby’s stated effect metric, after the established donor/region balancing.','',table(top[['version','measured_lipids','concordant_up','concordant_down','concordant_total','discordant_total','concordance_percent']]),'',
        '## Full measured-significant panels: no prediction-significance requirement','',table(reference_thresholds),'',
        '## Both significant: no added predicted fold cutoff','',table(thresholds[thresholds.minimum_predicted_fold.eq(1)][display]),'',
        '## All requested fold and significance thresholds','',
        'Raw p <0.05, BH <0.05 and BH <0.1 in both measured and predicted data; strict predicted absolute folds >1 (no additional effect cutoff), >1.1, >1.2, >1.5 and >2. Experimental effect log base is unconfirmed, so no experimental fold-magnitude cutoff is invented. All predicted folds are signed linear ratios: +2 means twice control; −2 means half control.','',table(thresholds[display]),'',
        '## Identical lipids passing all three versions','',
        'Every row uses the same lipid set across the three versions: source significance, prediction significance and predicted fold threshold must pass in all three. An empty subset has no interpretable concordance percentage.','',table(common),'',
        '## Previous saved results versus harmonized matching','',
        'The previous CP10k report used 132 matched experimental features. This fresh comparison preserves those identities and adds the approved CE(18:2)_POS and CE(20:3)_POS corrections, giving 134 matches for every version. The original predictions, all-202 p-values and folds reproduce unchanged. CE(20:3) adds one concordant significant lipid to the corrected CP10k version. Thus its old BH <0.1 value of 39/62 becomes 40/63 after matching correction; this change is not caused by a model refit. Previous result directories remain unchanged.','',
        table(historical[historical.minimum_predicted_fold.eq(1)][['model','statistic','p_cutoff','significant_in_both','concordant_total','discordant_total','concordance_percent']]),'',
        '## Source data and retained methods','',
        '- RNA: supplied GSE213001 counts, 15,065 genes × 139 profiles. Impute and retain every profile for every model/input representation; use 101 known-region IPF/control profiles in the established overall anatomical aggregation, yielding 20 IPF and 14 control donors. The two unknown-region profiles and ILD/CLAD profiles remain in full prediction exports. Donors are independent IPF/control groups, not paired patients.',
        '- Lipids: Dr. Clair’s `10_results_with_statistics.csv`, 544 measured features and 40 tissue-sample columns. Use supplied experimental raw p, BH and effect sign directly; no re-estimation of source statistics or unsupported source log-base assumption.',
        '- RNA library normalization uses all 15,065 source genes; aggregate duplicate mapped-gene counts before logging. CPM uses 1,000,000, CP10k uses 10,000, both with log2 and +1. Retain every one of 1,303 required RNA genes. Fill the approved eight unresolved inputs with gene-specific medians of the verified 45-profile training RNA table, identically in all versions.',
        '- Use original stored RNA scaler, lipid-specific genes, ElasticNetCV weights and inverse lipid target scaler for each model. Use `2**prediction` to obtain positive linear values without subtracting one. The prior output is its legacy processed coordinate, whereas the corrected output is native relative log2; this difference is explicit, not converted into an absolute concentration claim.',
        '- Average left/right within donor/anatomical region, then give available apex/base regions equal donor weight. Fold is the ratio of arithmetic group means of the balanced linear predictions. Welch raw tests use balanced donor log predictions; the source-only tables additionally show mean-log and arithmetic-abundance directions separately.',
        '- BH tests the full 202-lipid prediction family in every version, with no detection/abundance filter. All 202 outputs stay in the tables even where no unique corresponding experimental feature is established. No lipid-panel total normalization, disease-dependent alignment or new model fitting.',
        '- Threshold scans are descriptive. An increase in concordance at a selected threshold can reflect a smaller subset. Independent, unmatched bulk measurements support directional agreement; they do not validate patient-level abundance or unknown measured fold magnitudes.','',
        '## Every lipid in the primary bulk comparison','']
    for name,frame in detail.groupby('version',sort=False):
        lines += [f'### {name}','',table(frame[['model_lipid','measured_feature','status','supplied_log_effect','supplied_raw_p','supplied_adjusted_p','predicted_signed_fold','predicted_raw_p','predicted_BH']]),'']
    lines += ['## Explicit concordant and discordant lipid lists at every threshold','']
    for row in thresholds.itertuples():
        lines += [f'**{row.version}; {row.statistic} <{row.p_cutoff:g}; predicted absolute fold >{row.minimum_predicted_fold:g}**','',
            f'Concordant ({row.concordant_total}): '+(row.concordant_lipids if isinstance(row.concordant_lipids,str) else 'None.'),'',
            f'Discordant ({row.discordant_total}): '+(row.discordant_lipids if isinstance(row.discordant_lipids,str) else 'None.'),'']
    lines += ['## Verification and reproducibility','',
        'All saved predictions were independently recalculated from every stored estimator/scaler; donor aggregation, Welch p-values, complete-panel BH, arithmetic folds and threshold counts passed independent verification. All saved original CP10k and recent CPM predictions and overall differentials matched the fresh run within numerical tolerance. Both model hashes are unchanged.','',
        'Data: [all 202 lipids in all three versions](three_versions_all_202_lipids.csv), [all 45 threshold rows](three_versions_all_thresholds.csv), [same-lipid sets](same_lipids_significant_in_all_three.csv), [reference-only results](three_versions_reference_only.csv). Verification: [audit](THREE_VERSION_AUDIT.json) and [raw-estimator checks](full_panel_rerun/independent_verification.json).','',
        'Reproduce using `python -m altanalyze3.components.rna2lipid.compare_three_bulk_ipf_versions` from the repository environment. This invokes the existing source, authorization and identity gates.','']
    top_lines = ['## The exact 30 reference-selected lipids', '',
        'Predicted folds are signed arithmetic-abundance ratios. Agreement flags distinguish mean-log and arithmetic direction; predicted p-values do not determine this source-only selection.', '']
    for name, frame in top_detail.groupby('version', sort=False):
        top_lines += ['### '+name, '', table(frame[['model_lipid','measured_feature','supplied_log_effect','supplied_adjusted_p','predicted_signed_fold','predicted_raw_p','predicted_BH','mean_log_concordant','arithmetic_abundance_concordant']]), '']
    position = lines.index('## Full measured-significant panels: no prediction-significance requirement')
    lines[position:position] = top_lines
    (OUT/'BULK_IPF_THREE_VERSION_REPORT.md').write_text('\n'.join(lines))


if __name__=='__main__':
    with threadpool_limits(limits=1):
        main()
