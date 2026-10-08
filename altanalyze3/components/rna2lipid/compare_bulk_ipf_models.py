"""Authorized fixed-model, full-panel comparison on unmatched bulk IPF cohorts.

No fitting, feature intersection, outcome-directed calibration, model replacement,
or experimental log-base assumption. All input profiles and model outputs saved.
"""
from pathlib import Path
import hashlib
import json
import pickle

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

from .api import load_bundle
from .release import release_manifest, verify_release_bundle
from .ipf_utils import metadata, donor_means, fold_table

HERE = Path(__file__).resolve().parent
SOURCE = Path('/Users/saljh8/Downloads/Lipidomics/IPF')
PREVIOUS = HERE / 'artifacts/IPF_candidate_validation_20261001'
REFERENCE = HERE / 'artifacts/LungMAP_full202_native_log2_20261005'
OUT = HERE / 'artifacts/bulk_IPF_CPM_harmonized_20261008'
ANNOTATION = HERE.parent / 'fastCNV/resources/Hs_Ensembl_GRCh38_genes.tsv'
GENES_TO_FILL = ['ADORA3', 'CHGA', 'FABP7', 'FFAR1', 'GGT1', 'HBG1', 'ST6GALNAC6', 'STAR']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def signed_fold(log2fc):
    values = np.asarray(log2fc, float)
    return np.where(values > 0, np.exp2(values), np.where(values < 0, -np.exp2(-values), 0.))


def model_inputs(raw, published, annotation, genes, medians):
    """Keep all source rows in the CPM/CP10k denominators and all declared model genes."""
    if raw.index.has_duplicates or raw.columns.has_duplicates or published.index.has_duplicates:
        raise ValueError('Duplicate source gene/sample identities')
    if list(raw.index) != list(published.index) or list(raw.columns) != list(published.columns):
        raise ValueError('Raw/published RNA identity mismatch')
    if annotation.gene_id.duplicated().any():
        raise ValueError('Ambiguous supplied Ensembl annotation')
    if not np.isfinite(raw.to_numpy()).all() or (raw.to_numpy() < 0).any():
        raise ValueError('Invalid source counts')
    if not np.isfinite(published.to_numpy()).all():
        raise ValueError('Nonfinite publisher RNA values')
    totals = raw.sum(axis=0)
    if (totals <= 0).any():
        raise ValueError('Zero source library; no sample exclusion permitted')
    symbol = annotation.set_index('gene_id').gene.reindex(raw.index)
    matrices = {name: pd.DataFrame(index=raw.columns, columns=genes, dtype=float)
                for name in ('log2_CPM', 'log2_CP10k', 'publisher_logRPKM')}
    audit = []
    missing_genes = []
    for gene in genes:
        ids = raw.index[symbol.eq(gene)]
        if len(ids):
            counts = raw.loc[ids].sum(axis=0)
            matrices['log2_CPM'][gene] = np.log2(1 + counts / totals * 1e6)
            matrices['log2_CP10k'][gene] = np.log2(1 + counts / totals * 1e4)
            matrices['publisher_logRPKM'][gene] = np.log2(np.exp2(published.loc[ids]).sum(axis=0))
        else:
            missing_genes.append(gene)
            value = float(medians.loc[gene])
            if not np.isfinite(value):
                raise ValueError(f'No finite training median for {gene}')
            for matrix in matrices.values():
                matrix[gene] = value
        audit.append({'gene': gene, 'source_Ensembl_ids': ';'.join(ids),
                      'current_annotation_Ensembl_ids': ';'.join(annotation.loc[annotation.gene.eq(gene), 'gene_id']),
                      'source_match_count': len(ids), 'filled': not len(ids),
                      'training_median_fill': float(medians.loc[gene]) if not len(ids) else np.nan})
    if missing_genes != [g for g in genes if g in GENES_TO_FILL]:
        raise ValueError(f'New unresolved required RNA identities: {missing_genes}')
    for matrix in matrices.values():
        if list(matrix.columns) != list(genes) or not np.isfinite(matrix.to_numpy()).all():
            raise ValueError('Incomplete/nonfinite model input matrix')
    return matrices, pd.DataFrame(audit)


def thresholds(comparison, model, representation, region):
    valid = comparison.status.eq('matched')
    direction = comparison.supplied_log_effect * comparison.predicted_log2FC
    rows = []
    for statistic, cutoff in [('raw_p', .05), ('BH', .05), ('BH', .1)]:
        source_p = comparison.supplied_raw_p if statistic == 'raw_p' else comparison.supplied_adjusted_p
        predicted_p = comparison.predicted_raw_p if statistic == 'raw_p' else comparison.predicted_BH
        eligible_source = valid & source_p.lt(cutoff) & comparison.supplied_log_effect.ne(0)
        for fold in [1., 1.1, 1.2, 1.5, 2.]:
            overlap = eligible_source & predicted_p.lt(cutoff) & comparison.predicted_log2FC.abs().gt(np.log2(fold))
            concordant = overlap & direction.gt(0)
            discordant = overlap & direction.lt(0)
            n = int(overlap.sum())
            rows.append({'model': model, 'RNA_representation': representation, 'region': region,
                'statistic': statistic, 'p_cutoff': cutoff, 'minimum_predicted_fold': fold,
                'matched_experimental_significant': int(eligible_source.sum()), 'significant_in_both': n,
                'concordant_up': int((concordant & comparison.supplied_log_effect.gt(0)).sum()),
                'concordant_down': int((concordant & comparison.supplied_log_effect.lt(0)).sum()),
                'concordant_total': int(concordant.sum()), 'discordant_total': int(discordant.sum()),
                'concordance_percent': 100 * concordant.sum() / n if n else np.nan,
                'concordant_lipids': '; '.join(comparison.loc[concordant, 'measured_feature']),
                'discordant_lipids': '; '.join(comparison.loc[discordant, 'measured_feature'])})
    return rows



def verified_RNA_source(state):
    """Reject inference until the actual source answer and pinned code are verified."""
    q = next(q for q in state['questions'] if q['id'] == 'lung_RNA_training_matrix_preprocessing_20261007')
    expected = 'unmerged_clair_newRNA.csv = newrna_cell_clair_filtered_symbol.csv'
    if q['status'] != 'resolved' or q.get('resolved_by') != 'user' or q.get('user_message') != expected or not q.get('new_inference_allowed'):
        raise ValueError('Training RNA source identity unresolved; stop inference')
    source = q['provided_preprocessing_script']
    if sha(HERE.parents[1] / source['path']) != source['sha256']:
        raise ValueError('Verified RNA preprocessing source changed')
    if q['verified_script_normalization'] != 'log2(1 + df / column_sum(df) * 1000000)':
        raise ValueError('Verified training RNA transformation changed')
    return q


def corrected_matching(base, lipid_raw):
    """Correct the approved CE parser omission; preserve every previous identity."""
    fixed = base.copy()
    changes = []
    for gene in fixed.index[fixed.index.str.startswith('CE(')]:
        fixed.loc[gene, 'canonical'] = gene
        ion = fixed.loc[gene, 'ion_mode']
        suffix = {'P': '_POS', 'N': '_NEG'}.get(ion)
        candidates = lipid_raw.index[lipid_raw.Lipid.eq(gene) & lipid_raw.index.str.endswith(suffix)] if suffix else []
        if len(candidates) > 1:
            raise ValueError(f'Ambiguous CE mapping requires discussion: {gene}')
        if len(candidates) == 1:
            fixed.loc[gene, 'measured_feature'] = candidates[0]
            fixed.loc[gene, 'status'] = 'matched'
            changes.append({'model_lipid': gene, 'measured_feature': candidates[0], 'reason': 'CE omitted from original class parser'})
    pd.testing.assert_frame_equal(fixed.loc[base.status.eq('matched')], base.loc[base.status.eq('matched')])
    if {r['model_lipid'] for r in changes} != {'CE(18:2)', 'CE(20:3)'}:
        raise ValueError('CE correction differs from reviewed source; discuss before inference')
    if fixed.index.has_duplicates or fixed.loc[fixed.status.eq('matched'), 'measured_feature'].duplicated().any():
        raise ValueError('Nonunique frozen lipid mapping')
    return fixed, changes


def reference_only(comparison, model, representation, region):
    """Selection uses measured BH/effect only, never model p-values or outcome."""
    eligible = comparison[comparison.status.eq('matched') & comparison.supplied_adjusted_p.lt(.05) & comparison.supplied_log_effect.ne(0)]
    top = pd.concat([eligible[eligible.supplied_log_effect.gt(0)].sort_values(
        ['supplied_adjusted_p', 'supplied_log_effect'], ascending=[True, False]).head(15),
        eligible[eligible.supplied_log_effect.lt(0)].sort_values(
        ['supplied_adjusted_p', 'supplied_log_effect'], ascending=[True, True]).head(15)])
    summaries, details = [], []
    for selection, frame in [('top15_up_top15_down', top), ('all_measured_BH005', eligible)]:
        for estimand, field in [('mean_log_prediction', 'predicted_geometric_log2FC'), ('arithmetic_abundance', 'predicted_log2FC')]:
            agree = frame.supplied_log_effect * frame[field] > 0
            summaries.append({'model': model, 'RNA_representation': representation, 'region': region,
                'selection': selection, 'effect_estimand': estimand, 'measured_lipids': len(frame),
                'measured_up': int(frame.supplied_log_effect.gt(0).sum()), 'measured_down': int(frame.supplied_log_effect.lt(0).sum()),
                'concordant_up': int((agree & frame.supplied_log_effect.gt(0)).sum()),
                'concordant_down': int((agree & frame.supplied_log_effect.lt(0)).sum()),
                'concordant_total': int(agree.sum()), 'discordant_total': int((frame.supplied_log_effect * frame[field] < 0).sum()),
                'zero_predicted_effect': int(frame[field].eq(0).sum()),
                'concordance_percent': 100 * agree.mean() if len(frame) else np.nan})
        selected = frame.copy()
        selected['selection'] = selection
        selected['mean_log_concordant'] = selected.supplied_log_effect * selected.predicted_geometric_log2FC > 0
        selected['arithmetic_abundance_concordant'] = selected.supplied_log_effect * selected.predicted_log2FC > 0
        details.append(selected)
    return summaries, details

def main(output=OUT):
    OUT = Path(output)
    state_path = HERE / 'integrity/decision_state.json'
    state = json.loads(state_path.read_text())
    authorization = state['bulk_IPF_harmonized_plan_20261008']
    RNA_source_confirmation = verified_RNA_source(state)
    if authorization.get('granted_by') != 'user' or authorization.get('user_message') != 'proceed with this plan':
        raise ValueError('Bulk inference requires actual user authorization')
    decision = next(q for q in state['questions'] if q['id'] == 'bulk_IPF_eight_RNA_gene_mappings_20261007')
    if decision['status'] != 'resolved' or not decision.get('NA_filling_authorized'):
        raise ValueError('Required RNA mapping/filling decision unresolved')
    manifest = release_manifest()
    old_result = json.loads((PREVIOUS / 'IPF_candidate_validation_result.json').read_text())
    reviewed_hashes = {str(Path(path).resolve()): digest for path, digest in old_result['source_sha256'].items()}
    paths = [SOURCE / f for f in ['GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-raw_counts.csv.gz',
        'GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-logRPKMs-normalised.csv.gz',
        'GSE213001_series_matrix.txt.gz', '10_results_with_statistics.csv']] + [ANNOTATION]
    source_hashes = {}
    for path in paths:
        actual = sha(path)
        if actual != reviewed_hashes[str(path.resolve())]:
            raise ValueError(f'Previously reviewed source changed: {path}')
        source_hashes[str(path)] = actual
    raw = pd.read_csv(paths[0], index_col=0)
    published = pd.read_csv(paths[1], index_col=0)
    md = metadata(paths[2])
    if set(md.index) != set(raw.columns) or len(md) != 139:
        raise ValueError('Complete GEO sample roster mismatch')
    raw, published = raw.loc[:, md.index], published.loc[:, md.index]
    annotation = pd.read_csv(ANNOTATION, sep='\t')
    run_manifest = json.loads((REFERENCE / 'run_manifest.json').read_text())
    RNA_source = run_manifest['inputs']['RNA']
    if sha(RNA_source['path']) != RNA_source['sha256']:
        raise ValueError('Approved original training RNA source changed')
    training = pd.read_csv(RNA_source['path'], index_col=0).T
    training = training.T.groupby(level=0).mean().T
    training.index = training.index.str.strip()
    training.columns = training.columns.str.strip()
    training = training.loc[manifest['training_samples'], manifest['X_columns']]
    training = training.fillna(training.median())
    saved_training = pd.read_csv(REFERENCE / 'candidate_training_RNA.csv', index_col=0)
    if list(training.index) != list(saved_training.index) or list(training.columns) != list(saved_training.columns):
        raise ValueError('45-profile training-reference identities differ')
    np.testing.assert_allclose(training, saved_training, rtol=0, atol=1e-12)
    matrices, RNA_audit = model_inputs(raw, published, annotation, manifest['X_columns'], training.median())
    measured = pd.read_csv(REFERENCE / 'experimental_lung_IPF_lipids_all_544.csv', index_col=0)
    lipid_raw = pd.read_csv(paths[3], index_col=0)
    if list(measured.index) != list(lipid_raw.index) or len(measured) != 544:
        raise ValueError('Complete experimental lipid roster mismatch')
    np.testing.assert_array_equal(measured.supplied_log_effect, lipid_raw['log(IPF/Ctrl)'])
    np.testing.assert_array_equal(measured.supplied_raw_p, lipid_raw.IPF_vs_Ctrl_Ttest_p)
    np.testing.assert_array_equal(measured.supplied_adjusted_p, lipid_raw.IPF_vs_Ctrl_Ttest_padj)
    match = pd.read_csv(REFERENCE / 'lipid_matching_audit.csv')
    base_match = match[match.model.eq('candidate_native_log2_202')].set_index('gene')
    prior_match = match[match.model.eq('current_explicit_log2')].set_index('gene')
    if list(base_match.index) != manifest['Y_columns']:
        raise ValueError('Approved matching audit lacks required 202 identities')
    pd.testing.assert_frame_equal(base_match.drop(columns='model'), prior_match.drop(columns='model'))
    base_match, matching_changes = corrected_matching(base_match, lipid_raw)
    diagnoses = md.groupby('donorid').diseasegroup
    if (diagnoses.nunique() != 1).any():
        raise ValueError('Conflicting donor diagnoses')
    donor_groups = diagnoses.first()
    # Verify BOTH fixed algorithms before creating outputs or running either model.
    preflight_settings = None
    for label, filename, expected_hash in [('candidate', manifest['bundle_filename'], manifest['bundle_sha256']),
                                          ('prior', manifest['previous_bundle_filename'], manifest['previous_bundle_sha256'])]:
        path = HERE / filename
        if sha(path) != expected_hash:
            raise ValueError(f'{label} differs from the pinned fixed model')
        with path.open('rb') as handle:
            bundle = pickle.load(handle)
        if label == 'candidate':
            verify_release_bundle(path, bundle)
        if bundle['X_columns'] != manifest['X_columns'] or bundle['Y_columns'] != manifest['Y_columns']:
            raise ValueError(f'{label} full ordered panel mismatch')
        if list(bundle['models']) != manifest['Y_columns'] or {type(v['model']).__name__ for v in bundle['models'].values()} != {'ElasticNetCV'}:
            raise ValueError(f'{label} baseline estimators differ')
        settings = {key: bundle[key] for key in ['top_gene_options', 'l1_ratio_grid', 'alpha_grid', 'cv_folds', 'max_iter', 'sparsity_penalty', 'random_seed']}
        if preflight_settings is None:
            preflight_settings = settings
        else:
            for key in settings:
                np.testing.assert_array_equal(settings[key], preflight_settings[key])
    OUT.mkdir(parents=True, exist_ok=True)
    RNA_audit.to_csv(OUT / 'RNA_mapping_and_median_filling.csv', index=False)
    md.to_csv(OUT / 'all_139_RNA_sample_metadata.csv')
    base_match.to_csv(OUT / 'all_202_lipid_matching.csv')
    all_thresholds, all_comparisons, numerical = [], [], []
    reference_summaries, reference_details = [], []
    audit = {'authorization': authorization, 'RNA_source_confirmation': RNA_source_confirmation,
        'RNA_representations': list(matrices), 'matching_changes': matching_changes, 'source_sha256': source_hashes,
        'source_profiles': len(md), 'RNA_inputs': len(manifest['X_columns']), 'lipid_outputs': len(manifest['Y_columns']),
        'training_median_profiles': list(training.index), 'filled_RNA_genes': RNA_audit.loc[RNA_audit.filled, 'gene'].tolist(),
        'matched_experimental_features': int(base_match.measured_feature.notna().sum()),
        'experimental_log_base': 'unconfirmed; no experimental fold magnitude invented',
        'primary_RNA_representation': 'log2(1 + source_gene_counts / sum_all_15065_source_gene_counts * 1000000)',
        'sensitivity_RNA_representation': 'Previous CP10k and publisher logRPKM, unchanged except original duplicate-gene linear summation and authorized missing-gene medians',
        'study_offsets_added': False, 'training_executed': False, 'models_changed': False,
        'statistical_test': 'Established donor-level Welch test on averaged log predictions; BH over all 202 targets',
        'donor_aggregation': 'Established left/right mean within region, then equal apex/base weighting within donor; 101 known-region IPF/NDC profiles, all 34 IPF/NDC donors retained. All 139 profiles imputed and saved, including two unknown-region profiles and ILD/CLAD.',
        'reported_fold': 'Ratio of arithmetic means of positive 2**predictions, donor and region balanced; signed positive increase / negative decrease',
        'model_checks': {}}
    baseline_method = None
    for label, path, expected_hash in [('candidate', HERE / manifest['bundle_filename'], manifest['bundle_sha256']),
                                      ('prior', HERE / manifest['previous_bundle_filename'], manifest['previous_bundle_sha256'])]:
        if sha(path) != expected_hash:
            raise ValueError(f'{label} differs from approved bundle')
        with path.open('rb') as handle:
            original = pickle.load(handle)
        if label == 'candidate':
            verify_release_bundle(path, original)
        for key in ('X_columns', 'Y_columns'):
            if original[key] != manifest[key]:
                raise ValueError(f'{label} incomplete {key} roster')
        settings = {k: original[k] for k in ['top_gene_options', 'l1_ratio_grid', 'alpha_grid', 'cv_folds', 'max_iter', 'sparsity_penalty', 'random_seed']}
        if baseline_method is None:
            baseline_method = settings
        else:
            for k in settings:
                np.testing.assert_array_equal(settings[k], baseline_method[k])
        if {type(v['model']).__name__ for v in original['models'].values()} != {'ElasticNetCV'}:
            raise ValueError('Baseline estimator changed')
        model = load_bundle(path)
        audit['model_checks'][label] = {'path': str(path), 'sha256': expected_hash,
            'input_genes': len(model.input_genes), 'outputs': len(model.output_lipids),
            'architecture': model.architecture, 'target_scaling': model.target_scaling_mode}
        for representation, X in matrices.items():
            X.to_csv(OUT / f'{representation}_identical_RNA_inputs.csv')
            pred = model.predict_from_dataframe(X).predictions
            if list(pred.index) != list(md.index) or list(pred.columns) != manifest['Y_columns'] or not np.isfinite(pred.to_numpy()).all():
                raise ValueError('Incomplete/nonfinite predictions; do not drop/fill outputs')
            pred.to_csv(OUT / f'{label}_{representation}_all_139_predictions_log2.csv')
            linear = np.exp2(pred)
            if not np.isfinite(linear.to_numpy()).all() or not (linear.to_numpy() > 0).all():
                raise ValueError('Nonfinite/invertibility problem; do not clip predictions')
            for region in (None, 'Apex', 'Base'):
                region_name = region or 'donor_balanced'
                ld, logged = donor_means(linear, md, region), donor_means(pred, md, region)
                groups = donor_groups.loc[ld.index]
                # Region-specific cohorts have their own original rosters;
                # the overall 20/14 roster is not imposed on apex/base subsets.
                baseline_donors = pd.read_csv(PREVIOUS / 'supplied_bulk_MS1_47' /
                    f'publisher_logRPKM_{region_name}_donor_predictions_linear.csv', index_col=0).index
                baseline_counts = old_result['comparisons'][f'supplied_bulk_MS1_47/publisher_logRPKM/{region_name}']
                if list(ld.index) != list(baseline_donors):
                    raise ValueError(f'Original {region_name} donor identities/order changed')
                if groups.eq('IPF').sum() != baseline_counts['RNA_IPF_donors'] or groups.eq('NDC').sum() != baseline_counts['RNA_control_donors']:
                    raise ValueError(f'Original {region_name} donor diagnosis counts changed')
                stats = fold_table(ld, logged, groups, np.random.default_rng(20261007), 2000)
                stats['signed_fold'] = signed_fold(stats.log2FC)
                stats.to_csv(OUT / f'{label}_{representation}_{region_name}_all_202_differentials.csv')
                ld.to_csv(OUT / f'{label}_{representation}_{region_name}_donor_predictions_linear.csv')
                comparison = base_match.reset_index().rename(columns={'gene': 'model_lipid'})
                comparison['model'] = label
                comparison['RNA_representation'] = representation
                comparison['region'] = region_name
                for key, column in [('log2FC', 'predicted_log2FC'), ('signed_fold', 'predicted_signed_fold'),
                                    ('geometric_log2FC', 'predicted_geometric_log2FC'), ('pvalue', 'predicted_raw_p'), ('FDR', 'predicted_BH')]:
                    comparison[column] = comparison.model_lipid.map(stats[key])
                for column in ('supplied_log_effect', 'supplied_raw_p', 'supplied_adjusted_p'):
                    comparison[column] = comparison.measured_feature.map(measured[column])
                all_thresholds.extend(thresholds(comparison, label, representation, region_name))
                refs, ref_detail = reference_only(comparison, label, representation, region_name)
                reference_summaries.extend(refs)
                reference_details.extend(ref_detail)
                all_comparisons.append(comparison)
                matched = comparison[comparison.status.eq('matched')]
                concordant = matched.supplied_log_effect * matched.predicted_log2FC > 0
                numerical.append({'model': label, 'RNA_representation': representation, 'region': region_name,
                    'matched_lipids': len(matched), 'concordant_up_all': int((concordant & matched.supplied_log_effect.gt(0)).sum()),
                    'concordant_down_all': int((concordant & matched.supplied_log_effect.lt(0)).sum()),
                    'concordant_all': int(concordant.sum()), 'concordance_percent_all': 100 * concordant.mean(),
                    'Spearman_supplied_effect_vs_predicted_log2FC': float(spearmanr(matched.supplied_log_effect, matched.predicted_log2FC).statistic),
                    'median_absolute_predicted_fold_all202': float(np.median(np.exp2(stats.log2FC.abs()))),
                    'rawp005_predicted_all202': int(stats.pvalue.lt(.05).sum()),
                    'BH005_predicted_all202': int(stats.FDR.lt(.05).sum()), 'BH010_predicted_all202': int(stats.FDR.lt(.1).sum())})
            print(label, representation, 'complete 139 x 202 inference and three donor contrasts', flush=True)
    threshold_table = pd.DataFrame(all_thresholds)
    comparison_table = pd.concat(all_comparisons, ignore_index=True)
    threshold_table.to_csv(OUT / 'all_threshold_results.csv', index=False)
    comparison_table.to_csv(OUT / 'every_lipid_all_comparisons.csv', index=False)
    pd.DataFrame(reference_summaries).to_csv(OUT / 'reference_only_concordance.csv', index=False)
    pd.concat(reference_details, ignore_index=True).to_csv(OUT / 'reference_only_every_lipid.csv', index=False)
    pd.DataFrame(numerical).to_csv(OUT / 'overall_direction_and_effect_summary.csv', index=False)
    for label, values in audit['model_checks'].items():
        if sha(values['path']) != values['sha256']:
            raise ValueError(f'{label} bundle changed during inference')
    audit['inference_completed'] = True
    audit['script_sha256'] = sha(__file__)
    (OUT / 'comparison_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(threshold_table[(threshold_table.RNA_representation == 'log2_CPM') &
          (threshold_table.region == 'donor_balanced') & (threshold_table.minimum_predicted_fold == 1)]
          .drop(columns=['concordant_lipids', 'discordant_lipids']).to_string(index=False), flush=True)


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
