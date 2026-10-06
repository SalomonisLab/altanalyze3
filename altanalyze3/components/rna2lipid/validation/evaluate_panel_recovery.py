"""Diagnostic-only evaluation of complete-panel source and raw-MS recovery.

Never trains, changes a bundle, imputes corrected labels, or drops baseline rows.
Legacy target reconstruction is explicitly separate from corrected targets.
"""
from pathlib import Path
import argparse
import hashlib
import json
import pickle
import re
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import r2_score

HERE = Path(__file__).resolve().parents[1]
SOURCE = Path('/Users/saljh8/Downloads/Lipidomics')
DEFAULT_OUT = HERE / 'artifacts/full_panel_recovery_audit_20261005'
sys.path.insert(0, str(HERE))
from prepare_lungmap_targets import read_sheet
import massive_lipid_reanalysis as raw_tools
from candidate_integrity import require_diagnostic_authorization


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def base_name(name):
    return re.sub(r'_(?:POS|NEG|P|N)$', '', str(name))


def mode_name(name):
    return str(name).replace('_POS', '_P').replace('_NEG', '_N')


def run(out, scan_raw=False):
    require_diagnostic_authorization('full_panel_source_recovery')
    out.mkdir(parents=True, exist_ok=True)
    production = HERE / 'rna2lipid_hs_lung_lipidwise_bundle.pkl'
    before = digest(production)
    with production.open('rb') as handle:
        bundle = pickle.load(handle)
    lipids, profiles = list(bundle['Y_columns']), list(bundle['training_samples'])
    prep = SOURCE / 'Preprocessing Process'
    paths = {
        'production_bundle': production,
        'native_workbook': SOURCE / 'LMEX0000003692_data_table.xlsx',
        'unprocessed_cell': prep / 'Unprocessed_cell_type_clair_unnormalized.csv',
        'final_cell': prep / 'Final_cell_cell_lipids_cleaned_norm_median_286_log210.csv',
        'unprocessed_bulk': prep / 'unprocessed_bulk.csv',
        'final_bulk': prep / 'Final_bulk_Bulk_lipids_cleaned_normalized_median_527_log210.csv',
        'RNA': HERE / 'data/newrna_cell_clair_filtered_symbol.csv',
    }
    rows, formulas = read_sheet(paths['native_workbook'])
    sample_cols = [c for c, v in rows[3].items()
                   if isinstance(v, str) and re.fullmatch(r'Sample_\d+_\w+', v)][:45]
    lipid_rows = [r for i, r in rows.items() if i >= 4 and r.get('A')]
    native = pd.DataFrame([[r.get(c, np.nan) for r in lipid_rows] for c in sample_cols],
                          columns=[r['A'] for r in lipid_rows],
                          index=['D' + rows[3][c].split('_')[1].zfill(3) + '_' +
                                 rows[3][c].split('_')[2] for c in sample_cols], dtype=float)
    cached = pd.read_csv(HERE / 'validation/lungmap_scale_audit_20261001/targets_native_log2.csv', index_col=0)
    np.testing.assert_allclose(native.values, cached.values, atol=1e-12, rtol=0, equal_nan=True)
    cell = pd.read_csv(paths['unprocessed_cell'], index_col=0).T
    final = pd.read_csv(paths['final_cell'], index_col=0)
    bulk = pd.read_csv(paths['unprocessed_bulk'], index_col=0).T
    bulk.columns = [mode_name(c) for c in bulk.columns]
    final_bulk = pd.read_csv(paths['final_bulk'], index_col=0)
    broad = pd.read_csv(HERE / 'artifacts/reference_normalization_20261001/native_lipid_targets.csv', index_col=0)
    restricted = pd.read_csv(HERE / 'artifacts/MS1_reference_normalization_20261001/combined_reference_MS1_log2.csv', index_col=0)
    raw_dir = HERE / 'artifacts/MSV000081973_reanalysis_20261001'
    locked = pd.read_csv(raw_dir / 'locked_targets.tsv', sep='\t')
    qc = pd.read_csv(raw_dir / 'per_lipid_replication.csv')
    qc['lipid'] = qc.feature_id.str.split('|').str[0]
    qc['passes_prior_screen'] = (qc.isobaric_targets.isna() & (qc.ms2_support_scans > 0) &
        (qc.observed_profiles == 15) & (qc.rmse_log2 <= .75) & (qc.spearman >= .7))

    # Reconstruct the ORIGINAL target scale as a provenance test only.
    # Inversion zeros are missing, then sample medians use the entire 286-lipid
    # source table BEFORE selecting the production 202 outputs.
    recovered_linear = (np.exp2(final) - 1) / 10
    legacy_unfilled = np.log2(recovered_linear.where(recovered_linear > 0))
    legacy_filled = legacy_unfilled.T.fillna(legacy_unfilled.median(axis=1)).T
    legacy = legacy_filled.loc[profiles, lipids]
    legacy.to_csv(out / 'LEGACY_training_targets_reconstructed_NOT_corrected.csv')
    legacy_unfilled.isna().loc[profiles, lipids].to_csv(out / 'legacy_target_fill_mask.csv')
    x = pd.read_csv(paths['RNA'], index_col=0).T
    x = x.T.groupby(level=0).mean().T.loc[profiles, bundle['X_columns']]
    x = x.fillna(x.median())
    np.testing.assert_allclose(x.mean(), bundle['scaler_x'].mean_, atol=1e-12, rtol=0)
    np.testing.assert_allclose(legacy.mean(), bundle['scaler_y'].mean_, atol=1e-12, rtol=0)
    np.testing.assert_allclose(legacy.std(ddof=0), bundle['scaler_y'].scale_, atol=1e-12, rtol=0)
    xs = pd.DataFrame(bundle['scaler_x'].transform(x), index=x.index, columns=x.columns)
    ys = bundle['scaler_y'].transform(legacy)
    fit_summary = bundle['summary'].set_index('Lipid')
    fit_rows = []
    for j, lipid in enumerate(lipids):
        model = bundle['models'][lipid]
        pred = model['model'].predict(xs[model['genes']])
        r2 = r2_score(ys[:, j], pred)
        pearson = np.corrcoef(ys[:, j], pred)[0, 1] if np.std(pred) > 1e-14 else np.nan
        fit_rows.append({'lipid': lipid, 'replayed_train_R2': r2,
                         'saved_train_R2': fit_summary.loc[lipid, 'Train_R2_scaled'],
                         'R2_error': abs(r2 - fit_summary.loc[lipid, 'Train_R2_scaled']),
                         'replayed_train_Pearson': pearson,
                         'saved_train_Pearson': fit_summary.loc[lipid, 'Train_Pearson_scaled']})
    fit = pd.DataFrame(fit_rows)
    assert fit.R2_error.max() < 1e-12
    fit.to_csv(out / 'legacy_target_production_fit_replay.csv', index=False)

    crosswalk, source_matches = [], []
    for lipid in lipids:
        native_ids = [c for c in native if base_name(c) == lipid]
        assert native_ids, lipid
        for feature in native_ids:
            a = native[feature]
            b = cell.loc[native.index, lipid]
            mask = a.notna() & b.notna()
            source_matches.append({'lipid': lipid, 'native_feature': feature,
                'compared_values': int(mask.sum()), 'missing_native_values': int(a.isna().sum()),
                'max_error_vs_supplied_cell': float(abs(a[mask] - b[mask]).max())})
        merged = native[native_ids].fillna(0).max(axis=1)
        merge_error = float(abs(merged - cell.loc[native.index, lipid]).max())
        assert merge_error < 1e-8
        bulk_ids = [c for c in bulk if base_name(c) == lipid]
        shared_modes = [c for c in native_ids if c in bulk]
        in_broad = any(c in broad for c in native_ids)
        loss = ('retained_in_219' if in_broad else
                'removed_by_exact_mode_or_annotation_intersection' if not shared_modes else
                'removed_by_complete_case_filter')
        raw_rows = qc[qc.lipid == lipid]
        crosswalk.append({'production_lipid': lipid,
            'native_features': '; '.join(native_ids), 'native_mode_candidates': len(native_ids),
            'native_missing_entries_across_modes': int(native[native_ids].isna().sum().sum()),
            'native_profiles_with_no_mode_observed': int(native[native_ids].isna().all(axis=1).sum()),
            'supplied_cell_max_mode_after_zero_fill_error': merge_error,
            'bulk_exact_base_features': '; '.join(bulk_ids),
            'bulk_exact_mode_features': '; '.join(shared_modes),
            'in_final_cell': lipid in final, 'in_final_bulk': lipid in final_bulk,
            'intermediate_reference_status': loss,
            'in_restricted_47': any(c in restricted for c in native_ids),
            'prior_raw_exact_name_targets': '; '.join(raw_rows.feature_id),
            'prior_raw_pass_screen_any_mode': bool(raw_rows.passes_prior_screen.any()),
            'prior_raw_has_MS2_any_mode': bool((raw_rows.ms2_support_scans > 0).any()),
            'legacy_target_zero_fill_count': int(legacy_unfilled.loc[profiles, lipid].isna().sum()),
            'corrected_target_status': 'not_yet_reconciled_for_all_50_profiles'})
    crosswalk = pd.DataFrame(crosswalk)
    crosswalk.to_csv(out / 'all_202_lipid_reconciliation.csv', index=False)
    pd.DataFrame(source_matches).to_csv(out / 'all_native_mode_value_checks.csv', index=False)
    # Diagnose annotation suffixes numerically, without choosing a new merge
    # rule or treating near agreement as proof of biological equivalence.
    raw_bulk = pd.read_csv(paths['unprocessed_bulk'], index_col=0).T
    inverted_bulk = np.log2(((np.exp2(final_bulk)-1)/10).where(final_bulk > 0))
    bulk_alias_rows = []
    for lipid in lipids:
        source_ids = [c for c in raw_bulk
                      if re.sub(r'(?:_[A-D])?_(?:POS|NEG)$', '', c) == lipid]
        if not source_ids:
            continue
        samples = raw_bulk.index.intersection(inverted_bulk.index)
        z = raw_bulk.loc[samples, source_ids]
        for method, prediction in [('max', z.max(axis=1)), ('mean_log', z.mean(axis=1)), ('first', z.iloc[:, 0])]:
            error = abs(prediction - inverted_bulk.loc[samples, lipid])
            bulk_alias_rows.append({'lipid': lipid, 'source_ids': '; '.join(source_ids),
                'source_count': len(source_ids), 'method': method, 'max_abs_error': float(error.max()),
                'median_abs_error': float(error.median()), 'n': int(error.notna().sum())})
    bulk_alias = pd.DataFrame(bulk_alias_rows)
    bulk_alias.to_csv(out / 'bulk_source_alias_diagnostics.csv', index=False)
    metadata = []
    for sample in profiles:
        metadata.append({'sample': sample, 'donor': sample.split('_')[0],
                         'population': sample.split('_')[1], 'in_native_workbook': sample in native.index,
                         'in_unprocessed_cell': sample in cell.index, 'in_final_cell': sample in final.index,
                         'in_RNA': sample in x.index,
                         'legacy_nonzero_source_target_values': int(legacy_unfilled.loc[sample, lipids].notna().sum()),
                         'legacy_sample_median_filled_targets': int(legacy_unfilled.loc[sample, lipids].isna().sum()),
                         'legacy_target_values_after_documented_fill': int(legacy.loc[sample].notna().sum())})
    pd.DataFrame(metadata).to_csv(out / 'all_50_profile_reconciliation.csv', index=False)

    # Route B: retain every baseline identity in the audit. Only reuse adducts
    # recorded for the same lipid and ion mode; never guess missing adducts.
    rules_dir = HERE / 'validation/MSV000081973_reference/LIQUID_rules'
    cr = raw_tools.read_rules(rules_dir / 'DefaultCompositionRules.txt')
    fr = raw_tools.read_rules(rules_dir / 'DefaultFragmentationRules.txt')
    prior = raw_tools.targets_from_published(HERE / 'validation/MSV000081973_reference/published_lipid_profiles.csv', rules_dir)
    hypotheses, target_mapping, formula_rows = {}, [], []
    for lipid in lipids:
        formula = raw_tools.composition(lipid, cr)
        formula_text = ''.join(a + (str(n) if n != 1 else '') for a,n in formula.items() if n)
        formula_rows.append({'production_lipid': lipid, 'formula_from_annotation': formula_text,
                             'neutral_mass_from_annotation': raw_tools.mass(formula),
                             'interpretation': 'computed composition, not an observed identification'})
        native_ids = [c for c in native if base_name(c) == lipid]
        for feature in native_ids:
            positive = feature.endswith('_P')
            mode = 'positive' if positive else 'negative'
            existing = locked[(locked.lipid == lipid) & (locked.polarity == mode)]
            adducts = list(existing.adduct.unique())
            evidence = 'same_lipid_and_mode_in_published_raw_target_list'
            if not adducts:
                target_mapping.append({'production_lipid': lipid, 'native_feature': feature,
                                       'feature_id': None,
                                       'adduct_provenance': 'unresolved_source_annotation_required'})
            for adduct in adducts:
                shift = {'[M+H]+': raw_tools.ATOMS['H']-raw_tools.ELECTRON,
                         '[M-H]-': -raw_tools.ATOMS['H']+raw_tools.ELECTRON,
                         '[M+NH4]+': raw_tools.ATOMS['N']+4*raw_tools.ATOMS['H']-raw_tools.ELECTRON}[adduct]
                fid = lipid + '|' + adduct
                target = {'feature_id': fid, 'lipid': lipid, 'adduct': adduct, 'polarity': mode,
                          'formula': formula_text,
                          'mz': raw_tools.mass(formula)+shift, 'adduct_provenance': evidence}
                target['fragments'] = raw_tools.fragments(target, fr)
                hypotheses[fid] = target
                target_mapping.append({'production_lipid': lipid, 'native_feature': feature,
                                       'feature_id': fid, 'adduct_provenance': evidence})
    # Include competing existing annotations in the ambiguity assessment.
    competitors = {t['feature_id']: t for t in prior}
    competitors.update(hypotheses)
    for fid, t in hypotheses.items():
        t['isobaric_targets'] = ';'.join(k for k,c in competitors.items()
            if k != fid and c['polarity'] == t['polarity'] and abs(c['mz'] - t['mz']) <= t['mz'] * 10e-6)
    pd.DataFrame(target_mapping).to_csv(out / 'raw_target_hypotheses_all_202.csv', index=False)
    pd.DataFrame(formula_rows).to_csv(out / 'all_202_annotation_compositions.csv', index=False)
    pd.DataFrame([{k:v for k,v in t.items() if k != 'fragments'} for t in hypotheses.values()]).to_csv(
        out / 'raw_mass_hypotheses.tsv', sep='\t', index=False)
    files = sorted((SOURCE / 'MSV000081973/mzML').glob('*_L_*.mzML'))
    inventory = [{'file': str(p), 'bytes': p.stat().st_size,
                  'donor': 'D' + re.search(r'_D(\d+)_', p.name)[1].zfill(3),
                  'population': re.search(r'_D\d+_(\w+)_1_L_', p.name)[1]} for p in files]
    pd.DataFrame(inventory).to_csv(out / 'raw_file_inventory.csv', index=False)
    if scan_raw:
        sys.path.insert(0, '/Users/saljh8/Documents/GitHub/pyNeoQuant/src')
        from pyneoquant.quant.lipid_ms1 import iter_scans, mass_signal, quantify
        raw_out = out / 'expanded_raw_diagnostic'
        raw_out.mkdir(exist_ok=True)
        found, evidence, scans = raw_tools.discover(files, list(hypotheses.values()), iter_scans, mass_signal, ppm=5)
        pd.DataFrame(found).to_csv(raw_out / 'candidate_RT_windows.tsv', sep='\t', index=False)
        pd.DataFrame(evidence).to_csv(raw_out / 'MS2_candidate_evidence.csv', index=False)
        pd.DataFrame(scans).to_csv(raw_out / 'spectrum_inventory.csv', index=False)
        quant = []
        for p in files:
            mode = 'positive' if '_POS_' in p.name else 'negative'
            m = re.search(r'_D(\d+)_(\w+)_1_L_', p.name)
            sample = 'D' + m[1].zfill(3) + '_' + m[2]
            quant.extend(dict(row, sample=sample, file=p.name) for row in quantify(
                p, [t for t in found if t['polarity'] == mode], 5, 45, 3))
            print('diagnostic quantified', sample, mode, flush=True)
        pd.DataFrame(quant).to_csv(raw_out / 'candidate_peak_quantification.csv', index=False)
    summary = {
        'scope': 'Both recovery routes evaluated; no training or corrected target imputation authorized or performed',
        'production_lipids': len(lipids), 'production_profiles': len(profiles),
        'production_RNA_genes': len(bundle['X_columns']),
        'native_workbook_lipid_name_matches': int(crosswalk.native_features.ne('').sum()),
        'native_ion_mode_records_for_production_names': len(source_matches),
        'native_profiles': len(native),
        'native_sample_lipid_entries_missing_in_all_modes': int(crosswalk.native_profiles_with_no_mode_observed.sum()),
        'intermediate_reference_accounting': crosswalk.intermediate_reference_status.value_counts().to_dict(),
        'supplied_cell_reconstruction': 'native missing cells -> 0, then maximum across same-name ion modes',
        'supplied_cell_reconstruction_max_error': float(crosswalk.supplied_cell_max_mode_after_zero_fill_error.max()),
        'legacy_target_reconstruction': 'log2((2**final-1)/10); zeros become missing; fill each sample median across all 286 source lipids; select original 50x202',
        'legacy_target_mean_max_error': float(abs(legacy.mean().values-bundle['scaler_y'].mean_).max()),
        'legacy_target_sd_max_error': float(abs(legacy.std(ddof=0).values-bundle['scaler_y'].scale_).max()),
        'legacy_target_saved_R2_max_error': float(fit.R2_error.max()),
        'D071_nonzero_source_legacy_cells': int(legacy_unfilled.loc[[s for s in profiles if s.startswith('D071_')], lipids].notna().sum().sum()),
        'D071_filled_legacy_cells': int(legacy_unfilled.loc[[s for s in profiles if s.startswith('D071_')], lipids].isna().sum().sum()),
        'bulk_alias_candidate_name_coverage': int(bulk_alias.lipid.nunique()),
        'bulk_mean_log_rule_lipids_with_max_error_le_0101': int(((bulk_alias.method == 'mean_log') & (bulk_alias.max_abs_error <= .101)).sum()),
        'bulk_alias_unresolved_names': sorted(set(lipids) - set(bulk_alias.lipid)),
        'raw_prior_exact_name_lipids': int(crosswalk.prior_raw_exact_name_targets.ne('').sum()),
        'raw_prior_QC_pass_lipids_any_mode': int(crosswalk.prior_raw_pass_screen_any_mode.sum()),
        'raw_files': len(files), 'raw_donors': sorted({v['donor'] for v in inventory}),
        'raw_existing_exact_name_and_mode_lipids': len({t['lipid'] for t in hypotheses.values()}),
        'raw_existing_exact_name_and_mode_targets': len(hypotheses),
        'raw_computable_neutral_compositions': len(formula_rows),
        'raw_unresolved_native_mode_records': sum(row['feature_id'] is None for row in target_mapping),
        'raw_new_scan_completed': scan_raw,
        'remaining': ['D071 upstream original-scale provenance not yet reconciled; request original table and preprocessing script from user',
                      'Production target values and original relative workbook are different representations; no corrected 50x202 matrix established',
                      'Four production names have two native ion-mode records; historical max/zero rule must not be silently substituted for a corrected policy',
                      'Missing native values require explicit handling discussion, not feature removal',
                      'Raw target evidence must annotate, not determine, model-panel eligibility',
                      'Request identification/export table with lipid IDs, ion modes, adducts and retention times for unresolved raw targets; do not guess adducts'],
        'input_sha256': {k: {'path': str(p), 'sha256': digest(p)} for k,p in paths.items()},
        'production_unchanged': digest(production) == before,
    }
    assert summary['production_unchanged']
    (out / 'recovery_audit.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k:v for k,v in summary.items() if k != 'input_sha256'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--scan-raw', action='store_true')
    args = parser.parse_args()
    run(args.out, args.scan_raw)
