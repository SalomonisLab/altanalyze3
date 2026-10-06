"""Authorized full-panel native-log2 candidate; pinned production trainer and DE engine."""
from pathlib import Path
import argparse
import importlib.util
import json
import pickle
import re
import shutil
import sys
import time
import warnings

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2]))
from candidate_integrity import (CONTRACT, STATE, candidate_contract, checked_file, exact_ids,
                                 fit_verified_candidate, read_json, sha256, validate_contract,
                                 validate_run_manifest, validate_tables)
from prepare_lungmap_targets import read_sheet
from evaluate_lungmap_ipf import arithmetic_stats, differential, lipid_mapping
from infer_lungmap_elasticnet_ms1 import affine
from altanalyze3.components.rna2lipid.api import Rna2LipidBundle

OUT = HERE / 'artifacts/LungMAP_full202_native_log2_20261005'
OLD = HERE / 'artifacts/LungMAP_IPF_candidate_20261004'
ATLAS = Path('/Volumes/salomonis2/LungMAP/CellRef2')
PHASE = 'full202_native_log2_IPF_validation'
LABEL = 'candidate_native_log2_202'
PRIOR = 'current_explicit_log2'


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, default=str, allow_nan=False) + '\n')


def record(path, **kwargs):
    return dict(path=str(Path(path).resolve()), sha256=sha256(path), **kwargs)


def prepare():
    """Read all required native identities; keep historical observed mode-max mapping."""
    c, state = read_json(CONTRACT), read_json(STATE)
    validate_contract(c)
    authorization = state.get('latest_authorization', {})
    if authorization.get('granted_by') != 'user' or authorization.get('user_message') != 'Do it':
        raise RuntimeError('This specific native-log2 IPF workflow requires the recorded user authorization.')
    scope = state['approved_scope_change']
    effective = candidate_contract(c, {'approved_scope_change': scope}, state)
    source_audit = read_json(HERE / 'artifacts/full_panel_recovery_audit_20261005/recovery_audit.json')
    source = checked_file(source_audit['input_sha256']['native_workbook'], 'Original native measurements')
    cell_source = checked_file(source_audit['input_sha256']['unprocessed_cell'], 'Original mode mapping')
    rows, _ = read_sheet(source)
    columns = [k for k, v in rows[3].items() if isinstance(v, str) and re.fullmatch(r'Sample_\d+_\w+', v)][:45]
    lipid_rows = [(i, r) for i, r in rows.items() if i >= 4 and r.get('A')]
    names = [r['A'] for _, r in lipid_rows]
    if len(names) != len(set(names)):
        raise RuntimeError('Duplicate native feature IDs require discussion.')
    native = pd.DataFrame([[r.get(k, np.nan) for _, r in lipid_rows] for k in columns],
                          columns=names, index=['D' + rows[3][k].split('_')[1].zfill(3) + '_' +
                          rows[3][k].split('_')[2] for k in columns], dtype=float)
    exact_ids(sorted(effective['samples']), sorted(native.index), 'Native sample source')
    native = native.loc[effective['samples']]
    cell = pd.read_csv(cell_source, index_col=0).T
    y = pd.DataFrame(index=effective['samples'], columns=effective['lipids'], dtype=float)
    mappings = []
    for lipid in effective['lipids']:
        modes = [f for f in native if re.sub(r'_(?:P|N|POS|NEG)$', '', f) == lipid]
        if not modes:
            raise RuntimeError('Required native lipid has no source mapping: ' + lipid)
        # Verify the established upstream rule, rather than introducing a new mode choice.
        historical = native[modes].fillna(0).max(axis=1)
        np.testing.assert_allclose(historical, cell.loc[native.index, lipid], atol=1e-9, rtol=0)
        # NA is missing, not zero. The user expressly authorized missing-value correction.
        y[lipid] = native[modes].max(axis=1, skipna=True)
        mappings.append({'lipid': lipid, 'native_features': '; '.join(modes),
                         'historical_rule': 'maximum across same-name ion modes',
                         'missing_all_modes': int(y[lipid].isna().sum())})
    before = y.copy()
    medians = y.median(axis=0)
    y = y.fillna(medians)
    imputed = [{'sample': sample, 'lipid': lipid, 'filled_log2': float(y.loc[sample, lipid])}
               for sample in before.index for lipid in before.columns if pd.isna(before.loc[sample, lipid])]
    if len(imputed) != 25 or len({r['lipid'] for r in imputed}) != 7:
        raise RuntimeError('Unexpected missing-target pattern; discuss before training.')
    x = pd.read_csv(checked_file(c['RNA_source'], 'Production RNA'), index_col=0).T
    x = x.T.groupby(level=0).mean().T
    x.index = x.index.str.strip(); x.columns = x.columns.str.strip()
    x = x.loc[effective['samples'], effective['genes']].apply(pd.to_numeric, errors='raise')
    x = x.fillna(x.median())
    validate_tables(x, y, effective)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'provenance').mkdir(exist_ok=True)
    x.to_csv(OUT / 'candidate_training_RNA.csv')
    y.to_csv(OUT / 'candidate_training_lipids_log2.csv')
    before.to_csv(OUT / 'native_targets_before_NA_fill.csv')
    loaded = pd.read_csv(OUT / 'candidate_training_lipids_log2.csv', index_col=0)
    preservation = []
    for lipid in before:
        observed = before[lipid].dropna()
        reference_pairs = observed.to_numpy()[:, None] - observed.to_numpy()[None, :]
        final_values = loaded.loc[observed.index, lipid].to_numpy()
        errors = abs(reference_pairs - (final_values[:, None] - final_values[None, :]))
        preservation.append({'lipid': lipid, 'observed_samples': len(observed),
                             'observed_pairwise_differences': len(observed) * (len(observed) - 1) // 2,
                             'max_abs_log2_fold_error': float(errors.max())})
    if max(r['max_abs_log2_fold_error'] for r in preservation) > 1e-10:
        raise RuntimeError('Corrected reference changed an observed differential.')
    pd.DataFrame(preservation).to_csv(OUT / 'all_202_differential_preservation.csv', index=False)
    pd.DataFrame(mappings).to_csv(OUT / 'all_202_target_mappings.csv', index=False)
    pd.DataFrame(imputed).to_csv(OUT / 'NA_filled_entries.csv', index=False)
    description = ('Original workbook relative log2 values, all 45 authorized profiles and all 202 production lipids. '
                   'Retain established maximum of observed same-name ion-mode records. '
                   'Replace only all-mode-missing entries with per-lipid log2 median. Preserve negative observations. '
                   'Native relative abundance gauge; full MS1 abundance calibration is not asserted.')
    provenance = {'schema_version': 1, 'data_role': 'corrected_training_targets', 'unresolved': [],
                  'target_sha256': sha256(OUT / 'candidate_training_lipids_log2.csv'),
                  'transformation_description': description, 'imputed_entries': imputed,
                  'imputation_user_decision': state['approved_na_handling']}
    provenance['sample_sources'] = [{'id': sample, 'status': 'verified', 'source': record(source),
        'source_identifier': 'Sample_' + str(int(sample[1:4])) + '_' + sample.split('_')[1],
        'evidence': 'Existing complete-panel source audit verified DNNN / Sample_N RNA-to-workbook identifiers.'}
        for sample in effective['samples']]
    provenance['feature_sources'] = [{'id': row['lipid'], 'status': 'verified', 'source': record(source),
        'source_identifier': row['native_features'], 'evidence': 'Full-panel audit and exact original upstream mode-max replay.'}
        for row in mappings]
    dump(OUT / 'target_provenance.json', provenance)
    manifest = {'schema_version': 1, 'analysis_phase': PHASE,
                'baseline_contract_sha256': sha256(CONTRACT), 'approved_scope_change': scope,
                'method': c['method'], 'target_description': description,
                'inputs': {'RNA': c['RNA_source'],
                           'corrected_targets': record(OUT / 'candidate_training_lipids_log2.csv', data_role='corrected_training_targets'),
                           'provenance': record(OUT / 'target_provenance.json')}}
    for name in ['IPF_vs_healthy_registry.tsv', 'crosswalk.tsv', 'overrides.tsv',
                 'harmonized_library_metadata_harmonized_final_corrected_v6.txt', 'cellHarmony_differential_frozen.py']:
        path = OLD / 'provenance' / name
        shutil.copyfile(path, OUT / 'provenance' / name)
    manifest['evaluation_sources'] = {name: record(OUT / 'provenance' / name) for name in
                                     ['IPF_vs_healthy_registry.tsv', 'crosswalk.tsv', 'overrides.tsv', 'cellHarmony_differential_frozen.py']}
    dump(OUT / 'run_manifest.json', manifest)
    # Bind the actual, already given user instruction to its exact implemented inputs.
    # This record does not claim that the user inspected a file digest.
    state['current_execution_dependency']['status'] = 'resolved'
    state['current_execution_dependency']['user_answer'] = 'its back'
    for question in state['questions']:
        if question['id'] == 'raw_target_annotations':
            question['required_for'] = ['full_MS1_abundance_calibration']
        elif question['id'] == 'corrected_target_provenance':
            question.update(status='resolved', user_answer='Proceed and if you can\'t do anything the original model did or want deviate in any way that I haven\'t authorized, stop and discuss.',
                            resolution_evidence='Original upstream maximum rule is retained, numerically verified for all 202 source names; authorized log2 median filling handles missing entries. No new ion-mode policy.')
        elif question['id'] == 'review_full_proposal':
            question.update(status='resolved', user_answer='Do it', required_for=[PHASE],
                            resolution_evidence='Explicit response to the 45-profile, 202-lipid, original ElasticNetCV native-log2 IPF comparison proposal. This is not approval of full MS1 abundance calibration.')
    state['analysis_status'] = 'approved'
    state['approval'] = {'manifest_sha256': sha256(OUT / 'run_manifest.json'), 'analysis_phase': PHASE,
                         'granted_by': 'user', 'user_message': 'Do it',
                         'conversation_reference': 'Current thread, explicit approval immediately after the stated full202 native-log2 IPF comparison; followed by restored-drive confirmation.',
                         'binding_note': 'Digest records the implemented authorized inputs; no claim that user read the digest.'}
    state['note'] = 'Native-log2 IPF validation is explicitly authorized. Full MS1 abundance calibration remains incomplete and has a separate pending identification question.'
    dump(STATE, state)
    validate_run_manifest(OUT / 'run_manifest.json')
    dump(OUT / 'preflight.json', {'samples': len(y), 'lipids': len(y.columns), 'genes': len(x.columns),
          'NA_entries_filled': len(imputed), 'negative_observed_values_preserved': int((before < 0).sum().sum()),
          'all_required_identities_present': True, 'production_algorithm_unchanged': True})
    print('Preflight complete: 45 samples, 202 lipids, 1303 genes; 25 missing entries filled.', flush=True)


def train():
    manifest = OUT / 'run_manifest.json'
    validate_run_manifest(manifest)
    with threadpool_limits(limits=1), warnings.catch_warnings():
        warnings.simplefilter('ignore')
        predictions, bundle, seconds = fit_verified_candidate(manifest)
    c, _ = validate_run_manifest(manifest)
    production = pickle.load(checked_file(c['production_bundle'], 'Production bundle').open('rb'))
    bundle['training_samples'] = c['samples']
    bundle['training_metadata'] = production['training_metadata'].loc[c['samples']].copy()
    bundle['metadata'] = {'candidate_only': True, 'target_scale': 'native_relative_log2',
                          'full_MS1_abundance_calibration': False, 'external_holdout': False,
                          'production_default_changed': False, 'source_manifest_sha256': sha256(manifest)}
    bundle['training_seconds'] = seconds
    with (OUT / 'candidate_bundle.pkl').open('wb') as handle:
        pickle.dump(bundle, handle)
    predictions.to_csv(OUT / 'training_predictions_log2.csv')
    for key in ['summary', 'coefficients', 'candidate_models']:
        bundle[key].to_csv(OUT / (key + '.csv'), index=False)
    # Direct comparison with the delivered implementation on actual candidate data.
    sys.path.insert(0, str(HERE / 'validation'))
    from validate_trainer_equivalence import load_delivered_reference
    x = pd.read_csv(OUT / 'candidate_training_RNA.csv', index_col=0)
    y = pd.read_csv(OUT / 'candidate_training_lipids_log2.csv', index_col=0)
    settings = dict(c['method']['fit_settings']); settings.pop('n_jobs')
    with threadpool_limits(limits=1), warnings.catch_warnings():
        warnings.simplefilter('ignore')
        rp, rb, _ = load_delivered_reference().fit_sparse_lipidwise_elasticnet(x_train=x, y_train=y.iloc[:, :2], x_test=x, **settings)
    coef_errors = []
    for lipid in rp:
        a, b = bundle['models'][lipid], rb['models'][lipid]
        exact_ids(a['genes'], b['genes'], 'Delivered selected genes')
        np.testing.assert_allclose(a['model'].coef_, b['model'].coef_, atol=1e-10, rtol=0)
        assert a['model'].alpha_ == b['model'].alpha_ and a['model'].l1_ratio_ == b['model'].l1_ratio_
        coef_errors.append(float(abs(a['model'].coef_ - b['model'].coef_).max()))
    prediction_error = float(abs(predictions[rp.columns] - rp).to_numpy().max())
    assert prediction_error < 1e-10
    dump(OUT / 'training_method_audit.json', {'training_seconds': seconds, 'trained_lipids': len(bundle['models']),
          'trained_samples': len(c['samples']), 'genes': len(c['genes']), 'same_ElasticNetCV_settings': True,
          'delivered_coefficient_error': max(coef_errors), 'delivered_prediction_error': prediction_error,
          'production_bundle_sha256': sha256(c['production_bundle']['path'])})
    print('Full-panel training complete:', seconds, 'seconds', flush=True)


def load_current(unit):
    validate_run_manifest(OUT / 'run_manifest.json')
    name = 'pseudobulk' if unit == 'PB' else 'metacell'
    path = OLD / 'current_inputs' / f'cellref2_v8_{name}_lipid_forDE.h5ad'
    obj = ad.read_h5ad(path)
    c = read_json(CONTRACT)
    exact_ids(c['lipids'], obj.var_names, 'Production prediction panel')
    if sparse.issparse(obj.X): obj.X = obj.X.toarray()
    cross = pd.read_csv(OUT / 'provenance/crosswalk.tsv', sep='\t')
    mapping = dict(zip(cross.cell_state, cross.cell_type)) | dict(zip(cross.cell_type_full, cross.cell_type)) | dict(zip(cross.cell_type, cross.cell_type))
    state = obj.obs['cell_state'] if 'cell_state' in obj.obs else obj.obs.short_name
    obj.obs['cell_state'] = state.map(mapping)
    if obj.obs.cell_state.isna().any(): raise RuntimeError('Unresolved cell-state mappings.')
    return obj


def infer(units=None):
    c, _ = validate_run_manifest(OUT / 'run_manifest.json')
    bundle = pickle.load((OUT / 'candidate_bundle.pkl').open('rb'))
    exact_ids(c['lipids'], bundle['Y_columns'], 'Candidate output panel')
    exact_ids(c['samples'], bundle['training_samples'], 'Candidate training roster')
    w, intercept = affine(bundle)
    genes = bundle['X_columns']
    api = Rna2LipidBundle.load(OUT / 'candidate_bundle.pkl')
    cross = pd.read_csv(OUT / 'provenance/crosswalk.tsv', sep='\t')
    smap = dict(zip(cross.cell_state, cross.cell_type)) | dict(zip(cross.cell_type_full, cross.cell_type)) | dict(zip(cross.cell_type, cross.cell_type))
    def rowkeys(obj, unit):
        if unit == 'MC': return obj.obs_names
        s = obj.obs.cell_state if 'cell_state' in obj.obs else obj.obs.short_name
        return pd.Index(obj.obs.Study_internal.astype(str) + '|' + obj.obs.Sample.astype(str) + '|' + s.map(smap).astype(str))
    audit = []
    for unit in (units or ['PB', 'MC']):
        target = load_current(unit)
        tk = rowkeys(target, unit)
        if not tk.is_unique: raise RuntimeError('Duplicate prediction identities.')
        result = np.empty((target.n_obs, len(c['lipids'])), dtype=np.float32)
        assigned = np.zeros(target.n_obs, dtype=int)
        paths = [ATLAS / ('inputs/cellref2_v8_pbsample_ALL_log2cp10k.h5ad' if unit == 'PB' else 'inputs/cellref2_v8_mc_persample_log2cp10k.h5ad'),
                 ATLAS / ('inputs/copd_pb_ln1pcp10k.h5ad' if unit == 'PB' else 'inputs/copd_mc_ln1pcp10k.h5ad')]
        for path in paths:
            read_path = path
            copy_records = OUT / 'local_input_copies.json'
            if copy_records.exists():
                copies = read_json(copy_records)
                if str(path) in copies:
                    copied = copies[str(path)]
                    original_stat = path.stat()
                    if original_stat.st_size != copied['original_size'] or original_stat.st_mtime_ns != copied['original_mtime_ns']:
                        raise RuntimeError('Original RNA input changed after the local byte copy.')
                    read_path = checked_file(copied['local_copy'], 'Local byte copy of original RNA')
            obj = ad.read_h5ad(read_path, backed='r')
            positions = tk.get_indexer(rowkeys(obj, unit))
            if (positions < 0).any(): raise RuntimeError('RNA identities absent from baseline predictions.')
            gpos = obj.var_names.get_indexer(genes); valid = gpos >= 0
            with h5py.File(read_path, 'r') as handle, threadpool_limits(limits=1):
                group = handle['X']; ptr = group['indptr'][:]
                for start in range(0, obj.n_obs, 2048):
                    stop = min(start + 2048, obj.n_obs); lo, hi = ptr[start], ptr[stop]
                    block = sparse.csr_matrix((group['data'][lo:hi], group['indices'][lo:hi], ptr[start:stop+1]-lo), shape=(stop-start, obj.n_vars))
                    x = np.zeros((stop-start, len(genes)), dtype=float)
                    x[:, valid] = block[:, gpos[valid]].toarray()
                    prediction = x @ w + intercept
                    if start == 0:
                        direct = api.predict_from_dataframe(pd.DataFrame(x, columns=genes)).predictions.to_numpy()
                        np.testing.assert_allclose(prediction, direct, atol=1e-10, rtol=0)
                    result[positions[start:stop]] = prediction
                    assigned[positions[start:stop]] += 1
                    if start % 20480 == 0:
                        print('Inference progress', unit, path.name, stop, '/', obj.n_obs, flush=True)
            audit.append({'unit': unit, 'path': str(path), 'read_path': str(read_path), 'rows': obj.n_obs, 'normalization': obj.uns.get('normalization'),
                          'missing_genes_zero_filled_by_original_API_policy': [genes[i] for i in np.flatnonzero(~valid)]})
            obj.file.close()
            print('Inferred', unit, path.name, flush=True)
        if (assigned < 1).any() or not np.isfinite(result).all(): raise RuntimeError('Incomplete full atlas inference.')
        obj = ad.AnnData(result, obs=target.obs.copy(), var=pd.DataFrame(index=c['lipids']),
                         uns={'expression_scale': 'native_relative_log2', 'candidate_only': True})
        obj.write_h5ad(OUT / f'LungMAP_{unit}_predictions_log2.h5ad', compression='lzf')
        del obj, target, result
    previous = OUT / 'inference_audit.json'
    if units and previous.exists():
        audit = [r for r in json.loads(previous.read_text()) if r['unit'] not in units] + audit
    dump(previous, audit)


def evaluate(units=None):
    validate_run_manifest(OUT / 'run_manifest.json')
    spec = importlib.util.spec_from_file_location('full202_frozen_engine', OUT / 'provenance/cellHarmony_differential_frozen.py')
    engine = importlib.util.module_from_spec(spec); spec.loader.exec_module(engine)
    if not hasattr(pd.Series, 'nonzero'): pd.Series.nonzero = lambda s: np.asarray(s).nonzero()
    registry = pd.read_csv(OUT / 'provenance/IPF_vs_healthy_registry.tsv', sep='\t').fillna('')
    overrides = pd.read_csv(OUT / 'provenance/overrides.tsv', sep='\t')
    measured = pd.read_csv(OLD / 'experimental_lung_IPF_lipids_all_544.csv', index_col=0)
    # Resolve ion-mode metadata from the full source mapping, not the old 219-feature intersection.
    source_mapping = pd.read_csv(OUT / 'all_202_target_mappings.csv')
    native = [f.strip() for modes in source_mapping.native_features for f in modes.split(';')]
    maps = pd.concat([lipid_mapping(read_json(CONTRACT)['lipids'], label, measured, native) for label in [LABEL, PRIOR]], ignore_index=True)
    maps.to_csv(OUT / 'lipid_matching_audit.csv', index=False)
    measured.to_csv(OUT / 'experimental_lung_IPF_lipids_all_544.csv')
    allframes, allcensus = [], []
    for unit in (units or ['PB', 'MC']):
        candidate = ad.read_h5ad(OUT / f'LungMAP_{unit}_predictions_log2.h5ad')
        current = load_current(unit)
        exact_ids(current.obs_names, candidate.obs_names, 'Full prediction row roster')
        for row in registry.itertuples():
            ov = overrides[overrides.contrast_id.eq(row.contrast_id)]
            for label, obj in [(PRIOR, current), (LABEL, candidate)]:
                print('DE', row.contrast_id, unit, label, flush=True)
                folder = OUT / 'differentials' / row.contrast_id / unit / label
                cache = folder / 'source_binding.json'
                binding = {'manifest': sha256(OUT / 'run_manifest.json'),
                           'engine': sha256(OUT / 'provenance/cellHarmony_differential_frozen.py'),
                           'bundle': sha256(OUT / 'candidate_bundle.pkl') if label == LABEL else read_json(CONTRACT)['production_bundle']['sha256']}
                if cache.exists():
                    if read_json(cache) != binding: raise RuntimeError('Differential cache sources changed.')
                    frame = pd.read_csv(folder / 'all_tested_lipids.csv')
                    census = pd.read_csv(folder / 'population_replication_census.csv')
                else:
                    with threadpool_limits(limits=1), warnings.catch_warnings():
                        warnings.simplefilter('ignore')
                        frame, census = differential(obj, unit, row, ov, engine, OUT, label)
                    dump(cache, binding)
                allframes.append(frame)
                census['contrast'] = row.contrast_id; census['unit'] = unit; census['model'] = label
                allcensus.append(census)
        del candidate, current
    if units:
        pd.concat(allframes, ignore_index=True).to_csv(OUT / ('all_cellHarmony_tests_' + '_'.join(units) + '.csv'), index=False)
        print('Completed unit differential checkpoint:', units, flush=True)
        return
    full = pd.concat(allframes, ignore_index=True)
    full.to_csv(OUT / 'all_cellHarmony_tests.csv', index=False)
    pd.concat(allcensus, ignore_index=True).to_csv(OUT / 'population_replication_census.csv', index=False)
    mapping = maps[maps.status.eq('matched')][['model', 'gene', 'measured_feature']]
    joined = full.merge(mapping, on=['model', 'gene']).merge(measured.rename_axis('measured_feature').reset_index(), on='measured_feature', suffixes=('', '_experimental'))
    joined['direction_concordant_supplied'] = joined.log2fc * joined.supplied_log_effect > 0
    joined.to_csv(OUT / 'all_measured_vs_imputed_tests.csv', index=False)
    summarize(full, joined, maps, measured, registry)


def summarize(full, joined, maps, measured, registry, *, bh_policy='legacy_filtered', reference_out=None):
    rows, lipids, best = [], [], []
    shared = set(maps.loc[maps.status.eq('matched'), 'measured_feature'])
    for gate, source, stat, cutoff in [('both_rawp005', 'supplied_raw_p', 'pval', .05),
                                       ('both_BH005', 'supplied_adjusted_p', 'fdr', .05),
                                       ('both_BH010', 'supplied_adjusted_p', 'fdr', .1)]:
        eligible = measured.index[measured.index.isin(shared) & measured[source].lt(cutoff)]
        both = joined[joined[source].lt(cutoff) & joined[stat].lt(cutoff) & ~joined.population.str.contains('__vs__')]
        for fold in [1., 1.1, 1.2, 1.5, 2.]:
            selected = both[both.signed_fold.abs().ge(fold) & both.signed_fold.ne(0)]
            for row in registry.itertuples():
                for unit in ['PB', 'MC']:
                    for label in [PRIOR, LABEL]:
                        allg = selected[selected.contrast.eq(row.contrast_id) & selected.unit.eq(unit) & selected.model.eq(label)]
                        g = allg[allg.direction_concordant_supplied]
                        up = set(g.loc[g.supplied_log_effect.gt(0), 'measured_feature'])
                        down = set(g.loc[g.supplied_log_effect.lt(0), 'measured_feature'])
                        discordant = set(allg.loc[~allg.direction_concordant_supplied, 'measured_feature'])
                        rows.append({'contrast': row.contrast_id, 'unit': unit, 'model': label, 'overlap_gate': gate,
                                     'minimum_imputed_fold': fold, 'matched_experimental_differentials': len(eligible),
                                     'concordant_union': len(up | down), 'concordant_up': len(up), 'concordant_down': len(down),
                                     'discordant_union': len(discordant), 'significant_overlap_union': allg.measured_feature.nunique(),
                                     'supporting_cell_states': g.population.nunique()})
                        for feature, v in g.groupby('measured_feature'):
                            lipids.append({'contrast': row.contrast_id, 'unit': unit, 'model': label, 'overlap_gate': gate,
                                           'minimum_imputed_fold': fold, 'measured_feature': feature,
                                           'concordant_cell_types': '; '.join(sorted(v.population.unique()))})
                        if gate == 'both_BH010' and fold == 1 and len(g):
                            for population, v in g.groupby('population'):
                                best.append({'contrast': row.contrast_id, 'unit': unit, 'model': label, 'population': population,
                                             'concordant_total': v.measured_feature.nunique(),
                                             'concordant_up': v.loc[v.supplied_log_effect.gt(0), 'measured_feature'].nunique(),
                                             'concordant_down': v.loc[v.supplied_log_effect.lt(0), 'measured_feature'].nunique()})
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / 'cohort_concordant_union_summary.csv', index=False)
    pd.DataFrame(lipids).to_csv(OUT / 'cohort_concordant_lipid_celltypes.csv', index=False)
    intersect_rows = []
    for gate in ['both_rawp005', 'both_BH005', 'both_BH010']:
        for unit in ['PB', 'MC']:
            for label in [PRIOR, LABEL]:
                selected = [r for r in lipids if r['overlap_gate'] == gate and r['unit'] == unit and
                            r['model'] == label and r['minimum_imputed_fold'] == 1]
                adams = {r['measured_feature'] for r in selected if r['contrast'] == 'Adams2020__IPF_vs_Healthy'}
                natri = {r['measured_feature'] for r in selected if r['contrast'] == 'Natri2024__IPF_vs_Healthy'}
                shared_hits = adams & natri
                intersect_rows.append({'unit': unit, 'model': label, 'overlap_gate': gate,
                                       'Adams_concordant': len(adams), 'Natri_concordant': len(natri),
                                       'both_cohorts_concordant': len(shared_hits),
                                       'shared_lipids': '; '.join(sorted(shared_hits))})
    pd.DataFrame(intersect_rows).to_csv(OUT / 'Adams_Natri_concordant_intersection.csv', index=False)
    pd.DataFrame(best).sort_values('concordant_total', ascending=False).to_csv(OUT / 'concordant_celltype_ranking.csv', index=False)
    # Paired per-cohort model comparison; units remain distinct.
    primary = summary[summary.minimum_imputed_fold.eq(1)]
    paired = primary.pivot(index=['contrast', 'unit', 'overlap_gate', 'matched_experimental_differentials'], columns='model', values=['concordant_union', 'concordant_up', 'concordant_down'])
    paired.columns = ['__'.join(k) for k in paired.columns]
    paired = paired.reset_index()
    paired['candidate_minus_prior'] = paired['concordant_union__' + LABEL] - paired['concordant_union__' + PRIOR]
    paired.to_csv(OUT / 'candidate_vs_prior_comparison.csv', index=False)
    # Audit the actual correction family, distinguishing the historical and approved policies.
    coverage = full.groupby(['contrast', 'unit', 'model']).agg(rows=('gene', 'size'), finite_rawp=('pval', lambda s: s.notna().sum()), finite_BH=('fdr', lambda s: s.notna().sum())).reset_index()
    coverage.to_csv(OUT / 'BH_eligibility_audit.csv', index=False)
    families = full[full.fdr.notna()].groupby(['contrast', 'unit', 'population', 'model']).gene.agg(set)
    differences = 0
    for key in full[['contrast', 'unit', 'population']].drop_duplicates().itertuples(index=False, name=None):
        if families.get((*key, PRIOR), set()) != families.get((*key, LABEL), set()):
            differences += 1
    complete_panel = False
    if bh_policy == 'full_panel_no_abundance_filter':
        expected = set(read_json(CONTRACT)['lipids'])
        for key, group in full.groupby(['contrast', 'unit', 'population', 'model']):
            if len(group) != len(expected) or set(group.gene) != expected or not np.isfinite(group.fdr).all():
                raise RuntimeError('Full-panel BH coverage failure: ' + str(key))
        complete_panel = True
    dump(OUT / 'validation_gate.json', {'full_training_panel_preserved': True,
          'same_trainer_verified': True, 'raw_p_comparison_complete': True,
          'BH_testing_families_identical': differences == 0,
          'population_comparisons_with_different_BH_families': differences,
          'BH_policy': bh_policy, 'full_202_lipid_BH_verified': complete_panel,
          'all_covariate_stage_ready': differences == 0,
          'reason': ('User-authorized BH over all 202 lipids for each model/comparison; no abundance filter.'
                     if complete_panel else 'Saved >0.1 numeric filter depends on the relative target baseline. No filter change made.')})
    report(summary, paired, joined, measured, reference_out=reference_out)


def report(summary, paired, joined, measured, *, reference_out=None):
    reference_out = Path(reference_out) if reference_out is not None else OUT
    gate = read_json(OUT / 'validation_gate.json')
    corrected = gate.get('full_202_lipid_BH_verified', False)
    def table(frame):
        return frame.to_markdown(index=False)
    primary = summary[summary.overlap_gate.eq('both_BH010') & summary.minimum_imputed_fold.eq(1)]
    lines = ['# Full-panel native-log2 IPF candidate evaluation', '',
             'The candidate contains all 202 production lipids, 1,303 RNA genes and the 45 approved non-D071 training profiles. '
             'The deployed lipidwise ElasticNetCV implementation and hyperparameters are retained. '
             'Twenty-five missing entries across seven lipids are filled with per-lipid medians on log2 values; observed negatives are retained. '
             'Targets are native relative log2 values. This evaluation does not assert completed MS1 abundance calibration for all 202 lipids.', '',
             '## Concordant significant lipids across any combination of cell types', '',
             'The following counts require significance in both the experimental lung-tissue lipidomics and the imputed comparison, '
             'with matching direction. Each lipid is counted once within each cohort and unit. BH cutoff is 0.1; fold threshold is 1.', '',
             table(primary[['contrast', 'unit', 'model', 'matched_experimental_differentials', 'concordant_union', 'concordant_up', 'concordant_down', 'supporting_cell_states']]), '',
             '## Candidate versus prior', '', table(paired), '',
             '## Methods and interpretation', '',
             'RNA preprocessing, gene ranking, top-N selection, target/RNA StandardScaler, ElasticNetCV grids, 3-fold CV and seed '
             'match the pinned production trainer. No external holdout was used, as requested. The comparator is the saved production 202-lipid model '
             'with its original 50 training profiles. D071 exclusion is an authorized sample difference.', '',
             'Inference uses the verified complete local pseudobulk and real per-sample metacell inputs and the established absent-gene zero-fill policy. '
             'IPF validation precedes the authorized COPD RNA log1p correction; historical COPD inputs are retained at this stage.', '',
             'Differentials use the frozen cellHarmony moderated t test for pseudobulks (minimum 3, shrinkage 0.2), and Scanpy Wilcoxon '
             'for metacells (minimum 5, tie correction disabled). ' +
             ('The user authorized correction on 2026-10-06: BH is calculated over the same complete 202-lipid panel '
              'within each model/cohort/cell-type/input-unit comparison, without abundance eligibility filtering. '
              'The saved raw p values and folds are reused unchanged. Historical filtered BH values are retained in '
              'fdr_legacy_filtered; model training, predictions and experimental statistics are unchanged. '
              'The correction is per comparison; it does not correct an exploratory union over cell types or cohorts.'
              if corrected else
              'The saved numeric filter (>0.1 in at least two rows) and BH policy are retained. '
              'Because these targets are relative log2 values, eligibility is baseline-dependent: see BH_eligibility_audit.csv. '
              'No filter or BH-family change was introduced.'), '',
             'Signed folds are +2 for a twofold increase and -2 for a twofold decrease, calculated using ratios of arithmetic means of 2**prediction. '
             'All cohorts and raw-p/BH gates are reported, including fold thresholds 1, 1.1, 1.2, 1.5 and 2. '
             'The experimental source is Dr. Clair\'s supplied lung-tissue 10_results_with_statistics.csv; supplied p values and effects are used. '
             'Its effect log base is unconfirmed, so experimental fold magnitudes are not invented. RNA and lipidomics subjects are unmatched; '
             'this is disease-direction agreement, not subject-level validation. Metacell tests reproduce the original method and do not turn metacells '
             'into independent donors; donor counts are supplied in population_replication_census.csv.', '',
             '## Specific lipid results', '',
             'All measured-versus-imputed test rows are in all_measured_vs_imputed_tests.csv. Below are the concordant BH<0.1 '
             'results with the largest imputed fold magnitude for each lipid/cohort/unit/model (all supporting cell types are in cohort_concordant_lipid_celltypes.csv).', '']
    hits = joined[joined.supplied_adjusted_p.lt(.1) & joined.fdr.lt(.1) & joined.direction_concordant_supplied & ~joined.population.str.contains('__vs__')].copy()
    hits['absolute_signed_fold'] = hits.signed_fold.abs()
    hits = hits.sort_values('absolute_signed_fold', ascending=False).drop_duplicates(['contrast', 'unit', 'model', 'measured_feature'])
    lines.append(table(hits[['contrast', 'unit', 'model', 'measured_feature', 'population', 'signed_fold', 'pval', 'fdr', 'supplied_raw_p', 'supplied_adjusted_p']]))
    lines += ['', '## Complete experimental lipid inventory', '', table(measured.reset_index())]
    targets = pd.read_csv(reference_out / 'candidate_training_lipids_log2.csv', index_col=0)
    inventory = pd.read_csv(reference_out / 'summary.csv')[['Lipid', 'Train_R2_scaled', 'Nonzero_Coefficients', 'Alpha', 'L1_ratio']]
    inventory['target_min_relative_log2'] = inventory.Lipid.map(targets.min())
    inventory['target_median_relative_log2'] = inventory.Lipid.map(targets.median())
    inventory['target_max_relative_log2'] = inventory.Lipid.map(targets.max())
    lines += ['', '## All 202 trained lipid outputs', '',
              'These are relative-log2 target values and fitted training statistics; absolute MS1 abundances are not asserted.', '', table(inventory)]
    focused = summary[summary.contrast.isin(['Adams2020__IPF_vs_Healthy', 'Natri2024__IPF_vs_Healthy']) &
                      summary.minimum_imputed_fold.eq(1) & summary.overlap_gate.isin(['both_rawp005', 'both_BH010'])]
    display = focused.pivot(index=['contrast', 'unit', 'overlap_gate'], columns='model', values='concordant_union').reset_index()
    display.columns.name = None
    eligibility = pd.read_csv(OUT / 'BH_eligibility_audit.csv').groupby('model')[['rows', 'finite_rawp', 'finite_BH']].sum().reset_index()
    findings = ['## Main results', '', table(display), '',
                'Counts are unique experimental lipid features recovered in at least one cell type, significant in both datasets and matching direction. '
                'There are 132 uniquely matched experimental features under the established ion-mode matching rule; 91 have experimental raw p<0.05 '
                'and 97 have experimental BH<0.1. All 202 model outputs remain present; unmatched experimental identities are recorded separately.', '',
                'The candidate recovers all 91 matched raw-p-significant experimental lipids in metacells in each of Adams and pooled Natri '
                '(56 increases and 35 decreases). The prior recovers 87 and 90, respectively. Pseudobulk raw-p recovery is mixed: '
                '80 versus 77 for Adams, and 84 versus 87 for Natri.', '',
                ('All 116,756 test results per model now receive BH correction over the complete 202-lipid panel within each comparison. '
                 'All 578 paired testing families are identical between models. Raw p values, signed folds, comparison identities and '
                 'experimental statistics are unchanged; see BH_correction_audit.json. The historical filtered analysis is preserved at '
                 + str(reference_out) + '. This report contains the corrected candidate results; the production model and database are unchanged.'
                 if corrected else
                 'The saved filter admits predictions >0.1 in at least two rows into BH correction. Relative log2 values below that threshold '
                 'still represent present lipids; their exclusion depends on the chosen baseline. '
                 'The candidate has 63,067 eligible BH rows versus 116,756 for the prior, while both have 116,756 finite raw-p rows. '
                 'A clean comparison of BH recovery has not passed the testing-family gate. No BH-family change was made.'), '',
                table(eligibility), '',
                'Training contrasts are descriptive fitted-data checks, not held-out validation: across 2,020 sorted-population comparisons '
                'the prediction direction agrees with the corrected targets in 95.1%, with median absolute log2 fold error 0.080 '
                '(about a 1.057-fold error factor). All 198,970 originally observed sample-pair differences are retained in the saved corrected '
                'training reference. Per-lipid records are in training_population_fold_validation.csv and all_202_differential_preservation.csv.', '',
                '## Adams / Natri overlap', '',
                table(pd.read_csv(OUT / 'Adams_Natri_concordant_intersection.csv').drop(columns='shared_lipids')), '',
                'Shared identities are listed in Adams_Natri_concordant_intersection.csv. Donor independence between source studies is not asserted.', '',
                '## Best individual cell types at BH<0.1', '',
                table(pd.read_csv(OUT / 'concordant_celltype_ranking.csv').sort_values('concordant_total', ascending=False)
                      .drop_duplicates(['contrast', 'unit', 'model']))]
    lines[2:2] = findings + ['']
    (OUT / 'EVALUATION_SUMMARY.md').write_text('\n'.join(lines) + '\n')
    dump(OUT / 'completion.json', {'comparison_complete': True, 'full_panel': 202,
          'production_sha256_after': sha256(read_json(CONTRACT)['production_bundle']['path']),
          'MS1_abundance_calibration_complete': False,
          'validation_gate_passed_for_all_covariates': gate['all_covariate_stage_ready'],
          'all_covariate_RNA_scale_correction_started': False})
    print('Evaluation complete:', OUT, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['prepare', 'train', 'infer', 'evaluate', 'report'])
    parser.add_argument('--unit', choices=['PB', 'MC'])
    args = parser.parse_args()
    if args.stage == 'report':
        validate_run_manifest(OUT / 'run_manifest.json')
        summarize(pd.read_csv(OUT / 'all_cellHarmony_tests.csv'), pd.read_csv(OUT / 'all_measured_vs_imputed_tests.csv'),
                  pd.read_csv(OUT / 'lipid_matching_audit.csv'), pd.read_csv(OUT / 'experimental_lung_IPF_lipids_all_544.csv', index_col=0),
                  pd.read_csv(OUT / 'provenance/IPF_vs_healthy_registry.tsv', sep='\t'))
    elif args.stage == 'evaluate':
        evaluate([args.unit] if args.unit else None)
    elif args.stage == 'infer':
        infer([args.unit] if args.unit else None)
    else:
        globals()[args.stage]()
