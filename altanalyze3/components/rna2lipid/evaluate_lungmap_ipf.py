"""Isolated LungMAP candidate inference and frozen cellHarmony IPF evaluation.

Never writes to production bundles, imputed matrices, differential trees or databases.
Uses the saved cellHarmony statistical functions, including its EB and BH policies.
"""
from __future__ import annotations

try:
    from .candidate_integrity import reject_retired_workflow
except ImportError:
    from candidate_integrity import reject_retired_workflow

if __name__ == "__main__":
    reject_retired_workflow('evaluate_lungmap_ipf.py')


import argparse
import contextlib
import hashlib
import importlib.util
import json
import pickle
import shutil
import sys
import time
import warnings
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import logsumexp
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO.parent))
from reference_normalization import fit_bundle
from ipf_utils import canonical
from ipf_lung_targets import mode

TRANSFER = Path('/Users/saljh8/Dropbox/Transfer/CellRef2.0')
ATLAS = Path('/Volumes/salomonis2/LungMAP/CellRef2')
PRIOR = HERE / 'artifacts/IPF_candidate_validation_20261001'
CANDIDATES = {'candidate_MS1_47': 'supplied_bulk_MS1_47',
              'candidate_bulk_219_exploratory': 'bulk_gauge_219_exploratory'}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, default=str, allow_nan=False) + '\n')


def norm_counts(x):
    x = sparse.csr_matrix(x, dtype=np.float64)
    totals = np.asarray(x.sum(1)).ravel()
    scale = np.divide(1e4, totals, out=np.zeros_like(totals), where=totals > 0)
    x = sparse.diags(scale) @ x
    x.data = np.log2(1 + x.data)
    return x.tocsr()


def sum_groups(x, keys):
    names = pd.Index(sorted(set(keys)))
    rows = names.get_indexer(keys)
    d = sparse.csr_matrix((np.ones(len(rows)), (rows, np.arange(len(rows)))),
                          shape=(len(names), len(rows)))
    return d @ x, names, rows, d


def resolve_obs(a, meta, unit, states):
    obs = a.obs.copy()
    if unit == 'MC':
        info = meta.drop_duplicates(['study_acronym', 'Sample']).set_index(['study_acronym', 'Sample'])
        pairs = pd.MultiIndex.from_arrays([obs.study_acronym.astype(str), obs.Sample.astype(str)])
        if not pairs.isin(info.index).all():
            raise ValueError('Metacells have unresolved study/sample metadata')
        obs['Study_internal'] = info.Study_internal.reindex(pairs).to_numpy()
        obs['Donor'] = info.Donor.reindex(pairs).to_numpy()
        obs['cell_state'] = obs.short_name.map(states)
    else:
        if meta.Library.duplicated().any():
            raise ValueError('Ambiguous Library identifiers in metadata')
        info = meta.set_index('Library')
        if not obs.Library.astype(str).isin(info.index).all():
            raise ValueError('Library pseudobulks have unresolved metadata')
        for c in ('Sample', 'Donor', 'Study_internal', 'study_acronym'):
            obs[c] = info[c].reindex(obs.Library.astype(str)).to_numpy()
        obs['cell_state'] = obs.cell_state.map(states)
    if obs.cell_state.isna().any() or obs.Donor.isna().any():
        raise ValueError('Missing cell-state or donor annotation')
    return obs


def fit_candidates(genes, out):
    reject_retired_workflow('evaluate_lungmap_ipf.py:fit_candidates')
    bundles, info = {}, {}
    for label, prior in CANDIDATES.items():
        directory = out / label
        directory.mkdir(exist_ok=True)
        x = pd.read_csv(PRIOR / prior / 'candidate_training_RNA.csv', index_col=0)
        y = pd.read_csv(PRIOR / prior / 'candidate_training_lipids_log2.csv', index_col=0)
        missing = x.columns.difference(genes).tolist()
        used = x.columns.intersection(genes).tolist()
        x = x[used]
        fit_bundle(x, y, 'LungMAP_candidate_all_reference_shared_genes',
                   pd.Series(0., index=y.columns), 100., directory / 'candidate_bundle.pkl', [])
        b = pickle.load((directory / 'candidate_bundle.pkl').open('rb'))
        b['metadata'].update(candidate_only=True, production_default_changed=False,
                             donor_disjoint_validation=False, heldout_samples=0,
                             missing_LungMAP_genes_omitted=missing,
                             input_scale='log2(1 + CP10k), common healthy-study additive RNA alignment',
                             prediction_policy='Evaluation only; no disease-outcome tuning',
                             training_samples=len(x))
        x.to_csv(directory / 'candidate_training_RNA.csv')
        y.to_csv(directory / 'candidate_training_lipids_log2.csv')
        with (directory / 'candidate_bundle.pkl').open('wb') as f:
            pickle.dump(b, f)
        # Exact affine form of the frozen-alpha linear kernel ridge, checked below.
        w = b['model'].X_fit_.T @ b['model'].dual_coef_
        b['_coef'] = w * b['scaler_y'].scale_[None, :] / b['scaler_x'].scale_[:, None]
        b['_intercept'] = b['scaler_y'].mean_ - b['scaler_x'].mean_ @ b['_coef']
        fast = x.to_numpy() @ b['_coef'] + b['_intercept']
        original = b['scaler_y'].inverse_transform(b['model'].predict(b['scaler_x'].transform(x)))
        err = float(abs(fast - original).max())
        if err > 1e-8:
            raise ValueError('Affine ridge prediction differs from the fitted bundle')
        bundles[label] = b
        info[label] = {'training_profiles': len(x), 'input_genes': len(used),
                       'omitted_genes': missing, 'output_lipids': len(y.columns),
                       'alpha': 100, 'heldout_profiles': 0, 'affine_prediction_max_error': err}
    return bundles, info


def write_predictions(out, label, unit, values, obs, columns):
    directory = out / label
    path = directory / f'LungMAP_{unit}_predictions_log2.h5ad'
    a = ad.AnnData(np.asarray(values, dtype=np.float64), obs=obs.copy(),
                   var=pd.DataFrame(index=pd.Index(columns, name='lipid')))
    a.uns.update(expression_scale='log2_abundance', inverse='2**X; no pseudocount or clipping',
                 candidate_only=True, production_default_changed=False,
                 input_unit=unit, population_col='cell_state')
    if unit == 'PB':
        a.uns['pseudobulk_method'] = 'pseudobulk'
    a.write_h5ad(path, compression='lzf')
    linear = a.copy()
    linear.X = np.exp2(linear.X)
    if not np.isfinite(linear.X).all() or not (linear.X > 0).all():
        raise ValueError('Invalid positive candidate abundance')
    linear.uns['expression_scale'] = 'relative_linear_abundance'
    linear.write_h5ad(directory / f'LungMAP_{unit}_predictions_linear.h5ad', compression='lzf')
    return a


def infer(out, meta, cross, registry):
    reject_retired_workflow('evaluate_lungmap_ipf.py:infer')
    pb_source = TRANSFER / 'reference_v7_top500/pseudobulk/CellRef2.0_pseudobulk_library_x_cellstate_UNION.h5ad'
    mc_source = TRANSFER / 'reference_v7_top500/metacells/CellRef2.0_metacells_v8_persample.h5ad'
    p = ad.read_h5ad(pb_source)
    genes = p.var_names.astype(str)
    states = dict(zip(cross.cell_state, cross.cell_type)) | dict(zip(cross.cell_type_full, cross.cell_type)) | dict(zip(cross.cell_type, cross.cell_type))
    po = resolve_obs(p, meta, 'PB', states)
    keys = po.Study_internal.astype(str) + '|' + po.Sample.astype(str) + '|' + po.cell_state.astype(str)
    raw, names, rows, design = sum_groups(p.X, keys)
    first = po.loc[~keys.duplicated()].copy()
    first.index = keys.loc[~keys.duplicated()]
    obs = first.loc[names].copy()
    obs['n_cells'] = np.asarray(design @ po.n_cells.to_numpy()).ravel()
    obs['replicate_id'] = names
    obs['sample_uid'] = obs.Study_internal.astype(str) + '|' + obs.Sample.astype(str)
    obs = obs.drop(columns=['Library', 'counts_type'], errors='ignore')
    obs.index.name = 'observation_id'
    if obs.duplicated(['Study_internal', 'Sample', 'cell_state']).any():
        raise ValueError('Duplicate sample/state pseudobulk')
    bundles, training = fit_candidates(genes, out)
    needed = list(next(iter(bundles.values()))['X_columns'])
    if not all(list(b['X_columns']) == needed for b in bundles.values()):
        raise ValueError('Candidate RNA panels differ')
    positions = genes.get_indexer(needed)
    expression = norm_counts(raw)[:, positions].toarray()
    # Real whole-sample healthy RNA provides gene baselines. Sum raw counts
    # over states before normalization, then give each donor equal weight.
    sample_raw, snames, _, _ = sum_groups(raw, obs.sample_uid)
    sfirst = obs.loc[~obs.sample_uid.duplicated()].set_index('sample_uid').loc[snames]
    sx = pd.DataFrame(norm_counts(sample_raw)[:, positions].toarray(), index=snames, columns=needed)
    minfo = meta.drop_duplicates(['Study_internal', 'Sample']).copy()
    minfo.index = minfo.Study_internal.astype(str) + '|' + minfo.Sample.astype(str)
    normal = minfo.disease_group.eq('normal')
    normal_ids = sx.index.intersection(minfo.index[normal])
    donor = sfirst.loc[normal_ids, 'Study_internal'].astype(str) + '|' + sfirst.loc[normal_ids, 'Donor'].astype(str)
    healthy = sx.loc[normal_ids].groupby(donor).mean()
    studies = sfirst.loc[normal_ids, 'Study_internal'].groupby(donor).first()
    global_mean = healthy.mean()
    offsets = {}
    output = {}
    for label, b in bundles.items():
        md = pd.read_csv(PRIOR / CANDIDATES[label] / 'candidate_training_sample_metadata.csv', index_col=0)
        tx = pd.read_csv(out / label / 'candidate_training_RNA.csv', index_col=0)
        reference_mean = tx.loc[md.dataset.eq('bulk')].mean()
        offsets[label] = {study: (reference_mean - healthy.loc[studies.eq(study)].mean()
                                  if studies.eq(study).any() else reference_mean - global_mean)
                          for study in meta.Study_internal.unique()}
        offset_frame = pd.DataFrame(offsets[label]).T
        offset_frame.to_csv(out / label / 'healthy_study_RNA_offsets.csv')
        vals = expression @ b['_coef'] + b['_intercept']
        shifts = offset_frame.loc[obs.Study_internal].to_numpy() @ b['_coef']
        vals += shifts
        output[label] = {'PB': write_predictions(out, label, 'PB', vals, obs, b['Y_columns'])}
    del p, raw, sample_raw, expression
    m = ad.read_h5ad(mc_source, backed='r')
    mo = resolve_obs(m, meta, 'MC', states)
    mo['replicate_id'] = mo.index.astype(str)
    mo['sample_uid'] = mo.Study_internal.astype(str) + '|' + mo.Sample.astype(str)
    mo.index.name = 'observation_id'
    mcvals = {label: np.empty((m.n_obs, len(b['Y_columns'])), dtype=np.float64)
              for label, b in bundles.items()}
    reader = h5py.File(mc_source, 'r')
    node = reader['X']
    indptr = node['indptr'][:]
    for start in range(0, m.n_obs, 1024):
        stop = min(start + 1024, m.n_obs)
        left, right = int(indptr[start]), int(indptr[stop])
        raw_block = sparse.csr_matrix((node['data'][left:right], node['indices'][left:right],
                                       indptr[start:stop+1] - left), shape=(stop-start, m.n_vars))
        block = norm_counts(raw_block)[:, positions].toarray()
        for label, b in bundles.items():
            off = pd.DataFrame(offsets[label]).T.loc[mo.Study_internal.iloc[start:stop]].to_numpy()
            mcvals[label][start:stop] = (block + off) @ b['_coef'] + b['_intercept']
        if start % 20480 == 0:
            print(f'Metacell inference {stop:,}/{m.n_obs:,}', flush=True)
    reader.close()
    m.file.close()
    for label, b in bundles.items():
        output[label]['MC'] = write_predictions(out, label, 'MC', mcvals[label], mo, b['Y_columns'])
    for label, b in bundles.items():
        path = out / label / 'candidate_bundle.pkl'
        saved = pickle.load(path.open('rb'))
        saved['LungMAP_RNA_calibration'] = {'offsets_by_Study_internal': {
            study: values.to_dict() for study, values in offsets[label].items()},
            'policy': 'Healthy donor-balanced whole-sample RNA offsets; same offset in both disease arms',
            'fallback': 'All healthy donors for studies without measured normal donors'}
        with path.open('wb') as f:
            pickle.dump(saved, f)
    dump(out / 'training_and_input_audit.json', {'training': training,
        'PB_raw_library_rows': len(po), 'PB_sample_state_rows': len(obs), 'MC_rows': len(mo),
        'healthy_calibration_samples': len(normal_ids), 'healthy_calibration_donors': len(healthy),
        'control_donors_by_study': studies.value_counts().to_dict(),
        'RNA_normalization': 'log2(1 + raw_counts / full_transcriptome_count_sum * 10000)',
        'missing_RNA_policy': 'Refit all-reference candidates on 1293 shared genes; omit GPX1 and JMJD7-PLA2G4B, no zero fill',
        'source_PB': str(pb_source), 'source_MC': str(mc_source)})
    return output


def load_current(out, cross):
    result = {}
    mapping = dict(zip(cross.cell_state, cross.cell_type)) | dict(zip(cross.cell_type_full, cross.cell_type)) | dict(zip(cross.cell_type, cross.cell_type))
    for unit, name in [('PB', 'pseudobulk'), ('MC', 'metacell')]:
        source = ATLAS / f'imputed_v8/cellref2_v8_{name}_lipid_forDE.h5ad'
        dest = out / 'current_inputs' / source.name
        dest.parent.mkdir(exist_ok=True)
        if not dest.exists():
            shutil.copyfile(source, dest)
        a = ad.read_h5ad(dest)
        if sparse.issparse(a.X):
            a.X = a.X.toarray()
        if 'cell_state' not in a.obs:
            a.obs['cell_state'] = a.obs.short_name.map(mapping)
        else:
            a.obs['cell_state'] = a.obs.cell_state.map(mapping)
        a.obs.index.name = 'observation_id'
        result[unit] = a
    return result


def load_engine(out):
    source = ATLAS / 'code_v8/repair_codex/frozen/cellHarmony_differential.py'
    dest = out / 'provenance/cellHarmony_differential_frozen.py'
    dest.parent.mkdir(exist_ok=True)
    if not dest.exists():
        shutil.copyfile(source, dest)
    spec = importlib.util.spec_from_file_location('lungmap_frozen_cellharmony', dest)
    engine = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(engine)
    if not hasattr(pd.Series, 'nonzero'):
        pd.Series.nonzero = lambda self: np.asarray(self).nonzero()
    return engine


def pairs_for(obs, row, overrides, minimum):
    uid = obs.Study_internal.astype(str) + '|' + obs.Sample.astype(str)
    def keys(value):
        return {'|'.join(k.split('|')[-2:]) for k in str(value).split(';') if k and k != 'nan'}
    ck, nk = keys(row.case_keys), keys(row.control_keys)
    if ck & nk:
        raise ValueError('Case/control arms overlap')
    case, control = uid.isin(ck).to_numpy(), uid.isin(nk).to_numpy()
    rows = []
    for state in obs.cell_state.astype(str).unique():
        mask = obs.cell_state.astype(str).eq(state).to_numpy()
        rows.append((state, state, state, case & mask, control & mask, False))
    for ov in overrides.itertuples():
        d, r = ov.cell2_disease_state, ov.cell1_reference_state
        rows.append((d + '__vs__' + r, d, r,
                     case & obs.cell_state.eq(d).to_numpy(), control & obs.cell_state.eq(r).to_numpy(), True))
    accepted, census = [], []
    for label, d, r, ca, co, is_override in rows:
        status = 'tested' if ca.sum() >= minimum and co.sum() >= minimum else 'insufficient_replicates'
        census.append({'population': label, 'case_state': d, 'control_state': r,
                       'n_case': int(ca.sum()), 'n_control': int(co.sum()),
                       'n_case_donors': obs.loc[ca, ['Study_internal', 'Donor']].drop_duplicates().shape[0],
                       'n_control_donors': obs.loc[co, ['Study_internal', 'Donor']].drop_duplicates().shape[0],
                       'override': is_override, 'status': status})
        if status == 'tested':
            accepted.append((label, ca, co))
    return accepted, pd.DataFrame(census)


def arithmetic_stats(values, ncase, columns):
    values = np.asarray(values, dtype=np.float64)
    lc = logsumexp(values[:ncase] * np.log(2), axis=0) / np.log(2) - np.log2(ncase)
    ln = logsumexp(values[ncase:] * np.log(2), axis=0) / np.log(2) - np.log2(len(values) - ncase)
    return pd.DataFrame({'log2fc': lc - ln, 'case_mean_expr': lc, 'control_mean_expr': ln,
                         'geometric_log2fc': values[:ncase].mean(0) - values[ncase:].mean(0)}, index=columns)


def differential(a, unit, row, overrides, engine, out, label):
    folder = out / 'differentials' / row.contrast_id / unit / label
    folder.mkdir(parents=True, exist_ok=True)
    pairs, census = pairs_for(a.obs, row, overrides, 3 if unit == 'PB' else 5)
    census.to_csv(folder / 'population_replication_census.csv', index=False)
    all_rows, summary = [], []
    explicit_log = label != 'current_as_deployed'
    with (folder / 'cellHarmony.log').open('w') as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        for population, ca, co in pairs:
            selected = np.r_[np.flatnonzero(ca), np.flatnonzero(co)]
            values = np.asarray(a.X[selected])
            condition = pd.Categorical(['CASE'] * int(ca.sum()) + ['CONTROL'] * int(co.sum()))
            block = ad.AnnData(values if explicit_log else sparse.csr_matrix(values), obs=pd.DataFrame({'Condition': condition}, index=[str(i) for i in selected]),
                               var=a.var.copy())
            if explicit_log:
                # Stops Scanpy's raw-count heuristic. Its internal approximate FC
                # is discarded; compute the exact abundance ratio below.
                block.uns['log1p'] = {'base': 2.0}
            if unit == 'PB':
                frame, tested = engine._moderated_t_test(block, 'Condition', 'CASE', 'CONTROL', population)
                frame = frame.set_index('gene')
            else:
                names, fdr, _, pvals = engine._rank_genes_scanpy(block, 'Condition', 'CASE', 'CONTROL', 'wilcoxon')
                frame = pd.DataFrame({'pval': pvals, 'fdr': fdr}).reindex(names)
                tested = len(frame)
            exact = arithmetic_stats(values, int(ca.sum()), a.var_names)
            if explicit_log:
                frame = frame.drop(columns=['log2fc'], errors='ignore').join(exact)
            else:
                old = engine._compute_pseudobulk_log2fc(block, 'Condition', 'CASE', 'CONTROL')
                frame = frame.drop(columns=['log2fc'], errors='ignore').join(old[['log2fc', 'case_mean_log2', 'control_mean_log2']].rename(
                    columns={'case_mean_log2': 'case_mean_expr', 'control_mean_log2': 'control_mean_expr'}))
                frame['geometric_log2fc'] = exact.geometric_log2fc
                frame['correct_log2_abundance_arithmetic_log2fc'] = exact.log2fc
            frame.index.name = 'gene'
            frame['population'], frame['n_case'], frame['n_control'] = population, int(ca.sum()), int(co.sum())
            frame['signed_fold'] = np.where(frame.log2fc > 0, np.exp2(frame.log2fc), -np.exp2(-frame.log2fc))
            frame.loc[frame.log2fc.eq(0), 'signed_fold'] = 0.
            frame['significant_raw_p_005'] = frame.pval.lt(.05) & frame.log2fc.ne(0)
            frame['significant_BH_005'] = frame.fdr.lt(.05) & frame.log2fc.ne(0)
            frame['significant_BH_010'] = frame.fdr.lt(.1) & frame.log2fc.ne(0)
            all_rows.append(frame.reset_index())
            summary.append({'population': population, 'n_case': int(ca.sum()), 'n_control': int(co.sum()),
                            'tested_lipids': int(tested), 'raw_p_005_calls': int(frame.significant_raw_p_005.sum()),
                            'BH_005_calls': int(frame.significant_BH_005.sum()), 'BH_010_calls': int(frame.significant_BH_010.sum())})
    combined = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame(columns=['gene','population','log2fc','pval','fdr'])
    combined['contrast'], combined['unit'], combined['model'] = row.contrast_id, unit, label
    combined.to_csv(folder / 'all_tested_lipids.csv', index=False)
    combined.loc[combined.pval.lt(.05)].to_csv(folder / 'raw_p_005_differentials.csv', index=False)
    pd.DataFrame(summary).to_csv(folder / 'differential_counts.csv', index=False)
    dump(folder / 'parameters.json', {'engine': 'Frozen deployed cellHarmony statistical functions',
        'test': 'eBayes moderated t' if unit == 'PB' else 'Scanpy Wilcoxon, tie_correct=False',
        'alpha': .05, 'fc': 1, 'min_replicates': 3 if unit == 'PB' else 5,
        'metacell_cap': None, 'explicit_log2_abundance': explicit_log,
        'BH_scope': 'All independently filtered features in the model panel, per population/contrast'})
    return combined, census


def lipid_mapping(features, label, measured, native):
    rows = []
    for f in features:
        key = canonical(f)
        ion = mode(f)
        if ion is None:
            parents = [x for x in native if x in (f + '_P', f + '_N')]
            if len(parents) == 1:
                ion = mode(parents[0])
        matches = measured.index[measured.canonical.eq(key) & measured.ion_mode.eq(ion)] if key and ion else []
        rows.append({'model': label, 'gene': f, 'canonical': key, 'ion_mode': ion,
                     'measured_feature': matches[0] if len(matches) == 1 else None,
                     'status': 'matched' if len(matches) == 1 else 'absent_or_ambiguous'})
    return pd.DataFrame(rows)


def main():
    reject_retired_workflow('evaluate_lungmap_ipf.py:main')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', default=str(HERE / 'artifacts/LungMAP_IPF_candidate_20261004'))
    parser.add_argument('--reuse-imputation', action='store_true')
    parser.add_argument('--refresh-current', action='store_true', help='Recompute only the deployed baseline; reuse existing candidate DE tables')
    args = parser.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    protected = [HERE / 'rna2lipid_hs_lung_lipidwise_bundle.pkl', PRIOR / 'supplied_bulk_MS1_47/candidate_bundle.pkl']
    protected_before = {str(p): sha(p) for p in protected}
    source_meta = TRANSFER / 'pseudobulk/study_metadata/GPT5.6-sol/harmonized_library_metadata_harmonized_final_corrected_v6.txt'
    meta = pd.read_csv(source_meta, sep='\t', dtype=str).fillna('')
    frozen = ATLAS / 'code_v8/repair_codex/frozen'
    registry = pd.read_csv(frozen / 'registry.tsv', sep='\t').fillna('')
    # Only IPF-versus-healthy arms; other diseases and within-IPF contrasts
    # are not external validations of a whole-lung IPF/control experiment.
    registry = registry[registry.contrast_id.str.contains(r'__IPF(?:_.*)?_vs_Healthy')].copy()
    cross = pd.read_csv(frozen / 'crosswalk.tsv', sep='\t')
    overrides = pd.read_csv(frozen / 'overrides.tsv', sep='\t')
    (out / 'provenance').mkdir(exist_ok=True)
    registry.to_csv(out / 'provenance/IPF_vs_healthy_registry.tsv', sep='\t', index=False)
    for source in [source_meta, frozen / 'crosswalk.tsv', frozen / 'overrides.tsv']:
        shutil.copyfile(source, out / 'provenance' / source.name)
    start = time.time()
    with threadpool_limits(limits=2):
        if args.reuse_imputation:
            candidates = {label: {u: ad.read_h5ad(out / label / f'LungMAP_{u}_predictions_log2.h5ad') for u in ('PB','MC')}
                          for label in CANDIDATES}
        else:
            candidates = infer(out, meta, cross, registry)
        print('Candidate inference complete; loading deployed predictions', flush=True)
        current = load_current(out, cross)
        engine = load_engine(out)
        panels = {'current_as_deployed': current, 'current_explicit_log2': current, **candidates}
        measured = pd.read_csv(PRIOR / 'measured_lung_IPF_vs_control.csv', index_col=0)
        native = pd.read_csv(HERE / 'artifacts/reference_normalization_20261001/native_lipid_targets.csv', index_col=0).columns
        maps = pd.concat([lipid_mapping(a['PB'].var_names, label, measured, native) for label,a in panels.items()], ignore_index=True)
        maps.to_csv(out / 'lipid_matching_audit.csv', index=False)
        data, census, parity = [], [], []
        for row in registry.itertuples():
            ov = overrides[overrides.contrast_id.eq(row.contrast_id)]
            for unit in ('PB', 'MC'):
                for label, matrices in panels.items():
                    print(f'DE {row.contrast_id} {unit} {label}', flush=True)
                    if args.refresh_current and label != 'current_as_deployed':
                        directory = out / 'differentials' / row.contrast_id / unit / label
                        frame = pd.read_csv(directory / 'all_tested_lipids.csv')
                        counts = pd.read_csv(directory / 'population_replication_census.csv')
                    else:
                        frame, counts = differential(matrices[unit], unit, row, ov, engine, out, label)
                    data.append(frame)
                    counts['contrast'], counts['unit'], counts['model'] = row.contrast_id, unit, label
                    census.append(counts)
                    if label == 'current_as_deployed':
                        source = ATLAS / f'differential_v8/{row.contrast_id}/{unit}/lipid/CASE_vs_CONTROL/DEGs/DEG_detailed_CASE_vs_CONTROL.tsv'
                        if source.exists():
                            old = pd.read_csv(source, sep='\t')
                            dest = out / 'archived_current_differentials' / row.contrast_id / unit
                            dest.mkdir(parents=True, exist_ok=True); shutil.copyfile(source, dest / source.name)
                            joined = old.merge(frame, on=['gene','population'], suffixes=('_saved','_rerun'))
                            saved_calls = set(zip(old.gene, old.population))
                            called = frame[frame.pval.lt(.05)]
                            rerun_calls = set(zip(called.gene, called.population))
                            parity.append({'contrast': row.contrast_id, 'unit': unit, 'saved_rawp_calls': len(old),
                                'rerun_rawp_calls': len(rerun_calls), 'lost_saved_calls': len(saved_calls-rerun_calls),
                                'added_rerun_calls': len(rerun_calls-saved_calls),
                                'matched_rerun_rows': len(joined), 'pval_max_abs_error': float(abs(joined.pval_saved-joined.pval_rerun).max()) if len(joined) else None,
                                'log2fc_max_abs_error': float(abs(joined.log2fc_saved-joined.log2fc_rerun).max()) if len(joined) else None,
                                'n_case_disagreements': int(joined.n_case_saved.ne(joined.n_case_rerun).sum()),
                                'n_control_disagreements': int(joined.n_control_saved.ne(joined.n_control_rerun).sum())})
        full = pd.concat(data, ignore_index=True)
        full.to_csv(out / 'all_cellHarmony_tests.csv', index=False)
        pd.concat(census, ignore_index=True).to_csv(out / 'population_replication_census.csv', index=False)
        parity_frame = pd.DataFrame(parity)
        parity_frame.to_csv(out / 'current_saved_vs_rerun_parity.csv', index=False)
        parity_frame[['contrast','unit','saved_rawp_calls','rerun_rawp_calls','lost_saved_calls','added_rerun_calls']].rename(
            columns={'saved_rawp_calls':'saved_calls','rerun_rawp_calls':'rerun_calls','lost_saved_calls':'lost',
                     'added_rerun_calls':'added'}).to_csv(out / 'current_saved_call_set_audit.csv', index=False)
        match = maps[maps.status.eq('matched')]
        joined = full.merge(match[['model','gene','measured_feature']], on=['model','gene']).merge(
            measured.rename_axis('measured_feature').reset_index()[['measured_feature','original_annotation','supplied_log_effect','supplied_raw_p','supplied_adjusted_p','geometric_log2FC','pvalue','FDR']],
            on='measured_feature', suffixes=('', '_measured'))
        joined['direction_concordant_supplied'] = joined.log2fc * joined.supplied_log_effect > 0
        joined['direction_concordant_donor'] = joined.log2fc * joined.geometric_log2FC > 0
        joined.to_csv(out / 'all_measured_vs_imputed_tests.csv', index=False)
        summarize(out, joined, maps, measured, registry)
    protected = [HERE / 'rna2lipid_hs_lung_lipidwise_bundle.pkl', PRIOR / 'supplied_bulk_MS1_47/candidate_bundle.pkl']
    import scanpy, sklearn, scipy
    dump(out / 'run_manifest.json', {'candidate_only': True, 'official_model_replaced': False,
        'runtime': {'python':sys.version,'executable':sys.executable,'scanpy':scanpy.__version__,'sklearn':sklearn.__version__,'numpy':np.__version__,'scipy':scipy.__version__,'pandas':pd.__version__},
        'elapsed_seconds': round(time.time()-start, 1), 'contrasts': registry.contrast_id.tolist(),
        'source_metadata_sha256': sha(source_meta), 'engine_sha256': sha(out / 'provenance/cellHarmony_differential_frozen.py'),
        'protected_bundle_sha256_before': protected_before,
        'protected_bundle_sha256_after': {str(p):sha(p) for p in protected},
        'measured_source': '/Users/saljh8/Downloads/Lipidomics/IPF/10_results_with_statistics.csv',
        'source_table_sha256': sha('/Users/saljh8/Downloads/Lipidomics/IPF/10_results_with_statistics.csv')})
    if protected_before != {str(p): sha(p) for p in protected}:
        raise ValueError('Protected production/prior candidate bundle changed during evaluation')
    print(f'Completed: {out}', flush=True)


def summarize(out, joined, maps, measured, registry):
    summary = []
    scopes = [('supplied_rawp005', 'supplied_raw_p', .05), ('supplied_BH010','supplied_adjusted_p',.1),
              ('donor_rawp005','pvalue',.05), ('donor_BH010','FDR',.1)]
    matched_sets = {label:set(g.measured_feature.dropna()) for label,g in maps.groupby('model')}
    shared47 = matched_sets['current_as_deployed'] & matched_sets['candidate_MS1_47']
    shared219 = matched_sets['current_as_deployed'] & matched_sets['candidate_bulk_219_exploratory']
    selections = [('all_matched',None), ('shared_MS1_47',shared47), ('shared_bulk_219',shared219)]
    for contrast in registry.contrast_id:
      for unit in ('PB','MC'):
       for label in matched_sets:
        g = joined[joined.contrast.eq(contrast) & joined.unit.eq(unit) & joined.model.eq(label)]
        for source_name,col,cut in scopes:
            direction = 'direction_concordant_donor' if source_name.startswith('donor') else 'direction_concordant_supplied'
            source_sig = set(measured.index[measured[col].lt(cut)])
            for panel, allowed in selections:
                if panel == 'shared_MS1_47' and label == 'candidate_bulk_219_exploratory':continue
                if panel == 'shared_bulk_219' and label == 'candidate_MS1_47':continue
                available = matched_sets[label] if allowed is None else matched_sets[label] & allowed
                for include_overrides in [False,True]:
                    rows = g[g.measured_feature.isin(source_sig & available)]
                    if not include_overrides:
                        rows = rows[~rows.population.str.contains('__vs__')]
                    for gate,field,value in [('rawp005','pval',.05),('BH005','fdr',.05),('BH010','fdr',.1)]:
                        calls = rows[rows[field].lt(value)]
                        yes, no = calls[calls[direction]], calls[~calls[direction]]
                        matched_tested = set(rows.measured_feature)
                        up = set(yes.loc[yes.log2fc.gt(0),'measured_feature'])
                        down = set(yes.loc[yes.log2fc.lt(0),'measured_feature'])
                        summary.append({'contrast':contrast,'unit':unit,'model':label,'source_gate':source_name,
                            'panel':panel,'include_state_overrides':include_overrides,'prediction_gate':gate,
                            'experimentally_significant_total':len(source_sig),'eligible_matched_lipids':len(source_sig & available),
                            'matched_tested_lipids':len(matched_tested),'concordant_up_lipids':len(up),'concordant_down_lipids':len(down),
                            'concordant_unique_lipids':len(up|down),'discordant_unique_lipids':no.measured_feature.nunique(),
                            'both_direction_lipids':len((up|down)&set(no.measured_feature)),
                            'concordant_lipid_state_results':len(yes),'discordant_lipid_state_results':len(no),
                            'tested_population_comparisons':rows.population.nunique(),
                            'status': 'tested' if len(g) else 'no_testable_population; fold_only'})
    pd.DataFrame(summary).to_csv(out / 'directional_validation_summary.csv', index=False)
    measured.to_csv(out / 'experimental_lung_IPF_lipids_all_544.csv')


if __name__ == '__main__':
    warnings.filterwarnings('ignore', message='.*Trying to unpickle estimator.*')
    main()
