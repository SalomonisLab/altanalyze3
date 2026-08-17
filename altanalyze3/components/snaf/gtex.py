#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import h5py
import matplotlib.pyplot as plt
import anndata as ad
from scipy.optimize import minimize, minimize_scalar
from scipy import stats
from scipy.sparse import csr_matrix, find, hstack as sparse_hstack, issparse
from tqdm import tqdm
import re
import logging

logger = logging.getLogger(__name__)

# try:
#     import pymc3 as pm   # conda install -c conda-forge pymc3 mkl-service
#     import theano
#     import arviz as az
# except ImportError:
#     print('''
#         Optional package pymc3 is not installed, it is for calculating tumor specificity using hirerarchical bayesian model
#         For Linux: https://github.com/pymc-devs/pymc/wiki/Installation-Guide-(Linux)
#         For MacOS: https://github.com/pymc-devs/pymc/wiki/Installation-Guide-(MacOS)
#         For PC:    https://github.com/pymc-devs/pymc/wiki/Installation-Guide-(Windows)
#     ''')

# '''
# this script is to query the tumor specificity of the junction
# '''


def gtex_configuration(df,gtex_db,t_min_arg,n_max_arg,normal_cutoff_arg,tumor_cutoff_arg,normal_prevalance_cutoff_arg,tumor_prevalance_cutoff_arg,add_control=None,control_stats_path=None,use_summary='auto'):
    global adata_gtex
    global adata
    global t_min
    global n_max
    global normal_cutoff
    global tumor_cutoff
    global normal_prevalance_cutoff
    global tumor_prevalance_cutoff
    tested_junctions = set(df.index)

    # ---- Summary-backed control path (the precompute architecture) ----------------------
    # If a per-junction stats table exists for this control (built ONCE, offline, by
    # snaf-precompute-control), load ONLY that small table -- mean/std/mle/normal_prevalence
    # [+ BayesTS] -- and NEVER open the multi-GB count matrix. Sifting reads obs['mean'] and
    # obs['normal_prevalence']; tumor specificity reads obs['mle']; BayesTS is already done.
    # Skipped when add_control is given (those still need their full matrices to combine).
    if use_summary is not False and add_control is None:
        from .control_stats import default_stats_path, load_control_stats
        _sp = control_stats_path if control_stats_path is not None else default_stats_path(gtex_db)
        _stats = load_control_stats(_sp)
        if _stats is not None:
            keep = [u for u in _stats.index if u in tested_junctions]
            sub = _stats.loc[keep]
            obs = pd.DataFrame(index=sub.index)
            for c in ('mean','std','mle','normal_prevalence','bayests_sigma','bayests_percentile'):
                if c in sub.columns:
                    obs[c] = sub[c].values
            var = pd.DataFrame(index=['_summary']); var['tissue'] = 'summary'; var['total_count'] = 1.0
            adata_gtex = ad.AnnData(X=csr_matrix((sub.shape[0], 1), dtype=np.float32), obs=obs, var=var)
            adata_gtex.uns['snaf_summary_backed'] = True
            adata_gtex.uns['snaf_summary_cutoff'] = float(normal_cutoff_arg)
            adata = adata_gtex
            t_min = t_min_arg; n_max = n_max_arg
            normal_cutoff = normal_cutoff_arg; tumor_cutoff = tumor_cutoff_arg
            normal_prevalance_cutoff = normal_prevalance_cutoff_arg; tumor_prevalance_cutoff = tumor_prevalance_cutoff_arg
            print('Loaded control SUMMARY ({}/{} tested junctions) from {} -- full count matrix NOT loaded'.format(
                sub.shape[0], len(tested_junctions), _sp))
            return adata_gtex
    # -------------------------------------------------------------------------------------
    # Read the (potentially multi-GB) control DB in backed mode and subset to the
    # tested junctions ON DISK, so only the intersecting rows are pulled into RAM
    # (the whole point of passing `df` here). Falls back to a full read if the
    # h5ad cannot be opened backed.
    # Bulk-read the control into memory then subset in-memory. anndata's BACKED
    # fancy-indexing over ~1M labels does per-row h5ad seeks and is pathologically slow
    # (26+ min on a 1M-row intersection) whereas a full read + in-memory subset is ~4s.
    # Only fall back to on-disk backed subsetting if the control does not fit in RAM.
    try:
        adata = ad.read_h5ad(gtex_db)
        adata = adata[np.logical_not(adata.obs_names.duplicated()), :]
        keep = [o for o in adata.obs_names if o in tested_junctions]   # preserves DB order, deterministic
        adata = adata[keep, :].copy()
    except MemoryError:
        adata = ad.read_h5ad(gtex_db, backed='r')
        adata = adata[np.logical_not(adata.obs_names.duplicated()), :]
        keep = [o for o in adata.obs_names if o in tested_junctions]
        adata = adata[keep, :].to_memory()
    print('Current loaded gtex cohort with shape {}'.format(adata.shape))
    tissue_dict = adata.var['tissue'].to_dict()
    adata_gtex = adata   # already has mean and tissue variables
    if add_control is not None:
        left_X = csr_matrix(adata.X)
        left_rows = adata.obs_names.tolist()
        left_cols = adata.var_names.tolist()
        for id_, control in add_control.items():
            if isinstance(control,pd.DataFrame):
                assert len(set(control.columns).intersection(tissue_dict.keys())) == 0  # sample id can not be ambiguous
                control = control.loc[np.logical_not(control.index.duplicated()),:]
                control = control.loc[list(set(control.index).intersection(tested_junctions)),:]
                print('Adding cohort {} with shape {} to the database'.format(id_,control.shape))
                tissue_dict_right = {k:id_ for k in control.columns}
                right_X = csr_matrix(control.values.astype(np.float32))
                right_rows = control.index.tolist(); right_cols = control.columns.tolist()

            elif isinstance(control,ad.AnnData):
                assert len(set(control.var_names).intersection(tissue_dict.keys())) == 0
                control = control[np.logical_not(control.obs_names.duplicated()),:]
                control = control[list(set(control.obs_names).intersection(tested_junctions)),:]
                print('Adding cohort {} with shape {} to the database'.format(id_,control.shape))
                if 'tissue' in control.var.columns:   # if tissue is in var columns, it will be used
                    tissue_dict_right = control.var['tissue'].to_dict()
                else:
                    tissue_dict_right = {k:id_ for k in control.var_names}
                right_X = csr_matrix(control.X)
                right_rows = control.obs_names.tolist(); right_cols = control.var_names.tolist()

            else:
                raise Exception('control must be either in dataframe or anndata format')

            tissue_dict.update(tissue_dict_right)
            left_X, left_rows, left_cols = _combine_two_sparse(left_X, left_rows, left_cols, right_X, right_rows, right_cols)
            print('now the shape of control db is {}'.format((len(left_rows), len(left_cols))))

        adata = ad.AnnData(X=left_X, obs=pd.DataFrame(index=left_rows), var=pd.DataFrame(index=left_cols))
        adata.var['tissue'] = adata.var_names.map(tissue_dict).values
        adata.obs['mean'] = np.asarray(adata.X.mean(axis=1)).ravel()
        adata.var['total_count'] = np.asarray(adata.X.sum(axis=0)).ravel() / 1e6


    # ---- BayesTS: always available, never recomputed -----------------------------------
    # The bundled per-junction scores for GTEx and Tabula Sapiens ship with SNAF, so the
    # BayesTS tumour-specificity filter works in EVERY configuration -- including the
    # add_control path above, which skips the summary table and therefore used to leave
    # adata_gtex without a 'bayests_percentile' column, silently disabling
    # --max_bayests_percentile. Attaching it here costs one indexed read of a ~50 MB table.
    if 'bayests_percentile' not in adata_gtex.obs.columns:
        try:
            from .control_stats import load_bundled_bayests
            _bt = load_bundled_bayests(uids=set(adata_gtex.obs_names))
            if _bt is not None and _bt.shape[0] > 0:
                adata_gtex.obs['bayests_sigma'] = _bt['bayests_sigma'].reindex(adata_gtex.obs_names).values
                adata_gtex.obs['bayests_percentile'] = _bt['bayests_percentile'].reindex(adata_gtex.obs_names).values
                _n = int(np.isfinite(adata_gtex.obs['bayests_percentile'].values).sum())
                print('BayesTS scores attached from the bundled references: {} / {} control '
                      'junctions scored (junctions with no BayesTS score are treated as '
                      'maximally tumor-specific and are KEPT)'.format(_n, adata_gtex.n_obs))
            else:
                logger.warning('No bundled BayesTS reference was loaded; '
                               '--max_bayests_percentile will have NO effect this run.')
        except Exception as _e:
            logger.warning('Could not attach bundled BayesTS scores (%s); '
                           '--max_bayests_percentile will have NO effect this run.', _e)

    t_min = t_min_arg
    n_max = n_max_arg
    normal_cutoff = normal_cutoff_arg
    tumor_cutoff = tumor_cutoff_arg
    normal_prevalance_cutoff = normal_prevalance_cutoff_arg
    tumor_prevalance_cutoff = tumor_prevalance_cutoff_arg

    return adata


def get_all_normal_h5ad(uids,outdir,name):
    '''
    Get the normal tissue expression as a h5ad file, so that you can run BayesTS analysis (https://github.com/frankligy/BayesTS#interface-with-snaf)

    :param uids: list, all the uids you'd like to query
    :param outdir: string, the output directory
    :param name: string, the h5ad name

    Example::

        get_all_normal_h5ad(uids=uids,outdir='.',name='junction')
        
    '''
    uids = [u for u in uids if u in set(adata.obs_names)]
    adata_new = adata[uids,:].copy()
    print(adata_new)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    if not name.endswith('.h5ad'):
        name = '{}.h5ad'.format(name)
    adata_new.write(os.path.join(outdir, name))
    return adata_new


# Set SNAF_LEGACY_SIFTING=1 to force the original dense implementation (kept for
# numerical-equivalence testing). Default is the vectorized/sparse path, which lifts
# the ~100k-junction memory ceiling and is dramatically faster.
_USE_LEGACY_SIFTING = os.environ.get('SNAF_LEGACY_SIFTING', '0') == '1'


def multiple_crude_sifting(junction_count_matrix,add_control,dict_exonlist,outdir,filter_mode):
    if _USE_LEGACY_SIFTING:
        prevalance_fn, maxmin_fn = multiple_crude_sifting_prevalance, multiple_crude_sifting_maxmin
    else:
        prevalance_fn, maxmin_fn = multiple_crude_sifting_prevalance_fast, multiple_crude_sifting_maxmin_fast
    if filter_mode == 'prevalance':
        valid,invalid,cond_df = prevalance_fn(junction_count_matrix,add_control,dict_exonlist,outdir)
    elif filter_mode == 'maxmin':
        valid,invalid,cond_df = maxmin_fn(junction_count_matrix,add_control,dict_exonlist,outdir)
    return valid,invalid, cond_df

def multiple_crude_sifting_prevalance(junction_count_matrix,add_control=None,dict_exonlist=None,outdir='.'):

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    df_to_write = []
    df = pd.DataFrame(index=junction_count_matrix.index)
    prevalance_tumor = np.count_nonzero((junction_count_matrix > tumor_cutoff).values,axis=1) / junction_count_matrix.shape[1]
    df['prevalance_tumor'] = prevalance_tumor
    # consider gtex
    prevalance_normal = np.count_nonzero((adata_gtex.X > normal_cutoff).toarray(),axis=1) / adata_gtex.shape[1]
    prevalance_normal_dict = {j:v for j,v in zip(adata_gtex.obs_names,prevalance_normal)}
    df['prevalance_normal'] = df.index.map(prevalance_normal_dict).fillna(value=0)
    df['cond'] = (df['prevalance_tumor'] > tumor_prevalance_cutoff) & (df['prevalance_normal'] < normal_prevalance_cutoff)
    valid = df.loc[df['cond']].index.tolist()
    tmp = df.copy()
    df_to_write.append(tmp)
    print('reduce valid NeoJunction from {} to {} because they are present in GTEx'.format(df.shape[0],len(valid)))
    if dict_exonlist is not None:   # a valid junction can not be present in any ensembl documented transcript
        updated_valid = []
        for uid in tqdm(valid):
            ensg = uid.split(':')[0]
            exons = ':'.join(uid.split(':')[1:])
            if '_' in exons or 'U' in exons or 'ENSG' in exons or 'I' in exons:
                updated_valid.append(uid)
            else:
                exonlist = dict_exonlist[ensg]
                exonstring = '|'.join(exonlist)
                e1,e2 = exons.split('-')
                pattern1 = re.compile(r'^{}\|{}\|'.format(e1,e2))  # ^E1.1|E2.3|
                pattern2 = re.compile(r'\|{}\|{}$'.format(e1,e2))  # |E1.1|E2.3$
                pattern3 = re.compile(r'\|{}\|{}\|'.format(e1,e2)) # |E1.1|E2.3|
                if re.search(pattern3,exonstring) or re.search(pattern2,exonstring) or re.search(pattern1,exonstring):   # as long as match one pattern, should be eliminated
                    continue
                else:
                    updated_valid.append(uid)
        print('reduce valid Neojunction from {} to {} because they are present in Ensembl db'.format(len(valid),len(updated_valid)))
        valid = updated_valid
    # consider add_control
    if add_control is not None:
        for i,(id_,control) in enumerate(add_control.items()):
            n_previous_valid = len(valid)
            if isinstance(control,pd.DataFrame):
                prevalance_normal = np.count_nonzero((control > normal_cutoff).values,axis=1) / control.shape[1]
                prevalance_normal_dict = {j:v for j,v in zip(control.index, prevalance_normal)}
            elif isinstance(control,ad.AnnData):
                prevalance_normal = np.count_nonzero((control.X > normal_cutoff).toarray(),axis=1) / control.shape[1]
                prevalance_normal_dict = {j:v for j,v in zip(control.obs_names,prevalance_normal)}
            else:
                raise Exception('control must be either in dataframe or anndata format')
            df['prevalance_normal_add'] = df.index.map(prevalance_normal_dict).fillna(value=0)
            df['cond_add'] = (df['prevalance_tumor'] > tumor_prevalance_cutoff) & (df['prevalance_normal_add'] < normal_prevalance_cutoff)
            valid_add = df.loc[df['cond_add']].index.tolist()
            valid = list(set(valid).intersection(set(valid_add)))
            tmp = df.copy(); tmp.drop(columns=['prevalance_tumor','prevalance_normal','cond'],inplace=True); tmp.rename(columns=lambda x:x+'_{}'.format(id_),inplace=True)
            df_to_write.append(tmp)
            print('reduce valid Neojunction from {} to {} because they are present in added control {}'.format(n_previous_valid,len(valid),id_))
    invalid = list(set(junction_count_matrix.index).difference(set(valid)))
    # now, consider each entry
    t_min = tumor_cutoff
    valid_set = set(valid)
    cond_dict = {j:(True if j in valid_set else False) for j in junction_count_matrix.index}
    tmp = pd.DataFrame(index=junction_count_matrix.index,data={'placeholder':junction_count_matrix.index.map(cond_dict).values})
    first_half_cond_df = pd.concat([tmp]*junction_count_matrix.shape[1],axis=1)
    first_half_cond_df.columns = junction_count_matrix.columns
    cond_df = (first_half_cond_df) & (junction_count_matrix > t_min)
    # write the df
    df_to_write = pd.concat(df_to_write,axis=1)
    df_to_write.to_csv(os.path.join(outdir,'NeoJunction_statistics_prevalance.txt'),sep='\t')
    return valid,invalid,cond_df



def multiple_crude_sifting_maxmin(junction_count_matrix,add_control=None,dict_exonlist=None,outdir='.'):   # for JunctionCountMatrixQuery class, only consider gtex
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    df = pd.DataFrame(index=junction_count_matrix.index,data = {'max':junction_count_matrix.max(axis=1).values})
    df_to_write = []
    # consider gtex
    junction_to_mean = adata_gtex.obs.loc[adata_gtex.obs_names.isin(junction_count_matrix.index),'mean'].to_dict()
    df['mean'] = df.index.map(junction_to_mean).fillna(value=0)
    df['diff'] = df['max'] - df['mean']
    df['cond'] = (df['mean'] < n_max) & (df['diff'] > t_min)
    valid = df.loc[df['cond']].index.tolist()
    tmp = df.copy()
    df_to_write.append(tmp)
    print('reduce valid NeoJunction from {} to {} because they are present in GTEx'.format(df.shape[0],len(valid)))
    if dict_exonlist is not None:   # a valid junction can not be present in any ensembl documented transcript
        updated_valid = []
        for uid in tqdm(valid):
            ensg = uid.split(':')[0]
            exons = ':'.join(uid.split(':')[1:])
            if '_' in exons or 'U' in exons or 'ENSG' in exons or 'I' in exons:
                updated_valid.append(uid)
            else:
                exonlist = dict_exonlist[ensg]
                exonstring = '|'.join(exonlist)
                e1,e2 = exons.split('-')
                pattern1 = re.compile(r'^{}\|{}\|'.format(e1,e2))  # ^E1.1|E2.3|
                pattern2 = re.compile(r'\|{}\|{}$'.format(e1,e2))  # |E1.1|E2.3$
                pattern3 = re.compile(r'\|{}\|{}\|'.format(e1,e2)) # |E1.1|E2.3|
                if re.search(pattern3,exonstring) or re.search(pattern2,exonstring) or re.search(pattern1,exonstring):   # as long as match one pattern, should be eliminated
                    continue
                else:
                    updated_valid.append(uid)
        print('reduce valid Neojunction from {} to {} because they are present in Ensembl db'.format(len(valid),len(updated_valid)))
        valid = updated_valid
    # consider add_control
    mean_add_list = []
    if add_control is not None:
        for i,(id_,control) in enumerate(add_control.items()):
            n_previous_valid = len(valid)
            if isinstance(control,pd.DataFrame):
                junction_to_mean = control.mean(axis=1).to_dict()
            elif isinstance(control,ad.AnnData):
                junction_to_mean = control.to_df().mean(axis=1).to_dict()
            else:
                raise Exception('control must be either in dataframe or anndata format')
            df['mean_add'] = df.index.map(junction_to_mean).fillna(value=0)
            df['diff_add'] = df['max'] - df['mean_add']
            df['cond_add'] = (df['mean_add'] < n_max) & (df['diff_add'] > t_min)
            mean_add_list.append(df['mean_add'])
            valid_add = df.loc[df['cond_add']].index.tolist()
            valid = list(set(valid).intersection(set(valid_add)))
            tmp = df.copy(); tmp.drop(columns=['mean','diff','cond'],inplace=True); tmp.rename(columns=lambda x:x+'_{}'.format(id_),inplace=True)
            df_to_write.append(tmp)
            print('reduce valid Neojunction from {} to {} because they are present in added control {}'.format(n_previous_valid,len(valid),id_))
    invalid = list(set(junction_count_matrix.index).difference(set(valid)))
    # now, consider each entry
    gtex_df = pd.concat([df['mean']]*junction_count_matrix.shape[1],axis=1)
    gtex_df.columns = junction_count_matrix.columns
    diff_df_gtex = junction_count_matrix - gtex_df
    cond_df = (gtex_df < n_max) & (diff_df_gtex > t_min)
    if add_control is not None:
        for mean_add in mean_add_list:
            add_df = pd.concat([mean_add]*junction_count_matrix.shape[1],axis=1)    
            add_df.columns = junction_count_matrix.columns
            diff_df_add = junction_count_matrix - add_df
            cond_df = cond_df & (add_df < n_max) & (diff_df_add > t_min)
    # write the df
    df_to_write = pd.concat(df_to_write,axis=1)
    df_to_write.to_csv(os.path.join(outdir,'NeoJunction_statistics_maxmin.txt'),sep='\t')
    return valid,invalid,cond_df


# ---------------------------------------------------------------------------
# Vectorized / sparse re-implementations of the crude-sifting step.
#
# The original functions above tile a per-junction mean column into a full
# dense (n_junctions x n_samples) DataFrame with pd.concat([col]*n_samples) and
# then build several more full dense copies, and densify the entire sparse GTEx
# control matrix with .toarray(). That is the reason the workflow could not
# handle >~100k junctions. The _fast versions below produce byte-for-byte the
# same `valid`/`invalid`/`cond_df` using numpy broadcasting and scipy.sparse
# (no tiling, no densification of the control DB). Equivalence is verified in
# tests/test_sifting_equivalence.py.
# ---------------------------------------------------------------------------

def _build_gene_pair_sets(dict_exonlist):
    '''Precompute, per gene, the set of consecutive exon pairs (e1,e2) that occur
    in any documented transcript. Replaces the per-junction triple-regex scan
    (3 re.compile + 3 re.search per junction) in the legacy code with an O(1)
    set membership test. A junction ENSG:e1-e2 is "in Ensembl db" iff (e1,e2) is
    an adjacent pair in the pipe-joined exon string of some transcript, which is
    exactly what the ^e1|e2| / |e1|e2| / |e1|e2$ patterns matched.'''
    gene_pair_sets = {}
    for ensg, exonlist in dict_exonlist.items():
        pairs = set()
        for exonstring in exonlist:
            toks = exonstring.split('|')
            for a, b in zip(toks[:-1], toks[1:]):
                pairs.add((a, b))
        gene_pair_sets[ensg] = pairs
    return gene_pair_sets


def _filter_valid_not_in_db(valid, dict_exonlist):
    '''Keep only junctions NOT documented in any Ensembl transcript (same result
    as the legacy per-junction regex loop, but via precomputed pair sets).'''
    gene_pair_sets = _build_gene_pair_sets(dict_exonlist)
    updated_valid = []
    for uid in valid:
        ensg = uid.split(':')[0]
        exons = ':'.join(uid.split(':')[1:])
        if '_' in exons or 'U' in exons or 'ENSG' in exons or 'I' in exons:
            updated_valid.append(uid)
            continue
        e1, e2 = exons.split('-')
        pairs = gene_pair_sets.get(ensg)
        if pairs is not None and (e1, e2) in pairs:
            continue  # documented -> drop
        updated_valid.append(uid)
    return updated_valid


def _mean_vector_from_control(control, index):
    '''Per-junction mean across control samples, aligned to `index`, without
    densifying the whole control matrix. Works for DataFrame or AnnData.'''
    if isinstance(control, pd.DataFrame):
        s = control.mean(axis=1)
    elif isinstance(control, ad.AnnData):
        X = control.X
        if issparse(X):
            m = np.asarray(X.mean(axis=1)).ravel()
        else:
            m = np.asarray(X).mean(axis=1)
        s = pd.Series(m, index=control.obs_names)
    else:
        raise Exception('control must be either in dataframe or anndata format')
    return s.reindex(index).fillna(0.0).values.astype(np.float64)


def _prevalance_vector_from_control(control, index, cutoff):
    '''Per-junction fraction of control samples with count > cutoff, aligned to
    `index`, computed directly on the sparse matrix (no .toarray()).'''
    if isinstance(control, pd.DataFrame):
        frac = np.count_nonzero((control > cutoff).values, axis=1) / control.shape[1]
        s = pd.Series(frac, index=control.index)
    elif isinstance(control, ad.AnnData):
        # Summary-backed control: prevalence was precomputed (at snaf_summary_cutoff); use it
        # directly instead of scanning a matrix that isn't loaded.
        if control.uns.get('snaf_summary_backed') and 'normal_prevalence' in control.obs.columns:
            pc = control.uns.get('snaf_summary_cutoff')
            if pc is not None and float(pc) != float(cutoff):
                logger.warning('normal_cutoff=%s differs from precomputed summary cutoff=%s; '
                               'prevalence reflects the precompute cutoff.', cutoff, pc)
            return control.obs['normal_prevalence'].reindex(index).fillna(0.0).values.astype(np.float64)
        X = control.X
        n = control.shape[1]
        if issparse(X):
            # (X > cutoff) stays sparse for cutoff >= 0; count per row via sum
            frac = np.asarray((X > cutoff).sum(axis=1)).ravel() / n
        else:
            frac = np.count_nonzero(np.asarray(X) > cutoff, axis=1) / n
        s = pd.Series(frac, index=control.obs_names)
    else:
        raise Exception('control must be either in dataframe or anndata format')
    return s.reindex(index).fillna(0.0).values.astype(np.float64)


def _reindex_csr_rows(X, old_index, new_index):
    '''Return a CSR matrix whose rows follow `new_index`, filling missing rows with
    zeros. Pure sparse gather (permutation matmul) -- never densifies.'''
    X = csr_matrix(X)
    pos = pd.Index(old_index).get_indexer(pd.Index(new_index))   # -1 where missing
    present = pos >= 0
    rows = np.nonzero(present)[0]
    cols = pos[present]
    P = csr_matrix((np.ones(len(rows), dtype=X.dtype), (rows, cols)),
                   shape=(len(new_index), X.shape[0]))
    return P.dot(X).tocsr()


def _combine_two_sparse(left_X, left_rows, left_cols, right_X, right_rows, right_cols):
    '''Outer-join two sparse matrices on their row labels (union) and concatenate
    their columns -- the sparse equivalent of df_left.join(df_right,how="outer").
    Replaces adata.to_df()/control.to_df() dense materialization.'''
    left_idx = pd.Index(left_rows)
    right_idx = pd.Index(right_rows)
    right_only = right_idx.difference(left_idx, sort=False)
    new_rows = left_idx.append(right_only)
    L = _reindex_csr_rows(left_X, left_rows, new_rows)
    R = _reindex_csr_rows(right_X, right_rows, new_rows)
    X = sparse_hstack([L, R]).tocsr()
    return X, list(new_rows), list(left_cols) + list(right_cols)


def multiple_crude_sifting_maxmin_fast(junction_count_matrix, add_control=None, dict_exonlist=None, outdir='.'):
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    jc_index = junction_count_matrix.index
    jc_cols = junction_count_matrix.columns
    jc_vals = np.asarray(junction_count_matrix.values, dtype=np.float64)   # single dense copy of the input
    n_samples = jc_vals.shape[1]

    max_vec = jc_vals.max(axis=1)
    df = pd.DataFrame(index=jc_index, data={'max': max_vec})
    df_to_write = []

    # GTEx mean (already precomputed in adata_gtex.obs['mean']); align, missing -> 0
    mean_vec = adata_gtex.obs['mean'].reindex(jc_index).fillna(0.0).values.astype(np.float64)
    df['mean'] = mean_vec
    df['diff'] = max_vec - mean_vec
    df['cond'] = (mean_vec < n_max) & (df['diff'].values > t_min)
    valid = jc_index[df['cond'].values].tolist()
    df_to_write.append(df.copy())
    print('reduce valid NeoJunction from {} to {} because they are present in GTEx'.format(df.shape[0], len(valid)))

    if dict_exonlist is not None:
        updated_valid = _filter_valid_not_in_db(valid, dict_exonlist)
        print('reduce valid Neojunction from {} to {} because they are present in Ensembl db'.format(len(valid), len(updated_valid)))
        valid = updated_valid

    # per-cell condition matrix (single boolean array via broadcasting, no tiling)
    cond_arr = (mean_vec[:, None] < n_max) & ((jc_vals - mean_vec[:, None]) > t_min)

    if add_control is not None:
        for id_, control in add_control.items():
            n_previous_valid = len(valid)
            mean_add = _mean_vector_from_control(control, jc_index)
            diff_add = max_vec - mean_add
            cond_add = (mean_add < n_max) & (diff_add > t_min)
            valid_add = jc_index[cond_add].tolist()
            valid = list(set(valid).intersection(set(valid_add)))
            tmp = df.copy(); tmp.drop(columns=['mean', 'diff', 'cond'], inplace=True)
            tmp['mean_add'] = mean_add; tmp['diff_add'] = diff_add; tmp['cond_add'] = cond_add
            tmp.rename(columns=lambda x: x + '_{}'.format(id_), inplace=True)
            df_to_write.append(tmp)
            cond_arr &= (mean_add[:, None] < n_max) & ((jc_vals - mean_add[:, None]) > t_min)
            print('reduce valid Neojunction from {} to {} because they are present in added control {}'.format(n_previous_valid, len(valid), id_))

    invalid = list(set(jc_index).difference(set(valid)))
    cond_df = pd.DataFrame(cond_arr, index=jc_index, columns=jc_cols)
    df_to_write = pd.concat(df_to_write, axis=1)
    df_to_write.to_csv(os.path.join(outdir, 'NeoJunction_statistics_maxmin.txt'), sep='\t')
    return valid, invalid, cond_df


def multiple_crude_sifting_prevalance_fast(junction_count_matrix, add_control=None, dict_exonlist=None, outdir='.'):
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    jc_index = junction_count_matrix.index
    jc_cols = junction_count_matrix.columns
    jc_vals = np.asarray(junction_count_matrix.values, dtype=np.float64)

    df_to_write = []
    df = pd.DataFrame(index=jc_index)
    prevalance_tumor = np.count_nonzero(jc_vals > tumor_cutoff, axis=1) / jc_vals.shape[1]
    df['prevalance_tumor'] = prevalance_tumor
    # GTEx: count-over-cutoff directly on the sparse control matrix (no .toarray())
    prevalance_normal = _prevalance_vector_from_control(adata_gtex, jc_index, normal_cutoff)
    df['prevalance_normal'] = prevalance_normal
    df['cond'] = (prevalance_tumor > tumor_prevalance_cutoff) & (prevalance_normal < normal_prevalance_cutoff)
    valid = jc_index[df['cond'].values].tolist()
    df_to_write.append(df.copy())
    print('reduce valid NeoJunction from {} to {} because they are present in GTEx'.format(df.shape[0], len(valid)))

    if dict_exonlist is not None:
        updated_valid = _filter_valid_not_in_db(valid, dict_exonlist)
        print('reduce valid Neojunction from {} to {} because they are present in Ensembl db'.format(len(valid), len(updated_valid)))
        valid = updated_valid

    if add_control is not None:
        for id_, control in add_control.items():
            n_previous_valid = len(valid)
            prevalance_normal_add = _prevalance_vector_from_control(control, jc_index, normal_cutoff)
            cond_add = (prevalance_tumor > tumor_prevalance_cutoff) & (prevalance_normal_add < normal_prevalance_cutoff)
            valid_add = jc_index[cond_add].tolist()
            valid = list(set(valid).intersection(set(valid_add)))
            tmp = pd.DataFrame(index=jc_index)
            tmp['prevalance_normal_add'] = prevalance_normal_add
            tmp['cond_add'] = cond_add
            tmp.rename(columns=lambda x: x + '_{}'.format(id_), inplace=True)
            df_to_write.append(tmp)
            print('reduce valid Neojunction from {} to {} because they are present in added control {}'.format(n_previous_valid, len(valid), id_))

    invalid = list(set(jc_index).difference(set(valid)))
    # per-cell condition: valid-row gate AND per-sample tumor count > tumor_cutoff
    valid_set = set(valid)
    row_valid = np.array([j in valid_set for j in jc_index])
    cond_arr = row_valid[:, None] & (jc_vals > tumor_cutoff)
    cond_df = pd.DataFrame(cond_arr, index=jc_index, columns=jc_cols)
    df_to_write = pd.concat(df_to_write, axis=1)
    df_to_write.to_csv(os.path.join(outdir, 'NeoJunction_statistics_prevalance.txt'), sep='\t')
    return valid, invalid, cond_df


def crude_tumor_specificity(uid,count):    # for NeoJunction class, since we normally start from Jcmq with check_gtex=False, rarely being called.
    detail = ''
    if uid not in set(adata.obs_names):
        mean_value = 0
    else:
        mean_value = adata.obs.loc[uid,'mean']
    diff = count - mean_value
    if mean_value < n_max and diff >= t_min:
        identity = True
    else:
        identity = False
    return identity,mean_value


def mle_func(parameters,y):
    sigma = parameters
    ll = np.sum(stats.halfnorm.logpdf(y,0,sigma))
    neg_ll = -1 * ll
    return neg_ll

def split_df_to_chunks(df,cores=None):
    df_index = np.arange(df.shape[0])
    if cores is None:
        cores = mp.cpu_count()
    sub_indices = np.array_split(df_index,cores)
    sub_dfs = [df.iloc[sub_index,:] for sub_index in sub_indices]
    return sub_dfs


def split_array_to_chunks(array,cores=None):
    if not isinstance(array,list):
        raise Exception('split_array_to_chunks function works for list, not ndarray')
    array_index = np.arange(len(array))
    if cores is None:
        cores = mp.cpu_count()
    sub_indices = np.array_split(array_index,cores)
    sub_arrays = []
    for sub_index in sub_indices:
        item_in_group = []
        for i in sub_index:
            item_in_group.append(array[i])
        sub_arrays.append(item_in_group)
    return sub_arrays


def add_tumor_specificity_frequency_table(df,method='mean',remove_quote=True,cores=None,bayes_kwargs=None):
    '''
    add tumor specificty to each neoantigen-uid in the frequency table produced by SNAF T pipeline

    :param df: DataFrame, the frequency table produced by SNAF T pipeline
    :param method: string, either 'mean', or 'mle', or 'bayesian'
    :param remove quote: boolean, whether to remove the quotation or not, as one column in frequency table df is list, when loaded in memory using pandas, it will be added a quote, we can remove it
    :param cores: int, how many cpu cores to use for this computation, default None and use all the cpu the program detected
    :param bayes_kwargs: dict or None, extra options forwarded to BayesTS when method=='bayesian'
                         (e.g. {'mode':'XY','epoch':2000,'weights_dict':{...},'noise':3.0,'min_sample':10})

    :return new_df: a dataframe with one added column containing tumor specificity score

    Example::

        snaf.add_tumor_specificity_frequency_table(df,'mle',remove_quote=True)

    '''
    from ast import literal_eval
    if remove_quote:
        df['samples'] = [literal_eval(item) for item in df['samples']]

    all_unique_junctions = list(set([item.split(',')[1] for item in df.index]))

    if method != 'bayesian':
        # Fully vectorized batch score (byte-identical to the per-uid loop, per
        # tumor_specificity_batch's contract): mean = obs['mean'] lookup, mle = sparse
        # closed-form pass -- NO per-junction AnnData slicing / X.toarray(). This per-uid
        # loop was the dominant cost of generate_results (called 2 methods x 3 stages).
        all_score_dict = tumor_specificity_batch(all_unique_junctions, methods=(method,))[method].to_dict()
        col = [all_score_dict[item.split(',')[1]] for item in df.index]
    else:
        # Direct BayesTS integration (joint model over all queried junctions),
        # replacing the legacy per-junction pymc3 path. Cross-platform (CPU).
        try:
            from . import bayests
        except (ImportError, ValueError):   # allow standalone loading in tests
            import importlib.util as _ilu
            _spec = _ilu.spec_from_file_location('snaf_bayests', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bayests.py'))
            bayests = _ilu.module_from_spec(_spec); _spec.loader.exec_module(bayests)
        bk = {'mode': 'XY', 'epoch': 2000}
        if bayes_kwargs:
            bk.update(bayes_kwargs)
        res = bayests.compute_bayests_sigma(adata, uids=all_unique_junctions, **bk)
        score_dict = res['mean_sigma'].to_dict()
        col = [score_dict.get(item.split(',')[1], np.nan) for item in df.index]

    new_df = df.copy()
    new_df['tumor_specificity_{}'.format(method)] = col
    return new_df

def add_tumor_specificity_frequency_table_atomic_func(sub_array,method):
    uid_list = sub_array
    score_dict = {uid:tumor_specificity(uid,method) for uid in tqdm(uid_list,total=len(uid_list))}
    return score_dict


    

def tumor_specificity(uid,method,return_df=False):
    try:
        info = adata[[uid],:]
    except:
        logger.debug('%s not detected in gtex, impute as zero', uid)
        info = ad.AnnData(X=csr_matrix(np.full((1,adata.shape[1]),0)),obs=pd.DataFrame(data={'mean':[0]},index=[uid]),var=adata.var)  # weired , anndata 0.7.6 can not modify the X in place? anndata 0.7.2 can do that in scTriangulate
    df = pd.DataFrame(data={'value':info.X.toarray().squeeze(),'tissue':info.var['tissue'].values},index=info.var_names)
    if method == 'mean':
        try:
            sigma = adata.obs.loc[uid,'mean']
        except KeyError:
            sigma = 0
        if return_df:
            return sigma,df
        else:
            return sigma
    elif method == 'mle':
        scale_factor_dict = adata.var['total_count'].to_dict()
        df['value_cpm'] = df['value'].values / df.index.map(scale_factor_dict).values
        y = df['value_cpm'].values
        # Closed-form MLE of the half-normal scale is sqrt(mean(y^2)); the original
        # code used minimize_scalar over the bound (0,1), so clip identically. This
        # is exact (not an approximation) and removes a per-junction optimizer call.
        if len(y) > 0:
            sigma = float(min(np.sqrt(np.mean(np.square(y))), 1.0))
        else:
            sigma = 0.0
        if return_df:
            return sigma,df
        else:
            return sigma
    elif method == 'bayesian':
        # BayesTS is a JOINT model over all queried junctions -- scoring one uid in
        # isolation is not meaningful. Use the batch entry point instead, which runs
        # the real BayesTS (torch/pyro, cross-platform) once over the whole set:
        #     snaf.add_tumor_specificity_frequency_table(df, method='bayesian')
        # or, directly:
        #     from altanalyze3.components.snaf import bayests
        #     bayests.compute_bayests_sigma(adata, uids=[...], mode='XY')
        # (The legacy per-junction pymc3 model that lived here required pymc3/theano,
        # which does not install cleanly on macOS/Windows; it has been retired.)
        raise NotImplementedError(
            "Per-junction bayesian scoring is retired. BayesTS is a joint model; "
            "call add_tumor_specificity_frequency_table(df, method='bayesian') or "
            "bayests.compute_bayests_sigma(adata, uids=...) for the batch score.")


def tumor_specificity_batch(uids, methods=('mean', 'mle')):
    '''Vectorized equivalent of ``[tumor_specificity(u, m) for u in uids]`` for each
    method in ``methods``. Returns a DataFrame indexed by ``uids`` with one column per
    method ('mean' -> mean_gtex_count, 'mle' -> tumor_specificity_mle).

    Results are IDENTICAL to the per-uid path (this only removes the Python-level loop):
      - mean : ``adata.obs['mean']`` (missing junction -> 0)
      - mle  : ``min(sqrt(mean(cpm^2)), 1.0)`` where ``cpm = count / var['total_count']``,
               the same closed-form half-normal MLE, averaged over ALL control samples.

    Memory-safe: the MLE is computed with sparse row ops (no densification of the
    junctions x samples control matrix).'''
    uids = list(uids)
    present_set = set(adata.obs_names)
    out = pd.DataFrame(index=uids)

    if 'mean' in methods:
        out['mean'] = adata.obs['mean'].reindex(uids).fillna(0.0).values.astype(np.float64)

    if 'mle' in methods:
        # Summary-backed control: mle was precomputed per junction; look it up (no matrix).
        if adata.uns.get('snaf_summary_backed') and 'mle' in adata.obs.columns:
            out['mle'] = adata.obs['mle'].reindex(uids).fillna(0.0).values.astype(np.float64)
            return out
        sigma = np.zeros(len(uids), dtype=np.float64)
        present_uids = [u for u in uids if u in present_set]
        if present_uids:
            inv_scale = (1.0 / adata.var['total_count'].values.astype(np.float64))  # per-sample CPM factor
            pos = {u: i for i, u in enumerate(uids)}
            n_cols = adata.shape[1]
            # Process in chunks so the (chunk x samples) sparse intermediates stay small
            # -- keeps peak memory bounded on multi-hundred-k neojunction cohorts while
            # remaining fully vectorized. Row-wise math is unchanged, so results are identical.
            CHUNK = 50000
            for start in range(0, len(present_uids), CHUNK):
                batch = present_uids[start:start + CHUNK]
                X = adata[batch, :].X
                if not issparse(X):
                    X = csr_matrix(np.asarray(X))
                Xs = X.multiply(inv_scale[None, :])                        # cpm, sparse
                sq_sum = np.asarray(Xs.multiply(Xs).sum(axis=1)).ravel()   # sum of cpm^2 per junction
                s = np.minimum(np.sqrt(sq_sum / n_cols), 1.0)             # mean over ALL samples, clipped
                for u, val in zip(batch, s):
                    sigma[pos[u]] = val
        out['mle'] = sigma

    return out
            
















