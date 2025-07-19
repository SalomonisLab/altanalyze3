import os, sys
import time
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from anndata import AnnData, concat, read_h5ad
from joblib import Parallel, delayed
import multiprocessing

def extract_features_above_threshold(h5ad_files, min_reads=200):
    """Sum total reads across features and retain those above threshold."""
    feature_sums = {}
    for file in h5ad_files:
        adata = read_h5ad(file)
        counts = np.asarray(adata.X.sum(axis=0)).ravel()
        for feat, count in zip(adata.var_names, counts):
            feature_sums[feat] = feature_sums.get(feat, 0) + count
    retained = [feat for feat, total in feature_sums.items() if total >= min_reads]
    return sorted(retained)


def align_h5ad_to_union(h5ad_file, all_features, feature_to_idx):
    adata = read_h5ad(h5ad_file)
    sample_id = os.path.basename(h5ad_file).replace('.h5ad', '')
    adata.obs_names = [f"{sample_id}.{name}" for name in adata.obs_names]

    # Check for duplicates early
    if not adata.var_names.is_unique:
        duplicated = adata.var_names[adata.var_names.duplicated(keep=False)]
        first_duplicate = duplicated.unique()[0]

        print(f"\n🚨 Duplicate feature detected in file: {h5ad_file}")
        print(f"Feature name: {first_duplicate}")

        idxs = np.where(adata.var_names == first_duplicate)[0]

        for idx_num, idx in enumerate(idxs):
            feature_sum = adata.X[:, idx].sum()
            print(f"  Occurrence {idx_num+1}: Sum across cells = {feature_sum}")

        print(f"\nTotal number of duplicated feature names in file: {duplicated.nunique()}")
        sys.exit(1)

    old_features = adata.var_names
    old_to_new_idx = {i: feature_to_idx[feat] for i, feat in enumerate(old_features) if feat in feature_to_idx}

    X = adata.X.tocoo()
    rows, cols, data = [], [], []
    for i, j, v in zip(X.row, X.col, X.data):
        if j in old_to_new_idx:
            rows.append(i)
            cols.append(old_to_new_idx[j])
            data.append(v)
    new_X = csr_matrix((data, (rows, cols)), shape=(adata.n_obs, len(all_features)))

    new_var = pd.DataFrame(index=all_features)
    for col in adata.var.columns:
        dtype = float if pd.api.types.is_numeric_dtype(adata.var[col]) else object
        new_var[col] = pd.Series(index=all_features, dtype=dtype)
        shared_features = adata.var.index.intersection(new_var.index)
        new_var.loc[shared_features, col] = adata.var.loc[shared_features, col]

    adata.obs = adata.obs.copy()
    adata.obs['sample'] = sample_id
    return AnnData(X=new_X, obs=adata.obs, var=new_var, uns=adata.uns.copy())


def align_h5ad_to_union2(file, all_features, feature_to_idx):
    """Align a single .h5ad file to the global feature list."""
    adata = read_h5ad(file)
    sample_id = os.path.basename(file).replace('.h5ad', '')
    adata.obs_names = [f"{sample_id}.{name}" for name in adata.obs_names]

    old_features = adata.var_names
    X = adata.X.tocoo()

    rows, cols, data = [], [], []
    for i, j, v in zip(X.row, X.col, X.data):
        if old_features[j] in feature_to_idx:
            rows.append(i)
            cols.append(feature_to_idx[old_features[j]])
            data.append(v)
    new_X = csr_matrix((data, (rows, cols)), shape=(adata.n_obs, len(all_features)))

    new_var = pd.DataFrame(index=all_features)
    if adata.var.shape[1] > 0:
        for col in adata.var.columns:
            dtype = float if pd.api.types.is_numeric_dtype(adata.var[col]) else object
            new_var[col] = pd.Series(index=all_features, dtype=dtype)
            shared = adata.var.index.intersection(new_var.index)
            new_var.loc[shared, col] = adata.var.loc[shared, col]

    new_obs = adata.obs.copy()
    return AnnData(X=new_X, obs=new_obs, var=new_var, uns=adata.uns.copy())

def safe_parallel(func, iterable, n_jobs=4):
    """Safe parallel execution using multiprocessing spawn context."""
    ctx = multiprocessing.get_context('spawn')
    return Parallel(n_jobs=n_jobs, backend='loky')(delayed(func)(*args) for args in iterable)

def combine_h5ad_files_parallel(sample_h5ads_dict, output_file='combined.h5ad', n_jobs=4, min_total_reads=200):
    """Combine many h5ad files into a single AnnData, aligning features in parallel."""
    start_time = time.time()
    h5ad_files = sorted(set(sample_h5ads_dict.values()))

    print(f"Identifying union of features across {len(h5ad_files)} files...")
    all_features = extract_features_above_threshold(h5ad_files, min_total_reads)
    feature_to_idx = {feat: idx for idx, feat in enumerate(all_features)}

    print(f"Aligning .h5ad files in parallel with {n_jobs} workers...")
    aligned_adatas = safe_parallel(align_h5ad_to_union, [(file, all_features, feature_to_idx) for file in h5ad_files], n_jobs=n_jobs)

    print("Concatenating aligned AnnData objects...")
    combined_adata = concat(aligned_adatas, axis=0, join='outer')

    print(f"Writing combined file to: {output_file}")
    combined_adata.write_h5ad(output_file, compression='gzip')

    elapsed = time.time() - start_time
    print(f"Combination complete. Elapsed time: {elapsed:.2f} seconds.")
    return combined_adata


if __name__ == '__main__':

    folder = "/Users/saljh8/Dropbox/Revio/BoneMarrow/test"
    files = os.listdir(folder)
    h5ad_files = {}
    for file in files:
        if file.endswith(".h5ad"):
            h5ad_files[file[:-5]] = os.path.join(folder, file)

    combine_h5ad_files_parallel(
        h5ad_files,
        output_file="combined_samples.h5ad",
        n_jobs=1,
        min_total_reads=200
    )