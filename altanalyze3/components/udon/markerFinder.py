"""MarkerFinder CORRELATION CORE (Pearson of each feature to an idealized one-hot group
vector). This module does NOT render the AltAnalyze MarkerFinder heatmap.

For the full MarkerFinder that WRITES THE HEATMAP PDF (markers + heatmap TSV + figure), use
``components.cellHarmony.markerFinder`` -> ``find_markers_from_adata`` (AnnData in) or
``run_marker_finder_on_matrix(matrix_path, groups_path, output_dir)`` (matrix + groups file).
``render_heatmap.py`` here is only a lightweight raster imshow, not the canonical heatmap.
"""
import pandas as pd
import numpy as np
from scipy.stats import t


def unique_marker_finder(input_df, groups):

    # convert nmf_clusters to list object
    nmf_clusters_list = groups['cluster'].tolist()

    # ensure fold matrix has the same columns as rows/index of nmf_clusters.
    # Skip the reindex when the order already matches -- .loc copies the whole
    # cells x genes frame (2.17 GB at 30,000 x 18,100 float32) even when it is a no-op.
    if len(groups.index) == len(input_df.index) and input_df.index.equals(pd.Index(groups.index)):
        fold_matrix = input_df
    else:
        fold_matrix = input_df.loc[groups.index, :]

    corr_df, p_val_df = marker_finder(input_df=fold_matrix, groups=nmf_clusters_list)

    # find the top n markers for each cluster
    corr_arr = corr_df.to_numpy()
    best_col = np.argmax(corr_arr, axis=1)
    rows = np.arange(corr_arr.shape[0])
    top_r_vals = corr_arr[rows, best_col]
    top_cluster_vals = np.asarray(corr_df.columns)[best_col]

    # Gather the matching p-value by position. The previous form was a Python loop of
    # per-element .loc lookups, one per feature -- with the full gene pool that is
    # ~18,100 label-based lookups per MarkerFinder pass, twice per run.
    p_vals = p_val_df.to_numpy()[rows, best_col]

    # Assign each marker to its top scoring cluster
    markers_df = pd.DataFrame({"marker": corr_df.index.values,
                               "top_cluster": top_cluster_vals,
                               "pearson_r": top_r_vals,
                               "p_value": p_vals})

    # NUMERIC (natural) cluster order: P0,P1,...,P10 not lexicographic P0,P1,P10,...
    import re as _re
    _num = lambda x: int(_re.search(r"\d+", str(x)).group()) if _re.search(r"\d+", str(x)) else 0
    markers_df = (markers_df.assign(_k=markers_df["top_cluster"].map(_num))
                  .sort_values(by=["_k", "pearson_r"], ascending=[True, False]).drop(columns="_k"))

    return markers_df


def create_final_marker_heatmap(input_df, markers_df, groups_df):

    # Filter the input_df rows with the index of markers_df
    final_marker_heatmap = input_df.loc[:, markers_df['marker']]

    final_marker_heatmap = final_marker_heatmap.transpose()
    final_marker_heatmap = scale_final_marker_heatmap(final_marker_heatmap)

    # Create a new row 'column_clusters-flat' based on the 'cluster' values from groups_df
    column_clusters = groups_df.loc[final_marker_heatmap.columns].transpose()
    final_marker_heatmap = pd.concat([column_clusters, final_marker_heatmap])
    import re as _re, numpy as _np
    _num = lambda x: int(_re.search(r"\d+", str(x)).group()) if _re.search(r"\d+", str(x)) else 0
    _cc = final_marker_heatmap.iloc[0]                        # cluster per pseudobulk column
    sorted_columns = _cc.iloc[_np.argsort([_num(v) for v in _cc.values], kind="stable")].index
    final_marker_heatmap = final_marker_heatmap[sorted_columns]   # numeric cluster order

    # Create a new column 'row_clusters-flat' based on the 'top_cluster' values from markers_df
    row_clusters = markers_df.loc[final_marker_heatmap.index[1:], 'top_cluster']
    row_clusters = pd.concat([pd.Series('column_clusters-flat'), row_clusters])
    final_marker_heatmap.insert(0, 'row_clusters-flat', row_clusters)

    final_marker_heatmap.index.values[0] = 'column_clusters-flat'

    return final_marker_heatmap


def scale_final_marker_heatmap(final_marker_heatmap):

    # Drop the first row and column to get the matrix
    mat = final_marker_heatmap.values
    row_names = final_marker_heatmap.index
    col_names = final_marker_heatmap.columns

    # Compute row medians and perform row median centering
    row_medians = np.median(mat, axis=1)
    mat_n = mat - row_medians[:, None]

    # Create a dataframe for seaborn heatmap
    final_marker_heatmap_scaled = pd.DataFrame(mat_n, index=row_names, columns=col_names)

    return final_marker_heatmap_scaled


def marker_finder_wrapper(input_df, groups, top_n=60, rho_threshold=0.2, marker_finder_rho=0.3,
                          min_markers_per_cluster=3):

    # get unique markers
    markers_df_og = unique_marker_finder(input_df, groups)
    markers_df_og.index = markers_df_og['marker']
    markers_df = markers_df_og.copy()

    # count number of markers above rho_threshold per cluster
    markers_df = markers_df[markers_df['pearson_r'] >= rho_threshold]

    # clusters with too few markers above the supplied rho threshold are excluded.
    # Legacy RNA defaults used 3 here. Small feature panels such as ADT need a
    # lower configurable floor, while still requiring definitive markers.
    markers_df = markers_df.groupby('top_cluster').filter(lambda x: len(x) >= int(min_markers_per_cluster))

    # get Top n correlated marker for each cluster
    markers_df = markers_df.groupby('top_cluster').head(top_n)

    # keep only those markers whose pearson r is above rho_threshold
    markers_df = markers_df[markers_df['pearson_r'] >= marker_finder_rho]

    # create an output similar to FinalMarkerHeatmap file
    final_marker_heatmap = create_final_marker_heatmap(input_df, markers_df=markers_df, groups_df=groups)

    return markers_df_og, markers_df, final_marker_heatmap


def marker_finder(input_df, groups):
    """
    Function to find pearson correlation coefficient values and p-values for
    the given data and groups for groups to test. The function will perform a
    Pearson correlation of the input_df feature values to an "idealized"
    group specific expression vector, where each observation in a given group
    is set to a value of 1, and the observations in other groups are set to 0.

    Arguments
    ---------
    input_df : pandas.DataFrame
        Data Frame with observations as index and features as columns (Required)

    groups : list[str]
        List-like of specified groups corresponding to observations from the
        input_df. The order of groups should match the order in input_df.index
        (Required)

    Returns
    -------
    pandas.DataFrame, pandas.DataFrame
        The first item in the tuple is a pandas.DataFrame containing the pearson
        correlation coefficient values for each marker to the idealized vector
        for each cluster.
        The second item is also a pandas.DataFrame, but contains the p-values
        for each comparison.
    """
    ideal_vectors = pd.get_dummies(groups)
    ideal_vectors.index = input_df.index.values
    degrees_f = input_df.shape[0] - 2
    r_df = pearson_corr_df_to_df(input_df, ideal_vectors)
    r_df = r_df.dropna()
    # r == +/-1 makes 1 - r^2 zero and the t statistic infinite. Clip strictly inside
    # (-1, 1) so the divide is defined; sf of a huge |t| is 0 either way, so the p-value
    # is unchanged. This removes the RuntimeWarning storm rather than hiding it.
    r_arr = r_df.to_numpy(dtype=np.float64)
    r_safe = np.clip(r_arr, -1.0 + 1e-12, 1.0 - 1e-12)
    t_arr = r_safe * np.sqrt(degrees_f) / np.sqrt(1.0 - r_safe ** 2)
    p_vals = t.sf(np.abs(t_arr), df=degrees_f) * 2
    p_df = pd.DataFrame(p_vals, index=r_df.index, columns=r_df.columns)
    return r_df, p_df


# Computes pearson correlation coefficient for pairwise columns to columns
# of input DataFrames.
#
# Rewritten to run in float32 numpy rather than pandas. The pandas form built two
# full-size centred copies of df1 (cells x genes) and then a genes x clusters
# denominator via Series.apply, a Python-level loop. With the ICGS2-faithful gene
# pool df1 is 30,000 x 18,100, so each copy is 2.17 GB. Centring in place on a single
# float32 array removes one copy and halves the other. The arithmetic is identical.
def pearson_corr_df_to_df(df1, df2):
    cols1 = df1.columns
    cols2 = df2.columns
    a = np.asarray(df1.to_numpy(dtype=np.float32, copy=True))
    b = np.asarray(df2.to_numpy(dtype=np.float32, copy=True))
    a -= a.mean(axis=0, dtype=np.float32)
    b -= b.mean(axis=0, dtype=np.float32)
    sq1 = np.einsum("ij,ij->j", a, a)
    sq2 = np.einsum("ij,ij->j", b, b)
    num = a.T @ b
    denom = np.sqrt(np.outer(sq1, sq2))
    with np.errstate(divide="ignore", invalid="ignore"):
        r = num / denom
    # A zero-variance feature or group yields denom == 0 -> undefined correlation.
    # Leave it NaN so marker_finder's dropna removes it, rather than letting a filled
    # value rank as a real marker.
    r[~np.isfinite(r)] = np.nan
    return pd.DataFrame(r, index=cols1, columns=cols2)
