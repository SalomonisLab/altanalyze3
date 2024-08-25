import pandas as pd
import numpy as np
from scipy.stats import t


def unique_marker_finder(input_df, groups):

    # convert nmf_clusters to list object
    nmf_clusters_list = groups['cluster'].tolist()

    # ensure fold matrix has the same columns as rows/index of nmf_clusters
    fold_matrix = input_df.loc[groups.index, :]

    corr_df, p_val_df = marker_finder(input_df=fold_matrix, groups=nmf_clusters_list)

    # find the top n markers for each cluster
    top_cluster = corr_df.idxmax(axis=1)  # top cluster for each marker
    top_r = corr_df.max(axis=1)  # top r value for each marker

    # Assign each marker to its top scoring cluster
    markers_df = pd.DataFrame({"marker": top_cluster.index.values,
                               "top_cluster": top_cluster.values,
                               "pearson_r": top_r.values,
                               "p_value": [p_val_df.loc[i, j] for i, j in zip(top_cluster.index.values, top_cluster.values)]})

    markers_df = markers_df.sort_values(by=["top_cluster", "pearson_r"], ascending=[True, False])

    return markers_df


def create_final_marker_heatmap(input_df, markers_df, groups_df):

    # Filter the input_df rows with the index of markers_df
    final_marker_heatmap = input_df.loc[:, markers_df['marker']]

    final_marker_heatmap = final_marker_heatmap.transpose()
    final_marker_heatmap = scale_final_marker_heatmap(final_marker_heatmap)

    # Create a new row 'column_clusters-flat' based on the 'cluster' values from groups_df
    column_clusters = groups_df.loc[final_marker_heatmap.columns].transpose()
    final_marker_heatmap = pd.concat([column_clusters, final_marker_heatmap])
    sorted_columns = final_marker_heatmap.iloc[0].sort_values().index
    final_marker_heatmap = final_marker_heatmap[sorted_columns]

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


def marker_finder_wrapper(input_df, groups, top_n=60, rho_threshold=0.2, marker_finder_rho=0.3):

    # get unique markers
    markers_df_og = unique_marker_finder(input_df, groups)
    markers_df_og.index = markers_df_og['marker']
    markers_df = markers_df_og.copy()

    # count number of markers above rho_threshold per cluster
    markers_df = markers_df[markers_df['pearson_r'] >= rho_threshold]

    # clusters with fewer than two markers above the supplied rho threshold are excluded
    markers_df = markers_df.groupby('top_cluster').filter(lambda x: len(x) >= 3)

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
    t_df = r_df * np.sqrt(degrees_f) / np.sqrt(1 - (r_df ** 2))
    p_df = t_df.map(lambda x: t.sf(abs(x), df=degrees_f) * 2)
    return r_df, p_df


# Computes pearson correlation coefficient for pairwise columns to columns
# of input DataFrames
def pearson_corr_df_to_df(df1, df2):
    norm1 = df1 - df1.mean(axis=0)
    norm2 = df2 - df2.mean(axis=0)
    sqsum1 = (norm1**2).sum(axis=0)
    sqsum2 = (norm2**2).sum(axis=0)
    return (norm1.T @ norm2) / np.sqrt(sqsum1.apply(lambda x: x*sqsum2))
