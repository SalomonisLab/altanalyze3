import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import zscore
from matplotlib.colors import ListedColormap
from matplotlib import cm
from matplotlib.colors import to_hex

# Define a colormap to use for marker expression heatmap
N = 256
vals = np.ones((N*2, 4))
vals[:N, 0] = np.linspace(15/256, 0, N)
vals[:N, 1] = np.linspace(255/256, 0, N)
vals[:N, 2] = np.linspace(255/256, 0, N)
vals[N:, 0] = np.linspace(0, 255/256, N)
vals[N:, 1] = np.linspace(0, 243/256, N)
vals[N:, 2] = np.linspace(0, 15/256, N)
blue_black_yellow_cmap = ListedColormap(vals)


def assign_rainbow_colors_to_groups(groups):
    """ Creates a dictionary of cluster names to hexadecimal color strings

    This function takes a list of groups and assigns each unique item in the
    groups a color (using the matplotlib.cm.rainbow color-map) as a hexadecimal
    string value. This is useful for storing a single color scheme for clusters
    to be used with downstream visualizations.

    Arguments
    ---------
    groups : numpy.Array[str]
        List of cluster names. (Required)

    Returns
    -------
    dict {str:str}
        Dictionary of cluster-names to assigned colors (hexadecimal value)

    """
    unique_groups = np.unique(groups)
    groups_to_num = pd.Series(list(range(len(unique_groups))),
                              index=unique_groups)
    n = len(unique_groups)
    groups_to_color = pd.Series([to_hex(cm.rainbow(item / n)) for item in
                                 groups_to_num.values], index=groups_to_num.index.values).to_dict()
    return groups_to_color


def plot_markers_df(marker_heatmap, markers_df, clusters, path_to_save_figure):
    """ Plots a heatmap of the MarkerFinder results
    Arguments
    ---------
    marker_heatmap : pandas.DataFrame
        Data to plot. The rows must intersect with features in the \
        ordered_markers_df and the columns must intersect with the cells in \
        clusters.
    markers_df : pandas.DataFrame
        The markers_df where "marker" column are the genes and the "top_cluster" is the cluster to which the gene belongs to  (Required)
    clusters : pandas.DataFrame
        The cluster assignments of the samples, index of this data frame are the sample names (Required)
    path_to_save_figure : str
        The path to save the figure. (Required)

    Returns
    -------
    None
        Saves the heatmap to the file specified by path_to_save_figure

    """
    try:
        # get Dictionary of cluster-names to assigned colors (hexadecimal value)
        groups_to_colors = assign_rainbow_colors_to_groups(groups=np.array(clusters["cluster"]))

        # Filter the input data matrix to the markers and cells of interest
        input_df = marker_heatmap.loc[markers_df["marker"].values, clusters.index]
        input_df = input_df.astype(float)

        # Map the order of the groups_to_colors to make a cmap
        groups_to_order = pd.Series(list(range(len(groups_to_colors))),
                                    index=groups_to_colors.keys())

        # Build top heatmap to label clusters
        cell_labels_df = pd.DataFrame({"cluster": [groups_to_order[item] for item in clusters["cluster"].values]},
                                      index=clusters.index).T

        # Build row annotation dataframe
        row_labels_df = pd.DataFrame({"cluster": [groups_to_order[item] for item in markers_df["top_cluster"].values]},
                                     index=markers_df["marker"].values)

        # Get top markers for each cluster
        top_markers = markers_df.groupby('top_cluster')['marker'].first()

        # Build df to add cluster label ticks
        label_to_position = pd.pivot_table(pd.DataFrame({
            "cluster": clusters["cluster"].values,
            "position": list(range(clusters.shape[0]))}),
            index="cluster",
            values="position",
            aggfunc=np.mean)["position"].sort_values()

        # Build df to add marker label ticks
        label_to_position2 = pd.pivot_table(pd.DataFrame({
            "cluster": markers_df["top_cluster"].values,
            "position": list(range(markers_df.shape[0]))}),
            index="cluster",
            values="position",
            aggfunc=np.mean)["position"].sort_values()

        # Z-score normalize the expression matrix
        z_df = input_df.apply(lambda x: pd.Series(zscore(x.values),
                                                  index=x.index.values), axis=1)

        # Build the heatmap
        plt.close("all")
        fig = plt.figure(constrained_layout=True, figsize=(9, 5))
        ax = fig.add_gridspec(2, 2, width_ratios=(1, 20), height_ratios=(1, 20),
                      wspace=0.005, hspace=0.01)  # Adjusted to create space for row annotation and top markers
        ax2 = fig.add_subplot(ax[1, 1])
        ax1 = fig.add_subplot(ax[0, 1])
        ax3 = fig.add_subplot(ax[1, 0])  # Subplot for row annotation

        heat1 = sns.heatmap(cell_labels_df,
                            yticklabels=False,
                            xticklabels=False,
                            cmap=sns.color_palette(groups_to_colors.values()),
                            cbar=False,
                            ax=ax1)
        ax1.set_xticks(label_to_position.values)
        ax1.set_xticklabels(label_to_position.index.values, rotation=45, ha="left")
        ax1.xaxis.tick_top()

        heat2 = sns.heatmap(z_df,
                            vmin=-3,
                            vmax=3,
                            cmap=blue_black_yellow_cmap,  # Assuming you have blue_black_yellow_cmap defined
                            xticklabels=False,
                            yticklabels=False,  # Disable yticklabels here
                            cbar=True,
                            cbar_kws={"shrink": 0.5},
                            ax=ax2)
        ax2.collections[0].colorbar.set_label("Z-Score Normalized Expression")

        heat3 = sns.heatmap(row_labels_df,
                            xticklabels=False,
                            yticklabels=False,  # Disable yticklabels here as well
                            cmap=sns.color_palette(groups_to_colors.values()),
                            cbar=False,
                            ax=ax3)
        ax3.set_yticks(label_to_position2.values)
        ax3.set_yticklabels(label_to_position.index.values)
        # ax3.yaxis.tick_left()
        fig.suptitle('Marker Finder Heatmap', fontsize=16)
        # fig.tight_layout()
        plt.savefig(path_to_save_figure, dpi=10, bbox_inches='tight', format="pdf", pad_inches=0.9)

    except Exception as e:
        print(str(e))
        print("Warning! Failed to run plot_markers_df. See above Exception.")

# Usage example (this part should be replaced with actual data in your script)
# marker_heatmap = pd.DataFrame(...)
# markers_df = pd.DataFrame(...)
# clusters = pd.DataFrame(...)
# groups_to_colors = {'cluster1': '#1f77b4', 'cluster2': '#ff7f0e', ...}
# path_to_save_figure = 'heatmap.png'
# plot_markers_df(marker_heatmap, markers_df, clusters, groups_to_colors, path_to_save_figure)


def plot_metadata_heatmap(udon_clusters, udon_metadata, metadata_col, path_to_save_figure):
    try:
        # get Dictionary of cluster-names to assigned colors (hexadecimal value)
        groups_to_colors = assign_rainbow_colors_to_groups(groups=np.array(udon_clusters["cluster"]))

        # Map the order of the groups_to_colors to make a cmap
        groups_to_order = pd.Series(list(range(len(groups_to_colors))),
                                    index=groups_to_colors.keys())

        # Build top heatmap to label clusters
        cell_labels_df = pd.DataFrame({"cluster": [groups_to_order[item] for item in udon_clusters["cluster"].values]},
                                      index=udon_clusters.index).T


        # Build df to add cluster label ticks
        label_to_position = pd.pivot_table(pd.DataFrame({
            "cluster": udon_clusters["cluster"].values,
            "position": list(range(udon_clusters.shape[0]))}),
            index="cluster",
            values="position",
            aggfunc=np.mean)["position"].sort_values()

        # Get unique clusters and samples
        unique_clusters = udon_metadata[metadata_col].unique()
        samples = udon_metadata.index

        # Initialize the binary matrix with zeros
        binary_matrix = pd.DataFrame(0, index=unique_clusters, columns=samples)

        # Populate the binary matrix
        for sample in samples:
            cluster = udon_metadata.loc[sample, metadata_col]
            binary_matrix.loc[cluster, sample] = 1

        binary_matrix = binary_matrix.sort_index()
        print("sorted index")

        # Build df to add cluster label ticks
        label_to_position2 = pd.pivot_table(pd.DataFrame({
            "cluster": binary_matrix.index,
            "position": list(range(binary_matrix.shape[0]))}),
            index="cluster",
            values="position",
            aggfunc=np.mean)["position"].sort_values()

        # Build the heatmap
        plt.close("all")
        fig = plt.figure(constrained_layout=True, figsize=(9, 6))
        ax = fig.add_gridspec(2, 1, height_ratios=(1, 20),
                                  wspace=0.005,
                                  hspace=0.01)  # Adjusted to create space for row annotation and top markers
        ax2 = fig.add_subplot(ax[1, 0])
        ax1 = fig.add_subplot(ax[0, 0])

        heat1 = sns.heatmap(cell_labels_df,
                                yticklabels=False,
                                xticklabels=False,
                                cmap=sns.color_palette(groups_to_colors.values()),
                                cbar=False,
                                ax=ax1)
        ax1.set_xticks(label_to_position.values)
        ax1.set_xticklabels(label_to_position.index.values, rotation=45, ha="left")
        ax1.xaxis.tick_top()

        heat2 = sns.heatmap(binary_matrix,
                                vmin=0,
                                vmax=1,
                                cmap='YlGn',
                                xticklabels=False,
                                yticklabels=False,  # Disable yticklabels here
                                cbar=True,
                                cbar_kws={"shrink": 0.5},
                                ax=ax2)
        ax2.collections[0].colorbar.set_label("Z-Score Normalized Expression")
        ax2.set_yticks(label_to_position2.values)
        ax2.set_yticklabels(label_to_position2.index.values, ha="right")
        ax2.yaxis.tick_left()

        fig.suptitle('Metadata Heatmap', fontsize=16)
        # fig.tight_layout()
        plt.savefig(path_to_save_figure, dpi=10, bbox_inches='tight', format="pdf", pad_inches=0.9)

    except Exception as e:
        print(str(e))
        print("Warning! Failed to run plot_metadata_heatmap. See above Exception.")


def plot_markers_df_subplots(marker_heatmap, markers_df, clusters, groups_to_colors, path_to_save_figure):
    """
    Plots a heatmap of the MarkerFinder results
    Arguments
    ---------
    marker_heatmap : pandas.DataFrame
        Data to plot. The rows must intersect with features in the \
        ordered_markers_df and the columns must intersect with the cells in clusters.
    markers_df : pandas.DataFrame
        The markers_df where "marker" column are the genes and the "top_cluster" is the cluster to which the gene belongs to  (Required)
    clusters : pandas.DataFrame
        The cluster assignments of the samples, index of this data frame are the sample names (Required)
    groups_to_colors : dict {str:str}
        Dictionary of cluster-names to assigned colors (hexadecimal value) \
        The pyInfinityFlow.Plotting_Utilities.assign_rainbow_colors_to_groups \
        can be used to generate this dictionary from a list of clusters. \
        (Required)
    path_to_save_figure : str
        The path to save the figure. (Required)

    Returns
    -------
    None
        Saves the heatmap to the file specified by path_to_save_figure
    """
    try:
        # Filter the input data matrix to the markers and cells of interest
        input_df = marker_heatmap.loc[markers_df["marker"].values, clusters.index]
        input_df = input_df.astype(float)

        # Map the order of the groups_to_colors to make a cmap
        groups_to_order = pd.Series(list(range(len(groups_to_colors))),
                                    index=groups_to_colors.keys())

        # Build top heatmap to label clusters
        cell_labels_df = pd.DataFrame({"cluster": [groups_to_order[item] for item in clusters["cluster"].values]},
                                      index=clusters.index).T

        # Build row annotation dataframe
        row_labels_df = pd.DataFrame({"cluster": [groups_to_order[item] for item in markers_df["top_cluster"].values]},
                                     index=markers_df["marker"].values)

        # Get top markers for each cluster
        top_markers = markers_df.groupby('top_cluster')['marker'].first()

        # Build df to add cluster label ticks
        label_to_position = pd.pivot_table(pd.DataFrame({
            "cluster": clusters["cluster"].values,
            "position": list(range(clusters.shape[0]))}),
            index="cluster",
            values="position",
            aggfunc=np.mean)["position"].sort_values()

        # Build df to add marker label ticks
        label_to_position2 = pd.pivot_table(pd.DataFrame({
            "cluster": markers_df["top_cluster"].values,
            "position": list(range(markers_df.shape[0]))}),
            index="cluster",
            values="position",
            aggfunc=np.mean)["position"].sort_values()

        # Z-score normalize the expression matrix
        z_df = input_df.apply(lambda x: pd.Series(zscore(x.values),
                                                  index=x.index.values), axis=1)

        # Build the heatmap
        plt.close("all")
        fig, axs = plt.subplots(2, 2, figsize=(18, 10),
                                gridspec_kw={'width_ratios': [1, 20], 'height_ratios': [1, 20], 'wspace': 0.005,
                                             'hspace': 0.01})
        ax1, ax2, ax3 = axs[0, 1], axs[1, 1], axs[1, 0]  # ax3 remains the same

        heat1 = sns.heatmap(cell_labels_df,
                            yticklabels=False,
                            xticklabels=False,
                            cmap=sns.color_palette(groups_to_colors.values()),
                            cbar=False,
                            ax=ax1)
        ax1.set_xticks(label_to_position.values)
        ax1.set_xticklabels(label_to_position.index.values, rotation=45, ha="left")
        ax1.xaxis.tick_top()

        heat2 = sns.heatmap(z_df,
                            vmin=-3,
                            vmax=3,
                            cmap=blue_black_yellow_cmap,  # Assuming you have blue_black_yellow_cmap defined
                            xticklabels=False,
                            yticklabels=False,  # Disable yticklabels here
                            cbar=True,
                            cbar_kws={"shrink": 0.5},
                            ax=ax2)
        ax2.collections[0].colorbar.set_label("Z-Score Normalized Expression")

        heat3 = sns.heatmap(row_labels_df,
                            xticklabels=False,
                            yticklabels=False,  # Disable yticklabels here as well
                            cmap=sns.color_palette(groups_to_colors.values()),
                            cbar=False,
                            ax=ax3)
        ax3.set_yticks(label_to_position2.values)
        ax3.set_yticklabels(label_to_position.index.values)

        fig.suptitle('Marker Finder Heatmap', fontsize=16)

        # Save the figure with compression
        plt.savefig(path_to_save_figure, dpi=100, bbox_inches='tight', format="pdf", pad_inches=0.9)

    except Exception as e:
        print(str(e))
        print("Warning! Failed to run plot_markers_df. See above Exception.")
