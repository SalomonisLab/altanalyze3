import pandas as pd
from classificationFunctions import *
import matplotlib.pyplot as plt
from tqdm import tqdm
import warnings
import seaborn as sns


def determine_control_cells(adata, groups_col='cell_type', train_col='disease', train_id='train', test_id='test', hvg_matrix_field='top_n_hvg', save_obs_key='cosine_similarity',donor_key='donor_ids',disease_key='disease_type'):

    return adata


# find the top 500 highly variable features per cell type
def determine_cell_type_hvg(adata, groups_col='cell_type', train_col='disease', train_id='control', n_top_genes=500):
    adata_train = adata[adata.obs[train_col] == train_id]

    # Initialize an empty dictionary to store highly variable genes for each cell type
    hvg_dict = {}

    # Get unique cell types
    cell_types = adata_train.obs[groups_col].unique()

    # Loop through each cell type
    for cell_type in cell_types:
        print(cell_type)
        # Subset the AnnData object
        adata_subset = adata_train[adata_train.obs[groups_col] == cell_type]

        # Check the number of cells
        if adata_subset.shape[0] < 50:
            warnings.warn(
                f"Cell type '{cell_type}' has less than 50 cells ({adata_subset.shape[0]} cells). Highly variable gene computation may not be reliable.")

        if adata_subset.shape[0] < 10:
            hvg_dict[cell_type] = []
            continue

        # Compute highly variable genes
        sc.pp.highly_variable_genes(adata_subset, flavor='seurat_v3', n_top_genes=n_top_genes, layer='counts',span=0.8)

        # Get the top 500 highly variable genes
        hvg = adata_subset.var[adata_subset.var['highly_variable']].index.tolist()

        # Store the result in the dictionary
        hvg_dict[cell_type] = hvg

    # Convert the dictionary to a DataFrame
    hvg_df = pd.DataFrame.from_dict(hvg_dict, orient='index').transpose()

    adata.uns['top_n_hvg'] = hvg_df

    return adata


# find the degs per cell type
# find the degs per cell type
def determine_cell_type_deg(adata, groups_col='cell_type', train_col='train_test_col', train_id='control', test_id='disease', pval_cutoff=0.05, log2fc_min=1,
                                        log2fc_max=None, pt_min=0.10, save_key='deg_dict'):

    # Initialize an empty dictionary to store highly variable genes for each cell type
    deg_dict = {}
    full_deg_dict = {}

    # Get unique cell types
    cell_types = adata.obs[groups_col].unique()

    # Loop through each cell type
    for cell_type in cell_types:
        print(cell_type)
        # Subset the AnnData object
        adata_subset = adata[adata.obs[groups_col] == cell_type]

        # Check the number of cells
        if adata_subset.shape[0] < 50:
            warnings.warn(
                f"Cell type '{cell_type}' has less than 50 cells ({adata_subset.shape[0]} cells). DEG computation may not be reliable.")

        if adata_subset.shape[0] < 5:
            deg_dict[cell_type] = []
            continue

        # Compute DEGs using wilcoxon test
        sc.tl.rank_genes_groups(adata_subset, groupby=train_col, groups=[test_id], reference=train_id, method='wilcoxon', key_added="wilcoxon_deg", pts=True)

        # Get the DEGs
        deg_df = sc.get.rank_genes_groups_df(adata_subset, group=None, key='wilcoxon_deg', pval_cutoff=pval_cutoff, log2fc_min=log2fc_min,
                                        log2fc_max=log2fc_max)

        deg_df = deg_df.loc[deg_df['pct_nz_group'] >= pt_min, :]

        if deg_df.shape[0] < 20:
            deg_df = sc.get.rank_genes_groups_df(adata_subset, group=None, key='wilcoxon_deg')
            deg_df = deg_df.loc[deg_df['pvals'] < pval_cutoff]
            deg_df = deg_df.loc[deg_df['pct_nz_group'] >= pt_min]

        # Store the result in the dictionary
        full_deg_dict[cell_type] = deg_df

        degs = deg_df['names'].tolist()

        # Store the result in the dictionary
        deg_dict[cell_type] = degs

    full_df_key = save_key + "_full_df"
    adata.uns[full_df_key] = full_deg_dict

    deg_all_cell_types_df = pd.DataFrame.from_dict(deg_dict, orient='index').transpose()
    adata.uns[save_key] = deg_all_cell_types_df

    return adata


def determine_control_cells_cosine_similarity(adata, groups_col='cell_type', method='cosine_similarity', train_col='disease', train_id='train', test_id='test', hvg_matrix_field='top_n_hvg', save_obs_key='cosine_similarity',donor_key='donor_ids',disease_key='disease_type'):

    # split into train dataset and test dataset
    train_controls = adata[adata.obs[train_col] == train_id, :]
    test = adata[adata.obs[train_col] == test_id, :]

    print("finished train test split")

    cell_types = adata.obs[groups_col].unique()

    for cluster in tqdm(cell_types, desc="Predicting cell types with cosine similarity..."):

        print(cluster)
        top_500_var_feats_all_cell_types = adata.uns[hvg_matrix_field]
        var_feats_s = top_500_var_feats_all_cell_types.loc[:, cluster]

        # Remove missing values (NaNs)
        var_feats_s = var_feats_s.dropna()

        if len(var_feats_s) < 30:
            warnings.warn(
                f"Cell type '{cluster}' has less than 50 highly variable genes ({len(var_feats_s)} genes). not performing similarity metric.")
            continue

        # subset with the top 500 feats and the cells of that cluster
        train_controls_s = train_controls[train_controls.obs[groups_col] == cluster, var_feats_s]
        test_s = test[test.obs[groups_col] == cluster, var_feats_s]

        if np.shape(train_controls_s)[0] < 10:
            warnings.warn(
                f"Cell type '{cluster}' has less than 10 cells ({np.shape(train_controls_s)[0]} cells). not performing similarity metric.")
            continue

        if method == 'cosine_similarity':
            train_adata, test_adata, predicted_labels, row_mean_similarity = gene_cosine_similarity(train_adata=train_controls_s, test_adata=test_s)
        elif method == 'pearson_corr':
            train_adata, test_adata, predicted_labels, row_mean_similarity = gene_pearson_corr(train_adata=train_controls_s, test_adata=test_s)

        row_mean_similarity_df = pd.DataFrame(row_mean_similarity, index=test_s.obs_names,
                                              columns=["cosine_similarity"])
        row_mean_similarity_df['disease_type'] = test_s.obs[disease_key]

        # predicted_labels, row_mean_similarity = check_similarity_distribution(mean_threshold=0.6, count_in_range_threshold=200, train_controls_s=train_controls_s, test_s=test_s, row_mean_similarity_df=row_mean_similarity_df)

        adata.obs.loc[test_s.obs_names, (save_obs_key + '_predicted_labels')] = predicted_labels
        adata.obs.loc[test_s.obs_names, (save_obs_key + '_similarity_values')] = row_mean_similarity

        test_s.obs.loc[test_s.obs_names, (save_obs_key + '_predicted_labels')] = predicted_labels
        test_s.obs.loc[test_s.obs_names, (save_obs_key + '_similarity_values')] = row_mean_similarity

        plt.close("all")
        plt.figure()
        sc.pl.pca(train_controls_s, color=donor_key, legend_loc='right margin', dimensions=[(0, 1)])
        train_pca_file = f"{cluster.replace(' ', '_')}_train_pca_corr.svg"
        plt.savefig(train_pca_file, bbox_inches='tight', format='svg')
        plt.close("all")

        train_pca_coords = pd.DataFrame(
            train_controls_s.obsm["X_pca"],
            index=train_controls_s.obs.index,  # Cell names from AnnData
            columns=[f"PC{i+1}" for i in range(train_controls_s.obsm["X_pca"].shape[1])])  # Column names: PC1, PC2, ...

        train_pca_coord_file = f"{cluster.replace(' ', '_')}_train_pca_coordinates.txt"
        # Write to a text file (tab-separated)
        train_pca_coords.to_csv(train_pca_coord_file, sep="\t", index=True, header=True)

        plt.figure()
        sc.pl.pca(test_s, color=[donor_key, disease_key, (save_obs_key + '_similarity_values')],
                  dimensions=[(0, 1)])
        test_pca_file = f"{cluster.replace(' ', '_')}_test_pca_corr.svg"
        plt.savefig(test_pca_file, bbox_inches='tight', format='svg')
        plt.close("all")
        print("finished predicting")

        test_pca_coords = pd.DataFrame(
            test_s.obsm["X_pca"],
            index=test_s.obs.index,  # Cell names from AnnData
            columns=[f"PC{i+1}" for i in range(test_s.obsm["X_pca"].shape[1])])  # Column names: PC1, PC2, ...

        test_pca_coord_file = f"{cluster.replace(' ', '_')}_test_pca_coordinates.txt"
        # Write to a text file (tab-separated)
        test_pca_coords.to_csv(test_pca_coord_file, sep="\t", index=True, header=True)

        plt.figure()
        sns.displot(row_mean_similarity_df, x="cosine_similarity", hue="disease_type",
                    multiple="stack", binwidth=0.01)

        # Calculate the mean of cosine_similarity
        mean_value = row_mean_similarity_df['cosine_similarity'].mean()
        median_value = row_mean_similarity_df['cosine_similarity'].median()

        # Add a vertical line at the mean value
        plt.axvline(mean_value, color='red', linestyle='--')
        plt.axvline(median_value, color='blue', linestyle='--')

        # Annotate the asterisk at the mean value on the x-axis
        plt.text(mean_value, 0, '*', color='red', fontsize=15, ha='center', va='bottom')

        file = f"{cluster.replace(' ', '_')}_cosine_similarity_distribution.svg"
        plt.savefig(file, bbox_inches='tight', format='svg')
        plt.close("all")

    return adata


def determine_control_cells_one_class_svm(adata, groups_col='cell_type', train_col='disease', train_id='train', test_id='test', hvg_matrix_field='top_n_hvg', svm_matrix='pca', save_obs_key='one_class_svm', donor_key='donor_ids', disease_key='disease_type'):

    # split into train dataset and test dataset
    train_controls = adata[adata.obs[train_col] == train_id, :]
    test = adata[adata.obs[train_col] == test_id, :]

    print("finished train test split")

    adata.obs[(save_obs_key + '_predicted_labels')] = ''
    adata.obs[(save_obs_key + '_decision_values')] = ''
    cell_types = adata.obs[groups_col].unique()

    for cluster in tqdm(cell_types, desc="Predicting cell types with one-class SVM..."):

        print(cluster)
        top_500_var_feats_all_cell_types = adata.uns[hvg_matrix_field]
        var_feats_s = top_500_var_feats_all_cell_types.loc[:, cluster]


        # Remove missing values (NaNs)
        var_feats_s = var_feats_s.dropna()

        if len(var_feats_s) < 50:
            warnings.warn(
                f"Cell type '{cluster}' has less than 50 highly variable genes ({len(var_feats_s)} genes). not performing one-class SVM.")
            continue

        # subset with the top 500 feats and the cells of that cluster
        train_controls_s = train_controls[train_controls.obs[groups_col] == cluster, var_feats_s]
        test_s = test[test.obs[groups_col] == cluster, var_feats_s]

        # one class svm
        train_adata, test_adata, predicted_labels, decision_values = one_class_svm(train_adata=train_controls_s, test_adata=test_s, nu=0.1, matrix=svm_matrix)

        adata.obs.loc[test_s.obs_names, (save_obs_key + '_predicted_labels')] = predicted_labels
        adata.obs.loc[test_s.obs_names, (save_obs_key + '_decision_values')] = decision_values

        test_s.obs.loc[test_s.obs_names, (save_obs_key + '_predicted_labels')] = predicted_labels
        test_s.obs.loc[test_s.obs_names, (save_obs_key + '_decision_values')] = decision_values

        print("finished predicting")

        train_pca_file = f"{cluster.replace(' ', '_')}_train_pca_svm.pdf"
        sc.pl.pca(train_controls_s, color=donor_key, legend_loc='right margin', dimensions=[(0, 1)], save=train_pca_file)

        test_pca_file = f"{cluster.replace(' ', '_')}_test_pca_svm.pdf"
        sc.pl.pca(test_s, color=[donor_key, disease_key, (save_obs_key + '_decision_values')],
                  dimensions=[(0, 1)], save=test_pca_file)

    return adata


def check_similarity_distribution(mean_threshold, count_in_range_threshold, train_controls_s, test_s, row_mean_similarity_df):

    # Calculate the mean of cosine_similarity
    mean_value = row_mean_similarity_df['cosine_similarity'].mean()

    # Assuming you have already calculated the mean_value without outliers
    lower_bound_range = mean_value - 0.01
    upper_bound_range = mean_value + 0.01

    # Count the number of entries within the range
    count_in_range = row_mean_similarity_df[(row_mean_similarity_df['cosine_similarity'] >= lower_bound_range) &
                                            (row_mean_similarity_df['cosine_similarity'] <= upper_bound_range)].shape[0]

    n_genes = 500

    while (mean_value < mean_threshold and count_in_range < count_in_range_threshold) and n_genes > 150:

        # Compute highly variable genes
        try:
            sc.pp.highly_variable_genes(train_controls_s, flavor='seurat_v3', n_top_genes=n_genes, layer='counts', span=0.8)
        except:
            break

        # Get the top 500 highly variable genes
        hvg = train_controls_s.var[train_controls_s.var['highly_variable']].index.tolist()

        train_controls_s = train_controls_s[:, hvg]
        test_s = test_s[:, hvg]

        train_adata, test_adata, predicted_labels, row_mean_similarity = gene_cosine_similarity(train_adata=train_controls_s, test_adata=test_s)

        row_mean_similarity_df = pd.DataFrame(row_mean_similarity, index=test_s.obs_names,
                                              columns=["cosine_similarity"])
        row_mean_similarity_df['disease_type'] = test_s.obs['disease_type']

        # Calculate the mean of cosine_similarity
        mean_value = row_mean_similarity_df['cosine_similarity'].mean()

        # Assuming you have already calculated the mean_value without outliers
        lower_bound_range = mean_value - 0.01
        upper_bound_range = mean_value + 0.01

        # Count the number of entries within the range
        count_in_range = row_mean_similarity_df[(row_mean_similarity_df['cosine_similarity'] >= lower_bound_range) &
                                                (row_mean_similarity_df[
                                                     'cosine_similarity'] <= upper_bound_range)].shape[0]
        n_genes = n_genes - 100

    row_mean_similarity = row_mean_similarity_df['cosine_similarity']
    predicted_labels = np.zeros(np.shape(row_mean_similarity)[0])
    predicted_labels[row_mean_similarity < mean_value] = 1

    return predicted_labels, row_mean_similarity

