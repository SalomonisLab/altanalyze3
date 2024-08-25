import pandas as pd
import numpy as np
import scanpy as sc
from classificationFunctions import *
from sklearn import metrics
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity
from tqdm import tqdm
import warnings


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


def determine_control_cells_cosine_similarity(adata, groups_col='cell_type', train_col='disease', train_id='control', hvg_matrix_field='top_n_hvg', save_obs_key='cosine_similarity',donor_key='donor_ids',disease_key='disease_type'):

    # split into train dataset and test dataset
    train_controls = adata[adata.obs[train_col] == train_id, :]
    test = adata[~(adata.obs[train_col] == train_id), :]

    print("finished train test split")

    adata.obs[(save_obs_key + '_predicted_labels')] = ''
    adata.obs[(save_obs_key + '_decision_values')] = ''
    cell_types = adata.obs[groups_col].unique()

    for cluster in tqdm(cell_types, desc="Predicting cell types with cosine similarity..."):

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

        train_adata, test_adata, predicted_labels, row_mean_similarity = gene_cosine_similarity(train_adata=train_controls_s, test_adata=test_s)
        adata.obs.loc[test_s.obs_names, (save_obs_key + '_predicted_labels')] = predicted_labels
        adata.obs.loc[test_s.obs_names, (save_obs_key + '_similarity_values')] = row_mean_similarity

        test_s.obs.loc[test_s.obs_names, (save_obs_key + '_predicted_labels')] = predicted_labels
        test_s.obs.loc[test_s.obs_names, (save_obs_key + '_similarity_values')] = row_mean_similarity

        plt.close("all")
        plt.figure()
        sc.pl.pca(train_controls_s, color=donor_key, legend_loc='right margin', dimensions=[(0, 1)])
        train_pca_file = f"{cluster.replace(' ', '_')}_train_pca.png"
        plt.savefig(train_pca_file, bbox_inches='tight')
        plt.close("all")

        plt.figure()
        sc.pl.pca(test_s, color=[donor_key, disease_key, (save_obs_key + '_similarity_values')],
                  dimensions=[(0, 1)])
        test_pca_file = f"{cluster.replace(' ', '_')}_test_pca.png"
        plt.savefig(test_pca_file, bbox_inches='tight')
        plt.close()
        print("finished predicting")

    return adata


def determine_control_cells_one_class_svm(adata, groups_col='cell_type', train_col='disease', train_id='control', hvg_matrix_field='top_n_hvg', svm_matrix='pca', save_obs_key='one_class_svm', donor_key='donor_ids',disease_key='disease_type'):

    # split into train dataset and test dataset
    train_controls = adata[adata.obs[train_col] == train_id, :]
    test = adata[~(adata.obs[train_col] == train_id), :]

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

        plt.close("all")
        plt.figure()
        sc.pl.pca(train_controls_s, color=donor_key, legend_loc='right margin', dimensions=[(0, 1)])
        train_pca_file = f"{cluster.replace(' ', '_')}_train_pca.png"
        plt.savefig(train_pca_file, bbox_inches='tight')
        plt.close("all")

        plt.figure()
        sc.pl.pca(test_s, color=[donor_key, disease_key, (save_obs_key + '_decision_values')],
                  dimensions=[(0, 1)])
        test_pca_file = f"{cluster.replace(' ', '_')}_test_pca.png"
        plt.savefig(test_pca_file, bbox_inches='tight')
        plt.close()

    return adata

