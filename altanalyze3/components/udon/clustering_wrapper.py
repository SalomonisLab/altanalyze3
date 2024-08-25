import os
import time
from nmf import *
from markerFinder import *
from linearSVM import *
from feature_selection import *


def clustering_wrapper(adata, output_filename, rank=None, min_group_size=3, top_n_genes=60, rho_threshold=0.2, marker_finder_rho=0.2, write_files=True):

    if write_files is True:
        os.mkdir(output_filename)
        cd = os.getcwd()

    fold_matrix = adata.varm['pseudobulk_folds']

    # get var names for which 'pca_selected_genes' is True
    highly_variable_genes = adata.var.index[adata.var['correlated_genes'] == True]
    print("using corr genes...")
    # filter the full fold matrix with the highly variable events (currently comes from the file PCAbasedFeatureSelection.py ran before this function)
    fold_matrix_hvg = fold_matrix.loc[highly_variable_genes, :]

    # when rank is not 2 or not provided by the user (that they trust more), then automatically calculate number of clusters for NMF to determine (aka rank)
    # requires the determineNMFRank.py script
    if rank is None:
        est_k = determine_nmf_ranks(fold_matrix_hvg)
        rank = est_k

    # run and time the NMF analysis using the runNMF.py script
    start_time = time.time()
    print(("Running NMF analysis with rank set to " + str(rank) + "..."))
    basis_matrix, nmf_clusters = run_nmf(df=fold_matrix_hvg, rank=rank)
    print("--- %s seconds ---" % (time.time() - start_time))

    # remove clusters with less than 3 (min_groups_size) samples/pseudobulk folds in them
    # Count the number of samples in each cluster
    cluster_counts = nmf_clusters['cluster'].value_counts()

    # Filter clusters that have more than 3 samples
    large_clusters = cluster_counts[cluster_counts > min_group_size].index

    # Keep only samples that belong to the large clusters
    nmf_clusters = nmf_clusters[nmf_clusters['cluster'].isin(large_clusters)]

    fold_matrix = fold_matrix.loc[:, nmf_clusters.index]
    fold_matrix_hvg = fold_matrix_hvg.loc[:, nmf_clusters.index]

    start_time = time.time()
    print("Finding marker genes per cluster...")
    markers_df_og, markers_df, marker_heatmap = marker_finder_wrapper(input_df=fold_matrix.transpose(), groups=nmf_clusters, top_n=top_n_genes, rho_threshold=rho_threshold, marker_finder_rho=marker_finder_rho)
    print("--- %s seconds ---" % (time.time() - start_time))

    adata.uns['clusters_preSVM'] = nmf_clusters
    adata.uns['udon_marker_genes_full_preSVM'] = markers_df_og
    adata.uns['udon_marker_genes_top_n_preSVM'] = markers_df
    adata.uns['marker_heatmap_preSVM'] = marker_heatmap

    # remove the clusters with less than 10 markers
    # Count the number of markers in each cluster
    mf_gene_counts = markers_df['top_cluster'].value_counts()

    # Filter clusters that have more than 3 samples
    robust_clusters = mf_gene_counts.index

    # Keep only samples that belong to the large clusters
    markers_df = markers_df[markers_df['top_cluster'].isin(robust_clusters)]

    # filter the nmf_clusters with robust_clusters too
    nmf_clusters = nmf_clusters[nmf_clusters['cluster'].isin(robust_clusters)]

    start_time = time.time()
    print(("Determined " + str(len(markers_df)) + " unique marker genes across all clusters"))
    fold_matrix_with_mf_genes = fold_matrix.loc[markers_df['marker'], nmf_clusters.index]
    print("--- %s seconds ---" % (time.time() - start_time))

    # Linear SVM/SVC training dataset generation by creating centroids for the NMF clusters

    centroids = generate_train_data(fold_matrix_with_mf_genes=fold_matrix_with_mf_genes, groups=nmf_clusters)  # genes by clusters matrix
    centroids = centroids.dropna(axis=0)
    fold_matrix_with_mf_genes = fold_matrix_with_mf_genes.loc[centroids.index, :]

    start_time = time.time()
    print(("Running Linear SVM on final " + str(np.shape(centroids)[1]) + " clusters with..."))
    final_clusters = classify(train=centroids, fold_matrix_with_mf_genes=fold_matrix_with_mf_genes, groups=nmf_clusters)
    print("--- %s seconds ---" % (time.time() - start_time))

    fold_matrix = fold_matrix.loc[:, final_clusters.index]

    start_time = time.time()
    print("Finding final marker genes per cluster...")
    markers_df_og, markers_df, marker_heatmap = marker_finder_wrapper(input_df=fold_matrix.transpose(),
                                                                      groups=final_clusters, top_n=top_n_genes,
                                                                      rho_threshold=rho_threshold, marker_finder_rho=marker_finder_rho)
    print("--- %s seconds ---" % (time.time() - start_time))

    # clean the marker heatmap file for the final clusters in mf genes
    robust_clusters_final = markers_df['top_cluster'].value_counts(sort=False)
    robust_clusters_final = robust_clusters_final[robust_clusters_final > 15].index

    # Keep only samples that belong to the large clusters
    final_clusters_clean = final_clusters[final_clusters['cluster'].isin(robust_clusters_final)]
    markers_df_clean = markers_df[markers_df['top_cluster'].isin(robust_clusters_final)]

    final_samples = pd.concat([pd.Series('row_clusters-flat'), pd.Series(final_clusters_clean.index)])
    final_rows = pd.concat([pd.Series('column_clusters-flat'), pd.Series(markers_df_clean.index)])
    marker_heatmap = marker_heatmap.loc[final_rows, final_samples]

    adata.uns['udon_clusters'] = final_clusters_clean
    adata.uns['udon_marker_genes_full'] = markers_df_og
    adata.uns['udon_marker_genes_top_n'] = markers_df_clean
    adata.uns['marker_heatmap'] = marker_heatmap

    if write_files is True:
        mf_genes_path = os.path.join(cd, output_filename, "marker_genes_all.txt")
        mf_heatmap_path = os.path.join(cd, output_filename, "marker_heatmap.txt")
        final_clusters_path = os.path.join(cd, output_filename, "udon_clusters.txt")
        markers_df_og.to_csv(mf_genes_path, sep="\t")
        marker_heatmap.to_csv(mf_heatmap_path, sep="\t")
        final_clusters.to_csv(final_clusters_path, sep="\t")

    return adata


def find_udon_clusters(adata, output_filename, species, gtf_file_path=None, fold_threshold=1, samples_differing=3, intercorr_threshold=0.4, corr_n_events=5, pca_corr_threshold=0.4, n_components=30, rank=None, min_group_size=3, delta=0.1, p_val=0.05, min_differential_genes=100, top_n_differential_genes=150):

    print("beginning feature selection...")
    adata = feature_selection_wrapper(adata, species=species, gtf_file_path=gtf_file_path, fold_threshold=fold_threshold, samples_differing=samples_differing, intercorr_threshold=intercorr_threshold, corr_n_events=corr_n_events, pca_corr_threshold=pca_corr_threshold, n_components=n_components)

    print("performing NMF clustering...")
    final_clusters, marker_genes = clustering_wrapper(adata, output_filename=output_filename, rank=rank, min_group_size=min_group_size)

    adata.uns['udon_clusters'] = final_clusters
    adata.uns['udon_marker_genes'] = marker_genes

    return adata
