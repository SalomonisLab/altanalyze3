import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from itertools import chain


# filter out mitochondrial and ribosomal genes
def filter_out_non_coding_genes(adata):

    # remove the mitochondrial and ribosomal genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribosomal"] = adata.var_names.str.startswith(("RPS", "RPL"))

    adata = adata[:, adata.var['mt'] == False]
    adata = adata[:, adata.var['ribosomal'] == False]

    # remove genes that have a dot in their name
    adata = adata[:, ~adata.var_names.str.contains("\.")]

    # remove genes whose name starts with "Gm" or "GM"
    adata = adata[:, ~adata.var_names.str.startswith(("Gm", "GM"))]

    # remove genes whose last letter is "y" or "Y"
    adata = adata[:, ~adata.var_names.str.endswith(("y", "Y"))]

    # remove genes whose first 3 letters are HLA or Hla
    adata = adata[:, ~adata.var_names.str.startswith(("HLA", "Hla"))]

    # remove "xis" genes
    adata = adata[:, ~adata.var_names.str.startswith(("Xis", "XIS", "xis"))]

    # remove genes whose first 3 letters are "TSI" or "Tsi"
    adata = adata[:, ~adata.var_names.str.startswith(("TSI", "Tsi"))]

    # keep only the protein coding genes
    adata = adata[:, adata.var['protein_coding'] == True]

    return adata


# Variance based feature selection
def variance_based_feature_selection(adata, fold_threshold=0.5, samples_differing=3):

    adata.var['udon_hvg'] = True
    pseudobulk_folds_df = adata.varm['pseudobulk_folds']

    for gene in adata.var_names:
        gene_vector = pseudobulk_folds_df.loc[gene, :]
        gene_vector = gene_vector.sort_values(ascending=True)

        n_lowest = gene_vector.head(samples_differing).iloc[samples_differing-1]
        n_highest = gene_vector.tail(samples_differing).iloc[0]
        fold = n_highest - n_lowest
        adata.var.loc[gene, 'udon_hvg'] = fold > fold_threshold

    return adata


# Inter correlation step do this step if only less than 8000 to start with -- cut off for this 0.4 -- should bring it down to 1500
def intercorrelation_based_feature_selection(adata, corr_threshold=0.4, corr_n_events=5):

    adata_og = adata.copy()
    adata = adata[:, adata.var['udon_hvg'] == True]

    pseudobulk_folds_df = adata.varm['pseudobulk_folds']
    # Calculate correlation matrix
    pseudobulk_folds_df = pseudobulk_folds_df.transpose()

    corr_matrix = pseudobulk_folds_df.corr(method="pearson")

    # Get columns with correlation > 0.4
    correlated_features = (corr_matrix > corr_threshold).sum() > corr_n_events

    # Update adata.var with correlated features
    adata_og.var['correlated_genes'] = False
    adata_og.var.loc[correlated_features.index[correlated_features], 'correlated_genes'] = True

    return adata_og


# pca based variable feature selection
def pca_feature_selection(adata, corr_threshold=0.4, n_components=30):

    adata_og = adata.copy()
    adata = adata[:, adata.var['correlated_genes'] == True]

    exp = adata.varm['pseudobulk_folds']
    pca_obj = PCA(n_components, svd_solver='full')

    pca_obj.fit_transform(exp.transpose())
    loadings = pd.DataFrame(pca_obj.components_.T, index=exp.index)
    correlated_features_all_pc = []

    for pc in range(1, n_components + 1):
        # print(pc)
        pc_i = pd.Series(loadings.loc[:, pc - 1])  # get the genes and their loadings value from the given PC
        pc_i_sorted = pc_i.sort_values(ascending=False)

        pc_top200 = pc_i_sorted.head(200).index.tolist()  # top 200 features
        pc_bottom200 = pc_i_sorted.tail(200).index.tolist()  # bottom 200 features

        pc_400features = pc_top200 + pc_bottom200

        exp_400_features_pc = exp.loc[pc_400features, :]  # limit the expression to the above calculated 400 genes for that PC

        cor_400features_matrix = exp_400_features_pc.T.corr(method="pearson")
        cor_400features_matrix = np.tril(cor_400features_matrix, k=-1)
        idx_cor_400features_binary = abs(cor_400features_matrix) > corr_threshold
        idx_cor_400features_binary = idx_cor_400features_binary.astype(int)
        idx_cor_400features_binary = pd.DataFrame(idx_cor_400features_binary)
        idx_cor_400features_binary.index = exp_400_features_pc.index.tolist()
        idx_cor_400features_binary.columns = exp_400_features_pc.index.tolist()

        z = idx_cor_400features_binary.sum(axis=0)
        z = z[z > 0].index.tolist()
        z2 = idx_cor_400features_binary.sum(axis=1)
        z2 = z2[z2 > 0].index.tolist()

        correlated_events_pc = list(set(z) & set(z2))
        correlated_features_all_pc.append(correlated_events_pc)

    correlated_events_all_pc_unique = list(chain.from_iterable(correlated_features_all_pc))
    correlated_events_all_pc_unique = list(set(correlated_events_all_pc_unique))

    # Update adata.var with pca_selected_genes
    adata_og.var['pca_selected_genes'] = False
    adata_og.var.loc[correlated_events_all_pc_unique, 'pca_selected_genes'] = True

    return adata_og


def feature_selection_wrapper(adata, species, gtf_file_path=None, fold_threshold=1, samples_differing=3, intercorr_threshold=0.4, corr_n_events=5, pca_corr_threshold=0.4, n_components=30):

    # identify_protein_coding_genes
    print("Identifying protein coding genes...")
    adata = identify_protein_coding_genes(adata, species=species, gtf_file_path=gtf_file_path)

    print("before protein coding genes")
    print(adata.shape[1])
    # filter_out_non_coding_genes
    print("Filtering out non-coding genes...")
    adata = filter_out_non_coding_genes(adata)
    print("after protein coding genes")
    print(adata.shape[1])

    # variance_based_feature_selection
    print("Performing variance-based feature selection...")
    adata = variance_based_feature_selection(adata, fold_threshold=fold_threshold, samples_differing=samples_differing)

    print("after variance based feature selection")
    # get length of variables for which adata.var['udon_hvg'] is True
    print(adata.var['udon_hvg'].sum())

    # intercorrelation_based_feature_selection
    print("Performing intercorrelation-based feature selection...")
    adata = intercorrelation_based_feature_selection(adata, corr_threshold=intercorr_threshold, corr_n_events=corr_n_events)
    # filter the adata for only the correlated features/genes
    print("after intercorrelation based feature selection")
    # get length of variables for which adata.var['correlated_genes'] is True
    print(adata.var['correlated_genes'].sum())

    # pca_feature_selection
    print("Performing PCA-based feature selection...")
    adata = pca_feature_selection(adata, corr_threshold=pca_corr_threshold, n_components=n_components)
    print("after pca feature selection")
    # get length of variables for which adata.var['pca_selected_genes'] is True
    print(adata.var['pca_selected_genes'].sum())

    print("feature selection completed!")

    return adata


def identify_protein_coding_genes(adata, species, gtf_file_path=None):

    # read gtf file
    if gtf_file_path is not None:
        # read gtf file
        gtf_df = read_gtf_file(gtf_file_path)

        # Filter protein-coding genes
        protein_coding_genes = gtf_df[
            (gtf_df['feature'] == 'gene') & (gtf_df['attributes'].str.contains('protein_coding'))]
        protein_coding_gene_names = protein_coding_genes['attributes'].str.extract(r'gene_name "(.+?)";',
                                                                          expand=False).dropna().unique()

    else:
        gtf_df = pd.read_csv('ProteinCoding-Hs-Mm.txt', sep="\t", names=["GeneName"])
        # Filter based on species
        if species == "Hs":
            # Filter uppercase entries or those starting with "ENSG"
            gtf_df = gtf_df[(gtf_df['GeneName'].str.isupper()) | (gtf_df['GeneName'].str.startswith('ENSG'))]
        elif species == "Mm":
            # Filter uppercase entries or those starting with "ENSG"
            gtf_df = gtf_df[~(gtf_df['GeneName'].str.isupper()) | (gtf_df['GeneName'].str.startswith('ENSMU'))]

        protein_coding_gene_names = gtf_df['GeneName'].unique()

    gene_list = adata.var_names
    protein_coding_genes_in_list = [gene for gene in gene_list if gene.split("_")[0] in protein_coding_gene_names]

    # label protein_coding_genes_in_list as "protein_coding" in adata.var
    adata.var['protein_coding'] = False
    adata.var.loc[protein_coding_genes_in_list, 'protein_coding'] = True

    return adata


def read_gtf_file(gtf_file_path):

    columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
    df = pd.read_csv(gtf_file_path, sep='\t', comment='#', names=columns, dtype={"seqname": str})

    return df
