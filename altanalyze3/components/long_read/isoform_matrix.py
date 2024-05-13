import os
import sys
import csv
import anndata as ad
import pandas as pd
import numpy as np
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from scipy.sparse import csr_matrix
from scipy.io import mmread
from scipy.sparse import lil_matrix
import collections
import argparse
from tqdm import tqdm # progress bar
sys.path.insert(1, os.path.join(sys.path[0], '..'))

"""
Module for long-read sparse matrix single-cell processing.
"""

def mtx_to_adata(int_folder, gene_is_index, feature, feature_col, barcode, barcode_col, matrix, rev=False):
    """ Import isoform sparse matrices and key by isoform UID """
    feature_path = os.path.join(int_folder, feature)
    barcode_path = os.path.join(int_folder, barcode)
    matrix_path = os.path.join(int_folder, matrix)

    features = pd.read_csv(feature_path, sep='\t', header=None, index_col=feature_col)
    barcodes = pd.read_csv(barcode_path, sep='\t', header=None, index_col=barcode_col)
    mat = csr_matrix(mmread(matrix_path).tocsc())

    # If the barcodes are reverse complemented relative to standard 10x barcodes
    if rev:
        barcodes.index = barcodes.index.to_series().apply(reverse_complement_seq)

    # Remove gene suffix from isoform feature names if present
    features.index = features.index.to_series().apply(remove_feature_suffix)

    adata = ad.AnnData(X=mat.T, obs=pd.DataFrame(index=barcodes.index), var=pd.DataFrame(index=features.index))
    return adata

def append_sample_name(adata, sample_name):
    adata.obs_names = adata.obs_names.astype(str) + "_" + sample_name
    return adata

def import_barcode_clusters(barcode_cluster_dir):
    """ Import cell barcode to cluster associations and extract/key by sample name """
    barcode_sample_dict = {}
    df = pd.read_csv(barcode_cluster_dir, sep='\t', header=None, names=['barcode_cluster', 'cluster'])
    df[['barcode', 'sample_name']] = df['barcode_cluster'].str.split('.', expand=True)
    df['barcode'] = df['barcode'].apply(lambda x: x.split('.')[0])
    barcode_sample_dict = {sample: group.set_index('barcode')[['cluster']] for sample, group in df.groupby('sample_name')}
    return barcode_sample_dict

def reverse_complement_seq(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    nucleotide_part, suffix = seq.split('-')
    rev_comp_nucleotide_part = ''.join(complement[base] for base in reversed(nucleotide_part))
    return f"{rev_comp_nucleotide_part}-{suffix}"

def remove_feature_suffix(feature):
    return feature.split(':')[0]

def calculate_barcode_match_percentage(adata, barcode_clusters):
    # Verify that the cell barcodes maximally map to provided cell barcode annotations (no naming conflicts)
    matching_barcodes = set(adata.obs_names) & set(barcode_clusters.index)
    print(len(adata.obs_names),'matrix barcodes |',len(barcode_clusters.index),'annotation barcodes |',len(matching_barcodes),'overlap')

    matched_barcodes = set(adata.obs_names) & set(barcode_clusters.index)
    total_barcodes_in_clusters = len(barcode_clusters.index)
    percentage = (len(matched_barcodes) / total_barcodes_in_clusters) * 100
    print(f'Barcode mapping percentage: {percentage:.2f}%')

    adata = adata[sorted(matching_barcodes), :]
    adata.obs = adata.obs.join(barcode_clusters, how='inner')
    return adata

def pseudo_cluster_counts(combined_adata, cell_threshold=5, count_threshold=5, compute_tpm=False, tpm_threshold=1):
    # Convert sparse matrix to dense DataFrame for easier manipulation
    junction_counts_df = pd.DataFrame(combined_adata.X.toarray(), index=combined_adata.obs_names, columns=combined_adata.var_names)
    # Group by cluster
    grouped = junction_counts_df.groupby(combined_adata.obs['cluster'], observed=True)
    print("Cell counts per cluster:")
    print(grouped.size())
    # Filter for groups having at least 5 cells
    filtered_summed_groups = grouped.filter(lambda x: len(x) >= cell_threshold)
    # Sum the values post filtering
    summed_groups = filtered_summed_groups.groupby(combined_adata.obs['cluster'], observed=True).sum()
    # Filter out junctions with fewer than 5 reads
    summed_groups = summed_groups.loc[:, summed_groups.sum(axis=0) >= count_threshold]
    filtered_summed_groups_transposed = summed_groups.transpose()

    if compute_tpm:
        tpm_values = summed_groups.div(summed_groups.sum(axis=1), axis=0) * 1e6
        # Compute gene sums only once
        gene_sums = summed_groups.T.groupby(lambda x: x.split(":")[0]).transform('sum').T
        # Prepare isoform ratios efficiently
        isoform_ratios = summed_groups.div(gene_sums)
        # Filter isoform ratios where gene TPM < 1
        gene_tpms = tpm_values.T.groupby(lambda x: x.split(":")[0]).sum().T
        mask = gene_tpms >= tpm_threshold
        mask = mask.reindex(columns=tpm_values.columns, method='ffill')
        isoform_ratios = isoform_ratios.where(mask).transpose()
        print('Pseudobulk cluster junction counts exported')
        return filtered_summed_groups_transposed, tpm_values, isoform_ratios
    else:
        return filtered_summed_groups_transposed