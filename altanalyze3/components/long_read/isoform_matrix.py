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
from collections import defaultdict
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
    """Import cell barcode to cluster associations and extract/key by sample name."""
    if not isinstance(barcode_cluster_dir, list):
        barcode_cluster_dir = [barcode_cluster_dir]
    barcode_sample_dict_final = defaultdict(list)
    
    for dir in barcode_cluster_dir:
        df = pd.read_csv(dir, sep='\t', header=None, names=['barcode_cluster', 'cluster'])
        df[['barcode', 'sample_name']] = df['barcode_cluster'].str.split('.', expand=True)
        for sample, group in df.groupby('sample_name'):
            barcode_sample_dict_final[sample] = group.set_index('barcode')['cluster']
    return barcode_sample_dict_final

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

import pandas as pd
import anndata as ad

def tsv_to_adata(file_path: str):
    """Convert a TSV file with either (1) annot_transcript_id and counts or (2) isoform IDs and sample-level counts to an AnnData object."""
    df = pd.read_csv(file_path, sep='\t')
    
    if df.shape[1] < 2:
        raise ValueError("The TSV file must contain at least two columns.")

    if 'annot_transcript_id' in df.columns:
        # Handle the original format with 'annot_transcript_id' column
        transcript_ids = df['annot_transcript_id']
        counts = df.iloc[:, -1].values  # Assuming counts are in the last column
        count_df = pd.DataFrame(counts, index=transcript_ids, columns=['counts'])
        count_df = count_df.transpose()
    else:
        # Handle the standard count matrix format
        isoform_ids = df.iloc[:, 0]  # First column for isoform IDs
        counts = df.iloc[:, 1:].values  # All other columns for counts
        count_df = pd.DataFrame(counts, index=isoform_ids, columns=df.columns[1:])
        count_df = count_df.transpose()

    adata = ad.AnnData(X=count_df)
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
        tpm_values_transposed = tpm_values.transpose()
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
        return filtered_summed_groups_transposed, tpm_values_transposed, isoform_ratios
    else:
        return filtered_summed_groups_transposed

def pseudo_counts_no_cluster(combined_adata, count_threshold=5, compute_tpm=False, tpm_threshold=1):
    # Convert sparse matrix to dense DataFrame for easier manipulation
    junction_counts_df = pd.DataFrame(combined_adata.X.toarray(), index=combined_adata.obs_names, columns=combined_adata.var_names)
    
    # Filter out junctions with fewer than the specified count threshold in any sample
    filtered_counts = junction_counts_df.loc[:, (junction_counts_df >= count_threshold).any(axis=0)]
    
    if compute_tpm:
        # Compute TPM values for each sample
        total_counts = filtered_counts.sum(axis=1)
        tpm_values = (filtered_counts.T / total_counts).T * 1e6
        
        # Ensure the columns are correctly formatted as strings for grouping
        def extract_gene_id(x):
            return x.split(":")[0] if isinstance(x, str) else x
        
        # Calculate gene sums using transposed DataFrame
        gene_sums = filtered_counts.T.groupby(extract_gene_id).transform('sum').T
        isoform_ratios = filtered_counts / gene_sums
        
        # Filter isoform ratios where gene TPM < 1 for each sample using transposed DataFrame
        gene_tpms = tpm_values.T.groupby(extract_gene_id).transform('sum').T
        mask = gene_tpms >= tpm_threshold
        isoform_ratios = isoform_ratios.where(mask, np.nan)
        
        print('Junction counts and TPM values exported')
        return filtered_counts.T, tpm_values.T, isoform_ratios.T
    else:
        return filtered_counts.T


    def export_and_filter_pseudobulks(input_file, output_file, groups_file, cell_type_order_file=None):
        #Function to reduce a combined TPM, ratio or splicing matrix (with or without null values) to rows in which clusters have multiple samples with signal (>2). Assumes null = 0.

        # Step 1: Read in the optional cell type order file if provided
        if cell_type_order_file:
            with open(cell_type_order_file, 'r') as f:
                cell_type_order = f.read().splitlines()
        else:
            cell_type_order = None

        with open(input_file, 'r') as infile:
            header = infile.readline().strip().split('\t')
            samples = header[1:]  # Skip the first column (feature names)
            
            # Determine cell types from sample names
            cell_types = [sample.split('.')[0] for sample in samples]
            
            # If cell_type_order is provided, re-arrange columns
            if cell_type_order:
                reordered_columns = []
                reordered_cell_types = []
                for cell_type in cell_type_order:
                    cols = [sample for sample in samples if sample.startswith(cell_type+'.')]
                    reordered_columns.extend(cols)
                    reordered_cell_types.extend([cell_type] * len(cols))
            else:
                reordered_columns = samples
                reordered_cell_types = cell_types
            
            with open(output_file, 'w') as outfile:
                outfile.write('\t'.join(['Feature'] + reordered_columns) + '\n')

                valid_rows = []  # To collect all valid rows for later processing
                for line in infile:
                    fields = line.strip().split('\t')
                    feature = fields[0]
                    try:
                        values = list(map(float, fields[1:]))
                    except:
                        # for ratios
                        values = [float(value) if value != '' else 0 for value in fields[1:]]

                    # Step 3: Group values by cell type and check non-zero counts
                    cell_type_counts = {cell_type: 0 for cell_type in set(cell_types)}
                    
                    for value, cell_type in zip(values, cell_types):
                        if value > 0.1:
                            cell_type_counts[cell_type] += 1
                    
                    # If any cell type has at least 3 non-zero samples, keep the row
                    if any(count >= 3 for count in cell_type_counts.values()):
                        if cell_type_order:
                            reordered_values = [fields[samples.index(col) + 1] for col in reordered_columns]
                        else:
                            reordered_values = fields[1:]
                        valid_rows.append((feature, reordered_values))
                        outfile.write(f"{feature}\t" + "\t".join(map(str, reordered_values)) + "\n")
            
            # Step 6: Write the groups file
            with open(groups_file, 'w') as f:
                f.write('Original Column ID\tCell Type\n')
                for col, cell_type in zip(reordered_columns, reordered_cell_types):
                    f.write(f'{col}\t{cell_type}\n')

    def concatenate_h5ad_and_compute_pseudobulks(sample_files,collection_name = ''):

        # Initialize empty DataFrames to store results
        combined_pseudo_pdf = None
        combined_tpm = None
        combined_isoform_to_gene_ratio = None

        if collection_name!='':
            collection_name = '_'+collection_name

        # Loop through each sample file
        for sample_file in sample_files:
            sample_name = sample_file.split('.')[0]  # Extract the sample name from the file name
            print(f'Processing {sample_file}')
            
            # Load the data
            adata = ad.read_h5ad(sample_file)
            
            # Ensure indices are strings and observation names are unique
            adata.obs.index = adata.obs.index.astype(str)
            adata.obs_names_make_unique()

            # Generate pseudobulks for the current sample
            pseudo_pdf, tpm, isoform_to_gene_ratio = pseudo_cluster_counts(adata, cell_threshold=5, count_threshold=0, compute_tpm=True)
            
            # Rename the columns to include the sample name
            pseudo_pdf.columns = [f'{col}.{sample_name}' for col in pseudo_pdf.columns]
            tpm.columns = [f'{col}.{sample_name}' for col in tpm.columns]
            isoform_to_gene_ratio.columns = [f'{col}.{sample_name}' for col in isoform_to_gene_ratio.columns]
            
            # Ensure the original order of clusters is preserved
            if combined_pseudo_pdf is None:
                combined_pseudo_pdf = pseudo_pdf
                combined_tpm = tpm
                combined_isoform_to_gene_ratio = isoform_to_gene_ratio
            else:
                # Reindex the dataframes to match the combined dataframe's index (cluster names) and fill missing values with 0
                pseudo_pdf = pseudo_pdf.reindex(index=combined_pseudo_pdf.index, columns=pseudo_pdf.columns).fillna(0)
                tpm = tpm.reindex(index=combined_tpm.index, columns=tpm.columns).fillna(0)
                isoform_to_gene_ratio = isoform_to_gene_ratio.reindex(index=combined_isoform_to_gene_ratio.index, columns=isoform_to_gene_ratio.columns).fillna(0)
                
                # Concatenate the dataframes along columns
                combined_pseudo_pdf = pd.concat([combined_pseudo_pdf, pseudo_pdf], axis=1)
                combined_tpm = pd.concat([combined_tpm, tpm], axis=1)
                combined_isoform_to_gene_ratio = pd.concat([combined_isoform_to_gene_ratio, isoform_to_gene_ratio], axis=1)

        # Export to text file
        pseudo_counts = collection_name+'combined_pseudo_cluster_counts.txt'
        pseudo_tpm = collection_name+'combined_tpm.txt'
        pseudo_ratios = collection_name+'combined_isoform_to_gene_ratio.txt'
        combined_pseudo_pdf.to_csv(pseudo_counts, sep='\t')
        combined_tpm.to_csv(pseudo_tpm, sep='\t')
        combined_isoform_to_gene_ratio.to_csv(pseudo_ratios, sep='\t')
    return pseudo_counts,pseudo_tpm,pseudo_ratios
