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
from collections import defaultdict
import argparse
from tqdm import tqdm # progress bar
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from joblib import Parallel, delayed

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


def return_cluster_order(barcode_cluster_dir):
    from collections import Counter, OrderedDict
    pair_counts = Counter()
    import glob
    cluster_orders = {}

    if not isinstance(barcode_cluster_dir, list):
        barcode_cluster_dir = [barcode_cluster_dir]
    barcode_sample_dict_final = defaultdict(list)

    for file_path in barcode_cluster_dir:
        df = pd.read_csv(file_path, sep='\t', header=None, names=['barcode', 'cluster'])
        unique_clusters = []
        for cluster in df['cluster']:
            if cluster not in unique_clusters:
                unique_clusters.append(cluster)
        cluster_orders[file_path] = unique_clusters
        
    default_file = max(cluster_orders, key=lambda x: len(cluster_orders[x]))
    default_order = cluster_orders[default_file]
    for file_path, clusters in cluster_orders.items():
        if file_path == default_file:
            continue  # Skip the default file
        for i, cluster in enumerate(clusters):
            if cluster not in default_order:
                # Insert the cluster based on the preceding cluster (if available)
                if i > 0 and clusters[i - 1] in default_order:
                    pos = default_order.index(clusters[i - 1]) + 1
                    default_order.insert(pos, cluster)
                elif i == 0 or clusters[i - 1] not in default_order:
                    # If no preceding cluster is found in the default order, append at the end
                    default_order.append(cluster)
    return default_order

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


def pseudo_cluster_counts(sample, combined_adata, cell_threshold=5, count_threshold=5, compute_tpm=False, tpm_threshold=1, status=True):
    if compute_tpm:
        dataType = 'isoform'
    else:
        dataType = 'junction'
    # Convert sparse matrix to dense DataFrame for easier manipulation
    junction_counts_df = pd.DataFrame(combined_adata.X.toarray(), index=combined_adata.obs_names, columns=combined_adata.var_names)
    # Group by cluster
    grouped = junction_counts_df.groupby(combined_adata.obs['cluster'], observed=True)
    if status:
        print("Cell counts per cluster:")
        print(grouped.size())
    # Filter for groups having at least cell_threshold cells (with progress bar)
    with tqdm(total=len(grouped), desc="Filtering clusters by cell threshold") as pbar:
        filtered_summed_groups = grouped.filter(lambda x: len(x) >= cell_threshold)
        pbar.update(len(grouped))
    # Sum the values post filtering (with progress bar)
    with tqdm(total=len(filtered_summed_groups), desc="Summing filtered groups") as pbar:
        summed_groups = filtered_summed_groups.groupby(combined_adata.obs['cluster'], observed=True).sum()
        pbar.update(len(filtered_summed_groups))
    # Filter out junctions with fewer than count_threshold reads (with progress bar)
    with tqdm(total=len(summed_groups.columns), desc="Filtering junctions by read count") as pbar:
        summed_groups = summed_groups.loc[:, summed_groups.sum(axis=0) >= count_threshold]
        pbar.update(len(summed_groups.columns))
    # Transpose and save filtered summed groups
    filtered_summed_groups_transposed = summed_groups.transpose()
    filtered_summed_groups_transposed.to_csv(f"{sample}_{dataType}.txt", sep='\t')
    if compute_tpm:
        # Compute TPM values (with progress bar)
        with tqdm(total=summed_groups.shape[0], desc="Calculating TPM values") as pbar:
            tpm_values = summed_groups.div(summed_groups.sum(axis=1), axis=0) * 1e6
            pbar.update(summed_groups.shape[0])
        tpm_values_transposed = tpm_values.transpose()
        # Compute gene sums only once (with progress bar)
        with tqdm(total=summed_groups.shape[1], desc="Computing gene sums") as pbar:
            gene_sums = summed_groups.T.groupby(lambda x: x.split(":")[0]).transform('sum').T
            pbar.update(summed_groups.shape[1])
        # Prepare isoform ratios efficiently
        isoform_ratios = summed_groups.div(gene_sums)
        # Filter isoform ratios where gene TPM < 1 (with progress bar)
        with tqdm(total=tpm_values.shape[1], desc="Filtering isoform ratios by gene TPM") as pbar:
            gene_tpms = tpm_values.T.groupby(lambda x: x.split(":")[0]).sum().T
            mask = gene_tpms >= tpm_threshold
            mask = mask.reindex(columns=tpm_values.columns, method='ffill')
            isoform_ratios = isoform_ratios.where(mask).transpose()
            pbar.update(tpm_values.shape[1])
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


from tqdm import tqdm

def export_and_filter_pseudobulks(input_file, output_file, cell_type_order=None):
    # Function to reduce a combined TPM, ratio, or splicing matrix (with or without null values) to rows in which clusters have multiple samples with signal (>2). Assumes null = 0.

    # Step 1: Read in the optional cell type order file if provided
    if cell_type_order is None:
        cell_type_order = None

    with open(input_file, 'r') as infile:
        header = infile.readline().strip().split('\t')
        samples = header[1:]  # Skip the first column (feature names)

        # Determine cell types from sample names (prefix before '.')
        cell_types = [sample.split('.')[0] for sample in samples]

        # If cell_type_order is provided, re-arrange columns
        if cell_type_order:
            reordered_columns = []
            for cell_type in cell_type_order:
                cols = [sample for sample in samples if sample.startswith(cell_type + '.')]
                reordered_columns.extend(cols)
        else:
            reordered_columns = samples

        total_lines = sum(1 for _ in open(input_file)) - 1  # Count lines for progress tracking (excluding header)

        with open(output_file, 'w') as outfile:
            outfile.write('\t'.join(['Feature'] + reordered_columns) + '\n')

            with tqdm(total=total_lines, desc="Processing lines") as pbar:
                for line in infile:
                    fields = line.strip().split('\t')
                    feature = fields[0]
                    
                    # Convert to float where possible, fall back to 0 for non-convertible
                    try:
                        values = list(map(float, fields[1:]))
                    except ValueError:
                        values = [float(value) if value else 0 for value in fields[1:]]

                    # Group values by cell type and check non-zero counts
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
                        outfile.write(f"{feature}\t" + "\t".join(map(str, reordered_values)) + "\n")

                    pbar.update(1)

def concatenate_h5ad_and_compute_pseudobulks(sample_files,collection_name = ''):

    # Initialize empty DataFrames to store results
    combined_pseudo_pdf = None
    combined_tpm = None
    combined_isoform_to_gene_ratio = None

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
        if collection_name == 'isoform':
            pseudo_pdf, tpm, isoform_to_gene_ratio = pseudo_cluster_counts(sample_name, adata, cell_threshold=5, count_threshold=0, compute_tpm=True, status=False)
        
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
            pseudo_tpm = collection_name+'_combined_pseudo_cluster_tpm.txt'
            pseudo_ratios = collection_name+'_combined_pseudo_cluster_iso-gene_ratio.txt'
            combined_tpm.to_csv(pseudo_tpm, sep='\t')
            combined_isoform_to_gene_ratio.to_csv(pseudo_ratios, sep='\t')
        else:
            pseudo_pdf = pseudo_cluster_counts(sample_name, adata, cell_threshold=5, count_threshold=0, compute_tpm=False, status=False)
            pseudo_pdf.columns = [f'{col}.{sample_name}' for col in pseudo_pdf.columns]
            if combined_pseudo_pdf is None:
                combined_pseudo_pdf = pseudo_pdf
            else:
                # Reindex the dataframes to match the combined dataframe's index (cluster names) and fill missing values with 0
                pseudo_pdf = pseudo_pdf.reindex(index=combined_pseudo_pdf.index, columns=pseudo_pdf.columns).fillna(0)
                # Concatenate the dataframes along columns
                combined_pseudo_pdf = pd.concat([combined_pseudo_pdf, pseudo_pdf], axis=1)

    pseudo_counts = collection_name+'_combined_pseudo_cluster_counts.txt'
    combined_pseudo_pdf.to_csv(pseudo_counts, sep='\t')


def pseudo_cluster_counts_optimized(sample, combined_adata, cell_threshold=5, count_threshold=5, compute_tpm=False, tpm_threshold=1, status=True):
    dataType = 'isoform' if compute_tpm else 'junction'
    
    # Prepare the output file path
    if dataType == 'isoform':
        output_file = f"{sample}.txt"
    else:
        output_file = f"{sample}.txt"
    # Write the header first
    with open(output_file, 'w') as f_out:
        f_out.write('Feature\t' + '\t'.join(combined_adata.obs['cluster'].unique()) + '\n')
    
    chunk_size = 1000  # Number of rows to process per chunk
    num_chunks = int(np.ceil(combined_adata.shape[1] / chunk_size))
    
    # Prepare TPM and isoform ratio file paths if necessary
    if compute_tpm:
        tpm_output_file = f"{sample}_tpm.txt"
        isoform_ratio_file = f"{sample}_ratio.txt"
        with open(tpm_output_file, 'w') as tpm_out, open(isoform_ratio_file, 'w') as isoform_out:
            tpm_out.write('Feature\t' + '\t'.join(combined_adata.obs['cluster'].unique()) + '\n')
            isoform_out.write('Feature\t' + '\t'.join(combined_adata.obs['cluster'].unique()) + '\n')

    for chunk_idx in tqdm(range(num_chunks), desc="Processing chunks of data"):
        start = chunk_idx * chunk_size
        end = min((chunk_idx + 1) * chunk_size, combined_adata.shape[1])
        
        # Extract the chunk of data
        chunk_data = combined_adata[:, start:end].X.toarray()
        chunk_var_names = combined_adata.var_names[start:end]

        # Convert the chunk to a DataFrame for processing
        chunk_df = pd.DataFrame(chunk_data, index=combined_adata.obs_names, columns=chunk_var_names)
        
        # Group by clusters
        grouped = chunk_df.groupby(combined_adata.obs['cluster'], observed=True)
        
        # Filter clusters and sum values
        filtered_summed_groups = grouped.filter(lambda x: len(x) >= cell_threshold)
        summed_groups = filtered_summed_groups.groupby(combined_adata.obs['cluster'], observed=True).sum()

        # Filter out low counts
        filtered_summed_groups_transposed = summed_groups.loc[:, summed_groups.sum(axis=0) >= count_threshold].transpose()
        
        # Append the chunk results to the output file
        filtered_summed_groups_transposed.to_csv(output_file, sep='\t', mode='a', header=False)

        if compute_tpm:
            # Calculate TPM values
            tpm_values = summed_groups.div(summed_groups.sum(axis=1), axis=0) * 1e6
            tpm_values_transposed = tpm_values.transpose()
            tpm_values_transposed.to_csv(tpm_output_file, sep='\t', mode='a', header=False)

            # Compute gene sums and isoform ratios
            gene_sums = summed_groups.T.groupby(lambda x: x.split(":")[0]).transform('sum').T
            isoform_ratios = summed_groups.div(gene_sums)

            # Filter isoform ratios where gene TPM < 1
            gene_tpms = tpm_values.T.groupby(lambda x: x.split(":")[0]).sum().T
            mask = gene_tpms >= tpm_threshold
            mask = mask.reindex(columns=tpm_values.columns, method='ffill')
            isoform_ratios = isoform_ratios.where(mask).transpose()
            
            isoform_ratios.to_csv(isoform_ratio_file, sep='\t', mode='a', header=False)

    if compute_tpm:
        print(f'Finished processing and writing {output_file}, {tpm_output_file}, and {isoform_ratio_file}')
        return output_file, tpm_output_file, isoform_ratio_file
    else:
        print(f'Finished processing and writing {output_file}')
        return output_file

def concatenate_h5ad_and_compute_pseudobulks_optimized(sample_files, collection_name='junction', compute_tpm=False, tpm_threshold=1):
    combined_pseudo_pdf = None
    combined_tpm = None
    combined_isoform_to_gene_ratio = None
    
    # Process each sample file
    for sample_file in sample_files:
        sample_name = sample_file.split('.')[0]  # Extract sample name
        print(f'Processing {sample_file}')

        # Load the data
        adata = ad.read_h5ad(sample_file)

        # Ensure indices are strings and observation names are unique
        adata.obs.index = adata.obs.index.astype(str)
        adata.obs_names_make_unique()

        # Generate pseudobulks for the current sample
        if compute_tpm:
            pseudo_pdf_file, tpm_file, isoform_ratio_file = pseudo_cluster_counts_optimized(
                sample_name, adata, cell_threshold=5, count_threshold=0, compute_tpm=compute_tpm, tpm_threshold=tpm_threshold, status=False
            )
        else:
            pseudo_pdf_file = pseudo_cluster_counts_optimized(
                sample_name, adata, cell_threshold=5, count_threshold=0, compute_tpm=compute_tpm, tpm_threshold=tpm_threshold, status=False
            )

        # Load the pseudobulk file
        pseudo_pdf = pd.read_csv(pseudo_pdf_file, sep='\t', index_col=0)

        # Rename the columns to include the sample name
        pseudo_pdf.columns = [f'{col}.{sample_name}' for col in pseudo_pdf.columns]

        # Combine the pseudobulk counts
        if combined_pseudo_pdf is None:
            combined_pseudo_pdf = pseudo_pdf
        else:
            combined_pseudo_pdf = pd.concat([combined_pseudo_pdf, pseudo_pdf], axis=1).fillna(0)

        if compute_tpm:
            # Load the TPM and isoform ratio files
            tpm = pd.read_csv(tpm_file, sep='\t', index_col=0)
            isoform_to_gene_ratio = pd.read_csv(isoform_ratio_file, sep='\t', index_col=0)

            # Rename the columns to include the sample name
            tpm.columns = [f'{col}.{sample_name}' for col in tpm.columns]
            isoform_to_gene_ratio.columns = [f'{col}.{sample_name}' for col in isoform_to_gene_ratio.columns]

            # Combine the TPM and isoform-to-gene ratio dataframes
            if combined_tpm is None:
                combined_tpm = tpm
                combined_isoform_to_gene_ratio = isoform_to_gene_ratio
            else:
                combined_tpm = pd.concat([combined_tpm, tpm], axis=1).fillna(0)
                combined_isoform_to_gene_ratio = pd.concat([combined_isoform_to_gene_ratio, isoform_to_gene_ratio], axis=1).fillna(0)

    # Save the combined pseudobulk data to file
    combined_pseudo_file = f'{collection_name}_combined_pseudo_cluster_counts.txt'
    combined_pseudo_pdf.to_csv(combined_pseudo_file, sep='\t')

    if compute_tpm:
        # Save the combined TPM and isoform ratio data to files
        combined_tpm_file = f'{collection_name}_combined_pseudo_cluster_tpm.txt'
        combined_isoform_ratio_file = f'{collection_name}_combined_pseudo_cluster_ratio.txt'
        combined_tpm.to_csv(combined_tpm_file, sep='\t')
        combined_isoform_to_gene_ratio.to_csv(combined_isoform_ratio_file, sep='\t')

        return combined_pseudo_file, combined_tpm_file, combined_isoform_ratio_file

    return combined_pseudo_file

