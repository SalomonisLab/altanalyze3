import os, sys, time
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
from anndata import AnnData, concat, read_h5ad


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


def h5ad_to_adata(h5ad_path, rev=False):
    """Import isoform h5ad files and match mtx_to_adata conventions."""
    adata = read_h5ad(h5ad_path)
    def ensure_barcode_suffix(barcode):
        return barcode if '-' in barcode else f"{barcode}-1"

    def reverse_complement_barcode(barcode):
        if '-' in barcode:
            return reverse_complement_seq(barcode)
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        rev_comp = ''.join(complement.get(base, 'N') for base in reversed(barcode))
        return f"{rev_comp}-1"

    barcodes = adata.obs_names.to_series().astype(str)
    if rev:
        barcodes = barcodes.apply(reverse_complement_barcode)
    else:
        barcodes = barcodes.apply(ensure_barcode_suffix)
    adata.obs_names = barcodes.values
    adata.var_names = (
        adata.var_names.to_series()
        .astype(str)
        .apply(lambda value: value.split(':', 1)[1] if ':' in value else value)
        .values
    )
    return adata


def matrix_dir_to_adata(int_folder, gene_is_index, feature, feature_col,
                        barcode, barcode_col, matrix, rev=False):
    """Dispatch to h5ad or mtx based on the matrix_dir name."""
    path = str(int_folder)
    lower_path = path.lower()
    if lower_path.endswith('.h5ad') or lower_path.endswith('.h5ad.gz'):
        return h5ad_to_adata(path, rev=rev)
    return mtx_to_adata(
        path,
        gene_is_index,
        feature,
        feature_col,
        barcode,
        barcode_col,
        matrix,
        rev=rev
    )

def bulk_counts_to_adata(file_path: str):
    """Convert a text file with either (1) pbid and counts to an AnnData object."""
    df = pd.read_csv(file_path, sep='\t', comment='#')
    
    if df.shape[1] < 2:
        raise ValueError("The TSV file must contain at least two columns.")

    if 'pbid' in df.columns:
        # Handle the original format with 'annot_transcript_id' column
        transcript_ids = df['pbid']
        counts = df.iloc[:, -2].values  # Assuming counts are in the second to last column
        count_df = pd.DataFrame(counts, index=transcript_ids, columns=['counts'])
        count_df = count_df.transpose()

    adata = ad.AnnData(X=count_df)
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
    adata.obs.index.name = None
    adata.obs = adata.obs.join(barcode_clusters, how='inner')

    return adata

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


def pseudo_cluster_counts(sample, combined_adata, cell_threshold=0, count_threshold=0, compute_tpm=False, tpm_threshold=1, status=True):
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

def export_and_filter_pseudobulks(input_file, output_file, cell_type_order=None, min_group_size = 3):
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
                    if any(count >= min_group_size for count in cell_type_counts.values()):
                        if cell_type_order:
                            reordered_values = [fields[samples.index(col) + 1] for col in reordered_columns]
                        else:
                            reordered_values = fields[1:]
                        outfile.write(f"{feature}\t" + "\t".join(map(str, reordered_values)) + "\n")

                    pbar.update(1)

def export_and_filter_pseudobulk_chunks(input_file, output_file, cell_type_order=None, min_group_size=3, chunk_size=500000):
    # Read header and determine sample names
    with open(input_file, 'r') as infile:
        header = infile.readline().strip().split('\t')
        samples = header[1:]  # Skip the first column (feature names)

    # Precompute cell type mappings
    cell_types = np.array([sample.split('.')[0] for sample in samples])
    unique_cell_types = np.unique(cell_types)

    # If reordering is needed, create index mapping
    if cell_type_order:
        reordered_indices = [i for cell_type in cell_type_order for i, sample in enumerate(samples) if sample.startswith(cell_type + '.')]
    else:
        reordered_indices = np.arange(len(samples))  # No reordering needed

    # Open output file and write header
    with open(output_file, 'w') as outfile:
        outfile.write('\t'.join(['Feature'] + [samples[i] for i in reordered_indices]) + '\n')

        # Read in chunks instead of line-by-line processing
        chunk_iter = pd.read_csv(input_file, sep='\t', skiprows=1, header=None, names=['Feature'] + samples, chunksize=chunk_size)

        for chunk in tqdm(chunk_iter, desc="Processing Chunks"):
            # Convert values to NumPy array for speed
            values = chunk.iloc[:, 1:].values.astype(float)  # Convert data to float fast

            # Compute nonzero counts per cell type using NumPy fast operations
            cell_type_counts = {ct: np.sum(values[:, cell_types == ct] > 0.1, axis=1) for ct in unique_cell_types}

            # Identify rows meeting the min_group_size threshold
            mask = np.zeros(values.shape[0], dtype=bool)
            for ct in unique_cell_types:
                mask |= (cell_type_counts[ct] >= min_group_size)

            # Filter rows
            filtered_chunk = chunk[mask]

            # Reorder columns if needed
            if cell_type_order:
                filtered_chunk = filtered_chunk[['Feature'] + [samples[i] for i in reordered_indices]]

            # Write filtered data in one batch
            filtered_chunk.to_csv(outfile, sep='\t', index=False, header=False, mode='a')



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
            pseudo_pdf, tpm, isoform_to_gene_ratio = pseudo_cluster_counts(sample_name, adata, cell_threshold=0, count_threshold=0, compute_tpm=True, status=False)
        
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
            pseudo_pdf = pseudo_cluster_counts(sample_name, adata, cell_threshold=0, count_threshold=0, compute_tpm=False, status=False)
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

def pseudo_cluster_counts_optimized(sample, combined_adata, cell_threshold=0, count_threshold=0, compute_tpm=False, tpm_threshold=1, status=True):
    output_file = f"{sample}.txt"
    cluster_labels = combined_adata.obs['cluster'].astype(str)
    cluster_order = sorted(cluster_labels.unique().tolist())

    with open(output_file, 'w') as f_out:
        f_out.write('Feature\t' + '\t'.join(cluster_order) + '\n')

    if compute_tpm:
        tpm_output_file = f"{sample}_tpm.txt"
        isoform_ratio_file = f"{sample}_ratio.txt"
        with open(tpm_output_file, 'w') as tpm_out, open(isoform_ratio_file, 'w') as isoform_out:
            header = 'Feature\t' + '\t'.join(cluster_order) + '\n'
            tpm_out.write(header)
            isoform_out.write(header)

    chunk_size = 1000
    num_chunks = int(np.ceil(combined_adata.shape[1] / chunk_size))

    for chunk_idx in tqdm(range(num_chunks), desc="Processing chunks of data"):
        start = chunk_idx * chunk_size
        end = min((chunk_idx + 1) * chunk_size, combined_adata.shape[1])
        chunk_data = combined_adata[:, start:end].X.toarray()
        chunk_var_names = combined_adata.var_names[start:end]

        chunk_df = pd.DataFrame(chunk_data, index=combined_adata.obs_names, columns=chunk_var_names)
        chunk_df = chunk_df.join(cluster_labels.rename("cluster"))

        grouped = chunk_df.groupby("cluster", observed=True)
        valid_groups = grouped.filter(lambda x: len(x) >= cell_threshold)
        summed_groups = valid_groups.groupby("cluster", observed=True).sum()

        result = summed_groups.T.reindex(columns=cluster_order).fillna(0)
        result = result.loc[result.sum(axis=1) >= count_threshold]
        result.insert(0, "Feature", result.index)
        result.to_csv(output_file, sep="\t", index=False, mode="a", header=False)

        if compute_tpm:
            tpm = summed_groups.div(summed_groups.sum(axis=1), axis=0) * 1e6
            tpm_values = tpm.T.reindex(columns=cluster_order).fillna(0)
            tpm_values.insert(0, "Feature", tpm_values.index)
            tpm_values.to_csv(tpm_output_file, sep="\t", index=False, mode="a", header=False)

            isoform_cols = summed_groups.columns
            gene_ids = [col.split(":")[0] for col in isoform_cols]

            # TPM threshold mask (isoform level)
            mask = tpm >= tpm_threshold

            # Gene-level raw counts
            gene_counts = summed_groups.T.groupby(lambda x: x.split(":")[0]).sum().T

            
            gene_counts_broadcasted = pd.DataFrame(
                [[gene_counts.loc[idx, gid] if gid in gene_counts.columns else np.nan for gid in gene_ids] for idx in summed_groups.index],
                index=summed_groups.index,
                columns=isoform_cols
            )

            # Replace 0 with np.nan in denominator
            safe_denominator = gene_counts_broadcasted.replace(0, np.nan)

            # Calculate isoform-to-gene ratio
            ratios = summed_groups / safe_denominator
            ratios = ratios.where(mask)

            # Optional diagnostics
            if status:
                for gid in set(gene_ids):
                    if gene_ids.count(gid) > 1 and gid in gene_counts.columns:
                        isoform_subset = [col for col in isoform_cols if col.startswith(gid + ":")]
                        print(f"\n🔎 QC for gene: {gid}")
                        print("\nRaw counts:\n", summed_groups[isoform_subset])
                        print("\nTPM values:\n", tpm[isoform_subset])
                        print("\nTPM threshold mask:\n", mask[isoform_subset])
                        print("\nGene counts:\n", gene_counts[gid])
                        print("\nBroadcasted gene counts:\n", gene_counts_broadcasted[isoform_subset])
                        print("\nCalculated ratios:\n", ratios[isoform_subset])
                        break

            ratios = ratios.T.reindex(columns=cluster_order)
            ratios.insert(0, "Feature", ratios.index)
            ratios.to_csv(isoform_ratio_file, sep="\t", index=False, mode="a", header=False)

            #print(f"Final isoform_ratios non-zero entries: {(ratios.iloc[:, 1:] > 0).sum().sum()}")

    print(f"Done writing {output_file}" + (f", {tpm_output_file}, {isoform_ratio_file}" if compute_tpm else ''))

    if compute_tpm:
        return output_file, tpm_output_file, isoform_ratio_file
    else:
        return output_file

def just_return_dense_count_files(sample, compute_tpm=False):
    # Bypass the above function to regenerate the counts and TPM files

    output_file = f"{sample}.txt"
    if compute_tpm:
        tpm_output_file = f"{sample}_tpm.txt"
        isoform_ratio_file = f"{sample}_ratio.txt"
        return output_file, tpm_output_file, isoform_ratio_file
    else:
        return output_file

def concatenate_h5ad_and_compute_pseudobulks_optimized(sample_files, collection_name='junction', compute_tpm=False, tpm_threshold=1):
    
    combined_pseudo_pdf = None
    combined_tpm = None
    combined_isoform_to_gene_ratio = None
    only_import_existing = False

    # Process each sample file
    for sample_file in sample_files:
        sample_name = sample_file.split('.')[0]  # Extract sample name
        print(f'Processing {sample_file}')

        if only_import_existing:
            pseudo_pdf_file = just_return_dense_count_files(sample_name, compute_tpm=False)
            if compute_tpm:
                pseudo_pdf_file, tpm_file, isoform_ratio_file = pseudo_pdf_file
        else:
            # Load the data
            adata = ad.read_h5ad(sample_file)

            # Ensure indices are strings and observation names are unique
            adata.obs.index = adata.obs.index.astype(str)
            adata.obs_names_make_unique()

            # Generate pseudobulks for the current sample
            if compute_tpm:
                pseudo_pdf_file, tpm_file, isoform_ratio_file = pseudo_cluster_counts_optimized(
                    sample_name, adata, cell_threshold=0, count_threshold=0, compute_tpm=compute_tpm, tpm_threshold=tpm_threshold, status=False
                )
            else:
                pseudo_pdf_file = pseudo_cluster_counts_optimized(
                    sample_name, adata, cell_threshold=0, count_threshold=0, compute_tpm=compute_tpm, tpm_threshold=tpm_threshold, status=False
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


def extract_all_features_above_threshold(h5ad_files, min_reads):
    feature_sums = {}  # {feature: total_count_across_files}

    for file in h5ad_files:
        adata = read_h5ad(file)
        X = adata.X
        var_names = adata.var_names

        # Sum reads for this file
        counts = np.asarray(X.sum(axis=0)).ravel()
        for feat, count in zip(var_names, counts):
            if feat in feature_sums:
                feature_sums[feat] += count
            else:
                feature_sums[feat] = count

    # Keep only features with total reads above threshold
    retained_features = [feat for feat, total in feature_sums.items() if total >= min_reads]
    print(f"Retained {len(retained_features)} features with ≥{min_reads} total reads")

    return sorted(retained_features)

def align_adata_to_union(h5ad_file, all_features, feature_to_idx):
    adata = read_h5ad(h5ad_file)
    sample_id = os.path.basename(h5ad_file).replace('.h5ad', '')
    # Ensure unique obs_names across all files
    adata.obs_names = [f"{sample_id}.{name}" for name in adata.obs_names]

    old_features = adata.var_names

    old_to_new_idx = {i: feature_to_idx[feat] for i, feat in enumerate(old_features) if feat in feature_to_idx}

    X = adata.X.tocoo()
    rows, cols, data = [], [], []
    for i, j, v in zip(X.row, X.col, X.data):
        if j in old_to_new_idx:
            rows.append(i)
            cols.append(old_to_new_idx[j])
            data.append(v)
    new_X = csr_matrix((data, (rows, cols)), shape=(adata.n_obs, len(all_features)))

    new_var = pd.DataFrame(index=all_features)
    for col in adata.var.columns:
        dtype = float if pd.api.types.is_numeric_dtype(adata.var[col]) else object
        new_var[col] = pd.Series(index=all_features, dtype=dtype)
        shared_features = adata.var.index.intersection(new_var.index)
        new_var.loc[shared_features, col] = adata.var.loc[shared_features, col]

    adata.obs = adata.obs.copy()
    adata.obs['sample'] = os.path.basename(h5ad_file).replace('.h5ad', '')
    return AnnData(X=new_X, obs=adata.obs, var=new_var, uns=adata.uns.copy())

def combine_h5ad_files_parallel(sample_h5ads_dict, output_file='combined.h5ad', n_jobs=4, min_total_reads=200):
    start_time = time.time()
    h5ad_files = sorted(set(sample_h5ads_dict.values()))
    print(f"Identifying union of features across {len(h5ad_files)} files...")
    #all_features = extract_all_features(h5ad_files)
    all_features = extract_all_features_above_threshold(h5ad_files, min_total_reads)

    feature_to_idx = {feat: idx for idx, feat in enumerate(all_features)}

    print("Aligning features in parallel...")
    aligned_adatas = Parallel(n_jobs=n_jobs)(
        delayed(align_adata_to_union)(file, all_features, feature_to_idx)
        for file in h5ad_files
    )

    print("Concatenating aligned AnnData objects...")
    combined_adata = concat(aligned_adatas, axis=0, join='outer')
    print(f"Writing combined file to: {output_file}")
    combined_adata.write_h5ad(output_file, compression='gzip')
    end_time = time.time()  # End the timer
    print(f"Execution time: {end_time - start_time:.2f} seconds")
    return combined_adata

if __name__ == '__main__':
    sample_name = 'test'
    compute_tpm = True
    tpm_threshold = 1

    isoform_input = "WM34-isoform.h5ad"

    adata = ad.read_h5ad(isoform_input)

    # Ensure indices are strings and observation names are unique
    adata.obs.index = adata.obs.index.astype(str)
    adata.obs_names_make_unique()

    target_gene_id = "ENSG00000177674"
    isoform_vars = [v for v in adata.var_names if v.startswith(target_gene_id)]

    # ----- Extract counts matrix for selected isoforms -----
    X = adata[:, isoform_vars].copy()

    pseudo_pdf_file = pseudo_cluster_counts_optimized(
        sample_name, X, cell_threshold=0, count_threshold=0, compute_tpm=compute_tpm, 
        tpm_threshold=tpm_threshold, status=True)
