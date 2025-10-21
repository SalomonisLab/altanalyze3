import pandas as pd
import scipy.sparse as sp
from scipy.sparse import lil_matrix
import os, sys

def read_and_sparse(file_path):
    df = pd.read_csv(file_path, sep='\t', index_col=0)
    sparse_matrix = sp.csr_matrix(df.values)  # Convert the DataFrame to a sparse matrix
    return df.columns, df.index, sparse_matrix

def combine_sparse_matrices(barcode_lists, variant_lists, sparse_matrices):
    # Find the union of all barcodes and variants
    all_barcodes = sorted(set.union(*map(set, barcode_lists)))
    all_variants = sorted(set.union(*map(set, variant_lists)))

    # Create a dictionary to map barcodes and variants to indices
    barcode_to_index = {barcode: i for i, barcode in enumerate(all_barcodes)}
    variant_to_index = {variant: i for i, variant in enumerate(all_variants)}

    # Initialize a combined sparse matrix
    combined_matrix = lil_matrix((len(all_variants), len(all_barcodes)), dtype=int)

    # Fill in the combined matrix
    for barcodes, variants, sparse_matrix in zip(barcode_lists, variant_lists, sparse_matrices):
        barcodes_indices = [barcode_to_index[barcode] for barcode in barcodes]
        for variant_idx, variant in enumerate(variants):
            variant_index = variant_to_index[variant]
            # Extract the row corresponding to the variant from the sparse matrix
            variant_data = sparse_matrix[variant_idx, :].toarray().ravel()
            # Update the combined matrix
            combined_matrix[variant_index, barcodes_indices] = variant_data

    return all_barcodes, all_variants, combined_matrix.tocsr()

def import_txt_files(directory_path):
    directory_abs_path = os.path.abspath(directory_path)
    directory_items = os.listdir(directory_abs_path)
    txt_files = [item for item in directory_items if item.endswith('.txt')]
    full_paths = [os.path.join(directory_abs_path, file) for file in txt_files]
    return full_paths

if __name__ == "__main__":
    directory_path = '/Volumes/salomonis-archive/FASTQs/Grimes/RNA/scRNASeq/10X-Genomics/2-16-23PabioRevio/Sample34-2_FebRun_outputs/mutations/chr-results/'
    file_paths = import_txt_files(directory_path)
    barcode_lists = []
    variant_lists = []
    sparse_matrices = []

    for file_path in file_paths:
        barcodes, variants, sparse_matrix = read_and_sparse(file_path)
        barcode_lists.append(barcodes)
        variant_lists.append(variants)
        sparse_matrices.append(sparse_matrix)

    all_barcodes, all_variants, combined_matrix = combine_sparse_matrices(barcode_lists, variant_lists, sparse_matrices)
    print(f"Combined matrix shape: {combined_matrix.shape}")
