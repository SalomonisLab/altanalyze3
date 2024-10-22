import os,sys
import pandas as pd
import numpy as np
from copy import deepcopy
from tqdm import tqdm  # For the progress bar
import time
import multiprocessing
import argparse

# Original method for clique finding
def find_uid_in_clique(the_uid, region, strand, uid2coords):
    clique = {the_uid: region}
    region_start, region_end = map(int, region)

    for uid, coords in uid2coords.items():
        if uid != the_uid:
            coords_start, coords_end = map(int, coords)
            # Check overlap
            if strand == '+':
                if max(region_start, coords_start) - min(region_end, coords_end) < 0:
                    clique[uid] = coords
            elif strand == '-':
                if max(region_end, coords_end) - min(region_start, coords_start) < 0:
                    clique[uid] = coords

    return clique


# Early exit conditional checking to optimize
def is_valid(mat, uid):
    num_incl_events = np.count_nonzero(mat[0, :] >= 20)
    if num_incl_events <= 1:
        return False  # Early exit if inclusion events are too low

    num_excl_events = np.count_nonzero(mat[1:, :].sum(axis=0) >= 20)
    if num_excl_events <= 1:
        return False  # Early exit if exclusion events are too low

    total_number_junctions = np.count_nonzero(mat.sum(axis=0) >= 20)
    if total_number_junctions <= 2:
        return False  # Early exit if total junction count is too low

    return True


# Core PSI calculation with early exit conditions
def calculate_psi_core(clique, uid, count, sample_columns):
    sub_count = count.loc[list(clique.keys()), sample_columns]
    mat = sub_count.values
    with np.errstate(divide='ignore', invalid='ignore'):
        psi = mat[0, :] / mat.sum(axis=0)

    if sub_count.shape[0] > 1:
        bg_uid = sub_count.index.tolist()[np.argmax(mat[1:, :].sum(axis=1)) + 1]
        cond = is_valid(mat, uid)
    else:
        bg_uid = 'None'
        cond = False

    return psi, bg_uid, cond


# PSI calculation for a chunk of genes
def calculate_psi_per_gene_chunk(chunk_data):
    count, header, sample_columns, lookup_table = chunk_data
    result_data = []
    uid2coords = count.apply(lambda x: [x['start'], x['end']], axis=1, result_type='reduce').to_dict()

    for row in count.iterrows():
        uid = row[0]  # This should be the correct splice junction ID
        region = uid2coords[uid]
        strand = row[1]['strand']
        clique = find_uid_in_clique(uid, region, strand, uid2coords)
        psi_values, bg_uid, cond = calculate_psi_core(clique, uid, count, sample_columns)

        if np.nanmax(psi_values) - np.nanmin(psi_values) >= 0.1 and cond:
            # Ensure the `uid` splice junction IDs are retained by looking up in the lookup table
            uid_str = lookup_table.get(uid, uid)  # Use lookup table to get original string
            bg_uid_str = lookup_table.get(bg_uid, bg_uid)  # Map bg_uid to string
            data = (uid_str + '|' + bg_uid_str, *psi_values)
            result_data.append(data)

    if result_data:
        df = pd.DataFrame.from_records(result_data, columns=['uid'] + sample_columns)
        return df
    else:
        return None


# Write results to file synchronously (with a single write step after parallel processing)
def write_to_file(outdir, results, write_header=True):
    with open(outdir, 'a') as f:
        for df in results:
            if df is not None:
                df.to_csv(f, sep='\t', index=False, header=write_header)
                write_header = False  # Ensure the header is only written once


# Modify the function signature to accept a multiprocessing context
def process_junctions_in_chunks(junction_path, query_gene, outdir, total_lines, mp_context=None, num_cores=None):
    col_uid = []
    col_gene = []
    col_chrom = []
    col_start = []
    col_end = []
    col_strand = []
    data = []
    current_gene = None
    gene_chunks = []  # Store gene data chunks
    sample_columns = None
    lookup_table = {}  # Dictionary to map index to string IDs

    with tqdm(total=total_lines, desc="Processing junctions") as pbar:
        with open(junction_path, 'r') as f:
            header = next(f).strip().split('\t')
            sample_columns = header[1:]

            for line in f:
                pbar.update(1)
                row = line.strip().split('\t')
                uid, coords = row[0].split('=')  # This should be the correct splice junction ID
                chrom, coords = coords.split(':')
                start, end = coords.split('-')
                strand = '+' if start < end else '-'
                gene = uid.split(':')[0]

                # Store the mapping from index to original string ID in the lookup table
                index = len(lookup_table)  # Assuming index is the current length of the lookup table
                lookup_table[index] = uid  # Map the index to the string ID

                # If processing a specific gene, skip others
                if query_gene and gene != query_gene:
                    continue

                if current_gene and current_gene != gene:
                    count = pd.DataFrame(data, columns=header[1:])
                    for name, col in zip(['uid', 'gene', 'chrom', 'start', 'end', 'strand'],
                                         [col_uid, col_gene, col_chrom, col_start, col_end, col_strand]):
                        count[name] = col

                    # Collect chunks for multiprocessing
                    gene_chunks.append((count, header, sample_columns, lookup_table))

                    # Reset collections for the next gene
                    col_uid, col_gene, col_chrom, col_start, col_end, col_strand, data = [], [], [], [], [], [], []

                current_gene = gene
                col_uid.append(uid)  # Properly track splice junction IDs
                col_gene.append(gene)
                col_chrom.append(chrom)
                col_start.append(start)
                col_end.append(end)
                col_strand.append(strand)
                data.append([float(x) if x.replace('.', '', 1).isdigit() else np.nan for x in row[1:]])

            if data:
                count = pd.DataFrame(data, columns=header[1:])
                for name, col in zip(['uid', 'gene', 'chrom', 'start', 'end', 'strand'],
                                     [col_uid, col_gene, col_chrom, col_start, col_end, col_strand]):
                    count[name] = col
                gene_chunks.append((count, header, sample_columns, lookup_table))

    # Use the provided multiprocessing context to process gene chunks
    if mp_context and num_cores:
        with mp_context.Pool(num_cores) as pool:
            results = pool.map(calculate_psi_per_gene_chunk, gene_chunks)
    else:
        results = [calculate_psi_per_gene_chunk(chunk) for chunk in gene_chunks]

    # Write results in a single pass
    write_to_file(outdir, results, write_header=True)


def main(junction_path=None, query_gene=None, outdir=None, use_multiprocessing=True, mp_context=None, num_cores=None):
    print(f"Computing PSI from: {junction_path}")
    if junction_path is None or outdir is None:
        raise ValueError("Required parameters: junction_path and outdir must be provided.")

    # Count total lines for progress bar
    total_lines = sum(1 for line in open(junction_path))

    # Initialize the output file (overwrite if it exists)
    with open(outdir, 'w') as f:
        pass  # Simply open and close the file to overwrite/clear it

    # Process junctions in chunks with or without multiprocessing
    start_time = time.time()

    process_junctions_in_chunks(junction_path, query_gene, outdir, total_lines, mp_context, num_cores)

    end_time = time.time()
    print(f"Total time: {end_time - start_time:.2f} seconds")


# Protect multiprocessing entry point
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute PSI from count junction')
    parser.add_argument('--junction', type=str, default=None, help='the path to the junction count')
    parser.add_argument('--gene', type=str, default=None, help='gene you want to compute PSI')
    parser.add_argument('--outdir', type=str, default=None, help='output dir for the output file')
    args = parser.parse_args()

    # Create a multiprocessing context
    from multiprocessing import get_context, cpu_count
    mp_context = get_context('spawn')

    # Call main with parsed arguments and multiprocessing enabled
    main(junction_path=args.junction, query_gene=args.gene, outdir=args.outdir, use_multiprocessing=True, mp_context=mp_context, num_cores=cpu_count())
