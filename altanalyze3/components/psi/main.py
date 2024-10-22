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

def calculate_psi_per_gene_chunk(chunk_data):
    count, header, sample_columns, lookup_table = chunk_data
    result_data = []

    print("Processing chunk...")
    print("Header:", header)
    print("Sample columns:", sample_columns)
    print("First few rows of count DataFrame:\n", count.head())
    print("Index in count DataFrame:", count.index)

    # Create a mapping from UIDs to their coordinates
    try:
        uid2coords = count.apply(lambda x: [x['start'], x['end']], axis=1, result_type='reduce').to_dict()
    except Exception as e:
        print(f"Error creating UID to coordinate mapping: {str(e)}")
        sys.exit(1)

    # Iterate through each row in the DataFrame
    for uid, row in count.iterrows():
        print(f"Processing UID: {uid}")

        if uid not in uid2coords:
            print(f"Warning: Region not found for UID {uid}. Skipping.")
            continue

        try:
            region = uid2coords[uid]
            strand = row['strand']
        except KeyError as e:
            print(f"Error accessing region/strand for UID {uid}: {str(e)}")
            sys.exit(1)

        # Find the clique for this UID
        clique = find_uid_in_clique(uid, region, strand, uid2coords)
        print(f"Clique for UID {uid}: {clique}")

        # Use sub_count to filter by clique
        try:
            sub_count = count.loc[list(clique.keys()), sample_columns]
            print(f"sub_count DataFrame for UID {uid}:\n", sub_count.head())
        except Exception as e:
            print(f"Error selecting rows for UID {uid}: {str(e)}")
            sys.exit(1)

        # Calculate PSI values
        try:
            psi_values, bg_uid, cond = calculate_psi_core(clique, uid, sub_count, sample_columns)
            print(f"PSI values for {uid}: {psi_values}, Background UID: {bg_uid}, Condition: {cond}")
        except Exception as e:
            print(f"Error during PSI calculation for {uid}: {str(e)}")
            sys.exit(1)

        if np.nanmax(psi_values) - np.nanmin(psi_values) >= 0.1 and cond:
            try:
                uid_str = lookup_table.get(uid, uid)
                bg_uid_str = lookup_table.get(bg_uid, bg_uid)
                data = (uid_str + '|' + bg_uid_str, *psi_values)
                result_data.append(data)
                print(f"Added result: {data}")
            except Exception as e:
                print(f"Error while adding result for UID {uid}: {str(e)}")
                sys.exit(1)

    if result_data:
        df = pd.DataFrame.from_records(result_data, columns=['uid'] + sample_columns)
        print("Chunk processed successfully. Returning DataFrame...")
        return df
    else:
        print("No valid results found for this chunk.")
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
    col_uid, col_gene, col_chrom, col_start, col_end, col_strand = [], [], [], [], [], []
    data = []
    current_gene = None
    gene_chunks = []  # Store gene data chunks
    sample_columns = None
    lookup_table = {}  # Direct UID to UID mapping

    with tqdm(total=total_lines, desc="Processing junctions") as pbar:
        with open(junction_path, 'r') as f:
            header = next(f).strip().split('\t')
            sample_columns = header[1:]

            for line in f:
                pbar.update(1)
                row = line.strip().split('\t')
                uid, coords = row[0].split('=')
                chrom, coords = coords.split(':')
                start, end = coords.split('-')
                strand = '+' if int(start) < int(end) else '-'
                gene = uid.split(':')[0]

                # Store uid directly in the lookup table
                lookup_table[uid] = uid

                if query_gene and gene != query_gene:
                    continue

                if current_gene and current_gene != gene:
                    count = pd.DataFrame(data, columns=header[1:])
                    count.set_index('uid', inplace=True)
                    count['uid'] = col_uid
                    for name, col in zip(['uid', 'gene', 'chrom', 'start', 'end', 'strand'],
                                         [col_uid, col_gene, col_chrom, col_start, col_end, col_strand]):
                        count[name] = col

                    count.set_index('uid', inplace=True)

                    # Debugging print to confirm indexing worked
                    print("Index in count DataFrame:", count.index[:5])  # Optional: Remove after debugging


                    gene_chunks.append((count, header, sample_columns, lookup_table))
                    col_uid, col_gene, col_chrom, col_start, col_end, col_strand, data = [], [], [], [], [], [], []

                current_gene = gene
                col_uid.append(uid)
                col_gene.append(gene)
                col_chrom.append(chrom)
                col_start.append(start)
                col_end.append(end)
                col_strand.append(strand)
                data.append([float(x) if x.replace('.', '', 1).isdigit() else np.nan for x in row[1:]])

            # Ensure 'uid' becomes the index to align with expected DataFrame usage
            if data:
                # Create sample DataFrame from the collected data
                sample_df = pd.DataFrame(data, columns=sample_columns)

                # Create metadata DataFrame
                metadata_df = pd.DataFrame({
                    'uid': col_uid,
                    'gene': col_gene,
                    'chrom': col_chrom,
                    'start': col_start,
                    'end': col_end,
                    'strand': col_strand
                })

                # Ensure the lengths of the data match
                if len(metadata_df) != len(sample_df):
                    raise ValueError(f"Mismatch: metadata={len(metadata_df)}, sample={len(sample_df)}")

                # Merge the metadata and sample data along the columns
                count = pd.concat([metadata_df, sample_df], axis=1)

                # Set 'uid' as the index to enable proper lookups
                count.set_index('uid', inplace=True)  # This is critical

                # Add this chunk to the list of gene chunks for further processing
                gene_chunks.append((count, header, sample_columns, lookup_table))


                return gene_chunks

def main(junction_path=None, query_gene=None, outdir=None, use_multiprocessing=False, mp_context=None, num_cores=None):
    print(f"Computing PSI from: {junction_path}")
    if junction_path is None or outdir is None:
        raise ValueError("Required parameters: junction_path and outdir must be provided.")

    total_lines = sum(1 for _ in open(junction_path))

    with open(outdir, 'w') as f:
        pass

    gene_chunks = process_junctions_in_chunks(junction_path, query_gene, outdir, total_lines, mp_context, num_cores)

    results = [calculate_psi_per_gene_chunk(chunk) for chunk in gene_chunks]

    write_to_file(outdir, results, write_header=True)



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
    main(junction_path=args.junction, query_gene=args.gene, outdir=args.outdir, use_multiprocessing=False)
