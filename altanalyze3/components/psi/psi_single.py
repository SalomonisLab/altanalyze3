import os,sys
import pandas as pd
import numpy as np
from copy import deepcopy
import argparse
import warnings
import aiofiles  # Asynchronous file I/O
import asyncio   # For async handling
from tqdm import tqdm  # For the progress bar

# Chunk size and batch size
CHUNK_SIZE = 50  # Number of genes to process before writing
BATCH_SIZE = 10  # Number of genes per async file write

def initialize_output_file(file_path, header):
    with open(file_path, 'w') as f:
        f.write('\t'.join(header) + '\n')  # Write header once at the beginning

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
    min_reads = 20
    num_incl_events = np.count_nonzero(mat[0, :] >= min_reads)
    if num_incl_events <= 1:
        return False  # Early exit if inclusion events are too low

    num_excl_events = np.count_nonzero(mat[1:, :].sum(axis=0) >= min_reads)
    if num_excl_events <= 1:
        return False  # Early exit if exclusion events are too low

    total_number_junctions = np.count_nonzero(mat.sum(axis=0) >= min_reads)
    if total_number_junctions < 2:
        return False  # Early exit if total junction count is too low

    return True

# Core PSI calculation with early exit conditions
def calculate_psi_core(clique, uid, count, sample_columns):
    sub_count = count.loc[list(clique.keys()), sample_columns]
    mat = sub_count.values
    with np.errstate(divide='ignore', invalid='ignore'):
        psi = mat[0, :] / mat.sum(axis=0)
        psi[mat.sum(axis=0) < 5] = np.nan # Sets PSI to null where less than 5 reads
    if sub_count.shape[0] > 1:
        bg_uid = sub_count.index.tolist()[np.argmax(mat[1:, :].sum(axis=1)) + 1]
        cond = is_valid(mat, uid)
    else:
        bg_uid = 'None'
        cond = False

    return psi, bg_uid, cond

"""
### For debugging
def calculate_psi_core(clique, uid, count, sample_columns):

    TARGET_JUNC = "ENSG00000112062:E22.14-E26.1"

    sub_count = count.loc[list(clique.keys()), sample_columns]

    # enforce target junction as numerator (never return None)
    numerator_junc = sub_count.index[0]
    if numerator_junc != TARGET_JUNC:
        return np.full(len(sample_columns), np.nan), None, False

    mat = sub_count.values

    # ---- DEBUG: first sample only ----
    j = 0  # first sample encountered

    numerator = mat[0, j]
    denominator = mat[:, j].sum()

    print("PSI DEBUG (TARGET EVENT)")
    print("  Event UID:", uid)
    print("  Sample:", sample_columns[j])
    print("  Inclusion junction:", numerator_junc)
    print("  Exclusion junctions:", sub_count.index[1:].tolist())
    print("  Inclusion reads (numerator):", numerator)
    print("  Total reads (denominator):", denominator)
    print("  Per-junction counts:")
    for k, junc in enumerate(sub_count.index):
        print("   ", junc, "=", mat[k, j])
    print("-" * 60)

    with np.errstate(divide='ignore', invalid='ignore'):
        denom = mat.sum(axis=0)
        psi = mat[0, :] / denom
        psi[denom < 5] = np.nan  # existing rule

    if sub_count.shape[0] > 1:
        bg_uid = sub_count.index.tolist()[np.argmax(mat[1:, :].sum(axis=1)) + 1]
        cond = is_valid(mat, uid)
    else:
        bg_uid = 'None'
        cond = False

    return psi, bg_uid, cond
"""


# Asynchronously write to file
async def write_to_file(file_path, data, write_header):
    async with aiofiles.open(file_path, 'a') as f:
        header = '\t'.join(data.columns) + '\n' if write_header else ''
        await f.write(header)
        await f.write(data.to_csv(sep='\t', index=False, header=False))

# PSI calculation for a chunk of genes, writing in batches
def calculate_psi_per_gene(count, outdir, write_header=True, result_batch=None):
    return_data = []
    count = count.set_index('uid')
    uid2coords = count.apply(lambda x: [x['start'], x['end']], axis=1, result_type='reduce').to_dict()

    warnings.filterwarnings("ignore", category=RuntimeWarning, message="All-NaN slice encountered")
    for row in count.iterrows():
        uid = row[0]
        region = uid2coords[uid]
        strand = row[1]['strand']
        clique = find_uid_in_clique(uid, region, strand, uid2coords)
        index_list = deepcopy(row[1].index.tolist())
        for c in ['gene', 'chrom', 'start', 'end', 'strand']:
            index_list.remove(c)
        sample_columns = index_list
        psi_values, bg_uid, cond = calculate_psi_core(clique, uid, count, sample_columns)

        # Only keep PSI results that meet the range condition
        if np.nanmax(psi_values) - np.nanmin(psi_values) >= 0.1 and cond:
            data = (uid+'|'+bg_uid, *psi_values)
            return_data.append(data)

    if return_data:
        df = pd.DataFrame.from_records(return_data, columns=['uid-bg_uid'] + sample_columns)
        result_batch.append(df)


async def main(junction_path=None, query_gene=None, outdir=None):
    import time
    start_time = time.time()

    # Extract the header before initializing the output file
    with open(junction_path, 'r') as f:
        header = ['uid']+next(f).strip().split('\t')  # Read and split header from input file

    # Initialize output file with header
    initialize_output_file(outdir, header[1:])  # Call the function here

    col_uid = []
    col_gene = []
    col_chrom = []
    col_start = []
    col_end = []
    col_strand = []
    data = []

    current_gene = None
    write_header = True
    result_batch = []

    total_lines = sum(1 for line in open(junction_path))  # Count lines for progress bar
    with tqdm(total=total_lines, desc="Processing junctions") as pbar:
        with open(junction_path, 'r') as f:
            header = next(f).strip().split('\t')
            for line in f:
                pbar.update(1)  # Update progress bar
                row = line.strip().split('\t')
                uid, coords = row[0].split('=')
                chrom, coords = coords.split(':')
                start, end = coords.split('-')
                strand = '+' if start < end else '-'
                gene = uid.split(':')[0]

                # If processing a specific gene, skip others
                if query_gene and gene != query_gene:
                    continue

                if current_gene and current_gene != gene:
                    count = pd.DataFrame(data, columns=header[1:])
                    for name, col in zip(['uid', 'gene', 'chrom', 'start', 'end', 'strand'],
                                         [col_uid, col_gene, col_chrom, col_start, col_end, col_strand]):
                        count[name] = col

                    # Process chunk of genes and write to disk
                    calculate_psi_per_gene(count, outdir, write_header, result_batch)
                    if len(result_batch) >= BATCH_SIZE:
                        # Write results in batches
                        await write_to_file(outdir, pd.concat(result_batch), write_header)
                        result_batch.clear()

                    write_header = False

                    col_uid, col_gene, col_chrom, col_start, col_end, col_strand, data = [], [], [], [], [], [], []
                
                current_gene = gene
                col_uid.append(uid)
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
            calculate_psi_per_gene(count, outdir, write_header, result_batch)
            count={}

        if result_batch:
            await write_to_file(outdir, pd.concat(result_batch), write_header)
            
    end_time = time.time()
    print(f"Total time: {end_time - start_time:.2f} seconds")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Compute PSI from count junction')
    parser.add_argument('--junction', type=str, default=None, help='the path to the junction count')
    parser.add_argument('--gene', type=str, default=None, help='gene you want to compute PSI')
    parser.add_argument('--outdir', type=str, default=None, help='output dir for the output file')
    args = parser.parse_args()

    junction_path = args.junction
    query_gene = args.gene
    outdir = args.outdir

    asyncio.run(main(junction_path=junction_path, query_gene=None, outdir=outdir))
