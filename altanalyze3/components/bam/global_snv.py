import os,sys,string,h5py
import pysam
import multiprocessing
from collections import defaultdict
import argparse
import time
import pandas as pd
from pathlib import Path
import scipy.sparse as sp
from scipy.sparse import lil_matrix
import numpy as np

""" Note: code was partially implemented using ChatGPT prompts """

def analyze_chromosome(bam_file_path, chromosome, reference_genome_path, output_file, min_reads, min_percent):

    # Adjust chromosome name for BAM file (add 'chr' prefix if not present)
    try:
        coordinates = None
        bam_chromosome = chromosome if chromosome.startswith('chr') else 'chr' + chromosome
    except: 
        coordinates = chromosome[1]
        chromosome = chromosome[0] # for supplied genomic positions
        bam_chromosome = chromosome if chromosome.startswith('chr') else 'chr' + chromosome

    bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    ref_genome = pysam.FastaFile(reference_genome_path)
    mismatch_counts = defaultdict(lambda: {'total': 0, 'ref': 0, 'alt': 0})

    for position, counts in mismatch_counts.items():
        # Check if the position is in the specific_positions set
        if specific_positions and (chromosome, position) not in specific_positions:
            continue

    if chromosome == "chrM" or chromosome == "M":
        chromosome = "MT"

    if coordinates == None:
        for read in bam_file.fetch(bam_chromosome):  # use bam_chromosome here
            if read.is_unmapped:
                continue

            ref_seq = ref_genome.fetch(chromosome, read.reference_start, read.reference_end)  # use original chromosome name
                
            #if not read.is_secondary and not read.is_supplementary:
            for read_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
                if read_pos is not None and ref_pos is not None:
                    read_base = read.query_sequence[read_pos]
                    ref_base = ref_seq[ref_pos - read.reference_start]

                    mismatch_counts[ref_pos]['total'] += 1
                    if read_base == ref_base:
                        mismatch_counts[ref_pos]['ref'] += 1
                    else:
                        mismatch_counts[ref_pos]['alt'] += 1

        filtered_positions=[]
        with open(output_file, 'a') as f:
            for position, counts in mismatch_counts.items():
                if counts['total'] >= min_reads and counts['ref'] / counts['total'] >= min_percent / 100 and counts['alt'] / counts['total'] >= min_percent / 100:
                    f.write(f"{chromosome}\t{position}\t{counts['total']}\t{counts['ref']}\t{counts['alt']}\n")
                    filtered_positions.append(position)
    else:
        filtered_positions = set(coordinates)

    reference_dict = {}
    #print([bam_chromosome,chromosome],filtered_positions)
    for position in filtered_positions:
        for pileupcolumn in bam_file.pileup(bam_chromosome, position, position + 1, truncate=True):
            if pileupcolumn.pos == position:
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        try:
                            barcode = pileupread.alignment.get_tag('CB')
                            nt = pileupread.alignment.query_sequence[pileupread.query_position]
                            ref_seq = ref_genome.fetch(chromosome, pileupread.alignment.reference_start, pileupread.alignment.reference_end)
                            ref_base = ref_seq[position - pileupread.alignment.reference_start]
                            allele_type = "ref" if nt == ref_base else "alt"
                            key = f"{chromosome}:{position}|{nt}-{allele_type}"
                            if key not in reference_dict:
                                reference_dict[key] = {}
                            if barcode not in reference_dict[key]:
                                reference_dict[key][barcode] = 0
                            reference_dict[key][barcode] += 1
                        except Exception:
                            continue  
    bam_file.close()

    # Convert the reference_dict to a DataFrame and save it
    df_ref = pd.DataFrame.from_dict(reference_dict, orient='index').fillna(0).astype(int)
    ref_output_file_path = str(Path(output_file).parent) +'/chr-results/'+ f'alleles_{chromosome}.txt'
    try:
        df_ref.to_csv(ref_output_file_path, sep='\t')
    except:
        os.mkdir(str(Path(output_file).parent) +'/chr-results/')
        df_ref.to_csv(ref_output_file_path, sep='\t')


""" Functions to combine the chromosome specific SNV and export to text or h5 """

def combine_files(directory_path):
    directory = Path(directory_path)
    all_files = directory.glob('chr-results/*.txt')
    barcode_lists = []
    variant_lists = []
    sparse_matrices = []

    for file_path in all_files:
        barcodes, variants, sparse_matrix = read_and_sparse(file_path)
        barcode_lists.append(barcodes)
        variant_lists.append(variants)
        sparse_matrices.append(sparse_matrix)

    all_barcodes, all_variants, combined_matrix = combine_sparse_matrices(barcode_lists, variant_lists, sparse_matrices)
    print(f"Combined matrix shape: {combined_matrix.shape}")

    save_combined_matrix_h5(combined_matrix, all_barcodes, all_variants, directory / 'combined_matrix.h5')
    save_combined_matrix_tsv(combined_matrix, all_barcodes, all_variants, directory / 'combined_matrix.tsv')

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

def save_combined_matrix_tsv(combined_matrix, barcodes, variants, file_path):
    # Convert the combined sparse matrix to a dense format and save as TSV
    dense_matrix = combined_matrix.toarray()
    df_combined = pd.DataFrame(dense_matrix, index=variants, columns=barcodes)
    df_combined.to_csv(file_path, sep='\t')

def save_combined_matrix_h5(sparse_matrix, barcodes, variants, file_path):
    with h5py.File(file_path, 'w') as f:
        # Store the sparse matrix components with compression
        f.create_dataset('data', data=sparse_matrix.data, compression='gzip')
        f.create_dataset('indices', data=sparse_matrix.indices, compression='gzip')
        f.create_dataset('indptr', data=sparse_matrix.indptr, compression='gzip')
        f.create_dataset('shape', data=sparse_matrix.shape)
        # Convert barcodes and variants to strings and store
        barcodes_str = np.array(barcodes, dtype='S')
        variants_str = np.array(variants, dtype='S')
        f.create_dataset('barcodes', data=barcodes_str, compression='gzip')
        f.create_dataset('variants', data=variants_str, compression='gzip')

""" Multi-processing function """

def parallel_process(bam_file_path, reference_genome_path, output_file, min_reads, 
                        min_percent, test_chromosome=None, positions_file=None):

    initiate = time.time()
    bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    bam_chromosomes = set(bam_file.references)
    bam_file.close()

    ref_genome = pysam.FastaFile(reference_genome_path)
    ref_chromosomes = set(ref_genome.references)
    ref_genome.close()

    # Find matching chromosomes considering different naming conventions
    if test_chromosome:
        chromosomes = [test_chromosome]
    else:
        chromosomes = []
        for chrom in bam_chromosomes:
            if chrom in ref_chromosomes or ('chr' + chrom) in ref_chromosomes or chrom.replace('chr', '', 1) in ref_chromosomes:
                chromosomes.append(chrom.replace('chr',''))
            if chrom == 'chrM' or chrom == 'chrMT':
                chromosomes.append('M')

    with open(output_file, 'w') as f:
        f.write("Chromosome\tPosition\tTotal_Reads\tReference_Reads\tAlternative_Reads\n")

    processes = []
    max_concurrent_processes = 27

    # Read positions from file if provided
    if positions_file:
        chromosome_db = {}
        chromosomes=[]
        with open(positions_file, 'r') as file:
            for line in file:
                chrom, pos = line.strip().split(':')
                try: chromosome_db[chrom].append(int(pos))
                except: chromosome_db[chrom] = [int(pos)]
        for chrom in chromosome_db:
            chromosomes.append([chrom,chromosome_db[chrom]])

    for i, chromosome in enumerate(chromosomes):
        p = multiprocessing.Process(target=analyze_chromosome, args=(bam_file_path, chromosome, reference_genome_path, output_file, min_reads, min_percent))
        processes.append(p)
        p.start()

        # Wait for some processes to complete before starting more
        if i % max_concurrent_processes == 0:
            for p in processes:
                p.join()
            processes = []

    # Wait for any remaining processes to complete
    for p in processes:
        p.join()
    elapsed = time.time() - initiate

    combine_files(Path(output_file).parent)
    
    print(f"Processing completed in {elapsed:.2f} seconds.")

def main():

    # python3 isoform_snv.py MultiLin-2__RUNX1.5801-pre.bam SNV/genome.fa out/out.txt --test_chromosome chrM --min_reads 200 --min_percent 10
    # python3 isoform_snv.py MultiLin-2__RUNX1.5801-pre.bam SNV/genome.fa out/out.txt --min_reads 200 --min_percent 10
    parser = argparse.ArgumentParser(description='Analyze BAM file for mismatches.')
    parser.add_argument('bam_file', help='Path to the input BAM file')
    parser.add_argument('reference_genome', help='Path to the reference genome FASTA file')
    parser.add_argument('output_file', help='Path to the output file')
    parser.add_argument('--test_chromosome', help='Specify a single chromosome for testing', default=None)
    parser.add_argument('--min_reads', help='Minimum number of reads at a position', type=int, default=250)
    parser.add_argument('--min_percent', help='Minimum percentage of reference/alternative reads', type=float, default=10.0)
    parser.add_argument('--positions_file', help='Path to a file with specific genomic positions (formatted as chromosome:position)', default=None)

    args = parser.parse_args()

    parallel_process(args.bam_file, args.reference_genome, args.output_file, args.min_reads, args.min_percent, args.test_chromosome, args.positions_file)

if __name__ == "__main__":
    main()