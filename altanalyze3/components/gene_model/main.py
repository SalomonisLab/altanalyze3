import pandas as pd
import numpy as np
import os, sys
import re
from tqdm import tqdm
import argparse
import multiprocessing as mp

def build_gene_model(df):
    transcript_id = -1
    position_footprint = {}
    position_ense = {}
    boundary_set = set()
    gene_start, gene_end, strand, ensg, chromosome = None,None,None,None,None
    pat = re.compile(r'exon_id=(ENSE\d+\.\d+)')

    # phase1: build footprint and record all unique exon boundary positions
    for row in df.itertuples(index=False):
        if row.type == 'gene':
            gene_start, gene_end, strand, ensg, chromosome = row.start, row.end, row.strand, row.gene, row.chr
        elif row.type == 'transcript':
            transcript_id += 1
            continue
        elif row.type == 'exon':
            attrs = row.attrs
            ense_match = re.search(pat, attrs)
            if ense_match:
                ense = ense_match.group(1)
                for p in range(row.start, row.end):
                    position_footprint.setdefault(p, []).append(transcript_id)
                    position_ense.setdefault(p, []).append(ense)
                boundary_set.add(row.start)
                boundary_set.add(row.end)

    # phase2: sort all boundaries into ordered list
    ordered_boundaries = sorted(boundary_set, reverse=(strand == '-'))
    string_stream = ''
    block_index = 1
    segment_index = 1

    for i in range(len(ordered_boundaries) - 1):
        s = ordered_boundaries[i]
        e = ordered_boundaries[i + 1]
        start, end = min(s, e), max(s, e) # Corrects for order in negative strand positions
        mid_point = (start + end) // 2
        current_profile = position_footprint.get(mid_point)
        
        if current_profile is None:
            # This is an intron
            original_start, original_end = start, end
            if strand == '+':
                start += 1
                end -= 1
            else:
                start += 1
                end -= 1

            # If adjustment creates an invalid interval, fall back to original
            if end <= start:
                start, end = original_start, original_end
            subexon_identifier = f'I{block_index}.1'
            string_stream += f'{ensg}\t{subexon_identifier}\t{chromosome}\t{strand}\t{start}\t{end}\t\t\n'
            block_index += 1
            segment_index = 1
        else:
            # This is an exon segment
            ense_ids = []
            for p in range(start, end):
                ense_ids.extend(position_ense.get(p, []))
            ense_ids = list(dict.fromkeys(ense_ids))
            associated_ense = '|'.join(ense_ids)
            subexon_identifier = f'E{block_index}.{segment_index}'
            string_stream += f'{ensg}\t{subexon_identifier}\t{chromosome}\t{strand}\t{start}\t{end}\t\t{associated_ense}\n'
            segment_index += 1


    return string_stream



def split_array_to_chunks(array, cores=None):
    if not isinstance(array, list):
        raise Exception('split_array_to_chunks function works for list, not ndarray')
    array_index = np.arange(len(array))
    if cores is None:
        cores = mp.cpu_count()
    sub_indices = np.array_split(array_index, cores)
    return [[array[i] for i in sub_index] for sub_index in sub_indices]


def process_single_core(chunk):
    string_stream = ''
    for sub_df in tqdm(chunk, total=len(chunk)):
        string_stream += build_gene_model(sub_df)
    return string_stream


def main(args):
    gtf = args.gtf
    gene = args.gene
    outdir = args.outdir

    df = pd.read_csv(gtf, sep='\t', comment='#', header=None)
    df.columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attrs']
    df = df[df['type'].isin(['gene', 'transcript', 'exon'])]

    if gtf.endswith('.gff3') or gtf.endswith('.gff'):
        pat = re.compile(r'gene_id=(ENSG\d+)(?:\.\d+)?')
    else:
        pat = re.compile(r'gene_id "(ENSG\d+)"')
    df['gene'] = df['attrs'].apply(lambda x: re.search(pat, x).group(1) if re.search(pat, x) else None)

    if gene != 'all':
        df = df[df['gene'] == gene]
        string_stream = build_gene_model(df)
    else:
        sub_df_list = [group for _, group in df.groupby(by='gene')] 
        chunks = split_array_to_chunks(sub_df_list, mp.cpu_count())
        pool = mp.Pool(processes=mp.cpu_count())
        print(f'spawn {mp.cpu_count()} subprocesses')
        r = [pool.apply_async(func=process_single_core, args=(chunk,)) for chunk in chunks]
        pool.close()
        pool.join()
        string_stream = ''.join([res.get() for res in r])

    with open(os.path.join(outdir, f'gene_model_{gene}.tsv'), 'w') as f:
        f.write(string_stream)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build Gene Model from GTF')
    parser.add_argument('--gtf', type=str, required=True, help='the path to the gtf file')
    parser.add_argument('--gene', type=str, required=True, help='either all or stable ENSG ID')
    parser.add_argument('--outdir', type=str, required=True, help='output dir for the gene model txt file')
    args = parser.parse_args()
    main(args)
