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
    # GTF exon IDs are quoted; GFF3 uses key=value and custom IDs are valid.
    pat = re.compile(r'(?:exon_id[ =]+|ID=)[\"]?([^;\" ]+)')
    if not df.empty:
        gene_start, gene_end = int(df.start.min()), int(df.end.max())
        strand, ensg, chromosome = df.iloc[0].strand, df.iloc[0].gene, df.iloc[0].chr

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
            ense = ense_match.group(1) if ense_match else f'{row.gene}:{row.start}-{row.end}'
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
    df['type'] = df['type'].replace({'mRNA': 'transcript'})
    df = df[df['type'].isin(['gene', 'transcript', 'exon'])].copy()

    def attributes(text):
        fields = {}
        for field in text.split(';'):
            field = field.strip()
            if not field:
                continue
            key, value = field.split('=', 1) if '=' in field else field.split(None, 1)
            fields[key] = value.strip('" ')
        return fields
    parsed = [attributes(x) for x in df['attrs']]
    parent_gene = {}
    for row, attrs in zip(df.itertuples(index=False), parsed):
        if row.type == 'gene':
            ident = attrs.get('gene_id') or attrs.get('ID')
            if ident:
                parent_gene[attrs.get('ID', ident)] = ident.removeprefix('gene:')
    for row, attrs in zip(df.itertuples(index=False), parsed):
        if row.type == 'transcript':
            parent = attrs.get('gene_id') or parent_gene.get(attrs.get('Parent', ''), attrs.get('Parent', ''))
            if parent:
                parent_gene[attrs.get('ID', attrs.get('transcript_id', ''))] = parent.removeprefix('gene:')
    df['gene'] = [a.get('gene_id') or parent_gene.get(a.get('ID', '')) or parent_gene.get(a.get('Parent', '')) for a in parsed]
    if df.loc[df['type'] == 'exon', 'gene'].isna().any():
        raise ValueError('Exons require gene_id or resolvable GFF3 Parent attributes')

    if gene != 'all':
        df = df[df['gene'] == gene]
        string_stream = build_gene_model(df)
    else:
        sub_df_list = [group for _, group in df.groupby(by='gene')] 
        cores = max(1, min(getattr(args, 'cpus', 1), len(sub_df_list)))
        chunks = split_array_to_chunks(sub_df_list, cores)
        if cores == 1:
            string_stream = process_single_core(chunks[0])
        else:
            with mp.Pool(processes=cores) as pool:
                string_stream = ''.join(pool.map(process_single_core, chunks))
        if not string_stream:
            raise ValueError('Reference contains no usable exon models')

    with open(os.path.join(outdir, f'gene_model_{gene}.tsv'), 'w') as f:
        f.write(string_stream)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build Gene Model from GTF')
    parser.add_argument('--gtf', type=str, required=True, help='the path to the gtf file')
    parser.add_argument('--gene', type=str, required=True, help='either all or stable ENSG ID')
    parser.add_argument('--outdir', type=str, required=True, help='output dir for the gene model txt file')
    args = parser.parse_args()
    main(args)
