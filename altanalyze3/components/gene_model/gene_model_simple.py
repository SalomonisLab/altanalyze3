import pandas as pd
import numpy as np
import os
import sys
import re
import argparse
from tqdm import tqdm
import multiprocessing as mp

def build_gene_model(df):
    transcript_id = -1
    position_footprint = {}
    position_ense = {}
    gene_start, gene_end, strand, ensg, chromosome = None, None, None, None, None
    pat = re.compile(r'exon_id=(ENSE\d+\.\d+)')


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
                for p in range(row.start, row.end + 1):
                    position_footprint.setdefault(p, []).append(transcript_id)
                    position_ense.setdefault(p, []).append(ense)

    if gene_start is None:
        return ""

    if strand == '+':
        string_stream = ''
        running_profile = position_footprint[gene_start]
        block_index = 1
        intron_block_index = 0
        segment_index = 1
        anchor_position = gene_start
        p = gene_start
        while p <= gene_end:
            try:
                position_footprint[p]
            except KeyError:
                subexon_identifier = 'E'+str(block_index)+'.'+str(segment_index)
                associated_ense = '|'.join(position_ense[p-1]) if (p-1) in position_ense else ''
                string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ensg, chromosome, strand, subexon_identifier, anchor_position, p-1, associated_ense)
                pi = p
                while True:
                    try:
                        position_footprint[pi]
                    except KeyError:
                        pi += 1
                    else:
                        intron_block_index += 1
                        subexon_identifier = 'I'+str(intron_block_index)+'.'+str(1)
                        string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t\n'.format(ensg, chromosome, strand, subexon_identifier, p, pi-1)
                        p = pi
                        running_profile = position_footprint[p]
                        block_index += 1
                        segment_index = 1
                        anchor_position = p
                        break
            else:
                if position_footprint[p] != running_profile:
                    subexon_identifier = 'E'+str(block_index)+'.'+str(segment_index)
                    associated_ense = '|'.join(position_ense[p-1]) if (p-1) in position_ense else ''
                    string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ensg, chromosome, strand, subexon_identifier, anchor_position, p-1, associated_ense)
                    running_profile = position_footprint[p]
                    segment_index += 1
                    anchor_position = p
                    p += 1
                else:
                    p += 1
        subexon_identifier = 'E'+str(block_index)+'.'+str(segment_index)
        associated_ense = '|'.join(position_ense[p-1]) if (p-1) in position_ense else ''
        string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ensg, chromosome, strand, subexon_identifier, anchor_position, p-1, associated_ense)
    else:
        string_stream = ''
        running_profile = position_footprint[gene_end]
        block_index = 1
        intron_block_index = 0
        segment_index = 1
        anchor_position = gene_end
        p = gene_end
        while p >= gene_start:
            try:
                position_footprint[p]
            except KeyError:
                subexon_identifier = 'E'+str(block_index)+'.'+str(segment_index)
                associated_ense = '|'.join(position_ense[p+1]) if (p+1) in position_ense else ''
                string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ensg, chromosome, strand, subexon_identifier, p+1, anchor_position, associated_ense)
                pi = p
                while True:
                    try:
                        position_footprint[pi]
                    except KeyError:
                        pi -= 1
                    else:
                        intron_block_index += 1
                        subexon_identifier = 'I'+str(intron_block_index)+'.'+str(1)
                        string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t\n'.format(ensg, chromosome, strand, subexon_identifier, pi+1, p)
                        p = pi
                        running_profile = position_footprint[p]
                        block_index += 1
                        segment_index = 1
                        anchor_position = p
                        break
            else:
                if position_footprint[p] != running_profile:
                    subexon_identifier = 'E'+str(block_index)+'.'+str(segment_index)
                    associated_ense = '|'.join(position_ense[p+1]) if (p+1) in position_ense else ''
                    string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ensg, chromosome, strand, subexon_identifier, p+1, anchor_position, associated_ense)
                    running_profile = position_footprint[p]
                    segment_index += 1
                    anchor_position = p
                    p -= 1
                else:
                    p -= 1
        subexon_identifier = 'E'+str(block_index)+'.'+str(segment_index)
        associated_ense = '|'.join(position_ense[p+1]) if (p+1) in position_ense else ''
        string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ensg, chromosome, strand, subexon_identifier, p+1, anchor_position, associated_ense)
    return string_stream

def main(args):
    gtf = args.gtf
    outdir = args.outdir
    target_gene = args.gene

    os.makedirs(outdir, exist_ok=True)
    fout_path = os.path.join(outdir, f'gene_model_{target_gene}.txt')
    fout = open(fout_path, 'w')

    total_size = os.path.getsize(gtf)
    processed_bytes = 0
    update_every_n_lines = 100

    pbar = tqdm(total=total_size, unit='B', unit_scale=True, desc="Processing GTF/GFF3")

    current_gene_id = None
    current_gene_rows = []
    line_counter = 0

    with open(gtf, 'r') as f:
        for line in f:
            processed_bytes += len(line)
            line_counter += 1

            if line.startswith("#") or line.strip() == '':
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, source, feature_type, start, end, score, strand, phase, attrs = fields

            if feature_type not in {'gene', 'transcript', 'exon'}:
                continue

            if gtf.endswith('.gff3') or gtf.endswith('.gff'):
                pat = re.compile(r'gene_id=(ENSG\d+)(?:\.\d+)?')
            else:
                pat = re.compile(r'gene_id "(ENSG\d+)"')

            match = re.search(pat, attrs)
            if not match:
                continue
            gene_id = match.group(1)

            row_data = {
                'chr': chrom.replace('chr',''),
                'source': source,
                'type': feature_type,
                'start': int(start),
                'end': int(end),
                'score': score,
                'strand': strand,
                'phase': phase,
                'attrs': attrs,
                'gene': gene_id
            }

            if target_gene != 'all' and gene_id != target_gene:
                continue

            if current_gene_id is None:
                current_gene_id = gene_id

            if gene_id != current_gene_id:
                df = pd.DataFrame(current_gene_rows)
                string_stream = build_gene_model(df)
                fout.write(string_stream)
                current_gene_rows = []
                current_gene_id = gene_id

            current_gene_rows.append(row_data)

            if line_counter % update_every_n_lines == 0:
                pbar.update(processed_bytes - pbar.n)

    if current_gene_rows:
        df = pd.DataFrame(current_gene_rows)
        string_stream = build_gene_model(df)
        fout.write(string_stream)

    if pbar.n < processed_bytes:
        pbar.update(processed_bytes - pbar.n)
    pbar.close()
    fout.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build Gene Model from GTF/GFF3')
    parser.add_argument('--gtf', type=str, required=True, help='Input GTF or GFF3 file')
    parser.add_argument('--gene', type=str, required=True, help='Single ENSG id or "all"')
    parser.add_argument('--outdir', type=str, required=True, help='Output directory')
    args = parser.parse_args()
    main(args)
