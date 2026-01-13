#!/usr/bin/env python3
import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pysam
from scipy.sparse import coo_matrix

sys.path.insert(1, os.path.join(os.path.dirname(__file__), '..'))
from long_read import gff_process


def parse_gene_regions(exon_file):
    gene_regions = {}
    with open(exon_file, 'r') as handle:
        first = handle.readline()
        if not first:
            return gene_regions
        header = first.rstrip('\n').split('\t')
        has_header = header and header[0] == 'gene'
        if has_header:
            header_map = {name: idx for idx, name in enumerate(header)}
            idx_gene = header_map.get('gene', 0)
            idx_chr = header_map.get('chromosome', 2)
            idx_strand = header_map.get('strand', 3)
            idx_start = header_map.get('exon-region-start(s)', 4)
            idx_end = header_map.get('exon-region-stop(s)', 5)
        else:
            idx_gene, idx_chr, idx_strand, idx_start, idx_end = 0, 2, 3, 4, 5
            _update_gene_region(header, gene_regions, idx_gene, idx_chr, idx_strand, idx_start, idx_end)
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            _update_gene_region(parts, gene_regions, idx_gene, idx_chr, idx_strand, idx_start, idx_end)
    return gene_regions


def _update_gene_region(parts, gene_regions, idx_gene, idx_chr, idx_strand, idx_start, idx_end):
    if len(parts) <= max(idx_gene, idx_chr, idx_strand, idx_start, idx_end):
        return
    gene = parts[idx_gene].strip()
    chrom = parts[idx_chr].strip()
    strand = parts[idx_strand].strip()
    try:
        start = int(float(parts[idx_start]))
        end = int(float(parts[idx_end]))
    except ValueError:
        return
    if not gene or not chrom:
        return
    start, end = (start, end) if start <= end else (end, start)
    if gene in gene_regions:
        existing = gene_regions[gene]
        gene_regions[gene] = (
            existing[0],
            existing[1],
            min(existing[2], start),
            max(existing[3], end),
        )
    else:
        gene_regions[gene] = (chrom, strand, start, end)


def normalize_chrom(chrom, exon_coordinates):
    if not exon_coordinates:
        return chrom
    sample_chr = next(iter(exon_coordinates.keys()))[0]
    if sample_chr.startswith('chr'):
        return chrom if chrom.startswith('chr') else f'chr{chrom}'
    return chrom.replace('chr', '', 1) if chrom.startswith('chr') else chrom


def resolve_bam_chrom(chrom, bam_refs):
    if chrom in bam_refs:
        return chrom
    if chrom.startswith('chr'):
        alt = chrom.replace('chr', '', 1)
    else:
        alt = f'chr{chrom}'
    if alt in bam_refs:
        return alt
    return None


def get_read_strand(read):
    for tag in ('ts', 'XS', 'TS'):
        try:
            strand = read.get_tag(tag)
        except KeyError:
            continue
        if strand in ('+', '-'):
            return strand
    return '-' if read.is_reverse else '+'


def has_splice(cigartuples):
    if not cigartuples:
        return False
    return any(op == 3 for op, _ in cigartuples)


def extract_exons(read):
    exons = []
    for start, end in read.get_blocks():
        exon_start = start + 1
        exon_end = end
        if exon_end >= exon_start:
            exons.append((exon_start, exon_end))
    return exons


def has_known_splice_site(chrom, strand, exons, exon_coordinates):
    if len(exons) < 2:
        return False
    for idx in range(len(exons) - 1):
        donor = exons[idx][1]
        acceptor = exons[idx + 1][0]
        if ((chrom, donor, strand, 1) in exon_coordinates or
                (chrom, donor, strand, 2) in exon_coordinates or
                (chrom, acceptor, strand, 1) in exon_coordinates or
                (chrom, acceptor, strand, 2) in exon_coordinates):
            return True
    return False


def resolve_barcode(read, barcode_tags):
    for tag in barcode_tags:
        try:
            value = read.get_tag(tag)
        except KeyError:
            continue
        if value:
            return value
    return None


def resolve_molecule_id(read, molecule_tag):
    name = read.query_name or ''
    if molecule_tag:
        try:
            value = read.get_tag(molecule_tag)
        except KeyError:
            value = None
        if value is not None:
            tag_value = str(value)
            if name and name.endswith(tag_value) and '/' in name:
                return name
            return tag_value
    return name or 'unknown'


def sanitize_gff_value(value):
    if value is None:
        return ''
    return str(value).replace('"', '').replace(';', '').strip().replace(' ', '_')


def get_index(value, index_map, items):
    if value in index_map:
        return index_map[value]
    idx = len(items)
    index_map[value] = idx
    items.append(value)
    return idx


def write_gff_isoform(handle, chrom, strand, exons, isoform_id, gene, barcode, molecule_id, source):
    info_parts = []
    if gene:
        info_parts.append(f'gene_id "{gene}"')
    info_parts.append(f'transcript_id "{isoform_id}"')
    if molecule_id:
        info_parts.append(f'molecule_id "{molecule_id}"')
    if barcode:
        info_parts.append(f'cell_barcode "{barcode}"')
    info = ';'.join(info_parts) + ';'
    for exon_start, exon_end in exons:
        handle.write(
            f"{chrom}\t{source}\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t{info}\n"
        )
    tx_start = min(exon[0] for exon in exons)
    tx_end = max(exon[1] for exon in exons)
    handle.write(
        f"{chrom}\t{source}\ttranscript\t{tx_start}\t{tx_end}\t.\t{strand}\t.\t{info}\n"
    )


def extract_isoform_structures(bam_path, exon_file, output_prefix, target_gene=None,
                               min_mapq=1, barcode_tags=None, molecule_tag='zm',
                               source='bam', require_known_splice=True):
    if barcode_tags is None:
        barcode_tags = ['CB', 'CR', 'BC', 'BX']

    exon_coordinates, _, _ = gff_process.importEnsemblGenes(exon_file)
    gene_regions = parse_gene_regions(exon_file) if target_gene else {}

    bam = pysam.AlignmentFile(bam_path, 'rb')
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    gff_path = output_prefix.with_suffix('.gff')
    h5ad_path = output_prefix.with_suffix('.h5ad')

    rows = []
    cols = []
    data = []
    barcodes = []
    isoforms = []
    barcode_index = {}
    isoform_index = {}
    isoform_gene = {}
    stats = defaultdict(int)

    if target_gene:
        if target_gene not in gene_regions:
            raise ValueError(f"Gene {target_gene} not found in {exon_file}")
        gene_chr, _, gene_start, gene_end = gene_regions[target_gene]
        bam_chr = resolve_bam_chrom(gene_chr, bam.references)
        if bam_chr is None:
            raise ValueError(f"Chromosome {gene_chr} not found in BAM references")
        fetch_iter = bam.fetch(bam_chr, gene_start - 1, gene_end)
    else:
        fetch_iter = bam.fetch(until_eof=True)

    with open(gff_path, 'w') as gff_handle:
        for read in fetch_iter:
            stats['total_reads'] += 1
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            stats['mapped_primary_reads'] += 1
            if read.mapping_quality < min_mapq:
                continue
            if not has_splice(read.cigartuples):
                continue
            stats['spliced_reads'] += 1
            barcode = resolve_barcode(read, barcode_tags)
            if not barcode:
                continue
            stats['barcode_reads'] += 1

            exons = extract_exons(read)
            if len(exons) < 2:
                continue
            chrom = normalize_chrom(read.reference_name, exon_coordinates)
            strand = get_read_strand(read)
            gene, _, _, genes = gff_process.exonAnnotate(
                chrom, list(exons), strand, read.query_name
            )
            if not genes:
                continue
            if target_gene and target_gene not in genes and gene != target_gene:
                continue
            if target_gene and target_gene in genes:
                gene = target_gene
            stats['gene_assigned_spliced_reads'] += 1

            has_known = has_known_splice_site(chrom, strand, exons, exon_coordinates)
            if not has_known and require_known_splice:
                continue
            if has_known:
                stats['known_splice_reads'] += 1

            molecule_id = sanitize_gff_value(resolve_molecule_id(read, molecule_tag))
            isoform_id = sanitize_gff_value(f"{gene}:{molecule_id}")
            write_gff_isoform(
                gff_handle,
                chrom,
                strand,
                exons,
                isoform_id,
                gene,
                barcode,
                molecule_id,
                source,
            )

            row = get_index(barcode, barcode_index, barcodes)
            col = get_index(isoform_id, isoform_index, isoforms)
            rows.append(row)
            cols.append(col)
            data.append(1)
            isoform_gene[isoform_id] = gene
            stats['kept_reads'] += 1

    bam.close()

    if not isoforms or not barcodes:
        print("No isoforms or barcodes were collected. Skipping h5ad output.")
        return gff_path, None, stats

    matrix = coo_matrix((data, (rows, cols)),
                        shape=(len(barcodes), len(isoforms)),
                        dtype=np.int32).tocsr()
    adata = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=isoforms),
    )
    adata.var['gene'] = [isoform_gene.get(isoform, '') for isoform in isoforms]
    adata.write_h5ad(h5ad_path, compression='gzip')

    return gff_path, h5ad_path, stats


def main():
    parser = argparse.ArgumentParser(
        description="Extract spliced isoform structures from a long-read single-cell BAM."
    )
    parser.add_argument('bam_file', help='Input BAM file path.')
    parser.add_argument('--gene-model', dest='gene_model', required=True,
                        help='Gene model reference file.')
    parser.add_argument('--output-prefix', '--output_prefix', dest='output_prefix',
                        required=True, help='Output prefix for .gff and .h5ad files.')
    parser.add_argument('--gene', dest='gene', default=None,
                        help='Optional Ensembl gene ID for QC mode.')
    parser.add_argument('--min-mapq', dest='min_mapq', type=int, default=1,
                        help='Minimum mapping quality to retain a read.')
    parser.add_argument('--barcode-tags', dest='barcode_tags',
                        default='CB',
                        help='Comma-separated list of barcode tags to scan.')
    parser.add_argument('--molecule-tag', dest='molecule_tag', default='zm',
                        help='Tag to use for molecule ID (fallback: read name).')
    parser.add_argument('--source', dest='source', default='bam',
                        help='Source field to write in the GFF file.')
    parser.add_argument('--require-known-splice', dest='require_known_splice',
                        action=argparse.BooleanOptionalAction, default=True,
                        help='Require at least one known Ensembl splice site per read.')
    args = parser.parse_args()

    barcode_tags = [tag.strip() for tag in args.barcode_tags.split(',') if tag.strip()]
    gff_path, h5ad_path, stats = extract_isoform_structures(
        args.bam_file,
        args.gene_model,
        args.output_prefix,
        target_gene=args.gene,
        min_mapq=args.min_mapq,
        barcode_tags=barcode_tags,
        molecule_tag=args.molecule_tag,
        source=args.source,
        require_known_splice=args.require_known_splice,
    )

    print(f"GFF written to: {gff_path}")
    if h5ad_path:
        print(f"h5ad written to: {h5ad_path}")
    for key in (
        'total_reads',
        'mapped_primary_reads',
        'spliced_reads',
        'gene_assigned_spliced_reads',
        'barcode_reads',
        'known_splice_reads',
        'kept_reads',
    ):
        print(f"{key}: {stats.get(key, 0)}")


if __name__ == '__main__':
    main()
