#!/usr/bin/env python3
import argparse
import os
import sys
import time
from collections import defaultdict
from datetime import datetime
import multiprocessing
from pathlib import Path
import shutil
import csv
import gzip
import sqlite3

import anndata as ad
import numpy as np
import pandas as pd
import pysam
from scipy.sparse import coo_matrix

sys.path.insert(1, os.path.join(os.path.dirname(__file__), '..'))
from long_read import gff_process

_EXON_MODEL_PATH = None


def _ensure_exon_model(exon_file):
    global _EXON_MODEL_PATH
    exon_path = os.path.abspath(exon_file)
    if _EXON_MODEL_PATH != exon_path:
        gff_process.importEnsemblGenes(exon_file)
        _EXON_MODEL_PATH = exon_path
    return gff_process.exonCoordinates


def _init_chunk_worker(exon_file):
    _ensure_exon_model(exon_file)


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
    if not read.cigartuples:
        return exons
    ref_pos = read.reference_start
    exon_start = None
    for op, length in read.cigartuples:
        if op in (0, 7, 8, 2):  # M/=/X/D consume reference
            if exon_start is None:
                exon_start = ref_pos
            ref_pos += length
        elif op == 3:  # N splits exons
            if exon_start is not None:
                exon_end = ref_pos
                if exon_end >= exon_start + 1:
                    exons.append((exon_start + 1, exon_end))
                exon_start = None
            ref_pos += length
        else:
            # I/S/H/P do not consume reference
            continue
    if exon_start is not None:
        exon_end = ref_pos
        if exon_end >= exon_start + 1:
            exons.append((exon_start + 1, exon_end))
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
    transcript_id = isoform_id
    if isoform_id and ':' in isoform_id:
        transcript_id = isoform_id.split(':', 1)[1]
    if transcript_id and transcript_id.startswith('molecule/'):
        transcript_id = transcript_id.split('molecule/', 1)[1]
    info_parts.append(f'transcript_id "{transcript_id}"')
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
            molecule_core = molecule_id
            if molecule_core.startswith('molecule/'):
                molecule_core = molecule_core.split('molecule/', 1)[1]
            isoform_id = sanitize_gff_value(f"{gene}:{molecule_core}")
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


def _chunk_label(chrom, strand):
    safe_chrom = chrom.replace('/', '_')
    if strand is None:
        return f"{safe_chrom}_all"
    safe_strand = 'plus' if strand == '+' else 'minus'
    return f"{safe_chrom}_{safe_strand}"


def _write_chunk_stats(stats, stats_path):
    with open(stats_path, 'w') as handle:
        for key in (
            'total_reads',
            'mapped_primary_reads',
            'spliced_reads',
            'gene_assigned_spliced_reads',
            'barcode_reads',
            'known_splice_reads',
            'kept_reads',
        ):
            handle.write(f"{key}\t{stats.get(key, 0)}\n")


def extract_isoform_structures_chunk(bam_path, exon_file, output_prefix, chrom, strand,
                                     min_mapq=1, barcode_tags=None, molecule_tag='zm',
                                     source='bam', require_known_splice=True, chunk_dir=None):
    if barcode_tags is None:
        barcode_tags = ['CB', 'CR', 'BC', 'BX']
    if chunk_dir is None:
        chunk_dir = Path(output_prefix).parent / 'chr-results'
    chunk_dir = Path(chunk_dir)
    chunk_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = Path(output_prefix)
    chunk_label = _chunk_label(chrom, strand)
    gff_path = chunk_dir / f"{output_prefix.name}.{chunk_label}.gff"
    counts_path = chunk_dir / f"{output_prefix.name}.{chunk_label}.counts.tsv"
    stats_path = chunk_dir / f"{output_prefix.name}.{chunk_label}.stats.tsv"

    exon_coordinates = _ensure_exon_model(exon_file)
    bam = pysam.AlignmentFile(bam_path, 'rb')
    stats = defaultdict(int)
    counts_db_path = chunk_dir / f"{output_prefix.name}.{chunk_label}.counts.sqlite"
    conn = sqlite3.connect(counts_db_path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA locking_mode=EXCLUSIVE")
    conn.execute(
        "CREATE TABLE IF NOT EXISTS counts ("
        "barcode TEXT NOT NULL, "
        "isoform_id TEXT NOT NULL, "
        "count INTEGER NOT NULL, "
        "gene TEXT, "
        "PRIMARY KEY (barcode, isoform_id)"
        ")"
    )
    insert_sql = (
        "INSERT INTO counts (barcode, isoform_id, count, gene) "
        "VALUES (?, ?, ?, ?) "
        "ON CONFLICT(barcode, isoform_id) DO UPDATE SET "
        "count = count + excluded.count, "
        "gene = COALESCE(counts.gene, excluded.gene)"
    )
    counts_buffer = []
    buffer_limit = 10000
    counts_seen = 0
    conn.execute("BEGIN")

    with open(gff_path, 'w') as gff_handle:
        for read in bam.fetch(chrom):
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
            chrom_norm = normalize_chrom(read.reference_name, exon_coordinates)
            strand_call = get_read_strand(read)
            if strand is not None and strand_call != strand:
                continue
            gene, _, _, genes = gff_process.exonAnnotate(
                chrom_norm, list(exons), strand_call, read.query_name
            )
            if not genes:
                continue
            stats['gene_assigned_spliced_reads'] += 1

            has_known = has_known_splice_site(chrom_norm, strand_call, exons, exon_coordinates)
            if not has_known and require_known_splice:
                continue
            if has_known:
                stats['known_splice_reads'] += 1

            molecule_id = sanitize_gff_value(resolve_molecule_id(read, molecule_tag))
            molecule_core = molecule_id
            if molecule_core.startswith('molecule/'):
                molecule_core = molecule_core.split('molecule/', 1)[1]
            isoform_id = sanitize_gff_value(f"{gene}:{molecule_core}")
            write_gff_isoform(
                gff_handle,
                chrom_norm,
                strand_call,
                exons,
                isoform_id,
                gene,
                barcode,
                molecule_id,
                source,
            )
            gene_value = gene if gene else None
            counts_buffer.append((barcode, isoform_id, 1, gene_value))
            counts_seen += 1
            if len(counts_buffer) >= buffer_limit:
                conn.executemany(insert_sql, counts_buffer)
                counts_buffer.clear()
            stats['kept_reads'] += 1

    bam.close()
    if counts_buffer:
        conn.executemany(insert_sql, counts_buffer)
        counts_buffer.clear()
    conn.commit()
    if counts_seen:
        with open(counts_path, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['barcode', 'isoform_id', 'count', 'gene'])
            for barcode, isoform_id, count, gene in conn.execute(
                    "SELECT barcode, isoform_id, count, gene FROM counts"):
                writer.writerow([barcode, isoform_id, count, gene or ''])
    conn.close()
    try:
        counts_db_path.unlink()
    except FileNotFoundError:
        pass
    _write_chunk_stats(stats, stats_path)
    return str(gff_path), str(counts_path), stats


def _combine_chunk_gffs(chunk_dir, output_prefix):
    gff_path = Path(output_prefix).with_suffix('.gff')
    chunk_dir = Path(chunk_dir)
    with open(gff_path, 'w') as out_handle:
        for chunk in sorted(chunk_dir.glob(f"{Path(output_prefix).name}.*.gff")):
            with open(chunk, 'r') as in_handle:
                shutil.copyfileobj(in_handle, out_handle)
    return gff_path


def _combine_chunk_counts(chunk_dir, output_prefix):
    chunk_dir = Path(chunk_dir)
    counts_files = sorted(chunk_dir.glob(f"{Path(output_prefix).name}.*.counts.tsv"))
    barcodes = []
    isoforms = []
    barcode_index = {}
    isoform_index = {}
    isoform_gene = {}
    nnz = 0
    for counts_path in counts_files:
        with open(counts_path, 'r', newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            for row in reader:
                barcode = row.get('barcode')
                isoform_id = row.get('isoform_id')
                if not barcode or not isoform_id:
                    continue
                gene = row.get('gene', '')
                get_index(barcode, barcode_index, barcodes)
                get_index(isoform_id, isoform_index, isoforms)
                nnz += 1
                if gene:
                    isoform_gene[isoform_id] = gene
    if not isoforms or not barcodes:
        return None, None, {}
    if nnz == 0:
        return None, None, {}
    memmap_dir = chunk_dir / f"{Path(output_prefix).name}.memmap"
    memmap_dir.mkdir(parents=True, exist_ok=True)
    rows_path = memmap_dir / "rows.dat"
    cols_path = memmap_dir / "cols.dat"
    data_path = memmap_dir / "data.dat"
    rows = np.memmap(rows_path, dtype=np.int32, mode='w+', shape=(nnz,))
    cols = np.memmap(cols_path, dtype=np.int32, mode='w+', shape=(nnz,))
    data = np.memmap(data_path, dtype=np.int32, mode='w+', shape=(nnz,))
    idx = 0
    for counts_path in counts_files:
        with open(counts_path, 'r', newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            for row in reader:
                barcode = row.get('barcode')
                isoform_id = row.get('isoform_id')
                if not barcode or not isoform_id:
                    continue
                try:
                    count = int(row.get('count', 1))
                except ValueError:
                    count = 1
                rows[idx] = barcode_index[barcode]
                cols[idx] = isoform_index[isoform_id]
                data[idx] = count
                idx += 1
    matrix = coo_matrix((data, (rows, cols)),
                        shape=(len(barcodes), len(isoforms)),
                        dtype=np.int32).tocsr()
    h5ad_path = Path(output_prefix).with_suffix('.h5ad')
    adata = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=isoforms),
    )
    adata.var['gene'] = [isoform_gene.get(isoform, '') for isoform in isoforms]
    adata.write_h5ad(h5ad_path, compression='gzip')
    try:
        rows.flush()
        cols.flush()
        data.flush()
    except AttributeError:
        pass
    for path in (rows_path, cols_path, data_path):
        try:
            path.unlink()
        except FileNotFoundError:
            continue
    try:
        memmap_dir.rmdir()
    except OSError:
        pass
    return h5ad_path, isoform_gene, {}


def _cleanup_chunk_outputs(chunk_dir, output_prefix):
    chunk_dir = Path(chunk_dir)
    prefix = Path(output_prefix).name
    patterns = [
        f"{prefix}.*.gff",
        f"{prefix}.*.counts.tsv",
        f"{prefix}.*.counts.sqlite",
        f"{prefix}.*.stats.tsv",
    ]
    for pattern in patterns:
        for path in chunk_dir.glob(pattern):
            try:
                path.unlink()
            except FileNotFoundError:
                continue


def _gzip_file(path):
    path = Path(path)
    gz_path = path.with_suffix(path.suffix + '.gz')
    with open(path, 'rb') as source, gzip.open(gz_path, 'wb') as target:
        shutil.copyfileobj(source, target)
    path.unlink()
    return gz_path


def parallel_extract_isoform_structures(bam_path, exon_file, output_prefix, min_mapq=1,
                                        barcode_tags=None, molecule_tag='zm',
                                        source='bam', require_known_splice=True,
                                        max_processes=None):
    if barcode_tags is None:
        barcode_tags = ['CB', 'CR', 'BC', 'BX']
    bam = pysam.AlignmentFile(bam_path, 'rb')
    chromosomes = [ref for ref in bam.references if len(ref) < 6]
    bam.close()
    if not chromosomes:
        raise ValueError("No chromosomes found in BAM header.")

    cpu_total = multiprocessing.cpu_count()
    if max_processes is None:
        max_processes = max(1, cpu_total - 1)
    max_processes = max(1, min(max_processes, cpu_total))

    output_prefix = Path(output_prefix)
    chunk_dir = output_prefix.parent / 'chr-results'
    chunk_dir.mkdir(parents=True, exist_ok=True)
    _cleanup_chunk_outputs(chunk_dir, output_prefix)
    tasks = []
    for chrom in chromosomes:
        tasks.append((
            bam_path,
            exon_file,
            str(output_prefix),
            chrom,
            None,
            min_mapq,
            barcode_tags,
            molecule_tag,
            source,
            require_known_splice,
            str(chunk_dir),
        ))

    stats = defaultdict(int)
    if max_processes > 1:
        with multiprocessing.Pool(
                processes=max_processes,
                initializer=_init_chunk_worker,
                initargs=(exon_file,)) as pool:
            for _gff_path, _counts_path, chunk_stats in pool.starmap(extract_isoform_structures_chunk, tasks):
                for key, value in chunk_stats.items():
                    stats[key] += value
    else:
        _ensure_exon_model(exon_file)
        for task in tasks:
            _gff_path, _counts_path, chunk_stats = extract_isoform_structures_chunk(*task)
            for key, value in chunk_stats.items():
                stats[key] += value

    gff_path = _combine_chunk_gffs(chunk_dir, output_prefix)
    h5ad_path, _isoform_gene, _ = _combine_chunk_counts(chunk_dir, output_prefix)
    _cleanup_chunk_outputs(chunk_dir, output_prefix)
    gff_path = _gzip_file(gff_path)
    return gff_path, h5ad_path, stats, max_processes, len(tasks)


def chunked_extract_isoform_structures(bam_path, exon_file, output_prefix, min_mapq=1,
                                       barcode_tags=None, molecule_tag='zm',
                                       source='bam', require_known_splice=True):
    gff_path, h5ad_path, stats, _used_processes, _chunk_count = parallel_extract_isoform_structures(
        bam_path,
        exon_file,
        output_prefix,
        min_mapq=min_mapq,
        barcode_tags=barcode_tags,
        molecule_tag=molecule_tag,
        source=source,
        require_known_splice=require_known_splice,
        max_processes=1,
    )
    return gff_path, h5ad_path, stats


class _Tee:
    def __init__(self, *streams):
        self._streams = streams

    def write(self, message):
        for stream in self._streams:
            stream.write(message)
            stream.flush()

    def flush(self):
        for stream in self._streams:
            stream.flush()


def _init_logging(output_prefix):
    output_prefix = Path(output_prefix).expanduser()
    output_dir = output_prefix.parent if output_prefix.parent.as_posix() else Path('.')
    log_dir = output_dir / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_name = f"{output_prefix.name}_isoform_structure_extract_{timestamp}.log"
    log_path = log_dir / log_name
    log_handle = open(log_path, 'w')
    stdout_tee = _Tee(sys.stdout, log_handle)
    stderr_tee = _Tee(sys.stderr, log_handle)
    return log_path, log_handle, stdout_tee, stderr_tee


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
    parser.add_argument('--no-parallel', dest='no_parallel', action='store_true',
                        help='Disable parallel chromosome/strand processing when no gene is provided.')
    parser.add_argument('--max-processes', dest='max_processes', type=int, default=None,
                        help='Maximum parallel workers (default: cpu_count-1).')
    args = parser.parse_args()

    start_time = time.time()
    log_path, log_handle, stdout_tee, stderr_tee = _init_logging(args.output_prefix)
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    sys.stdout = stdout_tee
    sys.stderr = stderr_tee
    try:
        print(f"Run started: {datetime.now().isoformat(timespec='seconds')}")
        print("Command:", " ".join(sys.argv))
        print(f"Log file: {log_path}")
        barcode_tags = [tag.strip() for tag in args.barcode_tags.split(',') if tag.strip()]
        if args.gene is None:
            if args.no_parallel:
                gff_path, h5ad_path, stats, used_processes, chunk_count = parallel_extract_isoform_structures(
                    args.bam_file,
                    args.gene_model,
                    args.output_prefix,
                    min_mapq=args.min_mapq,
                    barcode_tags=barcode_tags,
                    molecule_tag=args.molecule_tag,
                    source=args.source,
                    require_known_splice=args.require_known_splice,
                    max_processes=1,
                )
                print(f"Processing mode: chunked sequential (chunks={chunk_count}, workers={used_processes})")
            else:
                gff_path, h5ad_path, stats, used_processes, chunk_count = parallel_extract_isoform_structures(
                    args.bam_file,
                    args.gene_model,
                    args.output_prefix,
                    min_mapq=args.min_mapq,
                    barcode_tags=barcode_tags,
                    molecule_tag=args.molecule_tag,
                    source=args.source,
                    require_known_splice=args.require_known_splice,
                    max_processes=args.max_processes,
                )
                print(f"Processing mode: parallel chrom (chunks={chunk_count}, workers={used_processes})")
        else:
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
            print("Processing mode: single-process (gene)")

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
        elapsed = time.time() - start_time
        print(f"Elapsed time (s): {elapsed:.2f}")
    finally:
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        log_handle.close()


if __name__ == '__main__':
    main()
