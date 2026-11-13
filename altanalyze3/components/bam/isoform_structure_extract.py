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
import heapq
import itertools
import json
import re

import anndata as ad
import numpy as np
import pandas as pd
import pysam
from scipy.sparse import coo_matrix

sys.path.insert(1, os.path.join(os.path.dirname(__file__), '..'))
from altanalyze3.components.long_read import gff_process

_EXON_MODEL_PATH = None
_COUNTS_AGGREGATOR = os.environ.get('ISOFORM_COUNTS_AGG', 'sort-merge')
_COUNTS_RUN_SIZE = int(os.environ.get('ISOFORM_COUNTS_RUN_SIZE', '200000'))
_MP_START_METHOD = os.environ.get('ISOFORM_MP_START', 'fork')
_MP_TASKS_PER_CHILD = int(os.environ.get('ISOFORM_MP_TASKS_PER_CHILD', '1'))
_PYSAM_VERBOSITY = int(os.environ.get('ISOFORM_PYSAM_VERBOSITY', '0'))


def _read_int_env(name):
    value = os.environ.get(name)
    if not value:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _resolve_max_processes(user_max):
    cpu_total = multiprocessing.cpu_count()
    env_candidates = (
        _read_int_env('ISOFORM_MAX_PROCESSES'),
        _read_int_env('LSB_DJOB_NUMPROC'),
        _read_int_env('SLURM_CPUS_PER_TASK'),
        _read_int_env('PBS_NP'),
        _read_int_env('NSLOTS'),
    )
    env_max = next((value for value in env_candidates if value and value > 0), None)
    if user_max is None:
        max_processes = max(1, min(cpu_total - 1, 8))
    else:
        max_processes = user_max
    if env_max is not None:
        max_processes = min(max_processes, env_max)
    return max(1, min(max_processes, cpu_total))


def _set_pysam_verbosity():
    try:
        pysam.set_verbosity(_PYSAM_VERBOSITY)
    except (AttributeError, TypeError, ValueError):
        pass


def _ensure_bam_index(bam_path):
    """Make sure a coordinate index (.bai/.csi) exists for ``bam_path``; build it if missing.

    Per-chromosome ``fetch`` requires an index, so this is called once before any worker reads the
    BAM. Indexing requires a coordinate-sorted BAM; if the BAM is not coordinate-sorted we sort it
    to a sibling ``*.sorted.bam`` and index that, returning the path actually usable for fetch. The
    common case (already indexed) returns immediately.
    """
    bam_path = str(bam_path)
    for idx in (bam_path + '.bai', bam_path + '.csi',
                bam_path[:-4] + '.bai' if bam_path.endswith('.bam') else bam_path + '.bai'):
        if os.path.exists(idx):
            return bam_path

    # Determine sort order from the header without reading records.
    try:
        with pysam.AlignmentFile(bam_path, 'rb', check_sq=False) as bam:
            sort_order = bam.header.to_dict().get('HD', {}).get('SO', '')
    except Exception as exc:
        raise ValueError(f"Unable to read BAM header for {bam_path}: {exc}") from exc

    if sort_order == 'coordinate':
        print(f"[bam] no index for {os.path.basename(bam_path)}; building .bai ...")
        try:
            pysam.index(bam_path)
            return bam_path
        except Exception as exc:
            # Index build failed despite a coordinate header (e.g. truncated/stale .bai, read-only
            # dir). Fall through to the sort+index path, which writes a fresh sibling copy.
            print(f"[bam] indexing failed ({type(exc).__name__}: {exc}); "
                  f"falling back to sort+index")

    # Not coordinate-sorted (or index build failed): sort then index. Idempotent -- reuse a prior
    # sorted+indexed copy when present.
    sorted_path = bam_path[:-4] + '.sorted.bam' if bam_path.endswith('.bam') else bam_path + '.sorted.bam'
    if os.path.exists(sorted_path) and os.path.exists(sorted_path + '.bai'):
        return sorted_path
    try:
        print(f"[bam] {os.path.basename(bam_path)} requires a sorted+indexed copy "
              f"(SO='{sort_order or 'unknown'}'); sorting -> {os.path.basename(sorted_path)} ...")
        pysam.sort('-o', sorted_path, bam_path)
        pysam.index(sorted_path)
        return sorted_path
    except Exception as exc:
        raise ValueError(
            f"Could not produce a usable coordinate index for {bam_path}. The BAM may be corrupt or "
            f"the output directory may not be writable. Underlying error: "
            f"{type(exc).__name__}: {exc}"
        ) from exc


def _ensure_exon_model(exon_file):
    global _EXON_MODEL_PATH
    exon_path = os.path.abspath(exon_file)
    if _EXON_MODEL_PATH != exon_path:
        gff_process.importEnsemblGenes(exon_file)
        _EXON_MODEL_PATH = exon_path
    return gff_process.exonCoordinates


def _init_chunk_worker(exon_file):
    _ensure_exon_model(exon_file)


def _parse_count_row(line):
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 2:
        return None
    barcode = parts[0]
    isoform_id = parts[1]
    gene = parts[2] if len(parts) > 2 else ''
    return barcode, isoform_id, gene


def _flush_sorted_run(records, run_path):
    records.sort(key=lambda item: (item[0], item[1]))
    with open(run_path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerows(records)


def _write_sorted_runs(raw_path, runs_dir, max_records):
    runs_dir.mkdir(parents=True, exist_ok=True)
    runs = []
    buffer_records = []
    with open(raw_path, 'r') as handle:
        for line in handle:
            parsed = _parse_count_row(line)
            if parsed is None:
                continue
            buffer_records.append(parsed)
            if len(buffer_records) >= max_records:
                run_path = runs_dir / f"run_{len(runs):05d}.tsv"
                _flush_sorted_run(buffer_records, run_path)
                runs.append(run_path)
                buffer_records.clear()
    if buffer_records:
        run_path = runs_dir / f"run_{len(runs):05d}.tsv"
        _flush_sorted_run(buffer_records, run_path)
        runs.append(run_path)
    return runs


def _merge_sorted_runs(run_paths, counts_path):
    handles = [open(path, 'r') for path in run_paths]
    heap = []
    for idx, handle in enumerate(handles):
        line = handle.readline()
        if line:
            parsed = _parse_count_row(line)
            if parsed is not None:
                barcode, isoform_id, gene = parsed
                heapq.heappush(heap, (barcode, isoform_id, gene, idx))
    with open(counts_path, 'w', newline='') as out_handle:
        writer = csv.writer(out_handle, delimiter='\t')
        writer.writerow(['barcode', 'isoform_id', 'count', 'gene'])
        current_key = None
        current_count = 0
        current_gene = ''
        while heap:
            barcode, isoform_id, gene, idx = heapq.heappop(heap)
            key = (barcode, isoform_id)
            if current_key is None:
                current_key = key
                current_count = 1
                current_gene = gene
            elif key == current_key:
                current_count += 1
                if not current_gene and gene:
                    current_gene = gene
            else:
                writer.writerow([current_key[0], current_key[1], current_count, current_gene])
                current_key = key
                current_count = 1
                current_gene = gene
            line = handles[idx].readline()
            if line:
                parsed = _parse_count_row(line)
                if parsed is not None:
                    b_next, i_next, g_next = parsed
                    heapq.heappush(heap, (b_next, i_next, g_next, idx))
        if current_key is not None:
            writer.writerow([current_key[0], current_key[1], current_count, current_gene])
    for handle in handles:
        handle.close()


def _aggregate_counts_sort_merge(raw_path, counts_path, runs_dir, max_records):
    run_paths = _write_sorted_runs(raw_path, runs_dir, max_records)
    if not run_paths:
        return False
    _merge_sorted_runs(run_paths, counts_path)
    for path in run_paths:
        try:
            path.unlink()
        except FileNotFoundError:
            continue
    try:
        runs_dir.rmdir()
    except OSError:
        pass
    return True


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


def get_read_strand(read, ts_read_relative=False):
    """Genomic strand of the transcript this read came from.

    ``exonAnnotate`` looks up ``exonCoordinates[(chr, pos, strand, 1|2)]``, and
    ``importEnsemblGenes`` registers every boundary under the GENE's own strand only. The strand
    returned here therefore gates gene assignment: get it wrong and the read finds no gene.

    ts_read_relative=False (default, unchanged): the ``ts``/``XS``/``TS`` tag is taken verbatim.

    ts_read_relative=True: ``ts`` is interpreted as minimap2 documents it, the transcript strand
    RELATIVE TO THE READ, so the genomic strand is the tag XOR the alignment orientation. Needed for
    BAMs where minimap2/pbmm2 emits ``ts`` this way. Measured on the ENCODE PC-3 long-read BAM
    (ENCFF044LIA), where every read carries ``ts:A:+``: of 20,000 chr7 spliced reads, the 9,352
    forward reads key into exonCoordinates under '+' (98.8%) and the 10,648 reverse reads key in
    under '-' (99.2%) and not under '+' (0.1%). Verified independently against the genome sequence:
    forward reads show GT..AG introns (24,443/25,008), reverse reads show CT..AC (10,860/11,064),
    which is the reverse complement, so their transcripts are on the minus strand.

    ``XS``/``TS`` are written by HISAT2/STAR/TopHat as GENOMIC strand, so they are always taken
    verbatim; only ``ts`` is subject to the read-relative interpretation.
    """
    for tag in ('ts', 'XS', 'TS'):
        try:
            strand = read.get_tag(tag)
        except KeyError:
            continue
        if strand in ('+', '-'):
            if tag == 'ts' and ts_read_relative and read.is_reverse:
                return '-' if strand == '+' else '+'
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


def resolve_barcode(read, barcode_tags, default_barcode=None):
    """Return the cell barcode from the first present tag, normalized to the standard single-cell
    convention with a ``-1`` lane suffix (e.g. ``AAACACCTGTTAGTGC-1``).

    Appending ``-1`` at extraction makes the suffix canonical at every downstream step (junction /
    isoform matrices, pseudobulk, comparisons) and keeps the barcodes compatible with cell
    annotations produced outside this workflow and with the reverse-complement code path. Barcodes
    that already carry a ``-<lane>`` suffix are left as-is.

    ``default_barcode`` is the BULK long-read path. Bulk libraries (PacBio CCS / ONT, e.g. the ENCODE
    long-read RNA-seq BAMs) carry no ``CB``/``CR``/``BC``/``BX`` tag, so every read would otherwise be
    discarded by the caller's ``if not barcode: continue``. Passing the library name here makes the
    whole BAM one pseudo-cell, which lets every downstream single-cell step (molecule h5ad, collapse,
    re-key, pseudobulk) run UNCHANGED and sum to a sample-level count. It is returned VERBATIM: no
    ``-1`` suffix is appended, because a library name is not a 10x barcode and the suffix would only
    break the barcode->annotation join. ``None`` (the default) preserves single-cell behaviour
    exactly: a read with no barcode tag is still dropped.
    """
    for tag in barcode_tags:
        try:
            value = read.get_tag(tag)
        except KeyError:
            continue
        if value:
            value = str(value)
            return value if '-' in value else f"{value}-1"
    return default_barcode


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
                               source='bam', require_known_splice=True, default_barcode=None,
                               ts_read_relative=False):
    if barcode_tags is None:
        barcode_tags = ['CB', 'CR', 'BC', 'BX']

    _set_pysam_verbosity()
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
            barcode = resolve_barcode(read, barcode_tags, default_barcode=default_barcode)
            if not barcode:
                continue
            stats['barcode_reads'] += 1

            exons = extract_exons(read)
            if len(exons) < 2:
                continue
            chrom = normalize_chrom(read.reference_name, exon_coordinates)
            strand = get_read_strand(read, ts_read_relative=ts_read_relative)
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
                                     source='bam', require_known_splice=True, chunk_dir=None,
                                     default_barcode=None, ts_read_relative=False):
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

    _set_pysam_verbosity()
    exon_coordinates = _ensure_exon_model(exon_file)
    bam = pysam.AlignmentFile(bam_path, 'rb')
    stats = defaultdict(int)
    counts_seen = 0
    use_sqlite = _COUNTS_AGGREGATOR == 'sqlite'
    counts_db_path = None
    conn = None
    insert_sql = None
    counts_buffer = []
    buffer_limit = 10000
    raw_counts_path = None
    raw_handle = None
    if use_sqlite:
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
        conn.execute("BEGIN")
    else:
        raw_counts_path = chunk_dir / f"{output_prefix.name}.{chunk_label}.counts.raw.tsv"
        raw_handle = open(raw_counts_path, 'w', newline='')

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
            barcode = resolve_barcode(read, barcode_tags, default_barcode=default_barcode)
            if not barcode:
                continue
            stats['barcode_reads'] += 1

            exons = extract_exons(read)
            if len(exons) < 2:
                continue
            chrom_norm = normalize_chrom(read.reference_name, exon_coordinates)
            strand_call = get_read_strand(read, ts_read_relative=ts_read_relative)
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
            gene_value = gene if gene else ''
            if use_sqlite:
                counts_buffer.append((barcode, isoform_id, 1, gene_value or None))
                if len(counts_buffer) >= buffer_limit:
                    conn.executemany(insert_sql, counts_buffer)
                    counts_buffer.clear()
            else:
                raw_handle.write(f"{barcode}\t{isoform_id}\t{gene_value}\n")
            counts_seen += 1
            stats['kept_reads'] += 1

    bam.close()
    if use_sqlite:
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
        if counts_db_path is not None:
            try:
                counts_db_path.unlink()
            except FileNotFoundError:
                pass
    else:
        raw_handle.close()
        if counts_seen:
            runs_dir = chunk_dir / f"{output_prefix.name}.{chunk_label}.counts.runs"
            _aggregate_counts_sort_merge(
                raw_counts_path,
                counts_path,
                runs_dir,
                _COUNTS_RUN_SIZE,
            )
        if raw_counts_path is not None:
            try:
                raw_counts_path.unlink()
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
        f"{prefix}.*.counts.raw.tsv",
        f"{prefix}.*.counts.sqlite",
        f"{prefix}.*.stats.tsv",
    ]
    for pattern in patterns:
        for path in chunk_dir.glob(pattern):
            try:
                path.unlink()
            except FileNotFoundError:
                continue
    for path in chunk_dir.glob(f"{prefix}.*.counts.runs"):
        if path.is_dir():
            shutil.rmtree(path, ignore_errors=True)


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
                                        max_processes=None, default_barcode=None,
                                        ts_read_relative=False):
    if barcode_tags is None:
        barcode_tags = ['CB', 'CR', 'BC', 'BX']
    _set_pysam_verbosity()
    # Per-chromosome fetch needs a coordinate index; build one (or a sorted+indexed copy) if missing.
    bam_path = _ensure_bam_index(bam_path)
    bam = pysam.AlignmentFile(bam_path, 'rb')
    chromosomes = [ref for ref in bam.references if len(ref) < 6]
    bam.close()
    if not chromosomes:
        raise ValueError("No chromosomes found in BAM header.")

    max_processes = _resolve_max_processes(max_processes)

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
            default_barcode,
            ts_read_relative,
        ))

    stats = defaultdict(int)
    if max_processes > 1:
        def _run_pool(ctx):
            with ctx.Pool(
                    processes=max_processes,
                    initializer=_init_chunk_worker,
                    initargs=(exon_file,),
                    maxtasksperchild=_MP_TASKS_PER_CHILD or None) as pool:
                for _gff_path, _counts_path, chunk_stats in pool.starmap(extract_isoform_structures_chunk, tasks):
                    for key, value in chunk_stats.items():
                        stats[key] += value

        ctx = multiprocessing.get_context(_MP_START_METHOD)
        try:
            _run_pool(ctx)
        except RuntimeError as exc:
            if _MP_START_METHOD != 'fork' and 'bootstrapping phase' in str(exc):
                ctx = multiprocessing.get_context('fork')
                _run_pool(ctx)
            else:
                raise
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
                                       source='bam', require_known_splice=True,
                                       default_barcode=None, ts_read_relative=False):
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
        default_barcode=default_barcode,
        ts_read_relative=ts_read_relative,
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


# ------------------------------------------------------------------------------------------------
# EXPERIMENTAL FORK: extra-chr-trans-splicing (2026-10-29). NOT validated altanalyze3 code.
# The validated extractors above skip supplementary alignments (is_supplementary -> continue), so a
# molecule whose alignment is split across chromosomes keeps only its primary piece. This mode
# rebuilds every aligned piece of such a molecule from the primary record and its SA tag, and writes
# ONLY cross-chromosome molecules to the read GFF and the molecule h5ad. Exon annotation, barcodes,
# molecule ids and the GFF/h5ad formats are the validated functions above.
# ------------------------------------------------------------------------------------------------
_CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')
_CIGAR_CODES = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
_COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def _cigar_from_string(cigar):
    parts = _CIGAR_RE.findall(cigar)
    if ''.join(n + op for n, op in parts) != cigar:
        raise ValueError(f'Unparsed SA CIGAR: {cigar}')
    return [(_CIGAR_CODES[op], int(n)) for n, op in parts]


def _exons_from_cigar(reference_start, cigartuples):
    """extract_exons(read) arithmetic for an alignment given as (0-based start, CIGAR tuples)."""
    exons = []
    ref_pos = reference_start
    exon_start = None
    for op, length in cigartuples:
        if op in (0, 7, 8, 2):
            if exon_start is None:
                exon_start = ref_pos
            ref_pos += length
        elif op == 3:
            if exon_start is not None:
                if ref_pos >= exon_start + 1:
                    exons.append((exon_start + 1, ref_pos))
                exon_start = None
            ref_pos += length
    if exon_start is not None and ref_pos >= exon_start + 1:
        exons.append((exon_start + 1, ref_pos))
    return exons


def _query_interval(cigartuples, is_reverse):
    """0-based half-open interval of the aligned query bases in the ORIGINAL read orientation."""
    lead = 0
    for op, length in cigartuples:
        if op in (4, 5):
            lead += length
        else:
            break
    aligned = sum(length for op, length in cigartuples if op in (0, 1, 7, 8))
    total = sum(length for op, length in cigartuples if op in (0, 1, 4, 5, 7, 8))
    start, end = lead, lead + aligned
    if is_reverse:
        start, end = total - end, total - start
    return start, end, total


def _alignment_segments(read, sa_entries, ts_read_relative, exon_coordinates):
    primary_exons = extract_exons(read)
    if _exons_from_cigar(read.reference_start, read.cigartuples) != primary_exons:
        raise ValueError(f'SA exon arithmetic differs from extract_exons for {read.query_name}')
    qs, qe, total = _query_interval(read.cigartuples, read.is_reverse)
    segments = [dict(chrom=normalize_chrom(read.reference_name, exon_coordinates),
                     strand=get_read_strand(read, ts_read_relative=ts_read_relative),
                     orientation='-' if read.is_reverse else '+', mapq=read.mapping_quality,
                     nm=read.get_tag('NM') if read.has_tag('NM') else '', cigar=read.cigarstring,
                     exons=primary_exons, primary=True, qstart=qs, qend=qe, read_length=total)]
    for entry in sa_entries:
        chrom, pos, strand, cigar, mapq, nm = entry.split(',')
        tuples = _cigar_from_string(cigar)
        qs, qe, total = _query_interval(tuples, strand == '-')
        segments.append(dict(chrom=normalize_chrom(chrom, exon_coordinates), strand=strand, orientation=strand,
                             mapq=int(mapq), nm=int(nm), cigar=cigar,
                             exons=_exons_from_cigar(int(pos) - 1, tuples), primary=False,
                             qstart=qs, qend=qe, read_length=total))
    segments.sort(key=lambda s: (s['qstart'], s['qend']))
    return segments


def _boundary_annotation(hit, gene_region_counts):
    if hit is None:
        return '', '', ''
    gene, exon, index = hit
    return f'{gene}:{exon}', index == 0, index == gene_region_counts[gene] - 1


def _segment_junction(a, b, exon_coordinates):
    """Breakpoint between read-adjacent segments a (5') and b (3'). When the two alignments overlap
    on the read by k bases, the k+1 possible split points are tested and the one with the most
    registered exon boundaries (donor = region end key 2, acceptor = region start key 1, the keys
    exonAnnotate uses) is kept; ties keep the smallest shift. Each shift assumes no indel inside the
    overlap."""
    gap = b['qstart'] - a['qend']
    overlap = max(0, -gap)
    a_last = a['exons'][-1] if a['strand'] == '+' else a['exons'][0]
    b_first = b['exons'][0] if b['strand'] == '+' else b['exons'][-1]
    donor0 = a_last[1] if a['strand'] == '+' else a_last[0]
    acceptor0 = b_first[0] if b['strand'] == '+' else b_first[1]
    step_a = 1 if a['strand'] == '+' else -1
    step_b = 1 if b['strand'] == '+' else -1
    best = None
    for shift in range(overlap + 1):
        if shift >= a_last[1] - a_last[0] + 1 or overlap - shift >= b_first[1] - b_first[0] + 1:
            continue
        donor = donor0 - shift * step_a
        acceptor = acceptor0 + (overlap - shift) * step_b
        donor_hit = exon_coordinates.get((a['chrom'], donor, a['strand'], 2))
        acceptor_hit = exon_coordinates.get((b['chrom'], acceptor, b['strand'], 1))
        score = (donor_hit is not None) + (acceptor_hit is not None)
        if best is None or score > best[0]:
            best = (score, shift, donor, acceptor, donor_hit, acceptor_hit)
    if best is None:
        return dict(gap=gap, overlap=overlap, shift='', donor=donor0, acceptor=acceptor0,
                    donor_hit=None, acceptor_hit=None, resolved=False)
    _, shift, donor, acceptor, donor_hit, acceptor_hit = best
    return dict(gap=gap, overlap=overlap, shift=shift, donor=donor, acceptor=acceptor,
                donor_hit=donor_hit, acceptor_hit=acceptor_hit, resolved=True)


def _trim_to_junction(a, b, junction):
    """Move the a 3' end to the donor and the b 5' end to the acceptor (removes the read overlap)."""
    if not junction['resolved']:
        return
    a_exons = list(a['exons']); b_exons = list(b['exons'])
    if a['strand'] == '+':
        a_exons[-1] = (a_exons[-1][0], junction['donor'])
    else:
        a_exons[0] = (junction['donor'], a_exons[0][1])
    if b['strand'] == '+':
        b_exons[0] = (junction['acceptor'], b_exons[0][1])
    else:
        b_exons[-1] = (b_exons[-1][0], junction['acceptor'])
    a['exons'] = a_exons; b['exons'] = b_exons


def _oriented_sequence(read, read_length):
    seq = read.query_sequence
    if not seq or len(seq) != read_length:
        return ''
    return seq.translate(_COMPLEMENT)[::-1] if read.is_reverse else seq


def _primary_for_supplementary(lookup, read):
    for entry in read.get_tag('SA').rstrip(';').split(';'):
        chrom, pos = entry.split(',')[:2]
        for candidate in lookup.fetch(chrom, int(pos) - 1, int(pos)):
            if (candidate.query_name == read.query_name and not candidate.is_supplementary
                    and not candidate.is_secondary):
                return candidate
    return None


def _sa_records_match(lookup, read, sa_entries):
    """Region-mode check: every SA entry must equal a real supplementary record of the same read."""
    matched = 0
    for entry in sa_entries:
        chrom, pos, strand, cigar = entry.split(',')[:4]
        for candidate in lookup.fetch(chrom, int(pos) - 1, int(pos)):
            if (candidate.query_name == read.query_name and candidate.is_supplementary
                    and candidate.reference_start == int(pos) - 1 and candidate.cigarstring == cigar
                    and ('-' if candidate.is_reverse else '+') == strand):
                matched += 1
                break
    return matched


def _write_cross_chromosome_gff(handle, segments, isoform_id, gene, source, first_in_file=False):
    """first_in_file: gff_process.consolidateLongReadGFFs uses the first non-comment line of a file
    only to locate transcript_id and never appends it as an exon (gff_process.py:878-888), so the
    first molecule written natively loses its first exon. Writing that molecule's transcript line
    once more at the top of the file keeps all of its exons; every reader ignores transcript lines
    except as the end-of-record marker."""
    transcript_id = isoform_id.split(':', 1)[1] if ':' in isoform_id else isoform_id
    if transcript_id.startswith('molecule/'):
        transcript_id = transcript_id.split('molecule/', 1)[1]
    lines = []
    for index, segment in enumerate(segments, 1):
        info = f'gene_id "{gene}";transcript_id "{transcript_id}";segment "{index}";'
        for exon_start, exon_end in sorted(segment['exons']):
            lines.append(f"{segment['chrom']}\t{source}\texon\t{exon_start}\t{exon_end}\t.\t{segment['strand']}\t.\t{info}\n")
    last = segments[-1]
    lines.append(f"{last['chrom']}\t{source}\ttranscript\t{min(s for s, _ in last['exons'])}\t"
                 f"{max(e for _, e in last['exons'])}\t.\t{last['strand']}\t.\t{info}\n")
    if first_in_file:
        lines.insert(0, lines[-1])
    handle.write(''.join(lines))
    return lines[-1]


_JUNCTION_COLUMNS = [
    'library', 'molecule_id', 'query_name', 'barcode', 'umi', 'status', 'isoform_id', 'read_length',
    'n_segments', 'pair_index', 'cross_chromosome', 'gene_a', 'gene_b',
    'a_chrom', 'a_strand', 'a_qstart', 'a_qend', 'a_mapq', 'a_nm', 'a_primary', 'a_cigar', 'a_exons',
    'b_chrom', 'b_strand', 'b_qstart', 'b_qend', 'b_mapq', 'b_nm', 'b_primary', 'b_cigar', 'b_exons',
    'query_gap', 'overlap_nt', 'shift_nt', 'donor', 'acceptor', 'donor_annotation',
    'donor_is_gene_first_region', 'donor_is_gene_last_region', 'acceptor_annotation',
    'acceptor_is_gene_first_region', 'acceptor_is_gene_last_region', 'junction_read_context',
    'structure_untrimmed', 'blacklist_membership', 'duplicated_sequence_near_partner']


def _fasta_name(chrom, names):
    bare = chrom[3:] if chrom.startswith('chr') else chrom
    if bare in names:
        return bare
    if chrom in names:
        return chrom
    return 'MT' if chrom == 'chrM' and 'MT' in names else chrom


def _hit_chrom(ctg):
    if ctg.startswith('chr'):
        return ctg
    return 'chrM' if ctg == 'MT' else 'chr' + ctg


def _duplicated_near_partner(segments, junctions, aligner, genome, window, min_identity, min_cover, cache):
    """Blacklist rule: a segment whose reference sequence also aligns, at >= min_identity over
    >= min_cover of its length, on the partner segment's chromosome within window bases of the
    partner. Such a joint can be one colinear transcript from duplicated sequence (paralog, orphon,
    segmental duplication) placed on two chromosomes. Returns a detail string or ''."""
    names = set(genome.references)
    for (a, b), junction in zip(zip(segments[:-1], segments[1:]), junctions):
        if not junction['cross']:
            continue
        for query, partner, label in ((b, a, 'b_near_a'), (a, b, 'a_near_b')):
            key = (query['chrom'], tuple(query['exons']))
            if key not in cache:
                c = _fasta_name(query['chrom'], names)
                seq = ''.join(genome.fetch(c, x - 1, y) for x, y in query['exons']) if c in names else ''
                cache[key] = (len(seq), [(_hit_chrom(h.ctg), h.r_st, h.r_en, h.mlen, h.blen, h.q_st, h.q_en)
                                         for h in aligner.map(seq)] if len(seq) >= 30 else [])
            qlen, hits = cache[key]
            p_lo = min(x for x, _ in partner['exons']) - window
            p_hi = max(y for _, y in partner['exons']) + window
            for ctg, r_st, r_en, mlen, blen, q_st, q_en in hits:
                identity, cover = mlen / blen if blen else 0.0, (q_en - q_st) / qlen if qlen else 0.0
                if ctg == partner['chrom'] and r_en >= p_lo and r_st <= p_hi and identity >= min_identity and cover >= min_cover:
                    return f'{label}:{ctg}:{r_st + 1}-{r_en}:identity={identity:.3f}:cover={cover:.3f}'
    return ''


def extract_cross_chromosome_isoforms(bam_path, exon_file, output_prefix, library, min_mapq=1,
                                      barcode_tags=None, molecule_tag='zm', source='bam',
                                      require_known_splice=True, default_barcode=None,
                                      ts_read_relative=False, regions=None, molecule_prefix='',
                                      blacklist_members=None, genome_index=None, genome_fasta=None,
                                      dup_window=2000000, dup_min_identity=0.9, dup_min_cover=0.5):
    """EXPERIMENTAL FORK. Writes <prefix>.gff.gz and <prefix>.h5ad holding ONLY molecules aligned to
    two or more chromosomes, plus <prefix>.cross_chromosome_junctions.tsv.gz (every cross-chromosome
    primary read, kept or not, one row per read-adjacent segment pair) and
    <prefix>.cross_chromosome_stats.json.

    require_known_splice mirrors the validated rule (at least one registered splice site) at the
    cross-chromosome joint: at least one side of one cross-chromosome joint must be a registered
    exon boundary. regions (list of (chrom, start0, end)) restricts the scan for demonstrations.
    molecule_prefix is prepended to every molecule id so ids stay unique across libraries: the reused
    collapse copies exemplar GFF lines with their bare molecule id into combined.gff, and
    gff_translate keys transcripts by that id, so a repeated id would merge two molecules.
    blacklist_members (reference/blacklist_members.tsv) and genome_index (minimap2 splice index of
    genome_fasta) enable the blacklist: chrM, unplaced contigs, rRNA/mitochondrial-RNA genes,
    IG/TR orphons, and duplicated sequence near the partner. A blacklisted molecule keeps its row in
    the junction table (status blacklist_*) and is not written to the GFF or h5ad."""
    if barcode_tags is None:
        barcode_tags = ['CB', 'CR', 'BC', 'BX']
    started = time.time()
    output_prefix = Path(output_prefix)
    gff_path = output_prefix.with_suffix('.gff')
    h5ad_path = output_prefix.with_suffix('.h5ad')
    table_path = Path(str(output_prefix) + '.cross_chromosome_junctions.tsv.gz')
    stats_path = Path(str(output_prefix) + '.cross_chromosome_stats.json')
    for path in (gff_path, Path(str(gff_path) + '.gz'), h5ad_path, table_path, stats_path):
        if path.exists():
            raise FileExistsError(f'Refusing to overwrite an existing output: {path}')
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    _set_pysam_verbosity()
    exon_coordinates, gene_data, _ = gff_process.importEnsemblGenes(exon_file)
    gene_region_counts = {gene: len(regions_) for gene, regions_ in gene_data.items()}
    blacklist_index = None
    if blacklist_members:
        from altanalyze3.components.bam import cross_chromosome_rules
        members = []
        with open(blacklist_members) as handle:
            next(handle)
            for line in handle:
                c, s, e, _gid, name, _gtype, reason = line.rstrip('\n').split('\t')
                members.append((c, int(s), int(e), (reason, name)))
        blacklist_index = cross_chromosome_rules.IntervalIndex(members)
    aligner = genome = None
    dup_cache = {}
    if genome_index:
        import mappy
        aligner = mappy.Aligner(fn_idx_in=str(genome_index), preset='splice', best_n=10)
        if not aligner:
            raise RuntimeError(f'Could not load the genome index {genome_index}')
        genome = pysam.FastaFile(str(genome_fasta))
    bam = pysam.AlignmentFile(bam_path, 'rb')
    lookup = pysam.AlignmentFile(bam_path, 'rb') if regions else None
    if regions:
        records = itertools.chain.from_iterable(bam.fetch(c, s, e) for c, s, e in regions)
    else:
        records = bam.fetch(until_eof=True)
    stats = defaultdict(int)
    rows, cols, data, barcodes, isoforms = [], [], [], [], []
    barcode_index, isoform_index, isoform_gene = {}, {}, {}
    seen_queries = set()
    with open(gff_path, 'w') as gff_handle, gzip.open(table_path, 'wt') as table:
        table.write('\t'.join(_JUNCTION_COLUMNS) + '\n')
        for read in records:
            stats['records_scanned'] += 1
            if read.is_unmapped or read.is_secondary:
                continue
            if read.is_supplementary:
                if not regions:
                    continue
                stats['region_supplementary_records'] += 1
                if not read.has_tag('SA'):
                    continue
                primary = _primary_for_supplementary(lookup, read)
                if primary is None:
                    stats['region_primary_not_found'] += 1
                    continue
                read = primary
            if regions:
                if read.query_name in seen_queries:
                    continue
                seen_queries.add(read.query_name)
            stats['primary_records'] += 1
            if not read.has_tag('SA'):
                continue
            stats['primary_records_with_SA'] += 1
            sa_entries = read.get_tag('SA').rstrip(';').split(';')
            if len({read.reference_name} | {e.split(',', 1)[0] for e in sa_entries}) < 2:
                stats['SA_same_chromosome_only'] += 1
                continue
            stats['cross_chromosome_reads'] += 1
            if regions:
                stats['SA_entries_checked'] += len(sa_entries)
                stats['SA_entries_matching_supplementary_record'] += _sa_records_match(lookup, read, sa_entries)
            segments = _alignment_segments(read, sa_entries, ts_read_relative, exon_coordinates)
            status = 'kept'
            primary_segment = next(s for s in segments if s['primary'])
            if primary_segment['strand'] != primary_segment['orientation']:
                status = 'strand_tag_differs_from_alignment'
            elif len({s['read_length'] for s in segments}) != 1:
                status = 'read_length_mismatch_between_alignments'
            barcode = resolve_barcode(read, barcode_tags, default_barcode=default_barcode)
            if status == 'kept' and not barcode:
                status = 'no_barcode'
            if status == 'kept' and min(s['mapq'] for s in segments) < min_mapq:
                status = 'mapq_below_min'
            untrimmed = [list(s['exons']) for s in segments]
            junctions = []
            for a, b in zip(segments[:-1], segments[1:]):
                junction = _segment_junction(a, b, exon_coordinates)
                junction['cross'] = a['chrom'] != b['chrom']
                junctions.append(junction)
            if status == 'kept' and require_known_splice and not any(
                    j['cross'] and (j['donor_hit'] is not None or j['acceptor_hit'] is not None)
                    for j in junctions):
                status = 'no_registered_boundary_at_cross_chromosome_joint'
            for (a, b), junction in zip(zip(segments[:-1], segments[1:]), junctions):
                _trim_to_junction(a, b, junction)
            membership = cross_chromosome_rules.segment_blacklist(segments, blacklist_index) if blacklist_index else None
            if membership:
                stats['blacklist_membership_all_reads:' + membership[0]] += 1
            if status == 'kept' and membership:
                status = membership[0]
            duplicated = ''
            if status == 'kept' and aligner is not None:
                duplicated = _duplicated_near_partner(segments, junctions, aligner, genome, dup_window,
                                                      dup_min_identity, dup_min_cover, dup_cache)
                if duplicated:
                    status = 'blacklist_duplicated_sequence_near_partner'
            gene, exon_ids, _, segment_genes = gff_process.annotate_cross_chromosome_segments(
                [(s['chrom'], s['strand'], s['exons']) for s in segments], read.query_name)
            _, raw_ids, _, _ = gff_process.annotate_cross_chromosome_segments(
                [(s['chrom'], s['strand'], e) for s, e in zip(segments, untrimmed)], read.query_name)
            molecule_id = sanitize_gff_value(resolve_molecule_id(read, molecule_tag))
            molecule_core = molecule_id.split('molecule/', 1)[1] if molecule_id.startswith('molecule/') else molecule_id
            isoform_id = sanitize_gff_value(f"{gene}:{molecule_prefix}{molecule_core}")
            if status == 'kept' and isoform_id in isoform_index:
                status = 'duplicate_isoform_id_in_library'
            stats['status:' + status] += 1
            if status == 'kept':
                _write_cross_chromosome_gff(gff_handle, segments, isoform_id, gene, source,
                                            first_in_file=not isoforms)
                rows.append(get_index(barcode, barcode_index, barcodes))
                cols.append(get_index(isoform_id, isoform_index, isoforms))
                data.append(1)
                isoform_gene[isoform_id] = gene
            sequence = _oriented_sequence(read, segments[0]['read_length'])
            for index, ((a, b), junction) in enumerate(zip(zip(segments[:-1], segments[1:]), junctions), 1):
                if junction['cross']:
                    stats['cross_chromosome_joints'] += 1
                d_ann, d_first, d_last = _boundary_annotation(junction['donor_hit'], gene_region_counts)
                a_ann, a_first, a_last = _boundary_annotation(junction['acceptor_hit'], gene_region_counts)
                context = ''
                if sequence:
                    left_end = a['qend'] - (junction['shift'] or 0) if junction['overlap'] else a['qend']
                    right_start = left_end if junction['overlap'] else b['qstart']
                    context = (sequence[max(0, left_end - 20):left_end] + '|' +
                               sequence[left_end:right_start].lower() + '|' + sequence[right_start:right_start + 20])
                values = [library, molecule_id, read.query_name, barcode or '',
                          read.get_tag('XM') if read.has_tag('XM') else '', status,
                          isoform_id if status == 'kept' else '', a['read_length'], len(segments), index,
                          junction['cross'], segment_genes[index - 1], segment_genes[index]]
                for s in (a, b):
                    values += [s['chrom'], s['strand'], s['qstart'], s['qend'], s['mapq'], s['nm'], s['primary'],
                               s['cigar'], ','.join(f'{x}-{y}' for x, y in s['exons'])]
                values += [junction['gap'], junction['overlap'], junction['shift'], junction['donor'],
                           junction['acceptor'], d_ann, d_first, d_last, a_ann, a_first, a_last, context,
                           '|'.join(raw_ids), f'{membership[0]}:{membership[1]}' if membership else '', duplicated]
                table.write('\t'.join(str(v) for v in values) + '\n')
    bam.close()
    if lookup is not None:
        lookup.close()
    stats['kept_molecules'] = len(isoforms)
    stats['kept_barcodes'] = len(barcodes)
    if isoforms:
        matrix = coo_matrix((data, (rows, cols)), shape=(len(barcodes), len(isoforms)), dtype=np.int32).tocsr()
        adata = ad.AnnData(X=matrix, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=isoforms))
        adata.var['gene'] = [isoform_gene.get(isoform, '') for isoform in isoforms]
        adata.write_h5ad(h5ad_path, compression='gzip')
    gff_gz = _gzip_file(gff_path)
    stats_out = dict(stats)
    stats_out.update(bam=str(bam_path), library=library, exon_file=str(exon_file), regions=regions,
                     min_mapq=min_mapq, molecule_tag=molecule_tag, barcode_tags=barcode_tags,
                     require_known_splice=require_known_splice, molecule_prefix=molecule_prefix,
                     blacklist_members=str(blacklist_members) if blacklist_members else None,
                     genome_index=str(genome_index) if genome_index else None,
                     dup_window=dup_window, dup_min_identity=dup_min_identity, dup_min_cover=dup_min_cover,
                     dup_segments_mapped=len(dup_cache), gff=str(gff_gz),
                     h5ad=str(h5ad_path) if isoforms else None, junction_table=str(table_path),
                     elapsed_seconds=round(time.time() - started, 1))
    stats_path.write_text(json.dumps(stats_out, indent=1, default=str) + '\n')
    return gff_gz, (h5ad_path if isoforms else None), stats_out


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
    parser.add_argument('--ts-read-relative', dest='ts_read_relative',
                        action='store_true',
                        help="Interpret the minimap2 'ts' tag as the transcript strand RELATIVE TO "
                             "THE READ (genomic strand = ts XOR alignment orientation). Off by "
                             "default, which keeps the tag verbatim.")
    parser.add_argument('--bulk-barcode', dest='bulk_barcode', default=None,
                        help='BULK long-read mode: barcode to assign when a read carries none of '
                             '--barcode-tags (bulk PacBio/ONT BAMs have no CB tag). Use the library '
                             'name; the whole BAM becomes one pseudo-cell. Omit for single-cell.')
    parser.add_argument('--no-parallel', dest='no_parallel', action='store_true',
                        help='Disable parallel chromosome/strand processing when no gene is provided.')
    parser.add_argument('--max-processes', dest='max_processes', type=int, default=None,
                        help='Maximum parallel workers (default: cpu_count-1).')
    parser.add_argument('--cross-chromosome-only', dest='cross_chromosome_only', action='store_true',
                        help='EXPERIMENTAL FORK: write only molecules aligned to two or more chromosomes '
                             '(primary + SA-tag alignments).')
    parser.add_argument('--library', dest='library', default=None,
                        help='Library name recorded in the cross-chromosome junction table.')
    parser.add_argument('--molecule-id-prefix', dest='molecule_prefix', default='',
                        help='Cross-chromosome mode: prefix added to molecule ids (unique across libraries).')
    parser.add_argument('--blacklist-members', dest='blacklist_members', default=None,
                        help='Cross-chromosome mode: reference/blacklist_members.tsv (enables the blacklist).')
    parser.add_argument('--genome-index', dest='genome_index', default=None,
                        help='Cross-chromosome mode: minimap2 splice index for the duplicated-sequence rule.')
    parser.add_argument('--genome-fasta', dest='genome_fasta', default=None,
                        help='Cross-chromosome mode: genome FASTA matching --genome-index.')
    parser.add_argument('--regions', dest='regions', default=None,
                        help='Cross-chromosome demo only: comma-separated chrom:start-end (1-based) windows.')
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
        if args.cross_chromosome_only:
            regions = None
            if args.regions:
                regions = []
                for window in args.regions.split(','):
                    chrom, span = window.rsplit(':', 1)
                    start, end = span.split('-')
                    regions.append((chrom, int(start) - 1, int(end)))
            gff_path, h5ad_path, stats = extract_cross_chromosome_isoforms(
                args.bam_file, args.gene_model, args.output_prefix,
                library=args.library or Path(args.output_prefix).name,
                min_mapq=args.min_mapq, barcode_tags=barcode_tags, molecule_tag=args.molecule_tag,
                source=args.source, require_known_splice=args.require_known_splice,
                default_barcode=args.bulk_barcode, ts_read_relative=args.ts_read_relative,
                regions=regions, molecule_prefix=args.molecule_prefix, blacklist_members=args.blacklist_members,
                genome_index=args.genome_index, genome_fasta=args.genome_fasta)
            print("Processing mode: EXPERIMENTAL cross-chromosome only (single process)")
            for key in sorted(stats):
                print(f"{key}: {stats[key]}")
        elif args.gene is None:
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
                    default_barcode=args.bulk_barcode,
                    ts_read_relative=args.ts_read_relative,
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
                    default_barcode=args.bulk_barcode,
                    ts_read_relative=args.ts_read_relative,
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
                default_barcode=args.bulk_barcode,
                ts_read_relative=args.ts_read_relative,
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
