#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from datetime import datetime
from collections import defaultdict, OrderedDict
import sqlite3

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['figure.facecolor'] = 'white'

sys.path.insert(1, os.path.join(os.path.dirname(__file__), '..'))
from long_read import isoform_matrix as iso

GENE_ID_PREFIXES = ('ENS',)
INDEX_DIR_NAME = "gene_indexes_v2"
_GROUPBY_CACHE = {}
_GROUPBY_ROWS_CACHE = {}
_H5AD_GENE_INDEX_CACHE = {}
_H5AD_COUNTS_CACHE = {}


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


def _init_logging(output_path):
    output_path = Path(output_path)
    output_dir = output_path.parent if output_path.suffix else output_path
    if not output_dir.as_posix():
        output_dir = Path('.')
    log_dir = output_dir / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_name = f"isoform_structure_view_{timestamp}.log"
    log_path = log_dir / log_name
    log_handle = open(log_path, 'w')
    stdout_tee = _Tee(sys.stdout, log_handle)
    stderr_tee = _Tee(sys.stderr, log_handle)
    return log_path, log_handle, stdout_tee, stderr_tee


def _normalize_group_value(value):
    if value is None:
        return set()
    if isinstance(value, (list, tuple, set)):
        values = value
    else:
        values = str(value).split(',')
    return {item.strip() for item in values if str(item).strip()}


def _matches_group(groups_value, target_group):
    if not target_group:
        return True
    return target_group in _normalize_group_value(groups_value)


def detect_gene_index(path):
    if os.path.isdir(path):
        if os.path.exists(os.path.join(path, "gene_model.db")) or os.path.exists(os.path.join(path, "transcripts.db")):
            return path
    base_dir = os.path.dirname(path)
    index_dir = os.path.join(base_dir, INDEX_DIR_NAME)
    if os.path.isdir(index_dir):
        if os.path.exists(os.path.join(index_dir, "gene_model.db")) or os.path.exists(os.path.join(index_dir, "transcripts.db")):
            return index_dir
    return None


def _is_index_stale(index_path, source_paths):
    if not os.path.exists(index_path):
        return True
    try:
        index_mtime = os.path.getmtime(index_path)
    except OSError:
        return True
    for source in source_paths:
        if not source:
            continue
        try:
            if os.path.getmtime(source) > index_mtime:
                return True
        except OSError:
            continue
    return False


def _transcript_index_dir(transcript_path):
    base_dir = os.path.dirname(transcript_path)
    index_dir = os.path.join(base_dir, INDEX_DIR_NAME)
    os.makedirs(index_dir, exist_ok=True)
    return index_dir


def _gene_model_index_dir(gene_model_path):
    base_dir = os.path.dirname(gene_model_path)
    index_dir = os.path.join(base_dir, INDEX_DIR_NAME)
    os.makedirs(index_dir, exist_ok=True)
    return index_dir


def _h5ad_gene_index_path(h5ad_path, _index_dir=None):
    base_dir = os.path.dirname(h5ad_path)
    index_dir = os.path.join(base_dir, INDEX_DIR_NAME)
    os.makedirs(index_dir, exist_ok=True)
    return os.path.join(index_dir, os.path.basename(h5ad_path) + ".gene_index.npz")

def _counts_cache_key(groupby, group_values, groupby_rev, groupby_sample):
    payload = {
        "groupby": groupby or "",
        "group_values": sorted(group_values) if group_values else [],
        "groupby_rev": bool(groupby_rev),
        "groupby_sample": groupby_sample or "",
    }
    digest = hashlib.md5(
        json.dumps(payload, sort_keys=True).encode("utf-8")
    ).hexdigest()[:12]
    return digest


def _h5ad_counts_cache_path(h5ad_path, groupby, group_values, groupby_rev, groupby_sample):
    base_dir = os.path.dirname(h5ad_path)
    index_dir = os.path.join(base_dir, INDEX_DIR_NAME)
    os.makedirs(index_dir, exist_ok=True)
    suffix = _counts_cache_key(groupby, group_values, groupby_rev, groupby_sample)
    return os.path.join(
        index_dir,
        os.path.basename(h5ad_path) + f".isoform_counts.{suffix}"
    )


def _load_isoform_counts_cache(cache_path):
    cached = _H5AD_COUNTS_CACHE.get(cache_path)
    if cached is not None:
        return cached
    counts_path = cache_path + ".counts.npy"
    var_path = cache_path + ".var_names.npy"
    if os.path.exists(counts_path) and os.path.exists(var_path):
        payload = {
            "counts": np.load(counts_path, mmap_mode='r'),
            "var_names": np.load(var_path, mmap_mode='r'),
        }
        _H5AD_COUNTS_CACHE[cache_path] = payload
        return payload
    legacy_path = cache_path + ".npz"
    if os.path.exists(legacy_path):
        try:
            cache = np.load(legacy_path, allow_pickle=True)
        except ValueError as exc:
            print(f"[warn] Failed to load isoform counts cache: {legacy_path} ({exc})")
            return None
        payload = {
            "counts": np.asarray(cache["counts"]),
            "var_names": np.asarray(cache["var_names"], dtype=str),
        }
        _H5AD_COUNTS_CACHE[cache_path] = payload
        return payload
    return None


def _build_isoform_counts_cache(adata, mask, cache_path):
    t_counts = time.time()
    if mask is not None:
        sub = adata[mask, :]
    else:
        sub = adata[:, :]
    counts = np.asarray(sub.X.sum(axis=0)).ravel().astype(np.float32)
    var_names = np.asarray(sub.var_names)
    var_names = np.asarray(var_names, dtype=str)
    max_len = max((len(name) for name in var_names), default=1)
    var_names_bytes = np.asarray(var_names, dtype=f"S{max_len}")
    counts_path = cache_path + ".counts.npy"
    var_path = cache_path + ".var_names.npy"
    np.save(counts_path, counts)
    np.save(var_path, var_names_bytes)
    _H5AD_COUNTS_CACHE[cache_path] = {
        "counts": counts,
        "var_names": var_names_bytes,
    }
    print(f"[index] isoform counts cache complete: {cache_path}")
    print(f"[timing] isoform counts cache build: {time.time() - t_counts:.2f}s")


def _records_from_counts_cache(cache, gene_idx, target_gene):
    counts = cache["counts"][gene_idx]
    var_names = cache["var_names"][gene_idx]
    if len(gene_idx) == 0:
        return []
    records = []
    for idx, var in enumerate(var_names):
        if isinstance(var, (bytes, np.bytes_)):
            var = var.decode()
        else:
            var = str(var)
        gene_id, isoform_id = parse_feature_identifier(var, gene_hint=target_gene)
        if gene_id is not None and gene_id != target_gene:
            continue
        records.append({
            'var': var,
            'gene_id': gene_id or target_gene,
            'isoform_id': isoform_id,
            'count': counts[idx]
        })
    return records


def _load_or_build_counts_cache(h5ad_path, adata, mask, groupby, group_values,
                                groupby_rev, groupby_sample, use_cache=True,
                                allow_build=True):
    if not use_cache:
        return None
    cache_path = _h5ad_counts_cache_path(
        h5ad_path, groupby, group_values, groupby_rev, groupby_sample
    )
    counts_path = cache_path + ".counts.npy"
    var_path = cache_path + ".var_names.npy"
    legacy_path = cache_path + ".npz"
    source_paths = [h5ad_path]
    if groupby and (groupby.endswith('.tsv') or groupby.endswith('.txt')):
        source_paths.append(groupby)
    if (os.path.exists(counts_path) and os.path.exists(var_path)
            and not _is_index_stale(counts_path, source_paths)):
        print(f"[index] Using isoform counts cache: {cache_path}")
        cached = _load_isoform_counts_cache(cache_path)
        if cached is not None:
            return cached
    if os.path.exists(legacy_path) and not _is_index_stale(legacy_path, source_paths):
        if allow_build and adata is not None:
            print(f"[index] Rebuilding legacy isoform counts cache: {legacy_path}")
            _build_isoform_counts_cache(adata, mask, cache_path)
            return _load_isoform_counts_cache(cache_path)
        print(f"[index] Using isoform counts cache: {cache_path}")
        cached = _load_isoform_counts_cache(cache_path)
        if cached is not None:
            return cached
    if not allow_build:
        return None
    if adata is None:
        return None
    _build_isoform_counts_cache(adata, mask, cache_path)
    return _load_isoform_counts_cache(cache_path)

def _decode_h5py_array(values):
    arr = np.asarray(values)
    if arr.dtype.kind in ('S', 'O'):
        decoded = []
        for value in arr:
            if isinstance(value, bytes):
                decoded.append(value.decode())
            else:
                decoded.append(str(value))
        return np.asarray(decoded, dtype=object)
    return arr.astype(str)


def _read_h5ad_gene_values_fast(h5ad_path):
    try:
        import h5py
    except ImportError:
        return None, None, "h5py-unavailable"
    try:
        with h5py.File(h5ad_path, 'r') as handle:
            var_group = handle.get('var')
            if var_group is None:
                return None, None, "missing-var-group"
            var_names = None
            if '_index' in var_group:
                var_names = _decode_h5py_array(var_group['_index'][()])
            gene_values = None
            if 'gene' in var_group:
                node = var_group['gene']
                if hasattr(node, 'shape') and len(node.shape) == 1:
                    gene_values = _decode_h5py_array(node[()])
                elif 'categories' in node and 'codes' in node:
                    categories = _decode_h5py_array(node['categories'][()])
                    codes = np.asarray(node['codes'][()])
                    gene_values = np.empty(len(codes), dtype=object)
                    for idx, code in enumerate(codes):
                        if code < 0:
                            gene_values[idx] = ''
                        else:
                            gene_values[idx] = categories[code]
            return gene_values, var_names, "h5py"
    except OSError:
        return None, None, "h5py-error"


def _build_h5ad_gene_index(h5ad_path, index_path):
    print(f"[index] Building h5ad gene index: {index_path}")
    gene_values, var_names, source = _read_h5ad_gene_values_fast(h5ad_path)
    if gene_values is None:
        adata = ad.read_h5ad(h5ad_path, backed='r')
        var_names = np.asarray(adata.var_names)
        if 'gene' in adata.var.columns:
            gene_values = np.asarray(adata.var['gene'])
        else:
            gene_values = np.array([
                name.split(':', 1)[0] if ':' in name else '' for name in var_names
            ])
        try:
            adata.file.close()
        except AttributeError:
            pass
        source = "anndata"
    elif gene_values is None and var_names is not None:
        gene_values = np.array([
            name.split(':', 1)[0] if ':' in name else '' for name in var_names
        ])
    print(f"[index] h5ad index source: {source}")
    gene_to_indices = defaultdict(list)
    for idx, gene in enumerate(gene_values):
        if not gene:
            continue
        gene_to_indices[str(gene)].append(idx)
    genes = np.array(sorted(gene_to_indices.keys()))
    indptr = np.zeros(len(genes) + 1, dtype=np.int64)
    indices = []
    for i, gene in enumerate(genes):
        idx_list = gene_to_indices[gene]
        indptr[i + 1] = indptr[i] + len(idx_list)
        indices.extend(idx_list)
    indices = np.asarray(indices, dtype=np.int64)
    np.savez_compressed(index_path, genes=genes, indptr=indptr, indices=indices)
    print(f"[index] h5ad gene index complete: {index_path}")


def _load_h5ad_gene_index(h5ad_path, index_dir=None):
    index_path = _h5ad_gene_index_path(h5ad_path, _index_dir=index_dir)
    cache_key = os.path.abspath(index_path)
    if cache_key in _H5AD_GENE_INDEX_CACHE:
        print(f"[index] Using cached h5ad gene index: {index_path}")
        return _H5AD_GENE_INDEX_CACHE[cache_key]
    if _is_index_stale(index_path, [h5ad_path]):
        _build_h5ad_gene_index(h5ad_path, index_path)
    else:
        print(f"[index] Using existing h5ad gene index: {index_path}")
    data = np.load(index_path, allow_pickle=False)
    genes = data['genes'].astype(str)
    indptr = data['indptr']
    indices = data['indices']
    gene_to_pos = {gene: idx for idx, gene in enumerate(genes)}
    payload = (gene_to_pos, indptr, indices)
    _H5AD_GENE_INDEX_CACHE[cache_key] = payload
    return payload


def looks_like_gene_id(value):
    if not value:
        return False
    if not value.startswith(GENE_ID_PREFIXES):
        return False
    # Basic heuristic for Ensembl-like IDs.
    return any(char.isdigit() for char in value)


def parse_feature_identifier(feature, gene_hint=None):
    if ':' in feature:
        left, right = feature.split(':', 1)
        if looks_like_gene_id(left):
            return left, right
        if looks_like_gene_id(right):
            return right, left
        return (gene_hint if looks_like_gene_id(gene_hint) else None), left
    return (gene_hint if looks_like_gene_id(gene_hint) else None), feature


def read_gene_model(path):
    gene_segments = defaultdict(list)
    exon_lookup = {}

    with open(path, 'r') as handle:
        first = handle.readline()
        if not first:
            return gene_segments, exon_lookup
        header = first.rstrip('\n').split('\t')
        has_header = header and header[0] == 'gene'
        if has_header:
            header_map = {name: idx for idx, name in enumerate(header)}
            idx_gene = header_map.get('gene', 0)
            idx_exon = header_map.get('exon-id', 1)
            idx_chr = header_map.get('chromosome', 2)
            idx_strand = header_map.get('strand', 3)
            idx_start = header_map.get('exon-region-start(s)', 4)
            idx_end = header_map.get('exon-region-stop(s)', 5)
        else:
            idx_gene, idx_exon, idx_chr, idx_strand, idx_start, idx_end = 0, 1, 2, 3, 4, 5
            parts = header
            _ingest_gene_model_row(parts, gene_segments, exon_lookup,
                                   idx_gene, idx_exon, idx_chr, idx_strand, idx_start, idx_end)

        for line in handle:
            parts = line.rstrip('\n').split('\t')
            _ingest_gene_model_row(parts, gene_segments, exon_lookup,
                                   idx_gene, idx_exon, idx_chr, idx_strand, idx_start, idx_end)

    return gene_segments, exon_lookup


def _build_gene_model_db(gene_model_path, db_path):
    print(f"[index] Building gene model index: {db_path}")
    conn = sqlite3.connect(db_path)
    try:
        cursor = conn.cursor()
        cursor.execute("PRAGMA journal_mode=OFF")
        cursor.execute("PRAGMA synchronous=OFF")
        cursor.execute("PRAGMA temp_store=MEMORY")
        cursor.execute("PRAGMA locking_mode=EXCLUSIVE")
        cursor.execute("DROP TABLE IF EXISTS gene_model")
        cursor.execute(
            "CREATE TABLE gene_model ("
            "gene TEXT, exon_id TEXT, chrom TEXT, strand TEXT, "
            "start INTEGER, end INTEGER, const TEXT, ens TEXT, "
            "splice_event TEXT, splice_junction TEXT)"
        )
        cursor.execute("CREATE INDEX idx_gene_model_gene ON gene_model(gene)")
        buffer = []
        buffer_limit = 50000
        with open(gene_model_path, 'r') as handle:
            first = handle.readline()
            if not first:
                return
            header = first.rstrip('\n').split('\t')
            has_header = header and header[0] == 'gene'
            if has_header:
                header_map = {name: idx for idx, name in enumerate(header)}
                idx_gene = header_map.get('gene')
                idx_exon = header_map.get('exon-id')
                idx_chr = header_map.get('chromosome')
                idx_strand = header_map.get('strand')
                idx_start = header_map.get('exon-region-start(s)')
                idx_end = header_map.get('exon-region-stop(s)')
                idx_const = header_map.get('constitutive')
                idx_ens = header_map.get('Ensembl-exon')
                idx_splice_event = header_map.get('splice-event')
                idx_splice_junc = header_map.get('splice-junction')
            else:
                idx_gene, idx_exon, idx_chr, idx_strand, idx_start, idx_end = 0, 1, 2, 3, 4, 5
                idx_const = idx_ens = idx_splice_event = idx_splice_junc = None
                parts = header
                if len(parts) > idx_end:
                    try:
                        start = int(float(parts[idx_start]))
                        end = int(float(parts[idx_end]))
                    except ValueError:
                        start = end = None
                    if start is not None and end is not None:
                        gene = parts[idx_gene].strip()
                        exon_id = parts[idx_exon].strip()
                        chrom = parts[idx_chr].strip()
                        strand = parts[idx_strand].strip()
                        const = parts[idx_const].strip() if idx_const is not None and idx_const < len(parts) else ''
                        ens = parts[idx_ens].strip() if idx_ens is not None and idx_ens < len(parts) else ''
                        splice_event = parts[idx_splice_event].strip() if idx_splice_event is not None and idx_splice_event < len(parts) else ''
                        splice_junc = parts[idx_splice_junc].strip() if idx_splice_junc is not None and idx_splice_junc < len(parts) else ''
                        start, end = (start, end) if start <= end else (end, start)
                        buffer.append((gene, exon_id, chrom, strand, start, end, const, ens, splice_event, splice_junc))
            for line in handle:
                parts = line.rstrip('\n').split('\t')
                if len(parts) <= idx_end:
                    continue
                gene = parts[idx_gene].strip()
                exon_id = parts[idx_exon].strip()
                chrom = parts[idx_chr].strip()
                strand = parts[idx_strand].strip()
                try:
                    start = int(float(parts[idx_start]))
                    end = int(float(parts[idx_end]))
                except ValueError:
                    continue
                if not gene or not exon_id:
                    continue
                start, end = (start, end) if start <= end else (end, start)
                const = parts[idx_const].strip() if idx_const is not None and idx_const < len(parts) else ''
                ens = parts[idx_ens].strip() if idx_ens is not None and idx_ens < len(parts) else ''
                splice_event = parts[idx_splice_event].strip() if idx_splice_event is not None and idx_splice_event < len(parts) else ''
                splice_junc = parts[idx_splice_junc].strip() if idx_splice_junc is not None and idx_splice_junc < len(parts) else ''
                buffer.append((gene, exon_id, chrom, strand, start, end, const, ens, splice_event, splice_junc))
                if len(buffer) >= buffer_limit:
                    cursor.executemany(
                        "INSERT INTO gene_model VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                        buffer
                    )
                    buffer.clear()
        if buffer:
            cursor.executemany(
                "INSERT INTO gene_model VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                buffer
            )
        conn.commit()
    finally:
        conn.close()
    print(f"[index] Gene model index complete: {db_path}")


def _ensure_gene_model_db(gene_model_path, index_dir=None):
    if not gene_model_path or not os.path.exists(gene_model_path):
        return None
    index_dir = index_dir or _gene_model_index_dir(gene_model_path)
    db_path = os.path.join(index_dir, "gene_model.db")
    if _is_index_stale(db_path, [gene_model_path]):
        _build_gene_model_db(gene_model_path, db_path)
    return db_path


def read_gene_model_from_sqlite(db_path, target_gene):
    gene_segments = defaultdict(list)
    exon_lookup = {}
    if not os.path.exists(db_path):
        return gene_segments, exon_lookup
    conn = sqlite3.connect(db_path)
    try:
        cursor = conn.cursor()
        cursor.execute(
            'SELECT * FROM gene_model WHERE gene = ?',
            (target_gene,)
        )
        for row in cursor.fetchall():
            gene, exon_id, chrom, strand, start, end, _const, _ens, _splice_ev, _splice_junc = row
            try:
                start = int(start)
                end = int(end)
            except (ValueError, TypeError):
                continue
            if not gene or not exon_id:
                continue
            start, end = (start, end) if start <= end else (end, start)
            feature_type = exon_id[0].upper()
            segment = {
                'gene': gene,
                'exon_id': exon_id,
                'chrom': chrom,
                'strand': strand,
                'start': start,
                'end': end,
                'type': feature_type
            }
            gene_segments[gene].append(segment)
            exon_lookup[(gene, exon_id)] = segment
    finally:
        conn.close()
    return gene_segments, exon_lookup


def load_transcript_associations_from_sqlite(db_path, target_gene):
    t_start = time.time()
    structures = {}
    gene_strand = {}
    if not os.path.exists(db_path):
        return structures, gene_strand
    conn = sqlite3.connect(db_path)
    try:
        cursor = conn.cursor()
        cursor.execute(
            'SELECT * FROM transcripts WHERE gene = ?',
            (target_gene,)
        )
        for row in cursor.fetchall():
            gene, strand, exon_sequence, transcript_id, _sample = row
            tokens = [t for t in exon_sequence.split('|') if t]
            if not tokens:
                continue
            structures[transcript_id] = tokens
            for alias in _iter_alias_tokens(transcript_id):
                if alias not in structures:
                    structures[alias] = tokens
            gene_strand[gene] = strand
    finally:
        conn.close()
    print(f"[timing] transcript_associations sqlite ({os.path.basename(db_path)}): {time.time() - t_start:.2f}s")
    return structures, gene_strand


def _build_transcripts_db(transcript_path, db_path):
    print(f"[index] Building transcript index: {db_path}")
    conn = sqlite3.connect(db_path)
    try:
        cursor = conn.cursor()
        cursor.execute("PRAGMA journal_mode=OFF")
        cursor.execute("PRAGMA synchronous=OFF")
        cursor.execute("PRAGMA temp_store=MEMORY")
        cursor.execute("PRAGMA locking_mode=EXCLUSIVE")
        cursor.execute("DROP TABLE IF EXISTS transcripts")
        cursor.execute(
            "CREATE TABLE transcripts ("
            "gene TEXT, strand TEXT, exon_sequence TEXT, "
            "transcript_id TEXT, sample TEXT)"
        )
        cursor.execute("CREATE INDEX idx_transcripts_gene ON transcripts(gene)")
        buffer = []
        buffer_limit = 50000
        with open(transcript_path, 'r') as handle:
            for line in handle:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 4:
                    continue
                gene, strand, exon_sequence, transcript_id = parts[:4]
                sample = parts[4] if len(parts) > 4 else ''
                buffer.append((gene, strand, exon_sequence, transcript_id, sample))
                if len(buffer) >= buffer_limit:
                    cursor.executemany(
                        "INSERT INTO transcripts VALUES (?, ?, ?, ?, ?)",
                        buffer
                    )
                    buffer.clear()
        if buffer:
            cursor.executemany(
                "INSERT INTO transcripts VALUES (?, ?, ?, ?, ?)",
                buffer
            )
        conn.commit()
    finally:
        conn.close()
    print(f"[index] Transcript index complete: {db_path}")


def _ensure_transcripts_db(transcript_path, index_dir=None):
    if not transcript_path or not os.path.exists(transcript_path):
        return None
    index_dir = index_dir or _transcript_index_dir(transcript_path)
    db_path = os.path.join(index_dir, "transcripts.db")
    if _is_index_stale(db_path, [transcript_path]):
        _build_transcripts_db(transcript_path, db_path)
    else:
        print(f"[index] Using existing transcript index: {db_path}")
    return db_path


def _delete_index(path, label):
    if not path or not os.path.exists(path):
        return
    try:
        os.remove(path)
        print(f"[index] Deleted {label}: {path}")
    except OSError as exc:
        print(f"[warn] Failed to delete {label} {path}: {exc}")


def _maybe_rebuild_indexes(rebuild, h5ad_paths, transcript_paths, gene_model_path):
    if not rebuild:
        return
    if gene_model_path:
        gene_dir = _gene_model_index_dir(gene_model_path)
        _delete_index(os.path.join(gene_dir, "gene_model.db"), "gene model index")
    for assoc_path in transcript_paths or []:
        trans_dir = _transcript_index_dir(assoc_path)
        _delete_index(os.path.join(trans_dir, "transcripts.db"), "transcript index")
    for h5ad_path in h5ad_paths or []:
        index_path = _h5ad_gene_index_path(h5ad_path)
        _delete_index(index_path, "h5ad gene index")
        _H5AD_GENE_INDEX_CACHE.pop(os.path.abspath(index_path), None)


def _ingest_gene_model_row(parts, gene_segments, exon_lookup,
                           idx_gene, idx_exon, idx_chr, idx_strand, idx_start, idx_end):
    if len(parts) <= max(idx_gene, idx_exon, idx_chr, idx_strand, idx_start, idx_end):
        return
    gene = parts[idx_gene].strip()
    exon_id = parts[idx_exon].strip()
    chrom = parts[idx_chr].strip()
    strand = parts[idx_strand].strip()
    try:
        start = int(float(parts[idx_start]))
        end = int(float(parts[idx_end]))
    except ValueError:
        return
    if not gene or not exon_id:
        return
    start, end = (start, end) if start <= end else (end, start)
    feature_type = exon_id[0].upper()
    segment = {
        'gene': gene,
        'exon_id': exon_id,
        'chrom': chrom,
        'strand': strand,
        'start': start,
        'end': end,
        'type': feature_type
    }
    gene_segments[gene].append(segment)
    exon_lookup[(gene, exon_id)] = segment


def build_gene_maps(gene_segments, intron_scale):
    gene_maps = {}
    for gene, segments in gene_segments.items():
        ordered = sorted(segments, key=lambda s: s['start'])
        display_pos = 0.0
        mapped_segments = []
        for seg in ordered:
            length = max(1, abs(seg['end'] - seg['start']) + 1)
            scale = intron_scale if seg['type'] == 'I' else 1.0
            display_start = display_pos
            display_end = display_pos + (length * scale)
            mapped_segments.append({
                **seg,
                'display_start': display_start,
                'display_end': display_end,
                'scale': scale
            })
            display_pos = display_end
        gene_maps[gene] = {
            'segments': mapped_segments,
            'total_length': display_pos
        }
    return gene_maps


def map_coord(gene_map, coord):
    segments = gene_map['segments']
    if not segments:
        return coord
    for seg in segments:
        if seg['start'] <= coord <= seg['end']:
            return seg['display_start'] + (coord - seg['start']) * seg['scale']
    if coord < segments[0]['start']:
        return segments[0]['display_start']
    return segments[-1]['display_end']


def _iter_alias_tokens(value, separators=('|', ':')):
    pending = [value]
    seen = set()
    while pending:
        token = pending.pop()
        if token in seen:
            continue
        seen.add(token)
        if not token:
            continue
        for sep in separators:
            if sep in token:
                pending.extend([part for part in token.split(sep) if part])
    return [t for t in seen if t]


def load_transcript_associations(path, target_gene):
    t_start = time.time()
    structures = {}
    gene_strand = {}
    with open(path, 'r') as handle:
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            gene, strand, exon_sequence, transcript_id = parts[:4]
            if gene != target_gene:
                continue
            tokens = [t for t in exon_sequence.split('|') if t]
            if not tokens:
                continue
            structures[transcript_id] = tokens
            for alias in _iter_alias_tokens(transcript_id):
                if alias not in structures:
                    structures[alias] = tokens
            gene_strand[gene] = strand
    print(f"[timing] transcript_associations tsv ({os.path.basename(path)}): {time.time() - t_start:.2f}s")
    return structures, gene_strand


def load_transcript_associations_auto(path, target_gene, index_dir=None):
    if path and path.endswith('.db') and os.path.exists(path):
        return load_transcript_associations_from_sqlite(path, target_gene)
    if path:
        trans_db = _ensure_transcripts_db(path, index_dir=index_dir)
        if trans_db and os.path.exists(trans_db):
            print(f"[index] Using transcripts.db: {trans_db}")
            structures, gene_strand = load_transcript_associations_from_sqlite(trans_db, target_gene)
            if structures:
                return structures, gene_strand
            print(
                f"[warn] No transcript structures in transcripts.db for {target_gene}; "
                f"falling back to TSV: {path}"
            )
        return load_transcript_associations(path, target_gene)
    raise FileNotFoundError("transcript_associations path is missing or invalid")


def sample_transcript_ids_from_associations(path, target_gene, limit=5):
    samples = []
    with open(path, 'r') as handle:
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            gene, _, _, transcript_id = parts[:4]
            if gene != target_gene:
                continue
            samples.append(transcript_id)
            if len(samples) >= limit:
                break
    return samples


def load_transcript_association_index(path, target_gene):
    transcript_ids = set()
    alias_to_transcript = {}
    with open(path, 'r') as handle:
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            gene, _, _, transcript_id = parts[:4]
            if gene != target_gene:
                continue
            transcript_ids.add(transcript_id)
            for alias in _iter_alias_tokens(transcript_id):
                if alias not in alias_to_transcript:
                    alias_to_transcript[alias] = transcript_id
    return transcript_ids, alias_to_transcript


def split_token(token, default_gene):
    gene_id = default_gene
    core = token
    if ':' in token:
        gene_id, core = token.split(':', 1)
    coord = None
    base = core
    if '_' in core:
        base, coord_str = core.rsplit('_', 1)
        try:
            coord = int(coord_str)
        except ValueError:
            base = core
            coord = None
    return gene_id, base, coord


def build_isoform_segments(tokens, exon_lookup, default_gene):
    parsed = []
    for token in tokens:
        gene_id, base, coord = split_token(token, default_gene)
        parsed.append({'gene': gene_id, 'base': base, 'coord': coord})

    grouped = []
    for item in parsed:
        if grouped and grouped[-1]['gene'] == item['gene'] and grouped[-1]['base'] == item['base']:
            if item['coord'] is not None:
                grouped[-1]['coords'].append(item['coord'])
        else:
            grouped.append({
                'gene': item['gene'],
                'base': item['base'],
                'coords': [item['coord']] if item['coord'] is not None else []
            })

    direction = None
    for i in range(len(grouped) - 1):
        seg1 = exon_lookup.get((grouped[i]['gene'], grouped[i]['base']))
        seg2 = exon_lookup.get((grouped[i + 1]['gene'], grouped[i + 1]['base']))
        if seg1 and seg2:
            if seg2['start'] > seg1['start']:
                direction = 1
            elif seg2['start'] < seg1['start']:
                direction = -1
            break
    if direction is None:
        direction = 1

    segments = []
    for idx, group in enumerate(grouped):
        seg = exon_lookup.get((group['gene'], group['base']))
        if not seg:
            continue
        start, end = seg['start'], seg['end']
        if start > end:
            start, end = end, start
        coords = [c for c in group['coords'] if c is not None]
        if coords:
            coord_min = min(coords)
            coord_max = max(coords)
            if idx == 0:
                if direction == 1:
                    start = max(start, coord_min)
                else:
                    end = min(end, coord_max)
            if idx == len(grouped) - 1:
                if direction == 1:
                    end = min(end, coord_max)
                else:
                    start = max(start, coord_min)
            if idx != 0 and idx != len(grouped) - 1 and len(coords) > 1:
                start = max(start, coord_min)
                end = min(end, coord_max)
        if end < start:
            start, end = end, start
        label = None
        if seg['type'] == 'E':
            label = normalize_token(group['base'], default_gene, 'block')
        segments.append({
            'gene': group['gene'],
            'start': start,
            'end': end,
            'type': seg['type'],
            'label': label
        })
    return segments


def normalize_token(token, default_gene, mode):
    gene_id, base, _ = split_token(token, default_gene)
    core = base
    if mode == 'full':
        return f"{gene_id}:{base}" if gene_id != default_gene else base
    if mode == 'base':
        return core
    if mode == 'block':
        if len(core) > 1 and '.' in core:
            prefix = core[0]
            number = core[1:].split('.', 1)[0]
            return f"{prefix}{number}"
        return core
    return core


def resolve_isoform_id(raw_id, transcript_structures):
    if raw_id in transcript_structures:
        return raw_id
    for candidate in _iter_alias_tokens(raw_id):
        if candidate in transcript_structures:
            return candidate
    if '_' in raw_id:
        tail = raw_id.split('_')[-1]
        if tail in transcript_structures:
            return tail
        head = raw_id.split('_')[0]
        if head in transcript_structures:
            return head
    return None


def resolve_isoform_id_fuzzy(raw_id, transcript_structures):
    resolved = resolve_isoform_id(raw_id, transcript_structures)
    if resolved:
        return resolved, []
    suffix_matches = [key for key in transcript_structures if key.endswith(raw_id)]
    if len(suffix_matches) == 1:
        return suffix_matches[0], []
    if not suffix_matches:
        contains_matches = [key for key in transcript_structures if raw_id in key]
        if len(contains_matches) == 1:
            return contains_matches[0], []
        return None, contains_matches
    return None, suffix_matches


def _load_adata(input_path):
    if input_path.endswith('.h5ad'):
        return ad.read_h5ad(input_path)
    if input_path.endswith('.mtx') or input_path.endswith('.mtx.gz') or os.path.isdir(input_path):
        mtx_dir = input_path if os.path.isdir(input_path) else os.path.dirname(input_path)
        feature = _find_mtx_file(mtx_dir, ['genes.tsv', 'features.tsv', 'genes.tsv.gz', 'features.tsv.gz'])
        barcode = _find_mtx_file(mtx_dir, ['barcodes.tsv', 'barcodes.tsv.gz'])
        matrix = _find_mtx_file(mtx_dir, ['matrix.mtx', 'matrix.mtx.gz'])
        if not feature or not barcode or not matrix:
            raise ValueError(f"Missing 10x mtx files in {mtx_dir}.")
        return iso.mtx_to_adata(
            int_folder=mtx_dir,
            gene_is_index=True,
            feature=os.path.basename(feature),
            feature_col=0,
            barcode=os.path.basename(barcode),
            barcode_col=0,
            matrix=os.path.basename(matrix),
            rev=False
        )
    raise ValueError(f"Unsupported input format: {input_path}")


def _infer_sample_name(input_path):
    base = os.path.basename(input_path)
    if base.endswith('.h5ad'):
        return base[:-5]
    if base.endswith('.h5ad.gz'):
        return base[:-8]
    return os.path.splitext(base)[0]


def _find_mtx_file(folder, candidates):
    for name in candidates:
        path = os.path.join(folder, name)
        if os.path.exists(path):
            return path
    return None


def _normalize_barcode(value):
    token = value.split('.', 1)[0]
    return token[:-2] if token.endswith('-1') else token


def _reverse_complement_barcode(seq):
    if '-' in seq:
        return iso.reverse_complement_seq(seq)
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def _split_barcode_sample(value):
    parts = value.split('.', 1)
    barcode = parts[0]
    sample = parts[1] if len(parts) > 1 else ''
    return barcode, sample


def _load_groupby_file(path, reverse_complement=True, sample_name=None):
    cache_key = (path, reverse_complement, sample_name)
    if cache_key in _GROUPBY_CACHE:
        return _GROUPBY_CACHE[cache_key]
    if path not in _GROUPBY_ROWS_CACHE:
        rows_all = []
        rows_by_sample = defaultdict(list)
        with open(path, newline='') as handle:
            reader = csv.reader(handle, delimiter='\t')
            for row in reader:
                if len(row) < 2:
                    continue
                barcode, group = row[0], row[1]
                if not barcode or not group:
                    continue
                rows_all.append((barcode, group))
                if '.' in barcode:
                    rows_by_sample[barcode.rsplit('.', 1)[1]].append((barcode, group))
        _GROUPBY_ROWS_CACHE[path] = (rows_all, rows_by_sample)
    rows_all, rows_by_sample = _GROUPBY_ROWS_CACHE[path]
    if sample_name:
        rows = rows_by_sample.get(str(sample_name), [])
    else:
        rows = rows_all
    if not rows:
        _GROUPBY_CACHE[cache_key] = {}
        return {}
    barcodes = [_normalize_barcode(row[0]) for row in rows]
    if reverse_complement:
        barcodes = [_reverse_complement_barcode(value) for value in barcodes]
    mapping = dict(zip(barcodes, [row[1] for row in rows]))
    _GROUPBY_CACHE[cache_key] = mapping
    return mapping


def _load_barcode_mask(adata, barcode_series, cell_states=None):
    if barcode_series is None:
        return np.ones(adata.n_obs, dtype=bool)
    if cell_states:
        barcode_series = barcode_series.loc[barcode_series.isin(cell_states)]
    barcodes = set(barcode_series.index.astype(str))
    obs_names = pd.Series(adata.obs_names.astype(str))
    mask = obs_names.isin(barcodes).values
    if mask.any():
        return mask
    norm_barcodes = {_normalize_barcode(value) for value in barcodes}
    norm_obs = obs_names.apply(_normalize_barcode)
    return norm_obs.isin(norm_barcodes).values


def _compute_groupby_mask(adata, groupby, group_values, groupby_rev=True, groupby_sample=None):
    if not groupby or not group_values:
        return None
    t_start = time.time()
    if groupby.endswith('.tsv') or groupby.endswith('.txt'):
        mapping = _load_groupby_file(groupby, reverse_complement=groupby_rev, sample_name=groupby_sample)
        obs_keys = adata.obs_names.astype(str).to_series().apply(_normalize_barcode)
        groups = obs_keys.map(mapping)
        mask = groups.isin(group_values).values
        print(f"[timing] groupby load+map: {time.time() - t_start:.2f}s")
        return mask
    if groupby in adata.obs:
        mask = adata.obs[groupby].isin(group_values).values
        print(f"[timing] groupby mask: {time.time() - t_start:.2f}s")
        return mask
    print(f"[warn] groupby field not found: {groupby}")
    return None


def load_isoform_counts(input_path, target_gene, groupby=None, group_values=None, groupby_rev=True, groupby_sample=None):
    adata = _load_adata(input_path)
    groupby_field = groupby
    if groupby and (groupby.endswith('.tsv') or groupby.endswith('.txt')):
        mapping = _load_groupby_file(groupby, reverse_complement=groupby_rev, sample_name=groupby_sample)
        obs_keys = adata.obs_names.astype(str).to_series().apply(_normalize_barcode)
        adata.obs['groupby_file'] = obs_keys.map(mapping)
        groupby_field = 'groupby_file'
    if groupby_field and groupby_field in adata.obs:
        if group_values:
            mask = adata.obs[groupby_field].isin(group_values)
            adata = adata[mask, :]
    counts = np.asarray(adata.X.sum(axis=0)).ravel()
    var_names = np.asarray(adata.var_names)
    gene_col = None
    if 'gene' in adata.var.columns:
        gene_col = np.asarray(adata.var['gene'])

    records = []
    for idx, var in enumerate(var_names):
        gene_hint = gene_col[idx] if gene_col is not None else None
        gene_id, isoform_id = parse_feature_identifier(str(var), gene_hint=gene_hint)
        if gene_id is not None and gene_id != target_gene:
            continue
        records.append({
            'var': var,
            'gene_id': gene_id,
            'isoform_id': isoform_id,
            'count': counts[idx]
        })
    return records


def load_isoform_counts_indexed(h5ad_path, target_gene, index_dir=None,
                                groupby=None, group_values=None,
                                groupby_rev=True, groupby_sample=None,
                                adata=None, gene_index=None, mask=None,
                                use_counts_cache=True, cache_build=True):
    t_start = time.time()
    if gene_index is None:
        t_index = time.time()
        gene_to_pos, indptr, indices = _load_h5ad_gene_index(h5ad_path, index_dir=index_dir)
        print(f"[timing] h5ad index load: {time.time() - t_index:.2f}s")
    else:
        gene_to_pos, indptr, indices = gene_index
    pos = gene_to_pos.get(target_gene)
    if pos is None:
        return []
    start = int(indptr[pos])
    end = int(indptr[pos + 1])
    if end <= start:
        return []
    gene_idx = np.asarray(indices[start:end], dtype=np.int64)
    if gene_idx.ndim == 0:
        gene_idx = gene_idx.reshape(1)

    cache = _load_or_build_counts_cache(
        h5ad_path,
        adata,
        mask,
        groupby,
        group_values,
        groupby_rev,
        groupby_sample,
        use_cache=use_counts_cache,
        allow_build=cache_build,
    )
    if cache is not None:
        t_cache = time.time()
        records = _records_from_counts_cache(cache, gene_idx, target_gene)
        print(f"[timing] isoform counts cache lookup: {time.time() - t_cache:.2f}s")
        print(f"[timing] total load_isoform_counts_indexed: {time.time() - t_start:.2f}s")
        return records

    close_adata = False
    if adata is None:
        t_open = time.time()
        adata = ad.read_h5ad(h5ad_path, backed='r')
        print(f"[timing] h5ad open: {time.time() - t_open:.2f}s")
        close_adata = True
    try:
        if mask is None:
            mask = _compute_groupby_mask(
                adata,
                groupby,
                group_values,
                groupby_rev=groupby_rev,
                groupby_sample=groupby_sample,
            )

        t_slice = time.time()
        if mask is not None:
            if isinstance(mask, tuple):
                print(f"[debug] groupby mask tuple: types={[type(item) for item in mask]}")
                raise ValueError("groupby mask is a tuple; expected 1D boolean array.")
            mask = np.asarray(mask, dtype=bool).ravel()
            if mask.ndim != 1:
                raise ValueError(f"groupby mask has unexpected shape: {mask.shape}")
            if gene_idx.ndim != 1:
                raise ValueError(f"gene_idx has unexpected shape: {gene_idx.shape}")
            print(f"[debug] groupby mask shape: {mask.shape} gene_idx shape: {gene_idx.shape}")
            row_idx = np.flatnonzero(mask)
            if row_idx.size == 0:
                return []
        else:
            row_idx = None

        sub_X = adata[:, gene_idx].X
        if row_idx is not None:
            sub_X = sub_X[row_idx, :]
        counts = np.asarray(sub_X.sum(axis=0)).ravel()
        print(f"[timing] h5ad slice+sum: {time.time() - t_slice:.2f}s")
        var_names = np.asarray(adata.var_names)[gene_idx]
        gene_col = None
        if 'gene' in adata.var.columns:
            gene_col = np.asarray(adata.var['gene'])[gene_idx]

        t_records = time.time()
        records = []
        for idx, var in enumerate(var_names):
            gene_hint = gene_col[idx] if gene_col is not None else None
            gene_id, isoform_id = parse_feature_identifier(str(var), gene_hint=gene_hint)
            if gene_id is not None and gene_id != target_gene:
                continue
            records.append({
                'var': var,
                'gene_id': gene_id,
                'isoform_id': isoform_id,
                'count': counts[idx]
            })
        print(f"[timing] record build: {time.time() - t_records:.2f}s")
        print(f"[timing] total load_isoform_counts_indexed: {time.time() - t_start:.2f}s")
        return records
    finally:
        if close_adata:
            try:
                adata.file.close()
            except AttributeError:
                pass

def plot_isoform_structures_by_conditions(sample_dict, conditions, cell_states, barcode_sample_dict,
                                          genes, gene_model, output_dir, index_dir=None,
                                          cluster_mode='block', cluster_features='tokens',
                                          cluster_strategy='substring',
                                          cluster_similarity_threshold=0.85, min_split_fraction=0.05,
                                          max_isoforms=200, min_count=1, intron_scale=0.2,
                                          row_height=0.0125, row_gap=0.0, group_gap=0.1,
                                          isoform_gap=None, label_mode='first',
                                          cluster_isoforms=True):
    start_time = time.time()
    if isinstance(genes, str):
        genes = [genes]
    if isinstance(conditions, str):
        conditions = [value.strip() for value in conditions.split(',') if value.strip()]
    if cell_states and isinstance(cell_states, str):
        cell_states = [state.strip() for state in cell_states.split(',') if state.strip()]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"[auto] Conditions: {conditions}")
    print(f"[auto] Cell states: {cell_states or ['all']}")
    print(f"[auto] Genes: {genes}")
    print(f"[auto] Output dir: {output_dir}")

    condition_counts = {
        condition: {gene: defaultdict(float) for gene in genes} for condition in conditions
    }
    transcript_union = {gene: {} for gene in genes}

    def infer_sample_name(sample_entry, h5ad_path):
        for key in ('library', 'gff_name'):
            value = sample_entry.get(key)
            if value:
                return str(value)
        return _infer_sample_name(h5ad_path)

    for uid, samples in sample_dict.items():
        for sample in samples:
            sample_conditions = [
                condition for condition in conditions
                if _matches_group(sample.get('groups'), condition)
            ]
            if not sample_conditions:
                continue
            h5ad_path = sample.get('matrix')
            if not h5ad_path:
                continue
            sample_name = infer_sample_name(sample, h5ad_path)
            print(f"[auto] Sample: {sample_name} | Conditions: {sample_conditions}")
            barcode_series = None
            if barcode_sample_dict and sample_name in barcode_sample_dict:
                barcode_series = barcode_sample_dict[sample_name]
            gene_to_pos, indptr, indices = _load_h5ad_gene_index(
                h5ad_path, index_dir=index_dir
            )
            cache_group_values = cell_states or ['all']
            counts_cache = _load_or_build_counts_cache(
                h5ad_path,
                None,
                None,
                "cell_states",
                cache_group_values,
                False,
                sample_name,
                use_cache=True,
                allow_build=False,
            )
            adata = None
            mask = None
            row_idx = None
            if counts_cache is None:
                adata = ad.read_h5ad(h5ad_path, backed='r')
                mask = _load_barcode_mask(adata, barcode_series, cell_states)
                if mask is not None and not mask.all():
                    row_idx = np.flatnonzero(mask)
                counts_cache = _load_or_build_counts_cache(
                    h5ad_path,
                    adata,
                    mask,
                    "cell_states",
                    cache_group_values,
                    False,
                    sample_name,
                    use_cache=True,
                    allow_build=True,
                )
            assoc_path = sample.get('transcript_associations')
            if not assoc_path:
                assoc_path = os.path.join(
                    os.path.dirname(h5ad_path),
                    'gff-output',
                    'transcript_associations.txt'
                )
                if not os.path.exists(assoc_path):
                    assoc_path = None
            for gene in genes:
                pos = gene_to_pos.get(gene)
                if pos is None:
                    continue
                start = int(indptr[pos])
                end = int(indptr[pos + 1])
                if end <= start:
                    continue
                raw_idx = None
                try:
                    raw_idx = indices[start:end]
                    if isinstance(raw_idx, tuple):
                        if len(raw_idx) == 1:
                            gene_idx = raw_idx[0]
                        else:
                            gene_idx = np.concatenate(
                                [np.asarray(part, dtype=np.int64).ravel() for part in raw_idx]
                            )
                    else:
                        gene_idx = np.asarray(raw_idx)
                    if gene_idx.ndim != 1 or gene_idx.dtype == object:
                        if 'gene' in adata.var.columns:
                            gene_idx = np.where(np.asarray(adata.var['gene']) == gene)[0]
                        else:
                            gene_idx = np.array([], dtype=np.int64)
                    gene_idx = np.asarray(gene_idx, dtype=np.int64).ravel()
                except Exception as exc:
                    print(
                        "[auto][error] gene_idx build failed "
                        f"gene={gene} sample={sample_name} h5ad={h5ad_path}: {exc}"
                    )
                    print(f"[auto][error] raw_idx type={type(raw_idx)} value={str(raw_idx)[:200]}")
                    raise
                if gene_idx.size == 0:
                    continue
                if counts_cache is not None:
                    try:
                        records = _records_from_counts_cache(counts_cache, gene_idx, gene)
                        isoform_counts = {rec['var']: rec['count'] for rec in records}
                    except Exception as exc:
                        print(
                            "[auto][error] counts cache lookup failed "
                            f"gene={gene} sample={sample_name}: {exc}"
                        )
                        raise
                else:
                    gene_idx = np.asarray(gene_idx, dtype=np.int64).ravel()
                    gene_idx_list = [int(value) for value in gene_idx]
                    try:
                        sub_X = adata[:, gene_idx_list].X
                    except Exception as exc:
                        print(
                            "[auto][error] h5ad slice failed "
                            f"gene={gene} sample={sample_name} idx_shape={getattr(gene_idx, 'shape', None)} "
                            f"idx_dtype={getattr(gene_idx, 'dtype', None)} adata_X={type(adata.X)}: {exc}"
                        )
                        print(f"[auto][error] gene_idx_list sample={gene_idx_list[:10]}")
                        raise
                    if row_idx is not None:
                        if row_idx.size == 0:
                            continue
                        sub_X = sub_X[row_idx, :]
                    counts = np.asarray(sub_X.sum(axis=0)).ravel()
                    isoforms = [str(value) for value in np.asarray(adata.var_names)[gene_idx]]
                    isoform_counts = dict(zip(isoforms, counts))
                for condition in sample_conditions:
                    store = condition_counts[condition][gene]
                    for iso_id, count in isoform_counts.items():
                        store[iso_id] += count
                if assoc_path:
                    structures, _strand_map = load_transcript_associations_auto(
                        assoc_path, gene, index_dir=None
                    )
                    if structures:
                        transcript_union[gene].update(structures)
            if adata is not None:
                try:
                    adata.file.close()
                except AttributeError:
                    pass

    gene_model_db = _ensure_gene_model_db(gene_model, index_dir=index_dir) if gene_model else None
    if gene_model_db and os.path.exists(gene_model_db):
        gene_segments = defaultdict(list)
        exon_lookup = {}
        for gene in genes:
            segments, lookup = read_gene_model_from_sqlite(gene_model_db, gene)
            for key, value in segments.items():
                gene_segments[key].extend(value)
            exon_lookup.update(lookup)
    else:
        gene_segments, exon_lookup = read_gene_model(gene_model)
    gene_maps = build_gene_maps(gene_segments, intron_scale)

    cell_state_label = 'all'
    if cell_states:
        cell_state_label = '+'.join(cell_states)

    for gene in genes:
        transcript_structures = transcript_union.get(gene)
        if not transcript_structures:
            print(f"[auto] No transcript structures for {gene}, skipping.")
            continue

        combined_counts = defaultdict(float)
        for condition in conditions:
            for iso_id, count in condition_counts[condition][gene].items():
                combined_counts[iso_id] += count
        print(f"[auto] Aggregated isoforms for {gene}: {len(combined_counts)}")

        combined_records = [
            {'isoform_id': iso_id, 'count': count, 'gene_id': gene}
            for iso_id, count in combined_counts.items()
        ]
        if not combined_records:
            print(f"[auto] No isoform counts for {gene}, skipping.")
            continue

        grouped_structures = None
        isoform_colors = None
        plotted, _, _ = build_plotted_isoforms(
            combined_records, transcript_structures, gene,
            cluster_mode, cluster_features, min_count, max_isoforms
        )
        if plotted:
            grouped_structures = group_structures(
                plotted, gene, cluster_features,
                cluster_strategy, cluster_similarity_threshold,
                min_split_fraction, include_introns=cluster_isoforms,
                report=cluster_isoforms
            )
            if cluster_isoforms:
                isoform_colors = assign_cluster_colors(grouped_structures, plt.get_cmap(DEFAULT_COLORMAP).colors)
            else:
                isoform_colors = assign_isoform_colors(grouped_structures, plt.get_cmap(DEFAULT_COLORMAP).colors)

        for condition in conditions:
            counts_map = condition_counts[condition][gene]
            if not counts_map:
                print(f"[auto] No counts for {gene} in {condition}, skipping.")
                continue
            isoform_records = [
                {'isoform_id': iso_id, 'count': count, 'gene_id': gene}
                for iso_id, count in counts_map.items()
            ]
            condition_dir = output_dir / condition / gene
            condition_dir.mkdir(parents=True, exist_ok=True)
            out_path = condition_dir / f"{condition}__{cell_state_label}__{gene}.pdf"
            print(f"[auto] Plotting {gene} for {condition} -> {out_path}")
            plotted, plot_cluster_label_map, plot_cluster_index_map = plot_isoform_structures(
                isoform_records,
                transcript_structures,
                exon_lookup,
                gene_maps,
                cluster_mode,
                gene,
                intron_scale,
                str(out_path),
                max_isoforms=max_isoforms,
                min_count=min_count,
                row_height=row_height,
                row_gap=row_gap,
                group_gap=group_gap,
                isoform_gap=isoform_gap,
                cluster_strategy=cluster_strategy,
                cluster_similarity_threshold=cluster_similarity_threshold,
                cluster_features=cluster_features,
                label_mode=label_mode,
                min_split_fraction=min_split_fraction,
                grouped_structures=grouped_structures,
                isoform_colors=isoform_colors,
                isoform_counts=counts_map,
                cluster_isoforms=cluster_isoforms,
                return_cluster_labels=True
            )
            export_plotted_isoforms(plotted, str(out_path), plot_cluster_label_map, plot_cluster_index_map)
    elapsed = time.time() - start_time
    print(f"[auto] Elapsed time (s): {elapsed:.2f}")

def filter_records_by_structures(records, transcript_structures):
    filtered = []
    for rec in records:
        if resolve_isoform_id(rec['isoform_id'], transcript_structures):
            filtered.append(rec)
    return filtered


def coerce_read_count(value):
    try:
        count = float(value)
    except (TypeError, ValueError):
        return 0
    if count <= 0:
        return 0
    rounded = int(round(count))
    if abs(count - rounded) < 1e-6:
        return rounded
    return max(1, rounded)


def build_feature_set(tokens, feature_mode):
    if feature_mode == 'junctions':
        return {f"{tokens[i]}->{tokens[i + 1]}" for i in range(len(tokens) - 1)}
    if feature_mode == 'both':
        features = {f"{tokens[i]}->{tokens[i + 1]}" for i in range(len(tokens) - 1)}
        features.update(tokens)
        return features
    return set(tokens)


def strip_terminal_coords(tokens, default_gene):
    if not tokens:
        return []
    trimmed = list(tokens)
    drop_indices = set()
    for idx in (0, len(trimmed) - 1):
        gene_id, base, coord = split_token(trimmed[idx], default_gene)
        if coord is not None and base.startswith('E'):
            drop_indices.add(idx)
    if drop_indices:
        return [tok for i, tok in enumerate(trimmed) if i not in drop_indices]
    return trimmed


def filter_exon_tokens(tokens, default_gene):
    exon_tokens = []
    for tok in tokens:
        gene_id, base, _ = split_token(tok, default_gene)
        if base.startswith('E'):
            if '_' in tok:
                continue
            if gene_id != default_gene and ':' not in tok:
                tok = f"{gene_id}:{base}"
            exon_tokens.append(tok)
    return exon_tokens


def filter_exon_intron_tokens(tokens, default_gene):
    kept_tokens = []
    for tok in tokens:
        gene_id, base, _ = split_token(tok, default_gene)
        if base.startswith(('E', 'I')):
            if base.startswith('E') and '_' in tok:
                continue
            if gene_id != default_gene and ':' not in tok:
                tok = f"{gene_id}:{base}"
            kept_tokens.append(tok)
    return kept_tokens


def cluster_by_feature(items, min_split_fraction=0.05):
    def split(group):
        if len(group) <= 1:
            return [group]
        total_weight = sum(item['weight'] for item in group)
        if total_weight <= 0:
            return [group]

        feature_weights = defaultdict(float)
        for item in group:
            for feat in item['feature_set']:
                feature_weights[feat] += item['weight']

        best_feat = None
        best_score = 0
        for feat, weight in feature_weights.items():
            if weight <= 0 or weight >= total_weight:
                continue
            score = min(weight, total_weight - weight)
            if score > best_score:
                best_score = score
                best_feat = feat

        if best_feat is None:
            return [group]

        split_fraction = best_score / total_weight
        if split_fraction < min_split_fraction:
            return [group]

        present = [item for item in group if best_feat in item['feature_set']]
        absent = [item for item in group if best_feat not in item['feature_set']]
        if not present or not absent:
            return [group]

        weight_present = sum(item['weight'] for item in present)
        weight_absent = sum(item['weight'] for item in absent)
        if weight_present >= weight_absent:
            ordered = [present, absent]
        else:
            ordered = [absent, present]
        return split(ordered[0]) + split(ordered[1])

    return split(items)


def cluster_by_jaccard(items, threshold=0.4):
    clusters = []

    def jaccard(a_set, b_set):
        union = a_set | b_set
        inter = a_set & b_set
        return (len(inter) / len(union)) if union else 0.0

    for item in sorted(items, key=lambda x: x['weight'], reverse=True):
        best_idx = None
        best_score = -1.0
        for idx, cluster in enumerate(clusters):
            score = jaccard(item['cluster_set'], cluster['rep_set'])
            if score > best_score:
                best_score = score
                best_idx = idx
        if best_idx is not None and best_score >= threshold:
            clusters[best_idx]['items'].append(item)
        else:
            clusters.append({'items': [item], 'rep_set': item['cluster_set']})
    return [cluster['items'] for cluster in clusters]


def longest_common_substring_length(seq_a, seq_b):
    if not seq_a or not seq_b:
        return 0
    if len(seq_a) > len(seq_b):
        seq_a, seq_b = seq_b, seq_a
    prev = [0] * (len(seq_a) + 1)
    best = 0
    for token in seq_b:
        current = [0]
        for idx, a_token in enumerate(seq_a, start=1):
            if token == a_token:
                val = prev[idx - 1] + 1
            else:
                val = 0
            current.append(val)
            if val > best:
                best = val
        prev = current
    return best


def substring_similarity(seq_a, seq_b):
    if not seq_a or not seq_b:
        return 0.0
    short = seq_a if len(seq_a) <= len(seq_b) else seq_b
    long = seq_b if short is seq_a else seq_a
    lcs = longest_common_substring_length(short, long)
    return lcs / float(len(short))


def cluster_by_substring(items, threshold=0.4):
    clusters = []

    for item in sorted(items, key=lambda x: x['weight'], reverse=True):
        best_idx = None
        best_score = -1.0
        for idx, cluster in enumerate(clusters):
            min_score = 1.0
            for other in cluster['items']:
                score = substring_similarity(item['cluster_seq'], other['cluster_seq'])
                if score < min_score:
                    min_score = score
                if min_score < threshold:
                    break
            if min_score >= threshold and min_score > best_score:
                best_score = min_score
                best_idx = idx
        if best_idx is not None:
            clusters[best_idx]['items'].append(item)
        else:
            clusters.append({'items': [item]})
    return [cluster['items'] for cluster in clusters]


def is_subsequence(short_seq, long_seq):
    if not short_seq:
        return True
    long_iter = iter(long_seq)
    for item in short_seq:
        for candidate in long_iter:
            if candidate == item:
                break
        else:
            return False
    return True


def build_plotted_isoforms(isoforms, transcript_structures, target_gene, cluster_mode,
                           cluster_features, min_count, max_isoforms):
    filtered = [i for i in isoforms if i['count'] >= min_count]
    missing_ids = []
    resolved_items = []
    for iso in filtered:
        resolved = resolve_isoform_id(iso['isoform_id'], transcript_structures)
        if not resolved:
            missing_ids.append(iso['isoform_id'])
            continue
        resolved_items.append((iso, resolved))

    resolved_items.sort(key=lambda x: x[0]['count'], reverse=True)
    resolved_items = resolved_items[:max_isoforms]

    plotted = []
    for iso, resolved in resolved_items:
        tokens = transcript_structures[resolved]
        tokens_trimmed = strip_terminal_coords(tokens, target_gene)
        cluster_tokens = filter_exon_tokens(tokens_trimmed, target_gene)
        cluster_set = set(cluster_tokens)
        tokens_norm = [normalize_token(tok, target_gene, cluster_mode) for tok in tokens]
        plotted.append({
            **iso,
            'resolved_id': resolved,
            'tokens': tokens,
            'tokens_trimmed': tokens_trimmed,
            'cluster_tokens': cluster_tokens,
            'cluster_set': cluster_set,
            'tokens_norm': tokens_norm,
            'feature_set': build_feature_set(tokens_trimmed, cluster_features)
        })
    return plotted, len(missing_ids), missing_ids


def build_plotted_rows(plotted, cluster_label_map=None, cluster_index_map=None):
    rows = []
    seen = set()
    for iso in plotted:
        isoform_id = iso.get('isoform_id') or iso.get('resolved_id')
        resolved_id = iso.get('resolved_id') or isoform_id
        key = (isoform_id, resolved_id)
        if key in seen:
            continue
        seen.add(key)
        row = {
            'isoform_id': isoform_id,
            'resolved_id': resolved_id,
            'count': iso.get('count', 0),
            'cluster_label': (cluster_label_map or {}).get(resolved_id)
                            or (cluster_label_map or {}).get(isoform_id)
        }
        if cluster_index_map:
            row['cluster_index'] = cluster_index_map.get(resolved_id)
        rows.append(row)
    return rows


def export_plotted_isoforms(plotted, output_path, cluster_label_map=None, cluster_index_map=None):
    if not plotted:
        return None, []
    rows = build_plotted_rows(plotted, cluster_label_map, cluster_index_map)
    output_base = os.path.splitext(output_path)[0]
    isoform_path = f"{output_base}_isoform_ids.tsv"
    pd.DataFrame(rows).to_csv(isoform_path, sep='\t', index=False)
    return isoform_path, rows


def sample_isoform_ids(records, limit=5):
    samples = []
    for rec in records:
        isoform_id = rec.get('isoform_id')
        if isoform_id:
            samples.append(isoform_id)
        if len(samples) >= limit:
            break
    return samples


def build_cluster_label_map(records, transcript_structures, target_gene, cluster_mode,
                            cluster_features, cluster_strategy,
                            cluster_similarity_threshold, min_split_fraction,
                            max_isoforms, min_count, cluster_isoforms,
                            inspect_isoforms=None):
    if cluster_strategy is None:
        cluster_strategy = 'substring' if cluster_isoforms else 'jaccard'
    if cluster_strategy in ('subsequence', 'overlap'):
        cluster_strategy = 'substring'
    if cluster_similarity_threshold is None:
        cluster_similarity_threshold = 0.85 if cluster_isoforms else 0.4
    plotted, _, _ = build_plotted_isoforms(
        records,
        transcript_structures,
        target_gene,
        cluster_mode,
        cluster_features,
        min_count,
        max_isoforms
    )
    plotted_ids = set()
    for iso in plotted:
        iso_key = iso.get('isoform_id', iso.get('resolved_id'))
        if iso_key:
            plotted_ids.add(iso_key)
        if iso.get('resolved_id'):
            plotted_ids.add(iso.get('resolved_id'))
    if not plotted:
        return {}, plotted_ids
    if inspect_isoforms:
        max_isoforms = max(max_isoforms, len(records))
        plotted, _, _ = build_plotted_isoforms(
            records,
            transcript_structures,
            target_gene,
            cluster_mode,
            cluster_features,
            min_count,
            max_isoforms
        )
    grouped_structures = group_structures(
        plotted,
        target_gene,
        cluster_features,
        cluster_strategy,
        cluster_similarity_threshold,
        min_split_fraction,
        include_introns=cluster_isoforms,
        report=False
    )
    if cluster_isoforms:
        for group in grouped_structures:
            group.sort(key=structure_max_length, reverse=True)
            for structure in group:
                structure['items'].sort(
                    key=lambda iso: (isoform_length(iso), iso.get('count', 0)),
                    reverse=True
                )
        grouped_structures.sort(
            key=lambda group: (group_max_length(group), sum(item['weight'] for item in group)),
            reverse=True
        )
    label_map = {}
    for group in grouped_structures:
        label = select_cluster_label(group) if cluster_isoforms else None
        for structure in group:
            for iso in structure['items']:
                iso_key = iso.get('isoform_id', iso.get('resolved_id'))
                if cluster_isoforms:
                    label_map[iso_key] = label
                    label_map[iso.get('resolved_id')] = label
                else:
                    label_map[iso_key] = iso.get('resolved_id')
    return label_map, plotted_ids


def report_isoform_structures(records, transcript_structures, target_gene, isoform_ids,
                              include_introns=False, cluster_label_map=None,
                              plotted_ids=None, plot_cluster_label_map=None,
                              plot_cluster_index_map=None):
    if not isoform_ids:
        return
    record_lookup = {}
    for rec in records:
        isoform_id = rec.get('isoform_id')
        var = rec.get('var')
        if isoform_id:
            record_lookup.setdefault(isoform_id, rec)
        if var:
            record_lookup.setdefault(var, rec)
    resolved_ids = []
    cluster_seqs = []
    print("Isoform structure inspection:")
    for raw_id in isoform_ids:
        resolved, matches = resolve_isoform_id_fuzzy(raw_id, transcript_structures)
        if not resolved:
            if matches:
                preview = ", ".join(matches[:3])
                print(f"  {raw_id}: multiple matches ({len(matches)}) {preview}")
            else:
                print(f"  {raw_id}: no match in transcript_associations")
            continue
        tokens = transcript_structures[resolved]
        trimmed = strip_terminal_coords(list(tokens), target_gene)
        if include_introns:
            cluster_seq = filter_exon_intron_tokens(list(trimmed), target_gene)
        else:
            cluster_seq = filter_exon_tokens(list(trimmed), target_gene)
        resolved_ids.append(resolved)
        cluster_seqs.append(cluster_seq)
        rec = record_lookup.get(raw_id) or record_lookup.get(resolved)
        label = resolved if raw_id == resolved else f"{raw_id} -> {resolved}"
        print(f"  {label}")
        if cluster_label_map:
            display_label = cluster_label_map.get(resolved) or cluster_label_map.get(raw_id)
            if display_label:
                print(f"    cluster_label={display_label}")
        if plot_cluster_label_map:
            plot_label = plot_cluster_label_map.get(resolved) or plot_cluster_label_map.get(raw_id)
            if plot_label:
                print(f"    plot_cluster_label={plot_label}")
        if plot_cluster_index_map:
            plot_index = plot_cluster_index_map.get(resolved) or plot_cluster_index_map.get(raw_id)
            if plot_index is not None:
                print(f"    plot_cluster_index={plot_index}")
        if plotted_ids is not None:
            plotted_flag = (resolved in plotted_ids) or (raw_id in plotted_ids)
            print(f"    plotted={str(plotted_flag).lower()}")
        if rec and rec.get('count') is not None:
            print(f"    count={rec.get('count')}")
        if rec and rec.get('var') and rec.get('var') not in (raw_id, resolved):
            print(f"    var={rec.get('var')}")
        print(f"    tokens={'|'.join(tokens)}")
        print(f"    trimmed={'|'.join(trimmed)}")
        print(f"    cluster_seq={'|'.join(cluster_seq)}")
    if len(resolved_ids) == 2:
        similarity = substring_similarity(cluster_seqs[0], cluster_seqs[1])
        print(f"  substring_similarity={similarity:.3f}")


def group_structures(plotted, target_gene, cluster_features, cluster_strategy,
                     cluster_similarity_threshold, min_split_fraction,
                     include_introns=False, report=False):
    structure_groups = OrderedDict()
    for iso in plotted:
        key = tuple(iso['tokens'])
        trimmed_key = tuple(strip_terminal_coords(list(key), target_gene))
        if include_introns:
            cluster_tokens = filter_exon_intron_tokens(list(trimmed_key), target_gene)
        else:
            cluster_tokens = filter_exon_tokens(list(trimmed_key), target_gene)
        cluster_set = set(cluster_tokens)
        if key not in structure_groups:
            structure_groups[key] = {
                'structure_key': key,
                'structure_key_trimmed': trimmed_key,
                'items': [],
                'weight': 0,
                'feature_set': build_feature_set(list(trimmed_key), cluster_features),
                'cluster_set': cluster_set,
                'cluster_seq': cluster_tokens
            }
        structure_groups[key]['items'].append(iso)
        structure_groups[key]['weight'] += iso['count']

    structure_items = list(structure_groups.values())
    if cluster_strategy == 'jaccard':
        grouped_structures = cluster_by_jaccard(
            structure_items,
            threshold=cluster_similarity_threshold
        )
    elif cluster_strategy in ('substring', 'subsequence', 'overlap'):
        grouped_structures = cluster_by_substring(
            structure_items,
            threshold=cluster_similarity_threshold
        )
    elif cluster_strategy == 'feature':
        grouped_structures = cluster_by_feature(
            structure_items,
            min_split_fraction=min_split_fraction
        )
    elif cluster_strategy == 'structure':
        grouped_structures = [[item] for item in structure_items]
    else:
        grouped_structures = [[item] for item in sorted(structure_items, key=lambda x: x['weight'], reverse=True)]

    if report:
        total_structures = len(structure_items)
        cluster_count = len(grouped_structures)
        merged_count = total_structures - cluster_count
        cluster_sizes = [len(group) for group in grouped_structures]
        mean_size = (sum(cluster_sizes) / len(cluster_sizes)) if cluster_sizes else 0
        max_size = max(cluster_sizes) if cluster_sizes else 0
        summary_parts = [
            "Cluster summary",
            f"strategy={cluster_strategy}",
            f"structures={total_structures}",
            f"clusters={cluster_count}",
            f"merged={merged_count}",
            f"mean_size={mean_size:.2f}",
            f"max_size={max_size}"
        ]
        if cluster_strategy in ('substring', 'subsequence', 'overlap', 'jaccard'):
            summary_parts.append(f"threshold={cluster_similarity_threshold}")
        else:
            summary_parts.append(f"min_split_fraction={min_split_fraction}")
        print(" ".join(summary_parts))

    for group in grouped_structures:
        rep = max(group, key=lambda x: x['weight'])
        rep_set = rep['cluster_set']

        def similarity_key(item):
            subset = item['cluster_set'].issubset(rep_set)
            union = rep_set | item['cluster_set']
            inter = rep_set & item['cluster_set']
            jaccard = len(inter) / len(union) if union else 0.0
            overlap = (len(inter) / min(len(rep_set), len(item['cluster_set']))) if rep_set and item['cluster_set'] else 0.0
            rep_rank = 0 if item is rep else 1
            return (rep_rank, -int(subset), -overlap, -jaccard, -item['weight'])

        group.sort(key=similarity_key)
        for entry in group:
            entry['items'].sort(key=lambda x: x['count'], reverse=True)

    return grouped_structures


def assign_isoform_colors(grouped_structures, palette):
    isoform_colors = {}
    color_idx = 0
    for group in grouped_structures:
        for structure in group:
            for iso in structure['items']:
                iso_key = iso.get('isoform_id', iso['resolved_id'])
                if iso_key not in isoform_colors:
                    isoform_colors[iso_key] = palette[color_idx % len(palette)]
                    color_idx += 1
    return isoform_colors


def assign_cluster_colors(grouped_structures, palette):
    isoform_colors = {}
    for group_idx, group in enumerate(grouped_structures):
        color = palette[group_idx % len(palette)]
        for structure in group:
            for iso in structure['items']:
                iso_key = iso.get('isoform_id', iso['resolved_id'])
                if iso_key not in isoform_colors:
                    isoform_colors[iso_key] = color
    return isoform_colors


def select_cluster_label(group):
    best_label = None
    best_score = None
    for structure in group:
        for iso in structure['items']:
            label = iso.get('resolved_id') or iso.get('isoform_id')
            if not label:
                continue
            tokens = iso.get('tokens_trimmed') or iso.get('tokens') or []
            score = (len(tokens), iso.get('count', 0), len(str(label)))
            if best_score is None or score > best_score:
                best_score = score
                best_label = label
    return best_label


def isoform_length(iso):
    tokens = iso.get('tokens_trimmed') or iso.get('tokens') or []
    return len(tokens)


def structure_max_length(structure):
    return max((isoform_length(iso) for iso in structure['items']), default=0)


def group_max_length(group):
    return max((structure_max_length(structure) for structure in group), default=0)


# Parameters:
# - isoforms: list of dicts with isoform IDs and read counts.
# - transcript_structures: transcript_id -> list of exon/intron tokens.
# - exon_lookup: (gene, exon_id) -> coordinate metadata.
# - gene_maps: gene -> display coordinate map (scaled introns).
# - cluster_mode: token normalization mode (full/base/block).
# - target_gene: gene ID used for filtering and token defaults.
# - intron_scale: display scale factor for intron lengths.
# - output_path: file path for the saved plot.
# - max_isoforms: cap on number of isoforms to render.
# - min_count: minimum read count to keep an isoform.
# - row_height: height of each read row in plot units.
# - row_gap: vertical gap between read rows.
# - group_gap: vertical gap between clustered isoform groups.
# - isoform_gap: vertical gap between isoform stacks (defaults to group_gap).
# - gene_gap: horizontal gap between genes in trans-splicing plots.
# - cluster_strategy: ordering method for isoform structures.
# - cluster_features: features used for clustering (tokens/junctions/both).
# - label_mode: label display mode (first/none/all).
# - min_split_fraction: minimum weighted split for feature clustering.
# - cluster_similarity_threshold: minimum similarity threshold for clustering.
# - grouped_structures: precomputed cluster order to reuse across files.
# - isoform_colors: precomputed isoform-to-color mapping.
# - isoform_counts: per-file isoform counts for consistent ordering.
# - cluster_isoforms: color isoforms by cluster and label with a representative isoform ID.
# - debug_transcripts: transcript IDs to print segment details for.
def plot_isoform_structures(isoforms, transcript_structures, exon_lookup, gene_maps,
                            cluster_mode, target_gene, intron_scale, output_path,
                            max_isoforms=50, min_count=1, row_height=0.0125,
                            row_gap=0.0, group_gap=0.0, isoform_gap=None, gene_gap=1000,
                            cluster_strategy='jaccard', cluster_features='tokens',
                            label_mode='first', min_split_fraction=0.05,
                            cluster_similarity_threshold=0.4,
                            debug_isoforms=None,
                            debug_transcripts=None,
                            grouped_structures=None,
                            isoform_colors=None,
                            isoform_counts=None,
                            transcript_associations_path=None,
                            cluster_isoforms=False,
                            return_cluster_labels=False):
    """Render stacked isoform read tracks for a single gene.

    isoforms: list of dicts with isoform IDs and read counts.
    transcript_structures: transcript_id -> list of exon/intron tokens.
    exon_lookup: (gene, exon_id) -> coordinate metadata.
    gene_maps: gene -> display coordinate map (scaled introns).
    cluster_mode: token normalization mode (full/base/block).
    target_gene: gene ID used for filtering and token defaults.
    intron_scale: display scale factor for intron lengths.
    output_path: file path for the saved plot.
    max_isoforms: cap on number of isoforms to render.
    min_count: minimum read count to keep an isoform.
    row_height: height of each read row in plot units.
    row_gap: vertical gap between read rows.
    group_gap: vertical gap between clustered isoform groups.
    isoform_gap: vertical gap between isoform stacks (defaults to group_gap).
    gene_gap: horizontal gap between genes in trans-splicing plots.
    cluster_strategy: ordering method for isoform structures.
    cluster_features: features used for clustering (tokens/junctions/both).
    label_mode: label display mode (first/none/all).
    min_split_fraction: minimum weighted split for feature clustering.
    cluster_similarity_threshold: minimum similarity threshold for clustering.
    grouped_structures: precomputed cluster order to reuse across files.
    isoform_colors: precomputed isoform-to-color mapping.
    isoform_counts: per-file isoform counts for consistent ordering.
    debug_transcripts: transcript IDs to print segment details for.
    """
    if cluster_strategy is None:
        cluster_strategy = 'substring' if cluster_isoforms else 'jaccard'
    if cluster_strategy in ('subsequence', 'overlap'):
        cluster_strategy = 'substring'
    if cluster_similarity_threshold is None:
        cluster_similarity_threshold = 0.85 if cluster_isoforms else 0.4
    if grouped_structures is None:
        t_build = time.time()
        plotted, missing_structures, missing_ids = build_plotted_isoforms(
            isoforms,
            transcript_structures,
            target_gene,
            cluster_mode,
            cluster_features,
            min_count,
            max_isoforms
        )
        print(f"[timing] build_plotted_isoforms: {time.time() - t_build:.2f}s")
        if not plotted:
            print("No isoforms with transcript structures available to plot.")
            missing_sample = sorted(set(missing_ids))[:5]
            if missing_sample:
                print("Sample isoform IDs missing from transcript_associations:", missing_sample)
            if transcript_associations_path:
                transcript_ids, alias_to_transcript = load_transcript_association_index(
                    transcript_associations_path, target_gene
                )
                matrix_ids = {rec['isoform_id'] for rec in isoforms if rec.get('isoform_id')}
                matched_map = {}
                for matrix_id in matrix_ids:
                    resolved = resolve_isoform_id(matrix_id, transcript_structures)
                    if not resolved:
                        continue
                    matched_map[matrix_id] = alias_to_transcript.get(resolved, resolved)
                matched_matrix_ids = set(matched_map.keys())
                matched_transcript_ids = set(matched_map.values())
                matrix_only = sorted(matrix_ids - matched_matrix_ids)
                transcript_only = sorted(transcript_ids - matched_transcript_ids)
                overlap_examples = sorted(matched_matrix_ids)[:5]
                print(f"QC matrix isoforms: {len(matrix_ids)}")
                print(f"QC transcript_associations isoforms: {len(transcript_ids)}")
                print(f"QC overlap (matrix->transcript): {len(matched_matrix_ids)}")
                print(f"QC matrix-only: {len(matrix_only)}")
                print(f"QC transcript-only: {len(transcript_only)}")
                if overlap_examples:
                    print("QC overlap examples:", overlap_examples)
                if matrix_only:
                    print("QC matrix-only examples:", matrix_only[:5])
                if transcript_only:
                    print("QC transcript-only examples:", transcript_only[:5])
            raise ValueError("No isoforms with transcript structures available to plot.")
        grouped_structures = group_structures(
            plotted,
            target_gene,
            cluster_features,
            cluster_strategy,
            cluster_similarity_threshold,
            min_split_fraction,
            include_introns=cluster_isoforms,
            report=cluster_isoforms
        )
    else:
        missing_structures = 0
        plotted = []
        seen = set()
        for group in grouped_structures:
            for structure in group:
                for iso in structure['items']:
                    iso_key = iso.get('isoform_id', iso.get('resolved_id'))
                    if iso_key in seen:
                        continue
                    seen.add(iso_key)
                    plotted.append(iso)

    if cluster_isoforms:
        for group in grouped_structures:
            group.sort(key=structure_max_length, reverse=True)
            for structure in group:
                structure['items'].sort(
                    key=lambda iso: (isoform_length(iso), iso.get('count', 0)),
                    reverse=True
                )
        grouped_structures.sort(
            key=lambda group: (group_max_length(group), sum(item['weight'] for item in group)),
            reverse=True
        )

    group_labels = {}
    if cluster_isoforms:
        for group_idx, group in enumerate(grouped_structures):
            group_labels[group_idx] = select_cluster_label(group)

    if debug_isoforms:
        debug_set = set(debug_isoforms)
        debug_hits = []
        debug_lookup = {}
        for structure_idx, group in enumerate(grouped_structures):
            for item in group:
                for iso in item['items']:
                    debug_lookup[iso['resolved_id']] = iso
                    if iso.get('var'):
                        debug_lookup[iso['var']] = iso
                    if iso['resolved_id'] in debug_set or iso.get('var') in debug_set:
                        debug_hits.append({
                            'iso_id': iso['resolved_id'],
                            'var': iso.get('var'),
                            'structure_key': item['structure_key'],
                            'group_index': structure_idx,
                            'weight': item['weight']
                        })
        if debug_hits:
            print("Debug isoform clustering:")
            for hit in debug_hits:
                print(
                    f"  {hit['iso_id']} (var={hit['var']}) group={hit['group_index']} "
                    f"weight={hit['weight']} key={'|'.join(hit['structure_key'])}"
                )
            if len(debug_isoforms) >= 2:
                iso_a = debug_lookup.get(debug_isoforms[0])
                iso_b = debug_lookup.get(debug_isoforms[1])
                if iso_a and iso_b:
                    def exon_only(tokens):
                        return [t for t in tokens if t.startswith('E')]

                    def junctions(tokens):
                        return [f"{tokens[i]}->{tokens[i + 1]}" for i in range(len(tokens) - 1)]

                    raw_a = iso_a['tokens']
                    raw_b = iso_b['tokens']
                    trim_a = strip_terminal_coords(raw_a, target_gene)
                    trim_b = strip_terminal_coords(raw_b, target_gene)
                    cluster_a = filter_exon_tokens(trim_a, target_gene)
                    cluster_b = filter_exon_tokens(trim_b, target_gene)
                    norm_a = iso_a['tokens_norm']
                    norm_b = iso_b['tokens_norm']
                    exon_a = exon_only(trim_a)
                    exon_b = exon_only(trim_b)
                    junc_a = junctions(trim_a)
                    junc_b = junctions(trim_b)

                    def jaccard(a_set, b_set):
                        union = a_set | b_set
                        inter = a_set & b_set
                        return (len(inter) / len(union)) if union else 0.0
                    def overlap(a_set, b_set):
                        if not a_set or not b_set:
                            return 0.0
                        inter = a_set & b_set
                        return len(inter) / min(len(a_set), len(b_set))

                    print("  Debug comparison:")
                    print(f"    A={iso_a['resolved_id']} raw={'|'.join(raw_a)}")
                    print(f"    B={iso_b['resolved_id']} raw={'|'.join(raw_b)}")
                    print(f"    A_trim={'|'.join(trim_a)}")
                    print(f"    B_trim={'|'.join(trim_b)}")
                    print(f"    A_cluster={'|'.join(cluster_a)}")
                    print(f"    B_cluster={'|'.join(cluster_b)}")
                    print(f"    A_exons={'|'.join(exon_a)}")
                    print(f"    B_exons={'|'.join(exon_b)}")
                    print(f"    A_juncs={'|'.join(junc_a)}")
                    print(f"    B_juncs={'|'.join(junc_b)}")
                    print(
                        "    Similarity:",
                        f"cluster_jaccard={jaccard(set(cluster_a), set(cluster_b)):.3f}",
                        f"cluster_overlap={overlap(set(cluster_a), set(cluster_b)):.3f}",
                        f"exon_jaccard={jaccard(set(exon_a), set(exon_b)):.3f}",
                        f"junc_jaccard={jaccard(set(junc_a), set(junc_b)):.3f}",
                        f"a_sub_b={is_subsequence(cluster_a, cluster_b)}",
                        f"b_sub_a={is_subsequence(cluster_b, cluster_a)}",
                    )
                else:
                    missing = [iso for iso in debug_isoforms[:2] if iso not in debug_lookup]
                    print(f"  Debug comparison skipped; missing: {missing}")

    # Track gene order for trans-splicing and compute panel offsets.
    genes_in_plot = OrderedDict()
    genes_in_plot[target_gene] = None
    for iso in plotted:
        for token in iso['tokens']:
            gene_id, _, _ = split_token(token, target_gene)
            if gene_id not in genes_in_plot:
                genes_in_plot[gene_id] = None

    gene_offsets = {}
    x_offset = 0.0
    for gene in genes_in_plot.keys():
        gene_map = gene_maps.get(gene)
        if not gene_map or not gene_map['segments']:
            continue
        gene_offsets[gene] = x_offset
        x_offset += gene_map['total_length'] + gene_gap
    total_width = max(1.0, x_offset - gene_gap)

    # Default isoform spacing to the group spacing when not explicitly provided.
    if isoform_gap is None:
        isoform_gap = group_gap
    if cluster_isoforms:
        isoform_gap = 0.0

    # Expand isoforms into per-read rows with per-isoform colors.
    colors = plt.get_cmap('tab20').colors
    if isoform_colors is None:
        if cluster_isoforms:
            isoform_colors = assign_cluster_colors(grouped_structures, colors)
        else:
            isoform_colors = assign_isoform_colors(grouped_structures, colors)
    rows = []
    y_positions = []
    y = 0.0
    for group_idx, group in enumerate(grouped_structures):
        for structure in group:
            for iso in structure['items']:
                iso_key = iso.get('isoform_id', iso['resolved_id'])
                iso_count = iso['count']
                if isoform_counts is not None:
                    iso_count = isoform_counts.get(iso_key, 0)
                read_count = coerce_read_count(iso_count)
                if read_count <= 0:
                    continue
                for row_idx in range(read_count):
                    rows.append({
                        'isoform': iso,
                        'color': isoform_colors.get(iso_key, colors[0]),
                        'group_idx': group_idx,
                        'group_label': group_labels.get(group_idx),
                        'row_idx': row_idx,
                        'read_count': read_count
                    })
                    y_positions.append(y)
                    y += row_height + row_gap
                y += isoform_gap
        y += max(0.0, group_gap - isoform_gap)
    if rows:
        y -= max(0.0, group_gap - isoform_gap)
        y -= isoform_gap
    total_height = max(1.0, y)

    # Scale figure height to the total stack height for consistent spacing.
    fig_height = max(1.5, total_height * 0.35)
    fig_width = max(10.0, total_width / 3000.0)
    t_render = time.time()
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    row_index = 0
    group_label_y = {}
    if cluster_isoforms and rows:
        group_min = {}
        group_max = {}
        for idx, row in enumerate(rows):
            center = y_positions[idx] + (row_height / 2.0)
            group_idx = row['group_idx']
            group_min[group_idx] = min(center, group_min.get(group_idx, center))
            group_max[group_idx] = max(center, group_max.get(group_idx, center))
        for group_idx in group_min:
            group_label_y[group_idx] = (group_min[group_idx] + group_max[group_idx]) / 2.0
    isoform_first_row = {}
    group_first_row = {}
    exon_labels_drawn = False
    debug_transcript_set = set(debug_transcripts or [])
    debug_transcript_seen = set()
    effective_label_mode = label_mode
    if cluster_isoforms and label_mode == 'all':
        effective_label_mode = 'first'
    for row in rows:
        iso = row['isoform']
        is_first_row = row_index == 0
        # Center the read within its row height.
        y_center = y_positions[row_index] + (row_height / 2.0)
        row_index += 1
        segments = build_isoform_segments(iso['tokens'], exon_lookup, target_gene)
        if iso['resolved_id'] in debug_transcript_set and iso['resolved_id'] not in debug_transcript_seen:
            debug_transcript_seen.add(iso['resolved_id'])
            print(f"Debug transcript {iso['resolved_id']}: tokens={'|'.join(iso['tokens'])}")
            for seg in segments:
                print(
                    "  segment",
                    seg['gene'],
                    seg['type'],
                    seg['start'],
                    seg['end'],
                    seg.get('label')
                )
        span_min = None
        span_max = None
        intron_spans = []
        exon_spans = []
        for seg in segments:
            gene = seg['gene']
            gene_map = gene_maps.get(gene)
            if not gene_map:
                continue
            x0 = map_coord(gene_map, seg['start']) + gene_offsets.get(gene, 0.0)
            x1 = map_coord(gene_map, seg['end']) + gene_offsets.get(gene, 0.0)
            if x1 < x0:
                x0, x1 = x1, x0
            if span_min is None or x0 < span_min:
                span_min = x0
            if span_max is None or x1 > span_max:
                span_max = x1
            if seg['type'] == 'I':
                intron_spans.append({
                    'x0': x0,
                    'x1': x1
                })
            else:
                exon_spans.append({
                    'x0': x0,
                    'x1': x1,
                    'label': seg['label']
                })
        # Draw a thin backbone across the read span.
        if span_min is not None and span_max is not None:
            ax.plot(
                [span_min, span_max],
                [y_center, y_center],
                color='lightgrey',
                linewidth=0.01,
                zorder=1
            )
        # Draw intron retention boxes behind exons.
        ir_face = row['color']
        ir_alpha = None
        for span in intron_spans:
            ax.add_patch(
                plt.Rectangle(
                    (span['x0'], y_center - (row_height / 2.0)),
                    max(1e-6, span['x1'] - span['x0']),
                    row_height,
                    facecolor=ir_face,
                    edgecolor='none',
                    linewidth=0,
                    alpha=ir_alpha,
                    zorder=1.5
                )
            )
        # Draw exon boxes on top of the backbone.
        for span in exon_spans:
            ax.add_patch(
                plt.Rectangle(
                    (span['x0'], y_center - (row_height / 2.0)),
                    max(1e-6, span['x1'] - span['x0']),
                    row_height,
                    facecolor=row['color'],
                    edgecolor='white',
                    linewidth=0.000000,
                    zorder=2
                )
            )
        # Add exon labels only once for the first isoform row shown.
        if is_first_row and not exon_labels_drawn:
            label_y = y_center - row_height * 0.7
            seen_labels = set()
            for span in exon_spans:
                label = span['label']
                if label and label not in seen_labels:
                    seen_labels.add(label)
                    ax.text(
                        (span['x0'] + span['x1']) / 2.0,
                        label_y,
                        label,
                        va='bottom',
                        ha='center',
                        fontsize=5,
                        color='black',
                        zorder=3
                    )
            exon_labels_drawn = True
        if effective_label_mode != 'none':
            if effective_label_mode == 'all':
                show_label = True
                label_y = y_center
            else:
                if cluster_isoforms:
                    show_label = row['row_idx'] == 0 and group_first_row.get(row['group_idx'], True)
                    if show_label:
                        label_y = group_label_y.get(row['group_idx'], y_center)
                else:
                    show_label = row['row_idx'] == 0 and isoform_first_row.get(iso['resolved_id'], True)
                    if show_label:
                        block_height = row_height * row['read_count'] + row_gap * (row['read_count'] - 1)
                        label_y = y_positions[row_index - 1] + (block_height / 2.0)
            if show_label:
                label_text = iso['resolved_id']
                if cluster_isoforms:
                    label_text = row.get('group_label') or label_text
                # Center the isoform label over its read stack.
                ax.text(-total_width * 0.01, label_y, label_text,
                        va='center', ha='right', fontsize=8, color=row['color'])
                if cluster_isoforms:
                    group_first_row[row['group_idx']] = False
                else:
                    isoform_first_row[iso['resolved_id']] = False

    for gene, offset in gene_offsets.items():
        gene_map = gene_maps.get(gene)
        if not gene_map:
            continue
        center = offset + gene_map['total_length'] / 2.0
        ax.text(center, -row_height * 4.0, gene, ha='center', va='bottom', fontsize=9)

    # Tight vertical framing around the read stack.
    ax.set_ylim(-row_height * 5.0, total_height)
    ax.set_xlim(-total_width * 0.1, total_width + total_width * 0.02)
    ax.invert_yaxis()
    ax.axis('off')

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"[timing] render+save plot: {time.time() - t_render:.2f}s")

    if missing_structures:
        print(f"Skipped {missing_structures} isoforms with no transcript_associations mapping.")
    if not return_cluster_labels:
        return plotted
    cluster_label_map = {}
    cluster_index_map = {}
    for group_idx, group in enumerate(grouped_structures):
        label = group_labels.get(group_idx)
        for structure in group:
            for iso in structure['items']:
                iso_key = iso.get('isoform_id', iso.get('resolved_id'))
                resolved_id = iso.get('resolved_id')
                if cluster_isoforms:
                    if iso_key:
                        cluster_label_map[iso_key] = label
                        cluster_index_map[iso_key] = group_idx
                    if resolved_id:
                        cluster_label_map[resolved_id] = label
                        cluster_index_map[resolved_id] = group_idx
                else:
                    if iso_key and resolved_id:
                        cluster_label_map[iso_key] = resolved_id
                        cluster_index_map[iso_key] = group_idx
    return plotted, cluster_label_map, cluster_index_map


def main():
    parser = argparse.ArgumentParser(
        description="Plot single-gene long-read isoform structures from an isoform matrix."
    )
    parser.add_argument('--h5ad', nargs='+', required=True,
                        help='Isoform h5ad file(s) or 10x mtx folder(s)')
    parser.add_argument('--gene', default=None, help='Target gene (e.g., ENSG...)')
    parser.add_argument('--genes', nargs='*', default=None,
                        help='Optional list of target genes (overrides --gene)')
    parser.add_argument('--genes-file', default=None,
                        help='Optional file with one gene ID per line')
    parser.add_argument('--transcript-associations', required=True, nargs='+',
                        help='Path(s) to gff-output/transcript_associations.txt')
    parser.add_argument('--gene-model', default=None,
                        help='Path to gene_model TSV reference')
    parser.add_argument('--gene-model-index', default=None,
                        help='Path to pre-built gene index directory (SQLite).')
    parser.add_argument('--rebuild-index', action='store_true',
                        help='Delete existing transcript/gene-model/h5ad indexes and rebuild')
    parser.add_argument('--out', default='isoform_structures.pdf', help='Output figure path')
    parser.add_argument('--groupby', default=None,
                        help='Cell annotation column or barcode-group TSV for filtering')
    parser.add_argument('--group', nargs='*', default=None,
                        help='Optional values in groupby column to keep')
    parser.add_argument('--groupby-rev', action=argparse.BooleanOptionalAction, default=True,
                        help='Reverse-complement barcodes from groupby TSV (default: True)')
    parser.add_argument('--groupby-sample', default=None,
                        help='Optional sample name suffix to match when groupby is a TSV/TXT file')
    parser.add_argument('--max-isoforms', type=int, default=200,
                        help='Maximum isoforms to display')
    parser.add_argument('--min-count', type=float, default=1,
                        help='Minimum summed count to include an isoform')
    parser.add_argument('--cluster-mode', choices=['full', 'base', 'block'], default='block',
                        help='How to group isoform structures for clustering')
    parser.add_argument('--cluster-strategy', choices=['substring', 'subsequence', 'overlap', 'jaccard', 'feature', 'structure', 'count'], default='subsequence',
                        help='How to order isoforms (substring/jaccard similarity, feature splits, structure grouping, or count)')
    parser.add_argument('--cluster-similarity-threshold', '--cluster-overlap-threshold',
                        dest='cluster_similarity_threshold', type=float, default=0.85,
                        help='Minimum similarity for clustering (subsequence when --cluster-isoforms)')
    parser.add_argument('--cluster-features', choices=['tokens', 'junctions', 'both'], default='tokens',
                        help='Features to use for clustering')
    parser.add_argument('--min-split-fraction', type=float, default=0.05,
                        help='Minimum fraction for a structural split (feature clustering)')
    parser.add_argument('--cluster-isoforms', action='store_true',
                        help='Color isoforms by cluster and label using the longest isoform in each cluster')
    parser.add_argument('--no-cluster-isoforms', action='store_false', dest='cluster_isoforms',
                        help='Disable isoform clustering and display isoforms independently')
    parser.add_argument('--inspect-isoforms', nargs='*', default=None,
                        help='Isoform IDs to print structure details for')
    parser.add_argument('--intron-scale', type=float, default=0.2,
                        help='Scale factor for intron lengths (0-1)')
    parser.add_argument('--row-height', type=float, default=0.0125,
                        help='Height of each read row')
    parser.add_argument('--row-gap', type=float, default=0.0,
                        help='Gap between read rows')
    parser.add_argument('--group-gap', type=float, default=0.1,
                        help='Gap between clustered groups')
    parser.add_argument('--isoform-gap', type=float, default=None,
                        help='Gap between isoform stacks (defaults to group-gap)')
    parser.add_argument('--label-mode', choices=['first', 'none', 'all'], default='first',
                        help='When to display isoform labels')
    parser.add_argument('--debug-transcripts', nargs='*', default=None,
                        help='Transcript IDs to print segment debugging for')
    parser.add_argument('--debug-isoforms', nargs='*', default=None,
                        help='Isoform IDs to print clustering diagnostics for')
    parser.set_defaults(cluster_isoforms=True)
    args = parser.parse_args()
    log_path, log_handle, stdout_tee, stderr_tee = _init_logging(args.out)
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    sys.stdout = stdout_tee
    sys.stderr = stderr_tee
    print(f"Log file: {log_path}")
    print(f"[debug] isoform_structure_view path: {__file__}")
    print("[debug] isoform_structure_view build: 20260206-mask-v2")
    start_time = time.time()
    try:
        def resolve_genes():
            gene_list = []
            if args.genes_file:
                with open(args.genes_file, 'r') as handle:
                    for line in handle:
                        token = line.strip()
                        if token:
                            gene_list.append(token)
            if args.genes:
                for token in args.genes:
                    if token and ',' in token:
                        gene_list.extend([item.strip() for item in token.split(',') if item.strip()])
                    elif token:
                        gene_list.append(token)
            if not gene_list and args.gene:
                gene_list = [args.gene]
            return gene_list

        genes = resolve_genes()
        if not genes:
            parser.error("Must provide --gene, --genes, or --genes-file")
        print(f"[run] Genes: {genes}")
        print(f"[run] Groupby: {args.groupby or 'none'} | Group values: {args.group or ['all']}")

        if args.group:
            expanded = []
            for value in args.group:
                if value and ',' in value:
                    expanded.extend([item.strip() for item in value.split(',') if item.strip()])
                elif value:
                    expanded.append(value)
            args.group = expanded or None
        if args.cluster_strategy is None:
            args.cluster_strategy = 'substring' if args.cluster_isoforms else 'jaccard'
        if args.cluster_strategy in ('subsequence', 'overlap'):
            args.cluster_strategy = 'substring'
        if args.cluster_similarity_threshold is None:
            args.cluster_similarity_threshold = 0.85 if args.cluster_isoforms else 0.4

        groupby_is_file = bool(args.groupby) and (args.groupby.endswith('.tsv') or args.groupby.endswith('.txt'))

        transcript_paths = args.transcript_associations
        if len(transcript_paths) == 1:
            transcript_paths = transcript_paths * len(args.h5ad)
        elif len(transcript_paths) != len(args.h5ad):
            parser.error("--transcript-associations must be one path or match --h5ad count")
        missing_transcripts = [path for path in transcript_paths if not os.path.exists(path)]
        if missing_transcripts:
            raise FileNotFoundError(
                "Missing transcript_associations file(s): "
                + ", ".join(missing_transcripts)
            )
        print(f"[run] Transcript associations: {transcript_paths}")

        _maybe_rebuild_indexes(args.rebuild_index, args.h5ad, transcript_paths, args.gene_model)

        index_dir = args.gene_model_index
        if not index_dir and args.gene_model:
            index_dir = detect_gene_index(args.gene_model) or _gene_model_index_dir(args.gene_model)

        if index_dir:
            print(f"[index] Using SQLite index: {index_dir}")
        else:
            print("[index] No SQLite index detected; reading transcript_associations text.")

        gene_segments = defaultdict(list)
        exon_lookup = {}
        gene_model_db = _ensure_gene_model_db(args.gene_model, index_dir=index_dir) if args.gene_model else None
        if gene_model_db and os.path.exists(gene_model_db):
            t_gene_model = time.time()
            for gene in genes:
                t_gene = time.time()
                segments, lookup = read_gene_model_from_sqlite(gene_model_db, gene)
                for key, value in segments.items():
                    gene_segments[key].extend(value)
                exon_lookup.update(lookup)
                print(f"[timing] gene_model sqlite {gene}: {time.time() - t_gene:.2f}s")
            print(f"[timing] gene_model sqlite total: {time.time() - t_gene_model:.2f}s")
        else:
            if not args.gene_model:
                raise ValueError("--gene-model required when not using SQLite index")
            t_gene_model = time.time()
            gene_segments, exon_lookup = read_gene_model(args.gene_model)
            print(f"[timing] gene_model tsv total: {time.time() - t_gene_model:.2f}s")
        gene_maps = build_gene_maps(gene_segments, args.intron_scale)

        transcript_structures_by_gene = {}
        for gene in genes:
            union_structures = {}
            for assoc_path in transcript_paths:
                t_assoc = time.time()
                structures, _ = load_transcript_associations_auto(assoc_path, gene, index_dir=None)
                print(f"[timing] transcript_associations load ({os.path.basename(assoc_path)}): {time.time() - t_assoc:.2f}s")
                if structures:
                    union_structures.update(structures)
            transcript_structures_by_gene[gene] = union_structures
            if not union_structures and len(genes) == 1:
                raise ValueError(
                    f"No transcript structures found for gene {gene}. "
                    f"Checked: {', '.join(transcript_paths)}"
                )

        def _infer_basename(input_path):
            if os.path.isdir(input_path):
                return os.path.basename(input_path)
            if input_path.endswith('.mtx') or input_path.endswith('.mtx.gz'):
                return os.path.basename(os.path.dirname(input_path))
            return os.path.splitext(os.path.basename(input_path))[0]

        def _resolve_output_path(input_path, gene):
            if len(genes) > 1:
                out_dir = args.out
                if args.out.endswith('.pdf'):
                    out_dir = os.path.dirname(args.out) or "."
                os.makedirs(out_dir, exist_ok=True)
                basename = _infer_basename(input_path)
                return os.path.join(out_dir, f"{basename}_{gene}_isoform_structures.pdf")
            if os.path.isdir(input_path) or input_path.endswith('.mtx') or input_path.endswith('.mtx.gz'):
                out_dir = args.out
                if args.out.endswith('.pdf'):
                    out_dir = os.path.dirname(args.out) or "."
                os.makedirs(out_dir, exist_ok=True)
                basename = _infer_basename(input_path)
                return os.path.join(out_dir, f"{basename}_isoform_structures.pdf")
            return args.out

        input_paths = args.h5ad
        use_counts_cache = True
        allow_cache_build = True
        if len(input_paths) == 1:
            path = input_paths[0]
            print(f"Processing single file: {path}")
            groupby_sample = args.groupby_sample
            if groupby_is_file and not groupby_sample:
                groupby_sample = _infer_sample_name(path)
            adata = None
            gene_index = None
            mask = None
            cache_path = _h5ad_counts_cache_path(
                path, args.groupby, args.group, args.groupby_rev, groupby_sample
            )
            cache_sources = [path]
            if groupby_is_file:
                cache_sources.append(args.groupby)
            cache_ready = use_counts_cache and not _is_index_stale(cache_path + ".counts.npy", cache_sources)
            try:
                if not cache_ready:
                    t_index = time.time()
                    gene_index = _load_h5ad_gene_index(path, index_dir=index_dir)
                    print(f"[timing] h5ad index load: {time.time() - t_index:.2f}s")
                    t_open = time.time()
                    adata = ad.read_h5ad(path, backed='r')
                    print(f"[timing] h5ad open: {time.time() - t_open:.2f}s")
                    mask = _compute_groupby_mask(
                        adata,
                        args.groupby,
                        args.group,
                        groupby_rev=args.groupby_rev,
                        groupby_sample=groupby_sample,
                    )
                for gene in genes:
                    transcript_structures = transcript_structures_by_gene.get(gene, {})
                    if not transcript_structures:
                        if len(genes) == 1:
                            raise ValueError(f"No transcript structures found for gene {gene}.")
                        print(f"No transcript structures found for gene {gene}; skipping.")
                        continue
                    isoform_records = load_isoform_counts_indexed(
                        path,
                        gene,
                        index_dir=index_dir,
                        groupby=args.groupby,
                        group_values=args.group,
                        groupby_rev=args.groupby_rev,
                        groupby_sample=groupby_sample,
                        adata=adata,
                        gene_index=gene_index,
                        mask=mask,
                        use_counts_cache=use_counts_cache,
                        cache_build=allow_cache_build and not cache_ready,
                    )
                    isoform_records = filter_records_by_structures(isoform_records, transcript_structures)
                    if not isoform_records:
                        if len(genes) == 1:
                            raise ValueError("No isoform counts found for the target gene.")
                        print(f"No isoform counts found for {gene}; skipping.")
                        continue

                    t_cluster_map = time.time()
                    cluster_label_map, _plotted_ids = build_cluster_label_map(
                        isoform_records,
                        transcript_structures,
                        gene,
                        args.cluster_mode,
                        args.cluster_features,
                        args.cluster_strategy,
                        args.cluster_similarity_threshold,
                        args.min_split_fraction,
                        args.max_isoforms,
                        args.min_count,
                        args.cluster_isoforms,
                        inspect_isoforms=args.inspect_isoforms
                    )
                    print(f"[timing] cluster_label_map: {time.time() - t_cluster_map:.2f}s")
                    output_path = _resolve_output_path(path, gene)
                    print(f"Generating plot for {gene}...")
                    plotted, plot_cluster_label_map, plot_cluster_index_map = plot_isoform_structures(
                        isoform_records,
                        transcript_structures,
                        exon_lookup,
                        gene_maps,
                        args.cluster_mode,
                        gene,
                        args.intron_scale,
                        output_path,
                        max_isoforms=args.max_isoforms,
                        min_count=args.min_count,
                        row_height=args.row_height,
                        row_gap=args.row_gap,
                        group_gap=args.group_gap,
                        isoform_gap=args.isoform_gap,
                        cluster_strategy=args.cluster_strategy,
                        cluster_similarity_threshold=args.cluster_similarity_threshold,
                        cluster_features=args.cluster_features,
                        label_mode=args.label_mode,
                        min_split_fraction=args.min_split_fraction,
                        debug_isoforms=args.debug_isoforms,
                        debug_transcripts=args.debug_transcripts,
                        transcript_associations_path=transcript_paths[0],
                        cluster_isoforms=args.cluster_isoforms,
                        return_cluster_labels=True
                    )
                    isoform_path, _ = export_plotted_isoforms(
                        plotted, output_path, plot_cluster_label_map, plot_cluster_index_map)
                    if isoform_path:
                        print(f"Saved plotted isoform IDs to: {isoform_path}")
            finally:
                if adata is not None:
                    try:
                        adata.file.close()
                    except AttributeError:
                        pass
        else:
            print(f"Processing {len(input_paths)} files...")
            per_gene_combined = {gene: defaultdict(float) for gene in genes}
            per_gene_record = {gene: {} for gene in genes}
            per_gene_per_file_counts = {gene: {} for gene in genes}

            for path, assoc_path in zip(input_paths, transcript_paths):
                print(f"  Loading {path}...")
                groupby_sample = args.groupby_sample
                if groupby_is_file and not groupby_sample:
                    groupby_sample = _infer_sample_name(path)
                adata = None
                gene_index = None
                mask = None
                cache_path = _h5ad_counts_cache_path(
                    path, args.groupby, args.group, args.groupby_rev, groupby_sample
                )
                cache_sources = [path]
                if groupby_is_file:
                    cache_sources.append(args.groupby)
                cache_ready = use_counts_cache and not _is_index_stale(cache_path + ".counts.npy", cache_sources)
                try:
                    if not cache_ready:
                        t_index = time.time()
                        gene_index = _load_h5ad_gene_index(path, index_dir=index_dir)
                        print(f"[timing] h5ad index load: {time.time() - t_index:.2f}s")
                        t_open = time.time()
                        adata = ad.read_h5ad(path, backed='r')
                        print(f"[timing] h5ad open: {time.time() - t_open:.2f}s")
                        mask = _compute_groupby_mask(
                            adata,
                            args.groupby,
                            args.group,
                            groupby_rev=args.groupby_rev,
                            groupby_sample=groupby_sample,
                        )
                    for gene in genes:
                        transcript_structures = transcript_structures_by_gene.get(gene, {})
                        if not transcript_structures:
                            continue
                        records = load_isoform_counts_indexed(
                            path,
                            gene,
                            index_dir=index_dir,
                            groupby=args.groupby,
                            group_values=args.group,
                            groupby_rev=args.groupby_rev,
                            groupby_sample=groupby_sample,
                            adata=adata,
                            gene_index=gene_index,
                            mask=mask,
                            use_counts_cache=use_counts_cache,
                            cache_build=allow_cache_build and not cache_ready,
                        )
                        records = filter_records_by_structures(records, transcript_structures)
                        per_gene_per_file_counts[gene][path] = {
                            rec['isoform_id']: rec['count'] for rec in records
                        }
                        for rec in records:
                            key = rec['isoform_id']
                            per_gene_combined[gene][key] += rec['count']
                            if key not in per_gene_record[gene]:
                                per_gene_record[gene][key] = rec
                finally:
                    if adata is not None:
                        try:
                            adata.file.close()
                        except AttributeError:
                            pass

            for gene in genes:
                transcript_structures = transcript_structures_by_gene.get(gene, {})
                if not transcript_structures:
                    if len(genes) == 1:
                        raise ValueError(f"No transcript structures found for gene {gene}.")
                    print(f"No transcript structures found for gene {gene}; skipping.")
                    continue

                record_by_isoform = per_gene_record.get(gene, {})
                combined_counts = per_gene_combined.get(gene, {})
                per_file_counts = per_gene_per_file_counts.get(gene, {})

                combined_records = []
                for key, rec in record_by_isoform.items():
                    merged = dict(rec)
                    merged['count'] = combined_counts.get(key, 0)
                    combined_records.append(merged)
                if not combined_records:
                    if len(genes) == 1:
                        raise ValueError("No isoform counts found for the target gene.")
                    print(f"No isoform counts found for {gene}; skipping.")
                    continue

                grouped_structures = None
                isoform_colors = None
                plotted, _, _ = build_plotted_isoforms(
                    combined_records, transcript_structures, gene,
                    args.cluster_mode, args.cluster_features,
                    args.min_count, args.max_isoforms
                )
                if plotted:
                    grouped_structures = group_structures(
                        plotted, gene, args.cluster_features,
                        args.cluster_strategy, args.cluster_similarity_threshold,
                        args.min_split_fraction, include_introns=args.cluster_isoforms,
                        report=args.cluster_isoforms
                    )
                    if args.cluster_isoforms:
                        isoform_colors = assign_cluster_colors(grouped_structures, plt.get_cmap(DEFAULT_COLORMAP).colors)
                    else:
                        isoform_colors = assign_isoform_colors(grouped_structures, plt.get_cmap(DEFAULT_COLORMAP).colors)

                print(f"Generating individual plots for {gene}...")
                isoform_path = None
                for path, assoc_path in zip(input_paths, transcript_paths):
                    out_path = _resolve_output_path(path, gene)
                    plotted, plot_cluster_label_map, plot_cluster_index_map = plot_isoform_structures(
                        record_by_isoform.values(),
                        transcript_structures,
                        exon_lookup,
                        gene_maps,
                        args.cluster_mode,
                        gene,
                        args.intron_scale,
                        out_path,
                        max_isoforms=args.max_isoforms,
                        min_count=args.min_count,
                        row_height=args.row_height,
                        row_gap=args.row_gap,
                        group_gap=args.group_gap,
                        isoform_gap=args.isoform_gap,
                        cluster_strategy=args.cluster_strategy,
                        cluster_similarity_threshold=args.cluster_similarity_threshold,
                        cluster_features=args.cluster_features,
                        label_mode=args.label_mode,
                        min_split_fraction=args.min_split_fraction,
                        debug_isoforms=args.debug_isoforms,
                        debug_transcripts=args.debug_transcripts,
                        grouped_structures=grouped_structures,
                        isoform_colors=isoform_colors,
                        isoform_counts=per_file_counts.get(path, {}),
                        transcript_associations_path=assoc_path,
                        cluster_isoforms=args.cluster_isoforms,
                        return_cluster_labels=True
                    )
                    isoform_path, _ = export_plotted_isoforms(
                        plotted, out_path, plot_cluster_label_map, plot_cluster_index_map)
                if isoform_path:
                    print(f"Saved plotted isoform IDs to: {isoform_path}")
    finally:
        elapsed = time.time() - start_time
        print(f"[run] Elapsed time (s): {elapsed:.2f}")
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        log_handle.close()


if __name__ == '__main__':
    main()
