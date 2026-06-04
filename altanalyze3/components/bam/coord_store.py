#!/usr/bin/env python3
"""Indexed isoform-structure -> genomic-coordinate store.

Built at the SAMPLE-LEVEL gff_process stage (``process_isoform`` in gff_process.py), which is the
first place a read's exon blocks are resolved into an *isoform structure string* (the canonical,
cross-sample key collapse decides on). At that moment ``chr, strand, exons`` (the genomic
coordinates) and the source GFF ``info`` line are all in hand, so the unique
``(gene, structure) -> genomic record`` map is materialized for free, deduplicated by structure
(thousands of molecules collapse to one structure -> one stored row).

At collapse time, ``stage_protein`` retrieves the genomic exon records for the KEPT exemplar
structures by an indexed point lookup -- turning what used to be a full sequential scan of every
per-sample read GFF (the dominant cost of writing combined.gff) into bounded
``SELECT ... WHERE (gene,structure) IN (...)`` queries, then writes combined.gff.gz.

Design (mirrors the existing per-chunk counts SQLite pattern in isoform_structure_extract.py):
  * one SQLite DB per sample, ``<sample>.struct_coords.sqlite``, beside the sample's
    gff-output/transcript_associations.txt
  * table ``structures(gene TEXT, structure TEXT, chrom TEXT, strand TEXT, source TEXT,
    blob BLOB, PRIMARY KEY(gene, structure))`` -- the composite PK gives the indexed lookup
  * ``blob`` = varint-delta-packed exon (start,end) pairs in the SAME order the source GFF lists
    them, so the emitted records are byte-identical to a GFF copy.

Tolerance: if the store is absent (older runs) or a structure is missing, the consumer falls back
to scanning the source GFF -- so this is a pure optimization, never a correctness dependency.
"""
import sqlite3
import struct
from pathlib import Path


_CREATE_SQL = (
    "CREATE TABLE IF NOT EXISTS structures ("
    "gene TEXT NOT NULL, "
    "structure TEXT NOT NULL, "
    "chrom TEXT NOT NULL, "
    "strand TEXT NOT NULL, "
    "source TEXT NOT NULL, "
    "blob BLOB NOT NULL, "
    "PRIMARY KEY (gene, structure)"
    ")"
)
_INSERT_SQL = (
    "INSERT INTO structures (gene, structure, chrom, strand, source, blob) "
    "VALUES (?, ?, ?, ?, ?, ?) "
    "ON CONFLICT(gene, structure) DO NOTHING"
)


def _pragmas(conn):
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA locking_mode=EXCLUSIVE")


# ----------------------------------------------------------------- varint pack
def _encode_varint(value, out):
    while True:
        b = value & 0x7F
        value >>= 7
        if value:
            out.append(b | 0x80)
        else:
            out.append(b)
            return


def _decode_varint(buf, pos):
    result = 0
    shift = 0
    while True:
        b = buf[pos]
        pos += 1
        result |= (b & 0x7F) << shift
        if not (b & 0x80):
            return result, pos
        shift += 7


def pack_exons(exons):
    """Varint-pack a list of (start, end) exon pairs, preserving input order (combined.gff copies
    exon lines in source order). mode 0 = non-negative delta chain from 0 (monotonic case);
    mode 1 = absolute off a base when any delta would be negative (out-of-order/overlapping)."""
    n = len(exons)
    monotonic = True
    prev_end = 0
    for s, e in exons:
        if s < prev_end or e < s:
            monotonic = False
            break
        prev_end = e
    out = bytearray()
    if monotonic:
        out.append(0)
        _encode_varint(n, out)
        prev_end = 0
        for s, e in exons:
            _encode_varint(s - prev_end, out)
            _encode_varint(e - s, out)
            prev_end = e
        return bytes(out)
    base = min(s for s, _ in exons)
    out.append(1)
    out += struct.pack('<q', base)
    _encode_varint(n, out)
    for s, e in exons:
        _encode_varint(s - base, out)
        _encode_varint(e - s, out)
    return bytes(out)


def unpack_exons(blob):
    mode = blob[0]
    pos = 1
    if mode == 0:
        n, pos = _decode_varint(blob, pos)
        exons = []
        prev_end = 0
        for _ in range(n):
            d, pos = _decode_varint(blob, pos)
            length, pos = _decode_varint(blob, pos)
            s = prev_end + d
            e = s + length
            exons.append((s, e))
            prev_end = e
        return exons
    elif mode == 1:
        (base,) = struct.unpack_from('<q', blob, pos)
        pos += 8
        n, pos = _decode_varint(blob, pos)
        exons = []
        for _ in range(n):
            ds, pos = _decode_varint(blob, pos)
            length, pos = _decode_varint(blob, pos)
            s = base + ds
            e = s + length
            exons.append((s, e))
        return exons
    raise ValueError(f"unknown coord blob mode {mode}")


# ----------------------------------------------------------------- writer
class StructCoordWriter:
    """Batched writer for one sample's structure->coords SQLite DB."""

    def __init__(self, db_path, buffer_limit=10000):
        self.db_path = str(db_path)
        self.conn = sqlite3.connect(self.db_path)
        _pragmas(self.conn)
        self.conn.execute(_CREATE_SQL)
        self.conn.execute("BEGIN")
        self._buf = []
        self._limit = buffer_limit

    def add(self, gene, structure, chrom, strand, source, exons):
        self._buf.append((gene, structure, chrom, strand, source, pack_exons(list(exons))))
        if len(self._buf) >= self._limit:
            self.conn.executemany(_INSERT_SQL, self._buf)
            self._buf.clear()

    def close(self):
        if self._buf:
            self.conn.executemany(_INSERT_SQL, self._buf)
            self._buf.clear()
        self.conn.commit()
        self.conn.close()


# ----------------------------------------------------------------- reader
def open_reader(db_path):
    """Open a structure-coords DB read-only; returns a connection or None if absent."""
    p = Path(db_path)
    if not p.exists():
        return None
    return sqlite3.connect(f"file:{p}?mode=ro", uri=True)


def _build_info(gene, transcript_id):
    """Reconstruct the exact col-9 string write_gff_isoform emits: gene_id then transcript_id."""
    parts = []
    if gene:
        parts.append(f'gene_id "{gene}"')
    parts.append(f'transcript_id "{transcript_id}"')
    return ';'.join(parts) + ';'


def format_records(chrom, source, strand, exons, info):
    """Build the multi-line GFF block (exon lines in stored order + a transcript summary line),
    byte-identical to isoform_structure_extract.write_gff_isoform."""
    lines = []
    for s, e in exons:
        lines.append(f"{chrom}\t{source}\texon\t{s}\t{e}\t.\t{strand}\t.\t{info}\n")
    tx_start = min(s for s, _ in exons)
    tx_end = max(e for _, e in exons)
    lines.append(f"{chrom}\t{source}\ttranscript\t{tx_start}\t{tx_end}\t.\t{strand}\t.\t{info}\n")
    return "".join(lines)


def fetch_structures(conn, gene_structure_pairs, transcript_id_for=None, chunk=None):
    """Yield (gene, structure, chrom, gff_block) for the requested (gene, structure) pairs.

    transcript_id_for: optional callable (gene, structure) -> transcript_id to stamp into the
    emitted col-9 (the kept exemplar id). If None, the transcript_id is set to the structure
    string. Missing pairs are simply not yielded (caller may fall back to a GFF scan).

    Implementation note: queries via a single-key point lookup ``WHERE gene=? AND structure=?``
    per pair, which the composite PRIMARY KEY index (sqlite_autoindex_structures_1) serves
    directly. The seemingly-natural ``WHERE (gene,structure) IN (VALUES ...)`` form does NOT use
    the index in SQLite -- it forces a full table scan per batch (verified via EXPLAIN QUERY PLAN),
    making it ~180x slower. ``chunk`` is accepted for backwards compatibility but ignored.
    """
    q = "SELECT chrom, strand, source, blob FROM structures WHERE gene=? AND structure=?"
    for gene, structure in gene_structure_pairs:
        row = conn.execute(q, (gene, structure)).fetchone()
        if row is None:
            continue
        chrom, strand, source, blob = row
        exons = unpack_exons(blob)
        tid = transcript_id_for(gene, structure) if transcript_id_for else structure
        info = _build_info(gene, tid)
        yield (gene, structure, chrom, format_records(chrom, source, strand, exons, info))
