"""Transparent compressed/uncompressed IO for long-read workflow artifacts.

The workflow writes several large text products per run: the per-sample molecule table, the collapse
catalog and its structure map, the protein summary, the coding regions, the per-sample counts
exports, and three sequence FASTAs. Uncompressed they total roughly 674 MB per sample; gzip takes
them to about 90 MB (measured ratios 6.1% to 18.2%).

Two rules make compressing them safe:

1. **Every reader accepts either form.** ``smart_open`` and ``resolve`` take a logical path and find
   whichever of ``<path>`` or ``<path>.gz`` exists. Runs produced before this change keep working,
   and a user who gunzips a file by hand keeps working.

2. **FASTA is compressed with BGZF, not plain gzip.** ``isv_web/data_api.py`` builds a
   transcript_id -> byte-offset index over the sequence FASTAs and calls ``fh.seek(start)`` to pull
   one record. A plain gzip stream cannot be seeked by byte offset, so plain gzip would silently
   break sequence lookup in the viewer. BGZF is gzip-compatible on read (``zcat``, ``gzip -d`` and
   ``gzip.open`` all work) and supports random access through ``pysam.BGZFile``. Measured cost over
   plain gzip: 1.8% more bytes.
"""

from __future__ import annotations

import os
import gzip
import shutil


#: Artifacts that are randomly accessed by byte offset and therefore need BGZF, not plain gzip.
SEEKABLE_SUFFIXES = ('.fasta', '.fa')


def resolve(path, missing_ok=False):
    """Return the path that actually exists, trying ``path`` then ``path + '.gz'``.

    A caller names the logical file (``protein_summary.txt``) and gets whichever form is on disk.
    Raises FileNotFoundError when neither exists, unless ``missing_ok``.
    """
    path = str(path)
    if os.path.exists(path):
        return path
    if path.endswith('.gz'):
        plain = path[:-3]
        if os.path.exists(plain):
            return plain
    else:
        gz = path + '.gz'
        if os.path.exists(gz):
            return gz
    if missing_ok:
        return None
    raise FileNotFoundError(f"Neither {path} nor its .gz counterpart exists")


def exists(path):
    return resolve(path, missing_ok=True) is not None


def smart_open(path, mode='rt', **kwargs):
    """Open ``path`` or ``path.gz``, decompressing transparently.

    BGZF files are valid gzip, so ``gzip.open`` reads them sequentially without any special case.
    Use ``open_seekable`` when you need random access instead.
    """
    resolved = resolve(path)
    if resolved.endswith('.gz'):
        if 'b' not in mode and 't' not in mode:
            mode = mode + 't'
        return gzip.open(resolved, mode, **kwargs)
    return open(resolved, mode.replace('t', '') or 'r', **kwargs)


def open_seekable(path):
    """Open for RANDOM ACCESS by byte offset, compressed or not.

    Returns a handle supporting ``seek``/``read``. For a BGZF file this is a ``pysam.BGZFile``,
    whose ``seek`` takes a virtual offset produced by ``tell`` on the same class, so an index built
    against a BGZF file must be built with this function too. For a plain file it is a normal
    handle and offsets are ordinary byte positions.
    """
    resolved = resolve(path)
    if resolved.endswith('.gz'):
        import pysam
        return pysam.BGZFile(resolved, 'rb')
    return open(resolved, 'rb')


def is_bgzf(path):
    """True when the file carries the BGZF extra-field signature in its gzip header."""
    resolved = resolve(path, missing_ok=True)
    if not resolved or not resolved.endswith('.gz'):
        return False
    try:
        with open(resolved, 'rb') as handle:
            head = handle.read(18)
        # gzip magic, FEXTRA set, and the BC subfield BGZF uses
        return len(head) >= 18 and head[:2] == b'\x1f\x8b' and (head[3] & 0x04) and head[12:14] == b'BC'
    except OSError:
        return False


#: Below this size gzip's header costs more than it saves; the unit test produced a 34-byte
#: file that grew to 62 bytes. Small artifacts are left alone.
MIN_COMPRESS_BYTES = 4096


def compress(path, level=6, remove_source=True, min_bytes=MIN_COMPRESS_BYTES, log=None):
    """Compress a freshly written file in place, returning the resulting path.

    FASTA gets BGZF so the viewer's offset index still works; everything else gets plain gzip.
    Already-compressed input is returned unchanged. A failure returns the original path rather than
    raising, because losing an output to a compression error is worse than leaving it uncompressed.
    """
    path = str(path)
    if path.endswith('.gz') or not os.path.exists(path):
        return path
    before = os.path.getsize(path)
    if before < min_bytes:
        return path
    target = path + '.gz'
    try:
        if path.endswith(SEEKABLE_SUFFIXES):
            import pysam
            pysam.tabix_compress(path, target, force=True)
        else:
            with open(path, 'rb') as src, gzip.open(target, 'wb', compresslevel=level) as dst:
                shutil.copyfileobj(src, dst, length=1024 * 1024)
    except Exception as exc:                      # noqa: BLE001 - never lose an output
        if log:
            log(f"[compress] {os.path.basename(path)} left uncompressed ({type(exc).__name__}: {exc})")
        if os.path.exists(target):
            try:
                os.replace(target, target + '.failed')
            except OSError:
                pass
        return path
    after = os.path.getsize(target)
    if remove_source:
        os.remove(path)
    if log:
        kind = 'bgzf' if path.endswith(SEEKABLE_SUFFIXES) else 'gzip'
        log(f"[compress] {os.path.basename(path)} {before:,} -> {after:,} bytes "
            f"({100.0 * after / before:.1f}%, {kind})")
    return target


def compress_many(paths, log=None):
    """Compress each path that exists. Returns the list of resulting paths."""
    out = []
    for p in paths:
        if p and os.path.exists(str(p)):
            out.append(compress(p, log=log))
    return out
