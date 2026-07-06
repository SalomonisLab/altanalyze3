"""Reference-data resolution and download for the SNAF full pipeline.

The full T-antigen pipeline needs the SNAF reference bundle (``Alt91_db/`` plus a
GTEx junction control h5ad under ``controls/``). This module *locates* that
bundle, *validates* it with actionable errors, and can *fetch* it directly from
the public mirror. The design goal is that a missing or partial reference fails
LOUDLY with the exact command to fix it -- never a silent skip or a cryptic
downstream ``KeyError``.

Public mirror (from the SNAF docs, ``docs/installation.rst``)::

    curl -o download.tar.gz http://altanalyze.org/SNAF/download.tar.gz
    tar -xzf download.tar.gz          # extracts to ./data/{Alt91_db,controls}/

Note the tarball nests everything under a top-level ``data/`` directory, so the
resolved reference root may be ``<db_dir>`` or ``<db_dir>/data`` -- both are
handled transparently here.
"""
import os
import logging
import tarfile
import urllib.request

logger = logging.getLogger(__name__)

REFERENCE_URL = 'http://altanalyze.org/SNAF/download.tar.gz'
REFERENCE_APPROX_GB = 2.7

# Files the pipeline actually reads, relative to the reference root.
_ALT91_REQUIRED = (
    'Hs_Ensembl_exon_add_col.txt',   # exon_table
    'mRNA-ExonIDs.txt',              # transcript_db
    'Hs_gene-seq-2000_flank.fa',     # gene-flank fasta
)
_CONTROL_BY_MODE = {
    'count': 'GTEx_junction_counts.h5ad',
    'psi':   'GTEx_junction_psi.h5ad',
}


class ReferenceError(RuntimeError):
    """Raised when the SNAF reference bundle is missing or incomplete."""


def control_filename(gtex_mode='count'):
    return _CONTROL_BY_MODE.get(gtex_mode, _CONTROL_BY_MODE['count'])


def _has_reference(root, gtex_mode='count'):
    alt = os.path.join(root, 'Alt91_db')
    if not all(os.path.exists(os.path.join(alt, f)) for f in _ALT91_REQUIRED):
        return False
    return os.path.exists(os.path.join(root, 'controls', control_filename(gtex_mode)))


def locate_reference_root(db_dir, gtex_mode='count'):
    """Return the abs path of the dir that actually holds ``Alt91_db/`` +
    ``controls/`` (tolerating the tarball's ``data/`` nesting), or ``None``."""
    for cand in (db_dir, os.path.join(db_dir, 'data')):
        if _has_reference(cand, gtex_mode):
            return os.path.abspath(cand)
    return None


def download_command(db_dir):
    """The exact shell command a user can run to obtain the reference."""
    return ('mkdir -p "{d}" && curl -L -o "{d}/download.tar.gz" {url} '
            '&& tar -xzf "{d}/download.tar.gz" -C "{d}"').format(d=db_dir, url=REFERENCE_URL)


def _download_stream(url, dest):
    """Stream ``url`` to ``dest``, logging progress every ~10%."""
    with urllib.request.urlopen(url) as resp:  # nosec - trusted public mirror
        total = int(resp.headers.get('Content-Length', 0))
        done = 0
        next_mark = 0.10
        chunk = 1 << 20  # 1 MiB
        with open(dest, 'wb') as fh:
            while True:
                buf = resp.read(chunk)
                if not buf:
                    break
                fh.write(buf)
                done += len(buf)
                if total and done / total >= next_mark:
                    logger.warning('  ...%d%% (%.2f/%.2f GB)',
                                   int(done / total * 100), done / 1e9, total / 1e9)
                    next_mark += 0.10
    return dest


def download_reference(db_dir, gtex_mode='count'):
    """Download + extract the reference bundle into ``db_dir``; return the resolved
    reference root. Raises :class:`ReferenceError` if it still can't be resolved."""
    os.makedirs(db_dir, exist_ok=True)
    dest = os.path.join(db_dir, 'download.tar.gz')
    logger.warning('Downloading SNAF reference (~%.1f GB) from %s', REFERENCE_APPROX_GB, REFERENCE_URL)
    _download_stream(REFERENCE_URL, dest)
    logger.warning('Extracting %s ...', dest)
    with tarfile.open(dest, 'r:gz') as tf:
        tf.extractall(db_dir)  # nosec - trusted public mirror
    root = locate_reference_root(db_dir, gtex_mode)
    if root is None:
        raise ReferenceError(
            'Downloaded and extracted the reference into {d} but the required files '
            'are still not present. The download may be corrupt; delete {dest} and '
            'retry.'.format(d=db_dir, dest=dest))
    logger.warning('SNAF reference ready at %s', root)
    return root


def ensure_reference(db_dir, gtex_mode='count', download=False):
    """Resolve and return the reference root under ``db_dir``.

    If the reference is missing: download it when ``download`` is True, otherwise
    raise a clear :class:`ReferenceError` that names the missing files and prints
    the exact command to fetch them. Never returns silently on a missing bundle.
    """
    if db_dir is None:
        raise ReferenceError('No --db_dir provided; the full SNAF pipeline needs the '
                             'reference bundle.\nGet it with:\n\n    {}\n'.format(
                                 download_command('snaf_reference')))
    root = locate_reference_root(db_dir, gtex_mode)
    if root is not None:
        return root
    if download:
        return download_reference(db_dir, gtex_mode)
    raise ReferenceError(
        "SNAF reference not found under '{d}'.\n"
        "Expected '{d}/Alt91_db/' (with {alt}) and '{d}/controls/{ctrl}'.\n"
        "Fetch it (~{gb} GB) with:\n\n    {cmd}\n\n"
        "or re-run with --download_ref to fetch it automatically.".format(
            d=db_dir, alt=', '.join(_ALT91_REQUIRED), ctrl=control_filename(gtex_mode),
            gb=REFERENCE_APPROX_GB, cmd=download_command(db_dir)))
