"""Custom cell-surface (surfaceome) gene database for SNAF-B.

SNAF-B decides which splicing events can encode a surface antigen with two gates
that both live inside ``Alt91_db``:

* ``uniprot_isoform_enhance.fasta`` -- ENSG -> {accession: protein sequence}. It is
  BOTH the gene whitelist (:func:`surface.main.filter_to_membrane_protein`) and the
  source of the reference protein used by the novelty alignment
  (``surface/alignment.py``) and by the length / extracellular comparison
  (``surface/main.py::_eval_sa_for_generate``). A gene with no record here cannot be
  scored -- the lookup raises ``KeyError``.
* ``human_membrane_proteins_acc2ens.txt`` -- the second membership test in
  :meth:`surface.main.SurfaceAntigen.is_membrane_protein`.

That built-in pair covers 2,810 ENSGs pulled from a UniProt subcellular-location
keyword search. A curated surfaceome such as SURFY (Bausch-Fluck et al., PNAS 2018)
disagrees with it on both sides, so this module lets a user supply their own list.

A *surface database* here is a directory holding three files:

    surface_genes.txt      TSV: ensg, symbol, has_reference, reference_source,
                           n_reference_isoforms  -- EVERY gene of the input list,
                           including the ones with no reference sequence (marked
                           has_reference=0), so nothing is hidden.
    surface_reference.fasta  the reference isoform sequences for the genes that have
                           one, in the exact header grammar SNAF already parses:
                           ``sp|ACC|ENTRY_NAME description|ENSG``. The canonical
                           accession is written FIRST per gene because the
                           downstream comparison takes the first record as THE
                           reference protein.
    surface_db_params.json  provenance: inputs, counts, mode, build time.

:func:`build_surface_db` writes one; :func:`load_surface_db` reads one. A bare
one- or two-column gene table is also accepted by :func:`load_surface_db`, in which
case the built-in reference sequences are used and genes without one are dropped
(and counted).
"""
import os
import io
import re
import csv
import json
import gzip
import hashlib
import logging
from collections import OrderedDict
from datetime import datetime

logger = logging.getLogger(__name__)

GENES_FILE = 'surface_genes.txt'
FASTA_FILE = 'surface_reference.fasta'
PARAMS_FILE = 'surface_db_params.json'

GENES_HEADER = ('ensg', 'symbol', 'has_reference', 'reference_source', 'n_reference_isoforms')

_ENSG_RE = re.compile(r'^ENSG\d+')


def _open_text(path):
    """Open a plain or gzipped text file, tolerating CRLF and a UTF-8 BOM."""
    if str(path).endswith('.gz'):
        return io.TextIOWrapper(gzip.open(path, 'rb'), encoding='utf-8-sig', newline='')
    return open(path, 'r', encoding='utf-8-sig', newline='')


def _rows(path):
    with _open_text(path) as f:
        for line in f:
            yield line.rstrip('\r\n')


# --------------------------------------------------------------------------- input parsing
def parse_gene_table(path):
    """Parse a user gene list into ``[(ensg, symbol), ...]`` plus a stats dict.

    Accepts any TSV/CSV with one column holding Ensembl gene IDs. The ENSG column is
    found by content (first column whose values match ``ENSG\\d+``), not by name, so
    column order does not matter. A second text column, when present, is taken as the
    gene symbol. Rows with no Ensembl ID are counted and reported, never silently
    dropped. Version suffixes (``ENSG00000123456.7``) are stripped.
    """
    lines = [ln for ln in _rows(path) if ln.strip() != '']
    if not lines:
        raise ValueError('surface gene table is empty: {}'.format(path))
    delim = '\t' if lines[0].count('\t') >= lines[0].count(',') else ','
    table = list(csv.reader(lines, delimiter=delim))

    # Which column holds the ENSG IDs? Pick the column with the most ENSG-looking values.
    ncol = max(len(r) for r in table)
    hits = [0] * ncol
    for r in table:
        for j, v in enumerate(r):
            if _ENSG_RE.match(v.strip()):
                hits[j] += 1
    if max(hits) == 0:
        raise ValueError(
            'no column of {} contains Ensembl gene IDs (ENSG...); a surface gene table '
            'must carry ENSG identifiers because SNAF-B keys every reference on ENSG'.format(path))
    ecol = hits.index(max(hits))
    scol = 0 if ecol != 0 else (1 if ncol > 1 else None)

    # A first row that carries no Ensembl ID is the header; drop it before counting.
    has_header = not any(_ENSG_RE.match(v.strip()) for v in table[0])
    body = table[1:] if has_header else table

    seen = OrderedDict()
    n_rows = 0
    n_no_ensg = 0
    n_dup = 0
    for r in body:
        n_rows += 1
        ensg = r[ecol].strip().split('.')[0] if len(r) > ecol else ''
        if not _ENSG_RE.match(ensg):
            n_no_ensg += 1
            continue
        sym = r[scol].strip() if (scol is not None and len(r) > scol) else ''
        if ensg in seen:
            n_dup += 1
            if not seen[ensg] and sym:
                seen[ensg] = sym
            continue
        seen[ensg] = sym
    stats = {'input_rows': n_rows, 'rows_without_ensembl_id': n_no_ensg,
             'duplicate_ensg_rows': n_dup, 'unique_ensg': len(seen),
             'has_header': has_header,
             'ensg_column_index': ecol, 'symbol_column_index': scol}
    return [(g, s) for g, s in seen.items()], stats


def read_reference_fasta_records(path):
    """``{ensg: [(acc, title, seq), ...]}`` from a ``sp|ACC|NAME desc|ENSG`` FASTA,
    keeping each record's ORIGINAL header title so a rebuild can re-emit it verbatim.

    Order is preserved: the FIRST record of a gene is the one SNAF-B uses as THE
    reference protein, so callers must write the canonical accession first.
    """
    out = OrderedDict()
    state = {'acc': None, 'ensg': None, 'title': None}
    chunks = []

    def flush():
        if state['ensg'] is None:
            return
        seq = ''.join(chunks)
        if not seq:
            return
        out.setdefault(state['ensg'], []).append((state['acc'], state['title'], seq))

    with _open_text(path) as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                flush()
                title = line[1:]
                parts = title.split('|')
                if len(parts) < 4:
                    raise ValueError(
                        "FASTA header must be 'db|ACC|ENTRY_NAME description|ENSG' "
                        "(SNAF reads field 1 as the accession and field 3 as the gene); "
                        "got: {}".format(title))
                state['acc'] = parts[1].strip()
                state['ensg'] = parts[3].strip().split()[0]
                state['title'] = title
                chunks = []
            else:
                chunks.append(line.strip())
        flush()
    return out


def read_reference_fasta(path):
    """``{ensg: OrderedDict(acc -> seq)}`` -- the shape SNAF-B's ``dict_uni_fa`` uses."""
    out = OrderedDict()
    for ensg, recs in read_reference_fasta_records(path).items():
        d = OrderedDict()
        for acc, _title, seq in recs:
            d[acc] = seq
        out[ensg] = d
    return out


# --------------------------------------------------------------------------- build
def _load_ens2uniprot(path):
    """``{ensg: [accession, ...]}`` from AltAnalyze's ``Hs_Ensembl-UniProt.txt``."""
    out = {}
    for i, line in enumerate(_rows(path)):
        if i == 0 or not line.strip():
            continue
        p = line.split('\t')
        if len(p) < 2:
            continue
        out.setdefault(p[0].strip(), []).append(p[1].strip())
    return out


def _load_uniprot_sequences(path):
    """``{accession: (entry_name, sequence)}`` from AltAnalyze's ``uniprot_sequence.txt``
    (``ENTRY_NAME <tab> acc1,acc2,... <tab> sequence``). This file carries CANONICAL
    sequences only -- it holds no ``ACC-2`` isoform records -- so a gene built from it
    gets exactly one reference isoform."""
    out = {}
    for line in _rows(path):
        p = line.split('\t')
        if len(p) < 3:
            continue
        name, accs, seq = p[0].strip(), p[1], p[2].strip()
        if not seq:
            continue
        for a in accs.split(','):
            a = a.strip()
            if a and a not in out:
                out[a] = (name, seq)
    return out


def _acc_sort_key(acc):
    """Canonical accession first, then isoforms by ascending index."""
    if '-' in acc:
        base, _, idx = acc.partition('-')
        try:
            return (base, 1, int(idx))
        except ValueError:
            return (base, 1, 10 ** 6)
    return (acc, 0, 0)


def build_surface_db(gene_table, db_dir, outdir, uniprot_dir=None, mode='replace',
                     name=None):
    """Format a user surfaceome list into a SNAF-B surface database.

    :param gene_table: the user's list (e.g. ``SURFY-Ensembl.txt``). Any TSV/CSV with an
        Ensembl-gene-ID column.
    :param db_dir: the SNAF reference root holding ``Alt91_db/``. Supplies the built-in
        reference sequences (``uniprot_isoform_enhance.fasta``) and the Ensembl-91 gene
        models used to report how many genes SNAF can actually model.
    :param outdir: directory to write the database into.
    :param uniprot_dir: optional AltAnalyze UniProt directory (e.g.
        ``.../EnsMart100/uniprot/Hs``) holding ``Hs_Ensembl-UniProt.txt`` and
        ``uniprot_sequence.txt``. When given, genes in the list that have no built-in
        reference sequence get one built from these tables. Without it, those genes are
        reported and excluded -- SNAF-B cannot score a gene with no reference protein.
    :param mode: ``replace`` (default) -- the database is exactly the user's list;
        ``union`` -- the user's list plus the built-in genes.
    :return: ``(outdir, params_dict)``. ``params_dict`` carries every count with its
        denominator and is also written to ``surface_db_params.json``.
    """
    mode = str(mode).lower()
    if mode not in ('replace', 'union'):
        raise ValueError("mode must be 'replace' or 'union', got {!r}".format(mode))
    alt = os.path.join(db_dir, 'Alt91_db')
    builtin_fa_path = os.path.join(alt, 'uniprot_isoform_enhance.fasta')
    exon_table = os.path.join(alt, 'Hs_Ensembl_exon_add_col.txt')
    for p in (gene_table, builtin_fa_path, exon_table):
        if not os.path.exists(p):
            raise FileNotFoundError('required input not found: {}'.format(p))
    os.makedirs(outdir, exist_ok=True)

    genes, parse_stats = parse_gene_table(gene_table)
    requested = OrderedDict(genes)
    logger.info('surface DB: %d unique ENSG from %s (%d rows, %d without an Ensembl ID)',
                parse_stats['unique_ensg'], gene_table, parse_stats['input_rows'],
                parse_stats['rows_without_ensembl_id'])

    builtin_fa = read_reference_fasta_records(builtin_fa_path)

    modelled = set()
    for i, line in enumerate(_rows(exon_table)):
        if i == 0:
            continue
        g = line.split('\t', 1)[0]
        if g:
            modelled.add(g)

    if mode == 'union':
        for g in builtin_fa:
            requested.setdefault(g, '')

    # ---- assemble the reference sequences -----------------------------------------
    ens2uni = seq_lookup = None
    if uniprot_dir:
        p_map = os.path.join(uniprot_dir, 'Hs_Ensembl-UniProt.txt')
        p_seq = os.path.join(uniprot_dir, 'uniprot_sequence.txt')
        for p in (p_map, p_seq):
            if not os.path.exists(p):
                raise FileNotFoundError('--uniprot_dir given but {} is missing'.format(p))
        ens2uni = _load_ens2uniprot(p_map)
        seq_lookup = _load_uniprot_sequences(p_seq)
        logger.info('surface DB: UniProt tables -> %d ENSG mapped, %d accessions with a sequence',
                    len(ens2uni), len(seq_lookup))

    records = OrderedDict()   # ensg -> [(acc, header_title, seq), ...]
    source = {}
    n_builtin = n_built = 0
    no_reference = []
    for ensg, symbol in requested.items():
        if ensg in builtin_fa:
            # re-emit the built-in records with their ORIGINAL headers, so a gene that
            # already worked keeps byte-identical reference sequences and accessions.
            records[ensg] = list(builtin_fa[ensg])
            source[ensg] = 'Alt91_db/uniprot_isoform_enhance.fasta'
            n_builtin += 1
            continue
        if ens2uni is not None:
            built = []
            used = set()
            for acc in sorted(set(ens2uni.get(ensg, ())), key=_acc_sort_key):
                base = acc.split('-')[0]
                hit = seq_lookup.get(acc) or seq_lookup.get(base)
                if hit is None:
                    continue
                # uniprot_sequence.txt is canonical-only, so label the record with the
                # accession the sequence actually belongs to, never the isoform id.
                real_acc = acc if acc in seq_lookup else base
                if real_acc in used:
                    continue
                used.add(real_acc)
                entry_name, seq = hit
                desc = 'reference isoform OS=Homo sapiens OX=9606'
                if symbol:
                    desc += ' GN={}'.format(symbol)
                title = 'sp|{}|{} {}|{}'.format(real_acc, entry_name or 'UNKNOWN_HUMAN', desc, ensg)
                built.append((real_acc, title, seq))
            if built:
                records[ensg] = built
                source[ensg] = 'uniprot_dir'
                n_built += 1
                continue
        no_reference.append(ensg)

    # ---- write surface_reference.fasta --------------------------------------------
    fasta_path = os.path.join(outdir, FASTA_FILE)
    n_seq = 0
    with open(fasta_path, 'w') as f:
        for ensg, recs in records.items():
            for _acc, title, seq in recs:
                f.write('>' + title + '\n')
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i + 60] + '\n')
                n_seq += 1

    # ---- write surface_genes.txt ---------------------------------------------------
    genes_path = os.path.join(outdir, GENES_FILE)
    with open(genes_path, 'w') as f:
        f.write('\t'.join(GENES_HEADER) + '\n')
        for ensg in requested:
            recs = records.get(ensg)
            f.write('{}\t{}\t{}\t{}\t{}\n'.format(
                ensg, requested[ensg], 1 if recs else 0,
                source.get(ensg, 'none'), len(recs) if recs else 0))

    params = {
        'built': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'mode': mode,
        'name': name or os.path.basename(os.path.normpath(outdir)),
        'inputs': {
            'gene_table': os.path.abspath(gene_table),
            'db_dir': os.path.abspath(db_dir),
            'builtin_reference_fasta': os.path.abspath(builtin_fa_path),
            'uniprot_dir': os.path.abspath(uniprot_dir) if uniprot_dir else None,
        },
        'gene_table_parse': parse_stats,
        'counts': {
            'genes_requested': len(requested),
            'genes_with_reference': len(records),
            'reference_from_builtin': n_builtin,
            'reference_built_from_uniprot': n_built,
            'genes_without_reference_excluded': len(no_reference),
            'reference_sequences_written': n_seq,
            'builtin_whitelist_size': len(builtin_fa),
            'overlap_with_builtin': len(set(requested) & set(builtin_fa)),
            'in_user_list_not_builtin': len(set(requested) - set(builtin_fa)),
            'in_builtin_not_user_list': len(set(builtin_fa) - set(requested)),
            'genes_with_ensembl91_model': len(set(records) & modelled),
            'genes_without_ensembl91_model': len(set(records) - modelled),
        },
        'genes_without_reference': sorted(no_reference),
        'outputs': {
            'genes': os.path.abspath(genes_path),
            'fasta': os.path.abspath(fasta_path),
        },
    }
    with open(os.path.join(outdir, PARAMS_FILE), 'w') as f:
        json.dump(params, f, indent=2)
    return outdir, params


# --------------------------------------------------------------------------- load
def load_surface_db(path, builtin_fasta=None):
    """Load a surface database for :func:`surface.main.initialize`.

    :param path: a directory written by :func:`build_surface_db`, or a bare gene table
        (TSV/CSV with an ENSG column). For a bare table the reference sequences come
        from ``builtin_fasta`` and genes without one are excluded and counted.
    :param builtin_fasta: path to ``Alt91_db/uniprot_isoform_enhance.fasta``; required
        when ``path`` is a bare gene table.
    :return: dict with ``genes`` (the usable ENSG set), ``fasta``
        (``{ensg: {acc: seq}}``), ``symbols``, ``signature`` (stable hash of the gene
        set, used to key cached indexes) and ``stats``.
    """
    symbols = {}
    if os.path.isdir(path):
        genes_path = os.path.join(path, GENES_FILE)
        fasta_path = os.path.join(path, FASTA_FILE)
        if not (os.path.exists(genes_path) and os.path.exists(fasta_path)):
            raise FileNotFoundError(
                "'{}' is a directory but is not a surface database: expected {} and {}. "
                "Build one with `altanalyze3 snaf-build-surface-db`.".format(
                    path, GENES_FILE, FASTA_FILE))
        fasta = read_reference_fasta(fasta_path)
        listed = 0
        for i, line in enumerate(_rows(genes_path)):
            if i == 0 or not line.strip():
                continue
            p = line.split('\t')
            listed += 1
            if len(p) > 1 and p[1]:
                symbols[p[0]] = p[1]
        genes = set(fasta)
        stats = {'source': 'surface_db_dir', 'genes_listed': listed,
                 'genes_with_reference': len(genes),
                 'reference_sequences': sum(len(v) for v in fasta.values())}
    else:
        if builtin_fasta is None:
            raise ValueError('a bare surface gene table needs builtin_fasta to supply '
                             'the reference protein sequences')
        pairs, parse_stats = parse_gene_table(path)
        symbols = {g: s for g, s in pairs if s}
        builtin = read_reference_fasta(builtin_fasta)
        requested = [g for g, _ in pairs]
        fasta = OrderedDict((g, builtin[g]) for g in requested if g in builtin)
        genes = set(fasta)
        stats = {'source': 'gene_table', 'genes_listed': len(requested),
                 'genes_with_reference': len(genes),
                 'genes_without_reference_excluded': len(requested) - len(genes),
                 'reference_sequences': sum(len(v) for v in fasta.values())}
        stats.update(parse_stats)

    if not genes:
        raise ValueError('surface database {} yields 0 usable genes'.format(path))
    sig = hashlib.sha1(('\n'.join(sorted(genes))).encode()).hexdigest()[:12]
    return {'genes': genes, 'fasta': fasta, 'symbols': symbols,
            'signature': sig, 'stats': stats, 'path': os.path.abspath(path)}
