"""Compare a reference transcriptome GTF/GFF against a BAM-derived collapsed isoform catalog.

Answers two questions with explicit denominators:
  1. Which reference transcripts have NO structure in our catalog?
  2. Which of our isoforms have NO structure in the reference?

and then, for the highest-expressed disagreement in each direction, says HOW the two structures
differ.

WHY BOTH FILES GO THROUGH ONE gff_process RUN
---------------------------------------------
The comparison key is the AltAnalyze exon-structure string (``E1.1|E2.1|I3.1_7654321|...``). Those
tokens are assigned by ``gff_process.consolidateLongReadGFFs`` against an Ensembl exon model. Known
exon boundaries get stable tokens, but a NOVEL boundary is named from what else that run has already
seen, so annotating two files SEPARATELY can give the same physical structure two different strings.
``isoform_collapse/reference.py`` records this for ENST00000262407 (ITGA2B). Both files are therefore
submitted to ONE ``consolidateLongReadGFFs`` call, as a list, exactly as the multi-GFF path is
already used elsewhere. ``verify_namespace_shared`` measures the size of that effect for the run at
hand, so the choice is evidenced rather than assumed.

The comparison is structure identity, not coordinate overlap: two transcripts match when they share
a gene AND an exon-structure string. That is the same key the collapse itself uses, so a match here
means the two pipelines built the same isoform.
"""

from __future__ import annotations

import os
import csv
import gzip
import shutil
import collections

from . import io_utils as _io

from . import gff_process
from . import isoform_collapse_utils as icu


# --------------------------------------------------------------------------- parsing

def read_transcript_associations(path):
    """Parse a gff_process transcript_associations.txt.

    Columns (gff_process.py:824): gene, strand, structure, transcript_id, source_label. The source
    label is the input file's basename truncated at the first '.g', so 'ENCFF801ZHP.gtf.gz' becomes
    'ENCFF801ZHP'.

    Returns (records, per_source_counts) where records is a list of dicts.
    """
    records = []
    per_source = collections.Counter()
    with _io.smart_open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            gene, strand, structure, transcript_id, source = parts[:5]
            if not structure:
                continue
            records.append({'gene': gene, 'strand': strand, 'structure': structure,
                            'transcript_id': transcript_id, 'source': source})
            per_source[source] += 1
    return records, per_source


def index_by_structure(records, source):
    """{(gene, structure): [transcript_id, ...]} for ONE source label.

    The key is the (gene, structure) PAIR. A bare structure string is not unique across genes, so
    keying on the structure alone would merge unrelated transcripts.
    """
    out = collections.defaultdict(list)
    for r in records:
        if r['source'] == source:
            out[(r['gene'], r['structure'])].append(r['transcript_id'])
    return out


# --------------------------------------------------------------------------- annotation

def annotate_together(gff_paths, exon_annot, work_dir):
    """Run ONE consolidateLongReadGFFs over every GFF so all structures share a token namespace.

    consolidateLongReadGFFs writes into ``os.getcwd()/gff-output`` when handed a LIST, so this
    chdirs into ``work_dir`` for the call. Without that it would overwrite the workflow's own
    gff-output/transcript_associations.txt, which the collapse catalog depends on.

    Returns the transcript_associations.txt path.
    """
    work_dir = os.path.abspath(str(work_dir))
    os.makedirs(work_dir, exist_ok=True)
    paths = [os.path.abspath(str(p)) for p in gff_paths]
    for p in paths:
        if not os.path.exists(p):
            raise FileNotFoundError(f"GFF to compare not found: {p}")
    ta = os.path.join(work_dir, 'gff-output', 'transcript_associations.txt')
    stamp = os.path.join(work_dir, 'gff-output', '.inputs')
    # realpath, not abspath: the reference is commonly reached through a symlink, and a lexical
    # path would make the cache key depend on which link was used to name the same file.
    signature = "\n".join(f"{os.path.realpath(p)}\t{os.path.getsize(p)}" for p in paths)
    # Reuse a previous joint annotation when it was built from EXACTLY these inputs. The run costs
    # minutes over millions of records, and re-deriving it changes nothing.
    if os.path.exists(ta) and os.path.getsize(ta) > 0 and os.path.exists(stamp):
        if open(stamp).read() == signature:
            log_n = sum(1 for _ in _io.smart_open(ta))
            print(f"[compare] reusing the cached joint annotation ({log_n:,} rows): {ta}")
            return ta
    previous = os.getcwd()
    try:
        os.chdir(work_dir)
        gff_process.consolidateLongReadGFFs(paths, str(exon_annot), mode='collapse')
    finally:
        os.chdir(previous)
    if not os.path.exists(ta):
        raise FileNotFoundError(f"consolidateLongReadGFFs produced no {ta}")
    with open(stamp, 'w') as handle:
        handle.write(signature)
    return ta


def annotate_alone(gff_path, exon_annot, work_dir):
    """Annotate ONE GFF on its own, for the namespace check. Same entry point, single-file mode."""
    work_dir = os.path.abspath(str(work_dir))
    os.makedirs(work_dir, exist_ok=True)
    local = os.path.join(work_dir, os.path.basename(str(gff_path)))
    if not os.path.exists(local):
        shutil.copy2(str(gff_path), local)
    gff_process.consolidateLongReadGFFs(local, str(exon_annot), mode='collapse')
    return os.path.join(work_dir, 'gff-output', 'transcript_associations.txt')


def source_label(gff_path):
    """The label consolidateLongReadGFFs writes in column 5: basename truncated at the first '.g'."""
    return os.path.basename(str(gff_path)).split('.g')[0]


def verify_namespace_shared(query_gff, exon_annot, combined_ta, work_dir, log=print):
    """Measure whether annotating the query GFF ALONE gives the same structures as annotating it
    WITH the sample GFFs.

    This is the evidence for running both files in one pass. It reports the fraction of the query's
    transcripts whose structure string is identical between the two runs. A fraction below 100%
    means a separate-run comparison would have produced false mismatches.
    """
    alone_dir = os.path.join(work_dir, 'namespace_check')
    alone_ta = annotate_alone(query_gff, exon_annot, alone_dir)
    label = source_label(query_gff)

    alone_records, _ = read_transcript_associations(alone_ta)
    combined_records, _ = read_transcript_associations(combined_ta)
    alone_by_tx = {r['transcript_id']: (r['gene'], r['structure'])
                   for r in alone_records if r['source'] == label}
    comb_by_tx = {r['transcript_id']: (r['gene'], r['structure'])
                  for r in combined_records if r['source'] == label}
    shared_tx = set(alone_by_tx) & set(comb_by_tx)
    same = sum(1 for t in shared_tx if alone_by_tx[t] == comb_by_tx[t])
    n = len(shared_tx)
    pct = (100.0 * same / n) if n else 0.0
    log(f"[namespace] {label}: {same:,}/{n:,} transcripts ({pct:.2f}%) keep an IDENTICAL "
        f"(gene, structure) when annotated alone vs alongside the sample GFF(s).")
    if n and same < n:
        log(f"[namespace] {n - same:,} transcripts ({100.0 - pct:.2f}%) would have been scored as "
            f"false mismatches by a separate-run comparison. One combined run is required.")
    return {'transcripts_compared': n, 'identical': same, 'pct_identical': pct,
            'alone_ta': alone_ta}


# --------------------------------------------------------------------------- expression

def load_catalog(path):
    """FINAL_isoform_catalog.tsv -> {final_isoform_id: row dict}.

    Columns: gene, final_isoform_id, exon_blocks, total_reads, bin, known. ``total_reads`` is the
    collapsed read count across every sample, so it is the expression measure for 'most highly
    expressed' when no per-sample matrix is supplied.
    """
    out = {}
    with _io.smart_open(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            row['total_reads'] = int(float(row.get('total_reads') or 0))
            row['exon_blocks'] = int(float(row.get('exon_blocks') or 0))
            out[row['final_isoform_id']] = row
    return out


def catalog_expression_by_transcript(catalog):
    """Map the transcript ids that appear in our collapsed GFF to their catalog read totals.

    A NOVEL final_isoform_id is ``<molecule>.<sample>`` while the GFF stamps the BARE ``<molecule>``
    as transcript_id (isoform_collapse/pipeline.py stage_protein). A KNOWN final_isoform_id is a
    version-stripped ENST while the reference GFF carries ``ENST.<version>``. Register BOTH spellings
    so either lookup resolves.
    """
    out = {}
    for final_id, row in catalog.items():
        reads = row['total_reads']
        out[final_id] = (final_id, reads)
        if '.' in final_id:
            bare = final_id.rsplit('.', 1)[0]
            out.setdefault(bare, (final_id, reads))
        else:
            out.setdefault(final_id, (final_id, reads))
    return out


def load_gene_symbols(path):
    """ENSG -> HGNC symbol from an Ensembl annotation table (column 0 -> column 1).

    The bundled ``Hs_Ensembl-annotations.txt`` is headerless: gene id, symbol, description, blank.
    Returns {} when the file is absent, so a missing table degrades the report rather than failing it.
    """
    out = {}
    if not path or not _io.exists(str(path)):
        return out
    with _io.smart_open(str(path)) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] and parts[1]:
                out.setdefault(parts[0], parts[1])
    return out


def load_reference_counts(path, id_col='annot_transcript_id', count_col=None):
    """Load a reference transcript quantification table (e.g. an ENCODE TALON TSV).

    Returns (counts, meta, count_col_used). ``count_col`` defaults to the LAST column, which is where
    the ENCODE TALON abundance files put the per-replicate count.
    """
    counts = {}
    meta = {}
    with _io.smart_open(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames or []
        if id_col not in fields:
            raise ValueError(f"{path} has no '{id_col}' column. Columns: {fields}")
        col = count_col or fields[-1]
        if col not in fields:
            raise ValueError(f"{path} has no '{col}' column. Columns: {fields}")
        for row in reader:
            tid = row[id_col]
            try:
                counts[tid] = float(row[col])
            except (TypeError, ValueError):
                counts[tid] = 0.0
            meta[tid] = {'novelty': row.get('transcript_novelty', ''),
                         'n_exons': int(float(row.get('n_exons') or 0)),
                         'gene_name': row.get('annot_gene_name', ''),
                         'transcript_name': row.get('annot_transcript_name', '')}
    return counts, meta, col


# --------------------------------------------------------------------------- structural diff

def structure_tokens(structure):
    """'E1.1|E2.1|I3.1' -> ['E1.1','E2.1','I3.1']."""
    return [t for t in str(structure).split('|') if t]


def describe_difference(structure, other):
    """Describe how ``structure`` differs from ``other`` in one line, at the exon-token level."""
    a = structure_tokens(structure)
    b = structure_tokens(other)
    sa, sb = set(a), set(b)
    only_a = [t for t in a if t not in sb]
    only_b = [t for t in b if t not in sa]
    try:
        contained = icu.is_contiguous_subsequence(a, b)
    except Exception:
        contained = False
    try:
        contains = icu.is_contiguous_subsequence(b, a)
    except Exception:
        contains = False
    if contained:
        relation = 'query is a contiguous SUBSTRING of the catalog isoform (shorter, same path)'
    elif contains:
        relation = 'query CONTAINS the catalog isoform (longer, same path)'
    elif not only_a and not only_b:
        relation = 'same exon set, DIFFERENT order'
    else:
        relation = 'different exon path'
    return {
        'relation': relation,
        'n_tokens_query': len(a),
        'n_tokens_match': len(b),
        'tokens_only_in_query': ','.join(only_a) if only_a else '',
        'tokens_only_in_match': ','.join(only_b) if only_b else '',
    }


def nearest_same_gene(structure, candidate_structures):
    """Pick the structure in the same gene sharing the most exon tokens. Returns (structure, jaccard)
    or (None, 0.0) when the gene has no candidate."""
    a = set(structure_tokens(structure))
    best, best_j = None, -1.0
    for cand in candidate_structures:
        b = set(structure_tokens(cand))
        union = a | b
        j = (len(a & b) / len(union)) if union else 0.0
        if j > best_j:
            best, best_j = cand, j
    return best, (best_j if best is not None else 0.0)


# --------------------------------------------------------------------------- driver

def _strip_version(tid):
    """ENST00000262407.6 -> ENST00000262407. Reference GFFs carry a version; the collapse catalog
    stores the version-stripped id."""
    return str(tid).split('.')[0]


def _stream_ta(path):
    """Yield (gene, strand, structure, transcript_id, source) from a transcript_associations.txt."""
    with _io.smart_open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5 and parts[2]:
                yield parts[0], parts[1], parts[2], parts[3], parts[4]


def _molecule_of(final_isoform_id, sample_labels):
    """A NOVEL final_isoform_id is '<molecule>.<library>'. Strip the library suffix to recover the
    molecule id the read GFF carries. A KNOWN id is an ENST and has no molecule."""
    for lab in sample_labels:
        suffix = '.' + lab
        if final_isoform_id.endswith(suffix):
            return final_isoform_id[:-len(suffix)]
    return None


def compare(query_gff, sample_gffs, exon_annot, catalog, final_ta, sample_tas,
            ref_gff=None, counts=None, reference_counts=None, out_dir='.', top_n=10,
            restrict_to_observed=True, gene_symbols=None, log=print):
    """Compare ``query_gff`` (a reference annotation) against the BAM-derived collapsed catalog.

    sample_gffs: the per-library READ GFFs (``<library>.gff.gz``). These are used, NOT the collapsed
      ``combined.gff.gz``, because combined.gff.gz interleaves two attribute syntaxes (novel records
      in GTF style, reference records copied verbatim from GENCODE GFF3) and gff_process fixes its
      attribute delimiter from the first record of a file, so that file cannot be re-parsed by the
      same reader.
    final_ta:   gff-output/transcript_associations.txt written by the collapse -- one row per FINAL
      isoform: gene, strand, representative structure, final_isoform_id, source.
    sample_tas: the per-library molecule tables (<library gff dir>/gff-output/transcript_associations.txt).
    ref_gff:    the reference GFF the collapse injected (e.g. gencode.v45.annotation.gff3). REQUIRED
      whenever the catalog contains KNOWN (ENST) isoforms. A known isoform absorbs TRUNCATED reads,
      so no molecule carries its full-length structure; its structure exists only in the reference
      annotation. Including the reference in the same joint run puts those structures in the shared
      namespace so they can be compared. Measured on the PC-3 catalog: 27,792 of 29,775 known
      isoforms (93.3%) have no molecule carrying their representative structure.

    HOW THE TWO SIDES ARE MADE COMPARABLE
    Structures are compared only within ONE joint gff_process run, so the exon-token namespace is
    shared by construction. The collapse catalog, however, was built in a SEPARATE (sample-only) run.
    The two are bridged through MOLECULE IDS, which are invariant across runs: each final isoform's
    representative structure identifies a molecule in the sample run, and that same molecule carries
    a structure in the joint run. No assumption about namespace stability is made; the drift is
    measured and reported.
    """
    out_dir = os.path.abspath(str(out_dir))
    os.makedirs(out_dir, exist_ok=True)
    work_dir = os.path.join(out_dir, 'annotation')

    if gene_symbols is None:
        gene_symbols = os.path.join(os.path.dirname(os.path.abspath(str(exon_annot))),
                                    'Hs_Ensembl-annotations.txt')
    symbols = load_gene_symbols(gene_symbols)
    log(f"[compare] gene symbols loaded: {len(symbols):,}")

    sample_gffs = [str(p) for p in sample_gffs]
    sample_tas = [str(p) for p in sample_tas]
    q_label = source_label(query_gff)
    s_labels = [source_label(p) for p in sample_gffs]
    log(f"[compare] query source '{q_label}', sample source(s) {s_labels}")

    # ---- 1. final isoforms -> representative structure (sample-run namespace)
    final_rep = {}            # final_isoform_id -> (gene, structure)
    for gene, _strand, struct, final_id, _src in _stream_ta(final_ta):
        final_rep[final_id] = (gene, struct)
    log(f"[compare] final isoforms: {len(final_rep):,} (from {final_ta})")

    wanted = set(final_rep.values())

    # ---- 2. representative structure -> a molecule that carries it (sample-run namespace)
    rep_molecule = {}         # (gene, structure) -> molecule id
    n_sample_rows = 0
    for ta in sample_tas:
        for gene, _strand, struct, mol, _src in _stream_ta(ta):
            n_sample_rows += 1
            key = (gene, struct)
            if key in wanted and key not in rep_molecule:
                rep_molecule[key] = mol
    log(f"[compare] molecule rows scanned: {n_sample_rows:,}; representative molecules resolved: "
        f"{len(rep_molecule):,}/{len(wanted):,} distinct representative structures")

    final_molecule = {}       # final_isoform_id -> molecule id  (NOVEL isoforms)
    final_enst = {}           # final_isoform_id -> ENST id       (KNOWN isoforms)
    for final_id, key in final_rep.items():
        mol = rep_molecule.get(key)
        if mol is not None:
            final_molecule[final_id] = mol
        elif final_id.startswith(('ENST', 'NM_', 'NR_', 'XM_', 'XR_')):
            # A known isoform absorbs truncated reads, so its own full-length structure need not be
            # carried by any molecule. Bridge it by transcript id through the reference GFF instead.
            final_enst[final_id] = _strip_version(final_id)
    unresolved = len(final_rep) - len(final_molecule) - len(final_enst)
    log(f"[compare] bridge: {len(final_molecule):,} isoforms via a representative molecule, "
        f"{len(final_enst):,} known isoforms via their reference transcript id, "
        f"{unresolved:,} unresolved")
    if final_enst and not ref_gff:
        raise ValueError(
            f"{len(final_enst):,} catalog isoforms are KNOWN (ENST) and have no molecule carrying "
            f"their full-length structure, so they can only be placed in the joint namespace "
            f"through the reference annotation. Pass --ref_gff (the same reference the collapse "
            f"used, e.g. gencode.v45.annotation.gff3).")
    needed_molecules = set(final_molecule.values())
    needed_enst = set(final_enst.values())

    # ---- 3. ONE joint annotation run: read GFF(s) + the reference + the query annotation.
    joint_inputs = list(sample_gffs)
    r_label = None
    if ref_gff:
        joint_inputs.append(str(ref_gff))
        r_label = source_label(ref_gff)
        log(f"[compare] reference source '{r_label}' included in the joint run")
    joint_inputs.append(str(query_gff))
    joint_ta = annotate_together(joint_inputs, exon_annot, work_dir)

    # ---- 4. stream the joint table
    mol_joint = {}                       # molecule -> (gene, structure)  [only the ones we need]
    enst_joint = {}                      # ENST (version stripped) -> (gene, structure)
    query_idx = collections.defaultdict(list)   # (gene, structure) -> [query transcript ids]
    drift_same = drift_total = 0
    sample_struct_of_mol = {}
    for ta in sample_tas:                # molecule -> sample-run structure, for the drift measure
        for gene, _s, struct, mol, _src in _stream_ta(ta):
            if mol in needed_molecules:
                sample_struct_of_mol[mol] = (gene, struct)
    for gene, _strand, struct, tid, src in _stream_ta(joint_ta):
        if src == q_label:
            query_idx[(gene, struct)].append(tid)
        elif r_label and src == r_label:
            bare = _strip_version(tid)
            if bare in needed_enst and bare not in enst_joint:
                enst_joint[bare] = (gene, struct)
        elif src in s_labels:
            if tid in needed_molecules:
                mol_joint[tid] = (gene, struct)
                prior = sample_struct_of_mol.get(tid)
                if prior is not None:
                    drift_total += 1
                    drift_same += (prior == (gene, struct))
    pct_stable = (100.0 * drift_same / drift_total) if drift_total else 0.0
    log(f"[namespace] {drift_same:,}/{drift_total:,} representative molecules ({pct_stable:.2f}%) "
        f"keep an IDENTICAL (gene, structure) between the sample-only run and the joint run.")
    if drift_total and drift_same < drift_total:
        log(f"[namespace] {drift_total - drift_same:,} molecules ({100.0 - pct_stable:.2f}%) shift "
            f"token namespace. The comparison uses the JOINT structures for both sides, so this "
            f"shift does not create false mismatches; it is reported as evidence that a "
            f"separate-run comparison would have.")

    # ---- 5. our catalog in the joint namespace
    ours_idx = collections.defaultdict(list)    # (gene, structure) -> [final_isoform_id]
    placed_mol = placed_enst = 0
    for final_id, mol in final_molecule.items():
        key = mol_joint.get(mol)
        if key is not None:
            ours_idx[key].append(final_id)
            placed_mol += 1
    for final_id, enst in final_enst.items():
        key = enst_joint.get(enst)
        if key is not None:
            ours_idx[key].append(final_id)
            placed_enst += 1
    placed = placed_mol + placed_enst
    log(f"[compare] catalog isoforms placed in the joint namespace: {placed:,}/{len(final_rep):,} "
        f"({placed_mol:,} via molecules, {placed_enst:,} via reference transcript ids)")
    if placed < len(final_rep):
        log(f"[compare] {len(final_rep) - placed:,} catalog isoforms could not be placed and are "
            f"EXCLUDED from the 'ours' denominator. Report this fraction.")

    # ---- 6. expression
    cat = load_catalog(catalog)
    reads_of_final = {fid: row['total_reads'] for fid, row in cat.items()}
    ref_counts, ref_meta, ref_col = ({}, {}, None)
    if reference_counts:
        ref_counts, ref_meta, ref_col = load_reference_counts(reference_counts)
        log(f"[compare] reference quantification: {len(ref_counts):,} transcripts, column "
            f"'{ref_col}', total {int(sum(ref_counts.values())):,} reads")

    def query_observed(tx_ids):
        if not (restrict_to_observed and ref_counts):
            return True
        return any(ref_counts.get(t, 0) > 0 for t in tx_ids)

    q_keys_all = set(query_idx)
    q_keys = {k for k in q_keys_all if query_observed(query_idx[k])}
    o_keys = set(ours_idx)
    shared = q_keys & o_keys
    q_only = q_keys - o_keys
    o_only = o_keys - q_keys

    def our_reads(key):
        return sum(reads_of_final.get(f, 0) for f in ours_idx.get(key, []))

    def query_reads(key):
        return sum(ref_counts.get(t, 0.0) for t in query_idx.get(key, []))

    def query_novelty(key):
        vals = [ref_meta.get(t, {}).get('novelty', '') for t in query_idx.get(key, [])]
        vals = [v for v in vals if v]
        return ','.join(sorted(set(vals)))

    summary = {
        'query_structures_in_file': len(q_keys_all),
        'query_structures_observed': len(q_keys),
        'our_structures': len(o_keys),
        'our_catalog_isoforms': len(final_rep),
        'our_isoforms_placed': placed,
        'shared': len(shared),
        'query_only': len(q_only),
        'our_only': len(o_only),
        'pct_query_found_in_ours': (100.0 * len(shared) / len(q_keys)) if q_keys else 0.0,
        'pct_ours_found_in_query': (100.0 * len(shared) / len(o_keys)) if o_keys else 0.0,
        'namespace_pct_identical': round(pct_stable, 2),
        'namespace_molecules_compared': drift_total,
    }

    q_by_gene = collections.defaultdict(list)
    for g, st in q_keys:
        q_by_gene[g].append(st)
    o_by_gene = collections.defaultdict(list)
    for g, st in o_keys:
        o_by_gene[g].append(st)

    diff_header = ['nearest_structure_other_side', 'jaccard_exon_tokens', 'relation',
                   'n_exon_tokens_nearest', 'tokens_only_here', 'tokens_only_in_nearest']

    def write_keys(keys, path, reads_fn, tx_idx, other_by_gene, with_diff, extra_cols=None):
        rows = sorted(keys, key=lambda k: -reads_fn(k))
        header = ['gene', 'symbol', 'reads', 'n_exon_tokens', 'ids', 'structure']
        header += (extra_cols or [])
        header += (diff_header if with_diff else [])
        with open(path, 'w') as handle:
            handle.write('\t'.join(header) + '\n')
            for (g, st) in rows:
                row = [g, symbols.get(g, ''), f"{reads_fn((g, st)):.0f}",
                       str(len(structure_tokens(st))), ','.join(tx_idx.get((g, st), [])), st]
                if extra_cols:
                    row.append(query_novelty((g, st)))
                if with_diff:
                    near, jac = nearest_same_gene(st, other_by_gene.get(g, []))
                    if near is None:
                        row += ['', '0.000', 'no isoform of this gene on the other side', '0', '', '']
                    else:
                        d = describe_difference(st, near)
                        row += [near, f"{jac:.3f}", d['relation'], str(d['n_tokens_match']),
                                d['tokens_only_in_query'], d['tokens_only_in_match']]
                handle.write('\t'.join(row) + '\n')
        return len(rows)

    p_shared = os.path.join(out_dir, 'shared_structures.txt')
    p_qonly = os.path.join(out_dir, 'reference_only_structures.txt')
    p_oonly = os.path.join(out_dir, 'bam_only_structures.txt')

    # The shared table carries BOTH id sets. Ours alone cannot be joined back to the reference
    # metadata (novelty class, exon count), which is exactly what a reader needs to interpret it.
    with open(p_shared, 'w') as handle:
        handle.write('\t'.join(['gene', 'symbol', 'our_reads', 'reference_reads', 'n_exon_tokens',
                                'our_ids', 'reference_ids', 'reference_novelty', 'structure']) + '\n')
        for (g, st) in sorted(shared, key=lambda k: -our_reads(k)):
            handle.write('\t'.join([
                g, symbols.get(g, ''), f"{our_reads((g, st)):.0f}", f"{query_reads((g, st)):.0f}",
                str(len(structure_tokens(st))), ','.join(ours_idx.get((g, st), [])),
                ','.join(query_idx.get((g, st), [])), query_novelty((g, st)), st]) + '\n')
    write_keys(q_only, p_qonly, query_reads, query_idx, o_by_gene, True, ['novelty'])
    write_keys(o_only, p_oonly, our_reads, ours_idx, q_by_gene, True)

    def top_table(keys, path, reads_fn, tx_idx, other_by_gene, title, with_novelty=False):
        rows = sorted(keys, key=lambda k: -reads_fn(k))[:top_n]
        header = ['rank', 'gene', 'symbol', 'reads', 'ids', 'n_exon_tokens', 'structure']
        if with_novelty:
            header.append('novelty')
        header += ['nearest_structure_other_side', 'jaccard_exon_tokens', 'relation',
                   'n_exon_tokens_nearest', 'tokens_only_here', 'tokens_only_in_nearest']
        with open(path, 'w') as handle:
            handle.write(f"# {title}\n")
            handle.write('\t'.join(header) + '\n')
            for i, (g, st) in enumerate(rows, 1):
                near, jac = nearest_same_gene(st, other_by_gene.get(g, []))
                if near is None:
                    d = {'relation': 'no isoform of this gene on the other side',
                         'n_tokens_match': 0, 'tokens_only_in_query': '', 'tokens_only_in_match': ''}
                    near = ''
                else:
                    d = describe_difference(st, near)
                row = [str(i), g, symbols.get(g, ''), f"{reads_fn((g, st)):.0f}",
                       ','.join(tx_idx.get((g, st), [])), str(len(structure_tokens(st))), st]
                if with_novelty:
                    row.append(query_novelty((g, st)))
                row += [near, f"{jac:.3f}", d['relation'], str(d['n_tokens_match']),
                        d['tokens_only_in_query'], d['tokens_only_in_match']]
                handle.write('\t'.join(row) + '\n')
        return rows

    p_top_q = os.path.join(out_dir, f'top{top_n}_reference_only.txt')
    p_top_o = os.path.join(out_dir, f'top{top_n}_bam_only.txt')
    top_table(q_only, p_top_q, query_reads, query_idx, o_by_gene,
              f'Top {top_n} reference structures by reference read count with NO match in our catalog',
              with_novelty=True)
    top_table(o_only, p_top_o, our_reads, ours_idx, q_by_gene,
              f'Top {top_n} catalog isoforms by our read count with NO match in the reference')

    p_sum = os.path.join(out_dir, 'structure_concordance_summary.txt')
    with open(p_sum, 'w') as handle:
        handle.write("metric\tvalue\n")
        for k, v in summary.items():
            handle.write(f"{k}\t{v}\n")
        handle.write(f"query_source_label\t{q_label}\n")
        handle.write(f"catalog_source_labels\t{','.join(s_labels)}\n")
        handle.write(f"reference_count_column\t{ref_col or ''}\n")
        handle.write(f"restricted_to_observed\t{bool(restrict_to_observed and ref_counts)}\n")

    log("")
    log("STRUCTURE CONCORDANCE (key = gene + AltAnalyze exon-structure string, one joint run)")
    log(f"  reference structures in the GTF            : {summary['query_structures_in_file']:,}")
    log(f"  reference structures with >=1 read         : {summary['query_structures_observed']:,}")
    log(f"  our catalog structures                     : {summary['our_structures']:,}")
    log(f"  shared                                     : {summary['shared']:,}")
    log(f"  in reference, MISSING from ours            : {summary['query_only']:,} "
        f"({100.0 - summary['pct_query_found_in_ours']:.1f}% of {summary['query_structures_observed']:,})")
    log(f"  in ours, MISSING from reference            : {summary['our_only']:,} "
        f"({100.0 - summary['pct_ours_found_in_query']:.1f}% of {summary['our_structures']:,})")
    log("")
    for p in (p_sum, p_shared, p_qonly, p_oonly, p_top_q, p_top_o):
        log(f"  wrote {p}")
    return summary
