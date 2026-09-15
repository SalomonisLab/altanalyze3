"""Compare two isoform annotations on SPLICE-JUNCTION COMPOSITION alone.

This is a stricter and more portable comparison than matching AltAnalyze exon-structure strings.
A junction is a pair of genomic coordinates, so the comparison is independent of:

  * the Ensembl release used to name exons (an exon called ``E3.13`` in one build and ``E3.1`` in
    another is the same coordinate either way);
  * the transcript identifier (two pipelines may name the same junction chain differently);
  * the transcript's first and last exon boundaries, which are TSS/TES calls rather than splicing.

Two further rules make the comparison answer the biological question rather than a bookkeeping one:

  1. Only transcripts assigned to a KNOWN Ensembl gene are considered. Neither pipeline can be
     credited or blamed for a locus one of them does not call a gene.
  2. A junction chain that is a CONTIGUOUS SUBSTRING of a chain on the other side is REDUNDANT, not
     missing. It is the same splice path observed over a shorter span, which is what a 5'/3'
     truncated read produces and what the collapse deliberately merges.

Mono-exonic transcripts have no junction and are excluded from both sides.
"""

from __future__ import annotations

import os
import re
import gzip
import collections

from . import io_utils as _io


# --------------------------------------------------------------------------- parsing

_ATTR_GTF = re.compile(r'(\w+)\s+"([^"]*)"')
_ATTR_GFF3 = re.compile(r'(\w+)=([^;]*)')


def parse_attributes(field):
    """Parse a GFF/GTF attribute column in EITHER convention.

    A file may mix them: ``isoform_collapse.pipeline.stage_protein`` writes novel records in GTF
    style (``gene_id "X";``) and copies reference records verbatim from a GENCODE GFF3
    (``gene_id=X;``). Both are read here so such a file parses correctly.
    """
    out = dict(_ATTR_GTF.findall(field))
    if 'transcript_id' not in out or 'gene_id' not in out:
        for key, value in _ATTR_GFF3.findall(field):
            out.setdefault(key, value)
    return out


def strip_version(value):
    return str(value).split('.')[0]


def parse_gff_junctions(path, log=print):
    """Read a GFF/GTF and return {transcript_id: (chrom, strand, gene, junctions)}.

    ``junctions`` is a tuple of (donor, acceptor) genomic coordinate pairs, built from the exon
    blocks sorted by genomic START. Building them in genomic order makes the chain independent of
    whether the file lists exons transcriptomically or genomically, so a '+' and a '-' strand
    transcript are handled the same way.

    Exon records are grouped BY TRANSCRIPT ID rather than by file order, so a file that lists the
    transcript line before its exons and a file that lists it after both parse correctly.
    """
    opener = gzip.open if str(path).endswith('.gz') else open
    exons = collections.defaultdict(list)
    info = {}
    n_lines = 0
    with opener(str(path), 'rt') as handle:
        for line in handle:
            if not line or line[0] == '#':
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 9 or f[2] != 'exon':
                continue
            n_lines += 1
            attrs = parse_attributes(f[8])
            tid = attrs.get('transcript_id')
            if not tid:
                continue
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                continue
            exons[tid].append((start, end))
            if tid not in info:
                info[tid] = (f[0], f[6], strip_version(attrs.get('gene_id', '')))

    out = {}
    for tid, blocks in exons.items():
        blocks.sort()
        junctions = tuple((blocks[i][1], blocks[i + 1][0]) for i in range(len(blocks) - 1))
        chrom, strand, gene = info[tid]
        out[tid] = (chrom, strand, gene, junctions)
    log(f"[junctions] {os.path.basename(str(path))}: {n_lines:,} exon records -> {len(out):,} "
        f"transcripts")
    return out


# --------------------------------------------------------------------------- filters

def has_known_splice_site(chrom, strand, junctions, exon_coordinates):
    """True when at least ONE junction has a donor or acceptor that Ensembl already annotates.

    This mirrors ``bam/isoform_structure_extract.has_known_splice_site``, the gate the BAM
    extractor applies, so the same requirement is imposed on the annotation being compared. One
    known site is enough; the extractor does not require all of them.
    """
    for donor, acceptor in junctions:
        if ((chrom, donor, strand, 1) in exon_coordinates or
                (chrom, donor, strand, 2) in exon_coordinates or
                (chrom, acceptor, strand, 1) in exon_coordinates or
                (chrom, acceptor, strand, 2) in exon_coordinates):
            return True
    return False


def load_exon_coordinates(exon_annot):
    from . import gff_process
    gff_process.importEnsemblGenes(str(exon_annot))
    return gff_process.exonCoordinates


# --------------------------------------------------------------------------- matching

def is_contiguous_sublist(short, long_):
    """True iff ``short`` appears as a consecutive run inside ``long_``.

    Contiguity is the point. A truncated read gives a consecutive run of the parent's junctions; an
    exon-skipping isoform gives an ordered but NON-consecutive subset and must not be folded in.
    """
    n, m = len(short), len(long_)
    if n == 0 or n > m:
        return False
    first = short[0]
    for i in range(m - n + 1):
        if long_[i] == first and long_[i:i + n] == short:
            return True
    return False


class JunctionIndex:
    """Index a set of junction chains so superstring lookup is fast.

    A chain that CONTAINS the query must contain every one of the query's junctions, so the
    candidate set is the intersection of the per-junction posting lists. Only those candidates are
    then tested for contiguity.
    """

    def __init__(self, chains):
        self.chains = list(chains)              # list of (id, junction tuple)
        self.exact = collections.defaultdict(list)
        self.posting = collections.defaultdict(set)
        for i, (tid, junctions) in enumerate(self.chains):
            self.exact[junctions].append(tid)
            for j in junctions:
                self.posting[j].add(i)

    def find_exact(self, junctions):
        return self.exact.get(junctions, [])

    def find_superstrings(self, junctions, limit=50):
        """Ids of indexed chains that contain ``junctions`` as a contiguous run."""
        if not junctions:
            return []
        candidates = None
        for j in junctions:
            posting = self.posting.get(j)
            if not posting:
                return []
            candidates = posting if candidates is None else (candidates & posting)
            if not candidates:
                return []
        hits = []
        for i in candidates:
            tid, cand = self.chains[i]
            if len(cand) > len(junctions) and is_contiguous_sublist(junctions, cand):
                hits.append(tid)
                if len(hits) >= limit:
                    break
        return hits


def classify(query, target_index):
    """Classify one junction chain against an indexed set.

    Returns one of:
      identical                 the same chain exists on the other side
      redundant (substring)     the query is a contiguous run inside a LONGER chain on the other
                                side, so the other side already represents this splice path
      unique                    neither
    """
    if not query:
        return 'no junctions', []
    exact = target_index.find_exact(query)
    if exact:
        return 'identical', exact
    supers = target_index.find_superstrings(query)
    if supers:
        return 'redundant (substring)', supers
    return 'unique', []


# --------------------------------------------------------------------------- driver

def compare(our_gff, query_gff, exon_annot, out_dir, our_catalog=None, reference_counts=None,
            gene_symbols=None,
            require_known_gene=True, require_known_splice=True, require_observed=True,
            top_n=10, log=print):
    """Answer two questions on junction composition alone.

    Q1  How many of OUR final collapsed isoforms have a junction chain the reference does not
        contain, even as part of a longer chain?
    Q2  How many reference isoforms that MEET OUR SEARCH REQUIREMENTS have a junction chain we do
        not contain, even as part of a longer chain?

    our_gff:   the final derived GFF (gff-output/combined.gff.gz).
    query_gff: the reference annotation (an ENCODE TALON GTF).
    """
    import csv
    out_dir = os.path.abspath(str(out_dir))
    os.makedirs(out_dir, exist_ok=True)

    exon_coordinates = load_exon_coordinates(exon_annot)
    log(f"[junctions] Ensembl exon boundaries loaded: {len(exon_coordinates):,}")

    symbols = {}
    if gene_symbols and _io.exists(str(gene_symbols)):
        with _io.smart_open(str(gene_symbols)) as handle:
            for line in handle:
                p = line.rstrip('\n').split('\t')
                if len(p) >= 2 and p[0] and p[1]:
                    symbols.setdefault(p[0], p[1])

    counts = {}
    meta = {}
    if reference_counts:
        with _io.smart_open(str(reference_counts)) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            col = (reader.fieldnames or [])[-1]
            for row in reader:
                tid = row['annot_transcript_id']
                try:
                    counts[tid] = float(row[col])
                except (TypeError, ValueError):
                    counts[tid] = 0.0
                meta[tid] = {'novelty': row.get('transcript_novelty', ''),
                             'n_exons': int(float(row.get('n_exons') or 0))}
        log(f"[junctions] reference quantification: {len(counts):,} transcripts, column '{col}'")

    # Our per-isoform read totals, so Q1 can be weighted by expression rather than counted alone.
    # FINAL_isoform_catalog.tsv keys a NOVEL isoform as '<molecule>.<library>' while the GFF stamps
    # the bare '<molecule>' as transcript_id; a KNOWN isoform is a version-stripped ENST. Register
    # both spellings so either lookup resolves.
    our_reads = {}
    if our_catalog and _io.exists(str(our_catalog)):
        with _io.smart_open(str(our_catalog)) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            for row in reader:
                fid = row['final_isoform_id']
                try:
                    rd = float(row.get('total_reads') or 0)
                except (TypeError, ValueError):
                    rd = 0.0
                our_reads[fid] = rd
                if '.' in fid:
                    our_reads.setdefault(fid.rsplit('.', 1)[0], rd)
        log(f"[junctions] our catalog read totals: {len(our_reads):,} keys, "
            f"{int(sum(v for k, v in our_reads.items() if '.' in k or k.startswith('ENST'))):,} reads")

    ours = parse_gff_junctions(our_gff, log=log)
    theirs = parse_gff_junctions(query_gff, log=log)

    # ---- funnel: apply OUR search requirements to the reference annotation
    funnel = collections.OrderedDict()
    funnel['transcripts in the file'] = len(theirs)
    keep_theirs = {}
    n_multi = n_gene = n_splice = n_obs = 0
    for tid, (chrom, strand, gene, junctions) in theirs.items():
        if not junctions:
            continue
        n_multi += 1
        if require_known_gene and not gene.startswith('ENSG'):
            continue
        n_gene += 1
        if require_known_splice and not has_known_splice_site(chrom, strand, junctions,
                                                              exon_coordinates):
            continue
        n_splice += 1
        if require_observed and counts and counts.get(tid, 0) <= 0:
            continue
        n_obs += 1
        keep_theirs[tid] = (chrom, strand, gene, junctions)
    funnel['spliced (>=1 junction)'] = n_multi
    funnel['in a known Ensembl gene'] = n_gene
    funnel['>=1 known Ensembl splice site'] = n_splice
    funnel['observed (>=1 read)'] = n_obs

    keep_ours = {t: v for t, v in ours.items() if v[3]}
    log("")
    log("REFERENCE FUNNEL: applying our own search requirements to the reference annotation")
    prev = None
    for label, n in funnel.items():
        pct = f"  ({100.0*n/funnel['transcripts in the file']:.1f}% of file)"
        drop = '' if prev is None else f"   dropped {prev-n:,}"
        log(f"  {label:<32}{n:>10,}{pct}{drop}")
        prev = n
    log(f"  {'our final isoforms (spliced)':<32}{len(keep_ours):>10,}")

    # ---- distinct junction chains on each side
    ours_by_chain = collections.defaultdict(list)
    for tid, (_c, _s, gene, junctions) in keep_ours.items():
        ours_by_chain[junctions].append(tid)
    theirs_by_chain = collections.defaultdict(list)
    for tid, (_c, _s, gene, junctions) in keep_theirs.items():
        theirs_by_chain[junctions].append(tid)
    log("")
    log(f"  distinct junction chains, ours      : {len(ours_by_chain):,} "
        f"(from {len(keep_ours):,} isoforms)")
    log(f"  distinct junction chains, reference : {len(theirs_by_chain):,} "
        f"(from {len(keep_theirs):,} transcripts)")

    # Q1 and Q2 need DIFFERENT target sets, and conflating them is wrong.
    #
    # Q2 asks "what did the reference SEE that we missed", so its query set is the reference reduced
    # by our own search requirements, including "observed".
    #
    # Q1 asks "is this chain IN the reference GFF". That is a question about the file, not about
    # what the reference quantified, so its target must be EVERY spliced transcript in the GFF. The
    # observed filter must not apply. Measured consequence of getting this wrong: FTH1
    # ENST00000530019 is present in the ENCODE GTF with identical exon coordinates but carries zero
    # reads in the TALON abundance table, so an observed-filtered target set scored our 19,145-read
    # chain as absent from a file that plainly contains it.
    all_theirs_by_chain = collections.defaultdict(list)
    for tid, (_c, _s, _g, junctions) in theirs.items():
        if junctions:
            all_theirs_by_chain[junctions].append(tid)
    log(f"  distinct junction chains, reference FILE (Q1 target): {len(all_theirs_by_chain):,} "
        f"(from {sum(1 for v in theirs.values() if v[3]):,} spliced transcripts, no observed filter)")

    idx_theirs = JunctionIndex([(t[0], c) for c, t in all_theirs_by_chain.items()])
    idx_ours = JunctionIndex([(t[0], c) for c, t in ours_by_chain.items()])

    def run_side(by_chain, index, gene_of, reads_of, label, out_path, unique_path):
        res = collections.Counter()
        reads = collections.Counter()
        rows_unique = []
        with open(out_path, 'w') as handle:
            handle.write('\t'.join(['gene', 'symbol', 'n_junctions', 'reads', 'verdict',
                                    'matched_on_other_side', 'ids', 'junctions']) + '\n')
            for chain, ids in by_chain.items():
                verdict, hits = classify(chain, index)
                rd = sum(reads_of(i) for i in ids)
                res[verdict] += 1
                reads[verdict] += rd
                gene = gene_of(ids[0])
                jstr = ','.join(f"{a}-{b}" for a, b in chain)
                handle.write('\t'.join([gene, symbols.get(gene, ''), str(len(chain)),
                                        f"{rd:.0f}", verdict, ','.join(hits[:5]),
                                        ','.join(ids[:5]), jstr]) + '\n')
                if verdict == 'unique':
                    rows_unique.append((rd, gene, symbols.get(gene, ''), len(chain), ids, jstr))
        rows_unique.sort(reverse=True)
        with open(unique_path, 'w') as handle:
            handle.write('\t'.join(['rank', 'gene', 'symbol', 'reads', 'n_junctions', 'ids',
                                    'junctions']) + '\n')
            for i, (rd, gene, sym, nj, ids, jstr) in enumerate(rows_unique[:top_n], 1):
                handle.write('\t'.join([str(i), gene, sym, f"{rd:.0f}", str(nj),
                                        ','.join(ids[:5]), jstr]) + '\n')
        total = sum(res.values())
        total_reads = sum(reads.values())
        log("")
        log(label)
        log(f"  {'verdict':<26}{'chains':>10}{'% chains':>11}{'reads':>13}{'% reads':>10}")
        for v in ('identical', 'redundant (substring)', 'unique'):
            if v not in res:
                continue
            log(f"  {v:<26}{res[v]:>10,}{100.0*res[v]/total:>10.1f}%"
                f"{int(reads[v]):>13,}{100.0*reads[v]/total_reads if total_reads else 0:>9.1f}%")
        log(f"  {'TOTAL':<26}{total:>10,}{100.0:>10.1f}%{int(total_reads):>13,}{100.0:>9.1f}%")
        log(f"  wrote {out_path}")
        log(f"  wrote {unique_path}")
        return res, reads, total, total_reads, rows_unique

    q1 = run_side(ours_by_chain, idx_theirs,
                  lambda t: keep_ours[t][2],
                  lambda t: our_reads.get(t, our_reads.get(strip_version(t), 0.0)),
                  "Q1  OUR final isoforms, tested against EVERY spliced transcript in the reference GFF",
                  os.path.join(out_dir, 'ours_vs_reference_junctions.txt'),
                  os.path.join(out_dir, f'top{top_n}_ours_unique.txt'))
    q2 = run_side(theirs_by_chain, idx_ours,
                  lambda t: keep_theirs[t][2], lambda t: counts.get(t, 0.0),
                  "Q2  REFERENCE transcripts meeting our requirements, tested against ours",
                  os.path.join(out_dir, 'reference_vs_ours_junctions.txt'),
                  os.path.join(out_dir, f'top{top_n}_reference_missed.txt'))

    summary = os.path.join(out_dir, 'junction_concordance_summary.txt')
    with open(summary, 'w') as handle:
        handle.write("metric\tvalue\n")
        for label, n in funnel.items():
            handle.write(f"reference_funnel: {label}\t{n}\n")
        handle.write(f"our_isoforms_spliced\t{len(keep_ours)}\n")
        handle.write(f"our_distinct_chains\t{len(ours_by_chain)}\n")
        handle.write(f"reference_distinct_chains_observed_Q2\t{len(theirs_by_chain)}\n")
        handle.write(f"reference_distinct_chains_whole_file_Q1\t{len(all_theirs_by_chain)}\n")
        for name, (res, reads, total, total_reads, _u) in (('Q1_ours', q1), ('Q2_reference', q2)):
            for v in ('identical', 'redundant (substring)', 'unique'):
                handle.write(f"{name}: {v}\t{res.get(v, 0)}\n")
                handle.write(f"{name}: {v} reads\t{int(reads.get(v, 0))}\n")
            handle.write(f"{name}: total\t{total}\n")
    log("")
    log(f"  wrote {summary}")
    return {'funnel': funnel, 'q1': q1[0], 'q2': q2[0]}
