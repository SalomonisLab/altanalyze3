"""EXPERIMENTAL FORK (extra-chr-trans-splicing, added 2026-11-05). NOT validated altanalyze3 code.

Same-chromosome chimeric molecules in long-read single-cell BAMs: circRNA-like back-splice joints,
lariat-like intronic joints, distal forward joints and inverted joints.

The validated extractor skips supplementary alignments. The fork's cross-chromosome mode keeps only
molecules whose pieces lie on two or more chromosomes. A read across a circRNA back-splice junction
aligns as two or more pieces on ONE chromosome and strand, and the 3' piece of the read lies upstream
of the 5' piece in transcript order. This module reads every primary record with an SA tag, rebuilds
its pieces with the fork helpers in isoform_structure_extract (_alignment_segments, _segment_junction,
_trim_to_junction, unchanged), and writes one row per read-adjacent pair of pieces on one chromosome.

Joint geometry, transcript orientation (donor = 3' end of the 5' piece, acceptor = 5' start of the
3' piece, after the fork's overlap resolution):
  back      same strand, acceptor at or upstream of the donor (non-colinear)
  forward   same strand, acceptor downstream of the donor (colinear, split beyond one alignment)
  inverted  opposite strands

Molecule evidence recorded for a back joint (no threshold is applied here):
  circle span            genomic interval acceptor..donor
  outside_circle_nt      aligned reference bases of the WHOLE molecule outside the circle span; a piece
                         on another chromosome or strand counts in full. 0 = every aligned base lies in
                         the span, as for a circular template. A linear molecule with a non-colinear
                         joint (exon repetition, template switch) usually extends past the span.
  crossings              adjacent pairs of this molecule at the same back joint (>= 2 = rolling RT)
  pair_overlap_nt        reference bases covered by both pieces (> 0 = more than one pass)
  lariat_like            acceptor is the first base of an annotated intron and the donor lies inside
                         that intron (5' splice site joined to a branch-point region)
Splice evidence: registered boundaries (Ensembl exon-model keys, as in the fork), genome dinucleotides
and junction homology (arithmetic of scripts/report_stream.py dinucleotides/junction_homology), and the
polyA/T concatemer signature (cross_chromosome_rules.concatemer_signature).
"""
import argparse
import gzip
import itertools
import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import pysam

sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from altanalyze3.components.bam import isoform_structure_extract as ise  # noqa: E402
from altanalyze3.components.bam import cross_chromosome_rules as R  # noqa: E402
from altanalyze3.components.long_read import gff_process  # noqa: E402

_HERE = Path(__file__).resolve().parent
for _module in (ise, R):
    if Path(_module.__file__).resolve().parent != _HERE:
        raise ImportError(f'{_module.__name__} resolved outside the fork: {_module.__file__}')

BIN = 10_000

COLUMNS = [
    'library', 'molecule_id', 'query_name', 'barcode', 'umi', 'status', 'read_length', 'n_segments',
    'n_chromosomes', 'pair_index', 'joint_type', 'joint_class', 'chrom',
    'a_strand', 'a_qstart', 'a_qend', 'a_mapq', 'a_nm', 'a_primary', 'a_cigar', 'a_exons',
    'b_strand', 'b_qstart', 'b_qend', 'b_mapq', 'b_nm', 'b_primary', 'b_cigar', 'b_exons',
    'query_gap', 'overlap_nt', 'shift_nt', 'aligned_donor', 'aligned_acceptor', 'slide_left_nt', 'slide_right_nt',
    'slide_nt', 'donor', 'acceptor', 'resolved',
    'donor_annotation', 'donor_is_gene_first_region', 'donor_is_gene_last_region',
    'acceptor_annotation', 'acceptor_is_gene_first_region', 'acceptor_is_gene_last_region',
    'donor_region', 'acceptor_region', 'donor_genes', 'acceptor_genes', 'same_gene',
    'donor_dinucleotide', 'acceptor_dinucleotide', 'canonical_dinucleotides', 'junction_homology_nt',
    'polyA_at_joint', 'junction_read_context', 'genomic_distance', 'pair_span_overlap', 'circle_start', 'circle_end',
    'circle_span_nt', 'outside_circle_nt', 'crossings', 'pair_overlap_nt', 'lariat_like',
    'donor_to_intron_end_nt', 'unaligned_read_nt']


class BinIndex:
    """Intervals per chromosome in fixed 10-kb bins; point query returns payloads containing pos.
    The fork's IntervalIndex scans back by the longest interval, which is slow for megabase introns."""
    def __init__(self, rows, size=BIN):
        self.size = size
        self.bins = defaultdict(list)
        for chrom, lo, hi, payload in rows:
            for b in range(lo // size, hi // size + 1):
                self.bins[(chrom, b)].append((lo, hi, payload))

    def at(self, chrom, pos):
        return [p for lo, hi, p in self.bins.get((chrom, pos // self.size), ()) if lo <= pos <= hi]


class Genome:
    """Genome lookups with the arithmetic of scripts/report_stream.py (fa_name, dinucleotides,
    oriented, junction_homology), so values match the cross-chromosome report."""
    def __init__(self, fasta_path):
        self.fasta = pysam.FastaFile(str(fasta_path))
        self.names = set(self.fasta.references)
        self.cache = {}

    def name(self, c):
        bare = c[3:] if c.startswith('chr') else c
        if bare in self.names:
            return bare
        return 'MT' if c in ('chrM', 'M') else c

    def oriented(self, chrom, start0, end0, strand):
        c = self.name(chrom)
        if c not in self.names:
            return ''
        piece = self.fasta.fetch(c, max(0, start0), max(0, end0)).upper()
        return piece if strand == '+' else R.revcomp(piece)

    def dinucleotides(self, chrom, strand, pos, kind):
        """Two intron bases next to a splice site, transcript-oriented. '' (undefined, counted by the
        caller) when the contig is absent or the two bases fall outside it, e.g. a joint at the chrM
        origin; report_stream.py raises there."""
        c = self.name(chrom)
        if c not in self.names:
            return ''
        before = (kind == 'donor') != (strand == '+')
        start = pos - 3 if before else pos
        if start < 0 or start + 2 > self.fasta.get_reference_length(c):
            return ''
        s = self.fasta.fetch(c, start, start + 2)
        return (s if strand == '+' else R.revcomp(s)).upper()

    def slide(self, ca, sa, d, cb, sb, a, span=20):
        """(left, right): bases the breakpoint can move upstream / downstream over genome sequence
        shared by both partners, at most span each. left + right = junction_homology."""
        a_exon_end = self.oriented(ca, d - span, d, '+') if sa == '+' else self.oriented(ca, d - 1, d - 1 + span, '-')
        a_intron = self.oriented(ca, d, d + span, '+') if sa == '+' else self.oriented(ca, d - 1 - span, d - 1, '-')
        b_intron = self.oriented(cb, a - 1 - span, a - 1, '+') if sb == '+' else self.oriented(cb, a, a + span, '-')
        b_exon = self.oriented(cb, a - 1, a - 1 + span, '+') if sb == '+' else self.oriented(cb, a - span, a, '-')
        return R.common_suffix(a_exon_end, b_intron), R.common_prefix(a_intron, b_exon)

    def junction_homology(self, ca, sa, d, cb, sb, a, span=20):
        return sum(self.slide(ca, sa, d, cb, sb, a, span))

    def canonical(self, chrom, a_strand, donor, b_strand, acceptor):
        key = ('canonical', chrom, a_strand, donor, b_strand, acceptor)
        if key not in self.cache:
            self.cache[key] = R.canonical_splice(self.dinucleotides(chrom, a_strand, donor, 'donor'),
                                                 self.dinucleotides(chrom, b_strand, acceptor, 'acceptor'))
        return self.cache[key]

    def joint(self, chrom, a_strand, donor, b_strand, acceptor):
        key = (chrom, a_strand, donor, b_strand, acceptor)
        if key not in self.cache:
            d = self.dinucleotides(chrom, a_strand, donor, 'donor')
            a = self.dinucleotides(chrom, b_strand, acceptor, 'acceptor')
            self.cache[key] = (d, a, R.canonical_splice(d, a),
                               self.junction_homology(chrom, a_strand, donor, chrom, b_strand, acceptor))
        return self.cache[key]


class Annotation:
    """Ensembl exon model as loaded by the validated gff_process.importEnsemblGenes."""
    def __init__(self, exon_file):
        self.exon_coordinates, gene_data, strand_data = gff_process.importEnsemblGenes(exon_file)
        gene_chr = gff_process.gene_chr
        self.gene_region_counts = {g: len(r) for g, r in gene_data.items()}
        self.intron_starts = {}
        rows = []
        for gene, regions in gene_data.items():
            chrom, strand = gene_chr[gene], strand_data[gene]
            run = None
            for start, stop, region in regions:
                rows.append((chrom, min(start, stop), max(start, stop), (gene, strand, region)))
                if region.startswith('I'):
                    # consecutive intron sub-regions (transcript order) form one intron
                    run = [gene, region, start, stop] if run is None else run[:3] + [stop]
                    self.intron_starts[(chrom, run[2], strand)] = tuple(run)
                else:
                    run = None
        self.regions = BinIndex(rows)
        self.symbols = {}
        annotations = Path(exon_file).parent / Path(exon_file).name.replace('_exon.txt', '-annotations.txt')
        if annotations.exists():
            with open(annotations) as handle:
                for line in handle:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) > 1 and parts[0].startswith('ENSG') and parts[1]:
                        self.symbols[parts[0]] = parts[1]
        self.cache = {}

    def at(self, chrom, strand, pos):
        key = (chrom, strand, pos)
        if key not in self.cache:
            hits = sorted({(g, r) for g, s, r in self.regions.at(chrom, pos) if s == strand})
            self.cache[key] = ('|'.join(f'{self.symbols.get(g, g)}:{r}' for g, r in hits),
                               sorted({g for g, _ in hits}))
        return self.cache[key]


def _normalize_breakpoint(a, b, junction, annotation, genome, span=20):
    """The aligner places a breakpoint anywhere inside genome sequence shared by both partners. For
    pieces that meet on the read (query gap <= 0 after overlap resolution), test every position in that
    window and keep the one with the most registered boundaries, then canonical dinucleotides, then the
    smallest move (the aligner's own placement wins ties). Positions that would empty a piece's terminal
    block are not tested. Records the window and the move in junction['slide']."""
    chrom = a['chrom']
    d0, a0 = junction['donor'], junction['acceptor']
    left, right = genome.slide(chrom, a['strand'], d0, chrom, b['strand'], a0, span)
    junction['slide'] = (left, right, 0)
    if junction['gap'] > 0 or not (left or right):
        return
    step_a = 1 if a['strand'] == '+' else -1
    step_b = 1 if b['strand'] == '+' else -1
    a_last = a['exons'][-1] if a['strand'] == '+' else a['exons'][0]
    b_first = b['exons'][0] if b['strand'] == '+' else b['exons'][-1]
    a_room = d0 - a_last[0] if step_a == 1 else a_last[1] - d0
    b_room = b_first[1] - a0 if step_b == 1 else a0 - b_first[0]
    best = None
    for s in range(-min(left, a_room), min(right, b_room) + 1):
        d, acc = d0 + s * step_a, a0 + s * step_b
        dh = annotation.exon_coordinates.get((chrom, d, a['strand'], 2))
        ah = annotation.exon_coordinates.get((chrom, acc, b['strand'], 1))
        score = ((dh is not None) + (ah is not None), genome.canonical(chrom, a['strand'], d, b['strand'], acc), -abs(s))
        if best is None or score > best[0]:
            best = (score, s, d, acc, dh, ah)
    _, s, d, acc, dh, ah = best
    junction.update(donor=d, acceptor=acc, donor_hit=dh, acceptor_hit=ah, slide=(left, right, s))


def _overlap_nt(intervals_a, intervals_b):
    return sum(max(0, min(y1, y2) - max(x1, x2) + 1) for x1, y1 in intervals_a for x2, y2 in intervals_b)


def _outside_nt(segments, chrom, strand, lo, hi):
    out = 0
    for s in segments:
        for x, y in s['exons']:
            if s['chrom'] != chrom or s['strand'] != strand:
                out += y - x + 1
            else:
                out += max(0, min(y, lo - 1) - x + 1) + max(0, y - max(x, hi + 1) + 1)
    return out


def _unaligned_nt(segments):
    covered, end = 0, -1
    for qs, qe in sorted((s['qstart'], s['qend']) for s in segments):
        if qe > end:
            covered += qe - max(qs, end)
            end = qe
    return segments[0]['read_length'] - covered


def joint_class(row):
    """Descriptive class of one same-chromosome pair. chrM: mitochondrial_<type>. back: lariat_like; donor and acceptor in different
    genes (backjoint_different_genes, '_registered_canonical' when both registered and canonical); then
    registered donor AND acceptor (canonical or not), one registered side, none (the last two split by
    canonical dinucleotides, '_canonical'). forward: same gene, different genes, or
    a side outside any gene. inverted: pieces overlapping on the genome or separate."""
    t = row['joint_type']
    if row['chrom'] in ('chrM', 'MT', 'M'):
        # circular mitochondrial genome and polycistronic, non-spliceosomal transcripts (fork blacklist rule)
        return 'mitochondrial_' + t
    if t == 'back':
        if row['lariat_like']:
            return 'lariat_like'
        registered = bool(row['donor_annotation'] and row['acceptor_annotation'])
        if row['donor_genes'] and row['acceptor_genes'] and not row['same_gene']:
            return 'backjoint_different_genes' + ('_registered_canonical' if registered and row['canonical_dinucleotides'] else '')
        if registered:
            return 'backsplice_registered_canonical' if row['canonical_dinucleotides'] else 'backsplice_registered_noncanonical'
        suffix = '_canonical' if row['canonical_dinucleotides'] else ''
        if row['donor_annotation'] or row['acceptor_annotation']:
            return 'backjoint_one_side_registered' + suffix
        return 'backjoint_unregistered' + suffix
    if t == 'forward':
        if row['same_gene']:
            return 'forward_same_gene'
        if row['donor_genes'] and row['acceptor_genes']:
            return 'forward_different_genes'
        return 'forward_side_outside_genes'
    return 'inverted_overlapping' if row['pair_span_overlap'] else 'inverted_separate'


def analyze_read(read, annotation, genome, library, barcode_tags=('CB',), molecule_tag='zm', min_mapq=1,
                 ts_read_relative=False, default_barcode=None):
    """Rows (dicts, COLUMNS order) for one primary record with an SA tag; [] when no read-adjacent
    pair of pieces lies on one chromosome."""
    sa_entries = read.get_tag('SA').rstrip(';').split(';')
    segments = ise._alignment_segments(read, sa_entries, ts_read_relative, annotation.exon_coordinates)
    pairs = list(zip(segments[:-1], segments[1:]))
    if not any(a['chrom'] == b['chrom'] for a, b in pairs):
        return []
    status = 'kept'
    primary = next(s for s in segments if s['primary'])
    if primary['strand'] != primary['orientation']:
        status = 'strand_tag_differs_from_alignment'
    elif len({s['read_length'] for s in segments}) != 1:
        status = 'read_length_mismatch_between_alignments'
    barcode = ise.resolve_barcode(read, list(barcode_tags), default_barcode=default_barcode)
    if status == 'kept' and not barcode:
        status = 'no_barcode'
    if status == 'kept' and min(s['mapq'] for s in segments) < min_mapq:
        status = 'mapq_below_min'
    junctions = [ise._segment_junction(a, b, annotation.exon_coordinates) for a, b in pairs]
    for (a, b), junction in zip(pairs, junctions):
        junction['aligned'] = (junction['donor'], junction['acceptor'])
        junction['slide'] = ('', '', 0)
        if a['chrom'] == b['chrom'] and junction['resolved']:
            _normalize_breakpoint(a, b, junction, annotation, genome)
    for (a, b), junction in zip(pairs, junctions):
        ise._trim_to_junction(a, b, junction)
    sequence = ise._oriented_sequence(read, segments[0]['read_length'])
    molecule_id = ise.sanitize_gff_value(ise.resolve_molecule_id(read, molecule_tag))
    umi = read.get_tag('XM') if read.has_tag('XM') else ''
    n_chrom = len({s['chrom'] for s in segments})
    unaligned = _unaligned_nt(segments)
    rows = []
    for index, ((a, b), junction) in enumerate(zip(pairs, junctions), 1):
        if a['chrom'] != b['chrom']:
            continue
        chrom = a['chrom']
        donor, acceptor = junction['donor'], junction['acceptor']
        d_ann, d_first, d_last = ise._boundary_annotation(junction['donor_hit'], annotation.gene_region_counts)
        a_ann, a_first, a_last = ise._boundary_annotation(junction['acceptor_hit'], annotation.gene_region_counts)
        d_region, d_genes = annotation.at(chrom, a['strand'], donor)
        a_region, a_genes = annotation.at(chrom, b['strand'], acceptor)
        d_dinuc, a_dinuc, canonical, homology = genome.joint(chrom, a['strand'], donor, b['strand'], acceptor)
        context = ''
        moved = junction['slide'][2]
        if sequence:
            left_end = (a['qend'] - (junction['shift'] or 0) if junction['overlap'] else a['qend']) + moved
            right_start = left_end if junction['overlap'] else b['qstart'] + moved
            context = (sequence[max(0, left_end - 20):left_end] + '|' +
                       sequence[left_end:right_start].lower() + '|' + sequence[right_start:right_start + 20])
        a_span = (min(x for x, _ in a['exons']), max(y for _, y in a['exons']))
        b_span = (min(x for x, _ in b['exons']), max(y for _, y in b['exons']))
        row = dict(library=library, molecule_id=molecule_id, query_name=read.query_name, barcode=barcode or '',
                   umi=umi, status=status, read_length=a['read_length'], n_segments=len(segments),
                   n_chromosomes=n_chrom, pair_index=index, chrom=chrom)
        for tag, s in (('a', a), ('b', b)):
            row.update({f'{tag}_strand': s['strand'], f'{tag}_qstart': s['qstart'], f'{tag}_qend': s['qend'],
                        f'{tag}_mapq': s['mapq'], f'{tag}_nm': s['nm'], f'{tag}_primary': s['primary'],
                        f'{tag}_cigar': s['cigar'], f'{tag}_exons': ','.join(f'{x}-{y}' for x, y in s['exons'])})
        row.update(query_gap=junction['gap'], overlap_nt=junction['overlap'], shift_nt=junction['shift'],
                   aligned_donor=junction['aligned'][0], aligned_acceptor=junction['aligned'][1],
                   slide_left_nt=junction['slide'][0], slide_right_nt=junction['slide'][1], slide_nt=moved,
                   donor=donor, acceptor=acceptor, resolved=junction['resolved'],
                   donor_annotation=d_ann, donor_is_gene_first_region=d_first, donor_is_gene_last_region=d_last,
                   acceptor_annotation=a_ann, acceptor_is_gene_first_region=a_first,
                   acceptor_is_gene_last_region=a_last, donor_region=d_region, acceptor_region=a_region,
                   donor_genes='|'.join(d_genes), acceptor_genes='|'.join(a_genes),
                   same_gene=bool(set(d_genes) & set(a_genes)), donor_dinucleotide=d_dinuc,
                   acceptor_dinucleotide=a_dinuc, canonical_dinucleotides=canonical,
                   junction_homology_nt=homology, polyA_at_joint=R.concatemer_signature(context) if context else '',
                   junction_read_context=context, unaligned_read_nt=unaligned,
                   pair_span_overlap=a_span[0] <= b_span[1] and b_span[0] <= a_span[1],
                   circle_start='', circle_end='', circle_span_nt='', outside_circle_nt='', crossings='',
                   pair_overlap_nt='', lariat_like=False, donor_to_intron_end_nt='')
        if a['strand'] != b['strand']:
            row['joint_type'] = 'inverted'
            row['genomic_distance'] = ''
        else:
            strand = a['strand']
            sign = 1 if strand == '+' else -1
            if sign * (acceptor - donor) > 0:
                row['joint_type'] = 'forward'
                row['genomic_distance'] = sign * (acceptor - donor) - 1
            else:
                row['joint_type'] = 'back'
                row['genomic_distance'] = sign * (acceptor - donor) - 1
                lo, hi = min(acceptor, donor), max(acceptor, donor)
                row.update(circle_start=lo, circle_end=hi, circle_span_nt=hi - lo + 1,
                           outside_circle_nt=_outside_nt(segments, chrom, strand, lo, hi),
                           pair_overlap_nt=_overlap_nt(a['exons'], b['exons']))
                # lariat 5' end = first intron base after a registered 5' splice site, never itself a registered
                # exon start (Ensembl 91 holds overlapping intron regions, e.g. TET2 I6.2 starts at exon E6.1)
                intron = annotation.intron_starts.get((chrom, acceptor, strand))
                if (intron is not None and junction['acceptor_hit'] is None
                        and (chrom, acceptor - sign, strand, 2) in annotation.exon_coordinates):
                    _, _, ts_start, ts_end = intron
                    if sign * (donor - ts_start) >= 0 and sign * (ts_end - donor) >= 0:
                        row['lariat_like'] = True
                        row['donor_to_intron_end_nt'] = sign * (ts_end - donor)
        row['joint_class'] = joint_class(row)
        rows.append(row)
    back = [(r['chrom'], r['a_strand'], r['donor'], r['acceptor']) for r in rows if r['joint_type'] == 'back']
    for r in rows:
        if r['joint_type'] == 'back':
            r['crossings'] = back.count((r['chrom'], r['a_strand'], r['donor'], r['acceptor']))
    return rows


def extract_noncolinear(bam_path, exon_file, genome_fasta, output_prefix, library, min_mapq=1,
                        barcode_tags=('CB',), molecule_tag='zm', ts_read_relative=False, default_barcode=None,
                        regions=None):
    """Writes <prefix>.noncolinear_joints.tsv.gz (one row per read-adjacent same-chromosome pair of
    every primary record with an SA tag) and <prefix>.noncolinear_stats.json. regions (list of
    (chrom, start0, end)) restricts the scan for tests; supplementary records found in a region are
    traced to their primary record, as in the fork's cross-chromosome region mode."""
    started = time.time()
    output_prefix = Path(output_prefix)
    table_path = Path(str(output_prefix) + '.noncolinear_joints.tsv.gz')
    stats_path = Path(str(output_prefix) + '.noncolinear_stats.json')
    denominator_path = Path(str(output_prefix) + '.barcode_denominators.tsv.gz')
    for path in (table_path, stats_path, denominator_path):
        if path.exists():
            raise FileExistsError(f'Refusing to overwrite an existing output: {path}')
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    ise._set_pysam_verbosity()
    annotation = Annotation(exon_file)
    genome = Genome(genome_fasta)
    bam = pysam.AlignmentFile(str(bam_path), 'rb')
    lookup = pysam.AlignmentFile(str(bam_path), 'rb') if regions else None
    records = itertools.chain.from_iterable(bam.fetch(c, s, e) for c, s, e in regions) if regions else bam.fetch(until_eof=True)
    stats = defaultdict(int)
    # Additive evidence only: count every mapped, nonsecondary primary before the SA
    # branch. Missing barcodes remain an explicit empty-barcode bucket. These are
    # read denominators, not UMI/molecule counts or an approved normalization.
    denominators = defaultdict(lambda: defaultdict(int))
    seen = set()
    with gzip.open(table_path, 'wt') as table:
        table.write('\t'.join(COLUMNS) + '\n')
        for read in records:
            stats['records_scanned'] += 1
            if read.is_unmapped or read.is_secondary:
                continue
            if read.is_supplementary:
                if not regions or not read.has_tag('SA'):
                    continue
                read = ise._primary_for_supplementary(lookup, read)
                if read is None:
                    stats['region_primary_not_found'] += 1
                    continue
            if regions:
                if read.query_name in seen:
                    continue
                seen.add(read.query_name)
            stats['primary_records'] += 1
            denominator_barcode = ise.resolve_barcode(read, list(barcode_tags), default_barcode=default_barcode) or ''
            denominator = denominators[denominator_barcode]
            denominator['mapped_primary_reads'] += 1
            denominator['primary_reads_with_SA'] += int(read.has_tag('SA'))
            denominator['primary_reads_mapq_ge_min'] += int(read.mapping_quality >= min_mapq)
            if not read.has_tag('SA'):
                continue
            stats['primary_records_with_SA'] += 1
            rows = analyze_read(read, annotation, genome, library, barcode_tags=barcode_tags,
                                molecule_tag=molecule_tag, min_mapq=min_mapq, ts_read_relative=ts_read_relative,
                                default_barcode=default_barcode)
            if not rows:
                stats['reads_without_same_chromosome_pair'] += 1
                continue
            stats['reads_with_same_chromosome_pair'] += 1
            stats['reads_status:' + rows[0]['status']] += 1
            for row in rows:
                stats['pairs'] += 1
                if not row['donor_dinucleotide'] or not row['acceptor_dinucleotide']:
                    stats['pairs_dinucleotide_undefined'] += 1
                stats[f"pairs:{row['status']}:{row['joint_class']}"] += 1
                table.write('\t'.join(str(row[c]) for c in COLUMNS) + '\n')
    bam.close()
    if lookup is not None:
        lookup.close()
    if sum(d['mapped_primary_reads'] for d in denominators.values()) != stats.get('primary_records', 0):
        raise AssertionError('barcode denominator primary-count partition failed')
    if sum(d['primary_reads_with_SA'] for d in denominators.values()) != stats.get('primary_records_with_SA', 0):
        raise AssertionError('barcode denominator SA-count partition failed')
    with gzip.open(denominator_path, 'wt') as handle:
        handle.write('library\tbarcode\tbarcode_missing\tmapped_primary_reads\tprimary_reads_with_SA\tprimary_reads_mapq_ge_min\tmin_mapq\n')
        for barcode, counts in sorted(denominators.items()):
            handle.write('\t'.join(map(str, [library, barcode, not bool(barcode), counts['mapped_primary_reads'],
                counts['primary_reads_with_SA'], counts['primary_reads_mapq_ge_min'], min_mapq])) + '\n')
    out = dict(sorted(stats.items()))
    out.update(bam=str(bam_path), library=library, exon_file=str(exon_file), genome_fasta=str(genome_fasta),
               regions=regions, min_mapq=min_mapq, molecule_tag=molecule_tag, barcode_tags=list(barcode_tags),
               table=str(table_path), elapsed_seconds=round(time.time() - started, 1),
               module=str(Path(__file__).resolve()))
    out['barcode_denominators'] = str(denominator_path)
    stats_path.write_text(json.dumps(out, indent=1, default=str) + '\n')
    return table_path, out


def _parse_regions(text):
    if not text:
        return None
    regions = []
    for item in text.split(','):
        chrom, span = item.rsplit(':', 1)
        start, end = span.replace('_', '').split('-')
        regions.append((chrom, int(start) - 1, int(end)))
    return regions


def main():
    p = argparse.ArgumentParser(description='EXPERIMENTAL FORK: same-chromosome chimeric / back-splice joints')
    p.add_argument('--bam', required=True)
    p.add_argument('--gene-model', required=True, help='Ensembl exon model (Hs_Ensembl_exon.txt)')
    p.add_argument('--genome-fasta', required=True)
    p.add_argument('--output-prefix', required=True)
    p.add_argument('--library', required=True)
    p.add_argument('--min-mapq', type=int, default=1)
    p.add_argument('--barcode-tags', default='CB')
    p.add_argument('--molecule-tag', default='zm')
    p.add_argument('--ts-read-relative', action='store_true')
    p.add_argument('--regions', default=None, help='tests only: comma-separated chrom:start-end (1-based)')
    args = p.parse_args()
    _, stats = extract_noncolinear(args.bam, args.gene_model, args.genome_fasta, args.output_prefix, args.library,
                                   min_mapq=args.min_mapq, barcode_tags=[t for t in args.barcode_tags.split(',') if t],
                                   molecule_tag=args.molecule_tag, ts_read_relative=args.ts_read_relative,
                                   regions=_parse_regions(args.regions))
    print(json.dumps(stats, indent=1, default=str))


if __name__ == '__main__':
    main()
