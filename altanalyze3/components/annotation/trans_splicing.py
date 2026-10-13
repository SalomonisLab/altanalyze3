"""Physical splice geometry for separating distant gene joins from annotation overlap.

Report distinct-locus joins as trans-splicing candidates, not proof of mechanism.
The same classifier is shared by event annotation and complete catalog audits.
"""
from __future__ import annotations
import bisect,csv,gzip,re
from collections import defaultdict
from pathlib import Path
from ..long_read import io_utils as IO
ATTR = re.compile(r'(?:^|;)\s*([\w]+)(?:=|\s+")([^";]+)')

def opened(path):
    return IO.smart_open(path)

def gene_id(value):
    return value.split('.', 1)[0] if value.startswith('ENSG') else value

def near(values, position, tolerance):
    index = bisect.bisect_left(values, position - tolerance)
    return index < len(values) and values[index] <= position + tolerance

def merged(intervals):
    out = []
    for start, end in sorted(intervals):
        if start > end:
            raise ValueError('Ascending inclusive genomic exon coordinates required')
        if out and start <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], end))
        else:
            out.append((start, end))
    return out

def token_parts(token, primary):
    owner, label = token.rsplit(':', 1) if ':' in token else (primary, token)
    base, coord = label.rsplit('_', 1) if '_' in label else (label, None)
    return gene_id(owner), base, int(coord) if coord is not None else None

def reference(exon_path, reference_gff):
    regions = defaultdict(dict)
    loci = defaultdict(dict)
    exonic = defaultdict(lambda: defaultdict(list))
    sites = defaultdict(lambda: defaultdict(lambda: {'donor': set(), 'acceptor': set()}))
    with opened(exon_path) as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            gene, label = gene_id(row['gene']), row['exon-id']
            chrom, strand = row['chromosome'], row['strand']
            start, end = sorted((int(row['exon-region-start(s)']), int(row['exon-region-stop(s)'])))
            locus = (chrom, strand)
            if label in regions[gene] and regions[gene][label] != (chrom, strand, start, end):
                raise ValueError('Required exon-region identity is ambiguous: ' + gene + ':' + label)
            regions[gene][label] = (chrom, strand, start, end)
            old = loci[gene].get(locus, (start, end))
            loci[gene][locus] = (min(old[0], start), max(old[1], end))
            if label.startswith('E'):
                exonic[gene][locus].append((start, end))
                sites[gene][locus]['donor'].add(end if strand == '+' else start)
                sites[gene][locus]['acceptor'].add(start if strand == '+' else end)
    # The same original GENCODE reference supplies explicit full gene spans.
    # It extends locus geometry, never replaces the native novelty reference.
    with opened(reference_gff) if reference_gff else _empty_reference() as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            row = line.rstrip('\n').split('\t')
            if len(row) != 9 or row[2] != 'gene':
                continue
            attr = dict(ATTR.findall(row[8]))
            gene = gene_id(attr.get('gene_id', attr.get('ID', '').removeprefix('gene:')))
            if not gene:
                raise ValueError('Original reference gene identity missing')
            locus = row[0], row[6]
            start, end = int(row[3]), int(row[4])
            old = loci[gene].get(locus, (start, end))
            loci[gene][locus] = (min(old[0], start), max(old[1], end))
    for gene in exonic:
        for locus, intervals in exonic[gene].items():
            exonic[gene][locus] = merged(intervals)
    for gene in sites:
        for locus in sites[gene]:
            for role in ('donor', 'acceptor'):
                sites[gene][locus][role] = sorted(sites[gene][locus][role])
    return regions, loci, exonic, sites

def token_boundary(token, primary, role, regions):
    gene, base, coord = token_parts(token, primary)
    record = regions.get(gene, {}).get(base)
    if record is None:
        return gene, None
    chrom, strand, start, end = record
    if coord is None:
        coord = end if (role == 'donor') == (strand == '+') else start
    return gene, (chrom, strand, coord)

def physical_blocks(exons):
    by_locus = defaultdict(list)
    for chrom, strand, start, end in exons:
        by_locus[chrom, strand].append((start, end))
    # Native exported exon rows are genomic-coordinate ordered. Restore the
    # transcript direction for either strand; do not create a join across contigs.
    if len(by_locus) != 1:
        return None
    (chrom, strand), intervals = next(iter(by_locus.items()))
    if strand not in ('+', '-'):
        raise ValueError('Physical model strand unresolved')
    return [(chrom, strand, start, end) for start, end in
            sorted(merged(intervals), reverse=strand == '-')]

def classify_trans(primary, structure, blocks, regions, loci, exonic, tolerance):
    tokens = [t for t in structure.split('|') if t]
    transitions = []
    for left, right in zip(tokens, tokens[1:]):
        a, ca = token_boundary(left, primary, 'donor', regions)
        b, cb = token_boundary(right, primary, 'acceptor', regions)
        if a != b:
            transitions.append((a, ca, b, cb, left, right))
    if not transitions:
        return {'status': 'single_gene', 'candidate': False, 'events': []}
    if blocks is None:
        return {'status': 'unresolved_multicontig_geometry', 'candidate': None, 'events': []}
    junctions = [(left[:2], left[3] if left[1] == '+' else left[2],
                  right[2] if right[1] == '+' else right[3])
                 for left, right in zip(blocks, blocks[1:])]
    events = []
    for a, ca, b, cb, left, right in transitions:
        event = {'gene_a': a, 'gene_b': b, 'token_a': left, 'token_b': right}
        if ca is None or cb is None:
            event['status'] = 'unresolved_token_region'
            events.append(event)
            continue
        la, lb = loci.get(a, {}).get(ca[:2]), loci.get(b, {}).get(cb[:2])
        event.update(boundary_a=ca, boundary_b=cb, gene_locus_a=la, gene_locus_b=lb)
        if la is None or lb is None:
            event['status'] = 'unresolved_gene_locus'
        elif ca[0] != cb[0]:
            # Source collapse supports single-chromosome primary models. Never
            # call a foreign chromosome label physical trans-splicing evidence.
            event['status'] = 'foreign_locus_label_without_physical_join'
        elif not (la[0] <= ca[2] <= la[1] and lb[0] <= cb[2] <= lb[1]):
            event['status'] = 'gene_label_outside_own_locus'
        else:
            lo, hi = max(la[0], lb[0]), min(la[1], lb[1])
            actual = any(locus == ca[:2] and abs(donor-ca[2]) <= tolerance
                         and abs(acceptor-cb[2]) <= tolerance
                         for locus, donor, acceptor in junctions)
            event['actual_splice_join'] = actual
            event['gene_locus_gap_nt'] = max(0, max(la[0], lb[0])-min(la[1], lb[1])-1)
            shared_exonic_locations = all(any(start-tolerance <= position <= end+tolerance
                                              for locus, intervals in exonic.get(owner, {}).items()
                                              if locus[0] == ca[0]
                                              for start, end in intervals)
                                           for owner in (a, b) for position in (ca[2], cb[2]))
            event['both_boundaries_exonic_in_both_gene_annotations'] = shared_exonic_locations
            if lo <= hi and lo-tolerance <= ca[2] <= hi+tolerance and lo-tolerance <= cb[2] <= hi+tolerance and (shared_exonic_locations or not actual):
                event['status'] = 'overlapping_locus_annotation_artifact'
            elif lo <= hi:
                event['status'] = 'overlapping_gene_loci_ambiguous'
            elif actual:
                event['status'] = 'distinct_locus_spliced_candidate'
            else:
                event['status'] = 'distinct_locus_labels_without_physical_splice'
        events.append(event)
    # Token order can include interior annotation blocks and therefore need not
    # reproduce physical exon-junction boundaries. Independently check the
    # actual joins against the explicitly represented owners' gene loci.
    owners = {token_parts(t, primary)[0] for t in tokens} | {primary}
    for locus, donor, acceptor in junctions:
        for a in sorted(owners):
            la = loci.get(a, {}).get(locus)
            if la is None or not la[0] <= donor <= la[1]:
                continue
            for b in sorted(owners - {a}):
                lb = loci.get(b, {}).get(locus)
                if lb is None or not lb[0] <= acceptor <= lb[1]:
                    continue
                if max(la[0], lb[0]) <= min(la[1], lb[1]):
                    continue
                event = {'gene_a': a, 'gene_b': b,
                         'boundary_a': [*locus, donor], 'boundary_b': [*locus, acceptor],
                         'gene_locus_a': la, 'gene_locus_b': lb,
                         'actual_splice_join': True,
                         'gene_locus_gap_nt': max(la[0], lb[0])-min(la[1], lb[1])-1,
                         'status': 'distinct_locus_spliced_candidate',
                         'evidence': 'physical_junction_between_explicitly_represented_gene_loci'}
                if not any(e.get('gene_a') == a and e.get('gene_b') == b
                           and e.get('boundary_a') == event['boundary_a']
                           and e.get('boundary_b') == event['boundary_b'] for e in events):
                    events.append(event)
    statuses = {e['status'] for e in events}
    if 'distinct_locus_spliced_candidate' in statuses:
        result = 'distinct_locus_spliced_candidate'
        candidate = True
    elif statuses == {'overlapping_locus_annotation_artifact'}:
        result = 'overlapping_locus_annotation_artifact'
        candidate = False
    elif any(not loci.get(owner) for owner in owners):
        result = 'ambiguous_or_unsupported_multigene_labels'
        candidate = None
    else:
        # These full source/annotation records lack a join between separate
        # loci. Overlap and unsupported labels are audit outcomes, not missing
        # samples or models, and are not counted as distinct-locus candidates.
        chrom = blocks[0][0]
        spans = [[span for locus, span in loci[owner].items() if locus[0] == chrom]
                 for owner in sorted(owners)]
        all_overlap = (all(len(group) == 1 for group in spans)
                       and all(max(a[0][0], b[0][0]) <= min(a[0][1], b[0][1])
                               for i, a in enumerate(spans) for b in spans[i+1:]))
        result = ('overlapping_gene_labels_not_distinct_loci'
                  if all_overlap or statuses <= {'overlapping_gene_loci_ambiguous', 'overlapping_locus_annotation_artifact'}
                  else 'unsupported_multigene_annotation_labels')
        candidate = False
    return {'status': result, 'candidate': candidate, 'events': events}

def characterize(primary, structure, exons, reference_data, tolerance=3):
    regions, loci, exonic, sites = reference_data
    owners = {token_parts(t, primary)[0] for t in structure.split('|') if t} | {primary}
    blocks = physical_blocks(exons)
    trans = classify_trans(primary, structure, blocks, regions, loci, exonic, tolerance)
    if blocks is None:
        return None, None, trans, [], []
    # Retain explicit uncertainty rather than silently omit models/partners.
    missing = [g for g in owners if g not in regions]
    if missing:
        trans['missing_native_gene_reference'] = sorted(missing)
        return None, None, trans, [], []
    new_sites, new_exons = [], []
    for left, right in zip(blocks, blocks[1:]):
        locus = left[:2]
        for role, pos in (('donor', left[3] if left[1] == '+' else left[2]),
                          ('acceptor', right[2] if right[1] == '+' else right[3])):
            if not any(near(sites.get(g, {}).get(locus, {}).get(role, []), pos, tolerance) for g in owners):
                new_sites.append({'chrom': locus[0], 'strand': locus[1], 'role': role, 'position': pos})
    for chrom, strand, start, end in blocks:
        intervals = [interval for g in owners for interval in exonic.get(g, {}).get((chrom, strand), [])]
        if not any(start <= stop and begin <= end for begin, stop in intervals):
            new_exons.append((chrom, strand, start, end))
    return bool(new_sites), bool(new_exons), trans, new_sites, new_exons


def _empty_reference():
    import io
    return io.StringIO('')


class JunctionGeometry:
    def __init__(self, exon_reference, reference_gff=None, tolerance=3):
        self.data=reference(exon_reference,reference_gff)
        self.tolerance=tolerance

    def classify(self,junction,coordinates,strand):
        primary,body=junction.split(':',1)
        left,right=body.split('-',1)
        owners={token_parts(left,primary)[0],token_parts(right,primary)[0]}
        if len(owners)==1:return {'status':'single_gene','candidate':False,'events':[]}
        match=re.fullmatch(r'([^:]+):(\d+)-(\d+)',str(coordinates))
        if not match:return {'status':'unresolved_physical_junction_coordinates','candidate':None,'events':[]}
        chrom,a,b=match.groups();a,b=int(a),int(b)
        if a==b:return {'status':'unsupported_zero_length_junction','candidate':False,'events':[]}
        # One-base endpoint blocks encode only the observed splice edge, not
        # invented transcript extent. All model-level classification uses GFF.
        blocks=[(chrom,strand,a,a),(chrom,strand,b,b)]
        regions,loci,exonic,_=self.data
        return classify_trans(primary,left+'|'+right,blocks,regions,loci,exonic,self.tolerance)
