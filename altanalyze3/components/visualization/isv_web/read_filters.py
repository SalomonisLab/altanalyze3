"""Optional read-display filters; scientific assignments and source counts are untouched."""
from collections import defaultdict
from copy import deepcopy
import math


def _chrom(value):
    value = str(value or '').upper().removeprefix('CHR')
    return 'M' if value == 'MT' else value


def blocks_by_chromosome(segments, models, default_gene):
    grouped = defaultdict(list)
    for segment in segments or []:
        gene = segment.get('gene') or default_gene
        chrom = _chrom(segment.get('chrom') or models.get(gene, {}).get('chrom'))
        if not chrom:
            raise ValueError('Verified chromosome required for read-display filtering: ' + str(gene))
        a, b = sorted((segment['start'], segment['end']))
        if not math.isfinite(a) or not math.isfinite(b):
            raise ValueError('Finite original genomic coordinates required')
        if 'genomic_start' in segment:
            raise ValueError('Apply the bleeding-read filter BEFORE display-coordinate projection')
        grouped[chrom].append((a, b))
    result = {}
    for chrom, values in grouped.items():
        merged = []
        for a, b in sorted(values):
            if merged and a <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
            else:
                merged.append((a, b))
        result[chrom] = merged
    return result


def bleeding_intervals(read, parent, tolerance=0):
    """An aligned block crosses a model splice edge into its removed intron.

    No minimum retention length. Transcript outer ends and wholly intronic
    novel exons are not splice-edge extensions. Inclusive genomic coordinates.
    """
    if not isinstance(tolerance, int) or tolerance < 0:
        raise ValueError('Bleed tolerance must be a nonnegative integer')
    hits = []
    for chrom, blocks in read.items():
        model = parent.get(chrom, [])
        for left, right in zip(model, model[1:]):
            donor, acceptor = left[1], right[0]
            for a, b in blocks:
                lo, hi = max(a, donor + 1 + tolerance), min(b, acceptor - 1 - tolerance)
                crosses_edge = (a <= donor < b) or (a < acceptor <= b)
                if crosses_edge and lo <= hi:
                    hits.append(dict(chrom=chrom, model_splice=[donor, acceptor],
                                     aligned_block=[a, b], intronic_extension=[lo, hi]))
    return hits


def filter_bleeding_reads(response, *, enabled=True, tolerance=0):
    """Return a drawing copy, retaining an audit of every hidden source row.

    Run on the FULL eligible read pool, before capping/downsampling. All model
    rows and unfiltered matrix-derived structure_reads are preserved.
    """
    if not enabled:
        return response
    models = {m['gene']: m for m in response.get('gene_models', [])}
    if response.get('gene_model', {}).get('gene'):
        models.setdefault(response['gene_model']['gene'], response['gene_model'])
    default_gene = response.get('gene')
    references = {m.get('final_isoform_id') or m.get('isoform_id'): m
                  for m in response.get('reference_isoforms', [])}
    parents = {}
    decisions = []
    for panel in response.get('panels', []):
        rows = []
        for row in panel.get('molecules', []):
            tid = row.get('final_isoform_id')
            if not tid:
                raise ValueError('Assigned isoform identity required for the bleeding-read filter')
            parent = row.get('registered_exon_segments') or references.get(tid, {}).get('registered_exon_segments')
            if not parent:
                raise ValueError('Exact registered model required for the bleeding-read filter: ' + tid)
            if tid not in parents:
                parents[tid] = blocks_by_chromosome(parent, models, default_gene)
            read = blocks_by_chromosome(row.get('filter_alignment_segments') or row.get('exon_segments'), models, default_gene)
            if not read:
                raise ValueError('Aligned read coordinates required: ' + str(row.get('isoform_id')))
            rows.append(bleeding_intervals(read, parents[tid], tolerance))
        decisions.append(rows)
    result = deepcopy(response)
    result['hidden_bleeding_reads'] = []
    total = dict(input_rows=0, input_RNA=0, hidden_rows=0, hidden_RNA=0, eligible_rows=0, eligible_RNA=0)
    for panel, choices in zip(result.get('panels', []), decisions):
        kept = []; stats = defaultdict(lambda: dict(input_rows=0, input_RNA=0, hidden_rows=0, hidden_RNA=0, eligible_rows=0, eligible_RNA=0))
        for row, hits in zip(panel.get('molecules', []), choices):
            s = stats[row['final_isoform_id']]; count = float(row.get('count', 1))
            for obj in (s, total):obj['input_rows'] += 1; obj['input_RNA'] += count
            if hits:
                result['hidden_bleeding_reads'].append(dict(condition=panel.get('condition'),
                    isoform_id=row['isoform_id'], final_isoform_id=row['final_isoform_id'], count=count,
                    detections=row.get('detections', []), bleeding_intervals=hits))
                for obj in (s, total):obj['hidden_rows'] += 1; obj['hidden_RNA'] += count
            else:
                kept.append(row)
                for obj in (s, total):obj['eligible_rows'] += 1; obj['eligible_RNA'] += count
        panel['unfiltered_n_molecules'] = panel.get('n_molecules')
        panel['molecules'] = kept
        panel['n_molecules'] = sum(float(m.get('count', 1)) for m in kept)
        panel['bleeding_filter_by_isoform'] = dict(stats)
    result['read_filter'] = dict(enabled=True, kind='splice_edge_extension', tolerance_nt=tolerance,
                                 filtering_precedes_sampling=True, **total)
    result['total_molecules'] = sum(p['n_molecules'] for p in result.get('panels', []))
    return result


def cap_reads_per_isoform(response, cap=200, seed=0):
    """Seeded sampling from each isoform's full eligible pool, AFTER filtering."""
    import numpy as np
    if not isinstance(cap, int) or cap < 1:
        raise ValueError('Positive read cap per isoform required')
    result = deepcopy(response)
    for pi, panel in enumerate(result.get('panels', [])):
        groups = defaultdict(list)
        for m in panel.get('molecules', []):groups[m['final_isoform_id']].append(m)
        selected = []
        for gi, (tid, rows) in enumerate(sorted(groups.items())):
            rows.sort(key=lambda m: (str(m['isoform_id']), str(m.get('detections', []))))
            counts = np.array([int(round(float(m.get('count', 1)))) for m in rows], dtype=np.int64)
            if any(float(m.get('count', 1)) != int(n) or n < 0 for m, n in zip(rows, counts)):
                raise ValueError('Nonnegative integer RNA counts required for display sampling')
            drawn = (np.random.default_rng([seed, pi, gi]).multivariate_hypergeometric(counts, cap)
                     if counts.sum() > cap else counts)
            for row, n in zip(rows, drawn):
                if n:
                    row['count'] = float(n); selected.append(row)
            if tid in panel.get('bleeding_filter_by_isoform', {}):
                panel['bleeding_filter_by_isoform'][tid]['shown_RNA'] = int(drawn.sum())
        panel['molecules'] = selected; panel['n_molecules'] = sum(m['count'] for m in selected)
    result['total_molecules'] = sum(p['n_molecules'] for p in result.get('panels', []))
    result['display_read_cap_per_isoform'] = cap
    result['display_sampling_seed'] = seed
    return result
