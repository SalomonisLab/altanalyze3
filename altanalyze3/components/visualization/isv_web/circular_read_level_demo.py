"""Show actual read-pair alignments in the native two-locus ISV layout.

Each blue row is one source query. Its query-ordered alignment blocks are
assigned to two drawing copies of the same registered gene, with a red join
at the measured donor/acceptor. No physical circle is inferred by this view.
"""

from copy import deepcopy

from .circular_repeat_isv import build_two_copy_path
from .isv_web_remote import render


def _exons(text, model, chrom):
    result = []
    for token in text.split(','):
        a, b = (int(x) for x in token.split('-'))
        if a > b:
            raise ValueError('Source exon interval is reversed')
        result.append(dict(gene=model['gene'], chrom=model['chrom'],
                           strand=model['strand'], start=a, end=b,
                           type='E', label=None))
    return result


def _admit(e, donor, acceptor):
    try:
        if e['joint_type'] != 'back' or e['status'] != 'kept':
            return False
        if e['query_gap'] != '0' or int(e['outside_circle_nt']) != 0:
            return False
        if int(e['a_mapq']) < 50 or int(e['b_mapq']) < 50:
            return False
        if int(e['a_qstart']) != 0 or int(e['a_qend']) != int(e['b_qstart']):
            return False
        if int(e['b_qend']) < int(e['read_length']) - 2:
            return False
        if int(e['a_qend']) < 20 or int(e['b_qend']) - int(e['b_qstart']) < 20:
            return False
        # Source exon intervals are genomic ascending, while a minus-strand
        # query traverses them in the reverse direction.
        minus = e.get('a_strand', '-' if donor < acceptor else '+') == '-'
        a = e['a_exons'].split(',')[0 if minus else -1]
        b = e['b_exons'].split(',')[-1 if minus else 0]
        a_boundary = int(a.split('-')[0 if minus else 1])
        b_boundary = int(b.split('-')[1 if minus else 0])
        if a_boundary != donor or b_boundary != acceptor:
            return False
        return True
    except (KeyError, ValueError, IndexError):
        return False


def build_read_level_demo(native_payload, event, per_condition=2):
    if per_condition < 1:
        raise ValueError('per_condition must be positive')
    seed = build_two_copy_path(native_payload, event)
    first, second = seed['gene_models']
    donor, acceptor = (int(x) for x in event.split('|')[4:6])
    lo, hi = min(donor, acceptor), max(donor, acceptor)
    panels, ledger, comparator_ledger = [], [], []
    for source_panel in native_payload['panels']:
        blue, gold = [], []
        for original in source_panel['molecules']:
            if original.get('display_color') == '#4884BC':
                for e in original.get('circular_evidence', []):
                    if _admit(e, donor, acceptor):
                        blue.append((original, e))
            elif original.get('display_color') == '#E9B940':
                gold.append(original)
        # Deterministic choice for a compact demo, not an inference filter.
        blue.sort(key=lambda item: (-int(item[1]['b_qend']), -int(item[1]['a_qend']),
                                    item[1]['library'], item[1]['query_name']))
        molecules = []
        for index, (original, e) in enumerate(blue[:per_condition]):
            a = _exons(e['a_exons'], first, first['chrom'])
            b = _exons(e['b_exons'], second, second['chrom'])
            minus = e.get('a_strand', '-' if donor < acceptor else '+') == '-'
            if ((a[0]['start'] if minus else a[-1]['end']) != donor or
                    (b[-1]['end'] if minus else b[0]['start']) != acceptor):
                raise ValueError('Selected read does not touch both BSJ endpoints')
            # Keep source labels legible in a dense native pileup. The full
            # library, query ID and query intervals remain in the audit ledger.
            name = '%s:%s' % (e['library'][:10],
                              e['query_name'].rsplit('/', 1)[-1][-8:])
            molecules.append(dict(isoform_id=name, final_isoform_id=name,
                                  source_id=e['library'], molecule_id=e['query_name'],
                                  count=1, cluster_index=index, color_index=index,
                                  display_color='#4884BC', known=False,
                                  use_track_color_label=True,
                                  exon_segments=a + b,
                                  read_bsj_link=dict(event=event, pair_index=e['pair_index'],
                                                     a_qstart=e['a_qstart'], a_qend=e['a_qend'],
                                                     b_qstart=e['b_qstart'], b_qend=e['b_qend'])))
            ledger.append(dict(condition=source_panel['condition'], library=e['library'],
                               query_name=e['query_name'], read_length=e['read_length'],
                               a_qstart=e['a_qstart'], a_qend=e['a_qend'],
                               b_qstart=e['b_qstart'], b_qend=e['b_qend'],
                               a_exons=e['a_exons'], b_exons=e['b_exons'],
                               donor=e['donor'], acceptor=e['acceptor'],
                               source_pair_index=e['pair_index'],
                               source_status=e['status'],
                               a_cigar=e['a_cigar'], b_cigar=e['b_cigar'],
                               a_mapq=e['a_mapq'], b_mapq=e['b_mapq'],
                               query_gap=e['query_gap'],
                               outside_circle_nt=e['outside_circle_nt']))
        def overlap(m):
            return sum(max(0, min(hi, int(s['end'])) - max(lo, int(s['start'])) + 1)
                       for s in m['exon_segments'] if s['type'] == 'E')
        # Retain multiple actual BSJ-negative source reads. A single maximal-
        # overlap comparator falsely suggests that the linear evidence is one
        # molecule or one structure when the frozen source panel has many.
        gold.sort(key=lambda m: (-overlap(m), str(m.get('source_id')),
                                 str(m.get('molecule_id'))))
        for gold_index, original in enumerate(gold[:per_condition]):
            rows = []
            for s in original['exon_segments']:
                a, b = max(lo, int(s['start'])), min(hi, int(s['end']))
                if a <= b and s['type'] == 'E':
                    rows.append(dict(gene=first['gene'], chrom=first['chrom'],
                                     strand=first['strand'], start=a, end=b,
                                     type='E', label=s.get('label')))
            if rows:
                name = 'linear %s:%s' % (str(original.get('source_id', ''))[:10],
                                          str(original.get('molecule_id', ''))[-8:])
                molecules.append(dict(isoform_id=name, final_isoform_id=name,
                                      source_id=original.get('source_id'),
                                      molecule_id=original.get('molecule_id'),
                                      count=1, cluster_index=per_condition + gold_index,
                                      color_index=per_condition + gold_index,
                                      display_color='#E9B940',
                                      known=False, use_track_color_label=True,
                                      exon_segments=rows))
                comparator_ledger.append(dict(
                    condition=source_panel['condition'],
                    library=original.get('source_id'),
                    molecule_id=original.get('molecule_id'),
                    original_exon_segments=[
                        dict(start=s['start'], end=s['end'], type=s['type'])
                        for s in original['exon_segments']],
                    displayed_exon_segments=[
                        dict(start=s['start'], end=s['end'], type=s['type'])
                        for s in rows]))
        if molecules:
            condition = source_panel['condition'].split(' · ')[0] + ' · selected source reads'
            panels.append(dict(condition=condition, molecules=molecules,
                               n_molecules=len(molecules)))
    if not ledger:
        raise ValueError('No exact source pair meets the displayed junction-geometry rule')
    seed['panels'] = panels
    seed['total_molecules'] = sum(p['n_molecules'] for p in panels)
    seed['read_level_demo'] = dict(event=event, source_blue_pairs=len(ledger),
                                   source_gold_reads=len(comparator_ledger),
                                   source_conditions=len(panels),
                                   selection=('zero query gap, >=50 MAPQ on both arms, '
                                              '<=2 unaligned terminal bases, no outside-circle alignment; '
                                              'up to %d displayed source queries per condition' % per_condition))
    seed['comparator_ledger'] = comparator_ledger
    return seed, ledger


def render_read_level_demo(native_payload, event, out_path, per_condition=2):
    data, ledger = build_read_level_demo(deepcopy(native_payload), event,
                                         per_condition=per_condition)
    caption = ('Each blue row is one actual source query, with the first query interval on copy 1 and '
               'the next query interval on copy 2; the red line marks its measured back-ordered '
               'alignment junction. The two reference rows duplicate the same annotated gene with '
               'independent native locus offsets. Each gold row, where available, '
               'is a separate BSJ-negative source read windowed to this locus. '
               'Query-interval continuity and source identities are in the '
               'companion ledger. These alignments do not establish a physical RNA circle; '
               'sequence-equivalent tandem and other alternative models remain under review.')
    # The native renderer reserves 10 px per unique read label even when its
    # dense pileup is capped at MAXPANEL. At 250+250 rows this separates labels
    # from their tracks and produces a >20,000-point page. Keep every source
    # row and its full identity in the companion ledger, suppress only labels.
    ok = render(data, out_path, width=1500, intron_scale=0.1, isoform_track=False,
                font_scale=1.1, show_ruler=False,
                show_read_labels=per_condition <= 50,
                show_read_counts=False, dense_exon_labels=True,
                structure_labels_above=True, label_protein_length=False,
                title_gene=data['symbol'],
                title_event='Read-level back-splice alignment across two vertically offset copies',
                caption=caption, iso_order='given')
    if not ok:
        raise RuntimeError('Native ISV did not render')
    return data, ledger
