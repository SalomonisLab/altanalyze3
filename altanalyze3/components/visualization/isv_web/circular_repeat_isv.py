"""Native trans-splicing layout for two traversals of one proposed circle.

The two reference copies are deliberately separate vertically stacked locus
rows. This is an annotated path hypothesis, never an observed full read.
"""

from __future__ import annotations

from copy import deepcopy

from .isv_web_remote import render


def build_two_copy_path(native_payload, event):
    gm = native_payload['gene_model']
    symbol = native_payload.get('symbol') or native_payload.get('gene') or 'gene'
    parts = event.split('|')
    if len(parts) != 6 or parts[0] != 'BSJ' or parts[2] != parts[3]:
        raise ValueError('Expected same-strand BSJ event')
    chrom, strand = parts[1], parts[2]
    donor, acceptor = int(parts[4]), int(parts[5])
    if (strand == '+' and donor <= acceptor) or (strand == '-' and donor >= acceptor):
        raise ValueError('BSJ order conflicts with annotated strand')
    if gm.get('chrom') != chrom or gm.get('strand') != strand:
        raise ValueError('Native registered gene model conflicts with event')
    lo, hi = sorted((donor, acceptor))
    blocks = []
    for b in gm['blocks']:
        a, z = max(lo, int(b['start'])), min(hi, int(b['end']))
        if a < z:
            blocks.append(dict(start=a, end=z, exon_id=str(b['exon_id']), type='E'))
    if not blocks:
        raise ValueError('No registered exon blocks within proposed BSJ interval')
    ordered = sorted(blocks, key=lambda b: b['start'], reverse=strand == '-')
    genes = [gm['gene'] + '__path_copy_1', gm['gene'] + '__path_copy_2']
    # Distinct drawing chromosomes force separate native locus frames.
    # Original chromosome and coordinates remain in source metadata.
    models = []
    for i, gene in enumerate(genes, 1):
        models.append(dict(gene=gene, symbol='%s · copy %d' % (symbol, i),
                           chrom='%s / copy %d' % (chrom, i), strand=strand,
                           extent=[lo, hi], blocks=deepcopy(blocks)))

    def segments(copy):
        gene = genes[copy - 1]
        return [dict(gene=gene, chrom=models[copy - 1]['chrom'], strand=strand,
                     start=b['start'], end=b['end'], type='E', label=b['exon_id'])
                for b in ordered]

    def mol(name, color, copies, cluster):
        rows = [s for copy in copies for s in segments(copy)]
        return dict(isoform_id=name, final_isoform_id=name,
                    source_id='annotated-path', molecule_id=name, count=1,
                    known=False, cluster_index=cluster, color_index=cluster,
                    display_color=color, exon_segments=rows)

    first = mol('Traversal 1', '#4884BC', (1,), 0)
    second = mol('Traversal 2', '#4884BC', (2,), 1)
    comparator = mol('Linear comparator', '#E9B940', (1,), 2)
    data = dict(gene=genes[0], symbol=symbol, gene_model=models[0],
                gene_models=models, has_multiple_loci=True,
                circular_repeat_view=True,
                panels=[dict(condition='Proposed back-splice path · schematic',
                             molecules=[first, second, comparator], n_molecules=3,
                             repeat_path_link=dict(donor_gene=genes[0],
                                                   acceptor_gene=genes[1],
                                                   donor=donor, acceptor=acceptor,
                                                   donor_row='Traversal 1',
                                                   acceptor_row='Traversal 2'))],
                total_molecules=3)
    data['repeat_path_metadata'] = dict(event=event, original_chrom=chrom,
                                        acceptor_exon=ordered[0]['exon_id'],
                                        donor_exon=ordered[-1]['exon_id'],
                                        copy1_blocks=len(ordered),
                                        copy2_blocks=len(ordered),
                                        stacked_reference_rows=True,
                                        schematic_only=True)
    return data


def render_two_copy_path(native_payload, event, out_path):
    data = build_two_copy_path(deepcopy(native_payload), event)
    caption = ('Annotated path hypothesis using the native trans-splicing layout. '
               'The two vertically stacked reference rows are the same gene in its original genomic coordinates, '
               'shown twice as separate drawing loci and blue traversal rows. The red link runs from the donor '
               'on copy 1 to the acceptor on copy 2; gold shows one linear traversal. '
               'This does not show two observed reads or prove '
               'physical circular topology. See the companion empirical ISV and read-level diagnostics.')
    ok = render(data, out_path, width=1200, intron_scale=0.1, isoform_track=False,
                font_scale=1.2, show_ruler=False, show_read_labels=True,
                show_read_counts=False, dense_exon_labels=True,
                structure_labels_above=True, title_gene=data['symbol'],
                title_event='Same gene on two vertically offset native reference tracks',
                caption=caption, iso_order='given')
    if not ok:
        raise RuntimeError('Native ISV did not render')
    return data['repeat_path_metadata']
