from altanalyze3.components.visualization.isv_web.circular_read_level_demo import build_read_level_demo


def test_source_query_crosses_two_separate_native_reference_copies():
    event = 'BSJ|chr1|+|+|249|100'
    evidence = dict(joint_type='back', status='kept', query_gap='0',
                    outside_circle_nt='0', a_mapq='60', b_mapq='60',
                    a_qstart='0', a_qend='101', b_qstart='101', b_qend='151',
                    read_length='151', a_exons='100-150,200-249',
                    b_exons='100-149', library='lib1',
                    query_name='molecule/42', pair_index='1',
                    donor='249', acceptor='100', a_cigar='51=49N50=50S',
                    b_cigar='101S50=')
    native = dict(symbol='GENE', gene_model=dict(gene='ENSG_TEST',
                  chrom='chr1', strand='+', extent=[100, 249],
                  blocks=[dict(exon_id='E2', start=100, end=150),
                          dict(exon_id='E3', start=200, end=249)]),
                  panels=[dict(condition='young · source', molecules=[
                      dict(display_color='#4884BC', circular_evidence=[evidence]),
                      dict(display_color='#E9B940', source_id='lib1',
                           molecule_id='linear1', isoform_id='linear1',
                           exon_segments=[dict(start=100, end=150, type='E', label='E2')]),
                      dict(display_color='#E9B940', source_id='lib2',
                           molecule_id='linear2', isoform_id='linear2',
                           exon_segments=[dict(start=110, end=149, type='E', label='E2')])
                  ])])
    result, ledger = build_read_level_demo(native, event)
    assert len(ledger) == 1
    assert result['has_multiple_loci'] is True
    assert len(result['gene_models']) == 2
    row = result['panels'][0]['molecules'][0]
    assert row['read_bsj_link']['a_qend'] == row['read_bsj_link']['b_qstart']
    assert [s['gene'] for s in row['exon_segments']] == [
        result['gene_models'][0]['gene']] * 2 + [result['gene_models'][1]['gene']]
    gold = [m for m in result['panels'][0]['molecules']
            if m['display_color'] == '#E9B940']
    assert len(gold) == 2
    assert {m['molecule_id'] for m in gold} == {'linear1', 'linear2'}
    assert len(result['comparator_ledger']) == 2


def test_minus_strand_uses_query_direction_endpoint_coordinates():
    event = 'BSJ|chr4|-|-|100|249'
    evidence = dict(joint_type='back', status='kept', query_gap='0',
                    outside_circle_nt='0', a_mapq='60', b_mapq='60',
                    a_strand='-', a_qstart='0', a_qend='101',
                    b_qstart='101', b_qend='151', read_length='151',
                    a_exons='100-150,200-249', b_exons='200-249',
                    library='lib1', query_name='molecule/43',
                    pair_index='1', donor='100', acceptor='249',
                    a_cigar='51=49N50=50S', b_cigar='101S50=')
    native = dict(symbol='GENE', gene_model=dict(gene='ENSG_TEST',
                  chrom='chr4', strand='-', extent=[100, 249],
                  blocks=[dict(exon_id='E3', start=100, end=150),
                          dict(exon_id='E2', start=200, end=249)]),
                  panels=[dict(condition='young · source',
                               molecules=[dict(display_color='#4884BC',
                                               circular_evidence=[evidence])])])
    result, ledger = build_read_level_demo(native, event)
    assert len(ledger) == 1
    row = result['panels'][0]['molecules'][0]
    assert row['exon_segments'][0]['start'] == 100
    assert row['exon_segments'][-1]['end'] == 249
    assert len({s['gene'] for s in row['exon_segments']}) == 2
