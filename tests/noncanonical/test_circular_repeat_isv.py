import pytest

from altanalyze3.components.visualization.isv_web.circular_repeat_isv import build_two_copy_path


@pytest.mark.parametrize('strand,donor,acceptor,expected', [
    ('+', 399, 100, ['E2', 'E3', 'E4']),
    ('-', 100, 399, ['E4', 'E3', 'E2']),
])
def test_two_copy_path_respects_transcript_orientation(strand, donor, acceptor, expected):
    native = {'symbol':'GENE','gene_model':{'gene':'ENSG_TEST','chrom':'chr1',
                                          'strand':strand,'blocks':[
        {'exon_id':'E2','start':100,'end':150},
        {'exon_id':'E3','start':200,'end':250},
        {'exon_id':'E4','start':350,'end':399}]}}
    event='BSJ|chr1|{0}|{0}|{1}|{2}'.format(strand,donor,acceptor)
    data=build_two_copy_path(native,event)
    first=data['panels'][0]['molecules'][0]['exon_segments']
    second=data['panels'][0]['molecules'][1]['exon_segments']
    assert [x['label'] for x in first]==expected
    assert [x['label'] for x in second]==expected
    assert data['has_multiple_loci'] is True
    assert len(data['gene_models']) == 2
    assert data['gene_models'][0]['gene'] != data['gene_models'][1]['gene']
    assert data['gene_models'][0]['chrom'] != data['gene_models'][1]['chrom']
    assert data['gene_models'][0]['blocks'] == data['gene_models'][1]['blocks']
    assert len({x['gene'] for x in first}) == 1
    assert len({x['gene'] for x in second}) == 1
    assert first[0]['gene'] != second[0]['gene']
    assert data['panels'][0]['repeat_path_link']['donor_row'] == 'Traversal 1'
    assert data['repeat_path_metadata']['schematic_only'] is True
    assert data['repeat_path_metadata']['stacked_reference_rows'] is True
    assert len(data['panels'][0]['molecules'][2]['exon_segments'])==3


def test_linear_order_rejected():
    native={'gene_model':{'gene':'ENSG_TEST','chrom':'chr1','strand':'+','blocks':[
        {'exon_id':'E2','start':100,'end':150}]}}
    with pytest.raises(ValueError,match='order'):
        build_two_copy_path(native,'BSJ|chr1|+|+|100|399')
