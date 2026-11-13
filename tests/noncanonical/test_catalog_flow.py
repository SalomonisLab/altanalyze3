"""Unit tests for the catalog flow summary stage and its renderer (synthetic inputs)."""
import gzip
from pathlib import Path

import pytest

from altanalyze3.components.long_read.noncanonical import catalog_flow as CF
from altanalyze3.components.visualization import catalog_flow as VF


def _gz(path, lines):
    with gzip.open(path, 'wt') as handle:
        handle.write(''.join(l + '\n' for l in lines))


def test_eligibility_matches_source_scan():
    assert CF.eligible(['ENSG1', '+', 'E1.1|E2.1', '5', 'lib'])
    assert not CF.eligible(['UNK1', '+', 'E1.1', '5', 'lib'])
    assert not CF.eligible(['ENSG1', '+', '', '5', 'lib'])


def test_protein_rule_and_gene_longest(tmp_path):
    p = tmp_path / 'protein_summary.txt.gz'
    _gz(p, ['Gene ID\tTranscript ID\tProtein Length\tNMD Status\tIntron Retention\tLongest Isoform Length',
            'G1\tI1\t400\tNot-NMD\tFalse\t400',
            'G1\tI2\t90\tNot-NMD\tFalse\t400',
            'G1;(NMD)\tI3\t300\tPotential-NMD\tFalse\t400',
            'G2\tI4\t50\tNMD\tFalse\t10'])
    classes, longest = CF.protein_classes(p)
    assert classes[('G1', 'I1')] == 'protein_coding'
    assert classes[('G1', 'I2')] == 'truncated_lt_25pct_longest'
    assert classes[('G1', 'I3')] == 'likely_NMD'
    assert classes[('G2', 'I4')] == 'likely_NMD'
    assert longest == {'G1': 400}


def test_protein_rule_refuses_conflicting_gene_longest(tmp_path):
    p = tmp_path / 'protein_summary.txt.gz'
    _gz(p, ['Gene ID\tTranscript ID\tProtein Length\tNMD Status\tIntron Retention\tLongest Isoform Length',
            'G1\tI1\t400\tNot-NMD\tFalse\t400', 'G1\tI2\t300\tNot-NMD\tFalse\t350'])
    with pytest.raises(ValueError):
        CF.protein_classes(p)


def test_transcript_id_both_attribute_styles():
    gtf = 'gene_id "ENSG1";transcript_id "123.lib";'
    gff3 = 'ID=exon:ENST1:1;Parent=ENST1;gene_id=ENSG2.6;transcript_id=ENST1;exon_number=1'
    assert CF.TRANSCRIPT_ID.search(gtf).group(1) == '123.lib'
    assert CF.TRANSCRIPT_ID.search(gff3).group(1) == 'ENST1'


def test_interchromosomal_exemplar_protein():
    longest = {'GA': 400}
    assert CF.exemplar_protein({'gene_a': 'GA', 'catalog_exemplars': ''}, longest) == 'no_catalog_isoform'
    assert CF.exemplar_protein({'gene_a': 'GA', 'catalog_exemplars': 'x.s1 (42 aa, Potential-NMD) x10'}, longest) == 'likely_NMD'
    assert CF.exemplar_protein({'gene_a': 'GA', 'catalog_exemplars': 'x.s1 (90 aa, Not-NMD) x6; y (300 aa, Not-NMD) x2'}, longest) == 'truncated_lt_25pct_longest'
    assert CF.exemplar_protein({'gene_a': 'GA', 'catalog_exemplars': 'x.s1 (113 aa, Not-NMD) x6'}, {'GA': 300}) == 'protein_coding'


def _trace_inputs(tmp_path):
    ta = tmp_path / 'ta.txt.gz'
    _gz(ta, ['G1\t+\tE1.1|E2.1\t1\ts', 'G1\t+\tE1.1|E2.1\t2\ts', 'G1\t+\tE1.1|E3.1\t3\ts',
             'G2\t+\tE1.1|G3:E2.1\t4\ts', 'UNK1\t+\tE1.1\t5\ts'])
    mol = tmp_path / 's.mol2struct.tsv'
    mol.write_text('molecule_id\tgene\tstructure\n1\tG1\tc1\n2\tG1\tc1\n3\tG1\tc2\n4\tG2\tc3\n')
    fmap = {('G1', 'c1'): 'iA', ('G1', 'c2'): 'iA', ('G2', 'c3'): 'iTS'}
    pmap = {('G1', 'c1'): 'pA', ('G1', 'c2'): 'pB', ('G2', 'c3'): 'pTS'}
    return ta, mol, fmap, pmap


def test_trace_source_partitions_linear_and_chimeric(tmp_path):
    ta, mol, fmap, pmap = _trace_inputs(tmp_path)
    res = CF.trace_source('s', ta, mol, fmap, pmap, ts_isoforms={('G2', 'iTS')}, chimeric_split={'3'})
    assert res['eligible'] == 4                      # the UNK row is not eligible
    assert res['fin'][('G1', 'iA')] == [3, 2]        # read 3 is a split-alignment read
    assert res['fin'][('G2', 'iTS')] == [1, 0]       # trans-splicing isoform reads are chimeric
    assert res['pre'][('G1', 'pB')] == [1, 0]        # pB loses its only read
    assert res['ts_reads'] == 1 and res['ts_and_split'] == 0
    assert res['structures'] == 3 and res['linear_structures'] == 1 and res['intra_structures'] == 1


def test_trace_source_refuses_missing_prefilter_class(tmp_path):
    ta, mol, fmap, pmap = _trace_inputs(tmp_path)
    del pmap[('G1', 'c2')]
    with pytest.raises(ValueError):
        CF.trace_source('s', ta, mol, fmap, pmap, ts_isoforms=set(), chimeric_split=set())


def test_trace_source_refuses_read_count_mismatch(tmp_path):
    ta, mol, fmap, pmap = _trace_inputs(tmp_path)
    mol.write_text('molecule_id\tgene\tstructure\n1\tG1\tc1\n2\tG1\tc1\n3\tG1\tc2\n')
    with pytest.raises(ValueError):
        CF.trace_source('s', ta, mol, fmap, pmap, ts_isoforms=set(), chimeric_split=set())


VALUES = dict(umi_label='UMI in n=2', umi_total=100, linear_structures=40, linear_collapsed=20, known=5, novel=7,
              known_protein_coding=3, novel_protein_coding=4, chimeric_structures=9, chimeric_collapsed=6,
              intra_trans_splice=2, inter_trans_splice=1, circrna=3, intra_trans_splice_protein_coding=1,
              inter_trans_splice_protein_coding=0)


def test_renderer_draws_every_value(tmp_path):
    result = VF.render(dict(VALUES), tmp_path / 'flow')
    assert Path(result['pdf']).stat().st_size > 0 and Path(result['png']).stat().st_size > 0
    for key in VF.REQUIRED:
        assert f'{VALUES[key]:,}' in result['drawn_text']
    assert 'UMI in n=2' in result['drawn_text']
    pdf = Path(result['pdf']).read_bytes()
    assert b'/Subtype /Image' not in pdf


def test_renderer_range_and_refusals(tmp_path):
    v = dict(VALUES); del v['linear_collapsed']; v['linear_collapsed_low'] = 18; v['linear_collapsed_high'] = 20
    assert VF.collapsed_label(v) == '18–20'
    with pytest.raises(ValueError):
        VF.render({k: x for k, x in VALUES.items() if k != 'circrna'}, tmp_path / 'x')
    bad = dict(VALUES); bad['known'] = 'many'
    with pytest.raises(ValueError):
        VF.render(bad, tmp_path / 'y')


def test_values_table_roundtrip(tmp_path):
    path = tmp_path / 'values.tsv'
    path.write_text('box\tvalue\n' + ''.join(f'{k}\t{v}\n' for k, v in VALUES.items()))
    assert VF.read_values(path) == VALUES
    path.write_text('box\tvalue\nknown\t1\nknown\t2\n')
    with pytest.raises(ValueError):
        VF.read_values(path)
