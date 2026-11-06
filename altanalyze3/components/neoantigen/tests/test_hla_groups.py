import pytest
from altanalyze3.components.neoantigen.io import read_table
from altanalyze3.components.neoantigen.hla import combine_hla_groups, read_hla_table


def _hla(tmp_path, name, sample, alleles):
    path = tmp_path / f'{name}.tsv'
    path.write_text(f'sample\thla\n{sample}\t{alleles}\n')
    return path


def _groups(tmp_path, pairs):
    path = tmp_path / 'groups.tsv'
    path.write_text('sample\tgroup\n' + ''.join(f'{s}\t{g}\n' for s, g in pairs))
    return path


def test_concordant_group_called_conflict_dropped_and_reported(tmp_path):
    inputs = [
        _hla(tmp_path, 'a1', 'D1_pre', 'HLA-A*02:01,HLA-A*24:02,HLA-B*07:02,HLA-B*44:02,HLA-C*07:02,HLA-C*05:01'),
        _hla(tmp_path, 'a2', 'D1_post', 'HLA-A*24:02,HLA-A*02:01,HLA-B*44:02,HLA-B*07:02'),
        _hla(tmp_path, 'b1', 'D2_pre', 'HLA-A*01:01,HLA-B*08:01,HLA-C*07:01'),
        _hla(tmp_path, 'b2', 'D2_post', 'HLA-A*01:01,HLA-B*08:02,HLA-C*07:01'),
        _hla(tmp_path, 'c1', 'D3', 'HLA-A*03:01,HLA-B*35:01,HLA-C*04:01'),
    ]
    groups = _groups(tmp_path, [('D1_pre', 'D1'), ('D1_post', 'D1'), ('D1_mouse', 'D1'),
                                ('D2_pre', 'D2'), ('D2_post', 'D2'), ('D3', 'D3')])
    out, report = tmp_path / 'hla_groups.tsv', tmp_path / 'report.tsv'
    called = combine_hla_groups(inputs, groups, out, report)
    table = read_hla_table(out)
    assert sorted(table) == ['D1', 'D3']
    assert sorted(table['D1']) == sorted(['HLA-A*02:01', 'HLA-A*24:02', 'HLA-B*07:02', 'HLA-B*44:02',
                                          'HLA-C*07:02', 'HLA-C*05:01'])
    assert called['D3'] == ['HLA-A*03:01', 'HLA-B*35:01', 'HLA-C*04:01']
    rows, _ = read_table(report)
    status = {(r['group'], r['gene']): r for r in rows}
    assert len(rows) == 9
    assert status[('D2', 'HLA-B')]['status'] == 'conflict'
    assert status[('D2', 'HLA-A')]['status'] == 'agreed'
    assert status[('D1', 'HLA-C')]['status'] == 'agreed'
    assert status[('D1', 'HLA-C')]['n_samples_called'] == '1'
    assert status[('D1', 'HLA-A')]['not_typed'] == 'D1_mouse'
    assert status[('D1', 'HLA-A')]['n_samples_typed'] == '2'


def test_homozygous_vs_heterozygous_is_a_conflict(tmp_path):
    inputs = [_hla(tmp_path, 'x1', 'S1', 'HLA-A*02:01,HLA-B*07:02'),
              _hla(tmp_path, 'x2', 'S2', 'HLA-A*02:01,HLA-A*01:01,HLA-B*07:02'),
              _hla(tmp_path, 'y1', 'S3', 'HLA-A*11:01')]
    groups = _groups(tmp_path, [('S1', 'D'), ('S2', 'D'), ('S3', 'E')])
    combine_hla_groups(inputs, groups, tmp_path / 'o.tsv', tmp_path / 'r.tsv')
    assert sorted(read_hla_table(tmp_path / 'o.tsv')) == ['E']


def test_refuses_unmapped_sample_and_all_conflicting(tmp_path):
    a = _hla(tmp_path, 'a', 'S1', 'HLA-A*02:01')
    b = _hla(tmp_path, 'b', 'S2', 'HLA-A*01:01')
    with pytest.raises(ValueError, match='has no row'):
        combine_hla_groups([a, b], _groups(tmp_path, [('S1', 'D')]), tmp_path / 'o.tsv', tmp_path / 'r.tsv')
    with pytest.raises(ValueError, match='No group has a concordant'):
        combine_hla_groups([a, b], _groups(tmp_path, [('S1', 'D'), ('S2', 'D')]),
                           tmp_path / 'o.tsv', tmp_path / 'r.tsv')
    assert (tmp_path / 'r.tsv').exists()
