import csv
import hashlib
import json

import pytest

from altanalyze3.components.visualization.psi_boxplot import donor_psi_records, native_psi_records, main


def source(tmp_path):
    path = tmp_path / 'counts.tsv'
    rows = [('healthy', 'H', 'I', 0), ('healthy', 'H', 'S', 0),
            ('MDS', 'D1', 'I', 6), ('MDS', 'D1', 'S', 2),
            ('MDS', 'D2', 'I', 1), ('MDS', 'D2', 'S', 1)]
    with path.open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['cohort', 'donor', 'isoform_id', 'assigned_RNA'])
        writer.writerows(rows)
    roster = {'healthy': ['H'], 'MDS': ['D1', 'D2']}
    roster_path = tmp_path / 'roster.json'
    roster_path.write_text(json.dumps(roster))
    return path, roster, roster_path


def test_zero_support_is_undefined_and_every_roster_identity_is_preserved(tmp_path):
    path, roster, _ = source(tmp_path)
    records = donor_psi_records(path, roster, ['healthy', 'MDS'], ['I'], ['S'])
    assert [(r['condition'], r['donor'], r['PSI']) for r in records] == [
        ('healthy', 'H', None), ('MDS', 'D1', .75), ('MDS', 'D2', .5)]
    assert records[0]['status'] == 'UNDEFINED_ZERO_INFORMATIVE_READS_RETAINED'


@pytest.mark.parametrize('failure', ['missing', 'duplicate', 'unknown', 'nan', 'negative'])
def test_bad_count_records_fail_without_reducing_the_roster(tmp_path, failure):
    path, roster, _ = source(tmp_path)
    text = path.read_text()
    if failure == 'missing':
        text = text.replace('MDS\tD2\tS\t1\n', '')
    elif failure == 'duplicate':
        text += 'MDS\tD2\tS\t1\n'
    elif failure == 'unknown':
        text += 'MDS\tUNEXPECTED\tS\t1\n'
    else:
        text = text.replace('MDS\tD2\tS\t1\n', 'MDS\tD2\tS\t'+('nan' if failure == 'nan' else '-1')+'\n')
    path.write_text(text)
    with pytest.raises(ValueError):
        donor_psi_records(path, roster, ['healthy', 'MDS'], ['I'], ['S'])


def test_all_conditions_and_disjoint_isoform_roles_are_required(tmp_path):
    path, roster, _ = source(tmp_path)
    with pytest.raises(ValueError, match='complete independent'):
        donor_psi_records(path, roster, ['MDS'], ['I'], ['S'])
    with pytest.raises(ValueError, match='disjoint'):
        donor_psi_records(path, roster, ['healthy', 'MDS'], ['I'], ['I'])


def test_plot_cli_retains_zero_support_and_records_exact_source_hash(tmp_path):
    path, _, roster_path = source(tmp_path)
    before = path.read_bytes()
    out = tmp_path / 'psi.pdf'
    argv = ['--counts', str(path), '--donor-roster', str(roster_path),
            '--inclusion-isoforms', 'I', '--exclusion-isoforms', 'S',
            '--cohort-order', 'healthy', 'MDS', '--y-min', '.5', '--out', str(out)]
    result = main(argv)
    assert out.read_bytes().startswith(b'%PDF-') and path.read_bytes() == before
    assert len(result['records']) == 3 and result['records'][0]['PSI'] is None
    assert result['counts_sha256'] == hashlib.sha256(before).hexdigest()
    assert result['all_defined_points_in_view'] and result['y_limits'] == [.5, 1.03]
    with pytest.raises(FileExistsError):
        main(argv)


def test_cached_isv_cli_forwards_caption_and_checks_input_hash(tmp_path, monkeypatch):
    from altanalyze3.components.visualization.isv_web import isv_web_remote as renderer
    response = tmp_path / 'input.json'
    response.write_text(json.dumps({'gene': 'G', 'symbol': 'G', 'panels': []}))
    caption = tmp_path / 'caption.txt'
    caption.write_text('Pinned scientific caption.')
    calls = []
    monkeypatch.setattr(renderer, 'render', lambda *a, **kw: calls.append(kw) or True)
    argv = ['--run_dir', str(tmp_path), '--gene', 'G', '--render-input', str(response),
            '--caption-file', str(caption), '--render-input-sha256', hashlib.sha256(response.read_bytes()).hexdigest()]
    with pytest.raises(SystemExit) as exc:
        renderer.main(argv)
    assert exc.value.code == 0 and calls[0]['caption'] == caption.read_text()
    with pytest.raises(SystemExit, match='SHA256 mismatch'):
        renderer.main(argv[:-1]+['0'*64])
    assert len(calls) == 1


def test_root_visualization_commands_are_registered():
    from altanalyze3.utilities.parser import ArgsParser
    parser = ArgsParser.__new__(ArgsParser).get_parser()
    args = parser.parse_args(['visualize', 'isv-domains', '--run_dir', '/tmp', '--gene', 'G'])
    assert args.func.__name__ == 'run_isv'
    args = parser.parse_args(['visualize', 'psi-boxplot', '--counts', 'c.tsv', '--donor-roster', 'r.json',
                             '--inclusion-isoforms', 'I', '--exclusion-isoforms', 'S',
                             '--cohort-order', 'A', '--out', 'plot.pdf'])
    assert args.func.__name__ == 'run_psi'
    with pytest.raises(SystemExit) as exc:
        ArgsParser(['visualize', 'isv-domains', '--run_dir', '/tmp', '--gene', 'G', '--typo-option'])
    assert exc.value.code == 2


def test_reference_only_tests_keep_subset_comparisons_out_of_family():
    import numpy as np
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests
    from altanalyze3.components.visualization.violin import _pairwise_group_tests
    groups = [dict(label='young', values=[.95,.9,.98]), dict(label='aged', values=[.8,.85,.82,.78,.81]),
              dict(label='all_MDS', values=[.6,.65,.7,.55,.62,.64,.58]),
              dict(label='SRSF2_subset', values=[.6,.65,.7])]
    rows = _pairwise_group_tests(groups, test='mannwhitney', correction='holm', reference='aged')
    assert len(rows)==3 and all('aged' in (r['group_1'],r['group_2']) for r in rows)
    p = [mannwhitneyu(groups[r['i']]['values'], groups[r['j']]['values'], alternative='two-sided').pvalue
         for r in rows]
    assert np.allclose([r['p'] for r in rows],p)
    assert np.allclose([r['p_adjusted'] for r in rows],multipletests(p,method='holm')[1])


def test_reference_with_zero_support_cannot_silently_be_tested(tmp_path):
    path, _, roster_path = source(tmp_path)
    with pytest.raises(ValueError,match='no defined measurements'):
        main(['--counts',str(path),'--donor-roster',str(roster_path),
              '--inclusion-isoforms','I','--exclusion-isoforms','S',
              '--cohort-order','healthy','MDS','--reference-condition','healthy','--out',str(tmp_path/'s.pdf')])
    assert not (tmp_path/'s.pdf').exists()


def test_native_psi_input_keeps_unmeasured_donor_and_literal_event(tmp_path):
    path=tmp_path/'psi.tsv'
    path.write_text('cohort\tdonor\tevent_uid\tPSI\naged\tA\tG:E11-E12|G:E10-E12\t0.8\n'
                    'MDS\tD\tG:E11-E12|G:E10-E12\t0.6\nMDS\tD0\tG:E11-E12|G:E10-E12\t\n')
    roster={'aged':['A'],'MDS':['D','D0']};event='G:E11-E12|G:E10-E12'
    rows=native_psi_records(path,roster,['aged','MDS'],event)
    assert [r['PSI'] for r in rows]==[.8,.6,None]
    with pytest.raises(ValueError,match='identity mismatch'):
        native_psi_records(path,roster,['aged','MDS'],'G:E10-E12|G:E11-E12')
    with pytest.raises(ValueError,match='complete donor roster'):
        native_psi_records(path,{'aged':['A'],'MDS':['D','D0','missing']},['aged','MDS'],event)


def test_existing_native_statistics_are_annotated_without_retesting(tmp_path):
    path=tmp_path/'psi.tsv';event='G:E11-E12|G:E10-E12'
    path.write_text('cohort\tdonor\tevent_uid\tPSI\naged\tA\t'+event+'\t0.8\nMDS\tD\t'+event+'\t0.6\n')
    roster=tmp_path/'roster.json';roster.write_text(json.dumps({'aged':['A'],'MDS':['D']}))
    stats=tmp_path/'native.tsv';stats.write_text('reference\tcomparison\tp\ttest\tcorrection\tnative_FDR\naged\tMDS\t0.012\tnative moderated test\traw P\t1.0\n')
    result=main(['--psi-table',str(path),'--event-uid',event,'--donor-roster',str(roster),
                 '--cohort-order','aged','MDS','--annotations-table',str(stats),'--out',str(tmp_path/'a.pdf')])
    assert result['statistics'][0]['p']==.012
    assert result['statistics'][0]['label']=='*' and result['statistics'][0]['native_FDR']=='1.0'
    assert result['existing_statistics_sha256']==hashlib.sha256(stats.read_bytes()).hexdigest()


def test_saved_pseudobulk_display_tests_each_donor_once(tmp_path):
    path=tmp_path/'native.tsv';event='G:E11-E12|G:E10-E12'
    rows=[('aged','A1.s1','A1',.8),('aged','A1.s2','A1',.9),
          ('aged','A2.s1','A2',.7),('aged','A2.s2','A2',.8),
          ('MDS','B1.s1','B1',.4),('MDS','B1.s2','B1',.5),
          ('MDS','B2.s1','B2',.5),('MDS','B2.s2','B2',.6),
          ('MDS','B3.s1','B3',None)]
    with path.open('w') as h:
        w=csv.writer(h,delimiter='\t');w.writerow(['cohort','pseudobulk','donor','event_uid','PSI'])
        w.writerows((c,p,d,event,v) for c,p,d,v in rows)
    roster={c:[p for group,p,d,v in rows if group==c] for c in ['aged','MDS']}
    rp=tmp_path/'roster.json';rp.write_text(json.dumps(roster))
    result=main(['--psi-table',str(path),'--event-uid',event,'--donor-roster',str(rp),
                 '--donor-column','pseudobulk','--stat-donor-column','donor',
                 '--cohort-order','aged','MDS','--reference-condition','aged',
                 '--unit-label','pseudobulks','--show-clipped-points','--y-min','.5','--out',str(tmp_path/'p.pdf')])
    assert len(result['records'])==9
    assert result['statistics'][0]['n_1']==result['statistics'][0]['n_2']==2
    assert result['records'][4]['PSI']==.4 and not result['all_defined_points_in_view']
    assert result['clipped_points_marked']
    with (tmp_path/'p.statistical_units.tsv').open() as h: units=list(csv.DictReader(h,delimiter='\t'))
    assert len(units)==5 and units[-1]['biological_donor']=='B3'
    assert units[-1]['mean_saved_PSI']=='' and units[-1]['status']=='NO_DEFINED_NATIVE_PSI_RETAINED'
