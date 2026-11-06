import hashlib,json
from pathlib import Path
import pytest
from pypdf import PdfReader
from altanalyze3.components.visualization.psi_boxplot_grid import main


def manifest(tmp_path):
    roster=tmp_path/'roster.json';roster.write_text(json.dumps({'healthy':['Y','A'],'MDS-other':['O'],'MDS-SRSF2':['S']}))
    markers=tmp_path/'markers.tsv';markers.write_text('cohort\tdonor\tmarker\tlabel\nhealthy\tY\t^\tYoung\nhealthy\tA\to\tAged / MDS\nMDS-other\tO\to\tAged / MDS\nMDS-SRSF2\tS\to\tAged / MDS\n')
    panels=[]
    for gene in ['G1','G2']:
        event=gene+':E1-E3|'+gene+':E2-E3';path=tmp_path/(gene+'.tsv')
        path.write_text('cohort\tdonor\tevent_uid\tPSI\n'+'\n'.join(c+'\t'+d+'\t'+event+'\t'+v for c,d,v in [('healthy','Y','0'),('healthy','A',''),('MDS-other','O','.1'),('MDS-SRSF2','S','.7')])+'\n')
        panels.append(dict(gene=gene,event_uid=event,event_label=gene+':E1-E3|E2-E3',native_argv=['--psi-table',str(path),'--psi-sha256',hashlib.sha256(path.read_bytes()).hexdigest(),'--event-uid',event,'--donor-roster',str(roster),'--cohort-order','healthy','MDS-other','MDS-SRSF2','--marker-table',str(markers),'--out',str(tmp_path/(gene+'_never_written.pdf'))]))
    p=tmp_path/'panels.json';p.write_text(json.dumps(dict(required_genes=['G1','G2'],panels=panels,title='Native MPP-1 PSI',caption='Original values; missing remains missing.')))
    return p


def test_one_page_native_panels_retain_every_value_and_full_event_header(tmp_path):
    p=manifest(tmp_path);out=tmp_path/'all.pdf';r=main(['--panel-manifest',str(p),'--columns','2','--width','12','--height','8','--out',str(out)])
    assert len(PdfReader(out).pages)==1
    text=PdfReader(out).pages[0].extract_text().replace('\n','')
    assert 'G1:E1-E3|E2-E3' in text and 'G2:E1-E3|E2-E3' in text
    assert len(r['panels'])==2 and all([x['PSI'] for x in q['records']]==[0,None,.1,.7] for q in r['panels'])
    assert all(q['records'][0]['point_marker']=='^' for q in r['panels'])
    assert not list(tmp_path.glob('*never_written*'))


def test_composite_rejects_incorrect_event_header_before_writing_pdf(tmp_path):
    p=manifest(tmp_path);spec=json.loads(p.read_text());spec['panels'][1]['event_label']='G2:E2-E3|E1-E3';p.write_text(json.dumps(spec))
    with pytest.raises(ValueError,match='header'):
        main(['--panel-manifest',str(p),'--columns','2','--width','12','--height','8','--out',str(tmp_path/'bad.pdf')])
    assert not (tmp_path/'bad.pdf').exists()


def test_axis_ranges_round_outward_and_keep_zero_missing_and_flat_values():
    from altanalyze3.components.visualization.psi_boxplot_grid import observed_psi_limits
    assert observed_psi_limits([{'PSI':None},{'PSI':.123},{'PSI':.276}],.1)==pytest.approx((.1,.3))
    assert observed_psi_limits([{'PSI':0},{'PSI':None}],.1)==pytest.approx((0,.1))
    assert observed_psi_limits([{'PSI':1}],.1)==pytest.approx((.9,1))
    with pytest.raises(ValueError,match='No observed'):
        observed_psi_limits([{'PSI':None}],.1)


def test_adaptive_axes_keep_literal_values_with_headers_and_no_annotation_headroom(tmp_path):
    p=manifest(tmp_path);out=tmp_path/'adaptive.pdf'
    r=main(['--panel-manifest',str(p),'--columns','2','--width','12','--height','8',
            '--panel-range-step','.1','--out',str(out)])
    assert len(PdfReader(out).pages)==1 and r['panel_range_step']==.1
    for panel in r['panels']:
        assert panel['y_limits']==pytest.approx([0,.7])
        assert [x['PSI'] for x in panel['records']]==[0,None,.1,.7]
        assert panel['all_defined_points_in_view']


def test_float32_decimal_boundary_keeps_actual_point_without_extra_tenth():
    import numpy as np
    from altanalyze3.components.visualization.psi_boxplot_grid import observed_psi_limits
    lo=float(np.float32(.7));hi=float(np.float32(.8))
    bounds=observed_psi_limits([{'PSI':lo},{'PSI':hi}],.1)
    assert [round(x,1) for x in bounds]==[.7,.8]
    assert bounds[0]<=lo and bounds[1]>=hi
    flat=observed_psi_limits([{'PSI':lo}],.1)
    assert flat[1]-flat[0]==pytest.approx(.1)
    assert flat[0]<=lo<=flat[1]


def test_selected_isoform_count_grid_calculates_ratios_and_labels_exact_ids(tmp_path):
    roster=tmp_path/'roster.json';roster.write_text(json.dumps({'healthy':['H1','H0'],'MDS':['M1']}))
    counts=tmp_path/'counts.tsv'
    counts.write_text('cohort\tdonor\tisoform_id\tassigned_RNA\n'
                      'healthy\tH1\tBLUE\t1\nhealthy\tH1\tGOLD\t3\n'
                      'healthy\tH0\tBLUE\t0\nhealthy\tH0\tGOLD\t0\n'
                      'MDS\tM1\tBLUE\t8\nMDS\tM1\tGOLD\t2\n')
    uid='G:E1-E3|G:E2-E3'
    panel=dict(gene='G',event_uid=uid,event_label='G:E1-E3|E2-E3',
               psi_source='selected_isoform_counts',isoform_definition=dict(numerator=['BLUE'],other=['GOLD']),
               native_argv=['--counts',str(counts),'--counts-sha256',hashlib.sha256(counts.read_bytes()).hexdigest(),
                            '--event-uid',uid,'--donor-roster',str(roster),'--cohort-order','healthy','MDS',
                            '--inclusion-isoforms','BLUE','--exclusion-isoforms','GOLD',
                            '--ylabel','Selected-isoform PSI','--out',str(tmp_path/'unused.pdf')])
    path=tmp_path/'panels.json';path.write_text(json.dumps(dict(required_genes=['G'],panels=[panel],title='Selected isoforms',caption='PSI = numerator/(numerator+other); no imputation.')))
    out=tmp_path/'counts_grid.pdf'
    result=main(['--panel-manifest',str(path),'--columns','2','--width','12','--height','8',
                 '--panel-range-step','.1','--out',str(out)])
    p=result['panels'][0]
    assert [r['PSI'] for r in p['records']]==[.25,None,.8]
    assert p['isoform_definition']=={'numerator':['BLUE'],'other':['GOLD']}
    assert p['counts_sha256']==hashlib.sha256(counts.read_bytes()).hexdigest()
    assert all(r['event_uid']==uid for r in p['records'])
    text=PdfReader(out).pages[0].extract_text()
    assert 'Numerator: BLUE' in text and 'Other: GOLD' in text
    panel['isoform_definition']['numerator']=['OTHER_ID']
    path.write_text(json.dumps(dict(required_genes=['G'],panels=[panel],title='Bad',caption='')))
    with pytest.raises(ValueError,match='exact explicit selected-isoform'):
        main(['--panel-manifest',str(path),'--out',str(tmp_path/'bad_counts.pdf')])
