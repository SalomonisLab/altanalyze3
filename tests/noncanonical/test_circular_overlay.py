from copy import deepcopy
import json
import pytest
from altanalyze3.components.visualization.isv_web.circular_overlay import marker_specs, draw_circular_markers, validate_circular_payload
from altanalyze3.components.visualization.isv_web import isv_web_remote as renderer


def molecule(strand='+'):
    return dict(count=1,isoform_id='read1',final_isoform_id='candidate',cluster_index=0,
                exon_segments=[dict(start=100,end=200,type='E',gene='G')],
                circular_junctions=[dict(donor=200 if strand=='+' else 100,
                    acceptor=100 if strand=='+' else 200,strand=strand,
                    coordinate_system='1-based-inclusive',validation_status='candidate')])


@pytest.mark.parametrize('strand', ['+', '-'])
def test_endpoint_identity_not_left_right_identity(strand):
    m=molecule(strand); before=deepcopy(m)
    specs=marker_specs(m,100,200)
    assert [(s['role'],s['position']) for s in specs]==[('donor',200 if strand=='+' else 100),('acceptor',100 if strand=='+' else 200)]
    assert all(s['linestyle']=='solid' for s in specs)
    assert m==before


def test_clipping_is_display_only():
    m=molecule(); before=deepcopy(m)
    assert [s['position'] for s in marker_specs(m,150,210)]==[200]
    assert m==before


@pytest.mark.parametrize('field,value', [('coordinate_system','0-based'),('strand','?'),('donor',-1),('donor',True),('donor',90)])
def test_invalid_endpoint_contract_rejected(field,value):
    m=molecule();m['circular_junctions'][0][field]=value
    with pytest.raises(ValueError):marker_specs(m,1,1000)


def test_validation_claim_requires_independent_source():
    m=molecule();j=m['circular_junctions'][0];j['validation_status']='orthogonally_validated'
    with pytest.raises(ValueError,match='source'):marker_specs(m,1,1000)
    j['validation_source']='fixture only: test orthogonal evidence'
    assert all(s['linestyle']=='solid' for s in marker_specs(m,1,1000))


def test_equal_diameter_and_correct_coordinates():
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots();m=molecule()
    artists=draw_circular_markers(ax,m,lambda x: 2*x,15,8,100,200)
    assert [a.center for a in artists]==[(400,15),(200,15)]
    assert [a.radius for a in artists]==[4,4]
    plt.close(fig)


def test_no_silent_overflow_or_multilocus_reprojection():
    with pytest.raises(ValueError,match='paginate'):
        validate_circular_payload(dict(panels=[dict(molecules=[molecule()]*4001)]))
    with pytest.raises(ValueError,match='single genomic'):
        validate_circular_payload(dict(has_multiple_loci=True,panels=[dict(molecules=[molecule()])]))


def test_native_renderer_preserves_evidence_and_overlay_count(tmp_path,monkeypatch):
    m=molecule();m['circular_evidence']=[dict(aligned_donor=197,donor=200,slide_nt=3,raw_cigar='100=')]
    original_evidence=deepcopy(m['circular_evidence']); called=[]
    import altanalyze3.components.visualization.isv_web.circular_overlay as overlay
    original=overlay.draw_circular_markers
    def spy(*args,**kwargs):
        artists=original(*args,**kwargs);called.extend(artists);return artists
    monkeypatch.setattr(overlay,'draw_circular_markers',spy)
    payload=dict(gene_model=dict(chrom='chr1',strand='+',extent=[90,210],blocks=[]),
                 panels=[dict(condition='test',n_molecules=1,molecules=[m])])
    assert renderer.render(payload,tmp_path/'native.pdf',show_read_labels=False,caption='Synthetic renderer test; not biological evidence.')
    assert (tmp_path/'native.pdf').stat().st_size>1000
    assert len(called)==2
    assert m['circular_evidence']==original_evidence
    assert m['count']==1


def test_linear_render_is_pixel_identical_to_saved_baseline(tmp_path,monkeypatch):
    import importlib.util
    from pathlib import Path
    import numpy as np
    from matplotlib.figure import Figure
    baseline_path=Path(__file__).parent/'fixtures/isv_web_remote_20261105.py'
    spec=importlib.util.spec_from_file_location('altanalyze3.components.visualization.isv_web._baseline_render_test',baseline_path)
    baseline=importlib.util.module_from_spec(spec);spec.loader.exec_module(baseline)
    frames=[]
    def capture(self,*args,**kwargs):
        self.canvas.draw();frames.append(np.asarray(self.canvas.buffer_rgba()).copy())
    monkeypatch.setattr(Figure,'savefig',capture)
    m=molecule();del m['circular_junctions']
    payload=dict(gene_model=dict(chrom='chr1',strand='+',extent=[90,210],blocks=[]),
                 panels=[dict(condition='linear control',n_molecules=1,molecules=[m])])
    baseline.render(deepcopy(payload),tmp_path/'baseline.pdf',show_read_labels=False)
    renderer.render(deepcopy(payload),tmp_path/'patched.pdf',show_read_labels=False)
    assert len(frames)==2
    np.testing.assert_array_equal(frames[0],frames[1])
