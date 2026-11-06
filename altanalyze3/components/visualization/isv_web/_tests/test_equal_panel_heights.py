"""Cohort layout must preserve all RNA while avoiding sparse-panel inflation."""
from copy import deepcopy
from unittest.mock import patch
import pytest
from altanalyze3.components.visualization.isv_web import isv_web_remote as remote


def response():
    panels=[]
    for label,n in [('Healthy',66),('Pre',41),('Post',45),('AML',45)]:
        rows=[]
        for i in range(n):
            rows.append(dict(isoform_id=f'{label}_{i}',final_isoform_id='A' if i%2 else 'B',
                known=True,count=1,cluster_index=i%2,color_index=i%2,
                exon_segments=[dict(start=100,end=200,type='E'),dict(start=300,end=400,type='E')]))
        panels.append(dict(condition=label,n_molecules=n,molecules=rows))
    return dict(gene='G',symbol='G',gene_model=dict(chrom='chr1',strand='+',extent=[100,400],
        blocks=[dict(start=100,end=200,exon_id='E1'),dict(start=300,end=400,exon_id='E2')]),panels=panels)


def panel_boxes(fig):
    return [p.get_height() for p in fig.axes[0].patches
            if p.get_x()==6 and p.get_width()==1088 and p.get_facecolor()[:3]==(1.,1.,1.)]


def test_equal_height_render_keeps_every_read_and_input_unchanged(tmp_path):
    source=response();before=deepcopy(source);figures=[]
    with patch.object(remote.plt,'close',side_effect=lambda fig: figures.append(fig)):
        assert remote.render(source,str(tmp_path/'equal.pdf'),show_read_labels=False,equal_panel_heights=True)
    assert source==before
    boxes=panel_boxes(figures[-1]);assert len(boxes)==4 and min(boxes)>0
    assert boxes==pytest.approx([boxes[0]]*4)
    # Header counts remain real depth; no downsampling or cohort removal.
    labels=[t.get_text() for t in figures[-1].axes[0].texts]
    for count in [66,41,45]:assert any(f'{count} reads' in t for t in labels)


def test_default_layout_and_cli_are_backward_compatible(tmp_path):
    source=response();figures=[]
    with patch.object(remote.plt,'close',side_effect=lambda fig: figures.append(fig)):
        remote.render(source,str(tmp_path/'default.pdf'),show_read_labels=False)
        remote.render(source,str(tmp_path/'explicit_default.pdf'),show_read_labels=False,equal_panel_heights=False)
    assert panel_boxes(figures[-2])==panel_boxes(figures[-1])
    args=remote.build_arg_parser().parse_args(['--run_dir',str(tmp_path),'--gene','G','--equal-panel-heights'])
    with patch.object(remote,'render',return_value=True) as draw:
        assert remote._render_selected(args,source,str(tmp_path/'cli.pdf'),None)
    assert draw.call_args.kwargs['equal_panel_heights'] is True
    assert draw.call_args.args[0] is source
