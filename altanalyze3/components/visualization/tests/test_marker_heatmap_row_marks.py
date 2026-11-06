import tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pytest
from altanalyze3.components.visualization import marker_heatmap_h5ad as mh


def test_marks_are_exact_oriented_rows_and_preserve_axes():
    a='ENSG00000000001:E1.1-E2.1|ENSG00000000001:E1.1-E3.1'
    b='|'.join(reversed(a.split('|')))
    fig,ax=plt.subplots();ax.imshow(np.zeros((2,3)));before=(ax.get_xlim(),ax.get_ylim())
    rows=mh._draw_event_row_marks(ax,['Healthy:'+a,'Mutant:'+b],[a])
    assert rows==[0]
    assert len(ax.lines)==1 and ax.lines[0].get_color()=='#d62728'
    np.testing.assert_array_equal(ax.lines[0].get_ydata(),[0,0])
    assert (ax.get_xlim(),ax.get_ylim())==before
    with pytest.raises(ValueError):mh._draw_event_row_marks(ax,[a],['ENSG999:E1-E2|ENSG999:E1-E3'])
    plt.close(fig)


def test_tight_export_callouts_keep_readable_font_and_offset():
    fig,ax=plt.subplots(figsize=(3,4));ax.imshow(np.zeros((100,3)))
    before=ax.get_ylim()
    mh._draw_row_callouts(ax,['FXR1:'+('E1-E2|'*20)+' [1: direct binding]'],[25],100,
                          x_offset=.05,fit_to_page=False)
    assert len(ax.texts)==1 and ax.texts[0].get_fontsize()==6.
    assert ax.texts[0].get_position()[0]==pytest.approx(1.085)
    assert ax.get_ylim()==before
    plt.close(fig)


def test_renderer_exports_tick_rows_and_exact_literature_callouts_without_changing_matrix():
    a='ENSG00000000001:E1.1-E2.1|ENSG00000000001:E1.1-E3.1'
    b='|'.join(reversed(a.split('|')))
    frame=pd.DataFrame([[-1.,1.],[1.,-1.]],index=[a,b],columns=['A','B']);before=frame.copy()
    with tempfile.TemporaryDirectory() as folder:
        output=str(Path(folder)/'marks.pdf')
        with patch.object(mh,'_save_figure'):
            mh._plot_heatmap(frame,output,[('A',1),('B',1)],['A','B'],['A','B'],['A','B'],
                label_events=[],row_marks=[a],row_mark_label='Same-direction S34 reference',
                callout_label_overrides={a:'Exact event [1]'},show_pattern_labels=True)
        marks=pd.read_csv(output.replace('.pdf','_row_marks.tsv'),sep='\t')
        labels=pd.read_csv(output.replace('.pdf','_label_events.tsv'),sep='\t')
        assert marks.event.tolist()==[a] and marks.heatmap_row.tolist()==[1]
        assert labels.event.tolist()==[a] and labels.label.tolist()==['Exact event [1]']
        assert labels.match.tolist()==['explicit literature callout']
    pd.testing.assert_frame_equal(frame,before,check_exact=True)
