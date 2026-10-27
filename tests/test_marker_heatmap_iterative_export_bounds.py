"""Long event callouts must survive iterative PDF/SVG page export."""
from altanalyze3.components.visualization import marker_heatmap_h5ad as mh


def test_iterative_export_retains_long_callouts_in_pdf_and_svg(tmp_path):
    import matplotlib.pyplot as plt
    from pypdf import PdfReader
    fig,ax=plt.subplots(figsize=(2,2))
    label='RPL17:E2.10-RPL17-C18orf32:E2.2|RPL17-C18orf32:E2.2-E3.2'
    ax.text(1.1,.5,label,transform=ax.transAxes,fontsize=8)
    out=tmp_path/'long.pdf'
    mh._save_figure(fig,str(out),dpi=100,bbox_inches='tight')
    plt.close(fig)
    page=PdfReader(out).pages[0]
    assert label in page.extract_text()
    assert float(page.mediabox.width)>144
    assert label in out.with_suffix('.svg').read_text()
