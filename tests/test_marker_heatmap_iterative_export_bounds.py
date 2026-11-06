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


def test_small_native_hook_at_2400dpi_never_allocates_a_full_page_rgba(tmp_path,monkeypatch):
    import matplotlib.backends.backend_agg as agg
    import anndata as ad
    import numpy as np
    import pandas as pd
    from pypdf import PdfReader
    monkeypatch.setenv('MARKER_HEATMAP_DPI','2400')
    init=agg.RendererAgg.__init__
    def guard(self,width,height,dpi):
        assert width*height*4 < 200*1024**2, 'Tiny heatmap attempted a >200MiB RGBA canvas'
        init(self,width,height,dpi)
    monkeypatch.setattr(agg.RendererAgg,'__init__',guard)
    x=np.array([[.1,.9,np.nan],[.2,.8,np.nan],[.15,.85,np.nan],
                [.8,.2,np.nan],[.9,.1,np.nan],[np.nan,.15,np.nan]])
    a=ad.AnnData(x,obs=pd.DataFrame({'state':['A']*3+['B']*3},index=list('abcdef')),
                var=pd.DataFrame(index=['up','down','unmeasured']))
    out=tmp_path/'small.pdf'
    result=mh.generate_marker_heatmap_from_adata(a,cluster_key='state',out=str(out),
        marker_method='markerfinder',validate_scaling=False,centroid_method='mean',
        min_rho=.4,pval_threshold=1.,max_missing_fraction=1.,return_statistics=True,cells_per_cluster=0)
    assert result['rho'].loc['unmeasured'].isna().all()
    assert a.shape==(6,3) and np.array_equal(a.X,x,equal_nan=True)
    images=list(PdfReader(out).pages[0].images)
    assert images and all(image.image.width<1000 and image.image.height<1000 for image in images)
    assert 'A' in PdfReader(out).pages[0].extract_text()
    assert out.with_suffix('.svg').is_file()


def test_large_matrix_keeps_high_dpi_raster_path(tmp_path,monkeypatch):
    import numpy as np
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots()
    image=ax.imshow(np.zeros((1001,1000)),interpolation='none')
    image.set_rasterized(True)
    calls=[]
    monkeypatch.setattr(fig,'savefig',lambda *args,**kw:calls.append((image.get_rasterized(),kw['dpi'])))
    try:
        mh._save_figure(fig,str(tmp_path/'large.pdf'),dpi=2400)
        assert calls==[(True,2400),(True,2400)]
        assert image.get_rasterized() is True
    finally:plt.close(fig)
