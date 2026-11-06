"""Small heatmaps must not embed one-pixel-per-cell images that viewers smooth."""
import tempfile
import subprocess
from PIL import Image
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pypdf import PdfReader
from altanalyze3.components.visualization import marker_heatmap_h5ad as mh

def test_small_pdf_has_page_resolution_flat_cells_and_retains_raster_flag():
    with tempfile.TemporaryDirectory() as directory:
        fig, ax = plt.subplots(figsize=(2, 2))
        artist = ax.imshow(np.array([[0., 1.], [1., 0.]]), cmap='gray',
                           interpolation='none', vmin=0, vmax=1)
        artist.set_rasterized(True)
        ax.set_axis_off()
        output = Path(directory)/'discrete.pdf'
        mh._save_figure(fig, str(output), dpi=100, write_svg=False)
        plt.close(fig)
        assert artist.get_rasterized()
        images = list(PdfReader(output).pages[0].images)
        assert len(images) == 1
        assert min(images[0].image.size) > 100  # not the original 2x2 image
        # Inspect a rendered page: packed 1-bit PDF image extraction can use a
        # different scanline layout than its correctly rendered PDF content.
        subprocess.run(['pdftoppm','-scale-to','400','-png','-singlefile',
                        str(output),str(Path(directory)/'page')],check=True)
        im = Image.open(Path(directory)/'page.png').convert('RGB')
        im = im.crop((51,48,359,356))
        quadrants = [im.getpixel((int(x*im.width), int(y*im.height)))
                     for x,y in [(.25,.25),(.75,.25),(.25,.75),(.75,.75)]]
        # PDF image storage may invert Y through its drawing transform.
        assert set(quadrants) == {(0,0,0), (255,255,255)}
        assert quadrants[0] == quadrants[3] and quadrants[1] == quadrants[2]
        pixels = np.asarray(im)
        # Ignore the two-pixel antialias border around actual cell edges.
        for x,y in [(0,0),(1,0),(0,1),(1,1)]:
            tile = pixels[y*154+5:(y+1)*154-5,x*154+5:(x+1)*154-5]
            assert np.unique(tile.reshape(-1,3),axis=0).shape == (1,3)
