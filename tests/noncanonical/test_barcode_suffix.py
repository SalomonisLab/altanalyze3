"""Standard 10x lane suffixes must not suppress the requested orientation."""
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from altanalyze3.components.long_read import isoform_matrix as I

@pytest.mark.parametrize('barcode', ['AAGCAGGTGTGTACGT', 'AAGCAGGTGTGTACGT-1', 'AAGCAGGTGTGTACGT-2'])
def test_recognizes_cell_barcode_with_optional_lane(barcode):
    assert I.is_cell_barcode(barcode)

@pytest.mark.parametrize('barcode', ['WF83', 'ACGT', 'source006', 'AAGCAGGTGTGTACGT-library'])
def test_preserves_non_cell_identity(barcode):
    assert not I.is_cell_barcode(barcode)

@pytest.mark.parametrize('reverse', [False, True])
def test_h5ad_preserves_counts_and_transforms_only_cell_barcodes(tmp_path, reverse):
    raw=['AAGCAGGTGTGTACGT-1','AAGCAGGTGTGTACGT-2','WF83']
    matrix=np.array([[7,0],[0,3],[2,5]])
    p=tmp_path/'counts.h5ad'
    ad.AnnData(matrix,obs=pd.DataFrame(index=raw),var=pd.DataFrame(index=['gene:iso1','gene:iso2'])).write_h5ad(p)
    result=I.h5ad_to_adata(p,rev=reverse)
    expected=['ACGTACACACCTGCTT-1','ACGTACACACCTGCTT-2','WF83'] if reverse else raw
    assert list(result.obs_names)==expected
    np.testing.assert_array_equal(result.X,matrix)
    assert list(result.var_names)==['iso1','iso2']
