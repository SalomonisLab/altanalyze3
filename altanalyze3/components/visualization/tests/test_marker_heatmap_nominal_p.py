"""Explicit nominal-P eligibility must preserve native scores and BH statistics."""
import tempfile
from pathlib import Path
from unittest.mock import patch
import anndata as ad
import numpy as np
import pandas as pd
from altanalyze3.components.visualization import marker_heatmap_h5ad as mh

def test_nominal_gate_changes_assignment_not_statistics_or_inputs():
    a=ad.AnnData(np.array([[.8,.7,.6],[.9,.8,.6],[.1,.2,.6],[.2,.1,.6]]),
        obs=pd.DataFrame({'group':['A','A','B','B']},index=list('abcd')),
        var=pd.DataFrame(index=['x','y','z']))
    p=pd.DataFrame([[.001,.001],[.04,.04],[.8,.8]],index=a.var_names,columns=['A','B'])
    r=pd.DataFrame([[.9,-.9],[.7,-.7],[.1,-.1]],index=a.var_names,columns=['A','B'])
    before=a.X.copy()
    with tempfile.TemporaryDirectory() as d,patch.object(mh,'_markerfinder_stats',return_value=(p,r)):
        kw=dict(cluster_key='group',marker_method='markerfinder',min_rho=.25,
            pval_threshold=.01,min_measured_per_group=2,centroid_method='mean',
            validate_scaling=False,cells_per_cluster=0,render_heatmap=False,
            write_expression_tsv=False,return_statistics=True)
        base=mh.generate_marker_heatmap_from_adata(a.copy(),out=str(Path(d)/'fdr.pdf'),**kw)
        nominal=mh.generate_marker_heatmap_from_adata(a.copy(),out=str(Path(d)/'raw.pdf'),raw_pval_threshold=.05,**kw)
        assert base['assignments'].gene.tolist()==['x']
        assert nominal['assignments'].gene.tolist()==['x','y']
        for metric in ['rho','pval','fdr']:
            pd.testing.assert_frame_equal(base[metric],nominal[metric],check_exact=True)
        table=pd.read_csv(nominal['markers_tsv'],sep='\t').set_index('Gene')
        assert table.loc['y','Raw p-value']==.04
        assert table.loc['y','FDR p-value']==.06
        assert nominal['sample_ids']==a.obs_names.tolist()
        np.testing.assert_array_equal(a.X,before)
