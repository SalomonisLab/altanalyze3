"""Native Wilcoxon score-only PSI route retains full source identities."""
import tempfile
from pathlib import Path
import unittest
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import ranksums
from altanalyze3.components.visualization import marker_heatmap_h5ad as mh
from altanalyze3.components.cellHarmony.iterative_marker_finder import (
    legacy_strongest_patterns, collapse_discrete_pattern_logic, Scheme)


class NativeWilcoxonScore(unittest.TestCase):
    def test_vectorized_fill_exactly_matches_pandas_baseline(self):
        rng=np.random.default_rng(2710)
        values=rng.random((401,16))
        values[rng.random(values.shape)<.3]=np.nan
        values[0,:]=np.nan
        values[1,:]=0
        frame=pd.DataFrame(values,index=['event'+str(i) for i in range(401)],
                           columns=['donor'+str(i) for i in range(16)])
        original=frame.copy(deep=True)
        expected=frame.T.fillna(frame.median(axis=1)).T
        result,mask=mh.fill_PSI_row_median(frame)
        pd.testing.assert_frame_equal(result,expected)
        pd.testing.assert_frame_equal(frame,original)
        pd.testing.assert_frame_equal(mask,frame.isna() & expected.notna())

    def test_float32_fill_preserves_baseline_reduction_layout_and_zscores(self):
        values=np.random.default_rng(2710).random((16,401)).astype('float32')
        values[values<.3]=np.nan
        frame=pd.DataFrame(values.T)
        expected=frame.T.fillna(frame.median(axis=1)).T
        result,_=mh.fill_PSI_row_median(frame)
        pd.testing.assert_frame_equal(result,expected)
        self.assertEqual(result.to_numpy().flags.c_contiguous,expected.to_numpy().flags.c_contiguous)
        pd.testing.assert_frame_equal(mh._zscore_rows(result),mh._zscore_rows(expected),check_exact=True)

    def test_filled_input_native_scores_match_scipy_rank_sum(self):
        raw=np.array([[.1,.8,np.nan,.4],[.2,.9,np.nan,.4],[.1,.8,np.nan,.4],
                      [.8,.1,np.nan,.4],[.9,.2,np.nan,.4],[np.nan,.1,np.nan,.4]])
        data=ad.AnnData(X=raw.copy(),obs=pd.DataFrame({'group':['Young']*3+['Shared']*3},
            index=list('abcdef')),var=pd.DataFrame(index=['up','down','all_missing','constant']))
        data.uns['lineage_order']=['Young','Shared']
        with tempfile.TemporaryDirectory() as folder:
            result=mh.generate_marker_heatmap_from_adata(data,cluster_key='group',
                out=str(Path(folder)/'rank.pdf'),marker_method='scanpy',method='wilcoxon',
                marker_ranking='statistic',pval_threshold=None,missing_psi='row_median',
                validate_scaling=False,centroid_method='mean',cells_per_cluster=0,
                top_n=100,return_statistics=True,render_heatmap=False,write_expression_tsv=False)
            np.testing.assert_array_equal(data.X,raw)
            self.assertEqual(result['scores'].index.tolist(),data.var_names.tolist())
            self.assertNotIn('rho',result)
            self.assertTrue(result['scores'].loc['all_missing'].isna().all())
            self.assertTrue(result['scores'].loc['constant'].isna().all())
            filled=mh.fill_PSI_row_median(pd.DataFrame(raw.T))[0].to_numpy().T
            for event in ['up','down']:
                j=data.var_names.get_loc(event)
                expected=ranksums(filled[3:,j],filled[:3,j]).statistic
                self.assertAlmostEqual(result['scores'].loc[event,'Shared'],expected,places=6)
            winners=legacy_strongest_patterns([('binary',result)],['Young','Shared'],
                rho_cutoff=0,score_key='scores')
            self.assertEqual(set(winners.index),{'up','down'})
            self.assertEqual(winners.loc['up','pattern'],'up-Shared')
            self.assertNotIn('rho',winners.columns)
            self.assertIn('score',winners.columns)
            collapsed,_=collapse_discrete_pattern_logic(winners,[Scheme('binary',data,'group',['Young','Shared'])])
            self.assertEqual(set(collapsed.index),{'up','down'})
            self.assertIn('score',collapsed.columns)
            self.assertIn('wilcoxon_z',pd.read_csv(result['markers_tsv'],sep='\t').columns)
            self.assertTrue(ad.read_h5ad(result['heatmap_cache']).layers['imputed_PSI'].sum()>0)

    def test_explicit_statistic_route_rejects_P_filter(self):
        data=ad.AnnData(np.ones((4,1)),obs=pd.DataFrame({'group':['A','A','B','B']},index=list('abcd')))
        with self.assertRaises(ValueError):
            mh.generate_marker_heatmap_from_adata(data,cluster_key='group',out='unused.pdf',
                marker_method='scanpy',method='wilcoxon',marker_ranking='statistic',pval_threshold=.05)


if __name__=='__main__':unittest.main()
