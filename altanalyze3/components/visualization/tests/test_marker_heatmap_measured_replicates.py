"""PSI observations remain distinct from missing data and assignment eligibility."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from altanalyze3.components.visualization import marker_heatmap_h5ad as mh


class MeasuredReplicateMarkerHeatmap(unittest.TestCase):
    def data(self):
        values = np.array([[np.nan, np.nan, .8, 0.], [.8, .8, .85, 0.],
                           [.8, np.nan, .9, np.nan], [np.nan, np.nan, .1, 0.],
                           [0., 0., .12, 0.], [0., 0., .14, 0.]])
        events = [f'ENSG{i:011d}:E1.1-E2.1|ENSG{i:011d}:E1.1-E3.1'
                  for i in range(1, 5)]
        data = ad.AnnData(values, obs=pd.DataFrame({'group': ['pre']*3+['post']*3},
                         index=list('abcdef')), var=pd.DataFrame(index=events))
        data.uns['lineage_order'] = ['pre', 'post']
        return data

    def run_hook(self, data, directory, minimum=None, ranking='fdr'):
        return mh.generate_marker_heatmap_from_adata(
            data.copy(), cluster_key='group', out=str(Path(directory)/'markers.pdf'),
            marker_method='markerfinder', min_rho=.25, max_missing_fraction=.55,
            pval_threshold=.01 if ranking == 'fdr' else None,
            marker_ranking=ranking, validate_scaling=False, centroid_method='mean',
            cells_per_cluster=0, top_n=100, return_statistics=True,
            min_measured_per_group=minimum, render_heatmap=False,
            write_expression_tsv=False)

    def test_assignment_gate_preserves_all_scores_samples_and_source_PSI(self):
        data = self.data()
        original = data.X.copy()
        for ranking in ['fdr', 'rho']:
            with self.subTest(ranking=ranking), tempfile.TemporaryDirectory() as first, \
                    tempfile.TemporaryDirectory() as second:
                baseline = self.run_hook(data, first, ranking=ranking)
                gated = self.run_hook(data, second, minimum=2, ranking=ranking)
                for key in ['rho', 'pval', 'fdr']:
                    pd.testing.assert_frame_equal(baseline[key], gated[key], check_exact=True)
                    self.assertEqual(gated[key].index.tolist(), data.var_names.tolist())
                self.assertEqual(gated['sample_ids'], data.obs_names.tolist())
                self.assertIn(data.var_names[1], baseline['assignments'].gene.tolist())
                self.assertNotIn(data.var_names[1], gated['assignments'].gene.tolist())
                self.assertEqual(set(gated['assignments'].gene), set(data.var_names[[0, 2]]))
                audit = pd.read_csv(gated['measurement_eligibility_tsv'], sep='\t')
                self.assertEqual(len(audit), data.n_vars*2)
                self.assertFalse(audit.loc[audit.event == data.var_names[1],
                                           'eligible_for_assignment'].any())
                cache = ad.read_h5ad(gated['heatmap_cache'])
                raw = data[cache.obs_names, cache.var.gene.tolist()].X
                np.testing.assert_array_equal(cache.layers['source_PSI'], raw.astype(np.float32))
                np.testing.assert_array_equal(np.isnan(cache.X), np.isnan(raw))
                j = cache.var.gene.tolist().index(data.var_names[0])
                self.assertTrue(np.isfinite(cache.X[4:, j]).all())
                self.assertTrue((cache.X[4:, j] < 0).all(), 'Measured PSI=0 must not become missing')
        np.testing.assert_array_equal(data.X, original)

    def test_sparse_implicit_zeros_are_measured(self):
        data = self.data()
        data.X = sparse.csr_matrix(data.X)
        eligible, measured, reference = mh._measured_group_eligibility(
            data, 'group', False, None, 2, ['pre', 'post'])
        self.assertEqual(measured.loc[data.var_names[0]].tolist(), [2, 2])
        self.assertEqual(reference.loc[data.var_names[0]].tolist(), [2, 2])
        self.assertTrue(eligible.loc[data.var_names[0]].all())
        self.assertFalse(eligible.loc[data.var_names[1]].any())

    def test_non_PSI_default_display_path_is_unchanged(self):
        data = self.data()
        data.var_names = ['A', 'B', 'C', 'D']
        matrix, _, _, order = mh._build_heatmap(data, 'group', data.var_names.tolist(),
                                               ['pre', 'post'], False, None)
        source = pd.DataFrame(data[order, :].X.T, index=data.var_names, columns=order)
        pd.testing.assert_frame_equal(matrix, mh._zscore_rows(source))

    def test_missing_color_is_passed_to_the_existing_renderer(self):
        data = self.data()
        with tempfile.TemporaryDirectory() as directory:
            result = self.run_hook(data, directory, minimum=2)
            with patch.object(mh, '_plot_heatmap') as plot:
                mh.render_heatmap_from_cache(result['heatmap_cache'],
                    str(Path(directory)/'review.pdf'), load_go_terms=False,
                    missing_color='#d9d9d9')
            self.assertEqual(plot.call_args.kwargs['missing_color'], '#d9d9d9')
            self.assertTrue(plot.call_args.args[0].isna().to_numpy().any())


if __name__ == '__main__':
    unittest.main()
