"""Opt-in rho-only discovery and source-preserving event-median PSI filling."""
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from altanalyze3.components.visualization import marker_heatmap_h5ad as mh


class RhoOnlyMarkerHeatmap(unittest.TestCase):
    def test_assignment_and_order_ignore_conflicting_or_missing_pvalues(self):
        groups = ["Young", "Shared"]
        rho = pd.DataFrame([[.61, .9], [.8, -.8], [.7, .71]],
                           index=["a", "b", "c"], columns=groups)
        p = pd.DataFrame([[.001, 1.], [np.nan, .0001], [.00001, .9]],
                         index=rho.index, columns=groups)
        result = mh._select_unique_markers(p, groups, None, effect_df=rho,
                                          min_rho=.6, pval_threshold=None, ranking="rho")
        self.assertEqual(result.gene.tolist(), ["b", "a", "c"])
        self.assertEqual(result.cluster.astype(str).tolist(), ["Young", "Shared", "Shared"])
        np.testing.assert_array_equal(result.effect, [.8, .9, .71])
        with self.assertRaises(ValueError):
            mh._select_unique_markers(p, groups, None, effect_df=rho, ranking="rho")

    def test_default_selection_keeps_existing_fdr_rule(self):
        rho = pd.DataFrame([[.7, .9], [.9, -.9]], index=["a", "b"], columns=["A", "B"])
        p = pd.DataFrame([[.001, .02], [.5, .5]], index=rho.index, columns=rho.columns)
        result = mh._select_unique_markers(p, ["A", "B"], None, effect_df=rho,
                                          min_rho=.6, pval_threshold=.01)
        self.assertEqual(result.gene.tolist(), ["a"])
        self.assertEqual(result.cluster.astype(str).tolist(), ["A"])

    def test_hook_preserves_source_and_records_median_filling(self):
        values = np.array([[.2, .8, np.nan, .4], [.2, .8, np.nan, .4],
                           [.2, .8, np.nan, .4], [.9, .1, np.nan, .4],
                           [.9, .1, np.nan, .4], [np.nan, .1, np.nan, .4]])
        data = ad.AnnData(X=values.copy(),
                          obs=pd.DataFrame({"group": ["Young"]*3+["Shared"]*3},
                                           index=["a", "b", "c", "d", "e", "f"]),
                          var=pd.DataFrame(index=["event1", "event2", "all_missing", "constant"]))
        data.uns["lineage_order"] = ["Young", "Shared"]
        with tempfile.TemporaryDirectory() as directory:
            result = mh.generate_marker_heatmap_from_adata(
                data, cluster_key="group", out=str(Path(directory)/"test.pdf"),
                marker_method="markerfinder", marker_ranking="rho", pval_threshold=None,
                missing_psi="row_median", min_rho=np.nextafter(.6, 1.),
                max_missing_fraction=1., validate_scaling=False, centroid_method="mean",
                top_n=100, cells_per_cluster=0, return_statistics=True,
                render_heatmap=False, write_expression_tsv=False)
            np.testing.assert_array_equal(data.X, values)
            self.assertEqual(result["rho"].index.tolist(), data.var_names.tolist())
            self.assertTrue(result["rho"].loc["all_missing"].isna().all())
            self.assertTrue(result["rho"].loc["constant"].isna().all())
            self.assertEqual(result["missing_PSI_handling"]["filled_cells"], 1)
            self.assertEqual(result["missing_PSI_handling"]["all_missing_events"], 1)
            cache = ad.read_h5ad(result["heatmap_cache"])
            self.assertEqual(int(cache.layers["imputed_PSI"].sum()), 1)
            j = cache.var.gene.tolist().index("event1")
            self.assertTrue(np.isnan(cache.layers["source_PSI"][-1, j]))
            filled = mh.fill_PSI_row_median(pd.DataFrame(values.T))[0]
            self.assertEqual(filled.iloc[0, -1], .2)
            expected = (filled.iloc[0] - filled.iloc[0].mean())/filled.iloc[0].std(ddof=0)
            np.testing.assert_allclose(cache.X[:, j], expected, atol=1e-6)


if __name__ == "__main__":
    unittest.main()
