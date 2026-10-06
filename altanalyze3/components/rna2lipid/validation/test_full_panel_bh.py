"""Meaningful coverage and scale-invariance checks for the approved BH correction."""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from full_panel_bh import correct_full_panel_bh


class FullPanelBHTests(unittest.TestCase):
    def data(self):
        return pd.DataFrame([{'contrast': 'IPF', 'unit': 'PB', 'population': 'AT2',
                              'model': model, 'gene': lipid, 'pval': p, 'fdr': np.nan,
                              'log2fc': fold, 'case_mean_expr': -5., 'control_mean_expr': -6.,
                              'significant_BH_005': False, 'significant_BH_010': False}
                             for model in ['prior', 'candidate']
                             for lipid, p, fold in [('a', .001, 1.), ('b', .02, -1.), ('c', .8, .1)]])

    def test_matches_independent_BH_reference_and_retains_raw_statistics(self):
        source = self.data()
        result = correct_full_panel_bh(source, ['a', 'b', 'c'], ['prior', 'candidate'])
        np.testing.assert_allclose(result.fdr, np.tile(multipletests([.001, .02, .8], method='fdr_bh')[1], 2))
        pd.testing.assert_frame_equal(source[['gene', 'pval', 'log2fc']], result[['gene', 'pval', 'log2fc']])
        self.assertTrue(result.fdr_legacy_filtered.isna().all())

    def test_negative_or_shifted_means_do_not_change_eligibility(self):
        source = self.data()
        result = correct_full_panel_bh(source, ['a', 'b', 'c'], ['prior', 'candidate'])
        shifted = source.copy()
        shifted[['case_mean_expr', 'control_mean_expr']] += 100
        other = correct_full_panel_bh(shifted, ['a', 'b', 'c'], ['prior', 'candidate'])
        np.testing.assert_array_equal(result.fdr, other.fdr)
        self.assertTrue(result.fdr.notna().all())

    def test_missing_duplicate_or_swapped_identity_blocks(self):
        source = self.data()
        for bad in [source.iloc[1:], pd.concat([source, source.iloc[[0]]], ignore_index=True),
                    source.assign(gene=source.gene.replace('c', 'unexpected'))]:
            with self.assertRaises(ValueError):
                correct_full_panel_bh(bad, ['a', 'b', 'c'], ['prior', 'candidate'])

    def test_nonfinite_pvalue_blocks_instead_of_reducing_panel(self):
        source = self.data()
        source.loc[0, 'pval'] = np.nan
        with self.assertRaises(ValueError):
            correct_full_panel_bh(source, ['a', 'b', 'c'], ['prior', 'candidate'])

    def test_unpaired_comparison_blocks(self):
        source = self.data()
        source.loc[source.model.eq('candidate'), 'population'] = 'AM'
        with self.assertRaises(ValueError):
            correct_full_panel_bh(source, ['a', 'b', 'c'], ['prior', 'candidate'])


if __name__ == '__main__':
    unittest.main()
