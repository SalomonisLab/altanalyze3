"""Authorized full-panel BH for imputed lipids; preserve every raw test and fold."""
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

KEYS = ['contrast', 'unit', 'population', 'model']


def correct_full_panel_bh(frame, lipids, models):
    """Correct each complete model/comparison over the exact instructed lipid panel.

    A missing feature, duplicate identity, nonfinite p-value or unmatched model
    comparison is an error. Neither statistical filtering nor intersection can
    reduce the correction family. This function never fits a model or calculates
    a replacement raw test.
    """
    expected = set(lipids)
    if len(expected) != len(lipids) or not expected:
        raise ValueError('Required lipid identities must be unique and nonempty.')
    if set(frame.model) != set(models):
        raise ValueError('Unexpected or missing comparator model.')
    if frame[KEYS + ['gene']].isna().any().any():
        raise ValueError('Missing comparison or lipid identity.')
    if frame.duplicated(KEYS + ['gene']).any():
        raise ValueError('Duplicate lipid test identity.')
    p = frame.pval.to_numpy(dtype=float)
    if not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError('Invalid raw p-value: stop rather than shrink the panel.')
    comparisons = {}
    result = frame.copy()
    if 'fdr_legacy_filtered' not in result:
        result['fdr_legacy_filtered'] = result.fdr
    corrected = np.full(len(frame), np.nan)
    for key, positions in frame.groupby(KEYS, sort=False, observed=True).indices.items():
        group = frame.iloc[positions]
        if len(group) != len(lipids) or set(group.gene) != expected:
            raise ValueError(f'Incomplete required lipid panel: {key}')
        comparisons.setdefault(key[-1], set()).add(key[:-1])
        corrected[positions] = multipletests(group.pval.to_numpy(), method='fdr_bh')[1]
    if any(comparisons[m] != comparisons[models[0]] for m in models):
        raise ValueError('The models do not have the same comparison identities.')
    result['fdr'] = corrected
    result['significant_BH_005'] = result.fdr.lt(.05) & result.log2fc.ne(0)
    result['significant_BH_010'] = result.fdr.lt(.1) & result.log2fc.ne(0)
    result['BH_panel_size'] = len(lipids)
    result['BH_policy'] = 'full_panel_no_abundance_filter'
    assert np.isfinite(result.fdr).all()
    changed = {'fdr', 'significant_BH_005', 'significant_BH_010'}
    pd.testing.assert_frame_equal(frame.drop(columns=list(changed), errors='ignore'),
                                  result[frame.columns].drop(columns=list(changed), errors='ignore'))
    return result
