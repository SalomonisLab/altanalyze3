"""Requested figure patterns never change the native fit or full input roster."""
import hashlib
import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.cellHarmony.iterative_marker_finder import (
    Scheme, review_legacy_markerfinder_fit, run_legacy_iterative_marker_finder,
)


def test_display_patterns_retain_other_assignments_and_report_empty_design(tmp_path):
    groups = np.repeat(['C', 'P', 'R'], 3)
    jitter = np.tile([-.02, 0, .02], 3)
    x = np.column_stack([np.where(groups == 'C', .8, .2) + jitter,
                         np.where(groups == 'P', .8, .2) + jitter,
                         np.ones(9) * .5, np.full(9, np.nan)])
    a = ad.AnnData(x, obs=pd.DataFrame({'group': groups}, index=[str(i) for i in range(9)]),
                   var=pd.DataFrame(index=['C_event', 'P_event', 'constant', 'unmeasured']))
    a.uns['lineage_order'] = ['C', 'P', 'R']
    source = tmp_path / 'source.h5ad'
    a.write_h5ad(source)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    fit = tmp_path / 'fit'
    assignments, parameters = run_legacy_iterative_marker_finder(
        [Scheme('clinical', a, 'group', ['C', 'P', 'R'])], a, 'group', fit,
        rho_cutoff=.6, top_n=2, heatmap_dpi=100)
    assert parameters['source_features'] == a.var_names.tolist()
    before = (fit / 'Consolidated-MarkerFinder.assignments.tsv').read_bytes()
    manifest = dict(display_h5ad=str(source), display_group_key='group',
                    schemes=[dict(name='clinical', h5ad=str(source), group_key='group',
                                  level_order=['C', 'P', 'R'])])
    out = tmp_path / 'review'
    result = review_legacy_markerfinder_fit(
        manifest, fit, out, input_sha256=digest, rho_cutoff=.7, top_n=1,
        min_observed=3, min_fraction=1., strict_top_n=True,
        display_patterns=['up-P', 'up-R'], dpi=100)
    cache = ad.read_h5ad(result['heatmap_cache'])
    assert cache.var.gene.tolist() == ['P_event']
    audit = pd.read_csv(out / 'all_original_assignments_and_review_gates.tsv', sep='\t', index_col=0)
    assert set(audit.index) == set(assignments.index)
    assert not audit.loc['C_event', 'display_pattern_requested']
    review = json.loads((out / 'native_review_parameters.json').read_text())
    assert review['requested_pattern_eligible_counts'] == {'up-P': 1, 'up-R': 0}
    assert (fit / 'Consolidated-MarkerFinder.assignments.tsv').read_bytes() == before
    assert hashlib.sha256(source.read_bytes()).hexdigest() == digest
    assert np.array_equal(ad.read_h5ad(source).X, x, equal_nan=True)
    with pytest.raises(ValueError, match='no source design'):
        review_legacy_markerfinder_fit(manifest, fit, tmp_path / 'bad', input_sha256=digest,
            rho_cutoff=.7, strict_top_n=True, display_patterns=['up-unknown'], dpi=100)
