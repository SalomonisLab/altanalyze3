import pandas as pd
import pytest

from altanalyze3.components.cellHarmony.webapp.feature_lookup import FeatureLookup


def test_ensembl_versions_case_and_supplied_symbols_preserve_primary_features():
    names = ['ENSG00000000001.2', 'ENSMUSG00000000002', 'CXCL12']
    var = pd.DataFrame({'feature_name': pd.Categorical(['SFTPC', 'Lyz2', 'CXCL12'])}, index=names)
    lookup = FeatureLookup(names, var)
    for request in ['ENSG00000000001.2', 'ENSG00000000001', 'sftpc']:
        assert lookup.resolve(request) == names[0]
    assert lookup.resolve('LYZ2') == names[1]
    assert lookup.resolve('cxcl12') == names[2]
    assert lookup.resolve('unresolved') is None
    assert set(names).issubset(lookup.suggestions)


def test_shared_symbol_is_explicitly_ambiguous_both_primary_features_remain():
    names = ['ENSG00000000001', 'ENSG00000000002']
    lookup = FeatureLookup(names, pd.DataFrame({'gene_symbols': ['SHARED', 'SHARED']}, index=names))
    with pytest.raises(ValueError, match='multiple primary IDs'):
        lookup.resolve('SHARED')
    assert lookup.resolve(names[0]) == names[0] and lookup.resolve(names[1]) == names[1]
    assert lookup.suggestions == names


def test_version_collisions_require_an_exact_native_id():
    names = ['ENSG00000000001.1', 'ENSG00000000001.2']
    lookup = FeatureLookup(names)
    with pytest.raises(ValueError, match='multiple primary IDs'):
        lookup.resolve('ENSG00000000001')
    assert lookup.resolve(names[0]) == names[0] and lookup.resolve(names[1]) == names[1]
