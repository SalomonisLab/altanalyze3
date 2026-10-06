import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from scipy.stats import hypergeom

from altanalyze3.components.clustering import ICGS
from altanalyze3.components.cellHarmony.scalable_discover.pipeline import attach_serving_gene_symbols


@pytest.fixture
def catalog(tmp_path):
    path = tmp_path / 'catalog.tsv'
    pd.DataFrame({'Ensembl': [f'ENSG0000000{i:04}' for i in [1, 2, 3, 4, 3, 4, 5]],
                  'Gene': ['A', 'B', 'C', 'D', 'C', 'D', 'E'],
                  'Term': ['Lung Type One'] * 4 + ['Lung Type Two'] * 3}).to_csv(path, sep='\t', index=False)
    return path


def run_enrichment(tmp_path, catalog, genes, background):
    config = ICGS.ICGS3Config(input_paths=[], output_dir=str(tmp_path), biomarker_file=str(catalog),
                            minimal_outputs=True)
    markers = pd.DataFrame({'marker': genes, 'top_cluster': ['C18'] * len(genes)})
    labels = ICGS.biomarker_enrichment(markers, background, config, str(tmp_path))
    evidence = pd.read_csv(tmp_path / 'GO-Elite/icgs3_biomarker_enrichment.tsv', sep='\t')
    recorded_bg = pd.read_csv(tmp_path / 'GO-Elite/icgs3_biomarker_background.tsv', sep='\t')
    assert list(recorded_bg.feature) == background
    return labels, evidence


@pytest.mark.parametrize('versions', [False, True])
def test_ensembl_and_symbol_inputs_have_identical_statistics_full_background(tmp_path, catalog, versions):
    symbols = ['A', 'B', 'C', 'D', 'E'] + [f'OTHER{i}' for i in range(95)]
    ids = [f'ENSG0000000{i:04}' + ('.12' if versions else '') for i in range(1, 6)] + symbols[5:]
    left, symbol_evidence = run_enrichment(tmp_path / 'symbols', catalog, symbols[:3], symbols)
    right, id_evidence = run_enrichment(tmp_path / 'ids', catalog, ids[:3], ids)
    pd.testing.assert_frame_equal(left, right)
    cols = ['term_name', 'p_value', 'fdr', 'overlap', 'query_size', 'term_size']
    pd.testing.assert_frame_equal(symbol_evidence[cols], id_evidence[cols])
    best = id_evidence[id_evidence.term_name == 'Lung Type One'].iloc[0]
    assert best.p_value == pytest.approx(hypergeom(100, 4, 3).sf(2))
    assert id_evidence.overlap_genes.iloc[0] == ','.join(ids[:3])


def test_unresolved_queries_fail_instead_of_silent_unk(tmp_path, catalog):
    config = ICGS.ICGS3Config(input_paths=[], output_dir=str(tmp_path), biomarker_file=str(catalog))
    markers = pd.DataFrame({'marker': ['unresolved'], 'top_cluster': ['C18']})
    with pytest.raises(ValueError, match='absent from its background'):
        ICGS.biomarker_enrichment(markers, ['A', 'B', 'C'], config, str(tmp_path))
    with pytest.raises(ValueError, match='No input identifiers match'):
        ICGS.biomarker_enrichment(markers, ['unresolved'], config, str(tmp_path))


def test_unknown_labels_are_lowercase_and_biological_labels_preserved():
    labels = pd.DataFrame({'cluster': ['C18', 'c19', '20'], 'term_name': ['UNK', '', 'Lung Type One']})
    actual = ICGS.clean_biomarker_prediction_labels(labels)
    assert list(actual.cell_type_prediction) == ['c18', 'c19', 'Lung Type One_c20']


def test_supplied_symbols_resolve_ensembl_without_changing_features_cells_or_expression():
    ids = ['ENSG00000000001', 'ENSG00000000002.3', 'CXCL12', 'ENSG00000000004']
    var = pd.DataFrame({'gene_symbols': ids, 'feature_name': pd.Categorical(['FGR', 'LYZ', 'CXCL12', None])},
                       index=ids)
    x = np.arange(12, dtype=np.float32).reshape(3, 4)
    a = ad.AnnData(sp.csr_matrix(x), obs=pd.DataFrame(index=['cell1', 'cell2', 'cell3']), var=var)
    a.layers['counts'] = a.X.copy()
    info = attach_serving_gene_symbols(a)
    assert list(a.var.gene_symbols) == ['FGR', 'LYZ', 'CXCL12', 'ENSG00000000004']
    assert list(a.var_names) == ids and list(a.obs_names) == ['cell1', 'cell2', 'cell3']
    np.testing.assert_array_equal(a.X.toarray(), x)
    np.testing.assert_array_equal(a.layers['counts'].toarray(), x)
    assert info['primary_ensembl_ids'] == 3
    assert info['symbol_sources'] == {'feature_name': 2, 'gene_symbols': 1, 'var_names': 1}
