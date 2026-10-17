"""The indexed reference must preserve the native exhaustive lookup results."""
import copy

import pytest

from altanalyze3.components.visualization import isoform_structure_view as S


def reference():
    # Interleave genes, reversed bounds, introns and adjacent exon subregions.
    return {
        ('B', 'E2.1'): {'start': 80, 'end': 70, 'type': 'E'},
        ('A', 'E2.1'): {'start': 30, 'end': 40, 'type': 'E'},
        ('B', 'I1.1'): {'start': 31, 'end': 69, 'type': 'I'},
        ('A', 'E1.1'): {'start': 10, 'end': 29, 'type': 'E'},
        ('C', 'E1.1'): {'start': 100, 'end': 120},
        ('B', 'E1.1'): {'start': 10, 'end': 30, 'type': 'E'},
    }


def assert_native(lookup):
    plain = dict(lookup)
    for genes in ({'A'}, {'B'}, {'C'}, {'absent'}, {'C', 'B', 'A'}, set()):
        expected = S._gene_feature_index(plain, genes)
        actual = S._gene_feature_index(lookup, genes)
        assert actual == expected
        assert list(actual[0]) == list(expected[0])
        assert list(actual[1]) == list(expected[1])


def test_preserves_all_features_values_and_cross_gene_order():
    original = reference()
    indexed = S._GeneIndexedExonLookup(original)
    assert list(indexed.items()) == list(original.items())
    assert all(indexed[key] is value for key, value in original.items())
    assert_native(indexed)
    assert list(S._gene_feature_index(indexed, {'A', 'B', 'C'})[0]) == ['B', 'A', 'C']
    assert S._gene_feature_index(indexed, {'A'})[1] == {'A': [(10, 40)]}


@pytest.mark.parametrize('mutation', [
    lambda x: x.__setitem__(('A', 'E3.1'), {'start': 50, 'end': 60}),
    lambda x: x.__setitem__(('A', 'E2.1'), {'start': 45, 'end': 50}),
    lambda x: x.__delitem__(('B', 'E2.1')),
    lambda x: x.pop(('A', 'E1.1')),
    lambda x: x.popitem(),
    lambda x: x.clear(),
    lambda x: x.setdefault(('A', 'E3.1'), {'start': 70, 'end': 80}),
    lambda x: x.setdefault(('A', 'E1.1'), {'start': 700, 'end': 800}),
    lambda x: x.update({('C', 'I2.1'): {'start': 121, 'end': 140, 'type': 'I'}}),
    lambda x: x.__ior__({('D', 'E1.1'): {'start': 1, 'end': 3}}),
    lambda x: x[('A', 'E1.1')].update(start=1, end=50),
])
def test_mutations_preserve_exhaustive_lookup(mutation):
    indexed = S._GeneIndexedExonLookup(copy.deepcopy(reference()))
    assert_native(indexed)  # Populate the cache before each mutation.
    mutation(indexed)
    assert_native(indexed)
    assert S._gene_feature_index(indexed, {'D'}) == S._gene_feature_index(dict(indexed), {'D'})


def test_delete_and_reinsert_preserves_new_dictionary_order():
    indexed = S._GeneIndexedExonLookup(reference())
    assert_native(indexed)
    for key in [('B', 'E2.1'), ('B', 'I1.1'), ('B', 'E1.1')]:
        value = indexed.pop(key)
        indexed[key] = value
    assert_native(indexed)
    assert list(S._gene_feature_index(indexed, {'A', 'B', 'C'})[0]) == ['A', 'C', 'B']


def test_nonset_membership_uses_native_iteration():
    indexed = S._GeneIndexedExonLookup(reference())
    assert S._gene_feature_index(indexed, ['B', 'B', 'A']) == S._gene_feature_index(dict(indexed), ['B', 'B', 'A'])
    assert S._gene_feature_index(indexed, iter(['B', 'A'])) == S._gene_feature_index(dict(indexed), iter(['B', 'A']))


def test_reader_retains_complete_reference(tmp_path):
    path = tmp_path / 'reference.tsv'
    path.write_text('gene\texon-id\tchromosome\tstrand\texon-region-start(s)\texon-region-stop(s)\n'
                    'A\tE1.1\tchr1\t+\t10\t20\n'
                    'B\tE1.1\tchr2\t-\t30\t40\n')
    genes, lookup = S.read_gene_model(str(path))
    assert list(genes) == ['A', 'B']
    assert isinstance(lookup, S._GeneIndexedExonLookup)
    assert list(lookup) == [('A', 'E1.1'), ('B', 'E1.1')]
    assert_native(lookup)
