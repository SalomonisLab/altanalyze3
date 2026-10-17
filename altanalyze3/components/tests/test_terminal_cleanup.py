"""Approved terminal merges and adversarial controls for retained-catalog cleanup."""
import numpy as np
import pytest
from scipy import sparse

from altanalyze3.components.long_read.isoform_collapse.terminal_cleanup import (
    CatalogModel, junction_chain, plan_cleanup, remap_sparse_counts,
    terminal_merge_preserves_evidence,
)


def model(tid, exons, *, gene='gene', chrom='chr1', strand='+', known=False,
          total=3, own=3, key=None):
    return CatalogModel(gene, tid, chrom, strand, tuple(exons), total, known,
                        (-own, key or tid))


def plan(*rows):
    return plan_cleanup(rows, expected_identities=[r.identity for r in rows])


@pytest.mark.parametrize('gene,known,novel,ref,obs,counts', [
    ('ENSG00000159216', 'ENST00000455571', '63468836.BF71-N-0421_CD34',
     [(34880606,34880713),(34886843,34887096),(35048842,35048958),(35049168,35049298)],
     [(34880208,34880713),(34886843,34887096),(35048842,35048958),(35049168,35049459)], (339,3222)),
    ('ENSG00000119866', 'ENST00000643459', '17010213.AM72-N-0206_CD34',
     [(60461904,60462424),(60468732,60468833),(60545971,60546300),(60550868,60550907)],
     [(60460287,60462424),(60468732,60468833),(60545971,60546300),(60550868,60550907)], (16,3)),
    ('ENSG00000143995', 'ENST00000495021', '64775232.Y2982-E3',
     [(66439321,66439671),(66439843,66439984),(66440562,66440612),(66441414,66441464),
      (66442902,66443048),(66464109,66464220),(66512149,66512294),(66547943,66548019),
      (66567453,66567511),(66568667,66568756),(66569050,66569145),(66571246,66572759)],
     [(66439268,66439671),(66439843,66439984),(66440562,66440612),(66441414,66441464),
      (66442902,66443048),(66464109,66464220),(66512149,66512294),(66547943,66548019),
      (66567453,66567511),(66568667,66568756),(66569050,66569145),(66571246,66572762)], (393,25)),
])
def test_actual_three_approved_pairs_merge(gene, known, novel, ref, obs, counts):
    a = model(known, ref, gene=gene, known=True, total=counts[0])
    b = model(novel, obs, gene=gene, total=counts[1])
    r = plan(a, b)
    assert r['mapping'] == {a.identity:a.identity, b.identity:a.identity}
    assert r['totals'] == {a.identity:sum(counts)}
    assert a.exons == tuple(ref) and b.exons == tuple(obs)


def test_internal_retention_skipping_alternative_terminal_splices_and_fragments_stay_separate():
    base = model('ref', [(100,150),(200,250),(300,350),(400,450)], known=True)
    others = [
        model('retained', [(100,150),(200,350),(400,450)]),
        model('skipped', [(100,150),(300,350),(400,450)]),
        model('first_splice', [(100,160),(200,250),(300,350),(400,450)]),
        model('last_splice', [(100,150),(200,250),(300,355),(400,450)]),
        model('fragment', [(200,250),(300,350),(400,450)]),
    ]
    assert plan(base, *others)['merged_novel_identities'] == 0


def test_fragment_fitting_two_distinct_longer_chains_is_not_assigned_to_either():
    a = model('a', [(100,150),(200,250),(300,350)], known=True)
    b = model('b', [(10,50),(200,250),(300,350)], known=True)
    fragment = model('fragment', [(210,250),(300,340)])
    assert plan(a, b, fragment)['mapping'][fragment.identity] == fragment.identity


def test_one_base_site_shift_and_transitive_near_matches_never_merge():
    rows = [model(str(i), [(100,150),(200+i,250),(300,350)]) for i in range(7)]
    assert plan(*rows)['retained_identities'] == 7


def test_single_junction_merges_by_exact_sites_but_single_exons_never_share_empty_key():
    r = plan(model('one', [(100,150),(200,250)], known=True),
             model('two', [(110,150),(200,1000)]),
             model('mono1', [(100,1000)]), model('mono2', [(110,999)]))
    assert r['retained_identities'] == 3
    assert r['mapping']['gene','two'] == ('gene','one')


def test_locus_gene_and_strand_boundaries_and_continuous_annotation_blocks():
    exons = [(100,150),(200,250)]
    rows = [model('ref', exons, known=True), model('chrom', exons, chrom='chr2'),
            model('strand', exons, strand='-'), model('gene', exons, gene='other'),
            model('labels', [(100,120),(121,150),(200,220),(218,250)])]
    r = plan(*rows)
    assert r['merged_novel_identities'] == 1
    assert r['mapping']['gene','labels'] == ('gene','ref')


def test_known_ids_preserved_and_native_rank_not_catalog_abundance_selects_target():
    exons = [(100,150),(200,250)]
    a = model('ref1', exons, known=True, own=5, total=1)
    b = model('ref2', exons, known=True, own=9, total=100)
    n = model('novel', exons, own=900, total=900)
    r = plan(a, b, n)
    assert r['mapping'][a.identity] == a.identity
    assert r['mapping'][b.identity] == b.identity
    assert r['mapping'][n.identity] == b.identity
    small = model('own10', exons, own=10, total=10)
    big = model('own1', exons, own=1, total=500)
    r = plan(small, big)
    assert r['totals'] == {small.identity:510}
    # Reapplying to surviving original representatives is idempotent.
    again = plan(model('own10', exons, own=10, total=510))
    assert again['totals'] == r['totals']


def test_missing_duplicate_or_unresolved_ranking_rejects_without_dropping_identities():
    a = model('a', [(100,150),(200,250)])
    with pytest.raises(ValueError, match='coverage'):
        plan_cleanup([a], expected_identities=[a.identity, ('gene','missing')])
    with pytest.raises(ValueError, match='Duplicate'):
        plan(a, a)
    with pytest.raises(ValueError, match='ranking'):
        plan(a, model('b', a.exons, key='a'))
    with pytest.raises(ValueError, match='geometry'):
        plan(model('empty', []))


def test_sparse_sum_keeps_each_cell_and_integer_precision_and_all_source_features():
    x = sparse.csr_matrix(np.array([[16777216,1,0],[0,3,8],[0,0,0]], dtype=np.float32))
    crosswalk = {'a':'a', 'b':'a', 'c':'c'}
    y, ids = remap_sparse_counts(x, ['a','b','c'], crosswalk)
    assert ids == ['a','c']
    np.testing.assert_array_equal(y.toarray(), [[16777217,0],[3,8],[0,0]])
    np.testing.assert_array_equal(y.sum(axis=1), x.astype(np.int64).sum(axis=1))
    for bad in ({'a':'a','c':'c'}, {'a':'b','b':'c','c':'c'}):
        with pytest.raises(ValueError):
            remap_sparse_counts(x, ['a','b','c'], bad)
    for bad in (0.5, float('nan'), -1):
        with pytest.raises(ValueError, match='raw counts'):
            remap_sparse_counts(sparse.csr_matrix([[bad]]), ['a'], {'a':'a'})


def test_exact_chain_preserves_internal_geometry_on_both_strands():
    for strand in ('+', '-'):
        a = model('a', [(100,150),(200,250),(300,350)], strand=strand)
        b = model('b', [(1,150),(200,250),(300,9000)], strand=strand)
        r = plan(a, b)
        assert r['merged_novel_identities'] == 1
        assert junction_chain(a.exons) == junction_chain(b.exons)


def test_actual_hopx_terminal_coverage_at_splice_boundary_is_not_artificially_clipped():
    old = [(56648000,56648797),(56655857,56656012),(56656340,56656581)]
    new = [(56647998,56648797),(56655857,56656012),(56656340,56656522)]
    evidence = dict(protected_coverage=[(56656340,56656442),(56656526,56656668)],
                    protected_junctions=[],protected_splice_spans=[])
    assert terminal_merge_preserves_evidence(evidence, old, new, splice_slop=3)
    # Actual internal gap coverage remains forbidden, even when outer ends vary.
    assert not terminal_merge_preserves_evidence(
        dict(protected_coverage=[(56650000,56650100)]),old,new,splice_slop=3)
    # A source splice/continuous-span contradiction cannot be waived as an end.
    assert not terminal_merge_preserves_evidence(
        dict(protected_junctions=[(56648797,56655858+10)]),old,new,splice_slop=3)
    assert not terminal_merge_preserves_evidence(
        dict(protected_splice_spans=[(56648797,56655857)]),old,new,splice_slop=3)
    with pytest.raises(ValueError,match='same complete'):
        terminal_merge_preserves_evidence({},old,[(56648000,56656581)],splice_slop=3)
