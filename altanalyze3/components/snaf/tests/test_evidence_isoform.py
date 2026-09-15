"""Construction and leakage boundaries for evidence-driven synthetic isoforms."""
import numpy as np
import pytest

from altanalyze3.components.snaf.surface import evidence_isoform as E


@pytest.mark.parametrize('strand', ['+', '-', '1', '-1'])
def test_unknown_first_exon_is_250nt_and_translated_in_strand_direction(strand):
    P = E.P
    minus = strand in ('-', '-1')
    target = (500, 1000) if minus else (1000, 1500)
    chain = [(351, 500), (1700, 1800)] if minus else [(100, 200), (1500, 1649)]
    expected = [(351, 500), (1000, 1249)] if minus else [(751, 1000), (1500, 1649)]
    models = {'T': (strand, chain)}
    mrna = 'C' * 200 + 'ATG' + 'GCC' * 49 + 'TAA' + 'C' * 47
    genomic = P.reverse_complement(mrna) if minus else mrna
    dna = list('C' * 2000)
    offset = 0
    for s, e in expected:
        size = e - s + 1
        dna[s - 1:e] = genomic[offset:offset + size]
        offset += size
    fetch = lambda s, e: ''.join(dna[s - 1:e])
    ev = E.Evidence({target: [10]})
    pred = P.predict_isoform(models, *target, fetch, evidence=ev, first_exon=True)
    assert pred is not None
    assert pred.chain == expected
    assert pred.mrna == mrna
    assert pred.protein == 'M' + 'A' * 49
    assert P.prediction_row(pred)['first_exon_boundary_source'] == 'assumed_250nt'
    # Selection works with a reference-only hypothesis, too.
    assert P.predict_isoform(models, *target, fetch, first_exon=True).chain == expected


@pytest.mark.parametrize('strand,target,chain,expected', [
    ('+', (92323960, 92324060), [(92323605, 92323960), (92324060, 92324200)],
     [(92323605, 92323960), (92324060, 92324200)]),
    ('-', (500, 1000), [(351, 500), (1000, 1355)], [(351, 500), (1000, 1355)]),
])
def test_annotated_356nt_first_exon_is_not_replaced_by_fallback(strand, target, chain, expected):
    assert list(E.first_exon_seeds(chain, strand, target, {'T': (strand, chain)})) == [
        (expected, 'annotated')]


def test_first_exon_assumption_respects_chromosome_edges_and_splice_gap():
    models = {'T': ('+', [(400, 500), (1000, 1100)])}
    assert list(E.first_exon_seeds(models['T'][1], '+', (200, 1000), models)) == []
    with pytest.raises(ValueError, match='splice gap'):
        list(E.first_exon_seeds(models['T'][1], '+', (999, 1000), models))
    # A minus-strand fallback running off the far contig end is unresolved.
    minus_models = {'T': ('-', [(351, 500), (800, 900)])}
    assert E.generate_candidates(minus_models, (500, 1000),
        lambda s, e: 'C' * max(0, min(e, 1100) - s + 1),
        E.Evidence({(500, 1000): [10]}), first_exon=True) == []


def test_first_exon_selection_does_not_retain_an_upstream_junction(monkeypatch):
    models = {'T': ('+', [(100, 200), (1500, 1649)])}
    target, upstream = (1000, 1500), (800, 900)
    captured = []
    def translate(gene_models, *args, **kwargs):
        captured.extend(chain for _, chain in gene_models.values())
        return []
    monkeypatch.setattr(E.P, 'predict_isoform', translate)
    E.generate_candidates(models, target, lambda s, e: 'C' * (e - s + 1),
        E.Evidence({target: [10], upstream: [10]}), first_exon=True)
    assert captured
    assert all(E.introns(chain)[0] == target for chain in captured)


def test_competing_donors_and_crossing_introns_cannot_coexist():
    assert not E.compatible([(100, 300), (100, 400)])
    assert not E.compatible([(100, 300), (200, 400)])
    assert E.compatible([(100, 200), (250, 400)])


def test_pairwise_codetection_does_not_imply_joint_codetection():
    a, b, c = (100, 200), (300, 400), (500, 600)
    ev = E.Evidence({a: [5, 5, 0], b: [5, 0, 5], c: [0, 5, 5]})
    assert all(ev.supported(pair) for pair in [(a, b), (a, c), (b, c)])
    assert not ev.supported([a, b, c])
    assert all(len(js) < 3 for js in E.combinations(a, ev, set()))


def test_more_than_two_novel_junctions_are_reachable_in_one_sample():
    js = [(100, 200), (300, 400), (500, 600), (700, 800)]
    ev = E.Evidence({j: [5, 0] for j in js})
    assert tuple(js) in E.combinations(js[0], ev, set(), max_edits=4)


def test_sample_restriction_prevents_cross_sample_pairing():
    a, b = (100, 200), (250, 400)
    ev = E.Evidence({a: [5, 5], b: [0, 5]}, sample_index=0)
    assert E.combinations(a, ev, {(100, 400)}) == [(a,)]
    ev.sample_index = 1
    assert (a, b) in E.combinations(a, ev, {(100, 400)})


def test_slc7a5_minus_strand_paired_exon_has_177_bases():
    left, right = (87841155, 87843643), (87843819, 87851724)
    backbone = [(87841000, 87841155), (87851724, 87851900)]
    ev = E.Evidence({left: [8, 0], right: [9, 0]})
    assert E.paired_exon(left, right, {(87841155, 87851724)}, ev)
    ed = E.local_edit(backbone, [left, right], set(backbone))
    assert ed == [backbone[0], (87843643, 87843819), backbone[1]]
    assert ed[1][1] - ed[1][0] + 1 == 177


def test_mpzl1_skipping_backbone_cannot_swallow_upstream_exons():
    sparse = [(167721950, 167722242), (167765583, 167765749)]
    local = [(167721950, 167722242), (167734295, 167734986),
             (167741512, 167741725), (167742473, 167742605), (167765583, 167765749)]
    annotated = set(local)
    target = (167744987, 167765583)
    assert E.local_edit(sparse, [target], annotated) is None
    ed = E.local_edit(local, [target], annotated)
    assert ed[:-2] == local[:-2]
    assert ed[-2] == (167742473, 167744987)


def test_slc24a4_preserves_upstream_exons_when_inserting_141_base_exon():
    sparse = [(92322581, 92322635), (92325868, 92325978)]
    complete = [sparse[0], (92323193, 92323265), (92323808, 92323960), sparse[1]]
    pair = [(92323960, 92324060), (92324200, 92325868)]
    assert E.local_edit(sparse, pair, set(complete)) is None
    ed = E.local_edit(complete, pair, set(complete))
    assert ed == complete[:-1] + [(92324060, 92324200), complete[-1]]


def test_mpzl1_inserted_block_cannot_fuse_introns_four_and_five():
    ch = [(167721950, 167722242), (167765583, 167765749)]
    annotation = {(167734295, 167734986), (167741512, 167741725),
                  (167742473, 167742605)} | set(ch)
    hosts = {(167734986, 167741512), (167741725, 167742473)}
    wrong = [(167722242, 167733157), (167744987, 167765583)]
    assert E.local_edit(ch, wrong, annotation, host_introns=hosts) is None
    right = [(167722242, 167744903), (167744987, 167765583)]
    assert E.local_edit(ch, right, annotation, host_introns=hosts) == [
        ch[0], (167744903, 167744987), ch[1]]


def test_long_unbounded_extension_needs_another_measured_boundary():
    ch = [(100, 200), (3000, 3100)]
    assert E.local_edit(ch, [(1500, 3000)], set(ch), max_extension_nt=500) is None
    assert E.local_edit(ch, [(200, 1400), (1500, 3000)], set(ch), max_extension_nt=500) == [
        ch[0], (1400, 1500), ch[1]]


def test_long_novel_exon_needs_recurrence_and_same_host_intron():
    a, b = (100, 200), (900, 1000)
    ev = E.Evidence({a: [5, 0, 0], b: [5, 0, 0]})
    assert not E.paired_exon(a, b, {(100, 1000)}, ev)
    ev = E.Evidence({a: [5, 5, 5], b: [5, 5, 5]})
    assert E.paired_exon(a, b, {(100, 1000)}, ev)
    assert not E.paired_exon(a, b, {(100, 400), (600, 1000)}, ev)


def test_exon_skipping_is_not_mistaken_for_collateral_extension():
    ch = [(100, 200), (300, 400), (500, 600)]
    assert E.local_edit(ch, [(200, 500)], set(ch)) == [ch[0], ch[-1]]


def test_retention_marker_is_contiguous_not_a_fake_splice_junction():
    ch = [(100, 200), (300, 400)]
    ed = E.local_edit(ch, [(200, 201)], set(ch))
    assert ed == [(100, 400)]
    assert E.contains(ed, (200, 201))
    assert E.introns(ed) == ()


@pytest.mark.parametrize('counts', [ {(1, 2): [1, 2], (3, 4): [1]},
                                    {(1, 2): [-1, 2]}, {(1, 2): [np.nan]} ])
def test_invalid_evidence_is_rejected(counts):
    with pytest.raises(ValueError):
        E.Evidence(counts)


def test_retention_marker_cannot_supply_the_other_arm_of_a_novel_exon():
    a, b = (100, 101), (200, 400)
    ev = E.Evidence({a: [5, 5, 5], b: [5, 5, 5]})
    assert not E.paired_exon(a, b, {(100, 400)}, ev)


def test_pair_can_skip_known_exons_but_inner_sites_must_share_a_host():
    a, b = (100, 620), (700, 1000)
    ev = E.Evidence({a: [5], b: [5]})
    assert E.paired_exon(a, b, {(600, 800)}, ev)
    assert not E.paired_exon(a, b, {(600, 650), (680, 800)}, ev)


def test_existing_prediction_api_uses_the_supplied_learned_ranker():
    from altanalyze3.components.snaf.surface import predict_isoform as P
    dna = list('CAA' * 1000)
    dna[99:102] = 'ATG'
    dna[399:402] = 'ATG'
    dna[2101:2104] = 'TAA'
    dna = ''.join(dna)
    fetch = lambda s, e: dna[s - 1:e]
    models = {'T': ('+', [(100, 900), (2000, 2500)])}
    ev = E.Evidence({(900, 2000): [10]})
    class PreferShorter:
        def predict(self, x):
            return -np.asarray(x)[:, 0]
    candidates = P.predict_isoform(models, 900, 2000, fetch, evidence=ev, return_all=True)
    selected = P.predict_isoform(models, 900, 2000, fetch, evidence=ev,
                                 ranking_model=PreferShorter())
    assert len(candidates) >= 2
    assert len(selected.protein) == min(len(c.protein) for c in candidates)
    assert len(selected.protein) < max(len(c.protein) for c in candidates)


def test_annotation_loader_keeps_historical_mrna_exons(tmp_path):
    p = tmp_path / 'exons.tsv'
    p.write_text('gene\texon\tchr\tstrand\tstart\tend\tsource\n'
                 'G\tE4.1\tchr1\t+\t100\t200\tAF092424\n'
                 'G\tI4.1\tchr1\t+\t201\t299\t\n')
    annotation = E.load_exon_annotation(p, {'G'})
    assert annotation['G']['exons'] == {(100, 200)}
    assert annotation['G']['introns'] == {(200, 300)}
