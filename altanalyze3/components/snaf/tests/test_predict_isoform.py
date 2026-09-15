"""Contract for long-read-free isoform prediction.

Every sequence and coordinate here is synthetic. No cohort junction, isoform or gene is used.
"""
import random

import pytest

from altanalyze3.components.snaf.surface import predict_isoform as P
from altanalyze3.components.snaf.surface.orf_finder import transcript2orf



def build_locus(length=8000):
    """A locus whose background holds no T, so it contains no stop codon and no start codon.

    Features are then placed deliberately, which makes every reading frame in these tests
    predictable rather than a property of a random sequence.
    """
    seq = ('AAG' * (length // 3 + 2))[:length]
    return seq


def put(seq, pos, text):
    return seq[:pos - 1] + text + seq[pos - 1 + len(text):]


def genomic_of_offset(chain, minus, offset):
    n = 0
    for s, e in (reversed(chain) if minus else chain):
        span = e - s + 1
        if offset < n + span:
            k = offset - n
            return (e - k) if minus else (s + k)
        n += span
    raise IndexError(offset)


def place_orf(seq, chain, minus, start_offset, codons):
    """Put a start codon at a transcript offset and an in-frame stop `codons` later.

    The offsets are transcript coordinates of the EDITED chain, so a test can place a reading
    frame that deliberately spans, or deliberately misses, the junction.
    """
    seq = put(seq, genomic_of_offset(chain, minus, start_offset), 'ATG')
    seq = put(seq, genomic_of_offset(chain, minus, start_offset + 3 * codons), 'TAA')
    return seq


def fetcher(seq):
    def fetch(start, end):
        return seq[start - 1:end]
    return fetch


def make_locus(seed=0, length=4000):
    """A reproducible pseudo-random forward-strand sequence, 1-based coordinates."""
    rng = random.Random(seed)
    seq = ''.join(rng.choice('ACGT') for _ in range(length))
    def fetch(start, end):
        return seq[start - 1:end]
    return seq, fetch


def place(seq, pos, text):
    return seq[:pos - 1] + text + seq[pos - 1 + len(text):]


# ---------------------------------------------------------------- editing

def test_exon_skip():
    chain = [(100, 200), (300, 400), (500, 600)]
    assert P.apply_junction(chain, 200, 500) == [(100, 200), (500, 600)]


def test_novel_splice_site_inside_exons():
    chain = [(100, 200), (300, 400)]
    assert P.apply_junction(chain, 150, 350) == [(100, 150), (350, 400)]


def test_intron_retention_marker_fuses_the_flanking_exons():
    chain = [(100, 200), (300, 400)]
    edited = P.apply_junction(chain, 200, 201)
    assert edited == [(100, 200), (201, 400)]
    assert sum(e - s + 1 for s, e in edited) == 301


def test_transcript_that_cannot_reach_the_junction_is_rejected():
    assert P.apply_junction([(100, 200)], 300, 400) is None


def test_two_junctions_apply_together():
    chain = [(100, 200), (300, 400), (500, 600), (700, 800)]
    both = P.apply_junctions(chain, [(200, 500), (600, 700)])
    assert both == [(100, 200), (500, 600), (700, 800)]


def test_edit_labels():
    chain = [(100, 200), (300, 400), (500, 600)]
    assert P.describe_edit(chain, P.apply_junction(chain, 200, 500), 200, 500) == 'skip(1 exon)'
    assert P.describe_edit(chain, P.apply_junction(chain, 200, 201), 200, 201) == \
        'intron retention'
    assert 'novel donor in exon' in P.describe_edit(chain, P.apply_junction(chain, 150, 300),
                                                    150, 300)


# ---------------------------------------------------------------- ORF enumeration

@pytest.mark.parametrize('seed', [1, 2, 3, 4, 5])
def test_enumerate_orfs_matches_the_shipped_caller(seed):
    seq, _ = make_locus(seed=seed, length=3000)
    assert P.enumerate_orfs(seq) == transcript2orf(seq)


def test_enumerate_orfs_handles_a_sequence_with_no_stop_codon():
    cdna = 'ATG' + 'AAG' * 40
    assert P.enumerate_orfs(cdna) == transcript2orf(cdna)


# ---------------------------------------------------------------- NMD

def test_nmd_fires_when_the_stop_is_far_upstream_of_the_last_junction():
    chain = [(1, 300), (400, 700), (800, 1000)]
    assert P.nmd_by_50nt_rule(orf_start=0, orf_len=30, chain=chain, minus=False) is True


def test_nmd_does_not_fire_in_the_last_exon():
    chain = [(1, 300), (400, 700), (800, 1000)]
    assert P.nmd_by_50nt_rule(orf_start=0, orf_len=900, chain=chain, minus=False) is False


def test_single_exon_transcript_is_never_nmd():
    assert P.nmd_by_50nt_rule(0, 30, [(1, 900)], minus=False) is False


# ---------------------------------------------------------------- prediction

def test_reading_frame_must_span_the_junction():
    """Every returned candidate carries coding sequence on both sides of the junction."""
    seq = build_locus()
    chain = [(100, 900), (2000, 3000)]
    edited = P.apply_junction(chain, 900, 2000)
    seq = place_orf(seq, edited, False, 600, 200)                   # spans the junction at 800
    seq = put(seq, 2600, 'ATG')                                     # a start AFTER the junction
    fetch = fetcher(seq)
    cands = P.predict_isoform({'T1': ('+', chain)}, 900, 2000, fetch, return_all=True)
    assert cands
    for cand in cands:
        offsets = P.transcript_offsets(cand.chain, minus=False)
        mrna = ''.join(fetch(s, e) for s, e in cand.chain).upper()
        pos = mrna.find(cand.orf)
        assert pos <= offsets[900], 'an ORF starting after the junction must be rejected'
        assert pos + len(cand.orf) > offsets[2000], 'an ORF ending before it must be rejected'


def test_an_nmd_isoform_is_never_assigned():
    """A transcript nonsense-mediated decay destroys presents nothing, so it is not a prediction."""
    seq = build_locus()
    chain = [(100, 900), (2000, 2500), (4000, 4500)]
    edited = P.apply_junction(chain, 900, 2000)
    # opens at 600, spans the junction at 800, stops at 900. The last exon-exon junction of the
    # three-exon transcript sits at 1302, so the stop is 399 nt upstream of it.
    seq = place_orf(seq, edited, False, 600, 100)
    fetch = fetcher(seq)
    cands = P.predict_isoform({'T1': ('+', chain)}, 900, 2000, fetch, return_all=True)
    assert cands and all(c.nmd for c in cands), 'this fixture must produce only NMD candidates'
    assert P.predict_isoform({'T1': ('+', chain)}, 900, 2000, fetch) is None


def test_edited_isoform_is_named_after_its_parent_and_junctions():
    seq = build_locus()
    chain = [(100, 900), (2000, 2500), (4000, 4500)]
    edited = P.apply_junction(chain, 900, 2000)
    seq = place_orf(seq, edited, False, 600, 200)
    fetch = fetcher(seq)
    cands = P.predict_isoform({'PARENT1': ('+', chain)}, 900, 2000, fetch,
                              co_junctions=[(2500, 4000)], junction_label='J1',
                              co_junction_labels={(2500, 4000): 'J2'}, return_all=True)
    ids = {c.isoform_id for c in cands}
    assert 'PARENT1|J1' in ids
    assert 'PARENT1|J1+J2' in ids
    for c in cands:
        assert c.isoform_id.startswith('PARENT1|'), 'parentage must stay in the identifier'


def test_full_mrna_is_exposed_and_the_cds_indexes_into_it():
    seq = build_locus()
    chain = [(100, 900), (2000, 2500)]
    edited = P.apply_junction(chain, 900, 2000)
    # exon 1 is 801 nt, so the junction sits at transcript offset 800. The reading frame opens
    # at 600, leaving a 5' untranslated region, and closes at 900, past the junction.
    seq = place_orf(seq, edited, False, 600, 100)
    fetch = fetcher(seq)
    cands = P.predict_isoform({'T1': ('+', chain)}, 900, 2000, fetch, return_all=True)
    assert cands
    for c in cands:
        assert len(c.mrna) == c.transcript_len
        assert len(c.mrna) > len(c.orf), 'the mRNA must carry more than the coding sequence'
        assert c.mrna[c.cds_start:c.cds_end] == c.orf


def test_tsv_carries_ids_sequences_and_a_row_for_every_junction(tmp_path):
    seq = build_locus()
    chain = [(100, 900), (2000, 2500)]
    edited = P.apply_junction(chain, 900, 2000)
    seq = place_orf(seq, edited, False, 600, 200)
    fetch = fetcher(seq)
    pred = P.predict_isoform({'PARENT1': ('+', chain)}, 900, 2000, fetch, junction_label='J1')
    out = str(tmp_path / 'predicted.tsv')
    P.write_predictions_tsv(out, [('J1', pred), ('J2', None)])
    lines = open(out).read().rstrip('\n').split('\n')
    header = lines[0].split('\t')
    assert len(lines) == 3, 'a junction with no prediction still gets a row'
    for col in ('predicted_isoform_id', 'parent_isoform_id', 'mrna_sequence', 'cds_sequence',
                'protein_sequence', 'junctions_applied', 'nmd'):
        assert col in header
    row = dict(zip(header, lines[1].split('\t')))
    assert row['parent_isoform_id'] == 'PARENT1'
    assert row['predicted_isoform_id'] == 'PARENT1|J1'
    assert row['nmd'] == 'no'
    assert len(row['mrna_sequence']) > len(row['cds_sequence']) > 0


def test_all_reference_isoforms_are_edited_not_one():
    seq, fetch = make_locus(seed=3, length=5000)
    models = {
        'T1': ('+', [(100, 900), (2000, 2500)]),
        'T2': ('+', [(100, 900), (2000, 2500), (3000, 3500)]),
        'T3': ('+', [(50, 900), (2000, 2200)]),
    }
    # min_orf_nt=3 isolates the property under test: which isoforms get edited, not how long
    # their reading frames are. The length rule has its own test.
    cands = P.predict_isoform(models, 900, 2000, fetch, return_all=True, min_orf_nt=3)
    assert {c.backbone for c in cands} == {'T1', 'T2', 'T3'}, \
        'every reference isoform must be edited and offered'


def test_co_applied_junction_reaches_an_isoform_one_edit_cannot():
    """A second novel junction of the same gene must widen the candidate set."""
    seq = build_locus()
    chain = [(100, 900), (2000, 2500), (4000, 4500)]
    edited_one = P.apply_junction(chain, 900, 2000)
    seq = place_orf(seq, edited_one, False, 600, 200)
    fetch = fetcher(seq)
    models = {'T1': ('+', [(100, 900), (2000, 2500), (4000, 4500)])}
    one = P.predict_isoform(models, 900, 2000, fetch, return_all=True)
    two = P.predict_isoform(models, 900, 2000, fetch, co_junctions=[(2500, 4000)],
                            return_all=True)
    assert one, 'the single-junction edit must yield at least one spanning ORF'
    assert {c.co_applied for c in one} == {0}
    assert len(two) > len(one), 'the co-applied junction must add candidates'
    assert 1 in {c.co_applied for c in two}


def test_minus_strand_prediction_uses_the_reverse_complement():
    """On the minus strand the transcript is read from the other end."""
    seq = build_locus()
    chain = [(100, 900), (2000, 2500)]
    edited = P.apply_junction(chain, 900, 2000)
    # On the minus strand the 5' exon is the one with the higher coordinates, so the junction
    # sits at transcript offset 500. The reading frame opens at 300 and closes at 600, spanning it.
    g_start = genomic_of_offset(edited, True, 300)
    seq = put(seq, g_start - 2, 'CAT')                    # reverse complement of ATG
    g_stop = genomic_of_offset(edited, True, 600)
    seq = put(seq, g_stop - 2, 'TTA')                     # reverse complement of TAA
    fetch = fetcher(seq)
    minus = P.predict_isoform({'T1': ('-', chain)}, 900, 2000, fetch, return_all=True)
    plus = P.predict_isoform({'T1': ('+', chain)}, 900, 2000, fetch, return_all=True)
    assert minus, 'the minus-strand model must read the reverse complement'
    assert {c.protein for c in minus} != {c.protein for c in plus}


def test_no_candidate_returns_none():
    _, fetch = make_locus(seed=13, length=2000)
    assert P.predict_isoform({'T1': ('+', [(100, 200)])}, 900, 1000, fetch) is None


def test_predict_isoforms_offers_same_gene_junctions_to_each_other():
    seq, fetch = make_locus(seed=17, length=6000)
    models = {'G1': {'T1': ('+', [(100, 900), (2000, 2500), (4000, 4500)])}}
    junctions = [('J1', 'G1', 900, 2000), ('J2', 'G1', 2500, 4000)]
    out = P.predict_isoforms(junctions, models, fetch)
    assert set(out) == {'J1', 'J2'}


def test_cache_does_not_change_the_answer():
    """The per-gene cache is an optimization, so it must return exactly the same isoforms."""
    seq = build_locus()
    chain = [(100, 900), (2000, 2500), (4000, 4500)]
    edited = P.apply_junction(chain, 900, 2000)
    seq = place_orf(seq, edited, False, 600, 200)
    fetch = fetcher(seq)
    models = {'T1': ('+', chain), 'T2': ('+', [(100, 900), (2000, 2500)])}
    plain = P.predict_isoform(models, 900, 2000, fetch, co_junctions=[(2500, 4000)],
                              junction_label='J1', return_all=True)
    cache = {}
    cached = P.predict_isoform(models, 900, 2000, fetch, co_junctions=[(2500, 4000)],
                               junction_label='J1', return_all=True, cache=cache)
    assert cache, 'the cache must actually be populated'
    assert len(plain) == len(cached)
    for a, b in zip(plain, cached):
        assert (a.isoform_id, a.protein, a.orf, a.chain, a.nmd, a.mrna, a.cds_start) == \
               (b.isoform_id, b.protein, b.orf, b.chain, b.nmd, b.mrna, b.cds_start)
    one = P.predict_isoform(models, 900, 2000, fetch, junction_label='J1')
    two = P.predict_isoform(models, 900, 2000, fetch, junction_label='J1', cache={})
    assert (one is None) == (two is None)
    if one is not None:
        assert one.isoform_id == two.isoform_id and one.protein == two.protein


def test_cache_is_reused_across_junctions_of_one_gene():
    seq = build_locus()
    chain = [(100, 900), (2000, 2500), (4000, 4500)]
    edited = P.apply_junction(chain, 900, 2000)
    seq = place_orf(seq, edited, False, 600, 200)
    fetch = fetcher(seq)
    models = {'T1': ('+', chain)}
    cache = {}
    P.predict_isoform(models, 900, 2000, fetch, junction_label='J1', cache=cache)
    first = len(cache)
    P.predict_isoform(models, 900, 2000, fetch, junction_label='J1', cache=cache)
    assert len(cache) == first, 'a repeated chain must not add a new entry'


def test_short_reading_frames_are_refused_like_snaf_refuses_them():
    """orf_finder.prioritize_orf requires 30 codons; enumerating frames must not bypass that."""
    from altanalyze3.components.snaf.surface.orf_finder import prioritize_orf
    assert P.MIN_ORF_NT == 30 * 3
    seq = build_locus()
    chain = [(100, 900), (2000, 2500)]
    edited = P.apply_junction(chain, 900, 2000)
    seq = place_orf(seq, edited, False, 780, 8)          # 8 codons, spans the junction at 800
    fetch = fetcher(seq)
    assert P.predict_isoform({'T1': ('+', chain)}, 900, 2000, fetch) is None
    kept = P.predict_isoform({'T1': ('+', chain)}, 900, 2000, fetch, min_orf_nt=3)
    assert kept is not None and len(kept.protein) == 8
    assert prioritize_orf([kept.orf]) == '', 'SNAF would also refuse this reading frame'
