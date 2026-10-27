from altanalyze3.components.psi.event_matching import event_junctions, junction_matches

G = 'ENSG00000154134'


def test_single_opposed_arm_match_inverts_direction_and_retains_subexons():
    a = event_junctions(f'ROBO3:{G}:E17.12-E17.21|E17.17-E17.21')
    b = event_junctions(f'{G}:E17.2-E17.3|{G}:E17.12-E17.21')
    m = junction_matches(a, b)
    assert len(m) == 1 and m[0]['direction_multiplier'] == -1
    assert (m[0]['left_arm'], m[0]['right_arm']) == (1, 2)
    # User's two illustrated cell-line events share no complete junction.
    for s in ['E17.13-E17.14|E17.12-E17.14', 'E17.14-E17.15|E17.14-E17.16']:
        assert junction_matches(a, event_junctions(f'{G}:' + s)) == []


def test_all_four_arm_pairs_and_reversed_event_are_retained():
    a = event_junctions(f'{G}:E1.1-E2.1|{G}:E1.1-E3.1')
    same = junction_matches(a, a)
    reverse = junction_matches(a, tuple(reversed(a)))
    assert [(m['left_arm'], m['right_arm'], m['direction_multiplier']) for m in same] == [(1, 1, 1), (2, 2, 1)]
    assert [(m['left_arm'], m['right_arm'], m['direction_multiplier']) for m in reverse] == [(1, 2, -1), (2, 1, -1)]


def test_gene_intronic_and_novel_coordinate_specificity():
    a = event_junctions(f'{G}:E1.1-I1.1_123|{G}:E2.1-E3.1')
    for s in [f'{G}:E1.1-I1.1_124|E7.1-E8.1',
              'ENSG00000000001:E1.1-I1.1_123|E2.1-E3.1']:
        assert not junction_matches(a, event_junctions(s))
    trans = event_junctions(f'{G}:E1.1-ENSG00000000001:E2.1|{G}:E3.1-E4.1')
    assert trans[0] == (G + ':E1.1', 'ENSG00000000001:E2.1')
