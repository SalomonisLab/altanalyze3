"""Collapse-side terminal U evidence; no upstream annotation changes."""
import pytest
from altanalyze3.components.long_read.isoform_collapse.splice_evidence import (
    Policy, prepare_gene, terminal_u_events, recover_selected_terminal_tokens)


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('side', ['left', 'right'])
def test_recurrent_u_splice_pools_variable_rna_ends(strand, side):
    regions = {'E1.1': (300, 400), 'I1.1': (401, 499), 'E2.1': (500, 600)}
    rows = []
    for n in range(3):
        exons = [(100+n, 150), (300, 400)] if side == 'left' else [(500, 600), (700, 750+n)]
        raw = f'U1.1_{100+n}|U1.1_150|E1.1|E1.1_400' if side == 'left' else f'E2.1_500|E2.1|U2.1_700|U2.1_{750+n}'
        if strand == '-': raw = '|'.join(reversed(raw.split('|')))
        rows.append(dict(library=str(n), gene='G', molecule='m', strand=strand,
                         legacy_structure='U1.1|E1.1' if side == 'left' else 'E2.1|U2.1',
                         raw_structure=raw, source_exons=exons))
    def provider(row):
        return terminal_u_events(row['raw_structure'], row['source_exons'], regions, 'G', 3)
    mapping, candidates = prepare_gene(rows, [], Policy(3, 100, 3),
        terminal_event_provider=provider,
        annotation_corrector=lambda r,e,p:recover_selected_terminal_tokens(r,e,p,regions))
    assert len(set(mapping.values())) == 1
    candidate = next(iter(candidates.values()))
    assert candidate['count'] == 3 and len(candidate['protected_junctions']) == 1
    # Recurrence is molecule support, not an RNA-end classification.
    _, rare = prepare_gene(rows[:2], [], Policy(3, 100, 3), terminal_event_provider=provider)
    assert not next(iter(rare.values()))['protected_junctions']


def test_u_requires_own_gene_actual_splice_not_endpoint_or_interior():
    regions = {'E1.1': (300, 400), 'E2.1': (500, 600)}
    assert not terminal_u_events('U1.1_100', [(100,150),(300,400)], regions, 'G', 3)
    assert not terminal_u_events('OTHER:U1.1_150', [(100,150),(300,400)], regions, 'G', 3)
    assert not terminal_u_events('U1.1_400', [(300,400),(500,600)], regions, 'G', 3)
    assert terminal_u_events('G:U1.1_153', [(100,150),(300,400)], regions, 'G', 3)
    assert not terminal_u_events('G:U1.1_154', [(100,150),(300,400)], regions, 'G', 3)
    assert not terminal_u_events('U1.1_150', [(100,150)], regions, 'G', 3)


def test_u_recovery_rejects_unmeasured_label():
    row=dict(library='A',gene='G',molecule='m',strand='+',legacy_structure='E1.1',
             raw_structure='U1.1_100|E1.1',source_exons=[(100,150),(300,400)])
    with pytest.raises(ValueError,match='original raw/primary-alignment'):
        recover_selected_terminal_tokens(row,[('exon','U1.1',(150,300),(100,150))],
                                         Policy(3,100,3),{'E1.1':(300,400)})
