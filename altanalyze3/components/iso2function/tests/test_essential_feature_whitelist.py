"""Essential feature classes (TM, signal, kinase, DNA-binding, ...) are always drawn, labelled by
class, nearest the backbone and above anything they overlap. Rule set by Nathan Salomonis, 2026-10-23."""
import re

import pytest

from altanalyze3.components.iso2function.plotting import protein_features as pf
from altanalyze3.components.iso2function.plotting.isv_architecture import select_tracks, layout_tracks
from altanalyze3.components.iso2function.protein_isoform_audit import align_mapping, project_feature

REF = 'MACDEFGHIKLMNPQRSTVWYAVCLDPEFQGRHSIKTWLMYNVFDAPEGKRLHISCTQWMNDYKAVERFGPSLHITCQNW'  # 80 aa, no repeated block


def row(iso, seq, start, end, provider, name, description, fid, accession='', ref='P'):
    p = project_feature(REF, seq, start, end, align_mapping(REF, seq))
    return dict(gene_symbol='G', isoform_id=iso, provider=provider, feature_id=fid, name=name, accession=accession,
                description=description, reference_protein=ref, start=str(start), end=str(end), evidence='',
                sequence_present=str(REF[start - 1:end] in seq), **p)


def pair_rows(start, end, provider, name, description, fid, accession='', ref='P', alt=None):
    alt = REF[:60] if alt is None else alt
    return [row('A', REF, start, end, provider, name, description, fid, accession, ref),
            row('B', alt, start, end, provider, name, description, fid, accession, ref)]


def test_whitelist_file_parses_and_covers_the_required_classes():
    wl = pf.load_essential_whitelist()
    classes = {r['cls'] for r in wl['rules']}
    assert {'TM', 'SIGNAL', 'KINASE', 'PHOSPHATASE', 'DNA_BINDING'} <= classes
    assert all(re.fullmatch(r'#[0-9A-F]{6}', r['color'], re.I) for r in wl['rules'])
    assert min(r['priority'] for r in wl['rules'] if r['cls'] == 'TM') == 1


def test_uniprot_transmembrane_rows_are_matched_and_pfam_family_is_not():
    wl = pf.load_essential_whitelist()
    tm = dict(provider='uniprot_features', name='UniProt', accession='', description='TRANSMEM-Helica')
    assert pf.essential_rule(tm, wl)['cls'] == 'TM'
    kin = dict(provider='protein_features', name='PF00069', accession='IPR000719', description='Prot_kinase_dom')
    assert pf.essential_rule(kin, wl)['cls'] == 'KINASE'
    dbd = dict(provider='uniprot_features', name='UniProt', accession='', description='DNA_BIND-Homeobo')
    assert pf.essential_rule(dbd, wl)['cls'] == 'DNA_BINDING'
    fam = dict(provider='protein_features', name='PF13903', accession='IPR004031', description='PMP22/EMP/MP20/Claudin')
    assert pf.essential_rule(fam, wl) is None
    helix = dict(provider='uniprot_features', name='UniProt', accession='', description='HELIX-')
    assert pf.essential_rule(helix, wl) is None


def test_named_signaling_motifs_survive_pfam_only_selection_and_remain_distinct():
    rows = pair_rows(7, 10, 'uniprot_features', 'UniProt', 'MOTIF-ITA', 'uniprot_features:ita')
    rows += pair_rows(30, 35, 'uniprot_features', 'UniProt', 'MOTIF-ITIM', 'uniprot_features:itim')
    pair = dict(gene_symbol='G', isoform_A='A', isoform_B='B')
    features, audit, _ = pf.feature_plan(
        pair, rows, dict(A=REF, B=REF[:60]),
        feature_types=('PF00059',), sources=('protein_features',))
    assert {f['essential_class'] for f in features} == {'NAMED_MOTIF:ita', 'NAMED_MOTIF:itim'}
    assert {f['class_label'] for f in features} == {'ITA', 'ITIM'}
    assert all(r['display_role'] == 'displayed' for r in audit if r['source_isoform'] == 'A')
    assert all(f['A_start'] != '' and f['B_native_start'] != '' for f in features)
    assert pf.essential_rule(dict(provider='uniprot_features', description='MOTIF-'),
                             pf.load_essential_whitelist()) is None


def test_compact_tiny_coding_segment_retains_its_width_and_id_in_key():
    rows = pair_rows(7, 10, 'uniprot_features', 'UniProt', 'MOTIF-ITA', 'uniprot_features:ita')
    pair = dict(gene_symbol='G', isoform_A='A', isoform_B='B')
    features, _, events = pf.feature_plan(pair, rows, dict(A=REF, B=REF[:60]))
    models = {i: dict(assigned_reads='', nmd='Not-NMD', reported_aa=n,
                     coding_exons=[dict(aa_start0=1., aa_end0=1.01, color='#E04F60',
                                        exon_label='E12', exon_label_status='assigned')])
              for i, n in [('A',80), ('B',60)]}
    fig, _ = pf.render_compact_pair(pair, features, events, dict(A=REF,B=REF[:60]), models, 'context')
    ax = fig.axes[0]
    assert any(p.get_gid()=='coding_exon_A' and p.get_width()==pytest.approx(.01) for p in ax.patches)
    assert any('Small coding segments' in t.get_text() and 'E12 [1,1.01)' in t.get_text() for t in ax.texts)


def test_transmembrane_helices_are_displayed_despite_source_and_type_selection_and_numbered():
    rows = []
    for i, (s, e) in enumerate([(5, 25), (41, 61), (66, 80)], 1):
        rows += pair_rows(s, e, 'uniprot_features', 'UniProt', 'TRANSMEM-Helica', 'uniprot_features:' + str(i), ref='P' + str(i))
    rows += pair_rows(1, 70, 'protein_features', 'PF13903', 'PMP22/EMP/MP20/Claudin', 'protein_features:9', 'IPR004031')
    pair = dict(gene_symbol='G', isoform_A='A', isoform_B='B')
    feats, audit, _ = pf.feature_plan(pair, rows, dict(A=REF, B=REF[:60]), feature_types=('PF13903',), sources=('protein_features',))
    tm = [f for f in feats if f['essential_class'] == 'TM']
    assert [f['class_display'] for f in sorted(tm, key=lambda f: f['A_start'])] == ['TM1', 'TM2', 'TM3']
    assert all(pf.architecture_label(f).startswith('TM') for f in tm)
    assert {e['display_role'] for e in audit if e['essential_class'] == 'TM' and e['source_isoform'] == 'A'} == {'displayed'}
    # the family hit is still drawn, but the TM helices take the first lane and are drawn last
    fam = [f for f in feats if f['feature_type'] == 'PF13903']
    assert len(fam) == 1 and fam[0]['essential_class'] == ''
    placed, _ = pf.architecture_features(feats, 'A')
    assert placed


def test_transmembrane_is_displayed_even_when_only_reviewed_sources_are_selected():
    rows = pair_rows(5, 25, 'uniprot_features', 'UniProt', 'TRANSMEM-Helica', 'uniprot_features:1')
    pair = dict(gene_symbol='G', isoform_A='A', isoform_B='B')
    feats, audit, _ = pf.feature_plan(pair, rows, dict(A=REF, B=REF[:60]))
    assert len(feats) == 1 and feats[0]['class_display'] == 'TM'
    assert audit[0]['display_role'] == 'displayed' and audit[0]['essential_class'] == 'TM'


def test_overlapping_same_annotation_from_several_reference_proteins_collapses_to_one_bar():
    rows = []
    for i, (s, e) in enumerate([(3, 70), (5, 72), (2, 68)], 1):
        rows += pair_rows(s, e, 'protein_features', 'PF13903', 'PMP22/EMP/MP20/Claudin', 'protein_features:' + str(i), 'IPR004031', ref='P' + str(i))
    rows += pair_rows(75, 80, 'protein_features', 'PF13903', 'PMP22/EMP/MP20/Claudin', 'protein_features:4', 'IPR004031', ref='P4')
    pair = dict(gene_symbol='G', isoform_A='A', isoform_B='B')
    feats, audit, _ = pf.feature_plan(pair, rows, dict(A=REF, B=REF[:60]), feature_types=('PF13903',), sources=('protein_features',))
    assert len(feats) == 2                       # one overlapping chain + one separate interval
    rep = next(f for f in feats if f['A_start'] <= 5)
    assert (rep['A_start'], rep['A_end']) == (3, 70)   # the longest source interval, not a union
    assert sorted(rep['source_feature_ids']) == ['protein_features:1', 'protein_features:2', 'protein_features:3']
    roles = {e['feature_id']: e['display_role'] for e in audit if e['source_isoform'] == 'A'}
    assert roles['protein_features:1'] == 'displayed'
    assert roles['protein_features:2'] == roles['protein_features:3'] == 'coalesced_overlapping_same_annotation'
    assert roles['protein_features:4'] == 'displayed'


def test_compact_renderer_puts_essential_features_in_the_first_lane_and_on_top(tmp_path):
    rows = pair_rows(5, 25, 'uniprot_features', 'UniProt', 'TRANSMEM-Helica', 'uniprot_features:1')
    rows += pair_rows(1, 70, 'protein_features', 'PF13903', 'PMP22/EMP/MP20/Claudin', 'protein_features:9', 'IPR004031')
    pair = dict(gene_symbol='G', isoform_A='A', isoform_B='B')
    feats, _, events = pf.feature_plan(pair, rows, dict(A=REF, B=REF[:60]), feature_types=('PF13903',), sources=('protein_features',))
    models = {i: dict(assigned_reads='', nmd='Not-NMD', reported_aa=n) for i, n in (('A', 80), ('B', 60))}
    fig, layout = pf.render_compact_pair(pair, feats, events, dict(A=REF, B=REF[:60]), models, 'context')
    ax = fig.axes[0]
    tm_bars = [p for p in ax.patches if (p.get_gid() or '').startswith('feature_A_') and p.get_facecolor()[:3] == pytest.approx((0xD5 / 255, 0x5E / 255, 0))]
    fam_bars = [p for p in ax.patches if (p.get_gid() or '').startswith('feature_A_') and p not in tm_bars]
    assert tm_bars and fam_bars
    assert tm_bars[0].get_zorder() > fam_bars[0].get_zorder()
    assert tm_bars[0].get_y() < fam_bars[0].get_y()          # nearest the backbone
    labels = [t.get_text() for t in ax.texts if (t.get_gid() or '').startswith('feature_label_A_')]
    assert 'TM' in labels
    fig.savefig(tmp_path / 'x.pdf')


def test_isv_tracks_carry_class_labels_colours_and_lane_precedence():
    rows = pair_rows(5, 25, 'uniprot_features', 'UniProt', 'TRANSMEM-Helica', 'uniprot_features:1')
    rows += pair_rows(1, 70, 'protein_features', 'PF13903', 'PMP22/EMP/MP20/Claudin', 'protein_features:9', 'IPR004031')
    pair = dict(gene_symbol='G', isoform_A='A', isoform_B='B')
    feats, _, _ = pf.feature_plan(pair, rows, dict(A=REF, B=REF[:60]), feature_types=('PF13903',), sources=('protein_features',))
    models = {}
    for iid, n in (('A', 80), ('B', 60)):
        models[iid] = dict(gene_symbol='G', protein_aa=n, reported_aa=str(n), protein_sha256='x', source_nmd='Not-NMD',
                           coding_exons=[dict(start=0.0, end=float(n), color='#e04f60', label='E1')])
    bundle = dict(plans=[dict(pair=pair, features=feats)], models=models)
    tracks = select_tracks(bundle, ['A', 'B'], 'G')
    tm = next(f for f in tracks[0]['features'] if f['essential_class'] == 'TM')
    assert tm['label'] == 'TM' and tm['color'] == '#D55E00'
    layout = layout_tracks(tracks, 600)
    rects = [p for p in layout['primitives'] if p['kind'] == 'rect' and p.get('role') == 'protein_feature' and p['isoform_id'] == 'A']
    tm_rect = next(p for p in rects if p.get('essential'))
    fam_rect = next(p for p in rects if not p.get('essential'))
    assert rects.index(tm_rect) > rects.index(fam_rect)     # drawn after, so on top
    assert tm_rect['y'] > fam_rect['y']                      # lane nearest the backbone (y grows downward)
    assert any(p['kind'] == 'text' and p['text'] == 'TM' and p.get('weight') == 'bold' for p in layout['primitives'])


def test_overlapping_rows_of_one_class_from_several_providers_become_one_bar():
    rows = pair_rows(10, 60, 'protein_features', 'PF00069', 'Prot_kinase_dom', 'protein_features:1', 'IPR000719')
    rows += pair_rows(12, 58, 'uniprot_features', 'UniProt', 'DOMAIN-Protein kinas', 'uniprot_features:2')
    rows += pair_rows(70, 80, 'uniprot_features', 'UniProt', 'TRANSMEM-Helica', 'uniprot_features:3')
    pair = dict(gene_symbol='G', isoform_A='A', isoform_B='B')
    feats, audit, _ = pf.feature_plan(pair, rows, dict(A=REF, B=REF[:60]), feature_types=('x',), sources=('protein_features',))
    kin = [f for f in feats if f['essential_class'] == 'KINASE']
    # one bar; the curated UniProt interval wins the tie over the Pfam row
    assert len(kin) == 1 and (kin[0]['A_start'], kin[0]['A_end']) == (12, 58) and kin[0]['class_display'] == 'Kinase'
    assert 'protein_features:Prot_kinase_dom' in kin[0]['coalesced_providers']
    roles = {e['feature_id']: e['display_role'] for e in audit if e['source_isoform'] == 'A'}
    assert roles['protein_features:1'] == 'coalesced_same_essential_class'
    assert [f['class_display'] for f in feats if f['essential_class'] == 'TM'] == ['TM']


def test_two_zinc_fingers_survive_a_spanning_family_row():
    rows = pair_rows(10, 30, 'uniprot_features', 'UniProt', 'ZN_FING-GATA-typ', 'uniprot_features:1')
    rows += pair_rows(40, 60, 'uniprot_features', 'UniProt', 'ZN_FING-GATA-typ', 'uniprot_features:2')
    rows += pair_rows(8, 62, 'protein_features', 'PF00320', 'Znf_GATA', 'protein_features:3', 'IPR000679')
    pair = dict(gene_symbol='G', isoform_A='A', isoform_B='B')
    feats, audit, _ = pf.feature_plan(pair, rows, dict(A=REF, B=REF[:60]), feature_types=('x',), sources=('protein_features',))
    dbd = sorted((f for f in feats if f['essential_class'] == 'DNA_BINDING'), key=lambda f: f['A_start'])
    assert [(f['A_start'], f['A_end'], f['class_display']) for f in dbd] == [(10, 30, 'DNA-binding1'), (40, 60, 'DNA-binding2')]
    roles = {e['feature_id']: e['display_role'] for e in audit if e['source_isoform'] == 'A'}
    assert roles['protein_features:3'] == 'coalesced_same_essential_class'


def test_named_native_domains_survive_pfam_only_filters_and_reach_isv():
    # These are distinct, overlapping domains, not interchangeable members of
    # one generic "Domain" class. Neither the source nor reference filter may
    # silently remove them from an architecture figure.
    rows = pair_rows(5, 35, 'uniprot_features', '', 'DOMAIN-Alpha', 'native:1')
    rows += pair_rows(20, 50, 'uniprot_features', '', 'DOMAIN-Beta', 'native:2')
    pair = dict(gene_symbol='G', isoform_A='A', isoform_B='B')
    feats, audit, _ = pf.feature_plan(pair, rows, dict(A=REF, B=REF[:60]),
                                     feature_types=('PF00000',), sources=('protein_features',),
                                     reference_proteins={'OTHER'})
    assert {f['class_label'] for f in feats} == {'Alpha', 'Beta'}
    assert len({f['essential_class'] for f in feats}) == 2
    assert {e['display_role'] for e in audit if e['source_isoform'] == 'A'} == {'displayed'}
    models = {iid: dict(gene_symbol='G', protein_aa=n, reported_aa=str(n),
                       coding_exons=[dict(start=0., end=float(n), color='#e04f60', label='E1')])
              for iid, n in [('A', len(REF)), ('B', 60)]}
    tracks = select_tracks(dict(plans=[dict(pair=pair, features=feats)], models=models), ['A', 'B'])
    layout = layout_tracks(tracks, 1000)
    assert {'Alpha', 'Beta'} <= {p.get('text') for p in layout['primitives']}
    assert len([p for p in layout['primitives'] if p.get('role') == 'protein_feature']) == 4


def test_named_domain_fallback_keeps_literal_labels_and_does_not_promote_arbitrary_regions():
    wl = pf.load_essential_whitelist()
    assert pf.essential_rule(dict(provider='reviewed_uniprot_snapshot', name='Domain',
                                  description='Novel domain'), wl)['label'] == 'Novel domain'
    assert pf.essential_rule(dict(provider='uniprot_features', description='DOMAIN-Truncate'), wl)['label'] == 'Truncate'
    assert pf.essential_rule(dict(provider='uniprot_features', description='REGION-Disordered'), wl) is None
    assert pf.essential_rule(dict(provider='protein_features', name='Unrelated', description='Domain'), wl) is None


def test_fbar_wh2_and_separate_ef_hands_have_native_specific_labels():
    wl = pf.load_essential_whitelist()
    assert pf.essential_rule(dict(provider='protein_features', accession='IPR031160'), wl)['label'] == 'F-BAR'
    assert pf.essential_rule(dict(provider='protein_features', accession='IPR003124'), wl)['label'] == 'WH2'
    assert pf.essential_rule(dict(provider='uniprot_features', description='DOMAIN-WH'), wl)['label'] == 'WH'
    rows = pair_rows(5, 20, 'uniprot_features', '', 'DOMAIN-EF-hand', 'native:1')
    rows += pair_rows(30, 45, 'uniprot_features', '', 'DOMAIN-EF-hand', 'native:2')
    feats, _, _ = pf.feature_plan(dict(gene_symbol='G', isoform_A='A', isoform_B='B'),
                                 rows, dict(A=REF, B=REF[:60]), feature_types=('NONE',),
                                 sources=('protein_features',))
    assert [(f['A_start'], f['A_end'], f['class_display']) for f in feats] == [(5, 20, 'EF-hand1'), (30, 45, 'EF-hand2')]
