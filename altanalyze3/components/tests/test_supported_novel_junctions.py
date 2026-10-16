"""Physical junction identities must survive normalization of E-coordinate tokens."""
import pytest
from altanalyze3.components.long_read import isoform_collapse_utils as U
from altanalyze3.components.long_read.isoform_collapse import splice_evidence as S
from altanalyze3.components.long_read.isoform_collapse import scored_collapse as C
from altanalyze3.components.long_read.isoform_collapse.source_stream import (
    build_event_index, measured_splice_events)


def runtime():
    result = {k:getattr(U,k) for k in ('structure_tokens_for_containment',
        'junction_core_tokens','is_contiguous_subsequence','intronic_edge_compatible')}
    return dict(result,collapse_gene_jc=C.collapse_gene_jc)


@pytest.mark.parametrize('shift',[(6,0),(0,6),(6,6)])
@pytest.mark.parametrize('reverse',[False,True])
def test_nearby_sites_do_not_pool_different_known_reference_compatibility(shift,reverse):
    """A known-compatible minority RNA must not be forced into a novel class."""
    left,right=shift
    group=[dict(library='A',gene='G',molecule=str(i),legacy_structure='E1.1|E2.1',
                source_exons=[(1,100+(left//2 if i==3 else 0)),
                              (200+(right//2 if i==3 else 0),300)]) for i in range(4)]
    if reverse:group.reverse()
    refs={'E1.1|E2.1':dict(id='ENST_KNOWN',exons=[(1,100+left),(200+right,300)])}
    mapping,candidates=S.prepare_gene(group,[],S.Policy(3,100,3),
        parent_models={key:ref['exons'] for key,ref in refs.items()},
        protect_splice_geometry=True,
        required_identities=[(r['library'],'G',r['molecule']) for r in group],
        terminal_event_provider=lambda r:measured_splice_events(r['source_exons']))
    result,final=S.collapse_prepared(candidates,refs,'G',runtime(),splice_slop=3)
    assert result['enst'][final[mapping['A','G','3']]]=='ENST_KNOWN'
    assert all(final[mapping['A','G',str(i)]] not in result['enst'] for i in range(3))
    assert len(mapping)==4 and sum(c['count'] for c in candidates.values())==4


@pytest.mark.parametrize('offset',[0,6])
def test_uniform_reference_compatibility_preserves_existing_physical_classes(offset):
    """Uniform eligibility adds no measured contrast or reason to relabel a class."""
    group=[dict(library='A',gene='G',molecule=str(i),legacy_structure='E1.1|E2.1',
                source_exons=[(1,100),(200,300)]) for i in range(3)]
    def run(parents):
        return S.prepare_gene(group,[],S.Policy(3,100,3),parent_models=parents,
            protect_splice_geometry=True,
            terminal_event_provider=lambda r:measured_splice_events(r['source_exons']))
    assert run({})==run({'known':[(1,100+offset),(200+offset,300)]})


def rows(counts=(8,3,1)):
    result=[]
    for acceptor,count in zip((513,643,700),counts):
        for i in range(count):
            result.append(dict(library='A',gene='G',molecule=f'{acceptor}_{i}',
                legacy_structure=f'E2.1|E3.1_{acceptor}|E3.1|E3.2|E4.1',
                source_exons=[(100+i,296),(acceptor,1093),(1637,1800+i)]))
    return result


def prepare(group):
    index=build_event_index({'E2.1':(100,296),'E3.1':(400,1093),'E4.1':(1637,1800)},
                            {'ref':[(100,296),(513,1093),(1637,1800)]})
    return S.prepare_gene(group,[],S.Policy(3,100,3),
        required_identities=[(r['library'],'G',r['molecule']) for r in group],
        terminal_event_provider=lambda r:measured_splice_events(r['source_exons']))


def test_supported_exonic_acceptor_survives_while_rare_site_keeps_native_fallback():
    group=rows();mapping,candidates=prepare(group)
    _,final=S.collapse_prepared(candidates,{},'G',runtime(),splice_slop=3)
    parent={a:final[mapping['A','G',f'{a}_0']] for a in (513,643,700)}
    assert parent[513]!=parent[643]
    assert parent[700]==parent[513]
    assert len(mapping)==12 and sum(c['count'] for c in candidates.values())==12
    assert (296,643) in candidates[mapping['A','G','643_0']]['protected_junctions']


def test_known_preference_cannot_override_supported_alternative_acceptor():
    mapping,candidates=prepare(rows((8,3,0)))
    ref={'E2.1|E3.1|E3.2|E4.1':dict(id='ENST1',exons=[(100,296),(513,1093),(1637,1800)])}
    result,final=S.collapse_prepared(candidates,ref,'G',runtime(),splice_slop=3)
    assert final[mapping['A','G','513_0']] in result['enst']
    assert final[mapping['A','G','643_0']] not in result['enst']


def test_duplicate_physical_source_aliases_do_not_manufacture_recurrence():
    group=[]
    for alias in ('A','B','C'):
        group.append(dict(library=alias,gene='G',molecule='same',
            support_source=('same_raw_file','G','same'),legacy_structure='E1.1|E2.1',
            source_exons=[(1,100),(300,400)]))
    _,candidates=prepare(group)
    assert sum(c['count'] for c in candidates.values())==3
    assert all(not c['protected_junctions'] for c in candidates.values())


def test_actual_splice_within_one_exon_block_is_not_lost():
    index=build_event_index({'E3.1':(100,900)}, {})
    assert measured_splice_events([(100,300),(500,900)])==[
        ('splice','physical_splice',(300,500),(300,500))]
    assert measured_splice_events([(100,900)])==[]


@pytest.mark.parametrize('offset,novel',[(0,False),(3,False),(4,True),(-3,False),(-4,True)])
def test_reference_compatibility_uses_exactly_existing_three_nt(offset,novel):
    group=[dict(library='A',gene='G',molecule=str(i),
        legacy_structure=f'E1.1|E2.1_{500+offset}|E2.1',
        source_exons=[(100,300),(500+offset,900)]) for i in range(3)]
    mapping,candidates=S.prepare_gene(group,[],S.Policy(3,100,3),
        terminal_event_provider=lambda r:measured_splice_events(r['source_exons']))
    refs={'E1.1|E2.1':dict(id='ENST1',exons=[(100,300),(500,900)])}
    result,final=S.collapse_prepared(candidates,refs,'G',runtime(),splice_slop=3)
    assert (final[mapping['A','G','0']] not in result['enst'])==novel


def test_terminal_annotation_recovery_receives_no_invented_region():
    seen=[]
    group=rows((0,3,0));index=build_event_index({}, {})
    def correct(row,events,policy):
        seen.extend(events)
        return row['legacy_structure']
    S.prepare_gene(group,[],S.Policy(3,100,3),annotation_corrector=correct,
        terminal_event_provider=lambda r:measured_splice_events(r['source_exons']))
    assert seen==[]


def test_public_pipeline_prepares_and_persists_complete_physical_evidence(tmp_path):
    from altanalyze3.components.tests.test_intronic_collapse_integration import source_fixture
    from altanalyze3.components.long_read.isoform_collapse import pipeline as P
    import json
    samples,annot,cache,refgff=source_fixture(tmp_path)
    out=tmp_path/'physical_pipeline'
    result=P.run_pipeline(samples,str(out),nproc=1,min_total=3,ref=str(annot),
        enst_cache=str(cache),ref_gff=str(refgff),write_h5ad=False,
        tier1_dir=str(tmp_path/'physical_tier1'),
        protect_intronic_evidence=True,protect_novel_junctions=True)
    assert sorted((row[1],row[3]) for row in result['catalog'])==[('p0.A',3),('r0.A',4)]
    assert len(P.load_molecule_candidates(out,'A'))==7
    scan=json.loads((out/'intronic_source_evidence/sample_0000/manifest.json').read_text())
    assert scan['complete'] and scan['protect_novel_junctions']
    assert scan['novel_junction_scope']=='all_physical_junctions_and_continuity_v4'
    assert scan['native_eligible_reads']==7
    assert scan['possible_evidence_reads']==7


def test_public_pipeline_rejects_novel_guard_without_complete_source_preparation(tmp_path):
    from altanalyze3.components.long_read.isoform_collapse import pipeline as P
    with pytest.raises(ValueError,match='complete source-evidence preparation'):
        P.run_pipeline([],str(tmp_path/'bad'),protect_novel_junctions=True)


def test_gencode_known_site_inside_older_native_region_requires_physical_guard():
    refs={
        'E1.1|E2.1_500|E2.1|E2.2|E3.1':dict(id='ENST_A',exons=[(20,100),(500,800),(1000,1200)]),
        'E1.1|E2.1_630|E2.1|E2.2|E3.1':dict(id='ENST_B',exons=[(20,100),(630,800),(1000,1200)])}
    index=build_event_index({'E1.1':(20,100),'E2.1':(400,800),'E3.1':(1000,1200)},
                            {k:r['exons'] for k,r in refs.items()})
    group=[dict(library='A',gene='G',molecule=str(i),
        legacy_structure='E1.0|E1.1|E2.1_630|E2.1|E2.2|E3.1|E3.2',
        source_exons=[(10,100),(630,800),(1000,1300)]) for i in range(3)]
    mapping,candidates=S.prepare_gene(group,[],S.Policy(3,100,3),
        terminal_event_provider=lambda r:measured_splice_events(r['source_exons']))
    result,final=S.collapse_prepared(candidates,refs,'G',runtime(),splice_slop=3)
    assert (100,630) in next(iter(candidates.values()))['protected_junctions']
    assert result['enst'][final[mapping['A','G','0']]]=='ENST_B'


def test_canonical_source_cannot_be_redirected_to_coordinate_normalized_alternative_reference():
    refs={
        'E1.1|E2.1|E2.2|E3.1':dict(id='ENST_CANONICAL',exons=[(20,100),(400,800),(1000,1200)]),
        'E1.1|E2.1_630|E2.1|E2.2|E3.1':dict(id='ENST_ALTERNATIVE',exons=[(20,100),(630,800),(1000,1200)])}
    group=[dict(library='A',gene='G',molecule=str(i),
        legacy_structure='E1.0|E1.1|E2.1|E2.2|E3.1|E3.2',
        source_exons=[(10,100),(400,800),(1000,1300)]) for i in range(3)]
    mapping,candidates=S.prepare_gene(group,[],S.Policy(3,100,3),
        terminal_event_provider=lambda r:measured_splice_events(r['source_exons']))
    result,final=S.collapse_prepared(candidates,refs,'G',runtime(),splice_slop=3)
    assert result['enst'][final[mapping['A','G','0']]]=='ENST_CANONICAL'


def continuity_rows(spliced=8,continuous=3):
    result=[]
    for label,count,exons in (
            ('spliced',spliced,[(20,100),(400,500),(600,800),(1000,1200)]),
            ('continuous',continuous,[(20,100),(400,800),(1000,1200)])):
        for i in range(count):
            result.append(dict(library='A',gene='G',molecule=f'{label}_{i}',
                legacy_structure=('E1.0|E1.1|E2.1_500|E2.1_600|E2.1|E2.2|E3.1|E3.2'
                    if label=='spliced' else 'E1.0|E1.1|E2.1|E2.2|E3.1|E3.2'),
                source_exons=exons))
    return result


def prepare_continuity(group,refs):
    return S.prepare_gene(group,[],S.Policy(3,100,3),
        parent_models={k:r['exons'] for k,r in refs.items()},protect_splice_geometry=True,
        terminal_event_provider=lambda r:measured_splice_events(r['source_exons']))


def test_continuous_canonical_exon_cannot_acquire_an_internal_reference_splice():
    refs={
        'E1.1|E2.1|E2.2|E3.1':dict(id='ENST_CANONICAL',exons=[(20,100),(400,800),(1000,1200)]),
        'E1.1|E2.1_500|E2.1_600|E2.1|E2.2|E3.1':dict(id='ENST_SPLICED',exons=[(20,100),(400,500),(600,800),(1000,1200)])}
    mapping,candidates=prepare_continuity(continuity_rows(0,3),refs)
    c=candidates[mapping['A','G','continuous_0']]
    assert (500,600) in c['protected_splice_spans']
    result,final=S.collapse_prepared(candidates,refs,'G',runtime(),splice_slop=3)
    assert result['enst'][final[mapping['A','G','continuous_0']]]=='ENST_CANONICAL'


def test_recurrent_internal_splice_cannot_absorb_positive_continuity():
    mapping,candidates=prepare_continuity(continuity_rows(),{})
    _,final=S.collapse_prepared(candidates,{},'G',runtime(),splice_slop=3)
    assert final[mapping['A','G','continuous_0']]!=final[mapping['A','G','spliced_0']]
    assert (500,600) in candidates[mapping['A','G','continuous_0']]['protected_splice_spans']


def test_unrecurrent_novel_gap_does_not_create_continuity_classes():
    mapping,candidates=prepare_continuity(continuity_rows(2,8),{})
    assert not candidates[mapping['A','G','continuous_0']]['protected_splice_spans']


def test_continuity_is_a_splice_contradiction_not_an_endpoint_constraint():
    assert S.splice_spans_compatible({(500,600)},[(400,550)],3)
    assert S.splice_spans_compatible({(500,600)},[(650,800)],3)
    assert S.splice_spans_compatible({(500,600)},[(400,800)],3)
    assert not S.splice_spans_compatible({(500,600)},[(400,500),(600,800)],3)
    assert not S.splice_spans_compatible({(500,600)},[(400,503),(603,800)],3)
    # Spanning one nested reference gap does not establish both flanks of a
    # different larger parent gap. Partial terminal coverage keeps its 100-nt
    # rule; it must not acquire a new event through interval overlap alone.
    assert S.splice_spans_compatible({(500,600)},[(100,400),(800,1000)],3)


def test_positive_continuity_requires_both_source_flanks():
    refs={'E2.1':dict(id='ENST1',exons=[(400,500),(600,800)])}
    group=[dict(library='A',gene='G',molecule=str(i),legacy_structure='E2.1',
        source_exons=[exons]) for i,exons in enumerate(((400,550),(550,800),(497,603)))]
    mapping,candidates=prepare_continuity(group,refs)
    assert candidates[mapping['A','G','2']]['protected_splice_spans']=={(500,600)}
    assert not candidates[mapping['A','G','0']]['protected_splice_spans']
    assert not candidates[mapping['A','G','1']]['protected_splice_spans']


def test_eventless_monoexon_geometry_is_protected_before_exact_aggregation(tmp_path):
    from altanalyze3.components.tests.test_intronic_collapse_integration import source_fixture
    from altanalyze3.components.long_read.isoform_collapse import source_stream as SS,stream_context as SC
    samples,annot,cache,refgff=source_fixture(tmp_path)
    # Keep all seven original identities, but give them continuous source
    # geometry across a reference splice too short for terminal continuation.
    sample=samples[0]
    from pathlib import Path
    for p in (Path(sample[2]),tmp_path/'transcript_associations_raw.txt'):
        rows=[line.split('\t') for line in p.read_text().splitlines()]
        for row in rows:row[2]='E1.1'
        p.write_text(''.join('\t'.join(row)+'\n' for row in rows))
    Path(sample[3]).write_text(''.join(
        f'chr1\tsource\texon\t1\t100\t.\t+\t.\tgene_id "G"; transcript_id "{m}";\n'
        for m in ['p0','p1','p2','r0','r1','r2','r3']))
    regions={'G':{'E1.1':(1,100)}}
    refs={'G':{'E1.1':dict(id='ENST1',exons=[(1,50),(70,100)],strand='+')}}
    p=S.Policy(3,100,3)
    scan=SS.scan_sample(sample,regions,{'G':SS.build_event_index(regions['G'],{'ref':refs['G']['E1.1']['exons']})},tmp_path/'scan',p,protect_novel_junctions=True)
    assert scan['native_eligible_reads']==scan['possible_evidence_reads']==7
    result=SC.build_contexts([scan],regions,refs,tmp_path/'contexts',p)
    ctx=SC.load_context(result['gene_evidence']['G'])
    assert sum(c['count'] for c in ctx['candidates'].values())==7
    assert all(c['protected_splice_spans']=={(50,70)} for c in ctx['candidates'].values())


def test_decoded_value_sharing_preserves_exact_values_and_is_immutable():
    import json,zlib
    from altanalyze3.components.bam.coord_store import pack_exons
    from altanalyze3.components.long_read.isoform_collapse.stream_context import make_source_decoders
    text,raw,exons,events=make_source_decoders()
    value='E1.1|E2.1|E3.1'
    assert text(value)==value and text(value) is text(value.encode().decode())
    blob=zlib.compress(value.encode())
    assert raw(blob)==value and raw(blob) is raw(bytes(bytearray(blob)))
    geometry=pack_exons([(1,100),(300,400)])
    assert exons(geometry)==((1,100),(300,400)) and exons(geometry) is exons(geometry)
    encoded=json.dumps([['splice','physical_splice',[100,300],[100,300]]])
    result=events(encoded)
    assert result==(('splice','physical_splice',(100,300),(100,300)),)
    assert result is events(encoded.encode().decode())
    with pytest.raises(TypeError):result[0][2][0]=101
    with pytest.raises(TypeError):exons(geometry)[0][0]=2


def test_frozen_terminal_projection_cannot_reuse_physical_event_cache():
    import json
    from altanalyze3.components.long_read.isoform_collapse.stream_context import make_source_decoders
    encoded=json.dumps([['splice','physical_splice',[100,300],[100,300]],
                        ['continuation','I1.1',None,[100,300]]])
    full=make_source_decoders()[-1]
    terminal=make_source_decoders({'continuation'})[-1]
    assert len(full(encoded))==2
    assert terminal(encoded)==(('continuation','I1.1',None,(100,300)),)


def test_shared_class_values_keep_every_required_molecule_assignment(tmp_path):
    from altanalyze3.components.long_read.isoform_collapse import pipeline as P
    expected={('G','a'):('E1.1|E2.1','E1.1|E2.1~intronic=abc'),
              ('G','b'):('E1.1|E2.1','E1.1|E2.1~intronic=abc')}
    P.persist_molecule_candidates(tmp_path,{'A':expected})
    actual=P.load_molecule_candidates(tmp_path,'A')
    assert actual==expected and len(actual)==2
    assert actual['G','a'] is actual['G','b']
