"""ISV consumes pinned domain placements without changing native identities."""
from copy import deepcopy
import csv
import hashlib
import json
import pytest
from altanalyze3.components.iso2function.plotting.isv_architecture import load_architectures, select_tracks, layout_tracks


def bundle(tmp_path):
    feature=dict(feature_type='Domain', description='DNA binding', provider='reviewed_uniprot_snapshot',
                 feature_id='P:1', reference_protein='P', display_id='F1', A_start=21, A_end=160,
                 B_native_segments=[[21, 100]], status='partial')
    plans=[dict(pair=dict(gene_symbol='G', isoform_A='A', isoform_B='B'),features=[feature])]
    (tmp_path/'plot_data.json').write_text(json.dumps(plans))
    rows=[]
    for iid,n in [('A',200),('B',100)]:
        for start,end,label,color in [(0,50,'E1','#e04f60'),(50,n,'E2','#4488bb')]:
            rows.append(dict(gene_symbol='G',isoform_id=iid,protein_aa=n,reported_aa=n,protein_sha256=iid*64,
                             source_nmd='Not-NMD',aa_start0=start,aa_end0=end,color=color,
                             exon_label=label,exon_label_status='assigned'))
    with (tmp_path/'coding_exon_coordinates.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    (tmp_path/'annotation_display_audit.tsv').write_text('complete source audit\n')
    manifest=dict(status='complete',plot_style='compact',coding_exon_map={'path':'pinned-map'},
                  source_manifest_sha256='f'*64,selection={},outputs={})
    for name in ['plot_data.json','coding_exon_coordinates.tsv','annotation_display_audit.tsv']:
        manifest['outputs'][name]=hashlib.sha256((tmp_path/name).read_bytes()).hexdigest()
    (tmp_path/'figure_manifest.json').write_text(json.dumps(manifest))
    return load_architectures(tmp_path)


def test_complete_native_roster_and_partial_domains_keep_original_coordinates(tmp_path):
    b=bundle(tmp_path);before=deepcopy(b);tracks=select_tracks(b,['B','A'],'G')
    assert [t['isoform_id'] for t in tracks]==['B','A'] and b==before
    assert tracks[0]['features'][0]['segments']==[(21,100)]
    assert tracks[0]['features'][0]['partial'] and not tracks[1]['features'][0]['partial']
    layout=layout_tracks(tracks,1000)
    exons=[p for p in layout['primitives'] if p.get('role')=='coding_exon']
    assert sum(p['width'] for p in exons if p['isoform_id']=='B')==500
    assert sum(p['width'] for p in exons if p['isoform_id']=='A')==1000
    assert any(p.get('text')=='DNA binding (partial)' for p in layout['primitives'])
    assert len([p for p in layout['primitives'] if p.get('role')=='exon_label'])==4


def test_unmatched_and_ambiguous_native_models_fail_without_shrinking_selection(tmp_path):
    b=bundle(tmp_path)
    with pytest.raises(ValueError,match='No audited'):select_tracks(b,['A','MISSING'])
    with pytest.raises(ValueError,match='unique'):select_tracks(b,['A','A'])
    with pytest.raises(ValueError,match='different gene'):select_tracks(b,['A'],'WRONG')
    alternate=deepcopy(b['plans'][0]);alternate['features'][0]['A_end']=150;b['plans'].append(alternate)
    with pytest.raises(ValueError,match='ambiguous'):select_tracks(b,['A','B'])


def test_modified_cache_or_unresolved_coding_exon_is_rejected(tmp_path):
    bundle(tmp_path)
    (tmp_path/'plot_data.json').write_text('[]')
    with pytest.raises(ValueError,match='unpinned'):load_architectures(tmp_path)


def test_tiny_exon_keeps_true_width_and_literal_label_in_a_legible_key(tmp_path):
    tracks=select_tracks(bundle(tmp_path),['A','B']);tracks[0]['coding_exons'][0]['end']=.01
    layout=layout_tracks(tracks,1000)
    block=next(p for p in layout['primitives'] if p.get('role')=='coding_exon' and p['isoform_id']=='A')
    assert block['width']==pytest.approx(.01*1000/max(t['protein_aa'] for t in tracks))
    assert any(p.get('role')=='small_exon_key' and tracks[0]['coding_exons'][0]['label'] in p['text']
               for p in layout['primitives'])


def test_unplaced_annotations_are_explicit_and_never_called_absent(tmp_path):
    b=bundle(tmp_path);b['plans'][0]['features'][0].update(A_start='',A_end='',status='ambiguous_placement')
    tracks=select_tracks(b,['A','B']);layout=layout_tracks(tracks,1000)
    assert all(t['unplaced']==['F1'] and not t['features'] for t in tracks)
    assert len([p for p in layout['primitives'] if p.get('text','').startswith('Unplaced annotations:')])==2
    assert not any('Not in this protein' in p.get('text','') for p in layout['primitives'])


def essential_bundle(tmp_path, placeable):
    """A panel whose only annotation is an essential DNA-binding class, placed on A or not placeable at all."""
    b=bundle(tmp_path)
    f=b['plans'][0]['features'][0]
    f.update(essential_class='DNA_BINDING',class_label='DNA-binding',class_display='DNA-binding',
             class_color='#1F78B4',class_priority='1',reference_start=117,reference_end=139)
    if not placeable:
        f.update(A_start='',A_end='',B_native_segments=[],status='unplaced_annotation')
    return b


def test_unplaceable_essential_class_is_drawn_on_a_reference_row(tmp_path):
    """Nathan, 2026-10-25: IKZF1 Pro-B-1 showed no DNA-binding bar because no isoform could place it."""
    tracks=select_tracks(essential_bundle(tmp_path,placeable=False),['A','B'],'G')
    assert [t['isoform_id'] for t in tracks]==['A','B','P']
    reference=tracks[-1]
    assert reference['reference'] and reference['coding_exons']==[] and reference['protein_aa']==139
    assert [f['segments'] for f in reference['features']]==[[(117,139)]]
    assert reference['features'][0]['essential_class']=='DNA_BINDING'
    layout=layout_tracks(tracks,1000)
    bars=[p for p in layout['primitives'] if p.get('role')=='protein_feature' and p.get('essential')]
    assert len(bars)==1 and bars[0]['isoform_id']=='P'
    assert any(p.get('role')=='reference_backbone' for p in layout['primitives'])
    assert any(p.get('text')=='DNA-binding' for p in layout['primitives'])
    assert any(p.get('text','').startswith('P (reference)') for p in layout['primitives'])


def test_placed_essential_class_adds_no_reference_row(tmp_path):
    """Panels that already draw their essential bars keep exactly the tracks they had."""
    tracks=select_tracks(essential_bundle(tmp_path,placeable=True),['A','B'],'G')
    assert [t['isoform_id'] for t in tracks]==['A','B']
    assert not any(t.get('reference') for t in tracks)


def test_one_reference_row_per_essential_class(tmp_path):
    """IKZF1 carried the same zinc fingers on ten reference proteins; one row per class is drawn."""
    from copy import deepcopy as _dc
    b=essential_bundle(tmp_path,placeable=False)
    first=b['plans'][0]['features'][0]
    second=_dc(first);second.update(feature_id='Q:1',display_id='F2',reference_protein='Q',
                                    reference_start=117,reference_end=139)
    third=_dc(first);third.update(feature_id='P:2',display_id='F3',reference_start=145,reference_end=167)
    b['plans'][0]['features']=[first,second,third]
    tracks=select_tracks(b,['A','B'],'G')
    references=[t for t in tracks if t.get('reference')]
    assert [t['isoform_id'] for t in references]==['P']
    assert sorted(f['segments'][0] for f in references[0]['features'])==[(117,139),(145,167)]
