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


def test_tiny_exon_labels_fail_legibly_instead_of_disappearing(tmp_path):
    tracks=select_tracks(bundle(tmp_path),['A','B']);tracks[0]['coding_exons'][0]['end']=.01
    with pytest.raises(ValueError,match='cannot fit'):layout_tracks(tracks,1000)


def test_unplaced_annotations_are_explicit_and_never_called_absent(tmp_path):
    b=bundle(tmp_path);b['plans'][0]['features'][0].update(A_start='',A_end='',status='ambiguous_placement')
    tracks=select_tracks(b,['A','B']);layout=layout_tracks(tracks,1000)
    assert all(t['unplaced']==['F1'] and not t['features'] for t in tracks)
    assert len([p for p in layout['primitives'] if p.get('text','').startswith('Unplaced annotations:')])==2
