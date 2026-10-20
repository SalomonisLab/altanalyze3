import argparse
import hashlib
import json
from pathlib import Path

import pytest

from altanalyze3.components.iso2function.cli import add_arguments
from altanalyze3.components.iso2function.protein_isoform_audit import align_mapping,project_feature
from altanalyze3.components.iso2function.plotting.protein_features import (
    runs,alignment_segments,feature_plan,selected_pairs,load_audit,render_pair,
    architecture_features,render_compact_pair,render_compact_figure,feature_reference_selection)


def rows_for(reference,alternative,start,end,feature_type='Domain'):
    rows=[]
    for name,seq in [('A',reference),('B',alternative)]:
        p=project_feature(reference,seq,start,end,align_mapping(reference,seq))
        rows.append(dict(gene_symbol='G',isoform_id=name,provider='reviewed_uniprot_snapshot',feature_id='P:1',
                         name=feature_type,description='Example feature',reference_protein='P',start=str(start),end=str(end),
                         evidence='',sequence_present=str(reference[start-1:end] in seq),**p))
    return rows


def plan(reference,alternative,start,end):
    p=dict(gene_symbol='G',isoform_A='A',isoform_B='B')
    return feature_plan(p,rows_for(reference,alternative,start,end),dict(A=reference,B=alternative))


def test_partial_domain_preserves_native_coordinates_and_remaining_residues():
    f,_,_=plan('MACDEFGHIKLMNPQRSTVWY','MACDEFGHIK',6,16)
    f=f[0]
    assert f['status']=='partial'
    assert f['identical_residues']==5 and f['absent_residues']==6
    assert f['identical_segments']==[(6,10)] and f['absent_segments']==[(11,16)]
    assert (f['B_native_start'],f['B_native_end'])==(6,10)


def test_n_terminal_truncation_is_not_left_shifted_on_reference_axis():
    reference='MACDEFGHIKLMNPQRSTVWY';alternative=reference[9:]
    f,_,events=plan(reference,alternative,11,18)
    assert f[0]['status']=='retained' and f[0]['B_native_start']==2
    assert f[0]['identical_segments']==[(11,18)]
    assert any(e['state']=='deleted' and e['A_start']==1 and e['A_end']==9 for e in events)


def test_substitution_and_insertion_are_never_called_complete_feature_loss():
    ref='MACDEFGHIKLMNPQRSTVWY'
    alt=ref[:10]+'YYY'+ref[10:]
    f,_,events=plan(ref,alt,6,16)
    assert f[0]['status']=='altered' and f[0]['inserted_residues']==3
    assert f[0]['identical_residues']==11 and not f[0]['absent_residues']
    assert any(e['state']=='inserted' and e['A_anchor']==10.5 for e in events)
    alt=ref[:10]+'Y'+ref[11:]
    f,_,_=plan(ref,alt,6,16)
    assert f[0]['status']=='altered' and f[0]['altered_residues']==1


def test_repeated_feature_is_not_a_loss_or_a_new_gain():
    f,_,_=plan('MACDEF','MACDEFACDEF',2,6)
    assert f[0]['status']=='ambiguous_placement'
    rows=rows_for('MACDEF','MACDEFACDEF',2,6)
    # Swap the plotted A/B orientation to exercise an ambiguous A placement.
    for r in rows:r['isoform_id']='B' if r['isoform_id']=='A' else 'A'
    f,_,_=feature_plan(dict(gene_symbol='G',isoform_A='A',isoform_B='B'),rows,dict(A='MACDEFACDEF',B='MACDEF'))
    assert f[0]['status']=='unplaced_annotation'


def test_feature_intact_only_in_B_stays_visible_without_fabricated_A_coordinates():
    reference='MACDEFGHIKLMNPQRSTVWY'
    rows=rows_for(reference,reference[:6],10,18)
    for r in rows:r['isoform_id']='B' if r['isoform_id']=='A' else 'A'
    f,audit,_=feature_plan(dict(gene_symbol='G',isoform_A='A',isoform_B='B'),rows,dict(A=reference[:6],B=reference))
    assert f[0]['status']=='intact_only_in_B_unplaced'
    assert f[0]['A_start']=='' and f[0]['B_native_start']==10
    assert len(audit)==2 and audit[1]['display_role']=='displayed_without_A_coordinate'


def test_display_audit_keeps_all_source_instances_and_low_read_pairs():
    rows=rows_for('MACDEFGHIKLM','MACDEF',3,5)
    other=dict(rows[0],provider='protein_features',feature_id='IPR:1')
    f,audit,_=feature_plan(dict(gene_symbol='G',isoform_A='A',isoform_B='B'),rows+[other],dict(A='MACDEFGHIKLM',B='MACDEF'))
    assert len(audit)==3 and any(a['display_role']=='source_not_in_view' for a in audit)
    pairs=[dict(gene_symbol='G',isoform_A='A',isoform_B='B',both_at_least_min_reads='False')]
    assert selected_pairs(pairs,[dict(gene_symbol='G')],'G')==pairs
    with pytest.raises(ValueError,match='Unresolved requested genes'):
        selected_pairs(pairs,[dict(gene_symbol='G')],'MISSING')


def test_hash_gate_rejects_changed_figure_inputs(tmp_path):
    for name in ['native_model_audit.tsv','pair_summary.tsv','feature_projections.tsv','native_proteins.fasta']:
        (tmp_path/name).write_text('original')
    manifest=dict(status='complete_with_explicit_annotation_limits',outputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in tmp_path.iterdir()})
    (tmp_path/'impact_manifest.json').write_text(json.dumps(manifest))
    (tmp_path/'feature_projections.tsv').write_text('changed')
    with pytest.raises(ValueError,match='Changed or unpinned'):
        load_audit(tmp_path)


def test_default_pair_selection_and_explicit_unpaired_gene():
    pairs=[dict(gene_symbol='G',isoform_A='A',isoform_B='B')]
    models=[dict(gene_symbol='G'),dict(gene_symbol='SINGLE')]
    assert selected_pairs(pairs,models)==pairs
    with pytest.raises(ValueError,match='without an audited pair'):
        selected_pairs(pairs,models,'SINGLE')


def test_pair_file_is_a_complete_explicit_selection(tmp_path):
    pairs=[dict(gene_symbol=g,isoform_A=g+'1',isoform_B=g+'2') for g in ['G','H']]
    models=[dict(gene_symbol=g) for g in ['G','H']]
    p=tmp_path/'pairs.tsv';p.write_text('gene_symbol\tisoform_A\tisoform_B\nH\tH1\tH2\n')
    assert selected_pairs(pairs,models,pair_file=p)==[pairs[1]]
    with pytest.raises(ValueError,match='not an audited native comparison'):
        selected_pairs(pairs,models,genes='G',pair_file=p)


def test_standard_renderer_has_valid_label_geometry_without_gene_rules():
    pytest.importorskip('matplotlib')
    ref='MACDEFGHIKLMNPQRSTVWY';alt=ref[:10]
    f,_,events=plan(ref,alt,6,16)
    pair=dict(gene_symbol='UNSEEN_GENE',isoform_A='A',isoform_B='B')
    # The gene label is presentation only; it does not change annotation logic.
    fig,layout=render_pair(pair,f,events,dict(A=ref,B=alt),{x:dict(assigned_reads='30',nmd='Not-NMD') for x in 'AB'},'RNA support: synthetic test fixture.')
    assert fig.axes[0].get_xlabel()=='Aligned amino-acid position in native A'
    assert layout and all(v['bounds_px'][2]>0 for v in layout)
    import matplotlib.pyplot as plt;plt.close(fig)


def test_shared_cli_registers_reproducible_plot_action():
    p=add_arguments(argparse.ArgumentParser())
    a=p.parse_args(['protein-impact-plots','--impact-results','audit','--genes','TOX,CITED2','--artifacts-dir','figures'])
    assert a.genes=='TOX,CITED2' and a.feature_type is None and a.annotation_source is None
    assert a.protein_plot_style=='compact'
    a=p.parse_args(['protein-impact-plots','--protein-plot-style','detailed'])
    assert a.protein_plot_style=='detailed'


def test_architecture_keeps_best_mapped_domain_despite_substitution():
    ref='MACDEFGHIKLMNPQRSTVWY';alt=ref[:10]+'Y'+ref[11:]
    f,_,_=plan(ref,alt,6,16)
    placed,unplaced=architecture_features(f,'B')
    assert not unplaced and len(placed)==1
    assert placed[0]['segments']==[(6,16)]
    assert not placed[0]['partial']


def test_architecture_uses_native_coordinates_after_n_terminal_shortening():
    ref='MACDEFGHIKLMNPQRSTVWY';alt=ref[9:]
    f,_,_=plan(ref,alt,11,18)
    placed,_=architecture_features(f,'B')
    assert placed[0]['segments']==[(2,9)]


def test_architecture_spans_internal_insertion_without_difference_marks():
    ref='MACDEFGHIKLMNPQRSTVWY';alt=ref[:10]+'YYY'+ref[10:]
    f,_,_=plan(ref,alt,6,16)
    placed,_=architecture_features(f,'B')
    assert placed[0]['segments']==[(6,19)] and not placed[0]['partial']


def test_architecture_retains_B_only_feature_from_text_table_coordinates():
    ref='MACDEFGHIKLMNPQRSTVWY'
    rows=rows_for(ref,ref[:6],10,18)
    for r in rows:r['isoform_id']='B' if r['isoform_id']=='A' else 'A'
    f,_,_=feature_plan(dict(gene_symbol='G',isoform_A='A',isoform_B='B'),rows,dict(A=ref[:6],B=ref))
    f[0]['B_native_start']='10';f[0]['B_native_end']='18'
    placed,unplaced=architecture_features(f,'B')
    assert not unplaced and placed[0]['segments']==[(10,18)]


def test_compact_renderer_has_two_tracks_without_caption_or_substitution_marks():
    ref='MACDEFGHIKLMNPQRSTVWY';alt=ref[:10]
    f,_,events=plan(ref,alt,6,16)
    pair=dict(gene_symbol='G',isoform_A='A',isoform_B='B')
    fig,layout=render_compact_pair(pair,f,events,dict(A=ref,B=alt),{x:dict(nmd='Not-NMD') for x in 'AB'},'SHOULD NOT APPEAR')
    texts=[t.get_text() for t in fig.texts+fig.axes[0].texts]
    assert not any('SHOULD NOT APPEAR' in t or 'AA change' in t or 'identical' in t for t in texts)
    assert 'G-1' in texts and 'G-2' in texts and any('(partial)' in t for t in texts)
    bars={p.get_gid():p for p in fig.axes[0].patches if p.get_gid()}
    assert bars['feature_A_F1'].get_facecolor()==bars['feature_B_F1'].get_facecolor()
    assert bars['feature_B_F1'].get_x()==5.5 and bars['feature_B_F1'].get_width()==5
    assert fig.get_figheight()<4 and layout
    import matplotlib.pyplot as plt;plt.close(fig)


def test_compact_coding_exon_blocks_use_the_cached_native_lengths_and_colors():
    ref='MACDEFGHIKLMNPQRSTVWY';alt=ref[:10]
    f,_,events=plan(ref,alt,6,16)
    pair=dict(gene_symbol='G',isoform_A='A',isoform_B='B')
    models={name:dict(nmd='Not-NMD',coding_exons=[dict(aa_start0=0,aa_end0=len(seq),color='#cc6699')])
            for name,seq in [('A',ref),('B',alt)]}
    fig,_=render_compact_pair(pair,f,events,dict(A=ref,B=alt),models,'')
    blocks=[p for p in fig.axes[0].patches if (p.get_gid() or '').startswith('coding_exon_')]
    assert len(blocks)==2
    assert [p.get_width() for p in blocks]==[len(ref),len(alt)]
    assert blocks[0].get_facecolor()==blocks[1].get_facecolor()
    assert all(p.get_x()==.5 for p in blocks)
    import matplotlib.pyplot as plt;plt.close(fig)


def test_compact_exon_labels_are_inside_blocks_and_keep_reference_numbers():
    ref='MACDEFGHIKLMNPQRSTVWY';alt=ref[:10]
    f,_,events=plan(ref,alt,6,16)
    pair=dict(gene_symbol='G',isoform_A='A',isoform_B='B')
    models={name:dict(nmd='Not-NMD',coding_exons=[dict(aa_start0=0,aa_end0=len(seq),color='#cc6699',
                         exon_group=1,exon_label='E9',exon_label_status='assigned')])
            for name,seq in [('A',ref),('B',alt)]}
    fig,_=render_compact_pair(pair,f,events,dict(A=ref,B=alt),models,'')
    ax=fig.axes[0];renderer=fig.canvas.get_renderer()
    labels=[t for t in ax.texts if (t.get_gid() or '').startswith('exon_label_')]
    blocks=[p for p in ax.patches if (p.get_gid() or '').startswith('coding_exon_')]
    assert len(labels)==2 and all(t.get_text()=='E9' for t in labels)
    for label,block in zip(labels,blocks):
        a=label.get_window_extent(renderer);b=block.get_window_extent(renderer)
        assert b.contains(a.x0,a.y0) and b.contains(a.x1,a.y1)
    import matplotlib.pyplot as plt;plt.close(fig)


def test_tiny_coding_fragment_has_legible_rotated_inside_label():
    pair=dict(gene_symbol='G',isoform_A='A',isoform_B='B')
    def block(lo,hi):return dict(aa_start0=lo,aa_end0=hi,color='#cc6699',exon_label='E4',exon_label_status='assigned')
    models=dict(A=dict(nmd='Not-NMD',coding_exons=[block(0,270)]),
                B=dict(nmd='Not-NMD',coding_exons=[block(0,392/3),block(392/3,135)]))
    fig,_=render_compact_pair(pair,[],[],dict(A='M'*270,B='M'*135),models,'')
    label=next(t for t in fig.axes[0].texts if t.get_gid()=='exon_label_B_2')
    assert label.get_rotation()==90 and label.get_fontsize()>=6.5
    models['B']['coding_exons'][1]['exon_label_status']='unresolved_reference_coverage'
    with pytest.raises(ValueError,match='Unresolved AltAnalyze exon label'):
        render_compact_pair(pair,[],[],dict(A='M'*270,B='M'*135),models,'')
    import matplotlib.pyplot as plt;plt.close('all')


def test_combined_figure_keeps_complete_vector_panels_on_one_page(tmp_path):
    ref='MACDEFGHIKLMNPQRSTVWY';alt=ref[:10]
    f,_,events=plan(ref,alt,6,16)
    panels=[dict(pair=dict(gene_symbol=g,isoform_A='A',isoform_B='B'),features=f,events=events)
            for g in ['FIRST','SECOND']]
    models={i:dict(nmd='Not-NMD') for i in ('A','B')}
    fig,layout=render_compact_figure(panels,dict(A=ref,B=alt),models,'')
    assert len(fig.subfigs)==len(layout)==2
    assert [p['gene_symbol'] for p in layout]==['FIRST','SECOND']
    assert fig.get_figwidth()==12 and fig.get_figheight()==pytest.approx(5.7)
    first,second=[p['bounds_px'] for p in layout]
    assert second[1]+second[3]<=first[1]+1
    for sf in fig.subfigs:
        assert len([p for p in sf.axes[0].patches if (p.get_gid() or '').startswith('feature_')])==2
        assert len([t for t in sf.axes[0].texts if (t.get_gid() or '').startswith('native_id_')])==2
    from matplotlib.backends.backend_pdf import PdfPages
    path=tmp_path/'combined.pdf'
    with PdfPages(path) as pdf:
        pdf.savefig(fig)
        assert pdf.get_pagecount()==1
    assert path.stat().st_size>1000
    with pytest.raises(ValueError,match='one positive height'):
        render_compact_figure(panels,dict(A=ref,B=alt),models,'',[1])
    import matplotlib.pyplot as plt;plt.close('all')


def test_standard_width_option_fits_a_six_residue_exon_without_distorting_geometry():
    def block(lo,hi,label):return dict(aa_start0=lo,aa_end0=hi,color='#cc6699',exon_label=label,exon_label_status='assigned')
    pair=dict(gene_symbol='G',isoform_A='A',isoform_B='B')
    models=dict(A=dict(nmd='Not-NMD',coding_exons=[block(0,6,'E16'),block(6,479,'E20')]),
                B=dict(nmd='Not-NMD',coding_exons=[block(0,6,'E16'),block(6,181,'E20')]))
    fig,_=render_compact_pair(pair,[],[],dict(A='M'*479,B='M'*181),models,'',figure_width=16)
    assert fig.get_figwidth()==16
    blocks=[p for p in fig.axes[0].patches if (p.get_gid() or '').startswith('coding_exon_')]
    assert [p.get_width() for p in blocks]==[6,473,6,175]
    import matplotlib.pyplot as plt;plt.close(fig)


def test_partial_alignment_and_conflicting_reported_lengths_stay_explicit():
    ref='MACDEFGHIKLMNPQRSTVWY';alt=ref[:10]
    f,_,events=plan(ref,alt,6,16)
    pair=dict(gene_symbol='G',isoform_A='A',isoform_B='B')
    models=dict(A=dict(nmd='Not-NMD',reported_aa=len(ref)),B=dict(nmd='Not-NMD',reported_aa=8))
    fig,_=render_compact_pair(pair,f,events,dict(A=ref,B=alt),models,'')
    texts=[t.get_text() for t in fig.texts]
    assert 'Partial: sequence alignment only' in texts and '10 aa (table: 8)' in texts
    import matplotlib.pyplot as plt;plt.close(fig)


def test_explicit_feature_reference_keeps_other_reference_records_in_audit():
    seq='MACDEFGHIKLMNPQRSTVWY';rows=rows_for(seq,seq[:10],3,8)
    other=[dict(r,reference_protein='OTHER',feature_id='OTHER:1') for r in rows]
    pair=dict(gene_symbol='G',isoform_A='A',isoform_B='B')
    selected=feature_reference_selection(['G=P'],[pair],rows+other)
    features,audit,_=feature_plan(pair,rows+other,dict(A=seq,B=seq[:10]),reference_proteins=selected['G'])
    assert len(features)==1 and features[0]['reference_protein']=='P'
    assert len(audit)==4 and any(r['feature_id']=='OTHER:1' and r['display_role']=='reference_not_in_view' for r in audit)


def test_feature_reference_selection_rejects_unresolved_or_duplicate_requests():
    pair=dict(gene_symbol='G',isoform_A='A',isoform_B='B');rows=rows_for('MACDEF','MACDEF',2,4)
    for values in [['MISSING=P'],['G=UNKNOWN'],['G=P','G=P'],['G']]:
        with pytest.raises(ValueError):feature_reference_selection(values,[pair],rows)


def test_unplaced_annotation_line_does_not_overlap_native_identifier():
    ref='MACDEFGHIKLMNPQRSTVWY';alt=ref[:10];f,_,events=plan(ref,alt,6,16)
    f[0]['status']='ambiguous_placement'
    pair=dict(gene_symbol='G',isoform_A='A',isoform_B='B')
    models={i:dict(nmd='Not-NMD') for i in ('A','B')}
    fig,_=render_compact_pair(pair,f,events,dict(A=ref,B=alt),models,'')
    ax=fig.axes[0];texts={t.get_gid():t for t in ax.texts};renderer=fig.canvas.get_renderer()
    assert not texts['native_id_B'].get_window_extent(renderer).overlaps(texts['unplaced_B'].get_window_extent(renderer))
    import matplotlib.pyplot as plt;plt.close(fig)


def test_long_key_uses_short_identifier_and_retains_complete_description():
    ref='MACDEFGHIKLMNPQRSTVWY';alt=ref[:10];f,_,events=plan(ref,alt,6,16)
    f[0]['description']='W'*18
    fig,_=render_compact_pair(dict(gene_symbol='G',isoform_A='A',isoform_B='B'),f,events,
                             dict(A=ref,B=alt),{i:dict(nmd='Not-NMD') for i in ('A','B')},'')
    assert any(t.get_text()==f[0]['display_id'] for t in fig.texts)
    assert any('W'*18 in t.get_text() for t in fig.texts)
    import matplotlib.pyplot as plt;plt.close(fig)
