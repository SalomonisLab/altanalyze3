"""Evidence classes must survive native aggregation, independent re-key and GFF export."""
import csv
import gzip
import json
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import pytest

from altanalyze3.components.long_read.isoform_collapse import pipeline as P
from altanalyze3.components.long_read.isoform_collapse.collapse_sample import collapse_sample
from altanalyze3.components.long_read.isoform_collapse.splice_evidence import Policy, prepare_gene


def fixture(tmp_path):
    rows=[]
    for i in range(3):
        rows.append(dict(library='A',gene='G',molecule=f'p{i}',legacy_structure='I1.1|E2.1|E3.1',
            source_exons=[(150+i,180),(400,450),(600,650)],strand='+'))
    for i in range(4):
        rows.append(dict(library='A',gene='G',molecule=f'r{i}',legacy_structure='I1.1|E2.1|E3.1',
            source_exons=[(150+i,450),(600,650)],strand='+'))
    mapping,candidates=prepare_gene(rows,[('I1.1',101,399)],Policy(3,100,3))
    overrides={('G',r['molecule']):(r['legacy_structure'],mapping['A','G',r['molecule']]) for r in rows}
    ta=tmp_path/'ta.txt'
    ta.write_text(''.join(f"G\t+\t{r['legacy_structure']}\t{r['molecule']}\tA\n" for r in rows))
    ref=tmp_path/'ref.tsv';ref.write_text('gene\texon-id\tchromosome\nG\tE1.1\tchr1\n')
    return rows,ta,ref,candidates,overrides


def test_native_split_and_persisted_rekey_conserve_each_cell(tmp_path):
    rows,ta,ref,candidates,overrides=fixture(tmp_path)
    gr,*_=P.stage1_collapse([('A',str(ta))],nproc=1,ref=ref,tier1_dir=tmp_path/'tier1',
        gene_evidence={'G':dict(candidates=candidates,references={},policy=Policy(3,100,3))},
        molecule_candidates={'A':overrides})
    catalog,kept=P.stage2_outliers(gr,min_total=3)
    assert len(catalog)==2
    assert sorted(r[3] for r in catalog)==[3,4]
    assert len(kept['G'])==2
    a=ad.AnnData(X=sp.csr_matrix(np.array([[1,2,3,4,5,6,7],[7,6,5,4,3,2,1]])),
        obs=pd.DataFrame(index=['cell1','cell2']),var=pd.DataFrame(index=['G:'+r['molecule'] for r in rows]))
    h5=tmp_path/'reads.h5ad';a.write_h5ad(h5)
    P.persist_molecule_candidates(tmp_path,{'A':overrides,'B':{}})
    loaded=P.load_molecule_candidates(tmp_path,'A')
    assert loaded==overrides
    assert P.load_molecule_candidates(tmp_path,'B')=={}
    (tmp_path/'out').mkdir()
    m2s = tmp_path/'viewer_tier1/A.mol2struct.tsv'
    out,*_=P.stage3_rekey_h5ad('A',h5,ta,kept,tmp_path/'out',molecule_candidates=loaded,
                              mol2struct_output=m2s)
    got=ad.read_h5ad(out)
    assert np.array_equal(np.asarray(got.X.sum(axis=1)),np.asarray(a.X.sum(axis=1)))
    for var,values in zip(got.var_names,got.X.toarray().T):
        if var=='G:p0.A':assert values.tolist()==[6,18]
        elif var=='G:r0.A':assert values.tolist()==[22,10]
        else:pytest.fail('Wrong source exemplar selected: '+var)
    # The existing viewer JOIN must see the corrected classes, too.
    from altanalyze3.components.visualization.isv_web import precompute_final_isoform_index as index
    import sqlite3
    final = tmp_path/'FINAL_structure_to_exemplar.tsv'
    final.write_text('gene\tstructure\tfinal_isoform_id\n'+
                     ''.join(f'G\t{s}\t{e}\n' for s,e in kept['G'].items()))
    index.build_mol2final_index({'A':str(m2s)},str(final),str(tmp_path/'viewer'))
    with sqlite3.connect(tmp_path/'viewer/A.mol2final.db') as con:
        resolved = dict(con.execute('SELECT mol,final FROM m2f'))
    assert resolved == {r['molecule']:kept['G'][overrides['G',r['molecule']][1]] for r in rows}


def test_mismatched_and_missing_required_molecules_rejected(tmp_path):
    _,ta,_,_,overrides=fixture(tmp_path)
    bad=dict(overrides);bad['G','p0']=('wrong','new')
    with pytest.raises(ValueError,match='mismatch'):
        collapse_sample('A',ta,tmp_path/'bad',molecule_candidates=bad)
    bad=dict(overrides);bad['G','missing']=('I1.1','new')
    with pytest.raises(ValueError,match='no subset accepted'):
        collapse_sample('A',ta,tmp_path/'missing',molecule_candidates=bad)
    with pytest.raises(ValueError,match='no subset accepted'):
        P.stage3_rekey_h5ad('A',tmp_path/'unused.h5ad',ta,{},tmp_path,molecule_candidates=bad)


def test_tampered_persisted_classes_and_missing_library_rejected(tmp_path):
    _,_,_,_,overrides=fixture(tmp_path)
    P.persist_molecule_candidates(tmp_path,{'A':overrides})
    with pytest.raises(ValueError,match='absent'):
        P.load_molecule_candidates(tmp_path,'other_library')
    manifest=json.loads((tmp_path/'intronic_evidence_manifest.json').read_text())
    path=tmp_path/manifest['samples']['A']['file']
    with path.open('ab') as f:f.write(b'changed')
    with pytest.raises(ValueError,match='Changed'):
        P.load_molecule_candidates(tmp_path,'A')


def test_exemplar_export_uses_exact_class_source_not_legacy_key(tmp_path,monkeypatch):
    rows,ta,ref,candidates,overrides=fixture(tmp_path)
    gr,*_=P.stage1_collapse([('A',str(ta))],nproc=1,ref=ref,tier1_dir=tmp_path/'tier1',
        gene_evidence={'G':dict(candidates=candidates,references={},policy=Policy(3,100,3))},
        molecule_candidates={'A':overrides})
    _,kept=P.stage2_outliers(gr,min_total=3)
    gff=tmp_path/'reads.gff'
    gff.write_text(''.join(f'chr1\tsource\texon\t{s}\t{e}\t.\t+\t.\tgene_id "G"; transcript_id "{r["molecule"]}";\n'
        for r in rows for s,e in r['source_exons']))
    monkeypatch.setattr(P,'translate_combined',lambda *a,**k:None)
    combined,out_ta=P.stage_protein(gr,kept,{'A':str(gff)},tmp_path/'export','unused')
    with gzip.open(combined,'rt') as f:lines=f.readlines()
    assert len(lines)==5
    assert all('p0' in l or 'r0' in l for l in lines)
    assert '~intronic=' not in open(out_ta).read()
    # Fail if a supposedly verified class exemplar changes in its primary file.
    gff.write_text(gff.read_text().replace('\t150\t180\t','\t151\t180\t'))
    with pytest.raises(ValueError,match='geometry differs'):
        P.stage_protein(gr,kept,{'A':str(gff)},tmp_path/'bad_export','unused')


def source_fixture(tmp_path):
    rows,ta,_,_,_=fixture(tmp_path)
    raw=tmp_path/'transcript_associations_raw.txt';raw.write_text(ta.read_text())
    gff=tmp_path/'reads.gff'
    gff.write_text(''.join(f'chr1\tsource\texon\t{s}\t{e}\t.\t+\t.\tgene_id "G"; transcript_id "{r["molecule"]}";\n'
        for r in rows for s,e in r['source_exons']))
    annot=tmp_path/'annotation.tsv'
    annot.write_text('gene\texon-id\tchromosome\tstrand\texon-region-start(s)\texon-region-stop(s)\n'+
        ''.join(f'G\t{k}\tchr1\t+\t{s}\t{e}\n' for k,s,e in
                [('E1.1',50,100),('I1.1',101,399),('E2.1',400,450),('I2.1',451,599),('E3.1',600,650)]))
    cache=tmp_path/'cache.tsv';cache.write_text('gene\tstructure\tENST\nG\tE2.1|E3.1\tENST001\n')
    refgff=tmp_path/'reference.gff3'
    refgff.write_text(''.join(f'chr1\tref\texon\t{s}\t{e}\t.\t+\t.\tgene_id=G.1;transcript_id=ENST001.2;\n'
                         for s,e in [(400,450),(600,650)]))
    return [('A','unused',str(ta),str(gff))],annot,cache,refgff


def test_general_source_preparation_preserves_full_roster(tmp_path):
    from altanalyze3.components.long_read.isoform_collapse.source_evidence import prepare_sources
    samples,annot,cache,refgff=source_fixture(tmp_path)
    evidence,overrides=prepare_sources(samples,annot,cache,refgff,tmp_path/'source_index',Policy(3,100,3))
    assert set(evidence)=={'G'}
    assert len(overrides['A'])==7
    assert sum(c['count'] for c in evidence['G']['candidates'].values())==7
    manifest=json.loads((tmp_path/'source_index/source_manifest.json').read_text())
    assert manifest['all_native_eligible_reads']==7
    assert manifest['complete_reference_identities']==1
    assert manifest['ordered_libraries']==['A']
    assert manifest['upstream_annotation_modified'] is False


def test_raw_source_identity_mismatch_rejects_preparation(tmp_path):
    from altanalyze3.components.long_read.isoform_collapse.source_evidence import prepare_sources
    samples,annot,cache,refgff=source_fixture(tmp_path)
    raw=tmp_path/'transcript_associations_raw.txt';raw.write_text(raw.read_text().replace('\tp0\t','\twrong\t'))
    with pytest.raises(ValueError,match='paired raw/native'):
        prepare_sources(samples,annot,cache,refgff,tmp_path/'source_index',Policy(3,100,3))


def test_missing_primary_reference_rejects_general_preparation(tmp_path):
    from altanalyze3.components.long_read.isoform_collapse.source_evidence import prepare_sources
    samples,annot,cache,refgff=source_fixture(tmp_path)
    cache.write_text(cache.read_text()+'G\tE1.1|E3.1\tENST_missing\n')
    with pytest.raises(ValueError,match='no subset accepted'):
        prepare_sources(samples,annot,cache,refgff,tmp_path/'source_index',Policy(3,100,3))


def test_missing_primary_read_rejects_general_preparation(tmp_path):
    from altanalyze3.components.long_read.isoform_collapse.source_evidence import prepare_sources
    samples,annot,cache,refgff=source_fixture(tmp_path)
    gff=tmp_path/'reads.gff';gff.write_text(''.join(l for l in gff.read_text().splitlines(True) if 'p0' not in l))
    with pytest.raises(ValueError,match='no subset accepted'):
        prepare_sources(samples,annot,cache,refgff,tmp_path/'source_index',Policy(3,100,3))


def test_characterized_retention_reference_with_lost_I_token_is_preserved(tmp_path):
    from altanalyze3.components.long_read.isoform_collapse.source_evidence import prepare_sources
    samples,annot,cache,refgff=source_fixture(tmp_path)
    refgff.write_text(refgff.read_text().replace('\t400\t450\t','\t250\t450\t'))
    evidence,overrides=prepare_sources(samples,annot,cache,refgff,tmp_path/'source_index',Policy(3,100,3))
    ref=evidence['G']['references']['E2.1|E3.1']
    assert ref['annotation_structure']=='I1.1|E2.1|E3.1'
    gr,*_=P.stage1_collapse([(s[0],s[2]) for s in samples],nproc=1,ref=annot,enst_cache=cache,
        tier1_dir=tmp_path/'tier1',gene_evidence=evidence,molecule_candidates=overrides)
    catalog,kept=P.stage2_outliers(gr,min_total=3)
    assert sorted((r[1],r[3],r[5]) for r in catalog)==[('ENST001',4,True),('p0.A',3,False)]
    assert gr['G']['exemplar_annotation_struct']['ENST001']=='I1.1|E2.1|E3.1'


def test_known_evidence_exemplar_GFF3_transcript_identity_matches_catalog(tmp_path,monkeypatch):
    from altanalyze3.components.long_read.isoform_collapse.source_evidence import prepare_sources
    samples,annot,cache,refgff=source_fixture(tmp_path)
    text=refgff.read_text().replace('\t400\t450\t','\t250\t450\t')
    refgff.write_text(text.replace('gene_id=','ID=exon:ENST001.2:1;Parent=ENST001.2;gene_id='))
    original=refgff.read_bytes()
    evidence,overrides=prepare_sources(samples,annot,cache,refgff,tmp_path/'source_index',Policy(3,100,3))
    gr,*_=P.stage1_collapse([(s[0],s[2]) for s in samples],nproc=1,ref=annot,enst_cache=cache,
        tier1_dir=tmp_path/'tier1',gene_evidence=evidence,molecule_candidates=overrides)
    _,kept=P.stage2_outliers(gr,min_total=3)
    monkeypatch.setattr(P,'translate_combined',lambda *a,**k:None)
    combined,ta=P.stage_protein(gr,kept,{'A':samples[0][3]},tmp_path/'export','unused',ref_gff=str(refgff))
    with gzip.open(combined,'rt') as f:known=[l for l in f if 'transcript_id=' in l]
    assert len(known)==2
    assert all('transcript_id=ENST001;' in l and 'Parent=ENST001;' in l for l in known)
    assert refgff.read_bytes()==original


def test_indexed_terminal_events_equal_full_primary_geometry_checks():
    from altanalyze3.components.long_read.isoform_collapse.source_stream import build_event_index,measured_terminal_events
    from altanalyze3.components.long_read.isoform_collapse.splice_evidence import terminal_events
    regions={'I1.1':(101,399),'I2.1':(451,599)}
    parents={'normal':[(50,100),(400,450),(600,650)],'alternative':[(50,180),(550,650)],
             'retention':[(80,450),(600,650)]}
    index=build_event_index(regions,parents)
    introns=[(k,*v) for k,v in regions.items()]
    for start in [50,80,100,101,150,300,399,400,420,450,451,550,599,600]:
        for end in [100,180,399,400,450,599,600,650]:
            if start>end:continue
            exons=[(start,end),(700,750)]
            assert set(measured_terminal_events(exons,index))==set(terminal_events(exons,introns,parents))


def test_streaming_source_scan_keeps_native_counts_and_is_identity_checked(tmp_path):
    from altanalyze3.components.long_read.isoform_collapse.source_stream import build_event_index,scan_sample
    samples,annot,cache,refgff=source_fixture(tmp_path)
    regions={'G':{'E1.1':(50,100),'I1.1':(101,399),'E2.1':(400,450),'I2.1':(451,599),'E3.1':(600,650)}}
    indexes={'G':build_event_index(regions['G'],{'ref':[(400,450),(600,650)]})}
    result=scan_sample(samples[0],regions,indexes,tmp_path/'scan',Policy(3,100,3))
    assert result['complete'] is True
    assert result['native_eligible_reads']==7
    assert result['native_structures']==1
    assert result['possible_evidence_reads']==7
    # The original aligned blocks, rather than the annotation's shared legacy
    # key, distinguish the isolated promoter and continuous intronic coverage.
    import sqlite3
    with sqlite3.connect(result['database']) as c:
        assert c.execute('SELECT SUM(count) FROM models').fetchone()[0]==7
        assert c.execute('SELECT COUNT(*) FROM events').fetchone()[0]==7


def test_stream_contexts_equal_complete_native_tier1_evidence(tmp_path):
    from altanalyze3.components.long_read.isoform_collapse.source_stream import build_event_index,scan_sample
    from altanalyze3.components.long_read.isoform_collapse.stream_context import build_contexts,load_context
    samples,annot,cache,refgff=source_fixture(tmp_path)
    regions={'G':{'E1.1':(50,100),'I1.1':(101,399),'E2.1':(400,450),'I2.1':(451,599),'E3.1':(600,650)}}
    refs={'G':{'E2.1|E3.1':dict(id='ENST001',exons=[(400,450),(600,650)],strand='+')}}
    indexes={'G':build_event_index(regions['G'],{s:r['exons'] for s,r in refs['G'].items()})}
    scan=scan_sample(samples[0],regions,indexes,tmp_path/'scan',Policy(3,100,3))
    result=build_contexts([scan],regions,refs,tmp_path/'contexts',Policy(3,100,3))
    assert result['all_native_reads']==7
    assert result['protected_reads']==7
    context=load_context(result['gene_evidence']['G'])
    assert sorted(c['count'] for c in context['candidates'].values())==[3,4]
    assert sum(c['count'] for c in context['candidates'].values())==7
    from altanalyze3.components.long_read.isoform_collapse import pipeline as pipeline
    _,gr,_=pipeline._collapse_chrom(('test',['G'],result['evidence_sample_paths'],str(cache),'wta',False,{'G':context}))
    catalog,_=pipeline.stage2_outliers(gr,min_total=3)
    assert sorted((r[1],r[3]) for r in catalog)==[('p0.A',3),('r0.A',4)]


def test_prepared_full_tables_and_primary_registry_emit_exact_models(tmp_path,monkeypatch):
    from altanalyze3.components.long_read.isoform_collapse.source_stream import build_event_index,scan_sample
    from altanalyze3.components.long_read.isoform_collapse.stream_context import build_contexts
    samples,annot,cache,refgff=source_fixture(tmp_path)
    regions={'G':{'E1.1':(50,100),'I1.1':(101,399),'E2.1':(400,450),'I2.1':(451,599),'E3.1':(600,650)}}
    refs={'G':{'E2.1|E3.1':dict(id='ENST001',exons=[(400,450),(600,650)],strand='+')}}
    indexes={'G':build_event_index(regions['G'],{s:r['exons'] for s,r in refs['G'].items()})}
    scan=scan_sample(samples[0],regions,indexes,tmp_path/'scan',Policy(3,100,3))
    manifest=build_contexts([scan],regions,refs,tmp_path/'contexts',Policy(3,100,3))
    gr,*_=P.stage1_collapse([(s[0],s[2]) for s in samples],nproc=1,ref=annot,enst_cache=cache,
        tier1_dir=tmp_path/'tier1',gene_evidence=manifest['gene_evidence'],
        prepared_tier1_paths=manifest['evidence_sample_paths'],source_registry={'A':scan['database']})
    _,kept=P.stage2_outliers(gr,min_total=3)
    monkeypatch.setattr(P,'translate_combined',lambda *a,**k:None)
    combined,ta=P.stage_protein(gr,kept,{'A':samples[0][3]},tmp_path/'export','unused',
                              ref_gff=str(refgff),canonical_transcript_ids=True)
    with gzip.open(combined,'rt') as f:records=f.read()
    assert 'transcript_id "p0.A"' in records
    assert 'transcript_id "r0.A"' in records
    assert '~intronic=' not in records+open(ta).read()
    # Baseline aggregation uses every original molecule and original first model.
    bg,*_=P.stage1_collapse([(s[0],s[2]) for s in samples],nproc=1,ref=annot,enst_cache=cache,
        tier1_dir=tmp_path/'baseline_tier1',prepared_tier1_paths=manifest['baseline_sample_paths'],
        source_registry={'A':scan['database']})
    assert sum(bg['G']['exemplar_hardtotal'].values())==7


def test_reference_locus_and_gene_namespace_reconcile_from_original_raw(tmp_path):
    from altanalyze3.components.long_read.isoform_collapse.source_evidence import load_reference_models
    cache=tmp_path/'cache.tsv';cache.write_text('gene\tstructure\tENST\nOLD_X\tE1.1\tENST1\nOLD_Y\tE2.1\tENST1\n')
    raw=tmp_path/'reference_raw.tsv';raw.write_text('OLD_X\t+\tE1.1\tENST1.2\tref\nOLD_Y\t+\tE2.1\tENST1.2_PAR_Y\tref\n')
    gff=tmp_path/'reference.gff3';gff.write_text(
        'chrX\tref\texon\t100\t200\t.\t+\t.\tgene_id=NEW_X.1;transcript_id=ENST1.2;\n'+
        'chrY\tref\texon\t500\t600\t.\t+\t.\tgene_id=NEW_Y.1;transcript_id=ENST1.2_PAR_Y;\n')
    models=load_reference_models(str(cache),str(gff),original_raw=str(raw))
    assert models['OLD_X']['E1.1']['exons']==[(100,200)]
    assert models['OLD_Y']['E2.1']['exons']==[(500,600)]
    assert models['OLD_Y']['E2.1']['primary_gene_id']=='NEW_Y'
    assert models['OLD_Y']['E2.1']['id']=='ENST1'


def test_public_option_runs_complete_streaming_path_and_persists_rekey(tmp_path):
    samples,annot,cache,refgff=source_fixture(tmp_path)
    out=tmp_path/'pipeline'
    result=P.run_pipeline(samples,str(out),nproc=1,min_total=3,ref=str(annot),enst_cache=str(cache),
        ref_gff=str(refgff),write_h5ad=False,tier1_dir=str(tmp_path/'tier1'),protect_intronic_evidence=True)
    assert sorted((r[1],r[3]) for r in result['catalog'])==[('p0.A',3),('r0.A',4)]
    assert len(P.load_molecule_candidates(str(out),'A'))==7
    assert (out/'intronic_source_evidence/pooled/complete_context_manifest.json').exists()


def test_duplicate_library_labels_preserve_both_physical_sources_and_P4_binding(tmp_path):
    import shutil
    samples,annot,cache,refgff=source_fixture(tmp_path)
    other=tmp_path/'second_source';other.mkdir()
    for name in ['ta.txt','transcript_associations_raw.txt','reads.gff']:
        shutil.copyfile(tmp_path/name,other/name)
    sources=samples+[('A','unused',str(other/'ta.txt'),str(other/'reads.gff'))]
    out=tmp_path/'pipeline'
    result=P.run_pipeline(sources,str(out),nproc=1,min_total=3,ref=str(annot),enst_cache=str(cache),
        ref_gff=str(refgff),write_h5ad=False,tier1_dir=str(tmp_path/'tier1'),protect_intronic_evidence=True)
    manifest=json.loads((out/'intronic_evidence_manifest.json').read_text())
    assert len(manifest['samples'])==2
    assert {b['library'] for b in manifest['source_bindings'].values()}=={'A'}
    assert sorted(r[3] for r in result['catalog'])==[6,8]
    bound=[P.bound_sample_name(str(out),'A',s[1],s[2]) for s in sources]
    assert len(set(bound))==2
    assert all(len(P.load_molecule_candidates(out,name))==7 for name in bound)
    with pytest.raises(ValueError,match='source binding unresolved'):
        P.bound_sample_name(str(out),'A','wrong.h5ad',samples[0][2])


def test_known_export_uses_exact_primary_PAR_locus_and_canonical_id(tmp_path, monkeypatch):
    import gzip
    ref = tmp_path/'reference.gff3'
    ref.write_text(
        'chrX\tref\texon\t100\t200\t.\t+\t.\tID=exon:ENST1.2;Parent=ENST1.2;gene_id=NEW_X.1;transcript_id=ENST1.2;\n'+
        'chrY\tref\texon\t500\t600\t.\t+\t.\tID=exon:ENST1.2_PAR_Y;Parent=ENST1.2_PAR_Y;gene_id=NEW_Y.1;transcript_id=ENST1.2_PAR_Y;\n')
    result={'OLD_Y':dict(struct2exemplar={'E2.1':'ENST1'},exemplar_total={'ENST1':3},
        exemplar_known={'ENST1':True},exemplar_sample={'ENST1':'ENSEMBL'})}
    models={'OLD_Y':{'E2.1':dict(id='ENST1',source_transcript_id='ENST1.2_PAR_Y',exons=[(500,600)])}}
    monkeypatch.setattr(P,'translate_combined',lambda *a,**k:None)
    combined,ta=P.stage_protein(result,{'OLD_Y':{'E2.1':'ENST1'}},{},tmp_path/'export','unused',
        ref_gff=str(ref),canonical_transcript_ids=True,reference_models=models)
    with gzip.open(combined,'rt') as f:text=f.read()
    assert text.startswith('chrY\t') and 'chrX' not in text
    assert 'transcript_id=ENST1;' in text and 'Parent=ENST1;' in text
    assert Path(ta).read_text().startswith('OLD_Y\t.\tE2.1\tENST1\t')


def test_known_novel_coding_NMD_summary_retains_every_model(tmp_path):
    import importlib.util
    from pathlib import Path
    path=Path(__file__).parents[1]/'long_read/isoform_collapse/validation/splice_evidence_20261011/summarize_known_novel_protein_classes.py'
    spec=importlib.util.spec_from_file_location('protein_class_validation',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    catalog=tmp_path/'catalog.tsv'
    catalog.write_text('gene\tfinal_isoform_id\texon_blocks\ttotal_reads\tbin\tknown\n'+
        'G\tknown1\t3\t100\tlong\tknown\nG\tnovel1\t2\t10\tother\tnovel\n'+
        'G\tnovel2\t1\t5\tother\tnovel\nG\tnovel3\t1\t4\tother\tnovel\n')
    protein=tmp_path/'protein.tsv'
    protein.write_text('Gene ID\tTranscript ID\tProtein Length\tNMD Status\tIntron Retention\tLongest Isoform Length\n'+
        'G\tknown1\t100\tNot-NMD\tFalse\t100\nG\tnovel1\t90\tPotential-NMD\tTrue\t100\n'+
        'G\tnovel2\t20\tNot-NMD\tFalse\t100\n')
    result=module.summarize(catalog,protein,tmp_path/'summary')
    assert result['all_catalog_models']==4
    counts={(r['known_novel'],r['protein_class']):r['isoforms'] for r in result['counts']}
    assert counts=={('known','protein_coding'):1,('novel','likely_NMD'):1,
        ('novel','truncated_lt_25pct_longest'):1,('novel','unclassified_missing_prediction'):1}


def test_spliced_source_cannot_fold_into_terminal_retention_of_its_gap(tmp_path):
    # Both retain a different C-terminal intron, and their coarse cores match.
    # The more abundant parent's N-terminal block retains 110 nt of a gap
    # explicitly spliced by the child. Checking only the child's C-tail misses it.
    rows=[]
    for prefix,count,exons in [('spliced',3,[(50,100),(310,350),(600,850)]),
                                ('retained',4,[(200,350),(600,850)])]:
        for i in range(count):
            rows.append(dict(library='A',gene='G',molecule=f'{prefix}{i}',strand='+',
                legacy_structure='E1.1|E2.1|I2.1',source_exons=exons))
    mapping,candidates=prepare_gene(rows,[('I2.1',700,999)],Policy(3,100,3),
        parent_models={'observed_spliced':[(50,100),(310,350),(600,699)]})
    overrides={('G',r['molecule']):(r['legacy_structure'],mapping['A','G',r['molecule']]) for r in rows}
    ta=tmp_path/'ta.tsv';ta.write_text(''.join(f"G\t+\t{r['legacy_structure']}\t{r['molecule']}\tA\n" for r in rows))
    ref=tmp_path/'annotation.tsv';ref.write_text('gene\texon-id\tchromosome\nG\tE1.1\tchr1\n')
    cache=tmp_path/'cache.tsv';cache.write_text('gene\tstructure\tENST\nG\tE1.1|E2.1\tENSTspliced\n')
    gr,*_=P.stage1_collapse([('A',str(ta))],nproc=1,ref=ref,tier1_dir=tmp_path/'tier1',
        enst_cache=str(cache),
        gene_evidence={'G':dict(candidates=candidates,references={
            'E1.1|E2.1':dict(id='ENSTspliced',exons=[(50,100),(310,350),(600,699)])},policy=Policy(3,100,3))},
        molecule_candidates={'A':overrides})
    catalog,kept=P.stage2_outliers(gr,min_total=3)
    assert sorted((r[1],r[3]) for r in catalog)==[('retained0.A',4),('spliced0.A',3)]
    assert kept['G'][mapping['A','G','spliced0']]=='spliced0.A'
    assert kept['G'][mapping['A','G','retained0']]=='retained0.A'


def test_unprotected_rare_intronic_splice_keeps_the_existing_retention_path(tmp_path):
    rows=[dict(library='A',gene='G',molecule='rare',strand='+',
        legacy_structure='I1.1_180|E2.1|E3.1',source_exons=[(150,180),(400,450),(600,650)])]
    rows += [dict(library='A',gene='G',molecule=f'ir{i}',strand='+',
        legacy_structure='E1.1|I1.1|E2.1|E3.1',source_exons=[(50,450),(600,650)]) for i in range(4)]
    mapping,candidates=prepare_gene(rows,[('I1.1',101,399)],Policy(3,100,3))
    assert mapping['A','G','rare']==rows[0]['legacy_structure']
    overrides={('G',r['molecule']):(r['legacy_structure'],mapping['A','G',r['molecule']]) for r in rows}
    ta=tmp_path/'ta.tsv';ta.write_text(''.join(f"G\t+\t{r['legacy_structure']}\t{r['molecule']}\tA\n" for r in rows))
    ref=tmp_path/'annotation.tsv';ref.write_text('gene\texon-id\tchromosome\nG\tE1.1\tchr1\n')
    gr,*_=P.stage1_collapse([('A',str(ta))],nproc=1,ref=ref,tier1_dir=tmp_path/'tier1',
        gene_evidence={'G':dict(candidates=candidates,references={},policy=Policy(3,100,3))},
        molecule_candidates={'A':overrides})
    catalog,kept=P.stage2_outliers(gr,min_total=3)
    assert [(r[1],r[3]) for r in catalog]==[('ir0.A',5)]
    assert kept['G'][mapping['A','G','rare']]=='ir0.A'


def test_selected_retention_preserves_original_overlap_coordinate_labels():
    from altanalyze3.components.long_read.isoform_collapse.splice_evidence import recover_selected_terminal_tokens
    path=Path(__file__).parents[1]/'long_read/isoform_collapse/validation/splice_evidence_20261011/overlapping_gene_recovery_fixture.json'
    fixture=json.loads(path.read_text());r=fixture['row'];before=json.dumps(r,sort_keys=True)
    result=recover_selected_terminal_tokens(r,fixture['selected_events'],Policy(3,100,3),fixture['regions'])
    # Measured boundary is preserved, even though the legacy annotator's
    # overlapping-gene prefix disagrees with the region namespace.
    assert 'ENSG00000011052:E5.2_51170013' in result.split('|')
    assert 'I5.1' in result.split('|')
    assert 'E5.1' not in result.split('|')  # truly interpolated outside every primary exon
    assert json.dumps(r,sort_keys=True)==before
    bad=dict(r,raw_structure=r['raw_structure'].replace('E5.2_51170013','E5.2'))
    with pytest.raises(ValueError,match='Required flanking annotation region unresolved: E5.2'):
        recover_selected_terminal_tokens(bad,fixture['selected_events'],Policy(3,100,3),fixture['regions'])


def test_parallel_gene_contexts_equal_serial_complete_outputs(tmp_path):
    from altanalyze3.components.long_read.isoform_collapse.source_stream import build_event_index,scan_sample
    from altanalyze3.components.long_read.isoform_collapse.stream_context import build_contexts,load_context
    from altanalyze3.components.long_read.isoform_collapse.stream_context_parallel import build_parallel_contexts
    import sqlite3
    samples,annot,cache,refgff=source_fixture(tmp_path)
    regions={'G':{'E1.1':(50,100),'I1.1':(101,399),'E2.1':(400,450),'I2.1':(451,599),'E3.1':(600,650)}}
    refs={'G':{'E2.1|E3.1':dict(id='ENST001',exons=[(400,450),(600,650)],strand='+')}}
    indexes={'G':build_event_index(regions['G'],{s:r['exons'] for s,r in refs['G'].items()})}
    scan=scan_sample(samples[0],regions,indexes,tmp_path/'scan',Policy(3,100,3))
    # Two disjoint mock genes, each with the verified fixture's seven records,
    # exercise both worker outputs, repeated headers and full ordered merging.
    with sqlite3.connect(scan['database']) as db:
        for table in ('events','models','ordinary_first'):
            fields=[r[1] for r in db.execute('PRAGMA table_info('+table+')')]
            selection=','.join("'H'" if f=='gene' else f for f in fields)
            db.execute('INSERT INTO '+table+' SELECT '+selection+' FROM '+table+" WHERE gene='G'")
    scan=dict(scan,native_eligible_reads=14)
    regions['H']=dict(regions['G']);refs['H']={k:dict(v,id='ENST002') for k,v in refs['G'].items()}
    serial=build_contexts([scan],regions,refs,tmp_path/'serial',Policy(3,100,3))
    parallel=build_parallel_contexts([scan],regions,refs,tmp_path/'parallel',Policy(3,100,3),workers=2)
    assert parallel['parallel_workers']==2
    assert parallel['complete_partition_gene_identities']==['G','H']
    for field in ('all_native_reads','protected_reads','evidence_genes','source_read_counts','override_counts'):
        assert serial[field]==parallel[field]
    for field in ('evidence_sample_paths','baseline_sample_paths'):
        assert Path(serial[field][0][1]).read_bytes()==Path(parallel[field][0][1]).read_bytes()
    with gzip.open(serial['override_paths']['A'],'rb') as a,gzip.open(parallel['override_paths']['A'],'rb') as b:
        assert a.read()==b.read()
    assert {g:load_context(e) for g,e in serial['gene_evidence'].items()}=={g:load_context(e) for g,e in parallel['gene_evidence'].items()}
    with pytest.raises(ValueError,match='Independent full source total differs'):
        build_parallel_contexts([dict(scan,native_eligible_reads=13)],regions,refs,tmp_path/'bad',Policy(3,100,3),workers=2)


def test_known_low_support_remains_in_full_protein_denominator(tmp_path):
    import importlib.util
    path=Path(__file__).parents[1]/'long_read/isoform_collapse/validation/splice_evidence_20261011/summarize_known_novel_protein_classes.py'
    spec=importlib.util.spec_from_file_location('support_threshold_summary',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    catalog=tmp_path/'catalog.tsv';catalog.write_text('gene\tfinal_isoform_id\ttotal_reads\tknown\nG\tENSTlow\t1\tknown\nG\tnovel\t3\tnovel\n')
    proteins=tmp_path/'protein.tsv';proteins.write_text('Gene ID\tTranscript ID\tProtein Length\tNMD Status\tLongest Isoform Length\nG\tENSTlow\t100\tNot-NMD\t100\nG\tnovel\t90\tPotential-NMD\t100\n')
    result=module.summarize(catalog,proteins,tmp_path/'out')
    assert result['all_catalog_models']==2
    assert result['models_below3_molecules']==1
    assert sum(r['isoforms'] for r in result['counts'])==2
    assert result['counts_ge3_molecules']==[dict(known_novel='novel',protein_class='likely_NMD',isoforms=1)]


def test_recovered_reference_never_overwrites_ordinary_source_annotation(tmp_path,monkeypatch):
    from altanalyze3.components.long_read.isoform_collapse.splice_evidence import recover_selected_terminal_tokens
    region={'E1.1':(50,100),'I1.1':(101,399),'E2.1':(400,450),'I2.1':(451,599),'E3.1':(600,650)}
    rows=[]
    for group,n,legacy,raw,exons in [
        ('ordinary',8,'E2.1|E3.1','E2.1_400|E2.1|E3.1|E3.1_650',[(400,450),(600,650)]),
        ('retention',3,'E2.1|E3.1','I1.1_150|E2.1|E3.1|E3.1_650',[(150,450),(600,650)]),
        ('full',4,'E1.1|E2.1|E3.1','E1.1_50|E1.1|E2.1|E3.1|E3.1_650',[(50,100),(400,450),(600,650)])]:
        for i in range(n):rows.append(dict(library='A',gene='G',strand='+',chromosome='chr1',gff_source='bam',molecule=f'{group}{i}',legacy_structure=legacy,raw_structure=raw,source_exons=exons))
    mapping,candidates=prepare_gene(rows,[(k,*v) for k,v in region.items() if k.startswith('I')],Policy(3,100,3),
        annotation_corrector=lambda r,e,p:recover_selected_terminal_tokens(r,e,p,region))
    before={k:c['annotation_structure'] for k,c in candidates.items()}
    refs={'E2.1|E3.1':dict(id='ENST2000',annotation_structure='I1.1|E2.1|E3.1',exons=[(150,450),(600,650)],strand='+',source_transcript_id='ENST2000.1'),
          'E1.1|E2.1|E3.1':dict(id='ENST1000',exons=[(50,100),(400,450),(600,650)],strand='+',source_transcript_id='ENST1000.1')}
    ta=tmp_path/'ta.txt';ta.write_text(''.join(f"G\t+\t{r['legacy_structure']}\t{r['molecule']}\tA\n" for r in rows))
    cache=tmp_path/'cache.tsv';cache.write_text('gene\tstructure\tENST\n'+''.join(f'G\t{s}\t{r["id"]}\n' for s,r in refs.items()))
    annotation=tmp_path/'annotation.tsv';annotation.write_text('gene\texon-id\tchromosome\nG\tE1.1\tchr1\n')
    overrides={(r['gene'],r['molecule']):(r['legacy_structure'],mapping['A','G',r['molecule']]) for r in rows}
    from altanalyze3.components.long_read.isoform_collapse.stream_context import dump_context,load_context
    roundtripped=load_context(dump_context(tmp_path/'context.json.gz',dict(candidates=candidates,references=refs,policy=Policy(3,100,3))))
    gr,*_=P.stage1_collapse([('A',str(ta))],nproc=1,ref=annotation,enst_cache=cache,tier1_dir=tmp_path/'tier1',
        gene_evidence={'G':roundtripped},molecule_candidates={'A':overrides})
    result=gr['G']
    assert result['struct2exemplar']['E2.1|E3.1']=='ENST1000'
    assert all(result['struct2exemplar'][mapping['A','G',f'retention{i}']]=='ENST2000' for i in range(3))
    assert result['exemplar_total']['ENST1000']==12
    assert result['exemplar_total']['ENST2000']==3
    assert result['exemplar_reference_structure']['ENST2000']=='E2.1|E3.1'
    assert before=={k:c['annotation_structure'] for k,c in candidates.items()}
    # Known GFF export must resolve the aliased representative back to its
    # original cache identity, not pretend the internal key is a source record.
    refgff=tmp_path/'reference.gff3';refgff.write_text(''.join(f'chr1\tref\texon\t{s}\t{e}\t.\t+\t.\tgene_id=G.1;transcript_id={r["source_transcript_id"]};\n' for r in refs.values() for s,e in r['exons']))
    monkeypatch.setattr(P,'translate_combined',lambda *a,**kw:None)
    _,kept=P.stage2_outliers(gr,min_total=3)
    out=tmp_path/'export';out.mkdir()
    path,_=P.stage_protein(gr,kept,{},str(out),'unused',ref_gff=str(refgff),canonical_transcript_ids=True,reference_models={'G':refs})
    with gzip.open(path,'rt') as f:text=f.read()
    assert 'transcript_id=ENST2000;' in text


def test_native_NMD_gene_marker_and_zero_ORF_IDs_reconcile_from_original_TA(tmp_path):
    import importlib.util
    path=Path(__file__).parents[1]/'long_read/isoform_collapse/validation/splice_evidence_20261011/summarize_known_novel_protein_classes.py'
    spec=importlib.util.spec_from_file_location('native_identifier_reconciliation',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    catalog=tmp_path/'catalog.tsv';catalog.write_text('gene\tfinal_isoform_id\ttotal_reads\tknown\nG\tknown\t3\tknown\nH\tzero\t3\tnovel\n')
    protein=tmp_path/'protein.tsv';protein.write_text('Gene ID\tTranscript ID\tProtein Length\tNMD Status\tLongest Isoform Length\nG;(NMD)\tknown\t100\tPotential-NMD\t200\nG;(NMD)\tzero\t0\tNMD\t200\n')
    ta=tmp_path/'ta.tsv';ta.write_text('G\t+\tE1\tknown\tENSEMBL\nH\t+\tE1\tzero\tA\n')
    result=module.summarize(catalog,protein,tmp_path/'out',ta)
    assert result['all_catalog_models']==2
    assert result['predictions_outside_catalog']==0
    assert result['native_prediction_gene_identifiers_reconciled']==2
    assert sum(r['isoforms'] for r in result['counts'] if r['protein_class']=='likely_NMD')==2
    protein.write_text(protein.read_text().replace('zero\t0','zero\t100'))
    with pytest.raises(ValueError,match='Nonzero native prediction gene/source identity unresolved'):
        module.summarize(catalog,protein,tmp_path/'bad',ta)


_REFERENCE_PRIORITY_CASES=json.loads((Path(__file__).parents[1]/'long_read/isoform_collapse/validation/splice_evidence_20261011/source_reference_priority_fixtures.json').read_text())['actual_cases']


def _priority_case_runtime():
    from altanalyze3.components.long_read import isoform_collapse_utils as U
    from altanalyze3.components.long_read.isoform_collapse import scored_collapse as S
    return U,S,dict(structure_tokens_for_containment=U.structure_tokens_for_containment,
        junction_core_tokens=U.junction_core_tokens,is_contiguous_subsequence=U.is_contiguous_subsequence,
        intronic_edge_compatible=U.intronic_edge_compatible,collapse_gene_jc=S.collapse_gene_jc)


def _priority_case_inputs(case):
    import copy
    c=copy.deepcopy(case['candidate'])
    for name in ('protected_junctions','protected_coverage','protected_terminal_blocks'):
        c[name]={tuple(v) for v in c.get(name,())}
    return {case['candidate_key']:c},copy.deepcopy(case['references'])


@pytest.mark.parametrize('case',_REFERENCE_PRIORITY_CASES,ids=[r['gene']+('_junction' if r['candidate']['protected_junctions'] else '_continuation') for r in _REFERENCE_PRIORITY_CASES])
def test_actual_known_compatible_RNA_keeps_native_matching_ranking_and_reference(case):
    from altanalyze3.components.long_read.isoform_collapse.splice_evidence import collapse_prepared
    U,S,runtime=_priority_case_runtime();candidates,refs=_priority_case_inputs(case);gene=case['gene'];key=case['candidate_key']
    native=key.rsplit('~intronic=',1)[0];annotation=candidates[key]['annotation_structure']
    assert U.junction_core_tokens(annotation,gene)!=U.junction_core_tokens(native,gene)
    def checked_native(sr,tok,core,flen,**kw):
        if kw.get('reference_targets_only'):
            # Native replay uses original keys, before evidence-class splitting.
            # It must preserve their counts and the original ranking tokens.
            assert sr[native] == sum(c['count'] for k,c in candidates.items()
                                   if k.rsplit('~intronic=',1)[0] == native)
            assert tok[native] == tuple(U.structure_tokens_for_containment(native,gene))
            assert core[native] == tuple(U.junction_core_tokens(native,gene))
            assert flen[native] == len(tok[native])
            assert 'edge_compatible' not in kw
            return S.collapse_gene_jc(sr,tok,core,flen,**kw)
        assert tok[key]==tuple(U.structure_tokens_for_containment(native,gene))
        assert core[key]==tuple(U.junction_core_tokens(native,gene))
        for s in refs:
            assert core[s]==tuple(U.junction_core_tokens(s,gene))
            assert flen[s]==len(U.structure_tokens_for_containment(s,gene))
        return S.collapse_gene_jc(sr,tok,core,flen,**kw)
    runtime['collapse_gene_jc']=checked_native
    result,final=collapse_prepared(candidates,refs,gene,runtime,splice_slop=3)
    assert result['enst'][final[key]]==case['expected_reference']
    assert candidates[key]['annotation_structure']==annotation
    assert candidates[key]['count']==case['candidate']['count']
    assert candidates[key]['representative_source']==case['candidate']['representative_source']


@pytest.mark.parametrize('case',_REFERENCE_PRIORITY_CASES,ids=[r['gene']+('_junction' if r['candidate']['protected_junctions'] else '_continuation') for r in _REFERENCE_PRIORITY_CASES])
def test_native_matching_view_cannot_bypass_actual_intronic_evidence(case):
    from altanalyze3.components.long_read.isoform_collapse.splice_evidence import collapse_prepared
    U,S,runtime=_priority_case_runtime();candidates,refs=_priority_case_inputs(case);key=case['candidate_key'];c=candidates[key]
    for r in refs.values():
        blocks=sorted(tuple(v) for v in r['exons'])
        if c['protected_junctions']:
            left,right=sorted(c['protected_junctions'])[0]
            index=next(i for i,(s,e) in enumerate(blocks) if abs(s-right)<=3)
            blocks[index]=(blocks[index][0]+20,blocks[index][1])
            assert blocks[index][0]<=blocks[index][1]
        else:
            split=False
            for a,b in sorted(c['protected_coverage']):
                for i,(s,e) in enumerate(blocks):
                    lo,hi=max(a,s),min(b,e)
                    if hi-lo>=40:
                        blocks[i:i+1]=[(s,lo+5),(lo+30,e)];split=True;break
                if split:break
            assert split
        r['exons']=blocks
        assert not U.intronic_edge_compatible(c['protected_junctions'],c['protected_coverage'],blocks,splice_slop=3)
    result,final=collapse_prepared(candidates,refs,case['gene'],runtime,splice_slop=3)
    assert final[key]==key
    assert not result['enst']


def test_recovered_reference_matching_view_requires_known_identity_and_evidence_guard():
    U,S,runtime=_priority_case_runtime()
    with pytest.raises(ValueError,match='known identities and an evidence guard'):
        S.collapse_gene_jc({'c':3},{'c':('E1',),'r':('E2',)}, {'c':('E1',),'r':('E2',)},
            {'c':1,'r':1},enst={'r':'ENST'},reference_evidence_cores={'r':[('E1',)]})


def test_actual_shared_source_reference_key_must_still_pass_retention_guard():
    from altanalyze3.components.long_read.isoform_collapse.splice_evidence import collapse_prepared
    from altanalyze3.components.long_read.isoform_collapse.stream_context import load_context
    import hashlib
    path=Path(__file__).parents[1]/'long_read/isoform_collapse/validation/splice_evidence_20261011/shared_reference_source_geometry_fixture.json.gz'
    ctx=load_context({'file':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    key='E8.1|ENSG00000255245:E3.1|ENSG00000255245:E4.1|ENSG00000255245:E5.1|ENSG00000255245:E6.1|ENSG00000255245:E7.1|E12.1|E12.2'
    assert key in ctx['candidates'] and key in ctx['references']
    assert ctx['candidates'][key]['annotation_structure']==ctx['references'][key].get('annotation_structure',key)
    U,S,runtime=_priority_case_runtime()
    result,final=collapse_prepared(ctx['candidates'],ctx['references'],'ENSG00000137726',runtime,splice_slop=3)
    assert final[key]!=key
    assert result.get('enst',{}).get(final[key])!='ENST00000583660'
    from altanalyze3.components.long_read.isoform_collapse.validation.splice_evidence_20261011.cluster.check_reciprocal_retention import known_junctions,check
    parent=final[key]
    exons=(result['resolved_references'][parent]['exons'] if parent in result['resolved_references'] else ctx['candidates'][parent]['parent_exons'])
    assert not check(ctx['candidates'][key]['parent_exons'],exons,result['parent_retention_coverage'].get(parent,()),known_junctions(ctx['references']),ctx['candidates'][key]['protected_junctions'])


def test_actual_conflicting_retention_molecules_survive_as_compatible_classes():
    from altanalyze3.components.long_read.isoform_collapse_utils import intronic_edge_compatible
    path=Path(__file__).parents[1]/'long_read/isoform_collapse/validation/splice_evidence_20261011/two_global_failure_diagnosis.json'
    data=json.loads(path.read_text());case=next(c for c in data['actual_cases'] if c['failure']['gene']=='ENSG00000099875')
    rows=[dict(library=r['library'],gene=case['failure']['gene'],molecule=r['molecule'],legacy_structure=r['legacy'],source_exons=r['exons'],measured_events=r['measured_events']) for r in case['matching_original_rows']]
    mapping,candidates=prepare_gene(rows,[],Policy(3,100,3),terminal_event_provider=lambda r:r['measured_events'])
    assert len(mapping)==2 and sum(c['count'] for c in candidates.values())==2
    assert len(candidates)==2
    assert {c['representative_source'][2] for c in candidates.values()}=={'54307838','56376487'}
    for c in candidates.values():
        assert intronic_edge_compatible(c['protected_junctions'],c['protected_coverage'],c['parent_exons'],splice_slop=3)
        assert not c['protected_junctions']  # Neither site gained invented >=3 support.


def test_exact_known_source_identity_survives_reference_namespace_separation():
    from altanalyze3.components.long_read.isoform_collapse.splice_evidence import collapse_prepared
    path=Path(__file__).parents[1]/'long_read/isoform_collapse/validation/splice_evidence_20261011/exact_known_source_priority_fixture.json'
    data=json.loads(path.read_text());row=data['control'];gene=row['gene'];key=row['structure']
    U,S,runtime=_priority_case_runtime();refs=data['references']
    tokens={s:tuple(U.structure_tokens_for_containment(s,gene)) for s in refs}
    cores={s:tuple(U.junction_core_tokens(s,gene)) for s in refs}
    native=S.collapse_gene_jc({key:int(row['count'])},tokens,cores,{s:len(t) for s,t in tokens.items()},
                            enst={s:r['id'] for s,r in refs.items()})
    assert native['enst'][native['assignment'].get(key,key)]==data['expected_baseline_parent']
    candidates={key:dict(annotation_structure=key,count=int(row['count']),
        representative_source=(row['library'],gene,row['rep_molecule_id']),
        parent_exons=row['source_exons'],protected_junctions=set(),protected_coverage=set(),
        protected_terminal_blocks=set(),minimum_continuation=100)}
    result,final=collapse_prepared(candidates,refs,gene,runtime,splice_slop=3)
    assert result['enst'][final[key]]==data['expected_baseline_parent']
    parent=result['resolved_references'][final[key]]['exons']
    observed=[(a[1],b[0]) for a,b in zip(row['source_exons'],row['source_exons'][1:])]
    parent_junctions=[(a[1],b[0]) for a,b in zip(parent,parent[1:])]
    assert all(any(max(abs(a-b) for a,b in zip(j,p))<=3 for p in parent_junctions) for j in observed)


@pytest.mark.parametrize('seed',range(20))
def test_without_selected_evidence_namespaced_collapse_replays_native_catalog_and_sources(seed):
    import random
    from collections import Counter
    from altanalyze3.components.long_read.isoform_collapse.splice_evidence import collapse_prepared
    U,S,runtime=_priority_case_runtime();rng=random.Random(seed);gene='NATIVE_REPLAY'
    structures=['|'.join(f'E{i}.1' for i in sorted(rng.sample(range(1,8),rng.randint(2,6))))
                for _ in range(30)]
    structures=list(dict.fromkeys(structures))
    counts={s:rng.randint(1,30) for s in structures[:20]}
    refs={s:dict(id=f'ENST_REPLAY_{i}',exons=[(10,50),(100,150)])
          for i,s in enumerate(structures[::3])}
    tokens={s:tuple(U.structure_tokens_for_containment(s,gene)) for s in set(counts)|set(refs)}
    cores={s:tuple(U.junction_core_tokens(s,gene)) for s in tokens}
    native=S.collapse_gene_jc(counts,tokens,cores,{s:len(t) for s,t in tokens.items()},
                            enst={s:r['id'] for s,r in refs.items()})
    candidates={s:dict(annotation_structure=s,count=n,representative_source=('fixture',gene,str(i)),
        parent_exons=[(10,50),(100,150)],protected_junctions=set(),protected_coverage=set(),
        protected_terminal_blocks=set(),minimum_continuation=100) for i,(s,n) in enumerate(counts.items())}
    result,final=collapse_prepared(candidates,refs,gene,runtime,splice_slop=3)
    canonical=lambda res,s:res['enst'].get(s,s)
    expected={s:canonical(native,native['assignment'].get(s,s)) for s in counts}
    observed={s:canonical(result,final[s]) for s in counts}
    assert observed==expected
    for name in ('long_reps','other_reps'):
        assert {canonical(result,s) for s in result[name]}=={canonical(native,s) for s in native[name]}
    for res in (native,result):
        res['rep_metadata']={canonical(res,s):(res['blocks'][s],res['score'][s])
                             for s in res['long_reps']+res['other_reps']}
    assert result['rep_metadata']==native['rep_metadata']
    assert Counter({s:sum(counts[k] for k,p in observed.items() if p==s) for s in set(observed.values())})==Counter({s:sum(counts[k] for k,p in expected.items() if p==s) for s in set(expected.values())})
