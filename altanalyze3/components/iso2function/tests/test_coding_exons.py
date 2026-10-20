import json
import sqlite3
from pathlib import Path

import pytest

from altanalyze3.components.iso2function.coding_exons import (
    peptide_segments,color_exons,build_coding_exon_map,load_coding_exon_map,assign_exon_labels,annotate_exon_map)
from altanalyze3.components.iso2function.protein_impact import sha,write


def test_split_codon_boundaries_and_terminal_stop_are_exact():
    rows,stop=peptide_segments([(10,14),(30,36)],10,36,'+',3)
    assert stop==3 and rows[0]['aa_end0']==5/3
    assert rows[1]['aa_start0']==5/3 and rows[1]['aa_end0']==3
    assert rows[1]['visible_genomic_end']==33


def test_negative_strand_is_in_transcription_order_with_stop_at_low_end():
    rows,stop=peptide_segments([(10,16),(30,34)],10,34,'-',3)
    assert stop==3 and rows[0]['genomic_start']==30
    assert rows[1]['visible_genomic_start']==13 and rows[1]['visible_genomic_end']==16
    assert rows[1]['aa_end0']==3


def test_unexplained_genomic_length_never_gets_silently_trimmed():
    with pytest.raises(ValueError,match='differs from verified peptide'):
        peptide_segments([(1,20)],1,20,'+',3)


def test_shared_exon_colors_do_not_depend_on_sequence_change_categories():
    a=dict(gene_symbol='G',chrom='1',strand='+',segments=[dict(genomic_start=10,genomic_end=20)])
    b=dict(gene_symbol='G',chrom='1',strand='+',segments=[dict(genomic_start=12,genomic_end=18)])
    color_exons([a,b])
    assert a['segments'][0]['color']==b['segments'][0]['color']
    assert a['segments'][0]['exon_group']==b['segments'][0]['exon_group']==1


def fixture(tmp_path,bad_index=False,chr_prefix=False):
    gff=tmp_path/'combined.gff.gz'
    import gzip
    with gzip.open(gff,'wt') as f:
        for lo,hi in [(10,14),(30,36)]:
            chrom='chr1' if chr_prefix else '1'
            f.write(f'{chrom}\tbam\texon\t{lo}\t{hi}\t.\t+\t.\tgene_id "ENSG1";transcript_id "N.source";\n')
    coding=tmp_path/'coding_regions.txt'
    write(coding,[{'Transcript ID':'N.source','Gene ID':'ENSG1','Chromosome':'1','Strand':'+',
                   'CDS Genomic Min':10,'CDS Genomic Max':36,'Protein Length':3,'NMD Status':'Not-NMD'}])
    native=tmp_path/'native_genomic_audit.tsv'
    write(native,[dict(isoform_id='N.source',native_exon_count=2,genomic_cds_nt=12,
                       coding_regions_aa=3,coding_regions_nmd='Not-NMD',protein_fasta_aa=3,source_table_nmd='Not-NMD')])
    db=tmp_path/'structures.db';c=sqlite3.connect(db)
    c.execute('create table structures(gene,tid,tid_base,chrom,strand,source,exons)')
    c.execute('insert into structures values (?,?,?,?,?,?,?)',
              ('ENSG1','N.source','N.source','chr1' if chr_prefix else '1','+','bam','10-15,31-36' if bad_index else '10-14,30-36'))
    c.commit();c.close()
    impact=dict(inputs={str(p):sha(p) for p in [gff,coding]},outputs={native.name:sha(native)})
    models=[dict(gene_symbol='G',isoform_id='N.source',ensg='ENSG1',nmd='Not-NMD',reported_aa='3',protein_sha256='fixed')]
    return impact,models,db


def test_existing_ISV_index_is_reused_and_attested(tmp_path):
    impact,models,db=fixture(tmp_path)
    rows,source=build_coding_exon_map(tmp_path,impact,models,{'N.source':'MAA'},db)
    assert source['indexed_rows_attested']==1 and source['index_path']==str(db)
    assert rows[0]['segments'][-1]['aa_end0']==3


def test_same_length_wrong_index_coordinates_fail_source_attestation(tmp_path):
    impact,models,db=fixture(tmp_path,bad_index=True)
    with pytest.raises(ValueError,match='differs from pinned native GFF'):
        build_coding_exon_map(tmp_path,impact,models,{'N.source':'MAA'},db)


def test_ISV_chromosome_normalization_reconciles_chr_prefix_without_dropping_model(tmp_path):
    impact,models,db=fixture(tmp_path,chr_prefix=True)
    rows,_=build_coding_exon_map(tmp_path,impact,models,{'N.source':'MAA'},db)
    assert len(rows)==1 and rows[0]['chrom']=='chr1'


def test_cached_map_never_needs_the_original_GFF_or_CDS(tmp_path):
    impact,models,db=fixture(tmp_path)
    rows,_=build_coding_exon_map(tmp_path,impact,models,{'N.source':'MAA'},db)
    (tmp_path/'impact_manifest.json').write_text('{}')
    root=tmp_path/'map';root.mkdir();p=root/'coding_exon_map.json';p.write_text(json.dumps(rows))
    (root/'exon_map_manifest.json').write_text(json.dumps(dict(status='complete',
        source_manifest_sha256=sha(tmp_path/'impact_manifest.json'),outputs={p.name:sha(p)})))
    for name in ['combined.gff.gz','coding_regions.txt','structures.db']:(tmp_path/name).unlink()
    models[0]['actual_aa']='3'
    result,_,_=load_coding_exon_map(str(root),tmp_path,models)
    assert list(result)==['N.source']


def region(name,start,end,strand='-'):
    return dict(gene='ENSG1',exon_id=name,chrom='1',strand=strand,start=start,end=end,type=name[0])


def test_reference_exon_numbers_are_not_color_group_ordinals():
    model=dict(chrom='chr1',strand='-',segments=[dict(visible_genomic_start=100,visible_genomic_end=149,exon_group=1)])
    assign_exon_labels(model,[region('E4.1',125,150),region('E4.2',99,125)])
    b=model['segments'][0]
    assert b['exon_label']=='E4' and b['altanalyze_region_ids']==['E4.1','E4.2']
    assert b['exon_group']==1 and b['exon_label_status']=='assigned'


def test_retained_intron_and_repeated_reference_exon_preserve_their_ids():
    model=dict(chrom='1',strand='-',segments=[dict(visible_genomic_start=90,visible_genomic_end=120),
                                             dict(visible_genomic_start=70,visible_genomic_end=80)])
    assign_exon_labels(model,[region('E3.1',100,120),region('I3.1',81,99),region('E3.2',70,80)])
    assert [b['exon_label'] for b in model['segments']]==['E3/I3','E3']


def test_uncovered_coding_sequence_remains_in_roster_without_invented_number():
    model=dict(chrom='1',strand='+',segments=[dict(visible_genomic_start=10,visible_genomic_end=40)])
    assign_exon_labels(model,[region('E8.1',10,20,'+'),region('E9.1',30,40,'+')])
    assert len(model['segments'])==1
    assert model['segments'][0]['exon_label']==''
    assert model['segments'][0]['exon_label_status']=='unresolved_reference_coverage'


def test_labels_use_existing_gene_model_sqlite_and_export_queried_rows(tmp_path):
    db=tmp_path/'gene_model.db';c=sqlite3.connect(db)
    c.execute('create table gene_model (gene,exon_id,chrom,strand,start,end,const,ens,splice_event,splice_junction)')
    c.execute('insert into gene_model values (?,?,?,?,?,?,?,?,?,?)',('ENSG1','E9.1','1','+',10,20,'','','',''))
    c.commit();c.close()
    models=[dict(ensg='ENSG1',chrom='1',strand='+',segments=[dict(visible_genomic_start=10,visible_genomic_end=20)])]
    source,rows=annotate_exon_map(models,db)
    assert source['reference_rows']==1 and source['unresolved_blocks']==0
    assert models[0]['segments'][0]['exon_label']=='E9'
    assert rows['ENSG1'][0]['exon_id']=='E9.1' and len(source['queried_gene_rows_sha256'])==64
