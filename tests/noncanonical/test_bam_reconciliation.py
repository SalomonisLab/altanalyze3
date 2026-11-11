import gzip,json,sqlite3
from pathlib import Path
import numpy as np
import pysam
import anndata as ad
from altanalyze3.components.bam.isoform_structure_extract import extract_isoform_structures
from altanalyze3.components.long_read.noncanonical.ledger import build_ledger
from altanalyze3.components.long_read.noncanonical.manifest import digest


def test_real_bam_native_matrix_and_ledger_conserve_all_queries(tmp_path):
    exon=tmp_path/'exons.tsv'
    exon.write_text('gene\texon\tchr\tstrand\tstart\tstop\nG\tE1.1\tchr1\t+\t100\t149\nG\tE2.1\tchr1\t+\t300\t349\n')
    bam=tmp_path/'reads.bam'
    header={'HD':{'VN':'1.6'},'SQ':[{'SN':'chr1','LN':1000}]}
    with pysam.AlignmentFile(bam,'wb',header=header) as handle:
        for i in range(3):
            read=pysam.AlignedSegment(handle.header);read.query_name=f'molecule/{i+101}'
            read.query_sequence='A'*100;read.reference_id=0;read.reference_start=99
            read.mapping_quality=60;read.cigarstring='50M150N50M';read.set_tag('zm',i+101);read.set_tag('ts','+')
            read.set_tag('CB' if i<2 else 'CR','AAAA-1');handle.write(read)
    pysam.index(str(bam));sources=[]
    # Same query names and identical BAM bytes across declared sources remain separate.
    for sid in ('sourceA','sourceB'):
        folder=tmp_path/'compatibility/samples'/sid;folder.mkdir(parents=True)
        extract_isoform_structures(str(bam),str(exon),str(folder/sid),barcode_tags=['CB'])
        sources.append(dict(source_id=sid,uid='sample',library='same_library',group='control',reverse=False,bam=str(bam),bam_sha256=digest(bam)))
    config=dict(output=str(tmp_path),sources=sources,exon_model=str(exon),modalities=['standard'],_sha256='SYNTHETIC_TEST_ONLY')
    result=build_ledger(config)
    assert result['unique_source_reads']==6
    assert result['actions']=={'retain_standard':4,'retain_unassigned':2}
    assert result['removed_standard_reads']==0
    for sid in ('sourceA','sourceB'):
        before=ad.read_h5ad(tmp_path/'compatibility/samples'/sid/(sid+'.h5ad'))
        after=ad.read_h5ad(tmp_path/'integrated/standard/samples'/sid/(sid+'.h5ad'))
        assert before.obs_names.equals(after.obs_names) and before.var_names.equals(after.var_names)
        np.testing.assert_array_equal(before.X.toarray(),after.X.toarray())
    db=sqlite3.connect(tmp_path/'reconciliation/evidence.sqlite')
    assert db.execute('SELECT COUNT(*) FROM segments').fetchone()[0]==6


def test_two_accepted_junctions_remove_one_native_read_contribution(tmp_path):
    from altanalyze3.components.long_read.noncanonical.manifest import write_tsv
    exon=tmp_path/'exons.tsv';exon.write_text('gene\texon\tchr\tstrand\tstart\tstop\nG\tE1.1\tchr1\t+\t100\t149\nG\tE2.1\tchr1\t+\t300\t349\n')
    bam=tmp_path/'reads.bam';header={'HD':{'VN':'1.6'},'SQ':[{'SN':'chr1','LN':1000},{'SN':'chr2','LN':1000}]}
    with pysam.AlignmentFile(bam,'wb',header=header) as h:
        r=pysam.AlignedSegment(h.header);r.query_name='molecule/1';r.query_sequence='A'*150
        r.reference_id=0;r.reference_start=99;r.mapping_quality=60;r.cigarstring='50M150N50M50S'
        r.set_tag('zm',1);r.set_tag('CB','AAAA-1');r.set_tag('ts','+');r.set_tag('SA','chr2,200,+,100S50M,60,0;');h.write(r)
    pysam.index(str(bam));sid='sourceA';folder=tmp_path/'compatibility/samples'/sid;folder.mkdir(parents=True)
    extract_isoform_structures(str(bam),str(exon),str(folder/sid),barcode_tags=['CB'])
    rows=[];classes=[]
    for d,a in [(349,200),(348,201)]:
        rows.append(dict(query_name='molecule/1',a_chrom='chr1',a_strand='+',donor=str(d),b_chrom='chr2',b_strand='+',acceptor=str(a),gene_a='G',gene_b='H',status='kept',cross_chromosome='True'))
        classes.append(dict(gene_a='G',gene_b='H',joint=f'chr1:{d}:+ > chr2:{a}:+',**{'class':'gene_fusion_candidate'}))
    # Synthetic classifier fixtures exercise accounting, not the biological detector.
    write_tsv(tmp_path/'trans_splicing/samples'/sid/(sid+'.cross_chromosome_junctions.tsv.gz'),list(rows[0]),rows)
    write_tsv(tmp_path/'trans_splicing/classification/joint_classification.tsv.gz',list(classes[0]),classes)
    source=dict(source_id=sid,uid='sample',library='lib',group='control',reverse=False,bam=str(bam),bam_sha256=digest(bam))
    config=dict(output=str(tmp_path),sources=[source],exon_model=str(exon),modalities=['standard','trans-splicing'],_sha256='SYNTHETIC_TEST_ONLY')
    result=build_ledger(config)
    assert result['removed_standard_reads']==1 and result['evidence_rows']==2
    after=ad.read_h5ad(tmp_path/'integrated/standard/samples'/sid/(sid+'.h5ad'))
    assert after.n_vars==1 and after.n_obs==1 and after.X.sum()==0 and after.layers['legacy_counts'].sum()==1
