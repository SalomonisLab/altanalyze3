import json
import sqlite3
from pathlib import Path
import pytest
from altanalyze3.components.long_read.noncanonical.ledger import initialize,support,route,event_id
from altanalyze3.components.long_read.noncanonical.manifest import load_config,atomic_json,receipt,verify_receipt

@pytest.fixture
def db():
    db=sqlite3.connect(':memory:');initialize(db)
    for sid in ('a','b'):
        db.execute('INSERT INTO reads VALUES(?,?,?,?,?,?,?,?,?)',(sid,'same_query','CB','zm1','UMI',1,'G:1','native_eligible','{}'))
    yield db
    db.close()

def test_same_query_in_two_bams_and_repeated_crossings_are_distinct(db):
    for sid in ('a','b'):
        for crossing in range(2):support(db,sid,'same_query','circrna',['chr1','+','+',200,100],'kept','candidate',{},crossing)
    assert db.execute('SELECT COUNT(*) FROM support').fetchone()[0]==4
    assert db.execute('SELECT COUNT(*) FROM (SELECT DISTINCT event_id,source_id,query_name FROM support)').fetchone()[0]==2

def test_unknown_source_identity_fails_instead_of_inner_join_loss(db):
    with pytest.raises(ValueError,match='Unknown source query'):support(db,'missing','same_query','circrna',[1],'kept','candidate',{})

def test_coordinate_orientation_is_part_of_event_identity():
    assert event_id('circrna',['chr1','+',200,100])!=event_id('circrna',['chr1','-',100,200])

def test_ambiguous_ptd_circle_preserves_standard_count():
    assert route('G:1',['circrna_candidate','ptd_circrna_ambiguous'],False)[0]=='retain_with_labels'

def test_itd_and_standard_are_compatible_identities():
    assert route('G:1',['itd_tandem_candidate'],False)[0]=='retain_with_labels'

def test_multiple_incompatible_labels_reallocate_once():
    assert route('G:1',['gene_fusion_candidate','trans_splicing_candidate'],True)[0]=='reallocate'

def test_unassigned_read_is_retained():
    assert route('',['gene_fusion_candidate'],True)[0]=='retain_unassigned'

def test_completed_output_must_match_its_hash(tmp_path):
    data=tmp_path/'counts';data.write_text('1')
    config={'_sha256':'pinned'};done=tmp_path/'done.json'
    receipt(done,config,{},[data]);verify_receipt(done,config)
    data.write_text('2')
    with pytest.raises(ValueError,match='changed'):verify_receipt(done,config)

def test_different_manifest_cannot_resume(tmp_path):
    data=tmp_path/'counts';data.write_text('1');done=tmp_path/'done.json'
    receipt(done,{'_sha256':'A'},{},[data])
    with pytest.raises(ValueError,match='another configuration'):verify_receipt(done,{'_sha256':'B'})

def test_duplicate_source_ids_rejected_before_execution(tmp_path):
    cfg=tmp_path/'config.json'
    atomic_json(cfg,{'output':str(tmp_path),'sources':[{'source_id':'s'},{'source_id':'s'}]})
    with pytest.raises(ValueError,match='unique source roster'):load_config(cfg)

def test_itd_cluster_verdict_uses_modal_length_and_cluster_median():
    from altanalyze3.components.long_read.noncanonical.onco_profile import itd_groups
    # A low individual tandem score belongs to a native accepted in-frame cluster.
    rows=[dict(uid='S',read=str(i),kind='insertion',ref_pos='100',length='21',tandem_8mer_fraction=str(v)) for i,v in enumerate([.1,.9,.9])]
    groups=list(itd_groups(rows))
    assert len(groups)==1 and groups[0]['accepted'] and len(groups[0]['rows'])==3
    # Individual in-frame reads must not override a non-in-frame modal cluster.
    rows[0]['length']='22';rows[1]['length']='22'
    assert not list(itd_groups(rows))[0]['accepted']


def test_ratio_preserves_gene_denominator_after_full_feature_reallocation():
    import numpy as np
    from scipy import sparse
    from altanalyze3.components.long_read.noncanonical.quantify import anchored_ratio
    from altanalyze3.components.long_read.isoform_matrix import cpm_and_ratio_from_counts
    baseline=sparse.csr_matrix([[8,2,0],[0,0,5]])
    current=sparse.csr_matrix([[0,2,0],[0,0,5]])
    ratio,mask,genes=anchored_ratio(current,baseline,['G:a','G:b','H:c'])
    np.testing.assert_allclose(ratio.toarray(),[[0,.2,0],[0,0,1]])
    assert mask.toarray().tolist()==[[1,0],[0,1]]
    same=anchored_ratio(baseline,baseline,['G:a','G:b','H:c'])[0]
    native=cpm_and_ratio_from_counts(baseline,['G:a','G:b','H:c'])[1]
    np.testing.assert_array_equal(same.toarray(),native.toarray())


def test_package_change_invalidates_receipt(tmp_path):
    path=tmp_path/'counts';path.write_text('1');done=tmp_path/'done.json'
    receipt(done,{'_sha256':'A','_code_sha256':'old'},{},[path])
    with pytest.raises(ValueError,match='changed package code'):
        verify_receipt(done,{'_sha256':'A','_code_sha256':'new'})

def test_native_ptd_alternative_exon2_is_retained_and_topology_unresolved():
    from altanalyze3.components.long_read.noncanonical.ptd_evidence import mane_exons,classify
    ex=mane_exons();a,b=ex[9]
    row=dict(same_gene_noncolinear='KMT2A',segments=f'chr11:{a}-{b}(+) q0-100 mq60 [KMT2A]; chr11:118439024-{ex[3][1]}(+) q100-199 mq60 [KMT2A]')
    result=classify(row)
    assert result['first_ptd_type_join']['acceptor_exon']=='alt-e2'
    assert result['first_ptd_type_join']['donor_exon']==9
    assert result['native_read_class'].startswith('uninformative')
    assert 'unresolved' in result['topology']


def test_untyped_same_gene_noncolinear_is_not_declared_ptd():
    from altanalyze3.components.long_read.noncanonical.ptd_evidence import mane_exons,classify
    ex=mane_exons();a,b=ex[4];c,d=ex[2]
    row=dict(same_gene_noncolinear='KMT2A',segments=f'chr11:{a}-{b}(+) q0-100 mq60 [KMT2A]; chr11:{c}-{d}(+) q100-200 mq60 [KMT2A]')
    assert classify(row) is None


def test_event_cells_retain_zero_read_annotation_roster_and_orientation():
    import pandas as pd
    from altanalyze3.components.long_read.noncanonical.event_counts import complete_cell_roster
    sources={sid:dict(uid='same_uid',library='same_label',group='g',reverse=reverse)
             for sid,reverse in [('a',False),('b',True)]}
    annotations={'a':pd.Series(['T','B'],index=['AAC-1','CCC-1']),
                 'b':pd.Series(['T','B'],index=['AAC-1','CCC-1'])}
    cells,rows=complete_cell_roster([('a','AAC-1'),('a','NNN-1'),('b','GTT-1'),('b','')],annotations,sources)
    assert set(cells)=={('a','AAC-1'),('a','CCC-1'),('a','NNN-1'),('b','GTT-1'),('b','GGG-1'),('b','')}
    by_id={(r['source_id'],r['barcode']):r for r in rows}
    assert by_id['b','GTT-1']['state']=='T'
    assert by_id['b','GGG-1']['state']=='B'
    assert not by_id['a','CCC-1']['observed_in_bam']
    assert not by_id['b','GGG-1']['observed_in_bam']
    assert by_id['a','NNN-1']['state']=='__UNMAPPED_BARCODE__'
    assert by_id['b','']['state']=='__MISSING_BARCODE__'


def test_duplicate_annotation_rejected_before_native_deduplication(tmp_path):
    from altanalyze3.components.long_read.noncanonical.manifest import verify_barcode_roster
    path=tmp_path/'barcodes.tsv'
    path.write_text('CellBarcode\tcluster\nAAC-1.a\tT\nAAC-1.a\tB\n')
    with pytest.raises(ValueError,match='Duplicate barcode'):verify_barcode_roster(path,['a'])
    path.write_text('AAC-1.a\tT\nAAC-1.b\tB\n')
    assert verify_barcode_roster(path,['a','b'])==2
    with pytest.raises(ValueError,match='Sources missing'):verify_barcode_roster(path,['a','b','c'])


def test_oncofusion_requires_pinned_native_samtools():
    from altanalyze3.components.long_read.noncanonical.manifest import verify_executables
    with pytest.raises(ValueError,match='Pinned native executable required'):
        verify_executables({'modalities':['oncofusion']})
    verify_executables({'modalities':['circrna']})
