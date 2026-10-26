from pathlib import Path
import csv,json,sqlite3
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from altanalyze3.components.annotation.splice_isoform_survey import healthy_counts,scan_events,carriers


def test_healthy_roster_and_state_counts_preserve_every_native_feature(tmp_path):
    import anndata as ad
    x=ad.AnnData(sparse.csr_matrix([[101,0,9],[0,102,8],[7,8,9]]),
        obs=pd.DataFrame(index=['HSC.d1-isoform','Mono.d2-isoform','HSC.d3-isoform']),
        var=pd.DataFrame(index=['G:A.full.id','G:B.full.id','G:C']))
    x.write_h5ad(tmp_path/'counts.h5ad')
    (tmp_path/'meta.tsv').write_text('uid\tgroups\nd1\thealthy\nd2\thealthy\nd3\tdisease\n')
    features,roster,states,donors,pooled,maximum=healthy_counts(tmp_path/'counts.h5ad',tmp_path/'meta.tsv','healthy')
    assert features==list(x.var_names) and len(roster)==3 and sum(r['selected'] for r in roster)==2
    np.testing.assert_array_equal(pooled.toarray(),[[101,0,9],[0,102,8]])
    np.testing.assert_array_equal(maximum,[101,102,9])
    (tmp_path/'meta.tsv').write_text('uid\tgroups\nd1\thealthy\nd2\thealthy\n')
    with pytest.raises(ValueError,match='Unresolved sample'):healthy_counts(tmp_path/'counts.h5ad',tmp_path/'meta.tsv','healthy')


def event_table(tmp_path):
    row={'Gene_Symbol':'GENE','Isoform_1|Length':'A|200','Isoform_2|Length':'B|100',
         'Feature':'G:E1-E3|G:E2-E3','Event_Type':'AltPromoter','Event_Direction':'proximal',
         'Coord_1':'chr1:100-300','Coord_2':'chr1:200-300','N_T':'2','N_Others':'5',
         'Mean_T':'.8','Mean_Others':'.2','DiffMeans_Others_m_T':'-.6','ebayesAdjPval':'.04'}
    p=tmp_path/'healthy-T-Others-T_stats-annotated.txt'
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(row),delimiter='\t');w.writeheader();w.writerow(row)
    return p


def test_original_fdr_and_reversed_direction_are_preserved(tmp_path):
    event_table(tmp_path);events,manifest=scan_events(tmp_path,['T'],.05)
    state=events['G:E1-E3|G:E2-E3']['states']['T']
    assert state['q']==.04 and abs(state['dpsi']-.6)<1e-10 and manifest[0]['negated']
    with pytest.raises(ValueError,match='roster differs'):scan_events(tmp_path,['T','UNMAPPED'],.05)


def test_strict_read_threshold_and_exact_dot_identifiers_use_indexed_native_junctions(tmp_path):
    p=tmp_path/'transcripts.db';c=sqlite3.connect(p)
    c.execute('CREATE TABLE transcripts (gene TEXT,strand TEXT,structure TEXT,tid TEXT,sample TEXT)')
    c.executemany('INSERT INTO transcripts VALUES (?,?,?,?,?)',[
        ('G','+','E1|E3','A.full.id','a'),('G','+','E2|E3','B.full.id','b'),
        ('G','+','E1|E3','C','c'),('G','+','E1|E3|E2|E3','BOTH','d')]);c.commit();c.close()
    features=['G:A.full.id','G:B.full.id','G:C','G:BOTH']
    links,structures,unresolved=carriers(['G:E1-E3|G:E2-E3'],features,np.array([101,102,100,103]),100,p)
    assert not unresolved and links['G:E1-E3']=={0,3} and links['G:E2-E3']=={1,3}
    assert set(structures)=={'G:A.full.id','G:B.full.id','G:BOTH'}
    _,_,unresolved=carriers(['G:E1-E3|G:E2-E3'],features+['G:MISSING'],np.array([101,102,100,103,104]),100,p)
    assert unresolved==[dict(feature='G:MISSING',reason='exact native ID absent from transcript index')]


def test_complete_export_retains_representative_mismatches_and_all_supported_pairs(tmp_path):
    import anndata as ad
    from types import SimpleNamespace
    from altanalyze3.components.annotation.splice_isoform_survey import run
    d=tmp_path/'events';d.mkdir();event_table(d)
    source=(d/'healthy-T-Others-T_stats-annotated.txt').read_text()
    (d/'healthy-U-Others-U_stats-annotated.txt').write_text(source.replace('N_T','N_U').replace('Mean_T','Mean_U').replace('Others_m_T','Others_m_U'))
    a=ad.AnnData(sparse.csr_matrix([[101,1],[2,102]]),obs=pd.DataFrame(index=['T.d1','U.d2']),var=pd.DataFrame(index=['G:A.full.id','G:B.full.id']))
    a.write_h5ad(tmp_path/'counts.h5ad');(tmp_path/'meta.tsv').write_text('uid\tgroups\nd1\thealthy\nd2\thealthy\n')
    c=sqlite3.connect(tmp_path/'transcripts.db');c.execute('CREATE TABLE transcripts (gene TEXT,strand TEXT,structure TEXT,tid TEXT,sample TEXT)')
    c.executemany('INSERT INTO transcripts VALUES (?,?,?,?,?)',[('G','+','E1|E3','A.full.id','a'),('G','+','E2|E3','B.full.id','b')]);c.commit();c.close()
    run(SimpleNamespace(counts=tmp_path/'counts.h5ad',metadata=tmp_path/'meta.tsv',group='healthy',events=d,
        transcript_index=tmp_path/'transcripts.db',min_reads=100,q=.05,outdir=tmp_path/'results'))
    m=json.loads((tmp_path/'results/manifest.json').read_text())
    assert m['status']=='complete' and m['native_isoform_pairs']==1 and m['input_feature_roster']==2
