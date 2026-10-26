"""A switch is a change in whole-gene leader, not a change between two minor isoforms."""
import pandas as pd
import pytest
from altanalyze3.components.long_read.isoform_state_switches import predominant_pair_context


def tables():
    # Independent, hand-specified usage distributions. gene minor's two minor
    # forms reverse while its85% leader is unchanged; close reverses by2points.
    values = {'switch': ([.6,.3,.1],[.2,.7,.1]),
              'minor': ([.05,.1,.85],[.1,.05,.85]),
              'close': ([.51,.49,0],[.49,.51,0]),
              'tie': ([.5,.5,0],[.2,.7,.1]),
              'one_context': ([.8,.2,0],None)}
    means=[];pred=[]
    for gene,contexts in values.items():
        for group,usages in zip(['A','B'],contexts):
            if usages is None:continue
            for i,value in enumerate(usages):
                if value:means.append(dict(state='S',group=group,gene=gene,feature=gene+':'+str(i),mean_usage=value,donors=3,assigned_reads=value*1000))
            i=usages.index(max(usages));pred.append(dict(state='S',group=group,gene=gene,predominant=gene+':'+str(i),mean_usage=usages[i],donors=3))
    return pd.DataFrame(pred),pd.DataFrame(means)


def test_whole_gene_maximum_and_twenty_point_margin_are_both_required():
    p,u=tables();before=u.copy(deep=True)
    result=predominant_pair_context(p,u,'S','A','S','B').set_index('gene')
    assert u.equals(before) and set(result.index)==set(p.gene)
    assert set(result.index[result.predominance_switch_at_margin])=={'switch'}
    assert result.loc['switch','margin_A']==pytest.approx(.3)
    assert result.loc['switch','margin_B']==pytest.approx(.5)
    assert result.loc['minor','predominant_A']==result.loc['minor','predominant_B']=='minor:2'
    assert result.loc['close','leader_changed'] and not result.loc['close','predominance_switch_at_margin']
    assert result.loc['tie','top_ties_A']==2 and not result.loc['tie','predominance_switch_at_margin']
    assert not result.loc['one_context','both_contexts_observed']
    assert pd.isna(result.loc['one_context','usage_B_in_B'])


def test_context_partitions_are_lossless_and_identical_to_native_queries():
    from altanalyze3.components.long_read.isoform_state_switches import predominant_context_partitions
    p,u=tables();partitions=predominant_context_partitions(p,u)
    for source,parts in zip([p,u],partitions):
        reconstructed=pd.concat(list(parts.values())).sort_index()
        pd.testing.assert_frame_equal(reconstructed,source)
    for a,b in [('A','B'),('B','A'),('A','UNAVAILABLE')]:
        baseline=predominant_pair_context(p,u,'S',a,'S',b)
        cached=predominant_pair_context(p,u,'S',a,'S',b,context_partitions=partitions)
        pd.testing.assert_frame_equal(baseline,cached)


def test_margin_can_pass_in_only_one_context_and_zero_is_not_missing():
    p=pd.DataFrame([dict(state='S1',group='H',gene='G',predominant='G:a',mean_usage=.9,donors=3),
                    dict(state='S2',group='H',gene='G',predominant='G:b',mean_usage=.51,donors=3)])
    u=pd.DataFrame([dict(state='S1',group='H',gene='G',feature='G:a',mean_usage=.9,donors=3,assigned_reads=90),
                    dict(state='S1',group='H',gene='G',feature='G:c',mean_usage=.1,donors=3,assigned_reads=10),
                    dict(state='S2',group='H',gene='G',feature='G:a',mean_usage=.49,donors=3,assigned_reads=49),
                    dict(state='S2',group='H',gene='G',feature='G:b',mean_usage=.51,donors=3,assigned_reads=51)])
    r=predominant_pair_context(p,u,'S1','H','S2','H').iloc[0]
    assert r.predominance_switch_at_margin and r.margin_A==.9
    assert r.margin_B==pytest.approx(.02) and r.usage_B_in_A==0 and r.reads_B_in_A==0


def test_native_context_means_use_all_gene_features_and_equal_donor_weights(tmp_path):
    import anndata as ad
    import numpy as np
    from scipy import sparse
    from altanalyze3.components.long_read.isoform_state_switches import predominant_isoforms
    a=ad.AnnData(sparse.csr_matrix([[9,1,1000000],[10,90,1000000],[1,9,1000000]]),
                 obs=pd.DataFrame(index=['S.D1-isoform','S.D2-isoform','S.D3-isoform']),
                 var=pd.DataFrame(index=['G:a','G:b','BACKGROUND:x']))
    counts=tmp_path/'counts.h5ad';a.write_h5ad(counts)
    meta=tmp_path/'groups.tsv';pd.DataFrame({'uid':['D1','D2','D3'],'groups':['H']*3}).to_csv(meta,sep='\t',index=False)
    default=predominant_isoforms(counts,meta,['H'],log=lambda *args:None)
    pred,u=predominant_isoforms(counts,meta,['H'],log=lambda *args:None,return_usage=True)
    pd.testing.assert_frame_equal(pred,default)
    g=u[u.gene.eq('G')].set_index('feature')
    assert g.loc['G:a','mean_usage']==pytest.approx((.9+.1+.1)/3)
    assert g.loc['G:b','mean_usage']==pytest.approx((.1+.9+.9)/3)
    assert g.loc['G:a','assigned_reads']==20 and g.loc['G:b','assigned_reads']==100
    assert set(u.feature)=={'G:a','G:b','BACKGROUND:x'} and set(u.donors)=={3}


def test_native_cli_keeps_unchanged_genes_and_coding_denominators(tmp_path):
    import json
    import anndata as ad
    from scipy import sparse
    from altanalyze3.components.long_read.isoform_state_switches import main
    obs=[s+'.'+d+'-isoform' for s in ['S1','S2'] for d in ['D1','D2','D3']]
    a=ad.AnnData(sparse.csr_matrix([[80,20,1000000]]*3+[[20,80,1000000]]*3),
                 obs=pd.DataFrame(index=obs),var=pd.DataFrame(index=['G:a','G:b','BG:x']))
    counts=tmp_path/'counts.h5ad';a.write_h5ad(counts)
    meta=tmp_path/'meta.tsv';pd.DataFrame({'uid':['D1','D2','D3'],'groups':['H']*3}).to_csv(meta,sep='\t',index=False)
    classes=tmp_path/'classes.tsv';pd.DataFrame([['G','a','coding',100,'Not-NMD'],['G','b','coding',80,'Not-NMD'],['BG','x','none',0,'']],
            columns=['gene','isoform','class','protein_length','nmd_status']).to_csv(classes,sep='\t',index=False)
    out=tmp_path/'out';main(['predominance','--counts',str(counts),'--groups',str(meta),'--group','H',
                           '--protein_classes',str(classes),'--outdir',str(out)])
    p=json.loads((out/'predominance_parameters.json').read_text());assert not p['statistical_test_performed'] and not p['differential_whitelist_used']
    row=p['contexts'][0];assert row['all_gene_roles']==2 and row['switches_at_margin']==row['well_covered_switches']==1
    r=pd.read_csv(row['output'],sep='\t').set_index('gene');assert set(r.index)=={'G','BG'}
    assert r.loc['G','usage_A_in_A']==pytest.approx(.8) and r.loc['G','usage_B_in_A']==pytest.approx(.2)
    c=pd.read_csv(row['predicted_coding_output'],sep='\t');assert len(c)==1 and c.iloc[0].well_covered_predominance_switch
    assert p['inputs']['protein_classes']['sha256']


def test_reads_below_native_cpm_conditioning_are_not_fabricated_as_zero(tmp_path):
    import anndata as ad
    from scipy import sparse
    from altanalyze3.components.long_read.isoform_state_switches import predominant_isoforms
    # In S2 the b-form's one read is confined to a donor where G is below CPM1.
    # G is measured in the other S2 donors, so b's conditional mean usage is0,
    # while its actual assigned-read support remains1.
    names=[s+'.'+d+'-isoform' for s in ['S1','S2'] for d in ['D1','D2','D3']]
    a=ad.AnnData(sparse.csr_matrix([[0,100,1000000]]*3+[[100,0,1000000],[0,1,1000000000],[100,0,1000000]]),
                 obs=pd.DataFrame(index=names),var=pd.DataFrame(index=['G:a','G:b','BG:x']))
    counts=tmp_path/'counts.h5ad';a.write_h5ad(counts)
    meta=tmp_path/'meta.tsv';pd.DataFrame({'uid':['D1','D2','D3'],'groups':['H']*3}).to_csv(meta,sep='\t',index=False)
    p,u=predominant_isoforms(counts,meta,['H'],log=lambda *args:None,return_usage=True)
    row=u[(u.state=='S2') & (u.feature=='G:b')].iloc[0]
    assert row.mean_usage==0 and row.assigned_reads==1 and row.donors==2
    r=predominant_pair_context(p,u,'S1','H','S2','H').set_index('gene').loc['G']
    assert r.predominant_A=='G:b' and r.usage_A_in_B==0 and r.reads_A_in_B==1
    from altanalyze3.components.long_read.isoform_state_switches import main
    classes=tmp_path/'classes.tsv'
    pd.DataFrame([['G','a','coding',100,'Not-NMD'],['G','b','coding',80,'Not-NMD'],['BG','x','none',0,'']],
                 columns=['gene','isoform','class','protein_length','nmd_status']).to_csv(classes,sep='\t',index=False)
    out=tmp_path/'coding_out'
    main(['predominance','--counts',str(counts),'--groups',str(meta),'--group','H',
          '--protein_classes',str(classes),'--outdir',str(out)])
    coding=pd.read_csv(out/'context_pairs/predicted_coding_all_context_pairs.tsv.gz',sep='\t').set_index('gene')
    assert coding.loc['G','predominant_A']=='G:b'
    assert coding.loc['G','usage_A_in_B']==0 and coding.loc['G','reads_A_in_B']==1


def test_explicit_pool_sums_counts_within_donors_and_preserves_original_contexts(tmp_path):
    import anndata as ad
    import numpy as np
    from scipy import sparse
    from altanalyze3.components.long_read.isoform_state_switches import predominant_isoforms
    a=ad.AnnData(sparse.csr_matrix(np.array([[9,1,1000000],[0,90,1000000],[80,20,1000000]],dtype=np.float32)),
                 obs=pd.DataFrame(index=['S1.D1-isoform','S2.D1-isoform','S1.D2-isoform']),
                 var=pd.DataFrame(index=['G:a','G:b','BG:x']))
    counts=tmp_path/'counts.h5ad';a.write_h5ad(counts)
    meta=tmp_path/'meta.tsv';pd.DataFrame({'uid':['D1','D2'],'groups':['H','H']}).to_csv(meta,sep='\t',index=False)
    original,u0=predominant_isoforms(counts,meta,['H'],return_usage=True,log=lambda *a:None)
    pred,u=predominant_isoforms(counts,meta,['H'],return_usage=True,state_pools={'Pool':['S1','S2']},log=lambda *a:None)
    pd.testing.assert_frame_equal(pred[pred.state!='Pool'].reset_index(drop=True),original,check_exact=True)
    pd.testing.assert_frame_equal(u[u.state!='Pool'].reset_index(drop=True),u0,check_exact=True)
    g=u[(u.state=='Pool') & (u.feature=='G:a')].iloc[0]
    assert g.donors==2 and g.assigned_reads==89
    assert g.mean_usage==pytest.approx((.09+.8)/2)
    with pytest.raises(ValueError,match='unresolved state pool'):
        predominant_isoforms(counts,meta,['H'],state_pools={'Pool':['S1','UNRESOLVED']},log=lambda *a:None)
