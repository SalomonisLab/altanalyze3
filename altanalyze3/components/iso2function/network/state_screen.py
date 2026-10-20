"""Current-catalog, cross-state isoform/PPI discovery through the iso2function harness.

Expression predominance is defined from every transcript of a TF gene, independently
of assay coverage. PPI differences use query.pair_targets (explicit co-tested calls).
Existing state-vs-other differential statistics are reused without refitting.
"""
from pathlib import Path
import collections,gzip,hashlib,itertools,json,math
import numpy as np
import pandas as pd
from .. import query,config
from . import build


def extract_expression(dataset, groups, clone_genes, out):
    """Audit native full counts/CPM and retain every observed healthy state/library row."""
    import h5py
    out=Path(out);out.mkdir(parents=True,exist_ok=True);dataset=Path(dataset)
    genes=set(clone_genes)
    cp=dataset/'isoform_combined_pseudo_cluster_counts.h5ad'
    pp=dataset/'isoform_combined_pseudo_cluster_cpm.h5ad'
    with h5py.File(cp,'r') as h,h5py.File(pp,'r') as p:
        obs=query._axis(h,'obs');var=query._axis(h,'var')
        assert query._axis(p,'obs')==obs and query._axis(p,'var')==var
        assert h['X'].attrs['encoding-type']=='csr_matrix' and p['X'].attrs['encoding-type']=='csr_matrix'
        features=[];select=[];fg=[]
        for i,v in enumerate(var):
            gene,tid=v.split(':',1)
            if gene in genes:select.append(i);features.append(v);fg.append(gene)
        to_small=np.full(len(var),-1,dtype=np.int32);to_small[select]=np.arange(len(select))
        selected=[];state_rows=collections.defaultdict(list)
        from ...long_read.comparisons import _resolve_uid
        for i,name in enumerate(obs):
            state,sample=name.split('.',1)
            group=groups.get(_resolve_uid(sample,groups,'isoform'),'unknown')
            if group=='healthy':selected.append(i);state_rows[state].append(len(selected)-1)
        if not selected:raise ValueError('No observed healthy libraries matched native metadata')
        counts=np.zeros((len(selected),len(features)));cpm=np.zeros_like(counts)
        row_records=[]
        for j,i in enumerate(selected):
            a,b=map(int,h['X/indptr'][i:i+2]);ix=h['X/indices'][a:b];v=h['X/data'][a:b].astype(float)
            c,d=map(int,p['X/indptr'][i:i+2]);pi=p['X/indices'][c:d];pv=p['X/data'][c:d].astype(float)
            assert len(ix)==len(np.unique(ix)) and np.isfinite(v).all() and (v>=0).all() and np.equal(v,np.rint(v)).all()
            assert np.array_equal(ix,pi) and v.sum()>0 and np.allclose(pv,v/v.sum()*1e6,rtol=1e-5,atol=1e-4)
            small=to_small[ix];mask=small>=0;counts[j,small[mask]]=v[mask];cpm[j,small[mask]]=pv[mask]
            state,sample=obs[i].split('.',1);row_records.append(dict(state=state,sample=sample,source_obs=obs[i],full_profile_counts=float(v.sum())))
            if (j+1)%100==0:print(f'Checked/extracted {j+1}/{len(selected)} healthy profiles',flush=True)
        proof=dict(counts_source=str(cp),cpm_source=str(pp),full_obs_sha256=query._roster_sha(obs),full_var_sha256=query._roster_sha(var),
                   full_features=len(var),full_profiles=len(obs),selected_profiles=len(selected),selected_gene_features=len(features),
                   counts_CPM_full_row_numeric_parity=True,observed_state_library_roster=[obs[i] for i in selected],
                   all_current_TF_gene_features_retained=True,missing_TF_gene_models=sorted(genes-set(fg)))
    np.savez_compressed(out/'EXPRESSION_CACHE.npz',counts=counts,cpm=cpm,features=np.array(features),genes=np.array(fg))
    pd.DataFrame(row_records).to_csv(out/'OBSERVED_HEALTHY_PROFILE_ROSTER.tsv',sep='\t',index=False)
    (out/'EXPRESSION_SOURCE_VERIFICATION.json').write_text(json.dumps(proof,indent=2)+'\n')
    return proof


def screen(source_dir,out):
    """Evaluate both user-authorized regulator rules at all three fold thresholds."""
    source_dir=Path(source_dir);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    contract=json.load(open(source_dir/'SCREEN_CONTRACT.json'))
    clones=json.load(open(source_dir/'ASSAY_CLONE_GENE_PROTEINS.json'));catalog={r['clone']:r for r in clones}
    match=json.load(open(source_dir/'LITERAL_CURRENT_PROTEIN_MATCHES.json'))['matches']
    raw=json.load(gzip.open(source_dir/'DE_STATE_DATA.json.gz','rt'));states=json.load(open(source_dir/'EXPECTED_HEALTHY_STATES.json'))
    assert len(states)==80 and len(raw)==160
    I={};G={};paths=[]
    for s in raw:
        assert '/isoform_final_20261016/analysis/healthy_diff/' in s['source']
        ori=1 if s['difference_column']==f"DiffMeans_{s['state']}_m_Others" else -1
        d={}
        for r in s['rows']:
            def n(x):return float(x) if x else np.nan
            d[r[0]]=dict(mean=n(r[1]),other=n(r[2]),difference=n(r[3])*ori,P=n(r[4]),q=n(r[5]),N=n(r[6]),N_others=n(r[7]))
        (I if s['kind']=='iso' else G)[s['state']]=d
        paths.append({k:s[k] for k in ['state','kind','source','sha256']})
    pd.DataFrame(paths).to_csv(out/'CURRENT_DIFFERENTIAL_SOURCE_MANIFEST.tsv',sep='\t',index=False)
    index=query.symbol_index(config.ENSEMBL91_ANNOT,config.MANE_SUMMARY)
    _,assays=build.build_graph(detected_only=False)
    assays=assays[assays.edge_type.eq('ppi')].copy()
    calls={};conflicts=[]
    for (c,t),a in assays.groupby(['source','target']):
        values={build.assay_call(r.detected) for r in a.itertuples() if build.assay_call(r.tested) is True}
        values.discard(None)
        if len(values)>1:conflicts.append((c,t))
        calls[c,t]=next(iter(values)) if len(values)==1 else None
    if conflicts:raise ValueError('Conflicting literal assay calls: '+repr(conflicts))
    by_gene=collections.defaultdict(list)
    for c in clones:by_gene[c['gene']].append(c['clone'])
    differences={};pair_audit=[]
    for gene,cs in by_gene.items():
        for ca,cb in itertools.combinations(sorted(cs),2):
            rr=query.pair_targets(ca,cb,assays,'ppi');differences[ca,cb]=rr
            for r in rr:pair_audit.append(dict(gene=gene,**r))
    pd.DataFrame(pair_audit).to_csv(out/'ALL_NATIVE_HARNESS_PAIR_TARGET_CALLS.tsv.gz',sep='\t',index=False,compression='gzip')
    z=np.load(source_dir/'EXPRESSION_CACHE.npz');counts=z['counts'];cpm=z['cpm'];features=z['features'];genes=z['genes']
    roster=pd.read_csv(source_dir/'OBSERVED_HEALTHY_PROFILE_ROSTER.tsv',sep='\t');assert len(roster)==len(counts)
    assert set(roster.state)==set(states)
    fi={f:i for i,f in enumerate(features)};by_ensg={g:np.flatnonzero(genes==g) for g in set(genes)}
    state_ix={s:np.flatnonzero(roster.state.to_numpy()==s) for s in states}
    means={s:cpm[ix].mean(axis=0) for s,ix in state_ix.items()}
    aliases=collections.defaultdict(list);links=[]
    for m in match:
        fid=m['ensg']+':'+m['isoform']
        if fid in fi:
            aliases[fid].append(m['clone']);links.append({**m,'feature':fid,'healthy_reads':float(counts[:,fi[fid]].sum())})
    pd.DataFrame(links).to_csv(out/'ALL_LITERAL_PROTEIN_LINKS_AND_READ_SUPPORT.tsv',sep='\t',index=False)
    # No assayed subset is used to decide predominance: unknown-protein RNA models can win.
    dominant=[];dom={}
    for gene,cs in by_gene.items():
        ensg=catalog[cs[0]]['ensg'];cols=by_ensg.get(ensg,[])
        for st in states:
            if not len(cols):dominant.append(dict(gene=gene,ensg=ensg,state=st,status='no current gene models'));continue
            a=means[st][cols];order=np.argsort(-a,kind='stable');top=cols[order[0]];value=float(a[order[0]]);tied=int(np.isclose(a,value,rtol=1e-12,atol=1e-12).sum())
            fid=features[top];ix=state_ix[st];total=float(a.sum());cl=sorted(aliases.get(fid,[]))
            support=int(np.sum(np.isclose(cpm[ix,top],cpm[ix][:,cols].max(axis=1))&(cpm[ix,top]>0)))
            r=dict(gene=gene,ensg=ensg,state=st,feature=fid,isoform=fid.split(':',1)[1],mean_CPM=value,
                   gene_mean_CPM=total,dominant_fraction=value/total if total else np.nan,top_ties=tied,
                   runner_up_CPM=float(a[order[1]]) if len(a)>1 else 0,clones=';'.join(cl),
                   state_profiles=len(ix),profiles_where_top=support,state_reads=float(counts[ix,top].sum()),
                   healthy_reads=float(counts[:,top].sum()),status='clear assayed dominant' if tied==1 and cl else 'unmatched/tied dominant')
            dominant.append(r)
            if tied==1 and cl and value>=contract['min_TF_mean_CPM'] and r['healthy_reads']>=30:dom[gene,st]=r
    pd.DataFrame(dominant).to_csv(out/'ALL_TF_ALL_STATE_PREDOMINANCE.tsv',sep='\t',index=False)
    # Resolve partners from both authoritative annotations, preserving every ambiguous/missing test.
    available_genes=set(next(iter(G.values())))
    mapping={}
    for target in set(assays.target):
        ids=sorted(index.get(target,set()));present=[g for g in ids if g in available_genes]
        mapping[target]=present[0] if len(present)==1 else ''
    pd.DataFrame([dict(target=t,candidate_genes=';'.join(sorted(index.get(t,set()))),resolved_gene=g,
                       status='resolved in current gene DE' if g else 'unresolved or unavailable current gene DE') for t,g in mapping.items()]).to_csv(out/'ALL_PPI_PARTNER_GENE_MAPPING.tsv',sep='\t',index=False)
    def up(r,cut):return r is not None and np.isfinite(r['q']) and r['q']<.05 and r['difference']>cut
    rows=[];edge_rows=[]
    for gene,cs in by_gene.items():
        good=[(s,dom[gene,s]) for s in states if (gene,s) in dom]
        for (sa,ra),(sb,rb) in itertools.combinations(good,2):
            for ca in ra['clones'].split(';'):
                for cb in rb['clones'].split(';'):
                    if ca==cb or catalog[ca]['aa']==catalog[cb]['aa']:continue
                    key=tuple(sorted((ca,cb)));pr=differences[key]
                    set_a={r['target'] for r in pr if r['comparison']==('A_only' if ca==key[0] else 'B_only')}
                    set_b={r['target'] for r in pr if r['comparison']==('B_only' if ca==key[0] else 'A_only')}
                    for fold in contract['fold_thresholds']:
                        cut=math.log2(fold)
                        ia=I[sa].get(ra['feature']);ib=I[sb].get(rb['feature'])
                        ga={t for t in set_a if mapping[t] and up(G[sa].get(mapping[t]),cut)}
                        gb={t for t in set_b if mapping[t] and up(G[sb].get(mapping[t]),cut)}
                        for mode in contract['modes']:
                            pass_iso=mode!='both_upregulated' or (up(ia,cut) and up(ib,cut))
                            ea=ga if pass_iso else set();eb=gb if pass_iso else set()
                            row=dict(gene=gene,state_A=sa,state_B=sb,clone_A=ca,clone_B=cb,isoform_A=ra['isoform'],isoform_B=rb['isoform'],
                                     aa_A=len(catalog[ca]['aa']),aa_B=len(catalog[cb]['aa']),mode=mode,fold=fold,log2_cutoff=cut,
                                     mean_CPM_A=ra['mean_CPM'],mean_CPM_B=rb['mean_CPM'],fraction_A=ra['dominant_fraction'],fraction_B=rb['dominant_fraction'],
                                     donor_top_A=ra['profiles_where_top'],donor_N_A=ra['state_profiles'],donor_top_B=rb['profiles_where_top'],donor_N_B=rb['state_profiles'],
                                     reads_A=ra['state_reads'],reads_B=rb['state_reads'],
                                     isoform_difference_A=np.nan if ia is None else ia['difference'],isoform_q_A=np.nan if ia is None else ia['q'],
                                     isoform_difference_B=np.nan if ib is None else ib['difference'],isoform_q_B=np.nan if ib is None else ib['q'],
                                     co_tested_A_only=len(set_a),co_tested_B_only=len(set_b),up_A_only=len(ea),up_B_only=len(eb),
                                     targets_A='|'.join(sorted(ea)),targets_B='|'.join(sorted(eb)),both_isoform_gates=pass_iso,
                                     bidirectional=bool(ea and eb),any_direction=bool(ea or eb))
                            rows.append(row)
                            if ea or eb:
                                for side,ts,st,clone,fid,other in [('A',ea,sa,ca,ra['isoform'],cb),('B',eb,sb,cb,rb['isoform'],ca)]:
                                    for t in sorted(ts):edge_rows.append(dict(gene=gene,state_A=sa,state_B=sb,clone_A=ca,clone_B=cb,isoform_A=ra['isoform'],isoform_B=rb['isoform'],mode=mode,fold=fold,
                                         side=side,state=st,clone=clone,isoform=fid,target=t,target_gene=mapping[t],call_selected=True,call_other=False,**G[st][mapping[t]]))
    pairs=pd.DataFrame(rows)
    if pairs.empty:raise ValueError('No different-protein dominant isoform state pairs; audits preserved')
    pairs['minimum_unique_partners']=pairs[['up_A_only','up_B_only']].min(axis=1);pairs['total_unique_partners']=pairs.up_A_only+pairs.up_B_only
    pairs['minimum_dominant_fraction']=pairs[['fraction_A','fraction_B']].min(axis=1)
    pairs=pairs.sort_values(['bidirectional','minimum_unique_partners','total_unique_partners','minimum_dominant_fraction'],ascending=False)
    pairs.to_csv(out/'ALL_STATE_PAIR_THRESHOLD_DECISIONS.tsv.gz',sep='\t',index=False,compression='gzip')
    pairs[pairs.any_direction].to_csv(out/'ALL_QUALIFYING_STATE_PAIR_THRESHOLD_CASES.tsv',sep='\t',index=False)
    pd.DataFrame(edge_rows).to_csv(out/'ALL_QUALIFYING_CO_TESTED_UNIQUE_PPI_GENES.tsv.gz',sep='\t',index=False,compression='gzip')
    pairgenes=set(pairs[pairs.any_direction].gene)
    coverage=[]
    for gene,cs in by_gene.items():
        coverage.append(dict(gene=gene,clones=len(cs),clones_with_PPI_tests=len(set(assays[assays.source.isin(cs)].source)),
                             clones_with_literal_matches=len({m['clone'] for m in match if m['gene']==gene}),
                             matched_healthy_transcripts=len({r['isoform'] for r in links if r['gene']==gene}),
                             states_with_clear_assayed_dominant=sum((gene,s) in dom for s in states),qualifying_switch=gene in pairgenes))
    coverage.append(dict(gene='MYC',clones=0,clones_with_PPI_tests=0,clones_with_literal_matches=0,matched_healthy_transcripts=0,states_with_clear_assayed_dominant=0,qualifying_switch=False))
    pd.DataFrame(coverage).to_csv(out/'ALL_TF_SOURCE_AND_RESULT_COVERAGE.tsv',sep='\t',index=False)
    summary=pairs.groupby(['mode','fold']).agg(state_pair_tests=('gene','size'),pairs_any=('any_direction','sum'),pairs_both=('bidirectional','sum')).reset_index()
    summary.to_csv(out/'THRESHOLD_SENSITIVITY.tsv',sep='\t',index=False)
    manifest=dict(status='complete',healthy_states=len(states),assayed_TF_genes=len(by_gene),assayed_clones=len(clones),
                  all_native_observed_profiles_retained=True,all_current_transcripts_used_for_predominance=True,
                  native_pair_call_function='altanalyze3.components.iso2function.query.pair_targets',
                  dominant_definition=contract['predominance'],min_TF_CPM=1,min_native_isoform_reads=30,
                  thresholds=contract['fold_thresholds'],modes=contract['modes'],statistics='original current healthy DE; BH q<0.05; log2(CPM+1) difference > log2(fold)',
                  state_pair_rows=len(pairs),qualifying_rows=int(pairs.any_direction.sum()),bidirectional_rows=int(pairs.bidirectional.sum()))
    (out/'SCREEN_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(summary.to_string(index=False));print(pairs[pairs.any_direction].head(12).to_string(index=False))
    return pairs
