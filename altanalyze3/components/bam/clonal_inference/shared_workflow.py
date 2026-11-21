"""Nested call-mask selection and full-roster reporting for shared genetic mixtures."""
from pathlib import Path
import itertools,json
import numpy as np
import pandas as pd
from .shared import fit,predict,scores,baseline,expectation
from .integrity import require


def mask_calls(X,events,seed,fraction=.2):
    mask=np.random.default_rng(seed).random(X.shape)<fraction
    groups=events['validation_group'].astype(str).to_numpy() if 'validation_group' in events else events.feature.to_numpy()
    for group in dict.fromkeys(groups):
        cols=np.flatnonzero(groups==group)
        mask[:,cols]=mask[:,cols[0],None]
    return mask&(X>=0)


def merge_profiles(f, observed_features):
    """Marginalize unobserved event states for reporting; retain original raw fits.

    Events are never removed from inference or exports. A globally unobserved
    event is displayed as unknown, so two otherwise identical components cannot
    become different reported clones solely because of an unobserved bit.
    """
    raw=f['G'].copy();raw[:,~observed_features]=-1
    G,inv=np.unique(raw,axis=0,return_inverse=True)
    Q=np.column_stack([f['posterior'][:,inv==i].sum(1) for i in range(len(G))])
    w=np.column_stack([f['w'][:,inv==i].sum(1) for i in range(len(G))])
    require(np.allclose(Q.sum(1),1) and np.allclose(w.sum(1),1),'Merged posterior lost mass')
    return G,Q,w


def label_genotype(g,labels):
    positive=[label for label,state in zip(labels,g) if state==1]
    uncertain=[label for label,state in zip(labels,g) if state==-1]
    main='|'.join(positive) if positive else ('No MUT resolved' if uncertain else 'Panel-WT prediction')
    return main+(' ; unknown: '+','.join(uncertain) if uncertain else '')


def run_shared(X,cells,events,output,args,families):
    require('patient' not in cells or cells.patient.nunique()==1,'Fit unrelated patients separately; shared clones must not cross donors')
    S=list(cells['sample'].unique());s=cells['sample'].map({u:i for i,u in enumerate(S)}).to_numpy()
    J=X.shape[1];seeds=list(range(args.seed,args.seed+args.starts));grid=sorted(set(args.k_grid))
    # K cannot exceed the number of distinct binary genotypes. This is a model
    # cardinality bound, not a reduction in cells/events, and is recorded.
    grid=[k for k in grid if k<=2**J]
    require(len(grid)>0,'No requested K fits the complete binary genotype space')
    require(args.starts>=2,'At least two starts are required for uncertainty reporting')
    observed=(X>=0).any(0);no_evidence=(X<0).all(1)
    outer=mask_calls(X,events,args.validation_seed)
    train=X.copy();train[outer]=-1
    inner=mask_calls(train,events,args.validation_seed+1)
    train_inner=train.copy();train_inner[inner]=-1
    require(inner.any() and outer.any(),'Insufficient observations for nested validation; no silent replacement or feature exclusion')
    np.savez_compressed(output/'validation_masks.npz',outer=outer,inner=inner,cell_ids=cells.cell_id.to_numpy(dtype=str),features=events.feature.to_numpy(dtype=str))
    rows=[]
    def fitted(data,k,seed):
        f=fit(data,s,k,'learned',seed,maxiter=args.max_iterations,fp_cap=args.fp_cap,fn_cap=args.fn_cap)
        require(f['converged'],f'EM did not converge: K={k}, seed={seed}. Preserve this run and increase the explicit budget.')
        require(np.isfinite(f['objective']) and np.allclose(f['posterior'].sum(1),1),'Invalid fitted posterior')
        return f
    for k in grid:
        for seed in seeds:
            f=fitted(train_inner,k,seed)
            rows.append(dict(K=k,seed=seed,objective=f['objective'],converged=f['converged'],**scores(predict(train_inner,s,f),train,inner)))
        pd.DataFrame(rows).to_csv(output/'inner_validation.tsv',sep='\t',index=False)
    means=pd.DataFrame(rows).groupby('K',as_index=False).agg(mean_logloss=('logloss','mean'),mean_brier=('brier','mean'))
    K=int(means.sort_values(['mean_logloss','K']).iloc[0].K);means.to_csv(output/'K_comparison.tsv',sep='\t',index=False)
    selection=dict(K=K,K_grid=grid,requested_K_grid=args.k_grid,selection_at_upper_boundary=K==max(grid) and K<2**J,
                   criterion='Mean inner-mask log loss across starts; exact ties favor smaller K',
                   baseline=scores(baseline(train,s),X,outer),inner_baseline=scores(baseline(train_inner,s),train,inner),
                   seeds=seeds,outer_mask_seed=args.validation_seed,inner_mask_seed=args.validation_seed+1,
                   positive_controls_used=False,known_all_features_retained=True)
    (output/'selection.json').write_text(json.dumps(selection,indent=2))
    full={};outer_scores=[];full_stats=[];maps={};profiles={};confidence={}
    for phase,data in [('outer',train),('full',X)]:
        for seed in seeds:
            f=fitted(data,K,seed);folder=output/f'{phase}_{seed}';folder.mkdir()
            np.savez_compressed(folder/'fit.npz',**{k:f[k] for k in ['G','p','w','errors','posterior']},cell_ids=cells.cell_id.to_numpy(dtype=str),features=events.feature.to_numpy(dtype=str))
            stats={k:f[k] for k in ['loglik','objective','converged','trace']};stats.update(seed=seed,phase=phase,K=K)
            if phase=='outer':
                pred=predict(data,s,f);stats.update(scores(pred,X,outer));outer_scores.append(stats)
                np.savez_compressed(folder/'prediction.npz',call_probability=pred,cell_ids=cells.cell_id.to_numpy(dtype=str),features=events.feature.to_numpy(dtype=str))
            else:
                G,Q,w=merge_profiles(f,observed);ix=Q.argmax(1);maps[seed]=G[ix];maps[seed][no_evidence]=-1
                profiles[seed]=G;confidence[seed]=Q.max(1);full[seed]=(f,G,Q,w,ix);full_stats.append(stats)
            (folder/'stats.json').write_text(json.dumps(stats,indent=2))
    pd.DataFrame([{k:v for k,v in r.items() if k!='trace'} for r in outer_scores]).to_csv(output/'heldout_scores.tsv',sep='\t',index=False)
    chosen=max(full_stats,key=lambda r:r['objective'])['seed'];f,G,Q,w,assignment=full[chosen]
    count=np.array([[np.sum((assignment==k)&(s==j)&~no_evidence) for j in range(len(S))] for k in range(len(G))])
    denominators=np.bincount(s,minlength=len(S));pct=100*count/denominators
    order=sorted(range(len(G)),key=lambda k:(-pct[k].max(),-count[k].sum(),tuple(G[k])))
    ids={k:f'SC{i+1}' for i,k in enumerate(order)};labels=events.event_label.astype(str).tolist() if 'event_label' in events else events.feature.tolist()
    names={k:label_genotype(G[k],labels) for k in order}
    full_agreement=np.mean(np.stack([(maps[seed]==maps[chosen]).all(1) for seed in seeds]),axis=0)
    event_agreement=np.mean(np.stack([maps[seed]==maps[chosen] for seed in seeds]),axis=0)
    ass=cells.copy();ass['clone']=[ids[k] if not ne else 'NE' for k,ne in zip(assignment,no_evidence)]
    ass['composition']=[names[k] if not ne else 'No evidence / unresolved' for k,ne in zip(assignment,no_evidence)]
    ass['assignment_probability']=Q.max(1);ass.loc[no_evidence,'assignment_probability']=np.nan
    ass['observed_events']=(X>=0).sum(1);ass['full_genotype_start_agreement']=full_agreement;ass.loc[no_evidence,'full_genotype_start_agreement']=np.nan
    ass['genotype']=[''.join('?' if v<0 else str(v) for v in g) for g in maps[chosen]]
    ass.to_csv(output/'cell_assignments.tsv.gz',sep='\t',index=False)
    frequency=[];support=[];genotype_rows=[];uncertainty=[]
    # Same genotype keys across all repeat starts; seed ranges are optimization
    # sensitivity, not biological confidence intervals.
    for k in order:
        for j,event in enumerate(events.feature):
            ix=(assignment==k)&~no_evidence
            genotype_rows.append(dict(clone=ids[k],composition=names[k],feature=event,event=labels[j],state=int(G[k,j]),MUT=int((X[ix,j]==1).sum()),WT=int((X[ix,j]==0).sum()),unknown=int((X[ix,j]<0).sum()),event_start_agreement=float(event_agreement[ix,j].mean()) if ix.any() else np.nan,genotype_unobserved_within_assigned_cells=bool(not (X[ix,j]>=0).any())))
        for j,u in enumerate(S):
            ix=(s==j);eq=[int(((maps[seed]==G[k]).all(1)&ix&~no_evidence).sum()) for seed in seeds]
            frequency.append(dict(clone=ids[k],composition=names[k],sample=u,cells=int(count[k,j]),all_cells=int(denominators[j]),percent_all=float(pct[k,j]),posterior_expected_percent_all=float(100*Q[ix&~no_evidence,k].sum()/denominators[j]),start_min_percent_all=100*min(eq)/denominators[j],start_max_percent_all=100*max(eq)/denominators[j]))
    for j,u in enumerate(S):
        n=int((no_evidence&(s==j)).sum());frequency.append(dict(clone='NE',composition='No evidence / unresolved',sample=u,cells=n,all_cells=int(denominators[j]),percent_all=100*n/denominators[j],posterior_expected_percent_all=100*n/denominators[j],start_min_percent_all=100*n/denominators[j],start_max_percent_all=100*n/denominators[j]))
    F=pd.DataFrame(frequency);require(np.allclose(F.groupby('sample').percent_all.sum(),100),'Lost frequency denominator')
    require(np.allclose(F.groupby('sample').posterior_expected_percent_all.sum(),100),'Lost posterior expected denominator')
    F.to_csv(output/'clone_frequencies.tsv',sep='\t',index=False);pd.DataFrame(genotype_rows).to_csv(output/'clone_genotypes.tsv',sep='\t',index=False)
    allfreq=[]
    for seed in seeds:
        patterns=np.array([''.join('?' if v<0 else str(v) for v in g) for g in maps[seed]])
        for j,u in enumerate(S):
            keys,n=np.unique(patterns[s==j],return_counts=True)
            for key,nn in zip(keys,n):allfreq.append(dict(seed=seed,sample=u,genotype=key,cells=int(nn),all_cells=int(denominators[j]),percent_all=100*nn/denominators[j]))
    pd.DataFrame(allfreq).to_csv(output/'all_start_genotype_frequencies.tsv',sep='\t',index=False)
    for a,b in itertools.combinations(seeds,2):
        for event,j in [('COMPLETE',None)]+list(zip(events.feature,range(J))):
            eq=(maps[a]==maps[b]).all(1) if j is None else maps[a][:,j]==maps[b][:,j]
            for u in ['ALL']+S:
                ix=~no_evidence if u=='ALL' else (~no_evidence)&(cells['sample']==u).to_numpy()
                uncertainty.append(dict(seed_a=a,seed_b=b,sample=u,event=event,informative_cells=int(ix.sum()),agreement_percent=float(100*eq[ix].mean()) if ix.any() else np.nan))
    pd.DataFrame(uncertainty).to_csv(output/'start_stability.tsv',sep='\t',index=False)
    errors=[]
    for seed,(fit_,gg,qq,ww,ii) in full.items():
        for j,event in enumerate(events.feature):errors.append(dict(seed=seed,feature=event,event=labels[j],FP=float(fit_['errors'][0,j]),FN=float(1-fit_['errors'][1,j]),FP_at_cap=bool(np.isclose(fit_['errors'][0,j],args.fp_cap)),FN_at_cap=bool(np.isclose(1-fit_['errors'][1,j],args.fn_cap)),observed_calls=int((X[:,j]>=0).sum()),mutant_components=int(fit_['G'][:,j].sum())))
    pd.DataFrame(errors).to_csv(output/'learned_errors.tsv',sep='\t',index=False)
    gp=np.einsum('ik,kj->ij',Q,(G==1).astype(float));gp[:,~observed]=np.nan
    np.savez_compressed(output/'cell_event_posteriors.npz',genotype_probability=gp,event_start_agreement=event_agreement,cell_ids=cells.cell_id.to_numpy(dtype=str),features=events.feature.to_numpy(dtype=str))
    selection.update(representative_seed=chosen,unique_reported_profiles=len(G),no_evidence_cells=int(no_evidence.sum()),all_features=J,globally_unobserved_features=events.feature[~observed].tolist(),interpretation='Operational default, not independently validated genetic truth. Start ranges measure optimization sensitivity, not confidence intervals. All unknown and WT cells retained in frequency denominators.')
    (output/'selection.json').write_text(json.dumps(selection,indent=2))
    family_rows=[]
    for name,definition in families.items():
        cols=[list(events.feature).index(e) for e in definition];take=(G[:,cols]==np.array(list(definition.values()))).all(1)
        for j,u in enumerate(S):family_rows.append(dict(family=name,sample=u,cells=int(count[take,j].sum()),all_cells=int(denominators[j]),percent_all=float(pct[take,j].sum())))
    if family_rows:pd.DataFrame(family_rows).to_csv(output/'family_frequencies.tsv',sep='\t',index=False)
    return selection
