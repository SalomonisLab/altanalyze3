"""Group accepted multi-event evidence without completing missing genotypes.

SC IDs denote identical positive event sets, not proven distinct biological clones.
Every cell retains its WT/unknown vector. No further genotype imputation occurs.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .integrity import require


def assignments(X):
    require(X.ndim==2 and X.shape[1]>0 and np.isin(X,[-1,0,1]).all(),'Invalid ternary matrix')
    # Full input vector retained separately; 0 here means not positive, never WT.
    signatures=np.array([''.join(map(str,row.astype(int))) for row in X==1],dtype=object)
    positive=(X==1).any(1)
    signatures[~positive&(X==0).all(1)]='WT'
    signatures[~positive&(X<0).any(1)&(X==0).any(1)]='WTobs'
    signatures[(X<0).all(1)]='NE'
    return signatures


def run_supported(X,cells,events,output):
    output=Path(output)
    require(X.shape==(len(cells),len(events)),'Roster shape mismatch')
    require(cells.cell_id.is_unique and events.feature.is_unique,'Duplicate identity')
    signatures=assignments(X);labels=events.event_label.astype(str).tolist() if 'event_label' in events else events.feature.tolist()
    counts=pd.crosstab(pd.Series(signatures,name='signature'),cells['sample'].reset_index(drop=True))
    totals=cells['sample'].value_counts();fractions=counts/totals[counts.columns]
    positives=[x for x in counts.index if x not in ['WT','WTobs','NE']]
    positives.sort(key=lambda x:(-float(fractions.loc[x].max()),-int(counts.loc[x].sum()),x))
    mapping={x:f'SC{i+1}' for i,x in enumerate(positives)};mapping.update({x:x for x in ['WT','WTobs','NE']})
    groups=[mapping[x] for x in signatures];profile_labels={}
    for s in positives:profile_labels[mapping[s]]='|'.join(labels[j] for j,b in enumerate(s) if b=='1')
    profile_labels.update(WT='Panel WT calls at all events',WTobs='WT-only calls; other events unknown',NE='No event calls')
    A=cells.copy();A['clone']=groups;A['positive_signature']=signatures
    A['input_states']=[''.join('?' if v<0 else str(v) for v in row) for row in X]
    A['composition']=[profile_labels[g] for g in groups]
    A['n_MUT']=(X==1).sum(1);A['n_WT']=(X==0).sum(1);A['n_unknown']=(X<0).sum(1)
    A['complete_panel_calls']=A.n_unknown==0
    A.to_csv(output/'cell_assignments.tsv.gz',sep='\t',index=False)
    genotypes=[];frequencies=[]
    order=[mapping[x] for x in positives]+[x for x in ['WT','WTobs','NE'] if x in counts.index]
    for g in order:
        ix=np.flatnonzero(np.array(groups)==g);q=X[ix]
        for j,e in events.iterrows():
            m=int((q[:,j]==1).sum());w=int((q[:,j]==0).sum());u=int((q[:,j]<0).sum())
            # state -1 means mixed WT/unknown or all unknown, never completed WT.
            state=1 if m==len(ix) else (0 if w==len(ix) else -1)
            genotypes.append(dict(clone=g,feature=e.feature,event_label=labels[j],state=state,MUT_cells=m,WT_cells=w,unknown_cells=u,cells=len(ix),MUT_fraction=m/len(ix),WT_fraction=w/len(ix),unknown_fraction=u/len(ix),role=e.get('role','')))
        for sample,n in totals.items():
            z=ix[cells.iloc[ix]['sample'].to_numpy()==sample];frequencies.append(dict(clone=g,composition=profile_labels[g],sample=sample,cells=len(z),all_cells=int(n),percent_all=100*len(z)/n,complete_panel_cells=int((X[z]>=0).all(1).sum()),interpretation='Supported mutation combination; unknown additional events are not WT or proven absent'))
    G=pd.DataFrame(genotypes);F=pd.DataFrame(frequencies)
    F.to_csv(output/'clone_frequencies.tsv',sep='\t',index=False);G.to_csv(output/'clone_genotypes.tsv',sep='\t',index=False)
    # Independent reconstruction verifies every event in every cell, not only totals.
    lookup={g:G[G.clone==g].set_index('feature').loc[events.feature,'state'].to_numpy()==1 for g in order}
    reported=np.stack([lookup[g] for g in groups]);require(np.array_equal(reported,X==1),'Unsupported positive or lost accepted MUT')
    require(np.allclose(F.groupby('sample').percent_all.sum(),100),'Cell mass lost')
    audit=[]
    for sample,ix in cells.groupby('sample',sort=False).groups.items():
        ix=np.asarray(ix)
        for j,e in events.iterrows():
            n=int((X[ix,j]==1).sum());assigned=int(reported[ix,j].sum());require(n==assigned,'Event prevalence changed')
            audit.append(dict(sample=sample,feature=e.feature,event_label=labels[j],all_cells=len(ix),input_MUT=n,reported_MUT=assigned,input_WT=int((X[ix,j]==0).sum()),unknown=int((X[ix,j]<0).sum()),percent_all=100*n/len(ix),unsupported_positive_cells=0,lost_positive_cells=0))
    pd.DataFrame(audit).to_csv(output/'event_frequency_audit.tsv',sep='\t',index=False)
    # Preserve matrix and ordered rosters for plotting/reconstruction.
    np.savez_compressed(output/'accepted_calls.npz',X=X,cell_ids=cells.cell_id.to_numpy(dtype=str),features=events.feature.to_numpy(dtype=str))
    (output/'evidence_gate.json').write_text(json.dumps(dict(cells=len(cells),events=len(events),samples=cells['sample'].nunique(),positive_groups=len(positives),unsupported_positive_cells=0,lost_positive_cells=0,unknown_to_WT=0,additional_imputation=False,all_cell_denominators=True,interpretation='Evidence-supported combinations, not fully resolved genetic clones'),indent=2))
