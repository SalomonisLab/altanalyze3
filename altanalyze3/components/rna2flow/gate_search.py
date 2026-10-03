"""Fast, auditable gates against annotation membership, with untouched test events.

These metrics validate approximation of labels, not stem-cell function or transfer truth.
Marker selection and thresholds use training events only; validation selects the strategy.
An exhaustive histogram rectangle search handles both positive and exclusion markers.
Short decision-tree paths provide sequential gates when two markers are insufficient.
"""
from itertools import combinations
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from .gating import rank_channels_for_population


def membership(X, features, conditions):
    mask = np.ones(len(X), bool)
    for c in conditions:
        if c['operator'] not in ('>', '<=') or not np.isfinite(c['threshold']):
            raise ValueError('Gate operators must be > or <= with finite thresholds')
        v = X[:, features.index(c['channel'])]
        mask &= v > c['threshold'] if c['operator'] == '>' else v <= c['threshold']
    return mask


def metrics(mask, target, beta=1):
    tp = int(np.sum(mask & target)); ng = int(mask.sum()); nt = int(target.sum())
    p = tp / ng if ng else 0.; r = tp / nt if nt else 0.
    return dict(n_events=len(target), n_target=nt, n_gated=ng, true_positive=tp,
                precision=p, recall=r, f1=2*p*r/(p+r) if p+r else 0.,
                fbeta=(1+beta**2)*p*r/(beta**2*p+r) if p+r else 0.)


def _rectangle(x, y, target, beta, bins=24):
    # Strict upper bounds and searchsorted(side='left') agree at tied detector values.
    ex = np.unique(np.quantile(x, np.linspace(0, 1, bins+1))[1:-1])
    ey = np.unique(np.quantile(y, np.linspace(0, 1, bins+1))[1:-1])
    ix=np.searchsorted(ex, x, side='left'); iy=np.searchsorted(ey, y, side='left')
    nx,ny=len(ex)+1,len(ey)+1
    total=np.zeros((nx,ny)); positive=np.zeros((nx,ny))
    np.add.at(total,(ix,iy),1);np.add.at(positive,(ix[target],iy[target]),1)
    def prefix(h):return np.pad(h.cumsum(0).cumsum(1),((1,0),(1,0)))
    ht,hp=prefix(total),prefix(positive)
    yl,yh=np.triu_indices(ny+1,1)
    best=(-1,None)
    for xl in range(nx):
        for xh in range(xl+1,nx+1):
            ng=ht[xh,yh]-ht[xl,yh]-ht[xh,yl]+ht[xl,yl]
            tp=hp[xh,yh]-hp[xl,yh]-hp[xh,yl]+hp[xl,yl]
            score=(1+beta**2)*tp / np.maximum(beta**2*target.sum()+ng,1)
            k=int(score.argmax())
            if score[k]>best[0]:best=(score[k],(xl,xh,int(yl[k]),int(yh[k])))
    xl,xh,yl,yh=best[1]
    return [None if xl==0 else float(ex[xl-1]),None if xh==nx else float(ex[xh-1]),
            None if yl==0 else float(ey[yl-1]),None if yh==ny else float(ey[yh-1])]


def optimize_population(X, labels, features, populations, *, beta=1., seed=42,
                        max_markers=6, bins=24, groups=None):
    X=np.asarray(X);labels=np.asarray(labels).astype(str);features=list(features)
    if X.ndim!=2 or X.shape!=(len(labels),len(features)) or not np.isfinite(X).all():
        raise ValueError('Expected finite events × named markers and one label per event')
    if not .25<=beta<=4:raise ValueError('beta must be between 0.25 and 4')
    populations=list(populations)
    missing=set(populations)-set(labels)
    if missing:raise ValueError('No events assigned to: '+', '.join(sorted(missing)))
    target=np.isin(labels,populations)
    if min(target.sum(),(~target).sum())<30:
        raise ValueError('Need at least 30 target and 30 background events for held-out gating')
    indices=np.arange(len(X))
    if groups is None:
        train,rest=train_test_split(indices,test_size=.4,stratify=target,random_state=seed)
        valid,test=train_test_split(rest,test_size=.5,stratify=target[rest],random_state=seed)
        unit='stratified events (60/20/20); biological replicate identity unavailable'
    else:
        from sklearn.model_selection import GroupShuffleSplit
        groups=np.asarray(groups)
        train,rest=next(GroupShuffleSplit(n_splits=1,test_size=.4,random_state=seed).split(X,target,groups))
        vi,ti=next(GroupShuffleSplit(n_splits=1,test_size=.5,random_state=seed).split(X[rest],target[rest],groups[rest]))
        valid,test=rest[vi],rest[ti];unit='held-out groups (approximately 60/20/20)'
    if any(min(target[z].sum(),(~target[z]).sum())<5 for z in [train,valid,test]):
        raise ValueError('Target/background insufficient in one split; change groups or target')
    # Binary target uses all background, rather than dropping small competing populations.
    binary=np.where(target,'target','background')
    rho=rank_channels_for_population(X[train],binary[train],features,'target',min_cells=5)
    ranked=rho.fillna(0).abs().sort_values(ascending=False).index.tolist()
    chosen=[c for c in ranked if np.ptp(X[train,features.index(c)])>0][:max_markers]
    if len(chosen)<2:raise ValueError('Need two informative channels')
    candidates=[]
    for a,b in combinations(chosen,2):
        box=_rectangle(X[train,features.index(a)],X[train,features.index(b)],target[train],beta,bins)
        conditions=[]
        for ch,lo,hi in [(a,box[0],box[1]),(b,box[2],box[3])]:
            if lo is not None:conditions.append(dict(channel=ch,operator='>',threshold=lo))
            if hi is not None:conditions.append(dict(channel=ch,operator='<=',threshold=hi))
        candidates.append(dict(method='rectangle',channels=[a,b],bounds=box,conditions=conditions))
    columns=[features.index(c) for c in chosen]
    for depth in [2,3,4]:
        tree=DecisionTreeClassifier(max_depth=depth,min_samples_leaf=max(10,len(train)//500),
                                    random_state=seed).fit(X[train][:,columns],target[train])
        leaves=[]
        def visit(node,path):
            tr=tree.tree_
            if tr.children_left[node]==tr.children_right[node]:leaves.append(path);return
            ch=chosen[tr.feature[node]];th=float(tr.threshold[node])
            visit(tr.children_left[node],path+[dict(channel=ch,operator='<=',threshold=th)])
            visit(tr.children_right[node],path+[dict(channel=ch,operator='>',threshold=th)])
        visit(0,[])
        # Select one path on TRAIN only; validation then compares depths and rectangles.
        path=max(leaves,key=lambda cs:metrics(membership(X[train],features,cs),target[train],beta)['fbeta'])
        candidates.append(dict(method='sequential_tree',max_depth=depth,
                               channels=list(dict.fromkeys(c['channel'] for c in path)),conditions=path))
    for c in candidates:
        c['train']=metrics(membership(X[train],features,c['conditions']),target[train],beta)
        c['validation']=metrics(membership(X[valid],features,c['conditions']),target[valid],beta)
    candidates.sort(key=lambda c:(-c['validation']['fbeta'],len(c['conditions'])))
    winner=dict(candidates[0]);mask=membership(X,features,winner['conditions'])
    winner['test']=metrics(mask[test],target[test],beta)
    winner['all_events']=metrics(mask,target,beta)
    null=np.random.default_rng(seed).permutation(target[test])
    winner['shuffled_test_labels']=metrics(mask[test],null,beta)
    return dict(populations=populations,objective='F-beta',beta=beta,seed=seed,
                split_unit=unit,n_candidates=len(candidates),best=winner,candidates=candidates,
                marker_ranking=[dict(channel=c,rho=float(rho[c])) for c in ranked],
                limitation='Capture and purity are against annotation membership. The test events are withheld from gate fitting and selection, but labels may have been inferred from these same antibodies. No functional stem-cell validation.',
                test_prevalence=float(target[test].mean()))


def optimize_signed_gate(X, labels, features, populations, signs, beta=1., seed=42):
    """Optimize a user-defined AND gate while preserving every marker's polarity.

    Different coordinate-update orders are fit on training data and selected on
    validation. Fixed signs are biological hypotheses, not inferred populations.
    """
    from copy import deepcopy
    X=np.asarray(X);features=list(features);target=np.isin(labels,populations)
    if min(target.sum(),(~target).sum())<30:raise ValueError('Need 30 target and background events')
    if not signs or any(c not in features or op not in ('>','<=') for c,op in signs.items()):
        raise ValueError('Signed gate needs named channels and > or <= operators')
    train,rest=train_test_split(np.arange(len(X)),test_size=.4,stratify=target,random_state=seed)
    valid,test=train_test_split(rest,test_size=.5,stratify=target[rest],random_state=seed)
    conditions=[dict(channel=c,operator=op,threshold=float(X[train,features.index(c)].min())-1e-7 if op=='>' else float(X[train,features.index(c)].max())) for c,op in signs.items()]
    thresholds={c:np.unique(np.quantile(X[train,features.index(c)],np.linspace(0,1,65))) for c in signs}
    rng=np.random.default_rng(seed);candidates=[]
    for attempt in range(5):
        cs=deepcopy(conditions);order=rng.permutation(len(cs))
        best=metrics(membership(X[train],features,cs),target[train],beta)['fbeta']
        for _ in range(8):
            improved=False
            for j in order:
                base=membership(X[train],features,[v for k,v in enumerate(cs) if k!=j]);v=X[train,features.index(cs[j]['channel'])];t=target[train]
                th=np.unique(np.append(thresholds[cs[j]['channel']],cs[j]['threshold']))
                values=np.sort(v[base]);positive=np.sort(v[base&t])
                ng=np.searchsorted(values,th,side='right');tp=np.searchsorted(positive,th,side='right')
                if cs[j]['operator']=='>':ng=len(values)-ng;tp=len(positive)-tp
                score=(1+beta**2)*tp/np.maximum(beta**2*t.sum()+ng,1)
                k=int(score.argmax())
                if score[k]>best+1e-9:
                    cs[j]['threshold']=float(th[k]);best=float(score[k]);improved=True
            if not improved:break
        candidates.append(dict(method='constrained_sequential',channels=list(signs),conditions=cs,
          update_order=[list(signs)[j] for j in order],train=metrics(membership(X[train],features,cs),target[train],beta),
          validation=metrics(membership(X[valid],features,cs),target[valid],beta)))
    candidates.sort(key=lambda c:-c['validation']['fbeta']);winner=deepcopy(candidates[0]);mask=membership(X,features,winner['conditions'])
    winner.update(test=metrics(mask[test],target[test],beta),all_events=metrics(mask,target,beta),
                  shuffled_test_labels=metrics(mask[test],rng.permutation(target[test]),beta))
    return dict(populations=list(populations),objective='F-beta with fixed marker polarities',beta=beta,seed=seed,
                split_unit='Stratified events (60/20/20); no biological replication',n_candidates=len(candidates),
                best=winner,candidates=candidates,test_prevalence=float(target[test].mean()),
                limitation='User-defined marker polarities, thresholds optimized against transferred labels. Gates and flow placement do not establish stem-cell function; low/negative thresholds require instrument controls for experimental sorting.')
