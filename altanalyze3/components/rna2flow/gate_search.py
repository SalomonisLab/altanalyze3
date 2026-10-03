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


def parent_population(X, features, target, conditions=None, mask=None):
    """Fixed measured gates; annotation exclusions are never used as sorting conditions."""
    parent=np.ones(len(X),bool);audit=[]
    if mask is not None:
        mask=np.asarray(mask)
        if mask.shape!=(len(X),) or mask.dtype!=bool:
            raise ValueError('Parent mask must be a boolean value for every event')
        parent &= mask
        audit.append(dict(stage='Existing parent gate',**metrics(parent,target)))
    for c in conditions or []:
        if c.get('channel') not in features:
            raise ValueError('Parent marker unavailable in this panel: '+str(c.get('channel')))
        parent &= membership(X,features,[c])
        audit.append(dict(stage=c.get('stage',c['channel']),condition=dict(c),**metrics(parent,target)))
    return parent,audit


def validation_frontier(candidates):
    """Indices not dominated on validation purity, recovery and marker count."""
    frontier=[]
    for i,c in enumerate(candidates):
        a=(c['validation']['precision'],c['validation']['recall'],-len({v['channel'] for v in c['conditions']}))
        dominated=False
        for j,d in enumerate(candidates):
            if i==j:continue
            b=(d['validation']['precision'],d['validation']['recall'],-len({v['channel'] for v in d['conditions']}))
            if all(x>=y for x,y in zip(b,a)) and any(x>y for x,y in zip(b,a)):
                dominated=True;break
        if not dominated:frontier.append(i)
    return frontier


def _greedy_candidates(X, target, features, beta, max_markers, bins=100):
    """Retention-constrained greedy threshold search, without reimplementing MarkerFinder.

    Scan all measured channels, select purity-improving cuts retaining 95% of
    current targets, and retain intermediate strategies for validation selection.
    """
    current=np.ones(len(X),bool);conditions=[];results=[]
    for step in range(12):
        nt=int(target[current].sum());ng=int(current.sum())
        if nt<5 or ng==0:break
        baseline=nt/ng;best=None;used={c['channel'] for c in conditions}
        for j,ch in enumerate(features):
            if ch not in used and len(used)>=max_markers:continue
            v=X[current,j];t=target[current];order=np.argsort(v,kind='stable')
            values=v[order];positives=np.cumsum(t[order]);cuts=np.unique(np.quantile(v,np.linspace(0,1,bins+1)))
            n=np.searchsorted(values,cuts,side='right');tp=np.where(n>0,positives[np.maximum(n-1,0)],0)
            for op,num,pos in [('<=',n,tp),('>',len(v)-n,nt-tp)]:
                purity=pos/np.maximum(num,1);eligible=(pos>=.95*nt)&(num<len(v))&(num>0)&(purity>baseline+1e-9)
                if not eligible.any():continue
                k=int(np.argmax(np.where(eligible,purity,-1)))
                score=(float(purity[k]),int(pos[k]),-int(num[k]))
                if best is None or score>best[0]:best=(score,dict(channel=ch,operator=op,threshold=float(cuts[k])))
        if best is None:break
        conditions.append(best[1]);current &= membership(X,features,[best[1]])
        results.append(dict(method='retention_greedy',channels=list(dict.fromkeys(c['channel'] for c in conditions)),
                            conditions=[dict(c) for c in conditions],per_step_target_retention=.95))
        if target[current].sum()/max(target.sum(),1)<.5:break
    return results


def _refine_candidates(X, target, features, seeds, beta, max_markers):
    """Train-only coordinate relaxation, pruning and forward expansion of AND gates.

    Relaxing earlier cuts can restore rare targets after later exclusions remove
    background. Both target and whole-population quantiles resolve rare-state tails.
    No validation or test values influence threshold proposals or path ordering.
    """
    thresholds={ch:np.unique(np.r_[np.quantile(X[:,j],np.linspace(0,1,65)),
                    np.quantile(X[target,j],np.linspace(0,1,65))]) for j,ch in enumerate(features)}
    def score(cs):return metrics(membership(X,features,cs),target,beta)['fbeta']
    ordered=sorted(seeds,key=lambda c:score(c['conditions']),reverse=True)
    starts=[[]];seen=set()
    for c in ordered:
        signature=tuple(sorted((v['channel'],v['operator']) for v in c['conditions']))
        if signature in seen:continue
        seen.add(signature);starts.append([dict(v) for v in c['conditions']])
        if len(starts)>=5:break
    out=[]
    for cs in starts:
        best=score(cs)
        for iteration in range(6):
            changed=False
            # Re-optimize existing cuts while all other cuts remain applied.
            for k in range(len(cs)-1,-1,-1):
                others=cs[:k]+cs[k+1:];base=membership(X,features,others)
                removal=score(others)
                if removal>=best-1e-10:
                    cs=others;best=removal;changed=True;continue
                c=cs[k];j=features.index(c['channel']);th=np.unique(np.r_[thresholds[c['channel']],c['threshold']])
                values=np.sort(X[base,j]);positive=np.sort(X[base&target,j]);ng=np.searchsorted(values,th,side='right');tp=np.searchsorted(positive,th,side='right')
                if c['operator']=='>':ng=len(values)-ng;tp=len(positive)-tp
                f=(1+beta**2)*tp/np.maximum(beta**2*target.sum()+ng,1);idx=int(f.argmax())
                if f[idx]>best+1e-9:c['threshold']=float(th[idx]);best=float(f[idx]);changed=True
            base=membership(X,features,cs);used={v['channel'] for v in cs};addition=None
            if len(cs)<2*max_markers:
                for j,ch in enumerate(features):
                    if ch not in used and len(used)>=max_markers:continue
                    values=np.sort(X[base,j]);positive=np.sort(X[base&target,j]);th=thresholds[ch]
                    ng=np.searchsorted(values,th,side='right');tp=np.searchsorted(positive,th,side='right')
                    for op,n,p in [('<=',ng,tp),('>',len(values)-ng,len(positive)-tp)]:
                        if any(c['channel']==ch and c['operator']==op for c in cs):continue
                        f=(1+beta**2)*p/np.maximum(beta**2*target.sum()+n,1);idx=int(f.argmax())
                        if f[idx]>best+1e-9 and (addition is None or f[idx]>addition[0]):
                            addition=(float(f[idx]),dict(channel=ch,operator=op,threshold=float(th[idx])))
            if addition:best=addition[0];cs.append(addition[1]);changed=True
            if cs:out.append(dict(method='coordinate_refined',channels=list(dict.fromkeys(c['channel'] for c in cs)),conditions=[dict(c) for c in cs]))
            if not changed:break
    return out


def serial_audit(X, features, target, conditions, parent_mask=None, beta=1.):
    """Descriptive serial gate counts; order changes intermediate, not final AND membership."""
    target=np.asarray(target,bool);current=np.ones(len(X),bool);rows=[]
    baseline=float(target.mean())
    def stage(name, condition=None):
        nonlocal current
        before=current.copy();nt=int((before&target).sum())
        if condition is not None:current &= membership(X,features,[condition])
        elif parent_mask is not None:current &= parent_mask
        row=dict(stage=name,condition=condition,**metrics(current,target,beta))
        row.update(n_entering=int(before.sum()),target_entering=nt,
                   target_lost=int((before&target&~current).sum()),
                   stage_recovery=float((current&target).sum()/nt) if nt else 0.,
                   enrichment=row['precision']/baseline if baseline else 0.)
        rows.append(row)
    rows.append(dict(stage='Starting population',condition=None,**metrics(current,target,beta),
                     n_entering=len(X),target_entering=int(target.sum()),target_lost=0,stage_recovery=1.,enrichment=1.))
    if parent_mask is not None:stage('Published / existing capture parent')
    for i,c in enumerate(conditions):stage(c.get('stage',f"{c['channel']} {c['operator']} {c['threshold']:.4g}"),c)
    return rows


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
                        max_markers=6, bins=24, groups=None, parent_conditions=None,
                        parent_mask=None, include_greedy=True, search_features=None, refine=True, include_hypergate=False, min_purity=0., min_recovery=0.):
    X=np.asarray(X);labels=np.asarray(labels).astype(str);features=list(features)
    if X.ndim!=2 or X.shape!=(len(labels),len(features)) or not np.isfinite(X).all():
        raise ValueError('Expected finite events × named markers and one label per event')
    if not .25<=beta<=4:raise ValueError('beta must be between 0.25 and 4')
    if not 0<=min_purity<=1 or not 0<=min_recovery<=1:raise ValueError('Purity and recovery constraints must be fractions from 0 to 1')
    populations=list(populations)
    missing=set(populations)-set(labels)
    if missing:raise ValueError('No events assigned to: '+', '.join(sorted(missing)))
    target=np.isin(labels,populations)
    if not isinstance(max_markers,int) or not 1<=max_markers<=20:
        raise ValueError('max_markers must be an integer between 1 and 20')
    parent,audit=parent_population(X,features,target,parent_conditions,parent_mask)
    search_features=features if search_features is None else list(search_features)
    if not set(search_features).issubset(features) or len(search_features)!=len(set(search_features)):
        raise ValueError('Search features must be unique available markers')
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
    full_valid=valid.copy();full_test=test.copy()
    train=train[parent[train]];valid=valid[parent[valid]];test=test[parent[test]]
    if any(min(target[z].sum(),(~target[z]).sum())<5 for z in [train,valid,test]):
        raise ValueError(f'Parent gate leaves insufficient target/background in a split; retains {int((target&parent).sum())}/{int(target.sum())} target events. Inspect target loss or change the parent')
    # Binary target uses all background, rather than dropping small competing populations.
    binary=np.where(target,'target','background')
    rho=rank_channels_for_population(X[train],binary[train],features,'target',min_cells=5)
    ranked=rho.fillna(0).abs().sort_values(ascending=False).index.tolist()
    chosen=[c for c in ranked if c in search_features and np.ptp(X[train,features.index(c)])>0][:max_markers]
    if len(chosen)<1:raise ValueError('Need an informative measured channel')
    candidates=[]
    if include_greedy:
        candidates.extend(_greedy_candidates(X[train][:,[features.index(c) for c in search_features]],target[train],search_features,beta,max_markers))
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
    if include_hypergate:
        from .hypergate import hypergate_conditions
        search_X=X[train][:,[features.index(c) for c in search_features]]
        cs=hypergate_conditions(search_X,target[train],search_features,beta)
        original_markers=len({c['channel'] for c in cs})
        # Honor the user's marker budget using training scores only.
        while len({c['channel'] for c in cs})>max_markers:
            options=[[c for c in cs if c['channel']!=ch] for ch in dict.fromkeys(c['channel'] for c in cs)]
            cs=max(options,key=lambda v:metrics(membership(search_X,search_features,v),target[train],beta)['fbeta'])
        candidates.append(dict(method='HyperGate_official' if original_markers<=max_markers else 'HyperGate_budget_pruned',
                               channels=list(dict.fromkeys(c['channel'] for c in cs)),conditions=cs,
                               original_marker_count=original_markers))
    if refine:
        candidates.extend(_refine_candidates(X[train][:,[features.index(c) for c in search_features]],target[train],search_features, candidates,beta,max_markers))
    for c in candidates:
        c['train']=metrics(membership(X[train],features,c['conditions']),target[train],beta)
        c['validation_within_parent']=metrics(membership(X[valid],features,c['conditions']),target[valid],beta)
        c['validation']=metrics(parent[full_valid]&membership(X[full_valid],features,c['conditions']),target[full_valid],beta)
        c['conditions']=[dict(v,fixed_parent=True) for v in (parent_conditions or [])]+c['conditions']
    def meets(c):return c['validation']['precision']>=min_purity and c['validation']['recall']>=min_recovery
    candidates.sort(key=lambda c:(not meets(c),-c['validation']['fbeta'],len(c['conditions'])))
    winner=dict(candidates[0]);mask=parent&membership(X,features,winner['conditions'])
    winner['test_within_parent']=metrics(mask[test],target[test],beta)
    winner['test']=metrics(mask[full_test],target[full_test],beta)
    winner['all_events']=metrics(mask,target,beta)
    null=np.random.default_rng(seed).permutation(target[full_test])
    winner['shuffled_test_labels']=metrics(mask[full_test],null,beta)
    return dict(populations=populations,objective='F-beta',beta=beta,seed=seed,
                split_unit=unit,n_candidates=len(candidates),best=winner,candidates=candidates,
                marker_ranking=[dict(channel=c,rho=float(rho[c])) for c in ranked],
                limitation='Capture and purity are against annotation membership. The test events are withheld from gate fitting and selection, but labels may have been inferred from these same antibodies. No functional stem-cell validation.',
                test_prevalence=float(target[full_test].mean()),parent_audit=audit,
                pareto_candidate_indices=validation_frontier(candidates),
                parent_summary=metrics(parent,target),
                parent_warning='Parent removes target cells; recovery includes these losses.' if (target&~parent).any() else '',
                constraints=dict(min_purity=min_purity,min_recovery=min_recovery,met=meets(winner)),
                validation_basis='Full starting population, including targets removed by fixed parents')


def optimize_signed_gate(X, labels, features, populations, signs, beta=1., seed=42,
                         parent_conditions=None,parent_mask=None,min_purity=0.,min_recovery=0.):
    """Optimize a user-defined AND gate while preserving every marker's polarity.

    Different coordinate-update orders are fit on training data and selected on
    validation. Fixed signs are biological hypotheses, not inferred populations.
    """
    from copy import deepcopy
    X=np.asarray(X);features=list(features);target=np.isin(labels,populations)
    parent,audit=parent_population(X,features,target,parent_conditions,parent_mask)
    if min(target.sum(),(~target).sum())<30:raise ValueError('Need 30 target and background events')
    if not signs or any(c not in features or op not in ('>','<=') for c,op in signs.items()):
        raise ValueError('Signed gate needs named channels and > or <= operators')
    train,rest=train_test_split(np.arange(len(X)),test_size=.4,stratify=target,random_state=seed)
    valid,test=train_test_split(rest,test_size=.5,stratify=target[rest],random_state=seed)
    full_valid=valid.copy();full_test=test.copy()
    train=train[parent[train]];valid=valid[parent[valid]];test=test[parent[test]]
    if any(min(target[z].sum(),(~target[z]).sum())<5 for z in [train,valid,test]):
        raise ValueError(f'Parent gate leaves insufficient target/background in a split; retains {int((target&parent).sum())}/{int(target.sum())} target events. Inspect target loss or change the parent')
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
          validation_within_parent=metrics(membership(X[valid],features,cs),target[valid],beta),
          validation=metrics(parent[full_valid]&membership(X[full_valid],features,cs),target[full_valid],beta)))
    for c in candidates:c['conditions']=[dict(v,fixed_parent=True) for v in (parent_conditions or [])]+c['conditions']
    def meets(c):return c['validation']['precision']>=min_purity and c['validation']['recall']>=min_recovery
    candidates.sort(key=lambda c:(not meets(c),-c['validation']['fbeta']));winner=deepcopy(candidates[0]);mask=parent&membership(X,features,winner['conditions'])
    winner.update(test=metrics(mask[full_test],target[full_test],beta),test_within_parent=metrics(mask[test],target[test],beta),all_events=metrics(mask,target,beta),
                  shuffled_test_labels=metrics(mask[full_test],rng.permutation(target[full_test]),beta))
    return dict(populations=list(populations),objective='F-beta with fixed marker polarities',beta=beta,seed=seed,
                split_unit='Stratified events (60/20/20); no biological replication',n_candidates=len(candidates),
                best=winner,candidates=candidates,test_prevalence=float(target[full_test].mean()),
                constraints=dict(min_purity=min_purity,min_recovery=min_recovery,met=meets(winner)),
                parent_audit=audit,parent_summary=metrics(parent,target),
                parent_warning='Parent removes target cells; recovery includes these losses.' if (target&~parent).any() else '',
                limitation='User-defined marker polarities, thresholds optimized against transferred labels. Gates and flow placement do not establish stem-cell function; low/negative thresholds require instrument controls for experimental sorting.')
