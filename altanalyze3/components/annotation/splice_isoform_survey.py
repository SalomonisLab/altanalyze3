"""Survey existing differential-splicing tests with native isoform read support.

Reuses the original statistics and ISV transcript index. No new P values, sequence
predictions or reference releases. Both annotated representatives and all expressed
exclusive junction carriers are retained; a gene-level match is not an isoform match.
"""
import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import sqlite3
import numpy as np
import pandas as pd
from scipy import sparse


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def dump(path,value):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'wt') as f:json.dump(value,f,indent=2,allow_nan=False)


def write(path,rows,fields=None):
    rows=iter(rows);first=next(rows,None)
    if first is None and not fields:raise ValueError('Empty table needs declared columns: '+str(path))
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'wt') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(first),delimiter='\t',lineterminator='\n');w.writeheader()
        if first is not None:w.writerow(first)
        w.writerows(rows)


def split_feature(feature):
    gene,sep,tid=feature.partition(':')
    if not sep or not gene or not tid:raise ValueError('Expected exact gene:isoform identity: '+feature)
    return gene,tid


def healthy_counts(counts,metadata,group):
    import anndata as ad
    from ..long_read.comparisons import _resolve_uid
    meta=pd.read_csv(metadata,sep='\t',dtype=str)
    pairs=meta[['uid','groups']].apply(lambda s:s.str.strip()).drop_duplicates()
    if pairs.uid.duplicated().any():raise ValueError('A donor UID belongs to multiple groups')
    mapping=dict(zip(pairs.uid,pairs.groups));a=ad.read_h5ad(counts)
    features=list(map(str,a.var_names));names=list(map(str,a.obs_names))
    if len(set(features))!=len(features) or len(set(names))!=len(names):raise ValueError('Duplicate count identity')
    roster=[]
    for name in names:
        state,sep,raw=name.partition('.')
        if not sep:raise ValueError('Pseudobulk identity needs state.sample: '+name)
        donor=_resolve_uid(raw,mapping,'isoform')
        if donor not in mapping:raise ValueError('Unresolved sample; retained and diagnosis required: '+name)
        roster.append(dict(pseudobulk=name,state=state,donor=donor,group=mapping[donor],selected=mapping[donor]==group))
    selected=[i for i,r in enumerate(roster) if r['selected']]
    if not selected:raise ValueError('No pseudobulks in requested group')
    x=sparse.csr_matrix(a.X[selected]);del a
    if not np.isfinite(x.data).all() or (x.data<0).any() or not np.equal(x.data,np.rint(x.data)).all():
        raise ValueError('Counts must be finite, nonnegative integer observations')
    states=sorted({roster[i]['state'] for i in selected});order={s:i for i,s in enumerate(states)}
    agg=sparse.csr_matrix((np.ones(len(selected)),([order[roster[i]['state']] for i in selected],range(len(selected)))),shape=(len(states),len(selected)))
    pooled=(agg@x).tocsr();maximum=np.asarray(pooled.max(axis=0).toarray()).ravel()
    donors={s:sorted({roster[i]['donor'] for i in selected if roster[i]['state']==s}) for s in states}
    return features,roster,states,donors,pooled,maximum


def scan_events(directory,states,q_threshold):
    pattern=re.compile(r'^healthy-(.+)-Others-\1_stats-annotated\.txt$')
    files={}
    for p in Path(directory).glob('*_stats-annotated.txt'):
        m=pattern.match(p.name)
        if not m:raise ValueError('Unexpected event filename: '+p.name)
        if m[1] in files:raise ValueError('Duplicate state event table')
        files[m[1]]=p
    if set(files)!=set(states):raise ValueError('Complete healthy state roster differs between event tests and counts: '+repr(set(files)^set(states)))
    events={};manifest=[]
    for state in states:
        p=files[state];rows=0;tested=0;significant=0;seen=set()
        with p.open() as f:
            reader=csv.DictReader(f,delimiter='\t');cols=reader.fieldnames
            fwd,rev='DiffMeans_'+state+'_m_Others','DiffMeans_Others_m_'+state
            if sum(c in cols for c in (fwd,rev))!=1:raise ValueError('Unresolved test orientation: '+str(p))
            effect=fwd if fwd in cols else rev
            for r in reader:
                rows+=1;event=r['Feature']
                if event in seen:raise ValueError('Duplicate event identity: '+event)
                seen.add(event)
                try:q=float(r['ebayesAdjPval']);raw=float(r[effect]);ms=float(r['Mean_'+state]);mo=float(r['Mean_Others'])
                except ValueError:continue  # original untested rows remain in their source file and file audit
                if not all(map(math.isfinite,[q,raw,ms,mo])):continue
                delta=ms-mo
                if not math.isclose(raw,delta if effect==fwd else -delta,rel_tol=1e-4,abs_tol=1e-5):
                    raise ValueError('Effect disagrees with source mean columns: '+str(p)+' '+event)
                tested+=1
                if q>=q_threshold:continue
                significant+=1
                entry=events.setdefault(event,dict(feature=event,gene=event.split(':',1)[0],gene_symbol=r['Gene_Symbol'],
                    event_type=r['Event_Type'],event_direction=r['Event_Direction'],coord_1=r['Coord_1'],coord_2=r['Coord_2'],
                    representatives=[],states={}))
                if (entry['event_type'],entry['coord_1'],entry['coord_2'])!=(r['Event_Type'],r['Coord_1'],r['Coord_2']):
                    raise ValueError('Conflicting source event annotation: '+event)
                reps=[r['Isoform_1|Length'],r['Isoform_2|Length']]
                if reps not in entry['representatives']:entry['representatives'].append(reps)
                entry['states'][state]=dict(dpsi=delta,q=q,mean_PSI=ms,mean_PSI_others=mo,
                    n_state=int(float(r['N_'+state])),n_others=int(float(r['N_Others'])))
        manifest.append(dict(state=state,path=str(p.resolve()),sha256=sha(p),rows=rows,tested=tested,
                             significant=significant,effect_column=effect,negated=effect==rev))
        print('[events]',state,rows,tested,significant,flush=True)
    return events,manifest


def carriers(events,features,maximum,threshold,index_path):
    wanted=defaultdict(dict)
    for i,feature in enumerate(features):
        if maximum[i]>threshold:
            gene,tid=split_feature(feature);wanted[gene][tid]=i
    needed=defaultdict(set)
    for event in events:
        for j in event.split('|'):
            gene,body=split_feature(j);needed[gene].add(body.split('=',1)[0])
    con=sqlite3.connect(Path(index_path).resolve().as_uri()+'?mode=ro',uri=True)
    result=defaultdict(set);structures={};unresolved=[]
    for gene,junctions in needed.items():
        found={}
        for _,strand,structure,tid,sample in con.execute('SELECT * FROM transcripts WHERE gene=?',(gene,)):
            if tid not in wanted.get(gene,{}):continue
            if tid in found and found[tid]!=structure:raise ValueError('Conflicting registered native token structure: '+gene+':'+tid)
            found[tid]=structure
        for tid,i in wanted.get(gene,{}).items():
            if tid not in found:unresolved.append(dict(feature=features[i],reason='exact native ID absent from transcript index'));continue
            tokens=found[tid].split('|');structures[features[i]]=tokens
            for a,b in zip(tokens,tokens[1:]):
                j=a+'-'+b
                if j in junctions:result[gene+':'+j].add(i)
    con.close()
    return result,structures,unresolved


def run(args):
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):raise ValueError('Use a fresh survey output directory')
    method=dict(statistics='Existing healthy state-versus-Others PSI tests; original per-file adjusted P values unchanged',
        q_threshold=args.q,read_threshold_strict_gt=args.min_reads,read_unit='Pooled assigned isoform reads per healthy cell state',
        pair_rule='Both native isoforms independently exceed threshold in at least one healthy state; states may differ',
        carrier_rule='Exact native transcript token adjacency from existing ISV index; opposite exclusive junction carriers',
        coding_filter=False,protein_length_filter=False,new_statistical_tests=False,
        pooled_dominance_reversal='Descriptive only; not a donor-paired or causal switch test')
    dump(out/'manifest.json',dict(status='incomplete',method=method))
    features,roster,states,donors,pooled,maximum=healthy_counts(args.counts,args.metadata,args.group)
    write(out/'pseudobulk_roster.tsv',roster);dump(out/'count_features.json.gz',features);dump(out/'states.json',states)
    sparse.save_npz(out/'healthy_state_counts.npz',pooled)
    # Indexed column access avoids rescanning the complete sparse matrix per link.
    columns=pooled.tocsc()
    events,sources=scan_events(args.events,states,args.q)
    links,structures,unresolved=carriers(events,features,maximum,args.min_reads,args.transcript_index)
    write(out/'unresolved_native_models.tsv',unresolved,fields=['feature','reason'])
    if unresolved:raise ValueError('Native structure reconciliation incomplete; see unresolved_native_models.tsv')
    findex={f:i for i,f in enumerate(features)};pairs=defaultdict(lambda:dict(events=set(),states=set()));event_rows=[]
    qualifying=set();linkrows=[]
    for event,e in events.items():
        js=event.split('|')
        if len(js)!=2:raise ValueError('Expected two competing junctions: '+event)
        one,two=(set(links.get(j.split('=',1)[0],set())) for j in js)
        left,right=one-two,two-one
        reps=[]
        for pair in e['representatives']:
            record=[]
            for rep in pair:
                tid=rep.split('|',1)[0];feature=e['gene']+':'+tid;i=findex.get(feature)
                record.append(dict(isoform_id=tid,count_feature=feature,present_in_counts=i is not None,
                    max_reads=None if i is None else int(maximum[i]),passes_read_rule=bool(i is not None and maximum[i]>args.min_reads)))
            reps.append(record)
        eligible=bool(left and right)
        if eligible:qualifying.add(event)
        best=max(e['states'],key=lambda s:abs(e['states'][s]['dpsi']))
        event_rows.append(dict(feature=event,gene=e['gene'],gene_symbol=e['gene_symbol'],event_type=e['event_type'],
            event_direction=e['event_direction'],coord_1=e['coord_1'],coord_2=e['coord_2'],
            expressed_exclusive_pair=eligible,n_side_1=len(left),n_side_2=len(right),n_both=len(one&two),
            min_q=min(v['q'] for v in e['states'].values()),strongest_state=best,dpsi=e['states'][best]['dpsi'],
            max_abs_dpsi_ge_0_1=abs(e['states'][best]['dpsi'])>=.1,n_significant_states=len(e['states']),
            source_representatives=json.dumps(reps,separators=(',',':')),state_tests=json.dumps(e['states'],separators=(',',':'))))
        for side,indices in [('1',left),('2',right),('both',one&two)]:
            for i in sorted(indices):linkrows.append(dict(event=event,gene_symbol=e['gene_symbol'],side=side,feature=features[i],
                max_reads=int(maximum[i]),max_read_state=states[int(columns[:,i].toarray().argmax())]))
        for i,j in itertools.product(left,right):
            if split_feature(features[i])[0]!=split_feature(features[j])[0]:
                raise ValueError('Competing event joins different native genes: '+event)
            key=tuple(sorted([i,j]));pairs[key]['events'].add(event);pairs[key]['states'].update(e['states'])
    write(out/'differential_event_catalog.tsv.gz',event_rows)
    write(out/'event_isoform_links.tsv.gz',linkrows,fields=['event','gene_symbol','side','feature','max_reads','max_read_state'])
    used=sorted({i for pair in pairs for i in pair});genes={split_feature(features[i])[0] for i in used}
    gene_cols=defaultdict(list)
    for i,f in enumerate(features):
        gene=split_feature(f)[0]
        if gene in genes:gene_cols[gene].append(i)
    totals={g:np.asarray(pooled[:,inds].sum(axis=1)).ravel() for g,inds in gene_cols.items()}
    counts={i:columns[:,i].toarray().ravel() for i in used};pairrows=[]
    for (i,j),p in pairs.items():
        gene=split_feature(features[i])[0];den=totals[gene]
        ui=np.divide(counts[i],den,out=np.full(len(states),np.nan),where=den>0)
        uj=np.divide(counts[j],den,out=np.full(len(states),np.nan),where=den>0);diff=ui-uj
        hi=int(np.nanargmax(diff));lo=int(np.nanargmin(diff))
        pairrows.append(dict(gene=gene,gene_symbol=events[next(iter(p['events']))]['gene_symbol'],isoform_1=features[i],isoform_2=features[j],
            max_reads_1=int(maximum[i]),max_reads_2=int(maximum[j]),max_read_state_1=states[int(counts[i].argmax())],
            max_read_state_2=states[int(counts[j].argmax())],dominance_reversal=bool(diff[hi]>0 and diff[lo]<0),
            state_favoring_1=states[hi],state_favoring_2=states[lo],usage_1_in_state_1=ui[hi],usage_2_in_state_1=uj[hi],
            usage_1_in_state_2=ui[lo],usage_2_in_state_2=uj[lo],reads_1_in_state_1=int(counts[i][hi]),reads_2_in_state_1=int(counts[j][hi]),
            reads_1_in_state_2=int(counts[i][lo]),reads_2_in_state_2=int(counts[j][lo]),
            events=';'.join(sorted(p['events'])),significant_event_states=';'.join(sorted(p['states']))))
    write(out/'native_isoform_pair_catalog.tsv.gz',pairrows)
    dump(out/'native_token_structures.json.gz',structures)
    write(out/'native_isoform_state_support.tsv.gz',(dict(feature=features[i],max_reads=int(maximum[i]),
        reads_by_state=json.dumps({s:int(v) for s,v in zip(states,counts[i])},separators=(',',':'))) for i in used))
    manifest=dict(status='complete',method=method,code_sha256=sha(__file__),states=states,donors_by_state=donors,
        healthy_donors=sorted({r['donor'] for r in roster if r['selected']}),healthy_pseudobulks=sum(r['selected'] for r in roster),
        input_feature_roster=len(features),native_high_read_features=int((maximum>args.min_reads).sum()),
        source_event_tables=sources,q_significant_events=len(events),events_with_read_supported_exclusive_pairs=len(qualifying),
        native_isoform_pairs=len(pairrows),genes_with_read_supported_events=len({events[e]['gene'] for e in qualifying}),
        inputs={str(Path(p).resolve()):sha(p) for p in [args.counts,args.metadata,args.transcript_index]},
        outputs={p.name:sha(p) for p in out.iterdir() if p.name!='manifest.json'})
    dump(out/'manifest.json',manifest);print(json.dumps({k:manifest[k] for k in ['status','healthy_pseudobulks','native_high_read_features','q_significant_events','events_with_read_supported_exclusive_pairs','native_isoform_pairs','genes_with_read_supported_events']},indent=2))


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--counts',required=True);p.add_argument('--metadata',required=True);p.add_argument('--group',default='healthy')
    p.add_argument('--events',required=True);p.add_argument('--transcript-index',required=True);p.add_argument('--outdir',required=True)
    p.add_argument('--min-reads',type=float,default=100);p.add_argument('--q',type=float,default=.05)
    args=p.parse_args(argv)
    if args.min_reads<0 or not 0<args.q<=1:raise ValueError('Invalid read/FDR threshold')
    run(args)


if __name__=='__main__':main()
