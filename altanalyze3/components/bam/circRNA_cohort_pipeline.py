"""Complete-roster circRNA abundance and native Figure 5 statistics.

Requires completed raw-scan coverage and evidence ledger. All BSJs remain in the
feature universe. Kept-status source-query counts are the baseline statistical
layer; all statuses, unknown cells and all50 sample-timepoints are also exported.
Read capture counts are not asserted to be independent original RNA molecules.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
import re
import itertools
from pathlib import Path
import time
import sqlite3
import resource
import sys
import numpy as np
import pandas as pd
from scipy import sparse
from .circRNA_native_statistics import native_differential, NORMALIZATION


def rc(barcode):
    core, sep, suffix = barcode.partition('-')
    return core.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1] + (sep + suffix if sep else '')


def barcode_state(barcode, annotation, orientation=None):
    if barcode == '__AMBIGUOUS_BARCODE__':return barcode
    direct, reverse = annotation.get(barcode), annotation.get(rc(barcode))
    if orientation is not None:
        chosen=reverse if orientation=='reverse' else direct
        if chosen:return chosen
        if direct or reverse:return '__OPPOSITE_ONLY_BARCODE__'
        return '__UNMAPPED_BARCODE__'
    if direct and reverse and direct != reverse:
        return '__AMBIGUOUS_BARCODE__'
    return direct or reverse or '__UNMAPPED_BARCODE__'


def source_group(row):
    """Literal mutation labels; never propagate mutation through unlabelled timepoints."""
    group = row['group']
    if 'iPSC' in group or 'MISTRG' in group or group == 'Anemia':
        return None
    if group in ('young','aged'):
        return group
    if group.startswith('MDS-'):
        return 'MDS-SRSF2' if 'SRSF2' in group else 'MDS-others'
    return 'AML-SRSF2' if group == 'AML-SRSF2' else None


def build_contrast_roster(samples, patients):
    groups = collections.defaultdict(set)
    for row in samples:
        group = source_group(row)
        if group:
            groups[group].add(patients[row['library']]['patient'])
    healthy = groups['young'] | groups['aged']
    contrasts = [('MDS-others_vs_healthy', groups['MDS-others'], healthy),
                 ('MDS-SRSF2_vs_healthy', groups['MDS-SRSF2'], healthy),
                 ('all-SRSF2_vs_healthy', groups['MDS-SRSF2'] | groups['AML-SRSF2'], healthy),
                 ('aged_vs_young', groups['aged'], groups['young'])]
    for name, cases, controls in contrasts:
        if cases & controls:
            raise ValueError('Overlapping biological donors in ' + name)
    return contrasts


def contrast_library_scope(name, samples):
    allowed = {'MDS-others_vs_healthy': {'MDS-others','young','aged'},
               'MDS-SRSF2_vs_healthy': {'MDS-SRSF2','young','aged'},
               'all-SRSF2_vs_healthy': {'MDS-SRSF2','AML-SRSF2','young','aged'},
               'aged_vs_young': {'aged','young'}}[name]
    return [r['library'] for r in samples if source_group(r) in allowed]


def write_tsv(path, rows, fields=None):
    path = Path(path)
    op = gzip.open if path.suffix == '.gz' else open
    with op(path, 'wt') as h:
        if fields is None:
            rows = list(rows)
            fields = list(rows[0]) if rows else []
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t'); w.writeheader();w.writerows(rows)


def plot_examples(out, selected, matrices, depths, contrasts, feature_meta, threshold):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for n,(state,contrast,event,why,statsrow) in enumerate(selected):
        if (state,contrast) not in matrices: continue
        values=matrices[state,contrast].loc[event]
        groups=depths[state,contrast]['groups']
        fig,ax=plt.subplots(figsize=(7,6.8))
        vals=[values[groups['grp']==g].dropna().values for g in ('control','case')]
        ax.boxplot([v if len(v) else [np.nan] for v in vals],whis=1.5,widths=.45,showfliers=False,
                   patch_artist=True,boxprops={'facecolor':'#eeeeee','edgecolor':'black'},medianprops={'color':'black'})
        for i,(g,v) in enumerate(zip(('control','case'),vals),1):
            ax.scatter(np.full(len(v),i)+np.linspace(-.065,.065,len(v)),v,s=36,c='black',zorder=3)
            ax.text(i,1.02,f'n={len(v)}/{int((groups["grp"]==g).sum())}',ha='center',transform=ax.get_xaxis_transform())
        ax.set_xticks([1,2],['Young' if contrast=='aged_vs_young' else 'Healthy',contrast.split('_vs_')[0]])
        ax.set_ylabel('log2(1 + BSJ support per million mapped primary records)')
        ax.margins(y=.20)
        md=feature_meta[event]
        genes=md.get('display_genes') or (md.get('donor_region','')+' / '+md.get('acceptor_region',''))[:80]
        state_label='Sample-wide' if state=='ALL_SAMPLE' else state
        ax.set_title(f'{genes}\n{state_label}\n{event}',fontsize=10)
        p=statsrow.get('ebayesPval');q=statsrow.get('ebayesAdjPval')
        text=f'Native P={p:.3g}; BH q={q:.3g}' if pd.notna(p) else 'Comparison unavailable'
        ax.text(.02,.98,text,transform=ax.transAxes,va='top',fontsize=9)
        ax.spines[['top','right']].set_visible(False)
        test_caption=('Full all-SRSF2 inference withheld: AML7 clinical independence is unresolved; no P/q generated. ' if contrast=='all-SRSF2_vs_healthy' else 'Native Figure5 limma-like two-sided test, fixed 20% variance shrinkage; BH over every testable BSJ in this state×contrast, no significance prefilter. ')
        caption=(f'{contrast}; {state}. Each point is one pooled primary-human donor; timepoints and technical sources '
                 f'pool counts and depth. Boxes median/IQR, whiskers 1.5×IQR. Kept-status distinct query/event within each source; '
                 f'all mapped primary records in matching cells form denominator. {test_caption}'
                 f'Selected {why}; label-blind display minimum {threshold} source-query reads; all candidates retained in tables. '
                 f'Healthy=young+aged except aging control=young; MDS-other=non-SRSF2-labelled primary MDS, all timepoints. No covariate adjustment. Undefined depth is missing, never zero. Capture/RT bias and unresolved fallback donor identities limit inference; '
                 f'backsplice support is not biochemical proof of circularity.')
        import textwrap
        fig.text(.07,.02,textwrap.fill(caption,100),ha='left',va='bottom',fontsize=7.5)
        fig.subplots_adjust(left=.15,right=.96,top=.80,bottom=.32)
        fig.savefig(out/f'example_{n+1:03d}.pdf'); fig.savefig(out/f'example_{n+1:03d}.png',dpi=140);plt.close(fig)


def run(args):
    started=time.time();root=Path(args.run_root); evidence=Path(args.evidence);out=Path(args.output)
    if out.exists() and any(out.iterdir()): raise FileExistsError('Refusing existing analysis outputs')
    samples=list(csv.DictReader(open(root/'config/samples50.tsv'),delimiter='\t'))
    if len(samples)!=50 or len({r['library'] for r in samples})!=50:raise ValueError('Independent all50 source roster required')
    libidx={r['library']:i for i,r in enumerate(samples)}
    validated=json.loads((root/'VALIDATED_EVIDENCE.json').read_text())
    evidence_manifest=json.loads((evidence/'MANIFEST.json').read_text())
    expected_libraries=[r['library'] for r in samples]
    if validated['expected_ordered_libraries']!=expected_libraries:raise ValueError('Validated exact50 library roster changed')
    if len(evidence_manifest['inputs'])!=50 or {Path(r['path']).parent.name for r in evidence_manifest['inputs']}!=set(expected_libraries):raise ValueError('Evidence aggregation does not cover all50 sources')
    if evidence_manifest['validated_evidence_sha256']!=hashlib.sha256((root/'VALIDATED_EVIDENCE.json').read_bytes()).hexdigest():raise ValueError('Evidence was aggregated under a different validation gate')
    def digest(path):
        h=hashlib.sha256()
        with open(path,'rb') as hh:
            for block in iter(lambda:hh.read(1024*1024),b''):h.update(block)
        return h.hexdigest()
    validated_libraries={r['library']:r for r in validated['libraries']}
    if len(validated['libraries'])!=50 or set(validated_libraries)!=set(expected_libraries):raise ValueError('Duplicate/missing validated libraries')
    for artifact in ['events.tsv.gz','read_support.tsv.gz','classes.tsv.gz']:
        if digest(evidence/artifact)!=evidence_manifest['outputs'][artifact]['sha256']:raise ValueError('Evidence changed since aggregation: '+artifact)
    for r in samples:
        d=root/'samples'/r['library']
        if not (d/'STAGE_NONCOLINEAR_COMPLETE.json').exists():
            # Exact scanner completion file may use the source stage naming.
            choices=list(d.glob('*COMPLETE.json'))
            if len(choices)!=1:raise ValueError('Missing/unresolved completion for '+r['library'])
        stats=d/(r['library']+'.noncolinear_stats.json')
        if not stats.exists():raise ValueError('Missing per-library stats '+r['library'])
        if not (d/(r['library']+'.barcode_denominators.tsv.gz')).exists():raise ValueError('Missing denominator '+r['library'])
        for filename,kind in [(r['library']+'.barcode_denominators.tsv.gz','denominators'),(r['library']+'.noncolinear_stats.json','stats'),('STAGE_NONCOLINEAR_COMPLETE.json','stage')]:
            if digest(d/filename)!=validated_libraries[r['library']]['files'][kind]['sha256']:raise ValueError('Validated input changed: '+r['library']+' '+kind)
    annotations=collections.defaultdict(dict)
    with open(args.annotation) as h:
        for line in h:
            key,state=line.rstrip('\n').split('\t')[:2];bc,lib=key.split('.',1)
            if lib not in libidx:raise ValueError('Annotation contains unexpected library '+lib)
            if bc in annotations[lib] and annotations[lib][bc]!=state:raise ValueError('Conflicting cell state '+key)
            annotations[lib][bc]=state
    realstates=sorted({s for m in annotations.values() for s in m.values()})
    early_decisions=json.loads(Path(args.decisions).read_text())
    if realstates!=early_decisions['expected_state_order']:raise ValueError('Exact independent state roster changed')
    for source,key in [(Path(args.annotation),'annotation_sha256'),(Path(args.patients),'patients_sha256'),(root/'config/samples50.tsv','samples50_sha256')]:
        if digest(source)!=early_decisions[key]:raise ValueError('Pinned source changed: '+key)
    if len(realstates)!=89:raise ValueError('Cell-state identity scope changed: review exact source')
    states=['ALL_SAMPLE']+realstates+['__UNMAPPED_BARCODE__','__AMBIGUOUS_BARCODE__','__OPPOSITE_ONLY_BARCODE__']
    stateidx={s:i for i,s in enumerate(states)}
    D=np.zeros((len(states),len(samples)),dtype=np.int64);mapping_rows=[];orientations={}
    for lib,i in libidx.items():
        path=root/'samples'/lib/(lib+'.barcode_denominators.tsv.gz')
        tally=collections.Counter()
        with gzip.open(path,'rt') as h:denominator_rows=list(csv.DictReader(h,delimiter='\t'))
        annotation=annotations[lib]
        forward=sum(r['barcode'] in annotation for r in denominator_rows)
        reverse=sum(rc(r['barcode']) in annotation for r in denominator_rows)
        orientations[lib]='reverse' if reverse>forward else 'forward'
        conflicting=sum(r['barcode'] in annotation and rc(r['barcode']) in annotation and annotation[r['barcode']]!=annotation[rc(r['barcode'])] for r in denominator_rows)
        for row in denominator_rows:
            if row['library']!=lib:raise ValueError('Denominator library identity mismatch')
            n=int(row['mapped_primary_reads']); s=barcode_state(row['barcode'],annotation,orientations[lib])
            D[0,i]+=n;D[stateidx[s],i]+=n;tally[s]+=n
        info=json.loads((root/'samples'/lib/(lib+'.noncolinear_stats.json')).read_text())
        if info.get('library')!=lib or info.get('bam')!=samples[i]['bam'] or info.get('regions') is not None:raise ValueError('Scan source identity/whole-BAM scope mismatch '+lib)
        # Counter is nested in scanner stats; support exact verified schema.
        expected=info.get('stats',{}).get('primary_records',info.get('primary_records'))
        if expected is None:expected=info.get('counters',{}).get('primary_records')
        if expected is None or D[0,i]!=expected:raise ValueError('Denominator completeness mismatch '+lib)
        mapping_rows.append(dict(library=lib,metadata_reverse=samples[i]['reverse'],chosen_orientation=orientations[lib],forward_matching_barcodes=forward,reverse_matching_barcodes=reverse,both_orientation_conflicting_states=conflicting,opposite_only_records=tally['__OPPOSITE_ONLY_BARCODE__'],all_mapped_primary_records=int(D[0,i]),annotated_mapped_primary_records=sum(v for k,v in tally.items() if not k.startswith('__')),unmapped_records=tally['__UNMAPPED_BARCODE__'],ambiguous_records=tally['__AMBIGUOUS_BARCODE__']))
    patient_rows=list(csv.DictReader(open(args.patients),delimiter='\t'))
    patients={r['library']:r for r in patient_rows}
    if len(patients)!=len(patient_rows):raise ValueError('Duplicate patient-map library identities')
    for r in samples:
        if r['library'] not in patients:
            if source_group(r) is not None:raise ValueError('Inferential library missing from verified patient map: '+r['library'])
            patients[r['library']]={'patient':r['original_uid'],'basis':'separate descriptive nonprimary-human/anemia sample; no independent-patient inference'}
    for uid in {r['original_uid'] for r in samples}:
        pt={patients[r['library']]['patient'] for r in samples if r['original_uid']==uid}
        if len(pt)!=1:raise ValueError('Patient identities conflict within original_uid '+uid)
    gene_symbols={}
    if getattr(args,'gene_annotations',None):
        with open(args.gene_annotations) as gh:
            for line in gh:
                parts=line.rstrip('\n').split('\t')
                if len(parts)>1 and parts[0].startswith('ENSG'):gene_symbols[parts[0]]=parts[1]
    feature_meta={}
    with gzip.open(evidence/'events.tsv.gz','rt') as h:
        for r in csv.DictReader(h,delimiter='\t'):
            if r['joint_type']=='back':
                if r['event'] in feature_meta:raise ValueError('Duplicate feature identifier')
                tokens=[]
                for region in (r.get('donor_region',''),r.get('acceptor_region','')):
                    for part in region.split('|'):
                        gene=part.split(':')[0]
                        if gene:tokens.append(gene_symbols.get(gene.split('.')[0],gene))
                r['display_genes']=' / '.join(dict.fromkeys(tokens))
                feature_meta[r['event']]=r
    class_memberships=collections.defaultdict(set)
    with gzip.open(evidence/'classes.tsv.gz','rt') as h:
        for r in csv.DictReader(h,delimiter='\t'):
            if r['event'] in feature_meta:class_memberships[r['event']].add(r['joint_class'])
    known_junctions=set()
    if getattr(args,'literature',None):
        for lit in csv.DictReader(open(args.literature),delimiter='\t'):
            try:
                lo=int(lit['GRCh38_start_1based']);hi=int(lit['GRCh38_end_1based']);chrom=lit['reported_start'].split(':')[0]
                donor,acceptor=(hi,lo) if lit['strand']=='+' else (lo,hi)
                known_junctions.add('BSJ|'+'|'.join(map(str,[chrom,lit['strand'],lit['strand'],donor,acceptor])))
            except (ValueError,KeyError):pass
    for event,meta in feature_meta.items():
        classes=class_memberships[event]
        meta['all_joint_classes']=';'.join(sorted(classes))
        other_origin=any('mitochondrial' in c or 'lariat' in c for c in classes)
        eligible=not other_origin and ('backsplice_registered_canonical' in classes or event in known_junctions)
        meta['circRNA_display_priority']=eligible
        meta['literature_exact_boundary_match']=event in known_junctions
        meta['evidence_label']='mitochondrial_or_lariat_like_backjunction' if other_origin else ('registered_canonical_or_literature_matched_candidate' if eligible else 'other_backjunction_candidate')
    features=list(feature_meta);fi={e:i for i,e in enumerate(features)}
    display_feature=np.array([feature_meta[e]['circRNA_display_priority'] for e in features],dtype=bool)
    # Sparse counters preserve all rows; absent numerical entries are observed zero
    # only because the exact50 completed source scope passed above.
    out.mkdir(parents=True,exist_ok=True)
    def memory_report(stage):
        peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        gib=peak/(1024**3) if sys.platform=='darwin' else peak/(1024**2)
        with open(out/'memory_profile.jsonl','a') as mh:mh.write(json.dumps(dict(stage=stage,peak_rss_GiB=gib,elapsed_seconds=time.time()-started,features=len(features)))+'\n')
        print('MEMORY',stage,round(gib,3),'GiB',flush=True)
    memory_report('features_loaded')
    counts={s:collections.Counter() for s in states};raw=collections.Counter();status=collections.Counter()
    flags=collections.defaultdict(collections.Counter);sample_counts=collections.Counter()
    querydb=sqlite3.connect(out/'cross_library_query_identity_diagnostic.sqlite')
    querydb.execute('CREATE TABLE ledger(event TEXT,query_name TEXT,library TEXT)')
    querybatch=[]
    uid_of={r['library']:r['original_uid'] for r in samples}
    with gzip.open(evidence/'read_support.tsv.gz','rt') as h:
        reader=csv.DictReader(h,delimiter='\t')
        for _,support_group in itertools.groupby(reader,key=lambda r:(r['event'],r['library'],r['status'],r['query_name'])):
            variants=list(support_group);r=dict(variants[0]);barcodes={v['barcode'] for v in variants}
            if len(barcodes)>1:r['barcode']='__AMBIGUOUS_BARCODE__'
            for flag in ('contained','multipass','polyA','homology_above_6','positive_query_gap','unresolved'):
                r[flag]=str(any(v.get(flag,'').lower() in ('true','1') for v in variants))
            e=r['event']
            if e not in fi:continue
            lib=r['library'];i=libidx[lib];f=fi[e];s=barcode_state(r['barcode'],annotations[lib],orientations[lib])
            raw[f,i]+=1;status[r['status']]+=1
            for flag in ('contained','multipass','polyA','homology_above_6','positive_query_gap','unresolved'):
                if r.get(flag,'').lower() in ('true','1'):
                    flags[e]['all_status_'+flag]+=1
                    if r['status']=='kept':flags[e]['kept_'+flag]+=1
            querybatch.append((e,r['query_name'],lib))
            if len(querybatch)>=10000:
                querydb.executemany('INSERT INTO ledger VALUES (?,?,?)',querybatch);querybatch=[]
            if r['status']=='kept':
                counts['ALL_SAMPLE'][f,i]+=1;counts[s][f,i]+=1
                sample_counts[f,uid_of[lib]]+=1
    out.mkdir(parents=True,exist_ok=True)
    write_tsv(out/'library_barcode_mapping.tsv',mapping_rows)
    write_tsv(out/'complete_library_state_depth.tsv.gz',(dict(library=r['library'],sample_uid=r['original_uid'],state=s,mapped_primary_records=int(D[stateidx[s],i])) for s in states for i,r in enumerate(samples)))
    write_tsv(out/'complete_feature_roster.tsv.gz',(feature_meta[e] for e in features))
    write_tsv(out/'complete_sample_roster.tsv',[dict(r,patient=patients[r['library']]['patient'],patient_basis=patients[r['library']]['basis']) for r in samples])
    querydb.executemany('INSERT INTO ledger VALUES (?,?,?)',querybatch);querydb.commit()
    write_tsv(out/'cross_library_query_identifier_overlap.tsv.gz',(dict(event=e,query_name=q,libraries=ls,interpretation='identifier overlap only; library-local IDs do not establish molecular duplication') for e,q,ls in querydb.execute('SELECT event,query_name,group_concat(DISTINCT library) FROM ledger GROUP BY event,query_name HAVING count(DISTINCT library)>1')),fields=['event','query_name','libraries','interpretation'])
    querydb.close()
    memory_report('support_ledger_aggregated')
    # Label-blind empirical display heuristic; never a validity or testing gate.
    observed=np.array(list(sample_counts.values()),dtype=float)
    threshold=max(3,int(math.ceil(np.quantile(observed,.75)))) if len(observed) else 3
    sensitivity=[]
    for t in sorted(set([1,2,3,5,10,20,50,threshold])):
        recurring=collections.Counter(f for (f,u),n in sample_counts.items() if n>=t)
        sensitivity.append(dict(min_reads_per_sample=t,events_at_least_one_sample=len(recurring),events_at_least_two_samples=sum(n>=2 for n in recurring.values())))
    write_tsv(out/'display_support_threshold_sensitivity.tsv',sensitivity)
    contrasts=build_contrast_roster(samples,patients)
    memberships=[]
    for name,cases,controls in contrasts:
        scope=set(contrast_library_scope(name,samples))
        for r in samples:
            pt=patients[r['library']]['patient']
            memberships.append(dict(contrast=name,library=r['library'],original_uid=r['original_uid'],patient=pt,literal_group=r['group'],role=('case' if pt in cases else 'control') if r['library'] in scope else 'outside_primary_inferential_scope'))
    write_tsv(out/'contrast_library_membership.tsv',memberships)
    def barcode_membership_rows():
        for r in samples:
            lib=r['library']
            with gzip.open(root/'samples'/lib/(lib+'.barcode_denominators.tsv.gz'),'rt') as bh:
                for br in csv.DictReader(bh,delimiter='\t'):
                    rawbc=br['barcode'];chosen=rc(rawbc) if orientations[lib]=='reverse' else rawbc
                    yield dict(library=lib,original_uid=r['original_uid'],raw_barcode=rawbc,annotation_barcode=chosen,state=barcode_state(rawbc,annotations[lib],orientations[lib]),orientation=orientations[lib],mapped_primary_records=br['mapped_primary_reads'])
    write_tsv(out/'barcode_state_membership.tsv.gz',barcode_membership_rows(),fields=['library','original_uid','raw_barcode','annotation_barcode','state','orientation','mapped_primary_records'])
    decision=json.loads(Path(args.decisions).read_text())
    selected=[];plot_matrices={};plot_depths={};best_examples={}
    summary_counts=[]
    identity_flags=[dict(library=r['library'],sample_uid=r['original_uid'],patient=patients[r['library']]['patient'],basis=patients[r['library']].get('basis','')) for r in samples if 'unverified' in patients[r['library']].get('basis','')]
    write_tsv(out/'donor_identity_flags.tsv',identity_flags)
    sample_uids=list(dict.fromkeys(r['original_uid'] for r in samples))
    sample_pool=np.zeros((50,len(sample_uids)),dtype=np.int8)
    for i,r in enumerate(samples):sample_pool[i,sample_uids.index(r['original_uid'])]=1
    sample_metadata=[]
    for uid in sample_uids:
        rows=[r for r in samples if r['original_uid']==uid]
        labels=sorted({r['group'] for r in rows})
        sample_metadata.append(dict(sample_uid=uid,groups=';'.join(labels),libraries=';'.join(r['library'] for r in rows),patient=patients[rows[0]['library']]['patient'],model_type='iPSC' if any('iPSC' in g for g in labels) else ('MISTRG' if any('MISTRG' in g for g in labels) else 'primary_human'),literal_SRSF2=any('SRSF2' in g for g in labels),healthy=labels in [['young'],['aged']]))
    write_tsv(out/'complete_sample_timepoint_roster.tsv',sample_metadata)
    inclusive_rows=[]
    for s in states:
        c=counts[s]
        if c:
            keys=list(c);M=sparse.coo_matrix(([c[k] for k in keys],([k[0] for k in keys],[k[1] for k in keys])),shape=(len(features),50)).tocsr()
        else:M=sparse.csr_matrix((len(features),50),dtype=np.int64)
        sparse.save_npz(out/f'counts_{stateidx[s]:03d}.npz',M)
        pooled_samples=M@sparse.csr_matrix(sample_pool)
        sparse.save_npz(out/f'original_uid_counts_{stateidx[s]:03d}.npz',pooled_samples)
        np.save(out/f'original_uid_depth_{stateidx[s]:03d}.npy',D[stateidx[s]]@sample_pool)
        if s=='ALL_SAMPLE':
            totals_for_display=np.asarray(M.sum(axis=1)).ravel()
            possible=np.flatnonzero(display_feature)
            top_idx=possible[np.argsort(totals_for_display[possible])[-10:][::-1]]
            sample_depth=D[0]@sample_pool
            topcounts=pooled_samples[top_idx].toarray()
            for k,f in enumerate(top_idx):
                for u,meta in enumerate(sample_metadata):
                    if meta['literal_SRSF2'] or meta['healthy']:
                        n=int(topcounts[k,u]);d=int(sample_depth[u])
                        inclusive_rows.append(dict(event=features[f],**meta,supporting_source_queries=n,mapped_primary_records=d,per_million=(n/d*1e6 if d else None),inference='descriptive only; model/treatment/donor confounding retained'))
        if s.startswith('__'):continue
        for name,case_donors,control_donors in contrasts:
            scope=contrast_library_scope(name,samples)
            units=sorted(control_donors)+sorted(case_donors)
            group=pd.DataFrame({'grp':['control' if u in control_donors else 'case' for u in units],'biological_donor':units},index=units)
            A=np.zeros((50,len(units)),dtype=np.int8)
            for lib in scope:A[libidx[lib],units.index(patients[lib]['patient'])]=1
            pooled=M@A; depth=pd.Series(D[stateidx[s]]@A,index=units)
            frame=pd.DataFrame(pooled,index=features,columns=units)
            unavailable_reason='UNAVAILABLE_UNRESOLVED_DONOR_IDENTITY' if name=='all-SRSF2_vs_healthy' else None
            result,normalized=native_differential(frame,depth,group,features,units,'case','control',decision,inference_unavailable_reason=unavailable_reason)
            support=(frame>=threshold).sum(axis=1)
            result['donors_with_display_support']=support
            result['case_donors_with_display_support']=(frame.loc[:,group['grp']=='case']>=threshold).sum(axis=1)
            result['control_donors_with_display_support']=(frame.loc[:,group['grp']=='control']>=threshold).sum(axis=1)
            result['total_kept_source_query_support']=frame.sum(axis=1)
            result['state']=s;result['contrast']=name
            result.index.name='event'
            numeric_columns=[c for c in result.columns if pd.api.types.is_numeric_dtype(result[c])]
            np.savez_compressed(out/f'stats_{stateidx[s]:03d}_{name}.npz',values=result[numeric_columns].to_numpy(dtype=float),columns=np.asarray(numeric_columns,dtype=str),log2_abundance=normalized.to_numpy(dtype=float),units=np.asarray(units,dtype=str),tested=result['ebayesPval'].notna().to_numpy(),unavailable_reason=np.asarray(unavailable_reason or ''),feature_order_sha256=np.asarray(hashlib.sha256('\n'.join(features).encode()).hexdigest()))
            # Full statistics and donor measurements retain feature order from the independent roster.
            # NPZ avoids hundreds of millions of repeated identifier and method-text strings.
            result.loc[(support>=2)].sort_values('ebayesPval').head(100).to_csv(out/f'discovery_preview_{stateidx[s]:03d}_{name}.tsv.gz',sep='\t',compression='gzip')
            pd.DataFrame({'unit':units,'group':group['grp'].values,'mapped_primary_records':depth.values}).to_csv(out/f'donor_depth_{stateidx[s]:03d}_{name}.tsv',sep='\t',index=False)
            v1=normalized.loc[:,group['grp']=='case'].var(axis=1,ddof=1);v2=normalized.loc[:,group['grp']=='control'].var(axis=1,ddof=1)
            pooled_variance=(v1+v2)/2
            testable=result['ebayesPval'].notna()
            prior=float(pooled_variance[testable].median()) if testable.any() else None
            summary_counts.append(dict(inference_status=unavailable_reason or 'NATIVE_TEST_AVAILABLE_WITH_RECORDED_LIMITATIONS',native_variance_prior=prior,zero_variance_fraction=float((pooled_variance[testable]==0).mean()) if testable.any() else None,native_sparse_calibration_warning=bool(prior==0),state=s,contrast=name,features=len(features),testable=int(result['ebayesPval'].notna().sum()),nominal_p_lt_0_05=int((result['ebayesPval']<.05).sum()),bh_q_lt_0_05=int((result['ebayesAdjPval']<.05).sum()),case_donors=len(case_donors),control_donors=len(control_donors)))
            eligible=result[(support>=2)&display_feature&(result['case_minus_control_mean_log2_CPM1'].notna() if unavailable_reason else result['ebayesPval'].notna())]
            if len(eligible):
                ranked=eligible['case_minus_control_mean_log2_CPM1'].abs().sort_values(ascending=False)
                e=ranked.index[0];score=float(ranked.iloc[0])
                old=best_examples.get(name)
                if old is None or score>old[0]:
                    best_examples[name]=(score,s,e,result.loc[e].to_dict(),normalized.loc[[e]].copy(),group.copy())
            if s=='ALL_SAMPLE':
                for e in result.loc[display_feature,'total_kept_source_query_support'].nlargest(2).index:
                    selected.append((s,name,e,'highest pooled kept source-query support among registered-canonical or exact-literature candidates',result.loc[e].to_dict()))
                plot_matrices[s,name]=normalized.loc[[e for ss,nn,e,*_ in selected if ss==s and nn==name]].copy()
                plot_depths[s,name]={'groups':group}
            del frame,pooled,result,normalized
        memory_report('state_complete:'+s)
        print('completed state',s,flush=True)
    for name,(score,s,e,row,values,group) in best_examples.items():
        why='largest absolute descriptive mean-log2-abundance difference; inference withheld' if name=='all-SRSF2_vs_healthy' else 'largest absolute native mean-log2-abundance difference'
        selected.append((s,name,e,why+' among recurrent supported registered-canonical or exact-literature candidates',row))
        if (s,name) in plot_matrices:
            plot_matrices[s,name]=pd.concat([plot_matrices[s,name],values]).loc[lambda x:~x.index.duplicated()]
        else:plot_matrices[s,name]=values
        plot_depths[s,name]={'groups':group}
    write_tsv(out/'comparison_summary.tsv',summary_counts)
    write_tsv(out/'all_SRSF2_inclusive_and_healthy_sampletime_examples.tsv',inclusive_rows,fields=['event','sample_uid','groups','libraries','patient','model_type','literal_SRSF2','healthy','supporting_source_queries','mapped_primary_records','per_million','inference'])
    totals=counts['ALL_SAMPLE'];ranks=[];sample_recurrence=collections.Counter(f for (f,uid),n in sample_counts.items() if n>0)
    for e,f in fi.items():
        n=sum(totals[f,i] for i in range(50));nr=sum(raw[f,i] for i in range(50))
        ranks.append(dict(feature_meta[e],kept_source_queries=n,all_status_source_queries=nr,pooled_per_million_mapped_primary=(n/float(D[0].sum())*1e6 if D[0].sum() else None),sample_timepoints_detected=sample_recurrence[f],source_libraries_detected=sum(totals[f,i]>0 for i in range(50)),**flags[e]))
    ranks.sort(key=lambda r:-r['kept_source_queries']);write_tsv(out/'all_candidates_ranked.tsv.gz',ranks,fields=list(dict.fromkeys(k for r in ranks for k in r)))
    permutation_diagnostics=[]
    for state,name,e,why,r in selected:
        if name=='all-SRSF2_vs_healthy':
            permutation_diagnostics.append(dict(state=state,contrast=name,event=e,exact_two_sided_mean_difference_P=None,permutations=0,interpretation='UNAVAILABLE_UNRESOLVED_DONOR_IDENTITY; no inferential test performed'))
            continue
        values=plot_matrices[state,name].loc[e];gg=plot_depths[state,name]['groups']['grp']
        values=values.dropna();gg=gg.loc[values.index];ncase=int((gg=='case').sum())
        vals=values.to_numpy();actual=float(vals[gg.to_numpy()=='case'].mean()-vals[gg.to_numpy()=='control'].mean()) if ncase and ncase<len(vals) else None
        if actual is not None:
            perms=[]
            for choice in itertools.combinations(range(len(vals)),ncase):
                mask=np.zeros(len(vals),dtype=bool);mask[list(choice)]=True
                perms.append(float(vals[mask].mean()-vals[~mask].mean()))
            pp=float(np.mean(np.abs(perms)>=abs(actual)-1e-12))
        else:pp=None;perms=[]
        permutation_diagnostics.append(dict(state=state,contrast=name,event=e,exact_two_sided_mean_difference_P=pp,permutations=len(perms),interpretation='selection-conditioned diagnostic on selected plotted features; no FDR or independent-validation claim'))
    write_tsv(out/'selected_example_permutation_diagnostics.tsv',permutation_diagnostics,fields=['state','contrast','event','exact_two_sided_mean_difference_P','permutations','interpretation'])
    plot_examples(out,selected,plot_matrices,plot_depths,contrasts,feature_meta,threshold)
    import matplotlib.pyplot as plt
    for number,event in enumerate(list(dict.fromkeys(r['event'] for r in inclusive_rows))[:2],1):
        subset=[r for r in inclusive_rows if r['event']==event]
        labels=['Healthy human','SRSF2 human','SRSF2 MISTRG','SRSF2 iPSC'];vals=[]
        for label in labels:
            vs=[math.log2(1+r['per_million']) for r in subset if r['per_million'] is not None and
                (('Healthy' in label and r['healthy']) or ('SRSF2' in label and r['literal_SRSF2'] and
                 r['model_type']==({'SRSF2 human':'primary_human','SRSF2 MISTRG':'MISTRG','SRSF2 iPSC':'iPSC'}[label])))]
            vals.append(vs)
        fig,ax=plt.subplots(figsize=(8,6))
        ax.boxplot([v if v else [np.nan] for v in vals],whis=1.5,showfliers=False)
        for i,v in enumerate(vals,1):ax.scatter(np.full(len(v),i)+np.linspace(-.07,.07,len(v)),v,s=36)
        ax.set_xticks(range(1,5),labels);ax.set_ylabel('log2(1 + supporting source queries per million primary records)')
        ax.set_title('All SRSF2-labelled sample timepoints, separated by model\n'+event,fontsize=10)
        fig.text(.08,.02,'Descriptive only. Points are sample timepoints, including repeated donors; these are not independent replicates.\nBoxes median/IQR, whiskers 1.5×IQR. Healthy includes young+aged. Model, therapy and capture effects are confounded.\nAll literal SRSF2 labels retained; no mutation inference for unlabelled timepoints. No pooled significance test.',fontsize=8)
        fig.subplots_adjust(bottom=.25,top=.84,left=.12);fig.savefig(out/f'all_SRSF2_inclusive_{number:02d}.pdf');fig.savefig(out/f'all_SRSF2_inclusive_{number:02d}.png',dpi=140);plt.close(fig)
    literature_path=Path(args.literature) if args.literature else None
    if literature_path:
        literature_rows=[]
        for lit in csv.DictReader(open(literature_path),delimiter='\t'):
            try:
                lo=int(lit['GRCh38_start_1based']);hi=int(lit['GRCh38_end_1based']);chrom=lit['reported_start'].split(':')[0]
                donor,acceptor=(hi,lo) if lit['strand']=='+' else (lo,hi)
                event='BSJ|'+'|'.join(map(str,[chrom,lit['strand'],lit['strand'],donor,acceptor]))
            except (ValueError,KeyError):
                literature_rows.append(dict(lit,event='',match_status='COORDINATES_UNRESOLVED_RETAINED'));continue
            for i,sample in enumerate(samples):
                f=fi.get(event);n=counts['ALL_SAMPLE'].get((f,i),0) if f is not None else 0
                alln=raw.get((f,i),0) if f is not None else 0
                literature_rows.append(dict(lit,event=event,match_status='EXACT_ORIENTED_BOUNDARY_MATCH' if f is not None else 'NOT_OBSERVED_IN_COMPLETE_SCANNED_BAMS',library=sample['library'],sample_uid=sample['original_uid'],kept_source_queries=n,all_status_source_queries=alln,mapped_primary_records=int(D[0,i])))
        write_tsv(out/'literature_exact_junction_all50.tsv',literature_rows,fields=list(dict.fromkeys(k for r in literature_rows for k in r)))
    write_tsv(out/'selected_ISV_example_candidates.tsv',[dict(state=s,contrast=n,event=e,selection=why,case_minus_control_mean_log2_CPM1=r.get('case_minus_control_mean_log2_CPM1'),P=r.get('ebayesPval'),BH_q=r.get('ebayesAdjPval')) for s,n,e,why,r in selected])
    manifest=dict(status='COMPLETE',libraries=50,features=len(features),states=states,normalization=NORMALIZATION,
                  display_threshold=threshold,display_threshold_rule='max(3,ceil(75th percentile of positive kept event-by-original_uid support)); label blind display only; not circularity/FDR calibration',
                  status_counts=dict(status),statistical_unit='primary-human donor pooled within each literal contrast; no genotype propagation',donor_identity_flags=identity_flags,identity_limitation='Existing patient map preserved; original_uid fallbacks remain flagged. Full all-SRSF2 inference is UNAVAILABLE_UNRESOLVED_DONOR_IDENTITY because AML7 clinical independence was not established; all its feature counts, units and descriptive plots retained without P/q.',
                  sparse_counts_semantics='feature rows follow complete_feature_roster.tsv.gz; columns follow complete_sample_roster.tsv; implicit zeros mean no kept supporting query after complete raw scan; missing depth yields missing abundance',
                  correction='BH all testable features per state per contrast',denominator_caveat='ALL_SAMPLE denominator includes all mapped-primary records, including missingCB/lowMAPQ; baseline kept numerator requires its recorded CB/MAPQ eligibility. This estimates captured retained-BSJ support per total mapped record, not absolute RNA abundance. Cell-state scopes use matching globally-oriented barcodes; opposite-only and conflicting identities remain explicit.',source_hashes={'annotation':digest(Path(args.annotation)),'patients':digest(Path(args.patients)),'samples50':digest(root/'config/samples50.tsv'),'decisions':digest(Path(args.decisions)),'pipeline':digest(Path(__file__))},elapsed_seconds=time.time()-started)
    (out/'COHORT_ANALYSIS_COMPLETE.json').write_text(json.dumps(manifest,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-root',required=True);p.add_argument('--evidence',required=True);p.add_argument('--output',required=True);p.add_argument('--annotation',required=True);p.add_argument('--patients',required=True);p.add_argument('--decisions',required=True);p.add_argument('--literature');p.add_argument('--gene-annotations')
    run(p.parse_args())
if __name__=='__main__':main()
