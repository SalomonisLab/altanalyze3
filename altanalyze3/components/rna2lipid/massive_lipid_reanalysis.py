"""Re-extract MSV000081973 MS1 apex intensities using published target annotations.

LIQUID composition/fragment rules define masses. RT locking uses spectra only,
never the published abundance patterns. Isobaric target annotations are excluded
from individual-species validation. This is a reproducible re-extraction, not a
claim that every original manual LIQUID identification has been reproduced.
"""
from __future__ import annotations

import argparse
import ast
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ATOMS={"C":12.0,"H":1.00782503223,"N":14.00307400443,
       "O":15.99491461957,"P":30.97376199842,"S":31.9720711744}
ELECTRON=0.000548579909065


def expression(text,x=0,y=0):
    if not text.strip():return 0
    text=re.sub(r"(\d)([XY])",r"\1*\2",text)
    def visit(node):
        if isinstance(node,ast.Constant) and isinstance(node.value,int):return node.value
        if isinstance(node,ast.Name) and node.id in ['X','Y']:return {'X':x,'Y':y}[node.id]
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return -visit(node.operand)
        if isinstance(node,ast.BinOp):
            a,b=visit(node.left),visit(node.right)
            if isinstance(node.op,ast.Add):return a+b
            if isinstance(node.op,ast.Sub):return a-b
            if isinstance(node.op,ast.Mult):return a*b
        raise ValueError('Unsupported composition expression '+text)
    return visit(ast.parse(text,mode='eval').body)


def chains(lipid):
    lipid=re.sub(r'_[A-Z]$','',lipid.split(';')[0].strip())
    if lipid=='CoQ10':return []
    body=lipid[lipid.index('(')+1:-1]
    result=[]
    for token in body.split('/'):
        m=re.fullmatch(r'(O-|P-|d|t|m)?(\d+):(\d+)',token.strip())
        if not m:raise ValueError('Unsupported chain annotation '+token)
        if int(m[2])==0:continue
        prefix=m[1] or '';kind={'O-':'Ether','P-':'Plasmalogen','d':'Dihydro','t':'Trihydro','m':'Monohydro'}.get(prefix,'Standard')
        result.append((int(m[2]),int(m[3]),kind))
    return result


def read_rules(path):
    with open(path) as handle:return list(csv.DictReader(handle,delimiter='\t'))


def composition(lipid,rules):
    if ';' in lipid:
        fs=[composition(l,rules) for l in lipid.split(';')]
        if any(f!=fs[0] for f in fs):raise ValueError('Ambiguous annotations have unequal masses: '+lipid)
        return fs[0]
    lipid=re.sub(r'_[A-Z]$','',lipid.strip())
    chain=chains(lipid);cls=lipid.split('(')[0]
    cls={'HexCer':'GlcCer','CoQ10':'Ubiquinone'}.get(cls,cls)
    chain_count=({'TG':3,'DG':2,'SM':2,'Cer':2,'GM3':2}.get(cls,len(chain))
                 if len(chain)==1 else len(chain))
    flags={'NumChains':chain_count,'ContainsEther':int(any(c[2]=='Ether' for c in chain)),
           'ContainsDiether':0,'ContainsPlasmalogen':int(any(c[2]=='Plasmalogen' for c in chain)),
           'ContainsLCB':int(any(c[2]=='Dihydro' for c in chain) or cls=='GM3'),
           'ContainsLCB+OH':int(any(c[2]=='Trihydro' for c in chain)),
           'ContainsLCB-OH':int(any(c[2]=='Monohydro' for c in chain)),
           'IsOxoCHO':0,'IsOxoCOOH':0,'NumOH':0,'ContainsOOH':0,'ContainsF2IsoP':0}
    selected=[]
    for r in rules:
        r={k.strip():v for k,v in r.items()}
        if r['LipidClass']!=cls:continue
        if lipid=='CoQ10' and r['Formula-NoAdduct']!='C59H90O4':continue
        if any(int(r[k] or 0)!=v for k,v in flags.items()):continue
        f={a:expression(r[a],sum(c[0] for c in chain),sum(c[1] for c in chain)) for a in ATOMS}
        if f not in selected:selected.append(f)
    if len(selected)!=1:raise ValueError(f'{lipid}: {len(selected)} distinct composition rules')
    return selected[0]


def mass(formula):return sum(ATOMS[a]*n for a,n in formula.items())


def constraint(text,value):
    text=(text or '').replace(' ','')
    if not text:return True
    if text.startswith('>'):return value>int(text[1:])
    if text.startswith('<'):return value<int(text[1:])
    return value==int(text)


def fragments(target,rules):
    chain=chains(target['lipid']);cls=target['lipid'].split('(')[0]
    cls={'HexCer':'GlcCer','CoQ10':'Ubiquinone'}.get(cls,cls)
    result=[]
    for r in rules:
        if r['lipidClass']!=cls or r['fragmentationMode'].lower()!=target['polarity']:continue
        if r['additional_element']:continue
        if not constraint(r['countOfChains'],len(chain)):continue
        if not constraint(r['countOfStandardAcylsChains'],sum(c[2]=='Standard' for c in chain)):continue
        if any(r.get(k) for k in ['containsHydroxy','sialic','acylChain.NumCarbons','acylChain.NumDoubleBonds','acylChain.HydroxyPosition','target_acylChains']):continue
        kind=r['acylChain.AcylChainType']
        sum_only=len(chain)==1 and cls in ['TG','DG','SM','GM3']
        if kind and sum_only:continue
        options=[None] if not kind else [c for c in chain if kind=='All' or c[2]==kind]
        for c in options:
            x,y=(c[0],c[1]) if c else (0,0)
            f={a:expression(r[a],x,y) for a in ATOMS}
            m=target['mz']-mass(f) if r['neutral_loss']=='1' else mass(f)+(-ELECTRON if target['polarity']=='positive' else ELECTRON)
            if m>0:result.append((m,r['desc'],f'{c[0]}:{c[1]}:{c[2]}' if c else '',r['diagnostic']=='1'))
    return result


def targets_from_published(published,rules_dir):
    cr=read_rules(rules_dir/'DefaultCompositionRules.txt');fr=read_rules(rules_dir/'DefaultFragmentationRules.txt')
    targets=[]
    for feature in pd.read_csv(published,nrows=0).columns[1:]:
        lipid,adduct=feature.split('|');f=composition(lipid,cr)
        shift={'[M+H]+':ATOMS['H']-ELECTRON,'[M-H]-':-ATOMS['H']+ELECTRON,
               '[M+NH4]+':ATOMS['N']+4*ATOMS['H']-ELECTRON}[adduct]
        row={'feature_id':feature,'lipid':lipid,'adduct':adduct,'polarity':'positive' if adduct.endswith('+') else 'negative',
             'formula':''.join(a+(str(n) if n!=1 else '') for a,n in f.items() if n),'mz':mass(f)+shift}
        row['fragments']=fragments(row,fr);targets.append(row)
    for r in targets:
        r['isobaric_targets']=';'.join(t['feature_id'] for t in targets if t is not r and t['polarity']==r['polarity'] and abs(t['mz']-r['mz'])<=r['mz']*10e-6)
        if ';' in r['lipid']:r['isobaric_targets']+=';composite_published_annotation'
    return targets


def discover(files,targets,scan_reader,mass_signal,ppm=5):
    bins=defaultdict(lambda:defaultdict(float));evidence=[];file_qc=[]
    for path in files:
        mode='positive' if '_POS_' in path.name else 'negative';ts=[t for t in targets if t['polarity']==mode]
        trace=defaultdict(list);n1=n2=0
        for scan in scan_reader(path):
            if scan.ms_level==1:
                n1+=1
                for t in ts:
                    value=mass_signal(scan,t['mz'],ppm)
                    if value>0:trace[t['feature_id']].append((scan.rt_seconds,value))
            elif scan.ms_level==2:
                n2+=1
                if scan.precursor_mz is None:continue
                for t in ts:
                    if abs(scan.precursor_mz-t['mz'])>t['mz']*10e-6:continue
                    found=[]
                    for fm,description,chain,diagnostic in t['fragments']:
                        # Original acquisition includes low-resolution CID as well as HCD.
                        signal=mass_signal(scan,fm,0.3/fm*1e6)
                        if signal>0 and signal>=max(scan.intensity,default=0)*.01:
                            found.append((description,chain,diagnostic))
                    distinct={f[1] for f in found if f[1]};diagnostics=sum(f[2] for f in found)
                    if found:
                        score=len(distinct)*2+int(diagnostics>0)
                        if score:bins[t['feature_id']][round(scan.rt_seconds/30)]+=score
                        evidence.append({'file':path.name,'feature_id':t['feature_id'],'scan_id':scan.scan_id,
                                         'rt_seconds':scan.rt_seconds,'precursor_mz':scan.precursor_mz,
                                         'matched_fragments':len(found),'matched_chain_count':len(distinct),
                                         'diagnostic_fragments':diagnostics,'score':score,
                                         'matched_labels':';'.join(sorted(set(f[0] for f in found)))})
        file_qc.append({'file':path.name,'ms1_scans':n1,'ms2_scans':n2})
        # Store spectrum-only RT candidates even when no diagnostic MS2 exists;
        # these remain marked unconfirmed and cannot establish identity.
        for fid,values in trace.items():
            rt,value=max(values,key=lambda p:p[1]);bins[fid][round(rt/30)]+=0.01
        print('scanned',path.name,n1,n2,flush=True)
    locked=[]
    for t in targets:
        b=bins[t['feature_id']]
        if not b:continue
        rtbin=max(b,key=lambda k:(b[k]+b.get(k-1,0)+b.get(k+1,0),-k))
        support=[r for r in evidence if r['feature_id']==t['feature_id'] and abs(r['rt_seconds']-rtbin*30)<=45]
        locked.append({k:v for k,v in t.items() if k!='fragments'}|{'rt_seconds':rtbin*30,
                      'ms2_support_scans':len(support),'max_matched_chains':max((r['matched_chain_count'] for r in support),default=0),
                      'annotation_status':'published_target_with_MS2_support' if support else 'accurate_mass_candidate_only'})
    return locked,evidence,file_qc


def run(args):
    sys.path.insert(0,str(Path(args.pyneoquant)/'src'))
    from pyneoquant.quant.lipid_ms1 import iter_scans,mass_signal,quantify
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    files=sorted(Path(args.mzml_dir).glob('*_L_*.mzML'))
    if len(files)!=30:raise ValueError(f'Expected 30 lipid runs, found {len(files)}')
    targets=targets_from_published(args.published,Path(args.rules_dir))
    locked,evidence,qc=discover(files,targets,iter_scans,mass_signal,args.ppm)
    pd.DataFrame(locked).to_csv(out/'locked_targets.tsv',sep='\t',index=False)
    pd.DataFrame(evidence).to_csv(out/'MS2_target_evidence.csv',index=False)
    pd.DataFrame(qc).to_csv(out/'spectrum_inventory.csv',index=False)
    rows=[]
    for path in files:
        m=re.search(r'D(\d+)_(END|EPI|MES|MIC|PMX)_',path.name);sample=f'D{int(m[1]):03d}_{m[2]}'
        mode='positive' if '_POS_' in path.name else 'negative'
        result=quantify(path,[r for r in locked if r['polarity']==mode],args.ppm,args.rt_window,3)
        rows.extend(dict(r,sample=sample,file=path.name) for r in result)
        print('quantified',sample,mode,flush=True)
    data=pd.DataFrame(rows);data.to_csv(out/'peak_quantification.csv',index=False)
    apex=data.pivot(index='sample',columns='feature_id',values='apex_intensity')
    apex.to_csv(out/'raw_apex_intensities.csv');log=np.log2(apex);log.to_csv(out/'raw_apex_log2.csv')
    # Combined ion-mode sample-median normalization follows the stated paper.
    normalized=log.sub(log.median(axis=1),axis=0)
    centered=normalized.sub(normalized.mean(axis=0),axis=1)
    centered.to_csv(out/'replicated_median_lipid_centered_log2.csv')
    pub=pd.read_csv(args.published,index_col=0)
    pub.index=[re.sub(r'^D0*(\d+)',lambda m:f'D{int(m[1]):03d}',s) for s in pub.index]
    shared=sorted(set(centered.columns)&set(pub.columns));samples=sorted(set(centered.index)&set(pub.index))
    diagnostics=[];contrasts=[]
    for fid in shared:
        target=next(r for r in locked if r['feature_id']==fid)
        valid=centered.loc[samples,fid].notna()&pub.loc[samples,fid].notna()
        a=centered.loc[samples,fid][valid];b=pub.loc[samples,fid][valid]
        diagnostics.append({'feature_id':fid,'observed_profiles':int(valid.sum()),
                            'rmse_log2':float(np.sqrt(np.mean((a-b)**2))),
                            'spearman':float(spearmanr(a,b).statistic) if len(a)>2 else None,
                            'isobaric_targets':target['isobaric_targets'],'ms2_support_scans':target['ms2_support_scans'],
                            'max_matched_chains':target['max_matched_chains']})
        for donor in ['D001','D008','D011']:
            for pop in ['END','EPI','MES','MIC']:
                s,p=donor+'_'+pop,donor+'_PMX'
                observed=float(centered.loc[s,fid]-centered.loc[p,fid]);published=float(pub.loc[s,fid]-pub.loc[p,fid])
                contrasts.append({'feature_id':fid,'sample':s,'raw_log2FC':observed,'published_log2FC':published,
                                  'unambiguous_published_target':not bool(target['isobaric_targets']),
                                  'ms2_supported':target['ms2_support_scans']>0})
    checks=pd.DataFrame(diagnostics);checks.to_csv(out/'per_lipid_replication.csv',index=False)
    c=pd.DataFrame(contrasts);c.to_csv(out/'published_contrast_replication.csv',index=False)
    selected=c[c['unambiguous_published_target']&c['ms2_supported']].dropna()
    strong=selected[selected['published_log2FC'].abs()>=.5]
    result={'status':'MS1_reextraction_completed_original_manual_identification_not_fully_replicated',
            'runs':len(files),'published_targets':len(targets),'locked_targets':len(locked),
            'individual_species_evaluation_excludes_isobaric_annotations':True,
            'evaluated_contrasts':len(selected),'evaluated_lipids':int(selected.feature_id.nunique()),
            'contrast_correlation':float(selected.raw_log2FC.corr(selected.published_log2FC)),
            'contrast_rmse_log2':float(np.sqrt(np.mean((selected.raw_log2FC-selected.published_log2FC)**2))),
            'strong_contrast_direction_accuracy':float(np.mean(np.sign(strong.raw_log2FC)==np.sign(strong.published_log2FC))),
            'MS2_fragment_tolerance_da':.3,'MS2_relative_intensity_threshold':.01,
            'ppm':args.ppm,'rt_half_window_seconds':args.rt_window,
            'not_established':['original manual RT/isomer assignments','absolute molar amounts','full-panel raw replication'],
            'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            'software':'pyNeoQuant targeted centroid MS1 extractor + LIQUID composition/fragment rules',
            'normalization':'log2 MS1 apex; subtract sample median; subtract lipid mean across 15 profiles',
            'RT_selection':'pooled MS2 fragment evidence; MS1-only fallback explicitly flagged; no published abundance pattern used'}
    (out/'replication_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='input_sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mzml-dir',required=True);p.add_argument('--published',required=True)
    p.add_argument('--rules-dir',required=True);p.add_argument('--out',required=True)
    p.add_argument('--pyneoquant',default='/Users/saljh8/Documents/GitHub/pyNeoQuant')
    p.add_argument('--ppm',type=float,default=5);p.add_argument('--rt-window',type=float,default=45)
    run(p.parse_args())
