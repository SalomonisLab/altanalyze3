"""Reproducible literature-screen inputs for a native healthy splice survey.

Scans existing usage statistics without new tests. Retains every high-read native
feature, all significant reciprocal usage pairs, and explicit literature-screen
coverage. Protein sequences are literal registered predictions, never reconstructed.
This endpoint supplies evidence for curation; it does not infer a literature match
from a gene symbol or a protein length.
"""
import argparse
from collections import defaultdict
import csv
import gzip
import itertools
import json
import math
from pathlib import Path
import re
import sqlite3
from scipy import sparse
import numpy as np
from .splice_isoform_survey import dump,write,sha,split_feature


def rows(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:
        yield from csv.DictReader(f,delimiter='\t')


def scan_usage(directory,states,high,q):
    pattern=re.compile(r'^healthy-(.+)-Others-\1_stats-annotated\.txt$')
    files={}
    for p in Path(directory).glob('*_stats-annotated.txt'):
        m=pattern.match(p.name)
        if not m:raise ValueError('Unexpected usage filename: '+p.name)
        files[m[1]]=p
    if set(files)!=set(states):raise ValueError('Usage source and complete healthy state roster disagree')
    stats=defaultdict(dict);symbols={};sources=[];unmatched=set()
    for state in states:
        p=files[state];n=tested=sig=0;seen=set()
        with open(p) as f:
            reader=csv.reader(f,delimiter='\t');header=next(reader);col={k:i for i,k in enumerate(header)}
            fwd,rev='DiffMeans_'+state+'_m_Others','DiffMeans_Others_m_'+state
            if sum(c in col for c in (fwd,rev))!=1:raise ValueError('Usage test orientation unresolved')
            effect=fwd if fwd in col else rev
            for r in reader:
                n+=1;feature=r[col['Feature']]
                if feature in seen:raise ValueError('Duplicate usage feature: '+feature)
                seen.add(feature);gene,tid=split_feature(feature);symbol=r[col['Gene_Symbol']]
                if gene in symbols and symbols[gene]!=symbol:raise ValueError('Conflicting source gene symbol')
                symbols[gene]=symbol
                if feature not in high:continue
                try:adj=float(r[col['ebayesAdjPval']]);raw=float(r[col[effect]]);ms=float(r[col['Mean_'+state]]);mo=float(r[col['Mean_Others']])
                except ValueError:continue
                if not all(map(math.isfinite,[adj,raw,ms,mo])):continue
                delta=ms-mo
                if not math.isclose(raw,delta if effect==fwd else -delta,rel_tol=1e-4,abs_tol=1e-5):raise ValueError('Usage effect disagrees with source means')
                tested+=1;sig+=adj<q
                stats[feature][state]=dict(dusage=delta,q=adj,mean_usage=ms,mean_usage_others=mo,
                    n_state=int(r[col['N_'+state]]),n_others=int(r[col['N_Others']]))
        sources.append(dict(state=state,path=str(p),sha256=sha(p),rows=n,high_read_tested=tested,high_read_significant=sig))
        print('[usage]',state,n,tested,sig,flush=True)
    pairs=defaultdict(dict);by=defaultdict(list)
    for feature,tests in stats.items():by[split_feature(feature)[0]].append(feature)
    for gene,features in by.items():
        for a,b in itertools.combinations(sorted(features),2):
            shared=set(stats[a])&set(stats[b]);both={}
            for state in shared:
                x,y=stats[a][state],stats[b][state]
                if x['q']<q and y['q']<q and x['dusage']*y['dusage']<0:both[state]=dict(isoform_1=x,isoform_2=y)
            if both:pairs[(a,b)]=both
    return stats,symbols,sources,pairs


def extract_registered(run,wanted,out):
    base=Path(run)/'isv/gff-output';summary={};coding={};sequences={};inputs={}
    for name,target in [('protein_summary.txt',summary),('coding_regions.txt',coding)]:
        p=base/name
        for r in rows(p):
            gene=r['Gene ID'].split(';')[0];feature=gene+':'+r['Transcript ID']
            if feature not in wanted:continue
            if feature in target and target[feature]!=r:raise ValueError('Conflicting registered '+name+' record: '+feature)
            target[feature]=r
        inputs[str(p)]=sha(p)
    p=base/'protein_sequences.fasta.gz';current=None;seq=[]
    def save():
        if current in wanted:
            s=''.join(seq)
            if current in sequences and sequences[current]!=s:raise ValueError('Conflicting literal peptide: '+current)
            sequences[current]=s
    with gzip.open(p,'rt') as f:
        for line in f:
            if line.startswith('>'):
                save();tid=line[1:].split(';')[0].strip();m=re.search(r'gene_id:([^;\s]+)',line)
                current=(m[1]+':'+tid) if m else None;seq=[]
            elif current in wanted:seq.append(line.strip())
        save()
    inputs[str(p)]=sha(p)
    write(out/'registered_protein_metadata.tsv',([dict(feature=k,**v) for k,v in sorted(summary.items())]),fields=['feature','Gene ID','Transcript ID','Protein Length','NMD Status','Intron Retention','Longest Isoform Length'])
    dump(out/'registered_coding_regions.json',coding)
    with open(out/'native_predicted_proteins.fasta','w') as f:
        for k,s in sorted(sequences.items()):f.write('>'+k+'\n'+s+'\n')
    audit=[]
    for feature in sorted(wanted):
        r=summary.get(feature);s=sequences.get(feature);c=coding.get(feature)
        audit.append(dict(feature=feature,summary_present=r is not None,peptide_present=s is not None,coding_present=c is not None,
            metadata_aa='' if r is None else r['Protein Length'],literal_peptide_aa='' if s is None else len(s),
            length_agrees='' if r is None or s is None else int(r['Protein Length'])==len(s),
            nmd='' if r is None else r['NMD Status']))
    write(out/'native_protein_source_audit.tsv',audit)
    return inputs


def run(args):
    survey=Path(args.survey);manifest=json.load(open(survey/'manifest.json'))
    if manifest['status']!='complete':raise ValueError('Primary survey is incomplete')
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):raise ValueError('Use a fresh review output directory')
    dump(out/'manifest.json',dict(status='incomplete'))
    features=json.load(gzip.open(survey/'count_features.json.gz','rt'));states=json.load(open(survey/'states.json'))
    counts=sparse.load_npz(survey/'healthy_state_counts.npz').tocsc();maxima=np.asarray(counts.max(axis=0).toarray()).ravel()
    high={f for f,v in zip(features,maxima) if v>manifest['method']['read_threshold_strict_gt']};index={f:i for i,f in enumerate(features)}
    stats,symbols,sources,pairs=scan_usage(args.usage,states,high,manifest['method']['q_threshold'])
    literature=json.load(open(args.literature));catalog={r['gene_symbol']:r for r in literature['genes']}
    selected={f for f in high if symbols.get(split_feature(f)[0]) in catalog}
    event_pairs=list(rows(survey/'native_isoform_pair_catalog.tsv.gz'))
    pairrows=[]
    for (a,b),tests in sorted(pairs.items()):
        gene=split_feature(a)[0];i,j=index[a],index[b];x=counts[:,i].toarray().ravel();y=counts[:,j].toarray().ravel()
        pairrows.append(dict(gene=gene,gene_symbol=symbols[gene],isoform_1=a,isoform_2=b,max_reads_1=int(maxima[i]),max_reads_2=int(maxima[j]),
            max_read_state_1=states[int(x.argmax())],max_read_state_2=states[int(y.argmax())],n_reciprocal_significant_states=len(tests),
            reciprocal_tests=json.dumps(tests,separators=(',',':'))))
    write(out/'reciprocal_usage_pair_catalog.tsv.gz',pairrows,fields=['gene','gene_symbol','isoform_1','isoform_2','max_reads_1','max_reads_2','max_read_state_1','max_read_state_2','n_reciprocal_significant_states','reciprocal_tests'])
    write(out/'high_read_usage_tests.tsv.gz',(dict(feature=f,state=s,**r) for f,ss in sorted(stats.items()) for s,r in sorted(ss.items())),fields=['feature','state','dusage','q','mean_usage','mean_usage_others','n_state','n_others'])
    allpairs={tuple(sorted([r['isoform_1'],r['isoform_2']])):dict(event_pair=True,usage_pair=False) for r in event_pairs}
    for pair in pairs:allpairs.setdefault(pair,dict(event_pair=False,usage_pair=False))['usage_pair']=True
    write(out/'literature_candidate_pairs.tsv',(dict(gene_symbol=symbols[split_feature(a)[0]],isoform_1=a,isoform_2=b,**v) for (a,b),v in sorted(allpairs.items()) if symbols[split_feature(a)[0]] in catalog),fields=['gene_symbol','isoform_1','isoform_2','event_pair','usage_pair'])
    write(out/'literature_screen_coverage.tsv',(dict(gene_symbol=g,genes=';'.join(k for k,v in symbols.items() if v==g),
        high_read_native_isoforms=sum(symbols.get(split_feature(f)[0])==g for f in high),
        junction_linked_pairs=sum(r['gene_symbol']==g for r in event_pairs),reciprocal_usage_pairs=sum(r['gene_symbol']==g for r in pairrows),
        published_mechanism=catalog[g]['published_mechanism'],primary_sources=';'.join(catalog[g]['primary_sources'])) for g in sorted(catalog)))
    con=sqlite3.connect(Path(args.transcript_index).resolve().as_uri()+'?mode=ro',uri=True);structures={}
    for gene in sorted({split_feature(f)[0] for f in selected}):
        for _,strand,tokens,tid,sample in con.execute('SELECT * FROM transcripts WHERE gene=?',(gene,)):
            feature=gene+':'+tid
            if feature in selected:
                if feature in structures and structures[feature]!=tokens:raise ValueError('Conflicting native structure')
                structures[feature]=tokens
    con.close()
    missing=selected-set(structures)
    write(out/'unresolved_selected_models.tsv',(dict(feature=f) for f in sorted(missing)),fields=['feature'])
    if missing:raise ValueError('Exact native models unresolved; no substitutions permitted')
    dump(out/'literature_native_token_structures.json',structures)
    write(out/'literature_native_read_support.tsv',(dict(gene_symbol=symbols[split_feature(f)[0]],feature=f,max_reads=int(maxima[index[f]]),
        max_read_state=states[int(counts[:,index[f]].toarray().argmax())],
        reads_by_state=json.dumps(dict(zip(states,map(int,counts[:,index[f]].toarray().ravel()))),separators=(',',':'))) for f in sorted(selected)))
    inputs=extract_registered(args.run,selected,out)
    dump(out/'manifest.json',dict(status='complete',code_sha256=sha(__file__),survey_manifest_sha256=sha(survey/'manifest.json'),
        literature_sha256=sha(args.literature),usage_sources=sources,inputs=inputs,selected_native_isoforms=len(selected),
        reciprocal_usage_pairs=len(pairrows),reciprocal_usage_genes=len({r['gene'] for r in pairrows}),
        method='Original per-state usage FDR below primary threshold in both opposite-direction isoforms; no new tests or effect filter; both independently pass strict pooled-read threshold',
        limitations='Event-linked and reciprocal-usage pairs are separate. Promoter identity and published mechanism require event/sequence verification. Literal predicted peptides do not demonstrate translation.',
        outputs={p.name:sha(p) for p in out.iterdir() if p.name!='manifest.json'}))


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['survey','run','usage','transcript-index','literature','outdir']:p.add_argument('--'+name,required=True)
    run(p.parse_args(argv))

if __name__=='__main__':main()
