"""Audit coding and registered-transcript differences for literature-screen pairs.

Uses AltAnalyze3's existing global protein alignment. Reports literal sequence gaps,
substitutions and transcript tokens without declaring published isoform identity.
All incomplete peptides remain explicit. No new sequence predictions or databases.
"""
import argparse
from collections import defaultdict
import json
from pathlib import Path
from Bio import SeqIO
from ..iso2function.protein_isoform_audit import align_mapping
from .splice_isoform_survey import dump,write,sha
from .splice_survey_review import rows


def ranges(indices):
    out=[]
    for i in sorted(indices):
        if out and i==out[-1][1]+1:out[-1][1]=i
        else:out.append([i,i])
    return [[a+1,b+1] for a,b in out]


def compare_peptides(a,b):
    mapping=align_mapping(a,b)
    prefix=0;suffix=0
    for x,y in zip(a,b):
        if x!=y:break
        prefix+=1
    for x,y in zip(a[::-1],b[::-1]):
        if x!=y:break
        suffix+=1
    deleted=set(range(len(a)))-set(mapping);inserted=set(range(len(b)))-set(mapping.values())
    return dict(aa_1=len(a),aa_2=len(b),common_prefix_aa=prefix,common_suffix_aa=suffix,identical_aligned_residues=sum(a[i]==b[j] for i,j in mapping.items()),
        substitutions=sum(a[i]!=b[j] for i,j in mapping.items()),
        unmapped_aa_1=ranges(deleted),unmapped_aa_2=ranges(inserted),
        missing_sequence_1=[dict(start=x,end=y,sequence=a[x-1:y]) for x,y in ranges(deleted)],
        missing_sequence_2=[dict(start=x,end=y,sequence=b[x-1:y]) for x,y in ranges(inserted)])


def run(args):
    review=Path(args.review);survey=Path(args.survey);out=Path(args.outdir)
    rm=json.load(open(review/'manifest.json'));sm=json.load(open(survey/'manifest.json'))
    if rm['status']!='complete' or sm['status']!='complete':raise ValueError('Incomplete source survey')
    if rm['survey_manifest_sha256']!=sha(survey/'manifest.json'):raise ValueError('Changed primary survey manifest')
    for name in ['native_predicted_proteins.fasta','literature_native_token_structures.json','literature_candidate_pairs.tsv']:
        if sha(review/name)!=rm['outputs'][name]:raise ValueError('Changed review source: '+name)
    out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):raise ValueError('Use a fresh correspondence directory')
    seq={r.id:str(r.seq) for r in SeqIO.parse(review/'native_predicted_proteins.fasta','fasta')}
    structures=json.load(open(review/'literature_native_token_structures.json'))
    events={r['feature']:r for r in rows(survey/'differential_event_catalog.tsv.gz')}
    ep={tuple(sorted([r['isoform_1'],r['isoform_2']])):r for r in rows(survey/'native_isoform_pair_catalog.tsv.gz')}
    up={tuple(sorted([r['isoform_1'],r['isoform_2']])):r for r in rows(review/'reciprocal_usage_pair_catalog.tsv.gz')}
    candidates={tuple(sorted([r['isoform_1'],r['isoform_2']])):dict(r,donor_pair=False) for r in rows(review/'literature_candidate_pairs.tsv')}
    donor_tests=defaultdict(list);donor_audit=[]
    support={r['feature']:r for r in rows(review/'literature_native_read_support.tsv')}
    if args.donor_switches:
        for r in rows(args.donor_switches):
            a,b=r['isoform_higher_in_A'],r['isoform_higher_in_B'];key=tuple(sorted([a,b]))
            # The complete source table remains audited; this join adds only
            # native models already independently meeting the read criterion.
            matched=a in support and b in support
            donor_audit.append(dict(r,literature_pair_with_both_read_supported=matched))
            if matched:
                candidates.setdefault(key,dict(gene_symbol=r['gene_symbol'],isoform_1=key[0],isoform_2=key[1],event_pair=False,usage_pair=False,donor_pair=False))['donor_pair']=True
                donor_tests[key].append(r)
        write(out/'donor_switch_source_audit.tsv.gz',donor_audit)
    diffs=[];links=[]
    for key,pair in sorted(candidates.items()):
        a,b=pair['isoform_1'],pair['isoform_2'];key=tuple(sorted([a,b]));e=ep.get(key);u=up.get(key)
        ta,tb=structures[a].split('|'),structures[b].split('|')
        difference=compare_peptides(seq[a],seq[b]) if a in seq and b in seq else None
        row=dict(pair,peptide_1_present=a in seq,peptide_2_present=b in seq,
            token_chain_equal=ta==tb,first_token_1=ta[0],first_token_2=tb[0],last_token_1=ta[-1],last_token_2=tb[-1],
            unique_tokens_1=json.dumps([t for t in ta if t not in tb]),unique_tokens_2=json.dumps([t for t in tb if t not in ta]),
            protein_comparison=json.dumps(difference,separators=(',',':')),
            reciprocal_usage_tests='' if u is None else u['reciprocal_tests'],
            donor_pair_tests=json.dumps(donor_tests.get(key,[]),separators=(',',':')),
            event_features='' if e is None else e['events'])
        diffs.append(row)
        if e:
            for name in e['events'].split(';'):
                event=events[name]
                for state,test in json.loads(event['state_tests']).items():
                    links.append(dict(gene_symbol=pair['gene_symbol'],isoform_1=a,isoform_2=b,event=name,state=state,
                        event_type=event['event_type'],coord_1=event['coord_1'],coord_2=event['coord_2'],**test))
    write(out/'native_pair_differences.tsv',diffs)
    write(out/'literature_pair_event_tests.tsv',links)
    dump(out/'manifest.json',dict(status='complete',code_sha256=sha(__file__),
        review_manifest_sha256=sha(review/'manifest.json'),survey_manifest_sha256=sha(survey/'manifest.json'),
        method='Existing iso2function.protein_isoform_audit.align_mapping global BLOSUM62 alignment, gap open -10 and extension -0.5; 1-based inclusive literal sequence ranges; exact registered ISV transcript tokens',
        disclaimer='Coding/transcript differences support manual literature-mechanism verification. Neither gene name, length nor alignment alone establishes a published switch or protein function.',
        donor_switch_source=None if not args.donor_switches else dict(path=args.donor_switches,sha256=sha(args.donor_switches),rows=len(donor_audit)),
        native_pairs=len(diffs),pairs_with_both_peptides=sum(r['peptide_1_present'] and r['peptide_2_present'] for r in diffs),
        outputs={p.name:sha(p) for p in out.iterdir()}))


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['review','survey','outdir']:p.add_argument('--'+name,required=True)
    p.add_argument('--donor-switches',help='Optional existing donor-paired isoform_state_switches.tsv; supplemental evidence retains original selection limits')
    run(p.parse_args(argv))

if __name__=='__main__':main()
