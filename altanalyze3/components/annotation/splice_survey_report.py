"""Validate curated published-mechanism matches and export a healthy splice survey.

Curation is an explicit, hashed input. Native sequence expectations fail closed;
all other qualifying pairs and every screened gene remain in separate full tables.
This does not generate new experimental evidence or new significance tests.
"""
import argparse
from collections import defaultdict
import gzip
from html import escape
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from Bio import SeqIO
from .splice_isoform_survey import dump,write,sha,split_feature
from .splice_survey_review import rows


def check_expectations(observed,expected,label):
    for key,value in expected.items():
        if observed.get(key)!=value:raise ValueError('Native literature-match expectation failed: '+label+' '+key+' '+repr(observed.get(key))+' != '+repr(value))


def shared_read_support(a,b,threshold):
    if a.keys()!=b.keys():raise ValueError('Native isoforms have different state rosters')
    shared=[s for s in sorted(a) if min(a[s],b[s])>threshold]
    best=max(shared,key=lambda s:min(a[s],b[s])) if shared else None
    return dict(n_states_both_gt_threshold=len(shared),same_state_read_supported=bool(shared),
        shared_high_read_state='' if best is None else best,
        reads_1_in_shared_state='' if best is None else a[best],reads_2_in_shared_state='' if best is None else b[best])


def html_table(records,fields):
    return '<table><thead><tr>'+''.join('<th>'+escape(k.replace('_',' '))+'</th>' for k in fields)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+escape(str(r.get(k,'')))+'</td>' for k in fields)+'</tr>' for r in records)+'</tbody></table>'


def validate_outputs(base,manifest):
    if manifest['status']!='complete':raise ValueError('Incomplete source result: '+str(base))
    for name,digest in manifest['outputs'].items():
        if sha(base/name)!=digest:raise ValueError('Changed source output: '+str(base/name))


def resolve_screen_symbols(screen,path):
    """Resolve missing symbols from a pinned existing annotation, never a zero fill."""
    missing=[r for r in screen if not r['genes']]
    found=defaultdict(set)
    if missing:
        if not path:raise ValueError('Missing literature gene mapping; supply the existing --gene-symbols annotation')
        wanted={r['gene_symbol'] for r in missing}
        with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:
            for line in f:
                p=line.rstrip('\n').split('\t')
                if len(p)>=2 and p[0].startswith('ENSG') and p[1] in wanted:found[p[1]].add(p[0])
    for r in screen:
        r['gene_mapping_source']='existing usage annotation' if r['genes'] else 'pinned bundled gene-symbol annotation; gene absent from usage source'
        if not r['genes']:
            if not found[r['gene_symbol']]:raise ValueError('Literature gene remains unresolved: '+r['gene_symbol'])
            r['genes']=';'.join(sorted(found[r['gene_symbol']]))
    return screen


def run(args):
    survey,review,correspondence=map(Path,[args.survey,args.review,args.correspondence]);out=Path(args.outdir)
    manifests=[json.load(open(b/'manifest.json')) for b in [survey,review,correspondence]]
    for base,m in zip([survey,review,correspondence],manifests):validate_outputs(base,m)
    sm,rm,cm=manifests
    if rm['survey_manifest_sha256']!=sha(survey/'manifest.json') or cm['review_manifest_sha256']!=sha(review/'manifest.json'):raise ValueError('Result provenance chain disagrees')
    curated=json.load(open(args.curation));out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):raise ValueError('Use a fresh report output directory')
    dump(out/'manifest.json',dict(status='incomplete'))
    differences={tuple(sorted([r['isoform_1'],r['isoform_2']])):r for r in rows(correspondence/'native_pair_differences.tsv')}
    support={r['feature']:r for r in rows(review/'literature_native_read_support.tsv')}
    screen=resolve_screen_symbols(list(rows(review/'literature_screen_coverage.tsv')),args.gene_symbols)
    sequences={r.id:str(r.seq) for r in SeqIO.parse(review/'native_predicted_proteins.fasta','fasta')}
    pair_events=defaultdict(list)
    for r in rows(correspondence/'literature_pair_event_tests.tsv'):pair_events[tuple(sorted([r['isoform_1'],r['isoform_2']]))].append(r)
    by_short={r['feature'].split(':',1)[1]:r['feature'] for r in support.values()}
    if len(by_short)!=len(support):raise ValueError('Ambiguous native isoform identity across genes')
    checked=[];case_by_pair={};detail=[]
    for c in curated['pairs']:
        a,b=by_short[c['isoform_id_1']],by_short[c['isoform_id_2']];key=tuple(sorted([a,b]));r=differences.get(key)
        if r is None:raise ValueError('Curated pair is absent from measured candidate evidence: '+str(key))
        if r['gene_symbol']!=c['gene_symbol']:raise ValueError('Curated gene identity mismatch')
        d=json.loads(r['protein_comparison'])
        if (a,b)!=(r['isoform_1'],r['isoform_2']):
            # Expectations are always in the explicit curation order.
            from .splice_survey_correspondence import compare_peptides
            d=compare_peptides(sequences[a],sequences[b])
        check_expectations(d,c['expected_protein'],c['gene_symbol'])
        sa,sb=support[a],support[b]
        if min(int(sa['max_reads']),int(sb['max_reads']))<=sm['method']['read_threshold_strict_gt']:raise ValueError('Curated pair fails strict read criterion')
        if c['gene_symbol']=='TCF4':
            if not sequences[a].startswith('MHHQQRMAALGTDKELSDLLDFSAMFSPPV') or not sequences[b].startswith('MYCAYTIPGMGGNSLMYYYNGKA'):raise ValueError('TCF4 published Figure S1 N-terminus check failed')
        tests=pair_events[key];strong=max(tests,key=lambda x:abs(float(x['dpsi']))) if tests else None
        row=dict(gene_symbol=c['gene_symbol'],published_mechanism=c['published_mechanism'],status=c['status'],
            isoform_1=a,isoform_2=b,predicted_aa_1=d['aa_1'],predicted_aa_2=d['aa_2'],
            max_reads_1=sa['max_reads'],max_read_state_1=sa['max_read_state'],max_reads_2=sb['max_reads'],max_read_state_2=sb['max_read_state'],
            strongest_event_state='' if strong is None else strong['state'],max_abs_source_dpsi='' if strong is None else abs(float(strong['dpsi'])),
            strongest_event_q='' if strong is None else strong['q'],n_donors_in_event_state='' if strong is None else strong['n_state'],
            strongest_event='' if strong is None else strong['event'],
            event_evidence=r['event_pair'],reciprocal_usage_evidence=r['usage_pair'],donor_pair_evidence=r['donor_pair'],
            interpretation=c['interpretation'],primary_sources=';'.join(c['primary_sources']),
            reciprocal_usage_tests=r['reciprocal_usage_tests'],donor_pair_tests=r['donor_pair_tests'])
        row.update(shared_read_support(json.loads(sa['reads_by_state']),json.loads(sb['reads_by_state']),sm['method']['read_threshold_strict_gt']))
        checked.append(row);case_by_pair[key]=row
        detail.append(dict(c,isoform_1=a,isoform_2=b,observed_protein=d,read_support_1=sa,read_support_2=sb,event_tests=tests,
            reciprocal_usage_tests=json.loads(r['reciprocal_usage_tests'] or '{}'),donor_pair_tests=json.loads(r['donor_pair_tests'] or '[]')))
    write(out/'verified_literature_correspondence.tsv',checked);dump(out/'verified_literature_correspondence.json',detail)
    assessments=[]
    for key,r in sorted(differences.items()):
        hit=case_by_pair.get(key)
        assessments.append(dict(r,curated_status='Not established for this native pair' if hit is None else hit['status'],
            curated_mechanism='' if hit is None else hit['published_mechanism'],
            primary_sources='' if hit is None else hit['primary_sources']))
    write(out/'all_literature_pair_assessments.tsv',assessments)
    features=json.load(gzip.open(survey/'count_features.json.gz','rt'));states=json.load(open(survey/'states.json'))
    matrix=sparse.load_npz(survey/'healthy_state_counts.npz').tocsc();maximum=np.asarray(matrix.max(axis=0).toarray()).ravel()
    if len(features)!=sm['input_feature_roster'] or int((maximum>sm['method']['read_threshold_strict_gt']).sum())!=sm['native_high_read_features']:raise ValueError('Complete native count roster/read-screen audit failed')
    gene_to_symbol={g:r['gene_symbol'] for r in screen for g in r['genes'].split(';') if g}
    all_native=[];gene_native=defaultdict(list)
    for i,f in enumerate(features):
        gene,tid=split_feature(f)
        if gene not in gene_to_symbol:continue
        value=int(maximum[i]);record=dict(gene_symbol=gene_to_symbol[gene],feature=f,max_reads=value,passes_strict_read_rule=value>sm['method']['read_threshold_strict_gt'])
        all_native.append(record);gene_native[record['gene_symbol']].append(record)
    write(out/'all_native_literature_gene_read_screen.tsv.gz',all_native)
    overview=[]
    for r in screen:
        g=r['gene_symbol'];rr=sorted(gene_native[g],key=lambda x:(-x['max_reads'],x['feature']));positive=[c for c in checked if c['gene_symbol']==g]
        actual_high=sum(x['passes_strict_read_rule'] for x in rr)
        if r['gene_mapping_source']=='existing usage annotation' and actual_high!=int(r['high_read_native_isoforms']):raise ValueError('Literature gene high-read roster disagrees: '+g)
        r['high_read_native_isoforms']=actual_high
        reason=('Resolved gene has no native feature in the complete count matrix; no expression value inferred' if not rr else 'Curated published correspondence' if positive else 'Fewer than two native isoforms independently exceed read threshold' if int(r['high_read_native_isoforms'])<2 else 'No qualifying differential pair in these existing test catalogs' if int(r['junction_linked_pairs'])+int(r['reciprocal_usage_pairs'])==0 else 'Differential native pair(s) present; published mechanism not established')
        overview.append(dict(r,registered_count_features=len(rr),curated_correspondence_pairs=len(positive),screen_conclusion=reason,
            highest_read_feature='' if not rr else rr[0]['feature'],highest_read_max='' if not rr else rr[0]['max_reads'],
            second_highest_read_feature='' if len(rr)<2 else rr[1]['feature'],second_highest_read_max='' if len(rr)<2 else rr[1]['max_reads']))
    write(out/'literature_gene_screen.tsv',overview)
    # Complete primary gene summary: no literature-based reduction of global results.
    global_pairs=list(rows(survey/'native_isoform_pair_catalog.tsv.gz'))
    needed={r[k] for r in global_pairs for k in ['isoform_1','isoform_2']}
    positions={f:i for i,f in enumerate(features) if f in needed}
    if set(positions)!=needed:raise ValueError('Incomplete global native pair roster')
    values={f:dict(zip(states,map(int,matrix[:,i].toarray().ravel()))) for f,i in positions.items()}
    global_genes=defaultdict(lambda:dict(events=set(),pairs=set(),states=set(),same_state_pairs=0))
    pair_read_audit=[]
    for r in global_pairs:
        shared=shared_read_support(values[r['isoform_1']],values[r['isoform_2']],sm['method']['read_threshold_strict_gt'])
        pair_read_audit.append(dict(r,**shared))
        g=(r['gene'],r['gene_symbol']);global_genes[g]['events'].update(r['events'].split(';'));global_genes[g]['pairs'].add((r['isoform_1'],r['isoform_2']));global_genes[g]['states'].update(r['significant_event_states'].split(';'));global_genes[g]['same_state_pairs']+=shared['same_state_read_supported']
    write(out/'global_read_supported_splice_pairs.tsv.gz',pair_read_audit)
    write(out/'global_read_supported_splice_genes.tsv',(dict(gene_id=g[0],gene_symbol=g[1],events=len(v['events']),native_pairs=len(v['pairs']),same_state_high_read_pairs=v['same_state_pairs'],significant_states=';'.join(sorted(v['states'])),literature_screened=g[1] in {r['gene_symbol'] for r in screen}) for g,v in sorted(global_genes.items())))
    same_state_pairs=sum(r['same_state_read_supported'] for r in pair_read_audit);same_state_genes=sum(v['same_state_pairs']>0 for v in global_genes.values())
    method='Healthy group: '+str(len(sm['healthy_donors']))+' donors, '+str(sm['healthy_pseudobulks'])+' donor-state pseudobulks, all '+str(len(states))+' recorded states (including the recorded stromal controls). Both native isoforms independently require >'+str(sm['method']['read_threshold_strict_gt'])+' pooled assigned reads in at least one healthy state; the states may differ. Original state-versus-Others adjusted P<'+str(sm['method']['q_threshold'])+' retained; no effect-size gate or new test. Source effects remain continuous. Pooled dominance is descriptive. These tests do not establish protein translation or function in these cells.'
    content='<!doctype html><html><head><meta charset="utf-8"><title>Healthy hematopoietic isoform survey</title><style>body{font:15px/1.5 Arial,sans-serif;margin:32px;color:#17212b}table{border-collapse:collapse;width:100%;font-size:13px}th,td{border:1px solid #c9d0d7;padding:8px;vertical-align:top}th{background:#edf1f5;text-align:left;position:sticky;top:0}h1{font-size:24px}h2{font-size:19px}a{color:#14589c}section{margin:28px 0}</style></head><body><h1>Healthy hematopoietic splice/promoter survey</h1><p>'+escape(method)+'</p><p>'+str(sm['events_with_read_supported_exclusive_pairs'])+' significant junction events link '+str(sm['native_isoform_pairs'])+' native pairs across '+str(sm['genes_with_read_supported_events'])+' genes. The stricter requirement that both partners exceed 100 reads together in a single state retains '+str(same_state_pairs)+' pairs across '+str(same_state_genes)+' genes; all broader pairs remain flagged in the full table. Existing reciprocal transcript-usage tests identify '+str(rm['reciprocal_usage_pairs'])+' pairs across '+str(rm['reciprocal_usage_genes'])+' genes (these can overlap the junction-linked pairs). The literature screen covers '+str(len(screen))+' genes. This is a genome-wide survey of the existing tests with targeted primary-literature mechanism verification; it is not proof that all published hematopoietic events have been indexed.</p><p>Downloads: <a href="verified_literature_correspondence.tsv">Verified correspondence and statistics</a>; <a href="verified_literature_correspondence.json">Complete per-state evidence</a>; <a href="all_literature_pair_assessments.tsv">Every screened differential pair</a>; <a href="literature_gene_screen.tsv">All screened genes and read failures</a>; <a href="global_read_supported_splice_genes.tsv">All genome-wide qualifying genes</a>.</p><h2>Published mechanisms and native evidence</h2>'+html_table(checked,['gene_symbol','published_mechanism','status','isoform_1','predicted_aa_1','max_reads_1','max_read_state_1','isoform_2','predicted_aa_2','max_reads_2','max_read_state_2','shared_high_read_state','reads_1_in_shared_state','reads_2_in_shared_state','strongest_event_state','max_abs_source_dpsi','strongest_event_q'])
    for r in checked:
        content+='<section><h2>'+escape(r['gene_symbol']+' — '+r['published_mechanism'])+'</h2><p>'+escape(r['interpretation'])+'</p><p>'+', '.join('<a href="'+escape(u,quote=True)+'">Primary source '+str(i+1)+'</a>' for i,u in enumerate(r['primary_sources'].split(';')))+'</p></section>'
    content+='<h2>Complete literature screen</h2>'+html_table(overview,['gene_symbol','high_read_native_isoforms','junction_linked_pairs','reciprocal_usage_pairs','curated_correspondence_pairs','second_highest_read_max','screen_conclusion'])+'</body></html>'
    (out/'healthy_hematopoietic_splice_survey.html').write_text(content)
    (out/'METHODS.txt').write_text(method+'\n\nThe supplementary donor-paired switches come from the existing coding-only, differential-expression-whitelisted, gene-read>200 run. They retain its donor pairing, moderated t statistics, BH correction and |usage change|>=0.2 gate. They supplement the global PSI survey and are not a complete independent test of all native isoforms.\n\nLiterature match statuses are explicit manual curation with executable expected native sequence checks. New molecule IDs are not classified as novel proteins. AltAnalyze exon-group labels are not assumed to equal publication exon numbers.\n\nSource annotation discrepancies remain in native_protein_source_audit.tsv. IGHM is reported only as an RNA-processing match because its metadata and literal peptides disagree.\n')
    dump(out/'manifest.json',dict(status='complete',code_sha256=sha(__file__),source_manifests={str(b):sha(b/'manifest.json') for b in [survey,review,correspondence]},
        gene_symbol_mapping=None if not args.gene_symbols else dict(path=args.gene_symbols,sha256=sha(args.gene_symbols)),
        curation=dict(path=str(Path(args.curation).resolve()),sha256=sha(args.curation)),verified_pairs=len(checked),verified_genes=len({r['gene_symbol'] for r in checked}),literature_screen_genes=len(screen),
        same_state_high_read_pairs=same_state_pairs,same_state_high_read_genes=same_state_genes,
        primary_gene_roster=len(global_genes),native_read_screen_records=len(all_native),outputs={p.name:sha(p) for p in out.iterdir() if p.name!='manifest.json'}))


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['survey','review','correspondence','curation','outdir']:p.add_argument('--'+name,required=True)
    p.add_argument('--gene-symbols',help='Pinned existing Ensembl gene-symbol annotation (plain or gzipped); required if a literature gene is absent from usage annotation')
    run(p.parse_args(argv))

if __name__=='__main__':main()
