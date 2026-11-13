"""Observed, reference-backed candidate sequence hypotheses; no classification.

Do not bridge unknown internal intervals. Keep every distinct compatible observed
chain and explicit unresolved construction outcomes. Full source paths are retained
elsewhere; clipping here defines a candidate-template hypothesis, never read geometry.
"""
from __future__ import annotations
import collections,gzip,hashlib,json,math,os,pathlib,re,subprocess,time
from .circRNA_full_query_validation_v2 import cigar_ops,contig,revcomp,reconcile_source,digest,event_key,parse_paf

def alignment_exons(a):
 t=int(a['start0']);start=t;blocks=[]
 for n,op in cigar_ops(a['cigar']):
  if op=='N':
   if t>start:blocks.append((start,t))
   t+=n;start=t
  elif op in 'MD=X':t+=n
 if t>start:blocks.append((start,t))
 return blocks

def clipped_chain(blocks,lo,hi):return [(max(a,lo),min(b,hi)) for a,b in blocks if max(a,lo)<min(b,hi)]
def union_blocks(blocks):
 out=[]
 for a,b in sorted(blocks):
  if out and a<=out[-1][1]:out[-1]=(out[-1][0],max(out[-1][1],b))
  else:out.append((a,b))
 return out

def compatible_overlap(a,b):
 """Two chains overlap or touch in reference and disagree at no splice gap."""
 if not any(max(x,u)<=min(y,v) for x,y in a for u,v in b):return False
 for first,second in [(a,b),(b,a)]:
  introns=[(x[1],y[0]) for x,y in zip(first,first[1:])]
  if any(max(x,u)<min(y,v) for x,y in introns for u,v in second):return False
 return True

def observed_model_chains(event,contexts):
 lo=min(int(event['donor']),int(event['acceptor']))-1;hi=max(int(event['donor']),int(event['acceptor']));strand=event['donor_strand']
 models={};outcomes=[]
 for r in contexts:
  if r['event']['event_id']!=event_key(event):raise ValueError('Event context mismatch')
  source=r.get('source',{});pairs=r['baseline_pairs'];v=reconcile_source(r['library'],r['query'],pairs,source)
  source_id=[r['library'],r['query']];result=[]
  if v['actual_query_identity']['status']!='verified_actual_query_full_sequence' or v['selected_arm_matching']['status']!='selected_actual_arms_verified':
   outcomes.append(dict(source_identity=source_id,status='unresolved_source_identity_or_selected_arms',issues=v['issues'],identity_axes={k:v[k] for k in ('actual_query_identity','baseline_descriptor_parity','selected_arm_matching')},models=[]));continue
  alignments=source['alignments'];checks={(str(x['pair_index']),x['arm']):x for x in v['selected_arm_checks']}
  for p in pairs:
   arms=[]
   for arm in ('a','b'):
    check=checks[(str(p['pair_index']),arm)];a=alignments[check['matching_actual_record_indices'][0]]
    arms.append(clipped_chain(alignment_exons(a),lo,hi))
   proposals=[]
   for a in arms:
    if a and a[0][0]==lo and a[-1][1]==hi:proposals.append((a,'single_actual_alignment_spans_candidate'))
   if all(arms) and compatible_overlap(*arms):
    combined=union_blocks(arms[0]+arms[1])
    if combined[0][0]==lo and combined[-1][1]==hi:proposals.append((combined,'overlapping_or_touching_splice_compatible_actual_arms'))
   for chain,basis in proposals:
    ordered=chain if strand=='+' else list(reversed(chain));key=tuple(ordered)
    entry=models.setdefault(key,dict(blocks=[dict(chrom=event['chrom'],start0=a,end0=b,strand=strand) for a,b in ordered],sources=[],construction_bases=set()))
    entry['sources'].append(dict(library=r['library'],query=r['query'],pair_index=p['pair_index'],baseline_status=p['status']));entry['construction_bases'].add(basis);result.append(key)
  outcomes.append(dict(source_identity=source_id,status='candidate_hypothesis_constructed' if result else 'unresolved_internal_chain_or_splice_conflict',identity_axes={k:v[k] for k in ('actual_query_identity','baseline_descriptor_parity','selected_arm_matching')},models=[list(x) for x in sorted(set(result))]))
 for m in models.values():m['construction_bases']=sorted(m['construction_bases'])
 return list(models.values()),outcomes

def block_sequence(blocks,fetch):
 return ''.join(fetch(b['chrom'],b['start0'],b['end0']).upper() if b['strand']=='+' else revcomp(fetch(b['chrom'],b['start0'],b['end0']).upper()) for b in blocks)

def build_event_models(event,contexts,fetch,contig_length,reference_sha256,external_chains=None):
 """Return all observed hypotheses plus explicitly sourced external chains.

 external_chains may contain curated exon chains or actual native full linear
 comparator structures, each with its own model_basis and source_pin. None is
 chosen as truth or silently substituted for an unresolved observed structure.
 """
 chains,outcomes=observed_model_chains(event,contexts);unresolved=[];templates=[]
 if len(reference_sha256)!=64:raise ValueError('Missing reference byte SHA')
 for m in external_chains or []:
  if not m.get('source_pin') or not m.get('model_basis'):raise ValueError('Unproven external model')
  if m.get('role','candidate')=='candidate':chains.append(dict(blocks=m['coordinate_blocks'],sources=[],construction_bases=[m['model_basis']],external_source_pin=m['source_pin']))
 maxread=max([len(r.get('source',{}).get('full_forward_sequence') or '') for r in contexts] or [0])
 for i,m in enumerate(chains):
  blocks=m['blocks'];seq=block_sequence(blocks,fetch)
  if not seq:raise ValueError('Empty candidate template')
  first=blocks[0]['start0']+1 if blocks[0]['strand']=='+' else blocks[0]['end0'];last=blocks[-1]['end0'] if blocks[-1]['strand']=='+' else blocks[-1]['start0']+1
  if (first,last)!=(int(event['acceptor']),int(event['donor'])):raise ValueError('External candidate model endpoints differ')
  copies=max(2,math.ceil(maxread/len(seq))+2);length=copies*len(seq);base='candidate_%04d'%i
  common=dict(reference_sha256=reference_sha256,coordinate_blocks=blocks,copy_length=len(seq),copy_number=copies,sequence_length=length,model_basis=m,opportunity_caveat='Equal template lengths do not equalize repeated target opportunities. Direct tandem is shared by circular-copy and linear-concatemer molecules.')
  direct=seq*copies
  templates.append(dict(model_id=base+'_direct',kind='exon_chain' if len(blocks)>1 else 'genomic_interval',sequence=direct,sequence_sha256=hashlib.sha256(direct.encode()).hexdigest(),boundaries=list(range(len(seq),length,len(seq))),**common))
  for label,pattern in [('alternating_orientation',seq+revcomp(seq)),('forward_forward_reverse_reverse',seq*2+revcomp(seq)*2)]:
   alt=(pattern*math.ceil(length/len(pattern)))[:length]
   templates.append(dict(model_id=base+'_'+label,kind='alternative_sequence',sequence=alt,sequence_sha256=hashlib.sha256(alt.encode()).hexdigest(),boundaries=[],**common))
  chrom=event['chrom'];n=contig_length(chrom);center=(min(int(event['donor']),int(event['acceptor']))+max(int(event['donor']),int(event['acceptor']))-1)//2
  if length>n:
   unresolved.append(dict(model_id=base+'_linear_genomic',reason='equal_length_context_exceeds_reference_contig',requested_length=length,available=n));continue
  start=max(0,min(n-length,center-length//2));linear=fetch(chrom,start,start+length).upper()
  if len(linear)!=length:raise ValueError('Reference fetch shorter than requested')
  templates.append(dict(model_id=base+'_linear_genomic',kind='alternative_genomic',sequence=linear,sequence_sha256=hashlib.sha256(linear.encode()).hexdigest(),reference_sha256=reference_sha256,coordinate_blocks=[dict(chrom=chrom,start0=start,end0=start+length,strand='+')],sequence_length=length,boundaries=[],model_basis={'type':'ordinary_genomic_context','candidate_model':base},opportunity_caveat='Ordinary genomic sequence; splicing assessed separately with whole-genome splice:hq. No circularity claim from tandem advantage.'))
 for i,m in enumerate(external_chains or []):
  if m.get('role') not in ('linear_comparator','alternative_query_path'):continue
  seq=block_sequence(m['coordinate_blocks'],fetch)
  templates.append(dict(model_id=('native_linear_' if m.get('role')=='linear_comparator' else 'reference_query_path_')+'%04d'%i,kind='alternative_sequence',sequence=seq,sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),reference_sha256=reference_sha256,coordinate_blocks=m['coordinate_blocks'],sequence_length=len(seq),boundaries=[],model_basis=m,opportunity_caveat='Actual reference-backed transcript/path length retained; not artificially repeated/padded. Alignment score/coverage and length must be reported jointly.'))
 if not any(m.get('role')=='linear_comparator' for m in external_chains or []):unresolved.append(dict(model_id='native_linear_comparator',reason='no_verified_native_full_transcript_model_supplied'))
 return dict(event=event,reference_sha256=reference_sha256,templates=templates,source_model_outcomes=outcomes,unresolved_models=unresolved,max_query_length=maxread,retention='All contexts and distinct chain hypotheses retained; models are sequence hypotheses, not inferred physical topology.')

def run_template_alignments(bundle,queries,outdir,threads=1,minimap2='minimap2',outcome_storage='inline'):
 """Identical map-hifi settings, each template separately; every query outcome retained.

 Genome remapping is run once for the union of queries by the parent batch runner,
 not once per event. All unmapped/failed/not-evaluated outcomes remain explicit.
 """
 if outcome_storage not in ('inline','external'):raise ValueError('Unknown outcome storage')
 model_ids=[m['model_id'] for m in bundle['templates']]
 if len(model_ids)!=len(set(model_ids)):raise ValueError('Duplicate template model ID')
 if any(not re.fullmatch('[A-Za-z0-9_.-]+',q['id']) for q in queries):raise ValueError('Missing/unsafe FASTA query ID')
 out=pathlib.Path(outdir)
 if out.exists():raise FileExistsError(out)
 out.mkdir(parents=True);qpath=out/'queries.fa';ids=set();rows=[]
 for r in queries:
  if r['id'] in ids:raise ValueError('Duplicate alignment query ID')
  ids.add(r['id']);seq=r.get('sequence')
  if seq:rows.append('>'+r['id']+'\n'+seq+'\n')
 qpath.write_text(''.join(rows));results=[];sequence_cache={};qlength={q['id']:len(q['sequence']) for q in queries if q.get('sequence')};version=subprocess.check_output([minimap2,'--version'],text=True).strip()
 for m in bundle['templates']:
  name=m['model_id']
  if not re.fullmatch('[A-Za-z0-9_.-]+',name):raise ValueError('Unsafe model name')
  fa=out/(name+'.fa');paf=out/(name+'.paf');err=out/(name+'.stderr')
  sequence_sha=hashlib.sha256(m['sequence'].encode()).hexdigest()
  if sequence_sha!=m['sequence_sha256']:raise ValueError('Model sequence SHA differs')
  if sequence_sha in sequence_cache:
   original=sequence_cache[sequence_sha]
   for suffix in ('.fa','.paf','.stderr'):os.link(out/(original['executed_model_id']+suffix),out/(name+suffix))
   results.append(dict(original,model_id=name,reused_exact_sequence=True));continue
  fa.write_text('>'+name+'\n'+m['sequence']+'\n')
  cmd=[minimap2,'-x','map-hifi','-c','--cs=long','--secondary=yes','-N','50','-p','0','-t',str(threads),str(fa),str(qpath)];start=time.time()
  with paf.open('w') as stdout,err.open('w') as stderr:r=subprocess.run(cmd,stdout=stdout,stderr=stderr)
  counts=collections.Counter();paf_hashes=collections.defaultdict(hashlib.sha256)
  if r.returncode==0:
   with paf.open() as handle:
    for line in handle:
     parsed=parse_paf(line);q=parsed['query']
     if q not in qlength or parsed['query_length']!=qlength[q]:raise ValueError('Unknown aligner query identity/length')
     counts[q]+=1;paf_hashes[q].update((json.dumps(parsed,sort_keys=True,separators=(',',':'))+'\n').encode())
  results.append(dict(model_id=name,executed_model_id=name,reused_exact_sequence=False,status='completed' if r.returncode==0 else 'failed',returncode=r.returncode,command=cmd,seconds=time.time()-start,template_sha256=digest(fa),paf_sha256=digest(paf),query_fasta_sha256=digest(qpath),model_sequence_sha256=m['sequence_sha256'],outcomes=[dict(id=q['id'],paf_records_sha256=paf_hashes[q['id']].hexdigest(),mapping_count=counts[q['id']] if r.returncode==0 and q.get('sequence') else None,status=('mapped' if counts[q['id']] else 'unmapped_retained') if r.returncode==0 and q.get('sequence') else 'unresolved_missing_sequence' if not q.get('sequence') else 'failed_retained') for q in queries]))
  if outcome_storage=='external':
   path=out/(name+'.outcomes.jsonl.gz')
   with gzip.open(path,'wt') as handle:
    for row in results[-1].pop('outcomes'):handle.write(json.dumps(row,separators=(',',':'))+'\n')
   results[-1]['outcomes_file']=path.name;results[-1]['outcomes_sha256']=digest(path);results[-1]['outcomes_count']=len(queries)
  if r.returncode==0:sequence_cache[sequence_sha]=results[-1]
 report=dict(status='DIAGNOSTIC_EXECUTION_ALL_OUTCOMES_RETAINED',minimap2_version=version,reference_sha256=bundle.get('reference_sha256'),query_sequence_sha256={q['id']:hashlib.sha256(q['sequence'].encode()).hexdigest() for q in queries if q.get('sequence')},models=results,limits=['Heuristic aligner with at most50 secondary chains per primary; not exhaustive alternatives.','map-hifi scores comparable only with identical scoring settings; splice:hq evaluated separately.','Tandem preference cannot distinguish circular-copy from linear concatemer.'])
 (out/'ALIGNMENT_EXECUTION.json').write_text(json.dumps(report,indent=2)+'\n');(out/'MODEL_BUNDLE.json').write_text(json.dumps(bundle,indent=2)+'\n');return report
