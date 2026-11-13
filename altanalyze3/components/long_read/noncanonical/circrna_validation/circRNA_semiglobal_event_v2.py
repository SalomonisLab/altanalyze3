"""Add exact, equally parameterized full-query evidence for every prepared model.

Minimap artifacts are only input/provenance and remain immutable. Distances are
not compared across backends. Every alias/model and every source query is retained.
"""
import argparse,collections,functools,gzip,hashlib,json,pathlib,time
from .circRNA_full_query_validation_v2 import digest,event_key
from .circRNA_semiglobal_v2 import align_model

def fasta(path):
 out={};name=None
 with open(path) as f:
  for line in f:
   if line.startswith('>'):
    name=line[1:].strip()
    if name in out:raise ValueError('Duplicate FASTA identity')
    out[name]=''
   else:out[name]+=line.strip()
 return out

def run(model_dir,expected_stage_sha,source_shard,expected_source_sha,output,admission_manifest,expected_admission_sha):
 d=pathlib.Path(model_dir);out=pathlib.Path(output);source=pathlib.Path(source_shard)
 if out.exists():raise FileExistsError(out)
 stagepath=d/'EVENT_STAGE_RESULT.json'
 if digest(stagepath)!=expected_stage_sha or digest(source)!=expected_source_sha:raise ValueError('Execution/source pin mismatch')
 admissionpath=pathlib.Path(admission_manifest)
 if digest(admissionpath)!=expected_admission_sha:raise ValueError('Independent admission manifest pin mismatch')
 admission=json.loads(admissionpath.read_text())
 if admission.get('status')!='INDEPENDENT_EVENT_MODEL_ARTIFACT_ADMISSION' or not admission.get('admission_basis'):raise ValueError('Independent output admission required')
 for rel,h in admission['model_artifacts'].items():
  if digest(d/rel)!=h:raise ValueError('Admitted model artifact changed '+rel)
 required={'EVENT_STAGE_RESULT.json','PREPARED_MODEL_BUNDLE.json','alignments/ALIGNMENT_EXECUTION.json','alignments/MODEL_BUNDLE.json'}
 if not required.issubset(admission['model_artifacts']) or admission['source_shard_sha256']!=expected_source_sha or admission['model_artifacts']['EVENT_STAGE_RESULT.json']!=expected_stage_sha:raise ValueError('Incomplete or inconsistent output admission')
 stage=json.loads(stagepath.read_text())
 if stage['status']!='EVENT_TEMPLATE_EXECUTION_COMPLETED_CHECK_PER_MODEL_OUTCOMES' or stage['input_sha256']!=expected_source_sha:raise ValueError('Prepared source stage incomplete')
 bundlepath=d/'PREPARED_MODEL_BUNDLE.json';bundle=json.loads(bundlepath.read_text());receiptpath=d/'alignments/ALIGNMENT_EXECUTION.json';receipt=json.loads(receiptpath.read_text())
 if json.loads((d/'alignments/MODEL_BUNDLE.json').read_text())!=bundle:raise ValueError('Prepared/executed model bundle differs')
 models=bundle['templates'];mids=[m['model_id'] for m in models]
 if len(mids)!=len(set(mids)) or len(models)!=stage['template_count']:raise ValueError('Model roster mismatch')
 execution={m['model_id']:m for m in receipt['models']}
 if len(execution)!=len(receipt['models']) or set(execution)!=set(mids):raise ValueError('Missing/duplicate minimap model outcome')
 queries=fasta(d/'alignments/queries.fa');expected={};source_keys=[];paircount=0
 with source.open() as f:
  for line in f:
   r=json.loads(line);qid='Q'+hashlib.sha256(json.dumps([r['library'],r['query']],separators=(',',':')).encode()).hexdigest()
   if r.get('alignment_query_id')!=qid or qid in expected or event_key(r['event'])!=stage['event']:raise ValueError('Source query/event identity mismatch')
   expected[qid]=r['source'].get('full_forward_sequence');source_keys.append(dict(id=qid,library=r['library'],query=r['query']));paircount+=len(r['baseline_pairs'])
 if len(expected)!=stage['source_contexts'] or paircount!=stage['raw_pair_rows'] or queries!={q:s for q,s in expected.items() if s}:raise ValueError('Exact source sequences/roster differ')
 if receipt['query_sequence_sha256']!={q:hashlib.sha256(s.encode()).hexdigest() for q,s in expected.items() if s}:raise ValueError('Prior query sequence proof mismatch')
 groups=collections.OrderedDict()
 for m in models:
  seqhash=hashlib.sha256(m['sequence'].encode()).hexdigest();e=execution[m['model_id']]
  if seqhash!=m['sequence_sha256'] or e['model_sequence_sha256']!=seqhash:raise ValueError('Model sequence proof mismatch')
  key=(seqhash,tuple(m.get('boundaries',[])))
  groups.setdefault(key,[]).append(m)
 out.mkdir();(out/'QUERY_IDENTITIES.json').write_text(json.dumps(source_keys,indent=2)+'\n');start=time.time();model_results=[]
 for index,aliases in enumerate(groups.values(),1):
  model=aliases[0];path=out/('model_%05d.outcomes.jsonl.gz'%index);counts=collections.Counter();sequence_hashes=set();n=0
  @functools.lru_cache(maxsize=128)
  def cached_sequence_evidence(sequence):return align_model(sequence,model)
  with gzip.open(path,'wt') as dst:
   for q,seq in expected.items():
    seqkey=hashlib.sha256(seq.encode()).hexdigest() if seq else None
    sequence_hashes.add(seqkey);evidence=cached_sequence_evidence(seq);counts[evidence['status']]+=1
    dst.write(json.dumps(dict(alignment_query_id=q,canonical_model_id=model['model_id'],evidence=evidence),separators=(',',':'))+'\n');n+=1
  if n!=len(expected):raise ValueError('Query conservation failure')
  model_results.append(dict(canonical_model_id=model['model_id'],all_model_aliases=[m['model_id'] for m in aliases],sequence_sha256=model['sequence_sha256'],boundaries=model.get('boundaries',[]),query_outcomes=n,status_counts=dict(counts),path=path.name,sha256=digest(path),unique_query_sequences=len(sequence_hashes),alignment_evaluations=cached_sequence_evidence.cache_info().misses,cache_max_entries=128))
  print(json.dumps(dict(completed_model_sequence_boundary_groups=index,total_groups=len(groups),model_aliases_completed=sum(len(r['all_model_aliases']) for r in model_results),queries_per_group=len(expected),seconds=time.time()-start)),flush=True)
 if {m for r in model_results for m in r['all_model_aliases']}!=set(mids):raise ValueError('Final model alias coverage mismatch')
 if digest(source)!=expected_source_sha or digest(stagepath)!=expected_stage_sha or digest(admissionpath)!=expected_admission_sha:raise ValueError('Source/admission changed during diagnostic execution')
 for rel,h in admission['model_artifacts'].items():
  if digest(d/rel)!=h:raise ValueError('Admitted model artifact changed during execution '+rel)
 report=dict(status='ALL_EVENT_QUERIES_ALL_PREPARED_MODELS_EXACT_SEMIGLOBAL_EVIDENCE_RETAINED',event=stage['event'],source_contexts=len(expected),raw_pair_rows=paircount,all_model_count=len(models),distinct_sequence_boundary_groups=len(groups),models=model_results,seconds=time.time()-start,source_pins={'independent_admission':expected_admission_sha,'source_shard':expected_source_sha,'event_stage':expected_stage_sha,'prepared_model_bundle':digest(bundlepath),'minimap_execution':digest(receiptpath),'query_identities':digest(out/'QUERY_IDENTITIES.json')},code_sha256=digest(__file__),backend_module_sha256=digest(pathlib.Path(__file__).with_name('circRNA_semiglobal_v2.py')),limits=['Evidence for all prepared competitors, not exhaustive biological mechanisms.','Minimap and exact edit-distance outputs remain separate.','Forced boundary alignment alone is not support; all competitor distances and ties retained.','Short reads and all errors retained without threshold; no physical-circle classification.','Repeated templates have more alignment opportunities; all model lengths/copy numbers remain in pinned model bundle.'])
 (out/'SEMIGLOBAL_EVENT_COMPLETE.json').write_text(json.dumps(report,indent=2)+'\n');return report
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--model-dir',required=True);p.add_argument('--model-stage-sha256',required=True);p.add_argument('--source-shard',required=True);p.add_argument('--source-sha256',required=True);p.add_argument('--output',required=True);p.add_argument('--admission-manifest',required=True);p.add_argument('--admission-sha256',required=True);a=p.parse_args();run(a.model_dir,a.model_stage_sha256,a.source_shard,a.source_sha256,a.output,a.admission_manifest,a.admission_sha256)
