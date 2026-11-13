"""Pinned, additive model-group partition of the unchanged exact read evaluator.

This runs whole source-query rosters for specified model sequence/boundary groups.
It does not infer physical topology or drop unresolved sources. Its outputs can
be admitted only after independent parity with the serial frozen evaluator.
"""
import argparse,collections,functools,gzip,hashlib,json,pathlib,time
from .circRNA_full_query_validation_v2 import digest,event_key
from .circRNA_semiglobal_event_v2 import fasta
from .circRNA_semiglobal_v2 import align_model

def run(a):
 d=pathlib.Path(a.model_dir);source=pathlib.Path(a.source_shard);admissionpath=pathlib.Path(a.admission_manifest);out=pathlib.Path(a.output)
 if out.exists():raise FileExistsError(out)
 stagepath=d/'EVENT_STAGE_RESULT.json';bundlepath=d/'PREPARED_MODEL_BUNDLE.json';receiptpath=d/'alignments/ALIGNMENT_EXECUTION.json'
 if digest(stagepath)!=a.model_stage_sha256 or digest(source)!=a.source_sha256 or digest(admissionpath)!=a.admission_sha256:raise ValueError('Pinned stage/source/admission changed')
 admission=json.loads(admissionpath.read_text());stage=json.loads(stagepath.read_text());bundle=json.loads(bundlepath.read_text());receipt=json.loads(receiptpath.read_text())
 if admission['status']!='INDEPENDENT_EVENT_MODEL_ARTIFACT_ADMISSION' or admission['source_shard_sha256']!=a.source_sha256 or admission['model_artifacts']['EVENT_STAGE_RESULT.json']!=a.model_stage_sha256:raise ValueError('Independent admission mismatch')
 for rel,h in admission['model_artifacts'].items():
  if digest(d/rel)!=h:raise ValueError('Admitted model artifact changed: '+rel)
 if stage['status']!='EVENT_TEMPLATE_EXECUTION_COMPLETED_CHECK_PER_MODEL_OUTCOMES' or stage['input_sha256']!=a.source_sha256:raise ValueError('Model stage incomplete')
 if json.loads((d/'alignments/MODEL_BUNDLE.json').read_text())!=bundle:raise ValueError('Model bundle execution mismatch')
 models=bundle['templates'];mids=[m['model_id'] for m in models]
 if len(mids)!=len(set(mids)) or len(models)!=stage['template_count']:raise ValueError('Model roster mismatch')
 execution={m['model_id']:m for m in receipt['models']}
 if len(execution)!=len(receipt['models']) or set(execution)!=set(mids):raise ValueError('Minimap execution roster mismatch')
 expected={};paircount=0
 with source.open() as f:
  for line in f:
   r=json.loads(line);qid='Q'+hashlib.sha256(json.dumps([r['library'],r['query']],separators=(',',':')).encode()).hexdigest()
   if r.get('alignment_query_id')!=qid or qid in expected or event_key(r['event'])!=stage['event']:raise ValueError('Source query identity mismatch')
   expected[qid]=r['source'].get('full_forward_sequence');paircount+=len(r['baseline_pairs'])
 if len(expected)!=stage['source_contexts'] or paircount!=stage['raw_pair_rows'] or fasta(d/'alignments/queries.fa')!={q:s for q,s in expected.items() if s}:raise ValueError('Source roster mismatch')
 if receipt['query_sequence_sha256']!={q:hashlib.sha256(s.encode()).hexdigest() for q,s in expected.items() if s}:raise ValueError('Prior query sequence proof mismatch')
 groups=collections.OrderedDict()
 for m in models:
  seqhash=hashlib.sha256(m['sequence'].encode()).hexdigest();e=execution[m['model_id']]
  if seqhash!=m['sequence_sha256'] or e['model_sequence_sha256']!=seqhash:raise ValueError('Model sequence proof mismatch')
  groups.setdefault((seqhash,tuple(m.get('boundaries',[]))),[]).append(m)
 total=len(groups)
 if not 1<=a.start_group<=a.end_group<=total:raise ValueError('Invalid requested group range')
 out.mkdir();start=time.time();results=[]
 for index in range(a.start_group,a.end_group+1):
  aliases=list(groups.values())[index-1];model=aliases[0];path=out/('model_%05d.outcomes.jsonl.gz'%index);counts=collections.Counter();sequence_hashes=set();n=0
  @functools.lru_cache(maxsize=128)
  def cached(sequence):return align_model(sequence,model)
  with gzip.open(path,'wt') as dst:
   for q,seq in expected.items():
    sequence_hashes.add(hashlib.sha256(seq.encode()).hexdigest() if seq else None);evidence=cached(seq);counts[evidence['status']]+=1
    dst.write(json.dumps(dict(alignment_query_id=q,canonical_model_id=model['model_id'],evidence=evidence),separators=(',',':'))+'\n');n+=1
  if n!=len(expected):raise ValueError('Query conservation failure')
  results.append(dict(group_index=index,canonical_model_id=model['model_id'],all_model_aliases=[m['model_id'] for m in aliases],sequence_sha256=model['sequence_sha256'],boundaries=model.get('boundaries',[]),query_outcomes=n,status_counts=dict(counts),path=path.name,sha256=digest(path),unique_query_sequences=len(sequence_hashes),alignment_evaluations=cached.cache_info().misses,cache_max_entries=128))
  print(json.dumps(dict(completed_group=index,group_end=a.end_group,total_groups=total,queries_per_group=len(expected),seconds=time.time()-start)),flush=True)
 if digest(stagepath)!=a.model_stage_sha256 or digest(source)!=a.source_sha256 or digest(admissionpath)!=a.admission_sha256:raise ValueError('Pinned source changed during execution')
 for rel,h in admission['model_artifacts'].items():
  if digest(d/rel)!=h:raise ValueError('Admitted artifact changed during execution')
 report=dict(status='PARTITION_ALL_SELECTED_MODELS_ALL_SOURCE_QUERIES_EXACT',event=stage['event'],source_contexts=len(expected),raw_pair_rows=paircount,all_model_count=len(models),distinct_sequence_boundary_groups=total,requested_group_start=a.start_group,requested_group_end=a.end_group,models=results,seconds=time.time()-start,source_pins=dict(independent_admission=a.admission_sha256,source_shard=a.source_sha256,event_stage=a.model_stage_sha256,prepared_model_bundle=digest(bundlepath),minimap_execution=digest(receiptpath)),backend_module_sha256=digest(pathlib.Path(__file__).with_name('circRNA_semiglobal_v2.py')),serial_dispatch_module_sha256=digest(pathlib.Path(__file__).with_name('circRNA_semiglobal_event_v2.py')),code_sha256=digest(__file__),limits=['Partitioned execution uses the unchanged exact align_model; independent parity with serial group outputs is required before merge.','All source query contexts and all requested model groups are retained.','No physical topology or calibrated circle support is inferred.'])
 (out/'GROUP_PART_COMPLETE.json').write_text(json.dumps(report,indent=2)+'\n');return report
if __name__=='__main__':
 p=argparse.ArgumentParser()
 for x in ['model_dir','model_stage_sha256','source_shard','source_sha256','output','admission_manifest','admission_sha256']:p.add_argument('--'+x.replace('_','-'),required=True)
 p.add_argument('--start-group',type=int,required=True);p.add_argument('--end-group',type=int,required=True)
 r=run(p.parse_args());print(json.dumps({k:r[k] for k in ('status','event','source_contexts','distinct_sequence_boundary_groups','requested_group_start','requested_group_end','seconds')}),flush=True)
