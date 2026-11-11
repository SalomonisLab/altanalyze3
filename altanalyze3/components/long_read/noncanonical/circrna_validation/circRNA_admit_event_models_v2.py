"""Independently reconstruct pinned event models and verify all execution outcomes.
No alignment rerun and no admission based on self-consistent model edits alone.
"""
import argparse,collections,gzip,hashlib,json,pathlib,subprocess,sys,tempfile
from .circRNA_full_query_validation_v2 import digest,parse_paf

def main():
 p=argparse.ArgumentParser();p.add_argument('--frozen-bundle',required=True);p.add_argument('--frozen-manifest-sha256',required=True);p.add_argument('--shard-manifest',required=True);p.add_argument('--shard-manifest-sha256',required=True);p.add_argument('--array-index',type=int,required=True);p.add_argument('--reference-manifest',required=True);p.add_argument('--reference-manifest-sha256',required=True);p.add_argument('--native-directory',required=True);p.add_argument('--native-complete-sha256',required=True);p.add_argument('--model-dir',required=True);p.add_argument('--output',required=True);a=p.parse_args()
 h=pathlib.Path(a.frozen_bundle);d=pathlib.Path(a.model_dir);out=pathlib.Path(a.output);native=pathlib.Path(a.native_directory)
 tracked=['EVENT_STAGE_RESULT.json','PREPARED_MODEL_BUNDLE.json','QUERY_IDENTITIES.json','BASELINE_ROSTER.json','alignments/ALIGNMENT_EXECUTION.json','alignments/MODEL_BUNDLE.json']
 initial_hashes={name:digest(d/name) for name in tracked}
 if out.exists():raise FileExistsError(out)
 if digest(h/'SOURCE_HASHES.json')!=a.frozen_manifest_sha256 or digest(a.shard_manifest)!=a.shard_manifest_sha256 or digest(a.reference_manifest)!=a.reference_manifest_sha256 or digest(native/'NATIVE_MODELS_COMPLETE.json')!=a.native_complete_sha256:raise ValueError('Independent input pin mismatch')
 for name,expected in json.loads((h/'SOURCE_HASHES.json').read_text()).items():
  if digest(h/name)!=expected:raise ValueError('Frozen reconstruction source differs '+name)
 shards=json.loads(pathlib.Path(a.shard_manifest).read_text())
 if not 1<=a.array_index<=108:raise ValueError('Array index outside full event roster')
 selected=shards['shards'][a.array_index-1];source=pathlib.Path(selected['path'])
 if not 1<=a.array_index<=108 or digest(source)!=selected['sha256']:raise ValueError('Source shard identity differs')
 nc=json.loads((native/'NATIVE_MODELS_COMPLETE.json').read_text());external_sha=nc['outputs']['external_chains_by_event.json']['sha256']
 with tempfile.TemporaryDirectory(prefix='circ-model-admission-') as tmp:
  rebuilt=pathlib.Path(tmp)/'rebuilt';cmd=[sys.executable,str(h/'circRNA_validate_event_v2.py'),'--input',str(source),'--input-sha256',selected['sha256'],'--shard-manifest',a.shard_manifest,'--shard-manifest-sha256',a.shard_manifest_sha256,'--reference-manifest',a.reference_manifest,'--reference-manifest-sha256',a.reference_manifest_sha256,'--external-models',str(native/'external_chains_by_event.json'),'--external-models-sha256',external_sha,'--native-complete',str(native/'NATIVE_MODELS_COMPLETE.json'),'--native-complete-sha256',a.native_complete_sha256,'--native-catalog-rows',str(native/'matching_native_catalog_rows.jsonl.gz'),'--output',str(rebuilt),'--prepare-only']
  subprocess.run(cmd,check=True)
  for name in ('PREPARED_MODEL_BUNDLE.json','QUERY_IDENTITIES.json','BASELINE_ROSTER.json'):
   if digest(rebuilt/name)!=digest(d/name):raise ValueError('Reference/source reconstructed model differs '+name)
 stage=json.loads((d/'EVENT_STAGE_RESULT.json').read_text());bundle=json.loads((d/'PREPARED_MODEL_BUNDLE.json').read_text());receipt=json.loads((d/'alignments/ALIGNMENT_EXECUTION.json').read_text());models={m['model_id']:m for m in bundle['templates']};runs={r['model_id']:r for r in receipt['models']}
 if len(models)!=len(bundle['templates']) or len(runs)!=len(receipt['models']) or set(models)!=set(runs):raise ValueError('Model/execution identity mismatch')
 if stage['status']!='EVENT_TEMPLATE_EXECUTION_COMPLETED_CHECK_PER_MODEL_OUTCOMES' or stage['template_count']!=len(models) or stage['source_contexts']!=selected['contexts'] or stage['raw_pair_rows']!=selected['pair_rows']:raise ValueError('Incomplete/inconsistent event execution stage')
 if stage['input_sha256']!=selected['sha256'] or stage['event']!=selected['event'] or stage['code_sha256']!=digest(h/'circRNA_validate_event_v2.py'):raise ValueError('Execution source identity differs')
 if json.loads((d/'alignments/MODEL_BUNDLE.json').read_text())!=bundle:raise ValueError('Executed model bundle differs')
 query_sequences={}
 with source.open() as f:
  for line in f:
   r=json.loads(line);q=r['alignment_query_id']
   if q in query_sequences:raise ValueError('Duplicate source query')
   query_sequences[q]=r['source'].get('full_forward_sequence')
 query_hashes={q:hashlib.sha256(seq.encode()).hexdigest() for q,seq in query_sequences.items() if seq}
 if receipt['query_sequence_sha256']!=query_hashes or receipt['reference_sha256']!=bundle['reference_sha256']:raise ValueError('Execution sequence/reference identity mismatch')
 qfa_sha=digest(d/'alignments/queries.fa');expected_fa=hashlib.sha256()
 for q,seq in query_sequences.items():
  if seq:expected_fa.update(('>'+q+'\n'+seq+'\n').encode())
 if qfa_sha!=expected_fa.hexdigest():raise ValueError('Executed query FASTA differs from exact source sequences')
 checked={};proofs=[]
 for ident,run in runs.items():
  model=models[ident];canonical=run.get('executed_model_id',ident)
  if canonical not in runs or models[canonical]['sequence']!=model['sequence']:raise ValueError('Invalid exact-sequence execution alias')
  for suffix,field in [('.fa','template_sha256'),('.paf','paf_sha256')]:
   if digest(d/'alignments'/(ident+suffix))!=run[field]:raise ValueError('Execution artifact differs '+ident+suffix)
  if qfa_sha!=run['query_fasta_sha256'] or run['model_sequence_sha256']!=model['sequence_sha256']:raise ValueError('Execution input sequence proof differs')
  expected_command=['minimap2','-x','map-hifi','-c','--cs=long','--secondary=yes','-N','50','-p','0','-t','1',str(d/'alignments'/(canonical+'.fa')),str(d/'alignments/queries.fa')]
  if run['command']!=expected_command:raise ValueError('Alignment parameters/path differ from frozen single-thread event execution')
  if canonical in checked:
   previous=checked[canonical]
   for field in ('paf_sha256','template_sha256','model_sequence_sha256','outcomes_sha256','outcomes_file','outcomes_count','status','returncode','command'):
    if run[field]!=previous[field]:raise ValueError('Alias execution receipt differs '+field)
   continue
  fa=(d/'alignments'/(ident+'.fa')).read_text().splitlines()
  if fa[0]!='>'+canonical or ''.join(fa[1:])!=model['sequence']:raise ValueError('Executed template bytes differ from reconstructed reference')
  counts=collections.Counter();hashes=collections.defaultdict(hashlib.sha256)
  if run['status']=='completed':
   if run['returncode']!=0:raise ValueError('Successful execution returncode differs')
   with (d/'alignments'/(ident+'.paf')).open() as f:
    for line in f:
     parsed=parse_paf(line);q=parsed['query']
     if q not in query_hashes or parsed['query_length']!=len(query_sequences[q]) or parsed['target']!=canonical or parsed['target_length']!=len(model['sequence']):raise ValueError('PAF source/model identity differs')
     counts[q]+=1;hashes[q].update((json.dumps(parsed,sort_keys=True,separators=(',',':'))+'\n').encode())
  elif run['status']!='failed':raise ValueError('Unknown model execution status')
  op=d/'alignments'/run['outcomes_file']
  if digest(op)!=run['outcomes_sha256']:raise ValueError('Outcome artifact pin mismatch')
  seen=set()
  with gzip.open(op,'rt') as f:
   for line in f:
    row=json.loads(line);q=row['id']
    if q not in query_sequences or q in seen:raise ValueError('Unknown/duplicate query outcome')
    seen.add(q);available=bool(query_sequences[q]);expected_count=counts[q] if run['status']=='completed' and available else None
    expected_status=('mapped' if counts[q] else 'unmapped_retained') if run['status']=='completed' and available else 'unresolved_missing_sequence' if not available else 'failed_retained'
    if row['mapping_count']!=expected_count or row['status']!=expected_status or row['paf_records_sha256']!=hashes[q].hexdigest():raise ValueError('Per-query outcome differs from actual execution')
  if seen!=set(query_sequences) or run['outcomes_count']!=len(seen):raise ValueError('Incomplete query outcome universe')
  checked[canonical]=run;proofs.append(dict(executed_model_id=canonical,query_outcomes=len(seen),mapping_records=sum(counts.values()),status=run['status']))
 artifacts={name:digest(d/name) for name in tracked}
 if artifacts!=initial_hashes or digest(source)!=selected['sha256'] or digest(a.shard_manifest)!=a.shard_manifest_sha256 or digest(a.reference_manifest)!=a.reference_manifest_sha256 or digest(native/'NATIVE_MODELS_COMPLETE.json')!=a.native_complete_sha256 or digest(h/'SOURCE_HASHES.json')!=a.frozen_manifest_sha256:raise ValueError('Source/model admission inputs changed during verification')
 report=dict(status='INDEPENDENT_EVENT_MODEL_ARTIFACT_ADMISSION',event=selected['event'],admission_basis='Exact rebuilt model bundle/rosters from independently pinned source shards, native catalog and reference; all actual PAF/template/outcome hashes and query identities verified. No alignment rerun.',source_shard_sha256=selected['sha256'],frozen_source_manifest_sha256=a.frozen_manifest_sha256,shard_manifest_sha256=a.shard_manifest_sha256,reference_manifest_sha256=a.reference_manifest_sha256,native_complete_sha256=a.native_complete_sha256,model_artifacts=artifacts,all_models=len(models),executed_sequence_groups=len(checked),queries=len(query_sequences),verified_executions=proofs,code_sha256=digest(__file__))
 out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2)+'\n');print('ADMITTED',selected['event'],len(models),len(query_sequences),flush=True)
if __name__=='__main__':main()
