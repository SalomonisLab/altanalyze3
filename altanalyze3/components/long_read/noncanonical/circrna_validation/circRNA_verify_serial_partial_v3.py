"""Additive admission of bounded complete serial exact-model groups.

Never alters the frozen serial calculation or its output. Used to decide whether
redundant partition work can be stopped safely after a serial memory failure.
"""
import argparse,collections,datetime,gzip,hashlib,json
from pathlib import Path
from .circRNA_merge_exact_partitions_v3 import sha,query_id,read_json


def main():
 p=argparse.ArgumentParser()
 p.add_argument('--model-dir',required=True)
 p.add_argument('--source-shard',required=True)
 p.add_argument('--admission-manifest',required=True)
 p.add_argument('--serial-dir',required=True)
 p.add_argument('--completed-through',type=int,required=True)
 p.add_argument('--output',required=True)
 a=p.parse_args()
 model_dir=Path(a.model_dir);source=Path(a.source_shard);admission=read_json(Path(a.admission_manifest))
 stage=read_json(model_dir/'EVENT_STAGE_RESULT.json');bundle=read_json(model_dir/'PREPARED_MODEL_BUNDLE.json')
 if admission['status']!='INDEPENDENT_EVENT_MODEL_ARTIFACT_ADMISSION':raise ValueError('Missing model admission')
 for rel,digest in admission['model_artifacts'].items():
  if sha(model_dir/rel)!=digest:raise ValueError('Model artifact hash mismatch '+rel)
 source_hash=sha(source)
 if source_hash!=admission['source_shard_sha256'] or source_hash!=stage['input_sha256']:raise ValueError('Source hash mismatch')
 groups=collections.OrderedDict()
 for model in bundle['templates']:
  seqhash=hashlib.sha256(model['sequence'].encode()).hexdigest()
  if seqhash!=model['sequence_sha256']:raise ValueError('Model sequence hash mismatch')
  groups.setdefault((seqhash,tuple(model.get('boundaries',[]))),[]).append(model)
 if not 0<a.completed_through<=len(groups):raise ValueError('Invalid completed bound')
 ids=set();n=0;pair_rows=0
 with source.open() as f:
  for line in f:
   row=json.loads(line);qid=query_id(row['library'],row['query'])
   if qid!=row['alignment_query_id'] or qid in ids:raise ValueError('Source query identity mismatch or duplicate')
   ids.add(qid);n+=1;pair_rows+=len(row['baseline_pairs'])
 if n!=stage['source_contexts'] or pair_rows!=stage['raw_pair_rows']:raise ValueError('Source count mismatch')
 receipts=[]
 for index in range(1,a.completed_through+1):
  model=list(groups.values())[index-1][0]
  path=Path(a.serial_dir)/f'model_{index:05d}.outcomes.jsonl.gz'
  seen=set();statuses=collections.Counter()
  with gzip.open(path,'rt') as f:
   for line in f:
    x=json.loads(line);qid=x['alignment_query_id']
    if qid not in ids or qid in seen or x['canonical_model_id']!=model['model_id']:raise ValueError(f'Group {index} identity/roster mismatch')
    seen.add(qid);statuses[x['evidence']['status']]+=1
  if seen!=ids:raise ValueError(f'Incomplete serial group {index}: {len(seen)} of {n}')
  receipts.append({'group_index':index,'model_id':model['model_id'],'payload_sha256':sha(path),'query_count':len(seen),'status_counts':dict(statuses)})
 out={'status':'PASS_PINNED_SERIAL_PARTIAL_ALL_QUERIES_PER_GROUP','checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':source_hash,'admission_sha256':sha(Path(a.admission_manifest)),'model_bundle_sha256':sha(model_dir/'PREPARED_MODEL_BUNDLE.json'),'source_contexts':n,'raw_pair_rows':pair_rows,'all_model_groups':len(groups),'completed_through':a.completed_through,'verified_groups':receipts,'interpretation':'These exact model-group files are complete and may be reused by the full merge; no circle topology inference.'}
 target=Path(a.output)
 if target.exists():raise FileExistsError(target)
 target.parent.mkdir(parents=True,exist_ok=True)
 target.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':out['status'],'source_contexts':n,'completed_through':a.completed_through,'all_model_groups':len(groups)}))
if __name__=='__main__':main()
