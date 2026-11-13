"""One exhaustive event shard: prepare all model hypotheses and optionally align.

Whole-genome mappings are deliberately a separately pinned, once-per-query union
job. Their absence is NOT a negative result; outcome records explicitly retain it.
"""
import argparse,collections,hashlib,json,pathlib,sys
from .circRNA_full_query_validation_v2 import digest,event_key,contig,validate_model
from .circRNA_candidate_models_v2 import build_event_models,run_template_alignments

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--input-sha256',required=True);ap.add_argument('--shard-manifest',required=True);ap.add_argument('--shard-manifest-sha256',required=True);ap.add_argument('--reference-manifest',required=True);ap.add_argument('--reference-manifest-sha256',required=True);ap.add_argument('--external-models',required=True);ap.add_argument('--external-models-sha256',required=True);ap.add_argument('--native-complete',required=True);ap.add_argument('--native-complete-sha256',required=True);ap.add_argument('--native-catalog-rows',required=True);ap.add_argument('--output',required=True);ap.add_argument('--prepare-only',action='store_true');ap.add_argument('--threads',type=int,default=1);a=ap.parse_args()
 import pysam
 src=pathlib.Path(a.input);out=pathlib.Path(a.output)
 if out.exists():raise FileExistsError(out)
 for path,h in [(src,a.input_sha256),(pathlib.Path(a.reference_manifest),a.reference_manifest_sha256),(pathlib.Path(a.external_models),a.external_models_sha256),(pathlib.Path(a.shard_manifest),a.shard_manifest_sha256),(pathlib.Path(a.native_complete),a.native_complete_sha256)]:
  if digest(path)!=h:raise ValueError('Pinned input differs '+str(path))
 contexts=[json.loads(l) for l in src.open()];ids=[(r['event']['event_id'],r['library'],r['query']) for r in contexts]
 if not contexts or len(ids)!=len(set(ids)):raise ValueError('Missing/duplicate source context')
 events={r['event']['event_id'] for r in contexts}
 if len(events)!=1:raise ValueError('Shard must contain exactly one event')
 event=contexts[0]['event'];eid=event_key(event)
 if events!={eid}:raise ValueError('Event coordinates mismatch')
 shards=json.loads(pathlib.Path(a.shard_manifest).read_text())
 if shards['status']!='PASS_ALL50_ALL108_EXACT_BASELINE_AND_SOURCE_PARTITION' or shards['events']!=108:raise ValueError('Complete source shard roster required')
 matches=[r for r in shards['shards'] if r['event']==eid]
 if len(matches)!=1:raise ValueError('Missing/duplicate expected event shard')
 expected=matches[0]
 if expected['sha256']!=a.input_sha256 or pathlib.Path(expected['path']).name!=src.name or expected['contexts']!=len(contexts) or expected['pair_rows']!=sum(len(r['baseline_pairs']) for r in contexts):raise ValueError('Independent event roster/pair partition differs')

 ref=json.loads(pathlib.Path(a.reference_manifest).read_text());rp=pathlib.Path(ref['path']);stat=rp.stat()
 if (stat.st_size,stat.st_mtime_ns)!=(ref['size'],ref['mtime_ns']) or digest(str(rp)+'.fai')!=ref['fai_sha256']:raise ValueError('Prehashed reference metadata/FAI changed')
 if not ref.get('sha256_verified_before_dispatch'):raise ValueError('Reference full-byte hash not independently verified before dispatch')
 fa=pysam.FastaFile(str(rp));lengths=dict(zip(fa.references,fa.lengths))
 def cname(c):
  matches=[n for n in lengths if contig(n)==contig(c)]
  if len(matches)!=1:raise ValueError('Missing/ambiguous reference contig '+c)
  return matches[0]
 def fetch(c,s,e):
  name=cname(c)
  if not 0<=s<e<=lengths[name]:raise ValueError('Reference bounds')
  return fa.fetch(name,s,e)
 native=json.loads(pathlib.Path(a.native_complete).read_text())
 if native['status']!='COMPLETE_NATIVE_CATALOG_ENUMERATION_NO_EXCLUSIONS' or native['events']!=108:raise ValueError('Incomplete native model enumeration')
 if native['outputs']['external_chains_by_event.json']['sha256']!=a.external_models_sha256 or digest(a.native_catalog_rows)!=native['outputs']['matching_native_catalog_rows.jsonl.gz']['sha256']:raise ValueError('Native model/source binding differs')
 all_external=json.loads(pathlib.Path(a.external_models).read_text());matches=[r for r in all_external['events'] if r['event']['event_id']==eid]
 if len(matches)!=1:raise ValueError('External event model missing/duplicated')
 external=matches[0]
 if external['status'] not in ('native_candidate_hypotheses_available','no_compatible_native_exonic_candidate_model'):raise ValueError('External model coverage status required')
 for m in external['external_chains']:
  if m['source_pin']!={'artifact':'matching_native_catalog_rows','sha256':native['outputs']['matching_native_catalog_rows.jsonl.gz']['sha256']}:raise ValueError('External chain source pin differs')
 for r in contexts:
  r['source']['source_details']=r['source'].get('source_pin');r['source']['source_pin']={'artifact':str(src),'sha256':a.input_sha256}
 bundle=build_event_models(event,contexts,fetch,lambda c:lengths[cname(c)],ref['sha256'],external.get('external_chains',[]))
 for m in bundle['templates']:
  test=dict(m,alignment_parameters=['prepared_not_executed'],paf_records=[]);validate_model(test,event)
 queries=[]
 for r in contexts:
  qid='Q'+hashlib.sha256(json.dumps([r['library'],r['query']],separators=(',',':')).encode()).hexdigest()
  if r.get('alignment_query_id',qid)!=qid:raise ValueError('Global alignment identity differs')
  queries.append(dict(id=qid,library=r['library'],query=r['query'],sequence=r['source'].get('full_forward_sequence')))
 if len({q['id'] for q in queries})!=len(queries):raise ValueError('Alignment identity collision')
 out.mkdir(parents=True);(out/'PREPARED_MODEL_BUNDLE.json').write_text(json.dumps(bundle,indent=2)+'\n');(out/'QUERY_IDENTITIES.json').write_text(json.dumps([{k:v for k,v in q.items() if k!='sequence'} for q in queries],indent=2)+'\n')
 (out/'BASELINE_ROSTER.json').write_text(json.dumps(dict(event=eid,ordered_source_keys=ids,source_contexts=len(contexts),raw_pair_rows=sum(len(r['baseline_pairs']) for r in contexts),baseline_status_pair_counts=dict(collections.Counter(p['status'] for r in contexts for p in r['baseline_pairs']))),indent=2)+'\n')
 if not a.prepare_only:run_template_alignments(bundle,queries,out/'alignments',a.threads,outcome_storage='external')
 if digest(src)!=a.input_sha256:raise ValueError('Source changed during model construction')
 if (rp.stat().st_size,rp.stat().st_mtime_ns)!=(ref['size'],ref['mtime_ns']):raise ValueError('Reference metadata changed')
 summary=dict(status='PREPARED_NOT_ALIGNED' if a.prepare_only else 'EVENT_TEMPLATE_EXECUTION_COMPLETED_CHECK_PER_MODEL_OUTCOMES',event=eid,source_contexts=len(contexts),raw_pair_rows=sum(len(r['baseline_pairs']) for r in contexts),template_count=len(bundle['templates']),unresolved_source_models=sum(x['status']!='candidate_hypothesis_constructed' for x in bundle['source_model_outcomes']),unresolved_models=bundle['unresolved_models'],whole_genome_alternatives='not_evaluated_by_this_event_job',input_sha256=a.input_sha256,reference_manifest_sha256=a.reference_manifest_sha256,external_model_manifest_sha256=a.external_models_sha256,code_sha256=digest(__file__),physical_circularity='not_established')
 (out/'EVENT_STAGE_RESULT.json').write_text(json.dumps(summary,indent=2)+'\n');fa.close();print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
