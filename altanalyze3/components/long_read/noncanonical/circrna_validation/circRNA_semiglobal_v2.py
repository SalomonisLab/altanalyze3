"""Exact all-query/end-free-target edit-distance evidence, separate from minimap.

Edlib HW reports optimal locations; NW supplies one optimal traceback per reported
location. This is not exhaustive enumeration of all equally optimal tracebacks.
No distance/anchor cutoff, topology classifier or feature exclusion is applied.
"""
import hashlib,importlib.metadata
from .circRNA_full_query_validation_v2 import cigar_boundary_evidence,cigar_ops

def rc(s):return s.translate(str.maketrans('ACGTNacgtn','TGCANtgcan'))[::-1]
def sha(s):return hashlib.sha256(s.encode()).hexdigest()
def align_model(query,model):
 import edlib
 target=model['sequence'];boundaries=model.get('boundaries',[])
 if not query:return {'status':'unresolved_missing_query_sequence','model_id':model['model_id'],'orientations':[]}
 if not target:return {'status':'unresolved_missing_model_sequence','model_id':model['model_id'],'orientations':[]}
 if model.get('sequence_sha256',sha(target))!=sha(target):raise ValueError('Model sequence hash mismatch')
 if any(not isinstance(b,int) or not 0<b<len(target) for b in boundaries) or len(boundaries)!=len(set(boundaries)):raise ValueError('Invalid model boundaries')
 orientations=[]
 for strand,seq in [('+',query),('-',rc(query))]:
  result=edlib.align(seq,target,mode='HW',task='locations',k=-1)
  if result['editDistance']<0:raise ValueError('Unlimited exact alignment unexpectedly failed')
  placements=[]
  for start,last in result['locations']:
   start=0 if start is None else start;end=last+1
   if not 0<=start<=end<=len(target):raise ValueError('Invalid exact alignment location')
   segment=target[start:end]
   traced=edlib.align(seq,segment,mode='NW',task='path',k=-1)
   if traced['editDistance']!=result['editDistance']:raise ValueError('HW/NW optimum mismatch')
   cigar=traced['cigar'] if segment else str(len(seq))+'I'
   ops=cigar_ops(cigar)
   if sum(n for n,o in ops if o in 'MI=X')!=len(seq) or sum(n for n,o in ops if o in 'MDN=X')!=len(segment):raise ValueError('Exact traceback spans differ')
   edit=sum(n for n,o in ops if o in 'IDX')
   if edit!=result['editDistance']:raise ValueError('Traceback edit distance differs')
   p=dict(query_length=len(query),query_start=0,query_end=len(query),strand=strand,target_start=start,target_end=end,tags={'cg':cigar})
   geometry=cigar_boundary_evidence(p,boundaries)
   # Edlib emits =/X runs. Preserve original run-based fields, and separately
   # measure uninterrupted aligned flanks across mismatches, stopping at I/D.
   merged=[]
   for n,op in ops:
    op='M' if op in '=XM' else op
    if merged and merged[-1][1]==op:merged[-1]=(merged[-1][0]+n,op)
    else:merged.append((n,op))
   merged_p=dict(p,tags={'cg':''.join(str(n)+op for n,op in merged)})
   uninterrupted=cigar_boundary_evidence(merged_p,boundaries)
   if len(uninterrupted['observations'])!=len(geometry['observations']):raise ValueError('Boundary geometry representation differs')
   for original,coalesced in zip(geometry['observations'],uninterrupted['observations']):
    original['uninterrupted_aligned_flank_left_including_mismatches']=coalesced['adjacent_aligned_bases_left']
    original['uninterrupted_aligned_flank_right_including_mismatches']=coalesced['adjacent_aligned_bases_right']
   placements.append(dict(target_start0=start,target_end0=end,cigar=cigar,boundaries=geometry['observations'],matched_bases=sum(n for n,o in ops if o=='='),mismatch_bases=sum(n for n,o in ops if o=='X'),inserted_bases=sum(n for n,o in ops if o=='I'),deleted_bases=sum(n for n,o in ops if o=='D')))
  orientations.append(dict(strand=strand,edit_distance=result['editDistance'],distance_per_query_base=result['editDistance']/len(query),locations=placements))
 return dict(status='exact_semiglobal_evidence_measured',model_id=model['model_id'],query_sequence_sha256=sha(query),model_sequence_sha256=sha(target),query_length=len(query),target_length=len(target),best_edit_distance=min(x['edit_distance'] for x in orientations),orientations=orientations,backend='edlib',backend_version=importlib.metadata.version('edlib'),parameters=dict(mode='HW',task='locations',k=-1,traceback='NW each reported location',both_query_orientations=True),limits=['All query bases participate; target end gaps are free.','Reported locations are Edlib optimal locations, with one NW traceback each, not every optimal traceback.','Edit distance has equal unit substitution/insertion/deletion costs; not comparable to minimap alignment scores.','Best alignment exists even for unrelated sequences; no support threshold or physical topology classification.'])

def align_models(query,models):
 ids=[m['model_id'] for m in models]
 if len(ids)!=len(set(ids)):raise ValueError('Duplicate model IDs')
 return [align_model(query,m) for m in models]
