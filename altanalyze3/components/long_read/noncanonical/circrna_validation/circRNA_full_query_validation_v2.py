"""Additive full-query circRNA evidence, schema 2; no candidate filtering/classifier.

Frozen scanner identities/statuses remain authoritative baseline records. This module
adds measured sequence/path evidence and unresolved states. It never changes counts,
normalizes endpoints, substitutes models, or claims physical RNA circularity.
Coordinates: genomic event endpoints 1-based inclusive; query/model 0-based half-open.
"""
from __future__ import annotations
import argparse,collections,gzip,hashlib,json,pathlib,re,sqlite3,tempfile

CIGAR=re.compile(r'(\d+)([MIDNSHP=X])')
RC=str.maketrans('ACGTNacgtn','TGCANtgcan')

def revcomp(s):return s.translate(RC)[::-1]
def digest(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def cigar_ops(s):
 a=[(int(n),op) for n,op in CIGAR.findall(s)]
 if not a or any(n<=0 for n,_ in a) or ''.join(str(n)+op for n,op in a)!=s:raise ValueError('Malformed CIGAR')
 return a
def contig(c):return c[3:] if c.startswith('chr') else c

def event_key(e):
 return 'BSJ|{}|{}|{}|{}|{}'.format(e['chrom'],e['donor_strand'],e['acceptor_strand'],int(e['donor']),int(e['acceptor']))
def pair_event(p):
 if p['joint_type']!='back':raise ValueError('Baseline row is not a back-junction')
 return 'BSJ|{}|{}|{}|{}|{}'.format(p['chrom'],p['a_strand'],p['b_strand'],int(p['donor']),int(p['acceptor']))
def alignment_key(a):return (contig(a['chrom']),int(a['start0']),'-' if a['reverse'] else '+',a['cigar'])

def terminal_clip_signature(cigar):
 """Equivalent terminal H/S descriptors retain exact aligned operation path."""
 ops=cigar_ops(cigar);lead=0;trail=0
 while ops and ops[0][1] in 'SH':lead+=ops.pop(0)[0]
 while ops and ops[-1][1] in 'SH':trail+=ops.pop()[0]
 if not ops or any(op in 'SH' for n,op in ops):raise ValueError('Internal/empty clipped alignment')
 return lead,tuple(ops),trail

def equivalent_alignment_key(a):
 return (contig(a['chrom']),int(a['start0']),'-' if a['reverse'] else '+',terminal_clip_signature(a['cigar']))

def hard_clipped_source_slice(a,sequence):
 """Verify actual stored bases against the full original-orientation sequence."""
 ops=cigar_ops(a['cigar']);lead_h=0;tail_h=0
 for n,op in ops:
  if op not in 'SH':break
  if op=='H':lead_h+=n
 for n,op in reversed(ops):
  if op not in 'SH':break
  if op=='H':tail_h+=n
 if not lead_h and not tail_h:return {'status':'not_hard_clipped'}
 if not sequence or not a.get('forward_sequence'):return {'status':'unresolved_missing_hardclip_sequence'}
 lo,hi=(tail_h,len(sequence)-lead_h) if a['reverse'] else (lead_h,len(sequence)-tail_h)
 exact=0<=lo<=hi<=len(sequence) and a['forward_sequence'].upper()==sequence[lo:hi].upper()
 return dict(status='verified_full_sequence_slice' if exact else 'hardclip_sequence_slice_mismatch',original_query_start0=lo,original_query_end0=hi,actual_stored_bases=len(a['forward_sequence']))

def alignment_geometry(a):
 ops=cigar_ops(a['cigar']);n=sum(k for k,o in ops if o in 'MIS H=X'.replace(' ',''));lead=0
 for k,o in ops:
  if o in 'SH':lead+=k
  else:break
 aligned=sum(k for k,o in ops if o in 'MI=X');qs,qe=lead,lead+aligned
 if a['reverse']:qs,qe=n-qe,n-qs
 end=int(a['start0'])+sum(k for k,o in ops if o in 'MDN=X')
 return dict(query_start=qs,query_end=qe,full_query_length=n,reference_start0=int(a['start0']),reference_end0=end,hard_clip_nt=sum(k for k,o in ops if o=='H'),reference_skip_nt=sum(k for k,o in ops if o=='N'))

def reference_exons_1based(a):
 """Match frozen scanner _exons_from_cigar: D extends exon; N splits it."""
 ref=int(a['start0']);start=None;exons=[]
 for n,op in cigar_ops(a['cigar']):
  if op in 'M=X D'.replace(' ',''):
   if start is None:start=ref
   ref+=n
  elif op=='N':
   if start is not None and ref>=start+1:exons.append((start+1,ref))
   start=None;ref+=n
 if start is not None and ref>=start+1:exons.append((start+1,ref))
 return exons

def query_cut_reference_context(a,cut):
 """Map an original-orientation query cut through actual CIGAR, including indels.

 No linear overlap subtraction. Multiple reference positions at deletion/skips or
 cuts inside insertion are reported explicitly instead of choosing a breakpoint.
 """
 g=alignment_geometry(a);n=g['full_query_length'];qc=n-int(cut) if a['reverse'] else int(cut)
 if not 0<=qc<=n:raise ValueError('Query cut outside read')
 q=0;t=int(a['start0']);hits=[]
 for length,op in cigar_ops(a['cigar']):
  qend=q+length if op in 'MISH=X' else q;tend=t+length if op in 'MDN=X' else t
  if op in 'M=X' and q<=qc<=qend:hits.append(dict(operation=op,reference_cut0=t+qc-q,query_offset_in_operation=qc-q))
  elif op=='I' and q<=qc<=qend:hits.append(dict(operation=op,reference_cut0=t,query_offset_in_operation=qc-q,insertion_length=length,inside_insertion=q<qc<qend))
  elif op in 'DN' and qc==q:hits.append(dict(operation=op,reference_cut_interval0=[t,tend],reference_gap_length=length))
  q,t=qend,tend
 return dict(original_query_cut0=cut,reference_candidates=hits,inside_insertion=any(h.get('inside_insertion') for h in hits),reference_gap_at_cut=any(h['operation'] in 'DN' for h in hits),status='mapped_cut' if hits else 'unresolved_cut_outside_aligned_blocks')

def original_context(seq,p):
 overlap=int(p.get('overlap_nt') or 0);shift=int(p.get('shift_nt') or 0);slide=int(p.get('slide_nt') or 0)
 left=int(p['a_qend'])-(shift if overlap else 0)+slide
 right=left if overlap else int(p['b_qstart'])+slide
 if not (0<=left<=right<=len(seq)):return None
 return seq[max(0,left-20):left]+'|'+seq[left:right].lower()+'|'+seq[right:right+20]

def reconcile_source(library,query,pairs,source):
 """Match actual arms using CIGAR, original-query interval and unnormalized endpoint.

 Stored exons may have been trimmed at both boundaries, so their first coordinate is
 deliberately NOT treated as an original BAM start. Original aligned_donor/acceptor
 plus the frozen overlap shift recover the untrimmed arm endpoint.
 """
 issues=[];alignments=source.get('alignments',[]);seq=source.get('full_forward_sequence');checks=[]
 if not source.get('source_pin'):issues.append('missing_source_provenance_pin')
 if not alignments:issues.append('missing_actual_alignments')
 identities=[(equivalent_alignment_key(a),bool(a['secondary']),bool(a['supplementary'])) for a in alignments]
 if len(identities)!=len(set(identities)):issues.append('duplicate_actual_alignment_identity')
 if sum(not a['secondary'] and not a['supplementary'] for a in alignments)!=1:issues.append('primary_missing_or_multiple')
 lengths={int(p['read_length']) for p in pairs}
 if len(lengths)!=1:issues.append('baseline_read_length_conflict')
 n=next(iter(lengths)) if len(lengths)==1 else None
 seqs={a['forward_sequence'].upper() for a in alignments if a.get('forward_sequence') and len(a['forward_sequence'])==n}
 if not seq or len(seq)!=n or seq.upper() not in seqs or len(seqs)!=1:issues.append('full_query_sequence_missing_or_conflicting')
 geometry=[];hardclip_checks=[]
 for a in alignments:
  g=alignment_geometry(a);geometry.append(g);hc=hard_clipped_source_slice(a,seq);hardclip_checks.append(hc)
  if hc['status'] not in ('not_hard_clipped','verified_full_sequence_slice'):issues.append(hc['status'])
  if a.get('query')!=query:issues.append('actual_record_query_mismatch')
  if g['full_query_length']!=n:issues.append('actual_cigar_query_length_mismatch')
 primary_records=[a for a in alignments if not a['secondary'] and not a['supplementary']]
 if len(primary_records)==1:
  tags=primary_records[0].get('tags',{});rule=source.get('barcode_rule',{'tags':['CB']});barcode=None
  for tag in rule.get('tags',['CB']):
   if tags.get(tag):barcode=tags[tag];break
  if barcode is None and rule.get('default_barcode') and rule.get('provenance_pin'):barcode=rule['default_barcode']
  for p in pairs:
   if p.get('barcode') and barcode is None:issues.append('primary_barcode_not_verified_missing_tag')
   elif barcode is not None and (barcode if '-' in barcode else barcode+'-1')!=p.get('barcode'):issues.append('primary_barcode_identity_mismatch')
   if p.get('umi') and 'XM' not in tags:issues.append('primary_UMI_not_verified_missing_tag')
 for p in pairs:
  if seq and original_context(seq,p)!=p.get('junction_read_context'):issues.append('scanner_context_replay_mismatch')
  for arm in ('a','b'):
   strand=p[arm+'_strand'];step=1 if strand=='+' else -1
   shift=int(p.get('shift_nt') or 0);overlap=int(p.get('overlap_nt') or 0)
   field='aligned_donor' if arm=='a' else 'aligned_acceptor'
   if p.get(field) in ('',None):issues.append('missing_original_aligned_endpoint');continue
   endpoint=int(p[field])+(shift*step if arm=='a' else -(overlap-shift)*step)
   matched=[]
   for i,(a,g) in enumerate(zip(alignments,geometry)):
    exons=reference_exons_1based(a)
    selected_exon=None if not exons else ((exons[-1] if strand=='+' else exons[0]) if arm=='a' else (exons[0] if strand=='+' else exons[-1]))
    actual_endpoint=(selected_exon[1] if strand=='+' else selected_exon[0]) if arm=='a' else (selected_exon[0] if strand=='+' else selected_exon[1]) if selected_exon else None
    if (not a['secondary'] and contig(a['chrom'])==contig(p['chrom']) and ('-' if a['reverse'] else '+')==strand and terminal_clip_signature(a['cigar'])==terminal_clip_signature(p[arm+'_cigar']) and (g['query_start'],g['query_end'])==(int(p[arm+'_qstart']),int(p[arm+'_qend'])) and (not a['supplementary'])==(str(p[arm+'_primary'])=='True') and actual_endpoint==endpoint):matched.append(i)
   overlap_cut=int(p['a_qend'])-shift if overlap else (int(p['a_qend']) if arm=='a' else int(p['b_qstart']))
   checks.append(dict(pair_index=p['pair_index'],arm=arm,unaltered_endpoint=endpoint,matching_actual_record_indices=matched,actual_pair_arm_mapqs=[alignments[i]['mapq'] for i in matched],terminal_clip_equivalence_used=[alignments[i]['cigar']!=p[arm+'_cigar'] for i in matched],baseline_pair_arm_mapq=p.get(arm+'_mapq'),baseline_split_query_cut0=overlap_cut,actual_cigar_cut_geometry=[query_cut_reference_context(alignments[i],overlap_cut) for i in matched]))
   if len(matched)!=1:issues.append('selected_arm_missing_or_ambiguous')
  for a in alignments:
   tags=a.get('tags',{});cb=tags.get('CB')
   if cb and (cb if '-' in cb else cb+'-1')!=p.get('barcode'):issues.append('barcode_identity_mismatch')
   if 'XM' in tags and tags['XM']!=p.get('umi',''):issues.append('UMI_identity_mismatch')
 missing=[];score_differences=[]
 expected=set()
 for a in alignments:
  for entry in a.get('tags',{}).get('SA','').rstrip(';').split(';'):
   if not entry:continue
   c,pos,s,cg,mq,nm=entry.split(',');expected.add((contig(c),int(pos)-1,s,cg,int(mq),int(nm)))
 for s in sorted(expected):
  found=[a for a in alignments if equivalent_alignment_key(a)==(s[0],s[1],s[2],terminal_clip_signature(s[3]))]
  if not found:missing.append(s)
  for a in found:
   if (a['mapq'],a.get('tags',{}).get('NM'))!=s[4:]:score_differences.append(dict(identity=s[:4],SA_mapq=s[4],SA_NM=s[5],actual_mapq=a['mapq'],actual_NM=a.get('tags',{}).get('NM')))
 if missing:issues.append('missing_actual_SA_alignment')
 descriptor_issue_names={'missing_actual_SA_alignment','scanner_context_replay_mismatch'}
 arm_issue_names={'missing_original_aligned_endpoint','selected_arm_missing_or_ambiguous'}
 actual_issues=sorted(set(issues)-descriptor_issue_names-arm_issue_names)
 descriptor_issues=sorted(set(issues)&descriptor_issue_names)
 arm_issues=sorted(set(issues)&arm_issue_names)
 identity_axes=dict(actual_query_identity={'status':'verified_actual_query_full_sequence' if not actual_issues else 'unresolved_actual_query_identity','issues':actual_issues},baseline_descriptor_parity={'status':'baseline_descriptors_concordant' if not descriptor_issues else 'unresolved_baseline_descriptor_parity','issues':descriptor_issues,'score_discrepancies_nonblocking':len(score_differences)},selected_arm_matching={'status':'selected_actual_arms_verified' if not arm_issues else 'unresolved_selected_arm_matching','issues':arm_issues})
 return dict(**identity_axes,status='resolved_source_identity' if not issues else 'unresolved_source_identity',issues=sorted(set(issues)),selected_arm_checks=checks,missing_SA=missing,SA_score_discrepancies=score_differences,all_actual_alignment_geometry=geometry,hardclip_sequence_checks=hardclip_checks,limitation='SA-linked retrieval is not an exhaustive search for unreported alignments.')

def parse_paf(line):
 a=line.rstrip('\n').split('\t')
 if len(a)<12:raise ValueError('Malformed PAF')
 tags={}
 for t in a[12:]:
  k,typ,v=t.split(':',2)
  if k in tags:raise ValueError('Duplicate PAF tag')
  tags[k]=int(v) if typ=='i' else float(v) if typ=='f' else v
 r=dict(query=a[0],query_length=int(a[1]),query_start=int(a[2]),query_end=int(a[3]),strand=a[4],target=a[5],target_length=int(a[6]),target_start=int(a[7]),target_end=int(a[8]),matching_bases=int(a[9]),block_length=int(a[10]),mapq=int(a[11]),tags=tags)
 if r['strand'] not in ('+','-') or not 0<=r['query_start']<r['query_end']<=r['query_length'] or not 0<=r['target_start']<r['target_end']<=r['target_length']:raise ValueError('Invalid PAF intervals')
 return r

def cigar_boundary_evidence(paf,boundaries):
 """Every template boundary spanned by the alignment, without an anchor cutoff.

 M means aligned bases, NOT exact matches. I/D may explain original positive query
 gaps. A deletion spanning the boundary leaves precise breakpoint unresolved.
 Report counts/operations; no arbitrary pass/fail threshold is applied.
 """
 if 'cg' not in paf['tags']:return dict(status='unresolved_missing_alignment_cigar',observations=[])
 q=paf['query_start'] if paf['strand']=='+' else paf['query_length']-paf['query_end'];t=paf['target_start'];blocks=[]
 for k,o in cigar_ops(paf['tags']['cg']):
  if o in 'SH':raise ValueError('PAF cg must describe aligned interval without clipping')
  tq=t+k if o in 'MDN=X' else t;qq=q+k if o in 'MI=X' else q
  blocks.append(dict(op=o,length=k,target_start=t,target_end=tq,query_start=q,query_end=qq));t,q=tq,qq
 if t!=paf['target_end'] or q!=(paf['query_end'] if paf['strand']=='+' else paf['query_length']-paf['query_start']):raise ValueError('PAF CIGAR span mismatch')
 out=[]
 for boundary in boundaries:
  b=int(boundary)
  if not paf['target_start']<b<paf['target_end']:continue
  left=sum(max(0,min(b,z['target_end'])-z['target_start']) for z in blocks if z['op'] in 'M=X' and z['target_start']<b)
  right=sum(max(0,z['target_end']-max(b,z['target_start'])) for z in blocks if z['op'] in 'M=X' and z['target_end']>b)
  crossing=[z for z in blocks if z['target_start']<b<z['target_end']]
  deleted=any(z['op'] in 'DN' for z in crossing)
  insertion=[z for z in blocks if z['op']=='I' and z['target_start']==b]
  contiguous_left=max([b-z['target_start'] for z in blocks if z['op'] in 'M=X' and z['target_start']<b<=z['target_end']] or [0])
  contiguous_right=max([z['target_end']-b for z in blocks if z['op'] in 'M=X' and z['target_start']<=b<z['target_end']] or [0])
  qpositions=[]
  for z in blocks:
   if z['op'] in 'M=X' and z['target_start']<=b<=z['target_end']:qpositions.append(z['query_start']+b-z['target_start'])
   elif z['op'] in 'DN' and z['target_start']<b<z['target_end']:qpositions.append(z['query_start'])
  if paf['strand']=='-':qpositions=[paf['query_length']-v for v in qpositions]
  out.append(dict(target_boundary0=b,query_boundary_positions0=sorted(set(qpositions)),aligned_bases_left=left,aligned_bases_right=right,adjacent_aligned_bases_left=contiguous_left,adjacent_aligned_bases_right=contiguous_right,insertion_at_boundary_nt=sum(z['length'] for z in insertion),reference_deletion_or_skip_spans_boundary=deleted,aligned_bases_both_sides=left>0 and right>0,adjacent_alignment_reaches_both_sides=contiguous_left>0 and contiguous_right>0,exact_match_claim=False))
 return dict(status='measured_cigar_geometry',observations=out,blocks=blocks)

def query_path(alignments,sequence):
 primary=[];secondary=[]
 for i,a in enumerate(alignments):
  row=dict(actual_record_index=i,alignment_key=alignment_key(a),**alignment_geometry(a))
  (secondary if a['secondary'] else primary).append(row)
 primary.sort(key=lambda x:(x['query_start'],x['query_end'],x['alignment_key']))
 adjacent=[]
 for a,b in zip(primary,primary[1:]):
  gap=b['query_start']-a['query_end'];s=sequence[a['query_end']:b['query_start']] if sequence and gap>0 else None
  adjacent.append(dict(left_record=a['actual_record_index'],right_record=b['actual_record_index'],query_gap_nt=gap,query_overlap_nt=max(0,-gap),gap_sequence=s))
 covered=sorted((a['query_start'],a['query_end']) for a in primary);union=[]
 for s,e in covered:
  if union and s<=union[-1][1]:union[-1][1]=max(e,union[-1][1])
  else:union.append([s,e])
 return dict(primary_supplementary_path=primary,secondary_alternatives=secondary,adjacent_pair_relations=adjacent,union_aligned_query_intervals=union,union_aligned_query_bases=sum(e-s for s,e in union),query_length=len(sequence) if sequence else None,limitation='Overlapping alignments retained; sorted adjacency is not a uniquely resolved molecular path.')

def gap_exact_evidence(gap,reference):
 """Longest exact local match at every gap start, in each orientation.

 No length cutoff. Repeated placements are counted; first placement is explicitly
 representative only. Original sequence/reference remain available for enumeration.
 This is a sequence diagnostic, not an alignment/classification substitute.
 """
 if not gap:return dict(status='no_positive_query_gap',matches=[])
 if reference is None:return dict(status='unresolved_missing_local_reference',matches=[])
 gap=gap.upper();reference=reference.upper();matches=[]
 for orientation,ref in [('+',reference),('-',revcomp(reference))]:
  for i in range(len(gap)):
   stop=i
   while stop<len(gap) and gap[stop] in 'ACGT':stop+=1
   lo,hi=0,stop-i
   while lo<hi:
    k=(lo+hi+1)//2
    if gap[i:i+k] in ref:lo=k
    else:hi=k-1
   if not lo:continue
   pattern=gap[i:i+lo];pos=ref.find(pattern);first=pos;count=0
   while pos>=0:count+=1;pos=ref.find(pattern,pos+1)
   rs,re_=first,first+lo
   if orientation=='-':rs,re_=len(ref)-re_,len(ref)-rs
   matches.append(dict(query_start=i,query_end=i+lo,reference_start=rs,reference_end=re_,orientation=orientation,length=lo,exact_placement_count=count,reference_position_is_representative=True))
 return dict(status='measured_exact_gap_matches',matches=matches,longest_forward=max([m['length'] for m in matches if m['orientation']=='+'] or [0]),longest_reverse=max([m['length'] for m in matches if m['orientation']=='-'] or [0]),limitation='Longest exact match per gap start and orientation; no minimum match length. Short repeats and palindromes are ambiguous; these matches do not establish mechanism or physical circularity.')

def validate_model(m,event):
 required={'model_id','kind','sequence_sha256','reference_sha256','coordinate_blocks','alignment_parameters','paf_records','boundaries','model_basis'}
 if not required.issubset(m):raise ValueError('Missing model provenance fields')
 if m['kind'] not in ('exon_chain','genomic_interval','alternative_genomic','alternative_sequence'):raise ValueError('Unknown model kind')
 if not all(m[x] for x in ('model_id','alignment_parameters','model_basis')):raise ValueError('Empty model provenance')
 for name in ('sequence_sha256','reference_sha256'):
  if not re.fullmatch('[a-f0-9]{64}',m[name]):raise ValueError('Malformed model sequence/reference SHA256')
 if len(m['boundaries'])!=len(set(m['boundaries'])):raise ValueError('Duplicate model boundary')
 if m['kind'] in ('alternative_genomic','alternative_sequence'):
  if m['boundaries']:raise ValueError('Ordinary genomic model cannot declare circular-copy boundaries')
  return
 blocks=m['coordinate_blocks']
 if not blocks or (m['kind']=='genomic_interval' and len(blocks)!=1):raise ValueError('Missing/invalid candidate coordinate blocks')
 strand=event['donor_strand']
 if event['acceptor_strand']!=strand:raise ValueError('Candidate model cannot silently project discordant strands')
 for b in blocks:
  if not (isinstance(b['start0'],int) and isinstance(b['end0'],int) and 0<=b['start0']<b['end0']):raise ValueError('Invalid coordinate block')
  if contig(b['chrom'])!=contig(event['chrom']) or b['strand']!=strand:raise ValueError('Candidate model contig/strand mismatch')
 first=blocks[0]['start0']+1 if strand=='+' else blocks[0]['end0'];last=blocks[-1]['end0'] if strand=='+' else blocks[-1]['start0']+1
 if (first,last)!=(int(event['acceptor']),int(event['donor'])):raise ValueError('Model endpoints differ from candidate')
 for a,b in zip(blocks,blocks[1:]):
  if (a['end0']>b['start0'] if strand=='+' else b['end0']>a['start0']):raise ValueError('Candidate model blocks not disjoint transcription order')
 length=sum(b['end0']-b['start0'] for b in blocks)
 if m.get('copy_length')!=length or not isinstance(m.get('copy_number'),int) or m['copy_number']<2:raise ValueError('Missing/inconsistent explicit tandem copy number')
 if m.get('sequence_length')!=length*m['copy_number']:raise ValueError('Template length/copy advantage not explicit')
 if sorted(m['boundaries'])!=list(range(length,m['sequence_length'],length)):raise ValueError('Missing/invented template boundary')
 if 'sequence' in m and (len(m['sequence'])!=m['sequence_length'] or hashlib.sha256(m['sequence'].encode()).hexdigest()!=m['sequence_sha256']):raise ValueError('Template sequence hash/length differs')


def pin_verified(pin,verified_artifacts):
 return isinstance(pin,dict) and isinstance(pin.get('artifact'),str) and re.fullmatch('[a-f0-9]{64}',str(pin.get('sha256',''))) is not None and verified_artifacts.get(pin['artifact'])==pin['sha256']

class DiskOutcomeLookup:
 def __init__(self,connection,model):self.connection=connection;self.model=model
 def get(self,key):
  r=self.connection.execute('SELECT payload FROM outcomes WHERE model=? AND query=?',(self.model,key)).fetchone()
  return json.loads(r[0]) if r else None

class VerifiedArtifacts(dict):
 def __init__(self):super().__init__();self.documents={};self._scratch=tempfile.TemporaryDirectory(prefix='circrna-outcomes-');self._db=sqlite3.connect(str(pathlib.Path(self._scratch.name)/'outcomes.sqlite'));self._db.execute('CREATE TABLE outcomes(model TEXT,query TEXT,payload TEXT,PRIMARY KEY(model,query)) WITHOUT ROWID')
 def close(self):self._db.close();self._scratch.cleanup()


def verify_artifact_manifest(path):
 """Hash declared extraction/reference/template/execution artifacts once, not per read.

 BAM full-byte identity versus metadata-only concordance must be stated by the
 upstream extraction manifest; this gate establishes these artifact bytes only.
 """
 manifest=json.loads(pathlib.Path(path).read_text());items=manifest['artifacts'];seen=set();verified=VerifiedArtifacts()
 for item in items:
  aid=item['artifact']
  if aid in seen:raise ValueError('Duplicate artifact identity')
  seen.add(aid)
  if digest(item['path'])!=item['sha256']:raise ValueError('Artifact hash differs '+aid)
  verified[aid]=item['sha256']
  if item.get('kind')=='execution_receipt':
   document=json.loads(pathlib.Path(item['path']).read_text());model_lookup={}
   loaded={}
   for run in document['models']:
    if run['model_id'] in model_lookup:raise ValueError('Duplicate execution model')
    if 'outcomes_file' in run:
     op=pathlib.Path(item['path']).parent/run['outcomes_file']
     if digest(op)!=run['outcomes_sha256']:raise ValueError('External outcome artifact differs')
     token=item['artifact']+'|'+run.get('executed_model_id',run['model_id'])
     if token not in loaded:
      count=0
      with gzip.open(op,'rt') as handle:
       for line in handle:
        o=json.loads(line);verified._db.execute('INSERT INTO outcomes VALUES (?,?,?)',(token,o['id'],line));count+=1
      if count!=run['outcomes_count']:raise ValueError('External outcome count differs')
      verified._db.commit();loaded[token]=run['outcomes_sha256']
     elif loaded[token]!=run['outcomes_sha256']:raise ValueError('Alias outcomes differ')
     outcomes=DiskOutcomeLookup(verified._db,token)
    else:
     outcomes={o['id']:o for o in run['outcomes']}
     if len(outcomes)!=len(run['outcomes']):raise ValueError('Duplicate execution query outcome')
    run['_outcomes_by_id']=outcomes;model_lookup[run['model_id']]=run
   document['_models_by_id']=model_lookup;verified.documents[aid]=document
 return verified

def validate_record(record,verified_artifacts=None,expected_model_ids=None,gap_cache=None):
 verified_artifacts=verified_artifacts or {}
 if record.get('schema_version')!=2:raise ValueError('Schema must be 2')
 e=record['event'];event=event_key(e)
 if e['event_id']!=event:raise ValueError('Event identity differs from coordinates')
 lib,query=record['library'],record['query'];pairs=record['baseline_pairs']
 if not lib or not query or not pairs:raise ValueError('Missing baseline source identity')
 for p in pairs:
  if p['library']!=lib or p['query_name']!=query or pair_event(p)!=event:raise ValueError('Baseline event/library/query mismatch')
 source=record.get('source',{});identity=reconcile_source(lib,query,pairs,source);seq=source.get('full_forward_sequence');path=query_path(source.get('alignments',[]),seq)
 exact=set();shifted=set()
 for p in pairs:
  positions=(int(p['a_qend']),int(p['b_qstart']))
  # Frozen normalized event occurrences and original aligned coordinates differ.
  exact.add(positions)
  shifted.add((int(p.get('aligned_donor') or p['donor']),int(p.get('aligned_acceptor') or p['acceptor']),positions))
 models=[];unresolved=list(identity['issues'])
 if not pin_verified(source.get('source_pin'),verified_artifacts):
  identity['issues'].append('source_artifact_not_verified');identity['status']='unresolved_source_identity';unresolved.append('source_artifact_not_verified')
 model_ids=[m['model_id'] for m in record.get('models',[])]
 if len(model_ids)!=len(set(model_ids)):raise ValueError('Duplicate model identity')
 if expected_model_ids is None:unresolved.append('expected_model_roster_not_pinned');expected_model_ids=[]
 if len(expected_model_ids)!=len(set(expected_model_ids)):raise ValueError('Duplicate expected model identity')
 if set(model_ids)-set(expected_model_ids):raise ValueError('Unexpected model identity')
 for missing in sorted(set(expected_model_ids)-set(model_ids)):
  models.append(dict(model_id=missing,status='not_evaluated',mapping_count=None,reason='expected_model_absent_retained'));unresolved.append('expected_model_not_evaluated:'+missing)
 for m in record.get('models',[]):
  validate_model(m,e)
  execution=m.get('execution',{});execution_status=execution.get('status','not_evaluated')
  if execution_status not in ('completed','failed','not_evaluated'):raise ValueError('Unknown model execution status')
  verified_execution=(execution_status=='completed' and pin_verified(execution.get('manifest_pin'),verified_artifacts) and pin_verified(execution.get('paf_pin'),verified_artifacts) and pin_verified(m.get('template_pin'),verified_artifacts))
  if not verified_execution:
   if m['paf_records']:raise ValueError('Alignment records supplied without verified successful execution')
   models.append(dict(model_id=m['model_id'],status=execution_status if execution_status!='completed' else 'unresolved_execution_provenance',mapping_count=None,model_provenance=m));unresolved.append('model_not_successfully_evaluated:'+m['model_id']);continue
  receipt=getattr(verified_artifacts,'documents',{}).get(execution['manifest_pin']['artifact'])
  if receipt is None:raise ValueError('Verified execution receipt contents unavailable')
  run=receipt['_models_by_id'].get(m['model_id'])
  if run is None:raise ValueError('Execution model missing')
  if run['status']!='completed' or run['returncode']!=0 or run['paf_sha256']!=execution['paf_pin']['sha256'] or run['template_sha256']!=m['template_pin']['sha256'] or run['model_sequence_sha256']!=m['sequence_sha256'] or run['command']!=m['alignment_parameters']:raise ValueError('Execution receipt does not bind model/parameters/artifacts')
  if receipt.get('reference_sha256')!=m['reference_sha256']:raise ValueError('Execution reference identity differs')
  qid=record.get('alignment_query_id',query)
  if seq and receipt['query_sequence_sha256'].get(qid)!=hashlib.sha256(seq.encode()).hexdigest():raise ValueError('Execution query sequence differs from source')
  expected_mapping=run['_outcomes_by_id'].get(qid)
  if expected_mapping is None:raise ValueError('Execution query outcome missing')
  if expected_mapping['mapping_count'] is None:
   if m['paf_records']:raise ValueError('Mappings exist for unavailable query')
   models.append(dict(model_id=m['model_id'],status=expected_mapping['status'],mapping_count=None));unresolved.append('model_query_not_evaluated:'+m['model_id']);continue
  if expected_mapping['mapping_count']!=len(m['paf_records']):raise ValueError('Execution query mapping count differs')
  observations=[];parsed_records_hash=hashlib.sha256()
  for a in m['paf_records']:
   p=parse_paf(a) if isinstance(a,str) else a
   parsed_records_hash.update((json.dumps(p,sort_keys=True,separators=(',',':'))+'\n').encode())
   if p['query']!=record.get('alignment_query_id',query) or (seq and p['query_length']!=len(seq)):raise ValueError('PAF query/length mismatch')
   if m['kind'] in ('exon_chain','genomic_interval') and (p['target']!=m.get('executed_model_id',m['model_id']) or p['target_length']!=m['sequence_length']):raise ValueError('PAF target model mismatch')
   v=cigar_boundary_evidence(p,m['boundaries']);observations.append(dict(alignment=p,boundary_evidence=v))
  if parsed_records_hash.hexdigest()!=expected_mapping['paf_records_sha256']:raise ValueError('PAF subset content differs from successful execution')
  models.append(dict(model_id=m['model_id'],kind=m['kind'],mapping_count=len(observations),status='mapped' if observations else 'unmapped_retained',alignments=observations,model_provenance={k:v for k,v in m.items() if k!='paf_records'}))
 if not expected_model_ids:unresolved.append('alternative_and_junction_models_not_evaluated')
 local=record.get('local_reference');gap_results=[]
 if local is not None:
  if contig(local['chrom'])!=contig(e['chrom']):raise ValueError('Local reference contig differs from event')
  if not pin_verified(local.get('reference_manifest_pin'),verified_artifacts):raise ValueError('Local reference artifact not verified')
  if not all(local.get(k) for k in ('reference_sha256','sequence_sha256','sequence','chrom')) or not isinstance(local.get('start0'),int) or local['start0']<0:raise ValueError('Missing local reference provenance')
  if hashlib.sha256(local['sequence'].encode()).hexdigest()!=local['sequence_sha256']:raise ValueError('Local reference sequence hash differs')
 for p in pairs:
  start,end=int(p['a_qend']),int(p['b_qstart']);gap=seq[start:end] if seq and end>start else ''
  if not seq:ge=dict(status='unresolved_missing_query_sequence',matches=[])
  elif gap_cache is None:ge=gap_exact_evidence(gap,local['sequence'] if local else None)
  else:
   cache_id=hashlib.sha256((gap+'\0'+(local['sequence_sha256'] if local else 'missing_reference')).encode()).hexdigest()
   if cache_id not in gap_cache:gap_cache[cache_id]=dict(gap_sequence=gap,reference_sequence_sha256=local['sequence_sha256'] if local else None,**gap_exact_evidence(gap,local['sequence'] if local else None))
   cached=gap_cache[cache_id];ge=dict(status=cached['status'],shared_gap_evidence_id=cache_id,longest_forward=cached.get('longest_forward'),longest_reverse=cached.get('longest_reverse'))
  gap_results.append(dict(pair_index=p['pair_index'],baseline_query_gap=int(p['query_gap']),query_interval=[start,end],sequence=gap if seq else None,evidence=ge))
 return dict(schema_version=2,event=e,library=lib,query=query,baseline=dict(pair_rows=len(pairs),status_pair_counts=dict(collections.Counter(p['status'] for p in pairs)),any_kept=any(p['status']=='kept' for p in pairs),original_pairs=pairs),revised_evidence=dict(source_identity=identity,full_query_path=path,gap_sequence_evidence=gap_results,candidate_specific_normalized_pair_occurrences=len(exact),original_aligned_coordinate_occurrences=[dict(donor=d,acceptor=a,query_positions=list(q)) for d,a,q in sorted(shifted)],reference_pair_overlap_nt=[int(p.get('pair_overlap_nt') or 0) for p in pairs],reported_crossings=[int(p.get('crossings') or 0) for p in pairs],model_alignments=models,unresolved_reasons=sorted(set(unresolved))),source_evidence=source,status='ADDITIVE_EVIDENCE_RETAINED',physical_circularity='not_established',limitation='No feature/sample exclusion or statistical replacement. Model preference and repeated junction sequence cannot distinguish circular-copy from linear-concatemer origins.')

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--expected-roster',required=True);ap.add_argument('--output',required=True);ap.add_argument('--artifact-manifest',required=True);args=ap.parse_args()
 src=pathlib.Path(args.input);roster=pathlib.Path(args.expected_roster);out=pathlib.Path(args.output)
 if out.exists():raise FileExistsError(out)
 verified=verify_artifact_manifest(args.artifact_manifest)
 expected=json.loads(roster.read_text());
 if digest(args.artifact_manifest)!=expected['artifact_manifest_sha256']:raise ValueError('Artifact manifest differs from independent roster')
 keys=[tuple(x) for x in expected['event_library_query_keys']];keyset=set(keys)
 if len(keys)!=len(set(keys)):raise ValueError('Duplicate expected roster identity')
 if digest(src)!=expected['input_sha256']:raise ValueError('Input hash differs')
 pending=out.with_suffix(out.suffix+'.pending');seen=set();totals=collections.Counter()
 if pending.exists():raise FileExistsError(pending)
 with pending.open('x') as f:
  for line in src.open():
   r=json.loads(line);key=(r['event']['event_id'],r['library'],r['query'])
   if key in seen or key not in keyset:raise ValueError('Duplicate/unexpected input identity')
   seen.add(key);v=validate_record(r,verified,expected['expected_model_ids_by_event'][key[0]]);f.write(json.dumps(v,separators=(',',':'))+'\n');totals['source_queries']+=1;totals['baseline_pair_rows']+=len(r['baseline_pairs'])
 if seen!=keyset:raise ValueError('Missing expected identities; partial pending file retained')
 if digest(src)!=expected['input_sha256']:raise ValueError('Input changed during evaluation')
 pending.rename(out)
 out.with_suffix(out.suffix+'.complete.json').write_text(json.dumps(dict(status='COMPLETE_EXACT_ROSTER_NO_EXCLUSIONS',input_sha256=digest(src),expected_roster_sha256=digest(roster),output_sha256=digest(out),code_sha256=digest(__file__),totals=dict(totals)),indent=2)+'\n')
if __name__=='__main__':main()
