"""Collapse-only evidence-preserving input path for the existing WTA collapse.

An evidence class is formed BEFORE exact aggregation. Unprotected reads retain
their original key. Protected classes carry actual source junctions/coverage;
variable transcript endpoints do not become independent classes. Policy values are explicit; this module does not replace collapse_gene_jc
or change upstream annotation.
"""
from __future__ import annotations
from collections import Counter, defaultdict
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json


@dataclass(frozen=True)
class Policy:
    splice_slop: int
    minimum_continuation: int
    minimum_junction_support: int

    def __post_init__(self):
        if self.splice_slop < 0 or self.minimum_continuation < 1 or self.minimum_junction_support < 1:
            raise ValueError('Explicit valid splice tolerance, continuation and recurrence settings are required')


def merge_exons(exons):
    result = []
    for s,e in sorted(exons):
        if s > e:
            raise ValueError('Source exons must use ascending genomic coordinates')
        if result and s <= result[-1][1]+1:
            result[-1] = (result[-1][0],max(e,result[-1][1]))
        else:
            result.append((s,e))
    if not result:
        raise ValueError('A source exon record is required for every input read')
    return result


def splice_spans_compatible(spans, parent_exons, splice_slop):
    """A parent cannot splice out sequence observed continuously across a trusted gap.

    A shorter terminal endpoint is not a splice contradiction. Spans are selected
    before aggregation only where a source exon covers both flanks of an actual
    reference or recurrent source junction. They are not inferred intron labels.
    """
    if not spans:
        return True
    blocks = merge_exons(parent_exons)
    gaps = [(a[1],b[0]) for a,b in zip(blocks,blocks[1:])]
    # Require the actual target splice, not merely an overlapping interval.
    # A short source fragment spanning a different nested gap does not prove
    # both flanks of this parent gap were observed. Partial terminal coverage
    # remains governed by the existing minimum-continuation policy.
    return not any(max(abs(left-a),abs(right-b)) <= splice_slop
                   for left,right in spans for a,b in gaps)


def terminal_parent_gap_coverage(source_exons, parent_exons):
    """Measure terminal continuation into true parent gaps, including E labels."""
    source_exons=merge_exons(source_exons)
    parent_exons=merge_exons(parent_exons)
    result=[]
    for i in sorted({0,len(source_exons)-1}):
        s,e=source_exons[i]
        for left,right in zip(parent_exons,parent_exons[1:]):
            a,b=left[1]+1,right[0]-1
            lo,hi=max(s,a),min(e,b)
            if lo<=hi and (s<a<=e or s<=b<e):
                result.append(dict(parent_gap=(a,b),covered=(lo,hi),
                    measured_nt=hi-lo+1,terminal_block=(s,e)))
    return result


def terminal_events(exons, introns, parent_models=None):
    """Measured continuous intronic coverage and wholly intronic terminal exons.

    introns contains (region_id, inclusive_start, inclusive_end). Coordinates
    are genomic, so both strands and either transcript end use the same checks.
    Only terminal exon blocks are considered; interior events are unchanged.
    """
    exons = merge_exons(exons)
    result = []
    for i in sorted({0,len(exons)-1}):
        start,end = exons[i]
        for region,lo,hi in introns:
            if lo <= start <= end <= hi and len(exons)>1:
                junction = (end,exons[i+1][0]) if i==0 else (exons[i-1][1],start)
                result.append(('exon',region,junction,(start,end)))
            elif (start<lo<=end or start<=hi<end):
                a,b = max(start,lo),min(end,hi)
                if a<=b:
                    result.append(('continuation',region,None,(a,b)))
    if parent_models:
        seen=set()
        for parent_exons in parent_models.values():
            for event in terminal_parent_gap_coverage(exons,parent_exons):
                key=event['parent_gap'],event['terminal_block']
                if key not in seen:
                    seen.add(key)
                    a,b=event['parent_gap']
                    result.append(('parent_gap_continuation',f'gap:{a}-{b}',
                                   event['terminal_block'],event['covered']))
    return result


def terminal_u_events(raw_structure, exons, regions, gene, splice_slop):
    """Measure terminal U splice sites outside the annotated gene body.

    U labels come from the original raw read annotation, not fabricated region
    bounds. Require its coordinate to identify the actual terminal splice site;
    an RNA endpoint alone is not splice evidence. Recurrence is selected later
    by prepare_gene using the same molecule and tolerance policy as I exons.
    """
    blocks = merge_exons(exons)
    bounds = [v for k, v in regions.items() if k.startswith(('E', 'I'))]
    if len(blocks) < 2 or not bounds:
        return []
    lo = min(v[0] for v in bounds)
    hi = max(v[1] for v in bounds)
    labels = []
    for token in raw_structure.split('|'):
        if ':' in token and token.rsplit(':', 1)[0] != gene:
            continue
        label = token.rsplit(':', 1)[-1]
        if '_' not in label:
            continue
        base, coordinate = label.rsplit('_', 1)
        if base.startswith('U'):
            labels.append((base, int(coordinate)))
    events = []
    for i in (0, len(blocks)-1):
        start, end = blocks[i]
        if not (end < lo or start > hi):
            continue
        junction = (end, blocks[i+1][0]) if i == 0 else (blocks[i-1][1], start)
        boundary = end if i == 0 else start
        for label in sorted({base for base, coordinate in labels
                             if abs(coordinate-boundary) <= splice_slop}):
            events.append(('exon', label, junction, (start, end)))
    return events


def junction_anchors(records, slop):
    """Support-first bounded grouping of equivalent sites; no transitive chaining.

    A point belongs to an anchor only when BOTH endpoints are within slop of
    that observed anchor. Every source molecule is counted once per site.
    """
    molecules = defaultdict(set)
    for region,junction,source in records:
        molecules[region,tuple(junction)].add(source)
    anchor_for = {}
    groups = defaultdict(set)
    ordered = sorted(molecules,key=lambda p:(-len(molecules[p]),p))
    anchors = defaultdict(list)
    for region,junction in ordered:
        compatible = [a for a in anchors[region]
                      if all(abs(x-y)<=slop for x,y in zip(junction,a))]
        anchor = compatible[0] if compatible else junction
        if not compatible:
            anchors[region].append(anchor)
        anchor_for[region,junction] = anchor
        groups[region,anchor].update(molecules[region,junction])
    return anchor_for,{k:len(v) for k,v in groups.items()}


def prepare_gene(rows, introns, policy, *, annotation_corrector=None, parent_models=None,
                 required_identities=None, terminal_event_provider=None,
                 protect_splice_geometry=False):
    """Return complete read -> candidate mappings and pooled candidate evidence.

    Each row must contain library, gene, molecule, legacy_structure and actual
    source_exons. The caller supplies one complete gene roster across all input
    libraries. No row is filtered, even if baseline abundance later removes its
    representative. An internal class suffix is separate from the annotation
    structure used by WTA; consumers must use annotation_structure below.
    """
    events = []
    junction_records = []
    identity = set()
    genes = set()
    if required_identities is not None:
        expected=list(required_identities)
        if len(set(expected))!=len(expected):
            raise ValueError('Duplicate identity in the independent source roster')
        observed=[(r['library'],r['gene'],r['molecule']) for r in rows]
        if set(observed)!=set(expected) or len(observed)!=len(expected):
            raise ValueError('Complete source identity roster must match before evidence preparation')
    for row in rows:
        source = row['library'],row['gene'],row['molecule']
        if source in identity:
            raise ValueError(f'Duplicate source read identity: {source}')
        identity.add(source)
        genes.add(row['gene'])
        ev = (terminal_events(row['source_exons'],introns,parent_models) if terminal_event_provider is None
              else terminal_event_provider(row))
        events.append(ev)
        support_source=tuple(row.get('support_source',source))
        junction_records.extend((region,junction,support_source) for kind,region,junction,coverage in ev if kind in ('exon','splice'))
    if len(genes)!=1:
        raise ValueError('prepare_gene requires one complete gene and its reference regions')
    anchor_for,support = junction_anchors(junction_records,policy.splice_slop)
    # Selection is complete. Do not retain this temporary per-observation
    # index while allocating the final source-to-class mappings.
    del junction_records
    trusted = set()
    if protect_splice_geometry:
        for exons in (parent_models or {}).values():
            blocks=merge_exons(exons)
            trusted.update((a[1],b[0]) for a,b in zip(blocks,blocks[1:]))
        trusted.update(anchor for (region,anchor),count in support.items()
                       if region=='physical_splice' and count>=policy.minimum_junction_support)
    trusted = sorted(trusted)
    trusted_left = [a for a,b in trusted]
    from ..isoform_collapse_utils import intronic_edge_compatible
    candidates = {}
    mapping = {}
    prepared_rows = []
    first_known_targets = {}
    mixed_known_targets = set()
    @lru_cache(maxsize=4096)
    def evidence_key(base, signature):
        # JSON encodes tuples/lists identically. Cache only the deterministic
        # class name, never observations, counts or source membership.
        digest = hashlib.sha256(json.dumps(signature).encode()).hexdigest()[:20]
        return base+'~intronic='+digest
    @lru_cache(maxsize=65536)
    def prepared_evidence(junctions, coverage, terminal_blocks, splice_spans, selected):
        # Share equal immutable values, never RNA identities or observations.
        # Highly expressed genes can otherwise allocate the same selected
        # reference/recurrent splice evidence millions of times.
        return junctions, coverage, terminal_blocks, splice_spans, selected
    @lru_cache(maxsize=65536)
    def compatible_known_references(junctions, coverage, splice_spans):
        # Bounded site anchors can pool RNAs with different feasible references:
        # sites 100 and 103 are within 3 nt, but only 103 fits reference 106.
        # Partition before pooling so known preference remains applicable to
        # each compatible RNA. This does not rank or select any parent.
        return tuple(sorted(name for name,exons in (parent_models or {}).items()
                            if intronic_edge_compatible(junctions,coverage,exons,
                                                       splice_slop=policy.splice_slop)
                            and splice_spans_compatible(splice_spans,exons,policy.splice_slop)))
    @lru_cache(maxsize=65536)
    def span_evidence(spans):
        return (spans,tuple(('splice_span',left,right) for left,right in spans),
                tuple(('splice_span','physical_splice',None,(left,right)) for left,right in spans))
    @lru_cache(maxsize=65536)
    def source_splice_spans(exons):
        # Trusted gaps and policy are fixed for this complete gene. Reuse only
        # the pure geometric result; every RNA still follows the full mapping
        # and count path. Intern equal outputs across different RNA endpoints.
        spans=set()
        for start,end in merge_exons(exons):
            lo=bisect_left(trusted_left,start+policy.splice_slop)
            hi=bisect_right(trusted_left,end-policy.splice_slop)
            for left,right in trusted[lo:hi]:
                if right<=end-policy.splice_slop:spans.add((left,right))
        return span_evidence(tuple(sorted(spans)))
    for row,ev in zip(rows,events):
        protected_junctions = []
        protected_coverage = []
        protected_terminal_blocks = set()
        protected_splice_spans = set()
        selected_events = []
        signature = []
        for kind,region,junction,coverage in ev:
            if kind in ('exon','splice'):
                anchor = anchor_for[region,tuple(junction)]
                if support[region,anchor] >= policy.minimum_junction_support:
                    protected_junctions.append(tuple(junction))
                    signature.append((kind,region,*anchor))
                    selected_events.append((kind,region,junction,coverage))
            elif kind=='parent_gap_continuation' and coverage[1]-coverage[0]+1>=policy.minimum_continuation:
                protected_terminal_blocks.add(tuple(junction))
                protected_coverage.append(tuple(coverage))
                signature.append(('parent_gap_continuation',region))
                selected_events.append((kind,region,junction,coverage))
            elif kind=='continuation' and coverage[1]-coverage[0]+1 >= policy.minimum_continuation:
                protected_coverage.append(tuple(coverage))
                # End extent is evidence, not a new isoform identity per offset.
                signature.append(('continuation',region))
                selected_events.append((kind,region,junction,coverage))
        if protect_splice_geometry:
            protected_splice_spans,span_signature,span_events=source_splice_spans(
                tuple(tuple(exon) for exon in row['source_exons']))
            signature.extend(span_signature)
            selected_events.extend(span_events)
        selected_value=tuple((kind,region,tuple(junction) if junction is not None else None,
                              tuple(coverage)) for kind,region,junction,coverage in selected_events)
        prepared=prepared_evidence(tuple(protected_junctions),tuple(protected_coverage),
                                   tuple(protected_terminal_blocks),tuple(protected_splice_spans),selected_value)
        known_targets=(compatible_known_references(prepared[0],prepared[1],prepared[3])
                       if protect_splice_geometry and signature and parent_models else None)
        base = row['legacy_structure']
        # Evidence must not leave the corrected candidate annotated as purely
        # exonic when the original terminal-strip bug removed its I region.
        regions = {t.rsplit(':',1)[-1].split('_')[0] for t in base.split('|')
                   if ':' not in t or t.rsplit(':',1)[0]==row['gene']}
        missing = {s[1] for s in signature if s[0]=='continuation'}-regions
        if missing and annotation_corrector is None:
            raise ValueError(f'Protected continuation requires collapse-side evidence recovery: {missing}')
        terminal_selected = [ev for ev in selected_events if ev[0] not in ('splice','splice_span')]
        if terminal_selected and annotation_corrector is not None:
            annotation = annotation_corrector(row,terminal_selected,policy)
            present = {t.rsplit(':',1)[-1].split('_')[0] for t in annotation.split('|')
                       if ':' not in t or t.rsplit(':',1)[0]==row['gene']}
            if missing-present:
                raise ValueError('Collapse-side recovery did not retain every protected intronic region')
        else:
            annotation = base
        # Recovery represents the chosen source model; it is NOT an identity
        # key per transcript endpoint. Variable ends remain in one class.
        key = base
        if signature:
            key = evidence_key(base,tuple(sorted(signature)))
        source = row['library'],row['gene'],row['molecule']
        mapping[source] = key
        if known_targets is not None:
            if key in first_known_targets and first_known_targets[key]!=known_targets:
                mixed_known_targets.add(key)
            else:
                first_known_targets.setdefault(key,known_targets)
        if key not in candidates:
            candidates[key] = dict(annotation_structure=annotation,count=0,
                representative_source=source,parent_exons=merge_exons(row['source_exons']),
                protected_junctions=set(),protected_coverage=set(),protected_terminal_blocks=set(),
                protected_splice_spans=set(),
                minimum_continuation=policy.minimum_continuation)
        item = candidates[key]
        item['count'] += 1
        item['protected_junctions'].update(protected_junctions)
        item['protected_coverage'].update(protected_coverage)
        item['protected_terminal_blocks'].update(protected_terminal_blocks)
        item['protected_splice_spans'].update(protected_splice_spans)
        # Refinement only iterates these values; avoid three per-molecule set
        # allocations (including empty sets) across complete source rosters.
        prepared_rows.append((row,key,annotation,*prepared))
    # Exact legacy strings can also pool incompatible continuation boundaries.
    # Refine ONLY a class whose first physical representative contradicts its
    # selected evidence. Observed terminal splice sites use the existing +/-3
    # tolerance; variable RNA endpoints never enter the splice signature.
    bad={key for key,c in candidates.items() if not intronic_edge_compatible(
        c['protected_junctions'],c['protected_coverage'],c['parent_exons'],splice_slop=policy.splice_slop)
        or not splice_spans_compatible(c['protected_splice_spans'],c['parent_exons'],policy.splice_slop)}
    # Preserve the exact existing class keys and native tie-breaking when
    # reference compatibility is homogeneous. Refine only mixed classes.
    bad.update(mixed_known_targets)
    if bad:
        intron_bounds={name:(a,b) for name,a,b in introns}
        groups=defaultdict(list)
        for row,key,annotation,pj,pc,pt,ps,selected in prepared_rows:
            if key not in bad:continue
            exons=merge_exons(row['source_exons'])
            sites=(tuple((a[1],b[0]) for a,b in zip(exons,exons[1:]))
                   if any(ev[0] in ('splice','splice_span') for ev in selected) else
                   (() if len(exons)==1 else ((exons[0][1],exons[1][0]),(exons[-2][1],exons[-1][0]))))
            sides=[]
            if not sites:
                for kind,region,junction,coverage in selected:
                    if kind=='continuation':a,b=intron_bounds[region]
                    elif kind=='parent_gap_continuation':a,b=map(int,region.split(':',1)[1].split('-'))
                    else:continue
                    s,e=exons[0]
                    sides.append((region,s<a<=e,s<=b<e))
            shape=(sites,tuple(sorted(sides)))
            known_targets=(compatible_known_references(pj,pc,ps)
                           if key in mixed_known_targets else None)
            target=next((group for group in groups[key]
                         if len(group['sites'])==len(sites) and group['sides']==shape[1]
                         and group['known_targets']==known_targets
                         and all(max(abs(a-b) for a,b in zip(j,k))<=policy.splice_slop
                                 for j,k in zip(sites,group['sites']))),None)
            if target is None:
                partition_shape=(shape+('known_reference_compatibility',known_targets)
                                 if key in mixed_known_targets else shape)
                digest=hashlib.sha256(json.dumps(['compatible_terminal_evidence',key,partition_shape]).encode()).hexdigest()[:20]
                new_key=key.rsplit('~intronic=',1)[0]+'~intronic='+digest
                if new_key in candidates:raise ValueError('Refined source evidence key collision')
                target=dict(key=new_key,sites=sites,sides=shape[1],known_targets=known_targets)
                groups[key].append(target)
                candidates[new_key]=dict(annotation_structure=annotation,count=0,
                    representative_source=(row['library'],row['gene'],row['molecule']),parent_exons=exons,
                    protected_junctions=set(),protected_coverage=set(),protected_terminal_blocks=set(),
                    protected_splice_spans=set(),
                    minimum_continuation=policy.minimum_continuation)
            c=candidates[target['key']];c['count']+=1
            c['protected_junctions'].update(pj);c['protected_coverage'].update(pc);c['protected_terminal_blocks'].update(pt)
            c['protected_splice_spans'].update(ps)
            mapping[row['library'],row['gene'],row['molecule']]=target['key']
        for key in bad:del candidates[key]
        if any(not intronic_edge_compatible(c['protected_junctions'],c['protected_coverage'],
                c['parent_exons'],splice_slop=policy.splice_slop) or not splice_spans_compatible(
                    c['protected_splice_spans'],c['parent_exons'],policy.splice_slop) for c in candidates.values()):
            raise ValueError('A refined source class still contradicts its measured intronic evidence')
    if len(mapping)!=len(rows) or sum(c['count'] for c in candidates.values())!=len(rows):
        raise ValueError('Read identity/count conservation failed before collapse')
    return mapping,candidates


def collapse_prepared(candidates, references, gene, runtime, *, splice_slop):
    """Call the existing parent selection with actual evidence at every edge.

    references: {structure: {'id': ENST, 'exons': actual_exon_blocks}}.
    This has no independent parent-ranking implementation. The original WTA
    computes bins, parent scores, reference priority and hard assignment.
    """
    # A corrected reference and an ordinary source class can share the OLD
    # cache key while having different annotations. Never overwrite the read's
    # observed structure with the reference's restored intronic evidence.
    # Parent matching repeatedly checks the same immutable geometry within a
    # gene. The index implements the same closed +/-slop predicate, with caches
    # released on return; public scan-based checks remain independent.
    from .splice_span_index import make_checker
    spans_compatible = make_checker(merge_exons)
    resolved_references={}
    for s,r in references.items():
        item=dict(r,native_cache_structure=s)
        key=s
        if s in candidates:
            digest=hashlib.sha256((s+'\0'+r['id']).encode()).hexdigest()[:20]
            key=s+'~reference='+digest
            if key in candidates or key in references:raise ValueError('Reference/source namespace collision')
        resolved_references[key]=item
    references=resolved_references
    annotation = {k:v['annotation_structure'] for k,v in candidates.items()}
    annotation.update({s:r.get('annotation_structure',s) for s,r in references.items()})
    if any(annotation[k]!=c['annotation_structure'] for k,c in candidates.items()):
        raise ValueError('Reference annotation overwrote an observed source class')
    sr = {k:v['count'] for k,v in candidates.items()}
    # Keep the instructed native matching/ranking representation. A recovered
    # terminal I endpoint describes evidence within a continuous exon; feeding
    # it back through junction_core_tokens invents an E/I transition and can
    # displace a fully compatible Ensembl or higher-support novel parent.
    # Evidence classes and actual-genome edge checks distinguish the errors;
    # annotation recovery is retained separately for final annotation/export.
    native_annotation={k:k.rsplit('~intronic=',1)[0] for k in candidates}
    native_annotation.update({s:r['native_cache_structure'] for s,r in references.items()})
    tok = {k:tuple(runtime['structure_tokens_for_containment'](s,gene)) for k,s in native_annotation.items()}
    core = {k:tuple(runtime['junction_core_tokens'](s,gene)) for k,s in native_annotation.items()}
    enst = {s:r['id'] for s,r in references.items()}
    # Replay only native reference selection on the complete original pooled
    # structures. Evidence splitting must not change which compatible known
    # model ordinary reads prefer merely because another reference now survives.
    native_reads = Counter()
    for key, count in sr.items():
        native_reads[native_annotation[key]] += count
    native_enst = {r['native_cache_structure']: r['id'] for r in references.values()}
    original_structures = set(native_reads) | set(native_enst)
    original_tokens = {s: tuple(runtime['structure_tokens_for_containment'](s, gene))
                       for s in original_structures}
    original_cores = {s: tuple(runtime['junction_core_tokens'](s, gene))
                      for s in original_structures}
    native_choices = runtime['collapse_gene_jc'](
        native_reads, original_tokens, original_cores,
        {s: len(t) for s, t in original_tokens.items()},
        enst=native_enst, reference_targets_only=True)['source_to_reference']
    reference_keys = {r['native_cache_structure']: s for s, r in references.items()}
    preferred_references = {s: reference_keys[native_choices[native_annotation[s]]]
                            for s in sr if native_annotation[s] in native_choices}
    # A characterized retention model can itself have lost a terminal I token
    # in the old cache. Its recovered view is an additional matching view of
    # the same verified model, without changing native reference ranking.
    reference_evidence_cores={s:[tuple(runtime['junction_core_tokens'](annotation[s],gene))]
                              for s in references
                              if core[s]!=tuple(runtime['junction_core_tokens'](annotation[s],gene))}
    parent_exons = {k:v['parent_exons'] for k,v in candidates.items()}
    parent_exons.update({s:r['exons'] for s,r in references.items()})
    thresholds = {c['minimum_continuation'] for c in candidates.values()}
    if len(thresholds) != 1:
        raise ValueError('Every candidate must use the same explicit continuation policy')
    minimum_continuation = next(iter(thresholds))
    reference_junctions = set()
    for model in references.values():
        blocks = merge_exons(model['exons'])
        reference_junctions.update((a[1],b[0]) for a,b in zip(blocks,blocks[1:]))

    def compatible(evidence,parent):
        if evidence is None:
            return True
        if not spans_compatible(evidence.get('protected_splice_spans',()),
                                       parent_exons[parent],splice_slop):
            return False
        if not runtime['intronic_edge_compatible'](
            evidence['protected_junctions'],evidence['protected_coverage'],
            parent_exons[parent],splice_slop=splice_slop):
            return False
        return True

    # A valid shorter Ensembl retention model must not disappear into a longer
    # spliced reference before its supporting reads can be assigned. Propagate
    # observed compatible constraints to the reference-containment check.
    reference_constraints={}
    for ref in references:
        reference_constraints[ref]=[(k,e) for k,e in candidates.items()
            if (e['protected_junctions'] or e['protected_coverage'] or e['protected_terminal_blocks'] or e.get('protected_splice_spans'))
            and any(runtime['is_contiguous_subsequence'](core[k],t)
                    for t in [core[ref]]+reference_evidence_cores.get(ref,[]))
            and compatible(e,ref)]

    parent_coverage = {k: c['protected_coverage'] for k,c in candidates.items()}
    for ref, constraints in reference_constraints.items():
        parent_coverage[ref] = {cv for _,e in constraints for cv in e['protected_coverage']}
    for key, coverage in parent_coverage.items():
        if any(b-a+1 < minimum_continuation for a,b in coverage):
            raise ValueError('Unselected short coverage entered the retention guard')
        # Reverse contradictions use ANY covered sequence, so a genomic union
        # is exactly equivalent and avoids repeated variable-end intervals.
        parent_coverage[key] = merge_exons(coverage) if coverage else []

    def reverse_compatible(child, parent):
        # Use ONLY already-selected >=100-nt intronic coverage. A rare novel
        # gap inside an ordinary exon must not turn the canonical exon into IR.
        # That would defeat the intentional normalisation of noisy splice sites.
        coverage = parent_coverage.get(parent, ())
        if not coverage:return True
        source_exons = (candidates[child]['parent_exons'] if child in candidates
                        else references[child]['exons'])
        source_blocks = merge_exons(source_exons)
        selected = candidates.get(child, {}).get('protected_junctions', ())
        trusted_gaps = []
        for left,right in zip(source_blocks,source_blocks[1:]):
            junction = (left[1],right[0])
            if (junction in selected or any(max(abs(a-b) for a,b in zip(junction,j)) <= splice_slop
                                           for j in reference_junctions)):
                trusted_gaps.append((left[1]+1+splice_slop,right[0]-1-splice_slop))
        for start,end in coverage:
            if end-start+1 < minimum_continuation:
                raise ValueError('Unselected short coverage entered the retention guard')
            for a,b in parent_exons[parent]:
                lo,hi = max(start,a),min(end,b)
                if lo>hi:continue
                if any(max(lo,a) <= min(hi,b) for a,b in trusted_gaps):
                    return False
        return True

    reference_constraints = {ref: [(k,e) for k,e in constraints if reverse_compatible(k,ref)]
                             for ref,constraints in reference_constraints.items()}

    def edge(child,parent):
        evidence=candidates.get(child)
        if not reverse_compatible(child,parent):return False
        if evidence is not None:return compatible(evidence,parent)
        return all(compatible(e,parent) and reverse_compatible(k,parent)
                   for k,e in reference_constraints.get(child,[]))

    result = runtime['collapse_gene_jc'](sr,tok,core,{k:len(v) for k,v in tok.items()},
                                       enst=enst,edge_compatible=edge,
                                       reference_evidence_cores=reference_evidence_cores,
                                       native_keys=native_annotation,
                                       preferred_reference_targets=preferred_references)
    result['parent_retention_coverage'] = parent_coverage
    result['resolved_references'] = references
    reps = set(result['long_reps'])|set(result['other_reps'])
    final = {s:s for s in reps}
    final.update(result['assignment'])
    if set(final)!=set(candidates)|set(result.get('enst',{})):
        raise ValueError('Collapsed candidate/reference identity conservation failed')
    for child,parent in final.items():
        if child!=parent and not edge(child,parent):
            raise ValueError('An emitted assignment contradicts source evidence')
    return result,final


def token_base(token):
    return token.rsplit(':',1)[-1].split('_',1)[0]


def recover_selected_terminal_tokens(row,events,policy,regions):
    """Return the internal annotation for one evidence-protected read.

regions maps region IDs to inclusive genomic (start,end) coordinates. It must
come from the instructed annotation source. This neither runs nor modifies
gff_process.exonAnnotate; the original raw tokens and aligned blocks are pinned
inputs. Recovery is limited to selected terminal intronic events.
"""
    cache=row.get('reference_cache_structure')
    raw=row.get('raw_structure') if cache is None else cache
    if cache is not None and row.get('raw_structure') is not None:
        raise ValueError('Reference cache and raw read annotation are distinct source types')
    if raw is None:
        raise ValueError(f"Required original raw annotation unresolved: {row['library']}:{row['gene']}:{row['molecule']}")
    tokens=raw.split('|') if isinstance(raw,str) else list(raw)
    if not tokens:raise ValueError('Original raw annotation must not be empty')
    if row.get('strand') not in ('+','-'):
        raise ValueError('Verified source strand is required for terminal evidence recovery')
    selected=[]
    for kind,region,junction,covered in events:
        measured_u = kind == 'exon' and region.startswith('U')
        if measured_u and (kind,region,junction,covered) not in terminal_u_events(
                raw, row['source_exons'], regions, row['gene'], policy.splice_slop):
            raise ValueError('Selected U splice site lacks original raw/primary-alignment evidence')
        if kind!='parent_gap_continuation' and not measured_u and region not in regions:
            raise ValueError(f'Required annotated intron region unresolved: {region}')
        if kind in ('continuation','parent_gap_continuation') and covered[1]-covered[0]+1<policy.minimum_continuation:
            continue
        selected.append((kind,region,junction,covered))
    if not selected:return row['legacy_structure']
    selected_regions={ev[1] for ev in selected}
    protected_blocks=[ev[2] for ev in selected if ev[0]=='parent_gap_continuation']
    # Apply the original terminal trimming to all unprotected token types.
    # An I endpoint survives only for the explicitly selected evidence.
    for i in (() if cache is not None else (0,-1)):
        if tokens and '_' in tokens[i] and token_base(tokens[i]) not in selected_regions:
            try:coordinate=int(tokens[i].rsplit('_',1)[1])
            except ValueError:raise ValueError('Required raw terminal coordinate is invalid')
            if not any(min(abs(coordinate-s),abs(coordinate-e))<=policy.splice_slop for s,e in protected_blocks):
                del tokens[i]
    for kind,region,junction,covered in selected:
        if kind!='continuation':continue
        lo,hi=regions[region]
        matching=[ex for ex in merge_exons(row['source_exons']) if
                  ex[0]<=covered[0]<=covered[1]<=ex[1] and
                  (ex[0]<lo<=ex[1] or ex[0]<=hi<ex[1])]
        if len(matching)!=1:
            raise ValueError('Selected continuation does not have one matching original continuous exon block')
        s,e=matching[0]
        # The legacy index bug can include an adjacent E region outside this
        # actual terminal block. Only remove such tokens from the same flanking
        # exon block as this selected intron, preserving the complete originals.
        intron_number=region[1:].split('.',1)[0]
        flanks={f'E{intron_number}',f'E{int(intron_number)+1}'}
        cleaned=[]
        for token in tokens:
            # Coordinate-bearing labels mark original measured boundaries;
            # they are not interpolated flanking regions. In overlapping genes
            # the legacy annotator can prefix such a label with the endpoint's
            # gene rather than the gene that supplied the region ID. Preserve
            # that input verbatim: actual primary blocks drive compatibility.
            # This cleanup only repairs bare, interpolated flank tokens.
            if '_' in token:
                cleaned.append(token)
                continue
            if ':' in token and token.rsplit(':',1)[0]!=row['gene']:
                cleaned.append(token)
                continue
            base=token_base(token)
            if base.split('.',1)[0] in flanks:
                if base not in regions:raise ValueError(f'Required flanking annotation region unresolved: {base}')
                a,b=regions[base]
                if max(s,a)>min(e,b)+policy.splice_slop:
                    # Do not remove another actual source exon's region.
                    if not any(max(x,a)<=min(y,b)+policy.splice_slop for x,y in row['source_exons']):
                        continue
            cleaned.append(token)
        tokens=cleaned
        if region not in {token_base(t) for t in tokens if '_' not in t
                          and (':' not in t or t.rsplit(':',1)[0]==row['gene'])}:
            # Insert alongside the retained endpoint where available. This
            # restores RUNX1's missing bare I region within collapse alone.
            indices=[i for i,t in enumerate(tokens) if token_base(t)==region
                     and (':' not in t or t.rsplit(':',1)[0]==row['gene'])]
            if indices:
                i=indices[0]
                terminal_before=(s>=lo and e>hi) == (row['strand']=='+')
                tokens.insert(i+1 if terminal_before else i,region)
            else:
                terminal_before=(s>=lo and e>hi) == (row['strand']=='+')
                if terminal_before:tokens.insert(0,region)
                else:tokens.append(region)
    return '|'.join(tokens)
