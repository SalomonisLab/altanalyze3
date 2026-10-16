"""Complete pooled evidence contexts and native Tier-1 tables from verified scans."""
from collections import defaultdict,Counter
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import zlib
from .splice_evidence import (Policy,prepare_gene,recover_selected_terminal_tokens,
                              terminal_events,merge_exons)
from ...bam.coord_store import unpack_exons


def make_source_decoders(event_kinds=None,*,cache_values=True):
    """Share immutable repeated values, never RNA identities or observations.

    Highly expressed genes contain millions of distinct source RNAs with the
    same junction list. Keeping a separate decoded copy for every RNA wastes
    memory. Caches are bounded, scoped to one gene, and preserve exact values.
    """
    kinds=None if event_kinds is None else frozenset(event_kinds)
    def text(value):return value
    def raw(blob):return zlib.decompress(blob).decode()
    def exons(blob):return tuple(unpack_exons(blob))
    def events(value):
        return tuple((k,a,tuple(j) if j is not None else None,tuple(cv))
                     for k,a,j,cv in json.loads(value) if kinds is None or k in kinds)
    decoders=(text,raw,exons,events)
    return tuple(lru_cache(maxsize=65536)(f) for f in decoders) if cache_values else decoders


def _record(mol,ordinal,strand,chrom,source,blob,library,gene,*,exon_unpacker=unpack_exons):
    return dict(library=library,gene=gene,molecule=mol,ordinal=ordinal,strand=strand,
                chromosome=chrom,gff_source=source,source_exons=exon_unpacker(blob))


def dump_context(path,context):
    payload={'policy':context['policy'].__dict__,'references':context['references'],'candidates':{}}
    for k,c in context['candidates'].items():
        item=dict(c)
        for field in ('protected_junctions','protected_coverage','protected_terminal_blocks','protected_splice_spans'):
            item[field]=sorted(item.get(field,()))
        payload['candidates'][k]=item
    with gzip.open(path,'wt') as f:json.dump(payload,f,separators=(',',':'))
    h=hashlib.sha256(Path(path).read_bytes()).hexdigest()
    return dict(file=str(path),sha256=h)


def load_context(entry):
    path=Path(entry['file'])
    if hashlib.sha256(path.read_bytes()).hexdigest()!=entry['sha256']:
        raise ValueError('Changed complete gene evidence context')
    with gzip.open(path,'rt') as f:result=json.load(f)
    result['policy']=Policy(**result['policy'])
    # JSON preserves values but changes tuple/set evidence to lists. Restore
    # the preparation API's types before pooled reference constraints use them.
    for candidate in result['candidates'].values():
        for field in ('protected_junctions','protected_coverage','protected_terminal_blocks','protected_splice_spans'):
            candidate[field]={tuple(value) for value in candidate.get(field,())}
    return result


def build_contexts(scans,regions,references,outdir,policy,*,log=print,
                   _gene_partition=None, _partition_read_counts=None,
                   event_kinds=None):
    """Preserve all native counts and actual representative source identities.

    references contains the COMPLETE instructed cache for observed genes, each
    with verified primary GFF blocks and strand. Ordinary reads aggregate using
    their actual primary first-source record; no synthetic molecule is introduced.
    """
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    libs=[s['library'] for s in scans]
    if len(set(libs))!=len(libs) or not all(s['complete'] for s in scans):
        raise ValueError('Complete unique physical-source scans required; no benchmark subset accepted')
    physical_scopes={s.get('novel_junction_scope','terminal_only') for s in scans}
    if len(physical_scopes)!=1:
        raise ValueError('All physical-source scans must use the same evidence scope')
    physical_geometry=(event_kinds is None and physical_scopes=={'all_physical_junctions_and_continuity_v4'})
    full_geometry_records=physical_scopes=={'all_physical_junctions_and_continuity_v4'}
    physical_files=[list(s['sources'])[-1] for s in scans]
    conns=[sqlite3.connect(f"file:{s['database']}?mode=ro",uri=True) for s in scans]
    handles=[];baseline_handles=[];override_handles=[]
    sample_paths=[];baseline_paths=[];override_paths={};override_counts=Counter();contexts={}
    input_count=output_count=0
    try:
        genes=set()
        if _gene_partition is None:
            if _partition_read_counts is not None:
                raise ValueError('Partition counts require explicit gene identities')
            for c in conns:genes.update(r[0] for r in c.execute('SELECT DISTINCT gene FROM models'))
        else:
            genes=set(_gene_partition)
            if len(genes)!=len(_gene_partition) or len(_partition_read_counts or [])!=len(scans):
                raise ValueError('Complete unique partition identities and source counts required')
        per_source_input=Counter()
        for i,lib in enumerate(libs):
            folder=out/f'sample_{i:04d}';folder.mkdir()
            ep=folder/'evidence.exemplars.tsv';bp=folder/'baseline.exemplars.tsv';op=folder/'molecule_classes.tsv.gz'
            eh=ep.open('w');bh=bp.open('w');oh=gzip.open(op,'wt')
            for f in (eh,bh):f.write('gene\tstrand\tstructure\tcount\trep_molecule_id\n')
            oh.write('gene\tmolecule\toriginal_structure\tcandidate_class\n')
            handles.append(eh);baseline_handles.append(bh);override_handles.append(oh)
            sample_paths.append((lib,str(ep)));baseline_paths.append((lib,str(bp)));override_paths[lib]=str(op)
        for gene in sorted(genes):
            text_value,raw_value,exon_value,event_value=make_source_decoders(
                event_kinds,cache_values=full_geometry_records)
            rows=[];event_rows_by_lib=[];models_by_lib=[];ordinary_by_lib=[]
            for i,c in enumerate(conns):
                evrows=[]
                for mol,ordinal,legacy,raw,strand,chrom,source,blob,ev in c.execute(
                    'SELECT mol,ordinal,legacy,raw,strand,chrom,source,exons,events FROM events WHERE gene=? ORDER BY ordinal',(gene,)):
                    r=_record(mol,ordinal,text_value(strand),text_value(chrom),text_value(source),
                              blob,libs[i],gene,exon_unpacker=exon_value)
                    r.update(legacy_structure=text_value(legacy),raw_structure=raw_value(raw),
                             measured_events=event_value(ev),
                             support_source=(physical_files[i],gene,mol))
                    rows.append(r);evrows.append(r)
                event_rows_by_lib.append(evrows)
                models={}
                for legacy,count,mol,ordinal,strand,chrom,source,blob,has_event in c.execute(
                    'SELECT legacy,count,mol,ordinal,strand,chrom,source,exons,has_event FROM models WHERE gene=?',(gene,)):
                    models[legacy]=dict(count=count,has_event=has_event,record=_record(mol,ordinal,strand,chrom,source,blob,libs[i],gene))
                    input_count+=count
                    per_source_input[i]+=count
                    baseline_handles[i].write(f'{gene}\t{strand}\t{legacy}\t{count}\t{mol}\n')
                models_by_lib.append(models)
                ordinary={}
                for legacy,mol,ordinal,strand,chrom,source,blob in c.execute(
                    'SELECT legacy,mol,ordinal,strand,chrom,source,exons FROM ordinary_first WHERE gene=?',(gene,)):
                    ordinary[legacy]=_record(mol,ordinal,strand,chrom,source,blob,libs[i],gene)
                ordinary_by_lib.append(ordinary)
            if not rows:
                for i,models in enumerate(models_by_lib):
                    for legacy,m in models.items():
                        r=m['record'];handles[i].write(f"{gene}\t{r['strand']}\t{legacy}\t{m['count']}\t{r['molecule']}\n")
                        output_count+=m['count']
                continue
            introns=[(k,*v) for k,v in regions[gene].items() if k.startswith('I')]
            ref=references.get(gene,{})
            mapping,candidates=prepare_gene(rows,introns,policy,
                protect_splice_geometry=physical_geometry,
                parent_models={s:r['exons'] for s,r in ref.items()},
                annotation_corrector=lambda r,e,p:recover_selected_terminal_tokens(r,e,p,regions[gene]),
                required_identities=[(r['library'],gene,r['molecule']) for r in rows],
                terminal_event_provider=lambda r:r['measured_events'])
            original_counts=Counter(r['legacy_structure'] for r in rows)
            # Per-library native class totals and the actual FIRST unprotected
            # molecule after splitting must be resolved before pooled collapse.
            per_sample=[];changed=0
            for i,models in enumerate(models_by_lib):
                class_counts=Counter();first={};protected=Counter();unprotected_event=defaultdict(list)
                for r in event_rows_by_lib[i]:
                    old=r['legacy_structure'];key=mapping[libs[i],gene,r['molecule']]
                    if key!=old:
                        protected[old]+=1;class_counts[key]+=1;changed+=1;override_counts[libs[i]]+=1
                        override_handles[i].write(f"{gene}\t{r['molecule']}\t{old}\t{key}\n")
                        first.setdefault(key,r)
                    else:unprotected_event[old].append(r)
                for legacy,m in models.items():
                    remaining=m['count']-protected[legacy]
                    if remaining<0:raise ValueError('Protected reads exceed complete original structure count')
                    if not remaining:continue
                    if not m['has_event']:r=m['record']
                    else:
                        possibilities=unprotected_event[legacy]+([ordinary_by_lib[i][legacy]] if legacy in ordinary_by_lib[i] else [])
                        if not possibilities:raise ValueError('Required first unprotected primary representative unresolved')
                        r=min(possibilities,key=lambda r:r['ordinal'])
                    class_counts[legacy]=remaining;first[legacy]=r
                if sum(class_counts.values())!=sum(m['count'] for m in models.values()):
                    raise ValueError('Full per-library read count conservation failed')
                per_sample.append((class_counts,first))
            if changed:
                # Reset the sparse possible-event counts to COMPLETE pooled
                # class counts, then pin each representative to native library
                # order / first appearance within that library.
                for c in candidates.values():c['count']=0
                for i,(counts,first) in enumerate(per_sample):
                    for key,count in counts.items():
                        r=first[key]
                        if key not in candidates:
                            candidates[key]=dict(annotation_structure=key,count=0,
                                protected_junctions=set(),protected_coverage=set(),protected_terminal_blocks=set(),
                                minimum_continuation=policy.minimum_continuation)
                        c=candidates[key]
                        if c['count']==0:
                            c.update(representative_source=(libs[i],gene,r['molecule']),parent_exons=merge_exons(r['source_exons']),
                                     source_record=dict(chrom=r['chromosome'],strand=r['strand'],source=r['gff_source'],exons=r['source_exons']))
                        c['count']+=count
                # No sparse possible-event-only candidate can survive with zero
                # reads; this is an accounting failure, never a filter.
                if any(c['count']<=0 for c in candidates.values()):raise ValueError('Incomplete class read coverage')
                if sum(c['count'] for c in candidates.values())!=sum(m['count'] for models in models_by_lib for m in models.values()):
                    raise ValueError('Complete pooled gene read count conservation failed')
                refs={s:dict(r) for s,r in ref.items()}
                for s,r in refs.items():
                    ev=[e for e in terminal_events(r['exons'],introns) if e[0]=='continuation'
                        and e[3][1]-e[3][0]+1>=policy.minimum_continuation]
                    if ev:r['annotation_structure']=recover_selected_terminal_tokens(dict(library='ENSEMBL',gene=gene,
                        molecule=r['id'],strand=r['strand'],source_exons=r['exons'],legacy_structure=s,
                        reference_cache_structure=s),ev,policy,regions[gene])
                path=out/f'gene_{len(contexts):06d}.json.gz'
                contexts[gene]=dump_context(path,dict(candidates=candidates,references=refs,policy=policy))
            for i,(counts,first) in enumerate(per_sample):
                for key in sorted(counts,key=lambda k:first[k]['ordinal']):
                    r=first[key];handles[i].write(f"{gene}\t{r['strand']}\t{key}\t{counts[key]}\t{r['molecule']}\n")
                    output_count+=counts[key]
            if len(contexts)%100==0 and changed:log(f'[pooled-source-context] {len(contexts):,} evidence genes; {gene}')
        expected_counts=(_partition_read_counts if _gene_partition is not None
                         else [s['native_eligible_reads'] for s in scans])
        if ([per_source_input[i] for i in range(len(scans))]!=list(expected_counts)
                or input_count!=output_count or input_count!=sum(expected_counts)):
            raise ValueError('Full cohort/native read count conservation failed')
    finally:
        for c in conns:c.close()
        for f in handles+baseline_handles+override_handles:f.close()
    result=dict(ordered_libraries=libs,all_native_reads=input_count,evidence_genes=len(contexts),
                protected_reads=sum(override_counts.values()),policy=policy.__dict__,gene_evidence=contexts,
                evidence_sample_paths=sample_paths,baseline_sample_paths=baseline_paths,
                override_paths=override_paths,override_counts=dict(override_counts))
    result['source_read_counts']=[per_source_input[i] for i in range(len(scans))]
    if _gene_partition is not None:
        result['partition_only']=True
        result['partition_genes']=sorted(genes)
    (out/'complete_context_manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
