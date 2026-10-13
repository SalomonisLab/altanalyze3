"""Prepare the opt-in terminal-evidence classes from unmodified source artifacts.

All molecules in native Tier 1's input roster are reconciled with raw annotation
and primary alignment exons. A disk index bounds ingestion memory. Nothing here
annotates reads, substitutes GFF models, or changes native WTA ranking.
"""
from __future__ import annotations
import collections
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import re


import sqlite3

from .. import io_utils as IO
from .splice_evidence import (Policy, prepare_gene, recover_selected_terminal_tokens,
                             merge_exons, terminal_events, terminal_u_events)

ATTR=re.compile(r'(?:^|;)\s*([A-Za-z_]+)(?:=|\s+")([^";]+)')

def distinct_source_samples(samples):
    """Preserve every input row while disambiguating repeated library labels."""
    counts = collections.Counter(s[0] for s in samples)
    result = []
    labels = {}
    for sample in samples:
        name = sample[0]
        if counts[name] > 1:
            physical = (str(Path(sample[1]).resolve()), str(Path(IO.resolve(str(sample[2]))).resolve()))
            name += '__source' + hashlib.sha256(json.dumps(physical).encode()).hexdigest()[:12]
        if name in labels:
            raise ValueError('Repeated identical source row requires explicit source identity reconciliation')
        labels[name] = sample[0]
        result.append((name,) + tuple(sample[1:]))
    return result, labels


def attributes(info):
    return dict(ATTR.findall(info))


def _sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def _rows(path):
    with IO.smart_open(path) as f:
        for line in f:
            p=line.rstrip('\n').split('\t')
            if len(p)<5:raise ValueError(f'Malformed required annotation source: {path}')
            yield p


def prepare_sources(samples, exon_reference, enst_cache, ref_gff, workdir, policy, *, log=print):
    """Return gene evidence and library-scoped overrides for complete input samples.

    samples are native (library, h5ad, TA, primary_read_GFF) tuples. Every source
    is explicit except raw TA, the sibling artifact produced with TA by native
    annotation. A missing or conflicting identity rejects preparation; callers
    must resolve the source rather than accepting a subset.
    """
    if not isinstance(policy,Policy):raise ValueError('An explicit evidence policy is required')
    libraries=[s[0] for s in samples]
    if not libraries or len(set(libraries))!=len(libraries):
        raise ValueError('Complete unique ordered source library roster is required')
    inputs={}
    resolved=[]
    for sample in samples:
        if len(sample)<4:raise ValueError(f'Primary read GFF required for {sample[0]}')
        lib,_,ta,gff=sample[:4]
        ta=IO.resolve(str(ta));gff=IO.resolve(str(gff))
        raw=IO.resolve(str(Path(ta).parent/'transcript_associations_raw.txt'))
        for path in (ta,raw,gff):
            if not os.path.isfile(path):raise ValueError(f'Required original source unresolved: {path}')
            inputs[path]=None
        resolved.append((lib,ta,raw,gff))
    exon_reference=IO.resolve(str(exon_reference));enst_cache=IO.resolve(str(enst_cache))
    ref_gff=IO.resolve(str(ref_gff))
    for path in (exon_reference,enst_cache,ref_gff):
        if not os.path.isfile(path):raise ValueError(f'Required original source unresolved: {path}')
        inputs[path]=None
    workdir=Path(workdir);workdir.mkdir(parents=True,exist_ok=True)
    db=workdir/'source_reads.sqlite'
    if db.exists():raise ValueError('Fresh evidence work directory required; cached source classes cannot be reused')
    conn=sqlite3.connect(db)
    conn.execute('PRAGMA journal_mode=WAL')
    conn.execute('CREATE TABLE reads (lib INTEGER,gene TEXT,mol TEXT,ordinal INTEGER,strand TEXT,legacy TEXT,raw TEXT,PRIMARY KEY(lib,gene,mol))')
    conn.execute('CREATE INDEX read_molecule ON reads(lib,mol)')
    conn.execute('CREATE INDEX read_gene ON reads(gene,lib,ordinal)')
    conn.execute('CREATE TABLE exons (lib INTEGER,mol TEXT,start INTEGER,end INTEGER,chrom TEXT,strand TEXT)')
    conn.execute('CREATE INDEX exon_molecule ON exons(lib,mol)')
    totals={}
    try:
        for lib_index,(lib,ta,raw,gff) in enumerate(resolved):
            batch=[];valid=0;native_ignored=0
            for ordinal,pair in enumerate(itertools.zip_longest(_rows(ta),_rows(raw))):
                p,r=pair
                if p is None or r is None or (p[0],p[1],p[3],p[4])!=(r[0],r[1],r[3],r[4]):
                    raise ValueError(f'Complete paired raw/native annotation identity check failed: {lib}:{ordinal}')
                if not p[2] or 'UNK' in p[0]:
                    # Exactly native Tier 1's pre-existing eligibility rule.
                    native_ignored+=1;continue
                if p[1] not in ('+','-'):raise ValueError(f'Unresolved source strand: {lib}:{p[3]}')
                batch.append((lib_index,p[0],p[3],ordinal,p[1],p[2],r[2]));valid+=1
                if len(batch)>=20000:
                    conn.executemany('INSERT INTO reads VALUES (?,?,?,?,?,?,?)',batch);batch=[]
            conn.executemany('INSERT INTO reads VALUES (?,?,?,?,?,?,?)',batch);conn.commit()
            duplicate=conn.execute('SELECT mol FROM reads WHERE lib=? GROUP BY mol HAVING COUNT(*)>1 LIMIT 1',(lib_index,)).fetchone()
            if duplicate:raise ValueError(f'Primary source molecule is assigned to multiple genes: {lib}:{duplicate[0]}')
            batch=[];source_exons=0
            with IO.smart_open(gff) as f:
                for line in f:
                    if line.startswith('#'):continue
                    p=line.rstrip('\n').split('\t')
                    if len(p)!=9 or p[2]!='exon':continue
                    mol=attributes(p[8]).get('transcript_id')
                    if mol is None:raise ValueError(f'Primary exon transcript identity unresolved in {gff}')
                    a,b=int(p[3]),int(p[4])
                    if a>b:raise ValueError('Primary genomic exon coordinates must be ascending')
                    batch.append((lib_index,mol,a,b,p[0],p[6]));source_exons+=1
                    if len(batch)>=20000:
                        conn.executemany('INSERT INTO exons VALUES (?,?,?,?,?,?)',batch);batch=[]
            conn.executemany('INSERT INTO exons VALUES (?,?,?,?,?,?)',batch);conn.commit()
            missing=conn.execute('SELECT r.gene,r.mol FROM reads r WHERE r.lib=? AND NOT EXISTS (SELECT 1 FROM exons e WHERE e.lib=r.lib AND e.mol=r.mol) LIMIT 1',(lib_index,)).fetchone()
            if missing:raise ValueError(f'Required primary source alignment unresolved: {lib}:{missing}; no subset accepted')
            mismatch=conn.execute('SELECT r.gene,r.mol FROM reads r JOIN exons e ON e.lib=r.lib AND e.mol=r.mol WHERE r.lib=? AND e.strand!=r.strand LIMIT 1',(lib_index,)).fetchone()
            if mismatch:raise ValueError(f'Primary alignment/annotation strand mismatch: {lib}:{mismatch}')
            totals[lib]=dict(native_eligible_reads=valid,native_ineligible_reads=native_ignored,primary_exon_records=source_exons)
            log(f'[intronic-evidence:{lib}] all {valid:,} native read identities resolved')
        genes={r[0] for r in conn.execute('SELECT DISTINCT gene FROM reads')}
        regions=collections.defaultdict(dict)
        with IO.smart_open(exon_reference) as f:
            for r in csv.DictReader(f,delimiter='\t'):
                if r['gene'] not in genes:continue
                a,b=int(r['exon-region-start(s)']),int(r['exon-region-stop(s)'])
                value=min(a,b),max(a,b)
                prev=regions[r['gene']].get(r['exon-id'])
                if prev is not None and prev!=value:raise ValueError('Conflicting required annotation region')
                regions[r['gene']][r['exon-id']]=value
        if genes-set(regions):raise ValueError(f'Required annotated source gene regions unresolved: {sorted(genes-set(regions))}')
        verified_references=load_reference_models(enst_cache,ref_gff,genes=genes)
        wanted={(g,r['id']) for g,m in verified_references.items() for r in m.values()}
        evidence={};overrides={lib:{} for lib in libraries};prepared=0
        for gene in sorted(genes):
            rows=[]
            for i,mol,ordinal,strand,legacy,raw in conn.execute('SELECT lib,mol,ordinal,strand,legacy,raw FROM reads WHERE gene=? ORDER BY lib,ordinal',(gene,)):
                blocks=list(conn.execute('SELECT start,end,chrom FROM exons WHERE lib=? AND mol=?',(i,mol)))
                if len({b[2] for b in blocks})!=1:raise ValueError(f'Required single-chromosome source unresolved: {libraries[i]}:{gene}:{mol}')
                rows.append(dict(library=libraries[i],gene=gene,molecule=mol,strand=strand,
                    legacy_structure=legacy,raw_structure=raw,source_exons=[b[:2] for b in blocks]))
            refs={s:dict(r) for s,r in verified_references.get(gene,{}).items()}
            introns=[(k,*v) for k,v in regions[gene].items() if k.startswith('I')]
            # The original reference cache can suffer the same terminal I-token
            # omission as the reads. Recover only verified >=100-nt continuation
            # so a valid characterized retention model can still receive reads.
            # This is not a second terminal-trimming pass over the saved cache.
            for s,r in refs.items():
                selected=[e for e in terminal_events(r['exons'],introns)
                          if e[0]=='continuation' and e[3][1]-e[3][0]+1>=policy.minimum_continuation]
                if selected:
                    r['annotation_structure']=recover_selected_terminal_tokens(dict(
                        library='ENSEMBL',gene=gene,molecule=r['id'],legacy_structure=s,
                        reference_cache_structure=s,source_exons=r['exons'],
                        strand=r['strand']),selected,policy,regions[gene])
            parent_models={s:r['exons'] for s,r in refs.items()}
            mapping,candidates=prepare_gene(rows,introns,policy,
                annotation_corrector=lambda r,e,p:recover_selected_terminal_tokens(r,e,p,regions[gene]),
                parent_models=parent_models,
                terminal_event_provider=lambda r:terminal_events(r['source_exons'],introns,parent_models)
                    + terminal_u_events(r['raw_structure'],r['source_exons'],regions[gene],gene,policy.splice_slop),
                required_identities=[(libraries[i],gene,mol) for i,mol in conn.execute('SELECT lib,mol FROM reads WHERE gene=?',(gene,))])
            changed=[r for r in rows if mapping[r['library'],gene,r['molecule']]!=r['legacy_structure']]
            if changed:
                evidence[gene]=dict(candidates=candidates,references=refs,policy=policy)
                for r in changed:
                    overrides[r['library']][gene,r['molecule']]=(r['legacy_structure'],mapping[r['library'],gene,r['molecule']])
            prepared+=len(rows)
        if prepared!=sum(r['native_eligible_reads'] for r in totals.values()):
            raise ValueError('Complete read identity/count conservation failed')
        for path in inputs:inputs[path]=_sha(path)
        manifest=dict(status='all_original_sources_verified',ordered_libraries=libraries,
            policy=policy.__dict__,source_sha256=inputs,source_counts=totals,all_native_eligible_reads=prepared,
            complete_reference_identities=len(wanted),protected_genes=len(evidence),
            protected_reads=sum(len(m) for m in overrides.values()),upstream_annotation_modified=False)
        (workdir/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        return evidence,overrides
    finally:conn.close()


def load_reference_models(enst_cache,ref_gff,*,genes=None,original_raw=None):
    """Reconcile native cache identities to exact measured reference transcripts.

    The reference builder strips transcript versions for display and assigns the
    Ensembl91 gene namespace through exonAnnotate. Primary GENCODE gene labels
    can differ. The original raw reference rows retain the full transcript ID,
    including PAR loci, so use those bindings without merging separate loci.
    """
    wanted=collections.defaultdict(dict)
    with IO.smart_open(enst_cache) as f:
        next(f)
        for line in f:
            g,s,tid=line.rstrip('\n').split('\t')[:3]
            if genes is None or g in genes:wanted[g][s]=tid
    ids={t for m in wanted.values() for t in m.values()}
    bindings={}
    original_raw=original_raw or str(Path(enst_cache).parent/'transcript_associations_raw.txt')
    if IO.exists(str(original_raw)):
        with IO.smart_open(str(original_raw)) as f:
            for line in f:
                p=line.rstrip('\n').split('\t')
                if len(p)<4 or p[3].split('.')[0] not in ids:continue
                ts=p[2].split('|')
                if ts and '_' in ts[0]:del ts[0]
                if ts and '_' in ts[-1]:del ts[-1]
                s='|'.join(ts);tid=p[3].split('.')[0]
                if wanted.get(p[0],{}).get(s)==tid:
                    # Same last-row semantics as reference._load_from_transcript_associations.
                    bindings[p[0],s,tid]=p[3]
    primary={};variants=collections.defaultdict(set)
    with IO.smart_open(ref_gff) as f:
        for line in f:
            if line.startswith('#'):continue
            p=line.rstrip('\n').split('\t')
            if len(p)!=9 or p[2]!='exon':continue
            attr=attributes(p[8]);tid=attr.get('transcript_id');gene=attr.get('gene_id')
            if tid is None or gene is None:raise ValueError('Required original reference exon identifiers unresolved')
            if tid.split('.')[0] not in ids:continue
            if tid not in primary:
                primary[tid]=dict(exons=[],strand=p[6],chromosome=p[0],primary_gene_id=gene.split('.')[0],source_transcript_id=tid)
            r=primary[tid]
            if (r['strand'],r['chromosome'])!=(p[6],p[0]):raise ValueError('One full reference identity has conflicting primary loci')
            r['exons'].append((int(p[3]),int(p[4])));variants[tid.split('.')[0]].add(tid)
    result={};missing=[]
    for gene,m in wanted.items():
        result[gene]={}
        for s,tid in m.items():
            full=bindings.get((gene,s,tid))
            if full is None:
                candidates=variants[tid]
                if len(candidates)!=1:
                    missing.append((gene,s,tid,'full original reference identity unresolved'));continue
                full=next(iter(candidates))
            if full not in primary:
                missing.append((gene,s,tid,full));continue
            result[gene][s]=dict(primary[full],id=tid,exons=merge_exons(primary[full]['exons']))
    if missing:raise ValueError(f'Complete original reference coverage failed: {missing}; no subset accepted')
    return result
