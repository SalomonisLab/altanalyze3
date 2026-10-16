"""Streaming primary-alignment input for cohort-scale collapse evidence.

Read annotation and GFF are emitted in the same transcript order by the native
annotator. Reconcile every native annotation identity with its actual aligned
block, aggregate ordinary source models immediately, and retain per-molecule
records only for possible selected terminal events. All input molecules still
contribute to native structure counts. This avoids indexing billions of ordinary
exon records in a second database.
"""
from __future__ import annotations
from bisect import bisect_left,bisect_right
from collections import defaultdict
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sqlite3
import zlib
from .. import io_utils as IO
from ...bam.coord_store import pack_exons,unpack_exons
from .source_evidence import attributes
from .splice_evidence import merge_exons,Policy,terminal_u_events


def source_blocks(path):
    """Yield original consecutive primary transcript blocks, in source order."""
    info=None;cols=None;exons=[]
    with IO.smart_open(path) as f:
        for line in f:
            if line.startswith('#'):continue
            p=line.rstrip('\n').split('\t')
            if len(p)!=9 or p[2] not in ('exon','transcript'):continue
            if p[8]!=info:
                if exons:
                    attr=attributes(info)
                    mol=attr.get('transcript_id')
                    if not mol:raise ValueError(f'Primary transcript identity unresolved in {path}')
                    yield mol,cols[0],cols[6],cols[1],exons
                info=p[8];cols=p;exons=[]
            if p[2]=='exon':exons.append((int(p[3]),int(p[4])))
        if exons:
            mol=attributes(info).get('transcript_id')
            if not mol:raise ValueError(f'Primary transcript identity unresolved in {path}')
            yield mol,cols[0],cols[6],cols[1],exons


def build_event_index(regions,parent_models):
    """Precompute actual reference gap boundaries; retain every annotated intron."""
    introns=sorted((a,b,k) for k,(a,b) in regions.items() if k.startswith('I'))
    nonoverlap=all(a[1]<b[0] for a,b in zip(introns,introns[1:]))
    gaps=set()
    for exons in parent_models.values():
        merged=merge_exons(exons)
        gaps.update((a[1]+1,b[0]-1) for a,b in zip(merged,merged[1:]))
    starts=sorted((a,b) for a,b in gaps)
    ends=sorted((b,a) for a,b in gaps)
    return dict(introns=introns,intron_starts=[a for a,b,k in introns],
                intron_ends=[b for a,b,k in introns],nonoverlap=nonoverlap,
                gap_starts=starts,gap_start_coords=[a for a,b in starts],
                gap_ends=ends,gap_end_coords=[b for b,a in ends])


def measured_splice_events(exons):
    """Measure every physical splice gap, independent of annotation membership.

    Canonical sources can also be misassigned to coordinate-normalized alternative
    references. Neither reference membership nor a native region edge exempts an
    observed junction from the same pooled recurrence/compatibility policy.
    RNA endpoints and adjacent subregions are never treated as splice junctions.
    """
    blocks = merge_exons(exons)
    return [('splice', 'physical_splice', (left[1], right[0]), (left[1], right[0]))
            for left, right in zip(blocks, blocks[1:])]


def measured_terminal_events(exons,index):
    """Equivalent to terminal_events, querying only overlapping/boundary gaps."""
    exons=merge_exons(exons);events=[]
    for i in sorted({0,len(exons)-1}):
        s,e=exons[i]
        if index['nonoverlap']:
            lo=bisect_left(index['intron_ends'],s)
            hi=bisect_right(index['intron_starts'],e)
            introns=index['introns'][lo:hi]
        else:introns=index['introns']
        for a,b,region in introns:
            if a<=s<=e<=b and len(exons)>1:
                j=(e,exons[i+1][0]) if i==0 else (exons[i-1][1],s)
                events.append(('exon',region,j,(s,e)))
            elif (s<a<=e or s<=b<e) and max(s,a)<=min(e,b):
                events.append(('continuation',region,None,(max(s,a),min(e,b))))
        gaps=set(index['gap_starts'][bisect_right(index['gap_start_coords'],s):
                                     bisect_right(index['gap_start_coords'],e)])
        gaps.update((a,b) for b,a in index['gap_ends'][bisect_left(index['gap_end_coords'],s):
                                        bisect_left(index['gap_end_coords'],e)])
        for a,b in sorted(gaps):
            lo,hi=max(s,a),min(e,b)
            if lo<=hi and (s<a<=e or s<=b<e):
                events.append(('parent_gap_continuation',f'gap:{a}-{b}',(s,e),(lo,hi)))
    return list(dict.fromkeys(events))


def scan_sample(sample,regions,event_indexes,workdir,policy,*,limit=None,log=print,
                protect_novel_junctions=False):
    """Reconcile sources and retain all counts before any optional evidence fold.

    limit is exclusively a bounded diagnostic/benchmark; such outputs are marked
    incomplete and must never be accepted by the full-cohort execution gate.
    """
    lib,_,ta,gff=sample[:4]
    ta=IO.resolve(str(ta));gff=IO.resolve(str(gff))
    raw=IO.resolve(str(Path(ta).parent/'transcript_associations_raw.txt'))
    out=Path(workdir);out.mkdir(parents=True,exist_ok=True)
    path=out/'source_scan.sqlite'
    if path.exists():raise ValueError('A fresh source-scan directory is required')
    conn=sqlite3.connect(path)
    conn.execute('PRAGMA journal_mode=OFF');conn.execute('PRAGMA synchronous=OFF')
    conn.execute('CREATE TABLE events (gene TEXT,mol TEXT,ordinal INTEGER,legacy TEXT,raw BLOB,strand TEXT,chrom TEXT,source TEXT,exons BLOB,events TEXT,PRIMARY KEY(gene,mol))')
    conn.execute('CREATE TABLE models (gene TEXT,legacy TEXT,count INTEGER,mol TEXT,ordinal INTEGER,strand TEXT,chrom TEXT,source TEXT,exons BLOB,has_event INTEGER,PRIMARY KEY(gene,legacy))')
    conn.execute('CREATE TABLE ordinary_first (gene TEXT,legacy TEXT,mol TEXT,ordinal INTEGER,strand TEXT,chrom TEXT,source TEXT,exons BLOB,PRIMARY KEY(gene,legacy))')
    models={};ordinary={};event_batch=[];total=eligible=possible=extra=0
    giter=iter(source_blocks(gff));block=next(giter,None)
    stats={p:(Path(p).stat().st_size,Path(p).stat().st_mtime_ns) for p in (ta,raw,gff)}
    try:
        with IO.smart_open(ta) as a,IO.smart_open(raw) as b:
            for ordinal,pair in enumerate(itertools.zip_longest(a,b)):
                if limit is not None and ordinal>=limit:break
                x,y=pair
                if x is None or y is None:raise ValueError(f'Complete raw/native roster mismatch for {lib}')
                p=x.rstrip('\n').split('\t');r=y.rstrip('\n').split('\t')
                if len(p)<5 or len(r)<5 or (p[0],p[1],p[3],p[4])!=(r[0],r[1],r[3],r[4]):
                    raise ValueError(f'Original annotation identity mismatch: {lib}:{ordinal}')
                total+=1
                while block is not None and block[0]!=p[3]:
                    extra+=1;block=next(giter,None)
                if block is None:raise ValueError(f'Required primary alignment unresolved: {lib}:{p[0]}:{p[3]}')
                mol,chrom,strand,source,exons=block;block=next(giter,None)
                if strand!=p[1]:raise ValueError(f'Original source strand mismatch: {lib}:{mol}')
                if not p[2] or 'UNK' in p[0]:continue # native pre-existing eligibility only
                gene,legacy=p[0],p[2];eligible+=1
                if gene not in regions:raise ValueError(f'Required gene annotation unresolved: {gene}')
                ev=measured_terminal_events(exons,event_indexes[gene])
                ev.extend(terminal_u_events(r[2],exons,regions[gene],gene,policy.splice_slop))
                if protect_novel_junctions:
                    ev.extend(measured_splice_events(exons))
                    # Continuous mono-exon records can contradict a reference
                    # or recurrent alternative splice. Preserve their geometry
                    # before exact aggregation, even without a terminal event.
                ev=[e for e in ev if e[0] in ('exon','splice') or e[3][1]-e[3][0]+1>=policy.minimum_continuation]
                if protect_novel_junctions and not ev:
                    ev.append(('geometry','source_geometry',None,(0,0)))
                key=gene,legacy
                if key not in models:
                    models[key]=[1,mol,ordinal,strand,chrom,source,pack_exons(exons),bool(ev)]
                else:models[key][0]+=1
                if not ev and models[key][-1] and key not in ordinary:
                    ordinary[key]=[mol,ordinal,strand,chrom,source,pack_exons(exons)]
                if ev:
                    possible+=1
                    event_batch.append((gene,mol,ordinal,legacy,zlib.compress(r[2].encode()),strand,chrom,
                                        source,pack_exons(exons),json.dumps(ev,separators=(',',':'))))
                    if len(event_batch)>=10000:
                        conn.executemany('INSERT INTO events VALUES (?,?,?,?,?,?,?,?,?,?)',event_batch);event_batch=[]
                if total%1000000==0:log(f'[source-scan:{lib}] {total:,} original identities; {possible:,} possible evidence molecules')
        conn.executemany('INSERT INTO events VALUES (?,?,?,?,?,?,?,?,?,?)',event_batch)
        conn.executemany('INSERT INTO models VALUES (?,?,?,?,?,?,?,?,?,?)',
                         ((g,s,*m) for (g,s),m in models.items()))
        conn.executemany('INSERT INTO ordinary_first VALUES (?,?,?,?,?,?,?,?)',
                         ((g,s,*m) for (g,s),m in ordinary.items()))
        conn.commit()
        if sum(m[0] for m in models.values())!=eligible:raise ValueError('Native read count conservation failed')
        if protect_novel_junctions and (possible!=eligible or
                conn.execute('SELECT COUNT(*) FROM events').fetchone()[0]!=eligible):
            raise ValueError('Every eligible source geometry must survive before aggregation')
        for p,before in stats.items():
            after=Path(p).stat()
            if before!=(after.st_size,after.st_mtime_ns):raise ValueError(f'Original source changed during scan: {p}')
        manifest=dict(library=lib,sources={p:dict(bytes=n,mtime_ns=mt) for p,(n,mt) in stats.items()},
            total_annotation_rows=total,native_eligible_reads=eligible,native_structures=len(models),
            possible_evidence_reads=possible,primary_blocks_outside_native_annotation=extra,
            complete=limit is None,policy=policy.__dict__,database=str(path),
            protect_novel_junctions=protect_novel_junctions,
            novel_junction_scope=('all_physical_junctions_and_continuity_v4'
                                  if protect_novel_junctions else 'terminal_only'))
        (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        log(f'[source-scan:{lib}] complete={manifest["complete"]}; all {eligible:,} native reads counted')
        return manifest
    finally:conn.close()
