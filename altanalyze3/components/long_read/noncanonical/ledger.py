"""Lossless read/event tables and auditable, non-destructive count reconciliation."""
from __future__ import annotations
import collections
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
from .manifest import atomic_json, digest, read_tsv, write_tsv, verify_receipt, receipt


def event_id(kind, geometry):
    payload=json.dumps(geometry,separators=(',',':'),sort_keys=True)
    return kind+'|'+hashlib.sha256(payload.encode()).hexdigest()[:32]

def initialize(db):
    db.executescript('''
    CREATE TABLE sources(source_id TEXT PRIMARY KEY, uid TEXT, library TEXT, group_name TEXT, original_json TEXT);
    CREATE TABLE reads(source_id TEXT, query_name TEXT, barcode TEXT, molecule_id TEXT,
      umi TEXT, primary_seen INT, legacy_feature TEXT, legacy_status TEXT, metadata_json TEXT,
      PRIMARY KEY(source_id,query_name));
    CREATE TABLE segments(source_id TEXT, query_name TEXT, alignment_index INT,
      chrom TEXT, start0 INT, end0 INT, reverse INT, flag INT, mapq INT, cigar TEXT, sa TEXT,
      PRIMARY KEY(source_id,query_name,alignment_index));
    CREATE TABLE events(event_id TEXT PRIMARY KEY, modality TEXT, geometry_json TEXT);
    CREATE TABLE support(event_id TEXT, source_id TEXT, query_name TEXT, status TEXT,
      label TEXT, row_index INT, evidence_json TEXT);
    CREATE INDEX support_identity ON support(source_id,query_name);
    CREATE TABLE assignments(source_id TEXT, query_name TEXT, legacy_feature TEXT,
      labels_json TEXT, action TEXT, reason TEXT, PRIMARY KEY(source_id,query_name));
    ''')

def support(db, sid, query, kind, geometry, status, label, evidence, index=0):
    if db.execute('SELECT 1 FROM reads WHERE source_id=? AND query_name=?',(sid,query)).fetchone() is None:
        raise ValueError(f'Unknown source query in detector output: {sid} / {query}')
    eid=event_id(kind,geometry)
    db.execute('INSERT OR IGNORE INTO events VALUES(?,?,?)',(eid,kind,json.dumps(geometry,sort_keys=True)))
    db.execute('INSERT INTO support VALUES(?,?,?,?,?,?,?)',
               (eid,sid,query,status,label,index,json.dumps(evidence,sort_keys=True)))


def export_query(db,path,sql):
    cursor=db.execute(sql)
    with gzip.open(path,'wt') as handle:
        import csv
        writer=csv.writer(handle,delimiter='\t')
        writer.writerow([r[0] for r in cursor.description]);writer.writerows(cursor)


def route(legacy_feature, labels, eligible_crosschrom):
    """The accepted first release policy; interpretation ambiguity is never a vetoed read.

    Native candidate-class cross-chromosome evidence contradicts a single-locus
    whole-read model. PTD/circle ambiguity and variant annotations alone do not.
    """
    if legacy_feature and eligible_crosschrom:
        return 'reallocate','native cross-chromosome candidate contradicts single-locus whole-read identity'
    if legacy_feature and labels:
        return 'retain_with_labels','compatible or unresolved additional identity; no automatic subtraction'
    if legacy_feature:
        return 'retain_standard','no incompatible accepted path evidence'
    return 'retain_unassigned','no original conventional contribution; evidence retained'


def _register_reads(db,config,source):
    import pysam
    from altanalyze3.components.bam import isoform_structure_extract as ise
    sid=source['source_id'];root=Path(config['output'])
    matrix=root/'compatibility/samples'/sid/(sid+'.h5ad')
    # Identity accounting follows the native extractor exactly. Genomic annotation
    # is evaluated in input order, including its original eligibility conditions.
    exon_coordinates,_,_=ise.gff_process.importEnsemblGenes(config['exon_model'])
    has_matrix=matrix.exists()
    legacy_counts=collections.Counter()
    alignment_indices=collections.Counter()
    with pysam.AlignmentFile(source['bam'],'rb') as bam:
        for read in bam.fetch(until_eof=True):
            q=read.query_name
            if not q:raise ValueError('BAM record has no query name')
            idx=alignment_indices[q];alignment_indices[q]+=1
            db.execute('INSERT INTO segments VALUES(?,?,?,?,?,?,?,?,?,?,?)',
                (sid,q,idx,read.reference_name,read.reference_start,read.reference_end,int(read.is_reverse),
                 read.flag,read.mapping_quality,read.cigarstring,read.get_tag('SA') if read.has_tag('SA') else ''))
            bc=ise.resolve_barcode(read,['CB']) or ''
            molecule=ise.sanitize_gff_value(ise.resolve_molecule_id(read,'zm'))
            umi=next((str(read.get_tag(t)) for t in ('UB','XM','UR') if read.has_tag(t)),'')
            primary=not (read.is_secondary or read.is_supplementary or read.is_unmapped)
            old=db.execute('SELECT primary_seen FROM reads WHERE source_id=? AND query_name=?',(sid,q)).fetchone()
            if old and old[0] and primary:raise ValueError(f'Duplicate primary query: {sid}/{q}')
            feature='';status='nonprimary_or_unmapped'
            if primary:
                status='not_native_eligible'
                if has_matrix and read.mapping_quality>=1 and ise.has_splice(read.cigartuples) and bc:
                    exons=ise.extract_exons(read)
                    if len(exons)>=2:
                        chrom=ise.normalize_chrom(read.reference_name,exon_coordinates)
                        strand=ise.get_read_strand(read)
                        gene,_,_,genes=ise.gff_process.exonAnnotate(chrom,list(exons),strand,q)
                        if genes and ise.has_known_splice_site(chrom,strand,exons,exon_coordinates):
                            core=molecule.split('molecule/',1)[1] if molecule.startswith('molecule/') else molecule
                            feature=ise.sanitize_gff_value(f'{gene}:{core}')
                            legacy_counts[bc,feature]+=1;status='native_eligible'
            metadata=dict(flag=read.flag,reference=read.reference_name,reference_start=read.reference_start,
                          cigar=read.cigarstring,mapq=read.mapping_quality)
            if old is None:
                db.execute('INSERT INTO reads VALUES(?,?,?,?,?,?,?,?,?)',
                           (sid,q,bc,molecule,umi,int(primary),feature,status,json.dumps(metadata)))
            elif primary:
                db.execute('UPDATE reads SET barcode=?,molecule_id=?,umi=?,primary_seen=1,legacy_feature=?,legacy_status=?,metadata_json=? WHERE source_id=? AND query_name=?',
                           (bc,molecule,umi,feature,status,json.dumps(metadata),sid,q))
    db.commit()
    if has_matrix:
        import anndata as ad
        import numpy as np
        a=ad.read_h5ad(matrix);coo=a.X.tocoo()
        observed={(str(a.obs_names[i]),str(a.var_names[j])):int(v) for i,j,v in zip(coo.row,coo.col,coo.data)}
        if dict(legacy_counts)!=observed:
            differences=[(k,legacy_counts[k],observed.get(k,0)) for k in legacy_counts.keys()|observed.keys() if legacy_counts[k]!=observed.get(k,0)]
            atomic_json(root/'reconciliation'/f'{sid}.identity_mismatch.json',differences[:100])
            raise ValueError(f'{sid}: native read-to-matrix identity mismatch; retained diagnostic, no subtraction performed')
    return matrix


def build_ledger(config):
    root=Path(config['output']);out=root/'reconciliation';out.mkdir(parents=True,exist_ok=True)
    done=out/'COMPLETE.json'
    if done.exists():return verify_receipt(done,config)
    database=out/'evidence.sqlite'
    if database.exists():raise FileExistsError('Incomplete ledger retained; use a fresh run directory: '+str(database))
    db=sqlite3.connect(database);initialize(db)
    for source in config['sources']:
        db.execute('INSERT INTO sources VALUES(?,?,?,?,?)',(source['source_id'],source['uid'],source['library'],source['group'],json.dumps(source)))
        _register_reads(db,config,source)
    classes={}
    if 'trans-splicing' in config['modalities']:
        for r in read_tsv(root/'trans_splicing/classification/joint_classification.tsv.gz'):
            classes[r['gene_a'],r['gene_b'],r['joint']]=r
    documented={}
    if 'oncofusion' in config['modalities']:
        for r in read_tsv(root/'oncofusion/native_profile/AML_MDS_FUSION_DOCUMENTED_PAIRS.tsv'):
            documented[tuple(sorted((r['gene_1'],r['gene_2'])))]=r
    for source in config['sources']:
        sid=source['source_id']
        if 'circrna' in config['modalities']:
            for i,r in enumerate(read_tsv(root/'circrna/samples'/sid/(sid+'.noncolinear_joints.tsv.gz'))):
                geometry=[r[k] for k in ('chrom','a_strand','b_strand','donor','acceptor')]
                label='circrna_candidate' if r['joint_type']=='back' else 'noncolinear_'+r['joint_type']
                support(db,sid,r['query_name'],'circrna',geometry,r['status'],label,r,i)
        if 'trans-splicing' in config['modalities']:
            for i,r in enumerate(read_tsv(root/'trans_splicing/samples'/sid/(sid+'.cross_chromosome_junctions.tsv.gz'))):
                geometry=[r[k] for k in ('a_chrom','a_strand','donor','b_chrom','b_strand','acceptor')]
                joint=f"{r['a_chrom']}:{r['donor']}:{r['a_strand']} > {r['b_chrom']}:{r['acceptor']}:{r['b_strand']}"
                ga='UNK' if r['gene_a'].startswith('UNK') else r['gene_a']
                gb='UNK' if r['gene_b'].startswith('UNK') else r['gene_b']
                classified=classes.get((ga,gb,joint))
                r['classification']=classified or {'class':'unresolved','reason':'outside legacy detailed-report scope'}
                label=r['classification']['class']
                support(db,sid,r['query_name'],'trans-splicing' if r['cross_chromosome']=='True' else 'cis_edge_in_noncanonical_path',geometry,r['status'],label,r,i)
        if 'oncofusion' in config['modalities']:
            for i,r in enumerate(read_tsv(root/'oncofusion/samples'/sid/(sid+'.fusion_reads.tsv'))):
                for pair in filter(None,r['gene_pairs'].split(';')):
                    genes,kind=pair.rsplit(':',1);a,b=genes.split('>');ann=documented.get(tuple(sorted((a,b))))
                    label='documented_fusion_pair_candidate' if ann else 'fusion_pair_candidate'
                    evidence=dict(r,annotation=ann or {},identity_level='gene_pair_and_observed_path',biological_mechanism='unresolved')
                    support(db,sid,r['read'],'oncofusion',{'gene_pair':sorted((a,b)),'resolution':'gene_pair'},'observed',label,evidence,i)
                from .ptd_evidence import classify as classify_ptd
                ptd=classify_ptd(r)
                for gene in filter(None,r['same_gene_noncolinear'].split(';')):
                    evidence=dict(r,kmt2a_ptd_evidence=ptd) if gene=='KMT2A' else r
                    label='ptd_type_topology_unresolved' if gene=='KMT2A' and ptd else 'same_gene_noncolinear_observed'
                    geometry={k:ptd['first_ptd_type_join'][k] for k in ('chrom','strand','donor','acceptor')} if gene=='KMT2A' and ptd else dict(resolution='observed_noncolinear_path',segments=r['segments'])
                    support(db,sid,r['read'],'ptd',geometry,'observed',label,evidence,i)
            events=root/'oncofusion/native_profile/FLT3_exon14_15_events_all_reads.tsv'
            original_uid=source['uid']+'__'+source['library']
            own={(r['read'],r['kind'],r['ref_pos'],r['length'],r['seq']) for r in read_tsv(root/'oncofusion/samples'/sid/(sid+'.flt3_events.tsv'))}
            all_events=read_tsv(events) if events.stat().st_size else []
            from .onco_profile import itd_groups
            accepted={(r['uid'],r['read'],r['ref_pos'],r['length'],r['seq']) for g in itd_groups(all_events) if g['accepted'] for r in g['rows']}
            seen=set()
            for i,r in enumerate(all_events):
                key=(r['read'],r['kind'],r['ref_pos'],r['length'],r['seq'])
                if r['uid']!=original_uid or key not in own or key in seen:continue
                seen.add(key)
                tandem=r['kind']=='insertion' and float(r.get('tandem_8mer_fraction') or 0)>=0.8
                label='itd_native_group_candidate' if (r['uid'],r['read'],r['ref_pos'],r['length'],r['seq']) in accepted else ('tandem_sequence_match_observed' if tandem else 'insertion_or_softclip_observed')
                support(db,sid,r['read'],'itd',['chr13',r['ref_pos'],r['kind'],r['length'],r['seq']],
                        'observed',label,r,i)
        db.commit()
    for sid,q,feature in db.execute('SELECT source_id,query_name,legacy_feature FROM reads').fetchall():
        evidence=db.execute('SELECT label,status FROM support WHERE source_id=? AND query_name=?',(sid,q)).fetchall()
        labels=sorted({r[0] for r in evidence})
        eligible=any(label in ('trans_splicing_candidate','gene_fusion_candidate') and status=='kept' for label,status in evidence)
        action,reason=route(feature,labels,eligible)
        db.execute('INSERT INTO assignments VALUES(?,?,?,?,?,?)',(sid,q,feature,json.dumps(labels),action,reason))
    db.commit()
    # Preserve feature and cell roster, including all-zero features after reallocation.
    deltas=[]
    for source in config['sources']:
        sid=source['source_id'];path=root/'compatibility/samples'/sid/(sid+'.h5ad')
        if not path.exists():continue
        import anndata as ad
        import numpy as np
        from scipy import sparse
        a=ad.read_h5ad(path);before=float(a.X.sum());x=a.X.tolil(copy=True)
        oi={str(v):i for i,v in enumerate(a.obs_names)};vi={str(v):i for i,v in enumerate(a.var_names)}
        for bc,feature,n in db.execute('''SELECT r.barcode,a.legacy_feature,COUNT(*) FROM assignments a
                JOIN reads r USING(source_id,query_name) WHERE a.source_id=? AND a.action='reallocate'
                GROUP BY r.barcode,a.legacy_feature''',(sid,)):
            i,j=oi[bc],vi[feature]
            if x[i,j]<n:raise ValueError('Reallocation exceeds original contribution')
            x[i,j]-=n;deltas.append(dict(source_id=sid,barcode=bc,legacy_feature=feature,removed_reads=n))
        a.X=x.tocsr();a.layers['legacy_counts']=ad.read_h5ad(path).X.copy()
        a.uns['assignment_view']='integrated; original counts retained in legacy_counts'
        if before-float(a.X.sum())!=sum(d['removed_reads'] for d in deltas if d['source_id']==sid):raise ValueError('Count conservation failed')
        target=root/'integrated/standard/samples'/sid/(sid+'.h5ad');target.parent.mkdir(parents=True,exist_ok=True);a.write_h5ad(target)
    write_tsv(out/'count_deltas.tsv',['source_id','barcode','legacy_feature','removed_reads'],deltas)
    export_query(db,out/'reads.tsv.gz','SELECT * FROM reads ORDER BY source_id,query_name')
    export_query(db,out/'segments.tsv.gz','SELECT * FROM segments ORDER BY source_id,query_name,alignment_index')
    export_query(db,out/'events.tsv.gz','SELECT * FROM events ORDER BY event_id')
    export_query(db,out/'read_event_support.tsv.gz','SELECT * FROM support ORDER BY source_id,query_name,event_id,row_index')
    export_query(db,out/'assignments.tsv.gz','SELECT * FROM assignments ORDER BY source_id,query_name')
    export_query(db,out/'translation_caveats.tsv.gz', '''SELECT source_id,query_name,legacy_feature,labels_json,
        'reference exon-chain translation does not resolve the observed insertion, duplicated path or circular topology' AS caveat
        FROM assignments WHERE labels_json!='[]' ORDER BY source_id,query_name''')
    export_query(db,out/'event_source_counts.tsv.gz', '''SELECT e.modality,s.event_id,s.source_id,s.status,
        COUNT(DISTINCT s.query_name) AS distinct_source_reads,COUNT(*) AS evidence_rows
        FROM support s JOIN events e USING(event_id) GROUP BY e.modality,s.event_id,s.source_id,s.status
        ORDER BY e.modality,s.event_id,s.source_id,s.status''')
    export_query(db,out/'event_cell_counts.tsv.gz', '''SELECT e.modality,s.event_id,s.source_id,r.barcode,s.status,
        COUNT(DISTINCT s.query_name) AS distinct_source_reads FROM support s JOIN events e USING(event_id)
        JOIN reads r USING(source_id,query_name) GROUP BY e.modality,s.event_id,s.source_id,r.barcode,s.status
        ORDER BY e.modality,s.event_id,s.source_id,r.barcode,s.status''')
    summary=dict(source_count=len(config['sources']),unique_source_reads=db.execute('SELECT COUNT(*) FROM reads').fetchone()[0],
                 events=db.execute('SELECT COUNT(*) FROM events').fetchone()[0],evidence_rows=db.execute('SELECT COUNT(*) FROM support').fetchone()[0],
                 actions=dict(db.execute('SELECT action,COUNT(*) FROM assignments GROUP BY action')),
                 removed_standard_reads=sum(d['removed_reads'] for d in deltas),
                 counting_unit='distinct source_id + query_name within each event; overlapping modalities must not be summed',
                 topology='candidate labels do not establish physical circularity or genomic rearrangement')
    if sum(summary['actions'].values())!=summary['unique_source_reads']:raise ValueError('Read partition failed')
    atomic_json(out/'SUMMARY.json',summary);db.close()
    files=[p for p in out.iterdir() if p.is_file()]+list((root/'integrated/standard/samples').glob('*/*.h5ad'))
    receipt(done,config,{},files,status='COMPLETE',summary=summary)
    return summary
