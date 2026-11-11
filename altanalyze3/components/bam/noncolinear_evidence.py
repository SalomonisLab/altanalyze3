#!/usr/bin/env python3
"""Lossless count correction for noncolinear_rna pair tables.

Raw pair tables remain authoritative and unchanged. Every class and status is retained.
Stable event ID uses the baseline normalized oriented donor/acceptor coordinates.
Counts distinguish pair observations, query names, declared molecule IDs and cell/UMI
identities. A molecule ID is not asserted to be an original RNA molecule. CB/UMI
counts are descriptive only, not a substituted scientific counting baseline.
"""
import argparse, csv, gzip, hashlib, json, pathlib, shutil, sqlite3, tempfile

EVENT=('chrom','a_strand','b_strand','donor','acceptor')
FLAGS=('contained','multipass','polyA','homology_above_6','positive_query_gap','unresolved')

def event_id(r):
    return 'BSJ|'+'|'.join(str(r[k]) for k in EVENT) if r['joint_type']=='back' else 'NCJ|'+'|'.join(str(r[k]) for k in EVENT)

def flags(r):
    back=r['joint_type']=='back'
    return (int(back and r['outside_circle_nt']=='0'),int(back and (int(r['pair_overlap_nt'] or 0)>0 or int(r['crossings'] or 0)>=2)),int(r['polyA_at_joint']=='True'),int(int(r['junction_homology_nt'])>6),int(int(r['query_gap'])>0),int(r['resolved']!='True'))

def initialize(db):
    db.execute('CREATE TABLE evidence(event TEXT, library TEXT,status TEXT,query_name TEXT,molecule_id TEXT,barcode TEXT,umi TEXT,contained INT,multipass INT,polyA INT,homology_above_6 INT,positive_query_gap INT,unresolved INT)')
    db.execute('CREATE TABLE events(event TEXT PRIMARY KEY, chrom TEXT, donor_strand TEXT,acceptor_strand TEXT,donor INT,acceptor INT,joint_type TEXT,donor_annotation TEXT,acceptor_annotation TEXT,donor_region TEXT,acceptor_region TEXT)')
    db.execute('CREATE TABLE classes(event TEXT,library TEXT,status TEXT,joint_class TEXT,pair_rows INT,PRIMARY KEY(event,library,status,joint_class))')

def ingest(db,path,expected_library=None):
    n=0
    with gzip.open(path,'rt') as fh:
      for r in csv.DictReader(fh,delimiter='\t'):
        if expected_library is not None and r['library']!=expected_library: raise ValueError('library identity mismatch: '+r['library'])
        eid=event_id(r)
        db.execute('INSERT OR IGNORE INTO events VALUES (?,?,?,?,?,?,?,?,?,?,?)',(eid,r['chrom'],r['a_strand'],r['b_strand'],int(r['donor']),int(r['acceptor']),r['joint_type'],r['donor_annotation'],r['acceptor_annotation'],r['donor_region'],r['acceptor_region']))
        db.execute('INSERT INTO evidence VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)',(eid,r['library'],r['status'],r['query_name'],r['molecule_id'],r['barcode'],r['umi'],*flags(r)))
        db.execute('INSERT INTO classes VALUES (?,?,?,?,1) ON CONFLICT(event,library,status,joint_class) DO UPDATE SET pair_rows=pair_rows+1',(eid,r['library'],r['status'],r['joint_class']))
        n+=1
        if n%100000==0: db.commit()
    db.commit();return n

def write_query(db,path,sql):
    cur=db.execute(sql)
    with gzip.open(path,'wt') as fh:
      writer=csv.writer(fh,delimiter='\t');writer.writerow([x[0] for x in cur.description]);writer.writerows(cur)

def finalize(db,out):
    db.execute('CREATE INDEX evidence_read ON evidence(event,library,status,query_name)')
    write_query(db,out/'events.tsv.gz','SELECT * FROM events ORDER BY event')
    write_query(db,out/'classes.tsv.gz','SELECT * FROM classes ORDER BY event,library,status,joint_class')
    # Never sum per-class counts to obtain unique reads. Multiple rows/crossings may support one event.
    write_query(db,out/'event_library_status_counts.tsv.gz', '''SELECT event,library,status,COUNT(*) AS raw_pair_rows,
      COUNT(DISTINCT query_name) AS unique_read_count,
      COUNT(DISTINCT NULLIF(molecule_id,'')) AS unique_declared_molecule_id_count,
      COUNT(DISTINCT CASE WHEN barcode<>'' AND umi<>'' THEN json_array(barcode,umi) END) AS unique_cell_umi_count,
      COUNT(DISTINCT NULLIF(barcode,'')) AS barcode_count,
      COUNT(DISTINCT CASE WHEN umi='' THEN query_name END) AS reads_missing_umi
      FROM evidence GROUP BY event,library,status ORDER BY event,library,status''')
    # The support ledger retains original query/molecule/barcode/UMI identities, all statuses,
    # multiplicity and QC evidence; query-name conflicts remain separate rows rather than disappearing.
    write_query(db,out/'read_support.tsv.gz','SELECT event,library,status,query_name,molecule_id,barcode,umi,COUNT(*) AS raw_pair_rows,'+
      ','.join('MAX('+f+') AS '+f for f in FLAGS)+' FROM evidence GROUP BY event,library,status,query_name,molecule_id,barcode,umi ORDER BY event,library,status,query_name')
    write_query(db,out/'read_identity_conflicts.tsv.gz', '''SELECT library,query_name,COUNT(DISTINCT json_array(molecule_id,barcode,umi)) AS identity_variants
      FROM evidence GROUP BY library,query_name HAVING identity_variants>1''')

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as fh:
      for chunk in iter(lambda:fh.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=pathlib.Path,required=True);p.add_argument('--out',type=pathlib.Path,required=True)
    p.add_argument('--scratch',type=pathlib.Path,help='Node-local parent for SQLite and sort/export artifacts; copied to out only after successful checks')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    paths=sorted(a.root.glob('samples/*/*.noncolinear_joints.tsv.gz'))
    if not paths:raise ValueError('No input tables')
    validation=None
    if (a.root/'RUN_MANIFEST.json').exists():
        # Production all-source runs cannot silently fall back to a completed subset.
        validation=json.loads((a.root/'VALIDATED_EVIDENCE.json').read_text())
        run=json.loads((a.root/'RUN_MANIFEST.json').read_text())
        expected=run['expected_ordered_libraries']
        if validation['manifest_sha256']!=sha256(a.root/'RUN_MANIFEST.json'):raise ValueError('Validation/run version mismatch')
        if validation['expected_ordered_libraries']!=expected:raise ValueError('Validation roster mismatch')
        if sorted(path.parent.name for path in paths)!=sorted(expected):raise ValueError('Exact full raw library roster mismatch')
        validated={r['library']:r for r in validation['libraries']}
        if len(validated)!=len(expected) or set(validated)!=set(expected):raise ValueError('Duplicate or missing validation libraries')
        paths=[a.root/'samples'/lib/(lib+'.noncolinear_joints.tsv.gz') for lib in expected]
    scratch=pathlib.Path(tempfile.mkdtemp(prefix='circrna_evidence_',dir=a.scratch)) if a.scratch else a.out
    db=sqlite3.connect(scratch/'evidence.sqlite');db.execute('PRAGMA temp_store=FILE');initialize(db);manifest=[]
    for path in paths:
        lib=path.parent.name;before=(path.stat().st_size,path.stat().st_mtime_ns);hashed=sha256(path)
        if validation and hashed!=validated[lib]['files']['raw']['sha256']:raise ValueError('Raw evidence changed after validation '+lib)
        n=ingest(db,path,lib)
        if validation and n!=validated[lib]['pair_rows']:raise ValueError('Raw pair retention failed '+lib)
        if (path.stat().st_size,path.stat().st_mtime_ns)!=before:raise ValueError('Raw evidence changed during aggregation '+lib)
        manifest.append(dict(path=str(path),sha256=hashed,rows=n));print('Aggregated',lib,n,flush=True)
    if db.execute('SELECT COUNT(*) FROM evidence').fetchone()[0]!=sum(x['rows'] for x in manifest):raise AssertionError('Raw pair storage conservation failed')
    finalize(db,scratch);db.close()
    artifacts={}
    for path in sorted(scratch.iterdir()):
        hashed=sha256(path);artifacts[path.name]=dict(sha256=hashed,size=path.stat().st_size)
        if scratch!=a.out:
            target=a.out/path.name;shutil.copy2(path,target)
            if sha256(target)!=hashed:raise AssertionError('Output copy checksum mismatch '+path.name)
    if validation and sha256(a.root/'RUN_MANIFEST.json')!=validation['manifest_sha256']:raise ValueError('Run manifest changed during aggregation')
    (a.out/'MANIFEST.json').write_text(json.dumps(dict(inputs=manifest,outputs=artifacts,code_sha256=sha256(__file__),
      validated_evidence_sha256=sha256(a.root/'VALIDATED_EVIDENCE.json') if validation else None,
      counting_unit='Distinct query name within event + library + status; separate declared molecule-id and CB/UMI descriptive counts',
      method_change='Correct pair-row/molecule label bug; no input rows/classes/statuses excluded',
      scientific_status='Descriptive evidence only; biological validation and statistical baseline gates remain external'),indent=2)+'\n')
    if scratch!=a.out:shutil.rmtree(scratch)

if __name__=='__main__':main()
