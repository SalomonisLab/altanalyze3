"""Independent gene partitions with complete-source gates and ordered merging."""
from collections import Counter
from pathlib import Path
import gzip,json,multiprocessing,shutil,sqlite3
from .stream_context import build_contexts

_STATE=None

def _partition(task):
    index,genes,counts,out=task
    scans,regions,references,policy=_STATE
    return build_contexts(scans,regions,references,Path(out)/f'partition_{index:03d}',policy,
                          _gene_partition=genes,_partition_read_counts=counts)


def build_parallel_contexts(scans,regions,references,outdir,policy,*,workers=8):
    """Run identical per-gene calculations, retaining every complete input row."""
    global _STATE
    if workers<1:raise ValueError('Positive worker count required')
    if not all(s['complete'] for s in scans):raise ValueError('Complete source scans required')
    libs=[s['library'] for s in scans]
    if len(set(libs))!=len(libs):raise ValueError('Unique physical source identities required')
    # Independent full SQL totals establish each partition's required coverage.
    gene_counts={}
    for i,s in enumerate(scans):
        c=sqlite3.connect(f"file:{s['database']}?mode=ro",uri=True)
        try:
            total=0
            for gene,count in c.execute('SELECT gene,SUM(count) FROM models GROUP BY gene'):
                gene_counts.setdefault(gene,[0]*len(scans))[i]=count;total+=count
            if total!=s['native_eligible_reads']:raise ValueError('Independent full source total differs')
        finally:c.close()
    ordered=sorted(gene_counts);parts=[];counts=[];chunk=[];current=[0]*len(scans)
    target=sum(s['native_eligible_reads'] for s in scans)/workers
    for gene in ordered:
        if chunk and sum(current)>=target and len(parts)<workers-1:
            parts.append(chunk);counts.append(current);chunk=[];current=[0]*len(scans)
        chunk.append(gene);current=[a+b for a,b in zip(current,gene_counts[gene])]
    if chunk:parts.append(chunk);counts.append(current)
    flat=[g for p in parts for g in p]
    if flat!=ordered or len(set(flat))!=len(flat):raise ValueError('Complete disjoint gene partition gate failed')
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    _STATE=scans,regions,references,policy
    print('FULL_CONTEXT_PARTITIONS_BEGIN',len(parts),'workers',len(ordered),'genes',flush=True)
    try:
        with multiprocessing.get_context('fork').Pool(min(workers,len(parts))) as pool:
            results=pool.map(_partition,[(i,p,n,str(out)) for i,(p,n) in enumerate(zip(parts,counts))],chunksize=1)
    finally:_STATE=None
    if [g for r in results for g in r['partition_genes']]!=ordered:raise ValueError('Returned partition gene roster differs')
    complete_counts=[sum(r['source_read_counts'][i] for r in results) for i in range(len(scans))]
    if complete_counts!=[s['native_eligible_reads'] for s in scans]:raise ValueError('Complete ordered per-source count gate failed')
    contexts={};override_counts=Counter()
    for r in results:
        if r['ordered_libraries']!=libs or set(contexts)&set(r['gene_evidence']):raise ValueError('Partition identity collision')
        contexts.update(r['gene_evidence']);override_counts.update(r['override_counts'])
    paths=[];baseline=[];overrides={}
    for i,lib in enumerate(libs):
        folder=out/f'sample_{i:04d}';folder.mkdir()
        for field,name in [('evidence_sample_paths','evidence.exemplars.tsv'),('baseline_sample_paths','baseline.exemplars.tsv')]:
            path=folder/name
            with path.open('wb') as f:
                for j,r in enumerate(results):
                    if r[field][i][0]!=lib:raise ValueError('Physical source order changed')
                    with open(r[field][i][1],'rb') as src:
                        header=src.readline()
                        if header!=b'gene\tstrand\tstructure\tcount\trep_molecule_id\n':raise ValueError('Invalid partition table header')
                        if j==0:f.write(header)
                        shutil.copyfileobj(src,f,16*1024**2)
            (paths if field=='evidence_sample_paths' else baseline).append((lib,str(path)))
        path=folder/'molecule_classes.tsv.gz'
        # Concatenate gzip members; all readers transparently read all members.
        # Strip subsequent headers by streaming recompression.
        with gzip.open(path,'wb') as f:
            for j,r in enumerate(results):
                with gzip.open(r['override_paths'][lib],'rb') as src:
                    header=src.readline()
                    if header!=b'gene\tmolecule\toriginal_structure\tcandidate_class\n':raise ValueError('Invalid partition override header')
                    if j==0:f.write(header)
                    shutil.copyfileobj(src,f,16*1024**2)
        overrides[lib]=str(path)
    result=dict(ordered_libraries=libs,all_native_reads=sum(complete_counts),evidence_genes=len(contexts),
                protected_reads=sum(override_counts.values()),policy=policy.__dict__,gene_evidence=contexts,
                evidence_sample_paths=paths,baseline_sample_paths=baseline,override_paths=overrides,
                override_counts=dict(override_counts),source_read_counts=complete_counts,
                complete_partition_gene_identities=ordered,parallel_workers=workers,
                all_partition_identities_and_source_counts_verified=True)
    if result['protected_reads']!=sum(r['protected_reads'] for r in results):raise ValueError('Protected count merge changed')
    (out/'complete_context_manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print('FULL_CONTEXT_PARTITIONS_COMPLETE',len(ordered),'genes',result['all_native_reads'],'reads',flush=True)
    return result
