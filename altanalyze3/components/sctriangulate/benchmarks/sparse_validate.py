"""Sparse-only compute path vs the densified path: labels, ARI, time, memory."""
import os, sys, json, time, shutil, resource, subprocess
HERE=os.path.dirname(os.path.abspath(__file__)); REPO=os.path.abspath(os.path.join(HERE,'..','..','..','..'))
PY=sys.executable
CHILD='''
import os,sys,json,time,resource
def main():
    mode,h5ad,outdir,keys=sys.argv[1:5]
    sys.path.insert(0,{repo!r})
    import matplotlib; matplotlib.use('Agg')
    import scanpy as sc
    from altanalyze3.components.sctriangulate import ScTriangulate
    from altanalyze3.components.sctriangulate import main_class as MC
    adata=sc.read(h5ad); t0=time.perf_counter()
    s=ScTriangulate(dir=outdir,adata=adata,query=keys.split(','),predict_doublet=False)
    s.sparse_compute=(mode=='sparse')
    s.lazy_run(compute_metrics_parallel=False,compute_shapley_parallel=False)
    obs=s.adata.obs
    json.dump({{'mode':mode,'seconds':time.perf_counter()-t0,
        'peak_gb':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9,
        'pruned':list(obs['pruned'].astype(str)),'raw':list(obs['raw'].astype(str)),
        'final_annotation':list(obs['final_annotation'].astype(str))}},
        open(os.path.join(outdir,'r.json'),'w'))
    print('@@OK@@')
if __name__=='__main__': main()
'''
def run(mode,h5ad,keys):
    outdir=os.path.join(HERE,'out_sparse_'+mode)
    if os.path.isdir(outdir): shutil.rmtree(outdir)
    os.makedirs(outdir)
    sp=os.path.join(outdir,'_c.py'); open(sp,'w').write(CHILD.format(repo=REPO))
    env=dict(os.environ,PYTHONHASHSEED='1000'); env.setdefault('OMP_NUM_THREADS','8')
    o=subprocess.run([PY,sp,mode,h5ad,outdir,','.join(keys)],capture_output=True,text=True,env=env,cwd=HERE)
    if '@@OK@@' not in o.stdout: print(o.stdout[-2500:]);print(o.stderr[-2500:]); raise SystemExit(1)
    return json.load(open(os.path.join(outdir,'r.json')))
if __name__=='__main__':
    from sklearn.metrics import adjusted_rand_score as ari, normalized_mutual_info_score as nmi
    h5ad=os.path.abspath(sys.argv[1]); keys=sys.argv[2].split(',')
    d=run('dense',h5ad,keys); s=run('sparse',h5ad,keys)
    print('  dense  path: %7.1f s  peak %5.2f GB'%(d['seconds'],d['peak_gb']))
    print('  sparse path: %7.1f s  peak %5.2f GB   (%.2fx time, %.2fx memory)'%(
        s['seconds'],s['peak_gb'],d['seconds']/s['seconds'],d['peak_gb']/s['peak_gb']))
    print('  %-20s %14s %9s %9s'%('COLUMN','EXACT','ARI','NMI'))
    for c in ['final_annotation','raw','pruned']:
        n=sum(1 for x,y in zip(d[c],s[c]) if x==y)
        print('  %-20s %7d/%-6d %9.4f %9.4f'%(c,n,len(d[c]),ari(d[c],s[c]),nmi(d[c],s[c])))
