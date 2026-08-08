"""Marker prefilter: how much smaller, and does it change the answer (ARI)?"""
import os, sys, json, time, shutil, subprocess
HERE=os.path.dirname(os.path.abspath(__file__)); REPO=os.path.abspath(os.path.join(HERE,'..','..','..','..'))
PY=sys.executable
CHILD='''
import os,sys,json,time,resource
def main():
    n_per,h5ad,outdir,keys=sys.argv[1:5]
    n_per=int(n_per); sys.path.insert(0,{repo!r})
    import matplotlib; matplotlib.use('Agg')
    import scanpy as sc
    from altanalyze3.components.sctriangulate import ScTriangulate
    adata=sc.read(h5ad); t0=time.perf_counter()
    s=ScTriangulate(dir=outdir,adata=adata,query=keys.split(','),predict_doublet=False)
    info=None
    if n_per>0: info=s.prefilter_marker_genes(n_per_cluster=n_per,max_cells_per_cluster=1000)
    s.lazy_run(compute_metrics_parallel=False,compute_shapley_parallel=False)
    obs=s.adata.obs
    json.dump({{'n_per':n_per,'seconds':time.perf_counter()-t0,
        'peak_gb':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9,'info':info,
        'pruned':list(obs['pruned'].astype(str)),'raw':list(obs['raw'].astype(str)),
        'final_annotation':list(obs['final_annotation'].astype(str))}},
        open(os.path.join(outdir,'r.json'),'w'))
    print('@@OK@@')
if __name__=='__main__': main()
'''
def run(n_per,h5ad,keys):
    outdir=os.path.join(HERE,'out_prefilter_%d'%n_per)
    if os.path.isdir(outdir): shutil.rmtree(outdir)
    os.makedirs(outdir)
    sp=os.path.join(outdir,'_c.py'); open(sp,'w').write(CHILD.format(repo=REPO))
    env=dict(os.environ,PYTHONHASHSEED='1000'); env.setdefault('OMP_NUM_THREADS','8')
    o=subprocess.run([PY,sp,str(n_per),h5ad,outdir,','.join(keys)],capture_output=True,text=True,env=env,cwd=HERE)
    if '@@OK@@' not in o.stdout: print(o.stdout[-2500:]);print(o.stderr[-2500:]); raise SystemExit(1)
    return json.load(open(os.path.join(outdir,'r.json')))
if __name__=='__main__':
    from sklearn.metrics import adjusted_rand_score as ari
    h5ad=os.path.abspath(sys.argv[1]); keys=sys.argv[2].split(',')
    base=run(0,h5ad,keys)
    print('  %-12s %8s %9s %8s %9s %9s %9s'%('top/cluster','genes','seconds','peak GB','ARI final','ARI raw','ARI pruned'))
    print('  '+'-'*72)
    print('  %-12s %8d %9.1f %8.2f %9s %9s %9s'%('all (none)',len(base['pruned']) and 32738,base['seconds'],base['peak_gb'],'-','-','-'))
    rows=[]
    for n in [int(x) for x in (sys.argv[3].split(',') if len(sys.argv)>3 else ['500','200','50'])]:
        r=run(n,h5ad,keys); i=r['info']
        a=[ari(base[c],r[c]) for c in ['final_annotation','raw','pruned']]
        print('  %-12d %8d %9.1f %8.2f %9.4f %9.4f %9.4f'%(n,i['genes_after'],r['seconds'],r['peak_gb'],a[0],a[1],a[2]))
        rows.append({'n_per_cluster':n,'genes':i['genes_after'],'seconds':r['seconds'],
                     'peak_gb':r['peak_gb'],'ari_final':a[0],'ari_raw':a[1],'ari_pruned':a[2]})
    json.dump({'baseline':{k:v for k,v in base.items() if k not in ('pruned','raw','final_annotation')},'rows':rows},
              open(os.path.join(HERE,'prefilter_validate.json'),'w'),indent=1)
