"""Compare only UMAP projection search, with explicit input and method gates.

Use --input for a fully hashed, independently verified prepared marker matrix,
not an unverified or reconstructed visitor input. --synthetic is a computational
stress test. Fit/landmarks/features/metric/transform optimization are identical.
Run in a fresh process with NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4.
"""
from __future__ import annotations
import argparse, json, time, resource, sys, hashlib
from pathlib import Path
from importlib.metadata import version
import numpy as np
from threadpoolctl import threadpool_limits
from altanalyze3.components.clustering.umap_fit import select_landmarks, projection_search
from altanalyze3.components.clustering.umap_neighbors import ExactCorrelationIndex
from altanalyze3.components.clustering.benchmarking.benchmark_landmark_umap import load_verified_input
from altanalyze3.components.clustering.benchmarking.benchmark_umap_separation import separation


def main():
    parser=argparse.ArgumentParser()
    source=parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--input',type=Path)
    source.add_argument('--synthetic',type=int,help='Total computational-fixture cells')
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--baseline-query-cells',type=int,default=0,help='0 compares the full roster')
    args=parser.parse_args()
    # This is the actual human decision, not a generic proceed or self-approval.
    authorization=json.loads((Path(__file__).parent/'discover_umap_transform_profile_20261007.json').read_text())['decision']
    if authorization['state']!='approved' or 'Yes, do that' not in authorization['source']:
        raise ValueError('Explicit approval of the proposed search-only comparison is required.')
    threadpool_limits(4)
    import numba
    numba.set_num_threads(4)
    from umap import UMAP
    parameters={'n_neighbors':15,'min_dist':.75,'metric':'correlation','random_state':0}
    if args.input:
        x,cells,features,states,manifest=load_verified_input(args.input)
        parameters=manifest['umap_parameters']
        provenance={'kind':'verified prepared input','sha256':manifest['sha256'],'source_sha256':manifest['source_sha256'],'parameters':parameters}
    else:
        rng=np.random.default_rng(20261007)
        features=[f'synthetic_feature_{i}' for i in range(827)]
        cells=[f'synthetic_cell_{i}' for i in range(args.synthetic)]
        states=(np.arange(args.synthetic)%51).astype(str)
        centers=rng.random((51,len(features)),dtype=np.float32);centers[centers<.85]=0
        x=np.empty((len(cells),len(features)),dtype=np.float32)
        for start in range(0,len(x),4096):
            block=centers[np.arange(start,min(start+4096,len(x)))%51].copy()
            block+=rng.random(block.shape,dtype=np.float32)*.2
            block[rng.random(block.shape,dtype=np.float32)<.65]=0
            np.log1p(block*20,out=block);x[start:start+len(block)]=block
        provenance={'kind':'synthetic computational fixture','seed':20261007,'parameters':parameters,'matrix_sha256':hashlib.sha256(x.tobytes()).hexdigest()}
    assert x.shape==(len(cells),len(features)) and len(states)==len(cells)
    selected=select_landmarks(states,30000,0,200)
    remaining=np.setdiff1d(np.arange(len(x)),selected)
    assert all(np.sum(states[selected]==s)>=min(200,np.sum(states==s)) for s in np.unique(states))
    model=UMAP(**parameters)
    t=time.perf_counter();fit=model.fit_transform(x[selected]);fit_seconds=time.perf_counter()-t
    original_index=model._knn_search_index
    print(json.dumps({'fit_seconds':fit_seconds,'cells':len(x),'features':len(features),'landmarks':len(selected)}),flush=True)
    matched=remaining if not args.baseline_query_cells else remaining[:args.baseline_query_cells]
    out={};outputs={};baseline_coords=None
    for backend in ('umap','exact_correlation'):
        rows=matched if backend=='umap' else remaining
        coords=np.empty((len(x),2),dtype=np.float32);coords[selected]=fit
        started=time.perf_counter();batches=[]
        with projection_search(model,backend) as effective:
            for block in np.array_split(rows,max(1,int(np.ceil(len(rows)/50000)))):
                t=time.perf_counter();coords[block]=model.transform(x[block]);batches.append({'cells':len(block),'seconds':time.perf_counter()-t})
                print(json.dumps({'backend':backend,'batch':batches[-1]}),flush=True)
        assert model._knn_search_index is original_index
        assert np.array_equal(coords[selected],fit) and np.isfinite(coords[np.r_[selected,rows]]).all()
        result={'query_cells':len(rows),'seconds':time.perf_counter()-started,'batches':batches,'fit_coordinates_unchanged':True,'finite_coordinates':len(selected)+len(rows),'search_backend':effective}
        if len(rows)==len(remaining):
            result['separation']=separation(coords,states)
        if backend=='umap':
            baseline_coords=coords if len(rows)==len(remaining) else None
        out[backend]=result;outputs[backend]=coords
        print(json.dumps({backend:{k:v for k,v in result.items() if k!='separation'}}),flush=True)
    # Exact independent float64 Pearson reference uses every feature, all fit
    # cells and a diagnostic anchor from every state. No analytical panel shrinks.
    anchors=np.concatenate([np.flatnonzero(states[remaining]==s)[:5] for s in np.unique(states[remaining])])
    anchors=remaining[anchors]
    train=np.asarray(x[selected],dtype=np.float64);train-=train.mean(axis=1,keepdims=True);norm=np.linalg.norm(train,axis=1);tc=norm==0;train/=np.where(tc,1,norm)[:,None]
    query=np.asarray(x[anchors],dtype=np.float64);query-=query.mean(axis=1,keepdims=True);norm=np.linalg.norm(query,axis=1);qc=norm==0;query/=np.where(qc,1,norm)[:,None]
    distances=1-query@train.T;distances[np.ix_(qc,tc)]=0
    k=int(parameters['n_neighbors']);reference=np.argsort(distances,axis=1)[:,:k]
    baseline_neighbors=original_index.query(x[anchors],k=k,epsilon=.12)[0]
    with projection_search(model,'exact_correlation'):
        neighbors,d=model._knn_search_index.query(x[anchors],k=k)
    recall=lambda ix: float(np.mean([len(set(a)&set(b))/k for a,b in zip(ix,reference)]))
    out['neighbor_validation']={'anchors':len(anchors),'all_features':len(features),'baseline_exact_recall':recall(baseline_neighbors),'candidate_exact_recall':recall(neighbors),'max_distance_error':float(np.max(np.abs(d-np.take_along_axis(distances,neighbors,axis=1))))}
    out.update(cells=len(cells),features=len(features),states=len(np.unique(states)),fit_seconds=fit_seconds,landmarks=len(selected),authorization=authorization,provenance=provenance,
               peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**3 if sys.platform=='darwin' else 1024**2),
               versions={p:version(p) for p in ('numpy','scipy','numba','umap-learn','pynndescent','scikit-learn')},
               scope='UMAP-only benchmark; does not establish full-workflow or visitor-job performance')
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'report.json').write_text(json.dumps(out,indent=2)+'\n')
    for backend,coords in outputs.items():
        if backend!='umap' or baseline_coords is not None:np.save(args.output/f'{backend}.npy',coords)
    np.save(args.output/'selected.npy',selected)
    print(json.dumps({k:v for k,v in out.items() if k not in ('umap','exact_correlation','provenance')}),flush=True)

if __name__=='__main__':main()
