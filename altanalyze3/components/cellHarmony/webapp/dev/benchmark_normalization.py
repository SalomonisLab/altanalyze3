"""Compare normalization candidates in isolated processes on a saved RNA matrix.

The baseline function can be snapshotted before editing cellHarmony_lite.py.
Both modes load counts as X, retain an independent counts layer, and hash all
CSR buffers after normalization. This is a normalization benchmark, not a full
analysis timing. Pass --baseline-source for the original function source file.
"""
import argparse
import ast
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import anndata as ad
import h5py


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('h5ad',type=Path)
    parser.add_argument('--baseline-source',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--rounds',type=int,default=3)
    args=parser.parse_args()
    from altanalyze3.components.cellHarmony import cellHarmony_lite as module
    try:
        from anndata.io import read_elem
    except ImportError:
        from anndata.experimental import read_elem
    with h5py.File(args.h5ad,'r') as fh:
        matrix=read_elem(fh['layers/counts'] if 'counts' in fh.get('layers',{}) else fh['X'])
    obj=ad.AnnData(matrix)
    obj.layers['counts']=matrix.copy()
    normalize=module.normalize_adata
    if args.baseline_source:
        node=next(n for n in ast.parse(args.baseline_source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='normalize_adata')
        namespace=dict(vars(module))
        exec(compile(ast.Module(body=[node],type_ignores=[]),str(args.baseline_source),'exec'),namespace)
        normalize=namespace['normalize_adata']
    before=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    times=[]
    for trial in range(args.rounds):
        if trial:
            obj.X.data[:] = obj.layers['counts'].data
            obj.uns.pop('log1p', None)
        start=time.perf_counter()
        normalize(obj)
        times.append(time.perf_counter()-start)
    elapsed=sum(times)/len(times)
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    hashes={}
    for name,matrix in [('X',obj.X),('counts',obj.layers['counts'])]:
        digest=hashlib.sha256()
        for buffer in (matrix.data,matrix.indices,matrix.indptr):
            digest.update(memoryview(buffer).cast('B'))
        hashes[name]=digest.hexdigest()
    divisor=1024**3 if sys.platform=='darwin' else 1024**2
    result=dict(mode='baseline' if args.baseline_source else 'bounded',shape=obj.shape,nnz=obj.X.nnz,
                elapsed_seconds=elapsed,round_seconds=times,peak_rss_gib=peak/divisor,pre_normalization_peak_gib=before/divisor,hashes=hashes)
    args.output.write_text(json.dumps(result,indent=2))
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    main()
