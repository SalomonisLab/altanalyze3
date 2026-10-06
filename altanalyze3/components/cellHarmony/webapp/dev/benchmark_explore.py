"""Read-only, fresh-process reproduction of large saved-job Explore requests.

Baseline emulates the float64 bundle rejection plus legacy all-view payloads.
Optimized uses the shared routes and lossless compact transport. No analysis is rerun.
"""
import argparse
import importlib
import json
import resource
import time
from pathlib import Path

from fastapi.testclient import TestClient
import h5py
try:
    from anndata.io import read_elem
except ImportError:
    from anndata.experimental import read_elem


def main():
    p = argparse.ArgumentParser()
    p.add_argument('job_json', type=Path)
    p.add_argument('--mode', choices=['baseline','optimized'], required=True)
    p.add_argument('--report', type=Path, required=True)
    args = p.parse_args()
    W = importlib.import_module('altanalyze3.components.cellHarmony.webapp.app')
    meta = json.loads(args.job_json.read_text())
    with h5py.File(meta['artifacts']['combined_h5ad']) as f:
        expected = read_elem(f['obs']).index.astype(str).tolist()
        genes = read_elem(f['var']).index.astype(str)
    gene = 'SFTPC' if 'SFTPC' in genes else str(genes[0])
    app = W.create_app({'JOB_STORAGE': str(args.job_json.parent.parent)})
    if args.mode == 'baseline': W._job_bundle.view = lambda *a, **k: None
    report = {'mode':args.mode, 'job_id':meta['job_id'], 'cells':len(expected), 'requests':[]}
    try:
        client = TestClient(app)
        for iteration in range(2):
            for endpoint in ('umap','expression'):
                params = {} if endpoint == 'umap' else {'gene':gene}
                if args.mode == 'optimized':
                    params['compact'] = 'true'
                    if endpoint == 'expression': params['view'] = 'umap'
                start = time.perf_counter()
                r = client.get(f"/api/jobs/{meta['job_id']}/{endpoint}",params=params)
                elapsed = time.perf_counter()-start
                assert r.status_code == 200, r.text[:1000]
                payload = r.json();points = payload['query' if endpoint=='umap' else 'umap']
                cells = points['columns']['barcode'] if isinstance(points,dict) else [x['barcode'] for x in points]
                assert len(cells) == len(expected) and set(cells) == set(expected)
                record = dict(iteration=iteration, endpoint=endpoint, seconds=elapsed, bytes=len(r.content), all_cell_ids_verified=True)
                report['requests'].append(record);print(record,flush=True)
                del r,payload,points,cells
        import sys
        report['peak_rss_gib'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**3 if sys.platform=='darwin' else 1024**2)
        args.report.write_text(json.dumps(report,indent=2)+'\n')
    finally:
        app.state.job_runner.executor.shutdown(wait=True)


if __name__=='__main__': main()
