"""Exercise a real all-modality upload, cell-level differentials and serving routes.

Replicas of real uploads test computational scaling; they are not independent
biological samples. Results are for performance validation only. Use a scratch
root: the source job is read-only and all outputs are written under that root.
"""
import argparse
import json
import os
from pathlib import Path
import threading
import time

import psutil


def monitor(root, limit_gib):
    process = psutil.Process()
    state = dict(peak_rss_gib=0.0, idle_rss_gib=process.memory_info().rss / 1024**3)
    stop = threading.Event()
    def watch():
        while not stop.wait(0.2):
            children = process.children(recursive=True)
            total = 0
            for proc in [process, *children]:
                try:
                    total += proc.memory_info().rss
                except psutil.NoSuchProcess:
                    pass
            state['peak_rss_gib'] = max(state['peak_rss_gib'], total / 1024**3)
            if limit_gib > 0 and total > limit_gib * 1024**3:
                for proc in reversed(children):
                    try:
                        proc.kill()
                    except psutil.NoSuchProcess:
                        pass
                (root / 'memory_guard.json').write_text(json.dumps(state))
                os._exit(70)
    threading.Thread(target=watch, daemon=True).start()
    return state, stop


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['prepare', 'prepare-saved', 'pipeline', 'differential', 'serve', 'serve-concurrent', 'differential-views'])
    parser.add_argument('root', type=Path)
    parser.add_argument('--source-job', type=Path)
    parser.add_argument('--registry', type=Path, default=Path(__file__).resolve().parents[2] / 'flask/reference_config.json')
    parser.add_argument('--replicas', type=int, default=3)
    parser.add_argument('--rss-limit-gib', type=float, default=28)
    parser.add_argument('--modality', default='rna')
    parser.add_argument('--rounds', type=int, default=3)
    args = parser.parse_args()
    args.root.mkdir(parents=True, exist_ok=True)
    from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
    from altanalyze3.components.cellHarmony.flask.tasks import JobRunner
    store = JobStore(args.root / 'jobs')
    fixture = args.root / 'fixture.json'
    if args.mode == 'prepare-saved':
        if fixture.exists() or args.source_job is None or args.replicas < 1:
            parser.error('prepare-saved needs --source-job, positive replicas and a fresh root')
        source = json.loads((args.source_job / 'job.json').read_text())
        ids = []
        for _ in range(args.replicas):
            job = store.create_job(source['species'], source['reference'], source.get('ambient_option'), [])['job_id']
            snapshot = dict(source)
            for key in ('job_id', 'created_at', 'updated_at'):
                snapshot.pop(key, None)
            store.update_job(job, **snapshot)
            ids.append(job)
        fixture.write_text(json.dumps(dict(job_id=ids[0], job_ids=ids,
                                           source=str(args.source_job.resolve()), replicas=args.replicas)))
        print(fixture.read_text(), flush=True)
        return
    if args.mode == 'prepare':
        if fixture.exists() or args.source_job is None or args.replicas < 1:
            parser.error('prepare needs --source-job, positive replicas and a fresh root')
        source = json.loads((args.source_job / 'job.json').read_text())
        files = []
        for replica in range(args.replicas):
            for entry in source['files']:
                files.append(dict(entry, filename=f"replica{replica}_{entry['filename']}",
                                  sample_name=f"replica{replica}_{entry['sample_name']}"))
        meta = store.create_job(source['species'], source['reference'], source.get('ambient_option'), files)
        for replica in range(args.replicas):
            for entry in source['files']:
                src = args.source_job / 'uploads' / entry['filename']
                (store.uploads_dir(meta['job_id']) / f"replica{replica}_{entry['filename']}").symlink_to(src.resolve())
        store.update_job(meta['job_id'], qc=dict(source['qc'], impute_modalities=['all']))
        fixture.write_text(json.dumps(dict(job_id=meta['job_id'], source=str(args.source_job.resolve()),
                                           replicas=args.replicas, input_files=len(files))))
        print(fixture.read_text(), flush=True)
        return
    job_id = json.loads(fixture.read_text())['job_id']
    state, stop = monitor(args.root, args.rss_limit_gib)
    started = time.perf_counter()
    runner = JobRunner(store, args.registry, isolate_jobs=True, export_approx_pdfs=False)
    try:
        if args.mode == 'pipeline':
            runner.submit_pipeline(job_id)
            runner._futures[f'pipeline:{job_id}'].result()
            meta = store.get_job(job_id)
            assert meta['status'] == 'completed', meta.get('message')
            state.update(status=meta['status'], modalities=list(meta['modality_artifacts']),
                         bundle_status=meta.get('bundle', {}).get('status'))
        elif args.mode == 'differential':
            meta = store.get_job(job_id)
            names = [entry['sample_name'] for entry in meta['files']]
            cut = max(1, len(names) // 2)
            cfg = dict(modality=args.modality, population_col=meta['cluster_key'], sample_field='Library',
                       group1_samples=names[:cut], group2_samples=names[cut:], comparison_type='cells')
            store.update_job(job_id, differential=dict(status='queued', config=cfg, artifacts={}))
            runner.submit_differential(job_id)
            runner._futures[f'differential:{job_id}'].result()
            diff = store.get_job(job_id)['differential']
            assert diff['status'] == 'completed', diff.get('message')
            state.update(modality=args.modality, status=diff['status'], config=cfg)
        elif args.mode == 'serve':
            from fastapi.testclient import TestClient
            from altanalyze3.components.cellHarmony.webapp.app import create_app, _get_expression_cache
            app = create_app(dict(JOB_STORAGE=str(store.root), REFERENCE_REGISTRY=str(args.registry)))
            client = TestClient(app)
            meta = store.get_job(job_id)
            outcomes = []
            for cycle in range(args.rounds):
                ids = json.loads(fixture.read_text()).get('job_ids', [job_id])
                job_id = ids[cycle % len(ids)]
                meta = store.get_job(job_id)
                for modality in meta['modality_artifacts']:
                    cache = _get_expression_cache(app, meta, modality)
                    genes = list(map(str, cache['var_names'][:3]))
                    del cache
                    for route, params in [('umap', {}), ('expression', {'gene': genes[0]}),
                                          ('dotplot', {'genes': ','.join(genes)}),
                                          ('combplot', {'genes': ','.join(genes)}), ('genes', {})]:
                        response = client.get(f'/api/jobs/{job_id}/{route}', params=dict(params, modality=modality))
                        assert response.status_code == 200, (route, modality, response.text[:500])
                        outcomes.append([cycle, modality, route, response.status_code])
                        del response
                for question in ['Which cell states are most abundant?', 'What is the best modality marker of HSC-1?',
                                 'Where is CD34 most expressed?']:
                    response = client.post(f'/api/jobs/{job_id}/chat', json={'question': question})
                    assert response.status_code == 200, response.text[:500]
                    outcomes.append([cycle, 'chat', question, response.json().get('status')])
                    del response
                print(json.dumps(dict(cycle=cycle, job_id=job_id,
                                      rss_gib=psutil.Process().memory_info().rss / 1024**3,
                                      retained=app.state.result_cache_budget.snapshot())), flush=True)
            state.update(requests=outcomes, retained=app.state.result_cache_budget.snapshot())
            app.state.job_runner.executor.shutdown(wait=True)
            client.close()
        elif args.mode == 'serve-concurrent':
            from concurrent.futures import ThreadPoolExecutor
            from fastapi.testclient import TestClient
            from altanalyze3.components.cellHarmony.webapp.app import create_app, _get_expression_cache
            app = create_app(dict(JOB_STORAGE=str(store.root), REFERENCE_REGISTRY=str(args.registry)))
            meta = store.get_job(job_id)
            features = {m: str(_get_expression_cache(app, meta, m)['var_names'][0])
                        for m in meta['modality_artifacts']}
            def visit(visitor):
                outcomes = []
                with TestClient(app) as client:
                    for cycle in range(args.rounds):
                        for modality, gene in features.items():
                            for route, params in [('umap', {}), ('expression', {'gene':gene}),
                                                  ('dotplot', {'genes':gene}), ('combplot', {'genes':gene})]:
                                response = client.get(f'/api/jobs/{job_id}/{route}', params=dict(params, modality=modality))
                                assert response.status_code == 200, response.text[:500]
                                outcomes.append([visitor, cycle, modality, route, 200])
                                del response
                        response = client.post(f'/api/jobs/{job_id}/chat',
                                               json={'question':'What is the best modality marker of HSC-1?'})
                        assert response.status_code == 200 and response.json().get('status') == 'ok'
                        outcomes.append([visitor, cycle, 'chat', 200])
                return outcomes
            try:
                with ThreadPoolExecutor(max_workers=2) as visitors:
                    results = list(visitors.map(visit, range(2)))
                state.update(requests=[r for rows in results for r in rows],
                             retained=app.state.result_cache_budget.snapshot())
            finally:
                app.state.job_runner.executor.shutdown(wait=True)
        else:
            from fastapi.testclient import TestClient
            from altanalyze3.components.cellHarmony.webapp.app import create_app, _get_differential_detail_table
            app = create_app(dict(JOB_STORAGE=str(store.root), REFERENCE_REGISTRY=str(args.registry)))
            outcomes = []
            try:
                with TestClient(app) as client:
                    history = store.get_job(job_id).get('differential_history', {})
                    assert len(history) >= 6, 'Run all modality differentials first.'
                    for run_id, run in history.items():
                        response = client.post(f'/api/jobs/{job_id}/differential/select', params={'contrast': run_id})
                        assert response.status_code == 200, response.text[:500]
                        payload = response.json()
                        for view, populations in payload.get('visualization_populations', {}).items():
                            if not populations or view not in {'summary', 'heatmap', 'volcano', 'network', 'go', 'table'}:
                                continue
                            response = client.get(f'/api/jobs/{job_id}/differential/interactive/{view}',
                                                  params={'population': populations[0]})
                            assert response.status_code == 200, (view, run['config']['modality'], response.text[:500])
                            outcomes.append([run['config']['modality'], view, 200])
                        meta = store.get_job(job_id)
                        has_detail = any(k.startswith('DEG_detailed_') for k in meta['differential'].get('artifacts', {}))
                        frame = _get_differential_detail_table(app, meta) if has_detail else None
                        if frame is not None and not frame.empty:
                            row = frame.iloc[0]
                            response = client.get(f'/api/jobs/{job_id}/differential/interactive/gene',
                                                  params={'population': str(row.population), 'gene': str(row.gene)})
                            assert response.status_code == 200, response.text[:500]
                            outcomes.append([run['config']['modality'], 'gene', 200])
                        response = client.get(f'/api/jobs/{job_id}/integrated/cross-pathways',
                                              params={'cell_state':'HSC-1', 'source':'differential', 'contrast':run_id})
                        assert response.status_code == 200, response.text[:500]
                        outcomes.append([run['config']['modality'], 'cross-pathways', 200])
                state.update(requests=outcomes, retained=app.state.result_cache_budget.snapshot())
            finally:
                app.state.job_runner.executor.shutdown(wait=True)
    finally:
        runner.executor.shutdown(wait=True)
        stop.set()
        state.update(elapsed_sec=time.perf_counter()-started,
                     final_rss_gib=psutil.Process().memory_info().rss / 1024**3)
        path = args.root / f'{args.mode}_{args.modality}.json'
        path.write_text(json.dumps(state, indent=2))
        print(json.dumps(state), flush=True)


if __name__ == '__main__':
    main()
