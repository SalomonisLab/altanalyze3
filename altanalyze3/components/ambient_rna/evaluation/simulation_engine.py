"""Execution engine for abundance-weighted RNA-release evaluation."""
import gc
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, io

from . import simulation_analysis as analysis

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'benchmarking/ND20_167_release_20261003'
SCENARIOS = [('heterogeneous', rho, seed) for rho in (.20, .30)
             for seed in (101, 102, 103)] + [('depth_proportional', rho, 101) for rho in (.20, .30)]
METHODS = ['contaminated', 'python_auto', 'SoupX_known', 'SoupX_auto']
def generate(base, model, rho, seed, directory):
    cache = directory / 'contaminated.h5ad'
    if cache.exists():
        return ad.read_h5ad(cache)
    p = weights(base.var_names, seed)
    matrices, truths = [], []
    for index, cap in enumerate(('HSC', 'MPP')):
        mask = base.obs.Library.to_numpy() == cap
        clean = base.layers['counts'][mask]
        total = np.asarray(clean.sum(1, dtype=np.float64)).ravel()
        rng = np.random.default_rng(np.random.SeedSequence([92431, seed, index, int(rho * 100)]))
        if model == 'heterogeneous':
            exposure = rng.gamma(shape=4, scale=.25, size=len(total))
            # Independent of cellular RNA content, scaled only to a library total.
            lam = rho / (1-rho) * total.sum() * exposure / exposure.sum()
        else:
            lam = rho / (1-rho) * total
        added = analysis.poisson_matrix(lam, p[cap], rng)
        sparse.save_npz(directory / f'added_{cap}.npz', added)
        matrices.append((clean + added).tocsr())
        pd.DataFrame({'gene': base.var_names, 'probability': p[cap]}).to_csv(
            directory / f'true_profile_{cap}.tsv', sep='\t', index=False)
        realized = np.asarray(added.sum(1, dtype=np.float64)).ravel()
        pd.DataFrame({'cell': base.obs_names[mask], 'expected_added': lam,
                      'realized_added': realized,
                      'realized_rho': realized / (total+realized)}).to_csv(
            directory / f'cell_loading_{cap}.tsv', sep='\t', index=False)
        erng = np.random.default_rng(np.random.SeedSequence([70253, seed, index, int(rho*100)]))
        empty_exposure = erng.gamma(4, 7.5, 50000)
        io.mmwrite(directory / f'empty_{cap}.mtx', analysis.poisson_matrix(empty_exposure, p[cap], erng).T)
        truths.append(dict(capture=cap, baseline_counts=total.sum(), added_counts=realized.sum(),
            observed_counts=total.sum()+realized.sum(), realized_rho=realized.sum()/(total.sum()+realized.sum()),
            median_cell_rho=np.median(realized/(total+realized)),
            cell_rho_p10=np.quantile(realized/(total+realized),.1),
            cell_rho_p90=np.quantile(realized/(total+realized),.9)))
    assert (base.obs.Library.iloc[:len(matrices[0].indptr)-1] == 'HSC').all()
    result = ad.AnnData(sparse.vstack(matrices, format='csr'), obs=base.obs.copy(), var=base.var.copy())
    result.write_h5ad(cache)
    pd.DataFrame(truths).to_csv(directory / 'simulation_truth.tsv', sep='\t', index=False)
    return result


def main(design_override=None):
    OUT.mkdir(exist_ok=True)
    analysis.OUT = OUT
    for name in ('baseline_full.h5ad', 'baseline', 'population_counts.tsv', 'baseline_cells.tsv', 'soupx_reference'):
        assert (OUT/name).exists(), name
    assert design_override is not None
    design = design_override
    dp=OUT/'design.json'
    if dp.exists():
        assert json.loads(dp.read_text()) == json.loads(json.dumps(design))
    else:
        dp.write_text(json.dumps(design,indent=2)+'\n')
    base=ad.read_h5ad(OUT/'baseline_full.h5ad')
    populations=pd.read_csv(OUT/'population_counts.tsv',sep='\t',index_col=0)
    primary=populations.index[populations.min(axis=1)>=100].tolist()
    shared=populations.index[populations.min(axis=1)>=20].tolist()
    selected=os.environ.get('STRUCTURED_CASE')
    for model,rho,seed in SCENARIOS:
        scenario=f'{model}_rho{rho:.2f}_seed{seed}'
        if selected and selected != scenario: continue
        directory=OUT/scenario; directory.mkdir(exist_ok=True)
        records,truths=[],[]
        observed=generate(base,model,rho,seed,directory)
        contaminated,_=analysis.align(observed.copy(),directory/'contaminated')
        for cap in ('HSC','MPP'):
            rd=directory/f'R_{cap}';rd.mkdir(exist_ok=True)
            mask=base.obs.Library.to_numpy()==cap
            if not (rd/'contaminated.mtx').exists():
                io.mmwrite(rd/'contaminated.mtx',observed.X[mask].T)
                os.link(directory/f'empty_{cap}.mtx',rd/'empty.mtx')
                (rd/'genes.txt').write_text('\n'.join(base.var_names)+'\n')
                pd.DataFrame({'cell':base.obs_names[mask],
                    'cluster':contaminated.obs.loc[mask,analysis.KEY].astype(str).to_numpy()}).to_csv(rd/'cells.tsv',sep='\t',index=False)
            if not all((rd/f'{m}_status.tsv').exists() for m in ('SoupX_known','SoupX_auto')):
                with (rd/'run.log').open('a') as log:
                    subprocess.run([analysis.RSCRIPT,str(ROOT/'evaluation/simulation_soupx.R'),str(rd),
                        str(OUT/'soupx_reference'),str(rho),str(seed)],stdout=log,stderr=subprocess.STDOUT,check=True)
        for method in METHODS:
            dest=directory/method;dest.mkdir(exist_ok=True)
            if method=='contaminated':
                corrected=contaminated;seconds=0.
            elif method=='python_auto':
                corrected,seconds=analysis.align(observed.copy(),dest,'auto')
            else:
                statuses=[pd.read_csv(directory/f'R_{cap}'/f'{method}_status.tsv',sep='\t').iloc[0] for cap in ('HSC','MPP')]
                if any(s.status!='success' for s in statuses): continue
                x=sparse.vstack([io.mmread(directory/f'R_{cap}'/f'{method}.mtx').T.tocsr().astype(np.float32)
                    for cap in ('HSC','MPP')],format='csr')
                corrected,_=analysis.align(ad.AnnData(x,obs=base.obs.copy(),var=base.var.copy()),dest)
                seconds=sum(s.correction_seconds+s.profile_seconds for s in statuses)
            signature=analysis.fingerprint(corrected.layers['counts'])
            (dest/'counts.sha256').write_text(signature+'\n')
            sparse.save_npz(dest/'corrected_counts.npz',corrected.layers['counts'])
            metrics=analysis.evaluate(corrected,base,method,dest,shared,primary)
            metrics.update(scenario=scenario,model=model,rho=rho,simulation_seed=seed)
            records.append(metrics)
            truth=analysis.truth_metrics(corrected,base,observed,method,seconds)
            for row in truth: row.update(scenario=scenario,model=model,rho=rho,simulation_seed=seed)
            truths.extend(truth)
            analysis.save_table(records,directory/'method_metrics.tsv')
            analysis.save_table(truths,directory/'truth_recovery.tsv')
            print('DONE',scenario,method,metrics,flush=True)
            if method!='contaminated': del corrected
            gc.collect()
        (directory/'scenario_completed.json').write_text(json.dumps(dict(scenario=scenario,methods=len(records)),indent=2)+'\n')
        del observed,contaminated;gc.collect()
    (OUT/f'provenance_{selected or "all"}.json').write_text(json.dumps(dict(
        source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            (Path(__file__),ROOT/'ambient_subtract.py',ROOT/'evaluation/simulation_soupx.R',ROOT/'evaluation/simulation_analysis.py')},
        baseline_count_hash=analysis.fingerprint(base.layers['counts'])),indent=2)+'\n')


if __name__=='__main__': main()
