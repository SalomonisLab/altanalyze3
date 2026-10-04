"""Population-size and endogenous-abundance weighted RNA release benchmark.

Uses corrected baseline expression as a cell-lysis model, never measured original
ambient RNA. This revision follows the user's explicit release-composition rule.
"""
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from . import simulation_engine as engine

OUT=Path(__file__).resolve().parent.parent/'benchmarking/ND20_167_release_20261003'
BASELINE_ROOT=OUT
SCENARIOS=engine.SCENARIOS
METHODS=engine.METHODS


def release_profiles(base):
    profiles,contributors={},[]
    for cap in ('HSC','MPP'):
        capmask=base.obs.Library.to_numpy()==cap
        burden=np.zeros(base.n_vars,dtype=np.float64)
        for state in sorted(base.obs.loc[capmask,'fixed_population'].unique()):
            mask=capmask&(base.obs.fixed_population.to_numpy()==state)
            n=int(mask.sum())
            # N_k * mean raw baseline count per cell exactly equals summed RNA.
            mean=np.asarray(base.layers['counts'][mask].sum(0,dtype=np.float64)).ravel()/n
            contribution=n*mean;burden+=contribution
            contributors.append(dict(capture=cap,population=state,cells=n,
                mean_RNA_per_cell=float(mean.sum()),released_RNA_weight=float(contribution.sum())))
        pooled=np.asarray(base.layers['counts'][capmask].sum(0,dtype=np.float64)).ravel()
        assert np.allclose(burden,pooled,rtol=1e-12,atol=1e-8)
        profiles[cap]=burden/burden.sum()
    return profiles,pd.DataFrame(contributors)


def main():
    OUT.mkdir(exist_ok=True)
    base=ad.read_h5ad(OUT/'baseline_full.h5ad')
    profiles,contributions=release_profiles(base)
    contributions.to_csv(OUT/'release_contributors.tsv',sep='\t',index=False)
    pd.DataFrame({'gene':base.var_names,**profiles}).to_csv(OUT/'release_profiles.tsv',sep='\t',index=False)
    profile_hash=hashlib.sha256(np.column_stack([profiles[c] for c in ('HSC','MPP')]).tobytes()).hexdigest()
    del base
    engine.OUT=OUT
    engine.weights=lambda genes,seed: {cap:p.copy() for cap,p in profiles.items()}
    design=dict(baseline=str(OUT/'baseline_full.h5ad'),
        profile_rule='p_cg proportional to sum_k N_ck * mean_corrected_baseline_count_ckg; equivalent to pooled corrected cell counts',
        release_assumption='Equal lysis propensity per cell; per-cell RNA mass preserved. All filtered baseline cells contribute.',
        models={'heterogeneous':'A_jg~Poisson(lambda_j*p_cg), lambda_j proportional to independent Gamma(4,0.25) exposure; expected library contamination .20 or .30',
                'depth_proportional':'A_jg~Poisson(rho/(1-rho)*baseline_total_j*p_cg)'},
        scenarios=SCENARIOS,methods=METHODS,umap_seeds=[0,1,2],
        empty_rule='50000 droplets/library, independent Gamma(4,7.5) exposure and Poisson(exposure*p)',
        profile_seed_rule='Profiles fixed by baseline cell abundance; seeds vary independent exposure and Poisson count sampling only',
        original_uncorrected_input_used=False,original_ambient_profile_used=False,
        corrected_baseline_expression_used_for_release=True,reference_expression_used=False,
        background_profile_hash=profile_hash,
        known_fraction_arm='Nominal library fraction, not cell-specific truth or true profile',
        design_scope='User-specified abundance/population weighting; fixed before release-model outcomes')
    engine.main(design_override=design)


if __name__=='__main__': main()
