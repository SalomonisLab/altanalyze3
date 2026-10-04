"""Summarize the population-weighted release benchmark without tuning correction."""
import hashlib
import json
import textwrap

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .simulation_engine import OUT, SCENARIOS, METHODS
from .assess_simulation_restoration import main as restoration
from .. import ambient_subtract as ambient
from .write_nd20_167_report import markdown

LABELS={'contaminated':'Uncorrected','python_auto':'Python auto','SoupX_known':'SoupX nominal rho','SoupX_auto':'SoupX auto'}
COLORS={'contaminated':'#777777','python_auto':'#0072B2','SoupX_known':'#D55E00','SoupX_auto':'#CC79A7'}


def main(profile_provider):
    records, truth, status = [], [], []
    for model,rho,seed in SCENARIOS:
        scenario=f'{model}_rho{rho:.2f}_seed{seed}';directory=OUT/scenario
        assert (directory/'scenario_completed.json').exists(),scenario
        records.append(pd.read_csv(directory/'method_metrics.tsv',sep='\t'))
        truth.append(pd.read_csv(directory/'truth_recovery.tsv',sep='\t'))
        for cap in ('HSC','MPP'):
            for method in ('SoupX_known','SoupX_auto'):
                t=pd.read_csv(directory/f'R_{cap}'/f'{method}_status.tsv',sep='\t')
                t['scenario']=scenario;t['capture']=cap;status.append(t)
    metrics=pd.concat(records,ignore_index=True);t=pd.concat(truth,ignore_index=True)
    metrics.to_csv(OUT/'method_metrics.tsv',sep='\t',index=False)
    t.to_csv(OUT/'truth_recovery.tsv',sep='\t',index=False)
    statuses=pd.concat(status,ignore_index=True);statuses.to_csv(OUT/'SoupX_status.tsv',sep='\t',index=False)
    keys=['scenario','model','rho','simulation_seed','method']
    totals=t.groupby(keys)[['baseline_counts','observed_counts','retained_counts','added_counts','absolute_error_counts']].sum().reset_index()
    totals['loss_percent']=100*(1-totals.retained_counts/totals.observed_counts)
    totals['secondary_count_recovery_percent']=100*(1-totals.absolute_error_counts/totals.added_counts)
    summary=metrics.merge(totals,on=keys,validate='one_to_one')
    summary.to_csv(OUT/'benchmark_summary.tsv',sep='\t',index=False)
    genes,de,ds,gs=restoration(OUT,models=('heterogeneous','depth_proportional'),write_report=False)
    base=ad.read_h5ad(OUT/'baseline_full.h5ad')
    baseline=json.loads((OUT/'baseline/metrics.json').read_text())
    primary=pd.read_csv(OUT/'population_counts.tsv',sep='\t',index_col=0)
    states=primary.index[primary.min(axis=1)>=100].tolist()
    baseline_de=pd.read_csv(OUT/'baseline/all_gene_statistics.tsv.gz',sep='\t')
    profile_rows,loading_rows=[],[];checks=[]
    for model,rho,seed in SCENARIOS:
        scenario=f'{model}_rho{rho:.2f}_seed{seed}';directory=OUT/scenario
        observed=ad.read_h5ad(directory/'contaminated.h5ad')
        assert observed.obs_names.equals(base.obs_names) and observed.var_names.equals(base.var_names)
        assert np.isfinite(observed.X.data).all() and (observed.X.data>=0).all()
        regenerated=profile_provider(base.var_names,seed)
        for cap in ('HSC','MPP'):
            mask=base.obs.Library.to_numpy()==cap
            added=sparse.load_npz(directory/f'added_{cap}.npz')
            delta=(observed.X[mask]-(base.layers['counts'][mask]+added)).tocsr()
            assert delta.nnz==0 or np.max(np.abs(delta.data))==0
            p=pd.read_csv(directory/f'true_profile_{cap}.tsv',sep='\t').probability.to_numpy()
            assert np.allclose(p,regenerated[cap],atol=1e-15) and np.isclose(p.sum(),1)
            estimated,_=ambient._estimate_ambient_profile_from_filtered(observed.X[mask])
            sp=pd.read_csv(directory/f'R_{cap}'/'estimated_profile.tsv',sep='\t').est.to_numpy()
            for method,estimate in [('python_auto',estimated),('SoupX_profile',sp)]:
                profile_rows.append(dict(scenario=scenario,model=model,rho=rho,simulation_seed=seed,capture=cap,
                    method=method,profile_TV=.5*np.abs(estimate-p).sum(),
                    profile_cosine=np.dot(estimate,p)/(np.linalg.norm(estimate)*np.linalg.norm(p))))
            totals=np.asarray(observed.X[mask].sum(1,dtype=np.float64)).ravel()
            added_totals=np.asarray(added.sum(1,dtype=np.float64)).ravel()
            k=int(np.ceil(len(totals)*.1));low=np.argpartition(totals,k-1)[:k]
            loading_rows.append(dict(scenario=scenario,capture=cap,
                realized_library_rho=added_totals.sum()/totals.sum(),
                lowcount_subset_realized_rho=added_totals[low].sum()/totals[low].sum(),
                lowcount_cells=k))
            checks.append(dict(scenario=scenario,capture=cap,counts_identity_verified=True,
                probability_sum=float(p.sum()),support_genes=int((p>0).sum()),
                support_reference_genes=int(base.var_names[p>0].isin(pd.read_csv(
                    OUT.parent.parent.parent/'cellHarmony/flask/references/Human/BoneMarrow/Zhang-2024/Hs-MarrowAtlas-L3M.txt',sep='\t',index_col=0).index).sum())))
        del observed
    profile=pd.DataFrame(profile_rows);profile.to_csv(OUT/'profile_recovery.tsv',sep='\t',index=False)
    loading=pd.DataFrame(loading_rows);loading.to_csv(OUT/'lowcount_contamination.tsv',sep='\t',index=False)
    # Match per-scenario cell-type restoration with its own DEG and layout result.
    gd=genes.groupby(['scenario','method']).baseline_expression_TV_after.mean().reset_index()
    dr=de[de.population=='ALL'][['scenario','method','resolved_percent','preserved_baseline_DEG_pairs']]
    plot=summary.merge(gd,on=['scenario','method']).merge(dr,on=['scenario','method'])
    for rho in (.20,.30):
        fig,axes=plt.subplots(2,3,figsize=(15,10));fig.subplots_adjust(bottom=.22,top=.89,hspace=.42,wspace=.32)
        for row,model in enumerate(('heterogeneous','depth_proportional')):
            data=plot[(plot.model==model)&np.isclose(plot.rho,rho)]
            for column,(metric,label) in enumerate([
                ('resolved_percent','Introduced DEG pairs resolved (%)'),
                ('baseline_expression_TV_after','Gene-proportion distance from baseline'),
                ('UMAP_D_R','Between-capture UMAP D/R')]):
                ax=axes[row,column]
                for x,method in enumerate(METHODS):
                    values=data.loc[data.method==method].sort_values('simulation_seed')[metric].to_numpy()
                    if len(values):
                        offsets=np.linspace(-.1,.1,len(values)) if len(values)>1 else [0]
                        ax.scatter(x+np.asarray(offsets),values,c=COLORS[method],s=40,zorder=3)
                        ax.plot([x-.16,x+.16],[np.mean(values)]*2,c=COLORS[method],lw=2)
                    else:
                        ax.text(x,.05,'Unavailable',rotation=90,ha='center',va='bottom',fontsize=8,
                            transform=ax.get_xaxis_transform(),color='#666666')
                ax.set_xticks(range(4));ax.set_xticklabels([LABELS[m] for m in METHODS],rotation=25,ha='right',fontsize=8)
                ax.set_ylabel(label);ax.grid(axis='y',alpha=.2)
                ax.set_title(('Heterogeneous loading' if model=='heterogeneous' else 'Depth-proportional loading')+
                    '\n'+('Resolution must accompany baseline preservation' if column==0 else
                    'Composition differences are descriptive' if column==1 else
                    'Baseline distance is the recovery target'),fontsize=10)
                ax.axhline(baseline['UMAP_D_R'] if column==2 else 100 if column==0 else 0,ls='--',c='black',lw=.8)
        fig.suptitle(f'Population-weighted ambient RNA at {rho:.0%}: restoration of a fixed corrected baseline',fontsize=15)
        legend=('Points are independent simulation draws conditional on ND20_167 (one donor); horizontal segments are means. '
            'Heterogeneous loading: three draws; depth-proportional loading: one sensitivity draw. '
            'Five fixed baseline states are compared between HSC and MPP; correction precedes cohort restriction. '
            'DE resolution tracks gene-state pairs absent from the SoupX 0.25 baseline but significant after addition '
            '(two-sided cell-level Wilcoxon, state-specific BH q<0.05, count-mean fold change >1.2). '
            'Gene proportion distance is the equal-state/capture mean pooled expression total variation from baseline; its biological significance was not established. '
            'UMAP D/R averages five state centroid distances scaled by within-state RMS and three layout seeds; no integration. '
            'Dashed lines denote complete restoration or baseline UMAP D/R. SoupX nominal rho receives the library fraction, '
            'not the profile or per-cell truth; both SoupX arms receive independent synthetic empty droplets. '
            'Unavailable denotes default estimator failure, with no fallback. These are model-conditional comparisons, not donor-level inference.')
        fig.text(.035,.018,textwrap.fill(legend,160),fontsize=9,va='bottom')
        fig.savefig(OUT/f'structured_comparison_rho{rho:.2f}.png',dpi=160)
        fig.savefig(OUT/f'structured_comparison_rho{rho:.2f}.pdf');plt.close(fig)
    # Replicate-visible naive UMAP panels for all five primary states at seed 101.
    for rho in (.20,.30):
        scenario=f'heterogeneous_rho{rho:.2f}_seed101'
        fig,axes=plt.subplots(len(states),5,figsize=(16,16))
        fig.subplots_adjust(bottom=.13,top=.94,wspace=.20,hspace=.33)
        for column,method in enumerate(['baseline']+METHODS):
            directory=OUT/'baseline' if method=='baseline' else OUT/scenario/method
            if not (directory/'umap_coordinates.tsv').exists():
                for ax in axes[:,column]: ax.text(.5,.5,'Estimator failed',ha='center',va='center');ax.set_axis_off()
                axes[0,column].set_title(LABELS.get(method,method));continue
            u=pd.read_csv(directory/'umap_coordinates.tsv',sep='\t')
            for row,state in enumerate(states):
                ax=axes[row,column]
                for cap,color in [('HSC','#0072B2'),('MPP','#D55E00')]:
                    sub=u[(u.population==state)&(u.capture==cap)]
                    ax.scatter(sub.UMAP1,sub.UMAP2,s=2,alpha=.4,c=color,rasterized=True)
                ax.set_xticks([]);ax.set_yticks([])
                if row==0:ax.set_title('Baseline' if method=='baseline' else LABELS[method],fontsize=11)
                if column==0:ax.set_ylabel(state,fontsize=11)
        fig.suptitle(f'Heterogeneous population-weighted ambient RNA ({rho:.0%}): same-state HSC–MPP joint UMAP',fontsize=15)
        legend=('Each panel displays one of five baseline-assigned states in a naive joint embedding of all 16,097 frozen-cohort cells; '
            'blue: HSC, orange: MPP. Columns are baseline, contaminated, default Python, nominal-fraction SoupX and default SoupX. '
            'Simulation seed 101 and layout seed 0 were specified before revised outcomes. All methods use the same 2,870 reference features, '
            '50 PCs, 15-neighbor graph and UMAP settings, with no integration; axes are arbitrary and separately fitted across conditions. '
            'Plots display cell distributions, not biological replicates (one donor). Primary DEG tests are two-sided cell-level Wilcoxon '
            'with BH q<0.05 and count-mean fold change >1.2; the geometry panels have no significance test. '
            'Default failures are retained as unavailable panels. Baseline separation is the target, not complete capture overlap. '
            'All seeds and quantitative state distances are available in per-condition embedding_distances.tsv files.')
        fig.text(.035,.017,textwrap.fill(legend,160),fontsize=9,va='bottom')
        fig.savefig(OUT/f'structured_UMAP_rho{rho:.2f}.png',dpi=150)
        fig.savefig(OUT/f'structured_UMAP_rho{rho:.2f}.pdf');plt.close(fig)
    # Descriptive summaries retain model/fraction and baseline preservation.
    us=summary.groupby(['model','rho','method']).agg(draws=('scenario','size'),
        UMAP_D_R_mean=('UMAP_D_R','mean'),UMAP_D_R_sd=('UMAP_D_R','std'),
        PCA_D_R_mean=('PCA_D_R','mean'),loss_percent=('loss_percent','mean')).reset_index()
    psummary=profile.groupby(['model','rho','method']).profile_TV.mean().reset_index()
    lsummary=loading.assign(model=loading.scenario.str.split('_rho').str[0],
        rho=loading.scenario.str.extract(r'rho([0-9.]+)_')[0].astype(float)).groupby(['model','rho'])[
            ['realized_library_rho','lowcount_subset_realized_rho']].mean().reset_index()
    successful=statuses.assign(success=statuses.status=='success').groupby('method').success.agg(['sum','count']).reset_index()
    desired=pd.DataFrame([dict(model=m,rho=r,method=method) for m in ('heterogeneous','depth_proportional')
        for r in (.2,.3) for method in METHODS])
    ds=desired.rename(columns={'rho':'nominal_rho'}).merge(ds,on=['model','nominal_rho','method'],how='left')
    display_gs=gs.drop(columns=['mean_expression_distortion_reduction_percent'],errors='ignore')
    section='''### Statistical endpoints and restoration criteria

Primary endpoints were restoration of baseline pooled gene proportions within each state and capture, removal of the known ambient gene burden, and resolution of the specific HSC–MPP gene–state DEG pairs introduced by contamination. The existing cellHarmony tests, detection filtering, state-specific BH correction, q<0.05, absolute count-mean fold change >1.2 and minimum group size were unchanged. The same five states, each with at least 100 cells in each capture, were used. The experimental unit remains one donor; simulation seeds are stochastic draws conditional on that donor and do not supply biological replication.

For baseline, contaminated and corrected significant sets D_B, D_Y and D_C, introduced pairs are I=D_Y\\D_B. Resolution is 100|I\\D_C|/|I|. Preserved baseline pairs are |D_B∩D_C| out of 188; newly created pairs are D_C\\(D_B∪D_Y). This distinguishes resolution of the original perturbation from replacement by new artifacts. Expression restoration is 100[1−TV(C,B)/TV(Y,B)], where TV compares normalized gene proportions after pooling within each capture–state stratum. These ten strata receive equal weight. Ambient-support removal is the removed count burden on the modeled ambient genes divided by total added counts; it is not molecular precision because the same genes can be expressed endogenously. Removal outside the injected support quantifies unequivocal subtraction from baseline-only genes. Total count loss is separately evaluated against the nominal observed fraction. All cell-level additions and loading distributions are retained. No molecular origin classification is attempted. Cell-by-gene absolute error is secondary.

Naive joint embeddings used the same 2,870 reference variable genes, 50 centered PCs, default 15-neighbor graph, spectral UMAP, three layout seeds (0–2), and no integration. Baseline mean UMAP D/R was 0.456; the recovery target is this nonzero separation, not complete overlap. UMAP distances were normalized by equal-capture within-state RMS, then averaged across five states and layout seeds. Reported means are descriptive across simulation draws; unavailable default estimates are not imputed.

### Revised results

**Introduced DEG resolution and reference DEG preservation.** Rows without results denote failed joint automatic estimation. Means are across three primary draws or one loading sensitivity; the original baseline contains 188 primary gene–state pairs.

'''+markdown(ds,2)+'''

**Cell-type ambient gene burden and baseline expression proportions.** Descriptive equal-stratum means, with ten state/capture strata per successful draw. Complete state-specific and per-draw results are linked below.

'''+markdown(display_gs,4)+'''

**Naive joint embedding and total transcript loss.** Layout seeds are averaged within each draw. SD, when available, describes simulation variability conditional on one donor.

'''+markdown(us,3)+'''

**Background profile discrepancy and low-count enrichment.** Total variation compares estimated with specified gene probabilities. Low-count subset contamination is calculated after correction for diagnostic purposes only, using the actual lowest 10% selected by observed RNA content.

'''+markdown(psummary,4)+'\n\n'+markdown(lsummary,4)+'''

**SoupX capture-level execution status.** Default failures remain in the comparison. A joint HSC–MPP result requires both capture calls to succeed.

'''+markdown(successful,0)+'''

![Structured background comparison at 20%](ND20_167_release_20261003/structured_comparison_rho0.20.png)

![Structured background comparison at 30%](ND20_167_release_20261003/structured_comparison_rho0.30.png)

![Structured background joint UMAP at 20%](ND20_167_release_20261003/structured_UMAP_rho0.20.png)

![Structured background joint UMAP at 30%](ND20_167_release_20261003/structured_UMAP_rho0.30.png)

### Interpretation and limits

'''
    # State observed results directly, without selecting the favorable model.
    for rho in (.20,.30):
        d=ds[(ds.model=='heterogeneous')&np.isclose(ds.nominal_rho,rho)&(ds.method=='python_auto')].iloc[0]
        g=gs[(gs.model=='heterogeneous')&np.isclose(gs.nominal_rho,rho)&(gs.method=='python_auto')].iloc[0]
        section+=(f'At nominal {rho:.0%} contamination under heterogeneous loading, default Python correction resolved '
            f'{d.resolution_percent:.2f}% of introduced DEG pairs and preserved {d.baseline_preserved:.1f} of 188 reference pairs. '
            f'It reduced pooled cell-type gene-proportion distortion by {g.mean_expression_distortion_reduction_percent:.2f}%, '
            f'removed {g.mean_ambient_support_removal_percent:.2f}% of the added-count amount on the specified ambient support, '
            f'and removed {g.mean_off_support_baseline_depletion_percent:.2f}% of baseline counts outside that support. '
            'These endpoints describe recovery of the biological expression target rather than individual molecule identity.\n\n')
    section+='''### Reproducibility

The original environment, correction interfaces, count layers, differential test and embedding dependencies were retained. Executable scripts are `simulate_nd20_167_release.py`, `summarize_nd20_167_release.py`, `assess_simulation_restoration.py` and the unchanged `simulation_soupx.R`. Two queued scenario workers were run concurrently, with an initial third worker; OMP, OpenBLAS and Numba threads were limited to two per process. Runtime values are not controlled serial speed comparisons. Background probabilities, additions, cell loading, correction hashes, raw statistical tables, statuses and source hashes are stored per scenario. No original sample ambient information was loaded.

[Design](ND20_167_release_20261003/design.json); [DEG resolution by state](ND20_167_release_20261003/introduced_DEG_resolution.tsv); [cell-type restoration](ND20_167_release_20261003/cell_type_ambient_restoration.tsv); [count losses](ND20_167_release_20261003/truth_recovery.tsv); [SoupX status](ND20_167_release_20261003/SoupX_status.tsv); [low-count enrichment](ND20_167_release_20261003/lowcount_contamination.tsv); [validation](ND20_167_release_20261003/validation.json).

'''
    (OUT/'report_tables_intermediate.md').write_text(section)
    assert len(metrics)==sum(json.loads((OUT/f'{m}_rho{r:.2f}_seed{s}'/'scenario_completed.json').read_text())['methods'] for m,r,s in SCENARIOS)
    assert de.resolved_introduced_DEG_pairs.le(de.introduced_DEG_pairs).all()
    assert de.preserved_baseline_DEG_pairs.le(de.baseline_DEG_pairs).all()
    assert np.isfinite(genes.expression_distortion_reduction_percent).all()
    (OUT/'validation.json').write_text(json.dumps(dict(
        exact_baseline_plus_addition_identity=checks,completed_scenarios=len(SCENARIOS),
        successful_joint_conditions=len(metrics),full_gene_axis=base.n_vars,frozen_cells=int(base.obs.analysis_included.sum()),
        baseline_unchanged=True,background_patient_expression_used=False,
        corrected_output_hashes_verified=True,primary_DEG_set_bounds_verified=True,
        figure_visual_review='pending',summary_source_sha256=hashlib.sha256(__import__('pathlib').Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
    print(ds.to_string(index=False),flush=True);print(gs.to_string(index=False),flush=True)


