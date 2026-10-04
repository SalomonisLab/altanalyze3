"""Export benchmark tables, interpretable figures, and Markdown after all methods finish."""
from pathlib import Path
import hashlib
import json
import platform
import subprocess
import textwrap

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .benchmark_nd20_167 import OUT, ROOT, REF, METHODS, BASE
from altanalyze3.components.cellHarmony import cellHarmony_differential as diff

def figure_save(fig, name, legend):
    fig.text(.04,.015,textwrap.fill(legend,155),ha='left',va='bottom',fontsize=9)
    fig.savefig(OUT/(name+'.png'),dpi=180)
    fig.savefig(OUT/(name+'.pdf'))
    plt.close(fig)

def main():
    metrics = pd.read_csv(OUT/'method_metrics.tsv',sep='\t').set_index('method').loc[METHODS]
    loss = pd.read_csv(OUT/'transcript_loss.tsv',sep='\t')
    distances = pd.read_csv(OUT/'embedding_distances.tsv',sep='\t')
    populations = pd.read_csv(OUT/'population_counts.tsv',sep='\t',index_col=0)
    predominant = populations.index[populations[['HSC','MPP']].min(axis=1)>=100].tolist()
    shared = populations.index[populations[['HSC','MPP']].min(axis=1)>=20].tolist()
    d = distances[distances.population.isin(predominant)]
    med = d.groupby('method').agg(UMAP_distance=('UMAP_distance','mean'),
          UMAP_within_scaled=('UMAP_distance_within_scaled','mean'),
          UMAP_global_scaled=('UMAP_distance_global_scaled','mean'),
          PCA_within_scaled=('PCA_distance_within_scaled','mean'),
          PCA_mixing=('PCA_neighbor_mixing','mean')).reindex(METHODS)
    summary = metrics.join(med)
    summary['DEG_reduction_percent'] = 100*(1-summary.predominant_DEG_gene_state_pairs/summary.loc['uncorrected','predominant_DEG_gene_state_pairs'])
    summary['UMAP_reduction_percent'] = 100*(1-summary.UMAP_within_scaled/summary.loc['uncorrected','UMAP_within_scaled'])
    summary.to_csv(OUT/'benchmark_summary.tsv',sep='\t')
    de_rows = []; secondary = []; common_family = []; baseline_eligible = {}
    for method in METHODS:
        table = pd.read_csv(OUT/method/'DE_summary.tsv',sep='\t').set_index('population')
        complete = pd.read_csv(OUT/method/'all_gene_statistics.tsv.gz',sep='\t')
        for pop in shared:
            s = complete[complete.population==pop]
            assert len(s)>0
            n = int(((s.fdr<.05)&(s.log2fc.abs()>np.log2(1.2))).sum())
            if pop in table.index: assert n == int(table.loc[pop,'num_DEG'])
            else: assert n == 0
            de_rows.append(dict(method=method,population=pop,DEGs=n))
        # Secondary effect-size check: same Wilcoxon p/FDR, mean CP10k rather than UMI means.
        a = ad.read_h5ad(OUT/method/'aligned_fixed_cohort.h5ad')
        stats = complete
        fixed = pd.read_csv(OUT/'fixed_population_assignments.tsv',sep='\t',index_col=0).iloc[:,0]
        a.obs['population'] = fixed.reindex(a.obs_names).values
        records = []
        normalized = a.X.copy(); normalized.data = np.expm1(normalized.data)
        for pop in predominant:
            mask = a.obs.population.to_numpy() == pop
            case = mask & (a.obs.Library.to_numpy() == 'MPP')
            control = mask & (a.obs.Library.to_numpy() == 'HSC')
            cm = np.asarray(normalized[case].mean(0)).ravel()
            km = np.asarray(normalized[control].mean(0)).ravel()
            fc = pd.Series(np.log2((cm+1)/(km+1)),index=a.var_names)
            s = stats[stats.population==pop].copy()
            eligible = pd.Series(np.asarray((a.X[mask] > .1).sum(0)).ravel() >= 2,index=a.var_names)
            if method == 'uncorrected': baseline_eligible[pop] = eligible
            use = s.gene.map(baseline_eligible[pop]).to_numpy(dtype=bool)
            assert np.isfinite(s.loc[use,'pval']).all()
            common_fdr = np.ones(len(s))
            common_fdr[use] = diff._bh_fdr(s.loc[use,'pval'].to_numpy())
            common_family.append(dict(method=method,population=pop,
                    method_BH_eligible_genes=int(eligible.sum()),
                    fixed_BH_eligible_genes=int(use.sum()),
                    fixed_family_DEGs=int(((common_fdr<.05)&(s.log2fc.abs()>np.log2(1.2))).sum())))
            s['CP10k_log2FC'] = s.gene.map(fc)
            records.append(s)
            secondary.append(dict(method=method,population=pop,
                 CP10k_DEGs=int(((s.fdr<.05)&(s.CP10k_log2FC.abs()>np.log2(1.2))).sum()),
                 median_abs_CP10k_log2FC=float(s.CP10k_log2FC.abs().median())))
        pd.concat(records).to_csv(OUT/method/'normalized_effect_sensitivity.tsv.gz',sep='\t',index=False)
        del a, normalized, stats
    de = pd.DataFrame(de_rows).pivot(index='method',columns='population',values='DEGs').reindex(METHODS)
    de.to_csv(OUT/'DEGs_by_state.tsv',sep='\t')
    secondary = pd.DataFrame(secondary)
    secondary.to_csv(OUT/'normalized_effect_sensitivity.tsv',sep='\t',index=False)
    common_family = pd.DataFrame(common_family)
    common_family.to_csv(OUT/'fixed_BH_family_sensitivity.tsv',sep='\t',index=False)
    fig,axes = plt.subplots(1,3,figsize=(15,7))
    xs = np.arange(len(METHODS))
    for i,pop in enumerate(predominant):
        color = ['#0072B2','#D55E00'][i]
        axes[0].plot(xs,de[pop],marker='o',color=color,label=pop)
        for seed in [0,1,2]:
            points = d[(d.population==pop)&(d.seed==seed)].set_index('method').loc[METHODS]
            axes[1].scatter(xs+(seed-1)*.10,points.UMAP_distance_within_scaled,
                            marker=['o','s','^'][seed],color=color,alpha=.7,s=27,
                            label=pop if seed==0 else None)
    for cap in ['HSC','MPP']:
        axes[2].plot(xs,loss[loss.capture==cap].set_index('method').loc[METHODS].loss_percent,marker='o',label=cap)
    for ax,title,ylabel in zip(axes,['Capture-associated DEGs','Same-state UMAP separation','Transcript removal'],
                         ['DEG gene–state pairs','Centroid distance / within-state RMS radius','Removed counts (%)']):
        ax.set_title(title); ax.set_ylabel(ylabel); ax.set_xticks(xs,METHODS,rotation=60,ha='right',fontsize=8)
        ax.legend(fontsize=9); ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Do filtered-only Python corrections reduce capture differences in ND20_167?',fontsize=14)
    fig.subplots_adjust(left=.07,right=.98,top=.86,bottom=.37,wspace=.35)
    figure_save(fig,'benchmark_comparison',
      'Human bone marrow, patient ND20_167; HSC and MPP are distinct sorted capture gates, one donor, '
      'not independent biological replicates. Left: cellHarmony two-sided Wilcoxon, BH within each '
      'state over features detected >0.1 log(CP10k+1) in at least two cells; FDR <0.05 and '
      '|log2((mean corrected UMI MPP +1)/(mean corrected UMI HSC +1))| >log2(1.2), matching scALABLE. '
      'HSC-2: 4,060 HSC/2,280 MPP cells; MPP-1: 2,738/3,780 cells. Center: three layout seeds '
      'per state (circle/square/triangle = seeds 0/1/2), no confidence intervals or donor replication; '
      'lower is less separation. Blue/orange denote HSC-2/MPP-1 at left and center, HSC/MPP capture at right. '
      'Right: total count loss across all filtered cells, measured input/corrected counts. '
      'All methods share fixed baseline labels and QC cells; all predominant states are shown. '
      'Gate biology and expression-defined labels prevent attributing all differences to ambient RNA.')
    fixed = pd.read_csv(OUT/'fixed_population_assignments.tsv',sep='\t',index_col=0).iloc[:,0]
    for pop in predominant:
        fig,axes=plt.subplots(3,3,figsize=(13,13))
        for ax,method in zip(axes.ravel(),METHODS):
            coords = pd.read_csv(OUT/method/'umap_coordinates.tsv',sep='\t')
            ax.scatter(coords.UMAP1,coords.UMAP2,c='#dddddd',s=.4,rasterized=True)
            means=[]
            for cap,color in [('HSC','#0072B2'),('MPP','#D55E00')]:
                part=coords[(coords.population==pop)&(coords.capture==cap)]
                ax.scatter(part.UMAP1,part.UMAP2,c=color,s=1.3,alpha=.5,rasterized=True,label=cap)
                means.append(part[['UMAP1','UMAP2']].mean().values)
                ax.scatter(*means[-1],c=color,s=65,marker='X',edgecolor='black',linewidth=.5)
            ax.plot(*np.asarray(means).T,c='black',lw=1)
            ax.set_title(method); ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
            ax.legend(markerscale=4,fontsize=8,loc='best')
        fig.suptitle(f'{pop}: do HSC- and MPP-gated cells occupy the same expression space?',fontsize=14)
        fig.subplots_adjust(left=.07,right=.97,top=.92,bottom=.16,hspace=.4,wspace=.3)
        figure_save(fig,'UMAP_'+pop,
          f'ND20_167 human marrow; {pop}: HSC n={populations.loc[pop,"HSC"]:,}, MPP '
          f'n={populations.loc[pop,"MPP"]:,} cells, one donor and two capture libraries. '
          'Blue/orange are HSC/MPP cells assigned to this state in uncorrected cellHarmony; gray '
          'shows every other QC cell. X marks capture centroids and the line joins them. '
          'Each method is embedded separately with the same reference genes, log(CP10k+1), '
          '50 PCs, 15 neighbors, spectral UMAP, seed 0; no integration. Coordinates have arbitrary '
          'rotation and scale; compare within-map overlap rather than absolute positions across '
          'panels. Corrected expression is inferred; labels are expression-derived. No '
          'significance tests or confidence intervals are shown. Capture biology and one donor '
          'limit ambient-specific conclusions; all tested methods are displayed.')
    source_paths = [REF, ROOT/'ambient_subtract.py',ROOT.parent/'cellHarmony/cellHarmony_lite.py',
                    ROOT.parent/'cellHarmony/cellHarmony_differential.py',ROOT.parent/'cellHarmony/flask/pipeline.py',
                    ROOT/'evaluation/benchmark_nd20_167.py',ROOT/'evaluation/summarize_nd20_167.py']
    provenance = dict(python=sys_version(),scanpy=sc.__version__,anndata=ad.__version__,numpy=np.__version__,
        platform=platform.platform(),git_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
        raw_source_files=[str(BASE/f'ND20_167_{cap}/outs/filtered_feature_bc_matrix.h5') for cap in ['HSC','MPP']],
        soup_source_files=[str(BASE/f'ND20_167_{cap}/outs/soupX-contamination-fraction-{rho}/matrix.mtx')
                           for cap in ['HSC','MPP'] for rho in ['0.10','0.15','0.20','0.25']])
    # Preserve the initial analysis manifest when summaries are regenerated.
    if not (OUT/'provenance.json').exists():
        (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    transitions = []; state_stability = []
    for method in METHODS:
        assigned = pd.read_csv(OUT/method/'cellHarmony_lite_assignments.txt',sep='\t').set_index('CellBarcode').iloc[:,0].reindex(fixed.index)
        assert assigned.notna().all()
        for pop in predominant:
            chosen = fixed == pop
            state_stability.append(dict(method=method,population=pop,reassigned_percent=100*float((assigned[chosen]!=pop).mean())))
        table = pd.DataFrame({'baseline':fixed,'assigned':assigned,'capture':fixed.index.str.rsplit('.',n=1).str[-1]})
        transitions.extend(table.groupby(['capture','baseline','assigned']).size().rename('cells').reset_index().assign(method=method).to_dict('records'))
    pd.DataFrame(transitions).to_csv(OUT/'population_transitions.tsv',sep='\t',index=False)
    state_stability = pd.DataFrame(state_stability).pivot(index='method',columns='population',values='reassigned_percent').reindex(METHODS)
    state_stability.to_csv(OUT/'primary_state_reassignment.tsv',sep='\t')
    primary_table = summary[['predominant_DEG_gene_state_pairs','predominant_unique_DEGs',
             'DEG_reduction_percent','UMAP_within_scaled','UMAP_reduction_percent','PCA_within_scaled','PCA_mixing']].reset_index()
    from .write_nd20_167_report import main as write_report
    write_report()
    print(primary_table.to_string(index=False))

def sys_version():
    import sys
    return sys.version

if __name__ == '__main__': main()
