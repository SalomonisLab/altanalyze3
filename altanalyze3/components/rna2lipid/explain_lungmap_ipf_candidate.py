"""Exact gene decomposition of selected candidate mean-log lipid differences."""
from pathlib import Path
import pickle
import anndata as ad
import numpy as np
import pandas as pd
from evaluate_lungmap_ipf import TRANSFER, resolve_obs, sum_groups, norm_counts, pairs_for
from rank_lungmap_ipf_comparisons import OUT


def explain(out=OUT):
    label='candidate_MS1_47'
    ranks=pd.read_csv(out/'ranked_dataset_celltype_comparisons.csv')
    eligible=ranks[(ranks.model==label)&(ranks.unit=='PB')&(ranks.panel=='shared_MS1_47')&(ranks.prediction_gate=='BH010')&ranks.donor_supported]
    selected=eligible.groupby('contrast',sort=False).head(2)
    source=ad.read_h5ad(TRANSFER/'reference_v7_top500/pseudobulk/CellRef2.0_pseudobulk_library_x_cellstate_UNION.h5ad')
    meta=pd.read_csv(out/'provenance/harmonized_library_metadata_harmonized_final_corrected_v6.txt',sep='\t')
    cross=pd.read_csv(out/'provenance/crosswalk.tsv',sep='\t')
    mapping={}
    for c in ['cell_state','cell_type_full','cell_type']:mapping.update(dict(zip(cross[c],cross.cell_type)))
    obs=resolve_obs(source,meta,'PB',mapping)
    keys=obs.Study_internal.astype(str)+'|'+obs.Sample.astype(str)+'|'+obs.cell_state.astype(str)
    raw,names,_,_=sum_groups(source.X,keys)
    b=pickle.load((out/label/'candidate_bundle.pkl').open('rb'))
    genes=list(b['X_columns']);lipids=list(b['Y_columns'])
    x=norm_counts(raw)[:,source.var_names.get_indexer(genes)].toarray()
    pred=ad.read_h5ad(out/label/'LungMAP_PB_predictions_log2.h5ad')
    if not names.equals(pred.obs_names):raise ValueError('RNA/prediction alignment mismatch')
    w=(b['model'].X_fit_.T@b['model'].dual_coef_)*b['scaler_y'].scale_[None,:]/b['scaler_x'].scale_[:,None]
    registry=pd.read_csv(out/'provenance/IPF_vs_healthy_registry.tsv',sep='\t').set_index('contrast_id')
    matches=pd.read_csv(out/'all_measured_vs_imputed_tests.csv')
    checks=[];drivers=[];allcontrib=[]
    for row in selected.itertuples():
        r=registry.loc[row.contrast].copy();r['contrast_id']=row.contrast
        pairs,_=pairs_for(pred.obs,r,pd.DataFrame(columns=['cell2_disease_state','cell1_reference_state']),3)
        _,ca,co=next(p for p in pairs if p[0]==row.population)
        delta=x[ca].mean(0)-x[co].mean(0)
        g=matches[(matches.contrast==row.contrast)&(matches.unit=='PB')&(matches.model==label)&(matches.population==row.population)&(matches.fdr<.1)&matches.direction_concordant_supplied]
        for lipid in g.itertuples():
            j=lipids.index(lipid.gene);contrib=delta*w[:,j]
            actual=float(pred.X[ca,j].mean()-pred.X[co,j].mean())
            err=abs(contrib.sum()-actual)
            if err>1e-10:raise ValueError('Gene decomposition does not reproduce prediction')
            checks.append(dict(contrast=row.contrast,population=row.population,lipid=lipid.gene,
                               measured_feature=lipid.measured_feature,arithmetic_signed_fold=lipid.signed_fold,
                               geometric_signed_fold=2**actual if actual>0 else -(2**(-actual)),
                               pval=lipid.pval,fdr=lipid.fdr,max_reconstruction_error=err))
            frame=pd.DataFrame(dict(contrast=row.contrast,population=row.population,lipid=lipid.gene,
                                    RNA_gene=genes,mean_log_RNA_difference=delta,model_coefficient=w[:,j],
                                    log2_prediction_contribution=contrib,fold_factor=2**contrib))
            allcontrib.append(frame)
            drivers.append(frame.iloc[np.argsort(-abs(contrib))[:10]])
    checks=pd.DataFrame(checks);drivers=pd.concat(drivers,ignore_index=True)
    checks.to_csv(out/'gene_explanation_reconstruction_audit.csv',index=False)
    pd.concat(allcontrib,ignore_index=True).to_csv(out/'selected_lipid_all_gene_contributions.csv',index=False)
    drivers.to_csv(out/'selected_lipid_top_gene_drivers.csv',index=False)
    lines=['# RNA drivers of selected lipid predictions','',
           'The linear candidate permits an exact algebraic explanation: each gene’s mean log-RNA difference between IPF and control is multiplied by its fitted lipid coefficient. Contributions sum to the difference between mean predicted log2 lipid values. Study offsets and the intercept cancel within each contrast. Reconstruction is checked against the full saved predictions.','',
           'This decomposition explains the geometric fold change. The differential tables use ratios of arithmetic abundance means; these two fold definitions can differ. Fold factors multiply to the unsigned geometric ratio. Coefficients are model associations, not causal effects, and correlated genes can share prediction contributions. Selected dataset/state/lipid combinations are exploratory and were chosen using validation agreement, not held-out outcomes.','',
           '## Reconstruction and original arithmetic folds','',checks.to_markdown(index=False,floatfmt='.5g'),'',
           '## Ten largest absolute gene contributions per selected lipid','',drivers.to_markdown(index=False,floatfmt='.5g'),'']
    (out/'RNA_GENE_EXPLANATIONS.md').write_text('\n'.join(lines))
    print(f'Explained {len(checks)} lipid/state comparisons; maximum reconstruction error {checks.max_reconstruction_error.max():.3g}')

if __name__=='__main__':explain()
