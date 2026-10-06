"""Donor-level whole-lung RNA aggregation for the tissue lipid validation scope.

Sum actual RNA counts across all observed states/libraries, then within donor.
Apply the evaluation bundles and the unchanged production bundle to those same
RNA profiles. This complements, rather than substitutes for, cell-state DE.
"""

try:
    from .candidate_integrity import reject_retired_workflow
except ImportError:
    from candidate_integrity import reject_retired_workflow

if __name__ == "__main__":
    reject_retired_workflow('evaluate_lungmap_whole_lung.py')

import argparse
import contextlib
import json
import pickle
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits

from evaluate_lungmap_ipf import (HERE, TRANSFER, CANDIDATES, norm_counts,
    sum_groups, resolve_obs, load_engine, arithmetic_stats, lipid_mapping, dump, sha)
from altanalyze3.components.rna2lipid.api import Rna2LipidBundle


def run(out):
    reject_retired_workflow('evaluate_lungmap_whole_lung.py:run')
    registry = pd.read_csv(out/'provenance/IPF_vs_healthy_registry.tsv',sep='\t')
    meta = pd.read_csv(out/'provenance/harmonized_library_metadata_harmonized_final_corrected_v6.txt',sep='\t',dtype=str).fillna('')
    cross=pd.read_csv(out/'provenance/crosswalk.tsv',sep='\t')
    mapping=dict(zip(cross.cell_state,cross.cell_type))|dict(zip(cross.cell_type_full,cross.cell_type))|dict(zip(cross.cell_type,cross.cell_type))
    source=TRANSFER/'reference_v7_top500/pseudobulk/CellRef2.0_pseudobulk_library_x_cellstate_UNION.h5ad'
    a=ad.read_h5ad(source)
    obs=resolve_obs(a,meta,'PB',mapping)
    keys=obs.Study_internal.astype(str)+'|'+obs.Sample.astype(str)
    raw, names, _, design=sum_groups(a.X,keys)
    first=obs.loc[~keys.duplicated()].copy();first.index=keys.loc[~keys.duplicated()]
    smeta=first.loc[names].copy()
    smeta['n_cells']=np.asarray(design@obs.n_cells.to_numpy()).ravel()
    genes=a.var_names.copy();del a
    engine=load_engine(out)
    production=Rna2LipidBundle.load(HERE/'rna2lipid_hs_lung_lipidwise_bundle.pkl')
    original_hash=sha(production.bundle_path)
    fitted={label:pickle.load((out/label/'candidate_bundle.pkl').open('rb')) for label in CANDIDATES}
    needed=list(production.input_genes)
    available=[g for g in needed if g in genes]
    positions=genes.get_indexer(available)
    measured=pd.read_csv(out/'experimental_lung_IPF_lipids_all_544.csv',index_col=0)
    native=pd.read_csv(HERE/'artifacts/reference_normalization_20261001/native_lipid_targets.csv',index_col=0).columns
    destination=out/'whole_lung_donor_validation';destination.mkdir(exist_ok=True)
    outputs,replicates,match_tables=[],[],[]
    for row in registry.itertuples():
        def select(value):return {'|'.join(k.split('|')[-2:]) for k in str(value).split(';') if k and k!='nan'}
        ca,co=names.isin(select(row.case_keys)),names.isin(select(row.control_keys))
        inds=np.flatnonzero(ca|co)
        donor=smeta.iloc[inds].Study_internal.astype(str)+'|'+smeta.iloc[inds].Donor.astype(str)
        # Disease arms must never merge into one donor profile.
        arm=pd.Series(np.where(ca[inds],'CASE','CONTROL'),index=smeta.index[inds])
        if arm.groupby(donor).nunique().max()>1:raise ValueError('A donor occurs in both IPF/control arms')
        donor_raw,dnames,_,_=sum_groups(raw[inds],donor)
        first=~donor.duplicated()
        md=smeta.iloc[inds].loc[first].copy();md.index=donor.loc[first];md=md.loc[dnames]
        md['Condition']=arm.groupby(donor).first().loc[dnames]
        expression=pd.DataFrame(norm_counts(donor_raw)[:,positions].toarray(),index=dnames,columns=available)
        directory=destination/row.contrast_id;directory.mkdir(exist_ok=True)
        expression.to_csv(directory/'donor_RNA_log2CP10k.csv')
        md.to_csv(directory/'donor_metadata.csv')
        ca=md.Condition.eq('CASE').to_numpy();co=md.Condition.eq('CONTROL').to_numpy()
        ncase,nctrl=int(ca.sum()),int(co.sum())
        print(f'Whole lung {row.contrast_id}: {ncase} IPF / {nctrl} control donors',flush=True)
        predictions={'current_explicit_log2':production.predict_from_dataframe(expression).predictions}
        for label,b in fitted.items():
            columns=b['X_columns']
            offsets=pd.read_csv(out/label/'healthy_study_RNA_offsets.csv',index_col=0).loc[md.Study_internal,columns]
            aligned=expression[columns].to_numpy()+offsets.to_numpy()
            vals=b['scaler_y'].inverse_transform(b['model'].predict(b['scaler_x'].transform(aligned)))
            predictions[label]=pd.DataFrame(vals,index=expression.index,columns=b['Y_columns'])
        for label,pred in predictions.items():
            pred.to_csv(directory/f'{label}_predictions_log2.csv')
            np.exp2(pred).to_csv(directory/f'{label}_predictions_linear.csv')
            # Ordering keeps CASE rows together without losing identities.
            order=np.r_[np.flatnonzero(ca),np.flatnonzero(co)]
            vals=pred.iloc[order].to_numpy()
            block=ad.AnnData(vals,obs=pd.DataFrame({'Condition':pd.Categorical(['CASE']*ncase+['CONTROL']*nctrl)},index=md.index[order]),
                            var=pd.DataFrame(index=pred.columns))
            block.uns['log1p']={'base':2.0}
            exact=arithmetic_stats(vals,ncase,pred.columns)
            if min(ncase,nctrl)>=3:
                with (directory/f'{label}_cellHarmony.log').open('w') as log,contextlib.redirect_stdout(log):
                    result,tested=engine._moderated_t_test(block,'Condition','CASE','CONTROL','whole_lung_RNA_aggregate')
                result=result.set_index('gene').drop(columns=['log2fc'],errors='ignore').join(exact)
                status='tested'
            else:
                result=exact.copy();result['pval']=np.nan;result['fdr']=np.nan;tested=0
                status='fold_only_below_three_donors_per_arm'
            result['signed_fold']=np.where(result.log2fc>0,np.exp2(result.log2fc),-np.exp2(-result.log2fc))
            result.loc[result.log2fc.eq(0),'signed_fold']=0
            result['n_case'],result['n_control'],result['model']=ncase,nctrl,label
            result['contrast'],result['status']=row.contrast_id,status
            result.index.name='gene'
            result.to_csv(directory/f'{label}_all_lipid_tests.csv')
            outputs.append(result.reset_index())
            replicates.append({'contrast':row.contrast_id,'model':label,'IPF_donors':ncase,'control_donors':nctrl,'status':status,'tested_features':tested})
            if row.contrast_id==registry.contrast_id.iloc[0]:match_tables.append(lipid_mapping(pred.columns,label,measured,native))
    full=pd.concat(outputs,ignore_index=True)
    full.to_csv(destination/'all_tests.csv',index=False)
    matches=pd.concat(match_tables,ignore_index=True)
    joined=full.merge(matches[matches.status.eq('matched')][['model','gene','measured_feature']],on=['model','gene']).merge(
        measured.rename_axis('measured_feature').reset_index()[['measured_feature','supplied_raw_p','supplied_adjusted_p','supplied_log_effect','pvalue','FDR','geometric_log2FC']],on='measured_feature')
    joined['direction_concordant_supplied']=joined.log2fc*joined.supplied_log_effect>0
    joined['direction_concordant_donor']=joined.log2fc*joined.geometric_log2FC>0
    joined.to_csv(destination/'measured_vs_imputed_all_tests.csv',index=False)
    shared47=set(matches.loc[matches.model.eq('candidate_MS1_47'),'measured_feature'].dropna()) & set(matches.loc[matches.model.eq('current_explicit_log2'),'measured_feature'].dropna())
    shared219=set(matches.loc[matches.model.eq('candidate_bulk_219_exploratory'),'measured_feature'].dropna()) & set(matches.loc[matches.model.eq('current_explicit_log2'),'measured_feature'].dropna())
    summary=[]
    for (contrast,label),group in joined.groupby(['contrast','model']):
        for source_gate,col,cut,direction in [('supplied_rawp005','supplied_raw_p',.05,'direction_concordant_supplied'),('supplied_BH010','supplied_adjusted_p',.1,'direction_concordant_supplied'),('donor_rawp005','pvalue',.05,'direction_concordant_donor'),('donor_BH010','FDR',.1,'direction_concordant_donor')]:
            for panel,allowed in [('all_matched',set(group.measured_feature)),('shared_MS1_47',shared47),('shared_bulk_219',shared219)]:
                if (panel=='shared_MS1_47' and label=='candidate_bulk_219_exploratory') or (panel=='shared_bulk_219' and label=='candidate_MS1_47'):continue
                select=group[group[col].lt(cut)&group.measured_feature.isin(allowed)]
                for gate,field,threshold in [('rawp005','pval',.05),('BH005','fdr',.05),('BH010','fdr',.1)]:
                    sig=select[select[field].lt(threshold)];yes=sig[sig[direction]];no=sig[~sig[direction]]
                    summary.append({'contrast':contrast,'model':label,'source_gate':source_gate,'panel':panel,'prediction_gate':gate,
                        'matched_experimental_differentials':len(select),'concordant_up':int(yes.log2fc.gt(0).sum()),'concordant_down':int(yes.log2fc.lt(0).sum()),
                        'concordant_total':len(yes),'discordant_total':len(no),'all_directions_concordant':int(select[direction].sum()),
                        'n_case_donors':int(group.n_case.iloc[0]),'n_control_donors':int(group.n_control.iloc[0]),'status':group.status.iloc[0]})
    pd.DataFrame(summary).to_csv(destination/'directional_validation_summary.csv',index=False)
    pd.DataFrame(replicates).to_csv(destination/'replication.csv',index=False)
    dump(destination/'method.json',{'RNA_aggregation':'Sum raw counts across measured lung cell states and libraries, then within Study_internal/Donor in each contrast arm, before CP10k/log2 normalization',
        'differential':'Frozen cellHarmony eBayes, minimum three donors/arm, alpha .05, no fold cutoff, same protocol for all models',
        'current_baseline':'Unchanged production bundle, existing alignment policy (1301/1303 genes; its two absent genes receive its default zero input), explicit correct log2 abundance for DE',
        'candidate_baseline':'Same 1293-gene candidate fits and healthy-study offsets as cell-state evaluation; no additional fitting or disease-outcome selection',
        'interpretation':'Whole-lung RNA sum of captured cells, not measured whole tissue composition; no matched RNA/lipid donors',
        'production_bundle_sha256_before':original_hash,'production_bundle_sha256_after':sha(production.bundle_path)})
    if original_hash!=sha(production.bundle_path):raise ValueError('Production bundle changed')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',default=str(HERE/'artifacts/LungMAP_IPF_candidate_20261004'))
    with threadpool_limits(limits=2):run(Path(p.parse_args().out))
