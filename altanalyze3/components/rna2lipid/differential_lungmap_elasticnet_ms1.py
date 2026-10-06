"""Frozen cellHarmony tests for same-algorithm MS1 candidate and matched-panel baseline."""

try:
    from .candidate_integrity import reject_retired_workflow
except ImportError:
    from candidate_integrity import reject_retired_workflow

if __name__ == "__main__":
    reject_retired_workflow('differential_lungmap_elasticnet_ms1.py')

import shutil
import time
import numpy as np
import pandas as pd
import anndata as ad
from threadpoolctl import threadpool_limits
from evaluate_lungmap_elasticnet_ms1 import OUT,OLD,LABEL,PRODUCTION,HERE
from evaluate_lungmap_ipf import load_engine,load_current,differential,lipid_mapping,sha,dump


def run(out=OUT):
    reject_retired_workflow('differential_lungmap_elasticnet_ms1.py:run')
    protected=sha(PRODUCTION)
    directory=out/'provenance';directory.mkdir(exist_ok=True)
    for file in ['IPF_vs_healthy_registry.tsv','crosswalk.tsv','overrides.tsv','harmonized_library_metadata_harmonized_final_corrected_v6.txt']:
        shutil.copyfile(OLD/'provenance'/file,directory/file)
    registry=pd.read_csv(directory/'IPF_vs_healthy_registry.tsv',sep='\t').fillna('')
    cross=pd.read_csv(directory/'crosswalk.tsv',sep='\t')
    overrides=pd.read_csv(directory/'overrides.tsv',sep='\t')
    current=load_current(out,cross)
    candidate={u:ad.read_h5ad(out/LABEL/f'LungMAP_{u}_predictions_log2.h5ad') for u in ['PB','MC']}
    lipids=list(candidate['PB'].var_names)
    matched={u:a[:,lipids].copy() for u,a in current.items()}
    for a in candidate.values():a.X=a.X.astype(np.float32) # deployed forDE storage precision
    measured=pd.read_csv(OLD/'experimental_lung_IPF_lipids_all_544.csv',index_col=0)
    native=pd.read_csv(HERE/'artifacts/reference_normalization_20261001/native_lipid_targets.csv',index_col=0).columns
    labels={'current_as_deployed':current,'current_explicit_log2':current,
            'current_matched47_explicit_log2':matched,LABEL:candidate}
    maps=pd.concat([lipid_mapping(a['PB'].var_names,label,measured,native) for label,a in labels.items()],ignore_index=True)
    maps.to_csv(out/'lipid_matching_audit.csv',index=False)
    engine=load_engine(out);frames=[];censuses=[]
    with threadpool_limits(limits=1):
        for row in registry.itertuples():
            ov=overrides[overrides.contrast_id.eq(row.contrast_id)]
            for unit in ['PB','MC']:
                for label,matrices in labels.items():
                    print('DE',row.contrast_id,unit,label,flush=True)
                    if label in ['current_as_deployed','current_explicit_log2']:
                        # Prior baseline was already replayed against the frozen engine and audited.
                        source=OLD/'differentials'/row.contrast_id/unit/label
                        dest=out/'differentials'/row.contrast_id/unit/label
                        shutil.copytree(source,dest,dirs_exist_ok=True)
                        frame=pd.read_csv(dest/'all_tested_lipids.csv')
                        census=pd.read_csv(dest/'population_replication_census.csv')
                    else:
                        frame,census=differential(matrices[unit],unit,row,ov,engine,out,label)
                    frames.append(frame);census['contrast']=row.contrast_id;census['unit']=unit;census['model']=label;censuses.append(census)
    full=pd.concat(frames,ignore_index=True);full.to_csv(out/'all_cellHarmony_tests.csv',index=False)
    census=pd.concat(censuses,ignore_index=True);census.to_csv(out/'population_replication_census.csv',index=False)
    mapping=maps[maps.status.eq('matched')][['model','gene','measured_feature']]
    join=full.merge(mapping,on=['model','gene']).merge(measured.rename_axis('measured_feature').reset_index(),on='measured_feature',suffixes=('','_experimental'))
    join['direction_concordant_supplied']=join.log2fc*join.supplied_log_effect>0
    join.to_csv(out/'all_measured_vs_imputed_tests.csv',index=False)
    measured.to_csv(out/'experimental_lung_IPF_lipids_all_544.csv')
    rows=[];supported=[];bestrows=[]
    shared=set(maps.loc[maps.model.eq(LABEL)&maps.status.eq('matched'),'measured_feature'])
    for gate,source,stat,cut in [('both_rawp005','supplied_raw_p','pval',.05),('both_BH005','supplied_adjusted_p','fdr',.05),('both_BH010','supplied_adjusted_p','fdr',.1)]:
        significant=join[join[source].lt(cut)&join[stat].lt(cut)&join.direction_concordant_supplied&~join.population.str.contains('__vs__')]
        for fold in [1.,1.1,1.2,1.5,2.]:
            selected=significant[significant.signed_fold.abs().ge(fold)&significant.signed_fold.ne(0)]
            for row in registry.itertuples():
                for unit in ['PB','MC']:
                    for label in labels:
                        g=selected[selected.contrast.eq(row.contrast_id)&selected.unit.eq(unit)&selected.model.eq(label)&selected.measured_feature.isin(shared)]
                        up=set(g.loc[g.supplied_log_effect.gt(0),'measured_feature']);down=set(g.loc[g.supplied_log_effect.lt(0),'measured_feature'])
                        rows.append(dict(contrast=row.contrast_id,unit=unit,model=label,overlap_gate=gate,minimum_imputed_fold=fold,
                                         matched_experimental_differentials=int(measured.loc[measured.index.isin(shared),source].lt(cut).sum()),
                                         concordant_union=len(up|down),concordant_up=len(up),concordant_down=len(down),
                                         supporting_cell_states=g.population.nunique()))
                        for feature,x in g.groupby('measured_feature'):
                            supported.append(dict(contrast=row.contrast_id,unit=unit,model=label,overlap_gate=gate,minimum_imputed_fold=fold,
                                                  measured_feature=feature,concordant_cell_types='; '.join(sorted(x.population.unique()))))
                        if fold==1 and gate=='both_BH010' and len(g):
                            statecounts=g.groupby('population').measured_feature.nunique()
                            for state in statecounts.index[statecounts.eq(statecounts.max())]:
                                x=g[g.population.eq(state)]
                                bestrows.append(dict(contrast=row.contrast_id,unit=unit,model=label,population=state,
                                                     concordant_total=int(statecounts[state]),concordant_up=int(x.supplied_log_effect.gt(0).sum()),
                                                     concordant_down=int(x.supplied_log_effect.lt(0).sum())))
    pd.DataFrame(rows).to_csv(out/'cohort_concordant_union_summary.csv',index=False)
    pd.DataFrame(supported).to_csv(out/'cohort_concordant_lipid_celltypes.csv',index=False)
    pd.DataFrame(bestrows).to_csv(out/'best_concordant_celltypes.csv',index=False)
    shutil.copyfile(OLD/'current_saved_vs_rerun_parity.csv',out/'current_saved_vs_rerun_parity.csv')
    after=sha(PRODUCTION)
    if protected!=after:raise ValueError('Production bundle changed')
    dump(out/'differential_method_audit.json',{'engine_sha256':sha(out/'provenance/cellHarmony_differential_frozen.py'),
                'production_sha256_before':protected,'production_sha256_after':after,'production_replaced':False,
                'matched_panel_lipids':lipids,'matched_panel_prior':'current_matched47_explicit_log2',
                'statistical_parameters':'frozen deployed functions; PB shrinkage .2 / minimum 3, MC Wilcoxon tie_correct False / minimum 5; no fc cutoff',
                'postprocessing':'explicit log2 for matched comparison; exact inverse 2**X arithmetic fold; no new abundance pseudocount',
                'deployed_baseline':'Copied previously audited replay, using same protected production bundle; residual MC p-value differences retained in parity audit',
                'training_sample_difference':'45 corrected profiles versus 50 production; five D071 profiles unavailable; no synthetic labels'})
    print('Completed ElasticNet differential comparison',flush=True)

if __name__=='__main__':run()
