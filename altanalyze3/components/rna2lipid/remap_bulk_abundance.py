"""Invert the supplied bulk transformation and map it to the sorted MS1 gauge."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from reference_normalization import preservation_gate, population_differentials


def invert_bulk(z):
    values=np.asarray(z,dtype=float)
    if not np.isfinite(values).all() or (values<=0).any():
        raise ValueError('Bulk log2(1+10*A) values must be finite and positive')
    return np.log2(np.expm1(values*np.log(2)))-np.log2(10)


def run(source,reference,ms1_reference,out):
    source,reference,ms1_reference,out=map(Path,[source,reference,ms1_reference,out]);out.mkdir(parents=True,exist_ok=True)
    final=pd.read_csv(source/'Final_bulk_Bulk_lipids_cleaned_normalized_median_527_log210.csv',index_col=0)
    inverted=pd.DataFrame(invert_bulk(final),index=final.index,columns=final.columns)
    inverted.to_csv(out/'supplied_bulk_all_features_native_log2.csv')
    np.exp2(inverted).to_csv(out/'supplied_bulk_all_features_native_linear.csv')
    original=pd.read_csv(reference/'legacy_preprocessing_audit/bulk_native_log2_mode_preserved.csv',index_col=0)
    offsets=pd.read_csv(ms1_reference/'calibration_offsets.csv',index_col=0).added_common_log2_offset
    rows=[];columns={}
    for lipid in offsets.index:
        name=lipid[:-2]
        if name not in inverted:
            rows.append({'lipid':lipid,'status':'absent_from_final_bulk'});continue
        candidates=[c for c in original if c[:-2]==name and c.endswith(('_P','_N'))]
        profiles=inverted.index.intersection(original.index)
        scores=[]
        for candidate in candidates:
            error=(inverted.loc[profiles,name]-original.loc[profiles,candidate]).abs()
            scores.append((float((error<=.101).mean()),float(error.mean()),candidate))
        accepted=[s for s in scores if s[0]>=.99]
        if len(accepted)!=1 or accepted[0][2]!=lipid:
            rows.append({'lipid':lipid,'status':'ion_mode_or_source_values_not_uniquely_verified',
                         'candidate_errors':str(scores)});continue
        columns[lipid]=name
        rows.append({'lipid':lipid,'status':'mapped','bulk_column':name,
                     'fraction_source_values_within_0_101_log2':accepted[0][0],
                     'mean_source_error_log2':accepted[0][1],'added_MS1_log2_offset':float(offsets[lipid])})
    if not columns:raise ValueError('No independently verified bulk annotation matches')
    pd.DataFrame(rows).to_csv(out/'mapping_and_exclusions.csv',index=False)
    recovered=pd.DataFrame({lipid:inverted[name] for lipid,name in columns.items()})
    mapped=recovered+offsets.loc[recovered.columns]
    recovered.to_csv(out/'supplied_bulk_inverted_native_log2.csv')
    mapped.to_csv(out/'supplied_bulk_MS1_log2.csv');np.exp2(mapped).to_csv(out/'supplied_bulk_MS1_linear.csv')
    # Combine the actual supplied final bulk profiles with the sorted reference,
    # keeping a donor's unsorted and sorted-PMX assays as distinct observations.
    existing=pd.read_csv(ms1_reference/'combined_reference_MS1_log2.csv',index_col=0)
    metadata=pd.read_csv(reference/'sample_metadata.csv',index_col=0)
    sorted_rows=metadata.index[metadata.dataset=='sorted']
    sorted_frame=existing.loc[sorted_rows,list(columns)]
    combined=pd.concat([sorted_frame,mapped.rename(index=lambda s:s+'_BULK')])
    combined.to_csv(out/'combined_sorted_and_supplied_bulk_MS1_log2.csv')
    np.exp2(combined).to_csv(out/'combined_sorted_and_supplied_bulk_MS1_linear.csv')
    md=metadata.loc[sorted_rows].copy()
    bmd=pd.DataFrame({'donor':mapped.index,'dataset':'bulk','population':'BULK'},index=mapped.index+'_BULK')
    md=pd.concat([md,bmd]);md.to_csv(out/'sample_metadata.csv')
    native=pd.concat([sorted_frame,recovered.rename(index=lambda s:s+'_BULK')])
    roundtrip=pd.read_csv(out/'combined_sorted_and_supplied_bulk_MS1_log2.csv',index_col=0)
    gate=preservation_gate(native,roundtrip,md)
    if not gate['pass']:raise ValueError('Supplied bulk remapping changed within-dataset fold differences')
    linear=pd.read_csv(out/'combined_sorted_and_supplied_bulk_MS1_linear.csv',index_col=0)
    if not (linear.to_numpy()>0).all() or not np.isfinite(linear.to_numpy()).all():raise ValueError('Invalid positive abundance export')
    np.testing.assert_allclose(np.log2(linear),roundtrip,atol=1e-12,rtol=0)
    result={'bulk_profiles':len(mapped),'sorted_profiles':len(sorted_frame),'combined_profiles':len(combined),
            'supplied_bulk_features':len(final.columns),'bulk_features_without_shared_MS1_calibration':len(final.columns)-len(columns),
            'mapped_lipids':len(columns),'excluded_candidate_lipids':len(offsets)-len(columns),
            'bulk_input_scale':'z = log2(1 + 10*A); independently checked against original bulk log2 values',
            'inverse':'x = log2((2**z - 1)/10)',
            'mapped_abundance':'A_MS1 = ((2**z - 1)/10) * 2**common_lipid_offset',
            'units':'reference-normalized MS1 peak intensity, same scale as sorted MS1 reference and predictions',
            'differential_gate':gate,'cross_dataset_interpretation':'Reference alignment assumes transferable PMX composition; not measured physical concentration equivalence',
            'source_exclusion_policy':'Require a unique ion-mode source with at least 99% of matched values within 0.101 log2; report all exclusions'}
    (out/'remapping_result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',required=True);p.add_argument('--reference',required=True)
    p.add_argument('--ms1-reference',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(a.source,a.reference,a.ms1_reference,a.out)
