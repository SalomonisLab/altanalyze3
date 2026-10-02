"""Test one-donor raw abundance anchoring against the other two MassIVE donors."""
import argparse
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd


def metrics(frame):
    d=frame.dropna(subset=['predicted_log2','measured_log2'])
    return {'observations':len(d),'rmse_log2':float(np.sqrt(np.mean((d.predicted_log2-d.measured_log2)**2))),
            'median_absolute_error_log2':float(np.median(np.abs(d.predicted_log2-d.measured_log2))),
            'correlation':float(d.predicted_log2.corr(d.measured_log2))}


def run(raw_dir,published_path,reference_dir):
    out=Path(raw_dir);pub=pd.read_csv(published_path,index_col=0)
    pub.index=[re.sub(r'^D0*(\d+)',lambda m:f'D{int(m[1]):03d}',s) for s in pub.index]
    raw=pd.read_csv(out/'raw_apex_log2.csv',index_col=0)
    targets=pd.read_csv(out/'locked_targets.tsv',sep='\t').set_index('feature_id')
    eligible=targets.index[targets.isobaric_targets.isna()&(targets.ms2_support_scans>0)]
    features=sorted(set(eligible)&set(pub.columns)&set(raw.columns))
    normal=raw.sub(raw.median(axis=1),axis=0)
    mode_normal=raw.copy()
    for mode in ['positive','negative']:
        cols=targets.index[targets.polarity==mode].intersection(raw.columns)
        mode_normal.loc[:,cols]=raw[cols].sub(raw[cols].median(axis=1),axis=0)
    scales={'raw_apex_log2':raw,'sample_median_normalized_log2':normal,
            'ion_mode_median_normalized_log2':mode_normal}
    results={};rows=[]
    for scale,measured in scales.items():
        for method in ['global_offset','lipid_common_offset','PMX_lipid_offset','population_specific_offset']:
            records=[];max_fc_change=0
            for donor in ['D001','D008','D011']:
                train=[s for s in pub.index if s.startswith(donor+'_')]
                held=[s for s in pub.index if not s.startswith(donor+'_')]
                residual=measured.loc[train,features]-pub.loc[train,features]
                if method=='global_offset':shift=pd.Series(np.nanmedian(residual.to_numpy()),index=features)
                elif method=='PMX_lipid_offset':shift=residual.loc[donor+'_PMX']
                else:shift=residual.median(axis=0)
                transformed=pub[features]+shift
                if method=='population_specific_offset':
                    for sample in transformed.index:
                        pop=sample.split('_')[1]
                        transformed.loc[sample]=pub.loc[sample,features]+residual.loc[donor+'_'+pop]
                for other in ['D001','D008','D011']:
                    for pop in ['END','EPI','MES','MIC']:
                        s,p=other+'_'+pop,other+'_PMX'
                        delta=(transformed.loc[s]-transformed.loc[p])-(pub.loc[s,features]-pub.loc[p,features])
                        max_fc_change=max(max_fc_change,float(delta.abs().max()))
                for sample in held:
                    for fid in features:
                        records.append({'scale':scale,'method':method,'anchor_donor':donor,'sample':sample,'lipid':fid,
                                        'predicted_log2':transformed.loc[sample,fid],'measured_log2':measured.loc[sample,fid]})
            frame=pd.DataFrame(records);rows.extend(records)
            results[scale+'/'+method]={**metrics(frame),'maximum_changed_log2FC':max_fc_change,
                                       'preserves_population_differentials':max_fc_change<=1e-9}
    pd.DataFrame(rows).to_csv(out/'one_donor_anchor_predictions.csv',index=False)
    # Independently compare the trained reference's donor-held-out RNA prediction
    # contrasts with MS1 re-extraction, on exact, unambiguous ion-mode annotations.
    reference=Path(reference_dir)
    model=pd.read_csv(reference/'MSV000081973_heldout_contrasts.csv')
    mapping={fid:fid.split('|')[0]+('_P' if fid.endswith('+') else '_N') for fid in features}
    inverse={v:k for k,v in mapping.items()};checks=[]
    for r in model.to_dict('records'):
        fid=inverse.get(r['lipid'])
        if fid is None:continue
        measured=normal.loc[r['sample'],fid]-normal.loc[r['reference'],fid]
        checks.append({**r,'MS1_reextracted_log2FC':measured})
    checks=pd.DataFrame(checks).dropna();checks.to_csv(out/'RNA_prediction_vs_reextracted_MS1.csv',index=False)
    strong=checks[checks.MS1_reextracted_log2FC.abs()>=.5]
    prediction={'contrasts':len(checks),'lipids':int(checks.lipid.nunique()),
                'correlation':float(checks.predicted_log2fc.corr(checks.MS1_reextracted_log2FC)),
                'rmse_log2FC':float(np.sqrt(np.mean((checks.predicted_log2fc-checks.MS1_reextracted_log2FC)**2))),
                'strong_contrast_direction_accuracy':float(np.mean(np.sign(strong.predicted_log2fc)==np.sign(strong.MS1_reextracted_log2FC)))}
    result={'eligible_published_targets':len(features),'anchor_tests':results,'RNA_prediction_vs_MS1':prediction,
            'anchor_holdout':'Each of the three donors anchors all five populations; the other two donors are held out',
            'interpretation':'Per-lipid offsets common across populations preserve fold changes; population-specific offsets must pass the same gate',
            'abundance_units':'MS1 peak intensity or reference-normalized intensity, not molar concentration',
            'identity_scope':'Published annotation plus fragment support; isobaric annotations excluded; original manual review not reproduced'}
    (out/'anchor_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir',required=True);parser.add_argument('--published',required=True)
    parser.add_argument('--reference-dir',required=True)
    a=parser.parse_args();run(a.raw_dir,a.published,a.reference_dir)
