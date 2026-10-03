"""Compare each CITE source's ML/CLP and DN annotations on the same flow events.

Writes full FlowSOM distributions, explicit absent labels, source-cell DN/marrow
cross-tabs, exact antibody mappings, and optional fixed-polarity gate evaluations.
"""
import argparse
import json
from pathlib import Path
import urllib.request
import urllib.error
import numpy as np
import pandas as pd
from ..visualization.flow_viewer.data_api import FlowBundle
from .crosswalk import build_crosswalk


def run(bundle, out, url=None):
    b=FlowBundle(bundle);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    flow_levels=np.asarray(b.levels('FlowSOM_J8DW'));flow_codes=b.labels('FlowSOM_J8DW')
    input_names=list(np.load(Path(bundle)/'validation_20261002/flow_rds_input.npz')['channels'])
    distributions=[];status=[];source_rows=[];maps=[];gates=[];flow_dn=[]
    targets=['ML-1a','ML-1b','CLP1-a','CLP1-b','CLP1-c']
    for space in ['cite_grimes','cite_chinese']:
        spec=b.space(space);features=b.features(space)
        pairs,missing,_=build_crosswalk(input_names,[f for f in features if not f.startswith('RNA:')])
        pairs=pairs.drop_duplicates('cite').drop_duplicates('flow')
        maps.extend(dict(source=space,**r) for r in pairs.to_dict('records'))
        semantic={}
        if 'cluster_name' in spec['labels']:
            cluster=np.asarray(b.levels('cluster_name',space))[b.labels('cluster_name',space)]
            pruned=np.asarray(b.levels('pruned',space))[b.labels('pruned',space)]
            for label in np.unique(pruned):
                names,counts=np.unique(cluster[pruned==label],return_counts=True)
                semantic[label]=str(names[counts.argmax()])
        for annotation,prefix in spec.get('transfer_links',{}).items():
            if annotation=='celltype':continue
            ref_levels=b.levels(annotation,space);ref_codes=b.labels(annotation,space)
            selections=targets if annotation=='Mm-MarrowAtlas-L4' else [x for x in ref_levels if 'DN' in semantic.get(x,x)]
            if annotation=='pruned':selections += [x for x in ref_levels if any(t in semantic.get(x,x) for t in targets)]
            for key in b.space('flow')['labels']:
                if not key.startswith(prefix):continue
                method=key[len(prefix):];levels=b.levels(key);codes=b.labels(key)
                for population in dict.fromkeys(selections):
                    nsource=int((ref_codes==ref_levels.index(population)).sum()) if population in ref_levels else 0
                    mask=codes==levels.index(population) if population in levels else np.zeros(b.n(),bool)
                    n=int(mask.sum());name=semantic.get(population,population) if annotation=='pruned' else population
                    state='absent_in_source' if not nsource else ('not_assigned_in_flow' if not n else 'assigned')
                    status.append(dict(source=space,annotation=annotation,method=method,label_set=key,population=population,display_name=name,n_source_cells=nsource,n_flow_events=n,status=state))
                    counts=np.bincount(flow_codes[mask & (flow_codes>=0)],minlength=len(flow_levels))
                    for cluster,count in zip(flow_levels,counts):
                        distributions.append(dict(source=space,annotation=annotation,method=method,population=population,display_name=name,FlowSOM_cluster=cluster,n_events=int(count),percent_of_population=100*count/n if n else np.nan))
                    if url and annotation=='Mm-MarrowAtlas-L4' and population in ['ML-1a','CLP1-a','CLP1-b','CLP1-c'] and method in ['kde_cellharmony','kde_k5_average']:
                        row=dict(source=space,method=method,population=population,status=state)
                        if n:
                            body=dict(space='flow',label_set=key,populations=[population],signs={'CD4':'<=','CD8':'<=','CD117':'>','CD25':'<=','CD11c':'<=','CD11b':'<=','CD27':'>','Sca-1':'>'})
                            try:
                                request=urllib.request.Request(url+'/api/optimize',data=json.dumps(body).encode(),headers={'Content-Type':'application/json'})
                                result=json.load(urllib.request.urlopen(request,timeout=120));row.update(result['best']['test']);row['conditions']=json.dumps(result['best']['conditions']);(out/(key+'_'+population+'_gate.json')).write_text(json.dumps(result,indent=2))
                            except urllib.error.HTTPError as e:row['error']=e.read().decode()
                        gates.append(row)
            if annotation in ['StJude','Author_celltype']:
                marrow=np.asarray(b.levels('Mm-MarrowAtlas-L4',space))[b.labels('Mm-MarrowAtlas-L4',space)]
                for population in selections:
                    mask=ref_codes==ref_levels.index(population);names,counts=np.unique(marrow[mask],return_counts=True)
                    for name,count in zip(names,counts):source_rows.append(dict(source=space,DN_label=population,marrow_label=name,n_CITE_cells=int(count),percent_of_DN=100*count/mask.sum()))
    for name,rows in [('population_status',status),('FlowSOM_distributions',distributions),('source_DN_marrow_associations',source_rows),('assignment_antibody_crosswalk',maps),('eight_marker_gate_validation',gates)]:
        pd.DataFrame(rows).to_csv(out/(name+'.tsv'),sep='\t',index=False)
    for space,annotation in [('cite_grimes','StJude'),('cite_chinese','Author_celltype')]:
        links=b.space(space)['transfer_links'];prefix=links[annotation]
        for key in b.space('flow')['labels']:
            if not key.startswith(prefix):continue
            method=key[len(prefix):];marrow_key=links['Mm-MarrowAtlas-L4']+method
            if marrow_key not in b.space('flow')['labels']:continue
            dn_codes=b.labels(key);marrow_codes=b.labels(marrow_key);marrow_levels=b.levels(marrow_key)
            for i,label in enumerate(b.levels(key)):
                if 'DN' not in label:continue
                mask=dn_codes==i;total=int(mask.sum())
                counts=np.bincount(marrow_codes[mask & (marrow_codes>=0)],minlength=len(marrow_levels))
                for population,count in zip(marrow_levels,counts):
                    flow_dn.append(dict(source=space,method=method,DN_label=label,marrow_label=population,n_flow_events=int(count),percent_of_DN=100*count/total if total else np.nan))
    pd.DataFrame(flow_dn).to_csv(out/'flow_DN_marrow_associations.tsv',sep='\t',index=False)
    df=pd.DataFrame(distributions);dom=df[df.n_events>0].sort_values('n_events',ascending=False).groupby(['source','annotation','method','population'],sort=False).head(1)
    dom.to_csv(out/'dominant_FlowSOM.tsv',sep='\t',index=False)
    (out/'interpretation.txt').write_text('Each source and annotation is evaluated separately on FlowSOM_J8DW. Full cluster distributions and absent labels are retained. Source DN/marrow associations refer to the same measured CITE cells; transferred DN/FlowSOM associations refer to separate flow events. Neither is proof of stem-cell function. Shared-marker mapping is from the RDS training matrix; CD4/CD8 display/gating channels were not used in these transfers. Gate precision is agreement with transferred labels, not experimental sorting purity. TCR_chain→TCRbeta is an existing assumed alias requiring panel/clone verification; integrin-7 is unmatched, not silently merged. Chinese values are DSB-normalized; Grimes are TotalVI. Raw cross-panel thresholds are not comparable.\n')
    print('Wrote source comparison:',out,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bundle',required=True);p.add_argument('--out',required=True);p.add_argument('--url');a=p.parse_args();run(a.bundle,a.out,a.url)

if __name__=='__main__':main()
