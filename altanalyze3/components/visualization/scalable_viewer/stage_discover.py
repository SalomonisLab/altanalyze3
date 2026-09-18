"""Stage a viewer release from exported pseudobulk calls without changing source bundles.

Copies only result tables and small metadata; expression sidecars remain symlinked
inside the viewer database. Unreported folds stay missing in the fold matrices.
"""
from pathlib import Path
import json
import sqlite3
import argparse
import pandas as pd


def stage(release, evidence, edge_index):
    release=Path(release);root=release.parent;evidence=Path(evidence)
    spec=json.loads(release.read_text());manifest=json.loads((evidence/'manifest.json').read_text())
    catalog=dict(spec);catalog['release']='copd-real-pseudobulk-integrated';catalog['datasets']=[];catalog['assets_root']='assets_integrated'
    assets=root/catalog['assets_root'];assets.mkdir(exist_ok=True)
    # Preserve Discover's own static per-state model scores, independent of each contrast.
    edges=json.loads(Path(edge_index).read_text())
    states=set(manifest['states']);edges={s:v for s,v in edges.items() if s in states}
    (evidence/'grn_edges_by_state.json').write_text(json.dumps(edges,separators=(',',':')))
    manifest['edge_scores']='grn_edges_by_state.json';manifest['source_edge_scores']=str(edge_index)
    (evidence/'manifest.json').write_text(json.dumps(manifest,indent=2))
    with sqlite3.connect(f'file:{evidence/manifest["differential_database"]}?mode=ro',uri=True) as con:
        for item in spec['datasets']:
            src=root/item['bundle_dir'];prefix=item['prefix'];dst=root/'bundles_integrated'/src.name;dst.mkdir(parents=True,exist_ok=True)
            for f in src.iterdir():
                if f.name in (prefix+'_metadata.json',prefix+'_deg_manifest.json',prefix+'_deg'):continue
                target=dst/f.name
                if not target.exists():target.symlink_to(f.resolve(),target_is_directory=f.is_dir())
            meta=json.loads((src/(prefix+'_metadata.json')).read_text())
            old=json.loads((src/(prefix+'_deg_manifest.json')).read_text())
            asset=json.loads((root/spec['assets_root']/(prefix+'_assets.json')).read_text())
            asset['bundle_dir']=str(dst);asset['integrated_root']=str(evidence)
            meta['scalable_viewer']['integrated_root']=str(evidence)
            comparisons=[];newassets={};template={}
            for entry in old['comparisons']:
                if entry.get('kind')!='per_cell_state':continue
                template.setdefault(entry['comparison'],entry)
            for comp,record in manifest['comparisons'].items():
                for mod,cov in record['modalities'].items():
                    if mod not in ('rna','adt','lipid','grn','grn_tf') or not cov['rows']:continue
                    frame=pd.read_sql_query('SELECT * FROM differential WHERE comparison=? AND modality=?',con,params=(comp,mod))
                    frame=frame.rename(columns={'case_mean':'case_mean_expr','control_mean':'control_mean_expr'})
                    frame['case_label']=record['case_label'];frame['control_label']=record['control_label'];frame['sig_metric']='pval'
                    base=template.get(comp,{'contrast':comp,'kind':'per_cell_state','comparison':comp})
                    cid=('' if mod=='rna' else mod+'::')+base.get('contrast',comp)+'::'+comp+'::per_cell_state'
                    rel=mod+'/'+comp+'.tsv';table=dst/(prefix+'_deg')/rel;table.parent.mkdir(exist_ok=True,parents=True)
                    frame.to_csv(table,sep='\t',index=False)
                    entry=dict(base,id=cid,modality=mod,file=rel,path=str(table),source=str(evidence/manifest['differential_database']),n_rows=len(frame),columns=list(frame),replicate_unit='pseudobulk',analysis_version='LungMAP v8',tested_populations=cov['states'])
                    comparisons.append(entry)
                    fold=evidence/'folds'/mod/(comp+'.tsv');fold.parent.mkdir(exist_ok=True,parents=True)
                    if not fold.exists():frame.pivot(index='gene',columns='population',values='log2fc').reindex(columns=cov['states']).to_csv(fold,sep='\t')
                    newassets[cid]={'modality':mod,'fold_matrix_tsv':str(fold),'networks':[],
                                    'provenance':'Existing LungMAP pseudobulk calls; unreported folds are missing, not zero.'}
            # Keep unrelated released modalities (for example cell communication).
            for entry in old['comparisons']:
                if entry.get('modality','rna') in ('rna','adt','lipid','grn','grn_tf'):continue
                retained=dict(entry)
                table=dst/(prefix+'_deg')/entry['file'];table.parent.mkdir(parents=True,exist_ok=True)
                if not table.exists():table.symlink_to(Path(entry['path']).resolve())
                retained['path']=str(table);comparisons.append(retained)
                if entry['id'] in asset['differential']:newassets[entry['id']]=asset['differential'][entry['id']]
            old['comparisons']=comparisons;meta['scalable_viewer']['deg']=old
            (dst/(prefix+'_metadata.json')).write_text(json.dumps(meta,indent=2))
            (dst/(prefix+'_deg_manifest.json')).write_text(json.dumps(old,indent=2))
            asset['differential']=newassets
            (assets/(prefix+'_assets.json')).write_text(json.dumps(asset,indent=2))
            catalog['datasets'].append(dict(item,bundle_dir=str(dst.relative_to(root))))
            print(prefix,len(comparisons),'pseudobulk comparisons staged',flush=True)
    catalog['publish']['tiers']+=['bundles_integrated','assets_integrated','integrated_pseudobulk']
    catalog['publish']['example']='Publish the release with all tiers listed above; preserve or dereference expression symlinks.'
    path=root/'viewer_release_integrated.json';path.write_text(json.dumps(catalog,indent=2));return path

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('release','evidence','edge-index'):p.add_argument('--'+name,required=True)
    print(stage(**vars(p.parse_args())))
