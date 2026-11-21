"""Manifest-driven clonal prediction and publication figures.

Defaults render saved validated fits, compact source heatmaps, and alluvials. Optional
ring and UMAP sections call their native entry points; AltAnalyze2 is the default UMAP
renderer. Re-fitting is explicit and should be submitted through the site scheduler.
All paths are relative to the manifest directory unless absolute. No shell evaluation.
"""
from pathlib import Path
import argparse
import json
import subprocess
import sys
from ..bam.clonal_inference.integrity import require, sha256


def resolve(base, value):
    p=Path(value)
    return p if p.is_absolute() else (base/p).resolve()


def run_manifest(manifest, output, patient=None, refit=False, section='all', include_unassigned=False):
    path=Path(manifest).resolve();cfg=json.loads(path.read_text());base=path.parent
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=True)
    for rel, digest in cfg.get('sha256',{}).items():
        require(sha256(resolve(base,rel))==digest,'Source changed: '+rel)
    report={'manifest':str(path),'manifest_sha256':sha256(path),'refit':refit,'patients':{},'section':section}
    patients=cfg['patients']
    require(patient is None or patient in patients,'Unknown patient: '+str(patient))
    if section in ['all','clones']:
        from ..bam.clonal_inference.detection_plotting import render
        for name,p in patients.items():
            if patient is not None and name!=patient:continue
            inputs=resolve(base,p['inputs']);result=resolve(base,p['result'])
            require(not (refit and p.get('reporting_workflow') == 'final-revised-1'), 'final-revised-1 renders saved lineage reports; refitting is a separate workflow')
            if refit:
                result=out/name/'fit'
                command=[sys.executable,'-m','altanalyze3.components.bam.clonal_inference','--matrix',str(inputs/'calls.npz'),'--cells',str(inputs/'cells.tsv'),'--events',str(inputs/'events.tsv'),'--contract',str(resolve(base,p['fit_contract'])),'--detector',str(inputs/'detector.npz'),'--sex',p['sex'],'--output',str(result)]
                subprocess.run(command,check=True)
            if p.get('reporting_workflow') == 'final-revised-1':
                require(not refit, 'final-revised-1 renders saved lineage reports; refitting is a separate workflow')
                from .clone_lineage_figures import render as render_lineage
                config=resolve(base,p['alluvial_config']) if p.get('alluvial_config') else None
                report['patients'][name]=render_lineage(inputs,result,out/name,config,include_unassigned=include_unassigned,review_panel=p.get('review_panel',True))
                continue
            kwargs={'include_unassigned':include_unassigned}
            if p.get('sensitivity_root'):kwargs['sensitivity_root']=resolve(base,p['sensitivity_root'])
            if p.get('alluvial_config'):kwargs['alluvial_config']=resolve(base,p['alluvial_config'])
            report['patients'][name]=render(result,inputs,out/name,**kwargs)
    if section in ['all','rings'] and cfg.get('rings'):
        from .variant_evidence_workflow import render
        args=dict(cfg['rings']);args['root']=resolve(base,args['root'])
        report['rings']=render(output=out/'evidence_rings',**args)
        if cfg.get('ring_highlight'):
            report['ring_highlight']=render(output=out/'evidence_rings_highlight',**dict(args,label_items=cfg['ring_highlight']))
    if section in ['all','umap'] and cfg.get('umap') and patient in [None,cfg['umap'].get('patient')]:
        from .clone_umap_mapping import prepare
        from .clone_umap_native import render as native
        from .clone_umap_panels import render as panels
        u=cfg['umap'];mapped=out/'umap_mapping'
        report['umap_mapping']=prepare(resolve(base,u['coordinates']),resolve(base,u['figures']),resolve(base,u['config']),mapped)
        report['umap_native']=native(mapped,resolve(base,u['expression']),resolve(base,u['coordinates']),resolve(base,u['altanalyze2']),u['python2'],out/'umap_native')
        report['umap_original_annotations']={}
        for name,groups in u.get('annotation_panels',{}).items():
            report['umap_original_annotations'][name]=native(mapped,resolve(base,u['expression']),resolve(base,u['coordinates']),resolve(base,u['altanalyze2']),u['python2'],out/('umap_native_'+name),resolve(base,groups))
        report['umap_panels']=panels(mapped,out/'umap_panels',resolve(base,u['config']))
    (out/('workflow_validation_'+section+('_'+patient if patient else '')+'.json')).write_text(json.dumps(report,indent=2))
    return report


def add_arguments(p):
    p.add_argument('--manifest',required=True);p.add_argument('--output',required=True);p.add_argument('--patient')
    p.add_argument('--section',choices=['all','clones','rings','umap'],default='all')
    p.add_argument('--include-unassigned',action='store_true',default=False,help='Include U in publication reporting')
    p.add_argument('--refit',action='store_true',help='Explicitly rerun the unchanged default inference; use bsub for real datasets')
def run(a):return run_manifest(a.manifest,a.output,a.patient,a.refit,a.section,a.include_unassigned)
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);add_arguments(p);return run(p.parse_args(argv))
if __name__=='__main__':main()
