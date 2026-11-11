"""Run existing long-read analysis functions in isolated, versioned views."""
from __future__ import annotations
from contextlib import contextmanager
import json
import os
from pathlib import Path
import shutil
from .manifest import atomic_json, receipt, verify_receipt, write_tsv, digest

@contextmanager
def working_directory(path):
    previous=Path.cwd();os.chdir(path)
    try:yield
    finally:os.chdir(previous)


def run_standard(config, view='compatibility'):
    """Same native catalog and statistical methods; explicit original biological UIDs.

    Compatibility is generated first. The integrated view reuses the complete
    catalog; reclassified reads remain in provenance and zero features remain in
    source matrices. Conventional junctions and gene expression retain the baseline
    physical-edge/gene counting view, independent of whole-read reassignment.
    """
    from altanalyze3.components.long_read import isoform_automate as A, isoform_matrix as I
    from altanalyze3.components.long_read import comparisons as C
    root=Path(config['output']);out=root/view if view=='compatibility' else root/'integrated/standard'
    out.mkdir(parents=True,exist_ok=True);done=out/'STANDARD_COMPLETE.json'
    if done.exists():return verify_receipt(done,config)
    if not config.get('standard_parameters'):
        raise ValueError('Explicit frozen standard_parameters required for complete native analyses')
    parameters=config['standard_parameters']
    required={'collapse_method','min_total','protect_intronic_evidence','protect_novel_junctions','analyses','method','conditions'}
    if not required<=parameters.keys():raise ValueError('Incomplete native baseline method configuration')
    rows=[]
    for s in config['sources']:
        sid=s['source_id'];folder=out/'samples'/sid;folder.mkdir(parents=True,exist_ok=True)
        raw=root/'compatibility/samples'/sid
        gff=raw/(sid+'.gff.gz')
        if not gff.exists():gff=raw/(sid+'.gff')
        matrix=folder/(sid+'.h5ad')
        if not matrix.exists():raise ValueError(f'No native source matrix for {sid}; no source may silently disappear')
        if view!='compatibility':
            target=folder/gff.name
            if not target.exists():target.symlink_to(gff)
            gff=target
        rows.append(dict(uid=s['uid'],library=sid,gff=str(gff),matrix=str(matrix),reverse=s['reverse'],groups=s['group']))
    meta=out/'metadata.tsv';write_tsv(meta,list(rows[0]),rows)
    bars=I.import_barcode_clusters([config['barcode_clusters']])
    with working_directory(out):
        samples=A.import_metadata(str(meta),include_hashed_samples=True)
        if view=='compatibility':
            for uid,libs in samples.items():
                A.export_sample_junctions(libs,config['exon_model'],bars,uid=uid)
                A.export_sample_gene(libs,bars,uid=uid)
            A.combine_junctions(str(meta),[config['barcode_clusters']],config['exon_model'])
            A.build_isoform_catalog(str(meta),config['exon_model'],config['reference_gff'],config['genome'],
                collapse_method=parameters['collapse_method'],min_total=parameters['min_total'],
                protect_intronic_evidence=parameters['protect_intronic_evidence'],
                protect_novel_junctions=parameters['protect_novel_junctions'])
        else:
            legacy=root/'compatibility'
            for p in legacy.iterdir():
                if p.name=='gff-output':
                    if not (out/p.name).exists():shutil.copytree(p,out/p.name)
                elif p.name.startswith(('junction_combined_','psi_combined_','gene_combined_')):
                    target=out/p.name
                    if not target.exists():target.symlink_to(p,target_is_directory=p.is_dir())
            # Per-sample annotation/catalog inputs are read-only bindings to the exact baseline.
            for s in config['sources']:
                sid=s['source_id'];raw=legacy/'samples'/sid;dest=out/'samples'/sid
                for p in raw.iterdir():
                    if p.name=='gff-output' or '-gene' in p.name or '-junction' in p.name:
                        if not (dest/p.name).exists():(dest/p.name).symlink_to(p,target_is_directory=p.is_dir())
        if view!='compatibility':
            bind_derived_catalog(config,out,root/'compatibility')
        for uid,libs in samples.items():
            A.export_sample_isoform(libs,str(out/'gff-output'),bars,uid=uid,
                                    collapse_method=parameters['collapse_method'],compute_cpm=True)
        A.combine_isoforms(str(meta),[config['barcode_clusters']])
        if 'gene' in parameters['analyses'] and view=='compatibility':
            A.combine_genes(str(meta),[config['barcode_clusters']])
        if view!='compatibility':
            from .quantify import write_anchored_ratios
            write_anchored_ratios(out,root/'compatibility')
        order=I.return_cluster_order([config['barcode_clusters']])
        if parameters['conditions']:
            C.compute_differentials(samples,parameters['conditions'],order,config['gene_symbols'],
                                    analyses=parameters['analyses'],method=parameters['method'])
    outputs=[p for p in out.rglob('*') if p.is_file() and not p.name.endswith('.log') and p.name!='STANDARD_COMPLETE.json']
    receipt(done,config,{},outputs,status='COMPLETE',assignment_view=view,standard_parameters=parameters)


def bind_derived_catalog(config,out,legacy):
    """Bind derived matrices only after verifying the complete reconciliation receipt.

    Catalog evidence and structure assignments remain byte-for-byte inherited; a
    new binding records both the original and derived inputs and their ledger.
    """
    root=Path(config['output'])
    proof=root/'reconciliation/COMPLETE.json'
    verify_receipt(proof,config)
    path=out/'gff-output/intronic_evidence_manifest.json'
    if not path.exists():return
    manifest=json.loads(path.read_text())
    bindings=manifest['source_bindings']
    changes=[]
    for source in config['sources']:
        sid=source['source_id'];original=legacy/'samples'/sid/(sid+'.h5ad')
        matches=[key for key,b in bindings.items() if Path(b['matrix']).resolve()==original.resolve()]
        if len(matches)!=1:raise ValueError('Exact catalog source cannot be rebound: '+sid)
        key=matches[0];derived=out/'samples'/sid/(sid+'.h5ad')
        changes.append(dict(source_id=sid,catalog_source=key,original_binding=dict(bindings[key]),
                            original_matrix_sha256=digest(original),derived_matrix=str(derived),
                            derived_matrix_sha256=digest(derived)))
        bindings[key]['matrix']=str(derived.resolve())
    manifest['derived_assignment_view']=dict(name='integrated',reconciliation_sha256=digest(proof),sources=changes)
    atomic_json(path,manifest)
