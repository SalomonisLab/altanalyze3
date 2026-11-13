"""Native cross-chromosome catalog, ordered-segment translation and quantification.

The validated stage-3/4 profile is WTA, minimum total 3, structure-coordinate cache
 disabled for multi-locus geometry. Full raw matrices and the complete source roster
remain available when the native catalog/annotated-cell view has no contribution.
"""
import argparse
import gzip
import json
import os
from pathlib import Path
from .manifest import load_config,write_tsv,atomic_json,receipt,verify_receipt


def run(config):
    os.environ['ISOFORM_STRUCT_COORDS']='0'
    from altanalyze3.components.long_read import gff_process as G,isoform_automate as A,isoform_matrix as I
    from .standard import working_directory
    import anndata as ad
    root=Path(config['output'])/'trans_splicing';done=root/'CATALOG_COMPLETE.json'
    if done.exists():return verify_receipt(done,config)
    rows=[];status=[];products=[]
    for source in config['sources']:
        sid=source['source_id'];folder=root/'samples'/sid
        verify_receipt(folder/'COMPLETE.json',config)
        matrix=folder/(sid+'.h5ad');gff=folder/(sid+'.gff.gz')
        if not matrix.exists():
            status.append(dict(source_id=sid,status='no_kept_cross_chromosome_molecules',raw_molecules=0,native_quant_status='not_applicable_empty'))
            continue
        G.consolidateLongReadGFFs(str(gff),config['exon_model'],mode='collapse')
        ta=folder/'gff-output/transcript_associations.txt'
        if not ta.exists():ta=Path(str(ta)+'.gz')
        opener=gzip.open if ta.suffix=='.gz' else open
        with opener(ta,'rt') as handle:annotated={f'{p[0]}:{p[3]}' for p in (line.rstrip('\n').split('\t') for line in handle if line.strip())}
        a=ad.read_h5ad(matrix);named=[v for v in a.var_names if not str(v).startswith('UNK')]
        if not set(named)<=annotated:raise ValueError('Named cross-chromosome source lacks transcript association: '+sid)
        status.append(dict(source_id=sid,status='raw_and_annotation_complete',raw_molecules=int(a.X.sum()),native_quant_status='pending'))
        rows.append(dict(uid=sid,library=sid,gff=str(gff),matrix=str(matrix),reverse=source['reverse'],groups=source['group']))
        products.extend(p for p in (folder/'gff-output').rglob('*') if p.is_file())
    metadata=root/'metadata.tsv';write_tsv(metadata,['uid','library','gff','matrix','reverse','groups'],rows)
    catalog_rows=0
    if rows:
        with working_directory(root):
            A.build_isoform_catalog(str(metadata),config['exon_model'],config['reference_gff'],config['genome'],
                                    collapse_method='wta',force_recollapse=True,min_total=3)
            catalog=root/'gff-output/FINAL_isoform_catalog.tsv'
            if not catalog.exists():catalog=Path(str(catalog)+'.gz')
            opener=gzip.open if catalog.suffix=='.gz' else open
            with opener(catalog,'rt') as handle:catalog_rows=max(0,sum(1 for _ in handle)-1)
            if catalog_rows:
                bars=I.import_barcode_clusters([config['barcode_clusters']]);samples=A.import_metadata(str(metadata),include_hashed_samples=True)
                statuses={r['source_id']:r for r in status}
                for sid,libs in samples.items():
                    try:
                        A.export_sample_isoform(libs,str(root/'gff-output'),bars,uid=sid,collapse_method='wta',compute_cpm=True)
                        statuses[sid]['native_quant_status']='complete'
                    except ValueError as error:
                        if 'barcodes matched cluster annotations' not in str(error):raise
                        statuses[sid]['native_quant_status']='no_annotated_cells_raw_matrix_retained'
                        statuses[sid]['detail']=str(error)
            else:
                for row in status:
                    if row['native_quant_status']=='pending':row['native_quant_status']='native_min_total_catalog_empty_raw_matrix_retained'
        atomic_json(root/'STAGE3_COMPLETE.json',dict(catalog_rows=catalog_rows,min_total=3,libraries_in_collapse=len(rows),libraries_without_cross_chromosome_molecules=[r['source_id'] for r in status if r['raw_molecules']==0]))
        atomic_json(root/'STAGE4_COMPLETE.json',dict(per_library={r['source_id']:r for r in status},libraries=len(status)))
        products.extend(p for p in (root/'gff-output').rglob('*') if p.is_file())
        products.extend([root/'STAGE3_COMPLETE.json',root/'STAGE4_COMPLETE.json'])
    atomic_json(root/'CATALOG_SOURCE_STATUS.json',dict(sources=status,catalog_rows=catalog_rows,
        parameters=dict(collapse_method='wta',min_total=3,structure_coordinate_cache=False),
        scope='native catalog threshold and annotated-cell view; every raw source/cell/query remains in evidence'))
    products += [metadata,root/'CATALOG_SOURCE_STATUS.json']
    products += [p for p in (root/'samples').rglob('*') if p.is_file() and ('-isoform' in p.name or p.name=='transcript_associations.txt.gz')]
    return receipt(done,config,{},products,status='COMPLETE',sources=status,catalog_rows=catalog_rows)


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);a=p.parse_args();run(load_config(a.config))
if __name__=='__main__':main()
