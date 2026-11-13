"""Complete source-read event matrices, with unresolved cells retained explicitly."""
from pathlib import Path
import sqlite3
import json
import collections
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
from .manifest import atomic_json,write_tsv,receipt,verify_receipt


def complete_cell_roster(observed, annotations, sources):
    """Retain the independent annotation roster plus every observed raw barcode."""
    from altanalyze3.components.bam.circRNA_cohort_pipeline import rc,barcode_state
    observed=set(observed);cells=set(observed);orientations={};annotation_maps={}
    for sid,source in sources.items():
        series=annotations.get(sid)
        if series is None:raise ValueError('Source missing from barcode annotation: '+sid)
        if series.index.has_duplicates:raise ValueError('Duplicate barcode annotation identity: '+sid)
        annotation_maps[sid]=series.to_dict()
        orientations[sid]='reverse' if str(source['reverse']).lower() in ('true','1','yes') else 'forward'
        cells.update((sid,rc(str(b)) if orientations[sid]=='reverse' else str(b)) for b in series.index)
    cells=sorted(cells);rows=[]
    for sid,barcode in cells:
        source=sources[sid];orientation=orientations[sid]
        state=barcode_state(barcode,annotation_maps[sid],orientation=orientation) if barcode else '__MISSING_BARCODE__'
        rows.append(dict(source_id=sid,barcode=barcode,normalized_barcode=rc(barcode) if orientation=='reverse' else barcode,
                         uid=source['uid'],library=source['library'],group=source['group'],state=state,
                         observed_in_bam=(sid,barcode) in observed))
    return cells,rows


def run(config):
    from altanalyze3.components.long_read.isoform_matrix import import_barcode_clusters
    root=Path(config['output']);out=root/'integrated/events';out.mkdir(parents=True,exist_ok=True)
    done=out/'COMPLETE.json'
    if done.exists():return verify_receipt(done,config)
    verify_receipt(root/'reconciliation/COMPLETE.json',config)
    db=sqlite3.connect(f'file:{root / "reconciliation/evidence.sqlite"}?mode=ro',uri=True)
    annotations=import_barcode_clusters([config['barcode_clusters']]);sources={s['source_id']:s for s in config['sources']}
    cells,cell_rows=complete_cell_roster(db.execute('SELECT DISTINCT source_id,barcode FROM reads'),annotations,sources)
    cell_index={key:i for i,key in enumerate(cells)}
    obs=pd.DataFrame(cell_rows,index=[json.dumps(k,separators=(',',':')) for k in cells])
    # Status is a separate matrix layer. A source query is counted once per event,
    # even when the raw table describes repeated crossings. Layers can overlap.
    summary={}
    modalities=['circrna','trans-splicing','oncofusion','ptd','itd']
    modalities += [r[0] for r in db.execute('SELECT DISTINCT modality FROM events ORDER BY modality') if r[0] not in modalities]
    for modality in modalities:
        branch='oncofusion' if modality in ('ptd','itd') else ('trans-splicing' if modality=='cis_edge_in_noncanonical_path' else modality)
        requested=branch in config['modalities']
        events=list(db.execute('SELECT event_id,geometry_json FROM events WHERE modality=? ORDER BY event_id',(modality,)))
        event_index={e[0]:j for j,e in enumerate(events)}
        var=pd.DataFrame({'geometry_json':[e[1] for e in events]},index=[e[0] for e in events])
        def matrix(status=None):
            condition='' if status is None else ' AND s.status=?'
            params=[modality] if status is None else [modality,status]
            rows=list(db.execute('''SELECT s.source_id,r.barcode,s.event_id,COUNT(DISTINCT s.query_name)
                FROM support s JOIN reads r USING(source_id,query_name) JOIN events e USING(event_id)
                WHERE e.modality=?'''+condition+''' GROUP BY s.source_id,r.barcode,s.event_id''',params))
            return sparse.csr_matrix(([r[3] for r in rows],([cell_index[r[0],r[1]] for r in rows],[event_index[r[2]] for r in rows])),shape=(len(cells),len(events)),dtype=np.int64)
        a=ad.AnnData(matrix(),obs=obs.copy(),var=var)
        statuses=[r[0] for r in db.execute('SELECT DISTINCT status FROM support JOIN events USING(event_id) WHERE modality=? ORDER BY status',(modality,))]
        for status in statuses:a.layers['status:'+status]=matrix(status)
        a.uns.update(count_unit='distinct source_id + query_name within event; not independent original molecules',
                     analysis_status='COMPLETE' if requested else 'NOT_REQUESTED_NOT_ASSESSED',
                     feature_universe='all observed event geometries including QC and unresolved evidence',
                     cross_modality_totals='overlap intentionally; do not sum as unique reads',
                     sequence_validation='not inferred from discovery status',barcode_universe='complete declared annotated cells plus all observed source/barcode identities, including unresolved')
        a.write_h5ad(out/(modality+'.cell_counts.h5ad'),compression='gzip')
        units=sorted({(r['uid'],r['state']) for r in cell_rows});ui={key:i for i,key in enumerate(units)}
        group=sparse.csr_matrix((np.ones(len(cells),dtype=np.int64),([ui[r['uid'],r['state']] for r in cell_rows],np.arange(len(cells)))),shape=(len(units),len(cells)))
        bulk=ad.AnnData(group@a.X,obs=pd.DataFrame({'uid':[u[0] for u in units],'state':[u[1] for u in units]},index=[json.dumps(u,separators=(',',':')) for u in units]),var=var.copy(),uns=dict(a.uns))
        for status in statuses:bulk.layers['status:'+status]=group@a.layers['status:'+status]
        if int(bulk.X.sum())!=int(a.X.sum()):raise ValueError('Event pseudobulk count conservation failed')
        bulk.write_h5ad(out/(modality+'.pseudobulk_counts.h5ad'),compression='gzip')
        summary[modality]=dict(features=len(events),source_cells=len(cells),uid_state_units=len(units),event_support=int(a.X.sum()),statuses=statuses,
                               analysis_status='COMPLETE' if requested else 'NOT_REQUESTED_NOT_ASSESSED',
                               annotated_cells_without_observed_reads=sum(not r['observed_in_bam'] for r in cell_rows))
    write_tsv(out/'cell_mapping.tsv',list(cell_rows[0]) if cell_rows else ['source_id','barcode'],cell_rows)
    atomic_json(out/'SUMMARY.json',summary);db.close()
    return receipt(done,config,{},[p for p in out.iterdir() if p.is_file()],status='COMPLETE',summary=summary)
