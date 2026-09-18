"""Shared integrated analysis endpoints for uploaded jobs and precomputed bundles."""
import csv
import io
from fastapi import HTTPException, Query
from fastapi.responses import PlainTextResponse
from altanalyze3.components.cellHarmony import discover_integration as integration
from .integration_data import integration_data


def install(app):
    def data(job_id):
        store=app.state.job_store
        if not store.job_exists(job_id):raise HTTPException(404,'Job not found.')
        return integration_data(app,store.get_job(job_id))

    @app.get('/api/jobs/{job_id}/integrated/network')
    def network(job_id:str,cell_state:str='',contrast:str='',features:str='',source:str=Query('differential',pattern='^(differential|marker)$'),
                limit:int=Query(50,ge=1,le=2000),min_fold:float=Query(1.2,ge=1),gene_fold:float=Query(1.2,ge=1),
                max_fdr:float=Query(.05,ge=0,le=1),min_expression:float=Query(.5,ge=0),
                significance:str=Query('reported',pattern='^(reported|fdr|pval)$'),min_score:float=Query(0,ge=0)):
        try:
            return integration.network(data(job_id),contrast=contrast,cell_state=cell_state,
                features=[x.strip() for x in features.split(',') if x.strip()],source=source,
                limit=limit,min_fold=min_fold,gene_fold=gene_fold,max_fdr=max_fdr,min_expression=min_expression,significance=significance,min_score=min_score)
        except (ValueError,KeyError,FileNotFoundError) as exc:
            return {'available':False,'nodes':[],'edges':[],'note':str(exc),'status':'unavailable'}

    @app.get('/api/jobs/{job_id}/integrated/network.tsv')
    def network_tsv(job_id:str,cell_state:str='',contrast:str='',features:str='',source:str=Query('differential',pattern='^(differential|marker)$'),
                    limit:int=Query(50,ge=1,le=2000),min_fold:float=Query(1.2,ge=1),gene_fold:float=Query(1.2,ge=1),
                    max_fdr:float=Query(.05,ge=0,le=1),min_expression:float=Query(.5,ge=0),
                    significance:str=Query('reported',pattern='^(reported|fdr|pval)$'),min_score:float=Query(0,ge=0)):
        result=network(job_id,cell_state,contrast,features,source,limit,min_fold,gene_fold,max_fdr,min_expression,significance,min_score)
        out=io.StringIO();writer=csv.writer(out,delimiter='\t')
        marker = source == 'marker'
        if marker:
            writer.writerow(['TF','target','edge_score','TF_expression','TF_marker_log2fc','target_marker_log2fc','cell_state','expression_scale'])
        else:
            writer.writerow(['TF','target','edge_score','edge_log2fc','edge_significance','significance_metric','TF_expression','TF_expression_log2fc','TF_expression_FDR','TF_activity_log2fc','TF_activity_FDR'])
        nodes={n['id']:n for n in result.get('nodes',[])}
        for e in result.get('edges',[]):
            n=nodes.get(e['source'],{})
            if marker:
                writer.writerow([e['source'],e['target'],e['score'],n.get('expression'),n.get('log2fc'),nodes.get(e['target'],{}).get('log2fc'),result.get('cell_state'),result.get('expression_scale')])
            else:
                writer.writerow([e.get(k) for k in ('source','target','score','log2fc','fdr','significance')]+[n.get(k) for k in ('expression','log2fc','fdr','activity_log2fc','activity_fdr')])
        if not result.get('edges'):out.write('# '+result.get('note','No edges pass the filters.')+'\n')
        return PlainTextResponse(out.getvalue(),media_type='text/tab-separated-values',headers={'Content-Disposition':'attachment; filename="regulatory_network.tsv"'})

    @app.get('/api/jobs/{job_id}/integrated/pathways')
    def pathways(job_id:str,cell_state:str='',contrast:str='',source:str=Query('differential',pattern='^(differential|marker)$'),features:str='',modality:str=Query('lipid',pattern='^(lipid|lipids|metabolite)$'),
                 max_fdr:float=Query(.05,ge=0,le=1),significance:str=Query('reported',pattern='^(reported|fdr|pval)$')):
        return integration.pathways(data(job_id),contrast=contrast,cell_state=cell_state,modality=modality,max_fdr=max_fdr,significance=significance,source=source,features=[x for x in features.split(',') if x])

    @app.get('/api/jobs/{job_id}/integrated/pathway')
    def pathway(job_id:str,id:str,cell_state:str='',contrast:str='',modality:str=Query('lipid',pattern='^(lipid|lipids|metabolite)$'),
                max_fdr:float=Query(.05,ge=0,le=1),significance:str=Query('reported',pattern='^(reported|fdr|pval)$')):
        return integration.pathway(data(job_id),id=id,contrast=contrast,cell_state=cell_state,modality=modality,max_fdr=max_fdr,significance=significance)

    @app.get('/api/jobs/{job_id}/integrated/cross-pathways')
    def cross_pathways(job_id: str, cell_state: str, source: str=Query('marker', pattern='^(marker|differential)$'), contrast: str=''):
        from .cross_pathways import collect
        data(job_id)  # Validate the job through the shared store.
        result = collect(app, app.state.job_store.get_job(job_id), cell_state, source, contrast)
        return dict(available=bool(result['rows']), pathways=[dict(r, name=r['pathway']) for r in result['rows']],
                    note='No retained features map to the bundled pathways.' if not result['rows'] else '')

    @app.get('/api/jobs/{job_id}/integrated/cross-pathway')
    def cross_pathway(job_id: str, id: str, cell_state: str, source: str=Query('marker', pattern='^(marker|differential)$'), contrast: str=''):
        from .cross_pathways import collect, diagram
        data(job_id)
        result = collect(app, app.state.job_store.get_job(job_id), cell_state, source, contrast)
        return diagram(result, id)
