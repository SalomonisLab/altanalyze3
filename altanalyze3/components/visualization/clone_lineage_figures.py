"""Publication adapters for saved minimum-supported-lineage reports; never refits.

U cells are excluded from reporting by default. WT remains; percentages are
renormalized over assigned cells. Complete original assignments are read-only.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .clone_source_heatmap import render as heatmap
from .clone_alluvial import render as alluvial
from ..bam.clonal_inference.detection_fishplot import clinical_fish_figure
from ..bam.clonal_inference.plotting import specimen_order


def render(inputs, figures, output, config=None, include_unassigned=False, review_panel=True):
    I,F,O=map(Path,[inputs,figures,output]);O.mkdir(parents=True,exist_ok=True)
    a=pd.read_csv(F/'cell_assignments.tsv.gz',sep='\t',dtype={'cell_id':str,'sample':str})
    c=pd.read_csv(I/'cells.tsv',sep='\t',dtype={'cell_id':str,'sample':str})
    if not a.cell_id.is_unique or not c.cell_id.is_unique or set(a.cell_id)!=set(c.cell_id):
        raise ValueError('Complete input and saved assignment identities must match')
    if not a.set_index('cell_id')['sample'].sort_index().equals(c.set_index('cell_id')['sample'].sort_index()):
        raise ValueError('Cell/sample identity mismatch')
    k=pd.read_csv(F/'clone_key.tsv',sep='\t');f=pd.read_csv(F/'clone_frequencies.tsv',sep='\t')
    observed=a.groupby(['sample','node']).size()
    for r in f.itertuples():
        if observed.get((r.sample,r.node),0)!=r.cells:raise ValueError('Saved frequency counts disagree')
    if int(f.cells.sum())!=len(a):raise ValueError('Saved frequency roster incomplete')
    if not set(a.node)<=set(k.node):raise ValueError('Unknown clone node')
    patient=str(c.patient.iloc[0]);meta=specimen_order(c)
    clinical=[s for s in meta.index if not bool(meta.loc[s,'is_xenograft'])]
    grafts=[s for s in meta.index if bool(meta.loc[s,'is_xenograft'])];samples=clinical+grafts
    d=a if include_unassigned else a.loc[a.node.ne(-1)].copy()
    key=k if include_unassigned else k.loc[k.node.ne(-1)].copy()
    denom=d.groupby('sample').size().reindex(samples,fill_value=0)
    if (denom==0).any():raise ValueError('No assigned cells remain in a specimen')
    counts=pd.crosstab(d.node,d['sample']).reindex(index=key.node,columns=samples,fill_value=0)
    freq=counts.div(denom,axis=1)*100
    if not np.allclose(freq.sum(),100):raise ValueError('Mass conservation failed')
    report=heatmap(I,F,O/'heatmap',patient,include_unassigned=include_unassigned)
    names=dict(zip(k.node,k.clone));colors=dict(zip(k.node,k.color));comp=dict(zip(k.node,k.mutation_combination))
    m=np.load(F/'model.npz')
    fig=clinical_fish_figure(patient,freq.to_numpy(),key.node.tolist(),m['parents'],colors,names,comp,meta,clinical,grafts,d,include_unassigned=include_unassigned)
    fig._suptitle.set_text(patient+' | Minimum supported lineage frequencies')
    legend=('One patient; direct and imputed RNA evidence. Minimum supported lineages are not proven exclusive genotypes. '
            +('All cells included. ' if include_unassigned else 'U excluded; percentages and n use assigned cells, including WT. ')
            +'Prediagnosis is schematic. Grafts are separate. No hypothesis test or confidence intervals.')
    fig.text(.035,.012,legend,fontsize=7)
    for ext in ['svg','png','pdf']:fig.savefig(O/('fishplot.'+ext),dpi=180)
    plt.close(fig)
    result={'include_unassigned':include_unassigned,'heatmap':report,'alluvial':alluvial(F,O/'alluvial',config,include_unassigned=include_unassigned),'inference_changed':False}
    rows=[]
    for s in samples:
        for node in key.node:rows.append({'sample':s,'node':int(node),'clone':names[node],'cells':int(counts.loc[node,s]),'denominator_cells':int(denom[s]),'percent':float(freq.loc[node,s])})
    pd.DataFrame(rows).to_csv(O/'displayed_clone_frequencies.tsv',sep='\t',index=False)
    accounting=pd.DataFrame({'input_cells':a.groupby('sample').size(),'plotted_cells':denom})
    accounting['excluded_unassigned_cells']=accounting.input_cells-accounting.plotted_cells
    accounting.to_csv(O/'displayed_cell_accounting.tsv',sep='\t')
    a.loc[~a.cell_id.isin(d.cell_id)].to_csv(O/'excluded_unassigned_cells.tsv.gz',sep='\t',index=False)
    import shutil
    for stem in ['clonal_tree.svg','clonal_tree.png','clonal_tree.pdf']:
        if (F/stem).exists() and (F/stem).resolve() != (O/stem).resolve():shutil.copy2(F/stem,O/stem)
    key.to_csv(O/'clone_key.tsv',sep='\t',index=False)
    (O/'display_reporting.json').write_text(json.dumps(result,indent=2))
    if review_panel:
        from .clone_review_panel import compose
        result['review_panel']=compose(O,O,patient)
    return result
