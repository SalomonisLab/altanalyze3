"""Editable vector clone composition, frequency, and source evidence figures.

No ancestry is invented. Grafts appear beside their matched clinical stage and
remain separate specimens/denominators. SVG text stays editable (no raster).
"""
from pathlib import Path
import argparse,json,textwrap,xml.etree.ElementTree as ET
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle,Patch
from .integrity import require,sha256


def specimen_order(cells):
    meta=cells.drop_duplicates('sample').set_index('sample')
    if 'stage_order' in meta:
        meta=meta.sort_values(['stage_order','is_xenograft'] if 'is_xenograft' in meta else ['stage_order'],kind='stable')
    return meta


def render(result,matrix,output=None,coverage_note='',max_main=12):
    manifest_path=Path(result)/'run_manifest.json'
    if manifest_path.exists() and json.loads(manifest_path.read_text()).get('algorithm',{}).get('method')=='detection-tree':
        from .detection_plotting import render as detection_render
        return detection_render(result,Path(matrix).parent,output or Path(result)/'figures')
    if manifest_path.exists() and json.loads(manifest_path.read_text()).get('arguments',{}).get('method')=='supported-combinations':
        from .supported_plotting import render as supported_render
        return supported_render(result,output,max_main)
    result=Path(result);output=Path(output) if output else result/'figures';output.mkdir(parents=True,exist_ok=True)
    require((result/'COMPLETE.json').exists(),'Fit is not complete')
    plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42,'font.family':'Arial','font.size':9,'axes.linewidth':.7,'svg.hashsalt':'altanalyze3-clones'})
    cells=pd.read_csv(result/'cell_roster.tsv',sep='\t',dtype={'cell_id':str,'sample':str});events=pd.read_csv(result/'event_roster.tsv',sep='\t')
    assignment=pd.read_csv(result/'cell_assignments.tsv.gz',sep='\t');F=pd.read_csv(result/'clone_frequencies.tsv',sep='\t');G=pd.read_csv(result/'clone_genotypes.tsv',sep='\t')
    selection=json.loads((result/'selection.json').read_text());meta=specimen_order(cells);samples=meta.index.tolist()
    manifest=json.loads((result/'run_manifest.json').read_text())
    require(sha256(matrix)==manifest['contract']['input_sha256']['matrix'],'Plot matrix hash differs from fitted input')
    z=np.load(matrix);X=z['X'];require(list(z['cell_ids'])==cells.cell_id.tolist()==assignment.cell_id.tolist(),'Plot cell identity/order mismatch');require(list(z['features'])==events.feature.tolist(),'Plot event identity mismatch')
    require(np.allclose(F.groupby('sample').percent_all.sum(),100),'Frequency denominator failure')
    actual=assignment.groupby(['sample','clone']).size()
    totals=cells.groupby('sample').size()
    for r in F.itertuples():
        require(r.cells==actual.get((r.sample,r.clone),0) and r.all_cells==totals[r.sample]
                and np.isclose(r.percent_all,100*r.cells/r.all_cells),'Frequency table disagrees with full assignment roster')
    patient=str(cells.patient.iloc[0]) if 'patient' in cells else 'Patient'
    frequencies=F.pivot(index='clone',columns='sample',values='percent_all').reindex(columns=samples)
    count=F.groupby('clone').cells.sum();order=sorted([k for k in frequencies.index if k!='NE' and count[k]>0],key=lambda k:int(k[2:]))
    names=F.drop_duplicates('clone').set_index('clone').composition.to_dict();names['Other']='Other predicted combinations';names['NE']='No evidence / unresolved'
    counts=cells['sample'].value_counts().reindex(samples)
    sample_labels=[str(meta.loc[u,'stage'])+'\nn='+format(int(counts[u]),',') for u in samples]
    informative=(X>=0).any(1);n=int(len(cells));J=len(events)
    event_labels=events.event_label.astype(str).tolist() if 'event_label' in events else events.feature.tolist()
    genotypes=G.pivot(index='clone',columns='feature',values='state').reindex(columns=events.feature)
    seed_agreement=assignment.groupby('clone').full_genotype_start_agreement.mean()
    palette=[plt.get_cmap('tab20').colors[i%20] for i in range(len(order))];color={k:palette[i] for i,k in enumerate(order)};color['Other']='#b4b4b4';color['NE']='#eeeeee'
    settings=json.loads((result/'run_manifest.json').read_text())['arguments']
    caption=(f'One patient ({patient}); {len(samples)} separate specimens, {n:,} cells, {J} retained events. Shared binary-genotype mixture with learned event errors capped at FP {100*settings["fp_cap"]:g}%, FN {100*settings["fn_cap"]:g}%; '
             f'component count selected by nested held-out call prediction; representative fit selected by training objective across {settings["starts"]} starts. '
             'Percentages use every profiled cell in each specimen, including predicted panel-WT and no-evidence cells. '
             'Direct and accepted imputed calls enter once per event; missing calls provide no evidence. Panel-WT is a model prediction, not proof of genome-wide normality. '
             'No significance test or biological-replicate confidence interval; starts measure optimization sensitivity, not biological replication. '
             'Clone genotypes and frequencies remain uncertain where observations are sparse or imputation is correlated. '+coverage_note)
    outputs=[]
    def save(fig,stem,legend):
        fig.canvas.draw();renderer=fig.canvas.get_renderer()
        for txt in fig.findobj(match=matplotlib.text.Text):
            if txt.get_visible() and txt.get_text():
                box=txt.get_window_extent(renderer)
                require(box.x0>=-1 and box.y0>=-1 and box.x1<=fig.bbox.width+1 and box.y1<=fig.bbox.height+1,
                        f'Figure text outside canvas in {stem}: {txt.get_text()[:100]}')
        for ext in ['svg','pdf','png']:
            path=output/(stem+'.'+ext);fig.savefig(path,dpi=180,metadata={'Creator':'AltAnalyze3 clonal inference'} if ext in ['svg','pdf'] else None)
            outputs.append(path.name)
            if ext=='svg':
                root=ET.parse(path).getroot();require(not any(e.tag.endswith('}image') for e in root.iter()),'SVG contains raster pixels');require(any(e.tag.endswith('}text') for e in root.iter()),'SVG text not editable')
        (output/(stem+'_legend.txt')).write_text(legend+'\n');plt.close(fig)
    def frequency_plot(shown,stem,pool=False,sensitivity=False):
        remaining=[k for k in order if k not in shown]
        rows=shown+(['Other'] if pool and remaining else [])+['NE']
        V=np.array([frequencies.loc[remaining].sum().to_numpy() if k=='Other' else frequencies.loc[k].to_numpy() for k in rows])
        gm=np.array([genotypes.loc[k].to_numpy() if k not in ['Other','NE'] else np.full(J,-1) for k in rows])
        if pool or not remaining:require(np.allclose(V.sum(0),100),'Displayed totals lost cells')
        wide=max(len(names[k]) for k in rows)>100
        labs=[textwrap.fill((k+' = '+names[k]).replace('|',' | '),95 if wide else 68,break_long_words=False) if k not in ['Other','NE'] else names[k] for k in rows]
        lines=max(label.count('\n')+1 for label in labs)
        height=max(9.8,3.3+.47*len(rows),len(rows)*(lines*9.6+4)/72/.57)
        fig=plt.figure(figsize=(22 if wide else 18,height))
        a=fig.add_axes([.43,.29,.25,.57] if wide else [.31,.29,.32,.57]);b=fig.add_axes([.72,.29,.26,.57] if wide else [.68,.29,.29,.57])
        cm=ListedColormap(['#d9d9d9','#e8eff6','#d94c3b']);norm=BoundaryNorm([-1.5,-.5,.5,1.5],3)
        a.pcolormesh(np.arange(J+1)-.5,np.arange(len(rows)+1)-.5,gm,cmap=cm,norm=norm,edgecolors='white',linewidth=.5,rasterized=False)
        a.set_ylim(len(rows)-.5,-.5);a.set_xticks(range(J),event_labels,rotation=75,ha='right',fontsize=8)
        a.set_yticks(range(len(rows)),labs,fontsize=8);a.tick_params(length=0);a.set_title('Predicted mutation combination',pad=18,fontsize=12)
        b.pcolormesh(np.arange(len(samples)+1)-.5,np.arange(len(rows)+1)-.5,V,cmap='Blues',vmin=0,vmax=100,edgecolors='white',linewidth=.5,rasterized=False)
        b.set_ylim(len(rows)-.5,-.5);b.set_yticks([]);b.set_xticks(range(len(samples)),[textwrap.fill(t.split('\n')[0],14)+'\n'+t.split('\n')[1] for t in sample_labels],fontsize=8);b.xaxis.tick_top();b.tick_params(length=0)
        low=F.pivot(index='clone',columns='sample',values='start_min_percent_all').reindex(columns=samples)
        high=F.pivot(index='clone',columns='sample',values='start_max_percent_all').reindex(columns=samples)
        for i,j in np.ndindex(V.shape):
            txt=f'{V[i,j]:.2f}'
            if sensitivity and rows[i]!='Other':txt+=f'\n[{low.loc[rows[i]].iloc[j]:.1f}, {high.loc[rows[i]].iloc[j]:.1f}]'
            if sensitivity and rows[i]=='Other':txt+='\n[—]'
            b.text(j,i,txt,ha='center',va='center',fontsize=7.5 if sensitivity else 9,color='white' if V[i,j]>55 else '#173047')
        local_unknown=G.pivot(index='clone',columns='feature',values='genotype_unobserved_within_assigned_cells').reindex(columns=events.feature)
        for i,k in enumerate(rows):
            if k not in ['Other','NE']:
                for j in range(J):
                    if bool(local_unknown.loc[k].iloc[j]):a.text(j,i,'○',ha='center',va='center',fontsize=10,color='#222')
        b.set_title('Frequency % [min, max across starts]' if sensitivity else 'Frequency (% of all specimen cells)',pad=48,fontsize=12)
        fig.suptitle(f'{patient} | mutation-defined clone composition across specimens',x=.025,y=.975,ha='left',fontsize=17)
        fig.text(.025,.925,'MISTRG specimens are displayed beside their matched patient stage; they are not additional chronological patient samples.' if ('is_xenograft' in meta and meta.is_xenograft.any()) else 'Columns follow the established specimen order; frequencies describe one patient, not biological replication.',fontsize=10)
        fig.legend(handles=[Patch(facecolor='#d94c3b',label='Predicted MUT'),Patch(facecolor='#e8eff6',label='Predicted WT'),Patch(facecolor='#d9d9d9',label='Unknown / pooled / no-evidence')],loc='lower left',bbox_to_anchor=(.025,.16),ncol=3,frameon=False,fontsize=9)
        legend=('Rows are complete predicted mutation combinations; columns at left are distinct events and columns at right are specimens. '
                + (f'The {max_main} profiles with highest maximum specimen frequency are shown; every remaining profile is pooled only in the display as Other and retained individually in the full SVG and tables. ' if pool else 'All nonzero assigned profiles are shown individually. ')
                +'White cell borders separate entries; blue intensity and printed values show percentages. Red means predicted MUT, pale blue predicted WT, gray unresolved or pooled. '
                +('Brackets show the minimum and maximum frequency of the exact genotype across starts, not confidence intervals. Other contains different profiles and its range is not combined. ' if sensitivity else '')
                +'A hollow circle marks an event without observations in the assigned cells; its state is model extrapolation. Full-genotype start agreement and per-event MUT/WT/unknown support are provided in companion tables. '+caption)
        fig.text(.025,.025,textwrap.fill(legend,225),va='bottom',fontsize=8);save(fig,stem,legend)
    frequency_plot(order[:max_main],patient+'_clone_composition_frequency',pool=True)
    frequency_plot(order,patient+'_all_clones',pool=False)
    frequency_plot(order[:max_main],patient+'_optimization_sensitivity',pool=True,sensitivity=True)
    # Stacked specimen bars keep grafts next to the matching stage without
    # claiming a lineage edge, a cell transition, or an elapsed-time scale.
    shown=order[:max_main];remaining=order[max_main:];rows=shown+(['Other'] if remaining else [])+['NE']
    fig,ax=plt.subplots(figsize=(15,10));fig.subplots_adjust(left=.07,right=.49,bottom=.25,top=.87)
    bottom=np.zeros(len(samples))
    for k in rows:
        vals=frequencies.loc[remaining].sum().to_numpy() if k=='Other' else frequencies.loc[k].to_numpy()
        bars=ax.bar(np.arange(len(samples)),vals,bottom=bottom,color=color[k],width=.78,edgecolor='white',linewidth=.3,label=k+' = '+names[k])
        for j,bar in enumerate(bars):
            if 'is_xenograft' in meta and bool(meta.loc[samples[j],'is_xenograft']):bar.set_hatch('///');bar.set_edgecolor('#555');bar.set_linewidth(.25)
        bottom+=vals
    require(np.allclose(bottom,100),'Stacked plot lost denominator')
    ax.set_ylim(0,100);ax.set_ylabel('% of all specimen cells');ax.set_xticks(range(len(samples)),[textwrap.fill(t.split('\n')[0],12)+'\n'+t.split('\n')[1] for t in sample_labels],fontsize=8)
    ax.spines[['top','right']].set_visible(False);ax.legend(loc='upper left',bbox_to_anchor=(1.03,1.02),labels=[textwrap.fill((k+' = '+names[k]).replace('|',' | '),65,break_long_words=False) for k in rows],frameon=False,fontsize=8)
    fig.suptitle(f'{patient} | clone frequencies across specimens',x=.03,ha='left',fontsize=15)
    legend='Each stacked bar is one specimen and totals 100%. Colors identify the same mutation combinations across this patient; hatched bars are MISTRG. No ancestry, cell transition, or interpolation is inferred. Up to 12 dominant profiles are shown; all others are pooled only for display. '+caption
    fig.text(.025,.025,textwrap.fill(legend,185),va='bottom',fontsize=8);save(fig,patient+'_specimen_composition',legend)
    # Source-resolved evidence is drawn as run-length vector rectangles. No
    # cell subsampling and no rasterization; source conflicts are visible.
    method_names=list(z['source_method_order']) if 'source_method_order' in z.files else []
    mut=z['imputed_MUT_sources'] if method_names else None;wt=z['imputed_WT_sources'] if method_names else None
    direct=z['direct'] if 'direct' in z.files else np.where(z['T']>0,X,-1)
    rank={k:i for i,k in enumerate(order+['NE'])};collected=[]
    indices=[]
    for u in samples:
        ix=np.flatnonzero(cells['sample'].to_numpy()==u);ix=ix[np.argsort([rank.get(k,len(order)) for k in assignment.clone.iloc[ix]],kind='stable')];indices.append(ix);collected.extend(ix.tolist())
    require(sorted(collected)==list(range(len(cells))),'Source plot excluded or duplicated a cell')
    source_colors=['#ffffff','#b4d7a8','#e2664a','#e4ad38']
    for page,start in enumerate(range(0,J,6),1):
        feats=list(range(start,min(start+6,J)));rows=[];labels=[]
        for j in feats:
            predicted=np.full(n,-1,int)
            for k in order:predicted[assignment.clone.to_numpy()==k]=int(genotypes.loc[k].iloc[j])
            for tag,values in [('MODEL',predicted),('CALL',X[:,j]),('READ',direct[:,j])]:rows.append(np.where(values<0,0,1+values));labels.append(event_labels[j]+' | '+tag)
            for bit,m in enumerate(method_names):
                pos=(mut[:,j]&(1<<bit))>0;neg=(wt[:,j]&(1<<bit))>0
                rows.append(np.where(pos&neg,3,np.where(pos,2,np.where(neg,1,0))));labels.append(event_labels[j]+' | '+m)
        mat=np.stack(rows);height=max(8,3+.15*len(rows));fig=plt.figure(figsize=(18,height))
        left=.28;right=.975;width=(right-left)/len(samples);axes=[]
        for si,(u,ix) in enumerate(zip(samples,indices)):
            ax=fig.add_axes([left+si*width,.20,width-.008,.64]);patches=[];colors=[]
            for ri,row in enumerate(mat[:,ix]):
                cuts=np.r_[0,np.flatnonzero(np.diff(row)!=0)+1,len(row)]
                for aa,bb in zip(cuts[:-1],cuts[1:]):
                    value=int(row[aa])
                    if value:patches.append(Rectangle((aa,ri-.5),bb-aa,1));colors.append(source_colors[value])
            collection=PatchCollection(patches,facecolors=colors,edgecolors='none',rasterized=False);collection.set_gid('source_evidence_'+str(si));ax.add_collection(collection)
            ax.set_xlim(0,len(ix));ax.set_ylim(len(rows)-.5,-.5);ax.set_yticks(range(len(rows)),labels if si==0 else ['']*len(rows),fontsize=6.8);ax.tick_params(length=0);ax.set_xticks([]);ax.set_xlabel('All cells, ordered by clone',fontsize=8)
            ax.set_title(textwrap.fill(str(meta.loc[u,'stage']),15)+'\nn='+format(len(ix),','),fontsize=9,pad=20)
            ca=assignment.clone.to_numpy()[ix];cuts=np.r_[0,np.flatnonzero(ca[1:]!=ca[:-1])+1,len(ix)]
            for aa,bb in zip(cuts[:-1],cuts[1:]):
                ax.axvline(aa,color='#666',lw=.3)
                if (bb-aa)/len(ix)>=.08:ax.text((aa+bb)/2,-1.6,ca[aa],fontsize=6.5,ha='center')
            block=3+len(method_names)
            for y in range(block,len(rows),block):ax.axhline(y-.5,color='#333',lw=.5)
        fig.suptitle(f'{patient} | serial clone assignments and source evidence (page {page})',x=.025,y=.975,ha='left',fontsize=15)
        fig.text(.025,.925,'Green = WT; red = MUT; white = unknown / no source call; amber = conflicting imputed MUT and WT.',fontsize=10)
        legend=('Every cell is shown once per specimen, grouped by the representative clone assignment. Each variant has MODEL (inferred genotype), CALL (consolidated input), READ (direct allele), and source-specific rows. '
                'PV=passenger variant; MPV=mitochondrial passenger; LOH=loss-of-heterozygosity classifier; CNV=copy-number association; DEG/DEI=gene/isoform expression classifier. '
                'Sources support the same target event and are not counted as independent mutation features. Multiple sources can mark the same cell; amber retains opposing source predictions. '
                'No evidence is inferred from white entries. '+caption)
        fig.text(.025,.025,textwrap.fill(legend,225),va='bottom',fontsize=8);save(fig,f'{patient}_serial_evidence_{page}',legend)
    from .fishplot import render as render_fishplots
    clinical_count=int((~meta.is_xenograft).sum()) if 'is_xenograft' in meta else len(meta)
    if clinical_count>=2:
        fish_review=render_fishplots(result,output,max_main=max_main)
        outputs.extend(fish_review['files'])
    review=dict(patient=patient,cells=n,events=J,samples=samples,MISTRG_specimens=int(meta.is_xenograft.sum()) if 'is_xenograft' in meta else 0,all_cells_in_source_heatmaps=True,denominator_includes_WT_and_unresolved=True,SVG_contains_no_raster=True,SVG_text_editable=True,files=outputs,coverage_note=coverage_note)
    (output/'figure_checks.json').write_text(json.dumps(review,indent=2));return review


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--result',required=True);p.add_argument('--matrix',required=True);p.add_argument('--output');p.add_argument('--coverage-note',default='');a=p.parse_args(argv)
    render(a.result,a.matrix,a.output,a.coverage_note)
if __name__=='__main__':main()
