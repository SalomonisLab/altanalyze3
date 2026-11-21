"""Fish-style composition trajectories and mutation-dot clone cartoons.

Uses saved assignments only. Smooth joins describe frequency interpolation,
not ancestry; graft specimens are separate bars beside the patient course.
"""
from pathlib import Path
import argparse,json,textwrap,xml.etree.ElementTree as ET
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Patch
from .integrity import require
from .plotting import specimen_order


def smooth_composition(values,points=70):
    """Convex cosine interpolation: exact knots, nonnegative, same total mass."""
    require(values.ndim==2 and values.shape[1]>=2,'Two patient specimens required')
    xs=[];ys=[]
    for j in range(values.shape[1]-1):
        t=np.linspace(0,1,points,endpoint=j==values.shape[1]-2);a=(1-np.cos(np.pi*t))/2
        xs.append(j+t);ys.append(values[:,j,None]*(1-a)+values[:,j+1,None]*a)
    return np.concatenate(xs),np.concatenate(ys,axis=1)


def render(result,output=None,max_main=12):
    result=Path(result);output=Path(output) if output else result/'figures'
    require((result/'COMPLETE.json').exists(),'Fit incomplete')
    C=pd.read_csv(result/'cell_roster.tsv',sep='\t');A=pd.read_csv(result/'cell_assignments.tsv.gz',sep='\t')
    F=pd.read_csv(result/'clone_frequencies.tsv',sep='\t');G=pd.read_csv(result/'clone_genotypes.tsv',sep='\t');E=pd.read_csv(result/'event_roster.tsv',sep='\t')
    require(C.cell_id.tolist()==A.cell_id.tolist(),'Assignment roster differs')
    counts=A.groupby(['sample','clone']).size();totals=C.groupby('sample').size()
    for r in F.itertuples():require(r.cells==counts.get((r.sample,r.clone),0) and r.all_cells==totals[r.sample] and np.isclose(r.percent_all,100*r.cells/r.all_cells),'Frequency/count mismatch')
    require(np.allclose(F.groupby('sample').percent_all.sum(),100),'Denominator failure')
    meta=specimen_order(C)
    if 'is_xenograft' not in meta:meta['is_xenograft']=False
    patients=meta.index[~meta.is_xenograft].tolist();grafts=meta.index[meta.is_xenograft].tolist()
    require(len(patients)>=2,'No patient time-course')
    patient=str(C.patient.iloc[0]);J=len(E);V=F.pivot(index='clone',columns='sample',values='percent_all')
    names=F.drop_duplicates('clone').set_index('clone').composition.to_dict();names.update(Other='Other predicted combinations',NE='No evidence / unresolved')
    order=sorted(F.loc[(F.clone!='NE')&(F.cells>0),'clone'].unique(),key=lambda k:int(k[2:]))
    shown=order[:max_main];remaining=order[max_main:];rows=shown+(['Other'] if remaining else [])+['NE']
    colors={k:plt.get_cmap('tab20').colors[i%20] for i,k in enumerate(order)};colors.update(Other='#b4b4b4',NE='#eeeeee')
    event_colors=[plt.get_cmap('tab20').colors[j] for j in range(J)] if J<=20 else [plt.get_cmap('hsv')(j/J) for j in range(J)]
    gp=G.pivot(index='clone',columns='feature',values='state').reindex(columns=E.feature)
    event_names=E.event_label.tolist() if 'event_label' in E else E.feature.tolist()
    plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42,'font.family':'Arial','font.size':9,'svg.hashsalt':'altanalyze3-fish-cartoons'})
    output.mkdir(parents=True,exist_ok=True);files=[]
    common=(f'One patient ({patient}), {len(meta)} specimens, {len(C):,} cells and {J} retained events. Existing shared-genotype MAP assignments; no refitting. '
            'All percentages include all cells, including WT-only observed and unresolved cells. Ring colors identify clone profiles; numbered dot colors identify distinct events. '
            'Filled dots are predicted MUT; empty positions are predicted WT (not proof of genomic WT); question marks would indicate unknown event states. '
            'Profiles combine direct and accepted imputed calls. No distinct complete panel-WT clone was recovered in these fits. '
            'No significance test or biological-replicate confidence interval. Exact profiles can change across starts; consult the initialization-sensitivity heatmap. '
            'Cartoons describe genotype composition; no ancestry edges are inferred.')
    def save(fig,stem,legend):
        fig.canvas.draw();rr=fig.canvas.get_renderer()
        for t in fig.findobj(match=matplotlib.text.Text):
            if t.get_visible() and t.get_text():
                b=t.get_window_extent(rr);require(b.x0>=-1 and b.y0>=-1 and b.x1<=fig.bbox.width+1 and b.y1<=fig.bbox.height+1,f'Clipped text: {stem}: {t.get_text()[:80]}')
        for ext in ['svg','pdf','png']:
            path=output/(stem+'.'+ext);fig.savefig(path,dpi=160);files.append(path.name)
            if ext=='svg':
                tags=[x.tag.rsplit('}',1)[-1] for x in ET.parse(path).getroot().iter()];require('image' not in tags and 'text' in tags,'Not editable vector SVG')
        (output/(stem+'_legend.txt')).write_text(legend+'\n');plt.close(fig)
    def cartoon(ax,k):
        ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
        # An inset with equal aspect keeps the cartoon circular at any panel size.
        ring=ax.inset_axes([.03,.41,.31,.56]);ring.set_aspect('equal');ring.set_xlim(-1.12,1.12);ring.set_ylim(-1.12,1.12);ring.axis('off')
        ring.add_patch(Circle((0,0),.97,facecolor='#f1f1f1',edgecolor=colors[k],lw=4))
        if k not in ['Other','NE']:
            cols=int(np.ceil(np.sqrt(J)));nrows=int(np.ceil(J/cols))
            for j,state in enumerate(gp.loc[k]):
                x=(j%cols-(cols-1)/2)*.38;y=((nrows-1)/2-j//cols)*.38
                if state==1:
                    ring.add_patch(Circle((x,y),.14,facecolor=event_colors[j],edgecolor='#333',lw=.6));ring.text(x,y,str(j+1),ha='center',va='center',fontsize=6,color='black')
                elif state<0:ring.text(x,y,'?',ha='center',va='center',fontsize=8)
        else:ring.text(0,0,'mixed' if k=='Other' else '?',ha='center',va='center',fontsize=10)
        ax.text(.37,.91,k,fontsize=12,fontweight='bold',va='top')
        ax.text(.37,.74,f'Max specimen frequency\n{(V.loc[remaining].sum().max() if k=="Other" else V.loc[k].max()):.2f}%',fontsize=8,va='top')
        label=textwrap.fill(names[k].replace('|',' | '),61,break_long_words=False)
        ax.text(.03,.35,label,fontsize=8,va='top')
    def event_key(fig,y):
        # Column key is independent of clone-ring colors.
        cols=4;nrows=int(np.ceil(J/cols))
        for j,label in enumerate(event_names):
            col=j//nrows;row=j%nrows
            fig.text(.035+col*.245,y-row*.021,f'{j+1}. {label}',fontsize=8,color='#222')
    def band_values(samples):
        return np.array([V.loc[remaining,samples].sum().to_numpy() if k=='Other' else V.loc[k,samples].to_numpy() for k in rows])
    def fish(ax):
        values=band_values(patients);x,smooth=smooth_composition(values)
        require(np.allclose(smooth.sum(0),100) and (smooth>=-1e-10).all(),'Interpolation changed proportions')
        unknown=smooth[-1];bottom=unknown/2
        ax.fill_between(x,0,bottom,color=colors['NE']);ax.fill_between(x,100-bottom,100,color=colors['NE'])
        for i,k in enumerate(rows[:-1]):
            top=bottom+smooth[i];ax.fill_between(x,bottom,top,color=colors[k],edgecolor='white',linewidth=.35)
            for j in range(len(patients)):
                if values[i,j]>=8:
                    y=values[-1,j]/2+values[:i,j].sum()+values[i,j]/2
                    ax.text(j,y,f'{k}\n{values[i,j]:.1f}%',ha='center',va='center',fontsize=8,bbox=dict(facecolor='white',alpha=.65,edgecolor='none',pad=1))
            bottom=top
        for j in range(len(patients)):ax.axvline(j,color='#444',lw=.5,alpha=.5)
        ax.set_xlim(-.06,len(patients)-1+.06);ax.set_ylim(0,100);ax.set_ylabel('Band thickness = % of all specimen cells')
        ax.set_xticks(range(len(patients)),[textwrap.fill(str(meta.loc[s,'stage']),16)+f'\nn={totals[s]:,}' for s in patients],fontsize=9)
        ax.set_title('Patient frequency course (schematic spacing)',loc='left',fontsize=13,pad=15)
        ax.spines[['top','right']].set_visible(False)
    def graft(ax):
        values=band_values(grafts);bottom=np.zeros(len(grafts))
        for i,k in enumerate(rows):
            ax.bar(range(len(grafts)),values[i],bottom=bottom,color=colors[k],edgecolor='#666',linewidth=.4,hatch='///',width=.65)
            for j in range(len(grafts)):
                if values[i,j]>=8:ax.text(j,bottom[j]+values[i,j]/2,f'{k}\n{values[i,j]:.1f}%',ha='center',va='center',fontsize=8,bbox=dict(facecolor='white',alpha=.75,edgecolor='none',pad=1))
            bottom+=values[i]
        require(np.allclose(bottom,100),'Graft denominator changed');ax.set_ylim(0,100);ax.set_xlim(-.65,len(grafts)-.35)
        ax.set_xticks(range(len(grafts)),[textwrap.fill(str(meta.loc[s,'stage']),16)+f'\nn={totals[s]:,}' for s in grafts],fontsize=9)
        ax.set_title('Matched MISTRG',fontsize=12,pad=15);ax.set_yticks([]);ax.spines[['top','right','left']].set_visible(False)
    legend=('Patient bands interpolate the measured specimen proportions with nonnegative cosine joins, at equally spaced stages; time between samples is not to scale and intervening abundance is not observed. '
            'No ancestral nesting or cellular tracking is asserted. Gray no-evidence area is split equally above and below the centered resolved profiles; its combined thickness is the unresolved percentage. '
            'Hatched bars are separate MISTRG specimens, not chronological patient continuations. Top 12 profiles are named; remaining profiles are pooled only as Other and retained in the full cartoon catalog and tables. '+common)
    # Requested combined panel: mutation-dot cartoons above a fish-style course.
    nrows=int(np.ceil(len(rows)/4));fig=plt.figure(figsize=(21,8.8+2.25*nrows))
    fig.suptitle(f'{patient} | clone composition and frequency course',x=.025,y=.98,ha='left',fontsize=19)
    top=.93;bottom=.51;hh=(top-bottom)/nrows
    for i,k in enumerate(rows):cartoon(fig.add_axes([.025+(i%4)*.245,top-(i//4+1)*hh,.238,hh*.97]),k)
    event_key(fig,.48)
    fish(fig.add_axes([.06,.17,.68 if grafts else .89,.20]))
    if grafts:graft(fig.add_axes([.81,.17,.14,.20]))
    fig.text(.025,.025,textwrap.fill(legend,230),va='bottom',fontsize=8)
    save(fig,patient+'_fishplot_cartoons',legend)
    # Standalone fish plot with complete named compact legend.
    fig=plt.figure(figsize=(20,11));fish(fig.add_axes([.06,.31,.55 if grafts else .68,.57]))
    if grafts:graft(fig.add_axes([.64,.31,.12,.57]))
    handles=[Patch(facecolor=colors[k],label=textwrap.fill((k+' = '+names[k]).replace('|',' | '),65,break_long_words=False)) for k in rows]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.78 if grafts else .76,.89),frameon=False,fontsize=7)
    fig.suptitle(f'{patient} | fish-style clone frequency trajectories',x=.025,y=.97,ha='left',fontsize=18)
    fig.text(.025,.025,textwrap.fill(legend,215),va='bottom',fontsize=8);save(fig,patient+'_fishplot',legend)
    # Complete unpooled catalog, including every rare MAP-assigned profile.
    for page,start in enumerate(range(0,len(order),12),1):
        subset=order[start:start+12];nrows=int(np.ceil(len(subset)/4));fig=plt.figure(figsize=(21,4+2.4*nrows));fig.suptitle(f'{patient} | complete mutation-dot clone catalog (page {page})',x=.025,y=.975,ha='left',fontsize=17)
        for i,k in enumerate(subset):cartoon(fig.add_axes([.025+(i%4)*.245,.91-(i//4+1)*(.62/nrows),.238,.60/nrows]),k)
        event_key(fig,.24);fig.text(.025,.025,textwrap.fill('Every nonzero assigned profile is shown without pooling. Dot positions and event numbers are fixed across all cartoons in this patient. '+common,230),va='bottom',fontsize=8)
        save(fig,f'{patient}_clone_cartoons_{page}',common)
    check=dict(patient=patient,cells=len(C),patient_specimens=patients,MISTRG_specimens=grafts,retained_profiles=order,files=files,no_refit=True,no_inferred_ancestry=True,all_cell_denominator=True,interpolation='convex cosine, exact observed knots, nonnegative, sum100',source_hashes={f:__import__('hashlib').sha256((result/f).read_bytes()).hexdigest() for f in ['clone_frequencies.tsv','clone_genotypes.tsv','cell_assignments.tsv.gz']})
    (output/'fishplot_checks.json').write_text(json.dumps(check,indent=2));return check


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--result',required=True);p.add_argument('--output');a=p.parse_args();render(a.result,a.output)
if __name__=='__main__':main()
