"""Editable longitudinal figures for supported mutation combinations."""
from pathlib import Path
import argparse,json,textwrap,xml.etree.ElementTree as ET
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Patch
from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap
from .integrity import require,sha256
from .fishplot import smooth_composition
from .plotting import specimen_order


def render(result,output=None,max_main=12):
    r=Path(result);o=Path(output) if output else r/'figures';o.mkdir(parents=True,exist_ok=True)
    require((r/'COMPLETE.json').exists(),'Incomplete result')
    C=pd.read_csv(r/'cell_roster.tsv',sep='\t');A=pd.read_csv(r/'cell_assignments.tsv.gz',sep='\t');E=pd.read_csv(r/'event_roster.tsv',sep='\t');F=pd.read_csv(r/'clone_frequencies.tsv',sep='\t');G=pd.read_csv(r/'clone_genotypes.tsv',sep='\t')
    z=np.load(r/'accepted_calls.npz');X=z['X'];require(z['cell_ids'].tolist()==C.cell_id.tolist()==A.cell_id.tolist(),'Cell order mismatch');require(z['features'].tolist()==E.feature.tolist(),'Event order mismatch')
    patient=str(C.patient.iloc[0]);meta=specimen_order(C)
    if 'is_xenograft' not in meta:meta['is_xenograft']=False
    clinical=meta.index[~meta.is_xenograft].tolist();grafts=meta.index[meta.is_xenograft].tolist();samples=clinical+grafts;N=C['sample'].value_counts()
    allgroups=sorted([k for k in F.clone.unique() if k.startswith('SC')],key=lambda x:int(x[2:]));special=['WT','WTobs','NE'];shown=allgroups[:max_main];rest=allgroups[max_main:];display=shown+(['Other'] if rest else [])+special
    V=F.pivot(index='clone',columns='sample',values='percent_all').reindex(index=allgroups+special,columns=samples,fill_value=0).fillna(0)
    names=F.drop_duplicates('clone').set_index('clone').composition.to_dict();names.update(WT='Panel WT calls',WTobs='WT-only calls; incomplete panel',NE='No event calls',Other='Other supported combinations')
    colors={g:plt.get_cmap('tab20').colors[i%20] for i,g in enumerate(allgroups)};colors.update(WT='#56A36F',WTobs='#B5D9BE',NE='#eeeeee',Other='#9c9c9c')
    uncertain=E['role'].astype(str).str.contains('inherited_status_uncertain') if 'role' in E else pd.Series(False,index=E.index)
    labels=[str(e.event_label)+(' [status uncertain]' if uncertain[j] else '') for j,e in E.iterrows()];J=len(E)
    plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42,'font.family':'Arial','font.size':9,'svg.hashsalt':'supported-clones-20261116'})
    files=[]
    common=(f'{patient}: one patient, {len(samples)} specimens ({len(grafts)} MISTRG), {len(C):,} cells and {J} retained events. '
      'SC groups share the same accepted MUT combination; they are not proven distinct complete genetic clones. Direct and accepted imputed calls are used equally; no additional genotype imputation. '
      'Unlisted events are unresolved or WT according to each cell, not automatically WT. All-cell denominators include WT, WTobs and NE. WT means calls at every panel event are WT, not whole-genome proof; WTobs means WT-only partial coverage; NE means no calls. '
      'No statistical test, confidence interval or inferred ancestry. Upstream imputation accuracy and unequal call availability limit biological prevalence interpretation.')
    if uncertain.any():common+=' DDX10::SKA3 is retained as an uncertain-inherited-status marker; somatic clone specificity is not established.'
    def save(fig,stem,legend):
        fig.canvas.draw();rr=fig.canvas.get_renderer()
        for t in fig.findobj(match=matplotlib.text.Text):
            if t.get_visible() and t.get_text():
                b=t.get_window_extent(rr);require(b.x0>=-1 and b.y0>=-1 and b.x1<=fig.bbox.width+1 and b.y1<=fig.bbox.height+1,f'Text clipping: {stem}: {t.get_text()[:80]}')
        for ext in ['svg','pdf','png']:
            p=o/(stem+'.'+ext);fig.savefig(p,dpi=140);files.append(p.name)
            if ext=='svg':
                tags=[x.tag.rsplit('}',1)[-1] for x in ET.parse(p).getroot().iter()];require('text' in tags and 'image' not in tags,'Noneditable SVG')
        (o/(stem+'_legend.txt')).write_text(legend+'\n');plt.close(fig)
    def footer(fig,text):fig.text(.025,.025,textwrap.fill(text,185),fontsize=8,va='bottom')
    def combo_name(k):return textwrap.fill(f'{k} = '+names[k].replace('|',' | '),53,break_long_words=False)
    def values(groups,samples):return np.array([V.loc[rest,samples].sum().to_numpy() if k=='Other' else V.loc[k,samples].to_numpy() for k in groups])
    def course(ax,groups):
        v=values(groups,clinical);x,y=smooth_composition(v);require(np.allclose(y.sum(0),100),'Mass not conserved')
        bottom=np.zeros(len(x))
        for i,k in enumerate(groups):
            ax.fill_between(x,bottom,bottom+y[i],color=colors[k],edgecolor='white',lw=.35)
            for j in range(len(clinical)):
                if v[i,j]>=8:ax.text(j,v[:i,j].sum()+v[i,j]/2,f'{k}\n{v[i,j]:.1f}%',ha='center',va='center',fontsize=8,bbox=dict(facecolor='white',alpha=.75,edgecolor='none',pad=1))
            bottom+=y[i]
        for j in range(len(clinical)):ax.axvline(j,color='#555',alpha=.35,lw=.5)
        ax.set_xlim(-.06,len(clinical)-.94);ax.set_ylim(0,100);ax.set_ylabel('% of all specimen cells');ax.set_xticks(range(len(clinical)),[textwrap.fill(str(meta.loc[s,'stage']),15)+f'\nn={N[s]:,}' for s in clinical]);ax.set_title('Supported-combination frequencies (schematic spacing)',loc='left');ax.spines[['top','right']].set_visible(False)
    def graft(ax,groups):
        v=values(groups,grafts);bottom=np.zeros(len(grafts))
        for i,k in enumerate(groups):ax.bar(range(len(grafts)),v[i],bottom=bottom,color=colors[k],hatch='//',edgecolor='white',lw=.4);bottom+=v[i]
        ax.set_ylim(0,100);ax.set_xticks(range(len(grafts)),[textwrap.fill(str(meta.loc[s,'stage']),15)+f'\nn={N[s]:,}' for s in grafts]);ax.set_title('Matched MISTRG',loc='left');ax.spines[['top','right']].set_visible(False)
    def card(ax,k):
        ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off');q=ax.inset_axes([.02,.40,.31,.56]);q.set_aspect('equal');q.axis('off');q.set_xlim(-1.1,1.1);q.set_ylim(-1.1,1.1);q.add_patch(Circle((0,0),.96,facecolor='#f3f3f3',edgecolor=colors[k],lw=4))
        if k in allgroups:
            gp=G[G.clone==k].set_index('feature').loc[E.feature,'state'].to_numpy();cols=int(np.ceil(np.sqrt(J)));nr=int(np.ceil(J/cols))
            for j,s in enumerate(gp):
                if s==1:
                    x=(j%cols-(cols-1)/2)*.37;y=((nr-1)/2-j//cols)*.37;q.add_patch(Circle((x,y),.14,facecolor=plt.get_cmap('tab20').colors[j%20],edgecolor='black',lw=.6));q.text(x,y,str(j+1),ha='center',va='center',fontsize=6)
        else:q.text(0,0,{'WT':'WT','WTobs':'WT?','NE':'?','Other':'mixed'}[k],ha='center',va='center')
        ax.text(.37,.92,k,fontsize=12,fontweight='bold',va='top');mx=float(V.loc[rest].sum().max()) if k=='Other' else float(V.loc[k].max());ax.text(.37,.73,f'Max specimen fraction\n{mx:.2f}%',fontsize=8,va='top');ax.text(.02,.34,textwrap.fill(names[k].replace('|',' | '),58,break_long_words=False),fontsize=8,va='top')
    def key(fig,y):
        nr=int(np.ceil(J/4))
        for j,label in enumerate(labels):fig.text(.03+(j//nr)*.24,y-(j%nr)*.020,f'{j+1}. {label}',fontsize=8)
    legend=common+' Colored dots are accepted MUT events; absent dots do not establish WT. Main panels show12 combinations ranked by maximum specimen fraction plus Other; catalogs and tables retain all groups. Curves interpolate sampled fractions without claiming tracked cells, ancestry, elapsed-time scaling or tumor growth. Hatched bars are separate matched MISTRG specimens.'
    nr=int(np.ceil(len(display)/4));height=9+nr*2.2
    fig=plt.figure(figsize=(21,height));fig.suptitle(f'{patient} | supported mutation combinations and frequencies',x=.025,ha='left',y=.985,fontsize=19)
    gs=fig.add_gridspec(nr,4,left=.025,right=.98,bottom=.47,top=.93,hspace=.16,wspace=.04)
    for i,k in enumerate(display):card(fig.add_subplot(gs[i//4,i%4]),k)
    key(fig,.44);course(fig.add_axes([.07,.17,.70 if grafts else .88,.18]),display)
    if grafts:graft(fig.add_axes([.82,.17,.15,.18]),display)
    footer(fig,legend);save(fig,patient+'_fishplot_cartoons',legend)
    fig=plt.figure(figsize=(20,10));fig.suptitle(f'{patient} | supported-combination frequency course',x=.025,ha='left',fontsize=18)
    course(fig.add_axes([.07,.51,.69 if grafts else .89,.36]),display)
    if grafts:graft(fig.add_axes([.82,.51,.15,.36]),display)
    handles=[Patch(facecolor=colors[k],label=combo_name(k)) for k in display];fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.03,.41),ncol=4,fontsize=8,frameon=False);footer(fig,legend);save(fig,patient+'_fishplot',legend)
    for page,start in enumerate(range(0,len(allgroups),12),1):
        ks=allgroups[start:start+12];nr=int(np.ceil(len(ks)/4));fig=plt.figure(figsize=(21,4+nr*2.3));fig.suptitle(f'{patient} | all supported combinations — page {page}',x=.025,ha='left',fontsize=18)
        gs=fig.add_gridspec(nr,4,left=.025,right=.98,bottom=.30,top=.91,hspace=.13,wspace=.04)
        for i,k in enumerate(ks):card(fig.add_subplot(gs[i//4,i%4]),k)
        key(fig,.25);footer(fig,common+' Every nonzero combination is included without pooling. Dot absence does not imply WT.');save(fig,f'{patient}_cartoons_{page}',common)
    # Composition and exact frequency tables, paginated without excluding profiles.
    for page,start in enumerate(range(0,len(allgroups)+3,18),1):
        ks=(allgroups+special)[start:start+18];h=4+len(ks)*.52;fig=plt.figure(figsize=(21,h));fig.suptitle(f'{patient} | accepted mutation combinations and all-cell percentages — page {page}',x=.025,ha='left',fontsize=17)
        ax=fig.add_axes([.36,.23,.26,.64]);fx=fig.add_axes([.68,.23,.29,.64]);mat=[]
        for k in ks:
            if k in allgroups:mat.append((G[G.clone==k].set_index('feature').loc[E.feature,'state'].to_numpy()==1).astype(int))
            else:mat.append(np.zeros(J))
        ax.pcolormesh(np.arange(J+1),np.arange(len(ks)+1),np.array(mat),cmap=ListedColormap(['#eeeeee','#e76f51']),vmin=0,vmax=1,edgecolor='white',linewidth=.6);ax.invert_yaxis();ax.set_yticks(np.arange(len(ks))+.5,[combo_name(k) for k in ks],fontsize=8);ax.set_xticks(np.arange(J)+.5,[str(i+1) for i in range(J)],fontsize=9);ax.set_title('Accepted MUT (orange); other states unresolved here',fontsize=10)
        v=V.loc[ks,samples].to_numpy();fx.pcolormesh(np.arange(len(samples)+1),np.arange(len(ks)+1),v,cmap='Blues',vmin=0,vmax=100,edgecolor='white',linewidth=.6);fx.invert_yaxis();fx.set_yticks([]);fx.set_xticks(np.arange(len(samples))+.5,[textwrap.fill(str(meta.loc[s,'stage']),12) for s in samples],fontsize=8);fx.set_title('% of all specimen cells',fontsize=11)
        for i in range(len(ks)):
            for j in range(len(samples)):fx.text(j+.5,i+.5,f'{v[i,j]:.2f}',ha='center',va='center',fontsize=9,color='white' if v[i,j]>50 else 'black')
        key(fig,.16);footer(fig,common+' Orange indicates membership in the supported MUT set; gray is not a WT assertion. Numbers are exact sampled percentages; MISTRG columns are independent specimens.');save(fig,f'{patient}_composition_frequency_{page}',common)
    # All cells and events, pure SVG vector run-length polygons; no subsampling.
    fig,axes=plt.subplots(1,len(samples),figsize=(22,7),sharey=True);axes=np.atleast_1d(axes);fig.subplots_adjust(left=.14,right=.985,top=.79,bottom=.29,wspace=.08);rank={k:i for i,k in enumerate(allgroups+special)}
    for ax,s in zip(axes,samples):
        ix=np.flatnonzero(C['sample'].to_numpy()==s);ix=np.array(sorted(ix,key=lambda i:(rank[A.clone.iloc[i]],A.input_states.iloc[i])));q=X[ix].T;verts=[];cs=[]
        for j,row in enumerate(q):
            edges=np.r_[0,np.flatnonzero(row[1:]!=row[:-1])+1,len(row)]
            for a,b in zip(edges[:-1],edges[1:]):verts.append([(a,j),(b,j),(b,j+1),(a,j+1)]);cs.append({-1:'#eeeeee',0:'#aad5ad',1:'#e76f51'}[int(row[a])])
        ax.add_collection(PolyCollection(verts,facecolors=cs,edgecolors='none',rasterized=False));ax.set_xlim(0,len(ix));ax.set_ylim(J,0);ax.set_xticks([]);ax.set_title(textwrap.fill(str(meta.loc[s,'stage']),15)+f'\nn={len(ix):,}',fontsize=10);ax.set_yticks(np.arange(J)+.5,labels,fontsize=8)
    fig.suptitle(f'{patient} | every cell’s accepted MUT / WT / unknown evidence',x=.025,ha='left',fontsize=18)
    fig.legend(handles=[Patch(color='#e76f51',label='MUT: direct or accepted imputed'),Patch(color='#aad5ad',label='WT call'),Patch(color='#eeeeee',label='Unknown')],loc='upper center',bbox_to_anchor=(.57,.91),ncol=3,frameon=False)
    footer(fig,common+' Rows are all retained events; columns are all cells, ordered by supported combination and then full ternary call pattern within each specimen. Every sample panel has equal display width; n gives its denominator.');save(fig,patient+'_serial_evidence',common)
    checks=dict(patient=patient,cells=len(C),events=J,samples=samples,MISTRG=len(grafts),groups=len(allgroups),files=files,all_cells_rendered=True,all_events_retained=True,no_fitted_genotype_completion=True,source_hashes={p:sha256(r/p) for p in ['cell_assignments.tsv.gz','clone_frequencies.tsv','accepted_calls.npz']})
    (o/'figure_checks.json').write_text(json.dumps(checks,indent=2));return checks

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--result',required=True);p.add_argument('--output');a=p.parse_args();render(a.result,a.output)
