"""Complete vector figures for the reviewed detection-tree model. No refitting."""
from pathlib import Path
import argparse, json, textwrap, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch, Circle, FancyArrowPatch
from matplotlib.path import Path as MplPath
from matplotlib.offsetbox import DrawingArea, AnnotationBbox
from matplotlib.text import Text
from matplotlib import font_manager
from .integrity import require, sha256
from .plotting import specimen_order
from .clone_fish import PALETTE, draw_clone_fish
from .detection_reporting import POLICY, DESCRIPTION, tables
from .detection_fishplot import clinical_fish_figure


def centered_mutation_dots(count, dot_radius=4.68, circle_radius=24.):
    """Center each row and the complete point cloud, using at most three rows."""
    if count == 0:
        return np.empty((0, 2))
    rows = min(3, int(np.ceil(count / 4)))
    sizes = [count // rows] * rows
    for row in sorted(range(rows), key=lambda r: abs(r-(rows-1)/2))[:count % rows]:
        sizes[row] += 1
    pitch = 2 * dot_radius + .64
    points = np.array([((col-(size-1)/2)*pitch, ((rows-1)/2-row)*pitch)
                       for row, size in enumerate(sizes) for col in range(size)])
    points -= points.mean(axis=0)
    require(np.max(np.linalg.norm(points, axis=1))+dot_radius < circle_radius,
            'Mutation dots do not fit inside clone circle')
    if count > 1:
        distances = np.linalg.norm(points[:, None]-points[None, :], axis=2)
        require(np.min(distances[np.triu_indices(count, 1)]) >= 2*dot_radius,
                'Mutation dots overlap')
    return points


def require_arial():
    """Use installed Arial or an explicit private font directory; never substitute."""
    directory = os.environ.get('ALTANALYZE_ARIAL_DIR')
    if directory:
        for path in Path(directory).glob('*.ttf'):
            font_manager.fontManager.addfont(str(path))
    for style in ['normal', 'italic']:
        font_manager.findfont(font_manager.FontProperties(family='Arial', style=style),
                              fallback_to_default=False)


# Approved 5801 presentation defaults, shared by every patient tree. Renderers
# must not resize the returned canvas or shrink its labels independently.
TREE_STYLE_NAME = "5801-approved-compact"
TREE_STYLE = {
    "clone_font_pt": 20, "mutation_font_pt": 16,
    "arrow_linewidth_pt": 0.5, "arrowhead_scale": 18,
    "depth_width_inches": 1.15, "circle_radius_pt": 24.0,
    "mutation_dot_radius_pt": 4.68,
}


def clonal_tree_figure(patient, G, parents, active, names, colors, labels,
                       mutation_colors, xy, leaf_count):
    """Publication styling only: preserve every fitted node, edge and inherited event."""
    require_arial()
    maxd = max(x for x, y in xy.values())
    def gain_label(label):
        if label == 'RUNX1-LOH':
            return '+RUNX1-\nLOH'
        if '::' in label:
            return '+'+label.replace('::', '::\n')
        if '_' in label:
            gene, detail = label.split('_', 1)
            return '+'+gene+'\n'+textwrap.fill(detail, 7)
        return textwrap.fill('+'+label, 9)
    multiline = any('_' in label and len(label) > 12 for label in labels)
    with plt.rc_context({'font.family': 'Arial', 'font.size': 18}):
        height = max(14, leaf_count*3+8) if multiline else max(14, leaf_count*1.9+7)
        fig, ax = plt.subplots(figsize=(max(11, maxd*TREE_STYLE["depth_width_inches"]+2), height))
        for n, (x, y) in xy.items():
            if n:
                px, py = xy[int(parents[n])]
                # Straight center-to-center connectors, trimmed in physical points.
                ax.add_patch(FancyArrowPatch((px, py), (x, y), arrowstyle='->',
                                            shrinkA=29, shrinkB=29, connectionstyle='arc3,rad=0',
                                            color='#444', linewidth=TREE_STYLE['arrow_linewidth_pt'],
                                            mutation_scale=TREE_STYLE['arrowhead_scale']))
            area = DrawingArea(64, 64, clip=False)
            area.add_artist(Circle((32, 32), 24, facecolor='#eeeeee',
                                  edgecolor=colors.get(n, '#888888'), lw=3,
                                  linestyle='-' if n in active else '--'))
            events = np.flatnonzero(G[n])
            for j, (dx, dy) in zip(events, centered_mutation_dots(len(events))):
                area.add_artist(Circle((32+dx, 32+dy), 4.68,
                                      facecolor=mutation_colors[j], edgecolor='black', lw=.4))
            ax.add_artist(AnnotationBbox(area, (x, y), frameon=False,
                                        box_alignment=(.5, .5), pad=0))
            ax.annotate(names.get(n, 'A'+str(n)), (x, y), xytext=(0, 34),
                        textcoords='offset points', ha='center', fontsize=TREE_STYLE['clone_font_pt'])
            if n:
                gains = np.flatnonzero(G[n]-G[int(parents[n])])
                gain_lines = [gain_label(labels[j]) for j in gains]
                shift = 13 if y < xy[int(parents[n])][1] else 0
                if any(xy[ch][1] > y for ch in xy if ch and parents[ch] == n):
                    shift = -13
                # Keep long HGVS labels centered so neighbouring gains cannot
                # collide after restoring the compact node spacing.
                if any('_' in labels[j] and len(labels[j]) > 12 for j in gains):
                    shift = 0
                ax.annotate('\n'.join(gain_lines),
                            (x, y), xytext=(shift, -34), textcoords='offset points',
                            ha='center', va='top', fontsize=TREE_STYLE['mutation_font_pt'], fontstyle='italic',
                            bbox=dict(facecolor='white', edgecolor='none', pad=.3))
        ax.set_xlim(-.65, maxd+.75)
        ax.set_ylim(leaf_count-.15, -1)
        ax.axis('off')
        fig.subplots_adjust(left=.035, right=.965, top=.9, bottom=.43)
        fig.suptitle(patient+('\n' if len(patient)>12 else ' | ')+'Inferred gain-only ancestry', fontsize=34, y=.975)
        # Every event is in the legend, including events retained only in latent nodes.
        legend_ax = fig.add_axes([.045, .25, .92, .14])
        legend_ax.axis('off')
        columns = 2 if max(map(len, labels), default=0) > 26 else 3
        rows = int(np.ceil(len(labels)/columns))
        for j, label in enumerate(labels):
            col, row = j % columns, j // columns
            x, y = col/columns, 1-(row+.5)/rows
            legend_ax.plot(x+.008, y, 'o', color=mutation_colors[j],
                           markeredgecolor='black', markeredgewidth=.4, markersize=9.36,
                           transform=legend_ax.transAxes, clip_on=False)
            legend_ax.text(x+.023, y, label, fontsize=TREE_STYLE['mutation_font_pt'], fontstyle='italic', va='center',
                           transform=legend_ax.transAxes)
    return fig


def interpolate_exclusive(F, points=121):
    """Convex cosine interpolation; preserve exclusive all-cell frequencies."""
    require(F.ndim==2 and F.shape[1]>=1 and np.isfinite(F).all() and (F>=-1e-9).all(), 'Invalid fish frequencies')
    require(np.allclose(F.sum(0),100), 'Fish denominators do not sum to 100')
    if F.shape[1]==1:return np.array([0.]),F.copy()
    xx=[];vv=[]
    for j in range(F.shape[1]-1):
        u=np.linspace(0,1,points);w=(1-np.cos(np.pi*u))/2
        sl=slice(None) if j==F.shape[1]-2 else slice(None,-1)
        xx.extend((j+u)[sl]);vv.append((F[:,j,None]*(1-w)+F[:,j+1,None]*w)[:,sl])
    V=np.concatenate(vv,axis=1)
    require(np.allclose(V.sum(0),100) and (V>=-1e-8).all(),'Interpolation changed mass')
    return np.array(xx),V


def render(result, inputs, output, sensitivity_root=None, alluvial_config=None, extended_figures=True, include_unassigned=False):
    result,inputs,output=map(Path,[result,inputs,output])
    require((result/'COMPLETE.json').exists(),'Incomplete fit')
    manifest=json.loads((result/'run_manifest.json').read_text())
    require(manifest['algorithm']['method']=='detection-tree','Wrong inference engine')
    for k,n in [('matrix','calls.npz'),('cells','cells.tsv'),('events','events.tsv')]:
        require(sha256(inputs/n)==manifest['contract']['input_sha256'][k],'Input hash mismatch: '+k)
    c=pd.read_csv(result/'cell_roster.tsv',sep='\t',dtype={'cell_id':str,'sample':str});e=pd.read_csv(result/'event_roster.tsv',sep='\t')
    a=pd.read_csv(result/'cell_assignments.tsv.gz',sep='\t',dtype={'cell_id':str,'sample':str});f=pd.read_csv(result/'clone_frequencies.tsv',sep='\t')
    pc=json.loads((inputs/'plot_source_contract.json').read_text())
    for name,digest in pc['sha256'].items():require(sha256(inputs/name)==digest,'Plot source hash mismatch: '+name)
    require(bool(pc['verified_source_provenance']),'Missing plot source provenance')
    z=np.load(result/'model.npz');src=np.load(inputs/'source_evidence.npz');er=pd.read_csv(inputs/'evidence_rows.tsv',sep='\t')
    require(np.array_equal(z['cell_ids'],src['cell_ids']) and np.array_equal(z['cell_ids'],c.cell_id) and np.array_equal(c.cell_id,a.cell_id),'Cell identities differ')
    require(np.array_equal(z['features'],src['features']) and np.array_equal(z['features'],e.feature),'Event identities differ')
    require(src['matrix'].shape==(len(er),len(c)),'Source shape mismatch')
    calls=np.load(inputs/'calls.npz')
    require(np.array_equal(calls['cell_ids'],z['cell_ids']) and np.array_equal(calls['features'],z['features']), 'Reporting input identity mismatch')
    X=calls['X']; active=z['active'];G=z['G'];parents=z['parents'];J=len(e)
    Q, reported_a, reported_f, audit = tables(X,c,e,z)
    if manifest['algorithm'].get('reporting_policy')==POLICY:
        require(np.array_equal(a.node,reported_a.node),'Saved reporting assignments differ')
        merged=f.merge(reported_f,on=['sample','node'],suffixes=('_saved','_expected'),validate='one_to_one')
        require(len(merged)==len(f)==len(reported_f) and np.allclose(merged.posterior_expected_percent_all_saved,merged.posterior_expected_percent_all_expected),'Saved reporting frequencies differ')
    else:
        # Verify legacy baseline before applying the current reporting policy.
        known=a.node.to_numpy()!=-1;nativeP=z['posterior']
        for row in f.itertuples():
            ix=c['sample'].eq(row.sample).to_numpy();n=ix.sum()
            expected=(~known&ix).sum() if row.node==-1 else nativeP[ix&known,list(active).index(row.node)].sum()
            require(row.all_cells==n and row.cells==((a.node==row.node)&ix).sum(),'Baseline assignment count mismatch')
            require(np.isclose(row.posterior_expected_percent_all,100*expected/n),'Baseline frequency mismatch')
    a,f,P=reported_a,reported_f,Q
    meta=specimen_order(c);clinical=[s for s in meta.index if not bool(meta.loc[s,'is_xenograft'])];grafts=[s for s in meta.index if bool(meta.loc[s,'is_xenograft'])];samples=clinical+grafts
    patient=str(c.patient.iloc[0]);pos=sorted([int(n) for n in active if n],key=lambda n:(G[n].sum(),tuple(G[n]),n));nodes=[0]+pos
    names={0:'WT',-1:'U',**{n:f'SC{i+1}' for i,n in enumerate(pos)}}
    colors={0:'#edcf91',-1:'#b9b9b9',**{n:PALETTE[i] if i<len(PALETTE) else matplotlib.colors.to_hex(plt.get_cmap('tab20')(i-len(PALETTE))) for i,n in enumerate(pos)}}
    labels=e.event_label.astype(str).str.replace('RUNX1_LOH','RUNX1-LOH').tolist()
    comp={n:('No direct or imputed mutant call' if n==0 else 'Unresolved: no target call' if n==-1 else '|'.join(labels[j] for j in np.flatnonzero(G[n]))) for n in nodes}
    F=f.pivot(index='node',columns='sample',values='posterior_expected_percent_all').reindex(index=nodes,columns=samples).to_numpy()
    require(np.isfinite(F).all() and np.allclose(F.sum(0),100),'Missing clone/specimen frequencies')
    output.mkdir(parents=True,exist_ok=True)
    a.to_csv(output/'cell_assignments.tsv.gz',sep='\t',index=False)
    audit.to_csv(output/'assignment_evidence_audit.tsv',sep='\t',index=False)
    np.savez_compressed(output/'reporting_memberships.npz',memberships=Q,nodes=active,cell_ids=z['cell_ids'],features=z['features'])
    (output/'reporting_policy.json').write_text(json.dumps(dict(policy=POLICY,description=DESCRIPTION),indent=2))
    key=pd.DataFrame([dict(clone=names[n],node=n,mutation_combination=comp[n],color=colors[n]) for n in nodes]);key.to_csv(output/'clone_key.tsv',sep='\t',index=False)
    f.merge(key,on='node',validate='many_to_one').to_csv(output/'clone_frequencies.tsv',sep='\t',index=False)
    samplelabels=[('MISTRG: ' if s in grafts else '')+str(meta.loc[s,'stage'])+'\n'+s+'\nn='+format(int(c['sample'].eq(s).sum()),',') for s in samples]
    plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42,'font.family':'DejaVu Sans','font.size':9,'svg.hashsalt':'detection-tree'})
    base=(f'{patient}: one patient, {len(clinical)} clinical and {len(grafts)} matched xenograft specimens; {len(c):,} cells and all {J} established events retained. '
      'Clone genotypes and gain-only ancestry are model predictions from direct RNA and accepted imputed evidence. Seed 222, pooling strength 2, RNA concentration 2; no fixed founder genes or expected frequencies. '
      'WT means no direct or imputed mutant call, including no-call and LOH-negative-only cells. Mutation-positive memberships are conditioned on the original malignant clones; percentages divide by ALL specimen cells. Source calls and unknown loci are unchanged. '
      'One patient is the biological unit; cells are not independent patients. No hypothesis test or biological confidence interval; sparse calibration, RNA detection and correlated surrogates limit interpretation. ')
    files=[]
    def save(fig,stem,extra):
        legend=base+extra
        if stem in ['clone_cartoons','fishplot_cartoons']:
            legend+=' Mutation dot key: '+'; '.join(str(j+1)+'='+label for j,label in enumerate(labels))+'.'
        is_tree = stem == 'clonal_tree'
        display_legend=legend
        if stem=='fishplot':
            display_legend=(f'{patient}: one patient; {len(clinical)} clinical and {len(grafts)} matched graft specimens. Areas show exclusive inferred clone frequencies (% of ALL cells, including WT). '
                'WT gradient = no direct or imputed MUT call, including no-call/LOH-negative-only cells. Positive cells retain relative malignant-clone memberships. Model uncertainty is not shown. '
                'Prediagnosis is schematic; other curves interpolate sampled stages. No hypothesis test. Clone keys show complete combinations; fitted ancestry and original calls are unchanged.')
        fig.text(.035,.015,textwrap.fill(display_legend,int(fig.get_figwidth()*(10 if is_tree else 19 if stem=='fishplot' else 9))),
                 fontsize=13 if is_tree else 7 if stem=='fishplot' else 8, va='bottom',
                 fontfamily='Arial' if is_tree else 'DejaVu Sans')
        fig.canvas.draw()
        for t in fig.findobj(match=matplotlib.text.Text):
            if t.get_visible() and t.get_text():
                b=t.get_window_extent(fig.canvas.get_renderer())
                require(b.x0>=-1 and b.y0>=-1 and b.x1<=fig.bbox.width+1 and b.y1<=fig.bbox.height+1,'Text outside canvas '+stem+': '+t.get_text()[:70])
        for ext in ['svg','pdf','png']:
            path=output/(stem+'.'+ext);fig.savefig(path,dpi=120);files.append(path.name)
            if ext=='svg':require('<image' not in path.read_text() and '<text' in path.read_text(),'Noneditable SVG')
        (output/(stem+'_legend.txt')).write_text(legend);plt.close(fig)
    def draw_key(ax):
        ax.axis('off');y=1.
        step=9*1.5/72/(ax.figure.get_figheight()*ax.get_position().height)
        for n in nodes:
            lines=textwrap.wrap(names[n]+' = '+comp[n].replace('|',' | '),105,break_long_words=False)
            ax.plot([0],[y-step*.4],marker='s',color=colors[n],markersize=8,transform=ax.transAxes,clip_on=False)
            ax.text(.02,y,'\n'.join(lines),va='top',fontsize=9,linespacing=1.5,transform=ax.transAxes)
            y-=step*len(lines)
    keylines=sum(len(textwrap.wrap(names[n]+' = '+comp[n].replace('|',' | '),105,break_long_words=False)) for n in nodes)
    K=len(nodes);height=max(12,7+keylines*.23)
    fig=plt.figure(figsize=(22,height));gs=fig.add_gridspec(2,2,left=.17,right=.97,top=.9,bottom=.20,width_ratios=[1,1.5],height_ratios=[max(4,keylines*.24),max(3,K*.33)],wspace=.2,hspace=.48)
    ax=fig.add_subplot(gs[0,0]);geno=np.column_stack([G[n] if n>=0 else np.full(J,2) for n in nodes]);ax.pcolormesh(geno,cmap=ListedColormap(['#eeeeee','#d95f02','#b9c0c7']),vmin=0,vmax=2,edgecolors='white');ax.set_ylim(J,0);ax.set_yticks(np.arange(J)+.5,labels,fontsize=8);ax.set_xticks(np.arange(K)+.5,[names[n] for n in nodes]);ax.set_title('A  Predicted mutation composition');ax.tick_params(length=0)
    ax=fig.add_subplot(gs[0,1]);draw_key(ax);ax.set_title('Clone identities: complete mutation combinations',loc='left')
    ax=fig.add_subplot(gs[1,:]);im=ax.pcolormesh(F,cmap='YlGnBu',vmin=0,vmax=100,edgecolors='white');ax.set_ylim(K,0);ax.set_yticks(np.arange(K)+.5,[names[n] for n in nodes]);ax.set_xticks(np.arange(len(samples))+.5,samplelabels,fontsize=8);ax.xaxis.tick_top();ax.set_title('B  Exclusive clone frequency (% of all cells)',pad=52)
    if grafts:ax.axvline(len(clinical),color='black',lw=2)
    for i,j in np.ndindex(F.shape):ax.text(j+.5,i+.5,f'{F[i,j]:.1f}',ha='center',va='center',fontsize=9,color='white' if F[i,j]>60 else 'black')
    cb=fig.colorbar(im,ax=ax,fraction=.017,pad=.015);cb.solids.set_rasterized(False)
    fig.suptitle(patient+' | Mutation-defined clones across specimens',fontsize=17,y=.97)
    save(fig,'clone_composition_frequency','Orange = predicted mutation; pale gray = predicted absence; blue-gray = unresolved. Frequency columns sum to 100% before rounding. Full combinations, including rare clones, are displayed. Matched grafts are separate columns.')
    # Complete source evidence, vector run lengths, every cell.
    R=len(er);height=max(15,7+R*.17+keylines*.17);fig,axs=plt.subplots(1,len(samples),squeeze=False,figsize=(27,height),gridspec_kw={'left':.19,'right':.99,'top':.84,'bottom':.35,'wspace':.04})
    palette=['white','#acd9a5','#f06b48','#9b59b6'];orders=[];widths=[]
    for si,s in enumerate(samples):
        ax=axs[0,si];ix=np.flatnonzero(c['sample'].eq(s));rank={n:i for i,n in enumerate(nodes)};ix=np.array(sorted(ix,key=lambda i:(rank[int(a.node.iloc[i])],str(c.cell_id.iloc[i]))));mat=src['matrix'][:,ix];verts=[];cc=[]
        for row,vals in enumerate(mat):
            br=np.r_[0,np.flatnonzero(vals[1:]!=vals[:-1])+1,len(vals)]
            for lo,hi in zip(br[:-1],br[1:]):
                if vals[lo]:verts.append([(lo,row),(hi,row),(hi,row+1),(lo,row+1)]);cc.append(palette[int(vals[lo])])
        ax.add_collection(PolyCollection(verts,facecolors=cc,edgecolors='none'));ax.set_xlim(0,len(ix));ax.set_ylim(R,-2);ax.set_xticks([]);ax.set_yticks(np.arange(R)+.5,er.label if si==0 else []);ax.tick_params(length=0,labelsize=8);ax.set_title(samplelabels[si],fontsize=10)
        for j in range(1,R):
            if er.event.iloc[j]!=er.event.iloc[j-1]:ax.axhline(j,color='#666',lw=.4)
        start=0
        for n in nodes:
            num=int((a.node.iloc[ix]==n).sum());widths.append(dict(sample=s,clone=names[n],MAP_cells=num,all_cells=len(ix),MAP_percent=100*num/len(ix)))
            if num:
                ax.add_patch(plt.Rectangle((start,-1.6),num,1.3,color=colors[n],lw=0));ax.axvline(start,color='black',lw=.5)
                if num/len(ix)>.045:ax.text(start+num/2,-.9,names[n],ha='center',va='center',fontsize=7)
            start+=num
        require(start==len(ix),'Lost heatmap cells')
        orders.extend(dict(sample=s,column=k,cell_id=c.cell_id.iloc[i],clone=names[int(a.node.iloc[i])]) for k,i in enumerate(ix))
    fig.suptitle(patient+' | Direct and imputed evidence for every cell',fontsize=18,y=.97)
    fig.legend(handles=[Patch(facecolor=palette[i],edgecolor='#aaa',label=t) for i,t in enumerate(['No call / unavailable','WT/reference','Mutation','MUT + WT conflict'])],loc='upper center',bbox_to_anchor=(.6,.93),ncol=4)
    ax=fig.add_axes([.04,.13,.92,.17]);draw_key(ax)
    save(fig,'cell_source_heatmaps','Each row is one event/source; all cells are ordered by MAP clone then barcode. Strip widths count maximum reporting memberships; frequency panels sum fractional reporting memberships. White is missing, never WT. PV/MPV: passenger inference; CNV: copy-number inference; DEG/DEI: expression/isoform inference; LOH: loss-of-heterozygosity inference. Strip labels below 4.5% are omitted for space; their cells and clone keys remain.')
    pd.DataFrame(orders).to_csv(output/'cell_order.tsv.gz',sep='\t',index=False);pd.DataFrame(widths).to_csv(output/'MAP_widths.tsv',sep='\t',index=False)
    # Numbered mutation dots and complete cartoons (all nodes, all events).
    mutation_colors=[plt.get_cmap('hsv')(j/max(J,1)) for j in range(J)]
    def cartoon(ax,n,x,y,r=.24):
        # Point-sized symbols stay circular even on deep, narrow event trees.
        area=DrawingArea(64,64,clip=False);radius=24.;cx=cy=32.
        area.add_artist(Circle((cx,cy),radius,facecolor='#eeeeee',edgecolor=colors.get(n,'#888888'),lw=3,linestyle='-' if n in active else '--'))
        js=np.flatnonzero(G[n]) if n>=0 else []
        for k,j in enumerate(js):
            cols=4;rr=int(np.ceil(len(js)/cols));dx=(k%cols-1.5)*radius*.38;dy=((rr-1)/2-k//cols)*radius*.38
            area.add_artist(Circle((cx+dx,cy+dy),radius*.15,facecolor=mutation_colors[j],edgecolor='black',lw=.4))
            area.add_artist(Text(cx+dx,cy+dy,str(j+1),ha='center',va='center',fontsize=5))
        ax.add_artist(AnnotationBbox(area,(x,y),frameon=False,box_alignment=(.5,.5),pad=0))
    fig,ax=plt.subplots(figsize=(20,max(11,K*1.05+3)));ax.set_xlim(0,10);ax.set_ylim(K+.8,-.6);ax.axis('off')
    for i,n in enumerate(nodes):
        cartoon(ax,n,.5,i,.28);ax.text(1,i,textwrap.fill(names[n]+' = '+comp[n].replace('|',' | '),125,break_long_words=False),va='center',fontsize=10)
    fig.subplots_adjust(left=.03,right=.98,top=.91,bottom=.17);fig.suptitle(patient+' | Complete clone mutation cartoons',fontsize=17)
    save(fig,'clone_cartoons','Each numbered dot identifies an event in event_number_key.tsv (the ordering of the composition heatmap). Circle outline colors match the fishplot and heatmap clone strips. All active mutant clones and WT are shown; cartoons encode inferred genotypes, not directly observed haplotypes.')
    pd.DataFrame({'number':np.arange(1,J+1),'event':labels,'feature':e.feature}).to_csv(output/'event_number_key.tsv',sep='\t',index=False)
    # Actual model tree, retaining latent ancestors rather than inventing shortcuts.
    children={n:[i for i,p in enumerate(parents) if i and p==n] for n in range(len(G))};xy={};cursor=[0]
    def layout(n,d):
        ys=[layout(ch,d+1) for ch in children[n]]
        y=float(np.mean(ys)) if ys else float(cursor[0]);cursor[0]+=int(not ys);xy[n]=(d,y);return y
    layout(0,0)
    fig=clonal_tree_figure(patient,G,parents,active,names,colors,labels,mutation_colors,xy,cursor[0])
    save(fig,'clonal_tree','Arrows are fitted event-gain relationships, not independently verified lineage tracing. A-prefixed dashed circles are latent ancestors with no sampled-state component; they have no invented frequency. Colored dots retain every inherited mutation, identified by the italicized mutation legend. Circle outline colors identify sampled clones consistently across companion plots; clone_key.tsv provides complete combinations.')
    # Rendering adapted from Xuan's Python Matplotlib implementation, not his results.
    def draw_fish(ax):
        return draw_clone_fish(ax,np.vstack([F[:,:len(clinical)],np.zeros(len(clinical))]),nodes+[-1],parents,colors,
            [str(meta.loc[s,'stage']).replace('Diagnosis / pre-HMA','Diagnosis /\npre-HMA')+'\nn='+format(int(c['sample'].eq(s).sum()),',') for s in clinical])
    def graftbars(ax):
        for k,s in enumerate(grafts):
            j=samples.index(s);bottom=0
            for i,n in enumerate(nodes):ax.bar(k,F[i,j],bottom=bottom,color=colors[n],width=.6);bottom+=F[i,j]
        ax.set_ylim(0,100);ax.set_xticks(range(len(grafts)),['MISTRG\n'+str(meta.loc[s,'stage'])+'\n'+s for s in grafts],fontsize=8);ax.set_title('Matched graft',fontsize=10,pad=45)
    fig=clinical_fish_figure(patient,F,nodes,parents,colors,names,comp,meta,clinical,grafts,c,include_unassigned=include_unassigned)
    save(fig,'fishplot','Approved WT reporting: all cells without accepted direct/imputed MUT are WT, including no-call and LOH-negative-only cells. Positive-cell memberships are conditioned on malignant clones. Rendering adapted from Xuan; no VAF-derived frequencies. Prediagnosis schematic; measured frequencies exact. Canvas 15.86 by 7.18 inches; all geometry and text editable.')
    fig=plt.figure(figsize=(23,max(13,K*.65+5)));gs=fig.add_gridspec(1,2,left=.04,right=.97,top=.88,bottom=.2,width_ratios=[1.2,1]);ax=fig.add_subplot(gs[0,0]);ax.set_xlim(0,10);ax.set_ylim(K,-1);ax.axis('off')
    for i,n in enumerate(nodes):cartoon(ax,n,.4,i,.25);ax.text(.9,i,textwrap.fill(names[n]+' = '+comp[n].replace('|',' | '),65,break_long_words=False),va='center',fontsize=8)
    draw_fish(fig.add_subplot(gs[0,1]));fig.suptitle(patient+' | Clone cartoons and clinical fishplot',fontsize=17)
    save(fig,'fishplot_cartoons','Left: complete clone combinations; numbered dots refer to the supplied event key. Right: same-colored nested fishplot, clinical specimens only. Matched graft bars and exact all-specimen percentages are in fishplot.svg and clone_composition_frequency.svg. Interpolated curves are schematic between sampling stages.')
    fig,ax=plt.subplots(figsize=(max(14,len(samples)*2),12));bottom=np.zeros(len(samples))
    for i,n in enumerate(nodes):ax.bar(range(len(samples)),F[i],bottom=bottom,color=colors[n],label=names[n]);bottom+=F[i]
    ax.set_ylim(0,100);ax.set_ylabel('% of all cells');ax.set_xticks(range(len(samples)),samplelabels,fontsize=9);ax.legend(ncol=5,loc='upper center',bbox_to_anchor=(.5,1.13));fig.subplots_adjust(top=.82,bottom=.28,left=.08,right=.97);fig.suptitle(patient+' | Specimen clonal composition',fontsize=17)
    save(fig,'specimen_composition','Stacked bars show exclusive reported clone frequencies for each separate specimen. Colors and complete mutation combinations are defined in clone_key.tsv and the companion composition heatmap. All WT cells remain in every denominator.')
    if sensitivity_root:
        rows=[]
        for seed in [222,111,333]:
            rr=Path(sensitivity_root)/f'{patient}_seed{seed}';require((rr/'COMPLETE.json').exists(),'Incomplete sensitivity run');zz=np.load(rr/'model.npz');_,_,ff,_=tables(X,c,e,zz);require(np.array_equal(zz['cell_ids'],z['cell_ids']) and np.array_equal(zz['features'],z['features']),'Sensitivity identities differ')
            for s in samples:
                sf=ff[ff['sample']==s].set_index('node')
                for j in range(J):rows.append(dict(seed=seed,sample=s,event=labels[j],percent=sum(sf.loc[int(n),'posterior_expected_percent_all']*zz['G'][n,j] for n in zz['active'])))
        sens=pd.DataFrame(rows);sens.to_csv(output/'initialization_sensitivity.tsv',sep='\t',index=False)
        fig,axs=plt.subplots(1,len(samples),squeeze=False,figsize=(max(17,len(samples)*3),max(10,J*.5+4)),sharey=True)
        for k,s in enumerate(samples):
            ax=axs[0,k];v=sens[sens['sample']==s].pivot(index='event',columns='seed',values='percent').reindex(labels);ax.hlines(range(J),v.min(1),v.max(1),color='#777',lw=3);ax.scatter(v[222],range(J),c='#0072B2',s=24,label='Reviewed seed 222');ax.set_xlim(-2,102);ax.set_xticks([0,50,100]);ax.set_title(samplelabels[k],fontsize=9);ax.set_xlabel('% of all cells');ax.set_yticks(range(J),labels,fontsize=8);ax.set_ylim(J-.5,-.5)
        fig.subplots_adjust(left=.23,right=.98,top=.83,bottom=.24,wspace=.08);fig.suptitle(patient+' | Initialization sensitivity of inferred event prevalence',fontsize=16)
        save(fig,'optimization_sensitivity','Blue dots: primary seed 222. Gray intervals: minimum to maximum event prevalence across seeds 111, 222 and 333, not confidence intervals. Event prevalence sums all clones containing that event, so rows need not sum to 100%. Additional starts do not replace the reviewed primary seed. Complete alternative trees and genotypes are retained.')
    checks=dict(reporting_policy=POLICY,reporting_module_sha256=sha256(Path(__file__).with_name('detection_reporting.py')),fish_layout_sha256=sha256(Path(__file__).with_name('detection_fishplot.py')),patient=patient,cells=len(c),events=J,source_rows=R,specimens=len(samples),grafts=len(grafts),all_cells_retained=True,all_events_retained=True,posterior_frequency_recomputed=True,MAP_counts_recomputed=True,fish_mass_verified=True,fish_prediagnosis_schematic=True,fish_canvas_inches=[15.86,7.18],fish_renderer_sha256=sha256(Path(__file__).with_name('clone_fish.py')),editable_vector_SVG=True,renderer_sha256=sha256(Path(__file__)),model_sha256=sha256(result/'model.npz'),plot_source_sha256=pc['sha256'],files=files)
    if extended_figures:
        from ...visualization.clone_source_heatmap import render as render_compact
        from ...visualization.clone_alluvial import render as render_alluvial
        checks['compact_source_heatmap'] = render_compact(inputs, output, output/'compact_source_heatmap', patient, include_unassigned=include_unassigned)
        checks['state_clone_alluvial'] = render_alluvial(output, output/'state_clone_alluvial', alluvial_config, include_unassigned=include_unassigned)
        checks['extended_defaults'] = True
    (output/'figure_checks.json').write_text(json.dumps(checks,indent=2));return checks


def main():
    p=argparse.ArgumentParser();p.add_argument('--result',required=True);p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--sensitivity-root');p.add_argument('--alluvial-config');p.add_argument('--legacy-vector-only',action='store_true',help='Explicitly replay the historical eight-figure suite');p.add_argument('--include-unassigned',action='store_true',default=False,help='Include U in publication reporting');a=p.parse_args();print(json.dumps(render(a.result,a.input,a.output,a.sensitivity_root,a.alluvial_config,not a.legacy_vector_only,a.include_unassigned),indent=2))
if __name__=='__main__':main()
