"""Approved longitudinal fishplot layout, shared by the default renderer."""
import textwrap
import numpy as np
import matplotlib.pyplot as plt
from .clone_fish import draw_clone_fish


def clinical_fish_figure(patient, F, nodes, parents, colors, names, comp,
                         meta, clinical, grafts, cells, include_unassigned=False):
    F = np.asarray(F, dtype=float).copy()
    nodes = list(nodes)
    if not include_unassigned and -1 in nodes:
        keep = [i for i,n in enumerate(nodes) if n != -1]
        F = F[keep]
        if np.any(F.sum(0) <= 0):
            raise ValueError('No assigned frequency mass remains')
        F = F / F.sum(0) * 100
        nodes = [nodes[i] for i in keep]
        if 'node' not in cells:
            raise ValueError('Cell assignments required to exclude U and report correct n')
        cells = cells.loc[cells.node.ne(-1)]
    with plt.rc_context({'font.family': 'Arial'}):
        fig = plt.figure(figsize=(15.86, 7.18))
        ax = fig.add_axes([.07, .35, .78 if grafts else .90, .45])
        wt = nodes.index(0)
        stage = [str(meta.loc[s, 'stage'])+'\nn='+format(int(cells['sample'].eq(s).sum()), ',')
                 +f'\nWT {F[wt,j]:.1f}%' for j,s in enumerate(clinical)]
        # The renderer retains an explicit zero-mass U slot for tree accounting.
        draw_clone_fish(ax, F[:, :len(clinical)] if -1 in nodes else np.vstack([F[:, :len(clinical)], np.zeros(len(clinical))]),
                        nodes if -1 in nodes else nodes+[-1], parents, colors, stage)
        ax.set_ylabel('% of all cells' if include_unassigned else '% of assigned cells (including WT)')
        ax.get_xticklabels()[-1].set_ha('right')
        fig.suptitle(patient+' | Clonal composition over time', fontsize=20, y=.98)
        for j,s in enumerate(clinical):
            if F[wt,j] >= 15:
                ax.text(j, F[wt,j]/4, f'WT {F[wt,j]:.1f}%',
                        ha='right' if j == len(clinical)-1 else 'center',
                        va='center', fontsize=10, color='#333333')
        if grafts:
            ga = fig.add_axes([.90,.35,.08,.45]); bottom = np.zeros(len(grafts))
            for i,n in enumerate(nodes):
                vals = F[i,len(clinical):]
                ga.bar(np.arange(len(grafts)), vals, bottom=bottom, color=colors[n], width=.8)
                bottom += vals
            ga.set_ylim(0,100)
            ga.set_xticks(np.arange(len(grafts)), ['MISTRG\n'+str(meta.loc[s,'stage'])
                          +f'\nWT {F[wt,len(clinical)+j]:.1f}%' for j,s in enumerate(grafts)], fontsize=8)
            ga.set_yticks([0,50,100]); ga.set_title('Matched graft', fontsize=10)
        for i,n in enumerate(nodes):
            label = 'WT = no direct or imputed mutant call' if n == 0 else names[n]+' = '+comp[n].replace('|',' | ')
            xx,yy = .035+(i%3)*.327, .275-(i//3)*.047
            fig.text(xx,yy,'■',color=colors[n],fontsize=12)
            fig.text(xx+.014,yy,textwrap.fill(label,65),fontsize=7.5,
                     fontstyle='italic' if n else 'normal',va='center')
    return fig
