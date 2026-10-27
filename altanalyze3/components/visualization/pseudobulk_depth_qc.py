"""Descriptive pseudobulk depth diagnostics; no differential testing or data replacement."""
from pathlib import Path
import numpy as np
from scipy.stats import hypergeom


def expected_detected_features(counts, target_reads):
    """Expected distinct features in a without-replacement draw from one pseudobulk.

    Sparse zero counts contribute zero. The original feature denominator and count
    vector are unchanged. A target above the available read count is unavailable,
    rather than an extrapolation or a normalized synthetic measurement.
    """
    c=np.asarray(counts)
    if c.ndim!=1 or not np.isfinite(c).all() or (c<0).any() or not np.equal(c,np.round(c)).all():
        raise ValueError('A finite nonnegative integer count vector is required')
    if not isinstance(target_reads,(int,np.integer)) or target_reads<0:
        raise ValueError('target_reads must be a nonnegative integer')
    total=int(c.sum())
    if target_reads>total:return float('nan')
    if not total or not target_reads:return 0.
    nonzero=c[c>0].astype(np.int64)
    return float((-np.expm1(hypergeom.logpmf(0,total,nonzero,target_reads))).sum())


def plot_depth_and_cell_numbers(table, output, *, state, caption, target_reads=7655):
    """Show all donor points, measured cells/reads and a depth-standardized expectation."""
    import matplotlib.pyplot as plt
    import textwrap
    d=table[table.state==state].copy()
    if d.uid.duplicated().any():raise ValueError('Expected one source pseudobulk per donor and state')
    d['display_group']=np.where(d.condition=='young','Young',np.where(d.condition=='aged','Aged',
        np.where(d.condition=='MDS-pre-SRSF2','MDS-pre SRSF2','MDS-pre other/UNK')))
    order=['Young','Aged','MDS-pre other/UNK','MDS-pre SRSF2']
    colors=['#4989bb','#777777','#8765aa','#cc4a49']
    fields=[('cells_with_long_read_input','Cells assigned to '+state),
            ('assigned_isoform_reads','Assigned isoform reads'),
            ('expected_detected_isoforms_at_'+str(target_reads),'Expected detected isoforms\nat '+format(target_reads,',')+' reads')]
    fig,axes=plt.subplots(1,3,figsize=(11.2,4.9));fig.subplots_adjust(bottom=.40,wspace=.4,top=.85)
    for ax,(field,ylabel) in zip(axes,fields):
        for j,(group,color) in enumerate(zip(order,colors)):
            part=d[d.display_group==group];v=part[field].to_numpy(float)
            offsets=np.linspace(-.13,.13,len(part))
            ax.scatter(j+offsets,v,color=color,s=24,zorder=3)
            for offset,(_,r) in zip(offsets,part.iterrows()):
                if np.isfinite(r[field]):ax.annotate(str(r.uid),(j+offset,r[field]),xytext=(3,2),textcoords='offset points',fontsize=4.8,rotation=45)
        ax.set_xticks(range(4),[g+'\nn='+str((d.display_group==g).sum()) for g in order],rotation=35,ha='right',fontsize=7)
        ax.set_ylabel(ylabel,fontsize=8);ax.tick_params(axis='y',labelsize=7)
        if field in ['cells_with_long_read_input','assigned_isoform_reads']:ax.set_yscale('log')
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle(state+': cell representation and detection-depth confounding',fontsize=12)
    fig.text(.04,.025,'\n'.join(textwrap.wrap(caption,width=170)),fontsize=7,va='bottom')
    Path(output).parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output);plt.close(fig)
