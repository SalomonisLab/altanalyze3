"""Explicit strict top-N review cannot add listed, but unqualified, source events."""
import hashlib,json
import anndata as ad
import numpy as np,pandas as pd
from altanalyze3.components.cellHarmony.iterative_marker_finder import review_legacy_markerfinder_fit


def test_strict_review_keeps_qualifying_winners_and_enforces_cache_cap(tmp_path):
    strong='ENSG1:E1.1-E2.1|ENSG1:E1.1-E3.1'
    weak='ENSG2:E1.1-E2.1|ENSG2:E1.1-E3.1'
    a=ad.AnnData(np.array([[.1,.2],[.15,.4],[.2,.3],[.8,.2],[.85,.4],[.9,.3]]),
        obs=pd.DataFrame({'group':['C']*3+['P']*3,'donor':list('abcdef')},index=list('abcdef')),
        var=pd.DataFrame(index=[strong,weak]))
    a.uns['lineage_order']=['C','P'];path=tmp_path/'source.h5ad';a.write_h5ad(path)
    fit=tmp_path/'fit';fit.mkdir()
    (fit/'parameters.json').write_text(json.dumps(dict(source_features=a.var_names.tolist(),
        rho_strictly_greater_than=.4,schemes=[dict(source_samples=a.obs_names.tolist(),scored_samples=a.obs_names.tolist())])))
    pd.DataFrame(dict(pattern=['up-P']*2,scheme=['clinical']*2,rho=[.95,.5]),index=a.var_names).to_csv(fit/'Consolidated-MarkerFinder.assignments.tsv',sep='\t')
    labels=tmp_path/'labels.tsv';labels.write_text('UID\n'+weak+'\n')
    out=review_legacy_markerfinder_fit(dict(display_h5ad=str(path),display_group_key='group',label_events=str(labels)),
        fit,tmp_path/'review',input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        rho_cutoff=.8,top_n=1,strict_top_n=True,collapse_pattern_logic=False,dpi=100)
    assert out['plotted']==1 and out['unassigned_callout_references']==0
    cache=ad.read_h5ad(out['heatmap_cache'])
    assert cache.var.gene.tolist()==[strong] and cache.var.cluster.tolist()==['up-P']
    assert np.array_equal(a.X,ad.read_h5ad(path).X)
