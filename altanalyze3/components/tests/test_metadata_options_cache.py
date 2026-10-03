import anndata as ad
import numpy as np
import pandas as pd
from altanalyze3.components.cellHarmony.webapp.metadata_options import cached_group_fields, _CACHE
from altanalyze3.components.cellHarmony.flask import pipeline


def test_group_menu_cache_preserves_values_and_avoids_rereads(tmp_path,monkeypatch):
    path=tmp_path/'query.h5ad'
    obj=ad.AnnData(np.zeros((8,3)),obs=pd.DataFrame({'Library':['a','b']*4,'condition':['case','control']*4},index=[f'c{i}' for i in range(8)]))
    obj.write_h5ad(path)
    reader=pipeline._read_h5ad_obs
    calls=[]
    def tracked(path):
        calls.append(path)
        return reader(path)
    monkeypatch.setattr(pipeline,'_read_h5ad_obs',tracked)
    expected=pipeline._candidate_group_fields(path,preferred=['Library'],max_categories=None)
    calls.clear()
    got=cached_group_fields(path,preferred=['Library'])
    assert got==expected
    got[1]['Library'].append('mutated')
    assert cached_group_fields(path,preferred=['Library'])==expected
    assert len(calls)==1
    assert _CACHE.budget.snapshot()['bytes']<=1024**2
    obj.obs['Library']=['new','other']*4
    obj.write_h5ad(path)
    assert cached_group_fields(path,preferred=['Library'])[1]['Library']==['new','other']
    assert len(calls)==2


def test_different_category_limits_and_missing_files(tmp_path):
    path=tmp_path/'query.h5ad'
    obj=ad.AnnData(np.zeros((8,2)),obs=pd.DataFrame({'donor':[f'd{i}' for i in range(8)]},index=[f'c{i}' for i in range(8)]))
    obj.write_h5ad(path)
    assert cached_group_fields(path,max_categories=3)==([], {})
    assert len(cached_group_fields(path,max_categories=None)[1]['donor'])==8
    path.unlink()
    assert cached_group_fields(path)==([], {})
