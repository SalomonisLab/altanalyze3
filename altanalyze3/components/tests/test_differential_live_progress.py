"""Native DE reports real stages to the web store, including zero-DE states."""
import anndata as ad
import numpy as np
import pandas as pd

from altanalyze3.components.cellHarmony import cellHarmony_differential as D
from altanalyze3.components.cellHarmony.flask.pipeline import _differential_de_progress
from altanalyze3.components.cellHarmony.flask.job_manager import JobStore


def test_web_grn_progress_covers_no_de_and_insufficient_states(tmp_path, monkeypatch):
    monkeypatch.setattr(D, 'diagnostic_report', False)
    store=JobStore(tmp_path/'jobs')
    job=store.create_job('human','test-reference',None,[])['job_id']
    store.update_job(job,differential={'status':'processing','progress':35,'config':{'modality':'grn'}})
    rng=np.random.default_rng(98)
    # Matched distributions give exactly no differential edges in A and B.
    a=rng.uniform(.1,2.,(24,80)).astype(np.float32)
    x=np.vstack([a,a,a,a,a[:1],a[:1]])
    obs=pd.DataFrame({'state':['A']*48+['B']*48+['C']*2,
        'condition':(['case']*24+['control']*24)*2+['case','control']},index=[f'c{i}' for i in range(len(x))])
    matrix=ad.AnnData(x,obs=obs,var=pd.DataFrame(index=[f'TF->g{i}' for i in range(80)]),
                     uns={'modality':'grn','expression_scale':'linear'})
    events=[];snapshots=[];ticks=[]
    update=_differential_de_progress(store,job,'grn')
    def report(**event):
        events.append(event)
        update(**event)
        snapshots.append(store.get_job(job)['differential'])
    result=D.run_de_for_comparisons(matrix,'state','condition','case','control','wilcoxon',.05,1.1,4,
                                   progress_callback=lambda:ticks.append(True),stage_callback=report)
    assert result['detailed_deg'].empty
    summary=result['summary_per_population'].set_index('population')
    assert summary.index.tolist()==['A','B','C']
    assert summary.num_DEG.eq(0).all()
    assert summary.loc['A','tested_genes']==80 and summary.loc['C','tested_genes']==0
    assert [e['label'] for e in events if e['phase']=='population']==['A','B','C']
    assert {'population','pooled','patterns','assembling'} <= {e['phase'] for e in events}
    assert len(ticks)==4  # Every population finishes, plus legacy final completion.
    assert [s['progress'] for s in snapshots]==sorted(s['progress'] for s in snapshots)
    assert all(s['config']=={'modality':'grn'} and s['status']=='processing' for s in snapshots)
    assert any('cell state 2 of 3: B' in s['message'] for s in snapshots)
    assert 'assembling' in snapshots[-1]['message'] and snapshots[-1]['progress']==59
