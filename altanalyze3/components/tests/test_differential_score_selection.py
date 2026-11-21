"""Native-score selection keeps identities/caps and survives row permutations."""
import json
import numpy as np
import pandas as pd
import anndata as ad
import pytest
from altanalyze3.components.cellHarmony.differential_sampling import select_cells, reference_scores, web_assignment_scores


def fixture():
    return pd.DataFrame({'state':['A']*4+['B']*2+['A'],
        'sample':['s1']*4+['s2']*2+['outside'],
        'AlignmentScore':[.1,.8,.8,-.5,.3,np.nan,1.]},index=['z','b','a','c','d','e','f'])


def test_highest_native_scores_and_barcode_ties_preserve_every_small_stratum():
    obs=fixture()
    picked=select_cells(obs,'state','sample',['s1','s2'],limit=2)
    assert obs.index[picked].tolist()==['a','b','d','e']
    assert obs.iloc[picked].groupby(['state','sample']).size().to_dict()=={('A','s1'):2,('B','s2'):2}
    obs.loc['z','AlignmentScore']=.9
    assert obs.index[select_cells(obs,'state','sample',['s1'],limit=2)].tolist()==['a','z']


def test_same_barcodes_and_analysis_order_after_reordering_cells_or_samples():
    obs=fixture();expected=obs.index[select_cells(obs,'state','sample',['s1','s2'],limit=2)]
    shuffled=obs.iloc[[6,5,3,2,0,4,1]]
    actual=shuffled.index[select_cells(shuffled,'state','sample',['s2','s1'],limit=2)]
    assert actual.equals(expected)
    scores=obs['AlignmentScore'].iloc[::-1]
    assert shuffled.index[select_cells(shuffled,'state','sample',['s2','s1'],limit=2,scores=scores)].equals(expected)


@pytest.mark.parametrize('missing',[np.nan,np.inf,-np.inf,'bad'])
def test_nonfinite_or_unmatched_capped_scores_stop_selection_without_dropping_cells(missing):
    obs=fixture();scores=obs.AlignmentScore.astype(object);scores.loc['a']=missing
    with pytest.raises(ValueError,match='missing or nonfinite'):
        select_cells(obs,'state','sample',['s1'],limit=2,scores=scores)
    assert len(obs)==7


def test_missing_score_source_does_not_fall_back_to_random_selection():
    obs=fixture().drop(columns='AlignmentScore')
    with pytest.raises(ValueError,match='No random selection'):
        select_cells(obs,'state','sample',['s1'],limit=2)
    assert obs.index[select_cells(obs,'state','sample',['s1'],limit=4)].tolist()==['a','b','c','z']


def test_both_scores_use_native_branch_assignments_and_below_cutoff_cells(tmp_path):
    job='bothfixture';output=tmp_path/'outputs';branch=output/'branches'/'supervised'/job
    branch.mkdir(parents=True)
    assignments=branch/'native_assignments.txt'
    pd.DataFrame({'CellBarcode':['c3','c1'],'AlignmentScore':[.8,.9]}).to_csv(assignments,sep='\t',index=False)
    below=branch/'below.h5ad'
    ad.AnnData(np.ones((1,1)),obs=pd.DataFrame({'AlignmentScore':[.1]},index=['c2'])).write_h5ad(below)
    (branch/'job.json').write_text(json.dumps({'artifacts':{'assignments':str(assignments),'unaligned_cells_h5ad':str(below)}}))
    meta={'analysis_mode':'both','job_id':job,'artifacts':{'combined_h5ad':str(output/'combined.h5ad'),'assignments':'labels_only.tsv'}}
    scores,sources=reference_scores(meta,['c2','c3','c1'])
    np.testing.assert_array_equal(scores,[.1,.8,.9])
    assert sources==[str(assignments),str(below)]


def test_duplicate_native_identities_are_not_averaged_or_dropped(tmp_path):
    assignments=tmp_path/'assignments.txt'
    assignments.write_text('CellBarcode\tAlignmentScore\nc1\t0.3\nc1\t0.9\n')
    with pytest.raises(ValueError,match='duplicate cell barcodes'):
        reference_scores({'artifacts':{'assignments':str(assignments)}},['c1'])


@pytest.mark.parametrize('column', ['SVM_score', 'ICGS3_SVM_score'])
def test_unsupervised_uses_saved_winning_score_by_barcode_without_reading_expression(tmp_path, monkeypatch, column):
    import altanalyze3.components.cellHarmony.differential_sampling as sampling
    obs=fixture().drop(columns='AlignmentScore')
    # Margins deliberately disagree: the winning score is the authorized metric.
    obs[column]=[.1,.8,.8,-.5,.3,np.nan,1.]
    obs['SVM_margin']=[99,0,0,99,0,0,0]
    path=tmp_path/'native.h5ad'
    ad.AnnData(np.ones((len(obs),2)),obs=obs).write_h5ad(path)
    real_read=sampling.read_elem
    def obs_only(elem):
        assert elem.name.startswith('/obs/')
        return real_read(elem)
    monkeypatch.setattr(sampling, 'read_elem', obs_only)
    # A modality sidecar can use a different row order without carrying scores.
    sidecar=obs.iloc[::-1].drop(columns=[column,'SVM_margin'])
    scores,sources,method,metric=web_assignment_scores(
        {'qc':{'analysis_mode':'unsupervised'},'artifacts':{'combined_h5ad':str(path)}},sidecar.index)
    rows=select_cells(sidecar,'state','sample',['s1','s2'],limit=2,scores=scores)
    assert sidecar.index[rows].tolist()==['a','b','d','e']
    assert sources==[str(path)] and method=='highest_icgs3_svm_score' and metric==column


@pytest.mark.parametrize('mode', ['supervised', 'both'])
def test_web_modes_keep_reference_policy_even_for_unsupervised_population_field(monkeypatch, mode):
    import altanalyze3.components.cellHarmony.differential_sampling as sampling
    expected=pd.Series([.9],index=['cell'])
    monkeypatch.setattr(sampling,'reference_scores',lambda meta,ids:(expected,['reference.tsv']))
    scores,sources,method,metric=web_assignment_scores({'analysis_mode':mode},['cell'])
    assert scores is expected and sources==['reference.tsv']
    assert (method,metric)==('highest_alignment_score','AlignmentScore')


@pytest.mark.parametrize('problem', ['missing', 'conflicting', 'duplicate'])
def test_unsupervised_missing_conflicting_or_duplicate_scores_never_substitute(tmp_path, problem):
    obs=pd.DataFrame({'SVM_margin':[9.,8.]},index=['a','b'])
    if problem=='conflicting':
        obs['SVM_score']=[.1,.2];obs['ICGS3_SVM_score']=[.2,.1]
    if problem=='duplicate':
        obs['SVM_score']=[.1,.2];obs.index=['a','a']
    path=tmp_path/'native.h5ad'
    ad.AnnData(np.ones((2,1)),obs=obs).write_h5ad(path)
    with pytest.raises(ValueError,match='missing|conflict|unique'):
        web_assignment_scores({'analysis_mode':'unsupervised','artifacts':{'combined_h5ad':str(path)}},['a','b'])
