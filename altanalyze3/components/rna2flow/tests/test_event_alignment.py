import numpy as np,pandas as pd,pytest
from altanalyze3.components.rna2flow.io import FlowData,align_event_labels,read_matrix_csv


def test_shuffled_label_rows_join_to_event_ids():
    f=FlowData(np.zeros((3,1)),['CD4'],'input',{'event_ids':np.array([10,20,30])})
    d=pd.DataFrame({'EventNumberDP':[30,10,20],'FlowSOM':['C','A','B']})
    assert align_event_labels(f,d).tolist()==['A','B','C']
    with pytest.raises(ValueError,match='missing'):align_event_labels(f,d.iloc[:2])
    with pytest.raises(ValueError,match='Duplicate'):align_event_labels(f,pd.concat([d,d.iloc[:1]]))


def test_csv_identifier_is_not_an_antibody(tmp_path):
    p=tmp_path/'flow.csv';pd.DataFrame({'EventNumberDP':[2,1],'CD4':[4,5],'CD8':[9,2]}).to_csv(p,index=False)
    f=read_matrix_csv(str(p));assert f.channels==['CD4','CD8'];assert f.meta['event_ids'].tolist()==[2,1]
