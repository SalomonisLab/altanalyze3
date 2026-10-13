"""Every isoform label stays readable when most reads occupy one cluster."""
import pytest
from altanalyze3.components.visualization.isv_web.isv_web_remote import read_label_positions

@pytest.mark.parametrize('n,height,gap', [(48,480,10),(47,470,10),(9,90,10),(2,100,10)])
def test_dense_tail_keeps_all_labels_ordered_and_inside(n,height,gap):
    centers=[height * .65]+[height * .8 + i*.1 for i in range(n-1)]
    pos=read_label_positions(centers,height,gap)
    assert len(pos)==n
    assert pos[0]>=gap/2 - 1e-9 and pos[-1]<=height-gap/2 + 1e-9
    assert all(b-a >= gap-1e-9 for a,b in zip(pos,pos[1:]))

def test_readable_centers_unchanged():
    assert read_label_positions([20,40,60],100,10)==[20,40,60]

def test_inadequate_panel_fails_without_omitting_labels():
    with pytest.raises(ValueError): read_label_positions([1,2,3],29,10)


def test_empty_labels():
    assert read_label_positions([],0,10)==[]
