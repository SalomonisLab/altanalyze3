from pathlib import Path
import importlib.util
import json
import numpy as np
import pytest

file = Path(__file__).resolve().parents[1]/'altanalyze3/components/clustering/benchmarking/benchmark_landmark_umap.py'
spec=importlib.util.spec_from_file_location('landmark_benchmark',file)
B=importlib.util.module_from_spec(spec);spec.loader.exec_module(B)


def fixture(root):
    np.save(root/'X.npy',np.ones((10,3),dtype=np.float32))
    np.save(root/'states.npy',np.array(['A']*8+['B']*2))
    (root/'cells.txt').write_text('\n'.join(f'C{i}' for i in range(10))+'\n')
    (root/'features.txt').write_text('G1\nG2\nG3\n')
    m={'baseline_verified':True, 'sha256':{f:B.digest(root/f) for f in ('X.npy','states.npy','cells.txt','features.txt')}}
    (root/'manifest.json').write_text(json.dumps(m));return m


def test_identity_hash_gate(tmp_path):
    m=fixture(tmp_path)
    x,c,g,s,_=B.load_verified_input(tmp_path)
    assert x.shape==(len(c),len(g)) and len(s)==len(c)
    (tmp_path/'features.txt').write_text('different\nG2\nG3\n')
    with pytest.raises(ValueError,match='hash'): B.load_verified_input(tmp_path)
    m['sha256']['features.txt']=B.digest(tmp_path/'features.txt')
    m['baseline_verified']=False
    (tmp_path/'manifest.json').write_text(json.dumps(m))
    with pytest.raises(ValueError,match='Baseline'): B.load_verified_input(tmp_path)


def test_landmarks_are_reproducible_and_rare_states_are_covered():
    states=np.array(['A']*1000+['B']*10+['C']*2)
    rows=B.select_landmarks(states,100,3,'state-stratified',min_per_state=20)
    assert len(rows)==len(set(rows))==100
    assert np.array_equal(rows,B.select_landmarks(states,100,3,'state-stratified',20))
    assert np.sum(states[rows]=='B')==10 and np.sum(states[rows]=='C')==2
    remaining=np.setdiff1d(np.arange(len(states)),rows)
    assert len(remaining)+len(rows)==len(states)
    expanded=B.select_landmarks(states,10,3,'state-stratified',20)
    assert len(expanded)==32
    assert np.sum(states[expanded]=='B')==10 and np.sum(states[expanded]=='C')==2


def test_neighbor_checks_do_not_assume_self_is_first():
    coords=np.array([[0.,0.],[0.,0.],[1.,0.],[2.,0.]])
    anchors=np.arange(4)
    rows=B.coordinate_neighbors(coords,anchors,2)
    assert rows.shape==(4,2)
    assert all(i not in n for i,n in enumerate(rows))


def test_missing_roster_hash_rejects_execution(tmp_path):
    m=fixture(tmp_path)
    del m['sha256']['features.txt']
    (tmp_path/'manifest.json').write_text(json.dumps(m))
    with pytest.raises(ValueError,match='Hashes'): B.load_verified_input(tmp_path)
