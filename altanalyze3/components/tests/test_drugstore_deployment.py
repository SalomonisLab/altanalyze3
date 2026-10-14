"""Deployment checks retain the complete native signature and gene rosters."""
import json
import numpy as np
import pytest
from altanalyze3.components.cellHarmony.discover_resources.verify_drugstore import verify
from altanalyze3.components.cellHarmony.discover_resources import drug_matching as dm


@pytest.fixture
def bundle(tmp_path):
    folder=tmp_path/'partitions'/'fixture'
    folder.mkdir(parents=True)
    (folder/'genes.json').write_text(json.dumps(['A','B']))
    (folder/'keys.json').write_text(json.dumps(['fixture:s0','fixture:s1']))
    records=[(json.dumps({'signature_key':key})+'\n').encode() for key in ['fixture:s0','fixture:s1']]
    (folder/'records.jsonl').write_bytes(b''.join(records))
    np.save(folder/'record_offsets.npy',np.array([0,len(records[0]),sum(map(len,records))]))
    np.save(folder/'values.npy',np.array([[1.,2.],[3.,4.]],dtype=np.float16))
    np.save(folder/'neglog10p.npy',np.ones((2,2),dtype=np.float16))
    np.save(folder/'is_landmark.npy',np.array([True,False]))
    (tmp_path/'manifest.json').write_text(json.dumps({'complete':True,'method_contract':dm.METHOD_VERSION,
        'store_version':'fixture','n_signatures':2,'local_effect_aggregation':'none','partitions':[{
            'partition_id':'fixture','n_genes':2,'n_signatures':2,'pvalue_available':True,
            'landmark_mask':'is_landmark.npy','n_landmark':1,'signature_type':'drug treatment'}]}))
    return tmp_path


def test_full_bundle_checks_and_env_path(bundle,monkeypatch):
    monkeypatch.setenv('CELLREF2_DRUGSTORE',str(bundle))
    result=verify()
    assert result['status']=='ready' and result['n_signatures']==2
    assert result['partitions'][0]['genes']==2 and result['partitions'][0]['shape']==[2,2]
    assert len(result['manifest_sha256'])==64


@pytest.mark.parametrize('problem', ['partial','count','mask','pvalues','offsets','duplicate_keys'])
def test_incomplete_bundle_never_subsets_or_silently_uses_other_sources(bundle,problem):
    folder=bundle/'partitions'/'fixture'
    manifest=json.loads((bundle/'manifest.json').read_text())
    if problem=='partial':manifest['complete']=False
    elif problem=='count':manifest['n_signatures']=3
    elif problem=='mask':(folder/'is_landmark.npy').unlink()
    elif problem=='pvalues':(folder/'neglog10p.npy').unlink()
    elif problem=='offsets':np.save(folder/'record_offsets.npy',np.array([0,1]))
    else:(folder/'keys.json').write_text(json.dumps(['fixture:s0','fixture:s0']))
    (bundle/'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError):verify(bundle)
