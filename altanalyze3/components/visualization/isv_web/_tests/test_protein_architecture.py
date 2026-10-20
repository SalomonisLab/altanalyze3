"""Optional native protein tracks preserve the ISV read data and defaults."""
from copy import deepcopy
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import pytest
from altanalyze3.components.iso2function.tests.test_isv_architecture import bundle
from altanalyze3.components.visualization.isv_web import isv_web_remote as remote


def response():
    rows=[]
    for i,(iid,n,end) in enumerate([('A',200,699),('B',100,399)]):
        exons=[dict(start=100,end=end,type='E',exon_id='E1.1',chrom='chr1')]
        rows.append(dict(isoform_id='read_'+iid,final_isoform_id=iid,known=True,count=1,
                         protein_length=n,nmd_status='Not-NMD',cluster_index=i,color_index=i,
                         exon_segments=exons,registered_exons=[[100,end]],cds_min=100,cds_max=end))
    return dict(gene='ENSG_G',symbol='G',gene_model=dict(chrom='chr1',strand='-',extent=[100,699],
                blocks=[dict(start=100,end=699,exon_id='E1.1')]),
                panels=[dict(condition='HSC-MPP',n_molecules=2,molecules=rows)])


def test_renderer_preserves_reads_order_coordinates_and_protein_sidecar(tmp_path):
    source=tmp_path/'source';source.mkdir();bundle(source)
    res=response();before=deepcopy(res);out=tmp_path/'plot.pdf'
    assert remote.render(res,str(out),isoform_track=True,protein_features=source,
                         show_read_labels=False,structure_labels_above=True,title_gene='G')
    assert res==before
    audit=json.loads(out.with_suffix('.protein_tracks.json').read_text())
    assert audit['RNA']['isoform_order']==['A','B']
    assert audit['RNA']['panels'][0]['reads']==2
    assert [t['protein_aa'] for t in audit['tracks']]==[200,100]
    import pymupdf
    with pymupdf.open(out) as doc:
        assert len(doc)==1 and not doc[0].get_images()
        assert 'DNA binding (partial)' in doc[0].get_text()


def test_default_renderer_has_no_protein_tracks_or_sidecar(tmp_path):
    out=tmp_path/'default.pdf';remote.render(response(),str(out),isoform_track=True)
    assert not out.with_suffix('.protein_tracks.json').exists()


def test_missing_native_model_and_conflicting_lengths_fail_before_output(tmp_path):
    source=tmp_path/'source';source.mkdir();bundle(source)
    res=response();res['panels'][0]['molecules'][1]['protein_length']=99
    with pytest.raises(ValueError,match='length differs'):remote.render(res,str(tmp_path/'bad.pdf'),isoform_track=True,protein_features=source)
    assert not (tmp_path/'bad.pdf').exists()
    res=response();res['panels'][0]['molecules'][1].pop('registered_exons')
    with pytest.raises(ValueError,match='registered native'):remote.render(res,str(tmp_path/'bad.pdf'),isoform_track=True,protein_features=source)


def test_saved_input_cli_does_not_load_run_context_or_allow_hidden_query_changes(tmp_path):
    source=tmp_path/'saved.json.gz'
    with gzip.open(source,'wt') as f:json.dump(response(),f)
    args=['--run_dir',str(tmp_path),'--gene','G','--render-input',str(source),'--out',str(tmp_path/'cache.pdf')]
    with patch.object(remote,'render',return_value=True) as draw:
        with pytest.raises(SystemExit) as done:remote.main(args)
        assert done.value.code==0 and draw.call_args.args[0]==response()
        assert draw.call_args.kwargs['protein_features'] is None
    with pytest.raises(SystemExit,match='query/filter options'):remote.main(args+['--conditions','young'])


def test_web_option_retains_rows_and_reports_missing_architecture(tmp_path,monkeypatch):
    from fastapi.testclient import TestClient
    from altanalyze3.components.visualization.isv_web import server
    ctx=SimpleNamespace(resolve_gene=lambda g:'ENSG_G',molecule_grouping='final_isoform')
    source=response();before=deepcopy(source)
    monkeypatch.setattr(server.da,'query_reads',lambda *a,**kw:source)
    client=TestClient(server.create_app(ctx))
    assert client.post('/api/reads',json=dict(gene='G')).json()==source
    r=client.post('/api/reads',json=dict(gene='G',protein_features=True));assert r.status_code==422
    bdir=tmp_path/'source';bdir.mkdir();ctx.protein_architectures=bundle(bdir)
    r=client.post('/api/reads',json=dict(gene='G',protein_features=True));assert r.status_code==200,r.text
    data=r.json();assert source==before
    assert len(data['protein_architecture']['tracks'])==2
    for p,old in zip(data['panels'],before['panels']):
        for row,orig in zip(p['molecules'],old['molecules']):
            assert {k:v for k,v in row.items() if k!='display_color'}==orig
    ctx.protein_architectures['models'].pop('B');client=TestClient(server.create_app(ctx))
    assert client.post('/api/reads',json=dict(gene='G',protein_features=True)).status_code==422
