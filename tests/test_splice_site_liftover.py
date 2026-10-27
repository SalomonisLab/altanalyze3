from altanalyze3.components.psi.splice_site_liftover import SpliceSiteLiftover


def test_actual_native_point_liftover_preserves_one_based_coordinates_and_chain(tmp_path):
    p=tmp_path/'synthetic.chain'
    p.write_text('chain 100 chrM 1000 + 99 199 chrH 1000 + 199 299 1\n100\n\n'
                 'chain 200 chrN 1000 + 99 199 chrH 1000 - 799 899 2\n100\n\n')
    lo=SpliceSiteLiftover(p)
    a=lo.junction('chrM:100-110')['same_chain_junction_mappings']
    assert len(a)==1 and (a[0]['first'],a[0]['second'],a[0]['chain_id'])==(200,210,'1')
    b=lo.junction('chrN:100-110')['same_chain_junction_mappings']
    assert len(b)==1 and (b[0]['first'],b[0]['second'],b[0]['chain_orientation'])==(201,191,'-')
    assert not lo.junction('chrM:100-210')['same_chain_junction_mappings']


def test_distinct_endpoint_chains_cannot_form_a_junction(tmp_path):
    p=tmp_path/'disjoint.chain'
    p.write_text('chain 100 chrM 1000 + 99 109 chrH 1000 + 199 209 1\n10\n\n'
                 'chain 100 chrM 1000 + 199 209 chrH 1000 + 299 309 2\n10\n\n')
    r=SpliceSiteLiftover(p).junction('chrM:100-200')
    assert len(r['first_point_mappings'])==len(r['second_point_mappings'])==1
    assert not r['same_chain_junction_mappings']
