"""Failure-mode tests for reference transfer, rather than assertions mirroring a report."""
import copy
import json
from types import SimpleNamespace

import pytest

from altanalyze3.components.annotation.domain_impact import (
    legacy_feature_fragment, feature_sequence_present, Domain, domain_coverage)
from altanalyze3.components.iso2function.protein_impact import (
    comparison, native_models, pinned_bundle, project_instance, sha, read_support)
from altanalyze3.components.iso2function.protein_isoform_audit import project_feature


def feature(ref="MACDEFGHIKLMNPQRSTVWY", a=3, b=12):
    return dict(reference_sequence=ref, start=a, end=b)


def test_invalid_reference_boundaries_do_not_become_intact_empty_strings():
    for a, b in [(0, 3), (3, 100), (8, 7), (50, 60)]:
        assert project_instance(feature(a=a, b=b), "MACDEF")["status"] == "invalid_reference_interval"
        with pytest.raises(ValueError, match="1-based"):
            project_feature("MACDEF", "MACDEF", a, b, {})
    with pytest.raises(ValueError, match="empty"):
        feature_sequence_present("", "MACDEF")


def test_unresolved_references_are_not_losses():
    p = project_instance(feature(ref=""), "MACDEF")
    assert p["status"] == "unresolved_reference_sequence"
    assert comparison(True, p["sequence_present"]) == "unresolved"


def test_repeats_are_present_but_their_placement_is_ambiguous():
    p = project_instance(feature("MACDEF", 2, 6), "MACDEFACDEF")
    assert p["sequence_present"] is True
    assert p["status"] == "ambiguous_repeated_sequence"
    assert comparison(True, p["sequence_present"]) == "shared_present"


def test_partial_domain_is_not_an_intact_domain_or_an_entirely_absent_protein():
    p = project_instance(feature(), "MACDEFGH")
    assert p["sequence_present"] is False
    assert p["status"] == "altered_or_partial"
    assert 0 < p["exact_fraction"] < 1


def test_both_short_forms_already_missing_feature_is_not_a_pair_loss():
    a = project_instance(feature(a=13, b=20), "MACDEF")
    b = project_instance(feature(a=13, b=20), "MACDEFG")
    assert comparison(a["sequence_present"], b["sequence_present"]) == "absent_from_both"


def test_legacy_baseline_preserves_historical_slice_and_short_feature_expansion():
    seq = "MACDEFGHIKLMNPQRSTVWY"
    # Frozen parity with FeatureAlignment.importEnsemblFTdata, including its indexing.
    assert legacy_feature_fragment(seq, 3, 12) == "DEFGHIKLM"
    assert legacy_feature_fragment(seq, 5, 10) == "EFGHIKL"
    assert legacy_feature_fragment(seq, 6, 6) == "DEFGHI"
    # Correct inclusive boundaries and historical boundaries can differ: expose both.
    p = project_instance(feature(seq, 3, 12), "CDEFGHIKLM")
    assert p["sequence_present"] is True
    assert p["legacy_sequence_present"] is True
    p = project_instance(feature(seq, 3, 12), "DEFGHIKLM")
    assert p["sequence_present"] is False
    assert p["legacy_sequence_present"] is True


def test_genomic_coverage_cannot_establish_feature_sequence_identity():
    d = Domain("P", "IPR1", "D", "D", 1, 6, 100, 117, ((100, 117),))
    assert domain_coverage([(100, 117)], d) == 1
    assert project_instance(feature("MACDEF", 1, 6), "MYYYYY")["sequence_present"] is False


def test_source_roster_contract_preserves_order_and_rejects_substitution(tmp_path):
    source = tmp_path / "native.tsv"
    source.write_text("gene_symbol\tisoform_id\tensg\tprotein_sequence\tcds_sequence\tcdna_sequence\nG\tnew.id\tENSG1\tMA\tATGGCT\tCCATGGCTTAA\nG\tOTHER\tENSG1\tMA\tATGGCT\tATGGCT\n")
    roster = tmp_path / "roster.json"
    c = dict(source_table_sha256=sha(source), required_roster=[dict(gene_symbol="G", isoform_id="new.id")])
    roster.write_text(json.dumps(c))
    rows, audit = native_models(source, roster)
    assert [r["isoform_id"] for r in rows] == ["new.id"]
    assert len(audit) == 2 and audit[1]["scope"] == "outside_independently_defined_roster"
    c["required_roster"][0]["isoform_id"] = "missing"
    roster.write_text(json.dumps(c))
    with pytest.raises(ValueError, match="Incomplete"):
        native_models(source, roster)
    source.write_text(source.read_text() + "\n")
    with pytest.raises(ValueError, match="hash"):
        native_models(source, roster)


def test_reference_release_and_content_are_pinned(tmp_path):
    data = tmp_path / "ref.txt"; data.write_text("reference")
    b = dict(schema="altanalyze-native-protein-impact-v1", release_decision="User chose test release",
             resources={k: dict(path="ref.txt", sha256=sha(data), source="test fixture", release="test") for k in
                        ["protein_relationships", "protein_features", "uniprot_features", "protein_coordinates", "protein_sequences", "ensembl_peptides"]})
    path = tmp_path / "bundle.json"; path.write_text(json.dumps(b))
    pinned_bundle(path)
    data.write_text("changed release")
    with pytest.raises(ValueError, match="hash"):
        pinned_bundle(path)


def test_reads_retain_literal_dotted_ids_and_zero_read_models(tmp_path):
    import h5py
    import numpy as np
    hpath = tmp_path / "counts.h5ad"
    with h5py.File(hpath, "w") as h:
        for axis, names in [("obs", ["HSC.L-isoform", "MPP.L-isoform"]),
                            ("var", ["ENSG1:a.source001", "ENSG1:b.source001"])]:
            g = h.create_group(axis); g.attrs["_index"] = "_index"
            g.create_dataset("_index", data=np.array(names, dtype=h5py.string_dtype()))
        h.create_dataset("X", data=[[10, 0], [20, 0]])
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("library\tgroups\toriginal_library\nL\tyoung\tL\n")
    models = [dict(ensg="ENSG1", isoform_id=v) for v in ["a.source001", "b.source001"]]
    counts, _ = read_support(hpath, metadata, ["HSC", "MPP"], ["young"], None, models)
    assert counts == {"ENSG1:a.source001": 30, "ENSG1:b.source001": 0}
    with pytest.raises(ValueError, match="Incomplete"):
        read_support(hpath, metadata, ["HSC", "ABSENT"], ["young"], None, models)
    models.append(dict(ensg="ENSG1", isoform_id="missing"))
    with pytest.raises(ValueError, match="missing from read matrix"):
        read_support(hpath, metadata, ["HSC"], ["young"], None, models)


def test_complete_cli_handler_repeats_byte_for_byte_and_fills_same_release_sequence_gap(tmp_path):
    from altanalyze3.components.iso2function.protein_impact import run_protein_impact
    gene = 'ENSG00000000001'
    texts = {
        'protein_relationships': 'Gene\tTrans\tProtein\n'+gene+'\tENST1\tENSP1\n',
        'protein_features': 'ID\tAA_Start\tAA_Stop\tStart\tStop\tName\tInterpro\tDescription\nENSP1\t2\t3\t4\t9\tD\tIPR1\tDomain\n',
        'uniprot_features': 'ensembl_prot\taa_start\taa_stop\tgenomic_start\tgenomic_stop\tname\tinterpro_id\tdescription\nENSP1\t2\t3\t4\t9\t\tUniProt\tMOTIF-example\n',
        'protein_coordinates': 'protienID\texonID\tAA_NT_Start\tAA_NT_Stop\tGenomic_Start\tGenomic_Stop\nENSP1\tE1\t1\t9\t1\t9\n',
        'protein_sequences': 'OTHER\tMAL\n',  # Distributed sequence file lacks the required ENSP.
        'feature_peptides': '>ENSP1.1 gene:'+gene+' transcript:ENST1\nMAL\n',
        'ensembl_peptides': '>ENSP1.1 gene:'+gene+' transcript:ENST1\nMAL\n',
    }
    resources = {}
    for key, text in texts.items():
        p = tmp_path / (key+'.txt'); p.write_text(text)
        resources[key] = dict(path=p.name, sha256=sha(p), source='synthetic test fixture', release='same-release-test')
    bundle = tmp_path / 'bundle.json'
    bundle.write_text(json.dumps(dict(schema='altanalyze-native-protein-impact-v1',release_decision='test fixture',resources=resources)))
    source = tmp_path/'source.tsv'
    source.write_text('gene_symbol\tisoform_id\tensg\tprotein_sequence\tcds_sequence\tcdna_sequence\tprotein_length\tannotated\tnmd_status\nG\ta\t'+gene+'\tMAL\tATGGCTCTG\tATGGCTCTG\t3\tyes\tNot-NMD\nG\tb.source001\t'+gene+'\tMA\tATGGCT\tATGGCT\t2\tno\tNot-NMD\n')
    roster=tmp_path/'roster.json'
    roster.write_text(json.dumps(dict(source_table_sha256=sha(source),required_roster=[dict(gene_symbol='G',isoform_id=i) for i in ['a','b.source001']])))
    args=SimpleNamespace(native_models=str(source),required_roster=str(roster),reference_bundle=str(bundle),
                         artifacts_dir=str(tmp_path/'one'),pseudobulk=None,metadata=None,state=None,
                         condition=None,donor=None,gff=None,coding_regions=None,min_isoform_reads=30,prior_feature_audit=None)
    assert run_protein_impact(args)==0
    args.artifacts_dir=str(tmp_path/'two')
    assert run_protein_impact(args)==0
    one,two=tmp_path/'one',tmp_path/'two'
    assert {p.name:p.read_bytes() for p in one.iterdir()}=={p.name:p.read_bytes() for p in two.iterdir()}
    import csv
    projections=list(csv.DictReader((one/'feature_projections.tsv').open(),delimiter='\t'))
    assert len(projections)==4
    assert {r['reference_sequence_source'] for r in projections}=={'same_release_feature_peptides'}
    assert all(r['status']!='unresolved_reference_sequence' for r in projections)
    pair=list(csv.DictReader((one/'pair_summary.tsv').open(),delimiter='\t'))[0]
    assert pair['builtin_lost_in_B']=='2'
    with pytest.raises(ValueError,match='empty output'):
        run_protein_impact(args)


def _two_gene_fixture(tmp_path, partner_column):
    """Gene A protein MAL (feature 2-3); partner gene B protein WYKV (feature 2-4); fusion MALWYKV."""
    ga, gb = 'ENSG00000000001', 'ENSG00000000002'
    texts = {
        'protein_relationships': 'Gene\tTrans\tProtein\n'+ga+'\tENST1\tENSP1\n'+gb+'\tENST2\tENSP2\n',
        'protein_features': 'ID\tAA_Start\tAA_Stop\tStart\tStop\tName\tInterpro\tDescription\n'
                            'ENSP1\t2\t3\t4\t9\tDA\tIPR1\tDomainA\nENSP2\t2\t4\t104\t112\tDB\tIPR2\tDomainB\n',
        'uniprot_features': 'ensembl_prot\taa_start\taa_stop\tgenomic_start\tgenomic_stop\tname\tinterpro_id\tdescription\n',
        'protein_coordinates': 'protienID\texonID\tAA_NT_Start\tAA_NT_Stop\tGenomic_Start\tGenomic_Stop\n'
                               'ENSP1\tE1\t1\t9\t1\t9\nENSP2\tE1\t1\t12\t101\t112\n',
        'protein_sequences': 'ENSP1\tMAL\nENSP2\tWYKV\n',
        'ensembl_peptides': '>ENSP1.1 gene:'+ga+' transcript:ENST1\nMAL\n>ENSP2.1 gene:'+gb+' transcript:ENST2\nWYKV\n',
    }
    resources = {}
    for key, text in texts.items():
        p = tmp_path / (key+'.txt'); p.write_text(text)
        resources[key] = dict(path=p.name, sha256=sha(p), source='synthetic test fixture', release='same-release-test')
    bundle = tmp_path / 'bundle.json'
    bundle.write_text(json.dumps(dict(schema='altanalyze-native-protein-impact-v1', release_decision='test fixture', resources=resources)))
    header = 'gene_symbol\tisoform_id\tensg\tprotein_sequence\tcds_sequence\tcdna_sequence\tprotein_length\tannotated\tnmd_status'
    fusion = 'GA\tf.source\t'+ga+'\tMALWYKV\tATGGCTCTGTGGTATAAAGTT\tATGGCTCTGTGGTATAAAGTT\t7\tno\tNot-NMD'
    single = 'GA\ta\t'+ga+'\tMAL\tATGGCTCTG\tATGGCTCTG\t3\tyes\tNot-NMD'
    if partner_column is not None:
        header += '\tpartner_ensg\tpartner_gene_symbol'
        fusion += '\t'+partner_column+'\tGB'
        single += '\t\t'
    source = tmp_path/'source.tsv'
    source.write_text('\n'.join([header, single, fusion])+'\n')
    roster = tmp_path/'roster.json'
    roster.write_text(json.dumps(dict(source_table_sha256=sha(source), required_roster=[dict(gene_symbol='GA', isoform_id=i) for i in ['a', 'f.source']])))
    return SimpleNamespace(native_models=str(source), required_roster=str(roster), reference_bundle=str(bundle),
                           artifacts_dir=str(tmp_path/'out'), pseudobulk=None, metadata=None, state=None,
                           condition=None, donor=None, gff=None, coding_regions=None, min_isoform_reads=30, prior_feature_audit=None)


def test_partner_gene_features_are_projected_and_compared_for_a_multi_gene_model(tmp_path):
    import csv
    from altanalyze3.components.iso2function.protein_impact import run_protein_impact
    args = _two_gene_fixture(tmp_path, 'ENSG00000000002')
    assert run_protein_impact(args) == 0
    out = tmp_path/'out'
    proj = list(csv.DictReader((out/'feature_projections.tsv').open(), delimiter='\t'))
    partner = {r['isoform_id']: r for r in proj if r['reference_protein'] == 'ENSP2'}
    assert set(partner) == {'a', 'f.source'}
    assert partner['f.source']['sequence_present'] == 'True' and partner['a']['sequence_present'] == 'False'
    pair = list(csv.DictReader((out/'pair_summary.tsv').open(), delimiter='\t'))[0]
    assert (pair['isoform_A'], pair['isoform_B']) == ('f.source', 'a') and pair['builtin_lost_in_B'] == '1'
    audit = {r['isoform_id']: r for r in csv.DictReader((out/'native_model_audit.tsv').open(), delimiter='\t')}
    assert audit['f.source']['partner_ensg'] == 'ENSG00000000002' and audit['f.source']['partner_gene_symbol'] == 'GB'
    manifest = json.loads((out/'impact_manifest.json').read_text())
    assert manifest['multi_gene_partner_features']['partners'] == {'GA': ['ENSG00000000002']}


def test_single_gene_models_ignore_partner_gene_reference_rows(tmp_path):
    import csv
    from altanalyze3.components.iso2function.protein_impact import run_protein_impact
    args = _two_gene_fixture(tmp_path, None)
    assert run_protein_impact(args) == 0
    out = tmp_path/'out'
    proj = list(csv.DictReader((out/'feature_projections.tsv').open(), delimiter='\t'))
    assert {r['reference_protein'] for r in proj} == {'ENSP1'}
    assert 'partner_ensg' not in next(csv.DictReader((out/'native_model_audit.tsv').open(), delimiter='\t'))
    assert 'multi_gene_partner_features' not in json.loads((out/'impact_manifest.json').read_text())


def test_invalid_partner_gene_is_rejected(tmp_path):
    from altanalyze3.components.iso2function.protein_impact import run_protein_impact
    args = _two_gene_fixture(tmp_path, 'ENSG00000000001')
    with pytest.raises(ValueError, match='Invalid partner_ensg'):
        run_protein_impact(args)
