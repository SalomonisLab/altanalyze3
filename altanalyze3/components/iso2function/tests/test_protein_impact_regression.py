"""Optional real-data audits against independent saved evidence, never runtime predictors."""
import json
import os
from pathlib import Path

import pytest
from altanalyze3.components.iso2function.query import _table
from altanalyze3.components.iso2function.protein_impact import sha


def results():
    value=os.environ.get('ISO2FUNCTION_VERIFIED_EVIDENCE')
    if not value:pytest.skip('Set ISO2FUNCTION_VERIFIED_EVIDENCE to the completed HSPC verification directory')
    return Path(value)


def test_complete_native_roster_and_counts_and_repeat_are_unchanged():
    base=results();out=base/'results';repeat=base/'results_repeat'
    source=Path(__file__).parents[3]/'output/hspc_tf_synthesis_20261020/ORIGINAL_PANEL_ANNOTATION_AND_COUNTS.tsv'
    expected=_table(source);native=_table(out/'native_model_audit.tsv')
    key=lambda r:(r['gene_symbol'],r['isoform_id'])
    assert [key(r) for r in native]==[key(r) for r in expected]
    assert len(native)==64 and len({r['gene_symbol'] for r in native})==33
    e={key(r):r for r in expected}
    for r in native:
        old=e[key(r)]
        assert r['protein_sha256']==old['protein_sha256']
        assert float(r['assigned_reads'])==float(old['healthy_HSC_MPP_assigned_reads'])
    assert {p.name:sha(p) for p in out.iterdir()}=={p.name:sha(p) for p in repeat.iterdir()}
    manifest=json.loads((out/'impact_manifest.json').read_text())
    assert manifest['prior_regression_changes']==0
    for name,digest in manifest['outputs'].items():assert sha(out/name)==digest
    hlf=[r for r in _table(out/'native_genomic_audit.tsv') if r['gene_symbol']=='HLF' and r['source_table_aa']=='144']
    assert len(hlf)==1 and hlf[0]['protein_fasta_aa']==hlf[0]['coding_regions_aa']=='169'
    assert hlf[0]['length_agrees']=='False'


def test_key_feature_claims_have_explicit_database_rows_and_qualifications():
    rows=_table(results()/'results/pair_feature_impacts.tsv')
    def find(gene,description,start,end,provider='uniprot_features'):
        found=[r for r in rows if r['gene_symbol']==gene and r['provider']==provider and r['description']==description
               and r['reference_start']==str(start) and r['reference_end']==str(end)]
        assert found
        return found[0]
    assert find('HOXB5','DNA_BIND-Homeobo',194,253)['AA_comparison']=='shared_present'
    assert find('HOXB5','MOTIF-Antp-type hexapeptid',176,181)['B_status']=='not_aligned'
    nr=find('NR4A2','DOMAIN-NR LB',360,595)
    assert nr['B_status']=='altered_or_partial' and float(nr['B_exact_fraction'])==95/236
    assert find('NR4A2','DNA_BIND-Nuclear recepto',260,335)['AA_comparison']=='shared_present'
    assert find('TOX','DNA_BIND-HMG bo',261,329)['B_status']=='not_aligned'
    assert find('LYL1','DOMAIN-bHL',150,202)['AA_comparison']=='shared_present'
    assert find('CITED2','CBP/p300-binding activation segment (literature construct)',216,269,'literature_construct_snapshot')['B_status']=='not_aligned'
    bcl=[r for r in rows if r['gene_symbol']=='BCL11A']
    assert all(r['AA_comparison'] not in ['lost_in_B','gained_in_B'] for r in bcl)


def test_substitution_query_keeps_every_baseline_pair_and_difference():
    b=results();out=b/'interaction_query';old=b.parent/'HSPC_PPI_PDI_CLI_20261020/results'
    pair=lambda r:(r['gene'],r['clone_A'],r['clone_B'])
    before=_table(old/'all_pair_decisions.tsv');after=_table(out/'all_pair_decisions.tsv')
    assert {pair(r) for r in before}=={pair(r) for r in after} and len(after)==803
    before_kept={pair(r) for r in _table(old/'retained_pairs.tsv')};after_kept={pair(r) for r in _table(out/'retained_pairs.tsv')}
    assert before_kept <= after_kept and len(after_kept)==35
    edge=lambda r:(*pair(r),r['modality'],r['target'],r['call_A'],r['call_B'])
    previous={edge(r) for r in _table(old/'divergent_interactions.tsv')}
    now=_table(out/'divergent_interactions.tsv')
    assert previous <= {edge(r) for r in now} and len(now)==556
    assert {r['gene'] for r in now if r['modality']=='ppi'} >= {'PBX1','MLX','NR4A1','NFYA','ZSCAN9','ZBTB44'}
    assert {r['target'] for r in now if r['gene']=='PBX1'}=={'PIN1','TMF1'}
    assert all(float(r['max_state_mean_cpm'])>1 for r in now if r['modality']=='ppi')
    assert all(r['expression_status']=='DNA_bait_not_gene_expression' and not r['max_state_mean_cpm'] for r in now if r['modality']=='pdi')
    assert all(float(r['max_single_native_reads_'+s])>=30 for r in _table(out/'retained_pairs.tsv') for s in ('A','B'))
    assert not any(r['gene']=='TCF4' and r['modality']=='pdi' for r in now)


def test_all_earlier_assay_comparisons_and_mismatch_cases_are_retained():
    out=results()/'assay_support'
    records=_table(out/'all_native_pair_evidence.tsv')
    assert len(records)==25
    assert {r['gene'] for r in records if r['both_native_proteins_exact']=='False'}=={'PBX1','NFYA','ZSCAN9','ZBTB44','NR4A1','MLX','ZNF655'}
    native=_table(out/'native_assayed_protein_comparison.tsv')
    assert len(native)==50 and sum(r['exact_sequence_match']=='False' for r in native)==9
    old=Path(__file__).parents[3]/'output/hspc_tf_synthesis_20261020/STRICT_PAIR_FUNCTION_SUMMARY.tsv'
    expected={(r['gene'],r['clone_A'],r['clone_B']):r for r in _table(old)}
    for r in records:
        before=expected[r['gene'],r['clone_A'],r['clone_B']]
        for modality in ['PPI','PDI']:
            for call in ['A_only','B_only']:
                assert r[modality.lower()+'_'+call]==before[modality+'_'+call]
    znf=next(r for r in records if r['gene']=='ZNF655')
    assert znf['M1H_A_call']==znf['M1H_B_call']=='neither'
