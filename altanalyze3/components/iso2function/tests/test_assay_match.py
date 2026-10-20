import argparse
import csv
import json

import pytest

from altanalyze3.components.iso2function.cli import add_arguments,run_iso2function
from altanalyze3.components.iso2function.assay_match import sequence_differences


def test_literal_substitutions_and_indels_have_explicit_residue_coordinates():
    assert sequence_differences('MACDEF','MACDEF') == []
    d=sequence_differences('MACDEF','MAYDEF')[0]
    assert (d['native_start'],d['native_end'],d['native_sequence'],d['assayed_sequence'])==(3,3,'C','Y')
    d=sequence_differences('MACDEF','MACYYDEF')[0]
    assert d['change']=='insert' and d['native_sequence']=='' and d['native_start']==''
    assert d['assayed_sequence']=='YY'


def test_standard_cli_keeps_mismatched_native_pair_and_uses_real_assay_query(tmp_path):
    clones=tmp_path/'clones.tsv'
    clones.write_text('clone_id\tgene_symbol\taa_seq\nG-1\tG\tMAL\nG-2\tG\tMM\n')
    native=tmp_path/'native.fa'
    native.write_text('>A gene_id:ENSG00000000001\nMAL\n>B.source001 gene_id:ENSG00000000001\nMA\n')
    pairs=tmp_path/'pairs.tsv'
    pairs.write_text('gene\tensg\tnative_A\tnative_B\tclone_A\tclone_B\nG\tENSG00000000001\tA\tB.source001\tG-1\tG-2\n')
    data=tmp_path/'data';data.mkdir()
    (data/'ppi_y2h.tsv').write_text('ad_clone_id\tad_gene_symbol\tdb_gene_symbol\tdetected\ttested\nG-1\tG\tX\tTrue\tTrue\nG-2\tG\tX\tFalse\tTrue\nG-1\tG\tY\tTrue\tTrue\nG-2\tG\tY\tTrue\tTrue\n')
    (data/'pdi_ey1h.tsv').write_text('clone_id\tgene_symbol\tbait_id\tdetected\ttested\n')
    (data/'activation_m1h.tsv').write_text('clone_id\tM1H_mean\tm1h_call\nG-1\t2\tactivator\nG-2\t-2\trepressor\n')
    parser=add_arguments(argparse.ArgumentParser())
    args=parser.parse_args(['assay-match','--pairs',str(pairs),'--protein-fasta',str(native),'--clonelist',str(clones),'--out-dir',str(data),'--artifacts-dir',str(tmp_path/'results')])
    assert run_iso2function(args)==0
    out=tmp_path/'results'
    rows=list(csv.DictReader((out/'all_native_pair_evidence.tsv').open(),delimiter='\t'))
    assert len(rows)==1
    r=rows[0]
    assert r['exact_native_A']=='True' and r['exact_native_B']=='False'
    assert r['both_native_proteins_exact']=='False'
    assert r['assay_applies_to_native_pair']=='related_protein_sequences_only'
    assert r['ppi_A_only']=='1' and r['ppi_A_only_targets']=='X'
    assert r['ppi_shared_positive']=='1' and r['ppi_shared_positive_targets']=='Y'
    assert (r['M1H_A_call'],r['M1H_B_call'])==('activator','repressor')
    assert json.load((out/'assay_match_manifest.json').open())['status']=='complete'
    assert (out/'assay_calls/query_manifest.json').exists()
