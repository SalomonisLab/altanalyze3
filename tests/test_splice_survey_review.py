import csv
import gzip
import json
from pathlib import Path
import pytest
from altanalyze3.components.annotation.splice_survey_review import scan_usage,extract_registered


def test_usage_keeps_source_fdr_and_requires_both_opposite_significant(tmp_path):
    fields=['Gene_Symbol','Feature','N_A','N_Others','Mean_A','Mean_Others','DiffMeans_Others_m_A','ebayesAdjPval']
    records=[['GENE','G:a',3,9,.8,.2,-.6,.01],['GENE','G:b',3,9,.1,.7,.6,.02],
             ['GENE','G:c',3,9,.7,.3,-.4,.03],['GENE','G:d',3,9,.05,.6,.55,.2],
             ['GENE','G:low',3,9,.1,.7,.6,.001]]
    with open(tmp_path/'healthy-A-Others-A_stats-annotated.txt','w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(fields);w.writerows(records)
    stats,symbols,sources,pairs=scan_usage(tmp_path,['A'],{'G:a','G:b','G:c','G:d'},.05)
    assert set(stats)=={'G:a','G:b','G:c','G:d'}
    assert set(pairs)=={('G:a','G:b'),('G:b','G:c')}
    assert stats['G:b']['A']['q']==.02
    assert stats['G:a']['A']['dusage']==pytest.approx(.6)
    assert sources[0]['rows']==5 and sources[0]['high_read_significant']==3
    with pytest.raises(ValueError,match='roster'):scan_usage(tmp_path,['B'],set(),.05)


def test_literal_peptides_keep_dotted_id_and_flag_metadata_discrepancy(tmp_path):
    base=tmp_path/'run/isv/gff-output';base.mkdir(parents=True);out=tmp_path/'out';out.mkdir()
    (base/'protein_summary.txt').write_text('Gene ID\tTranscript ID\tProtein Length\tNMD Status\tIntron Retention\tLongest Isoform Length\nG;(NMD)\ta.source001\t4\tNMD\tTrue\t10\nG\tb\t0\tNMD\tFalse\t10\n')
    (base/'coding_regions.txt').write_text('Transcript ID\tGene ID\tProtein Length\na.source001\tG\t3\n')
    with gzip.open(base/'protein_sequences.fasta.gz','wt') as f:f.write('>a.source001 ;Protein;gene_id:G;(NMD)\nMAS\n>other ;Protein;gene_id:Q\nXXXX\n')
    extract_registered(tmp_path/'run',{'G:a.source001','G:b'},out)
    audit=list(csv.DictReader(open(out/'native_protein_source_audit.tsv'),delimiter='\t'))
    assert audit[0]['literal_peptide_aa']=='3' and audit[0]['length_agrees']=='False'
    assert audit[1]['peptide_present']=='False' and audit[1]['literal_peptide_aa']==''
    assert (out/'native_predicted_proteins.fasta').read_text()=='>G:a.source001\nMAS\n'
