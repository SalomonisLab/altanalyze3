import csv,gzip,pathlib,sqlite3,tempfile,unittest
import sys
from pathlib import Path
from altanalyze3.components.bam import noncolinear_evidence as M
class Tests(unittest.TestCase):
 def test_crossings_identity_and_retention(self):
  with tempfile.TemporaryDirectory() as d:
   p=pathlib.Path(d);db=sqlite3.connect(':memory:');M.initialize(db)
   base=dict(chrom='chr1',a_strand='+',b_strand='+',donor='200',acceptor='100',joint_type='back',library='L1',status='kept',query_name='read1',molecule_id='zm1',barcode='CB',umi='UMI',outside_circle_nt='0',pair_overlap_nt='20',crossings='2',polyA_at_joint='False',junction_homology_nt='2',query_gap='0',resolved='True',donor_annotation='G:E2',acceptor_annotation='G:E1',donor_region='G:E2',acceptor_region='G:E1',joint_class='backsplice_registered_canonical')
   rows=[base,base,dict(base,query_name='read2'),dict(base,library='L2'),dict(base,status='mapq_below_min',joint_class='backjoint_unregistered')]
   with gzip.open(p/'input.gz','wt') as f:
    w=csv.DictWriter(f,list(base),delimiter='\t');w.writeheader();w.writerows(rows)
   self.assertEqual(M.ingest(db,p/'input.gz'),5);M.finalize(db,p)
   with gzip.open(p/'event_library_status_counts.tsv.gz','rt') as f:r=list(csv.DictReader(f,delimiter='\t'))
   self.assertEqual(len(r),3)
   x=next(x for x in r if x['library']=='L1' and x['status']=='kept')
   self.assertEqual((x['raw_pair_rows'],x['unique_read_count'],x['unique_declared_molecule_id_count'],x['unique_cell_umi_count']),('3','2','1','1'))
   with gzip.open(p/'read_support.tsv.gz','rt') as f:s=list(csv.DictReader(f,delimiter='\t'))
   self.assertEqual(len(s),4);self.assertEqual(sum(int(x['raw_pair_rows']) for x in s),5)
   self.assertEqual(sum(int(x['pair_rows']) for x in db.execute('SELECT pair_rows FROM classes')) if False else db.execute('SELECT SUM(pair_rows) FROM classes').fetchone()[0],5)
 def test_missing_umi_and_identity_gate(self):
  self.assertEqual(M.event_id(dict(chrom='chr21',a_strand='-',b_strand='-',donor='100',acceptor='200',joint_type='back')),'BSJ|chr21|-|-|100|200')
if __name__=='__main__':unittest.main()
