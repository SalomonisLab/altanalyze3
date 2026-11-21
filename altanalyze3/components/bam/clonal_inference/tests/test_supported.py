import unittest,tempfile,json
from pathlib import Path
import numpy as np,pandas as pd
from types import SimpleNamespace
from altanalyze3.components.bam.clonal_inference.supported import assignments,run_supported
from altanalyze3.components.bam.clonal_inference.__main__ import main,algorithm_spec
from altanalyze3.components.bam.clonal_inference.integrity import sha256

class SupportedTests(unittest.TestCase):
 def test_unknown_does_not_inherit_rare_mutation(self):
  X=np.array([[1,1,1]]+[[1,1,-1]]*99+[[0,0,0],[0,-1,0],[-1,-1,-1]],dtype=np.int8)
  with tempfile.TemporaryDirectory() as d:
   c=pd.DataFrame({'cell_id':[str(i) for i in range(len(X))],'sample':'s'});e=pd.DataFrame({'feature':['A','B','C']})
   run_supported(X,c,e,Path(d));a=pd.read_csv(Path(d)/'event_frequency_audit.tsv',sep='\t');self.assertEqual(a.reported_MUT.tolist(),[100,100,1]);self.assertTrue((a.unsupported_positive_cells==0).all())
   f=pd.read_csv(Path(d)/'clone_frequencies.tsv',sep='\t');self.assertAlmostEqual(f.percent_all.sum(),100);self.assertEqual(set(f.clone)&{'WT','WTobs','NE'},{'WT','WTobs','NE'})
 def test_WT_and_unknown_retained(self):
  X=np.array([[1,0,-1],[1,-1,0],[0,0,0],[0,-1,-1],[-1,-1,-1]])
  self.assertEqual(assignments(X).tolist(),['100','100','WT','WTobs','NE'])
  with self.assertRaises(ValueError):assignments(np.array([[2]]))
 def test_cell_order_does_not_change_frequency(self):
  X=np.array([[1,0],[1,-1],[0,1],[-1,-1]])
  self.assertEqual(sorted(assignments(X)),sorted(assignments(X[::-1])))
 def test_explicit_cli_and_contract_guard(self):
  with tempfile.TemporaryDirectory() as d:
   r=Path(d);X=np.array([[1,1,-1],[1,0,1],[0,-1,0]],dtype=np.int8);ids=np.array(['a','b','c']);features=np.array(['A','B','C'])
   np.savez(r/'calls.npz',X=X,cell_ids=ids,features=features);pd.DataFrame({'cell_id':ids,'sample':['s1','s2','s2']}).to_csv(r/'cells.tsv',sep='\t',index=False);pd.DataFrame({'feature':features}).to_csv(r/'events.tsv',sep='\t',index=False);(r/'independent').write_text('independent fixture')
   c=dict(schema_version=1,method='supported-combinations',algorithm=algorithm_spec(SimpleNamespace(method='supported-combinations',seed=0)),source_precedence='fixture accepted calls',expected_cell_ids=list(ids),expected_samples=['s1','s2','s2'],expected_features=list(features),provenance_sources=[dict(path='independent',sha256=sha256(r/'independent'))],input_sha256={k:sha256(r/f) for k,f in [('matrix','calls.npz'),('cells','cells.tsv'),('events','events.tsv')]});(r/'contract.json').write_text(json.dumps(c))
   argv=['--matrix',str(r/'calls.npz'),'--cells',str(r/'cells.tsv'),'--events',str(r/'events.tsv'),'--contract',str(r/'contract.json'),'--output',str(r/'out')]
   with self.assertRaisesRegex(ValueError,'detector'):main(argv)
   self.assertFalse((r/'out').exists())
   argv+=['--method','supported-combinations'];main(argv);self.assertTrue((r/'out/COMPLETE.json').exists())
   with self.assertRaisesRegex(ValueError,'overwrite'):main(argv)
   (r/'cells.tsv').write_text('bad')
   with self.assertRaisesRegex(ValueError,'Changed cells'):main(argv)
