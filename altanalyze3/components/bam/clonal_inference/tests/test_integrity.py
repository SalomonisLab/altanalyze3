import hashlib
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from altanalyze3.components.bam.clonal_inference.integrity import load_inputs, numerical_gate, sha256
from altanalyze3.components.bam.clonal_inference.reporting import export_point_assignments


class IntegrityTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.root = Path(self.tmp.name)
        self.X = np.array([[1,0],[-1,-1],[0,1]], dtype=np.int8)
        self.ids = np.array(['a','b','c']); self.features = np.array(['v1','v2'])
        np.savez(self.root/'calls.npz', X=self.X, cell_ids=self.ids, features=self.features)
        pd.DataFrame({'cell_id': self.ids, 'sample':['s1','s1','s2']}).to_csv(self.root/'cells.tsv', sep='\t', index=False)
        pd.DataFrame({'feature': self.features}).to_csv(self.root/'events.tsv', sep='\t', index=False)
        (self.root/'source.txt').write_text('independent primary roster a b c / v1 v2')
        self.contract = dict(schema_version=1, method='exhaustive', source_precedence='accepted calls unchanged',
                             expected_cell_ids=list(self.ids), expected_features=list(self.features), expected_samples=['s1','s1','s2'],
                             provenance_sources=[dict(path='source.txt',sha256=sha256(self.root/'source.txt'))],
                             input_sha256={k:sha256(self.root/f) for k,f in [('matrix','calls.npz'),('cells','cells.tsv'),('events','events.tsv')]})
        self.save_contract()

    def tearDown(self): self.tmp.cleanup()
    def save_contract(self): (self.root/'contract.json').write_text(json.dumps(self.contract))
    def load(self): return load_inputs(self.root/'calls.npz',self.root/'cells.tsv',self.root/'events.tsv',self.root/'contract.json','exhaustive')

    def test_complete_and_unknown(self):
        X,c,e,_=self.load(); self.assertTrue(np.array_equal(X,self.X))
        out=self.root/'out';out.mkdir()
        a=export_point_assignments(X,c,e,np.array(['1','1','2']),pd.DataFrame([[1,0],[0,1]],index=['1','2']),out)
        self.assertEqual(a.loc[1,'display_clone'],'Unresolved_no_evidence')
        f=pd.read_csv(out/'genotype_group_frequencies.tsv',sep='\t')
        self.assertEqual(f.groupby('sample').cells.sum().to_dict(),{'s1':2,'s2':1})

    def test_shuffled_identity_rejected_even_with_new_file_hash(self):
        np.savez(self.root/'calls.npz',X=self.X,cell_ids=self.ids[::-1],features=self.features)
        self.contract['input_sha256']['matrix']=sha256(self.root/'calls.npz');self.save_contract()
        with self.assertRaisesRegex(ValueError,'identity/order'):self.load()

    def test_feature_exclusion_rejected_even_with_new_file_hash(self):
        np.savez(self.root/'calls.npz',X=self.X[:,:1],cell_ids=self.ids,features=self.features[:1])
        self.contract['input_sha256']['matrix']=sha256(self.root/'calls.npz');self.save_contract()
        with self.assertRaisesRegex(ValueError,'Event identity'):self.load()

    def test_source_change_rejected(self):
        (self.root/'source.txt').write_text('changed')
        with self.assertRaisesRegex(ValueError,'Changed provenance'):self.load()

    def test_wrong_method_rejected(self):
        self.contract['method']='bnpc-native';self.save_contract()
        with self.assertRaisesRegex(ValueError,'Contract method'):self.load()

    def test_nonfinite_rejected(self):
        X=self.X.astype(float);X[0,0]=np.nan;np.savez(self.root/'calls.npz',X=X,cell_ids=self.ids,features=self.features)
        self.contract['input_sha256']['matrix']=sha256(self.root/'calls.npz');self.save_contract()
        with self.assertRaisesRegex(ValueError,'nonfinite'):self.load()

    def test_backend(self): numerical_gate()


if __name__=='__main__':unittest.main()
