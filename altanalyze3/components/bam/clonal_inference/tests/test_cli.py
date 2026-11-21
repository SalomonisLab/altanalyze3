import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
from altanalyze3.components.bam.clonal_inference.__main__ import main,algorithm_spec
from altanalyze3.components.bam.clonal_inference.integrity import sha256


class CliTests(unittest.TestCase):
    def test_complete_roundtrip_and_algorithm_gate(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);ids=np.array([f'c{i}' for i in range(9)]);features=np.array(['v1','v2'])
            X=np.array([[1,0],[1,-1],[0,1],[-1,-1],[0,0],[1,1],[0,-1],[1,0],[-1,1]])
            cells=pd.DataFrame({'cell_id':ids,'sample':['s1']*5+['s2']*4})
            np.savez(root/'calls.npz',X=X,cell_ids=ids,features=features)
            cells.to_csv(root/'cells.tsv',sep='\t',index=False);pd.DataFrame({'feature':features}).to_csv(root/'events.tsv',sep='\t',index=False)
            (root/'primary.txt').write_text('independent roster fixture')
            spec=algorithm_spec(SimpleNamespace(method='exhaustive',seed=0,false_positive_rate=.1,false_negative_rate=.1))
            contract=dict(schema_version=1,method='exhaustive',algorithm=spec,source_precedence='fixture',expected_cell_ids=list(ids),expected_features=list(features),expected_samples=list(cells['sample']),provenance_sources=[dict(path='primary.txt',sha256=sha256(root/'primary.txt'))],input_sha256={k:sha256(root/f) for k,f in [('matrix','calls.npz'),('cells','cells.tsv'),('events','events.tsv')]})
            (root/'contract.json').write_text(json.dumps(contract));(root/'families.json').write_text(json.dumps({'double':{'v1':1,'v2':1}}))
            argv=['--matrix',str(root/'calls.npz'),'--cells',str(root/'cells.tsv'),'--events',str(root/'events.tsv'),'--contract',str(root/'contract.json'),'--method','exhaustive','--false-positive-rate','.1','--false-negative-rate','.1','--families',str(root/'families.json'),'--output',str(root/'run')]
            main(argv)
            a=pd.read_csv(root/'run/cell_assignments.tsv.gz',sep='\t');self.assertEqual(list(a.cell_id),list(ids));self.assertEqual(a.display_config.iloc[3],'Unresolved_no_evidence')
            z=np.load(root/'run/sample_000.npz',allow_pickle=False);self.assertEqual(list(z['cell_ids']),list(ids[:5]))
            self.assertTrue((root/'run/COMPLETE.json').is_file())
            with self.assertRaisesRegex(ValueError,'overwrite'):main(argv)
            rejected=argv.copy(); rejected[-1]=str(root/'irrelevant_options')
            with self.assertRaisesRegex(ValueError,'BnpC-only options'):main(rejected+['--iterations','50'])
            self.assertFalse((root/'irrelevant_options').exists())
            with self.assertRaisesRegex(ValueError,'Shared-genotype-only'):main(rejected+['--fp-cap','.2'])
            contract['algorithm']['FP']=.2;(root/'contract.json').write_text(json.dumps(contract));argv[-1]=str(root/'rejected')
            with self.assertRaisesRegex(ValueError,'Algorithm differs'):main(argv)
            self.assertFalse((root/'rejected').exists())

    def test_retained_shared_engine_and_vector_export(self):
        from altanalyze3.components.bam.clonal_inference.plotting import render,specimen_order
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);X=np.repeat([[0,0],[1,0],[1,1]],40,axis=0);X=np.vstack([X,[-1,-1]])
            ids=np.array([f'c{i}' for i in range(len(X))]);features=np.array(['a','b'])
            cells=pd.DataFrame(dict(cell_id=ids,sample=['s']*len(X),patient=['p']*len(X),stage=['Diagnosis']*len(X),stage_order=[0]*len(X),is_xenograft=[False]*len(X)))
            events=pd.DataFrame(dict(feature=features,event_label=features))
            np.savez(root/'calls.npz',X=X,T=(X>=0).astype(int),cell_ids=ids,features=features)
            cells.to_csv(root/'cells.tsv',sep='\t',index=False);events.to_csv(root/'events.tsv',sep='\t',index=False)
            (root/'independent.txt').write_text('independently specified test roster')
            args=SimpleNamespace(method='shared-genotype',seed=7,fp_cap=.1,fn_cap=.3,starts=2,k_grid=[2,3,4],max_iterations=10000,validation_seed=17)
            contract=dict(schema_version=1,method='shared-genotype',algorithm=algorithm_spec(args),source_precedence='fixture',expected_cell_ids=list(ids),expected_features=list(features),expected_samples=list(cells['sample']),provenance_sources=[dict(path='independent.txt',sha256=sha256(root/'independent.txt'))],input_sha256={k:sha256(root/f) for k,f in [('matrix','calls.npz'),('cells','cells.tsv'),('events','events.tsv')]})
            (root/'contract.json').write_text(json.dumps(contract))
            argv=['--matrix',str(root/'calls.npz'),'--cells',str(root/'cells.tsv'),'--events',str(root/'events.tsv'),'--contract',str(root/'contract.json'),'--output',str(root/'run'),'--seed','7','--starts','2','--k-grid','2,3,4','--validation-seed','17']
            main(argv+['--method','shared-genotype'])  # Retain the prior engine as an explicit baseline.
            review=render(root/'run',root/'calls.npz');self.assertEqual(review['cells'],len(X));self.assertTrue(review['SVG_text_editable']);self.assertTrue(review['SVG_contains_no_raster'])
            f=pd.read_csv(root/'run/clone_frequencies.tsv',sep='\t');self.assertEqual(f.cells.sum(),len(X));self.assertEqual(f.loc[f.clone=='NE','cells'].iloc[0],1)
            X[0,0]=1;np.savez(root/'changed.npz',X=X,cell_ids=ids,features=features)
            with self.assertRaisesRegex(ValueError,'matrix hash'):render(root/'run',root/'changed.npz')
            shuffled=pd.DataFrame(dict(sample=['post','graft','pre'],stage_order=[1,0,0],is_xenograft=[False,True,False]))
            self.assertEqual(specimen_order(shuffled).index.tolist(),['pre','graft','post'])


if __name__=='__main__':unittest.main()
