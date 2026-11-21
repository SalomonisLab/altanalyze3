import unittest
import numpy as np
from altanalyze3.components.bam.clonal_inference.phylogeny_validation import validate_gain_tree, lineage_frequencies

class PhylogenyValidationTests(unittest.TestCase):
    def setUp(self):
        self.G=np.array([[0,0,0,0],[1,1,0,0],[1,1,1,0],[1,1,0,1]])
        self.parents=np.array([-1,0,1,1])
        self.ids=['normal','founder','BCOR','PHF6']
        self.events=['SRSF2','RUNX1','BCOR','PHF6']
    def check(self,G=None,parents=None):
        return validate_gain_tree(self.G if G is None else G,self.parents if parents is None else parents,self.ids,self.events,founder=['SRSF2','RUNX1'],malignant_clones=self.ids[1:])
    def test_multimutation_founder_and_sibling_branches(self):
        self.assertEqual(self.check()['acquired_events'][1].tolist(),[1,1,0,0])
    def test_rejects_mutation_combination_that_loses_founder(self):
        G=self.G.copy();G[2,0]=0
        with self.assertRaisesRegex(ValueError,'loses an ancestral'):self.check(G=G)
    def test_rejects_separate_malignant_root_without_founder(self):
        G=self.G.copy();G[3]=[0,0,0,1];p=self.parents.copy();p[3]=0
        with self.assertRaisesRegex(ValueError,'lacks required founder'):self.check(G,p)
    def test_rejects_recurrent_gain(self):
        G=self.G.copy();G[3,2]=1
        with self.assertRaisesRegex(ValueError,'more than once'):self.check(G)
    def test_rejects_cycle_unknown_latent_state_and_duplicate_clone(self):
        with self.assertRaisesRegex(ValueError,'cycle'):self.check(parents=np.array([-1,2,1,1]))
        G=self.G.copy();G[2,0]=-1
        with self.assertRaisesRegex(ValueError,'binary'):self.check(G)
        G=self.G.copy();G[2]=G[1]
        with self.assertRaisesRegex(ValueError,'zero-gain'):self.check(G)
    def test_exclusive_founder_can_be_smaller_than_descendant(self):
        out=lineage_frequencies([[20,5,50,15]],self.parents,[100])
        np.testing.assert_allclose(out['lineage_percent'],[[90,70,50,15]])
        np.testing.assert_allclose(out['exclusive_percent'],[[20,5,50,15]])
        np.testing.assert_allclose(out['unresolved_percent'],[10])
    def test_WT_unresolved_and_invalid_denominators(self):
        with self.assertRaisesRegex(ValueError,'exceed'):lineage_frequencies([[20,5,50,30]],self.parents,[100])
        with self.assertRaisesRegex(ValueError,'Nonfinite'):lineage_frequencies([[20,5,float('nan'),30]],self.parents,[100])

if __name__=='__main__':unittest.main()
