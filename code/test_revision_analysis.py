"""Numerical and leakage checks for the revised scientific analysis."""
import unittest
from itertools import combinations
import numpy as np
from scipy.spatial.distance import cdist

import revision_analysis as r


class DistanceTests(unittest.TestCase):
    def test_hand_calculated_distances_and_empty_vectors(self):
        z = np.array([[1, 1, 0, 0], [0, 0, 0, 0]])
        a = np.array([[1, 0, 1, 0], [0, 0, 0, 0]])
        ds = dict(r.distance_matrices(z, a))
        self.assertAlmostEqual(ds['jaccard'][0,0], 2/3)
        self.assertAlmostEqual(ds['dice'][0,0], 1/2)
        self.assertAlmostEqual(ds['hamming'][0,0], 1/2)
        for name,d in ds.items():
            self.assertEqual(d[1,1],0,name)
            self.assertTrue(np.isfinite(d).all(),name)

    def test_against_independent_library(self):
        rng = np.random.default_rng(42)
        z = rng.integers(0,2,(9,13)).astype(bool)
        a = rng.integers(0,2,(11,13)).astype(bool)
        ds = dict(r.distance_matrices(z,a))
        for metric in ['jaccard','dice','hamming']:
            np.testing.assert_allclose(ds[metric],cdist(z,a,metric=metric),atol=1e-12)
        w = np.log(len(a)/np.maximum(1,a.sum(axis=0)))
        expected = np.array([[np.sum(w*(u!=v))/np.sum(w*(u|v)) for v in a] for u in z])
        np.testing.assert_allclose(ds['idf_jaccard'],expected,atol=1e-12)

    def test_hamming_invariant_to_column_polarity(self):
        _,x,m,_,_ = r.load_design()
        a = x.loc[m.role.isin(['pronoun_anchor','compound_anchor'])].to_numpy()
        z = x.loc[list(r.TARGETS)].to_numpy()
        original = dict(r.distance_matrices(z,a))['hamming']
        a[:,[0,9,43]] = 1-a[:,[0,9,43]]
        z[:,[0,9,43]] = 1-z[:,[0,9,43]]
        np.testing.assert_allclose(dict(r.distance_matrices(z,a))['hamming'],original)

    def test_percentile_is_not_a_category_midpoint(self):
        # A positive contrast can still be at the middle of a biased reference.
        summary = r.rank_summary(.2,np.array([.1,.2,.3]))
        self.assertAlmostEqual(summary['percentile'],.5)
        self.assertAlmostEqual(summary['lower_tail'],.75)
        self.assertEqual(r.rank_summary(-1,np.array([0.,1.]))['lower_tail'],1/3)


class DesignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.protocol,cls.x,cls.m,cls.blocks,cls.excluded = r.load_design()

    def test_manifest_and_singleton_filter(self):
        self.assertEqual(self.x.shape,(138,154))
        self.assertEqual(self.excluded,['each_other'])
        self.assertEqual(self.m.role.value_counts()['pronoun_anchor'],41)
        compounds = set(self.m.index[self.m.role.eq('compound_anchor')])
        self.assertEqual(len(compounds),16)
        self.assertTrue({'anything','everything','nothing','something','anywhere'} <= compounds)
        self.assertFalse({'us_det','we_det','you_det'} & compounds)

    def test_no_target_or_control_self_inclusion(self):
        for item in self.m.index[self.m.role.ne('context')]:
            for pool in self.protocol['comparator_pools']:
                for holdout in self.protocol['holdouts']:
                    mask = r.comparator_mask(self.m,item,pool,holdout)
                    self.assertFalse(mask.loc[item])
                    self.assertFalse(mask.loc[list(r.TARGETS)].any())
                    if holdout=='family':
                        self.assertFalse((self.m.loc[mask,'family']==self.m.loc[item,'family']).any())
        # Explicitly test both case/reading variants and shared compound bases.
        she = r.comparator_mask(self.m,'she','full','family')
        self.assertFalse(she.loc[['she','herself','her_acc','her_dep','hers']].any())
        every = r.comparator_mask(self.m,'everyone','full','family')
        self.assertFalse(every.loc[['everyone','everybody','everything','everywhere']].any())

    def test_randomization_constraints_and_support(self):
        x = np.array([1,0,1,0,0,1])
        b = np.array(['a','a','b','b','b','b'])
        z = r.reference_profiles(x,b,'within_block',1000,np.random.default_rng(12))
        self.assertTrue((z[:,:2].sum(axis=1)==1).all())
        self.assertTrue((z[:,2:].sum(axis=1)==2).all())
        z = r.reference_profiles(x,b,'fixed_total',10000,np.random.default_rng(12))
        patterns,counts = np.unique(z,axis=0,return_counts=True)
        self.assertEqual(len(patterns),len(list(combinations(range(6),3))))
        self.assertTrue((counts>400).all() and (counts<600).all())
        for edge in [np.zeros(6),np.ones(6)]:
            np.testing.assert_array_equal(r.reference_profiles(edge,b,'within_block',2,np.random.default_rng(12)),np.tile(edge,(2,1)))


if __name__=='__main__':
    unittest.main()
