# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 14:06:01 2023

@author: fgrimald
"""

import pandas as pd
import numpy as np

# equality does not work with uncertainty
from uncertainties import unumpy, ufloat

import SEAT.Serpent2UncertaintyPropagation.MomentumPropagation.RepresentativitySuite as RS
import SEAT.Serpent2UncertaintyPropagation.MomentumPropagation.SensitivitySuite as SS
import SEAT.Serpent2UncertaintyPropagation.MomentumPropagation.CovarianceSuite as CS

from functions import check_ufloat_equality

class Test_Functions:
    idx = pd.MultiIndex.from_tuples([("U235", "(z,fission)", 1e-7),
                                     ("U235", "(z,fission)", 5.4e-7),
                                     ("U235", "(z,fission)", 4e-6)],
                                     names=['MAT', 'MT', 'E [MeV]'])
    s1 = pd.DataFrame([1, 2, 3], index=idx)
    s2 = pd.DataFrame([4, 5, 6], index=idx)
    cov = pd.DataFrame([[1, .1, .2],
                        [.1, 1, .3],
                        [.2, .3, 1]], index=idx, columns=idx)

    def test_sandwich_unsafe(self):
        assert RS.sandwich_unsafe(self.s1, self.cov, self.s2).values == 45

    def test_sandwich(self):
        idx_ = pd.MultiIndex.from_tuples([("U235", "(z,fission)", 1e-7),
                                         ("U235", "(z,fission)", 5.4e-7),
                                         ("U238", "(z,fission)", 4e-6)],
                                         names=['MAT', 'MT', 'E [eV]'])
        s1 = self.s1.copy()
        s1.index = idx_
        assert RS.sandwich(s1, self.cov, self.s2).values == 15.3 ## first 2
        assert RS.sandwich(self.s1, self.cov, self.s2).values == 45 ## all 3

    def test_stdev(self):
        assert np.round(RS.stdev(self.s1, self.cov).values[0,0], 10) == 4.3817804600

    def test_representativity(self):
        assert np.round(RS.representativity(self.s1, self.cov, self.s2).values[0,0], 10) == 0.9854781211


class Test_SCS:
    idx1 = pd.MultiIndex.from_tuples([("a", "U235", "(z,fission)", 1e-7),
                                     ("a", "U235", "(z,fission)", 5.4e-7),
                                     ("a", "U238", "(z,fission)", 4e-6)],
                                     names=['Observable', 'N', 'MT', 'E [MeV]'])
    idx2 = pd.MultiIndex.from_tuples([("a", "U235", "(z,fission)", 1e-7),
                                     ("a", "U235", "(z,fission)", 5.4e-7),
                                     ("a", "U235", "(z,fission)", 4e-6)],
                                     names=['Observable', 'N', 'MT', 'E [MeV]'])

    s1 = SS.Sensitivity.from_df(pd.DataFrame([[1, 0.01], [2, 0.01], [3, 0.01]],
                                             index=idx1, columns=["S", 'stat. err.']))
    s2 = SS.Sensitivity.from_df(pd.DataFrame([[4, 0.01], [5, 0.01], [6, 0.01]],
                                             index=idx2, columns=["S", 'stat. err.']))

    cov = CS.Covariance(pd.DataFrame([[1, .1, .2],
                                      [.1, 1, .3],
                                      [.2, .3, 1]],
                                     index=idx1.droplevel(0),
                                     columns=idx1.droplevel(0)))

    scs = RS.SCS(s1, cov, s2)

    shared_idx = pd.MultiIndex.from_tuples([("(z,fission)", 1e-7),
                                            ("(z,fission)", 5.4e-7)],
                                     names=['MT', 'E [MeV]'])

    # All the nuclides in the two systems for which a covariance entry exists are reported
    union_intersection_idx = pd.MultiIndex.from_tuples([("U235", "(z,fission)", 1e-7),
                                                        ("U235", "(z,fission)", 5.4e-7),
                                                        ('U238', '(z,fission)',   4e-06)],
                                                        names=['N', 'MT', 'E [MeV]'])

    def test_get_sharing(self):
        a, b, c = self.scs._get_sharing("U235")
        assert (a.index == self.shared_idx).all()
        assert np.array([check_ufloat_equality(i, j) for i, j
                 in zip(a.S, unumpy.uarray([1, 2], [0.01, 0.02]))]).all()
        assert ((b == pd.DataFrame([[1, .1],
                                   [.1, 1]],
                                  index=self.shared_idx,
                                  columns=self.shared_idx)).all()).all()
        assert (c.index == self.shared_idx).all()
        assert np.array([check_ufloat_equality(i, j) for i, j
                 in zip(c.S, unumpy.uarray([4, 5], [0.04, 0.05]))]).all()

    def test_idx(self):
        assert (self.scs.idx == self.union_intersection_idx).all()

    def test_s1_(self):
        assert (self.scs.s1_.index == self.union_intersection_idx).all()
        assert np.array([check_ufloat_equality(i, j) for i, j
                 in zip(self.scs.s1_,
                        unumpy.uarray([1, 2, 3], [0.01, 0.02, 0.03]))]).all()

    def test_s2_(self):
        assert (self.scs.s2_.index == self.union_intersection_idx).all()
        assert np.array([check_ufloat_equality(i, j) for i, j
                 in zip(self.scs.s2_,
                        unumpy.uarray([4, 5, 0], [0.04, 0.05, 0]))]).all()

    def test_cov_(self):
        assert ((self.scs.cov_ == pd.DataFrame([[1. , 0.1, 0.2],
                                                [0.1, 1. , 0.3],
                                                [0.2, 0.3, 1. ]],
                                                index=self.union_intersection_idx,
                                                columns=self.union_intersection_idx
                                                )).all()).all()

    def test_sandwich(self):
        assert check_ufloat_equality(self.scs.sandwich,
                                     ufloat(22.2, 0.214788267836025))

    def test_representativity(self):
        assert check_ufloat_equality(self.scs.representativity,
                                     ufloat(0.7552593373581465,
                                            0.002761092095039872))

    def test_reaction_wise(self):
        pass

    def test_reaction_wise_rep(self):
        pass

    def test_nuclide_wise(self):
        pass

    def test_nuclide_wise_rep(self):
        pass

    def test_stdev(self):
        assert check_ufloat_equality(self.scs.stdev(1),
                                     ufloat(4.381780460041329,
                                            0.02968585521759479))
        assert check_ufloat_equality(self.scs.stdev(2),
                                     ufloat(6.708203932499369,
                                            0.0483735464897913))

    def test_stdev_no_unc(self):
        assert self.scs.stdev_no_unc(1) == 4.3817804600413295
        assert self.scs.stdev_no_unc(2) == 6.708203932499369

    def test_nuclide_apportion(self):
        pass

    def test_reaction_apportion(self):
        pass
