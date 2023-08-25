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
        idx = pd.MultiIndex.from_tuples([("U235", "(z,fission)", 1e-7),
                                         ("U235", "(z,fission)", 5.4e-7)],
                                         names=['N', 'MT', 'E [MeV]'])
        assert (self.scs.idx == idx).all()

    def test_s1_(self):
        idx = pd.MultiIndex.from_tuples([("U235", "(z,fission)", 1e-7),
                                         ("U235", "(z,fission)", 5.4e-7)],
                                         names=['N', 'MT', 'E [MeV]'])
        assert (self.scs.s1_.index == idx).all()
        assert np.array([check_ufloat_equality(i, j) for i, j
                 in zip(self.scs.s1_,
                        unumpy.uarray([1, 2], [0.01, 0.02]))]).all()

    def test_s2_(self):
        idx = pd.MultiIndex.from_tuples([("U235", "(z,fission)", 1e-7),
                                         ("U235", "(z,fission)", 5.4e-7)],
                                         names=['N', 'MT', 'E [MeV]'])
        assert (self.scs.s2_.index == idx).all()
        assert np.array([check_ufloat_equality(i, j) for i, j
                 in zip(self.scs.s2_,
                        unumpy.uarray([4, 5], [0.04, 0.05]))]).all()

    def test_cov_(self):
        idx = pd.MultiIndex.from_tuples([("U235", "(z,fission)", 1e-7),
                                         ("U235", "(z,fission)", 5.4e-7)],
                                         names=['N', 'MT', 'E [MeV]'])
        assert ((self.scs.cov_ == pd.DataFrame([[1, .1],
                                              [.1, 1]],
                                             index=idx,
                                             columns=idx)).all()).all()

    def test_sandwich(self):
        assert check_ufloat_equality(self.scs.sandwich,
                                     ufloat(15.3, 0.1643715303816327))

    def test_representativity(self):
        assert check_ufloat_equality(self.scs.representativity,
                                     ufloat(0.9814954576223638,
                                            0.0015588325271614727))

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
                                     ufloat(2.32379000772445,
                                            0.01879716290649558))
        assert check_ufloat_equality(self.scs.stdev(2),
                                     ufloat(6.708203932499369,
                                            0.0483735464897913))

    def test_nuclide_apportion(self):
        pass

    def test_reaction_apportion(self):
        pass
