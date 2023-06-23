# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 14:06:01 2023

@author: fgrimald
"""

import pandas as pd
import numpy as np

import SEAT.Serpent2UncertaintyPropagation.MomentumPropagation.RepresentativitySuite as RS
import SEAT.Serpent2UncertaintyPropagation.MomentumPropagation.SensitivitySuite as SS
import SEAT.Serpent2UncertaintyPropagation.MomentumPropagation.CovarianceSuite as CS

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
        idx1 = pd.MultiIndex.from_tuples([("U235", "(z,fission)", 1e-7),
                                         ("U235", "(z,fission)", 5.4e-7),
                                         ("U238", "(z,fission)", 4e-6)],
                                         names=['MAT', 'MT', 'E [eV]'])
        idx2 = pd.MultiIndex.from_tuples([("U235", "(z,fission)", 1e-7),
                                         ("U235", "(z,fission)", 5.4e-7),
                                         ("U235", "(z,fission)", 4e-6)],
                                         names=['MAT', 'MT', 'E [MeV]'])
        s1 = self.s1.copy()
        s1.index = idx1
        s2 = self.s2.copy()
        s2.index = idx2
        assert RS.sandwich(s1, self.cov, s2).values == 15.3
        assert RS.sandwich(self.s1, self.cov, self.s2).values == 45

    def test_stdev(self):
        assert np.round(RS.stdev(self.s1, self.cov).values[0,0], 15) == 4.381780460041330

    def test_representativity(self):
        assert np.round(RS.representativity(self.s1, self.cov, self.s2).values[0,0], 15) == 0.985478121051818

    def test_u_sandwich(self):
        pass

    def test_u_representativity(self):
        pass


class Test_SCS:
    idx1 = pd.MultiIndex.from_tuples([("a", "U235", "(z,fission)", 1e-7),
                                     ("a", "U235", "(z,fission)", 5.4e-7),
                                     ("a", "U238", "(z,fission)", 4e-6)],
                                     names=['Observable', 'N', 'MT', 'E [MeV]'])
    idx2 = pd.MultiIndex.from_tuples([("a", "U235", "(z,fission)", 1e-7),
                                     ("a", "U235", "(z,fission)", 5.4e-7),
                                     ("a", "U235", "(z,fission)", 4e-6)],
                                     names=['Observable', 'N', 'MT', 'E [MeV]'])

    s1 = SS.Sensitivity.from_df(pd.DataFrame([1, 2, 3], index=idx1, columns=["S"]))
    us1 = pd.DataFrame([0, 0, 0], index=idx1)
    s2 = SS.Sensitivity.from_df(pd.DataFrame([4, 5, 6], index=idx2, columns=["S"]))
    us2 = pd.DataFrame([0, 0, 0], index=idx2)

    cov = CS.Covariance(pd.DataFrame([[1, .1, .2],
                                      [.1, 1, .3],
                                      [.2, .3, 1]],
                                     index=idx1.droplevel(0),
                                     columns=idx1.droplevel(0)))

    scs = RS.SCS(s1, us1, cov, s2, us2)

    shared_idx = pd.MultiIndex.from_tuples([("(z,fission)", 1e-7),
                                            ("(z,fission)", 5.4e-7)],
                                     names=['MT', 'E [MeV]'])

    def test_get_sharing(self):
        a, b, c = self.scs._get_sharing("U235")
        assert ((a == pd.DataFrame([1, 2],
                                  index=self.shared_idx, columns=["S"])).all()).all()
        assert ((b == pd.DataFrame([[1, .1],
                                   [.1, 1]],
                                  index=self.shared_idx,
                                  columns=self.shared_idx)).all()).all()
        assert ((c == pd.DataFrame([4, 5],
                                  index=self.shared_idx, columns=["S"])).all()).all()

    def test_idx(self):
        idx = pd.MultiIndex.from_tuples([("U235", "(z,fission)", 1e-7),
                                         ("U235", "(z,fission)", 5.4e-7)],
                                         names=['N', 'MT', 'E [MeV]'])
        assert (self.scs.idx == idx).all()

    def test_s1_(self):
        idx = pd.MultiIndex.from_tuples([("U235", "(z,fission)", 1e-7),
                                         ("U235", "(z,fission)", 5.4e-7)],
                                         names=['N', 'MT', 'E [MeV]'])
        assert (self.scs.s1_ == pd.Series([1, 2],
                                          index=idx, name="S")).all()

    def test_s2_(self):
        idx = pd.MultiIndex.from_tuples([("U235", "(z,fission)", 1e-7),
                                         ("U235", "(z,fission)", 5.4e-7)],
                                         names=['N', 'MT', 'E [MeV]'])
        assert (self.scs.s2_ == pd.Series([4, 5],
                                          index=idx, name="S")).all()

    def test_us1_(self):
        pass

    def test_us2_(self):
        pass

    def test_cov_(self):
        idx = pd.MultiIndex.from_tuples([("U235", "(z,fission)", 1e-7),
                                         ("U235", "(z,fission)", 5.4e-7)],
                                         names=['N', 'MT', 'E [MeV]'])
        assert ((self.scs.cov_ == pd.DataFrame([[1, .1],
                                              [.1, 1]],
                                             index=idx,
                                             columns=idx)).all()).all()

    def test_sandwich(self):
        assert self.scs.sandwich == 15.3

    def test_u_sandwich(self):
        pass

    def test_representativity(self):
        assert self.scs.representativity == 0.9814954576223638

    def test_u_representativity(self):
        pass

    def test_reaction_wise(self):
        pass

    def test_reaction_wise_rep(self):
        pass

    def test_nuclide_wise(self):
        pass

    def test_nuclide_wise_rep(self):
        pass

    def test_stdev(self):
        assert self.scs.stdev(1) == 2.32379000772445
        assert self.scs.stdev(2) == 6.708203932499369

    def test_nuclide_apportion(self):
        pass

    def test_reaction_apportion(self):
        pass
