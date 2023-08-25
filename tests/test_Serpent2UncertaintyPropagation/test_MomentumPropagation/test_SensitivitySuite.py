# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 14:05:51 2023

@author: fgrimald
"""
import pandas as pd
import numpy as np

# equality does not work with uncertainty
from uncertainties import unumpy

from functions import check_ufloat_equality

import SEAT.Serpent2UncertaintyPropagation.MomentumPropagation.SensitivitySuite as SS

class Test_Functions:
    def test_get_zais(self):
        pass

    def test_get_perts_mt(self):
        pass

    def test_get_perts_rea(self):
        pass

    def test_get_ene(self):
        pass

    def test_sts2df(self):
        pass

class Test_Sensitivity:
    idx = pd.MultiIndex.from_tuples([("a", "U235", "(z,fission)", 1e-7),
                                     ("a", "U235", "(z,fission)", 5.4e-7),
                                     ("a", "Pu239", "(z,fission)", 1e-7),
                                     ("a", "Pu239", "(z,fission)", 5.4e-7),
                                     ("b", "U235", "(z,fission)", 1e-7),
                                     ("b", "U235", "(z,fission)", 5.4e-7),
                                     ("b", "Pu239", "(z,fission)", 1e-7),
                                     ("b", "Pu239", "(z,fission)", 5.4e-7),],
                                     names=['Observable', 'N', 'MT', 'E [MeV]'])
    s_ = np.array([1, 2, 3, 4, 5, 6, 7, 8])
    u_ = np.array([.1, .1, .1, .1, .1, .1, .1, .1]) # absolute statistical uncertainty

    df = pd.DataFrame({"S": unumpy.uarray(s_, u_), "stat. err.": u_ / s_,
                       "nom_val": s_, "abs_sdev": u_}, index=idx)

    reference = SS.Sensitivity.from_df(df)

    def test_from_csv(self):
        pass

    def test_from_df(self):
        assert self.reference.file is None
        assert ((self.reference._df == self.df).all()).all()

    def test_idx(self):  # to be complemented
        assert (self.reference.idx == self.df.index).all()

    def test_data(self):
        assert ((self.reference._df == self.df).all()).all()

    def test_reader(self):
        pass

    def test_observables(self):
        pass

    def test_upper(self):
        up = self.df.S + self.df.S.apply(lambda x: x.nominal_value) * self.df["stat. err."]
        assert (self.reference.upper.data.index == up.index).all()
        assert np.array([check_ufloat_equality(i, j) for i, j
                         in zip(self.reference.upper.data.S, up)]).all()

    def test_lower(self):
        low = self.df.S - self.df.S.apply(lambda x: x.nominal_value) * self.df["stat. err."]
        assert (self.reference.lower.data.index == low.index).all()
        assert np.array([check_ufloat_equality(i, j) for i, j
                         in zip(self.reference.lower.data.S, low)]).all()

    def test_observe(self):
        assert ((self.reference.observe("a").data == self.df.iloc[:4]).all()).all()

    def test_query(self):
        assert ((self.reference.query("nom_val == 1").data == self.df.iloc[0]).all()).all()
