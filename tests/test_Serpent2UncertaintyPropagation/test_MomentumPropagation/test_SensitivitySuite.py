# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 14:05:51 2023

@author: fgrimald
"""
import pandas as pd

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
    df = pd.DataFrame({"S": [1, 2, 3, 4, 5, 6, 7, 8],
                       "stat. err.": [.1, .1, .1, .1, .1, .1, .1, .1]},
                      index=idx)

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
        assert ((self.reference.upper.data.S == self.df.S + self.df.S * self.df["stat. err."]).all()).all()

    def test_lower(self):
        assert ((self.reference.lower.data.S == self.df.S - self.df.S * self.df["stat. err."]).all()).all()

    def test_observe(self):
        assert ((self.reference.observe("a").data == self.df.iloc[:4]).all()).all()

    def test_query(self):
        assert ((self.reference.query("S == 1").data == self.df.iloc[0]).all()).all()
