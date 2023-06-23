# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 14:05:36 2023

@author: fgrimald
"""
import pandas as pd
import numpy as np

import SEAT.Serpent2UncertaintyPropagation.MomentumPropagation.CovarianceSuite as CS

class Test_Functions:
    idx = pd.MultiIndex.from_tuples([("9228", "18", "(0, 1e-1]"),
                                     ("9228", "18", "(1e-1, 5.4e-1]"),
                                     ("9228", "18", "(5.4e-1, 4]")],
                                    names=["A", "B", "C"])
    cov = pd.DataFrame(np.zeros((3,3)), index=idx, columns=idx)

    def test_get_cov(self):
        idx = pd.MultiIndex.from_tuples([(9228, 18, 1e-1),
                                         (9228, 18, 5.4e-1),
                                         (9228, 18, 4)],
                                        names=['MAT', 'MT', 'E [eV]'])
        assert (CS.get_cov(self.cov).index == idx).all()
        assert (CS.get_cov(self.cov).columns == idx).all()

    def test_reindex_cov(self):
        idx = pd.MultiIndex.from_tuples([("U235", "(z,fission)", 1e-7),
                                         ("U235", "(z,fission)", 5.4e-7),
                                         ("U235", "(z,fission)", 4e-6)],
                                        names=['N', 'MT', 'E [eV]'])
        assert (CS.reindex_cov(CS.get_cov(self.cov)).index == idx).all()
        assert (CS.reindex_cov(CS.get_cov(self.cov)).columns == idx).all()

class Test_Covariance:
    def test_from_csv(self):
        pass

    def test_idx(self):
        pass
