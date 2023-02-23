# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 10:55:19 2023

@author: fgrimald
"""

import SEAT.lattice_functions as lf

class TestFunctions:
    A, B, C, D, E = range(5)
    def test_NxNy_params(self):
        assert lf._NxNy_params(self.A) == '0 0'
        assert lf._NxNy_params(self.A, self.B, self.C, self.D, self.E) == '0 1 2 3 4'

    def test_square_params(self):
        assert lf.square_params(self.A) == '0 0'
        assert lf.square_params(self.A, self.B, self.C, self.D, self.E) == '0 1 2 3 4'