# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 17:32:33 2023

@author: fgrimald
"""

import SEAT.surface_functions as sf

class TestFunctions:
    A, B, C, D, E, F, G, H, I = range(9)
    def test_p_params(self):
        assert sf._p_params(self.A) == '0'

    def test_px_params(self):
        assert sf.px_params(self.A) == '0'

    def test_py_params(self):
        assert sf.py_params(self.A) == '0'

    def test_pz_params(self):
        assert sf.pz_params(self.A) == '0'

    def test_plane_params_equaiton(self):
        assert sf.plane_params(self.A, self.B, self.C, self.D) == '0 1 2 3'

    def test_plane_params_points(self):
        assert sf.plane_params(c1=(self.A, self.B, self.C),
                               c2=(self.D, self.E, self.F),
                               c3=(self.G, self.H, self.I)) == '0 1 2 3 4 5 6 7 8'

    def test_cylinder_params(self):
        assert sf._cylinder_params(self.A, self.B, self.C, self.D, self.E) == '1 2 0 3 4'
        assert sf._cylinder_params(self.A, self.B, self.C, self.D) == '1 2 0 3'

    def test_sqc_params(self):
        assert sf.sqc_params(self.A, self.B, self.C, self.D, self.E, self.F) == '1 2 0 3 4 5'

    def test_cyl_params(self):
        assert sf.cyl_params(self.A, self.B, self.C, self.D, self.E) == '1 2 0 3 4'

    def test_cylx_params(self):
        assert sf.cylx_params(self.A, self.B, self.C, self.D, self.E) == '1 2 0 3 4'

    def test_cyly_params(self):
        assert sf.cyly_params(self.A, self.B, self.C, self.D, self.E) == '1 2 0 3 4'

    def test_cylz_params(self):
        assert sf.cylz_params(self.A, self.B, self.C, self.D, self.E) == '1 2 0 3 4'
