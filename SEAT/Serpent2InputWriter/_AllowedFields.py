# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 10:05:04 2023

@author: fgrimald
"""

PREDICTOR_CORRECTOR = {'CE': 0,
                       'CELI': 1,
                       'LE': 2,
                       'LELI': 3,
                       'LEQI': 4,
                       'CECE': 6,
                       0: 0,
                       1: 1,
                       2: 2,
                       3: 3,
                       4: 4,
                       5: 5,
                       6: 6}
BC = {'vacuum': 1,
      'reflective': 2,
      'periodic': 3,
      'V': 1,
      'R': 2,
      'P': 3,
      1: 1,
      2: 2,
      3: 3}
XSCALC = {None: None,
          'NO': None,
          'DIRECT': 1,
          'SPECTRUM': 2,
          1: 1,
          2: 2}
BUMODE = {"CRAM": 2,
          "TTA": 1,
          2: 2,
          1: 1}
PLOT_KIND = {'yz': 1,
             'xz': 2,
             'xy': 3}