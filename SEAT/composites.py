# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 17:54:12 2023

Materials from:
    McConn, Ronald J, Gesh, Christopher J, Pagh, Richard T, Rucker, Robert A, and Williams, III, Robert.
    Compendium of Material Composition Data for Radiation Transport Modeling.
    United States: N. p., 2011.
    Web. doi:10.2172/1023125.

@author: fgrimald
"""

## Air: dry near seal level
Air = {'C': 0.000150, 'N': 0.784431, 'O': 0.210784, 'Ar': 0.004671}

## Stainless Steel
SS_HT9 = {'C': 0.002, 'Si': 0.004, 'P': 0.000534, 'S': 0.000344, 'Va': 0.003, 'Cr': 0.121971, 'Mn': 0.006023, 'Fe': 0.838897, 'Ni': 0.007698, 'Mo': 0.005748, 'W': 0.0015}

## Zircaloy
Zircaloy4 = {'O': 0.00679, 'Cr': 0.001741, 'Fe': 0.003242, 'Zr': 0.977549, 'Sn':0.010677}

## Stoichiometry
UO2 = {'U': 1/3, 'O': 2/3}
H2O = {'H': 2/3, 'O': 1/3}
B4C = {'B': 4/5, 'C': 1/5}

# Elements
Bi = {'Bi': 1}