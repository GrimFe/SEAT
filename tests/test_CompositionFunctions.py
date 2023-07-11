# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 18:51:24 2023

@author: fgrimald
"""
import SEAT.natural
from SEAT.composition_functions import *

class Test_Functions:
    def test_are_elements(self) -> None:
        assert are_elements({'Am', 'U'})
        assert not are_elements({'Am241m', 'U'})

    def test_are_nuclides(self) -> None:
        assert are_nuclides({'Am241m', 'U235'})
        assert not are_nuclides({'Am241m', 'U'})

    def test_unfold_composite(self) -> None:
        u = unfold_composite(SEAT.composites.H2O, atomic=True)
        assert u == {'H1': 0.66659, 'H2': 7.666666666666667e-05, 'O16': 0.3325233333333333, 'O17': 0.00012666666666666666, 'O18': 0.0006833333333333333}

    # def test_get_existing_xs(self) -> None:
    #     assert get_existing_xs('endfb71')[0] == 'Ba138'

    def test_pollute(self) -> None:
        p1 = 7e-2
        p2 = 3e-2
        test1 = pollute(SEAT.natural.U, {'Pu239': p1 + p2})
        assert test1 == {'Pu239': 0.1, 'U238': 0.893466, 'U235': 0.00648, 'U234': 4.5e-05}
        assert sum(test1.values()) == 0.999991
        test2 = pollute(SEAT.natural.U, {'Pu239': p1, 'Pu240': p2})
        assert test2 == {'Pu239': 0.07, 'Pu240': 0.03, 'U238': 0.893466, 'U235': 0.00648, 'U234': 4.5e-05}
        assert sum(test2.values()) == 0.999991

    def test_enrich(self) -> None:
        e = 3e-2
        test1 = enrich(SEAT.natural.U, {'U235': e})
        assert test1 == {'U238': 0.9699411482790922, 'U235': 0.03, 'U234': 4.885172090774484e-05}
        assert sum(test1.values()) == 0.9999899999999999
        test2 = enrich(SEAT.natural.U, {'U235': e, 'U234': 0})
        assert test2 == {'U238': 0.9699899999999999, 'U235': 0.03, 'U234': 0}
        assert sum(test2.values()) == 0.9999899999999999
        test3 = enrich(SEAT.natural.U, {'U235': e, 'U234': 1e-2})
        assert test3 == {'U238': 0.9599899999999999, 'U235': 0.03, 'U234': 0.01}
        assert sum(test3.values()) == 0.9999899999999999
        test4 = enrich({'U238': 0.5, 'U235': 0.3, 'H1': 0.2}, {'U235': 0})
        assert test4 == {'U238': 0.8, 'U235': 0, 'H1': 0.2}
        assert sum(test4.values()) == 1
        components = {'U238': 0.3309133333333333, 'U235': 0.0024, 'U234': 1.6666666666666667e-05, 'O16': 0.6650466666666666, 'O17': 0.00025333333333333333, 'O18': 0.0013666666666666666}
        target = {'U238': 0.30331472335539233, 'U235': e, 'U234': 1.5276644607620947e-05, 'O16': 0.6650466666666666, 'O17': 0.00025333333333333333, 'O18': 0.0013666666666666666}
        test5 = enrich(components, {'U235': 0.03})
        assert test5 == target
        assert test5['U235'] == e
        assert test5['U238'] < SEAT.natural.U['U238']
        assert test5['O16'] == components['O16']
