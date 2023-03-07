# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 10:02:22 2023

@author: fgrimald
"""
from SEAT.Serpent2InputWriter.base import Comment
from SEAT.Serpent2InputWriter.composition import MaterialRepresentation

EMPTY_COMMENT = Comment('')
COMMENT_BASE_DCT: dict[str, Comment] = {'Intro': EMPTY_COMMENT,
                                        'Geometry': EMPTY_COMMENT,
                                        'Steps': EMPTY_COMMENT,
                                        'Materials': EMPTY_COMMENT,
                                        'Others': EMPTY_COMMENT,
                                        'Concluding': EMPTY_COMMENT}
POPULATION: dict[str, int] = {'particles': 100000,
                              'generations': 250,
                              'inactive': 30,
                              'k_guess': 1,
                              'batch_interval': 1,
                              'parallel': 1}
URES: dict[str, int] = {"active": False,
                        "representation": None,
                        "dilution_cut": 1e-9}
