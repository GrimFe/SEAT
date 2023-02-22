import numpy as np

import TargetLatticesStrings
import UniverseNames

from SEAT.Serpent2InputWriter import geometry
from SEAT.Serpent2InputWriter.composition import Material
from SEAT.Serpent2InputWriter.base import Comment

TEST_NAME = "Test"
TEST_COMMENT = Comment("Test Comment")
TEST_MATERIAL = Material(name=TEST_NAME)

# Universes that are constants defined here for their shared use in tests of LatticeRepresentation and Lattice
UNIVERSE = geometry.Universe(name=UniverseNames.GU3_LIKE)
UNIVERSE1 = geometry.Universe(name=UniverseNames.GU3_LIKE1)
UNIVERSE2 = geometry.Universe(name=UniverseNames.GU3_LIKE2)


class Test_Universe:
    uni = geometry.Universe(name=UniverseNames.UNIVERSE, comment=TEST_COMMENT)

    def test_str(self) -> None:
        assert self.uni.__str__() == TEST_COMMENT.__str__() + UniverseNames.UNIVERSE

    def test_materials(self) -> None:
        assert self.uni.materials == [None]


class Test_NestedUniverse:
    uni = geometry.NestedUniverse(name=UniverseNames.NESTED_UNIVERSE, comment=TEST_COMMENT)

    def test_str(self) -> None:
        assert self.uni.__str__() == TEST_COMMENT.__str__() + UniverseNames.NESTED_UNIVERSE

    def test_materials(self) -> None:
        assert self.uni.materials == [None]

    def test_nest_universe(self) -> None:
        extra_universe = geometry.Universe(name="Extra")
        self.uni.nest_universe(extra_universe)
        assert self.uni.daughters == [extra_universe]


class Test_Cell:
    SURFACE_NAME = 'DELIMITER'

    father = geometry.NestedUniverse(name=UniverseNames.FATHER_UNIVERSE)
    filler = geometry.Universe(name=UniverseNames.DAUGHTER_UNIVERSE)

    parameters = [0, 0, 21.5 / 2]
    delimiter = geometry.Surface(name=SURFACE_NAME, parameters=parameters)

    def test_universe_cell(self) -> None:
        cell = geometry.Cell(name=UniverseNames.CELL_UNI, father=self.father, kind='fill', filler=self.filler,
                             delimiters=[self.delimiter])
        str_uni = f"cell {cell.name} {self.father.name} fill {self.filler.name}  {self.delimiter.name}\n"
        assert cell.__str__() == str_uni

    def test_material_cell(self) -> None:
        cell = geometry.Cell(name=UniverseNames.CELL_MAT, father=self.father, kind='material', material=TEST_MATERIAL,
                             delimiters=[self.delimiter])
        str_mat = f"cell {cell.name} {self.father.name} {TEST_MATERIAL.name}  {self.delimiter.name}\n"
        assert cell.__str__() == str_mat

    def test_outside_cell(self) -> None:
        cell = geometry.Cell(name=UniverseNames.CELL_OUT, father=self.father, delimiters=[self.delimiter])
        str_out = f"cell {cell.name} {self.father.name} outside  {self.delimiter.name}\n"
        assert cell.__str__() == str_out

    def test_flip(self) -> None:
        delimiter_flip = self.delimiter.copy().flip()
        cell = geometry.Cell(name=UniverseNames.CELL, father=self.father, kind='material', material=TEST_MATERIAL,
                             delimiters=[delimiter_flip])
        str_flip = f"cell {cell.name} {self.father.name} {TEST_MATERIAL.name} - {self.delimiter.name}\n"
        assert cell.__str__() == str_flip


class Test_Surface:
    parameters = {'x0': 0, 'y0': 0, 'r': 21.5 / 2}
    surface = geometry.Surface(name=TEST_NAME, parameters=parameters)

    def test_str(self) -> None:
        assert self.surface.__str__() == "surf " + TEST_NAME + ' ' + self.surface.kind + ' ' + ' '.join(
            str(i) for i in self.parameters.values()) + '\n'

    def test_flip(self) -> None:
        surface_ = self.surface.copy()
        surface_._operator = '-'
        assert self.surface.flip() == surface_

    def test_kind_kwarg(self) -> None:
        KIND = 'cylx'
        parameters_= {'y0': 0, 'z0': 0, 'r': 21.5 / 2}
        surface_ = geometry.Surface(name=TEST_NAME, parameters=parameters_, kind=KIND)
        assert surface_.__str__() == "surf " + TEST_NAME + ' ' + KIND + ' ' + ' '.join(
            str(i) for i in self.parameters.values()) + '\n'


class Test_Pin:
    r1, r2 = (0.2, 0.4)
    TEST_PIN_STRING = f"/* Test Comment */\npin {UniverseNames.PIN}\n{TEST_NAME} {r1}\n{TEST_NAME} {r2}\n{TEST_NAME} \n\n"""
    mat1 = mat2 = mat3 = TEST_MATERIAL
    pin = geometry.Pin(name=UniverseNames.PIN, comment=TEST_COMMENT, radi=[(mat1, r1), (mat2, r2), (mat3, 0)])

    def test_str(self) -> None:
        assert self.pin.__str__() == self.TEST_PIN_STRING

    def test_materials(self) -> None:
        assert self.pin.materials == [self.mat1, self.mat2, self.mat3]

    def test_radi_kwarg(self) -> None:
        pin_ = geometry.Pin(name=UniverseNames.PIN + '_', comment=TEST_COMMENT,
                            radi=[(self.mat1, self.r1), (self.mat2, self.r2), (self.mat3, None)])
        ref = f"/* Test Comment */\npin {pin_.name}\n{TEST_NAME} {self.r1}\n{TEST_NAME} {self.r2}\n{TEST_NAME} \n\n"""
        assert pin_.__str__() == ref


class Test_LatticeRepresentation:
    control_rod_coordinates = [("C", 3), ("F", 3), ("J", 3), ("M", 3),
                               ("E", 5), ("H", 5), ("K", 5),
                               ("C", 6), ("M", 6),
                               ("E", 8), ("K", 8),
                               ("C", 10), ("M", 10),
                               ("E", 11), ("H", 11), ("K", 11),
                               ("C", 13), ("F", 13), ("J", 13), ("M", 13)]
    sample_position = [(13, 7)]

    def test_from_cartesian(self) -> None:
        representation = geometry.LatticeRepresentation.from_cartesian(shape=(15, 15), filler=UNIVERSE, other=[
            (UNIVERSE1, self.control_rod_coordinates), (UNIVERSE2, self.sample_position)])
        assert representation.__str__() == TargetLatticesStrings.GU3_LIKE

    def test_from_cartesian_letters(self) -> None:
        sample_position_ = [('M', 7)]
        representation_ = geometry.LatticeRepresentation.from_cartesian(shape=(15, 15), filler=UNIVERSE, other=[
            (UNIVERSE1, self.control_rod_coordinates), (UNIVERSE2, sample_position_)])
        assert representation_.__str__() == TargetLatticesStrings.GU3_LIKE

    def test_from_rows(self) -> None:
        representation_rows = geometry.LatticeRepresentation.from_rows([UNIVERSE] * 15, 15)  # 15x15
        assert representation_rows.__str__() == TargetLatticesStrings.GU3_LIKE.replace('S', 'T').replace('G', 'T')

    def test_merge(self) -> None:
        rep1 = geometry.LatticeRepresentation.from_rows([UNIVERSE] * 2, 15)
        rep2 = geometry.LatticeRepresentation.from_cartesian(shape=(1, 15), filler=UNIVERSE,
                                                             other=[
                                                                 (UNIVERSE1, [(3, 1), (6, 1), (10, 1), (13, 1)])])
        rep3 = geometry.LatticeRepresentation.from_rows([UNIVERSE], 15)
        rep4 = geometry.LatticeRepresentation.from_cartesian(shape=(1, 15), filler=UNIVERSE,
                                                             other=[(UNIVERSE1, [(5, 1), (8, 1), (11, 1)])])
        rep5 = geometry.LatticeRepresentation.from_cartesian(shape=(1, 15), filler=UNIVERSE,
                                                             other=[(UNIVERSE1, [(3, 1), (13, 1)])])
        rep6 = geometry.LatticeRepresentation.from_cartesian(shape=(1, 15), filler=UNIVERSE,
                                                             other=[(UNIVERSE2, [(13, 1)])])
        rep7 = geometry.LatticeRepresentation.from_cartesian(shape=(1, 15), filler=UNIVERSE,
                                                             other=[(UNIVERSE1, [(5, 1), (11, 1)])])
        merged = geometry.LatticeRepresentation.merge(
            [rep1, rep2, rep3, rep4, rep5, rep6, rep7, rep3, rep5, rep4, rep3, rep2, rep1])
        assert merged.__str__() == TargetLatticesStrings.GU3_LIKE

    def test_flatten(self) -> None:
        representation_small = geometry.LatticeRepresentation.from_rows([UNIVERSE] * 2, 2)
        assert (representation_small.flatten == np.array([UNIVERSE] * 4)).all()


class Test_Lattice:
    simple_universe = geometry.Universe(name=UniverseNames.SUB_UNIVERSE, _materials=[TEST_MATERIAL],
                                        comment=TEST_COMMENT)
    REPRESENTATION = geometry.LatticeRepresentation.from_rows([simple_universe, simple_universe], 2)

    def test_str(self) -> None:
        control_rod_coordinates = [("C", 3), ("F", 3), ("J", 3), ("M", 3),
                                   ("E", 5), ("H", 5), ("K", 5),
                                   ("C", 6), ("M", 6),
                                   ("E", 8), ("K", 8),
                                   ("C", 10), ("M", 10),
                                   ("E", 11), ("H", 11), ("K", 11),
                                   ("C", 13), ("F", 13), ("J", 13), ("M", 13)]
        sample_position = [(13, 7)]
        representation = geometry.LatticeRepresentation.from_cartesian(shape=(15, 15), filler=UNIVERSE, other=[
            (UNIVERSE1, control_rod_coordinates),
            (UNIVERSE2, sample_position)])
        lattice = geometry.Lattice(name=UniverseNames.LATTICE, parameters=[0, 0, 15, 15, 1.43],
                                   representation=representation,
                                   kind=1)
        assert lattice.__str__() == "lat " + UniverseNames.LATTICE + " 1 0 0 15 15 1.43\n" + TargetLatticesStrings.GU3_LIKE

    def test_sub_universes(self) -> None:
        lattice = geometry.Lattice(name=UniverseNames.LATTICE1, representation=self.REPRESENTATION)
        assert lattice.sub_universes[0] == self.simple_universe
        assert np.alltrue(lattice.sub_universes == np.array([self.simple_universe] * 4))
        assert (lattice.sub_universes == np.array([self.simple_universe] * 4)).all()

    def test_materials(self) -> None:
        lattice = geometry.Lattice(name=UniverseNames.LATTICE2, representation=self.REPRESENTATION)
        assert (lattice.materials == np.array([TEST_MATERIAL] * 4)).all()


class Test_Geometry:
    father = geometry.NestedUniverse(name=UniverseNames.GEOMETRY_FATHER)
    material = Material(name='TEST_MATERIAL')
    external_material = Material(name='EXTERNAL')
    r1, r2 = 0.1, None
    pin = geometry.Pin(name=UniverseNames.GEOMETRY_PIN, radi=[(material, r1), (external_material, r2)])
    p1, p2, p3 = 0, 0, r1
    surface = geometry.Surface(name='SURFACE', parameters={'x0': p1, 'y0': p2, 'r': p3})
    cell = geometry.Cell(name=UniverseNames.GEOMETRY_CELL, father=father, delimiters=[surface], kind='material',
                         material=material)
    n_side_pins = 2
    representation = geometry.LatticeRepresentation.from_cartesian(shape=(n_side_pins, n_side_pins), filler=father)
    c1, c2 = 0, 0
    lattice = geometry.Lattice(name=UniverseNames.GEOMETRY_LATTICE, parameters=[0, 0, n_side_pins, n_side_pins, r1],
                               representation=representation)
    geom = geometry.Geometry(pins=[pin], surfaces=[surface], cells=[cell], lattices=[lattice])

    GEOMETRY_STRING = f"pin PIN\n{material.name} {r1}\n{external_material.name} \n\n" + \
                      f"surf {surface.name} sqc {p1} {p2} {p3}\n\n" + \
                      f"cell {cell.name} {father.name} {material.name}  {surface.name}\n\n" + \
                      f"lat {lattice.name} 1 {c1} {c2} {n_side_pins} {n_side_pins} {r1}\n" + \
                      f"{father.name} {father.name}\n{father.name} {father.name}\n\n"

    def test_string(self) -> None:
        assert self.geom.__str__() == self.GEOMETRY_STRING
