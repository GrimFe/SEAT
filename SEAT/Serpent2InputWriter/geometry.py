import warnings

import numpy as np
import copy as cp
from .base import Entity, reformat
from .composition import Material
from dataclasses import dataclass, field

__author__ = "Federico Grimaldi"
__all__ = [
    "Universe",
    "NestedUniverse",
    "Pin",
    "Surface",
    "Cell",
    "LatticeRepresentation",
    "Lattice",
    "Geometry"
]

UniversesIncluded = {}


class ExistingUniverse(Exception):
    pass


@dataclass(slots=True)
class Universe(Entity):
    """
    Handles:
    --------
    Handles the universe definition as well as operations on a universe.

    Methods:
    --------------
    * `get_materials()`: returns a list with the materials in the Universe.

    Inherits from:
    --------------
    Entity

    Required inherited parameters:
    ------------------------------
    * `name`: string or integer - is the identity of the Serpent 2 entity
    """
    materials: list[Material] = None

    def __post_init__(self):
        self.materials = self.materials if self.materials is not None else [None]
        global UniversesIncluded
        if self.name not in UniversesIncluded.keys():
            UniversesIncluded[self.name] = self
        else:
            raise ExistingUniverse(f'Universe named {self.name} already exists.')

    def __str__(self):
        string = self.comment.__str__() + f"{self.name}"
        return string

    def get_materials(self) -> list:
        # this method is here for consistency with the inheriting classes. It is called iteratively in the nested
        # universes
        return self.materials


@dataclass(slots=True)
class NestedUniverse(Universe):
    """
    Handles universes with sub-universes.

    Methods:
    --------------
    * `add_nested()`: adds a Universe-like object to the nested universe.
    * `get_materials()`: returns a list with the materials in the Universe.

    Inherits from:
    --------------
    Universe

    Takes:
    ------
    * `daughters`: list - list of sub-universes in the nested universe.

    Required inherited parameters:
    ------------------------------
    * `name`: string or integer - is the identity of the Serpent 2 entity
    """
    daughters: list[Universe] = field(default_factory=list)

    def __iter__(self):
        return self.daughters.__iter__()

    def __str__(self):
        string = self.comment.__str__() + f"{self.name}"
        return string

    def add_universe(self, uni: [Universe]):
        if isinstance(self.daughters, list):
            self.daughters.append(uni)
        else:
            self.daughters = [uni]

    def get_materials(self) -> list:
        out = []
        for uni in self:
            out.extend(uni.get_materials())
        return out


@dataclass(slots=True)
class Pin(Universe):
    """
    Handles:
    --------
    Handles the pin definition as well as operations on single pins.

    Methods:
    --------
    * `write()`: internal method to write on a file with proper formatting for Serpent 2 input cell definition

    Class methods:
    --------------
    * `from_dict()`: creates the Pin from a dictionary coupling Material and radius.
    * `get_materials()`: returns a list with the materials in the Pin.

    Inherits from:
    --------------
    Universe

    Takes:
    ------
    * `radi`: list - couples in tuples the Material object instance and its radius (0 or `None` for external material).
        It has Material instances as first elements and floats as second elements. Default is `None`

    Required inherited parameters:
    ------------------------------
    * `name`: string or integer - is the identity of the Serpent 2 universe
    """
    radi: list[tuple[Material, float]] = None

    def __str__(self):
        string = self.comment.__str__() + f"pin {self.name}" + self.inline_comment.__str__()
        for mat, r in self.radi:
            if r == 0 or r is None:
                r_ = ''
            else:
                r_ = r
            string += f"{mat.name} {r_}\n"
        string += '\n'
        return string

    @classmethod
    def from_dict(cls, radi: dict, *args, **kwargs):
        """
        DEPRECATED - UPDATE IS NEEDED NOW THAT THE DATACLASSES ARE USED TO MAKE MATERIALS HASABLE

        Creates the Pin from a dictionary coupling Material and radius.

        Takes:
        ------
        * `radi`: dictionary - couples in the Material object instance (key) and its radius (value).
                0 or None for external material.
        """
        return cls([(k, v) for k, v in radi.items()], *args, **kwargs)

    def get_materials(self) -> list:
        return [t[0] for t in self.radi]


@dataclass(slots=True)
class Surface(Entity):
    """
    Handles:
    --------
    Handles the surface definition as well as operations on single surfaces.

    Methods:
    --------
    * `write()`: internal method to write on a file with proper formatting for Serpent 2 input surface definition.
    * `flip()`: changes the operator of the surface flipping the direction of the normal direction to it.
                Applies '-' operator on the surface complementing it.
    * `copy()`: deep copies the object.

    Inherits from:
    --------------
    Entity

    Takes:
    ------
    * `parameters`: list - list of the input parameters required by the type expressed by `kind`. Default is `None`.
    * `kind`: string - identifies the surface type, default is `'sqc'`.

    Required inherited parameters:
    ------------------------------
    * `name`: string or integer - is the identity of the Serpent 2 entity

    Private parameters:
    -------------------
    * `_operator`: string - operator to the surface. Default is `''`.
    """
    parameters: list = None
    kind: str = 'sqc'
    _operator: str = ''

    def __str__(self):
        string = self.comment.__str__()
        string += f"surf {self.name} {self.kind}"
        for p in self.parameters:
            string += f" {p}"
        string += self.inline_comment.__str__()
        return string

    def flip(self, ret=True):
        """
        Method to apply '-' operator to the surface allowing the user to consider its complement.
        Flips the direction of the normal to the surface. It happens in place.

        Takes:
        ------
        * `ret`: boolean - flag to decide whether to return the flipped object. Default is `True`
        """
        self._operator = '-'
        if ret:
            return self

    def copy(self):
        """
        Copies the object instance to another memory allocation.

        Returns:
        --------
        Returns a variable pointing to the new memory allocation
        """
        return cp.deepcopy(self)


@dataclass(slots=True)
class Cell(Entity):
    """
    Handles:
    --------
    Handles the cell definition as well as operations on single cells.

    Methods:
    --------
    * `write()`: internal method to write on a file with proper formatting for Serpent 2 input cell definition
    * `get_materials()`: returns a list with the materials in the Cell.

    Inherits from:
    --------------
    Entity

    Takes:
    ------
    * `delimiters`: list - is a list of the Surface object instances encompassing the cell. Default is `None`.
    * `father`: NestedUniverse object instance - is the universe the cell belongs to.
        The nesting of the Cell in the universe is done right after the initialisation.
        Default is `None`.
    * `kind`: string - refers to the cell type. It can either be:
        - `'fill'` - to fill the cell with a universe
        - `'material'` - to fill the cell with a material - requires passing material to the Cell
        - `'outside'` - to define outside cells
        Default is `'outside'`.
    * `filler`: Universe-like object instance - is the universe the Cell belongs to in case `kind` is `'fill'`.
        Only one single filler universe is allowed.
        Default is None.
    * `material`: Material object instance - is the cell material in case `kind` is `'material'`. Default is None.

    Required inherited parameters:
    ------------------------------
    * `name`: string or integer - is the identity of the Serpent 2 universe
    * `father`: Universe object instance - is the universe the Cell belongs to
    """

    delimiters: list[Surface] = None
    father: NestedUniverse = None  # required
    kind: str = 'outside'
    filler: Universe | None = None  # should this be a nested universe instead?
    material: Material | None = None

    def __post_init__(self):
        self._nest_to_father()

    def __str__(self):
        string = self.comment.__str__()
        string += f"cell {self.name} {self.father.name}"
        if self.kind == 'fill':
            string += f' {self.kind} {self.filler.name}'
        elif self.kind == 'material':
            string += f' {self.material.name}'
        elif self.kind == 'outside':
            string += f" {self.kind}"
        for s in self.delimiters:
            string += f' {s._operator} {s.name}'
        string += self.inline_comment.__str__()
        return string

    def _nest_to_father(self):
        global UniversesIncluded
        # The father universe will always be included in UniversesIncluded as the user is asked to pass a Universe to
        # the Cell creation. When creating the Universe, it gets also added to the UniversesIncluded
        if self.kind == 'fill':
            UniversesIncluded[self.father.name].add_universe(self.filler)
        elif self.kind == 'material':
            sub_universe_number = len(UniversesIncluded[self.father.name].daughters)
            try:
                UniversesIncluded[self.father.name].add_universe(
                    Universe(name=f'{self.father.name}.{self.material.name}', materials=[self.material]))
            except ExistingUniverse:
                UniversesIncluded[self.father.name].add_universe(
                    Universe(name=f'{self.father.name}.{self.material.name}{sub_universe_number}',
                             materials=[self.material]))
        else:  # is the exception needed here as well?
            UniversesIncluded[self.father.name].add_universe(Universe(name=f'{self.father.name}.{None}'))


@dataclass(slots=True)
class LatticeRepresentation:
    """
    Handles:
    --------
    Handles the lattice representation.

    Methods:
    --------
    * `write()`: internal method to write on a file with proper formatting for Serpent 2 input lattice definition
    * `copy()`: copies the object instance to another memory allocation
    * `transpose()`: transposes the LatticeRepresentation instance modifying it
    * `rotate()`: rotates the LatticeRepresentation instance modifying it

    Class Methods:
    --------------
    * `from_cartesian()`: creates lattice from cartesian coordinate representation: se below for more
    * `from_rows()`: creates the lattice representation from the rows of matrix in which each element is the same cell.
    * `merge()`: creates the lattice representation by merging several representations.

    Takes:
    ------
    * `rep`: 2D iterable - is the matrix representing the lattice. Each element is a Universe object instance
        located where it should be in the lattice.
    """
    rep: list[list[Universe]]

    @property
    def as_array(self):
        return self.rep if self.rep.__class__.__name__ == 'ndarray' else np.array(self.rep)

    @property
    def flatten(self) -> np.array:
        return self.as_array.flatten()

    def __post_init__(self):
        self.rep = np.array([UniversesIncluded[uni.name] for uni in self.flatten]).reshape(self.as_array.shape)

    def __iter__(self):
        return self.rep.__iter__()

    def __str__(self):
        return (reformat(str(np.array([uni.name for uni in self.flatten]).reshape(self.as_array.shape)), "[]")
                .replace('\n ', '\n') + '\n').replace("'", '')

    @classmethod
    def from_cartesian(cls, shape: tuple[int, int], filler: Universe, other: list[tuple[Universe, list[tuple]]] = None):
        """
        Creates the lattice representation from a filled matrix in which some specific elements are substituted.

        Takes:
        ------
        * `shape`: tuple - it represents the lattice shape in terms of lattice elements
        * `filler`: Cell object instance - it is the universe filling the lattice or most of it
        * `other`:  dictionary - it has Cell object instances as keys and a list of tuples containing the
            coordinates of the corresponding universe in the lattice. Coordinates can be expressed as values
            starting from 1 or strings containing one single letter

        Returns:
        --------
        * LatticeRepresentation object instance
        """
        out = np.array([filler] * np.prod(shape)).reshape(shape)
        if other is None:
            other = {}
        for k, v in other:
            for t in v:
                pos = list(t)
                if isinstance(t[0], str):
                    pos[0] = ord(t[0].lower()) - 97
                else:
                    pos[0] -= 1
                if isinstance(t[1], str):
                    pos[1] = ord(t[1].lower()) - 97
                else:
                    pos[1] -= 1
                out[(pos[1], pos[0])] = k
        return cls(out)

    @classmethod
    def from_rows(cls, universes: list[Universe], length):
        """
        Creates the lattice representation from the rows of matrix in which each element is the same universe.

        Takes:
        ------
        * `universes`: list - represents the lattice row by the cell present in that row.
        * `length`: integer - length of each row in terms of number of cells.

        Returns:
        --------
        * LatticeRepresentation object instance
        """
        out = []
        for c in universes:
            out.append([c] * length)
        return cls(out)

    @classmethod
    def merge(cls, lst):
        """
        Creates the lattice representation by merging several representations.

        Takes:
        ------
        * `lst`: list - ordered list of the LatticeRepresentations to merge.

        Returns:
        --------
        * LatticeRepresentation object instance
        """
        rep = lst[0].as_array if len(np.array(lst).shape) == 1 else lst[0]
        for lattice in lst[1:]:
            rep = np.append(rep, lattice.rep, axis=0)
        return cls(rep)

    def write(self, file: str):
        """
        Internal method to write on a file with proper formatting for Serpent 2 input lattice representation definition

        Takes:
        ------
        * `file`: string - is the name of the file where to write
        """
        with open(file, 'a') as f:
            f.write(self.__str__())

    def copy(self):
        """
        Copies the object instance to another memory allocation.

        Returns:
        --------
        Returns a variable pointing to the new memory allocation
        """
        return cp.deepcopy(self)

    def transpose(self):
        pass

    def rotate(self, deg):
        pass


@dataclass(slots=True)
class Lattice(NestedUniverse):
    """
    Handles:
    --------
    Handles the lattice as a universe with specific geometrical arrangement.

    Methods:
    --------
    * `write()`: internal method to write on a file with proper formatting for Serpent 2 input lattice definition
    * `get_materials()`: returns a list with the materials in the Lattice.

    Takes:
    ------
    * `parameters`: list - is a list with the input parameters required by the type expressed by `typ.
        Default is `None`.
    * `representation`: LatticeRepresentation object instance - is the lattice visual representation in therms of
        composing universes.
        Default is `None`.
    * `typ`: integer - is the type of the lattice, default is `1`, can be:
        * `1`: - square lattice
        * `2`: - x-hexagonal lattice
        * `3': - y-hexagonal lattice

    Required inherited parameters:
    ------------------------------
    * `name`: string or integer - is the identity of the Serpent 2 universe
    """
    parameters: list = None
    representation: LatticeRepresentation = None
    typ: int = 1

    @property
    def sub_universes(self):
        return self.representation.flatten

    def __iter__(self):
        return self.sub_universes.__iter__()

    def __str__(self):
        string = self.comment.__str__()
        string += f"lat {self.name} {self.typ}"
        for p in self.parameters:
            string += f" {p}"
        string += self.inline_comment.__str__()
        string += self.representation.__str__()
        return string

    def get_materials(self) -> list:
        """
        Returns a list with the materials in the lattice, nested universe by nested universe.
        """
        out = []
        for uni in self:
            out.extend(uni.get_materials())
        return out


@dataclass(slots=True)
class Geometry:
    """
    Handles:
    --------
    Handles the geometry section of the Serpent 2 input file, composed of pins, surfaces, cells and lattice.

    Methods:
    --------
    * `write()`: internal method to write on a file the handled items. They are written in the following order:
        - pins
        - surfaces
        - cells
        - lattice
    * `copy()`: copies the object instance to another memory allocation

    Takes:
    --------
    * `pins`: list - is a list containing the pins in the geometry
    * `surfaces`: list - is a list containing the surfaces in the geometry
    * `cells`: list - is a list containing the cells in the geometry
    * `lattices`: list - is a list containing the lattices in the geometry
    """

    pins: list[Pin]
    surfaces: list[Surface]
    cells: list[Cell]
    lattices: list[Lattice]

    def __str__(self):
        string = ''
        for p in self.pins:
            string += p.__str__()
        for s in self.surfaces:
            string += s.__str__()
        string += '\n'
        for c in self.cells:
            string += c.__str__()
        string += '\n'
        for ll in self.lattices:
            string += ll.__str__()
        string += '\n'
        return string

    def write(self, file):
        """
        Internal method to write on a file the handled items. They are written in the following order:
        - pins
        - surfaces
        - cells
        - lattice

        Takes:
        ------
        * `file`: string - is the name of the file where to write
        """
        with open(file, 'a') as f:
            f.write(self.__str__())
