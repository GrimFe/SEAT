import warnings

import numpy as np
import copy as cp
from SEAT.Serpent2InputWriter.base import Entity, reformat
from SEAT.Serpent2InputWriter.composition import Material

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
    Handles the universe definition as well as operations on a universe.
    Inherits from `SEAT.Entity`.
    
    Attributes
    ----------
    name : str | int
        the identity of the Serpent 2 entity.
    comment : `SEAT.Comment`, optional
        the comment to the Serpent entity. The default is SEAT.Comment('').
    inline_comment : `SEAT.InlineComment`, optional
        the comment to be written on the same line as the Serpent 2 entity id.
        The default is SEAT.Comment('').
    materials : list[`SEAT.Material`], optional
        the materials in the universe. The default is [].

    Methods
    -------
    assess :
        prints the `SEAT.Universe` python id.
    write :
        writes the `SEAT.Universe` to a file.
    get_materials :
        lists the materials in the Universe.

    """
    materials: list[str] = field(default_factory=list)

    def __post_init__(self):
        self.materials = self.materials if self.materials else [None]  # is this really needed?
        global UniversesIncluded
        if self.name not in UniversesIncluded.keys():
            UniversesIncluded[self.name] = self
        else:
            raise ExistingUniverse(f'Universe named {self.name} already exists.')

    def __str__(self):
        string = self.comment.__str__() + f"{self.name}"
        return string

    def get_materials(self) -> list[Material]:
        """
        Lists the materials in the universe.

        Returns
        -------
        list[`SEAT.Material`]
            the materials in the universe.

        Note
        ----
        Called iteratively in the nested universes.

        """
        return self.materials


@dataclass(slots=True)
class NestedUniverse(Universe):
    """
    Handles universes with sub-universes.
    Inherits from Universe.

    Attributes
    ----------
    name : str | int
        the identity of the Serpent 2 entity.
    comment : `SEAT.Comment`, optional
        the comment to the Serpent entity. The default is SEAT.Comment('').
    inline_comment : `SEAT.InlineComment`, optional
        the comment to be written on the same line as the Serpent 2 entity id.
        The default is SEAT.Comment('').
    materials : list[`SEAT.Material`], optional
        the materials in the universe. The default is [].
    daughters : list[Universe], optional
        the universes nested in the universe.

    Methods
    -------
    assess :
        prints the `SEAT.NestedUniverse` python id.
    write :
        writes the `SEAT.NestedUniverse` to a file.
    get_materials :
        lists with the materials in the universe.
    nest_universe :
        adds a Universe-like object to the nested universe.

    """
    daughters: list[Universe] = field(default_factory=list)

    def __iter__(self):
        return self.daughters.__iter__()

    def __str__(self):
        string = self.comment.__str__() + f"{self.name}"
        return string

    def get_materials(self) -> list:
        """
        Lists the materials in the nested universe.

        Returns
        -------
        list[`SEAT.Material`]
            the materials in the nested universe.

        Note
        ----
        Called iteratively in the nested universes.

        """
        out = []
        for uni in self:
            out.extend(uni.get_materials())
        return out

    def nest_universe(self, uni: Universe):
        """
        Nests a universe in the nested universe.

        Parameters
        ----------
        uni : Universe
            the universe to nest.

        Returns
        -------
        None.

        """
        if isinstance(self.daughters, list):
            self.daughters.append(uni)
        else:
            self.daughters = [uni]


@dataclass(slots=True)
class Pin(Universe):
    """
    Handles the pin definition as well as operations on single pins.
    Inherits from Universe.

    Attributes
    ----------
    name : str | int
        the identity of the Serpent 2 entity.
    comment : `SEAT.Comment`, optional
        the comment to the Serpent entity. The default is SEAT.Comment('').
    inline_comment : `SEAT.InlineComment`, optional
        the comment to be written on the same line as the Serpent 2 entity id.
        The default is SEAT.Comment('').
    materials : list[`SEAT.Material`], optional
        the materials in the universe. The default is [].
    radi : list[tuple[`SEAT.Material`, float]], optional
        couples the Material and its radius (0 or `None` for external material).
        The default is [].

    Methods:
    --------
    assess :
        prints the `SEAT.Pin` python id.
    write :
        writes the `SEAT.Pin` to a file.
    get_materials :
        lists the materials in the Pin.

    Class methods:
    --------------
    from_dict :
        creates the Pin from a dictionary coupling Material and radius.

    """
    radi: list[tuple[Material, float]] = field(default_factory=list)

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

        Parameters
        ----------
        * `radi`: dictionary - couples in the Material object instance (key) and its radius (value).
                0 or None for external material.
        """
        return cls([(k, v) for k, v in radi.items()], *args, **kwargs)

    def get_materials(self) -> list[Material]:
        """
        Lists the materials in the pin.

        Returns
        -------
        list[`SEAT.Material`]
            the materials in the pin.

        Note
        ----
        Called iteratively in the nested universes.

        """
        return [t[0] for t in self.radi]


@dataclass(slots=True)
class Surface(Entity):
    """
    Handles the surface definition as well as operations on single surfaces.
    Inherits from `SEAT.Entity`.
    
    Attributes
    ----------
    name : str | int
        the identity of the Serpent 2 entity.
    comment : `SEAT.Comment`, optional
        the comment to the Serpent entity. The default is SEAT.Comment('').
    inline_comment : `SEAT.InlineComment`, optional
        the comment to be written on the same line as the Serpent 2 entity id.
        The default is SEAT.Comment('').
    parameters : list[float], optional
        the input parameters required by the type expressed by `kind`. The
        default is [].
    kind : str, optional
        identifies the surface type. The default is 'sqc'.
    _operator: str, optional
        the surface oprator. The default is ''.

    Methods
    -------
    assess :
        prints the `SEAT.Surface` python id.
    write :
        writes the `SEAT.Surface` to a file.
    flip :
        changes the operator of the surface to '-' flipping its normal direction.
    copy :
        copies the object instance to another memory allocation.

    """
    parameters: list[float] = field(default_factory=list)
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
        Method to apply '-' operator to the surface allowing the user to
        consider its complement. Flips the direction of the normal to the
        surface. It happens in place.

        Parameters
        ----------
        ret : bool
            defines whether the flipped object should be returned. The default
            is True.

        Returns
        -------
        `SEAT.Surface`
            teh flipped surface if `ret` is True, else None.

        """
        self._operator = '-'
        if ret:
            return self

    def copy(self):
        """
        Copies the object instance to another memory allocation.

        Returns
        -------
        Returns a variable pointing to the new memory allocation
        """
        return cp.deepcopy(self)


@dataclass(slots=True)
class Cell(Entity):
    """
    Handles the cell definition as well as operations on single cells.
    Inherits from `SEAT.Entity`.
    
    
    Attributes
    ----------
    name : str | int
        the identity of the Serpent 2 entity.
    comment : `SEAT.Comment`, optional
        the comment to the Serpent entity. The default is SEAT.Comment('').
    inline_comment : `SEAT.InlineComment`, optional
        the comment to be written on the same line as the Serpent 2 entity id.
        The default is SEAT.Comment('').
    delimiters : list[`SEAT.Surface`]
        the surfaces encompassing the cell. The default is [].
    father : `SEAT.NestedUniverse`
        the universe the cell belongs to. The nesting of the Cell in the
        universe is done right after the initialisation. The default is `None`.
    kind : str
        identifies the cell type.
        Allowed `kind` values are:
            - 'fill': to fill the cell with a universe.
            - 'material': to fill the cell with a material; requires passing material to the Cell.
            - 'outside': to define outside cells.
        The default is 'outside'.
    filler : `SEAT.Universe`
        the universe the Cell belongs to in case `kind` is 'fill'. Only one
        single filler universe is allowed. The default is None.
    material : `SEAT.Material`
        the cell material in case `kind` is 'material'. The default is None.

    Methods
    --------
    assess :
        prints the `SEAT.Entity` python id.
    write :
        writes the `SEAT.Entity` to a file.
    get_materials :
        lists the materials in the Cell.
    _nest_to_father :
        nest the Cell to the father nested universe.

    """

    delimiters: list[Surface] = field(default_factory=list)
    father: NestedUniverse = None  # required
    kind: str = 'outside'
    filler: Universe | None = None
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
        """
        Nests the Cell to the fater universe. This is done accessing and
        adapting the UniversesIncluded global variable.

        Returns
        -------
        None.

        """
        global UniversesIncluded
        # The father universe will always be included in UniversesIncluded as the user is asked to pass a Universe to
        # the Cell creation. When creating the Universe, it gets also added to the UniversesIncluded.
        if self.kind == 'fill':
            UniversesIncluded[self.father.name].nest_universe(self.filler)
        elif self.kind == 'material':
            sub_universe_number = len(UniversesIncluded[self.father.name].daughters)
            try:
                UniversesIncluded[self.father.name].nest_universe(
                    Universe(name=f'{self.father.name}.{self.material.name}',
                             materials=[self.material]))
            except ExistingUniverse:
                UniversesIncluded[self.father.name].nest_universe(
                    Universe(name=f'{self.father.name}.{self.material.name}{sub_universe_number}',
                             materials=[self.material]))
        else:  # is the exception needed here as well?
            UniversesIncluded[self.father.name].nest_universe(Universe(name=f'{self.father.name}.{None}'))


@dataclass(slots=True)
class LatticeRepresentation:
    """
    Handles the lattice representation.
    
    Attributes
    ----------
    rep: list[list[`SEAT.Universe`]]
        the matrix representing the lattice. Each element is a `SEAT.Universe`
        located where it should be in the lattice.

    Properties
    ----------
    as_array :
        the representation universes as `np.ndarray`.
    flatten :
        the represenattion universe as 1D `np.array`.

    Methods
    -------
    write :
        writes the `SEAT.Other` instance to a file.
    copy :
        copies the object instance to another memory allocation.
    transpose :
        transposes the LatticeRepresentation instance modifying it.
    rotate :
        rotates the LatticeRepresentation instance modifying it.

    Class Methods:
    --------------
    from_cartesian :
        creates lattice from cartesian coordinate representation.
    from_rows : 
        creates the lattice representation from matrix rows.
    merge :
        creates the lattice representation merging several representations.

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
    def from_cartesian(cls, shape: tuple[int, int], filler: Universe,
                       other: list[tuple[Universe, list[tuple[int | str]]]]=[]):
        """
        Creates the lattice representation from a full matrix in which some
        specific elements are substituted.

        Parameters
        ----------
        shape : tuple[int]
            represents the lattice shape in terms of lattice elements.
        filler : `SEAT.Universe`
            the universe filling the lattice or most of it.
        other: list[tuple[`SEAT.Universe`, list[tuple[int | str]]]], optional
            couples universes with a list of coordinate where those belong to
            in the lattice. The coordinates can be expressed as values starting
            from 1 or one-letter strings. The default is [].

        Returns
        -------
        `SEAT.LatticeRepresentation`
            the instance created by the classmethod.

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
        Creates the lattice representation from the rows of matrix in which
        each element is the same universe.

        Parameters
        ----------
        universes : list[`SEAT.Universe`]
            the lattice row by the universes present in that row.
        length : int
            length of each row in terms of number of universes.

        Returns
        -------
        `SEAT.LatticeRepresentation`
            the instance created by the classmethod.

        """
        out = []
        for c in universes:
            out.append([c] * length)
        return cls(out)

    @classmethod
    def merge(cls, lst: list):
        """
        Creates the lattice representation by merging several representations.

        Parameters
        ----------
        lst : list[`SEAT.LatticeRepresentation`]
            ordered list of the LatticeRepresentations to merge.

        Returns
        -------
        `SEAT.LatticeRepresentation`
            the instance created by the classmethod.

        """
        rep = lst[0].as_array if len(np.array(lst).shape) == 1 else lst[0]
        for lattice in lst[1:]:
            rep = np.append(rep, lattice.rep, axis=0)
        return cls(rep)

    def write(self, file: str):
        """
        Writes the `SEAT.LatticeRepresentation` instance to a file.

        Parameters
        ----------
        file : str
            name of the file to write to.
        mode : str, optional
            mode to open the file. The default is 'a'.

        Returns
        -------
        None.

        """
        with open(file, 'a') as f:
            f.write(self.__str__())

    def copy(self):
        """
        Copies the object instance to another memory allocation.

        Returns:
        --------
        SEAT.Other
            a copy of the `SEAT.LatticeRepresentation` instance.

        """
        return cp.deepcopy(self)

    def transpose(self):
        pass

    def rotate(self, deg):
        pass


@dataclass(slots=True)
class Lattice(NestedUniverse):
    """
    Handles the lattice as a universe with specific geometrical arrangement.
    Inhertis from `SEAT.NestedUniverse`.

    Attributes
    ----------
    name : str | int
        the identity of the Serpent 2 entity.
    comment : `SEAT.Comment`, optional
        the comment to the Serpent entity. The default is SEAT.Comment('').
    inline_comment : `SEAT.InlineComment`, optional
        the comment to be written on the same line as the Serpent 2 entity id.
        The default is SEAT.Comment('').
    materials : list[`SEAT.Material`], optional
        the materials in the universe. The default is [].
    daughters : list[Universe], optional
        the universes nested in the universe.
    parameters : list[float], optional
        the input parameters required by the type expressed by `kind`. The
        default is [].
    representation : `SEAT.LatticeRepresentation`, optional
        the lattice visual representation in terms of composing universes. The
        default is None.
    kind : int
        the type of the lattice.
        Allowed `kind` values are:
            - 1: square lattice
            - 2: x-hexagonal lattice
            - 3: y-hexagonal lattice
        The default value is 1.

    Properties
    ----------
    sub_universes :
        the universes composing the lattice.

    Methods
    -------
    assess :
        prints the `SEAT.Lattice` python id.
    write :
        writes the `SEAT.Lattice` to a file.
    get_materials :
        lists with the materials in the universe.
    nest_universe :
        adds a Universe-like object to the nested universe.

    """
    parameters: list[float] = field(default_factory=list)
    representation: LatticeRepresentation = None
    kind: int = 1

    @property
    def sub_universes(self) -> np.array:
        """
        the universes in the Lattice.

        Returns
        -------
        `np.array`
            the universes in the Lattice.

        """
        return self.representation.flatten

    def __iter__(self):
        return self.sub_universes.__iter__()

    def __str__(self):
        string = self.comment.__str__()
        string += f"lat {self.name} {self.kind}"
        for p in self.parameters:
            string += f" {p}"
        string += self.inline_comment.__str__()
        string += self.representation.__str__()
        return string

    def get_materials(self) -> list[Material]:
        """
        Lists the materials in the Lattice.

        Returns
        -------
        list[`SEAT.Material`]
            the materials in the Lattice.

        Note
        ----
        Called iteratively in the nested universes.

        """
        out = []
        for uni in self:
            out.extend(uni.get_materials())
        return out


@dataclass(slots=True)
class Geometry:
    """
    Handles the geometry section of the Serpent 2 input file, which composes of:
        * pins
        * surfaces
        * cells
        * lattice.

    Attributes
    ----------
    pins : list[`SEAT.Pin`], optional
        the pins in the geometry. The default is [].
    surfaces : list[`SEAT.Surface`], optional
        the surfaces in the geometry. The default is [].
    cells : list[`SEAT.Cell`], optional
        the cells in the geometry. The default is [].
    lattices : list[`SEAT.Lattice`], optional
        the lattices in the geometry. The default is [].

    Methods
    -------
    write :
        writes the `SEAT.Geometry` to a file.

    """
    pins: list[Pin] = field(default_factory=list)
    surfaces: list[Surface] = field(default_factory=list)
    cells: list[Cell] = field(default_factory=list)
    lattices: list[Lattice] = field(default_factory=list)

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
        Writes the `SEAT.Geoemtry.__str__()` to a file.

        Parameters
        ----------
        file : str
            the name of the file where to write.
        mode : str, optional
            mode to open the file. The default is 'a'.

        Returns
        -------
        None.

        """
        with open(file, 'a') as f:
            f.write(self.__str__())
