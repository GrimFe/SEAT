import copy as cp
import warnings
from dataclasses import dataclass

import numpy as np

from SEAT.nuclides import zam2za, nuclide2za
from .base import Entity, reformat


__author__ = "Federico Grimaldi"
__all__ = [
    "MaterialComposition",
    "MaterialRepresentation",
    "Material",
    "Composition",
]


@dataclass(slots=True)
class MaterialComposition:
    """
    Handles the nuclide composition of the materials in the simulation.

    Methods:
    --------
    * `to_file()`: writes to a file the components in the material composition

    Class methods:
    --------------
    * `from_file()`: reads the material composition from a file
    * `from_za()`: creates the material composition with nuclides represented through their ZA number
    * `from_nuclides()`: creates the material composition with nuclides represented through their string identifier
    * `write()`: internal method to write on a file with proper formatting for Serpent 2 input material composition
    * `from_output_file()`: creates the material composition from a Serpent2 output file - not implemented yet

    Takes:
    ------
    * `components`: str - is the nuclide components of a material
    """
    components: str

    def __str__(self):
        return self.components

    @classmethod
    def from_file(cls, file: str):
        """
        Creates material components from a file

        Takes:
        ------
        * `file`: string - name of the file from which the material composition should be read
        """
        with open(file, 'r') as f:
            components = f.read()
        return cls(components)

    @classmethod
    def from_za(cls, za: dict[int, float], atomic: bool = True):
        """
        Creates material composition from a dictionary

        Takes:
        ------
        * `za` is a dictionary having the ZA number as keys and the fraction (as positive floats) as values
        * `atomic` is a flag to select whether the passed fractions are atomic or mass fractions.
                    Default is True corresponding to atomic fractions.
                    Serpent2 uses negative fraction values to identify mass fractions.
        """
        za_ = za if atomic else {k: -v for k, v in za.items()}
        components = reformat(str(za_), "{}:").replace(', ', '\n')
        return cls(components)

    @classmethod
    def from_zam(cls, zam: dict[int, float], **kwargs):
        """
        Creates material composition from a dictionary

        Takes:
        ------
        * `zam` is a dictionary having the ZAM number as keys and the fraction (as positive floats) as values
        """
        components = dict(zip([zam2za(z)[0] for z in zam.keys()], zam.values()))
        return cls.from_za(components, **kwargs)

    @classmethod
    def from_nuclides(cls, nuclides: dict[str, float], **kwargs):
        """
        Creates material composition from a dictionary

        Takes:
        ------
        * `nuclides`: dictionary - has the nuclides as keys and the composition (as positive floats) as values.
            Nuclides shall be in the form Nn-###, e.g.: U-235, H-1, Am-241, Am-241m to have them represented as strings
            in the Serpent2 simulation input file, other options which would get converted to ZA are:
            * NN###
            * nn###
            * Nn###
            All keys in the dictionary should either have the `'-'` or not have it.
        """
        if '-' in list(nuclides.keys())[0]:
            components = reformat(str(nuclides), "{}:'").replace(', ', '\n')
            return cls(components)
        else:
            components = dict(zip([nuclide2za(i.lower().capitalize())[0] for i in nuclides.keys()],
                                  nuclides.values()))
            return cls.from_za(components, **kwargs)

    @classmethod
    def from_output_file(cls, file: str):  # interaction with serpentTools missing
        """
        Creates material composition from a Serpent 2 output file

        Takes:
        ------
        * `file`: string - name of the file from which the material components should be read
        """
        pass

    def to_file(self, file: str):
        """
        Writes the material composition to a file.

        Takes:
        ------
        * `file`: string - is the name of the file in which the material components should be witter
        """
        with open(file, 'w') as f:
            f.write(self.__str__())

    def copy(self):
        """
        Copies the object instance to another memory allocation.

        Returns:
        --------
        Returns an object instance pointing to the new memory allocation
        """
        return cp.deepcopy(self)


@dataclass(slots=True)
class MaterialRepresentation:
    """
    Handles the material representation coupling its components to the library identifier and the temperature at which
    nuclear data should be generated.

    Methods:
    --------
    * `write()`: internal method to write on a file with proper formatting for Serpent 2 input material composition

    Class methods:
    --------------
    * `from_string()`: creates material representation from a string where no information on the string nor on the
        library is included
    * `from_material_composition()`: creates material composition from a `MaterialComposition` object instance

    Takes:
    ------
    * `components`: dictionary - is the components of the material. It has fractions as values (positive values for
        nuclide fractions, negative values for mass fractions) and nuclide identifiers as keys. The keys can be in
        the form:
            - ZA number (integer)
            - Nn-###, e.g.: U-235, H-1, Am-241, Am-241m (string)
    * `tmp`: float - is the temperature at which nuclear data should be generated [K]
    * `data_type`: string - is the library identifier. Default is `'c'` for cross sections
    """
    components: dict[any, float]
    tmp: float
    data_type: str = 'c'

    def __str__(self):
        string = ''
        for k, v in zip(self.get_nuclide_representation(), self.get_values()):
            string += f"{k} {v}\n"
        return string

    @classmethod
    def from_material_composition(cls, components: MaterialComposition, *args, **kwargs):
        dct = {c.split()[0]: c.split()[1] for c in components.components.split('\n')}
        return cls(dct, *args, **kwargs)

    @classmethod
    def from_string(cls, string: str, *args, **kwargs):
        dct = {c.split()[0]: c.split()[1] for c in string.split('\n')}
        return cls(dct, *args, **kwargs)

    def get_temperature(self) -> str:
        if self.tmp < 1000:
            tmp_ = '0' + '{:.0f}'.format(self.tmp / 100)
        else:
            tmp_ = '{:.0f}'.format(self.tmp / 100)
        return tmp_

    def get_nuclide_representation(self):
        temperature, data_type = self.get_temperature(), self.data_type
        return [f"{k}.{temperature}{data_type}" for k in self.components.keys()]

    def get_values(self):
        return [v for v in self.components.values()]

    def write(self, file):
        """
        Internal method to write on a file with proper formatting for Serpent 2 input material definition

        Takes:
        ------
        * `file`: string - is the name of the file where to write
        """
        with open(file, 'a') as f:
            f.write(self.__str__())


@dataclass(slots=True)
class Material(Entity):
    """
    Handles:
    --------
    Handles the material definition and operations on single materials.

    Methods:
    --------
    * `write()`: internal method to write on a file with proper formatting for Serpent 2 input material definition
    * `divide()`: computes and creates the desired material divisions
    * `get_temperature()`: gets the temperature that is not none among the many that can be passed

    Class Methods:
    --------------
    * `mix()`: allows mixing materials to one that is the mixture to be written in the Serpent 2 input

    Inherits from:
    --------------
    Entity

    Takes:
    ------
    * `dens`: float - is the material density (positive values for nuclide density [cm-3] and negative values for
        mass densities [g/cm3]).
    * `representation`: MaterialComposition - is the material composition representation. Default is `None`
    * `tmp`: float - is the material temperature [K] for doppler preprocessing. Default is `None`
    * `tms`: float - is the material temperature [K] for on-the-fly temperature treatment. Default is `None`
    * `tft`: tuple - contains two floats for the material temperature limits [K]. Default is `None`
    * `rgb`: tuple - contains three integers in range 0-255 for color channels in the plot:
        red, green, blue respectively. Default is `None`
    * `vol`: float - material volume or cross-section (for 2D calculations) [cm3 (cm2)] calculations
    * `mass`: float - material mass [g]
    * `burn`: integer - indicates the regions to burn in the material. Default is 0, better use of divisions is handled
        by the `divide()` method
    * `fix`: tuple - indicates the library information for decay nuclides (i.e. without cross-section data).
        It is composed of:
            * string - library id
            * float - temperature [K]
        Default is `None`
    * `moder`: list - lists the moderators in the material and uses thermal scattering library for that nuclides.
        Each element in the list is a tuple composed as follows:
            - string - name of the thermal scattering data library defined using the therm card
            - integer - ZA number of the thermal scatter (1001 for H1)
        Default is `None`
    * `_mixed`: bool - identifies the material as a mixed one. Intended for internal use. Default is False.
    * `_mixture`: string - String of the mixed materials and their fraction. Default is `None`

    Required inherited parameters:
    ------------------------------
    * `name`: string or integer - is the identity of the Serpent 2 entity

    Default internal parameters:
    ----------------------------
    * `divisions`: list - list of material divisions. Default is `[]`.
    * `division_string`: string - string to set the material division. Default is ''
    """
    dens: float | None = None
    representation: MaterialRepresentation = None
    tmp: float = None
    tms: float = None
    tft: tuple[float, float] = None
    rgb: tuple[int, int, int] = None
    vol: float = None
    mass: float = None
    burn: int = 0
    fix: tuple = None
    moder: list = None
    _mixture: str = None
    _divisions: list = None
    _division_string: str = ''

    def __str__(self):
        string = self.comment.__str__()
        if self._mixture is not None:
            string += f"mix {self.name}" + self.inline_comment.__str__() + f"{self._mixture}"
        else:
            string += f"mat {self.name} {self.dens} {self.get_temperature_kind()} {self.get_temperature()}"
        # Check which cards work for the mix material and relate it to the transfer of kwargs in its definition
        if self.rgb is not None:
            string += f" rgb {self.rgb[0]} {self.rgb[1]} {self.rgb[2]}"
        if self.vol is not None:
            string += f" vol {self.vol}"
        if self.mass is not None:
            string += f" mass {self.mass}"
        if self.burn:
            string += " burn 1"
        if self.fix:
            string += f" fix {self.fix[0]} {self.fix[1]}"
        if self.moder is not None:
            string += reformat(f" moder {self.moder}", ",'()[]{}")
        string += self.inline_comment.__str__()
        if self._mixture is None:
            string += self.representation.__str__()
        string += "\n"
        return string

    @classmethod
    def mix(cls, name: str, materials: list[tuple], *args, **kwargs) -> object:
        """
        Creates material mixing other materials

        Takes:
        ------
        * `name`: string or integer - is the identity of the Serpent 2 entity
        * `materials`: list of tuples - having Material object instances as first elements and fractions following
        """
        a = []
        for k, v in materials:
            a.append(f"{k.name} {v}\n")
        return cls(name=name, _mixture=''.join(a), *args, **kwargs)

    def get_temperature(self) -> str:
        multiple = False
        temperature = None
        if self.tms is not None:
            temperature = str(self.tms)
        if self.tmp is not None:
            if temperature is not None:
                multiple = True
            temperature = str(self.tmp)
        if self.tft is not None:
            if temperature is not None:
                multiple = True
            temperature = reformat(str(self.tft), "(),'")
        if temperature is None or multiple:
            warnstring = f"No temperature or multiple temperatures given to {self.name}.\n" + \
                         f"Temperature will be written as if it was the same as the material representation temperature"
            warnings.warn(warnstring)
            temperature = self.representation.tmp
        return temperature

    def get_temperature_kind(self):
        kind = None
        if self.tms is not None:
            kind = 'tms'
        if self.tmp is not None:
            kind = 'tmp'
        if self.tft is not None:
            kind = 'tft'
        # if double:
        #     raise Exception(f'More than one temperature type was given to Material {self.name}')
        return kind

    def write(self, file: str):
        """
        Internal method to write on a file with proper formatting for Serpent 2 input material definition

        Takes:
        ------
        * `file`: string - is the name of the file where to write
        """
        with open(file, 'a') as f:
            f.write(self.__str__())

    def divide(self, lvl: int, kind: str, sub: tuple):
        """
        Method to create material divisions and compute sub-volumes.

        Takes:
        ------
        * `lvl`: integer - indicates the cell level of the division
            (0: no division, 1: last level division, 2: second last level division, ...)
        * `kind`: string - indicates the type of division to implement. Allowed values are:
            * `'x_cartesian'`
            * `'y_cartesian'`
            * `'z_cartesian'`
            * `'radial'`
            * `'sectorial'`
        * `sub`: tuple - it contains float data for the division (lengths in [cm]).
            It can be composed as follows:
            * for  cartesian and radial division:
                `sub` = (N, min, max), where N is the number of zones and min and max are respectively
                                the minimum an maximum coordinates in the direction defined by `kind`
            * for radial division it can also be:
                `sub` = (N, r0, r1, .., rn), where N is the number of zones and ri is the radial
                                coordinate of the i-th division
            * for sectorial division:
                `sub` = (N, s0), where N is the number of zones and s0 the angular coordinate of
                                the first in [deg]

        Global variables:
        -----------------
        * `mvol_counter`

        Notes:
        ------
        * mvol calculation not implemented for cartesian and sectorial division yet
        * Only 2D subdivision is implemented
        """
        self._division_string = f"div {self.name} sep {lvl} sub{kind[0]} " + reformat(str(sub), "(),")
        if self._divisions is None:
            self._divisions = []
        # mvol calculation
        if 'cartesian' in kind:
            area = []
        elif kind == 'radial' and len(sub) == 3:
            area = sub[-1] ** 2 * np.pi / sub[0] * np.ones(sub[0])
        elif kind == 'radial' and len(sub) > 3:
            area_out = np.asarray(sub[1:]) ** 2 * np.pi
            r_in = [0]
            r_in.extend(sub[1:-1])
            area_in = np.asarray(r_in) ** 2 * np.pi
            area = area_out - area_in
        elif kind == 'sectorial':
            area = []
        for a in area:
            self._divisions.insert(0, Division(self, a))


@dataclass(slots=True)
class Division:
    """
    Handles:
    --------
    Handles the cell division definition as well as operations on single divisions.

    Methods:
    --------
    * `write()`: internal method to write on a file with proper formatting for Serpent 2 input surface definition
    * `copy()`: deep copies the object.

    Takes:
    ------
    * `material`: Material object instance - is material of the division
    * `volume`: float - is the volume of the division
    """
    material: Material
    volume: float

    def __str__(self):
        string = f"{self.material.name} {self.volume}"
        return string

    def write(self, file: str):
        """
        Internal method to write on a file with proper formatting for Serpent 2 input division definition

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


@dataclass(slots=True)
class Composition:
    """
    Handles:
    --------
    Handles the composition section of the Serpent 2 input file, which composes of:
        * nuclear data libraries link
        * unresolved resonance probability table
        * thermal scattering data
        * restart files

    Methods:
    --------
    * `write()`: internal method to write on a file the handled items. They are written in the following order:
        - material
        - libraries
        - thermal scattering data
        - restart files
    * `get_subdivisions()`: gets the subdivisions of the materials in the composition

    Takes:
    ------
    * `materials`: list - list of Material object instances to introduce in the Composition
    * `libraries`: dictionary - has Serpent library identifiers as keys and raw string paths as values
    * `scattering`: dictionary - includes a tuple with the interpolation type and the name of the library as key
        and tuples as values. The interpolation type can either be `''` or `'stoch'`, for stochastic interpolation.
        The values are composed of a float indicating the temperature and of a list of strings identifying the
        libraries. A temperature value of 0 induces an interpolation of the libraries whereas a temperature set to
        `None` or negative values write no temperature.
        Default is `None`.
    * `to_restart`: string - the name of the binary restart file to which the composition should be written.
        Default is `None`.
    * `from_restart`: string - the name of the binary restart file from which the composition should be read.
        Default is `None`.
    """
    materials: list[Material]
    libraries: dict[str, str]
    scattering: dict[tuple[str, str], tuple[float, list]] = None
    to_restart: str = None
    from_restart: str = None

    def __str__(self):
        string = ''
        for k, v in self.libraries.items():
            string += f"{k} '{v}'\n"
        string += '\n'
        if self.scattering is not None:
            for k, v in self.scattering.items():
                type_ = k[0]
                name = k[1]
                interpolation = str(v[0]) if v[0] is not None and v[0] >= 0 else ''
                lib = reformat(str(v[1:]), "[](),'")
                string += f'therm{type_} {name} {interpolation} {lib}\n'
        if self.to_restart is not None or self.from_restart is not None:
            string += '\n'
            string += "/* Composition restart file definition */\n"
            string += f"rfw {self.to_restart}\n" if self.to_restart is not None else ''
            string += f"rfw {self.from_restart}\n" if self.from_restart is not None else ''
        string += '\n'
        for m in self.materials:
            string += m.__str__()
        string += '\n'
        return string

    def write(self, file: str):
        """
        internal method to write on a file the handled items. They are written in the following order:
        - material
        - libraries
        - thermal scattering data
        - restart files

        Takes:
        ------
        * `file`: string - is the name of the file where to write
        """
        with open(file, 'a') as f:
            f.write(self.__str__())

    def get_subdivisions(self) -> list:
        """
        Transforms the subdivision lists in the composition composing materials in one single list for the whole
        composition.

        Returns:
        --------
        `out`: list - contains the subdivisions
        """
        out = []
        for m in self.materials:
            out.extend(m._divisions)
        return out
