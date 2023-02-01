import copy as cp
import warnings
from dataclasses import dataclass, field

import numpy as np

import SEAT.natural
import SEAT.composites
from SEAT.composition_functions import *

from SEAT.nuclides import zam2nuclide, nuclide2zam, nuclide2za
from SEAT.Serpent2InputWriter.base import Entity, reformat

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

    Args:
    -----
    components : dict[str, float]
        * key: nuclide symbol (`'Nn###'`)
        * value: densities as positive floats.
    atomic : bool
        makes the material densities be interpreted as atomic or mass.
    _already_za : bool
        set to true when `components.keys()` was made of ZAs.

    Methods:
    --------
    write :
        writes the `SEAT.MaterialComposition` instance to a file.
    adjust :
        adjusts the material composition to a nuclear data library.

    Class methods:
    --------------
    parse :
        reads the material composition from a string.
    _from_list :
        reads the material composition from a list of strings.
    from_file :
        reads the material composition from a file.
    from_output_file :
        readsthe material composition from a Serpent2 output file - not implemented yet.
    enriched_custom :
        enriches a custom material composition.
    enriched :
        enriches a predefined material.
    polluted_custom :
        creates the material composition from a pollution criterion.
    polluted :
        creates the material composition polluting a natural element or one form `SEAT.composites`.
    natural :
        creates the material composition form natural abundances.
    predefined :
        creates the material composition form those predefined in `SEAT.composites`.
    composite :
        creates the material composition form elements.
    from_zam :
        creates the material composition giving the `components.keys()` as ZAMs.
    from_za :
        creates the material composition giving the `components.keys()` as ZAs.
    from_sep_nuclides :
        creates the material composition giving the `components.keys()` as (`Nn{sep}###`).

    """
    components: dict[str, float]
    atomic: bool
    _already_za: bool = field(init=False, default=False)

    def __str__(self):
        if not self._already_za:
            factor = 1 if self.atomic else -1
            nuclides = {nuclide2zam(k.lower().capitalize()): factor * v for k, v in self.components.items()}
            string = reformat(str(nuclides).replace(',', '\n'), '{} ').replace(':', ' ')
        else:
            za = self.components
            string = reformat(str(za).replace(',', '\n'), '{} ').replace(':', ' ')
        return string

    @classmethod
    def parse(cls, string: str, separator: str=None, **kwargs):
        """
        Reads the material composition from a string.

        Parameters
        ----------
        string : str
            the string to parse.
        separator : str, optional
            the separator of the nuclide symbol and the corresponding density.
            The default is None.
        **kwargs :
            - atomic : bool
                makes the material densities be interpreted as atomic or mass.

        Returns
        -------
        SEAT.MaterialComposition
            the instance created by he classmethod.

        Note
        ----
        The nuclides in `string` should either be in the form 'Nn###', or ZAM.

        """
        return cls._from_list(string.split('\n'), separator, **kwargs)

    @classmethod
    def _from_list(cls, lines: list[str], separator: str=None, **kwargs):
        """
        Reads the material composition from a list of strings.

        Parameters
        ----------
        lines : list[str]
            the list of strings with the composition and the density.
        separator : , optional
            the separator of the nuclide symbol and the corresponding density.
            The default is None.
        **kwargs :
            - atomic : bool
                makes the material densities be interpreted as atomic or mass.

        Returns
        -------
        SEAT.MaterialComposition
            the instance created by he classmethod.

        Note
        ----
        The nuclides in `lines` should either be in the form 'Nn###', or ZAM.

        """
        components = {}
        for line in lines:
            line = [line_.strip() for line_ in line.split(separator)]
            nuclide = line[0] if not line[0].isnumeric() else zam2nuclide(int(line[0]))
            components[nuclide] = float(line[1])
        return cls(components, **kwargs)

    @classmethod
    def from_file(cls, file: str, separator: str=None, **kwargs):
        """
        Reads the material composition from a file.

        Parameters:
        ------
        file : str
            name of the file from which the material composition should be read.
        separator: str
            the separator of the nuclide symbol and the corresponding density.
            The default is None.
        **kwargs :
            - atomic : bool
                makes the material densities be interpreted as atomic or mass.

        Returns
        -------
        SEAT.MaterialComposition
            the instance created by he classmethod.

        Note
        ----
        The nuclides in the file should either be in the form 'Nn###', or ZAM.

        """
        with open(file, 'r') as f:
            lines = f.readlines()
        return cls._from_list(lines, separator, **kwargs)

    @classmethod
    def enriched_custom(cls, abundance: dict[str, float], enrichers: dict[str, float], **kwargs):
        """
        Creates the material composition from a custom abundance and an
        enrichment criterion.

        Parameters
        ----------
        abundance : dict[str, float]
            * key: the symbol of the nuclides (`'Nn###'`) or elements (`''Nn`)
                    in the reference material.
            * value: the abundance of each nuclide or element.
        enrichers : dict[str, float]
            * key: the symbol (`'Nn###'`) of the nuclides enriching the material.
            * value: the target isotopic abundance of the enriching nuclides.
        **kwargs :
            - atomic : bool
                makes the material densities be interpreted as atomic or mass.

        Raises
        ------
        Exception
            `abundance.keys()` should either contain elements or nuclides.
        Exception
            `enrichers.keys()` should contain nuclides.

        Returns
        -------
        SEAT.MaterialComposition
            the instance created by he classmethod.

        """
        if are_elements(set(abundance.keys())):
            abundance_ = unfold_composite(abundance)
        elif are_nuclides(set(abundance.keys())):
            abundance_ = abundance
        else:
            raise Exception("The keys of the abundance should be consistent in being elements or nuclides.")
        if are_nuclides(set(enrichers.keys())):
            enrichers_ = enrichers
        else:
            raise Exception("The keys of the enrichers should be nuclides.")
        return cls(enrich(abundance_, enrichers_), **kwargs)

    @classmethod
    def enriched(cls, base: str, **kwargs):
        """
        Enriches a natural material or a `SEAT.composites` item according to an
        enriching criterion.

        Parameters
        ----------
        base : str
            the element symbol or `SEAT.composites` symbol of the material to
            enrich.
        **kwargs :
            - atomic : bool
                makes the material densities be interpreted as atomic or mass.
            - enrichers : dict[str, float]
                * key: the symbol (`'Nn###'`) of the nuclides enriching the
                        material.
                * value: the target abundance of the enriching nuclides.

        Raises
        ------
        Exception
            `base` should either be in `SEAT.natural` or in `SEAT.composites`.

        Returns
        -------
        SEAT.MaterialComposition
            the instance created by he classmethod.

        """
        if base in dir(SEAT.natural):
            return cls.enriched_custom(getattr(SEAT.natural, base), atomic=True, **kwargs)
        elif base in dir(SEAT.composites):
            return cls.enriched_custom(getattr(SEAT.composites, base), atomic=True, **kwargs)
        else:
            raise Exception(f"{base} not found in the SEAT database.")

    @classmethod
    def polluted_custom(cls, abundance: dict[str, float], pollutants: dict[str, float], **kwargs):
        """
        Creates the material composition from a custom abundance and a pollution
        criterion.

        Parameters
        ----------
        abundance : dict[str, float]
            * key: the symbol of the nuclides (`'Nn###'`) or elements (`''Nn`)
                    in the reference material.
            * value: the abundance of each nuclide or element.
        pollutants : dict[str, float]
            * key: the symbol (`'Nn###'`) of the nuclides pollutant.
            * value: the target isotopic abundance of the polluting nuclides.
        **kwargs :
            - atomic : bool
                makes the material densities be interpreted as atomic or mass.

        Raises
        ------
        Exception
            `abundance.keys()` should either contain elements or nuclides.
        Exception
            `enrichers.keys()` should contain nuclides.

        Returns
        -------
        SEAT.MaterialComposition
            the instance created by he classmethod.

        """
        if are_elements(set(abundance.keys())):
            abundance_ = unfold_composite(abundance)
        elif are_nuclides(set(abundance.keys())):
            abundance_ = abundance
        else:
            raise Exception("The keys of the abundance should be consistent in being elements or nuclides.")
        if are_nuclides(set(pollutants.keys())):
            pollutants_ = pollutants
        elif are_elements(pollutants.keys()):
            pollutants_ = unfold_composite(pollutants)
        else:
            raise Exception("The keys of the pollutants should be nuclides.")
        return cls(pollute(abundance_, pollutants_), **kwargs)

    @classmethod
    def polluted(cls, base: str, **kwargs):
        """
        Pollutes a natural material or a `SEAT.composites` item according to a
        pollution criterion.

        Parameters
        ----------
        base : str
            the element symbol or `SEAT.composites` symbol of the material to
            pollute.
        **kwargs :
            - atomic : bool
                makes the material densities be interpreted as atomic or mass.
            - pollutants : dict[str, float]
                * key: the symbol (`'Nn###'`) of the nuclides pollutant.
                * value: the target isotopic abundance of the polluting nuclides.

        Raises
        ------
        Exception
            `base` should either be in `SEAT.natural` or in `SEAT.composites`.

        Returns
        -------
        SEAT.MaterialComposition
            the instance created by he classmethod.

        """
        if base in dir(SEAT.natural):
            return cls.polluted_custom(getattr(SEAT.natural, base), atomic=True, **kwargs)
        elif base in dir(SEAT.composites):
            return cls.polluted_custom(getattr(SEAT.composites, base), atomic=True, **kwargs)
        else:
            raise Exception(f"{base} not found in the SEAT database.")

    @classmethod
    def natural(cls, element: str, **kwargs):
        """
        Creates the material composition form natural abundances.

        Parameters
        ----------
        element : str
            the element symbol of the material in the material compostition.
        **kwargs :
            - atomic : bool
                makes the material densities be interpreted as atomic or mass.

        Returns
        -------
        SEAT.MaterialComposition
            the instance created by he classmethod.

        """
        return cls(getattr(SEAT.natural, element), atomic=True, **kwargs)

    @classmethod
    def predefined(cls, composite: str, **kwargs):
        """
        Creates the material composition form those predefined in `SEAT.composites`.

        Parameters
        ----------
        composite : str
            the composite symbol of the material in the material composition.
        **kwargs :
            - atomic : bool
                makes the material densities be interpreted as atomic or mass.

        Returns
        -------
        SEAT.MaterialComposition
            the instance created by he classmethod.

        """
        
        return cls.composite(getattr(SEAT.composites, composite), atomic=True, **kwargs)

    @classmethod
    def composite(cls, elements: dict[str, float], **kwargs):
        """
        Creates the material composition form elements.

        Parameters
        ----------
        elements : dict[str, float]
            * key: the element symbol.
            * value: the share of that element in the target composite.
        **kwargs :
            - atomic : bool
                makes the material densities be interpreted as atomic or mass.

        Returns
        -------
        SEAT.MaterialComposition
            the instance created by he classmethod.

        """
        return cls(unfold_composite(elements), **kwargs)

    @classmethod
    def from_zam(cls, zam: dict[int, float], **kwargs):
        """
        Creates material composition from a dictionary of ZAMs.

        Parameters
        ----------
        zam : dict[int, float]
            * key: the ZAM number of the nuclides.
            * value: the corresponding fraction (as positive floats) as values.
        **kwargs :
            - atomic : bool
                makes the material densities be interpreted as atomic or mass.

        Returns
        -------
        SEAT.MaterialComposition
            the instance created by he classmethod.

        """
        components = dict(zip([zam2nuclide(z) for z in zam.keys()], zam.values()))
        return cls(components, **kwargs)

    @classmethod
    def from_za(cls, za: dict[int, float], **kwargs):
        """
        Creates material composition from a dictionary of ZAs.

        Parameters
        ----------
        za : dict[int, float]
            * key: the ZA number of the nuclides.
            * value: the corresponding fraction (as positive floats) as values.
        **kwargs :
            - atomic : bool
                makes the material densities be interpreted as atomic or mass.

        Returns
        -------
        SEAT.MaterialComposition
            the instance created by he classmethod.

        """
        instance = cls(za, **kwargs)
        instance._already_za = True
        return instance

    @classmethod
    def from_sep_nuclides(cls, nuclides: dict[str, float], sep: str='', **kwargs):
        """
        Creates the material composition giving the `components.keys()` as (`Nn{sep}###`).

        Parameters
        ----------
        nuclides : dict[int, float]
            * key: the nuclide symbols (`'Nn{sep}###'`).
            * value: the corresponding fraction (as positive floats) as values.
        sep: str
            the nuclide separator.
        **kwargs :
            - atomic : bool
                makes the material densities be interpreted as atomic or mass.

        Returns
        -------
        SEAT.MaterialComposition
            the instance created by he classmethod.

        """
        components = dict(zip([n.replace(sep, '') for n in nuclides.keys()], nuclides.values()))
        return cls(components, **kwargs)

    @classmethod
    def from_output_file(cls, file: str):  # interaction with serpentTools missing
        """
        TO BE IMPLEMENTED
        Creates material composition from a Serpent 2 output file

        Takes:
        ------
        * `file`: string - name of the file from which the material components should be read
        """
        pass

    def write(self, file: str, mode: str='w'):
        """
        Writes the `SEAT.MaterialComposition` instance to a file.

        Parameters
        ----------
        file : str
            name of the file to write to.
        mode : str, optional
            mode to open the file. The default is 'w'.

        Returns
        -------
        None.

        """
        with open(file, mode=mode) as f:
            f.write(self.__str__())

    def adjust(self, library: str, verbose: bool=True):
        """
        Adjusts the components of the MaterialComposition redistributing the
        nuclides which cross section data are not evaluated in a library to the
        other isotopes of the same element, according to their natural abundance.

        Parameters
        ----------
        library : str
            the library of which the available cross section set should be got.
        verbose: bool, optional
            prints the excluded nuclides. The default is False.

        Returns
        -------
        None

        """
        available = get_existing_xs(library)
        exclude = {k: 0 for k in self.components.keys() - available}
        if verbose: print(f"Excluding:\n{exclude}")
        if exclude: self.components = {k: v for k, v in enrich(self.components, exclude).items() if v != 0}
        return self

    def copy(self):
        """
        Copies the object instance to another memory allocation.

        Returns:
        --------
        SEAT.Other
            a copy of the `SEAT.MaterialComposition` instance.

        """
        return cp.deepcopy(self)


@dataclass(slots=True)
class MaterialRepresentation(MaterialComposition):
    """
    Handles the material representation coupling its components to the library
    identifier and the temperature at which nuclear data should be generated.
    Inherits from `SEAT.MaterialComposition`.

    Args:
    ------
    tmp : float
        the temperature at which nuclear data should be generated [K].
    data_type : str, optional
        the library identifier. The default is 'c' for cross sections.

    Methods:
    --------
    get_temperature : 
        formats the temperature according to Serpent 2 syntax.
    get_nuclide_representation :
        formats the nuclides in the composition as: {ZA}.{`tmp`}{`data_type`}.
    get_values :
        gets the density of the nuclides in the `SEAT.MaterialRepresentation`.
    write:
        writes the `SEAT.MaterialRepresentation` instance to a file.

    Class methods:
    --------------
    parse :
        reads the material composition from a string.
    _from_list :
        reads the material composition from a list of strings.
    from_file :
        reads the material composition from a file.
    from_output_file :
        readsthe material composition from a Serpent2 output file - not implemented yet.
    enriched_custom :
        enriches a custom material composition.
    enriched :
        enriches a predefined material.
    polluted_custom :
        creates the material composition from a pollution criterion.
    polluted :
        creates the material composition polluting a natural element or one form `SEAT.composites`.
    natural :
        creates the material composition form natural abundances.
    predefined :
        creates the material composition form those predefined in `SEAT.composites`.
    composite :
        creates the material composition form elements.
    from_zam :
        creates the material composition giving the `components.keys()` as ZAMs.
    from_za :
        creates the material composition giving the `components.keys()` as ZAs.
    from_sep_nuclides :
        creates the material composition giving the `components.keys()` as (`Nn{sep}###`).

    """
    tmp: float
    data_type: str = 'c'

    def __str__(self):
        string = ''
        for k, v in zip(self.get_nuclide_representation(), self.get_values()):
            string += f"{k} {v}\n"
        return string

    def get_temperature(self) -> str:
        """
        Gets the temperature of the `SEAT.MaterialRepresentation` formatted as
        in the Serpent 2 input file.

        Returns
        -------
        str
            the formatted temperature.

        """
        
        if self.tmp < 1000:
            tmp_ = '0' + '{:.0f}'.format(self.tmp / 100)
        else:
            tmp_ = '{:.0f}'.format(self.tmp / 100)
        return tmp_

    def get_nuclide_representation(self) -> list[str]:
        """
        Gets the representation of each nuclide in the composition as:
            {ZA}.{temperature}{data_type}.

        Returns
        -------
        list[str]
            with formatted nuclides as elements.

        """
        temperature, data_type = self.get_temperature(), self.data_type
        out = [f"{k}.{temperature}{data_type}" for k in self.components.keys()
               ] if self._already_za else [
                   f"{nuclide2za(k)[0]}.{temperature}{data_type}"
                   for k in self.components.keys()]
        return out

    def get_values(self) -> list[float]:
        """
        Gets the density of the nuclides in the `SEAT.MaterialRepresentation`.

        Returns
        -------
        list[float]
            with the densities as items.

        """
        return [v for v in self.components.values()]

    def write(self, file, mode: str='a'):
        """
        Writes the `SEAT.MaterialRepresentation` instance to a file.

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
        with open(file, mode=mode) as f:
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
    def mix(cls, name: str, materials: list[tuple], **kwargs) -> object:
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
        return cls(name=name, _mixture=''.join(a), **kwargs)

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
            warnstring = f"No temperature given to {self.name}.\n" if temperature is None else f"Multiple temperatures given to {self.name}.\n"
            warnstring += "Temperature will be written as if it was the same as the material representation temperature"
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
