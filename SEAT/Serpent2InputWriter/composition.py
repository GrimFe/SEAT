import copy as cp
import warnings
from dataclasses import dataclass, field

import numpy as np

import SEAT.natural
import SEAT.composites
from SEAT.composition_functions import *

from SEAT.nuclides import zam2nuclide, nuclide2zam, nuclide2za
from SEAT.Serpent2InputWriter.base import Entity, Comment, reformat

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

    Attributes
    ----------
    components : dict[str, float]
        * key: nuclide symbol (`'Nn###'`)
        * value: densities as positive floats.
    atomic : bool
        makes the nuclide densities be interpreted as atomic or mass.
    _already_za : bool
        set to True when `components.keys()` was made of ZAs.

    Methods
    --------
    write :
        writes the `SEAT.MaterialComposition` instance to a file.
    adjust :
        adjusts the material composition to a nuclear data library.

    Class methods
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
        creates the material composition polluting a natural element or one from `SEAT.composites`.
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
            the instance created by the classmethod.

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
            the instance created by the classmethod.

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

        Parameters
        ----------
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
            the instance created by the classmethod.

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
            the instance created by the classmethod.

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
            the instance created by the classmethod.

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
            the instance created by the classmethod.

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
            the instance created by the classmethod.

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
            the instance created by the classmethod.

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
            the instance created by the classmethod.

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
            the instance created by the classmethod.

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
            the instance created by the classmethod.

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
            the instance created by the classmethod.

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
            the instance created by the classmethod.

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
            prints the excluded nuclides. The default is True.

        Returns
        -------
        None.

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

    Attributes
    ----------
    tmp : float
        the temperature at which nuclear data should be generated [K].
    data_type : str, optional
        the library identifier. The default is 'c' for cross sections.

    Methods:
    --------
    write :
        writes the `SEAT.MaterialRepresentation` instance to a file.
    adjust :
        adjusts the material composition to a nuclear data library.
    get_temperature : 
        formats the temperature according to Serpent 2 syntax.
    get_nuclide_representation :
        formats the nuclides in the composition as: {ZA}.{`tmp`}{`data_type`}.
    get_values :
        gets the density of the nuclides in the `SEAT.MaterialRepresentation`.

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
    Handles the material definition and the operations on single materials.
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
    dens : float, optional
        the material density. The default is None.
    nuclide_density : bool, optional
        makes the material density be interpreted as nuclide [cm-3] or mass
        [g/cm3]. The default is True.
    vol : float, optional
        volume or cross-section (for 2D calculations) [cm3 (cm2)]. The default
        is None.
    mass : float, optional
        material mass [g]. The default is None.
    representation : `SEAT.MaterialRepresentation`
        the material composition representation. The default is None.
    tmp : float, optional
        the temperature [K] for doppler preprocessing. The efault is None.
    tms : float, optional
        the temperature [K] for on-the-fly temperature treatment. The default
        is None.
    tft : tuple[float], optional
        contains two floats for the material temperature limits [K]. The
        default is None.
    rgb : tuple[int], optional
        contains three integers in range 0-255 for color channels in the plot:
        red, green, blue respectively. The default is None.
    burn : int, optional
        the number of regions to burn the material. The default is 0.
        Better use of material divisions is handled by the `divide()` method.
    fix : dict[str, float], optional
        the library information for decay nuclides (i.e. without cross-section data).
        * key: library id.
        * valu: temperature [K].
        Only one key-value pair is allowed. The default is None.
    moder : dict[str, int], optional
        moderators and corresponding thermal scattering library (tsl).
        * key: name of the thermal scattering data library defined using the
                therm card.
        * value: ZA number of the thermal scatter (1001 for H1). The default is
                None.
    _mixture : str, optional
        the mixed materials and their fractions. The efault is None.
    _divisions: list[`SEAT.Division`], optional
        the material divisins. The default is None.
    _division_string: str, optional
        header string to the divisions. The default is ''.

    Methods:
    --------
    assess :
        prints the `SEAT.Material` python id.
    write :
        writes the `SEAT.Material` to a file.
    get_temperature :
        gets the temperature of the material.
    get_temperature_kind :
        gets the kind of temperature in the material definition.
    divide :
        computes and creates the desired material divisions.

    Class Methods:
    --------------
    mix :
        mixes materials in a composed one.

    """
    dens: float | None = None
    nuclide_density: bool = True
    vol: float = None
    mass: float = None
    representation: MaterialRepresentation = None
    tmp: float = None
    tms: float = None
    tft: tuple[float, float] = None
    rgb: tuple[int, int, int] = None
    burn: int = 0
    fix: dict[str, float] = None
    moder: dict[str, int] = None
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
            string += f" fix {self.fix.keys()[0]} {self.fix.values()[0]}"
        if self.moder is not None:
            string += reformat(f" moder {self.moder}", ",'()[]{}")
        string += self.inline_comment.__str__()
        if self._mixture is None:
            string += self.representation.__str__()
        string += "\n"
        return string

    @classmethod
    def mix(cls, name: str, materials: list[tuple], **kwargs):
        """
        Creates material mixing other materials.

        Parameters
        ----------
        name : str | int
            the identity of the Serpent 2 entity.
        materials : list[tuples[`SEAT.Material`, flaot]]
            each tuple contains:
                * a material.
                * the corresponging fraction.
        **kwargs :
            - .. missing

        Returns
        -------
        SEAT.Material
            the instance created by the classmethod.

        """
        a = []
        for k, v in materials:
            a.append(f"{k.name} {v}\n")
        return cls(name=name, _mixture=''.join(a), **kwargs)

    def get_temperature(self) -> str:
        """
        Gets the temperature of the material if any of `tms`, `tmp` or `tft`
        is defined, else gets the temperature of `SEAT.Material.representation`.

        Returns
        -------
        str
            the material temperature.
        
        Note
        ----
        If more no temperature or more than one temperatures are given to the
        material, this returns the representation temperature. In doing so it
        raises some warnings.

        """
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

    def get_temperature_kind(self) -> str:
        """
        Gets the kind of temperature in the material definition.

        Raises
        ------
        Exception
            If no temperature type or multiple temperature types were given.

        Returns
        -------
        kind : str
            the temperature type.

        """
        kind = None
        number = 0
        if self.tms is not None:
            kind = 'tms'
            number += 1
        if self.tmp is not None:
            kind = 'tmp'
            number += 1
        if self.tft is not None:
            kind = 'tft'
            number += 1
        string = "More than one temperature type" if number > 1 else "No temperature type"
        if number != 1:
            raise Exception(string + f' was given to Material {self.name}')
        return kind

    def divide(self, lvl: int, kind: str, sub: tuple[int | float]):
        """
        Method to create material divisions and compute sub-volumes.

        Parameters
        ----------
        lvl : int
            indicates the cell level of the division.
            (0: no division, 1: last level division, 2: second last level division, ...)
        kind : str
            indicates the type of division to implement.
            Allowed `kind` values are:
                - 'x_cartesian'
                - 'y_cartesian'
                - 'z_cartesian'
                - 'radial'
                - 'sectorial'
        sub : tuple[int | float]
            data for the division (lengths in [cm]).
            It can be composed as follows:
                - for cartesian and radial division:
                    `sub` = (N, min, max), where N is the number of zones and
                            min and max are respectively the minimum an maximum
                            coordinates in the direction defined by `kind`.
                - for radial division it can also be:
                    `sub` = (N, r0, r1, .., rn), where N is the number of zones
                            and ri is the radial coordinate of the i-th division.
                - for sectorial division:
                    `sub` = (N, s0), where N is the number of zones and s0 the
                            angular coordinate of the first in [deg].

        Notes:
        ------
        Division volume calculation not implemented for cartesian and sectorial
        division yet.
        Only 2D subdivision is implemented.

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
    Handles the cell division definition and operations on single divisions.

    Attributes:
    -----------
    material : `SEAT.Material`
        the material of the division.
    volume : float
        the volume of the division.

    Methods:
    --------
    write :
        writes the `SEAT.Division` to a file.
    copy :
        copies the object instance to another memory allocation.

    """
    material: Material
    volume: float

    def __str__(self):
        string = f"{self.material.name} {self.volume}"
        return string

    def write(self, file: str, mode: str='a'):
        """
        Writes the `SEAT.Division` instance to a file.

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

    def copy(self):
        """
        Copies the object instance to another memory allocation.

        Returns
        -------
        SEAT.Other
            a copy of the `SEAT.Other` instance.

        """
        return cp.deepcopy(self)


@dataclass(slots=True)
class Composition:
    """
    Handles the composition of the Serpent 2 input, which composes of:
        * nuclear data libraries link
        * unresolved resonance probability table
        * thermal scattering data
        * restart files
        * materials

    Attributes
    ----------
    materials : list[`SEAT.Material`]
        the materials of the simulation.
    libraries : dict[str, str]
        nuclear data libraries to use in the simulation.
        * keys: Serpent library identifiers (.. list options ..).
        * values: raw-string paths as values.
    scattering_name : str, optional
        thermal scattering data name. The default is None.
    scattering_type : str, optional
        interpolation method.
        Allowed `scattering_type` values are:
            - '': interpolation with the makxsf code methodology.
            - 'stoch': stochastic interpolation. Not allowed for on-the-fly
                    interpolation.
        The default is ''.
    scattering_tmp : float, optional
        interpolation temperature for the scattering data. The default is None.
    scattering_libs : list[str], optional
        items are the thermal scattering data identifiers. The default is [].
    to_restart : str, optional
        the name of the binary restart file the composition should be written
        to. The default is None.
    from_restart : str
        the name of the binary restart file the composition should be read
        from. The default is None.

    Methods
    -------
    write :
        writes the `SEAT.Composition` to a file.
    get_subdivisions :
        lists the subdivisions of the materials in the composition.

    
    """
    materials: list[Material]
    libraries: dict[str, str]
    scattering_name : str = None
    scattering_type : str = ''
    scattering_tmp : float = None
    scattering_libs : list[str] = []
    to_restart: str = None
    from_restart: str = None

    def __str__(self):
        string = ''
        for k, v in self.libraries.items():
            string += f"{k} '{v}'\n"
        string += '\n'
        if self.scattering is not None:
            libs = reformat(str(self.scattering_libs), "[],'")
            string += f'therm{self.scattering_type} {self.scattering_name} {self.scattering_tmp} {libs}\n'
        if self.to_restart is not None or self.from_restart is not None:
            string += '\n'
            string += Comment("Composition restart file definition").__str__()
            string += f"rfw {self.to_restart}\n" if self.to_restart is not None else ''
            string += f"rfw {self.from_restart}\n" if self.from_restart is not None else ''
        string += '\n'
        for m in self.materials:
            string += m.__str__()
        string += '\n'
        return string

    def write(self, file: str, mode: str='a'):
        """
        Writes the `SEAT.Composition.__str__()` to a file.

        Args:
        -----
        file : str
            the name of the file where to write.
        mode : str, optional
            mode to open the file. The default is 'a'.

        Returns
        -------
        None.

        """
        with open(file, mode=mode) as f:
            f.write(self.__str__())

    def get_subdivisions(self) -> list[Division]:
        """
        Lists of all the subdivisions in the `SEAT.Composition`.

        Returns:
        --------
        list[`SEAT.Division`]
            the subdivisions.

        """
        out = []
        for m in self.materials:
            out.extend(m._divisions)
        return out
