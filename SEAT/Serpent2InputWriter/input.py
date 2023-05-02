import copy as cp
import numpy as np
import warnings
import dataclasses
from dataclasses import dataclass, field

from SEAT.Serpent2InputWriter.composition import Composition, MaterialRepresentation, Division
from SEAT.Serpent2InputWriter.geometry import Geometry, Universe
from SEAT.Serpent2InputWriter.depletion import Depletion
from SEAT.Serpent2InputWriter.base import Comment, Other, reformat

from SEAT.nuclides import nuclide2zam

import SEAT.Serpent2InputWriter._DefaultOptions as DefaultOptions
import SEAT.Serpent2InputWriter._AllowedFields as AllowedFields

__author__ = "Federico Grimaldi"
__all__ = [
    "Simulation",
]


API_LINK = 'https://github.com/GrimFe/SEAT'
HEADER_COMMENT = f"""/*
#############################################################################################################################################
#                                                                                                                                           #
#                                            Serpent2 input file written with SEAT 0.0.1 python API                                         #
#                                                           https://github.com/GrimFe/SEAT                                                  #
#                                                                                                                                           #
#############################################################################################################################################
*/\n\n"""

mvol_counter = {}

def check_allowed_key(to_check, target) -> bool:
    """
    Checks whether a dictionary is a sub-dictionary of another.

    Parameters
    ----------
    to_check: dict
        the dictionary to check
    target: dict
        the target dictionary

    Returns
    -------
    test: bool
        True if the `to_check` is a sub-dictinary of `target`.

    """
    test = np.prod([k in target.keys() for k in to_check.keys()])
    return test


@dataclass(slots=True)
class DivisionWrapper:
    """
    Handles the labelling of the divisions.

    Parameters
    ----------
    div: list[`SEAT.Division`]
        the divisions to label and sort.

    """
    divisions: list[Division]

    def __iter__(self):
        return self.divisions.__iter__()

    def format_sort(self) -> str:
        """
        Labels and sorts the divisions in the wrapper writing them to a string.

        Returns
        -------
        str
            the formatted string of the divisions.

        """
        counter = {}
        ordered = [d.__str__() for d in self.divisions]
        ordered.sort()
        formatted = []
        for item in ordered:
            item_ = item.split()
            if item_[0] not in counter.keys():  # counter keys are division material names
                counter[item_[0]] = 1
            formatted.insert(0, f"{item_[0]} {counter[item_[0]]} {item_[1]}")
            counter[item_[0]] += 1
        return reformat(str(formatted), "[]'").replace(', ', '\n')


@dataclass(slots=True)
class DivisionWriter:
    """
    Handles the writing process of the divisions.

    Attributes
    ----------
    strings : list[str]
        the division strings.
    divisions : `SEAT.DivisionWrapper`
        properly wrapped divisions to write.

    """
    strings: list[str]
    divisions: DivisionWrapper

    def to_string(self) -> str:
        """
        Writes the wrapped divisions to a string properly headed.

        Returns
        -------
        str
            the division string to write to the simulation.

        """
        string = ''
        if not self.empty:
            string = reformat(str(set(self.strings)), "{}'").replace(', ', '\n')  # is the cast to set really needed? 'set' should also be removed form the string
            string += '\n'
            if None not in self.divisions:
                string += '\n'
                string += "set mvol\n\n"
                string += self.divisions.format_sort()
        return string

    @property
    def empty(self) -> bool:
        """
        Assess whether there are elements in the DivisionWriter.

        Returns
        -------
        bool
            True if the DivisionWriter is empty.

        """
        return not bool(self.strings)


@dataclass(slots=True)
class Plot:
    """
    Handles the various plot types.

    Attributes
    ----------
    xpix : int
        horizontal image size in pixels.
    ypix : int
        vertical image size in pixels.

    """
    xpix: int
    ypix: int


@dataclass(slots=True)
class GeometryPlot(Plot):
    """
    Geometry plot defition.

    Attributes
    ----------
    xpix : int
        horizontal image size in pixels.
    ypix : int
        vertical image size in pixels.
    kind : str
        the plot type can be `'xy'`, `'xz'` or `'yz'`.
    pos : float, optional
        location of the plot plane on the axis perpendicular to it. The default
        is None.
    min1 : float, optional
        lower boundary of the first plot direction. The default is None.
    min2 : float, optional
        lower boundary of the second plot direction. The default is None.
    max1 : float, optional
        upper boundary of the first plot direction. The default is None.
    max2 : float, optional
        upper boundary of the second plot direction. The default is None.

    """
    kind : str
    pos : float = None
    min1 : float = None
    min2 : float = None
    max1 : float = None
    max2 : float = None

    def __str__(self):
        string = str(AllowedFields.PLOT_KIND[self.kind]) + ' '
        if self._importance is not None:
            string += self._importance + ' '
        string += str(self.xpix) + ' '
        string += str(self.ypix) + ' '
        string += str(self.pos) + ' '
        string += str(self.min1) + ' '
        string += str(self.max1) + ' '
        string += str(self.min2) + ' '
        string += str(self.max2) + '\n'
        return string

    @property
    def _importance(self) -> str:
        """
        The formatted importance plot options.
        
        Returns
        -------
        None
            as the base geometry plot is not an importance map plot.
            Use SEAT.ImportancePlot instead.

        """
        return None

@dataclass(slots=True)
class ImportancePlot(GeometryPlot):
    """
    Importance map plot defition.

    Attributes
    ----------
    xpix : int
        horizontal image size in pixels.
    ypix : int
        vertical image size in pixels.
    kind : str
        the plot type can be `'xy'`, `'xz'` or `'yz'`.
    pos : float, optional
        location of the plot plane on the axis perpendicular to it. The default
        is None.
    min1 : float, optional
        lower boundary of the first plot direction. The default is None.
    min2 : float, optional
        lower boundary of the second plot direction. The default is None.
    max1 : float, optional
        upper boundary of the first plot direction. The default is None.
    max2 : float, optional
        upper boundary of the second plot direction. The default is None.
    fmin : float, optional
        minimum importance for importance map plots. The default is None.
    fmax : float, optional
        maximum importance for importance map plots. The default is None.
    energy : float.
        particle energy for importance map plots. The default is None.

    """
    fmin : float =  None
    fmax : float = None
    energy : float = None

    @property
    def _importance(self) -> str:
        """
        The formatted importance plot options.
        
        Returns
        -------
        str
            the formatted importance plot string.

        """
        warn = ''
        string = ''
        if self.fmin is not None:
            string += str(self.fmin) + ' '
        else:
            warn += 'fmin '
        if self.fmin is not None:
            string += str(self.fmax) + ' '
        else:
            warn += 'fmax '
        if self.fmin is not None:
            string += str(self.energy) + ' '
        else:
            warn += 'energy'
        if warn != '':
            w = f"The following required parameters are set to None: {warn}"
            warnings.warn(w)
        return string


@dataclass(slots=True)
class MeshPlot(Plot):
    pass


@dataclass(slots=True)
class Simulation:
    """
    Defines a static simulation.
    Hanldes the model as a whole allowing for the definition of simulation
    specifications. It is composed of sub-parts corresponding to the model
    sections.

    Attributes
    ----------
    geometry : `SEAT.Geometry`
        the geometry of the simulation.
    composition : `SEAT.Composition`
        the materials and nuclear data used in the simulation.
    others : list[`SEAT.Other`], optional
        for what is not implemented yet. The default is None
    comments : dict[str, `SEAT.Comment`], optional
        comments to the simulation sections.
        * keys: comment section identifier.
        * values: comment for the section.
        Allowed dictionary keys are:
            - 'Intro': introductory comments
            - 'Geometry': geometry comments
            - 'Materials': composition comments
            - 'Others': other comments
            - 'Concluding': concluding comments
            The default is None (all comments are empty).
    seed : int, optional
        seed for the random number generation. The default is None.
    dix : bool
        double indexing method for the cross section energy grid look-up. The
        default is False.
    population : dict[str, float]
        the particle population settings in the simulation.
        Allowed dictionary keys are:
            - 'particles': the number of neutrons per generation.
                            The default is 100000.
            - 'generations': the number of active generations.
                            The default is 250.
            - 'inactive': the number of inactive generations.
                            The default is 30.
            - 'k_guess': the guess value for k effective.
                            The default is 1.
            - 'batch_interval': the batching interval.
                            The default is 1.
            - 'parallel': the number of independent parallel eigenvalue
                            calculations. The default is 1.
    gcu : list[`SEAT.Universe`]
        the universes of the group constants. The default is None (no
        generation).
    ures : dict[str, any], optional
        Unresolved REsonance Probability table Sampling (ures) settings.
        Allowed dictionary keys are:
        - 'active': boolean to activate ures. The default is False.
        - 'representation': `SEAT.MaterialRepresentation` of the nuclides ures
                            should be applyed to. The default is None.
        - 'dilution_cut': float for the infinite dilution cut-off. The default
                        is 1e-9.
    bc : int | str | tuple[int|str], optional
        the boundary conditions to the simulation.
            The tuple is used to set separate bc for the X, Y, Z directions.
            Allowed `bc` values are:
                - 'vacumm'        or 'V'   i.e. 1: vacuum boundary conditions
                - 'reflective'    or 'R'   i.e. 2: reflective boundary conditions
                - 'periodic'      or 'P'   i.e. 3: periodic boundary conditions
            The default is 1, corresponding to 'vacuum'.
    albedo : float
        the value of the albedo to set for reflective and/or periodic boundary
        conditions. The default is 1.
    plots : list[`SAET.Plot`], optional
        the plots to include in the simulation. The default is None.
    _file : str, optional
        name of the file which the simulation is written to. The default is
        None.
    _written : bool, optional
        indicates whether the simulation has already been written. The
        default is False.
    _restart_filename : str, optional
        name of the binary file the simulation composition should be written
        to. The default is None.

    Methods:
    --------
    clear :
        clears what is written on the input file.
    write :
        writes the model to the input file. This process goes section by section.

    """
    geometry: Geometry
    composition: Composition
    others: list[Other] = None
    comments: dict[str, Comment] = field(
        default_factory=lambda: DefaultOptions.COMMENT_BASE_DCT)
    seed: int = None
    dix: bool = False
    population: dict[str, int] = field(
        default_factory=lambda: DefaultOptions.POPULATION)
    gcu: list[Universe] = None
    ures: dict[str, any] = field(default_factory=lambda: DefaultOptions.URES)
    bc: tuple[str | int] | str | int = 1
    albedo: float = 1
    plots: list[Plot] = None
    _file: str = None
    _written: bool = False

    def __str__(self):
        string = HEADER_COMMENT
        string += self.comments['Intro'].__str__()
        string += Comment("Global simulation settings").__str__()
        string += self._seed
        string += self._population
        string += self._bc
        string += '\n'
        string += Comment("Unresolved resonance probability table settings").__str__()
        string += self._ures
        string += Comment("Double indexing method settings").__str__()
        string += self._dix
        string += Comment("Group constant generation settings").__str__()
        string += self._gcu
        if self.plots is not None:
            for p in self.plots:
                string += p.__str__()
        string += '\n'
        string += self.comments['Geometry'].__str__()
        string += self.geometry.__str__()
        string += self.comments['Materials'].__str__()
        string += self.composition.__str__()
        if self.others is not None:
            string += self.comments['Others'].__str__()
            string += '\n'
            for other in self.others:
                string += other.__str__()
        string += '\n'
        if self._depletion is not None:
            string += self.comments['Depletion'].__str__()
            string += self._depletion
        string += self.comments['Concluding'].__str__()
        return string

    @property
    def _seed(self) -> str:
        """
        Serpent2-formatted seed. set seed card.

        Returns
        -------
        str
            the formatted seed.

        """
        return f"set seed {self.seed}\n" if self.seed is not None else ''

    @property
    def _dix(self) -> str:
        """
        Serpent2-formatted double indexing method. set dix card.

        Returns
        -------
        str
            the formatted dix.

        """
        return f'set dix {int(self.dix)}\n'

    @property
    def _gcu(self) -> str:
        """
        The formatted universe group constant generator. set gcu card.

        Returns
        -------
        str
            the formatted gcu.

        """
        if self.gcu is None:
            out = '-1'
        else:
            out = reformat(str([i.name for i in self.gcu]), "[],'")
        return 'set gcu ' + out + '\n\n'

    @property
    def _population(self) -> str:
        """
        Serpent2-formatted particle population. set pop card.

        Returns
        -------
        str
            the formatted particle population.

        """
        if not check_allowed_key(self.population, DefaultOptions.POPULATION):
            warnings.warn("Some of the population keys are not alowed.")
        pop = DefaultOptions.POPULATION | self.population
        pop_str = f"{pop['particles']} {pop['generations']} {pop['inactive']} "\
                  f"{pop['k_guess']} {pop['batch_interval']} {pop['parallel']}"
        return 'set pop ' + pop_str + "\n" 

    @property
    def _ures(self) -> str:
        """
        Serpent2-formatted unresolved resonance table sampling. set ures card.

        Returns
        -------
        str
            the formatted ures.

        """
        if not check_allowed_key(self.ures, DefaultOptions.URES):
            warnings.warn("Some of the ures keys are not alowed.")
        ures = DefaultOptions.URES | self.ures
        ures_str = 'set ures '
        ures_str += f'{int(ures["active"])}'
        if ures["active"]:
            ures_str += ' '
            ures_str += reformat(str(ures["representation"].get_nuclide_representation()), "[],'")
            ures_str += f' {ures["dilution_cut"]}'
        return ures_str + '\n\n'

    @property
    def _bc(self) -> str:
        """
        Serpent2-formatted model boundary conditions. set bc card.

        Returns
        -------
        str
            the formatted bc.

        """
        if isinstance(self.bc, tuple):
            mode_ = [AllowedFields.BC[bc] for bc in self.bc]
            mode = reformat(str(mode_), "[],")
            _albedo = f' {self.albedo}' if (AllowedFields.BC['R'] in mode_ or
                                            AllowedFields.BC['P'] in mode_) else ''
        else:
            mode = AllowedFields.BC[self.bc]
            _albedo = f' {self.albedo}' if (mode == AllowedFields.BC['R'] or
                                            mode == AllowedFields.BC['P']) else ''
        return f'set bc {mode}{_albedo}\n'

    @property
    def _depletion(self) -> str:
        """
        The formatted depletion options for the simulation.
        
        Returns
        -------
        None
            as the base simulation does not foresee any depletion.
            Use SEAT.DepletionSimulation instead.

        """
        return None

    def clear(self):
        """
        Clears what is written on the input file.

        Returns
        -------
        None.

        Note
        ----
            If the output file name is not selected raises a warning.

        """
        if self._file is not None:
            with open(self._file, 'w'):
                pass
        else:
            warnings.warn("No output file is set for the simulation object")

    def write(self, file: None):
        """
        Writes the Serpent 2 input file progressively appending section by section in the order:
            * Introductory comments
            * Comments on the geometry
            * Geometry
            * Comments on the depletion
            * Depletion
            * Comments on the materials
            * Materials
            * Comments on other elements
            * Other elements
            * Cells
            * Concluding comments

        Parameters
        ----------
        file : str, optional
            the name of the file to write the simulation to. The default is
            None.

        Returns
        -------
            None.

        """
        if self._file is not None and file is not None:
            warnings.warn("Ignoring a previously defined file was already defined for this simulation.")
        if file is not None:
            self._file = file            
        with open(self._file, 'w') as f:
            f.write(self.__str__())
        self._written = True


@dataclass(slots=True)
class DepletionSimulation(Simulation):
    """
    Defines a depletion simulation.

    Attributes
    ----------
    geometry : `SEAT.Geometry`
        the geometry of the simulation.
    composition : `SEAT.Composition`
        the materials and nuclear data used in the simulation.
    others : list[`SEAT.Other`], optional
        for what is not implemented yet. The default is None
    comments : dict[str, `SEAT.Comment`], optional
        comments to the simulation sections.
        * keys: comment section identifier.
        * values: comment for the section.
        Allowed dictionary keys are:
            - 'Intro': introductory comments
            - 'Geometry': geometry comments
            - 'Materials': composition comments
            - 'Others': other comments
            - 'Concluding': concluding comments
            The default is None (all comments are empty).
    seed : int, optional
        seed for the random number generation. The default is None.
    dix : bool
        double indexing method for the cross section energy grid look-up. The
        default is False.
    population : dict[str, float]
        the particle population settings in the simulation.
        Allowed dictionary keys are:
            - 'particles': the number of neutrons per generation.
                            The default is 100000.
            - 'generations': the number of active generations.
                            The default is 250.
            - 'inactive': the number of inactive generations.
                            The default is 30.
            - 'k_guess': the guess value for k effective.
                            The default is 1.
            - 'batch_interval': the batching interval.
                            The default is 1.
            - 'parallel': the number of independent parallel eigenvalue
                            calculations. The default is 1.
    gcu : list[`SEAT.Universe`]
        the universes of the group constants. The default is None (no
        generation).
    ures : dict[str, any], optional
        Unresolved REsonance Probability table Sampling (ures) settings.
        Allowed dictionary keys are:
        - 'active': boolean to activate ures. The default is False.
        - 'representation': `SEAT.MaterialRepresentation` of the nuclides ures
                            should be applyed to. The default is None.
        - 'dilution_cut': float for the infinite dilution cut-off. The default
                        is 1e-9.
    bc : int | str | tuple[int|str], optional
        the boundary conditions to the simulation.
            The tuple is used to set separate bc for the X, Y, Z directions.
            Allowed `bc` values are:
                - 'vacumm'        or 'V'   i.e. 1: vacuum boundary conditions
                - 'reflective'    or 'R'   i.e. 2: reflective boundary conditions
                - 'periodic'      or 'P'   i.e. 3: periodic boundary conditions
            The default is 1, corresponding to 'vacuum'.
    albedo : float
        the value of the albedo to set for reflective and/or periodic boundary
        conditions. The default is 1.
    _file : str, optional
        name of the file which the simulation is written to. The default is
        None.
    _written : bool, optional
        indicates whether the simulation has already been written. The
        default is False.
    _restart_filename : str, optional
        name of the binary file the simulation composition should be written
        to. The default is None.
    depletion : `SEAT.Depletion`, optional
        the neutron flux normalization and the time definition. The default is
        None.
    fpcut : float, optional
        the cumulative fission product yield cutoff in each mass chain. The
        default is 0.
    xscalc : int | str, optional
        control on how transmutation cross-sections are computed.
        Allowed `xscalc` vlaues are:
            - None: no calculation
            - 'NO': no calculation
            - 'DIRECT': direct tallies
            - 'SPECTRUM': spectrum collapse
            - 1: direct tallies
            - 2: spectrum collapse
        The default is None.
    printm_fraction : tuple[bool, float], optional
        control on the material composition printed to '.bumat' files.
        It composes of an on/off switch and of the minimum atomic fraction to
        print: 1 prints any decay nuclide, 0 prints those that don't have a
        transport cross section. The default is (False, 1).
        The use of resetart files is recommeded instead.
    restart_filename: str
        the name of the restart file for the material composition.
        The default is None (restart file name = restart+self._file in the
                             `write` method).
    bumode : str, optional
        burnup calculation mode.
        Allowed `bumode` options are:
            - 'TTA': Transmutation Trajectory Analysis
            - 'CRAM': Chebyshev Rational Approximation Method (parameters in
                    `cram_bumode`)
            - 1: Transmutation Trajectory Analysis
            - 2: Chebyshev Rational Approximation Method (parameters in
                    `cram_bumode`)
        The default is 'CRAM'.
    cram_param: tuple[int], optional
        CRAM burnup mode settings. The first item is the order, the second
        is the number of substeps.
        Allowed `cram_param` order values are:
            - 2
            - 4
            - 6
            - 8
            - 10
            - 12
            - 14
            - 16
            - -16
            - -48
         The default is None.
    pcc : int | str, optional
        predictor corrector scheme options.
        Allowed `pcc` valeus are:
                                  Predictor method            Corrector method
            - 'CE'   i.e. 0: constant extrapolation
            - 'CELI' i.e. 1: constant extrapolation       linear interpolation
            - 'LE'   i.e. 2: linear extrapolation
            - 'LELI' i.e. 3: linear extrapolation         linear interpolation
            - 'LEQI' i.e. 4: linear extrapolation         quadratic interpolation
            - 'CECE' i.e. 6: constant extrapolation       constant backwards extrapolation
            - None           : no predictor corrector scheme
        The default is 0, corresponding to 'CE'.
    pcc_param : tuple[int], optional
        the number of steps for the pcc. The first value is for the predictor,
        the second is for the corrector. The default is None.
    inventory : list | string | tuple, optional
        nuclides to report in the simulation output.
        It can be:
            - list of materials to include in 
                Isotope form `['U-235', 'Am-242']` or `['U235', 'Am242']`,
                ZAI form `[922350, 952421]`,
                Element in letters or mass number `['U', 'Am']` or `[92, 95]`
                All the items in the list should follow the same notation.
            - tuple of number of nuclides to include and parameter according to
              which choose them.
                Allowed parameters are:
                    ~ `'mass'`      - for contribution to mass fraction
                    ~ `'activity'`  - for activity
                    ~ `'sf'`        - for spontaneous fission
                    ~ `'gsrc'`      - for gamma emission rate
                    ~ `'dh'`        - for decay heat
                    ~ `'ingtox'`    - for ingestion toxicity
                    ~ `'inhtox'`    - for inhalation toxicity
            - string with Serpent 2 special entry.
                Allowed entries are:
                    ~ `'all'`               - all nuclides
                    ~ `'accident'`          - accident relevant nuclides
                    ~ `'actinides'`         - actinides
                    ~ `'burnupcredit'`      - burn-up credit relevant nuclides
                    ~ `'burnupindicators'`  - burn-up indicators
                    ~ `'cosi6'`             - COSI6 code input inventory
                    ~ `'longterm'`          - long term relevant nuclides
                    ~ `'minoractinides'`    - minor actinides
                    ~ `'fp'`                - fission products
                    ~ `'dp'`                - actinides decay products
                    ~ `'ng'`                - noble gases
            The default is `'all'`.

    Methods:
    --------
    clear :
        clears what is written on the input file.
    write :
        writes the model to the input file. This process goes section by section.
    composition_to_restart_file :
        writes the material composition to a binary restart file when te
        simulation is run.
    divisions :
        gets the divisions of the materials in the Simulation lattice(s)

    """
    depletion: Depletion = None
    fpcut: float = 0
    xscalc: int = None
    printm_fraction: tuple[bool, float] = field(default_factory=lambda: (False, 1))
    bumode: str = 'CRAM'
    cram_param: tuple[int, int] = None
    pcc: str | int = 0
    pcc_param: tuple[int, int] = None
    inventory: list | tuple | str = 'all'
    _restart_filename: str = None

    @property
    def _depletion(self) -> str:
        """
        The formatted depletion options for the simulation.
        
        Returns
        -------
        str
            the formatted depleition simulation options.

        """
        string = self._inventory
        string += self._fpcut
        string += self._bumode
        string += self._pcc
        string += self._xscalc
        string += self._printm_fraction
        if self.depletion is not None:
            string += self.comments['Steps'].__str__()
            string += self.depletion.__str__()
        else:
            warnings.warn("No Depletion included in the simulation.")
        string += "/* Material divisions */\n"
        string += self.divisions().to_string() + '\n'
        string += '\n'
        return string

    @property
    def _fpcut(self) -> str:
        """
        Serpent2-formatted fission product cut. set fpcut card.

        Returns
        -------
        str
            the formatted fpcut.

        """
        return f'set fpcut {self.fpcut}\n\n'

    @property
    def _xscalc(self) -> str:
        """
        Serpent2-formatted cross section calculation. set xscalc card.

        Returns
        -------
        str
            the formatted xscalc.

        """
        return f'set xscalc {self.xscalc}' if self.xscalc is not None else ''

    @property
    def _printm_fraction(self) -> str:
        """
        Serpent2-formatted material frcatcion to write to the bumat files. set
        bumat card.

        Returns
        -------
        str
            the formatted printm.

        """
        return f"set printm {int(self.printm_fraction[0])} {self.printm_fraction[1]}\n"

    @property
    def _pcc(self) -> str:
        """
        Serpent2-formatted predictor corrector. set pcc card.

        Returns
        -------
        str
            the formatted pcc.

        """
        pcc_ = AllowedFields.PREDICTOR_CORRECTOR[self.pcc]
        if self.pcc_param is not None:
            parameters = reformat(str(self.pcc_param), "(),")
        else:
            parameters = ''
        return f'set pcc {pcc_} {parameters}\n'

    @property
    def _bumode(self) -> str:
        """
        Serpent2-formatted burnup mode. bumode card.

        Returns
        -------
        str
            the formatted bumode.

        """
        params = reformat(str(self.cram_param), "(),") if self.cram_param is not None else ''
        return f'set bumode {AllowedFields.BUMODE[self.bumode]} ' + params + '\n'

    @property
    def _inventory(self) -> str:
        """
        Serpent2-formatted inventory to include in the simulation output. set
        inventory card.

        Raises
        ------
        TypeError
            If `self.inventory` is not a list, tuple or string.

        Returns
        -------
        str
            the formatted inventory.

        """
        if isinstance(self.inventory, str):
            inventory = self.inventory
        elif isinstance(self.inventory, list):
            if self.inventory[0].isalpha() or self.inventory[0].isnumeric():
                inventory = reformat(str(self.inventory), "[],'").replace(' ', '\n')
            else:
                inventory = reformat(str([nuclide2zam(i.replace('-', '')) for i in self.inventory]), "[],'").replace(
                    ' ', '\n')
        elif isinstance(self.inventory, tuple):
            inventory = 'top ' + reformat(str(self.inventory), "(),'")
        else:
            raise TypeError(f'Inventory should be list, tuple or string, {type(self.inventory)} was given instead.')
        return f'set inventory {inventory}\n'

    def composition_to_restart_file(self, file: str):
        """
        Sets the `restart_filename` parameter to what the user wants and writes
        in the simulation input if it is already written.

        Parameters
        ----------
        file : str
            the name of the target binary restart file.

        Returns
        -------
            None.

        """
        self._restart_filename = file
        if self._written:
            with open(self._file, 'a') as f:
                f.write("/* Composition restart file definition */\n")
                f.write(f"set rfw yes {self._restart_filename}\n")

    def divisions(self) -> DivisionWriter:
        """
        Divisions of the materials in the lattice(s).

        Returns
        -------
        DivisionWriter
            the formatted divisions.

        """
        strings = []
        divisions = []
        for lat in self.geometry.lattices:
            for mat in lat.materials:
                if mat is not None and mat._division_string != '':
                    strings.append(mat._division_string)
                    divisions.extend(mat._divisions)
        return DivisionWriter(strings, DivisionWrapper(divisions))
