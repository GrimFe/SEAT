import copy as cp
import warnings
import dataclasses
from dataclasses import dataclass, field

from SEAT.Serpent2InputWriter.composition import Composition, MaterialRepresentation, Division
from SEAT.Serpent2InputWriter.geometry import Geometry, Universe
from SEAT.Serpent2InputWriter.depletion import Depletion
from SEAT.Serpent2InputWriter.base import Comment, Other, reformat

from SEAT.nuclides import nuclide2zam

__author__ = "Federico Grimaldi"
__all__ = [
    "Simulation",
]

EMPTY_COMMENT = Comment('')
COMMENT_BASE_DCT: dict[str, Comment] = {'Intro': EMPTY_COMMENT,
                                        'Geometry': EMPTY_COMMENT,
                                        'Steps': EMPTY_COMMENT,
                                        'Materials': EMPTY_COMMENT,
                                        'Others': EMPTY_COMMENT,
                                        'Concluding': EMPTY_COMMENT}
POPULATION_DEFAULT: dict[str, int] = {'particles': 100000,
                                      'generations': 250,
                                      'inactive': 30,
                                      'k_guess': 1,
                                      'batch_interval': 1,
                                      'parallel': 1}
API_LINK = 'https://github.com/GrimFe/Serpent2InputWriter.git'
HEADER_COMMENT = f"""/*
#############################################################################################################################################
#                                                                                                                                           #
#                                            Serpent2 input file written with SOUP 0.0.1 python API                                         #
#                                               https://github.com/GrimFe/Serpent2InputWriter.git                                           #
#                                                                                                                                           #
#############################################################################################################################################
*/\n\n"""

PREDICTOR_CORRECTOR_MODES = {'CE': 0,
                             'CELI': 1,
                             'LE': 2,
                             'LELI': 3,
                             'LEQI': 4,
                             'CECE': 6,
                             0: 0,
                             1: 1,
                             2: 2,
                             3: 3,
                             4: 4,
                             5: 5,
                             6: 6}
BC_MODES = {'vacuum': 1,
            'reflective': 2,
            'periodic': 3,
            'V': 1,
            'R': 2,
            'P': 3,
            1: 1,
            2: 2,
            3: 3}

mvol_counter = {}


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
class Simulation:
    """
    Hanldes the model as a whole allowing for the definition of simulation
    specifications. It is composed of sub-parts corresponding to the model
    sections.

    Attributes
    ----------
    geometry : `SEAT.Geometry`, optional
        the geometry of the simulation. The default is None.
    composition : `SEAT.Composition`, optional
        the materials and nuclear data used in the simulation. The default is
        None.
    depletion : `SEAT.Depletion`, optional
        the neutron flux normalization and the time definition. The default is
        None.
    others : `SEAT.Others`, optional
        for what is not implemented yet. The default is None
    comments : dict[str, `SEAT.Comment`], optional
        comments to the simulation sections.    
        Allowed dictionary keys are:
            - 'Intro': introductory comments
            - 'Geometry': geometry comments
            - 'Steps': depletion comments
            - 'Materials': composition comments
            - 'Others': other comments
            - 'Concluding': concluding comments
            The default is None (all comments are empty).
    seed : int, optional
        seed for the random number generation. The default is None.
    dix : bool
        double indexing method settings. The default is False.
    fpcut : float, optional
        the cumulative fission product yield cutoff in each mass chain. The
        default is 1.
    xscalc : int, optional
        control on how transmutation cross-sections are computed.
        Allowed `xscalc` vlaues are:
            - None: no calculation
            - 1: direct tallies
            - 2: spectrum collapse
        The default is None.
    printm_fraction : tuple[bool, float], optional
        control on the material composition printed to '.bumat' files.
        It composes of an on/off switch and of the minimum atomic fraction to
        print: 1 prints any decay nuclide, 0 prints those that don't have a
        transport cross section. The default is (False, 1).
    restart_filename: str
        the name of the restart file for the material composition.
        The default is None (restart file name = restart+self._file in the
                             `write` method).
    population : dict[str, float]
        the particle population settings in the simulation.
        Allowed dictionary keys are:
            - 'particles': the number of neutrons per generation.
                            The default is 100000.
            - 'generations': the number of active generations.
                            The default is 250
            - 'inactive': the number of inactive generations.
                            The default is 30.
            - 'k_guess': the guess value for k effective.
                            The default is 1.
            - 'batch_interval': the batching interval.
                            The default is 1.
            - 'parallel': the number of independent parallel eigenvalue
                            calculations. The default is 1.
    gcu : list[`SEAT.Universe`]
        the universes of the group constants. The default is None (no generation).
    bumode : str, optional
        burnup calculation mode.
        Allowed `bumode` options are:
            - 'CRAM': Chebyshev Rational Approximation Method (parameters in
                    `cram_bumode`)
            - 'TTA': Transmutation Trajectory Analysis
        The default is 'CRAM'.
    cram_bumode: tuple[int]
        CRAM burnup mode settings. The first items is the order, the second
        is the number of substeps.
        Allowed CRAM order values are:
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
    ures : tuple[bool, `SEAT.MaterialRepresentation`, float], optional
        - bool: a flag to turn the Unresolved REsonance Probability table 
                Sampling (ures) on/off. The default is False = off.
        - MaterialRepresentation: the representation of the materials for
                which ures should be considered. The default is None.
        - float: the dilution cut. The default is 1e-9.
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
    inventory : list | string | tuple, optional
        nuclides to report in the simulation output.
        It can be:
            - list of materials to include in Isotope form `['U-235', 'Am-242']` or `['U235', 'Am242']`,
                ZAI form `[922350, 952421]`,
                Element in letters or mass number `['U', 'Am']` or `[92, 95]`
                All the items in the list should follow the same notation.
            - tuple of number of nuclides to include and parameter according to which choose them.
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
    set_BC :
        allows the user to impose the model boundary conditions.
    set_inventory:
        allows the user to set the simulation inventory.
    set_plot_parameters:
        allows the user to set the plot parameters.
    set_seed :
        allows the user to set a seed for simulation reproduction.
    composition_to_restart_file :
        writes the material composition to a binary restart file when te
        simulation is run.
    divisions :
        gets the divisions of the materials in the Simulation lattice(s)

    """
    geometry: Geometry = None
    composition: Composition = None
    depletion: Depletion = None
    others: Other = None
    comments: dict[str, Comment] = field(default_factory=lambda: COMMENT_BASE_DCT)
    seed: int = None
    dix: bool = False
    fpcut: float = 0
    xscalc: int = None
    printm_fraction: tuple[bool, float] = field(default_factory=lambda: (False, 1))
    population: dict[str, int] = field(default_factory=lambda: POPULATION_DEFAULT)
    gcu: list[Universe] = None
    bumode: str = 'CRAM'
    cram_param: tuple[int, int] = None
    pcc: str | int = 0
    pcc_param: tuple[int, int] = None
    ures: tuple[bool, MaterialRepresentation, float] = field(default_factory=lambda: (False, None, 0.000000001))
    inventory: list | tuple | str = 'all'
    bc: tuple[str | int] | str | int = field(default_factory=lambda: 1)
    albedo: float = 1
    _file: str = None
    _written: bool = False
    _restart_filename: str = None

    def __str__(self):
        string = HEADER_COMMENT
        string += self.comments['Intro'].__str__()
        string += Comment("Global simulation settings").__str__()
        string += self._seed
        string += self._population
        string += self._bc
        string += self._inventory
        string += '\n'
        string += Comment("Unresolved resonance probability table settings").__str__()
        string += self._ures
        string += Comment("Double indexing method settings").__str__()
        string += self._dix
        string += Comment("Fission products ut-off settings").__str__()
        string += self._fpcut
        string += Comment("Group constant generation settings").__str__()
        string += self._gcu
        string += Comment("Burnup calculation settings").__str__()
        string += self._bumode
        string += self._pcc
        string += self._xscalc
        string += self._printm_fraction
        string += '\n'
        if self.geometry is not None:
            string += self.comments['Geometry'].__str__()
            string += self.geometry.__str__()
        if self.depletion is not None:
            string += self.comments['Steps'].__str__()
            string += self.depletion.__str__()
        if self.composition is not None:
            string += self.comments['Materials'].__str__()
            string += self.composition.__str__()
        if self.others is not None:
            string += self.comments['Others'].__str__()
            string += '\n'
            string += self.others.__str__()
        string += '\n'
        string += "/* Material divisions */\n"
        string += self.divisions().to_string() + '\n'
        string += '\n'
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
        pcc_ = PREDICTOR_CORRECTOR_MODES[self.pcc]
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
        if self.bumode == 'TTA':
            out = '1'
        elif self.cram_param is not None:
            out = '2 ' + reformat(str(self.cram_param), "(),")
        else:
            out = '2'
        return 'set bumode ' + out + '\n'

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
        pop = POPULATION_DEFAULT | self.population
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
        set_ures = self.ures[0]
        if set_ures:
            if isinstance(self.ures[1], MaterialRepresentation):
                nuclides = ' '
                nuclides += reformat(str(self.ures[1].get_nuclide_representation()), "[],'")
                nuclides += ' '
            else:
                nuclides = ' '
            dilcut = self.ures[2]
            out = f"set ures {int(set_ures)}{nuclides}{dilcut}\n\n"
        else:
            out = f"set ures {int(set_ures)}\n\n"
        return out

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
            mode_ = [BC_MODES[bc] for bc in self.bc]
            mode = reformat(str(mode_), "[],")
            _albedo = f' {self.albedo}' if (BC_MODES['R'] in mode_ or BC_MODES['P'] in mode_) else ''
        else:
            mode = BC_MODES[self.bc]
            _albedo = f' {self.albedo}' if (mode == BC_MODES['R'] or mode == BC_MODES['P']) else ''
        return f'set bc {mode}{_albedo}\n'

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
