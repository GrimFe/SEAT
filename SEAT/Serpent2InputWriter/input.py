import copy as cp
import warnings
import dataclasses
from dataclasses import dataclass, field

from .composition import Composition, MaterialRepresentation, Division
from .geometry import Geometry, Universe
from .depletion import Depletion
from .base import Comment, Other, reformat

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
class Simulation:
    """
    Handles:
    --------
    Hanldes the whole model as an assembly. It is composed of sub-parts corresponding to the model sections.

    Methods:
    --------
    * `clear()`: clears what is written on the input file
    * `write()`: writes the model on the input file. This process goes section by section.
    * `set_BC()`: allows the user to impose the model boundary conditions
    * `set_inventory()`
    * `set_plot_parameters()`
    * `set_seed()`: allows the user to set a seed for simulation reproduction
    * `composition_to_restart_file()`: writes the material composition to a binary restart file
    * `get_divisions()`: gets the divisions of the materials in the Simulation lattice(s)

    Takes:
    ------
    * `geometry`: Geometry object instance - is the geometry of the simulation.
    * `depletion`: Depletion object instance - is the time description of the simulation and the relative power
        normalisation.
    * `composition`: Composition object instance - is the representation of material composition of the simulation.
    * `others`: Others object instance - is the shortcut for what is not implemented yet.

    * `comments`: dictionary of Comment object instances as values. The allowed dictionary keys are:
            - `'Intro'`: Introductory comments
            - `'Geometry'`: Geometry's comments
            - `'Steps'`: Steps' comments
            - `'Materials'`: Materials' comments
            - `'Others'`: Others' comments
            - `'Concluding'`: Concluding comments
            Default is `None`, meaning they are all empty.
    * `seed`: int - seed for the random number generation.
    * `dix`: bool - double indexing method settings. Default is `False`.
    * `fpcut`: float - fission product cumulative yield cutoff in each mass chain. (Serpent2 default is 1, reproduced here)
    * `xscalc`: int - control on how transmutation cross-sections are computed. Default is `None`, corresponding to
        no calculation.
        Other available options are:
            - 1: direct tallies
            - 2: spectrum collapse
    * `printm_fraction`: tuple - a tuple with:
            - boolean: a fleg to set set the printing on/off (Default is False = off)
            - float: minimum atomic fraction to print to the .bumat files. (Default is 1, printing any decay nuclide;
                                                                            0 will print all nuclides but those with transport cross-section)

    * `restart_filename`: string - is the name of the restart file for the material composition of the simulation.
        Default is `None`, which sets the restart file name to restart+self._file in the `write` method.

    * `population`: dictionary - neutron population settings. Default is `None`.
                        Available dictionary keys are:
                            - 'particles': the number of neutrons per generation (Default is 100000)
                            - 'generations': the number of active generations (Default is 250)
                            - 'inactive': the number of inactive generations (default is 30)
                            - 'k_guess': the guess value for k effective (Serpent2 default is 1, reproduced here)
                            - 'batch_interval': the batching interval (Serpent2 default is 1, reproduced here)
                            - 'parallel': the number of independent parallel eigenvalue calculations (Serpent2 default is 1, reproduced here)
    * `gcu`: list - setting of the universes of the group constants. Default is `None`, corresponding to no group
        constant generation.
    * `bumode`: string - burnup calculation mode. Default is `'CRAM'` (Chebyshev Rational Approximation Method),
        for which parameters are set in `cram_param` if not specified, Serpent2 default options will be used.
        Other allowed options are:
            - `'TTA'`: Transmutation Trajectory Analysis
    * `cram_bumode`: tuple - CRAM burnup mode settings. The first value is the order, the second value is the number
        of substeps. Default is `None`.
        Allowed CRAM orders are: 2, 4, 6, 8, 10, 12, 14, 16, -16, -48.
    * `pcc`: int or string - predictor corrector scheme options. Default is 0, corresponding to `'CE'`.
        Allowed options are:
                                  Predictor method            Corrector method
            - `'CE'`   i.e. 0: constant extrapolation
            - `'CELI'` i.e. 1: constant extrapolation       linear interpolation
            - `'LE'`   i.e. 2: linear extrapolation
            - `'LELI'` i.e. 3: linear extrapolation         linear interpolation
            - `'LEQI'` i.e. 4: linear extrapolation         quadratic interpolation
            - `'CECE'` i.e. 6: constant extrapolation       constant backwards extrapolation
        Setting this to None would result in no predictor corrector.
    * `pcc_param`: tuple - setting of the number of steps for the pcc. The first value is for the predictor, the
        second is for the corrector. Default is `None`.
    * `ures`: tuple - containing:
            - boolean: a flag to switch the Unresolved REsonance Probability table Sampling (ures card) on/off (Default is False = off)
            - MaterialRepresentation: the representation of the materials for which ures should be considered (Default is None)
            - float: the dilution cut (Serpent2 default is 1e-9, reproduced here).
    * `bc`: integer, string ot tuple of strings or integers setting the boundary conditions to the simulation.
            The tuple is used to set separate bc for the X, Y, Z directions. (Serpent2 default is 1, reproduced here)
            Options for the string are:
                - `'vacumm'`        or `'V'`   i.e. 1: vacuum boundary conditions
                - `'reflective'`    or `'R'`   i.e. 2: reflective boundary conditions
                - `'periodic'`      or `'P'`   i.e. 3: periodic boundary conditions
    * `albedo`: float - the value of the albedo to set for reflective and/or periodic bc. Default is 1
    * `inventory`: list, string or tuple - nuclides to report in the simulation output.
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
            Default is `'all'`.

    Default internal parameters relevant to the user:
    -------------------------------------------------
    * `_file`: string - name of the file to which the simulation is written. Default is `None`.
    * `_written`: bool - indicates whether the simulation has already been written. Default is `False`.
    """
    geometry: Geometry
    depletion: Depletion
    composition: Composition
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
        string += self.comments['Geometry'].__str__()
        string += self.geometry.__str__()
        string += self.comments['Steps'].__str__()
        string += self.depletion.__str__()
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
        return f"set seed {self.seed}\n" if self.seed is not None else ''

    @property
    def _dix(self) -> str:
        return f'set dix {int(self.dix)}\n'

    @property
    def _fpcut(self) -> str:
        return f'set fpcut {self.fpcut}\n\n'

    @property
    def _xscalc(self) -> str:
        return f'set xscalc {self.xscalc}' if self.xscalc is not None else ''

    @property
    def _printm_fraction(self) -> str:
        return f"set printm {int(self.printm_fraction[0])} {self.printm_fraction[1]}\n"

    @property
    def _pcc(self) -> str:
        pcc_ = PREDICTOR_CORRECTOR_MODES[self.pcc]
        if self.pcc_param is not None:
            parameters = reformat(str(self.pcc_param), "(),")
        else:
            parameters = ''
        return f'set pcc {pcc_} {parameters}\n'

    @property
    def _bumode(self) -> str:
        if self.bumode == 'TTA':
            out = '1'
        elif self.cram_param is not None:
            out = '2 ' + reformat(str(self.cram_param), "(),")
        else:
            out = '2'
        return 'set bumode ' + out + '\n'

    @property
    def _gcu(self) -> str:
        if self.gcu is None:
            out = '-1'
        else:
            out = reformat(str([i.name for i in self.gcu]), "[],'")
        return 'set gcu ' + out + '\n\n'

    @property
    def _population(self) -> str:
        """
        Takes:
        ------
        * `population`: dictionary.
                        Available dictionary keys are:
                            - 'particles': the number of neutrons per generation (Default is 100000)
                            - 'generations': the number of active generations (Default is 250)
                            - 'inactive': the number of inactive generations (Default is 30)
                            - 'k_guess': the guess value for k effective (Default is 1, reproduced here)
                            - 'batch_interval': the batching interval (Default is 1, reproduced here)
                            - 'parallel': the number of independent parallel eigenvalue calculations (Default is 1, reproduced here)
        """
        pop = f"{self.population['particles']} {self.population['generations']} {self.population['inactive']} " \
              f"{self.population['k_guess']} {self.population['batch_interval']} {self.population['parallel']}\n"
        return 'set pop ' + pop

    @property
    def _ures(self) -> str:
        """
        Takes:
        ------
        * `ures`: tuple containing:
                - boolean: a flag to switch the ures card on/off (Default is False = off)
                - MaterialRepresentation: the representation of the materials for which ures should be considered (Default is None)
                - float: the dilution cut (Serpent2 default is 1e-9, reproduced here).
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
        Model boundary conditions.
        Takes:
        ------
        * `bc`: integer, string ot tuple of strings or integers setting the boundary conditions to the simulation.
            The tuple is used to set separate bc for the X, Y, Z directions. (Serpent2 default is 1, reproduced here)
            Options for the string are:
                - `'vacumm'`        or `'V'`   i.e. 1: vacuum boundary conditions
                - `'reflective'`    or `'R'`   i.e. 2: reflective boundary conditions
                - `'periodic'`      or `'P'`   i.e. 3: periodic boundary conditions
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
        Nuclide inventory to include in simulation output.

        Takes:
        ------
        * `inventory`: can be:
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
            Default is `'all'`.
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
        Clears what is written on the input file
        """
        if self._file is not None:
            with open(self._file, 'w'):
                pass
        else:
            warnings.warn("No output file is set for the simulation object")

    def write(self, file):
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

        Takes:
        ------
        * `file`: string - name of the file to write the simulation to
        """
        if self._file is not None:
            warnings.warn("Ignoring a previously defined file was already defined for this simulation.")
        self._file = file
        with open(self._file, 'w') as f:
            f.write(self.__str__())
        self._written = True

    def composition_to_restart_file(self, file: str):
        """
        Sets the `restart_filename` parameter to what the user wants and writes in the simulation input if it is
        already written.

        Takes:
        ------
        * `file`: string - name of the target binary restart file.
        """
        self._restart_filename = file
        if self._written:
            with open(self._file, 'a') as f:
                f.write("/* Composition restart file definition */\n")
                f.write(f"set rfw yes {self._restart_filename}\n")

    def divisions(self):
        """
        Divisions of the materials in the lattice(s)
        """
        strings = []
        divisions = []
        for lat in self.geometry.lattices:
            for mat in lat.get_materials():
                if mat is not None and mat.division_string != '':
                    strings.append(mat.division_string)
                    divisions.extend(mat.divisions)
        return DivisionWriter(strings, DivisionWrapper(divisions))


@dataclass(slots=True)
class DivisionWrapper:
    """
    Handles the labelling of the divisions.

    Takes:
    ------
    * `div`: list - list of divisions to label and sort.
    """

    divisions: list[Division]

    def format_sort(self) -> str:
        """
        Labels and sorts the divisions in the wrapper writing them to a string.
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
    Handles the writing process of the divisions

    Takes:
    ------
    * `strings`: list - list of division strings.
    * `divisions`: DivisionWrapper object instance - the divisions properly wrapped.
    """
    strings: list[str]
    divisions: DivisionWrapper

    def to_string(self):
        string = reformat(str(set(self.strings)), "{}'").replace(', ', '\n') + '\n'
        if None not in {self.divisions}:
            string += '\n'
            string += "set mvol\n\n"
            string += self.divisions.format_sort()
        return string
