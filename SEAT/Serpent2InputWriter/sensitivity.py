import SEAT.Serpent2InputWriter._DefaultOptions as DefaultOptions
import warnings

from dataclasses import dataclass, field
from SEAT.Serpent2InputWriter.detectors import Detector
from SEAT.Serpent2InputWriter.composition import Material

__author__ = "Federico Grimaldi"
__all__ = [
    "Sensitivity"
    ]

@dataclass(slots=True)
class Detratio:
    """
    Class for the detector ratios.

    Attributes
    ----------
    numerator : `SEAT.Detector`
        the detector ratio numerator.
    denominator : `SEAT.Detector`
        the detector ratio denominator.
    wn : float, optional
        the weight of the detector ratio numerator. The default is 1.    
    wd : float, optional
        the weight of the detector ratio denominator. The default is 1.

    """
    numerator: Detector
    denominator: Detector
    wn: float = 1
    wd: float = 1

    def __str__(self):
        string = ''
        string += f'{self.numerator.name} '
        string += f'{self.denominator.name} '
        if self.wn != 1 or self.wd != 1:
            string += f'{self.wn} '
            string += f'{self.wd}'
        return string


@dataclass(slots=True)
class Response:
    """
    Class for the sensitivity response.

    Attributes
    ----------
    base : dict[str, bool], optional
        * keys: base perturbation identifier.
        * values: flag to activate it.
        Allowed dictionary keys are:
            - 'keff': effective multiplication factor
            - 'beff': effective delayed neutron fraction
            - 'leff': effective prompt generation time
            - 'lambda': delayed neutron precursor decay constants
            By default all base responses are set to False (suppressed).
    detratios : dict[str, `SEAT.Detratio`], optional
        * keys: name of the detector ratio response.
        * values: detector ratios `SEAT.Detratio` objects.
        The default is None.
    void : list[`SEAT.Material`], optional
        The materials to which voiding the sensitivity shall be computed. The
        default is None.
    lambda_total: bool, optional
        response specification for the delayed neutron precursor decay constants
        as total. The default is False.
    lambda_groupwise: bool, optional
        response specification for the delayed neutron precursor decay constants
        as groupwise. The default is False.

    """
    base: dict[str, bool] = field(default_factory=lambda: DefaultOptions.RESP)
    detratios: dict[str, Detratio] = None
    voids: list[Material] = None
    lambda_total: bool = False
    lambda_groupwise: bool = False

    def __str__(self):
        string = ''
        set_resp = 'sens resp '
        string += set_resp + 'keff\n' if self.base['keff'] else ''
        string += set_resp + 'beff\n' if self.base['beff'] else ''
        string += set_resp + 'leff\n' if self.base['leff'] else ''
        if self.base['lambda']:
            string += set_resp + 'lambda '
            string += f'{int(self.lambda_total)} {int(self.lambda_groupwise)}\n'
            if not self.lambda_total and not self.lambda_groupwise:
                warnstring = "No response specification was defined for the"\
                    "delayed neutron precursor decay constants. This therefore"\
                    "results to be suppressed."
                warnings.warn(warnstring)
        if self.detratios is not None:
            for dr in self.detratios:
                string += dr.__str__() + '\n'
        if self.voids is not None:
            for m in self.voids:
                string += f'{m.name}\n'
        return string


@dataclass(slots=True)
class _BasePerturbation:
    """
    Base class for the perturbations.

    Attributes
    ----------
    realist: list[str], optional
        the string identifiers for the cross section.
        Allowed values are:
            - 'ela': elastic scattering cross section
            - 'sab': subthermal scattering cross section
            - 'inl': inelastic scattering cross section (MTs 51-91)
            - 'capt': capture cross section
            - 'fiss': fission cross section
            - 'nxn': non fission reactions with more than one outgoing neutron
            - ['all']: sum reaction modes; exclusive
            - ['allmt']: partial reaction modes (MT wise perturbation); exclusive
            The defaul value is None.
        Alternative to `mtlist`.
    mtlist: list[int], optional
        the MT of the reactions to perturb. The default is None.
        Alternative to realist.
    zailist: list[str | int], optional
        the nuclides to perturb. Can be in the form of nuclide symblos (string)
        or of ZAM (integer). The default is 'all', perturbing all the nuclides.
    matlist: list[`SEAT.Material`], optional
        the materials where the perturbations are applied. The default is 'all',
        peturbing all the materials.

    Note
    ----
    `realist` and `mtlist` are alternative. "Sum reaction modes and partial
    reaction modes cannot currently be perturbed at the same time (in the same
    input). For inputs containing both perturbations (e.g. a sens pert xs
    realist card and a sens pert xs mtlist card) only the latter perturbation
    is used."

    """
    realist: list[str] = None
    mtlist: list[int] = None
    zailist: list[str | int] = 'all'
    matlist: list[Material] = 'all'


@dataclass(slots=True)
class Perturbation(_BasePerturbation):
    """
    Class for the sensitivity perturbation.

    Attributes
    ----------
    realist: list[str], optional
        the string identifiers for the cross section.
        Allowed values are:
            - 'ela': elastic scattering cross section
            - 'sab': subthermal scattering cross section
            - 'inl': inelastic scattering cross section (MTs 51-91)
            - 'capt': capture cross section
            - 'fiss': fission cross section
            - 'nxn': non fission reactions with more than one outgoing neutron
            - ['all']: sum reaction modes; exclusive
            - ['allmt']: partial reaction modes (MT wise perturbation); exclusive
            The defaul value is None.
        Alternative to `mtlist`.
    mtlist: list[int], optional
        the MT of the reactions to perturb. The default is None.
        Alternative to realist.
    zailist: list[str | int], optional
        the nuclides to perturb. Can be in the form of nuclide symblos (string)
        or of ZAM (integer). The default is 'all', perturbing all the nuclides.
    matlist: list[`SEAT.Material`], optional
        the materials where the perturbations are applied. The default is 'all',
        peturbing all the materials.
    base : dict[str, bool], optional
        * keys: base perturbation identifier.
        * values: flag to activate it.
        Allowed dictionary keys are:
            - 'xs': basic cross sections
            - 'chi': fission spectrum
            - 'nubar': fission nubar
            - 'elamu': elastic scattering cosine
            - 'inlmu': inelastic scattering cosine
            - 'eleg': Legendre moments of elastic scattering angular distribution
            - 'temperature': temperature
            By default all base perturbations are set to False (suppressed).
    nleg: int, optional
        specification to the eleg base perturbation. Number of Legendre moments
        to perturb. Can be up to 7. The default is 1.

    Note
    ----
    `realist` and `mtlist` are alternative. "Sum reaction modes and partial
    reaction modes cannot currently be perturbed at the same time (in the same
    input). For inputs containing both perturbations (e.g. a sens pert xs
    realist card and a sens pert xs mtlist card) only the latter perturbation
    is used."

    """
    base: dict[str, bool] = field(default_factory=lambda: DefaultOptions.PERT)
    nleg: int = 1

    def __str__(self):
        string = ''
        set_pert = 'sens pert '

        perturb_xs = False
        if self.realist is not None:
            if 'allmt' in self.realist:
                string += set_pert + 'xs allmt\n'
                perturb_xs = True
            elif 'all' in self.realist:
                string += set_pert + 'xs all\n'
                perturb_xs = True
            else:
                string += set_pert + f'xs ....\n'
                perturb_xs = True
        if self.mtlist is not None:
            string += set_pert + f'xs ....\n'
            perturb_xs = True

        if self.realist is not None and self.mtlist is not None:
            warnstring = "`realist` and `mtlist` are exclusive."
            warnings.warn(warnstring)
        if not self.base['xs'] and perturb_xs:
            warnstring = "A xs perturbation is imposed but it is set to false"\
                "the base perturbation."
            warnings.warn(warnstring)

        string += set_pert + 'chi\n' if self.base['chi'] else ''
        string += set_pert + 'nubar\n' if self.base['nubar'] else ''
        string += set_pert + 'elamu\n' if self.base['elamu'] else ''
        string += set_pert + 'inlmu\n' if self.base['inlmu'] else ''
        string += set_pert + f'eleg {self.nleg}\n' if self.base['eleg'] else ''
        string += set_pert + 'temperature\n' if self.base['temperature'] else ''

        if self.zailist == 'all':
            string += set_pert + 'zailist all\n'
        else:
            string += ' '.join(self.zailist) + '\n'

        if self.matlist == 'all':
            string += set_pert + 'matlist all\n'
        else:
            string += ' '.join([m.name for m in self.matlist]) + '\n'
        return string        


@dataclass(slots=True)
class CustomPerturbation(_BasePerturbation):
    """
    Class for custom sensitivity perturbations.
    
    Attributes
    ----------
    realist: list[str], optional
        the string identifiers for the cross section.
        Allowed values are:
            - 'ela': elastic scattering cross section
            - 'sab': subthermal scattering cross section
            - 'inl': inelastic scattering cross section (MTs 51-91)
            - 'capt': capture cross section
            - 'fiss': fission cross section
            - 'nxn': non fission reactions with more than one outgoing neutron
            - ['all']: sum reaction modes; exclusive
            - ['allmt']: partial reaction modes (MT wise perturbation); exclusive
            The defaul value is None.
        Alternative to `mtlist`.
    mtlist: list[int], optional
        the MT of the reactions to perturb. The default is None.
        Alternative to realist.
    zailist: list[str | int], optional
        the nuclides to perturb. Can be in the form of nuclide symblos (string)
        or of ZAM (integer). The default is 'all', perturbing all the nuclides.
    matlist: list[`SEAT.Material`], optional
        the materials where the perturbations are applied. The default is 'all',
        peturbing all the materials.
    name: str
        the custom perturbation name
    efunc: str
        the name of the file where the energy dependent perturbation is stored.

    Note
    ----
    `realist` and `mtlist` are alternative. "Sum reaction modes and partial
    reaction modes cannot currently be perturbed at the same time (in the same
    input). For inputs containing both perturbations (e.g. a sens pert xs
    realist card and a sens pert xs mtlist card) only the latter perturbation
    is used."

    `name` and `efunc` are required parameters.

    """
    name = None
    efunc = None

    def __str__(self):
        string = f'sens pert custom {self.name} {self.efunc}'
        string += ' zailist ' + ' '.join(self.zailist) + '\n' if self.zailist else ''
        string += ' matlist ' + ' '.join(self.matlist) + '\n' if self.matlist else ''
        string += ' mtlist ' + ' '.join([str(mt) for mt in self.matlist]) +\
            '\n' if self.matlist else ''
        string += ' realist ' + ' '.join(self.realist) if self.realist and\
            'all' not in self.realist and 'allmt' not in self.realist else ''
        return string + '\n'


@dataclass(slots=True)
class Sensitivity:
    """
    Class for the Sensitivity coupling `SEAT.Response` and `SEAT.Perturbation`.

    Attributes
    ----------
    response: `SEAT.Response`
        the response to monitor
    perturbation: `SEAT.Perturbation`, optional
        the perturbation to apply. The default is None.
    custom_perturbation: list[`SEAT.CustomPerturbation`], optional
        the custom perturbations to apply. The default is None.

    """
    response: Response
    perturbation: Perturbation = None
    custom_perturbation: CustomPerturbation = None

    def __str__(self):
        string = '--- Responses ---\n'
        string += self.response.__str__() + '\n'
        string += '--- Perturbations ---\n'
        if self.perturbation is None and self.custom_perturbation is None:
            warnstring = "No perturbation included"
            warnings.warn(warnstring)
        string += self.perturbation.__str__() if self.perturbation\
            is not None else ''
        string += self.custom_perturbation.__str__() if self.custom_perturbation\
            is not None else ''
        return string
