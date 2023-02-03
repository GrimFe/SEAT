import copy as cp
from .composition import Material
from dataclasses import dataclass

__author__ = "Federico Grimaldi"
__all__ = [
    "Normalization",
    "Interval",
    "Depletion",
]


@dataclass(slots=True)
class Normalization:
    """
    Class for the flux normalization.

    Attributes
    ----------
    value : float
        the power density (kW/g), the power (W) or anything needed for flux
        normalization.
    material : SEAT.Material
        the `SEAT.Material` instance in which the power is normalized.
        If `material == ''` then the normalization will be referred to the
        average power produced in the system.
    kind : str
        the type of normalization imposed. The default is 'powdens'.
        Allowed `kind` values are:
            - 'powdens': normalization to power density
            - 'power': normalization to total fission power
            - 'flux': normalization to total flux
            - 'genrate': normalization to neutron generation rate
            - 'fissrate': normalization to fission rate
            - 'absrate': normalization to absorption rate
            - 'lossrate': normalization to loss rate
            - 'srcrate': normalization to source rate
            - 'sfrate': normalization to spontaneous fission rate

    Methods:
    --------
    write :
        writes the `SEAT.Normalization` instance to a file.
    copy :
        copies the object instance to another memory allocation.

    """
    value: float
    material: Material
    kind: str = 'powdens'

    def __str__(self):
        string = f'set {self.kind} {self.value} {self.material.name}\n'
        return string

    def write(self, file: str, mode: str='a'):
        """
        Writes the `SEAT.Normalization` instance to a file.

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

        Returns:
        --------
        SEAT.Other
            a copy of the `SEAT.Other` instance.

        """
        return cp.deepcopy(self)


@dataclass(slots=True)
class Interval:
    """
    Class for the normalization interval.

    Attributes
    ----------
    span : float
        the duration of the interval.
    total : bool
        identifyes whether the step is integral. The default is False.
    kind : str
        the type of interval chosen, which also defines the `span` units.
        The default is 'day'.
        Allowed `kind` values are:
            - 'bu': burn-up interval (MWd/kgHM)
            - 'day': time interval (days)
            - 'dec': decay interval (days)
            - 'act': activation interval (days)

    Methods:
    --------
    write :
        writes the `SEAT.Interval` instance to a file.
    copy :
        copies the object instance to another memory allocation.

    """
    span: float
    total: bool = False
    span_type: str = 'day'

    @property
    def step(self) -> str:
        """
        The type of step chosen: total or step.

        Returns
        -------
        step_ : str
            Serpent 2 identifier for the selected step type.

        """
        step_ = 'tot' if self.total else 'step'
        return step_

    def __str__(self):
        string = f"dep {self.span_type}{self.step} {self.span}\n"
        return string

    def write(self, file: str, mode: str='a'):
        """
        Writes the `SEAT.Interval` instance to a file.

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

        Returns:
        --------
        SEAT.Other
            a copy of the `SEAT.Other` instance.

        """
        return cp.deepcopy(self)


@dataclass(slots=True)
class Depletion:
    """
    Class for the depletion, coupling `SEAT.Normalization` and `SEAT.Interval`.

    Attributes
    ----------
    steps : list[tuple]
        couples a `SEAT.Normalization` instance and a `SEAT.Interval` instance.

    Methods:
    --------
    write :
        writes the `SEAT.Interval` instance to a file.

    """
    steps: list[tuple]

    def __str__(self):
        string = ''
        for n, i in self.steps:
            string += n.__str__()
            string += i.__str__()
        string += "\n"
        return string

    def write(self, file: str, mode: str='a'):
        """
        Writes the `SEAT.Depletion` instance to a file.

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
