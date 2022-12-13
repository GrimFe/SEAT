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
    Handles:
    --------
    Handles the flux normalization option and relative value to which normalise

    Methods:
    --------
    * `write()`: internal method to write on a file with proper formatting for Serpent 2 input normalization condition
        definition
    * `copy()`: copies the object instance to another memory allocation

    Takes:
    ------
    * `value`: float - is either the power density (kW/g), the power (W) or anything needed for flux normalization
    * `material`: Material object instance - is the Material in which the normalization power is generated.
        `material == ''` then `power` will be referred to the average power produced in the system.
    * `kind`: string - is the type of normalization imposed. Default is `'powdens'`. It can be:
        - `'powdens'`: normalization to power density
        - `'power'`: normalization to total fission power
        - `'flux'`: normalization to total flux
        - `'genrate'`: normalization to neutron generation rate
        - `'fissrate'`: normalization to fission rate
        - `'absrate'`: normalization to absorption rate
        - `'lossrate'`: normalization to loss rate
        - `'srcrate'`: normalization to source rate
        - `'sfrate'`: normalization to spontaneous fission rate
    """
    value: float
    material: Material
    kind: str = 'powdens'

    def __str__(self):
        string = f'set {self.kind} {self.value} {self.material.name}\n'
        return string

    def write(self, file: str):
        """
        Internal method to write on a file with proper formatting for Serpent 2 input normalization definition

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
class Interval:
    """
    Handles:
    --------
    Handles the time interval over which a calculation step takes place

    Methods:
    --------
    * `write()`: internal method to write on a file with proper formatting for Serpent 2 input time interval definition
    * `copy()`: copies the object instance to another memory allocation

    Takes:
    ------
    * `span`: float - is the duration of the interval
    * `span_type`: string - is the type of the desired time interval, defining its units as well.
        It can either be:
            *`'bu'` for burn-up interval (MWd/kgHM)
            *`'day'` for time interval (days)
            *`'dec'` for decay interval (days)
            *`'act'` for activation interval (days)
    * `total`: bool - is the boolean identifying if the step is integral or not.
        Default value is `False`, corresponding to non-total step
    """
    span: float
    span_type: str = 'day'
    total: bool = False

    @property
    def step(self):
        step_ = 'tot' if self.total else 'step'
        return step_

    def __str__(self):
        string = f"dep {self.span_type}{self.step} {self.span}\n"
        return string

    def write(self, file: str):
        """
        Internal method to write on a file with proper formatting for Serpent 2 input syntax

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
class Depletion:
    """
    Handles:
    --------
    Handles the coupling of depletion time step and relative normalization

    Methods:
    --------
    * `write()`: internal method to write on a file with proper formatting for Serpent 2 input syntax

    Takes:
    ------
    * `steps`: list - is a list of tuples linking a Normalization object instance and an Interval object instance,
        in this order
    """
    steps: list

    def __str__(self):
        string = ''
        for n, i in self.steps:
            string += n.__str__()
            string += i.__str__()
        string += "\n"
        return string

    def write(self, file: str):
        """
        Internal method to write on a file with proper formatting for Serpent 2 input comment, inserted in '/*' '*/'

        Takes:
        ------
        * `file`: string - is the name of the file where to write
        """
        with open(file, 'a') as f:
            f.write(self.__str__())
