import Serpent2InputWriter.composition
from dataclasses import dataclass


@dataclass(slots=True)
class NaturalAbundances:
    """
    Stores the natural abundances of a number of commonly used elements.
    Stored as dictionaries of nuclide symbols to percentage abundance.
    """
    Pb = {'Pb204': -1.378e-2, 'Pb206': -23.956e-2, 'Pb207': -22.074e-2, 'Pb208': -52.592e-2}


@dataclass(slots=True)
class MaterialCompositions:
    """
    Stores predefined material compositions.
    """
    Pb = Serpent2InputWriter.composition.MaterialComposition.from_nuclides(NaturalAbundances().Pb)
