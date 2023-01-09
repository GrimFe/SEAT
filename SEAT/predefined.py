import Serpent2InputWriter.composition
from dataclasses import dataclass


@dataclass(slots=True)
class NaturalAbundances:
    """
    Stores the natural abundances of a number of commonly used elements.
    Stored as dictionaries of nuclide symbols to percentage abundance.
    """
    Pb = {'Pb204': 1.378e-2, 'Pb206': 23.956e-2, 'Pb207': 22.074e-2, 'Pb208': 52.592e-2}
    U = {'U238': 99.274e-2, 'U235': 0.720e-2, 'U234': 0.005e-2}


@dataclass(slots=True)
class MaterialCompositions:
    """
    Stores predefined material compositions.
    """
    Pb = Serpent2InputWriter.composition.MaterialComposition.from_nuclides(NaturalAbundances().Pb)
    U = Serpent2InputWriter.composition.MaterialComposition.from_nuclides(NaturalAbundances().U)

    @staticmethod
    def U_enriched(e: float | dict[str, float]) -> Serpent2InputWriter.composition.MaterialComposition:
        """
        Computes the material composition of enriched uranium.
        * e: int or dict. Is the uranium enrichment:
            - int: enrichment in U235 only (percentage value)
            - dict: dictionaries having uranium isotope symbols as keys and percentage enrichment as values
        """
        if type(e) == int:
            cmp = {'U235': e * 1e-2, 'U238': 1 - e * 1e-2}
        elif type(e) == dict:
            cmp = {k: v * 1e-2 for k, v in e.items()}
            cmp['U238'] = 1 - sum(e.values()) * 1e-2
        else:
            raise TypeError(f"The enrichment can be either int or dict, {type(e)} was given instead.")
        return Serpent2InputWriter.composition.MaterialComposition.from_nuclides(cmp)
