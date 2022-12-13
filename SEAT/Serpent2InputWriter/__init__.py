from SEAT.Serpent2InputWriter.base import *
from SEAT.Serpent2InputWriter.input import *
from SEAT.Serpent2InputWriter.geometry import *
from SEAT.Serpent2InputWriter.composition import *
from SEAT.Serpent2InputWriter.depletion import *

__doc__ = """
Usage:
    * Definition of the composition:
        - Create materials
        - Split materials if needed
        - Assemble composition
        
    * Definition of the Geometry:
        - Create pins
        - Create surfaces
        - Create universes
        - Create cells
        - Create representation
        - Create lattice
        - Assemble geometry
        
    * Definition of the depletion:
        - Create normalisation
        - Create time intervals
        - Assemble Depletion
    
    * Definition of Others
    
    * Definition of section comments
    
    * Definition of Simulation
    
    * Write simulation
"""

EMPTY_COMMENT = Comment("")
