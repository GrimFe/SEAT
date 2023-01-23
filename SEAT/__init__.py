import SEAT.Serpent2InputWriter
import SEAT.Serpent2InputParser
# import SOUP.MultipleOutputReader
from SEAT import nuclides
from SEAT import composites
from SEAT import natural
import SEAT.composition_functions

__all__ = [
    "Serpent2InputWriter",
    "Serpent2InputParser",
    "nuclides"
]

__doc__ = """
Usage:
"""

# when changing version also adapt the default string in the input.py file
__version__ = '0.0.1'


class SOUPException(Exception):
    pass
