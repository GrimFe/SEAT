import copy as cp
from dataclasses import dataclass

__author__ = "Federico Grimaldi"
__all__ = [
    "intersect",
    "unite",
    "surface_complement",
    # "cell_complement",
    "Comment",
    "StandaloneComment",
    "InlineComment",
    "Entity",
    "Other",
]


def reformat(string: str, reform: str) -> str:
    """
    Cleans a string from what is to be reformed.

    :param string: 'str' - the string to reformat
    :param reform: 'str' - the string with the items to remove from the string
    :return: reformatted string
    """
    null = dict(zip([ord(i) for i in reform], [None] * len(reform)))
    return string.translate(null)


def get_size(nested: list) -> int:
    """
    Function to get number of elements in a nested list

    Takes:
    ------
    * `nested`: list - the nested list to count the elements of
    """
    count = 0
    for elem in nested:
        if isinstance(elem, list):
            count += get_size(elem)
        else:
            count += 1
    return count


def flatten(nested: list) -> list:
    """
    Function to flatten a nested list or a list of lists to a 1D list.

    Takes:
    ------
    * `nested`: list - the nested list or list of list to flatten
    """
    flat = []
    for elem in nested:
        if isinstance(elem, list):
            flat.extend(flatten(elem))
        else:
            flat.append(elem)
    return flat


def intersect(surfaces: list) -> list:
    """
    Function to handle operators for surface intersection: operator is ''

    Takes:
    ------
    * `surfaces`: list - contains Surface object instances to intersect
    """
    flat = flatten(surfaces)
    out = [flat[0].copy()]
    for s in flat[1:]:
        new = s.copy()
        new._operator = '' + new._operator
        out.append(new)
    return out


def unite(surfaces: list) -> list:
    """
    Function to handle operators for surface union: operator is ':'

    Takes:
    ------
    * `surfaces`: list - contains Surface object instances to unite
    """
    flat = flatten(surfaces)
    out = [flat[0].copy()]
    for s in flat[1:]:
        new = s.copy()
        new._operator = ' : ' + new._operator
        out.append(new)
    return out


def surface_complement(surfaces: list) -> list:
    """
    Function to handle surface complement (flipping): operator is '-'

    Takes:
    ------
    * `surfaces`: list - contains Surface object instances to complement
    """
    flat = flatten(surfaces)
    out = []
    for s in flat:
        new = s.copy()
        new.flip()
        out.append(new)
    return out


# def cell_complement(surfaces: list) -> list:  # is it correct that this is a list of surfaces or is it a list of cells? --- it takes a Cell Pablo
#     """
#     Function to handle operators for cell complement: operator is '#'
#
#     Takes:
#     ------
#     * `surfaces` is a list of Surface object instances to complement
#     """
#     flat = flatten(surfaces)
#     out = [flat[0].copy()]
#     for s in flat[1:]:
#         new = s.copy()
#         new.operator = ' # ' + new.operator
#         out.append(new)
#     return out

@dataclass(slots=True)
class Comment:
    """
    Handles:
    --------
    Base class to handle the comments to the Serpent 2 input

    Methods:
    --------
    * `write()`: internal method to write on a file with proper formatting for Serpent 2 input time interval definition
    * `copy()`: copies the object instance to another memory allocation

    Takes:
    ------
    * `txt`: string - the comment text
    """
    txt: str

    def __str__(self):
        return f'/* {self.txt} */\n' if self.txt != '' else ''

    def __add__(self, comment2):
        return Comment(self.txt + '\n   ' + comment2.txt)

    def write(self, file: str):
        """
        Method to write on a file with proper formatting for Serpent 2 input comment, inserted in '/*' '*/'

        Takes:
        ------
        * `file`: string - is the name of the file where to write
        """
        with open(file=file, mode='a') as f:
            f.write(self.__str__())

    def copy(self):
        """
        Copies the object instance to another memory allocation.

        Returns:
        --------
        Returns a variable pointing to a new memory allocation
        """
        return cp.deepcopy(self)


@dataclass(slots=True)
class StandaloneComment(Comment):
    """
    Handles:
    --------
    Handles the comments to the Serpent 2 input standing alone

    Inherits from:
    --------------
    Comment

    Required inherited parameters:
    ------------------------------
    * `txt`: string - the comment text
    """

    def __str__(self):
        return f'/* {self.txt} */\n\n' if self.txt != '' else ''


@dataclass(slots=True)
class InlineComment(Comment):
    """
    Handles:
    --------
    Handles the in-line comments to the Serpent 2 input

    Inherits from:
    --------------
    Comment

    Required inherited parameters:
    ------------------------------
    * `txt`: string - the comment text
    """

    def __str__(self):
        txt = ''.join(self.txt.split('\n'))
        return f'  % {txt}\n' if txt != '' else '\n'


@dataclass(slots=True)
class Entity:
    """
    Handles:
    --------
    General attributes of Serpent 2 input entities: provides them with identity and comment

    Takes:
    ------
    * `name`: string or integer - is the identity of the Serpent 2 entity
    * `comment`: Comment object instance - is the comment to the Serpent entity. Default is empty comment
    * `inline_comment`: InlineComment object instance - is the comment to be written on the same line: ment for short
                        description of the specific instance. Default is empty comment
    """
    name: any
    comment: Comment = Comment('')
    inline_comment: InlineComment = InlineComment('')

    def __str__(self):
        string = self.comment.__str__() + f"""{self.name}"""
        return string

    def assess(self):
        print(id(self))


@dataclass(slots=True)
class Other:
    """
    Handles:
    --------
    Handles whatever Serpent 2 feature not implemented in this code yet

    Methods:
    --------
    * `write()`: internal method to write on a file with same formatting as in the instance definition
    * `copy()`: copies the object instance to another memory allocation

    Takes:
    ------
    * `string`: string - Serpent 2 properly formatted string for the desired implementation
    """
    string: str

    def __str__(self):
        return self.string

    def write(self, file: str):
        """
        Internal method to write on a file with same formatting as in the instance definition

        Takes:
        ------
        * `file`: string - is the name of the file where to write
        """
        with open(file, 'a') as f:
            f.write(self.string)

    def copy(self):
        """
        Copies the object instance to another memory allocation.

        Returns:
        --------
        Returns a variable pointing to a new memory allocation
        """
        return cp.deepcopy(self)
