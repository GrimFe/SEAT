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
    Removes from one string all characters present in another string.

    Parameters
    ----------
    string : str
        the string to remove characters from.
    reform : str
        the string of the characters to remove.

    Returns
    -------
    str
        reformatted string: the `string` without the characters in `reform`.

    """
    null = dict(zip([ord(i) for i in reform], [None] * len(reform)))
    return string.translate(null)


def get_size(nested: list) -> int:
    """
    Comutes the number of elements in a nested list.

    Parameters
    ----------
    nested : list
        the nested list to count the elements of.

    Returns
    -------
    int
        the number of elements in the nested list.

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
    Flattens a nested list to a 1D-iterable list.

    Parameters
    ----------
    nested : list
        the nested list or list of lists to flatten.

    Returns
    -------
    list
        the flattened 1D-iterable list.

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
    Handles the surface intersection operator ('' in Serpent 2).

    Parameters
    ----------
    surfaces : list
        the `SEAT.Surface` instances to intersect.

    Returns
    -------
    list
        the `SEAT.Surface` instances with operator modified to be also
        intersection.

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
    Handles the surface union operator (':' in Serpent 2).

    Parameters
    ----------
    surfaces : list
        the `SEAT.Surface` instances to intersect.

    Returns
    -------
    list
        the `SEAT.Surface` instances with operator modified to be also union.

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
    Handles the surface surface complement operator ('-' in Serpent 2).

    Parameters
    ----------
    surfaces : list
        the `SEAT.Surface` instances to complement.

    Returns
    -------
    list
        the `SEAT.Surface` instances with operator modified to be also
        complement.

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
    Base class to handle the comments to the Serpent 2 input.

    Attributes
    ----------
    txt : str
        the text of the comment.

    Methods
    -------
    write :
        writes the `SEAT.Comment` instance to a file.
    copy :
        copies the object instance to another memory allocation.

    """
    txt: str

    def __str__(self):
        return f'/* {self.txt} */\n' if self.txt != '' else ''

    def __add__(self, comment2):
        return Comment(self.txt + '\n   ' + comment2.txt)

    def write(self, file: str, mode: str='a') -> None:
        """
        Writes the `SEAT.Comment` instance to a file according to the Serpent 2
        syntax ('/* {comment} */').

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
        with open(file=file, mode=mode) as f:
            f.write(self.__str__())

    def copy(self):
        """
        Copies the object instance to another memory allocation.

        Returns
        -------
        SEAT.Comment
            a copy of the copied `SEAT.Comment` instance.

        """
        return cp.deepcopy(self)


@dataclass(slots=True)
class StandaloneComment(Comment):
    """
    Class for paragraphs of comments standing alone with respect to the text.
    Inherits from `SEAT.Comment`.

    Attributes
    ----------
    txt : str
        the text of the comment.

    """
    def __str__(self):
        return f'/* {self.txt} */\n\n' if self.txt != '' else ''


@dataclass(slots=True)
class InlineComment(Comment):
    """
    Class for inline comments standing by other text.
    Inherits from `SEAT.Comment`.

    Attributes
    ----------
    txt : str
        the text of the comment.

    """
    def __str__(self):
        txt = ''.join(self.txt.split('\n'))
        return f'  % {txt}\n' if txt != '' else '\n'


@dataclass(slots=True)
class Entity:
    """
    Base class to handle the comments to the Serpent 2 input. Provides them
    with an identity, comments and a write method.

    Attributes
    ----------
    name : str | int
        the identity of the Serpent 2 entity.
    comment : `SEAT.Comment`, optional
        the comment to the Serpent entity. The default is SEAT.Comment('').
    inline_comment : `SEAT.InlineComment`, optional
        the comment to be written on the same line as the Serpent 2 entity id.
        The default is SEAT.Comment('').

    Methods
    --------
    assess :
        prints the `SEAT.Entity` python id.
    write :
        writes the `SEAT.Entity` to a file.

    """
    name: str | int
    comment: Comment = Comment('')
    inline_comment: InlineComment = InlineComment('')

    def __str__(self):
        string = self.comment.__str__() + f"""{self.name}"""
        return string

    def assess(self):
        """
        Prints the `SEAT.Entity` python id.

        Returns
        -------
        None.

        """
        print(id(self))

    def write(self, file: str, mode: str='a'):
        """
        Writes the `SEAT.Entity.__str__()` to a file.

        Parameters
        ----------
        file : str
            the name of the file where to write.
        mode : str, optional
            mode to open the file. The default is 'a'.

        Returns
        -------
        None.

        """
        with open(file, mode=mode) as f:
            f.write(self.__str__())


@dataclass(slots=True)
class Other:
    """
    Object for not-yet-implemented Serpent 2 cards and entities.

    Attributes
    ----------
    string : str
        the Serpent 2 properly formatted string for the desired implementation.

    Methods
    -------
    write :
        writes the `SEAT.Other` instance to a file.
    copy :
        copies the object instance to another memory allocation.

    """
    string: str

    def __str__(self):
        return self.string

    def write(self, file: str, mode: str='a'):
        """
        Writes the `SEAT.Other` instance to a file.

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
            f.write(self.string)

    def copy(self):
        """
        Copies the object instance to another memory allocation.

        Returns:
        --------
        SEAT.Other
            a copy of the `SEAT.Other` instance.

        """
        return cp.deepcopy(self)
