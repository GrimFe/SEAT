import os

# import pandas as pd
#
# import sandy


## This module was inspired from sandy API by Luca Fiorito
__all__ = [
    "ELEMENTS",
    "METASTATES",
    "expand_za",
    "expand_zam",
    "za2latex",
    "zam2latex",
    "zam2za",
    "zam2nuclide",
    "za2zam",
    "nuclide2zam",
    "nuclide2za",
    "nuclide2element"
]

# pd.options.display.float_format = '{:.5e}'.format


ELEMENTS = {
    1: 'H',
    2: 'He',
    3: 'Li',
    4: 'Be',
    5: 'B',
    6: 'C',
    7: 'N',
    8: 'O',
    9: 'F',
    10: 'Ne',
    11: 'Na',
    12: 'Mg',
    13: 'Al',
    14: 'Si',
    15: 'P',
    16: 'S',
    17: 'Cl',
    18: 'Ar',
    19: 'K',
    20: 'Ca',
    21: 'Sc',
    22: 'Ti',
    23: 'V',
    24: 'Cr',
    25: 'Mn',
    26: 'Fe',
    27: 'Co',
    28: 'Ni',
    29: 'Cu',
    30: 'Zn',
    31: 'Ga',
    32: 'Ge',
    33: 'As',
    34: 'Se',
    35: 'Br',
    36: 'Kr',
    37: 'Rb',
    38: 'Sr',
    39: 'Y',
    40: 'Zr',
    41: 'Nb',
    42: 'Mo',
    43: 'Tc',
    44: 'Ru',
    45: 'Rh',
    46: 'Pd',
    47: 'Ag',
    48: 'Cd',
    49: 'In',
    50: 'Sn',
    51: 'Sb',
    52: 'Te',
    53: 'I',
    54: 'Xe',
    55: 'Cs',
    56: 'Ba',
    57: 'La',
    58: 'Ce',
    59: 'Pr',
    60: 'Nd',
    61: 'Pm',
    62: 'Sm',
    63: 'Eu',
    64: 'Gd',
    65: 'Tb',
    66: 'Dy',
    67: 'Ho',
    68: 'Er',
    69: 'Tm',
    70: 'Yb',
    71: 'Lu',
    72: 'Hf',
    73: 'Ta',
    74: 'W',
    75: 'Re',
    76: 'Os',
    77: 'Ir',
    78: 'Pt',
    79: 'Au',
    80: 'Hg',
    81: 'Tl',
    82: 'Pb',
    83: 'Bi',
    84: 'Po',
    85: 'At',
    86: 'Rn',
    87: 'Fr',
    88: 'Ra',
    89: 'Ac',
    90: 'Th',
    91: 'Pa',
    92: 'U',
    93: 'Np',
    94: 'Pu',
    95: 'Am',
    96: 'Cm',
    97: 'Bk',
    98: 'Cf',
    99: 'Es',
    100: 'Fm',
    101: 'Md',
    102: 'No',
    103: 'Lr',
    104: 'Rf',
    105: 'Db',
    106: 'Sg',
    107: 'Bh',
    108: 'Hs',
    109: 'Mt',
    110: 'Ds',
    111: 'Rg',
    112: 'Uub',
    113: 'Uut',
    114: 'Uuq',
    115: 'Uup',
    116: 'Uuh',
    117: 'Uus',
    118: 'UUp',
}

ATOMIC_NUMBERS = {v: k for k, v in ELEMENTS.items()}

METASTATES = {
    0: "g",
    1: "m",
    2: "n",
    3: "o",
}

METASTATES_FLIP = {v: k for k, v in METASTATES.items()}


def expand_za(za: int, method: str = "nndc", meta: int = 0):
    """
    Expands a ZA number to Z, A, M.

    Takes:
    ------
    * `za`: integer - the ZA to expand
    * `method`: string - the method to use for including the metastate. Default is `'nndc'`.
    * `meta`: integer - the metastate included

    Returns:
    --------
    * `z`: integer - atomic number of the nuclide.
    * `a`: integer - mass number of the nuclide.
    * `m`: integer - metastate number of the nuclide.
    """

    z = int(za // 1000)
    a = int(za - z * 1000)
    if method == "nndc":
        m = 0
        if a >= 300:
            m = 1
            a = a - 300 - m * 100
    else:
        m = int(meta)
    return z, a, m


def expand_zam(zam):
    """
    Transforms the ZAM to its components: Z, A, M.

    Takes:
    ------
    * `zam`: integer - ZAM number of the nuclide.

    Returns:
    ------
    * `z`: integer - atomic number of the nuclide.
    * `a`: integer - mass number of the nuclide.
    * `m`: integer - metastate number of the nuclide.
    """
    z = int(zam // 10000)
    a = int(zam - z * 10000) // 10
    m = int(zam - z * 10000 - a * 10)
    return z, a, m


def get_zam(z, a, m):
    """
    Composes the ZAM from its components: Z, A, M

    Takes:
    ------
    * `z`: integer - atomic number of the nuclide.
    * `a`: integer - mass number of the nuclide.
    * `m`: integer - metastate number of the nuclide.

    Returns:
    ------
    * `zam`: integer - ZAM number of the nuclide.
    """
    zam = z * 10000 + a * 10 + m
    return int(zam)


def get_za(z, a, m, method="nndc"):
    """
    Writes a nuclide ZA given its ZAM

    Takes:
    ------
    * `z`: integer - atomic number of the nuclide.
    * `a`: integer - mass number of the nuclide.
    * `m`: integer - metastate number of the nuclide.
    * `method`: string - method to use for the conversion, default is `'nndc'`.

    Returns:
    --------
    * `za`: integer - ZA of the nuclide.
    * `m`: integer - metastate of the nuclide.
    """
    if m != 0 and method == "nndc":
        za = z * 1000 + a + 300 + m * 100
    else:
        za = z * 1000 + a
    return int(za), m


def za2zam(za: int, method: str = "nndc", meta: int = 0):
    """
    Transforms the ZA to the corresponding ZAM.

    Takes:
    ------
    * `za`: integer - ZA of the nuclide
    * `method`: string - method to use in the ZA decomposition. Default is `'nndc'`
    * `meta`: integer - metastate number of the nuclide.

    Returns:
    ------
    * integer - ZAM number of the nuclide.
    """
    return get_zam(*expand_za(za, method=method, meta=meta))


def zam2za(zam: int, method: str = "nndc"):
    """
    Transforms the ZAM to the corresponding ZA.

    Takes:
    ------
    * `zam`: integer - ZAM of the nuclide
    * `method`: string - method to use in the ZA compositin. Default is `'nndc'`

    Returns:
    ------
    * integer - ZA number of the nuclide.
    """
    z, a, m = expand_zam(zam)
    return get_za(z, a, m, method=method)

def to_latex(z: int, a: int, m: int) -> str:
    A = f"{a}"
    m_ = get_meta_letter(m, skip_ground=True)
    M = "\text{" + m_ + "}" if m_ != '' else ''
    Z = ELEMENTS[z]
    return "$^{" + A + M + "}$" + Z

def za2latex(za):
    z, a, m = expand_za(za)
    return to_latex(z, a, m)


def zam2latex(zam):
    z, a, m = expand_zam(zam)
    return to_latex(z, a, m)


def zam2nuclide(zam, atomic_number=False, sep=""):
    """
    Convert ZAM to string such with symbol and mass, such as `922350` to
    `"U235"` or `952421` to `"Am242m"`.
    Parameters
    ----------
    zam : `int`
        nuclide ZAM indicator
    atomic_number : `bool`, optional, default is `False`
        flag to include the atomic number in the nuclide name
    sep : `str`, optional, default is `''`
        separation character(s) to place between the atomic number
        (if present), the element ID, and the mass number.
    Returns
    -------
    `string`
        nuclide expressed with symbol and mass
    Examples
    --------
    >>> zam2nuclide(922350)
    'U235'
    >>> zam2nuclide(922350, atomic_number=True)
    '92U235'
    >>> zam2nuclide(922350, atomic_number=True, sep="-")
    '92-U-235'
    >>> zam2nuclide(922350, atomic_number=False, sep="-")
    'U-235'
    >>> zam2nuclide(952420)
    'Am242'
    >>> zam2nuclide(952421)
    'Am242m'
    >>> zam2nuclide(952421, atomic_number=True, sep="_")
    '95_Am_242m'
    >>> zam2nuclide(952422)
    'Am242n'
    """
    z, a, m = expand_zam(zam)
    sym = ELEMENTS[z]
    meta = get_meta_letter(m, skip_ground=True)
    out = f"{sym}{sep}{a}{meta}"
    if atomic_number:
        out = f"{z}{sep}{out}"
    return out


def nuclide2zam(nuclide, atomic_number=False, sep=""):
    """
    Convert string with symbol and mass number to ZAM, such as `"U235"` to
    `922350` or `"Am242m"` to `952421`.
    Parameters
    ----------
    nuclide : `str`
        nuclide expressed with symbol and mass.
    atomic_number : `bool`, optional, default is `False`
        flag to pass a string with the atomic number in `nuclide`
    sep : `str`, optional, default is `''`
        separation character(s) placed between the atomic number
        (if present), the element ID, and the mass number.
    Returns
    -------
    zam : `int`
        nuclide ZAM indicator
    Examples
    --------
    >>> nuclide2zam('U235')
    922350
    >>> nuclide2zam('92U235', atomic_number=True)
    922350
    >>> nuclide2zam('92-U-235', atomic_number=True, sep="-")
    922350
    >>> nuclide2zam('Am242')
    952420
    >>> nuclide2zam('Am242m')
    952421
    >>> nuclide2zam('95_Am_242m', atomic_number=True, sep="_")
    952421
    >>> nuclide2zam('Am242n')
    952422
    """
    sym = ""
    num = ""
    if not nuclide[-1].isalpha():
        nuclide_ = nuclide + "g"
    else:
        nuclide_ = nuclide
    if atomic_number and nuclide_[:2].isnumeric():
        nuclide_ = nuclide_[2:]
    elif atomic_number and nuclide_[:1].isnumeric():
        nuclide_ = nuclide_[1:]
    if sep != "":
        nuclide_.replace(sep, "")
    for i in nuclide_[:-1]:
        if i.isalpha():
            sym += i
        if i.isnumeric():
            num += i
    z = ATOMIC_NUMBERS[sym]
    a = int(num)
    m = METASTATES_FLIP[nuclide_[-1]]
    out = get_zam(z, a, m)
    return out


def nuclide2za(nuclide):
    """
    Convert ZA to string with symbol and mass, such as 92235 to
    "U235" or 95642 to "Am242m".

    Takes:
    ------
    `za` : int
        nuclide ZA indicator
    `method` : str, optional, default is 'nndc'
        method of the representation of the metastable state in the ZA number
    `meta` : int, optional, default is 0
        metastable state of the nuclide
    `atomic_number` : bool, optional, default is False
        flag to include the atomic number in the nuclide name
    `sep` : str, optional, default is ''
        separation character(s) to place between the atomic number
        (if present), the element ID, and the mass number.

    Returns:
    --------
    `string`
        nuclide expressed with symbol and mass number

    """
    zam = nuclide2zam(nuclide)
    return zam2za(zam)


def get_meta_letter(m, skip_ground=False):
    meta = METASTATES[m]
    if skip_ground and m == 0:
        meta = ""
    return meta


def nuclide2element(nuclide: str, *args, **kwargs):
    return ELEMENTS[expand_zam(nuclide2zam(nuclide, *args, **kwargs))[0]]