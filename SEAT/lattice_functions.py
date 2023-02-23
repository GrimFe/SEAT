# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 10:17:34 2023

@author: fgrimald
"""

def _NxNy_params(x0: float=0, y0: float=0, Nx: int=None, Ny: int=None,
                 pitch: float=None) -> str:
    """
    Formats the parameters for a lattice with Nx elements in the X direction
    and Ny in the Y direction.

    Parameters
    ----------
    x0 : float, optional
        X coordinate of the lattice center. The default is 0.
    y0 : float, optional
        Y coordinate of the lattice center. The default is 0.
    Nx : int, optional
        Number of lattice elements in the X direction. The default is None.
    Ny : int, optional
        Number of lattice elements in the Y direction. The default is None.
    pitch : float, optional
        center-to-center distance of the lattice elements. The default is None.

    Returns
    -------
    str
        the ordered parameters as Serpent2 formatted string.

    """
    return " ". join([str(i) for i in [x0, y0, Nx, Ny, pitch] if i is not None])

def square_params(x0: float=0, y0: float=0, Nx: int=None, Ny: int=None,
                 pitch: float=None) -> str:
    """
    Formats the parameters for a square lattice with Nx elements in the X
    direction and Ny in the Y direction explicitly.

    Parameters
    ----------
    x0 : float, optional
        X coordinate of the lattice center. The default is 0.
    y0 : float, optional
        Y coordinate of the lattice center. The default is 0.
    Nx : int, optional
        Number of lattice elements in the X direction. The default is None.
    Ny : int, optional
        Number of lattice elements in the Y direction. The default is None.
    pitch : float, optional
        center-to-center distance of the lattice elements. The default is None.

    Returns
    -------
    str
        the ordered parameters as Serpent2 formatted string.

    """
    return _NxNy_params(x0, y0, Nx, Ny, pitch)