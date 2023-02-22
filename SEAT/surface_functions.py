# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 14:11:04 2023

@author: fgrimald
"""
# from SEAT.Serpent2InputWriter import reformat

def _p_params(offset: float) -> str:
    """
    Formats the parameters for a plane surface.

    Parameters
    ----------
    offset : float
        the plane offset [cm].

    Returns
    -------
    str
        the ordered parameters as Serpent2 formatted string.

    """
    return str(offset)

def px_params(x: float) -> str:
    """
    Formats the parameters for a X-parallel plane surface.

    Parameters
    ----------
    x : float
        the plane offset [cm].

    Returns
    -------
    str
        the ordered parameters as Serpent2 formatted string.

    """
    return _p_params(x)

def py_params(y: float) -> str:
    """
    Formats the parameters for a Y-parallel plane surface.

    Parameters
    ----------
    y : float
        the plane offset [cm].

    Returns
    -------
    str
        the ordered parameters as Serpent2 formatted string.

    """
    return _p_params(y)

def pz_params(z: float) -> str:
    """
    Formats the parameters for a Z-parallel plane surface.

    Parameters
    ----------
    z : float
        the plane offset [cm].

    Returns
    -------
    str
        the ordered parameters as Serpent2 formatted string.

    """
    return _p_params(z)

def plane_params(A: float=None, B: float=None, C: float=None, D: float=None,
                 c1: tuple[float]=None, c2: tuple[float]=None, c3: tuple[float]=None) -> str:
    """
    Formats the parameters for a custom plane surface.

    Parameters
    ----------
    A : float, optional
        A in equation Ax + By + Cz - D. The default is None.
    B : float, optional
        B in equation Ax + By + Cz - D. The default is None.
    C : float, optional
        C in equation Ax + By + Cz - D. The default is None.
    D : float, optional
        D in equation Ax + By + Cz - D. The default is None.
    c1 : tuple[float], optional
        the coordinates [cm] of the first point of the plane. The default is
        None.
    c2 : tuple[float], optional
        the coordinates [cm] of the second point of the plane. The default is
        None.
    c3 : tuple[float], optional
        the coordinates [cm] of the third point of the plane. The default is
        None.

    Returns
    -------
    str
        the ordered parameters as Serpent2 formatted string.

    """
    out = ''
    if c1 is None and c2 is None and c3 is None:
        # equation-defined plane
        out = f"{A} {B} {C} {D}"
    elif A is None and B is None and C is None and D is None:
        # point-defined plane
        if out != '': raise Exception("Excessive parameters passed.")
        out = f"{c1[0]} {c1[1]} {c1[2]} {c2[0]} {c2[1]} {c2[2]} {c3[0]} {c3[1]} {c3[2]}"
    else:
        msg = "Insufficient parameters passed."
        raise Exception(msg)
    return out

def _cylinder_params(r: float, i0: float=0, j0: float=0, k0: float=None,
                     k1: float=None) -> str:
    """
    Formats the parameters to cylindrical square surfaces

    Parameters
    ----------
    r : float
        the cylinder radus [cm].
    i0 : float, optional
        the centrar coordinate in the first direction. The default is 0.
    j0 : float, optional
        the centrar coordinate in the second direction. The default is 0.
    k0 : float, optional
        the lower bound of the cylinder [cm]. The default is None.
    k1 : float, optional
        the upper bound of the cylinder [cm]. The default is None.

    Returns
    -------
    str
        the ordered parameters as Serpent2 formatted string.

    """
    return " ". join([str(i) for i in [i0, j0, r, k0, k1] if i is not None])

def sqc_params(r: float=0, x0: float=0, y0: float=0, z0: float=None,
               z1: float=None, s: float=None) -> str:
    """
    Formats the parameters to cylindrical square surfaces oriented parallel to
    the Z-direction.

    Parameters
    ----------
    r : float
        the cylinder radus [cm].
    x0 : float, optional
        the X-coordinate of the center [cm]. The default is 0.
    y0 : float
        the Y-coordinate of the center [cm]. The default is 0.
    z0 : float, optional
        the lower bound of the cylinder [cm]. The default is None.
    z1 : float, optional
        the upper bound of the cylinder [cm]. The default is None.
    s : float, optional
        the radius of rounded corners [cm]. The default is None.

    Returns
    -------
    str
        the ordered parameters as Serpent2 formatted string.

    """
    s_ = f" {s}" if s is not None else ''
    return _cylinder_params(r, x0, y0, z0, z1) + s_

def cyl_params(r: float=0, x0: float=0, y0: float=0, z0: float=None,
               z1: float=None) -> str:
    """
    Formats the parameters to cylindrical surfaces oriented parallel to the
    Z-direction.

    Parameters
    ----------
    r : float
        the cylinder radus [cm].
    x0 : float, optional
        the X-coordinate of the center [cm]. The default is 0.
    y0 : float
        the Y-coordinate of the center [cm]. The default is 0.
    z0 : float, optional
        the lower bound of the cylinder [cm]. The default is None.
    z1 : float, optional
        the upper bound of the cylinder [cm]. The default is None.

    Returns
    -------
    str
        the ordered parameters as Serpent2 formatted string.

    """
    return _cylinder_params(r, x0, y0, z0, z1)

def cylx_params(r: float=0, y0: float=0, z0: float=0, x0: float=None,
               x1: float=None) -> str:
    """
    Formats the parameters to cylindrical surfaces oriented parallel to the
    X-direction.

    Parameters
    ----------
    r : float
        the cylinder radus [cm].
    y0 : float, optional
        the X-coordinate of the center [cm]. The default is 0.
    z0 : float
        the Y-coordinate of the center [cm]. The default is 0.
    x0 : float, optional
        the lower bound of the cylinder [cm]. The default is None.
    x1 : float, optional
        the upper bound of the cylinder [cm]. The default is None.

    Returns
    -------
    str
        the ordered parameters as Serpent2 formatted string.

    """
    return _cylinder_params(r, y0, z0, x0, x1)

def cyly_params(r: float=0, x0: float=0, z0: float=0, y0: float=None,
               y1: float=None) -> str:
    """
    Formats the parameters to cylindrical surfaces oriented parallel to the
    Y-direction.

    Parameters
    ----------
    r : float
        the cylinder radus [cm].
    x0 : float, optional
        the X-coordinate of the center [cm]. The default is 0.
    z0 : float
        the Y-coordinate of the center [cm]. The default is 0.
    y0 : float, optional
        the lower bound of the cylinder [cm]. The default is None.
    y1 : float, optional
        the upper bound of the cylinder [cm]. The default is None.

    Returns
    -------
    str
        the ordered parameters as Serpent2 formatted string.

    """
    return _cylinder_params(r, x0, z0, y0, y1)

def cylz_params(r: float=0, x0: float=0, y0: float=0, z0: float=None,
               z1: float=None) -> str:
    """
    Formats the parameters to cylindrical surfaces oriented parallel to the
    Z-direction.

    Parameters
    ----------
    r : float
        the cylinder radus [cm].
    x0 : float, optional
        the X-coordinate of the center [cm]. The default is 0.
    y0 : float
        the Y-coordinate of the center [cm]. The default is 0.
    z0 : float, optional
        the lower bound of the cylinder [cm]. The default is None.
    z1 : float, optional
        the upper bound of the cylinder [cm]. The default is None.

    Returns
    -------
    str
        the ordered parameters as Serpent2 formatted string.

    """
    return _cylinder_params(r, x0, y0, z0, z1)
