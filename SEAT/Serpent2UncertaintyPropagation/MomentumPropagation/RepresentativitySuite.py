import numpy as np
import pandas as pd

import SEAT.Serpent2UncertaintyPropagation.MomentumPropagation.SensitivitySuite as SS
import SEAT.Serpent2UncertaintyPropagation.MomentumPropagation.CovarianceSuite as CS

from functools import wraps
from uncertainties import unumpy

__all__ = [
    "SCS"
    ]

def sandwich_unsafe(s1: pd.DataFrame, cov: pd.DataFrame, s2: pd.DataFrame) -> pd.DataFrame:
    """
    Computes the sandwich rule for the covovariancove propagation.
    Computes s1' cov s2.

    Parameters
    ----------
    s1: pd.DataFrame
        the first MxN sensitivity vector
    cov: pd.DataFrame
        the MxM covariancove matrix
    s2: pd.DataFrame
        the second MxN sensitivity vector

    Returns
    -------
    pd.DataFrame of the new covariancove matrix

    """
    return s1.T @ cov @ s2

def sandwich(s1: pd.DataFrame, cov: pd.DataFrame, s2: pd.DataFrame) -> pd.DataFrame:
    """
    Computes the sandwich rule for the covariancove propagation.
    Computes s1' cov s2.

    Parameters
    ----------
    s1: pd.DataFrame
        the first MxN sensitivity vector
    cov: pd.DataFrame
        the MxM covariancove matrix
    s2: pd.DataFrame
        the second MxN sensitivity vector

    Returns
    -------
    pd.DataFrame of the new covariancove matrix

    """
    idx = cov.index.intersection(s1.index).intersection(s2.index)
    return sandwich_unsafe(s1.loc[idx], cov.loc[idx, idx], s2.loc[idx])

def stdev(s1: pd.DataFrame, cov: pd.DataFrame, check_idx: bool=True) -> pd.DataFrame:
    """
    Computes the standard deviation propagated from the covariance
    and a senstitivity vector.

    Parameters:
    s1: pd.DataFrame
        the sensitivity vector
    cov: pd.DataFrame
        the covariancove matrix
    check_idx: bool
        performs the operation over intersecting indices only. Default is True.

    Returns
    -------
    pd.DataFrame of the computed standard deviation

    Note:
    -----
    This is implemented with **0.5 rather than with np.sqrt() as that function
    only takes float arguments and now we work with ufloats.

    """
    if check_idx:
        return sandwich(s1, cov, s1)**0.5
    else:
        return sandwich_unsafe(s1, cov, s1)**0.5

def representativity(s1: pd.DataFrame, cov: pd.DataFrame, s2: pd.DataFrame,
                     check_idx: bool=True) -> pd.DataFrame:
    """
    Computes the representativity of one model (s1) to another (s2),
    which is the correlation of the two models induced by the
    covariance matrix cov.
    Computes s1' cov s2 / (sqrt(s1' cov s1) * sqrt(s2' cov s2)).

    Parameters
    ----------
    s1: pd.DataFrame
        the first MxN sensitivity vector
    cov: pd.DataFrame
        the MxM covariancove matrix
    s2: pd.DataFrame
        the second MxN sensitivity vector
    check_idx: bool
        performs the operation over intersecting indices only. Default is True.

    Returns
    -------
    pd.DataFrame of the representativity

    """
    n = sandwich(s1, cov, s2) if check_idx else sandwich_unsafe(s1, cov, s2)
    d = stdev(s1, cov, check_idx) * stdev(s2, cov, check_idx)
    return n / d

def drop_uncertainty_calc(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        user_uncertainty = self.uncertain  ## store this for later
        self.uncertain = False
        result = func(self, *args, **kwargs)
        self.uncertain = user_uncertainty  ## resetting the uncertainty
        return result
    return wrapper


class SCS:
    """
    Container of sensitivity and covariance information.

    Attributes
    ----------
    s1: `SEAT.Sensitivity`
        The first sensitivity vector queried for one observable only.
        Use `SEAT.Sensitivity.observe()` for the purpose.
    s2: `SEAT.Sensitivity`
        The second sensitivity vector queried for one observable only.
        Use `SEAT.Sensitivity.observe()` for the purpose.
    cov: `SEAT.Covariance`
        The covariance matrix to consider.
    uncertain: bool
        flag to track the uncertainties of the computed quantities.
        The default is True.

    Properties
    ----------
    idx: `pandas.MiltiIndex`
        union of the s1 and s2 indices intersected with the one of cov.
        Nuclides and reactions considered in the calculations.
    s1_:
        s1 reindexed and filled with zeros where needed.
    s2_:
        s2 reindexed and filled with zeros where needed.
    cov_:
        cov reindexed.
    sandwich:
        computes the bare sandwich formula using s1_, cov_ and s2_.
    representativity:
        computes the representativity for s1_, cov_ and s2_.
    reaction_wise:

    reaction_wise_rep:

    nuclide_wise:

    nuclide_wise_rep:

    Methods
    --------
    _get_sharing:

    _zero_filled_s:

    stdev:
        computes the standard deviation using s1 or s2 and cov.
        All the indices in the original data are considered.
    stdev_no_unc:
        computes the standard deviation using s1 or s2 and cov.
        All the indices in the original data are considered.
        Does not propagate the uncertanty to the output.
    reaction_apportion:

    nuclide_apportion:

    Note
    ----
    

    """
    def __init__(self, s1: SS.Sensitivity, cov: CS.Covariance,
                 s2: SS.Sensitivity, uncertain: bool=True):
        self.s1 = s1
        self.s2 = s2
        self.cov = cov
        self.uncertain = uncertain

    def _get_sharing(self, mat, mt=None):
        query = "N == @mat" if mt is None else "N == @mat and MT == @mt"
        loc = mat if mt is None else (mat, mt)
        s1 = self.s1_.reset_index().query(query).set_index(["N", "MT", "E [MeV]"]).loc[loc]
        s2 = self.s2_.reset_index().query(query).set_index(["N", "MT", "E [MeV]"]).loc[loc]
        cov = self.cov_.loc[loc, loc]
        return s1, cov, s2

    @property
    def idx(self) -> pd.MultiIndex:
        idx1 = self.s1.idx.droplevel("Observable")
        idx2 = self.s2.idx.droplevel("Observable")
        s_idx = idx1.union(idx2)  ## accounting for all the nuclides/reactions in the system
        return self.cov.index.intersection(s_idx)  ## if something is not in cov is out of the equation

    @property
    def s1_(self) -> pd.Series:
        s = self.s1.data.droplevel("Observable")
        return self._zero_filled_s(s)

    @property
    def s2_(self) -> pd.Series:
        s = self.s2.data.droplevel("Observable")
        return self._zero_filled_s(s)
        

    @property
    def cov_(self) -> pd.DataFrame:
        return self.cov.loc[self.idx, self.idx]

    def _zero_filled_s(self, s) -> pd.Series:
        old_idx = s.index
        dif_idx = self.idx.difference(old_idx, sort=False)
        zero = np.zeros(len(dif_idx))
        empty = pd.Series(unumpy.uarray(zero, zero),
                          index=dif_idx)
        if self.uncertain:
            out = pd.concat([s.S, empty])
            out.name = "S"
        else:
            out = pd.concat([s.nom_val, empty])
            out.name = "nom_val"
        return out.loc[self.idx]

    @property
    def sandwich(self) -> pd.DataFrame:
        ## indices are already handled in the x_ variables
        return sandwich_unsafe(self.s1_, self.cov_, self.s2_)

    @property
    def representativity(self) -> pd.DataFrame:
        ## indices are already handled in the x_ variables
        return representativity(self.s1_, self.cov_, self.s2_, check_idx=False)

    @property
    @drop_uncertainty_calc
    def nuclide_wise(self) -> pd.Series:
        """
        Computes the sandwich rule for each nuclide without propagating the
        uncertainty to it.

        """
        out = {}
        for mat, df in self.cov_.groupby(level=0):
            s1, cov, s2 = self._get_sharing(mat)
            ## indices are already handled in the x_ variables
            portion = sandwich_unsafe(s1, cov, s2).nom_val.values[0]  # portion must be scalar
            out[mat] = portion
        return pd.Series(out)

    @property
    @drop_uncertainty_calc
    def nuclide_wise_rep(self) -> pd.Series:
        """
        Computes the representativity for each nuclide.

        """
        out = {}
        for mat, df in self.cov_.groupby(level=0):
            s1, cov, s2 = self._get_sharing(mat)
            ## indices are already handled in the x_ variables
            portion = representativity(s1, cov, s2, check_idx=False).nom_val.values[0]  # portion must be scalar
            out[mat] = portion
        return pd.Series(out)

    @property
    @drop_uncertainty_calc
    def reaction_wise(self) -> pd.Series:
        """
        Computes the sandwich rule for each nuclide and reaction without
        propagating the uncertainty to it.

        """
        out = {}
        for (mat, mt), df in self.cov_.groupby(level=(0, 1)):
            s1, cov, s2 = self._get_sharing(mat, mt)
            ## indices are already handled in the x_ variables
            portion = sandwich_unsafe(s1, cov, s2).nom_val.values[0]  # portion must be scalar
            out[mat, mt] = portion
        return pd.Series(out)

    @property
    @drop_uncertainty_calc
    def reaction_wise_rep(self) -> pd.Series:
        """
        Computes the representativity for each nuclide and reaction without
        propagating the uncertainty to it.

        """
        out = {}
        for (mat, mt), df in self.cov_.groupby(level=(0, 1)):
            s1, cov, s2 = self._get_sharing(mat, mt)
            ## indices are already handled in the x_ variables
            portion = representativity(s1, cov, s2, check_idx=False).nom_val.values[0]  # portion must be scalar
            out[mat, mt] = portion
        return pd.Series(out)

    def stdev(self, side=1):
        ## all indices in the original data considered
        s = self.s1.data.S.copy() if side == 1 else self.s2.data.S.copy()
        s.index = s.index.droplevel("Observable")
        idx = self.cov.index.intersection(s.index)
        return stdev(s.loc[idx], self.cov.loc[idx, idx], check_idx=False)

    @drop_uncertainty_calc
    def stdev_no_unc(self, side=1):
        ## all indices in the original data considered
        s = self.s1.data.nom_val.copy() if side == 1 else self.s2.data.nom_val.copy()
        s.index = s.index.droplevel("Observable")
        idx = self.cov.index.intersection(s.index)
        return stdev(s.loc[idx], self.cov.loc[idx, idx], check_idx=False)

    @drop_uncertainty_calc
    def nuclide_apportion(self, denominator: float=None) -> pd.Series:
        """
        Apportions the nuclide-wise sandwich rule to a denominator.

        Parameters
        ----------
        denominator: float, optional
            the denominator the sandwich rule should be apportioned to.
            The default is None, normalizing over the full-set
            representativity.

        """
        # denominator_ must be float even if it is the representativity
        denominator_ = self.stdev(1) * self.stdev(2) if denominator is None else denominator
        return self.nuclide_wise / denominator_

    @drop_uncertainty_calc
    def reaction_apportion(self, denominator: float=None) -> pd.Series:
        """
        Apportions the nuclide-reaction-wise sandwich rule to a denominator.

        Parameters
        ----------
        denominator: float, optional
            the denominator the sandwich rule should be apportioned to.
            The default is None, normalizing over the full-set
            representativity.

        """
        # denominator_ must be float even if it is the representativity
        denominator_ = self.stdev(1) * self.stdev(2) if denominator is None else denominator
        return self.reaction_wise / denominator_
