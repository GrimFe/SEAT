import pandas as pd

import SEAT.Serpent2UncertaintyPropagation.MomentumPropagation.SensitivitySuite as SS
import SEAT.Serpent2UncertaintyPropagation.MomentumPropagation.CovarianceSuite as CS

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

def stdev(s1: pd.DataFrame, cov: pd.DataFrame) -> pd.DataFrame:
    """
    Computes the standard deviation propagated from the covariance
    and a senstitivity vector.

    Parameters:
    s1: pd.DataFrame
        the sensitivity vector
    cov: pd.DataFrame
        the covariancove matrix

    Returns
    -------
    pd.DataFrame of the computed standard deviation

    Note:
    -----
    This is implemented with **0.5 rather than with np.sqrt() as that function
    only takes float arguments and now we work with ufloats.

    """
    return sandwich(s1, cov, s1)**0.5

def representativity(s1: pd.DataFrame, cov: pd.DataFrame, s2: pd.DataFrame) -> pd.DataFrame:
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

    Returns
    -------
    pd.DataFrame of the representativity

    """
    n = sandwich(s1, cov, s2)
    d = stdev(s1, cov) * stdev(s2, cov)
    return n / d


class SCS:
    """
    Container of sensitivity and covariance information.
    """
    def __init__(self, s1: SS.Sensitivity, cov: pd.DataFrame,
                 s2: SS.Sensitivity):
        self.s1 = s1
        self.s2 = s2
        self.cov = cov

    def _get_sharing(self, mat, mt=None):
        query = "N == @mat" if mt is None else "N == @mat and MT == @mt"
        loc = mat if mt is None else (mat, mt)
        s1 = self.s1_.reset_index().query(query).set_index(["N", "MT", "E [MeV]"]).loc[loc]
        s2 = self.s2_.reset_index().query(query).set_index(["N", "MT", "E [MeV]"]).loc[loc]
        cov = self.cov_.loc[loc, loc]
        return s1, cov, s2

    @property
    def idx(self) -> pd.MultiIndex:
        idx1 = self.s1.idx.droplevel("Observable")  # removing the observable information
        idx2 = self.s2.idx.droplevel("Observable")  # removing the observable information
        return self.cov.index.intersection(idx1).intersection(idx2)

    @property
    def s1_(self) -> pd.DataFrame:
        s = self.s1.data.droplevel("Observable")
        return s.loc[self.idx].S

    @property
    def s2_(self) -> pd.DataFrame:
        s = self.s2.data.droplevel("Observable")
        return s.loc[self.idx].S

    @property
    def cov_(self) -> pd.DataFrame:
        return self.cov.loc[self.idx, self.idx]

    @property
    def sandwich(self) -> pd.DataFrame:
        return sandwich(self.s1_, self.cov_, self.s2_)

    @property
    def representativity(self) -> pd.DataFrame:
        return representativity(self.s1_, self.cov_, self.s2_)

    @property
    def reaction_wise(self) -> pd.Series:
        """
        Computes the sandwich rule for each nuclide and reaction.

        """
        out = {}
        for (mat, mt), df in self.cov_.groupby(level=(0, 1)):
            if mat in self.s1_.index.get_level_values(0) and mt in self.s1_.index.get_level_values(1) and\
            mat in self.s2_.index.get_level_values(0) and mt in self.s2_.index.get_level_values(1):  # might be not needed working with x_ variables
                s1, cov, s2 = self._get_sharing(mat, mt)
                portion = sandwich(s1, cov, s2).S.values[0]  # portion must be scalar
            else:
                portion = 0
                if 'total' not in mt:
                    print(f"Skip {mat} - {mt}")
            out[mat, mt] = portion
        return pd.Series(out)

    @property
    def reaction_wise_rep(self) -> pd.Series:
        """
        Computes the representativity for each nuclide and reaction.

        """
        out = {}
        for (mat, mt), df in self.cov_.groupby(level=(0, 1)):
            if mat in self.s1_.index.get_level_values(0) and mt in self.s1_.index.get_level_values(1) and\
            mat in self.s2_.index.get_level_values(0) and mt in self.s2_.index.get_level_values(1):  # might be not needed working with x_ variables
                s1, cov, s2 = self._get_sharing(mat, mt)
                portion = representativity(s1, cov, s2).S.values[0]  # portion must be scalar
            else:
                portion = 0
                if 'total' not in mt:
                    print(f"Skip {mat} - {mt}")
            out[mat, mt] = portion
        return pd.Series(out)

    @property
    def nuclide_wise(self) -> pd.Series:
        """
        Computes the sandwich rule for each nuclide.

        """
        out = {}
        for mat, df in self.cov_.groupby(level=0):
            if mat in self.s1_.index.get_level_values(0) and mat in self.s2_.index.get_level_values(0):  # might be not needed working with x_ variables
                s1, cov, s2 = self._get_sharing(mat)
                portion = sandwich(s1, cov, s2).S.values[0]  # portion must be scalar
            else:
                portion = 0
                print(f"Skip {mat}")                                                
            out[mat] = portion
        return pd.Series(out)

    @property
    def nuclide_wise_rep(self) -> pd.Series:
        """
        Computes the representativity for each nuclide.

        """
        out = {}
        for mat, df in self.cov_.groupby(level=0):
            if mat in self.s1_.index.get_level_values(0) and mat in self.s2_.index.get_level_values(0):  # might be not needed working with x_ variables
                s1, cov, s2 = self._get_sharing(mat)
                portion = representativity(s1, cov, s2).S.values[0]  # portion must be scalar
            else:
                portion = 0
                print(f"Skip {mat}")                                                
            out[mat] = portion
        return pd.Series(out)

    def stdev(self, side=1):
        s = self.s1_.copy() if side == 1 else self.s2_.copy()
        return stdev(s, self.cov_)

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
        denominator_ = self.representativity if denominator is None else denominator
        return self.nuclide_wise / denominator_

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
        denominator_ = self.representativity if denominator is None else denominator
        return self.reaction_wise / denominator_
