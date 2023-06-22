import pandas as pd
import numpy as np

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

    """
    return np.sqrt(sandwich(s1, cov, s1))

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

def u_sandwich(s1: pd.DataFrame, us1: pd.DataFrame, cov: pd.DataFrame,
               s2: pd.DataFrame, us2: pd.DataFrame):
    V1 = pd.concat([us1, us2], axis=0).loc[cov.index] **2
    SS = pd.concat([s1, s2], axis=1)
    V2 = SS.loc[cov.index].T.dot(cov).melt().set_index(SS.index.names) **2
    return np.sqrt(V1.T.dot(V2))

def u_representativity(s1: pd.DataFrame, us1: pd.DataFrame, cov: CS.Covariance,
                       s2: pd.DataFrame, us2: pd.DataFrame):
    A = s2 - (representativity(s1, cov, s2) * np.sqrt(sandwich(s2, cov, s2) / sandwich(s1, cov, s1))).S.values[0] * s1
    B = s1 - (representativity(s1, cov, s2) * np.sqrt(sandwich(s1, cov, s1) / sandwich(s2, cov, s2))).S.values[0] * s2
         
    Sr1 = 1/(stdev(s1, cov) * stdev(s2, cov)).S.values[0] * A.T @ cov
    Sr2 = 1/(stdev(s1, cov) * stdev(s2, cov)).S.values[0] * B.T @ cov

    # Src = C # not implemented yet
    # UC = pd.DataFrame(0, idx_M, idx_M) if UC is None else UC
    return np.sqrt((Sr1 **2).dot(us1 **2) + (Sr2 **2).dot(us2 **2))


class SCS:
    """
    Container of sensitivity and covariance information.
    """
    def __init__(self, s1: SS.Sensitivity, us1: pd.DataFrame,
                 cov: pd.DataFrame, s2:  SS.Sensitivity, us2: pd.DataFrame):
        self.s1 = s1
        self.us1 = us1
        self.s2 = s2
        self.us2 = us2
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
        idx1 = self.s1.idx.droplevel(0)  # removing the observable information
        idx2 = self.s2.idx.droplevel(0)  # removing the observable information
        return self.cov.index.intersection(idx1).intersection(idx2)#.intersection(self.us1.index).intersection(self.us2.index)

    @property
    def s1_(self) -> pd.DataFrame:
        s = self.s1.data.droplevel("Observable")
        return s.loc[self.idx].S

    @property
    def s2_(self) -> pd.DataFrame:
        s = self.s2.data.droplevel("Observable")
        return s.loc[self.idx].S

    @property
    def us1_(self) -> pd.DataFrame:
        s = self.us1.data.droplevel("Observable")
        return s.loc[self.idx]

    @property
    def us2_(self) -> pd.DataFrame:
        s = self.us2.data.droplevel("Observable")
        return s.loc[self.idx]

    @property
    def cov_(self) -> pd.DataFrame:
        return self.cov.loc[self.idx, self.idx]

    @property
    def sandwich(self) -> pd.DataFrame:
        return sandwich(self.s1_, self.cov_, self.s2_)

    @property
    def u_sandwich(self) -> pd.DataFrame:
        """
        Propagates the sensitivity statistical error to the sandwich rule.
        Assumes no correlation.

        """
        return u_sandwich(self.s1_, self.us1_, self.cov_, self.s2_, self.us2_)

    @property
    def representativity(self) -> pd.DataFrame:
        return representativity(self.s1_, self.cov_, self.s2_)

    @property
    def u_representativity(self) -> pd.DataFrame:
        """
        Propagates the sensitivity statistical error to the representativity.
        Assumes no correlation.

        """
        return u_representativity(self.s1_, self.us1_, self.cov_, self.s2_, self.us2_)

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
